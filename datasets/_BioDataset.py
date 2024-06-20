import os
import re
import sys
import bbi
import h5py
import logging
import pathlib
import subprocess

from logging import Logger

import numpy as np
import pandas as pd
from enum import Enum
from pathlib import Path
from abc import abstractmethod
from typing import Any, Callable, Dict, List, Optional, Tuple, Union

from torch.utils.data import Dataset

from mini_utils import enum_value_list, enum_name_list

class BioDataset(Dataset):

    source_list = ["https://hgdownload-test.gi.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeUwRepliSeq/wgEncodeUwRepliSeqBg02esWaveSignalRep1.bigWig",
                   "https://hgdownload-test.gi.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeUwRepliSeq/wgEncodeUwRepliSeqBjWaveSignalRep2.bigWig"]

    def __init__(
        self, 
        raw_path: Union[str, Path], 
        logger: Union[str, Logger] = logging.getLogger(),
        force_download: bool = False,
    ) -> None:
        

        self.logger = logging.getLogger(logger) if isinstance(logger, str) else logger
        self.logger.debug(f"init BioDataset start")

        if not hasattr(self, 'dataset_name') or self.dataset_name is None:
            self.dataset_name = "Bio"
        self.raw_path = Path(raw_path)
        
        self.force_download = force_download

        if force_download or not os.path.isdir(self.raw_path):
            self.download_rawdata()

        self.logger.debug("init BioDataset end.")

    def download_rawdata(self):
        """
        download rawdata listed in self.source_list to the directory self.raw_path
        """
        pathlib.Path.mkdir(self.raw_path, exist_ok=True, parents=True)
        self.logger.info(f"Initialize download directory: {self.raw_path}")
        self.failed_list = []

        for src in self.source_list :
            fname = src.split('/')[-1]
            tgt_f = self.raw_path.joinpath(fname)
            try:
                if not os.path.isfile(tgt_f):
                    self.logger.info(f"Starting download {fname}")
                    self.logger.info(f"{src}")
                    out = subprocess.run(['wget', '-c', '--no-check-certificate', src, '-P', self.raw_path], check=True, shell=False)
                else:
                    self.logger.info(f"Downloaded data file: {fname}")

            except Exception as e:
                self.logger.info(out.stdout)
                self.logger.error(out.stderr)
                self.logger.error(e)
                self.failed_list.append(src)

        if len(self.failed_list) > 0 :
            self.logger.error(f"Failed : \n{self.failed_list}")
        

class BioBigWigDataset(BioDataset):

    """
    BioBigWigDataset manage BigWig file type datasets. It is a compressed table looks like the following example :
    > import bbi
    > fd.fetch_intervals('chr1', 1, 100000)
    > bigwig_fd = bbi.open(bigwig_path)
    > bigwig_fd.fetch_intervals('chr1', 1, 10000)
        chrom	start	end	value
    0	chr1	1	10117	0.15904
    1	chr1	10117	10131	0.48031
    2	chr1	10131	10164	0.39292
    3	chr1	10164	10166	0.25747
    ....

    BioBigWigDataset building steps :
    1. download raw data source
    2. build h5 file for each sample file
    3. in each h5 sample file, each chromosome is present as group which contains different resolutions

    So that the h5 path layout should be like this : 
    \{h5_path}
    |- datafile1.h5
      \- group_chr1
       |- resolution_100_overlap_0
       |- resolution_1000_overlap_0
       |- resolution_10000_overlap_0
       ....
      \- group_chr2
       |- resolution_100_overlap_0
       |- resolution_1000_overlap_0
       |- resolution_10000_overlap_0
       ....
      \- group_chr3
       |- ...
      ....
    |- datafile2.h5
    |- datafile3.h5
    .....

    for more details about BigWig format, please refer to : https://genome.ucsc.edu/goldenPath/help/bigWig.html
    """

    class BigWigChm(Enum):
        chr1 = 1
        chr2 = 2
        chr3 = 3
        chr4 = 4
        chr5 = 5
        chr6 = 6
        chr7 = 7
        chr8 = 8
        chr9 = 9
        chr10 = 10
        chr11 = 11
        chr12 = 12
        chr13 = 13
        chr14 = 14
        chr15 = 15
        chr16 = 16
        chr17 = 17
        chr18 = 18
        chr19 = 19
        chr20 = 20
        chr21 = 21
        chr22 = 22
        chrX  = 23
        # chrY  = 24
        # chrM  = 25

    BigWigChromSizes = {
        BigWigChm(1): 249250621,
        BigWigChm(2): 243199373,
        BigWigChm(3): 198022430,
        BigWigChm(4): 191154276,
        BigWigChm(5): 180915260,
        BigWigChm(6): 171115067,
        BigWigChm(7): 159138663,
        BigWigChm(8): 146364022,
        BigWigChm(9): 141213431,
        BigWigChm(10): 135534747,
        BigWigChm(11): 135006516,
        BigWigChm(12): 133851895,
        BigWigChm(13): 115169878,
        BigWigChm(14): 107349540,
        BigWigChm(15): 102531392,
        BigWigChm(16): 90354753,
        BigWigChm(17): 81195210,
        BigWigChm(18): 78077248,
        BigWigChm(19): 59128983,
        BigWigChm(20): 63025520,
        BigWigChm(21): 48129895,
        BigWigChm(22): 51304566,
        BigWigChm(23): 155270560,
        # BigWigChm(Y): 59373566
        # BigWigChm(M): 16571,
    }

    class H5Attrs(Enum):
        COLUMNS   = 'columns'
        INDEX     = 'index'

    class BigWigSummary(Enum):
        mean = 'mean'
        max  = 'max'
        min  = 'min'
        std  = 'std'
        cov  = 'cov' # coverage

        chrom= 0
        start= 1
        end  = 2
        summary= 3
    
    def __init__(
        self, 
        h5_path: Union[str, Path], 
        raw_path: Union[str, Path], 
        resolutions: list[int],
        overlap : int = 0,
        h5_chunk_size: int = 100,
        summary : Union[str, List[str]] = 'mean',
        logger: Union[str, Logger] = logging.getLogger(),
        force_download:bool = False,
        rebuild_h5:bool = False,
        preprocess: Optional[Callable] = None, 
        transform:  Optional[Callable] = None, 
        lazy_load: bool = True
        ) -> None:

        if not hasattr(self, 'dataset_name') or self.dataset_name is None:
            self.dataset_name = "BioBigWig"
    
        super().__init__(raw_path = raw_path, logger = logger, force_download = force_download)

        self.logger.debug("init BioBigWigDataset start")
        self.resolutions = resolutions
        self.overlap = overlap
        self.h5_chunk_size= h5_chunk_size

        # will raise exception if summary function name is not in enum
        self.summary = [ self.BigWigSummary(s) for s in np.array([summary]).reshape(-1).tolist() ]

        self.h5_path = Path(h5_path)
        self.h5_list = []
        self.rebuild_h5 = rebuild_h5

        self.preprocess = preprocess
        self.transform  = transform

        list.sort(self.resolutions)
        chromsizes = np.array([self.BigWigChromSizes[chr] for chr in self.BigWigChromSizes.keys()])

        self.sample_nums = np.ceil(chromsizes/(self.resolutions[0] - self.overlap)/self.h5_chunk_size)
        self.sample_cum_nums = np.cumsum(self.sample_nums)

        for bigwig_src in self.source_list:
            bigwig_fname = Path(bigwig_src).name
            bigwig_fname = self.raw_path.joinpath(bigwig_fname)
            h5_fname     = self._h5_fname(bigwig_fname.name)
            self.h5_list.append(h5_fname)
            self.build_h5(bigwig = bigwig_fname, 
                          h5 = h5_fname, 
                          resolutions = self.resolutions, 
                          summary=self.summary) 
            
        self.summary_h5_fname = self.h5_path.joinpath(f"{self.dataset_name}.h5")
        self.build_h5_summary()

        if self.preprocess is not None:
            self.preprocess(self.summary_h5_fname)

        if os.path.isfile(self.summary_h5_fname):
            self.summary_h5_fd = h5py.File(self.summary_h5_fname, 'r')
        else:
            self.summary_h5_fd = None

        self.logger.debug("init BioBigWigDataset end.")

    def _h5_fname(self, bigwig_fname: str) -> Path:
        # replace bigwig by ignoring cases
        h5_fname = re.compile('bigwig', re.IGNORECASE).sub('h5', bigwig_fname)
        return self.h5_path.joinpath(h5_fname)

    def _best_cover(self, start_position: int, end_position: int, chm_size: int, resolution: int):
        """
        return the best centralized chunked region which cover the given interval
        """
        mid = np.mean([start_position, end_position])
        start = mid - resolution*self.h5_chunk_size /2 
        start = max(0,start)
        end   = start + resolution
        end_position   = int(np.rint(min(chm_size, end)))
        start_position = end_position - resolution
        return start_position, end_position

    def _build_position_encoding(self, chm: int, start_position: int, end_position: int, rslt: int, df_len: int):
        # add position encoding (chrom, start, end)
        assert (start_position-end_position)/rslt == df_len
        position_df = pd.DataFrame(np.zeros((df_len, 3)), columns=[self.BigWigSummary(i).name for i in range(3)])
        position_df[self.BigWigSummary.chrom.name] = chm
        position_df[self.BigWigSummary.start.name] = np.arange(start_position, end_position, rslt)
        position_df[self.BigWigSummary.end.name]   = np.arange(start_position + rslt, end_position + rslt, rslt)

    def _bigwig2df(self, 
                   bigwig: Path, 
                   chr: BigWigChm, 
                   resolution: int, 
                   summary: List[BigWigSummary]
        ) -> pd.DataFrame:

        self.logger.info(f"Open BigWig file: {bigwig}")
        with bbi.open(str(bigwig)) as epig_fd :
            # chr_size = epig_fd.chromsizes[chr]
            chr_size = self.BigWigChromSizes[chr]
            self.logger.info(f"{chr.name}, length {chr_size}")
            starts = np.arange(0, chr_size, resolution - self.overlap)
            ends   = starts + resolution
            chrs   = np.array([chr.name] * len(starts))
            self.logger.info(f"summary functions: {[ s.value for s in summary ]}")
            values = [ epig_fd.stackup(chrs, starts, ends, bins=1, summary = s.value).reshape(-1).tolist() for s in summary ]
            # bigwig always provides the same chromosome sizes, so within same resolution, 
            # to count from position 0 will always gives the same index
            # don't have to always save start and end position
            epig_df= pd.DataFrame(np.array(values).T, columns=[s.value for s in summary])

        return epig_df

    def _dataset_name(self, rslt: int, overlap: int) -> str:
        return f"{rslt}_{overlap}"

    def _dataset_fullname(self, chr: str, rslt: int, overlap: int) -> str:
        dataset_name = self._dataset_name(rslt=rslt, overlap=overlap)
        return f"{chr}/{dataset_name}"

    def build_h5(self, 
                 bigwig: Path, 
                 h5: Path, 
                 resolutions: list[int],
                 summary: List[BigWigSummary]
        ) -> None:

        mode = 'a'
        if not os.path.isfile(h5) or self.rebuild_h5:
            pathlib.Path.mkdir(h5.parent, exist_ok=True, parents=True)
            mode = 'w'

        with h5py.File(h5, mode) as h5fd:
            self.logger.debug(f"open h5 file {h5}")
            for rslt in resolutions:
                dataset_name = self._dataset_name(rslt=rslt, overlap=self.overlap)
                for chr in self.BigWigChm:
                    self.logger.info(f"Building chromosome {chr.name} in resolution {rslt}. ")
                    if self.rebuild_h5 or chr.name not in h5fd.keys() or dataset_name not in h5fd[chr.name].keys() :
                        dataset_fullname = self._dataset_fullname(chr=chr.name, rslt=rslt, overlap=self.overlap)
                        self.logger.info(f"create dataset {dataset_fullname} in the h5 file")
                        data_df = self._bigwig2df(bigwig, chr, rslt, summary)
                        h5fd.create_dataset(name = dataset_fullname, data = data_df)
                        h5fd[dataset_fullname].attrs[self.H5Attrs.COLUMNS.value] = data_df.columns.to_list()

    def _concat_summary_table(self, 
                              tgt_h5fd: h5py.File, 
                              src_h5fd_dict: Dict[str, h5py.File], 
                              chr: BigWigChm, 
                              rslt: int, 
                              overlap: int
        ) -> h5py.File :

        L = -1
        columns_count = 0

        self.logger.debug(f"{src_h5fd_dict}")

        for k in src_h5fd_dict:
            # check if the summary tables have same length
            self.logger.debug(self._dataset_fullname(chr.name, rslt, overlap))
            ds = src_h5fd_dict[k][self._dataset_fullname(chr.name, rslt, overlap)]
            self.logger.debug(f"shape of dataset: {ds}")
            if L < 0 :
                assert ds.shape[0] > 0
                L = ds.shape[0]
            else :
                assert L == ds.shape[0]
            columns_count += ds.shape[1]

        self.logger.debug(f"estimate the shape of entire dataset : {(L, columns_count)}")
        self.logger.debug(f"set the chunk size : {(self.h5_chunk_size, columns_count)}")
        # create chunked dataset
        tgt_h5fd.create_dataset(name = self._dataset_fullname(chr.name, rslt, overlap), 
                                shape = (L, columns_count), 
                                chunks= (self.h5_chunk_size, columns_count), 
                                dtype = float)

        # update to tgt_h5fd dataset one by one, since each one can be very large
        column_names = []
        columns_idx  = 0
        for k, fd in src_h5fd_dict.items():
            ds = fd[self._dataset_fullname(chr.name, rslt, overlap)]
            attrs = ds.attrs[self.H5Attrs.COLUMNS.value]
            self.logger.debug(f"attributes of dataset: {attrs}")
            column_names += [ f"{k}_{s}" for s in attrs ]
            self.logger.debug(f"column names: {column_names}")
            tgt_h5fd[self._dataset_fullname(chr.name,rslt,overlap)][:,columns_idx: columns_idx + ds.shape[1]] = ds[:]
            columns_idx += ds.shape[1]

        tgt_h5fd.attrs[self.H5Attrs.COLUMNS.value] = column_names
        return tgt_h5fd
    
    def build_h5_summary(self):
        print("You should implement the function: build_h5_summary")
    
    def __getitem__(self, index) -> Any:
        """
        index iterate over chunked dataframe
        """
        # find the index for the minimum resolution
        chm_idx = sum(self.sample_cum_nums < index)
        chm = self.BigWigChm(chm_idx+1)
        idx = index - self.sample_cum_nums[chm]
        start_position = idx * min(self.resolutions) * self.h5_chunk_size
        end_position   = (idx+1) * min(self.resolutions) * self.h5_chunk_size

        summary_list = []

        for rslt in self.resolutions:
            ds = self.summary_h5_fd[self._dataset_fullname(chm,rslt,self.overlap)]
            coll  = len(ds.attrs[self.H5Attrs.COLUMNS.value])
            start_position, end_position = self._best_cover(start_position, 
                                                            end_position,
                                                            self.BigWigChromSizes[chm.name],
                                                            rslt)
            
            values_df = pd.DataFrame(ds[(slice(int(start_position/rslt), int(end_position/rslt), 1),slice(0,coll,1))], 
                                     columns = ds.attrs[self.H5Attrs.COLUMNS.value])
            position_df = self._build_position_encoding(chm.value, start_position, end_position, rslt, values_df.shape[0])
            summary_df = pd.concat([position_df, values_df], axis=1)
            summary_list.append(summary_df)

        # 2. concat all the resolutions
        df = pd.concat(summary_list, axis=0)

        if self.transform is not None:
            df = self.transform(df)
        
        return df
    
    def __len__(self):
        # regarding to the minimum resolution
        return np.sum(self.sample_nums)

class BioMafDataset(BioDataset):
    def __init__(self, 
                 raw_path: str | Path, 
                 logger=logging.getLogger(os.getcwd()), 
                 force_download=False
        ) -> None:

        self.dataset_name = "BioMaf"
    
        super().__init__(raw_path = raw_path, logger = logger, force_download = force_download)

        # TODO: build h5 for maf files

class BioBigWigDatasetFolder(Dataset):
    def __init__(
        self, 
        root_path: Union[str, Path],
        resolution: Union[int, List[int]],
        preprocess: Callable, 
        dataset_cls: BioDataset,
        transform:  Optional[Callable] = None
    ) -> None:
        
        self.root_path = Path(root_path)

        if type(resolution) == int:
            self.resolution = [resolution]
        else:
            self.resolution = resolution

        self.datasets = [
            dataset_cls(
                raw_path = self.root_path.join(fn),
                resolution= r, 
                preprocess= preprocess, 
                transform = transform
            ) for fn in os.listdir(self) 
                for r in self.resolution 
                    if os.path.isfile(self.root_path.joinpath(fn))
        ]
