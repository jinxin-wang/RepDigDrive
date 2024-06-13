import os
import re
import sys
import bbi
import h5py
import logging
import pathlib
import subprocess

import numpy as np
import pandas as pd
from enum import Enum
from pathlib import Path
from typing import Any, Callable, Dict, List, Optional, Tuple, Union


from torch.utils.data import Dataset

class BioDataset(Dataset):

    source_list = []

    def __init__(
        self, 
        raw_path: Union[str, Path], 
        logger = logging.getLogger(os.getcwd()),
        force_download = False,
    ) -> None:
        self.dataset_name = "Bio"
        self.raw_path = Path(raw_path)
        self.logger = logger 
        self.force_download = force_download

        if force_download or not os.path.isdir(self.raw_path):
            self.download_rawdata()

    def download_rawdata(self):
        """
        download rawdata listed in self.source_list to the directory self.raw_path
        """
        pathlib.Path.mkdir(self.raw_path, exist_ok=True, parents=True)
        self.failed_list = []

        for src in self.source_list :
            fname = src.split('/')[-1]
            tgt_f = self.raw_path.joinpath(fname)
            self.logger.info(f"download {fname} to {self.raw_path}")
            self.logger.info(f"{src}")
            try:
                if os.path.isfile(tgt_f):
                    out = subprocess.run(['wget', '-c', '--no-check-certificate', src, '-P', self.raw_path], check=True, shell=False)

            except Exception as e:
                self.logger.info(out.stdout)
                self.logger.error(out.stderr)
                self.logger.error(e)
                self.failed_list.append(src)

        if len(self.failed_list) > 0 :
            self.logger.error(f"Failed : {"\n".join(self.failed_list)}")
        

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

    source_list = []

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
        'chr1': 249250621,        
        'chr2': 243199373,
        'chr3': 198022430,
        'chr4': 191154276,
        'chr5': 180915260,
        'chr6': 171115067,
        'chr7': 159138663,
        'chr8': 146364022,
        'chr9': 141213431,
        'chr10': 135534747,
        'chr11': 135006516,
        'chr12': 133851895,
        'chr13': 115169878,
        'chr14': 107349540,
        'chr15': 102531392,
        'chr16': 90354753,
        'chr17': 81195210,
        'chr18': 78077248,
        'chr19': 59128983,
        'chr20': 63025520,
        'chr21': 48129895,
        'chr22': 51304566,
        # 'chrM': 16571,
        'chrX': 155270560,
        # 'chrY': 59373566
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
        logger = logging.getLogger(os.getcwd()),
        force_download = False,
        rebuild_h5 = False,
        ) -> None:
        
        self.dataset_name = "BioBigWig"
    
        super().__init__(raw_path = raw_path, logger = logger, force_download = force_download)

        self.resolutions = resolutions
        self.overlap = overlap
        self.h5_chunk_size= h5_chunk_size

        # will raise exception if summary function name is not in enum
        self.summary = [ self.BigWigSummary(s).value for s in np.array([summary]).reshape(-1).tolist() ]

        self.h5_path = Path(h5_path)
        self.rebuild_h5 = rebuild_h5
        self.h5_list = []

        for bigwig_fname in os.listdir(self.raw_path):
            bigwig_fname = self.raw_path.joinpath(bigwig_fname)
            h5_fname     = self._h5_fname(bigwig_fname.name)
            self.h5_list.append(h5_fname)
            self.build_h5(bigwig = bigwig_fname, 
                          h5 = h5_fname, 
                          resolutions = self.resolutions) 

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
                   summary: List[str]
        ) -> pd.DataFrame:

        with bbi.open(bigwig) as epig_fd :
            # chr_size = epig_fd.chromsizes[chr]
            chr_size = self.BigWigChromSizes[chr]
            self.logger.info(f"{chr.name}, length {chr_size}")
            starts = np.arange(0, chr_size, resolution - self.overlap)
            ends   = starts + resolution
            chrs   = np.array([chr.name] * len(starts))
            values = [ epig_fd.stackup(chrs, starts, ends, bins=1, summary = s).reshape(-1).tolist() for s in summary ]

            # epig_df= pd.DataFrame(np.array([chrs, starts, ends, *values]).T, columns=['chrom','start','end', *summary])
            # epig_df= pd.DataFrame(np.array([starts, ends, *values]).T, columns=[self.BigWigSummary.start.name,self.BigWigSummary.end.name, *summary])

            # bigwig always provides the same chromosome sizes, so within same resolution, 
            # to count from position 0 will always gives the same index
            # don't have to always save start and end position
            epig_df= pd.DataFrame(np.array(values).T, columns=summary)

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
                 summary: List[str]
        ) -> None:
        
        self.logger.info(f"Building chromosome {chr.name} of {bigwig} to {h5} ")
        
        mode = 'a'
        if not os.path.isfile(h5) or self.rebuild_h5:
            mode = 'w'

        with h5py.File(self.h5, mode) as h5fd:
            for rslt in resolutions:
                dataset_name = self._dataset_name(rslt=rslt, overlap=self.overlap)
                for chr in self.BigWigChm:
                    if self.rebuild_h5 or chr.name not in h5fd.keys() or dataset_name not in h5fd[chr.name].keys() :
                        dataset_fullname = self._dataset_fullname(chr=chr.name, rslt=rslt, overlap=self.overlap)
                        data_df = self._bigwig2df(bigwig, chr, rslt, summary)
                        h5fd.create_dataset(name = dataset_fullname, data = data_df)
                        h5fd[dataset_fullname].attrs[self.H5Attrs.COLUMNS.value] = data_df.columns


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
