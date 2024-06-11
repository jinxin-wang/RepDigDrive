import os
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
    2. save the table to h5 file for each chromosome
    3. in each chromosome, each datafile is present in different resolutions

    So that the h5 path layout should be like this : 
    \{h5_path}
    |- chr1
      \- datafile1.h5
      |- datafile2.h5
      |- datafile3.h5
      ....
    |- chr2
    |- chr3
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
        chrY  = 24
        chrM  = 25

    def __init__(
        self, 
        h5_path: Union[str, Path], 
        raw_path: Union[str, Path], 
        resolutions: list[int],
        overlap : int = 0,
        logger = logging.getLogger(os.getcwd()),
        force_download = False,
        rebuild_h5 = False,
    ) -> None:
        
        self.dataset_name = "BioBigWig"
    
        super().__init__(raw_path = raw_path, logger = logger, force_download = force_download)

        self.resolutions = resolutions
        self.overlap = overlap
        self.h5_path = Path(h5_path)
        self.rebuild_h5 = rebuild_h5

        for bigwig_fname in os.listdir(self.raw_path):
            bigwig_fname = self.raw_path.joinpath(bigwig_fname)

            # don't +1 to the len, ignore chrM
            for i in range(1,len(self.BigWigChm)):
                chr = self.BigWigChm(i)
                h5_fname = self.h5_path.joinpath(chr.name).joinpath(bigwig_fname.name.replace('bigwig','h5'))
                if not os.path.isfile(h5_fname) or self.rebuild_h5:
                    self.build_h5(chr = chr,
                                  bigwig = bigwig_fname, 
                                  h5  = h5_fname,
                                  resolutions = self.resolutions)

    # def __getitem__(self, index) -> Any:
    #     return super().__getitem__(index)
    
    # def __len__(self):
    #     pass
    
    def _bigwig2df(self, 
                   bigwig: Path, 
                   chr: BigWigChm, 
                   resolution: int, 
                   summary: str = 'mean'
        ) -> pd.DataFrame:

        with bbi.open(bigwig) as epig_fd :
            chr_size = epig_fd.chromsizes[chr]
            self.logger.info(f"{chr.name}, length {chr_size}")
            if self.overlap <= 0:
                bin = np.ceil(chr_size/resolution).astype(int)
                self.logger.info(f"fetch table in {resolution} BP resolution. The number of bin is {bin}")
                epig_df = epig_fd.fetch(chr.name, 0, chr_size, bin=bin, summary = summary)
            else:
                starts = np.arange(0,chr_size, self.overlap)
                ends   = starts + resolution
                chrs   = np.array([chr.name] * len(starts))
                values = epig_fd.stackup(chrs, starts, ends, bins=1, summary = summary)
                values = np.reshape(values, -1)
                epig_df= pd.DataFrame([chrs, starts, ends, values], columns=['chrom','start','end', 'value'])
            return epig_df

    def build_h5(self, 
                 chr: BigWigChm, 
                 bigwig: Path, 
                 h5: Path, 
                 resolutions: list[int],
                 summary: str = 'mean') -> None:
        
        self.logger.info(f"Building chromosome {chr.name} of {bigwig} to {h5} ")
        
        mode = 'a'
        if not os.path.isfile(h5) or self.rebuild_h5:
            mode = 'w'

        with h5py.File(self.h5_path, mode) as h5fd:
            # h5fd.create
            for rslt in resolutions:
                dataset_name = f"{rslt}_{summary}_{self.overlap}"

                h5fd.create_dataset(name = dataset_name, 
                                    data = self._bigwig2df(bigwig, chr, rslt, summary))


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

    def __getitem__(self, index) -> Any:
        return super().__getitem__(index)
    
    def __len__(self):
        pass
