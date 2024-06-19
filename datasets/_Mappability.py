import os
import logging
import h5py

import numpy as np
import pandas as pd

from pathlib import Path
from _BioDataset import BioBigWigDataset

from enum import Enum
from typing import Any, Callable, Dict, List, Optional, Tuple, Union

class MappabilityDataset(BioBigWigDataset):

    mirror = "https://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeMapability"
    
    class MER(Enum):
        Align24mer = 24
        Align36mer = 36
        Align40mer = 40
        Align50mer = 50
        Align75mer = 75
        Align100mer= 100

    def __init__(
        self, 
        h5_path: Union[str, Path], 
        raw_path: Union[str, Path], 
        resolutions: list[int] = [10000, 100000],
        overlap : int = 0,
        h5_chunk_size: int = 100,
        logger = logging.getLogger(os.getcwd()),
        force_download = False,
        rebuild_h5 = False,
        design_mers: List[int] = [24, 36, 40, 50, 75, 100],
        preprocess: Optional[Callable] = None, 
        transform:  Optional[Callable] = None
    ) -> None:
        
        self.dataset_name = "Mappability"
        self.design_mers = design_mers
        self.preprocess = preprocess
        self.transform  = transform
        self.source_list  = [ f"{self.mirror}/{self._bigwig_fname(self.MER(num).name)}" for num in self.design_mers ] 

        list.sort(resolutions)

        super().__init__(h5_path = h5_path, 
                         raw_path=raw_path, 
                         resolutions = resolutions,
                         overlap = overlap,
                         h5_chunk_size = h5_chunk_size,
                         logger = logger,
                         force_download = force_download,
                         rebuild_h5 = rebuild_h5)
        
        self.mapp_h5_fname = self.h5_path.joinpath(f"{self.dataset_name}.h5")

        chromsizes = np.array([self.BigWigChromSizes[chr] for chr in self.BigWigChromSizes.keys()])

        self.sample_lens = np.ceil(chromsizes/(self.resolutions[0] - self.overlap)/self.h5_chunk_size)
        self.sample_cum_lens = np.cumsum(self.sample_lens)

        self.build_h5_summary()

        if self.preprocess is not None:
            self.mapp_h5_fname = self.preprocess(self.mapp_h5_fname)

        self.mapp_summary_h5_fd = h5py.File(self.mapp_h5_fname, 'r')

    def build_h5_summary(self):
        mode = 'a'
        if self.rebuild_h5:
            mode = 'w'

        with h5py.File(self.mapp_h5_fname, mode=mode) as h5fd:
            for rslt in self.resolutions:
                dataset_name = self._dataset_name(rslt, self.overlap) 
                for chr in self.BigWigChm:
                    if self.rebuild_h5 or chr.name not in h5fd.keys() or dataset_name not in h5fd[chr.name].keys() :
                        # open designed h5 files 
                        h5fd_dict = { self._h5_fname_to_mer(h5): h5py.File(h5, 'r') for h5 in self.h5_list if self._h5_fname_to_mer(h5).value in self.design_mers }
                        self._concat_summary_table(h5fd, h5fd_dict, chr, rslt, self.overlap)
                        for k in h5fd_dict:
                            h5fd_dict[k].close()

    def _bigwig_fname(self, key):
        return f"wgEncodeCrgMapability{key}.bigWig"

    def _h5_fname_to_mer(self, h5_fname: Path):
        """
        given filename, return mer enum element
        """
        return self.MER(int(h5_fname.name.split('mer')[0].replace('wgEncodeCrgMapabilityAlign','')))

    def _concat_summary_table(self, 
                              tgt_h5fd: h5py.File, 
                              src_h5fd_dict: dict[MER, h5py.File], 
                              chr: str, 
                              rslt: int, 
                              overlap: int
        ) -> h5py.File :

        L = -1

        for mer in src_h5fd_dict:
            # check if the summary tables have same length
            ds = src_h5fd_dict[mer][self._dataset_fullname(chr, rslt, overlap)]
            if L < 0 :
                assert ds.shape[0] > 0
                L = ds.shape[0]

            else :
                assert L == ds.shape[0]
                columns_count += ds.shape[1]
        
        # create chunked dataset
        tgt_h5fd.create_dataset(name=self._dataset_fullname(chr, rslt, overlap), 
                            shape=(L, columns_count), 
                            chunks=(self.h5_chunk_size, columns_count), 
                            dtype=float)

        # update to tgt_h5fd dataset one by one, since each one can be very large
        column_names = []
        columns_idx  = 0
        for mer, fd in src_h5fd_dict.items():
            ds = fd[self._dataset_fullname(chr, rslt, overlap)]
            attrs = fd.attrs[self.H5Attrs.COLUMNS.value]
            column_names += [ f"{mer}_{s}" for s in attrs ]
            tgt_h5fd[self._dataset_fullname(chr,rslt,overlap)][:,columns_idx: columns_idx + ds.shape[1]] = ds[:]
            columns_idx += ds.shape[1]

        tgt_h5fd.attrs[self.H5Attrs.COLUMNS.value] = column_names
        return tgt_h5fd
                
    def __getitem__(self, index) -> Any:
        """
        index iterate over chunked dataframe
        """
        # find the index for the minimum resolution
        chm_idx = sum(self.sample_cum_lens < index)
        chm = self.BigWigChm(chm_idx+1)
        idx = index - self.sample_cum_lens[chm]
        start_position = idx * min(self.resolutions) * self.h5_chunk_size
        end_position   = (idx+1) * min(self.resolutions) * self.h5_chunk_size

        summary_list = []

        for rslt in self.resolutions:
            ds = self.mapp_summary_h5_fd[self._dataset_fullname(chm,rslt,self.overlap)]
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
        return np.sum(self.sample_lens)


# addr_list =["https://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeMapability/wgEncodeCrgMapabilityAlign24mer.bigWig",
#             "https://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeMapability/wgEncodeCrgMapabilityAlign36mer.bigWig",
#             "https://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeMapability/wgEncodeCrgMapabilityAlign40mer.bigWig",
#             "https://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeMapability/wgEncodeCrgMapabilityAlign50mer.bigWig",
#             "https://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeMapability/wgEncodeCrgMapabilityAlign75mer.bigWig",
#             "https://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeMapability/wgEncodeCrgMapabilityAlign100mer.bigWig"]