import os
import logging
import h5py

from pathlib import Path
from _BioDataset import BioBigWigDataset

from enum import Enum
from typing import Any, Callable, Dict, List, Optional, Tuple, Union

class MappabilityDataset(BioBigWigDataset):

    mirror = "https://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeMapability"
    
    class WINSIZE(Enum):
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
        
        self.h5_chunk_size= h5_chunk_size
        self.dataset_name = "Mappability"
        self.design_mers = design_mers
        self.preprocess = preprocess
        self.transform  = transform
        self.source_list  = [ f"{self.mirror}/{self._bigwig_fname(self.WINSIZE(alg).name)}" for alg in self.design_mers ] 

        super().__init__(h5_path = h5_path, 
                         raw_path=raw_path, 
                         resolutions=resolutions,
                         overlap=overlap,
                         logger=logger,
                         force_download=force_download,
                         rebuild_h5=rebuild_h5)
        
        self.mapp_h5_fname = self.h5_path.joinpath(f"{self.dataset_name}.h5")

        mode = 'a'
        if rebuild_h5:
            mode = 'w'

        with h5py.File(self.mapp_h5_fname, mode=mode) as h5fd:
            for rslt in self.resolutions:
                dataset_name = f"{rslt}_{self.overlap}"
                for chr in self.BigWigChm:
                    if self.rebuild_h5 or chr.name not in h5fd.keys() or dataset_name not in h5fd[chr.name].keys() :
                        dataset_fullname = f"{chr.name}/{dataset_name}"
                        self._concat_summary_table(h5fd, self.design_mers, chr, rslt, overlap)

    def _bigwig_fname(self, key):
        return f"wgEncodeCrgMapability{key}.bigWig"

    def _h5_fname_to_mer(self, h5_fname: Path):
        """
        given filename, return mer enum element
        """
        return self.WINSIZE(int(h5_fname.name.split('mer')[0].replace('wgEncodeCrgMapabilityAlign','')))

    def _concat_summary_table(self, h5fd: h5py.File, design_mers: List[int], chr: str, rslt: int, overlap: int):
        # open h5 files and check summary tables 
        L = -1
        h5fd_dict = { self._h5_fname_to_mer(h5): h5py.File(h5, 'r') for h5 in self.h5_list }
        for mer,fd in h5fd_dict.items():
            # check if the summary tables have same length
            ds = fd[self._dataset_fullname(chr, rslt, overlap)]
            if L < 0 :
                L = ds.shape[0]
                start = ds[:,0]

            else :
                len_summary_table = ds.shape[0]
                assert L == len_summary_table

                num_same_start = sum(start == ds[:,self.BigWigSummary.start.value])
                assert L == num_same_start

                columns_count += (ds.shape[0]-2) # count without columns of start and end

                # TODO: remove start and end if they are always the same, keep only summary
                







# addr_list =["https://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeMapability/wgEncodeCrgMapabilityAlign24mer.bigWig",
#             "https://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeMapability/wgEncodeCrgMapabilityAlign36mer.bigWig",
#             "https://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeMapability/wgEncodeCrgMapabilityAlign40mer.bigWig",
#             "https://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeMapability/wgEncodeCrgMapabilityAlign50mer.bigWig",
#             "https://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeMapability/wgEncodeCrgMapabilityAlign75mer.bigWig",
#             "https://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeMapability/wgEncodeCrgMapabilityAlign100mer.bigWig"]