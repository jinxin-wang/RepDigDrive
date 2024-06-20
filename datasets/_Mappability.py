import os
import logging
import h5py

from pathlib import Path
# from _BioDataset 
from datasets import BioBigWigDataset

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
        self.design_mers = [ self.MER(s) for s in design_mers ]
        self.source_list  = [ f"{self.mirror}/{self._bigwig_fname(mer.name)}" for mer in self.design_mers ] 

        super().__init__(h5_path = h5_path, 
                         raw_path=raw_path, 
                         resolutions = resolutions,
                         overlap = overlap,
                         h5_chunk_size = h5_chunk_size,
                         logger = logger,
                         force_download = force_download,
                         rebuild_h5 = rebuild_h5,
                         preprocess = preprocess,
                         transform  = transform)

    def build_h5_summary(self):
        mode = 'a'
        if self.rebuild_h5:
            mode = 'w'

        # open designed h5 files 
        h5fd_dict = { self._h5_fname_to_mer(h5).name: h5py.File(h5, 'r') for h5 in self.h5_list if self._h5_fname_to_mer(h5) in self.design_mers }

        self.logger.info(f"start building summary: {self.summary_h5_fname}")
        with h5py.File(self.summary_h5_fname, mode=mode) as h5fd:
            for rslt in self.resolutions:
                dataset_name = self._dataset_name(rslt, self.overlap) 
                for chr in self.BigWigChm:
                    if self.rebuild_h5 or chr.name not in h5fd.keys() or dataset_name not in h5fd[chr.name].keys() :
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



# addr_list =["https://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeMapability/wgEncodeCrgMapabilityAlign24mer.bigWig",
#             "https://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeMapability/wgEncodeCrgMapabilityAlign36mer.bigWig",
#             "https://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeMapability/wgEncodeCrgMapabilityAlign40mer.bigWig",
#             "https://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeMapability/wgEncodeCrgMapabilityAlign50mer.bigWig",
#             "https://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeMapability/wgEncodeCrgMapabilityAlign75mer.bigWig",
#             "https://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeMapability/wgEncodeCrgMapabilityAlign100mer.bigWig"]