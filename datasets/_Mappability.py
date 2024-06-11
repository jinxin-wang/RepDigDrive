import os
import logging

from pathlib import Path
from _BioDataset import BioBigWigDataset

from enum import Enum
from typing import Any, Callable, Dict, List, Optional, Tuple, Union

class MappabilityDataset(BioBigWigDataset):

    addr_list =["https://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeMapability/wgEncodeCrgMapabilityAlign24mer.bigWig",
                "https://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeMapability/wgEncodeCrgMapabilityAlign36mer.bigWig",
                "https://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeMapability/wgEncodeCrgMapabilityAlign40mer.bigWig",
                "https://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeMapability/wgEncodeCrgMapabilityAlign50mer.bigWig",
                "https://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeMapability/wgEncodeCrgMapabilityAlign75mer.bigWig",
                "https://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeMapability/wgEncodeCrgMapabilityAlign100mer.bigWig"]
    
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
        resolutions: list[int],
        overlap : int = 0,
        logger = logging.getLogger(os.getcwd()),
        force_download = False,
        rebuild_h5 = False,
    ) -> None:
        
        self.dataset_name = "Roadmap Epigenomics"
    
        super().__init__(h5_path = h5_path, 
                         raw_path=raw_path, 
                         resolutions=resolutions,
                         overlap=overlap,
                         logger=logger,
                         force_download=force_download,
                         rebuild_h5=rebuild_h5)

    
