import logging
from logging import Logger
from pathlib import Path
import numpy as np

from typing import Any, Callable, Dict, List, Optional, Tuple

from ._BioDataset import BioDigDriverfDataset


class Dietlein(BioDigDriverfDataset):

    """
    designed_sets = 'Kidney', 'Esogastric', 'LungNSC', 
                    'Pancreas', 'Breast', 'Liver', 'Brain',
                    'HeadNeck', 'Bladder', 'Prostate'
    """

    mirror = "https://cb.csail.mit.edu/cb/DIG/downloads/mutation_files/megacohorts"

    subset_list =['Kidney', 'Esogastric', 'LungNSC', 'Pancreas', 'Breast', 'Liver', 'Brain',
                  'HeadNeck', 'Bladder', 'Prostate']

    def __init__(self, 
                 h5_path: str | Path, 
                 raw_path: str | Path, 
                 designed_subsets: List[str] | str = 'all',
                 resolutions: List[int] = [10000, 100000], 
                 overlap: int = 0, 
                 h5_chunk_size: int = 100, 
                 N_grams: List[int] | int = 3, 
                 logger: str | Logger = logging.getLogger(),
                 force_download: bool = False, 
                 concurrent_download: int = 0, 
                 rebuild_h5: bool = False, 
                 preprocess: Callable[..., Any] | None = None, 
                 transform: Callable[..., Any] | None = None, 
                 lazy_load: bool = True ) -> None:
        
        self.dataset_name = "PCAWG"

        self.resolutions = resolutions
        self.overlap = overlap
        self.h5_chunk_size = h5_chunk_size

        self.designed_subsets = []

        if type(designed_subsets) is str and designed_subsets == 'all':
            self.designed_subsets = self.subset_list
        else :
            for s in np.array([designed_subsets]).flatten():
                if s in self.subset_list:
                    self.designed_subsets.append(s)
                else:
                    logger.warning(f"{s} is not identified in subset list")

        self.source_list = [ f"{self.mirror}/{fn}_SNV.DEDUP.no_hypermut.annot.txt.gz" for fn in self.designed_subsets ]

        super().__init__(h5_path, raw_path, N_grams, logger, force_download, concurrent_download, rebuild_h5, preprocess, transform, lazy_load)

    def build_h5_summary(self):
        pass