import logging
from logging import Logger
from pathlib import Path
import numpy as np

from typing import Any, Callable, Dict, List, Optional, Tuple

from ._BioDataset import BioDigDriverfDataset


class PCAWG(BioDigDriverfDataset):

    """
    designed_sets = 'Sarcoma_tumors', 'Breast-AdenoCa', 'Breast-DCIS', 'Breast-LobularCa',
                'Carcinoma_tumors', 'Female_reproductive_system_tumors', 'Lymph-CLL',
                'Lymph-NOS', 'Myeloid-AML', 'Myeloid-MDS', 'Myeloid-MPN', 'Ovary-AdenoCA',
                'Panc-AdenoCA', 'Skin-Melanoma', 'CNS-PiloAstro', 'CNS_tumors', 'Squamous_tumors',
                'Stomach-AdenoCA', 'Kidney-RCC', 'Kidney_tumors', 'Liver-HCC', 'Lymph-BNHL',
                'Lymph_tumors', 'Myeloid_tumors', 'Pancan', 'Adenocarcinoma_tumors', 'Biliary-AdenoCA',
                'Eso-AdenoCa', 'Head-SCC', 'Hematopoietic_tumors',  'Bone-Cart', 'Bone-Epith',
                'Bone-Osteosarc','Breast_tumors', 'CNS-Medullo','Digestive_tract_tumors',
                'Glioma_tumors', 'Panc-Endocrine', 'Prost-AdenoCA'
    """

    mirror = "https://cb.csail.mit.edu/cb/DIG/downloads/mutation_files/PCAWG/ICGC_only"

    subset_list =['Sarcoma_tumors', 'Breast-AdenoCa', 'Breast-DCIS', 'Breast-LobularCa',
                'Carcinoma_tumors', 'Female_reproductive_system_tumors', 'Lymph-CLL',
                'Lymph-NOS', 'Myeloid-AML', 'Myeloid-MDS', 'Myeloid-MPN', 'Ovary-AdenoCA',
                'Panc-AdenoCA', 'Skin-Melanoma', 'CNS-PiloAstro', 'CNS_tumors', 'Squamous_tumors',
                'Stomach-AdenoCA', 'Kidney-RCC', 'Kidney_tumors', 'Liver-HCC', 'Lymph-BNHL',
                'Lymph_tumors', 'Myeloid_tumors', 'Pancan', 'Adenocarcinoma_tumors', 'Biliary-AdenoCA',
                'Eso-AdenoCa', 'Head-SCC', 'Hematopoietic_tumors',  'Bone-Cart', 'Bone-Epith',
                'Bone-Osteosarc','Breast_tumors', 'CNS-Medullo','Digestive_tract_tumors',
                'Glioma_tumors', 'Panc-Endocrine', 'Prost-AdenoCA']

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
        
        self.logger.debug("init PCAWG start")

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

        self.source_list = [ f"{self.mirror}/{fn}_SNV_MNV_INDEL.ICGC.annot.txt.gz" for fn in self.designed_subsets ]

        super().__init__(h5_path, raw_path, N_grams, logger, force_download, concurrent_download, rebuild_h5, preprocess, transform, lazy_load)

        self.logger.debug("init PCAWG end")

    def build_h5_summary(self):
        pass