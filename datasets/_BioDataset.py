import os
import h5py
import logging
from pathlib import Path
from typing import Any, Callable, Dict, List, Optional, Tuple, Union

from torch.utils.data import Dataset

class BioH5Dataset(Dataset):
    def __init__(
        self, 
        h5_path: Union[str, Path], 
        raw_path: Union[str, Path], 
        logger = logging.getLogger(os.getcwd()),
        force_download = False,
        rebuild_h5 = False,
    ) -> None:
        
        self.h5_path  = Path(h5_path)
        self.raw_path = Path(raw_path)
        self.logger   = logger 
    
        if force_download:
            self.download_rawdata()

        if (not os.path.isfile(self.h5_path)) or rebuild_h5: 
            self.logger.info(f"h5 dataset [{self.h5_path}] is not found")
            if not os.path.exists(self.raw_path):
                self.logger.info(f"raw dataset [{self.raw_path}] is not found. ")
                self.download_rawdata()
            self.build_h5()

    # def __getitem__(self, index) -> Any:
    #     return super().__getitem__(index)
    
    # def __len__(self):
    #     pass
    
    def build_h5(self):
        pass

    def download_rawdata(self):
        pass
        

class BioDatasetFolder(Dataset):
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
