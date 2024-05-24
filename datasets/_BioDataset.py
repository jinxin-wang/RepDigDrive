import os
from pathlib import Path
from typing import Any, Callable, Dict, List, Optional, Tuple, Union

from torch.utils.data import Dataset

class BioDataset(Dataset):
    def __init__(
        self, 
        raw_path: Union[str, Path], 
        resolution: int, 
        preprocess: Callable, 
        transform:  Optional[Callable] = None,
        transname: Optional[str] = None
    ) -> None:
        
        self.raw_path = Path(raw_path)
        self.resolution = resolution
        self.preprocessed_path = self.raw_path.parent.joinpath(f"{resolution}").joinpath(self.raw_path.name)
        if transname:
            self.transform_path = self.raw_path.parent.joinpath(f"{resolution}").joinpath(transname).joinpath(self.raw_path.name)

        if os.path.isdir():
            raise Exception("Expect path of file")
        
        if not os.path.isfile(self.preprocessed_path):
            preprocess(self.raw_path, self.preprocessed_path, self.resolution)

        if transform and not os.path.isfile(self.transform_path):
            transform(self.preprocessed_path, self.transform_path)


    def __getitem__(self, index) -> Any:
        return super().__getitem__(index)
    
    def __len__(self):
        # key = self.keys[]
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
