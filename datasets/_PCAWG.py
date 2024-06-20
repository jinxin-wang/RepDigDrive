import numpy as np
from config._log import DatasetsConfigFactory

from torch.utils.data import Dataset
from typing import Any, Callable, Dict, List, Optional, Tuple

class PCAWG(Dataset):
    def __init__(self) -> None:
        super().__init__()
        self.config = DatasetsConfigFactory.get_config()

    def __getitem__(self, index) -> Any:
        return super().__getitem__(index)
    
    def __len__(self):
        # key = self.keys[]
        pass

