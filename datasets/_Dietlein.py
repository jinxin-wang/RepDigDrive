import numpy as np
from torch.utils.data import Dataset

from torchvision import transforms, datasets

# class DatasetExample(Dataset):
#     def __init__(self, data_path, keys = ('x', 'y')) -> None:
#         super().__init__()

#     def __getitem__(self, index) -> Any:
#         return super().__getitem__(index)
    
#     def __len__(self):
#         pass

from torchvision import datasets