import torch 
import numpy as np
from typing import Any, Callable, Dict, List, Optional, Tuple, Union

def is_even(num: Union[List[int], Tuple[int], np.array[int], torch.Tensor[int]]) -> torch.Tensor:
    if type(num) is Dict:
        num = [ num[k] for k in num ]
    return (torch.Tensor(num) & 1) == 0

def is_odd(num: Union[List[int], Tuple[int], np.array[int], torch.Tensor[int]]) -> torch.Tensor:
    return is_even(num) == False

def quat2dec(n: str) -> int:
    return int(n, 4)

def dec2quat(n: int) -> str:
    return np.base_repr(n, base=10)

