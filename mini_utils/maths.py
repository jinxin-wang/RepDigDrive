from torch import Tensor
import numpy as np
from typing import Any, Callable, Dict, List, Optional, Tuple, Union

def is_even(num) -> bool:
    if type(num) is Dict:
        num = [ num[k] for k in num ]
    return (Tensor(num) & 1) == 0

def is_odd(num) -> bool:
    return is_even(num) == False

def quat2dec(n: str) -> int:
    return int(n, 4)

def dec2quat(n: int) -> str:
    return np.base_repr(n, base=10)

