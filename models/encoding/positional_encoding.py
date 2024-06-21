import numpy as np
from typing import Any, Callable, Dict, List, Optional, Tuple, Union

from numpy.core.multiarray import array as array
import torch 
from torch import nn
from abc import abstractmethod

class GenomicPositionalEncoding(nn.Module):

    def __init__(self, 
                 d_model: Union[int, List[int], Tuple[int], np.array, torch.Tensor, Dict[str,int]], 
                 dropout: float = -1., 
                 length: int = 100
        ) -> None:

        super().__init__()

        # there is no odd number in d_model
        assert sum(self._is_odd(d_model)) == 0

        self.d_model= d_model
        self.length = length

        if type(d_model) is int :
            _d_model = torch.Tensor([d_model])
        elif type(d_model) is Dict :
            _d_model = torch.Tensor([d_model[k] for k in d_model])
        else:
            _d_model = torch.Tensor(d_model)

        self._d_model= _d_model

        # initialize dropout
        if dropout > 0:
            self.dropout = nn.Dropout(p = dropout)
        else:
            self.dropout = None

    def _is_even(self, num: Union[List[int], Tuple[int], np.array[int], torch.Tensor[int]]) -> torch.Tensor:
        if type(num) is Dict:
            num = [ num[k] for k in num ]
        return (torch.Tensor(num) & 1) == 0

    def _is_odd(self,  num: Union[List[int], Tuple[int], np.array[int], torch.Tensor[int]]) -> torch.Tensor:
        return self._is_even(num) == False

    @abstractmethod
    def get_position_encoding(self, t: Union[List[int], np.array[int], torch.Tensor[int]]) -> torch.Tensor:
        raise NotImplementedError

class GenomicLinearPositionalEncoding(GenomicPositionalEncoding):
    def __init__(self, 
                 d_model: Union[List[int], Tuple[int], np.array[int], torch.Tensor[int], Dict[str,int]], 
                 dropout: float = -1, 
                 length: int = 100,
                 total_progression: int = 10000, 
        ) -> None:

        super().__init__(d_model, dropout, length)
        self.total_progression = total_progression

    def _lin_pos_formula(self, 
                         d: int,
                         t:Union[List[int], np.array[int], torch.Tensor[int]], 
                         k:Union[List[int], np.array[int], torch.Tensor[int]], 
                         total_progression:int = 10000
        ) -> torch.Tensor[float]:

        t = torch.Tensor(t)
        k = torch.Tensor(k)
        return t/torch.pow(total_progression, 2*k/d)
    
    def get_position_encoding(self, t:  Union[List[int], np.array[int], torch.Tensor[int]]) -> torch.Tensor:
        pe = [ self._lin_pos_formula(d, t, torch.arange(d).unsqueeze(0).view(-1,1) * torch.ones((d,t.shape[-1]))) for d in self._d_model ]

class TrigonometricPositionalEncoding(GenomicLinearPositionalEncoding):
    def __init__(self, 
                 d_model: Union[int, List[int], Tuple[int], np.array[int], torch.Tensor[int], Dict[str,int]], 
                 dropout: float = -1, 
                 length: int = 100, 
                 total_progression: int = 10000, 
        ) -> None:

        super().__init__(d_model, dropout, length, total_progression)

    def _tri_pos_formula(self, 
                         t: Union[int, List[int], np.array[int], torch.Tensor[int]],
                         k: Union[int, List[int], np.array[int], torch.Tensor[int]]
        ) -> np.array[float]:

        # https://kazemnejad.com/blog/transformer_architecture_positional_encoding/
        pe = np.zeros((self.dmodel,self.chr_len))
        pe[0::2,:] = np.sin(self._lin_pos_formula(t, k[0::2, :]))
        pe[1::2,:] = np.cos(self._lin_pos_formula(t, k[1::2, :]))
        return pe

class RelativePositionalEncoding(TrigonometricPositionalEncoding):
    def __init__(self, 
                 d_model: Union[int, List[int], Tuple[int], np.array[int], torch.Tensor[int], Dict[str,int]], 
                 dropout: float = -1, 
                 length: int = 100,
                 total_progression: int = 10000, 
        ) -> None:

        super().__init__(d_model, dropout, length, total_progression)


pe = GenomicLinearPositionalEncoding([3, 8, 8], 0.6, 100)