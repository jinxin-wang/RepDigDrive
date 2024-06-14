import numpy as np

class GenomicPositionEncoding(object):
    def __init__(self, chr: int, chr_len: int, dmodel: int) -> None:
        self.chr = chr
        self.chr_len= chr_len
        self.dmodel = dmodel
        self.position = self._positional_encoding()

    def get_position_encoding(self, start:int, end: int, step: int = 1):
        idx = np.arange(start, end, step)
        return self.position[idx, :]
    
    def _positional_encoding(self):
        raise NotImplementedError


class LinearPositionEncoding(GenomicPositionEncoding):
    def __init__(self, chr: int, chr_len: int, dmodel: int, a: int = 10000, b: int = 2) -> None:
        super().__init__(chr, chr_len, dmodel)
        self.a = a
        self.b = b

    def _lin_pos_formula(self, a, b, i):
        return np.arange(self.chr_len)/np.power(a, b*i/self.dmodel)
    
    def _positional_encoding(self):
        self.position = np.zeros((self.chr_len, self.dmodel))
        for i in np.arange(self.dmodel):
            self.position[:, i] = self._formula(self.a, self.b, i)


class TrigonometricPositionEncoding(LinearPositionEncoding):
    def __init__(self) -> None:
        super().__init__()

    def _positional_encoding(self):
        # https://kazemnejad.com/blog/transformer_architecture_positional_encoding/
        pass