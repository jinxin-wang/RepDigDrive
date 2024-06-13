import os
import logging
import h5py

from pathlib import Path
from _BioDataset import BioBigWigDataset

from enum import Enum
from typing import Any, Callable, Dict, List, Optional, Tuple, Union

class ReplicationTimingDataset(BioBigWigDataset):

    mirror = "https://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeMapability"
    
    class WINSIZE(Enum):
        Align24mer = 24
        Align36mer = 36
        Align40mer = 40
        Align50mer = 50
        Align75mer = 75
        Align100mer= 100

    def __init__(
        self, 
        h5_path: Union[str, Path], 
        raw_path: Union[str, Path], 
        resolutions: list[int] = [10000, 100000],
        overlap : int = 0,
        h5_chunk_size: int = 100,
        logger = logging.getLogger(os.getcwd()),
        force_download = False,
        rebuild_h5 = False,
        design_mers: List[int] = [24, 36, 40, 50, 75, 100],
        preprocess: Optional[Callable] = None, 
        transform:  Optional[Callable] = None
    ) -> None:
        
        self.h5_chunk_size= h5_chunk_size
        self.dataset_name = "Mappability"
        self.design_mers = design_mers
        self.preprocess = preprocess
        self.transform  = transform
        self.source_list  = [ f"{self.mirror}/{self._bigwig_fname(self.WINSIZE(alg).name)}" for alg in self.design_mers ] 

        super().__init__(h5_path = h5_path, 
                         raw_path=raw_path, 
                         resolutions=resolutions,
                         overlap=overlap,
                         logger=logger,
                         force_download=force_download,
                         rebuild_h5=rebuild_h5)
        
        self.mapp_h5_fname = self.h5_path.joinpath(f"{self.dataset_name}.h5")

        mode = 'a'
        if rebuild_h5:
            mode = 'w'

        with h5py.File(self.mapp_h5_fname, mode=mode) as h5fd:
            for rslt in self.resolutions:
                dataset_name = self._dataset_name(rslt, self.overlap) 
                for chr in self.BigWigChm:
                    if self.rebuild_h5 or chr.name not in h5fd.keys() or dataset_name not in h5fd[chr.name].keys() :
                        # open designed h5 files 
                        h5fd_dict = { self._h5_fname_to_mer(h5): h5py.File(h5, 'r') for h5 in self.h5_list if self._h5_fname_to_mer(h5).value in self.design_mers }
                        self._concat_summary_table(h5fd, h5fd_dict, chr, rslt, overlap)
        
        if self.preprocess is not None:
            self.preprocess(self.mapp_h5_fname)

    def _bigwig_fname(self, key):
        return f"wgEncodeCrgMapability{key}.bigWig"

    def _h5_fname_to_mer(self, h5_fname: Path):
        """
        given filename, return mer enum element
        """
        return self.WINSIZE(int(h5_fname.name.split('mer')[0].replace('wgEncodeCrgMapabilityAlign','')))

    def _concat_summary_table(self, 
                              tgt_h5fd: h5py.File, 
                              src_h5fd_dict: dict[WINSIZE, h5py.File], 
                              chr: str, 
                              rslt: int, 
                              overlap: int
        ) -> h5py.File :

        L = -1

        for mer in src_h5fd_dict:
            # check if the summary tables have same length
            ds = src_h5fd_dict[mer][self._dataset_fullname(chr, rslt, overlap)]
            if L < 0 :
                assert ds.shape[0] > 0
                L = ds.shape[0]

            else :
                assert L == ds.shape[0]
                columns_count += ds.shape[1]
        
        # create chunked dataset
        tgt_h5fd.create_dataset(name=self._dataset_fullname(chr, rslt, overlap), 
                            shape=(L, columns_count), 
                            chunks=(self.h5_chunk_size, columns_count), 
                            dtype=float)

        # update to tgt_h5fd dataset one by one, since each one can be very large
        column_names = []
        columns_idx  = 0
        for mer, fd in src_h5fd_dict.items():
            ds = fd[self._dataset_fullname(chr, rslt, overlap)]
            attrs = fd.attrs[self.H5Attrs.COLUMNS.value]
            column_names += [ f"{mer}_{s}" for s in attrs ]
            tgt_h5fd[self._dataset_fullname(chr,rslt,overlap)][:,columns_idx: columns_idx + ds.shape[1]] = ds[:]
            columns_idx += ds.shape[1]

        tgt_h5fd.attrs[self.H5Attrs.COLUMNS.value] = column_names
        return tgt_h5fd
        



# addr_list = ["https://hgdownload-test.gi.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeUwRepliSeq/wgEncodeUwRepliSeqBg02esWaveSignalRep1.bigWig",
#             "https://hgdownload-test.gi.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeUwRepliSeq/wgEncodeUwRepliSeqBjWaveSignalRep2.bigWig",
#             "https://hgdownload-test.gi.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeUwRepliSeq/wgEncodeUwRepliSeqGm12878WaveSignalRep1.bigWig",
#             "https://hgdownload-test.gi.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeUwRepliSeq/wgEncodeUwRepliSeqHelas3WaveSignalRep1.bigWig",
#             "https://hgdownload-test.gi.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeUwRepliSeq/wgEncodeUwRepliSeqHepg2WaveSignalRep1.bigWig",
#             "https://hgdownload-test.gi.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeUwRepliSeq/wgEncodeUwRepliSeqHuvecWaveSignalRep1.bigWig",
#             "https://hgdownload-test.gi.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeUwRepliSeq/wgEncodeUwRepliSeqImr90WaveSignalRep1.bigWig",
#             "https://hgdownload-test.gi.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeUwRepliSeq/wgEncodeUwRepliSeqMcf7WaveSignalRep1.bigWig",
#             "https://hgdownload-test.gi.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeUwRepliSeq/wgEncodeUwRepliSeqNhekWaveSignalRep1.bigWig",
#             "https://hgdownload-test.gi.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeUwRepliSeq/wgEncodeUwRepliSeqSknshWaveSignalRep1.bigWig"]