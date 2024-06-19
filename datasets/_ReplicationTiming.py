import os
import logging
import h5py

from pathlib import Path
from _BioDataset import BioBigWigDataset

from enum import Enum
from typing import Any, Callable, Dict, List, Optional, Tuple, Union

class ReplicationTimingDataset(BioBigWigDataset):

    """
    For each cell type, this track contains the following views:

    Percentage-normalized Signal: (view: PctSignal)
        Replication signal at 1 kb intervals as a percentage of normalized +/-25 kb tag densities for all cell cycle fractions (G1/G1b, S1, S2, S3, S4, G2).

    *Wavelet-smoothed Signal: (view: WaveSignal)
        Wavelet-smoothed transform of the six fraction profile that is a weighted average of the percentage-normalized signals such that earlier replication has higher values.

    Peaks: (view: Peaks)
        Local maxima in the wavelet-smoothed signal data corresponding to replication initiation (replication origin) zones.

    Valleys: (view: Valleys)
        Local minima in the wavelet-smoothed signal data corresponding to replication termination zones.

    *Summed Densities: (view: SumSignal)
        A measure of relative copy number at each genomic location that is the sum of normalized tag densities for each cell cycle fraction.

    More details about replication timing dataset, please refer to : https://genome.ucsc.edu/cgi-bin/hgTables?db=hg19&hgta_group=regulation&hgta_track=wgEncodeUwRepliSeq&hgta_table=wgEncodeUwRepliSeqImr90S3PctSignalRep1&hgta_doSchema=describe+table+schema
    """

    mirror = "https://hgdownload-test.gi.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeUwRepliSeq"
    
    class RepliSeqSignal(Enum):
        WaveSignal= 0
        SumSignal = 1

    class RepliSeqCell(Enum):
        BG02ES = 'Bg02es'
        BJ     = 'Bj'
        GM06990= 'Gm06990'
        GM12801= 'Gm12801'
        GM12812= 'Gm12812'
        GM12813= 'Gm12813'
        GM12878= 'Gm12878'
        HUVEC  = 'Huvec'
        HeLaS3 = 'Helas3'
        HepG2  = 'Hepg2'
        IMR90  = 'Imr90'
        K562   = 'K562'
        MCF7   = 'Mcf7'
        NHEK   = 'Nhek'
        SKNSH  = 'Sknsh'

    def __init__(
        self, 
        h5_path: Union[str, Path], 
        raw_path: Union[str, Path], 
        resolutions: List[int] = [10000, 100000],
        overlap : int = 0,
        h5_chunk_size: int = 100,
        logger = logging.getLogger(os.getcwd()),
        force_download = False,
        rebuild_h5 = False,
        design_signals: List[int] = [0,1],
        design_cells: List[str] | str = 'all',
        preprocess: Optional[Callable] = None, 
        transform:  Optional[Callable] = None
    ) -> None:
        
        self.dataset_name = "ReplicationTiming"
        self.design_signals = [ self.RepliSeqSignal(s) for s in design_signals ]
        if design_cells == 'all':
            self.design_cells = [ c for c in self.RepliSeqCell ]
        else:
            self.design_cells = [ self.RepliSeqCell(c) for c in design_cells ]
        self.preprocess = preprocess
        self.transform  = transform
        self.source_list  = [ f"{self.mirror}/{self._bigwig_fname(s.name, c.value)}" for c in self.RepliSeqCell for s in self.design_signals ] 

        list.sort(resolutions)

        super().__init__(h5_path = h5_path, 
                         raw_path=raw_path, 
                         resolutions = resolutions,
                         overlap = overlap,
                         h5_chunk_size = h5_chunk_size,
                         logger = logger,
                         force_download = force_download,
                         rebuild_h5 = rebuild_h5)
        
        self.rep_h5_fname = self.h5_path.joinpath(f"{self.dataset_name}.h5")

        mode = 'a'
        if rebuild_h5:
            mode = 'w'

        h5fd_dict = {}
        for signal in self.design_signals:
            for cell in self.design_cells:
                h5 = self._h5_fname(self._bigwig_fname(cell, signal))
                k  = self._bigwig_fname_key(cell, signal)
                h5fd_dict[k] = h5py.File(h5, 'r')

        with h5py.File(self.rep_h5_fname, mode=mode) as h5fd:
            for rslt in self.resolutions:
                dataset_name = self._dataset_name(rslt, self.overlap) 
                for chr in self.BigWigChm:
                    if self.rebuild_h5 or chr.name not in h5fd.keys() or dataset_name not in h5fd[chr.name].keys() :
                        self._concat_summary_table(h5fd, h5fd_dict, chr, rslt, overlap)
        
        for k in h5fd_dict:
            h5fd_dict[k].close()

        if self.preprocess is not None:
            self.preprocess(self.rep_h5_fname)

    def _bigwig_fname_key(self, cell: str, signal: str, replicate: int = 1):
        return f"{cell}{signal}Rep{replicate}"

    def _bigwig_fname(self, cell: str, signal: str, replicate: int = 1):
        return f"wgEncodeUwRepliSeq{self._bigwig_fname_key(cell, signal, replicate)}.bigWig"

    def _concat_summary_table(self, 
                              tgt_h5fd: h5py.File, 
                              src_h5fd_dict: Dict[str, h5py.File], 
                              chr: str, 
                              rslt: int, 
                              overlap: int
        ) -> h5py.File :

        L = -1

        for k in src_h5fd_dict:
            # check if the summary tables have same length
            ds = src_h5fd_dict[k][self._dataset_fullname(chr, rslt, overlap)]
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