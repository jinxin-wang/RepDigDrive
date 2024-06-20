import os
import logging
import h5py

from pathlib import Path
from datasets import BioBigWigDataset

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
        self.source_list  = [ f"{self.mirror}/{self._bigwig_fname(cell=c,signal=s)}" for c in self.design_cells for s in self.design_signals ] 

        logger.debug(f"sources list: {self.source_list}")

        super().__init__(h5_path = h5_path, 
                         raw_path=raw_path, 
                         resolutions = resolutions,
                         overlap = overlap,
                         h5_chunk_size = h5_chunk_size,
                         logger = logger,
                         force_download = force_download,
                         rebuild_h5 = rebuild_h5,
                         preprocess = preprocess,
                         transform  = transform)

    def build_h5_summary(self):
        mode = 'a'
        if self.rebuild_h5:
            mode = 'w'

        h5fd_dict = {}
        
        for signal in self.design_signals:
            for cell in self.design_cells:
                h5 = self._h5_fname(self._bigwig_fname(cell, signal))
                k  = self._bigwig_fname_key(cell, signal)
                h5fd_dict[k] = h5py.File(h5, 'r')

        self.logger.info(f"start building summary: {self.summary_h5_fname}")
        with h5py.File(self.summary_h5_fname, mode=mode) as h5fd:
            for rslt in self.resolutions:
                dataset_name = self._dataset_name(rslt, self.overlap) 
                for chr in self.BigWigChm:
                    if self.rebuild_h5 or chr.name not in h5fd.keys() or dataset_name not in h5fd[chr.name].keys() :
                        self._concat_summary_table(h5fd, h5fd_dict, chr, rslt, self.overlap)

        for k in h5fd_dict:
            h5fd_dict[k].close()

    def _bigwig_fname_key(self, cell: RepliSeqCell, signal: RepliSeqSignal, replicate: int = 1):
        return f"{cell.value}{signal.name}Rep{replicate}"

    def _bigwig_fname(self, cell: RepliSeqCell, signal: RepliSeqSignal, replicate: int = 1):
        return f"wgEncodeUwRepliSeq{self._bigwig_fname_key(cell, signal, replicate)}.bigWig"

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