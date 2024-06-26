from ._BioDataset import BioDataset, BioBigWigDataset, BioMafDataset, BioDigDriverfDataset

from ._Mappability import MappabilityDataset
from ._ReplicationTiming import ReplicationTimingDataset
from ._Epigenomics import RoadmapEpigenomicsDataset

from ._PCAWG import PCAWG

__all__ = (
    "BioDataset",
    "BioBigWigDataset", 
    "BioMafDataset", 
    "MappabilityDataset",
    "ReplicationTimingDataset",
    "RoadmapEpigenomicsDataset", 
    "BioDigDriverfDataset", "PCAWG",
)

# def __getattr__(name):
#     if name in ("module_name",):
#         from module import module_name 
#         return module_name
#     raise AttributeError(f"module {__name__!r} has no attribute {name!r}")

