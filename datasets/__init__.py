from ._BioDataset import BioDataset, BioBigWigDataset, BioMafDataset, BioDigDriverfDataset

from ._Mappability import MappabilityDataset
from ._ReplicationTiming import ReplicationTimingDataset
from ._Epigenomics import RoadmapEpigenomicsDataset

from ._PCAWG import PCAWGDataset

__all__ = (
    "BioDataset",
    "BioBigWigDataset", 
    "BioMafDataset", 
    "MappabilityDataset",
    "ReplicationTimingDataset",
    "RoadmapEpigenomicsDataset", 
    "BioDigDriverfDataset", "PCAWGDataset",
)

# def __getattr__(name):
#     if name in ("module_name",):
#         from module import module_name 
#         return module_name
#     raise AttributeError(f"module {__name__!r} has no attribute {name!r}")

