from ._BioDataset import BioDataset, BioBigWigDataset, BioMafDataset, BioBigWigDatasetFolder
from ._Mappability import MappabilityDataset
from ._ReplicationTiming import ReplicationTimingDataset

__all__ = (
    "BioDataset",
    "BioBigWigDataset", 
    "BioMafDataset", 
    "BioBigWigDatasetFolder",
    "MappabilityDataset",
    "ReplicationTimingDataset"
)

# def __getattr__(name):
#     if name in ("module_name",):
#         from module import module_name 
#         return module_name
#     raise AttributeError(f"module {__name__!r} has no attribute {name!r}")

