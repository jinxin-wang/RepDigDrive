from ._singleton import Singleton
from ._convert import enum_name_list, enum_value_list, enum_elt_list
from ._bio import Chm, Nucl, build_N_gram_nucl_enum, BigWigChromSizesDict, BASE_SUBSTITUTION_TYPES, BASE_SUBSTITUTION_CLASSES, build_N_ctx_mut_enum
from ._maths import is_even, is_odd, quat2dec, dec2quat


import _bio as bio
import _convert as convert
import _maths as maths
import _singleton as singleton

__all__ = (
    "bio", "convert", "maths", "singleton",
    "Singleton", 
    "enum_name_list", "enum_value_list", "enum_elt_list",
    "is_even", "is_odd", "quat2dec", "dec2quat",
    "Chm", "Nucl", "build_N_gram_nucl_enum", "build_N_ctx_mut_enum", "BigWigChromSizesDict","BASE_SUBSTITUTION_TYPES", "BASE_SUBSTITUTION_CLASSES"
)