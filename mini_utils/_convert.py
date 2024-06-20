from enum import Enum
from typing import Any, Callable, Dict, List, Optional, Tuple, Union

def enum_elt_list(enum_cls: Enum, func: Callable = lambda x: x):
    return [ func(e) for e in enum_cls ]

def enum_value_list(enum_cls: Enum, func: Callable = lambda x: x):
    return [ func(e.value) for e in enum_cls ]

def enum_name_list(enum_cls: Enum, func: Callable = lambda x: x):
    return [ func(e.name)  for e in enum_cls ]



# class H3K27ac(Enum):
#     E001 = 1
#     E002 = 2
#     E003 = 3


# # print(enum_value_list(H3K27ac))
# # print(enum_name_list(H3K27ac))

# print(H3K27ac(1).__class__.__name__)