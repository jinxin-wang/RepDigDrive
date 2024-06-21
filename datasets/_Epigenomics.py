import os
import logging

import h5py

from pathlib import Path
from logging import Logger
from _BioDataset import BioBigWigDataset

from enum import Enum
from typing import Any, Callable, Dict, List, Optional, Tuple, Union

from mini_utils import enum_elt_list, enum_value_list

class RoadmapEpigenomicsDataset(BioBigWigDataset):
    """
    Assemble Roadmap Epigenomics Dataset : 
    1. download 723 Roadmap Epigenomics bigwig files
    2. resolve the track in a specific resolution
    3. encapsulate all the track into a table
    4. compress the table to h5 file

    'DNase', 'H3K27ac', 'H3K27me3', 'H3K36me3', 'H3K4me1', 'H3K4me3', 'H3K9ac', 'H3K9me3'
    """

    mirror = "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval"

    class DNase(Enum):
        E001 = 1;  E002 = 2;  E003 = 3;  E004 = 4;  E005 = 5;  E006 = 6;  E007 = 7;  E008 = 8;  E009 = 9
        E010 = 10; E011 = 11; E012 = 12; E013 = 13; E014 = 14; E015 = 15; E016 = 16; E017 = 17; E018 = 18
        E019 = 19; E020 = 20; E021 = 21; E022 = 22; E023 = 23; E024 = 24; E025 = 25; E026 = 26; E027 = 27
        E028 = 28; E029 = 29; E030 = 30; E031 = 31; E032 = 32; E033 = 33; E034 = 34; E035 = 35; E036 = 36
        E037 = 37; E038 = 38; E039 = 39; E040 = 40; E041 = 41; E042 = 42; E043 = 43; E044 = 44; E045 = 45
        E046 = 46; E047 = 47; E048 = 48; E049 = 49; E050 = 50; E051 = 51; E052 = 52; E053 = 53; E054 = 54
        E055 = 55; E056 = 56; E057 = 57; E058 = 58; E059 = 59; E060 = 60; E061 = 61; E062 = 62; E063 = 63
        E064 = 64; E065 = 65; E066 = 66; E067 = 67; E068 = 68; E069 = 69; E070 = 70; E071 = 71; E072 = 72
        E073 = 73; E074 = 74; E075 = 75; E076 = 76; E077 = 77; E078 = 78; E079 = 79; E080 = 80; E081 = 81
        E082 = 82; E083 = 83; E084 = 84; E085 = 85; E086 = 86; E087 = 87; E088 = 88; E089 = 89; E090 = 90
        E091 = 91; E092 = 92; E093 = 93; E094 = 94; E095 = 95; E096 = 96; E097 = 97; E098 = 98; E099 = 99
        E100 = 100; E101 = 101; E102 = 102; E103 = 103; E104 = 104; E105 = 105; E106 = 106; E107 = 107; E108 = 108
        E109 = 109; E110 = 110; E111 = 111
        
    class H3K27ac(Enum):
        E001 = 1;  E002 = 2;  E003 = 3;  E004 = 4;  E005 = 5;  E006 = 6;  E007 = 7;  E008 = 8;  E009 = 9
        E010 = 10; E011 = 11; E012 = 12; E013 = 13; E014 = 14; E015 = 15; E016 = 16; E017 = 17; E018 = 18
        E019 = 19; E020 = 20; E021 = 21; E022 = 22; E023 = 23; E024 = 24; E025 = 25; E026 = 26; E027 = 27
        E028 = 28; E029 = 29; E030 = 30; E031 = 31; E032 = 32; E033 = 33; E034 = 34; E035 = 35; E036 = 36
        E037 = 37; E038 = 38; E039 = 39; E040 = 40; E041 = 41; E042 = 42; E043 = 43; E044 = 44; E045 = 45
        E046 = 46; E047 = 47; E048 = 48; E049 = 49; E050 = 50; E051 = 51; E052 = 52; E053 = 53; E054 = 54
        E055 = 55; E056 = 56; E057 = 57; E058 = 58; E059 = 59; E061 = 61; E062 = 62; E063 = 63; E065 = 65
        E066 = 66; E067 = 67; E068 = 68; E069 = 69; E070 = 70; E071 = 71; E072 = 72; E073 = 73; E074 = 74
        E075 = 75; E076 = 76; E077 = 77; E078 = 78; E079 = 79; E080 = 80; E081 = 81; E082 = 82; E083 = 83
        E084 = 84; E085 = 85; E086 = 86; E087 = 87; E088 = 88; E089 = 89; E090 = 90; E091 = 91; E092 = 92
        E093 = 93; E094 = 94; E095 = 95; E096 = 96; E097 = 97; E098 = 98; E099 = 99; E100 = 100; E101 = 101
        E102 = 102; E103 = 103; E104 = 104; E105 = 105; E106 = 106; E107 = 107; E108 = 108; E109 = 109; E110 = 110
        E111 = 111; E112 = 112; E113 = 113

    class H3K27me3(Enum):
        E001 = 1;  E002 = 2;  E003 = 3;  E004 = 4;  E005 = 5;  E006 = 6;  E007 = 7;  E008 = 8;  E009 = 9
        E010 = 10; E011 = 11; E012 = 12; E013 = 13; E014 = 14; E015 = 15; E016 = 16; E017 = 17; E018 = 18
        E019 = 19; E020 = 20; E021 = 21; E022 = 22; E023 = 23; E024 = 24; E025 = 25; E026 = 26; E027 = 27
        E028 = 28; E029 = 29; E030 = 30; E031 = 31; E032 = 32; E033 = 33; E034 = 34; E035 = 35; E036 = 36
        E037 = 37; E038 = 38; E039 = 39; E040 = 40; E041 = 41; E042 = 42; E043 = 43; E044 = 44; E045 = 45
        E046 = 46; E047 = 47; E048 = 48; E049 = 49; E050 = 50; E051 = 51; E052 = 52; E053 = 53; E054 = 54
        E055 = 55; E056 = 56; E057 = 57; E058 = 58; E059 = 59; E061 = 61; E062 = 62; E063 = 63; E065 = 65
        E066 = 66; E067 = 67; E068 = 68; E069 = 69; E070 = 70; E071 = 71; E072 = 72; E073 = 73; E074 = 74
        E075 = 75; E076 = 76; E077 = 77; E078 = 78; E079 = 79; E080 = 80; E081 = 81; E082 = 82; E083 = 83
        E084 = 84; E085 = 85; E086 = 86; E087 = 87; E088 = 88; E089 = 89; E090 = 90; E091 = 91; E092 = 92
        E093 = 93; E094 = 94; E095 = 95; E096 = 96; E097 = 97; E098 = 98; E099 = 99; E100 = 100; E101 = 101
        E102 = 102; E103 = 103; E104 = 104; E105 = 105; E106 = 106; E107 = 107; E108 = 108; E109 = 109; E110 = 110
        E111 = 111; E112 = 112; E113 = 113

    class H3K36me3(Enum):
        E001 = 1;  E002 = 2;  E003 = 3;  E004 = 4;  E005 = 5;  E006 = 6;  E007 = 7;  E008 = 8;  E009 = 9
        E010 = 10; E011 = 11; E012 = 12; E013 = 13; E014 = 14; E015 = 15; E016 = 16; E017 = 17; E018 = 18
        E019 = 19; E020 = 20; E021 = 21; E022 = 22; E023 = 23; E024 = 24; E025 = 25; E026 = 26; E027 = 27
        E028 = 28; E029 = 29; E030 = 30; E031 = 31; E032 = 32; E033 = 33; E034 = 34; E035 = 35; E036 = 36
        E037 = 37; E038 = 38; E039 = 39; E040 = 40; E041 = 41; E042 = 42; E043 = 43; E044 = 44; E045 = 45
        E046 = 46; E047 = 47; E048 = 48; E049 = 49; E050 = 50; E051 = 51; E052 = 52; E053 = 53; E054 = 54
        E055 = 55; E056 = 56; E057 = 57; E058 = 58; E059 = 59; E061 = 61; E062 = 62; E063 = 63; E065 = 65
        E066 = 66; E067 = 67; E068 = 68; E069 = 69; E070 = 70; E071 = 71; E072 = 72; E073 = 73; E074 = 74
        E075 = 75; E076 = 76; E077 = 77; E078 = 78; E079 = 79; E080 = 80; E081 = 81; E082 = 82; E083 = 83
        E084 = 84; E085 = 85; E086 = 86; E087 = 87; E088 = 88; E089 = 89; E090 = 90; E091 = 91; E092 = 92
        E093 = 93; E094 = 94; E095 = 95; E096 = 96; E097 = 97; E098 = 98; E099 = 99; E100 = 100; E101 = 101
        E102 = 102; E103 = 103; E104 = 104; E105 = 105; E106 = 106; E107 = 107; E108 = 108; E109 = 109; E110 = 110
        E111 = 111

    class H3K4me1(Enum):
        E001 = 1;  E002 = 2;  E003 = 3;  E004 = 4;  E005 = 5;  E006 = 6;  E007 = 7;  E008 = 8;  E009 = 9
        E010 = 10; E011 = 11; E012 = 12; E013 = 13; E014 = 14; E015 = 15; E016 = 16; E017 = 17; E018 = 18
        E019 = 19; E020 = 20; E021 = 21; E022 = 22; E023 = 23; E024 = 24; E025 = 25; E026 = 26; E027 = 27
        E028 = 28; E029 = 29; E030 = 30; E031 = 31; E032 = 32; E033 = 33; E034 = 34; E035 = 35; E036 = 36
        E037 = 37; E038 = 38; E039 = 39; E040 = 40; E041 = 41; E042 = 42; E043 = 43; E044 = 44; E045 = 45
        E046 = 46; E047 = 47; E048 = 48; E049 = 49; E050 = 50; E051 = 51; E052 = 52; E053 = 53; E054 = 54
        E055 = 55; E056 = 56; E057 = 57; E058 = 58; E059 = 59; E061 = 61; E062 = 62; E063 = 63; E065 = 65
        E066 = 66; E067 = 67; E068 = 68; E069 = 69; E070 = 70; E071 = 71; E072 = 72; E073 = 73; E074 = 74
        E075 = 75; E076 = 76; E077 = 77; E078 = 78; E079 = 79; E080 = 80; E081 = 81; E082 = 82; E083 = 83
        E084 = 84; E085 = 85; E086 = 86; E087 = 87; E088 = 88; E089 = 89; E090 = 90; E091 = 91; E092 = 92
        E093 = 93; E094 = 94; E095 = 95; E096 = 96; E097 = 97; E098 = 98; E099 = 99; E100 = 100; E101 = 101
        E102 = 102; E103 = 103; E104 = 104; E105 = 105; E106 = 106; E107 = 107; E108 = 108; E109 = 109; E110 = 110
        E111 = 111

    class H3K4me3(Enum):
        E001 = 1;  E002 = 2;  E003 = 3;  E004 = 4;  E005 = 5;  E006 = 6;  E007 = 7;  E008 = 8;  E009 = 9
        E010 = 10; E011 = 11; E012 = 12; E013 = 13; E014 = 14; E015 = 15; E016 = 16; E017 = 17; E018 = 18
        E019 = 19; E020 = 20; E021 = 21; E022 = 22; E023 = 23; E024 = 24; E025 = 25; E026 = 26; E027 = 27
        E028 = 28; E029 = 29; E030 = 30; E031 = 31; E032 = 32; E033 = 33; E034 = 34; E035 = 35; E036 = 36
        E037 = 37; E038 = 38; E039 = 39; E040 = 40; E041 = 41; E042 = 42; E043 = 43; E044 = 44; E045 = 45
        E046 = 46; E047 = 47; E048 = 48; E049 = 49; E050 = 50; E051 = 51; E052 = 52; E053 = 53; E054 = 54
        E055 = 55; E056 = 56; E057 = 57; E058 = 58; E059 = 59; E061 = 61; E062 = 62; E063 = 63; E065 = 65
        E066 = 66; E067 = 67; E068 = 68; E069 = 69; E070 = 70; E071 = 71; E072 = 72; E073 = 73; E074 = 74
        E075 = 75; E076 = 76; E077 = 77; E078 = 78; E079 = 79; E080 = 80; E081 = 81; E082 = 82; E083 = 83
        E084 = 84; E085 = 85; E086 = 86; E087 = 87; E088 = 88; E089 = 89; E090 = 90; E091 = 91; E092 = 92
        E093 = 93; E094 = 94; E095 = 95; E096 = 96; E097 = 97; E098 = 98; E099 = 99; E100 = 100; E101 = 101
        E102 = 102; E103 = 103; E104 = 104; E105 = 105; E106 = 106; E107 = 107; E108 = 108; E109 = 109; E110 = 110
        E111 = 111

    class H3K9ac(Enum):
        E001 = 1;  E002 = 2;  E003 = 3;  E004 = 4;  E005 = 5;  E006 = 6;  E007 = 7;  E008 = 8;  E009 = 9
        E010 = 10; E011 = 11; E012 = 12; E013 = 13; E014 = 14; E015 = 15; E016 = 16; E017 = 17; E018 = 18
        E019 = 19; E020 = 20; E023 = 23; E025 = 25; E026 = 26; E027 = 27; E038 = 38; E047 = 47; E049 = 49
        E052 = 52; E062 = 62; E063 = 63; E066 = 66; E067 = 67; E068 = 68; E069 = 69; E072 = 72; E073 = 73
        E074 = 74; E075 = 75; E076 = 76; E077 = 77; E083 = 83; E086 = 86; E087 = 87; E088 = 88; E101 = 101
        E102 = 102; E103 = 103; E107 = 107; E108 = 108; E110 = 110; E111 = 111

    class H3K9me3(Enum):
        E001 = 1;  E002 = 2;  E003 = 3;  E004 = 4;  E005 = 5;  E006 = 6;  E007 = 7;  E008 = 8;  E009 = 9
        E010 = 10; E011 = 11; E012 = 12; E013 = 13; E014 = 14; E015 = 15; E016 = 16; E017 = 17; E018 = 18
        E019 = 19; E020 = 20; E021 = 21; E022 = 22; E023 = 23; E024 = 24; E025 = 25; E026 = 26; E027 = 27
        E028 = 28; E029 = 29; E030 = 30; E031 = 31; E032 = 32; E033 = 33; E034 = 34; E035 = 35; E036 = 36
        E037 = 37; E038 = 38; E039 = 39; E040 = 40; E041 = 41; E042 = 42; E043 = 43; E044 = 44; E045 = 45
        E046 = 46; E047 = 47; E048 = 48; E049 = 49; E050 = 50; E051 = 51; E052 = 52; E053 = 53; E054 = 54
        E055 = 55; E056 = 56; E057 = 57; E058 = 58; E059 = 59; E061 = 61; E062 = 62; E063 = 63; E065 = 65
        E066 = 66; E067 = 67; E068 = 68; E069 = 69; E070 = 70; E071 = 71; E072 = 72; E073 = 73; E074 = 74
        E075 = 75; E076 = 76; E077 = 77; E078 = 78; E079 = 79; E080 = 80; E081 = 81; E082 = 82; E083 = 83
        E084 = 84; E085 = 85; E086 = 86; E087 = 87; E088 = 88; E089 = 89; E090 = 90; E091 = 91; E092 = 92
        E093 = 93; E094 = 94; E095 = 95; E096 = 96; E097 = 97; E098 = 98; E099 = 99; E100 = 100; E101 = 101
        E102 = 102; E103 = 103; E104 = 104; E105 = 105; E106 = 106; E107 = 107; E108 = 108; E109 = 109; E110 = 110
        E111 = 111

    class H4K20me1(Enum):
        E001 = 1;  E002 = 2;  E003 = 3;  E004 = 4;  E005 = 5;  E006 = 6;  E007 = 7;  E008 = 8;  E009 = 9
        E010 = 10; E011 = 11; E012 = 12; E013 = 13; E014 = 14; E015 = 15; E016 = 16; E017 = 17; E018 = 18
        E019 = 19; E020 = 20; E021 = 21; E022 = 22; E023 = 23; E024 = 24; E025 = 25; E026 = 26; E027 = 27
        E028 = 28; E029 = 29; E030 = 30; E031 = 31; E032 = 32; E033 = 33; E034 = 34; E035 = 35; E036 = 36
        E037 = 37; E038 = 38; E039 = 39; E040 = 40; E041 = 41; E042 = 42; E043 = 43; E044 = 44; E045 = 45
        E046 = 46; E047 = 47; E048 = 48; E049 = 49; E050 = 50; E051 = 51; E052 = 52; E053 = 53; E054 = 54
        E055 = 55; E056 = 56; E057 = 57; E058 = 58; E059 = 59; E061 = 61; E062 = 62; E063 = 63; E065 = 65
        E066 = 66; E067 = 67; E068 = 68; E069 = 69; E070 = 70; E071 = 71; E072 = 72; E073 = 73; E074 = 74
        E075 = 75; E076 = 76; E077 = 77; E078 = 78; E079 = 79; E080 = 80; E081 = 81; E082 = 82; E083 = 83
        E084 = 84; E085 = 85; E086 = 86; E087 = 87; E088 = 88; E089 = 89; E090 = 90; E091 = 91; E092 = 92
        E093 = 93; E094 = 94; E095 = 95; E096 = 96; E097 = 97; E098 = 98; E099 = 99; E100 = 100; E101 = 101
        E102 = 102; E103 = 103; E104 = 104; E105 = 105; E106 = 106; E107 = 107; E108 = 108; E109 = 109; E110 = 110
        E111 = 111

    EpigMod = {
        'DNase':    DNase,    'H3K27ac': H3K27ac, 'H3K27me3': H3K27me3, 
        'H3K36me3': H3K36me3, 'H3K4me1': H3K4me1, 'H3K4me3': H3K4me3, 
        'H3K9ac':   H3K9ac,   'H3K9me3': H3K9me3
    }

    def __init__(
        self, 
        h5_path: Union[str, Path], 
        raw_path: Union[str, Path], 
        resolutions: list[int] = [10000, 100000],
        overlap : int = 0,
        h5_chunk_size: int = 100,
        logger: Union[str, Logger] = logging.getLogger(),
        force_download = False,
        rebuild_h5 = False,
        design_epig_modi: List[str] | str = 'all',
        design_cell_line: List[int] | str = 'all',
        preprocess: Optional[Callable] = None, 
        transform:  Optional[Callable] = None,
        lazy_load: bool = True
    ) -> None:
        
        self.dataset_name = "Epigenomics"
        self.h5_chunk_size= h5_chunk_size
        if design_epig_modi == 'all':
            # ['DNase', 'H3K27ac', 'H3K27me3', 'H3K36me3', 'H3K4me1', 'H3K4me3', 'H3K9ac', 'H3K9me3']
            self.design_epig_modi = list(self.EpigMod.values())
        else:
            self.design_epig_modi = [ self.EpigMod[m] for m in design_epig_modi ]

        self.design_cell_line = design_cell_line

        self.source_list = []
        for epig_modi in self.design_epig_modi:
            if self.design_cell_line == 'all':
                src = enum_elt_list(epig_modi, lambda x: self._bigwig_src(epig_cell_line = x))
                self.source_list.extend(src)
            else:
                cell_line_list = enum_value_list(epig_modi)
                for cell_line_num in self.design_cell_line:
                    if cell_line_num in cell_line_list:
                        self.source_list.append(self._bigwig_src(epig_cell_line = epig_modi(cell_line_num)))

        super().__init__(h5_path = h5_path, 
                         raw_path=raw_path, 
                         resolutions = resolutions,
                         overlap = overlap,
                         h5_chunk_size = h5_chunk_size,
                         logger = logger,
                         force_download = force_download,
                         rebuild_h5 = rebuild_h5,
                         preprocess = preprocess,
                         transform  = transform,
                         lazy_load  = lazy_load)

    def build_h5_summary(self):
        mode = 'a'
        if self.rebuild_h5:
            mode = 'w'

        for epig_modi in self.design_epig_modi: 
            h5fd_dict = {}
            epig_modi_name = epig_modi.__class__.__name__

            if self.design_cell_line == 'all':
                epig_modi = enum_elt_list(epig_modi)
            else: 
                epig_modi = [ epig_modi(c)  for c in self.design_cell_line if c in enum_value_list(epig_modi) ]

            for epig_cell_line in epig_modi:
                bigwig_fname = self._bigwig_fname(epig_cell_line)
                bigwig_fname = self.raw_path.joinpath(bigwig_fname)
                h5_fname     = self._h5_fname(bigwig_fname.name)
                ds_key = f"{epig_modi_name}_{epig_cell_line.name}"
                h5fd_dict[ds_key] = h5py.File(h5_fname, 'r') 

            self.logger.info(f"start building summary: {self.summary_h5_fname}")
            with h5py.File(self.summary_h5_fname, mode=mode) as h5fd:
                for rslt in self.resolutions:
                    dataset_name = self._dataset_name(rslt, self.overlap) 
                    for chr in self.BigWigChm:
                        ############# TODO: Rewrite the logic 
                        if self.rebuild_h5 or chr.name not in h5fd.keys() or dataset_name not in h5fd[chr.name].keys() :
                            self._concat_summary_table(h5fd, h5fd_dict, chr, rslt, self.overlap) 

            for k in h5fd_dict:
                h5fd_dict[k].close()

    def _bigwig_track_key(self, cell_line: Enum, epig_modi: str):
        return f"{cell_line.name}-{epig_modi}"

    def _bigwig_fname(self,  epig_cell_line: Enum):
        epig_modi = epig_cell_line.__class__.__name__
        key = self._bigwig_track_key(epig_cell_line, epig_modi)
        return f"{key}.pval.signal.bigwig"
    
    def _bigwig_src(self, epig_cell_line: Enum):
        return f"{self.mirror}/{self._bigwig_fname(epig_cell_line=epig_cell_line)}"
    
# source_list =  ["https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E001-H3K27me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E001-H3K36me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E001-H3K4me1.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E001-H3K4me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E001-H3K9ac.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E001-H3K9me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E002-H3K27me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E002-H3K36me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E002-H3K4me1.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E002-H3K4me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E002-H3K9ac.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E002-H3K9me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E003-DNase.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E003-H3K27ac.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E003-H3K27me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E003-H3K36me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E003-H3K4me1.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E003-H3K4me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E003-H3K9ac.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E003-H3K9me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E004-DNase.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E004-H3K27ac.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E004-H3K27me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E004-H3K36me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E004-H3K4me1.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E004-H3K4me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E004-H3K9ac.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E004-H3K9me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E005-DNase.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E005-H3K27ac.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E005-H3K27me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E005-H3K36me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E005-H3K4me1.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E005-H3K4me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E005-H3K9ac.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E005-H3K9me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E006-DNase.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E006-H3K27ac.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E006-H3K27me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E006-H3K36me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E006-H3K4me1.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E006-H3K4me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E006-H3K9ac.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E006-H3K9me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E007-DNase.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E007-H3K27ac.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E007-H3K27me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E007-H3K36me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E007-H3K4me1.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E007-H3K4me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E007-H3K9ac.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E007-H3K9me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E008-DNase.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E008-H3K27ac.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E008-H3K27me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E008-H3K36me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E008-H3K4me1.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E008-H3K4me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E008-H3K9ac.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E008-H3K9me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E009-H3K27me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E009-H3K36me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E009-H3K4me1.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E009-H3K4me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E009-H3K9me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E010-H3K27me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E010-H3K36me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E010-H3K4me1.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E010-H3K4me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E010-H3K9me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E011-H3K27ac.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E011-H3K27me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E011-H3K36me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E011-H3K4me1.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E011-H3K4me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E011-H3K9ac.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E011-H3K9me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E012-H3K27ac.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E012-H3K27me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E012-H3K36me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E012-H3K4me1.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E012-H3K4me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E012-H3K9me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E013-H3K27ac.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E013-H3K27me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E013-H3K36me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E013-H3K4me1.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E013-H3K4me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E013-H3K9me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E014-H3K27ac.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E014-H3K27me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E014-H3K36me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E014-H3K4me1.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E014-H3K4me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E014-H3K9ac.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E014-H3K9me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E015-H3K27ac.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E015-H3K27me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E015-H3K36me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E015-H3K4me1.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E015-H3K4me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E015-H3K9ac.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E015-H3K9me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E016-H3K27ac.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E016-H3K27me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E016-H3K36me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E016-H3K4me1.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E016-H3K4me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E016-H3K9ac.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E016-H3K9me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E017-DNase.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E017-H3K27ac.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E017-H3K27me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E017-H3K36me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E017-H3K4me1.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E017-H3K4me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E017-H3K9ac.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E017-H3K9me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E018-H3K27me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E018-H3K36me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E018-H3K4me1.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E018-H3K4me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E018-H3K9ac.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E018-H3K9me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E019-H3K27ac.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E019-H3K27me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E019-H3K36me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E019-H3K4me1.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E019-H3K4me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E019-H3K9ac.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E019-H3K9me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E020-H3K27ac.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E020-H3K27me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E020-H3K36me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E020-H3K4me1.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E020-H3K4me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E020-H3K9ac.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E020-H3K9me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E021-DNase.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E021-H3K27ac.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E021-H3K27me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E021-H3K36me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E021-H3K4me1.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E021-H3K4me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E021-H3K9me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E022-DNase.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E022-H3K27ac.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E022-H3K27me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E022-H3K36me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E022-H3K4me1.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E022-H3K4me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E022-H3K9me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E023-H3K27me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E023-H3K36me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E023-H3K4me1.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E023-H3K4me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E023-H3K9ac.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E023-H3K9me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E024-H3K27me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E024-H3K36me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E024-H3K4me1.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E024-H3K4me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E024-H3K9me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E025-H3K27me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E025-H3K36me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E025-H3K4me1.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E025-H3K4me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E025-H3K9ac.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E025-H3K9me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E026-H3K27ac.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E026-H3K27me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E026-H3K36me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E026-H3K4me1.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E026-H3K4me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E026-H3K9ac.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E026-H3K9me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E027-H3K27me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E027-H3K36me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E027-H3K4me1.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E027-H3K4me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E027-H3K9ac.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E027-H3K9me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E028-DNase.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E028-H3K27me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E028-H3K36me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E028-H3K4me1.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E028-H3K4me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E028-H3K9me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E029-DNase.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E029-H3K27ac.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E029-H3K27me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E029-H3K36me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E029-H3K4me1.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E029-H3K4me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E029-H3K9me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E030-H3K27me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E030-H3K36me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E030-H3K4me1.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E030-H3K4me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E030-H3K9me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E031-H3K27me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E031-H3K36me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E031-H3K4me1.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E031-H3K4me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E031-H3K9me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E032-DNase.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E032-H3K27ac.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E032-H3K27me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E032-H3K36me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E032-H3K4me1.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E032-H3K4me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E032-H3K9me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E033-DNase.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E033-H3K27me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E033-H3K36me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E033-H3K4me1.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E033-H3K4me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E033-H3K9me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E034-DNase.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E034-H3K27ac.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E034-H3K27me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E034-H3K36me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E034-H3K4me1.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E034-H3K4me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E034-H3K9me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E035-H3K27me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E035-H3K36me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E035-H3K4me1.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E035-H3K4me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E035-H3K9me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E036-H3K27me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E036-H3K36me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E036-H3K4me1.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E036-H3K4me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E036-H3K9me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E037-H3K27ac.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E037-H3K27me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E037-H3K36me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E037-H3K4me1.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E037-H3K4me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E037-H3K9me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E038-H3K27ac.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E038-H3K27me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E038-H3K36me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E038-H3K4me1.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E038-H3K4me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E038-H3K9ac.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E038-H3K9me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E039-H3K27ac.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E039-H3K27me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E039-H3K36me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E039-H3K4me1.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E039-H3K4me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E039-H3K9me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E040-H3K27ac.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E040-H3K27me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E040-H3K36me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E040-H3K4me1.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E040-H3K4me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E040-H3K9me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E041-H3K27ac.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E041-H3K27me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E041-H3K36me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E041-H3K4me1.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E041-H3K4me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E041-H3K9me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E042-H3K27ac.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E042-H3K27me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E042-H3K36me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E042-H3K4me1.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E042-H3K4me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E042-H3K9me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E043-H3K27ac.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E043-H3K27me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E043-H3K36me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E043-H3K4me1.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E043-H3K4me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E043-H3K9me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E044-H3K27ac.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E044-H3K27me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E044-H3K36me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E044-H3K4me1.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E044-H3K4me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E044-H3K9me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E045-H3K27ac.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E045-H3K27me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E045-H3K36me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E045-H3K4me1.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E045-H3K4me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E045-H3K9me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E046-DNase.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E046-H3K27ac.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E046-H3K27me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E046-H3K36me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E046-H3K4me1.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E046-H3K4me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E046-H3K9me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E047-H3K27ac.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E047-H3K27me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E047-H3K36me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E047-H3K4me1.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E047-H3K4me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E047-H3K9ac.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E047-H3K9me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E048-H3K27ac.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E048-H3K27me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E048-H3K36me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E048-H3K4me1.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E048-H3K4me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E048-H3K9me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E049-H3K27ac.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E049-H3K27me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E049-H3K36me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E049-H3K4me1.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E049-H3K4me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E049-H3K9ac.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E049-H3K9me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E050-DNase.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E050-H3K27ac.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E050-H3K27me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E050-H3K36me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E050-H3K4me1.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E050-H3K4me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E050-H3K9me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E051-DNase.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E051-H3K27me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E051-H3K36me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E051-H3K4me1.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E051-H3K4me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E051-H3K9me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E052-H3K27me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E052-H3K36me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E052-H3K4me1.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E052-H3K4me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E052-H3K9ac.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E052-H3K9me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E053-H3K27me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E053-H3K36me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E053-H3K4me1.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E053-H3K4me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E053-H3K9me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E054-H3K27me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E054-H3K36me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E054-H3K4me1.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E054-H3K4me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E054-H3K9me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E055-DNase.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E055-H3K27ac.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E055-H3K27me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E055-H3K36me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E055-H3K4me1.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E055-H3K4me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E055-H3K9me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E056-DNase.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E056-H3K27ac.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E056-H3K27me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E056-H3K36me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E056-H3K4me1.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E056-H3K4me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E056-H3K9me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E057-DNase.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E057-H3K27me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E057-H3K36me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E057-H3K4me1.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E057-H3K4me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E057-H3K9me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E058-H3K27ac.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E058-H3K27me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E058-H3K36me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E058-H3K4me1.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E058-H3K4me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E058-H3K9me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E059-DNase.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E059-H3K27ac.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E059-H3K27me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E059-H3K36me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E059-H3K4me1.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E059-H3K4me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E059-H3K9me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E061-H3K27ac.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E061-H3K27me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E061-H3K36me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E061-H3K4me1.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E061-H3K4me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E061-H3K9me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E062-H3K27ac.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E062-H3K27me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E062-H3K36me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E062-H3K4me1.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E062-H3K4me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E062-H3K9ac.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E062-H3K9me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E063-H3K27ac.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E063-H3K27me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E063-H3K36me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E063-H3K4me1.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E063-H3K4me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E063-H3K9ac.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E063-H3K9me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E065-H3K27ac.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E065-H3K27me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E065-H3K36me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E065-H3K4me1.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E065-H3K4me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E065-H3K9me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E066-H3K27ac.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E066-H3K27me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E066-H3K36me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E066-H3K4me1.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E066-H3K4me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E066-H3K9ac.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E066-H3K9me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E067-H3K27ac.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E067-H3K27me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E067-H3K36me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E067-H3K4me1.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E067-H3K4me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E067-H3K9ac.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E067-H3K9me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E068-H3K27ac.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E068-H3K27me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E068-H3K36me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E068-H3K4me1.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E068-H3K4me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E068-H3K9ac.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E068-H3K9me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E069-H3K27ac.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E069-H3K27me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E069-H3K36me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E069-H3K4me1.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E069-H3K4me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E069-H3K9ac.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E069-H3K9me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E070-H3K27me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E070-H3K36me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E070-H3K4me1.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E070-H3K4me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E070-H3K9me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E071-H3K27ac.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E071-H3K27me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E071-H3K36me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E071-H3K4me1.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E071-H3K4me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E071-H3K9me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E072-H3K27ac.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E072-H3K27me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E072-H3K36me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E072-H3K4me1.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E072-H3K4me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E072-H3K9ac.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E072-H3K9me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E073-H3K27ac.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E073-H3K27me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E073-H3K36me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E073-H3K4me1.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E073-H3K4me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E073-H3K9ac.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E073-H3K9me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E074-H3K27ac.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E074-H3K27me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E074-H3K36me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E074-H3K4me1.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E074-H3K4me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E074-H3K9ac.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E074-H3K9me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E075-H3K27ac.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E075-H3K27me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E075-H3K36me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E075-H3K4me1.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E075-H3K4me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E075-H3K9ac.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E075-H3K9me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E076-H3K27ac.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E076-H3K27me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E076-H3K36me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E076-H3K4me1.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E076-H3K4me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E076-H3K9ac.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E076-H3K9me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E077-H3K27me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E077-H3K36me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E077-H3K4me1.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E077-H3K4me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E077-H3K9ac.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E077-H3K9me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E078-H3K27ac.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E078-H3K27me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E078-H3K36me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E078-H3K4me1.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E078-H3K4me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E078-H3K9me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E079-H3K27ac.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E079-H3K27me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E079-H3K36me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E079-H3K4me1.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E079-H3K4me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E079-H3K9me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E080-DNase.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E080-H3K27ac.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E080-H3K27me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E080-H3K36me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E080-H3K4me1.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E080-H3K4me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E080-H3K9me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E081-DNase.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E081-H3K27me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E081-H3K36me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E081-H3K4me1.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E081-H3K4me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E081-H3K9me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E082-DNase.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E082-H3K27me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E082-H3K36me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E082-H3K4me1.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E082-H3K4me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E082-H3K9me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E083-DNase.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E083-H3K27me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E083-H3K36me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E083-H3K4me1.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E083-H3K4me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E083-H3K9ac.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E083-H3K9me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E084-DNase.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E084-H3K27ac.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E084-H3K27me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E084-H3K36me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E084-H3K4me1.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E084-H3K4me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E084-H3K9me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E085-DNase.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E085-H3K27ac.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E085-H3K27me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E085-H3K36me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E085-H3K4me1.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E085-H3K4me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E085-H3K9me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E086-DNase.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E086-H3K27me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E086-H3K36me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E086-H3K4me1.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E086-H3K4me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E086-H3K9ac.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E086-H3K9me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E087-H3K27ac.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E087-H3K27me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E087-H3K36me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E087-H3K4me1.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E087-H3K4me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E087-H3K9ac.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E087-H3K9me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E088-DNase.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E088-H3K27me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E088-H3K36me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E088-H3K4me1.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E088-H3K4me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E088-H3K9ac.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E088-H3K9me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E089-DNase.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E089-H3K27ac.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E089-H3K27me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E089-H3K36me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E089-H3K4me1.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E089-H3K4me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E089-H3K9me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E090-DNase.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E090-H3K27ac.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E090-H3K27me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E090-H3K36me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E090-H3K4me1.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E090-H3K4me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E090-H3K9me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E091-DNase.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E091-H3K27ac.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E091-H3K27me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E091-H3K36me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E091-H3K4me1.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E091-H3K4me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E091-H3K9me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E092-DNase.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E092-H3K27ac.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E092-H3K27me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E092-H3K36me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E092-H3K4me1.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E092-H3K4me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E092-H3K9me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E093-DNase.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E093-H3K27ac.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E093-H3K27me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E093-H3K36me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E093-H3K4me1.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E093-H3K4me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E093-H3K9me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E094-DNase.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E094-H3K27ac.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E094-H3K27me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E094-H3K36me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E094-H3K4me1.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E094-H3K4me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E094-H3K9me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E095-H3K27ac.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E095-H3K27me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E095-H3K36me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E095-H3K4me1.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E095-H3K4me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E095-H3K9me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E096-H3K27ac.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E096-H3K27me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E096-H3K36me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E096-H3K4me1.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E096-H3K4me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E096-H3K9me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E097-DNase.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E097-H3K27ac.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E097-H3K27me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E097-H3K36me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E097-H3K4me1.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E097-H3K4me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E097-H3K9me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E098-DNase.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E098-H3K27ac.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E098-H3K27me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E098-H3K36me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E098-H3K4me1.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E098-H3K4me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E098-H3K9me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E099-H3K27ac.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E099-H3K27me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E099-H3K36me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E099-H3K4me1.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E099-H3K4me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E099-H3K9me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E100-DNase.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E100-H3K27ac.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E100-H3K27me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E100-H3K36me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E100-H3K4me1.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E100-H3K4me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E100-H3K9me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E101-H3K27ac.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E101-H3K27me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E101-H3K36me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E101-H3K4me1.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E101-H3K4me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E101-H3K9ac.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E101-H3K9me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E102-H3K27ac.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E102-H3K27me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E102-H3K36me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E102-H3K4me1.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E102-H3K4me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E102-H3K9ac.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E102-H3K9me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E103-H3K27ac.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E103-H3K27me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E103-H3K36me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E103-H3K4me1.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E103-H3K4me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E103-H3K9ac.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E103-H3K9me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E104-H3K27ac.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E104-H3K27me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E104-H3K36me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E104-H3K4me1.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E104-H3K4me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E104-H3K9me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E105-H3K27ac.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E105-H3K27me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E105-H3K36me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E105-H3K4me1.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E105-H3K4me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E105-H3K9me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E106-H3K27ac.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E106-H3K27me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E106-H3K36me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E106-H3K4me1.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E106-H3K4me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E106-H3K9me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E107-H3K27me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E107-H3K36me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E107-H3K4me1.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E107-H3K4me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E107-H3K9ac.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E107-H3K9me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E108-H3K27ac.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E108-H3K27me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E108-H3K36me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E108-H3K4me1.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E108-H3K4me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E108-H3K9ac.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E108-H3K9me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E109-DNase.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E109-H3K27ac.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E109-H3K27me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E109-H3K36me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E109-H3K4me1.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E109-H3K4me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E109-H3K9me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E110-H3K27me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E110-H3K36me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E110-H3K4me1.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E110-H3K4me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E110-H3K9ac.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E110-H3K9me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E111-H3K27ac.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E111-H3K27me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E111-H3K36me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E111-H3K4me1.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E111-H3K4me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E111-H3K9ac.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E111-H3K9me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E112-H3K27ac.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E112-H3K27me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E112-H3K36me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E112-H3K4me1.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E112-H3K4me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E112-H3K9me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E113-H3K27ac.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E113-H3K27me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E113-H3K36me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E113-H3K4me1.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E113-H3K4me3.pval.signal.bigwig",
#  "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E113-H3K9me3.pval.signal.bigwig"]

# # List of cell line for each epigenomic modification
# DNase   = ["E003", "E004", "E005", "E006", "E007", "E008", "E017", "E021", "E022",
#            "E028", "E029", "E032", "E033", "E034", "E046", "E050", "E051", "E055",
#            "E056", "E057", "E059", "E080", "E081", "E082", "E083", "E084", "E085",
#            "E086", "E088", "E089", "E090", "E091", "E092", "E093", "E094", "E097",
#            "E098", "E100", "E109"]

# H3K27ac = ["E003", "E004", "E005", "E006", "E007", "E008", "E011", "E012", "E013",
#            "E014", "E015", "E016", "E017", "E019", "E020", "E021", "E022", "E026",
#            "E029", "E032", "E034", "E037", "E038", "E039", "E040", "E041", "E042",
#            "E043", "E044", "E045", "E046", "E047", "E048", "E049", "E050", "E055",
#            "E056", "E058", "E059", "E061", "E062", "E063", "E065", "E066", "E067",
#            "E068", "E069", "E071", "E072", "E073", "E074", "E075", "E076", "E078",
#            "E079", "E080", "E084", "E085", "E087", "E089", "E090", "E091", "E092", 
#            "E093", "E094", "E095", "E096", "E097", "E098", "E099", "E100", "E101",
#            "E102", "E103", "E104", "E105", "E106", "E108", "E109", "E111", "E112",
#            "E113"]

# H3K27me3= ["E001", "E002", "E003", "E004", "E005", "E006", "E007", "E008", "E009", 
#            "E010", "E011", "E012", "E013", "E014", "E015", "E016", "E017", "E018",
#            "E019", "E020", "E021", "E022", "E023", "E024", "E025", "E026", "E027", 
#            "E028", "E029", "E030", "E031", "E032", "E033", "E034", "E035", "E036",
#            "E037", "E038", "E039", "E040", "E041", "E042", "E043", "E044", "E045",
#            "E046", "E047", "E048", "E049", "E050", "E051", "E052", "E053", "E054",
#            "E055", "E056", "E057", "E058", "E059", "E061", "E062", "E063", "E065", 
#            "E066", "E067", "E068", "E069", "E070", "E071", "E072", "E073", "E074",
#            "E075", "E076", "E077", "E078", "E079", "E080", "E081", "E082", "E083", 
#            "E084", "E085", "E086", "E087", "E088", "E089", "E090", "E091", "E092",
#            "E093", "E094", "E095", "E096", "E097", "E098", "E099", "E100", "E101", 
#            "E102", "E103", "E104", "E105", "E106", "E107", "E108", "E109", "E110",
#            "E111", "E112", "E113"]

# H3K36me3= ["E001", "E002", "E003", "E004", "E005", "E006", "E007", "E008", "E009",
#            "E010", "E011", "E012", "E013", "E014", "E015", "E016", "E017", "E018",
#            "E019", "E020", "E021", "E022", "E023", "E024", "E025", "E026", "E027",
#            "E028", "E029", "E030", "E031", "E032", "E033", "E034", "E035", "E036",
#            "E037", "E038", "E039", "E040", "E041", "E042", "E043", "E044", "E045",
#            "E046", "E047", "E048", "E049", "E050", "E051", "E052", "E053", "E054",
#            "E055", "E056", "E057", "E058", "E059", "E061", "E062", "E063", "E065",
#            "E066", "E067", "E068", "E069", "E070", "E071", "E072", "E073", "E074",
#            "E075", "E076", "E077", "E078", "E079", "E080", "E081", "E082", "E083", 
#            "E084", "E085", "E086", "E087", "E088", "E089", "E090", "E091", "E092",
#            "E093", "E094", "E095", "E096", "E097", "E098", "E099", "E100", "E101",
#            "E102", "E103", "E104",  "E105", "E106", "E107","E108", "E109", "E110",
#            "E111", "E112", "E113"]

# H3K4me1 = ["E001", "E002", "E003", "E004", "E005", "E006", "E007", "E008", "E009", 
#            "E010", "E011", "E012", "E013", "E014", "E015", "E016", "E017", "E018",
#            "E019", "E020", "E021", "E022", "E023", "E024", "E025", "E026", "E027",
#            "E028", "E029", "E030", "E031", "E032", "E033", "E034", "E035", "E036",
#            "E037", "E038", "E039", "E040", "E041", "E042", "E043", "E044", "E045", 
#            "E046", "E047", "E048", "E049", "E050", "E051", "E052", "E053", "E054", 
#            "E055", "E056", "E057", "E058", "E059", "E061", "E062", "E063", "E065", 
#            "E066", "E067", "E068", "E069", "E070", "E071", "E072", "E073", "E074", 
#            "E075", "E076", "E077", "E078", "E079", "E080", "E081", "E082", "E083", 
#            "E084", "E085", "E086", "E087", "E088", "E089", "E090", "E091", "E092", 
#            "E093", "E094", "E095", "E096", "E097", "E098", "E099", "E100", "E101", 
#            "E102", "E103", "E104", "E105", "E106", "E107", "E108", "E109", "E110",
#            "E111", "E112", "E113"]     

# H3K4me3 = ["E001", "E002", "E003", "E004", "E005", "E006", "E007", "E008", "E009", 
#            "E010", "E011", "E012", "E013", "E014", "E015", "E016", "E017", "E018", 
#            "E019", "E020", "E021", "E022", "E023", "E024", "E025", "E026", "E027", 
#            "E028", "E029", "E030", "E031", "E032", "E033", "E034", "E035", "E036", 
#            "E037", "E038", "E039", "E040", "E041", "E042", "E043", "E044", "E045", 
#            "E046", "E047", "E048", "E049", "E051", "E052", "E053", "E054", "E055", 
#            "E056", "E057", "E058", "E059", "E061", "E062", "E063", "E065", "E066", 
#            "E067", "E068", "E069", "E070", "E071", "E072", "E074", "E075", "E076", 
#            "E077", "E078", "E079", "E080", "E081", "E082", "E083", "E084", "E085", 
#            "E086", "E087", "E088", "E089", "E090", "E091", "E092", "E093", "E094", 
#            "E095", "E096", "E097", "E098", "E099", "E100", "E101", "E102", "E103", 
#            "E104", "E105", "E106", "E107", "E108", "E109", "E110", "E111", "E112", 
#            "E113"]  
    
# H3K9ac  = ["E001", "E002", "E003", "E004", "E005", "E006", "E007", "E008", "E011", 
#            "E014", "E015", "E016", "E017", "E018", "E019", "E020", "E023", "E025", 
#            "E026", "E027", "E038", "E047", "E049", "E052", "E062", "E063", "E066", 
#            "E067", "E068", "E069", "E072", "E073", "E074", "E075", "E076", "E077", 
#            "E083", "E086", "E087", "E088", "E101", "E102", "E103", "E107", "E108",
#            "E110", "E111"]

# H3K9me3 = ["E001", "E002", "E003", "E004", "E005", "E006", "E007", "E008", "E009", 
#            "E010", "E011", "E012", "E013", "E014", "E015", "E016", "E017", "E018", 
#            "E019", "E020", "E021", "E022", "E023", "E024", "E025", "E026", "E027", 
#            "E028", "E029", "E030", "E031", "E032", "E033", "E034", "E035", "E036", 
#            "E037", "E038", "E039", "E040", "E041", "E042", "E043", "E044", "E045", 
#            "E046", "E047", "E048", "E049", "E050", "E051", "E052", "E053", "E054", 
#            "E055", "E056", "E057", "E058", "E059", "E061", "E062", "E063", "E065", 
#            "E066", "E067", "E068", "E069", "E070", "E071", "E072", "E073", "E074", 
#            "E075", "E076", "E077", "E078", "E079", "E080", "E081", "E082", "E083"]
