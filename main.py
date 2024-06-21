# import os
# import sys
# import yaml
# import torch
# import argparse
# import trainer
# # from utils import metric
# from models import *
from datasets import *
# from torch.utils.data import DataLoader

from pathlib import Path

# def parse_args():
#     parser = argparse.ArgumentParser()
#     parser.add_argument('--model-config-path', dest='model_config_path', type=str, default=r"./config/train_config.yaml")
#     # .....
#     return parser.parse_args()

# def main(args):
#     print(PCAWG().config)

# if __name__ == '__main__':
#     args = parse_args()
#     main(args)

from config import LogSingletonFactory, DatasetConfig

logFactory = LogSingletonFactory()

logger = logFactory.getLogger('development')

# root_path = Path('/home/raf/Workspace/')
# bioDataset_path = root_path.joinpath("RepDigDriver/Test/Datasets/BioDataset")
# h5Dataset_path  = root_path.joinpath("RepDigDriver/Test/Datasets/BioDataset/h5")

bioDataset_path = DatasetConfig.getDatasetPath('Epigenomics', 'download')
h5Dataset_path  = DatasetConfig.getDatasetPath('Epigenomics', 'h5')

# BioDataset(raw_path=bioDataset_path, 
#            logger=logger,
#            force_download=True)

# from tests.datasets_testcases import BioBigWigDataset_TestCase

# BioBigWigDataset(raw_path=bioDataset_path, 
#                  h5_path =h5Dataset_path,
#                  resolutions=[10000, 100000],
#                  logger=logger,
#                  force_download=True)

# MappabilityDataset(raw_path=bioDataset_path, 
#                   h5_path =h5Dataset_path,
#                   resolutions=[100000],
#                   design_mers=[36],
#                   logger=logger,
#                   force_download=True)

# ReplicationTimingDataset(raw_path=bioDataset_path, 
#                     h5_path =h5Dataset_path,
#                     resolutions=[100000],
#                     design_signals=[0],
#                     design_cells=['Bj'],
#                     logger=logger,
#                     force_download=True)

RoadmapEpigenomicsDataset(raw_path=bioDataset_path, 
                    h5_path = h5Dataset_path, 
                    resolutions = [100000], 
                    design_epig_modi = 'H3K27ac', 
                    design_cell_line = [5,6], 
                    logger = logger,
                    force_download=True)