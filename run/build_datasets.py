from pathlib import Path

from datasets import *
from config import LogSingletonFactory, DatasetConfig

# logFactory = LogSingletonFactory()
# logger = logFactory.getLogger('development')
# datasetConfig = DatasetConfig()

def build_datasets(datasetConfig: DatasetConfig, logger):

    bioDataset_path = datasetConfig.getDatasetPath('Mappability', 'download')
    h5Dataset_path  = datasetConfig.getDatasetPath('Mappability', 'h5')

    MappabilityDataset(raw_path = bioDataset_path, 
                      h5_path = h5Dataset_path,
                      resolutions = [10000, 100000],
                      logger = logger,
                      force_download = True)

    bioDataset_path = datasetConfig.getDatasetPath('ReplicationTiming', 'download')
    h5Dataset_path  = datasetConfig.getDatasetPath('ReplicationTiming', 'h5')

    ReplicationTimingDataset(raw_path = bioDataset_path, 
                        h5_path = h5Dataset_path,
                        resolutions = [10000, 100000],
                        logger = logger,
                        force_download = True)

    bioDataset_path = datasetConfig.getDatasetPath('Epigenomics', 'download')
    h5Dataset_path  = datasetConfig.getDatasetPath('Epigenomics', 'h5')

    RoadmapEpigenomicsDataset(raw_path = bioDataset_path, 
                        h5_path = h5Dataset_path, 
                        resolutions = [10000, 100000], 
                        logger = logger,
                        force_download = True)
    

def build_datasets_test():

    logFactory = LogSingletonFactory()
    logger = logFactory.getLogger('development')

    # root_path = Path('/home/raf/Workspace/')
    # bioDataset_path = root_path.joinpath("RepDigDriver/Test/Datasets/BioDataset")
    # h5Dataset_path  = root_path.joinpath("RepDigDriver/Test/Datasets/BioDataset/h5")

    # MappabilityDataset(raw_path=bioDataset_path, 
    #                 h5_path =h5Dataset_path,
    #                 resolutions=[100000],
    #                 design_mers=[36],
    #                 logger=logger,
    #                 force_download=True)

    # ReplicationTimingDataset(raw_path=bioDataset_path, 
    #                     h5_path =h5Dataset_path,
    #                     resolutions=[100000],
    #                     design_signals=[0],
    #                     design_cells=['Bj'],
    #                     logger=logger,
    #                     force_download=True)

    # RoadmapEpigenomicsDataset(raw_path=bioDataset_path, 
    #                     h5_path = h5Dataset_path, 
    #                     resolutions = [100000], 
    #                     design_epig_modi = 'H3K27ac', 
    #                     design_cell_line = [5,6], 
    #                     logger = logger,
    #                     force_download=True)
    

    # root_path = Path('/home/raf/Workspace/')
    root_path = Path('D:/TMP/')
    dataset_path = root_path.joinpath("RepDigDriver/Test/Datasets/PCAWG")
    h5Dataset_path  = root_path.joinpath("RepDigDriver/Test/Datasets/PCAWG/h5")

    PCAWG(h5_path = h5Dataset_path, 
          raw_path = dataset_path,
          logger=logger,
          designed_subsets='Breast-DCIS')
