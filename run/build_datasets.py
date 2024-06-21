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