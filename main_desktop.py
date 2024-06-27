from config import LogSingletonFactory, DatasetConfig

from run import build_datasets, build_datasets_test

# logFactory = LogSingletonFactory()
# logger = logFactory.getLogger('development')
# datasetConfig = DatasetConfig()

# build_datasets(datasetConfig = datasetConfig, logger = logger)
build_datasets_test()