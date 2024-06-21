from datasets import *
from config import LogSingletonFactory, DatasetConfig

from run import build_datasets

logFactory = LogSingletonFactory()
logger = logFactory.getLogger('development')
datasetConfig = DatasetConfig()

build_datasets(datasetConfig = datasetConfig, logger = logger)