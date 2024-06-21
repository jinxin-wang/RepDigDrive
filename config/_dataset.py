import yaml

from mini_utils import Singleton

@Singleton
class DatasetConfig(object):

    def __init__(self, config_file: str = "config/datasets.yaml") -> None:
        self.config_file = config_file

        with open(self.config_file, 'r') as fd:
            self.config = yaml.safe_load(fd)

    def getConfigAttr(self, name):
        return self.config[name]
    
    def getDatasetNames(self):
        return list(self.config.keys())

