import yaml
from pathlib import Path

from mini_utils.singleton import Singleton

@Singleton
class DatasetConfig(object):

    def __init__(self, config_file: str = "config/datasets.yaml") -> None:
        self.config_file = config_file

        with open(self.config_file, 'r') as fd:
            self.config = yaml.safe_load(fd)

    def getDatasetRoot(self, sub):
        return Path(self.config["root"][sub])
    
    def getDatasetPath(self, name: str, sub: str):
        root = self.getDatasetRoot(sub)
        dspath = root.joinpath(self.config[name][sub])
        return dspath

    def getDatasetNames(self):
        return list(self.config.keys())

