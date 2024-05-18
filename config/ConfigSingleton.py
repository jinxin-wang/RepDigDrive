import yaml
from abc import ABC

class ABCConfigSingletonFactory(ABC):
    __config_singleton = None
    __config_path = None

    def __init__(self, config_path) -> None:
        super().__init__()
        self.__config_path = config_path
        self.set_config()

    @classmethod
    def _get_config(cls, config_path):
        if cls.__config_path is None:
            cls.__config_path = config_path

        if cls.__config_singleton is None:
            with open(cls.__config_path, 'r', encoding='utf-8') as fd:
                cls.__config_singleton = yaml.safe_load(fd)

        return cls.__config_singleton

class DatasetsConfigFactory(ABCConfigSingletonFactory):
    def __init__(self, config_path='config/datasets.yaml') -> None:
        super().__init__(config_path)

    @classmethod
    def get_config(cls, config_path = 'config/datasets.yaml'):
        return cls._get_config(config_path = config_path)

class ModelConfigFactory(ABCConfigSingletonFactory):
    def __init__(self, config_path = 'config/model.yaml') -> None:
        super().__init__(config_path)

    @classmethod
    def get_config(cls, config_path = None):
        return cls._get_config(config_path = 'config/model.yaml')

class TrainConfigFactory(ABCConfigSingletonFactory):
    def __init__(self, config_path = 'config/train.yaml') -> None:
        super().__init__(config_path)

    @classmethod
    def get_config(cls, config_path = None):
        return cls.__get_config(config_path='config/train.yaml')


