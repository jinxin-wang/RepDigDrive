import yaml
import logging.config

from mini_utils import Singleton

@Singleton
class LogSingletonFactory(object):

    log_dict = {}

    def __init__(self, config_file: str = 'config/logging_config.yaml') -> None:
        self.config_file = config_file

    def _createLogger(self, name: str, config_file: str):
        with open(config_file, 'rt') as f:
            # Load the config file
            config = yaml.safe_load(f.read())
        # Configure the logging module with the config file
        logging.config.dictConfig(config)
        return logging.getLogger(name)


    def getLogger(self, name: str = None) -> logging.Logger:
        if name in self.log_dict :
            return self.log_dict[name]
        else:
            # Get a logger object
            self.log_dict[name] = self._createLogger(name, self.config_file)
            return self.log_dict[name]

