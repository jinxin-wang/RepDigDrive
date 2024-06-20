import yaml
import logging
import unittest

from config import LogSingletonFactory
from unittest.mock import patch, mock_open, MagicMock

# logFactory = LogSingletonFactory()
# logger = logFactory.getLogger('development')

# # Log some messages
# logger.debug('This is a debug message')
# logger.info('This is an info message')
# logger.warning('This is a warning message')
# logger.error('This is an error message')
# logger.critical('This is a critical message')

class TestLogSingletonFactory(unittest.TestCase):

    def test_singleton_behavior(self):
        factory1 = LogSingletonFactory()
        factory2 = LogSingletonFactory()
        self.assertIs(factory1, factory2, "LogSingletonFactory is not a singleton!")

    @patch('builtins.open', new_callable=mock_open, read_data="{}")
    @patch('yaml.safe_load', return_value={})
    @patch('logging.config.dictConfig')
    def test_create_logger(self, mock_dictConfig, mock_safe_load, mock_open):
        factory = LogSingletonFactory()
        logger_name = 'development'
        
        # Create logger
        logger = factory.getLogger(logger_name)

        # Ensure logger is created
        self.assertIsInstance(logger, logging.Logger)

        # Ensure logger is cached
        self.assertIn(logger_name, factory.log_dict)
        self.assertIs(logger, factory.log_dict[logger_name])

        # Ensure configuration file was read
        mock_open.assert_called_once_with(factory.config_file, 'rt')
        mock_safe_load.assert_called_once()
        mock_dictConfig.assert_called_once_with({})

    @patch('builtins.open', new_callable=mock_open, read_data="{}")
    @patch('yaml.safe_load', return_value={})
    @patch('logging.config.dictConfig')
    def test_get_logger_caching(self, mock_dictConfig, mock_safe_load, mock_open):
        factory = LogSingletonFactory()
        logger_name = 'development'

        # Create logger
        logger1 = factory.getLogger(logger_name)
        logger2 = factory.getLogger(logger_name)

        # Ensure the same logger instance is returned
        self.assertIs(logger1, logger2)

        # Ensure the configuration file is read only once
        mock_open.assert_called_once_with(factory.config_file, 'rt')
        mock_safe_load.assert_called_once()
        mock_dictConfig.assert_called_once_with({})

    @patch('builtins.open', new_callable=mock_open, read_data="{}")
    @patch('yaml.safe_load', return_value={})
    @patch('logging.config.dictConfig')
    def test_logger_output_levels(self, mock_dictConfig, mock_safe_load, mock_open):
        factory = LogSingletonFactory()
        logger_name = 'development'

        # Create logger
        logger = factory.getLogger(logger_name)

        with self.assertLogs(logger, level='DEBUG') as log:
            logger.debug('Debug message')
            logger.info('Info message')
            logger.warning('Warning message')
            logger.error('Error message')
            logger.critical('Critical message')

        # Check the log output
        self.assertIn('DEBUG:testLogger:Debug message', log.output)
        self.assertIn('INFO:testLogger:Info message', log.output)
        self.assertIn('WARNING:testLogger:Warning message', log.output)
        self.assertIn('ERROR:testLogger:Error message', log.output)
        self.assertIn('CRITICAL:testLogger:Critical message', log.output)

def test_logger():
    unittest.main()
