from config import DatasetConfig

ds_config = DatasetConfig()
print(ds_config.getConfigAttr('root')['download'])

print(ds_config.getDatasetNames())