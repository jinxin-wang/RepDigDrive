import os
import sys
import yaml
import torch
import argparse
import trainer
from utils import metric
from models import *
from datasets import *
from torch.utils.data import DataLoader

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--model-config-path', dest='model_config_path', type=str, default=r"./config/train_config.yaml")
    # .....
    return parser.parse_args()

def load_config(yaml_path):
    with open(yaml_path, 'r', encoding='utf-8') as fd:
        config = yaml.safe_load(fd)
    return config

def main(args):
    print('Hello World!')
    model_config = load_config(args.model_config_path)
    print(model_config)

if __name__ == '__main__':
    args = parse_args()
    main(args)