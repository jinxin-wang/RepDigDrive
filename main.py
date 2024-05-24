import os
import sys
import yaml
import torch
import argparse
import trainer
# from utils import metric
from models import *
from datasets import *
from torch.utils.data import DataLoader

from datasets._PCAWG import PCAWG

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--model-config-path', dest='model_config_path', type=str, default=r"./config/train_config.yaml")
    # .....
    return parser.parse_args()

def main(args):
    print(PCAWG().config)

if __name__ == '__main__':
    args = parse_args()
    main(args)