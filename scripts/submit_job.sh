#!/usr/bin/bash

#SBATCH --job-name=wes_vae
#SBATCH --nodes=1
#SBATCH --mem=30gb
#SBATCH --time=7-00:00:00
#SBATCH --cpus-per-task=2
#SBATCH --partition=gpgpuq
#SBATCH --ntasks-per-node=3
#SBATCH --gres=gpu:a100:3
#SBATCH --error=./logs/wes_vae.err
#SBATCH --output=./logs/wes_vae.output

#### Notice: 
#### 1. types of gpu on flamingo
#### a100: gpgpuq
#### v100: gpgpuq
#### t4:   visual
#### 
#### partitions
#### gpgpuq       up 7-00:00:00      1    mix gpu03
#### gpgpuq       up 7-00:00:00      1   idle gpu02
#### visuq        up 7-00:00:00      1    mix gpu01
#### inter        up   infinite      1    mix gpu03
#### inter        up   infinite      1   idle gpu02
#### trainings    up   infinite      1    mix n25


#### Notice: 
#### 2. Number of GPUs
#### --ntasks-per-node=[num] should be equal to --gres=gpu:[num]
#### the number should be set in trainer as follow :
#### trainer = L.Trainer(
####    accelerator='gpu',
####    devices=[0,...,num], or device=[num],
####    ...
#### )

set -e ; 
trap 'exit' INT ; 

#### Notice : 
#### 3. conda environment
#### it needs to load .bashrc before activate the conda env 
source ~/.bashrc ; 
conda activate lightning ;

#### enable cuda debug
# CUDA_LAUNCH_BLOCKING=1

#### Notice :
#### 4. to start the training in SLURM env 
#### It needs to use srun to start 
srun python train.py ;
