#!/usr/bin/bash

set -e ;

CONDA_ENV_NAME="RepDigDrive"

# export env variables
source ~/.bashrc

conda create -n ${CONDA_ENV_NAME} "python=3.10" --yes ;
conda activate ${CONDA_ENV_NAME} ;

# Installing pytorch
# Stable (2.3.0)/Linux/CUDA 12.1
# check the install cmd for your system : https://pytorch.org/
# pip3 install torch torchvision torchaudio ;
pip3 install torch torchvision torchaudio --index-url https://download.pytorch.org/whl/cu118
# conda install --yes pytorch torchvision torchaudio pytorch-cuda=12.1 -c pytorch -c nvidia

pip3 install gpytorch ;

# Installing scikit-learn
# Linux/pip3
# check the install cmd for your system : https://scikit-learn.org/stable/install.html
pip3 install scikit-learn

# more stats packages
pip3 install numpy pandas matplotlib ;
pip3 install seaborn ;
pip3 install tqdm ;
pip3 install scipy statsmodels ;

# python packages
pip3 install PyYAML ;
pip3 install h5py ;
pip3 install pypickle ;
pip3 install tables ;
pip3 install tensorboardX ;

# bio packages :
pip3 install pyBigWig pybbi ;
pip3 install pysam pybedtools ; 

conda install --yes -c bioconda htslib ;

# jupyter
# conda install --yes -c anaconda ipykernel ;
pip3 install ipykernel
python -m ipykernel install --user --name=${CONDA_ENV_NAME} ;