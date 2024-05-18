#!/usr/bin/bash

LOCAL_PATH='./'

cd $LOCAL_PATH

#### download from the resources provided by paper 
wget -c -r --no-parent --progress --no-check-certificate https://cb.csail.mit.edu/cb/DIG/downloads/ 
mv ./cb.csail.mit.edu/cb/DIG .
rm -rf ./cb.csail.mit.edu/