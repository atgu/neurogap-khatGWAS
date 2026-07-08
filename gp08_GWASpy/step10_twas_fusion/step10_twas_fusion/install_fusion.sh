#!/bin/bash
set -euxo pipefail

export DEBIAN_FRONTEND=noninteractive

apt-get update
apt-get install -y \
  unzip \
  wget \
  git \
  r-base \
  r-base-dev \
  build-essential \
  libcurl4-openssl-dev \
  libssl-dev \
  libxml2-dev

mkdir -p /opt/fusion_setup
cd /opt/fusion_setup

# FUSION
wget -O fusion_master.zip https://github.com/gusevlab/fusion_twas/archive/master.zip
unzip -o fusion_master.zip

# Make a stable path for the code
ln -sfn /opt/fusion_setup/fusion_twas-master /opt/fusion_twas

# LD reference data
wget -O LDREF.tar.bz2 https://data.broadinstitute.org/alkesgroup/FUSION/LDREF.tar.bz2
tar xjvf LDREF.tar.bz2

# plink2R
wget -O plink2R_master.zip https://github.com/gabraham/plink2R/archive/master.zip
unzip -o plink2R_master.zip

# R packages
Rscript -e "install.packages(c('optparse','RColorBrewer'), repos='https://cloud.r-project.org')"
Rscript -e "install.packages('/opt/fusion_setup/plink2R-master/plink2R', repos=NULL, type='source')"
