FROM ubuntu:20.04

MAINTAINER yuhao<yqyuhao@outlook.com>

RUN sed -i 's/http:\/\/archive\.ubuntu\.com\/ubuntu\//http:\/\/mirrors\.aliyun\.com\/ubuntu\//g' /etc/apt/sources.list

# set timezone
RUN set -x \
&& export DEBIAN_FRONTEND=noninteractive \
&& apt-get update \
&& apt-get install -y tzdata \
&& ln -sf /usr/share/zoneinfo/Asia/Shanghai  /etc/localtime \
&& echo "Asia/Shanghai" > /etc/timezone

# install packages
RUN apt-get update \

&& apt-get install -y less curl apt-utils vim wget gcc-7 g++-7 make cmake git unzip dos2unix libncurses5 \

# lib
&& apt-get install -y zlib1g-dev libjpeg-dev libncurses5-dev libbz2-dev liblzma-dev libcurl4-gnutls-dev \
 
# python3 perl java r-base
&& apt-get install -y python3 python3-dev python3-pip python perl openjdk-8-jdk r-base r-base-dev

ENV software /yqyuhao

# create software folder

RUN mkdir -p $software/database $software/bin /data /data/analysis 

# install R packages
WORKDIR $software/source
RUN Rscript -e "install.packages(c('BiocManager','dplyr','readxl','flowCore','cowplot','ggplot2'));BiocManager::install('Seurat')"

# copy esssential files
#COPY flowcyto.R $software/bin/

# chown root:root
WORKDIR $software
RUN chown root:root -R $software

# mkdir fastq directory and analysis directory
WORKDIR /data/analysis
