FROM jupyter/datascience-notebook:lab-3.4.5

RUN pip install --upgrade pip
RUN pip install jupyterlab
RUN jupyter serverextension enable --py jupyterlab

WORKDIR $HOME
RUN mkdir code tools data out

EXPOSE 8888
VOLUME ["/home/jovyan/code"]

FROM rocker/tidyverse:4

RUN apt-get update
RUN apt-get install -y \
libcurl4-openssl-dev \
libssl-dev libjq-dev \
libprotobuf-dev \
protobuf-compiler \
make \
libgeos-dev \
libglpk40 \
libudunits2-dev \
libgdal-dev \
gdal-bin \
libproj-dev \
libv8-dev

CMD ['/bin/bash']
