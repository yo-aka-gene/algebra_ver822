#!/bin/bash

docker exec $1 apt-get update
docker exec $1 apt-get install -y \
libcurl4-openssl-dev libssl-dev libjq-dev \
libprotobuf-dev protobuf-compiler make \
libgeos-dev libglpk40 libudunits2-dev \
libgdal-dev gdal-bin libproj-dev libv8-dev
