# algebra_ver822
- transcriptional_algebraプロジェクトの開発用。公開用はまた別に作るつもり。
- algebra_devの暖簾分け
- factor analysisをいつ使うか、データセット同士の関係性が違う。
- データをなるべくsparseに扱うというところもrenewalしたいポイント

## How to start
### Initial settings
1. clone this repository with `git clone git@github.com:yo-aka-gene/algebra_dev.git`
2. install Docker to your local env
3. swich your working dir to `/algebra-dev/`
4. run `docker compose up -d`
5. run `make init`

### Run Rython codes
`*.ipynb` files are indexed with numbers. Please follow the numbering.<br>
***Notes***: some codes require that R scripts already done. Please check the description.
1. run `docker start algebra_dev-jupyterlab-1`
2. access `localhost:8080` via your local browser (edit `docker-compose.yml` to use different port)
3. run codes in jupyter notebook environment
4. run `make lib-py` or `make lib` to install packages in dependencies (preinstalled with `make init` command)

### Run R codes
***Notes***: some packages are only installed rstuio environment. Please use `algebra_dev-rstudio-1` container.
1. run `docker start algebra_dev-rstudio-1`
2. access `localhost:8787` via your local browser (edit `docker-compose.yml` to use different port)
3. run codes in rstudio environment
4. run `make lib-r` or `make lib` to install packages in dependencies (preinstalled with `make init` command)
---
## Copyright of the dataset
1. Allen Institute for Brain Science
    - human primary mortor cortex: [https://portal.brain.map.org/atlases-and-data/rnaseq/human-m1-10x](https://portal.brain.map.org/atlases-and-data/rnaseq/human-m1-10x)
    - run `make m1_10x` or `make all-data` or `make init` to install this data.
    - For the sake of the data size, we seplit the original files (make command includes this process as well)
---
## For developers
- please use `pip` instead of `conda` for installing Python packages
- run `make write-lib` to export dependencies (both Python and R)
