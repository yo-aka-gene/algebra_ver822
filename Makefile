.PHONY: clean clean-build clean-pyc all-data gse m1_10x lib lib-py lib-r write-lib init fmt_10x
.DEFAULT_GOAL := help

clean: clean-build clean-pyc

clean-build: ## remove build artifacts
	rm -fr build/
	rm -fr dist/
	rm -fr .eggs/
	find . -name '*.egg-info' -exec rm -fr {} +
	find . -name '*.egg' -exec rm -f {} +

clean-pyc: ## remove Python file artifacts
	find . -name '*.pyc' -exec rm -f {} +
	find . -name '*.pyo' -exec rm -f {} +
	find . -name '*~' -exec rm -f {} +
	find . -name '__pycache__' -exec rm -fr {} +

all-data: gse m1_10x

gse : ## download gse165388 data
	mkdir ./data/gse165388
	curl -o ./data/GSE165388_RAW.tar 'https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE165388&format=file'
	tar -xvf ./data/GSE165388_RAW.tar -C ./data/gse165388
	rm ./data/GSE165388_RAW.tar
	chmod 741 fmt_10x.sh
	./fmt_10x.sh ./data/gse165388


m1_10x: ## download m1_10x data and split them
	mkdir ./data/m1_10x
	curl -o ./data/m1_10x/matrix.csv https://idk-etl-prod-download-bucket.s3.amazonaws.com/aibs_human_m1_10x/matrix.csv
	split ./data/m1_10x/matrix.csv -l 2000
	curl -o ./data/m1_10x_meta/metadata.csv https://idk-etl-prod-download-bucket.s3.amazonaws.com/aibs_human_m1_10x/metadata.csv
	curl -o ./data/m1_10x_info/readme.txt https://brainmapportal-live-4cc80a57cd6e400d854-f7fdcae.divio-media.net/filer_public/70/32/70326830-e306-4743-a02c-a8da5bf9eb56/readme-m1-10.txt


lib: write-lib lib-py lib-r write-lib

lib-py: ## install required packages in Python
	docker exec algebra_ver822-jupyterlab-1 python -m pip install -r ./tools/requirements_py.txt

lib-r: ## install required packages in R
	docker exec  algebra_ver822-rstudio-1 Rscript ./home/rstudio/tools/install_deps.R

write-lib: ## export reqired packages
	docker exec algebra_ver822-jupyterlab-1 pip list --format=freeze > ./tools/requirements_py.txt
	docker exec  algebra_ver822-rstudio-1 Rscript ./home/rstudio/tools/export_deps.R

init:
	make all-data
	docker exec algebra_dev-jupyterlab-1 touch ./tools/requirements.txt
	make lib

