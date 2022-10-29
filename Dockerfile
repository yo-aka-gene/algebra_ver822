FROM jupyter/datascience-notebook:lab-3.4.5

RUN pip install --upgrade pip
RUN pip install jupyterlab
RUN jupyter serverextension enable --py jupyterlab

WORKDIR $HOME
RUN mkdir code tools data out

EXPOSE 8888
VOLUME ["/home/jovyan/code"]

CMD ['/bin/bash']
