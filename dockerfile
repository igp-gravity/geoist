FROM continuumio/miniconda3
COPY ./examples /home/
COPY ./notebooks /home/
COPY ./environment.yml /home/
WORKDIR /home
RUN conda env create -f environment.yml
RUN echo "source activate $(head -1 /home/environment.yml | cut -d' ' -f2)" > ~/.bashrc
ENV PATH /opt/conda/envs/$(head -1 /home/environment.yml | cut -d' ' -f2)/bin:$PATH
RUN conda activate $(head -1 /home/environment.yml | cut -d' ' -f2) && conda install gitlab && python geoist/setup.py
CMD jupyter lab --no-browser
