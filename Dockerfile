# Use python 3.7 
FROM python:3.7-slim

# install this way to fix paths in coverage report
ENV PYTHONPATH=$PYTHONPATH:/code/barycorrpy
ENV PYTHONHASHSEED=0

# change me to point to your coveralls view of the repo
ENV COVERALLS_REPO_TOKEN=pU74HZ9rvqc0u7iXr95LbzQ8NTzycWSLx

# turn off built-in Python multithreading
ENV MKL_NUM_THREADS=1
ENV NUMEXPR_NUM_THREADS=1
ENV OMP_NUM_THREADS=1

# setup the working directory
RUN mkdir /code && \
    mkdir /code/barycorrpy && \
    apt-get --yes update && \
    apt install build-essential -y --no-install-recommends && \
    apt-get install --yes git vim emacs && \
    /usr/local/bin/python -m pip install --upgrade pip && \
    cd /code

# Set the working directory to KPF-Pipeline
WORKDIR /code/barycorrpy
