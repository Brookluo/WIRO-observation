FROM python:3.10.6-bullseye

RUN apt -y update && apt install -y apt-utils

RUN apt install -y --no-install-recommends \
    gcc \
    g++ \
    make \
    git \
    file \
    pkg-config \
    wget \
    swig \
    netpbm \
    wcslib-dev \
    wcslib-tools \
    zlib1g-dev \
    libbz2-dev \
    libcairo2-dev \
    libcfitsio-dev \
    libcfitsio-bin \
    libgsl-dev \
    libjpeg-dev \
    libnetpbm10-dev \
    libpng-dev \
    libeigen3-dev \
    libgoogle-glog-dev \
    libceres-dev \
    postgresql-common \
    libpq-dev \
    curl \
    # Remove APT files
    && apt-get clean 


# Pip installs
RUN for pkg in \
    setuptools \
    wheel \
    cython \
    numpy \
    scipy \
    pillow \
    psycopg2 \
    fitsio \
    matplotlib \
    astropy \
    photutils \
    astroplan \
    ccdproc \
    ; do pip install $pkg; done

RUN mkdir /src
WORKDIR /src
RUN git clone https://github.com/Brookluo/WIRO-observation.git


ENTRYPOINT [ "/bin/bash" ]