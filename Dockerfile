# so3g
# A containerized so3g installation.

# Build on spt3g base image
FROM simonsobs/spt3g:0.3-289-g4bd3275

# Set locale
ENV LANG C.UTF-8

# Build tools needed for pixell; blas needed for so3g.
RUN apt update && apt install -y \
    build-essential \
    automake \
    gfortran \
    libopenblas-openmp-dev \
    libbz2-dev \
    python-is-python3 \
    libfftw3-dev
    libgoogle-glog-dev \
    libgflags-dev \
    libmetis-dev \
    libgtest-dev \
    libabsl-dev \
    libeigen3-dev

# Set the working directory
WORKDIR /app_lib/so3g

# Fetch and install ceres-solver
RUN git clone --depth 1 --branch 2.2.0 --recurse-submodules https://github.com/ceres-solver/ceres-solver

WORKDIR /app_lib/so3g/ceres-solver

RUN mkdir build \
    && cd build \
    && cmake .. -DBUILD_TESTING=OFF \
    && make -j$(nproc) \
    && make install

# Set the working directory back to so3g
WORKDIR /app_lib/so3g

# Copy the current directory contents into the container
ADD . /app_lib/so3g

# Install any needed packages specified in requirements.txt
RUN pip3 install -r requirements.txt
RUN pip3 install -r test-requirements.txt

# Build so3g
RUN /bin/bash /app_lib/so3g/docker/so3g-setup.sh
