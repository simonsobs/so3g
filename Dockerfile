# so3g
# A containerized so3g installation.

# Build on spt3g base image
FROM simonsobs/spt3g:0.3-16-g1341ea5

# Set locale
ENV LANG C.UTF-8

# Build tools needed for pixell; blas needed for so3g.
RUN apt install -y build-essential automake gfortran libopenblas-dev

# Set the working directory
WORKDIR /app_lib/so3g

# Copy the current directory contents into the container
ADD . /app_lib/so3g

# Install any needed packages specified in requirements.txt
RUN pip3 install -r requirements.txt

# Install qpoint
RUN /bin/bash /app_lib/so3g/docker/qpoint-setup.sh

# Build so3g
RUN /bin/bash /app_lib/so3g/docker/so3g-setup.sh
