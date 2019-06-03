# so3g
# A containerized so3g installation.

# Build on spt3g base image
FROM grumpy.physics.yale.edu/spt3g:ee03194

# Set the working directory
WORKDIR /app_lib/so3g

# Copy the current directory contents into the container
ADD . /app_lib/so3g

# Install any needed packages specified in requirements.txt
RUN pip3 install -r requirements.txt

# Build so3g
RUN /bin/bash /app_lib/so3g/docker/so3g-setup.sh
