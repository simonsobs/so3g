# so3g-jupyter
# Jupyter notebook environment with so3g installed.

# Use so3g base image
FROM so3g:latest

# Install python dependencies
COPY requirements.txt /tmp/
RUN pip3 install -r /tmp/requirements.txt

# Adds jupyter autocomplete extension
RUN jupyter contrib nbextension install \
    && jupyter nbextension enable hinterland/hinterland

# Adds OCS support
WORKDIR /app_lib/
RUN git clone https://github.com/simonsobs/ocs.git \
    && pip3 install -r ocs/requirements.txt \
    && pip3 install -e ocs

ENV OCS_CONFIG_DIR /config/
ENV JUPYTER_CONFIG_DIR /config/jupyter/

ENTRYPOINT ["jupyter", "notebook"]

# By default puts you into $OCS_CONFIG_DIR
CMD ["/config/"]
