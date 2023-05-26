# Start from the base Miniconda image.
FROM continuumio/miniconda3:latest

# Set working directory in destination.
WORKDIR /data-visualizer

# Use faster solver for Miniconda.
# RUN conda install -n base conda-libmamba-solver
# RUN conda config --set solver libmamba

# Install the required dependencies.
RUN conda install --channel conda-forge panel
RUN conda install --channel conda-forge geoviews
RUN conda install --channel conda-forge spatialpandas
RUN conda install --channel conda-forge "holoviews>=1.16.0"
RUN conda install --channel conda-forge "bokeh>=3.1.0"
RUN conda install --channel conda-forge dask-geopandas
RUN conda install --channel conda-forge geopandas
RUN conda install --channel conda-forge rioxarray
RUN conda install --channel conda-forge cartopy
RUN conda install --channel conda-forge jupyterlab

# Copy all relevant files into the container.
COPY ./ .

# Run `panel serve` to start the app.
# CMD panel serve --address="0.0.0.0" --port=$PORT app.ipynb --allow-websocket-origin=geospatial-data-viewer.herokuapp.com
CMD panel serve app.ipynb --show --allow-websocket-origin=localhost:8000