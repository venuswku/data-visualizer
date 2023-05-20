# Start from the base Miniconda image.
FROM continuumio/miniconda3:latest

# Set working directory in destination.
WORKDIR /data-visualizer

# Use faster solver for Miniconda.
# RUN conda install -n base conda-libmamba-solver
# RUN conda config --set solver libmamba

# Install the required dependencies.
# ^ failing environment has spatialpandas=0.4.7 holoviews=1.16.0 dask-geopandas=0.3.1 geopandas=0.13.0 rioxarray=0.14.1 jupyterlab=4.0.0
RUN conda install --channel conda-forge panel=0.14.0
RUN conda install --channel conda-forge geoviews=1.9.6
RUN conda install --channel conda-forge spatialpandas=0.4.6
RUN conda install --channel conda-forge holoviews=1.15.4
RUN conda install --channel conda-forge dask-geopandas=0.3.0
RUN conda install --channel conda-forge geopandas=0.12.2
RUN conda install --channel conda-forge rioxarray=0.14.0
RUN conda install --channel conda-forge cartopy=0.21.1
RUN conda install --channel conda-forge jupyterlab=3.5.0
RUN conda install --channel bokeh bokeh=2.4.3

# Copy all relevant files into the container.
COPY ./ .

# Run `panel serve` to start the app.
# CMD panel serve --address="0.0.0.0" --port=$PORT app.ipynb --allow-websocket-origin=geospatial-data-viewer.herokuapp.com
CMD panel serve app.ipynb --allow-websocket-origin=localhost:8000