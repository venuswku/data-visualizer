# # Start from the base Miniconda image.
FROM continuumio/miniconda3
# FROM osgeo/gdal

# Set working directory in destination.
WORKDIR /data-visualizer

# Install the required dependencies.
RUN conda install --channel conda-forge panel
RUN conda install --channel conda-forge geoviews
RUN conda install --channel conda-forge spatialpandas
RUN conda install --channel conda-forge holoviews
RUN conda install --channel conda-forge dask-geopandas
RUN conda install --channel conda-forge cartopy
RUN conda install --channel conda-forge geopandas
RUN conda install --channel conda-forge rioxarray

# Copy all relevant files into the container.
COPY ./ .

# Run `panel serve` to start the app.
CMD panel serve --address="0.0.0.0" --port=$PORT app.ipynb --allow-websocket-origin=geospatial-data-viewer.herokuapp.com