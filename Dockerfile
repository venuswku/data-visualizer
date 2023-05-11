# Start from base Anaconda image.
FROM continuumio/anaconda3
# Set working directory in destination.
WORKDIR /data-visualizer
# Install the required dependencies.
# COPY environment.yml .
# RUN conda env create --file environment.yml
RUN conda install -c pyviz panel
RUN conda install -c pyviz geoviews
RUN conda install -c pyviz spatialpandas
RUN conda install -c pyviz holoviews
RUN conda install -c conda-forge dask-geopandas
RUN conda install -c conda-forge cartopy
RUN pip install geopandas
RUN pip install rioxarray
# Copy all relevant files into the container.
COPY ./ .
# Run `panel serve` to start the app.
CMD panel serve --address="0.0.0.0" --port=$PORT app.ipynb --allow-websocket-origin=geospatial-data-viewer.herokuapp.com