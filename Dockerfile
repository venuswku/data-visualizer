# Start from base Anaconda image.
FROM continuumio/anaconda3
# Set working directory in destination.
WORKDIR /data-visualizer
# Install the required dependencies.
RUN conda install -c conda-forge panel
RUN conda install -c conda-forge geoviews
RUN conda install -c conda-forge spatialpandas
RUN conda install -c conda-forge holoviews
RUN conda install -c conda-forge dask-geopandas
RUN conda install -c conda-forge cartopy
RUN conda install -c conda-forge rioxarray
# Copy all relevant files into the container.
COPY . .
# Run `panel serve` to start the app.
CMD panel serve --address="0.0.0.0" --port=$PORT data_viewer_app.ipynb --allow-websocket-origin=geospatial-data-viewer.herokuapp.com