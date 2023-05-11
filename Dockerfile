# Start from base Anaconda image.
FROM continuumio/anaconda3
# Set working directory in destination.
WORKDIR /data-visualizer
# Install the required dependencies.
COPY environment.yml .
RUN conda env create --file environment.yml
# Copy all relevant files into the container.
COPY . .
# Run `panel serve` to start the app.
CMD panel serve --address="0.0.0.0" --port=$PORT data_viewer_app.ipynb --allow-websocket-origin=geospatial-data-viewer.herokuapp.com