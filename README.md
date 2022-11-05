# Data Visualizer
## Create an Environment with Anaconda
Create a new environment named `visualizer` with all the required packages by entering the following commands in succession into Anaconda Prompt (Windows) or Terminal (Mac/Linux):
```
conda create -n visualizer
conda activate visualizer
# Install Panel dependencies.
conda install -c bokeh ipywidgets_bokeh -y
conda install -c conda-forge panel -y
# Install other dependencies.
conda install -c conda-forge jupyterlab geoviews rioxarray pandas -y
```

## Launch Jupyter Notebook as a Web Server
- Make sure your Anaconda environment is activated by running `conda activate visualizer` in your terminal.
- Run the command `panel serve --show --autoreload app.ipynb` in your terminal.
- A webpage with the URL http://localhost:5006/app will display all Panel objects marked with `.servable()`.
- Any changes in the notebook will automatically be reflected on the webpage. Just in case, refresh the webpage to make sure you see your latest changes.

## Launch Jupyter Notebook
- Make sure your Anaconda environment is activated by running `conda activate visualizer` in your terminal.
- Run the command `jupyter notebook` in your terminal.
- Open the `app.ipynb` file when a webpage with the URL http://localhost:8888/tree appears.
  - Directly running `jupyter notebook app.ipynb` will skip this step of selecting a notebook to open.
- Run all the notebook cells from top to bottom. The Panel app will be outputted after the last cell is run.
- Reload the [`app.ipynb` webpage](http://localhost:8888/notebooks/app.ipynb) when you want to see your new changes.