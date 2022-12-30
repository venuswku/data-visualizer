# cd C:\Users\Venuxk\Projects\data-visualizer
# conda activate visualizer
# panel serve --show elwha.py

# Standard library imports

# External dependencies imports
import panel as pn
import geoviews.tile_sources as gts

# Import the data visualizer components.
from data_visualizer.components import (
	Application,
	DataMap,
	PopupModal
)
# from themes.DefaultCustomTheme import DefaultCustomTheme

# -------------------------------------------------- Constant Variables --------------------------------------------------

# Set the main color for the app.
app_main_color = "#2196f3"

# Set base path to data directories (contains category subfolders, which contain data files for each data category).
map_data_dir_path = "./data/Elwha"
time_series_data_dir_path = "./data/Elwha/Time-Series Data"

# Assign names for map's layer types.
topography_data = "Topography"
bathymetry_kayak_data = "Nearshore Bathymetry - Kayak"
bathymetry_watercraft_data = "Nearshore Bathymetry - Personal Watercraft"
grainsize_data = "Surface-Sediment Grain-Size Distributions"

elwha_data_types = [
	topography_data,
	bathymetry_kayak_data,
	bathymetry_watercraft_data,
	grainsize_data
]

data_type_colors = {
	topography_data: "red",
	bathymetry_kayak_data: "blue",
	bathymetry_watercraft_data: "green",
	grainsize_data: "#975411"
}

elwha_basemap_options = {
	"Default": gts.OSM,
	"Satellite": gts.EsriImagery,
	"Topographic": gts.OpenTopoMap,
	"Black & White": gts.StamenToner,
	"Dark": gts.CartoDark
}

all_latitude_col_names = topobathy_lat_cols = ["latitude", "Latitude"]
grainsize_lat_cols = ["Latitude (deg. N)", "Latitude"]
all_latitude_col_names.extend([col for col in grainsize_lat_cols if col not in topobathy_lat_cols])

all_longitude_col_names = topobathy_long_cols = ["longitude", "Longitude"]
grainsize_long_cols = ["Longitude (deg. E)", "Longitude"]
all_longitude_col_names.extend([col for col in grainsize_long_cols if col not in topobathy_long_cols])

all_datetime_col_names = topobathy_datetime_cols = ["Survey_Date", "datetime_utc", "Time_GMT"]
grainsize_datetime_cols = ["Date Collected"]
all_datetime_col_names.extend([col for col in grainsize_datetime_cols if col not in topobathy_datetime_cols])

all_ortho_height_col_names = ["Ortho_Ht_m", "Ortho_ht_m", "ortho_ht_m"]

all_weight_col_names = ["Wt. percent in -2.00 phi bin"]

# -------------------------------------------------- Initializing Data Visualization App --------------------------------------------------

# Instantiate the app's template.
template = pn.template.BootstrapTemplate(
	site = "Data Visualizer",
    title = "Elwha Topo-Bathy Data",
    header_background = app_main_color
	# theme = DefaultCustomTheme
)

# Instantiate the main components required by the Application.
data_map = DataMap(
	data_dir_path = map_data_dir_path,
  	latitude_col_names = all_latitude_col_names,
  	longitude_col_names = all_longitude_col_names,
	# colors = data_type_colors,
	basemap_options = elwha_basemap_options
)
popup_modal = PopupModal(
	data_converter = data_map,
	template = template,
	time_series_data_col_names = all_ortho_height_col_names + all_weight_col_names
)

# Create the application.
app = Application(
	data_map = data_map,
	popup_modal = popup_modal
)

# Populate the template with the sidebar, main, and modal layout.
template.sidebar.extend([
	*(data_map.param_widgets)
])
template.main.append(pn.panel(data_map.plot, loading_indicator = True))
template.modal.extend([
	pn.panel(popup_modal.content, loading_indicator = True)
])

# Use the Panel extension to load BokehJS, any pn.config variables, any custom models required, or optionally additional custom JS and CSS in Jupyter notebook environments.
pn.extension(loading_spinner = "dots", loading_color = app_main_color, sizing_mode = "stretch_width")

# Launch the app (`panel serve --show --autoreload elwha.py`).
template.servable()