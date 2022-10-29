# Standard library imports
import os
import datetime as dt

# External dependencies imports
import panel as pn
from ipyleaflet import basemaps
from ipywidgets import Button

# Import the data visualizer components.
from data_visualizer.components import (
    Application,
	DataMap
)

# Use the Panel extension to load BokehJS, any pn.config variables, any custom models required, or optionally additional custom JS and CSS in Jupyter notebook environments.
pn.extension(loading_spinner = "dots", loading_color = app_main_color, sizing_mode = "stretch_width")

# -------------------------------------------------- Constant Variables --------------------------------------------------

# Set the main color for the app.
app_main_color = "#2196f3"

# Set base path to data directories.
data_dir_path = "./data/Elwha"

# Assign names for map's layer types.
topography_data = "Topography"
bathymetry_kayak_data = "Nearshore Bathymetry - Kayak"
bathymetry_watercraft_data = "Nearshore Bathymetry - Personal Watercraft"
grainsize_data = "Surface-Sediment Grain-Size Distributions"
basemap_data = "Basemap"

elwha_data_types = [
  topography_data,
  bathymetry_kayak_data,
  bathymetry_watercraft_data,
  grainsize_data
]

# For all data types, get optional styling information ahead of time in order for GeoJSON data layers to appear on map.
data_type_styles = {}
for data_type in elwha_data_types:
  data_type_styles[data_type] = {}
  # Not assigning styles will keep ipyleaflet's default feature styling (marker for points).
  data_type_styles[data_type]["point_style"] = get_data_type_point_style(data_type)
  data_type_styles[data_type]["hover_style"] = get_data_type_hover_style(data_type)

data_type_colors = {
  topography_data: "red",
  bathymetry_kayak_data: "blue",
  bathymetry_watercraft_data: "green",
  grainsize_data: "#975411"
}

elwha_basemap_options = {
  "Default": basemaps.OpenStreetMap.Mapnik,
  "Satellite": basemaps.Esri.WorldImagery,
  "Topographic": basemaps.OpenTopoMap,
  "Black & White": basemaps.Stamen.Toner,
  "Dark": basemaps.CartoDB.DarkMatter
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

# -------------------------------------------------- Helper Functions --------------------------------------------------

# Gets the point_style based on the data type.
def get_data_type_point_style(data_type):
  color = data_type_colors[data_type]
  return {"color": color, "opacity": 0.5, "fillColor": color, "fillOpacity": 0.3, "radius": 8, "weight": 1, "dashArray": 2}

# Gets the hover_style based on the data type.
def get_data_type_hover_style(data_type):
  return {"color": app_main_color, "fillColor": app_main_color, "weight": 3}

# Gets info about the data file and creates a new GeoJSON layer with it.
def create_layer(file, data_type):
  # print("Loading data from " + file + "...")
  # Determine popup content based on different types of data.
  popup_info = {}
  if data_type == grainsize_data:
    popup_info = {
      "Date & Time Collected": [
        "Date Collected",
        [" "],
        {
          "Time (GMT)": "GMT",      # for grainsize data before July 2018
          "Time_GMT": "GMT"         # for grainsize data at and after July 2018
        }
      ],
      "Sample Type": ["Sample Type"],
      "Weight": ["Wt. percent in -2.00 phi bin", ["%"]],
      "Gravel": ["Percent Gravel", ["%"]],
      "Sand": ["Percent Sand", ["%"]],
      "Silt": ["Percent Silt", ["%"]],
      "Clay": ["Percent Clay", ["%"]],
      "Mud": ["Percent Mud", ["%"]]
    }
  elif (data_type == topography_data) or (data_type == bathymetry_kayak_data) or (data_type == bathymetry_watercraft_data):
    popup_info = {
      "Date & Time Collected": [
        {
          "Survey_Date": "",        # for topo-bathy data before July 2018
          "datetime_utc": "UTC"     # for topo-bathy data at and after July 2018
        }
      ],
      "Orthometric Height": [
        {
          "Ortho_Ht_m": "meters",
          "Ortho_ht_m": "meters",
          "ortho_ht_m": "meters"
        },
      ]
    }
  # Create and display GeoJSON layer on map.
  elwha.create_geojson(
    data_path = data_dir_path + "/" + data_type + "/" + file,
    name = file,
    popup_content = popup_info,
    longitude_col_names = all_longitude_col_names,
    latitude_col_names = all_latitude_col_names
  )

# Checks if the file contains data from the user's selected date range.
def data_within_date_range(filename):
  (selected_start_date, selected_end_date) = data_date_range_slider.value
  
  # Get the data's month and year from its file name.
  month_num = {"jan": 1, "feb": 2, "mar": 3, "apr": 4, "may": 5, "june": 6, "july": 7, "aug": 8, "sept": 9, "oct": 10, "nov": 11, "dec": 12}
  [month_name] = filter(lambda m: m in filename, month_num.keys())
  month = month_num[month_name]
  year = 2000 + int("".join(char for char in filename if char.isdigit()))
  file_date = dt.datetime(year, month, 1)
  
  return selected_start_date <= file_date <= selected_end_date

# -------------------------------------------------- Elwha Topo-Bathy Data Widgets --------------------------------------------------

# Instantiate the main components required by the Application.
map = DataMap(
	data_dir_path = data_dir_path,
	map_center = (48.148, -123.553),
	# category_styles = data_type_styles,
	data_details_button = see_data_point_details_button,
	basemap_options = elwha_basemap_options,
	legend_name = "Types of Data"
)
basemap_select = pn.widgets.Select(name="Basemap", options=list(elwha_basemap_options.keys()))
elwha_data_type_multi_choice = pn.widgets.MultiChoice(name="Type of Data", options=elwha_data_types, placeholder="Choose one or more types of data to display", solid=False)
data_date_range_slider = pn.widgets.DateRangeSlider(
	name = "Data Collection Range",
	start = dt.datetime(2010, 9, 5), end = dt.datetime.utcnow(),
	value = (dt.datetime(2018, 1, 1), dt.datetime(2019, 1, 1)),
	bar_color = app_main_color
)
see_data_point_details_button = Button(
	description = "See Details",
	tooltip = "Displays a time-series of data previously collected at this location. Because large datasets overcrowd the map, this data point is sampled from the original dataset. Click the button to also view the original dataset.",
	icon = "chart-line",
	button_style = "primary",
	style = dict(button_color = app_main_color)
)

# -------------------------------------------------- Callbacks & Reactive Functions --------------------------------------------------

# Update basemap whenever a different basemap value is selected.
basemap_select.param.watch(elwha.update_basemap, "value")

# Opens a modal containing a time-series plot and another plot with the original dataset that the selected data point was sampled from.
def display_data_point_details(event):
  # Open the app modal to display the scatter plots.
  elwha_data_visualizer.open_modal()

  # Show loading spinner while the data point's scatter plots are being created.
  elwha.plotter.plot_data_point_details(
    data = elwha.selected_geojson_data,
    category_latitude_cols = {
      topography_data: topobathy_lat_cols,
      bathymetry_kayak_data: topobathy_lat_cols,
      bathymetry_watercraft_data: topobathy_lat_cols,
      grainsize_data: grainsize_lat_cols
    },
    category_longitude_cols = {
      topography_data: topobathy_long_cols,
      bathymetry_kayak_data: topobathy_long_cols,
      bathymetry_watercraft_data: topobathy_long_cols,
      grainsize_data: grainsize_long_cols
    },
    category_datetime_cols = {
      topography_data: topobathy_datetime_cols,
      bathymetry_kayak_data: topobathy_datetime_cols,
      bathymetry_watercraft_data: topobathy_datetime_cols,
      grainsize_data: grainsize_datetime_cols
    },
    category_y_axis_cols = {
      topography_data: all_ortho_height_col_names,
      bathymetry_kayak_data: all_ortho_height_col_names,
      bathymetry_watercraft_data: all_ortho_height_col_names,
      grainsize_data: all_weight_col_names
    },
    category_y_axis_label = {
      topography_data: "Orthometric Height (meters)",
      bathymetry_kayak_data: "Orthometric Height (meters)",
      bathymetry_watercraft_data: "Orthometric Height (meters)",
      grainsize_data: "Weight Percentage in -2.00 phi bin"
    }
  )

# Display scatter plots in a modal whenever the user clicks on the button for viewing how a dataset changes over time.
see_data_point_details_button.on_click(display_data_point_details)

# Filters data based on what data type(s) and date range that the user selects.
def filter_data_on_map(event):
  selected_data_types = elwha_data_type_multi_choice.value
  for data_type in elwha_data_types:
    data_type_files = os.listdir(data_dir_path + "/" + data_type)
    for file in data_type_files:
      if (data_type in selected_data_types) and data_within_date_range(file):
        # Create and display the selected data if we never read the file before.
        if file not in elwha.geojsons:
          # print("create", file)
          create_layer(file, data_type)
        # Display the selected data if it isn't in map yet.
        else:
          # print("display", file)
          elwha.display_geojson(file)
      # Else hide the data if user didn't select to display it.
      else:
        # print("hide", file)
        elwha.hide_geojson(file)

# Filter data whenever the selected data type(s) or date range change.
elwha_data_type_multi_choice.param.watch(filter_data_on_map, "value")
data_date_range_slider.param.watch(filter_data_on_map, "value")

# -------------------------------------------------- Initializing Data Visualization App --------------------------------------------------

# Create the application.
app = Application(map = map)

# Populate the template with the sidebar and the main layout.
template = pn.template.BootstrapTemplate(
	site = "Data Visualizer",
    title = "Elwha Topo-Bathy Data",
    header_background = app_main_color,
    sidebar = [
		basemap_select,
		elwha_data_type_multi_choice,
		data_date_range_slider
	],
    main = [],
	modal = [

	]
)

# Launch the app (`panel serve --show --autoreload elwha.py`).
template.servable()