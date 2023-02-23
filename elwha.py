# cd C:\Users\Venuxk\Projects\data-visualizer
# conda activate visualizer
# panel serve --show elwha.py

# Standard library imports

# External dependencies imports
import panel as pn

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

# -------------------------------------------------- Initializing Data Visualization App --------------------------------------------------

# Instantiate the app's template.
template = pn.template.BootstrapTemplate(
	site = "Data Visualizer",
    title = "Elwha Topo-Bathy Data",
    header_background = app_main_color
	# theme = DefaultCustomTheme
)

# Instantiate the main components required by the Application.
data_map = DataMap()
popup_modal = PopupModal(
	data_map = data_map,
	template = template,
	time_series_data_col_names = ["Ortho_Ht_m", "Ortho_ht_m", "ortho_ht_m", "F-W Mean"]
)

# Create the application.
app = Application(
	data_map = data_map,
	popup_modal = popup_modal
)

# Populate the template with the sidebar, main, and modal layout.
template.sidebar.extend([
    # app.wiki_info_button,
	*(data_map.param_widgets),
    popup_modal.time_series_controls
])
template.main.append(pn.panel(data_map.plot, loading_indicator = True))
template.modal.extend([
	pn.panel(popup_modal.content, loading_indicator = True)
])

# Use the Panel extension to load BokehJS, any pn.config variables, any custom models required, or optionally additional custom JS and CSS in Jupyter notebook environments.
pn.extension(loading_spinner = "dots", loading_color = app_main_color, sizing_mode = "stretch_width")

# Launch the app (`panel serve --show --autoreload elwha.py`).
template.servable()