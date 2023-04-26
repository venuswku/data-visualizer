# Standard library imports

# External dependencies imports
import panel as pn
import param
from .DataMap import DataMap
from .PopupModal import PopupModal

class Application(param.Parameterized):
    # -------------------------------------------------- Main Components --------------------------------------------------
    data_map = param.ClassSelector(class_ = DataMap, is_instance = True)
    popup_modal = param.ClassSelector(class_ = PopupModal, is_instance = True)

    # -------------------------------------------------- Parameters --------------------------------------------------

    # -------------------------------------------------- Constructor --------------------------------------------------
    def __init__(self, template: pn.template, **params) -> None:
        """
        Creates a new instance of the Application class with its instance variables.

        Args:
            template (panel.template): Data visualizer app's template
        """
        super().__init__(**params)

        # -------------------------------------------------- Constants --------------------------------------------------
        self._app_template = template

        # _wiki_info_button = button that opens a tab to the GitHub Wiki page of the Data Visualizer app
        self._wiki_info_button = pn.widgets.Button(name = "\u2139", button_type = "light", width = 30)
        self._wiki_info_button.js_on_click(
            args = {"wiki_url": "https://github.com/venuswku/data-visualizer/wiki"},
            code = "window.open(wiki_url)"
        )
        
        # -------------------------------------------------- Internal Class Properties --------------------------------------------------
        # _sidebar_accordion = sidebar widget with expandable and collapsable sections that contain controls and information
        self._sidebar_accordion = pn.Accordion(objects = [], active = [], toggle = True, sizing_mode = "stretch_width")
        
    # -------------------------------------------------- Private Class Methods --------------------------------------------------
    @param.depends("data_map.clicked_transects_info", watch = True)
    def _update_clicked_transects_info(self) -> None:
        """
        Updates PopupModal's clicked_transects_info parameter, which stores information about the most recently clicked transect(s) or user-drawn transect from the data map,
        whenever DataMap's clicked_transects_info parameter changes because the user wants to view information about a different transect.
        """
        if self.data_map.clicked_transects_info:
            self.popup_modal.clicked_transects_info = self.data_map.clicked_transects_info
            # Open the app's modal to display info/error message about the selected transect(s).
            self._app_template.open_modal()

    @param.depends("data_map.collection", watch = True)
    def _update_time_series_collection_path(self) -> None:
        """
        Triggers event to update PopupModal's _collection_dir_path internal property and its related objects if DataMap's collection parameter changed.
        """
        self.popup_modal.update_collection_dir_path = True
    
    @param.depends("popup_modal.displayed_data_file", watch = True)
    def _update_selected_data_file(self) -> None:
        """
        Updates DataMap's data_file_paths parameter with one of the selected data files highlighted in PopupModal's MultiSelect widgets.
        """
        if self.popup_modal.displayed_data_file is not None:
            self.data_map.data_file_paths = [self.popup_modal.displayed_data_file]
        else:
            self.data_map.data_file_paths = []
    
    # -------------------------------------------------- Public Class Properties & Methods --------------------------------------------------
    @param.depends("data_map.collection")
    def sidebar_widgets(self) -> pn.Column:
        """
        Returns a column of widgets from both DataMap and PopupModal to display in the app's sidebar.
        """
        data_map_widgets = self.data_map.get_sidebar_widgets()
        poupup_modal_widgets = self.popup_modal.get_sidebar_widgets()
        return pn.Column(
            objects = data_map_widgets + poupup_modal_widgets
        )

    @param.depends("data_map.collection", "data_map.update_accordion_section", "popup_modal.update_accordion_section")
    def sidebar_accordion(self) -> pn.Accordion:
        """
        Returns an accordion widget with updated contents from both DataMap and PopupModal.
        """
        data_map_sections = self.data_map.get_accordion_sections()
        poupup_modal_sections = self.popup_modal.get_accordion_sections()
        new_accordion = pn.Accordion(
            objects = data_map_sections + poupup_modal_sections,
            active = self._sidebar_accordion.active,
            toggle = True, sizing_mode = "stretch_width"
        )
        self._sidebar_accordion = new_accordion
        return new_accordion

    @property
    def wiki_info_button(self) -> pn.widgets.Button:
        """
        Returns the button widget that opens a tab to the GitHub Wiki page of the Data Visualizer app.
        """
        return self._wiki_info_button