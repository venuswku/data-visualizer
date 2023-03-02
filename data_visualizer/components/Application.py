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
    def __init__(self, **params) -> None:
        """
        Creates a new instance of the Application class with its instance variables.
        """
        super().__init__(**params)

        # -------------------------------------------------- Constants --------------------------------------------------
        # _wiki_info_button = button that opens a tab to the GitHub Wiki page of the Data Visualizer app
        self._wiki_info_button = pn.widgets.Button(name = "\u2139", button_type = "light", width = 30)
        self._wiki_info_button.js_on_click(
            args = {"wiki_url": "https://github.com/venuswku/data-visualizer/wiki"},
            code = "window.open(wiki_url)"
        )
        
        # -------------------------------------------------- Internal Class Properties --------------------------------------------------
        
    # -------------------------------------------------- Private Class Methods --------------------------------------------------
    @param.depends("data_map.clicked_transects_info", watch = True)
    def _update_clicked_transects_info(self) -> None:
        """
        Updates PopupModal's clicked_transects_info parameter, which stores information about the most recently clicked transect(s) or user-drawn transect from the data map,
        whenever DataMap's clicked_transects_info parameter changes because the user wants to view information about a different transect.
        """
        if self.data_map.clicked_transects_info:
            self.popup_modal.clicked_transects_info = self.data_map.clicked_transects_info

    @param.depends("data_map.collection", watch = True)
    def _update_time_series_collection_path(self) -> None:
        """
        Triggers event to update PopupModal's _collection_dir_path internal property and its related objects if DataMap's collection parameter changed.
        """
        self.popup_modal.update_collection_dir_path = True
    
    @param.depends("popup_modal.plot_time_series_data", watch = True)
    def _update_last_selected_data_file(self) -> None:
        """
        Updates DataMap's data_file_paths parameter with the most recently selected data files from PopupModal's checkbox widgets.
        """
        self.data_map.data_file_paths = self.popup_modal.user_selected_data_files
    
    # -------------------------------------------------- Public Class Methods --------------------------------------------------
    @property
    def wiki_info_button(self) -> pn.widgets.Button:
        """
        Returns the button widget that opens a tab to the GitHub Wiki page of the Data Visualizer app.
        """
        return self._wiki_info_button
