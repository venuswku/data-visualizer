# Standard library importspn
import os

# External dependencies imports
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
        
        # -------------------------------------------------- Internal Class Properties --------------------------------------------------
        
    # -------------------------------------------------- Private Class Methods --------------------------------------------------
    @param.depends("data_map.clicked_transect_data", "popup_modal.user_selected_data_files", watch = True)
    def _update_clicked_transect_data(self) -> None:
        """
        Updates the pipe that stores information about the most recently clicked transect from the data map (whenever DataMap's clicked_transect_data parameter changes).
        Also triggers a new event whenever the PopupModal's user_selected_data_files parameter changes because a new time-series needs to be created when the user modifies what data should be displayed.
        """
        self.popup_modal.clicked_transect_pipe.event(data = self.data_map.clicked_transect_data)
