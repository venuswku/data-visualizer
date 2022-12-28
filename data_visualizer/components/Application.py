# Standard library imports

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
    @param.depends("data_map.clicked_transects_info", "popup_modal.user_selected_data_files", watch = True)
    def _update_clicked_transects_info(self) -> None:
        """
        Updates the pipe that stores information about the most recently clicked transect(s) from the data map whenever:
            - DataMap's clicked_transects_info parameter changes
            - PopupModal's user_selected_data_files parameter changes because the time-series needs to be updated
              when the user modifies what data should be displayed
            - PopupModal's clicked_transect_buffer parameter changes because the time-series and clicked transects map
              needs to be updated when the user modifies the search radius for extracting point data near a transect
        """
        if self.data_map.clicked_transects_info:
            self.popup_modal.clicked_transects_pipe.event(data = self.data_map.clicked_transects_info)
