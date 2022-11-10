import param
from panel.template.theme import Theme
import panel as pn

class DefaultCustomTheme(Theme):
    """
    DefaultCustomTheme provides the default custom color palette for the Data Visualizer app.
    """
    print()

class BootstrapDefaultCustomTheme(DefaultCustomTheme):
    # Add extra CSS to apply on the DefaultCustomTheme.
    css = param.Filename(default = "./default_custom.css")

    # Specify which implementation the DefaultCustomTheme should use.
    _template = pn.template.BootstrapTemplate