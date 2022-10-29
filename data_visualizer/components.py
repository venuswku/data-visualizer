# External dependencies imports
import param
import holoviews as hv

class DataPlotter(param.Parameterized):
  def __init__(self):
    self.original_dataset = hv.DynamicMap()

  @property
  def plot(self):
    return self.original_dataset