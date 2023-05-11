importScripts("https://cdn.jsdelivr.net/pyodide/v0.22.1/full/pyodide.js");

function sendPatch(patch, buffers, msg_id) {
  self.postMessage({
    type: 'patch',
    patch: patch,
    buffers: buffers
  })
}

async function startApplication() {
  console.log("Loading pyodide!");
  self.postMessage({type: 'status', msg: 'Loading pyodide'})
  self.pyodide = await loadPyodide();
  self.pyodide.globals.set("sendPatch", sendPatch);
  console.log("Loaded!");
  await self.pyodide.loadPackage("micropip");
  const env_spec = ['https://cdn.holoviz.org/panel/0.14.4/dist/wheels/bokeh-2.4.3-py3-none-any.whl', 'https://cdn.holoviz.org/panel/0.14.4/dist/wheels/panel-0.14.4-py3-none-any.whl', 'pyodide-http==0.1.0', 'data_visualizer']
  for (const pkg of env_spec) {
    let pkg_name;
    if (pkg.endsWith('.whl')) {
      pkg_name = pkg.split('/').slice(-1)[0].split('-')[0]
    } else {
      pkg_name = pkg
    }
    self.postMessage({type: 'status', msg: `Installing ${pkg_name}`})
    try {
      await self.pyodide.runPythonAsync(`
        import micropip
        await micropip.install('${pkg}');
      `);
    } catch(e) {
      console.log(e)
      self.postMessage({
	type: 'status',
	msg: `Error while installing ${pkg_name}`
      });
    }
  }
  console.log("Packages loaded!");
  self.postMessage({type: 'status', msg: 'Executing code'})
  const code = `
  
import asyncio

from panel.io.pyodide import init_doc, write_doc

init_doc()

#!/usr/bin/env python
# coding: utf-8

# In[ ]:


# cd C:\Users\Venuxk\Projects\data-visualizer
# conda activate visualizer
# panel serve --show data_viewer_app.ipynb

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
elwha_river_delta_collection_id = "5a01f6d0e4b0531197b72cfe"

# Map each data collection (name of folders in the root data directory) to a list of column names, which contains data used for the time-series.
collection_time_series_data = {
    elwha_river_delta_collection_id: ["Ortho_Ht_m", "Ortho_ht_m", "ortho_ht_m", "F-W Mean"]
}
# Map each data collection (name of folders in the root data directory) to a name for the data viewer app's title.
collection_app_title = {
    elwha_river_delta_collection_id: "Elwha Topo-Bathy Data Viewer"
}

# Set the main color for the app.
app_main_color = "#2196f3"

# -------------------------------------------------- Initializing Data Visualization App --------------------------------------------------

# Instantiate the app's template.
template = pn.template.BootstrapTemplate(
    title = collection_app_title[elwha_river_delta_collection_id],
    header_background = app_main_color
	# theme = DefaultCustomTheme
)

# Instantiate the main components required by the Application.
data_map = DataMap(
    time_series_data = collection_time_series_data
)
popup_modal = PopupModal(
	data_map = data_map,
    elwha_collection_name = elwha_river_delta_collection_id
)

# Create the application.
app = Application(
	data_map = data_map,
	popup_modal = popup_modal,
	template = template
)

# Populate the template with the sidebar, main, and modal layout.
template.sidebar.extend([
    # app.wiki_info_button,
	pn.panel(app.sidebar_widgets),
    pn.panel(app.sidebar_accordion)
])
template.main.append(pn.panel(data_map.plot, loading_indicator = True))
template.modal.extend([
	pn.panel(popup_modal.content, loading_indicator = True)
])

# Use the Panel extension to load BokehJS, any pn.config variables, any custom models required, or optionally additional custom JS and CSS in Jupyter notebook environments.
pn.extension(loading_spinner = "dots", loading_color = app_main_color, sizing_mode = "stretch_width")#, nthreads = 0

# Launch the app (\`panel serve --show --autoreload app.py\`).
template.servable()



await write_doc()
  `

  try {
    const [docs_json, render_items, root_ids] = await self.pyodide.runPythonAsync(code)
    self.postMessage({
      type: 'render',
      docs_json: docs_json,
      render_items: render_items,
      root_ids: root_ids
    })
  } catch(e) {
    const traceback = `${e}`
    const tblines = traceback.split('\n')
    self.postMessage({
      type: 'status',
      msg: tblines[tblines.length-2]
    });
    throw e
  }
}

self.onmessage = async (event) => {
  const msg = event.data
  if (msg.type === 'rendered') {
    self.pyodide.runPythonAsync(`
    from panel.io.state import state
    from panel.io.pyodide import _link_docs_worker

    _link_docs_worker(state.curdoc, sendPatch, setter='js')
    `)
  } else if (msg.type === 'patch') {
    self.pyodide.runPythonAsync(`
    import json

    state.curdoc.apply_json_patch(json.loads('${msg.patch}'), setter='js')
    `)
    self.postMessage({type: 'idle'})
  } else if (msg.type === 'location') {
    self.pyodide.runPythonAsync(`
    import json
    from panel.io.state import state
    from panel.util import edit_readonly
    if state.location:
        loc_data = json.loads("""${msg.location}""")
        with edit_readonly(state.location):
            state.location.param.update({
                k: v for k, v in loc_data.items() if k in state.location.param
            })
    `)
  }
}

startApplication()