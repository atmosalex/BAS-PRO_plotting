# Overview

BAS-PRO plotting is a set of Python scripts to plot model solutions generated by the BAS-PRO model.

# Getting started:

1. Clone this repository by typing ```git clone <repo>``` replacing <repo> with the URL (on Github/Gitlab, go to Clone > Clone with HTTPS and copy the link)

2. Optionally, create a new Python virtual environment. Install requirements by running, from the BAS-PRO plotting directory: ```pip install -r requirements.txt``` 

3. Generate a template configuration file for the plotting routines. Either using the Python interpreter (type ```python``` to start), or from a Python script, run:
```
import plot
plot.regenerate_template_config()
```

4. A configuration template will be generated in the ```basproplot_configurations``` directory. Rename it to ```example.json``` and modify the template based on the plotting requirements (see guide below).

5. Plots can now be generated from Python. Continuing from step 3., run:
```
plot.generate_from_config("full_path_to_solution_dir", "basproplot_configurations/example.json") 
```
replacing the paths as appropriate.

# This documentation is a work in progress