#!/bin/bash
## script changed from original by for bioconda-lefse by Nicholas Lund & Olivia Whitelaw 2025.12.03
## adapted from the original lefse run.sh available at:
## https://bitbucket.org/nsegata/lefse/src/54694b4b0d9e335ff1ecafff8db4f1e0cf7004da/example/run.sh?at=default&fileviewer=file-view-default
##
## Modified to use local LEfSe installation instead of conda
##

# Set up paths for local LEfSe installation
export PYTHONPATH="project/lefse:$PYTHONPATH"
cd project/example

# lefse_format_input.py convert the input data matrix to the format for LEfSe.
python3 ../lefse/lefse_format_input.py mouse_data.txt mouse_data.in -c 1 -s 2 -u 3 -o 1000000

# lefse_run.py performs the actual statistical analysis
python3 ../lefse/lefse_run.py mouse_data.in mouse_data.res -l 0.1

# lefse_plot_res.py visualizes the output
python3 ../lefse/lefse_plot_res.py mouse_data.res mouse_data.png

# lefse_plot_cladogram.py visualizes the output on a hierarchical tree
python3 ../lefse/lefse_plot_cladogram.py mouse_data.res mouse_data.cladogram.png --format png --dpi 600

# Create a directory for storing the raw-data representation of the discovered biomarkers
mkdir -p biomarkers_raw_images

# lefse_plot_features.py visualizes the raw-data features
python3 ../lefse/lefse_plot_features.py mouse_data.in mouse_data.res biomarkers_raw_images/