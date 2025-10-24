# Functional stability and recurrent STDP in Rhythmogenesis 

These scripts and functions generate the figures from the paper "Functional stability and recurrent STDP in Rhythmogenesis"

## Requirements
- MATLAB R2021a or newer
- Mathematica (ONLY for DISPLAYING figure 3a)
- MikTex/MaTex (for displaying titles and labels in Mathematica - not mandatory)

## Folder layout
- `src/` — reusable MATLAB functions (organized as `proj.*` packages)
- `scripts/` — runnable scripts for reproducing figures
- `figures/` — output figures and files 
- `setup.m` — adds the project to MATLAB path
- `README.md` — project overview

## Quick start

## Install path
In MATLAB run:
cd('path/to/this/repo') % Enter here the path of the repository
setup

%% Generate Figures 
1. In folder "scripts" find the folder of the figure that you want to generate and ENTER (important!) this folder.
2. Open the "FigureX_Simulation".
3. Run (Click F5).

%% Notes 
- Folders named "private" in some of the "scripts/figureX" folders contain auxiliary functions which are used only by their main script "FigureX_Simulation".
- Folders named "data" in some of the "scripts/figureX" folders contain some of the data generated and saved by the their main script "FigureX_Simulation". 