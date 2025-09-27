# Functional stability and recurrent STDP in Rhythmogenesis 

These scripts and functions generate the figures from the paper "Functional stability and recurrent STDP in Rhythmogenesis"

## Requirements
- MATLAB R2021a or newer
- Mathematica (ONLY for DISPLAYING figure 3a)

## Folder layout
- `src/` — reusable MATLAB functions (organized as `proj.*` packages)
- `scripts/` — runnable scripts for reproducing figures
- `figures/` — output figures and files (ignored by git)
- `setup.m` — adds the project to MATLAB path
- `README.md` — project overview

## Quick start
In MATLAB run:
```matlab
cd('path/to/this/repo') % Enter here the path to the repository
setup
run('scripts/YourMainScript.m')