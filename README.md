# AMY-HPC-Emotion

Code for analyses in the manuscript:  
**“Transient beta bursts coordinate amygdala–hippocampus activity during emotional processing in humans.”**

## Overview

This repository contains MATLAB code to reproduce the core analyses and figures. We analyze simultaneously recorded **single-unit activity** and **local field potentials (LFPs)** from the **amygdala (AMY)** and **hippocampus (HPC)** during an emotional image-rating task. The pipeline covers behavioral analyses, LFP preprocessing and time–frequency decomposition, burst detection, spike/LFP coupling, and population-level statistics.

> Figure/analysis mapping follows Figs. 1–6 in the manuscript.

## System Requirements

- **Operating System**: Windows 10/11, macOS, or Linux
- **MATLAB**: R2021a or later
- **Toolboxes**:
  - Statistics and Machine Learning Toolbox (for GLMs)
- **Additional Packages**:
  - [FieldTrip toolbox](https://www.fieldtriptoolbox.org/) (for preprocessing, time–frequency decomposition, cluster-based permutation tests)
  - [Wave_Clus3](https://github.com/csn-le/wave_clus) (for spike sorting)

A minimum of **16 GB RAM** is recommended for time–frequency and burst analyses, though smaller datasets can be processed with less.

## Installation Guide

1. **Clone the repository**
  
  ```bash
  git clone https://github.com/damisahlab/AMY-HPC-Emotion.git
  cd AMY-HPC-Emotion
  ```
  
2. **Set up MATLAB paths**
  
  - Download and install [FieldTrip](https://www.fieldtriptoolbox.org/download/).
  - Download and install [Wave_Clus3](https://github.com/csn-le/wave_clus).
  - Add both toolboxes and the repository to your MATLAB path:
    
    ```matlab
    addpath(genpath('path_to_fieldtrip'));
    addpath(genpath('path_to_waveclus'));
    addpath(genpath('path_to_AMY-HPC-Emotion'));
    ```
    
3. **Test your installation**  
  Run the behavioral analysis script to ensure dependencies are loaded correctly:
  
  ```matlab
  run figure_1.m
  ```
  
4. **Reproduce results**  
  Each script `figure_X.m` corresponds to a main figure in the manuscript. Execute them to replicate the analysis and visualization.
  

## Repository Structure

- `figure_1.m` — **Task & behavior**  
- `figure_2.m` — **Single-neuron encoding of valence**  
- `figure_3.m` — **Time–frequency & beta-burst pipeline**  
- `figure_4.m` — **Burst-locked LFP & spiking**  
- `figure_5.m` — **Cell-type contributions & AMY→HPC gating**  
- `figure_6.m` — **Beta–gamma phase–amplitude coupling**  
- `PlotExampleNeuronalActivity.m` — Helper for example unit **rasters**, **PSTHs**, and waveform insets.
- `PlotPopulationActivity.m` — Helper for **PSTH smoothing**, group concatenation, and wrapper functions for cluster-based tests.
- `PlotPopulationGLM.m` — Helper for **GLM** fitting and plotting of coefficients/effects.
  

## Usage

Each `figure_X.m` script reproduces the corresponding figure using processed data. For example:

```matlab
addpath(genpath('path_to_fieldtrip'));
addpath(genpath('path_to_waveclus'));
run figure_1.m
```

⚠️ **Note**: Raw patient data cannot be shared due to privacy restrictions. Instead, scripts are provided with placeholders and example structures that can be adapted to similar datasets.
