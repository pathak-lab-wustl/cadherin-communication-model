# Repository for organoid invasion with leader cell polarization based on Cdh1/Cdh3 functionality
This repository contains files to run simulations with CompuCell3D. The simulations explore effects of Cdh1 and Cdh3 communication function for leader cell polarization and organoid displacement

## Software
This code used CompuCell3D version 4.2.5. Binaries to install CompuCell3D for Windows or Mac are located on their [main site](https://compucell3d.org/SrcBin).

## Running simulations
To run a parameter scan, values for specific variables can be set in the "ParameterScanSpecs.json" file. The command that can be run in a terminal for a parameter scan is:
```
<path to CC3D paramScan.command> --input="<path to .cc3d file of Simulation folder>" --output-dir="<Path where simulation scans are saved>" --install-dir="<Path to installation folder of CC3D>" --gui --output-frequency=50 --screenshot-output-frequency=50
```


