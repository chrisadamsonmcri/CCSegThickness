# CCSegThickness


## Software Pipeline for Midsagittal Corpus Callosum Thickness Profile Processing

CCSegThickness performs fully automated pipeline for thickness profile evaluation and analysis of the human corpus callosum (CC) in 3D structural T1-weighted magnetic resonance images

The pipeline performs the following sequence of steps:

1. Midsagittal plane extraction
2. CC segmentation algorithm - Automated Segmentation
3. Quality control tool - Manual Editor
4. Thickness profile generation
5. Statistical analysis - Group-Wise Statistical Comparison
6. Results figure generator - Results Display

## Installation

Tested in Ubuntu 16.04 LTS, 18.04 LTS.

`git clone https://github.com/chrisadamsonmcri/CCSegThickness`

The following packages are dependencies and need to be installed via apt

`sudo apt install python3-dev python3-pip python3-numpy python3-scipy python3-nibabel python3-opencv python3-matplotlib libsuitesparse-dev cython3`

The following package needs to be installed in pip

`sudo -H pip3 install -U scikits.sparse`

These installations can be done with the script

`install_dependencies.sh`

Finally, make the streamline cython file

`make_cython.sh`

If you use the software in a publication, please cite the following article:

Adamson, C., Beare, R., Walterfang, M. et al. Neuroinform (2014) 12: 595. https://doi.org/10.1007/s12021-014-9236-3
