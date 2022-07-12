# README

This repository is meant to support two separate manuscripts:

1. [Plant spectra as integrative measures of plant phenotypes](https://ecoevorxiv.org/bfc5t/) by Kothari & Schweiger (2022), _EcoEvoRxiv_ DOI: 10.32942/osf.io/bfc5t (in press at _Journal of Ecology_)
2. [Predicting leaf traits across functional groups using reflectance spectroscopy](https://www.biorxiv.org/content/10.1101/2022.07.01.498461v2) by Kothari et al. (2022), _bioRxiv_ DOI: 10.1101/2022.07.01.498461

The repository is maintained by Shan Kothari (shan.kothari \[at\] umontreal \[dot\] ca). The current version of the repository is being archived at Zenodo (DOI forthcoming) and is meant to represent the analyses carried out in the version of manuscript (1) that will appear in _Journal of Ecology_. The code base may change further before manuscript (2) is published in the journal.

## Components

The repository contains scripts numbered as belonging to 10 distinct 'stages' of analysis, as well as stage 00 (useful functions that are called in later scripts).

The stages are:

1. Processing reflectance and transmittance spectra from the raw form they were downloaded in (see https://data.caboscience.org/leaf/). The spectra have already been resampled to 1 nm resolution, but here I interpolate over the sensor overlap region, average leaves for a sample, and apply a Savitzky-Golay filter. These steps are done per project.
2. Compiling spectral data from across projects.
3. Attaching trait data to the spectral data.
4. Dividing up training and testing data for PLSR analyses.
5. Calibrating models using different kinds of data (reflectance A + D, transmittance B, absorptance C, continuum-removed reflectance E, brightness-normalized reflectance F) and traits (area-based \[D\] vs mass-based \[all others\] chemical traits) 
6. Plotting models with (A) just reflectance or (B) comparing reflectance, transmittance, and absorptance.
7. Comparing trait distributions in the data to TRY data.
8. Evaluating transferability of models across functional groups.
9. Externally validating the models built in step 5 across LOPEX, ANGERS, and Dessain data.
10. Identifying the dimensionality and visualizing major dimensions of the spectral data, along with some miscellaneous analyses described in manuscript (1) above.

Manuscript (1) only involves scripts 1-3 and 10. Manuscript (2) involves all scripts. The processed spectral and trait data products produced at the end of script 3, and read in at the beginning of script 4, are archived elsewhere (see **Associated data** below).

## How to use

This repository is *not* a software package or any sort of user-oriented product that people can use without further modification. It *is* meant to be a reasonably well-documented and faithful record of the analyses carried out in the two manuscripts listed above. Some analyses should be easily reproducible, with some modification, given the scripts and the archived data (see **Associated data** below). However, we do not (for example) include the TRY data.

## Associated data

There are a few associated data products (DOIs forthcoming):

1. The main CABO dataset, including traits and spectra, are found [here](https://ecosis.org/package/cabo-2018-2019-leaf-level-spectra) at EcoSIS. Reflectance, transmittance, and absorptance spectra are available.
2. The complete Dessain project is not yet archived, but will be upon publication of manuscript (2).
3. [LOPEX](https://ecosis.org/package/leaf-optical-properties-experiment-database--lopex93-) and [ANGERS](https://ecosis.org/package/angers-leaf-optical-properties-database--2003-), used in the external validation, are available at those respective links
4. In script 10, I read in [fresh-](https://ecosis.org/package/fresh-leaf-cabo-spectra-from-herbarium-project), [pressed-](https://ecosis.org/package/pressed-leaf-cabo-spectra-from-herbarium-project), and [ground-](https://ecosis.org/package/pressed-leaf-cabo-spectra-from-herbarium-project)leaf spectra from another manuscript ([Kothari et al. 2022](https://www.biorxiv.org/content/10.1101/2021.04.21.440856), in press _Methods in Ecology and Evolution_) available at those respective links.
