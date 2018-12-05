# DC2 cosmology

A place for DC2 related analysis and pipeline code.

Everything that follows assuming you're running on NERSC.

## Overview

In terms of information flow, almost everything is set by the config files. The basic idea is that one YAML file represents one (galaxy) sample, with a particular selection function.

As an example, perhaps look at `config/cosmoDC2_source.yaml` in this repository.
The settings are grouped into blocks. From the top,

#### basic

Sets the top-level details to do with the sample and where to save outputs. 

The working directory is given by `workdir`. Note that this should be a full path. A directory at this location, and various sub-directories, will be created if they aren't already there.

`simulation` tells the code which parent cosmological simulation to load. The initial load step calls the GCRCatalogues package https://github.com/LSSTDESC/gcr-catalogs, which makes things much easier (the github link also has a list of currently available catalogues).

`columns` is a whitespace-delimited list of column names. See https://github.com/LSSTDESC/gcr-catalogs/blob/master/GCRCatalogs/SCHEMA.md

The selection mask is defined by `cuts`. You can specify as many cuts as you like, separated by spaces. Each one has the format column_name,value,upper/lower. So in the example we have two:

*`mag_i_lsst,24.1,-` places an upper bound on the i-band magnitude at i=24.1

*`redshift,0.2,+` places a lower bound on redshift at z=0.2

You can also choose to impose an additional (slightly more complicated) selection using the `colour_split` field. Values of either "red" or "blue" will define a cut in colour-magniutude space (as defined by true magnitudes). This bit is a work in progress.

Finally `sample` is the name of the population defined by the above choices. This is just used to decide on file names and internal bookkeeping; in the examples in this repo we've make unimaginative choices ("source" and "lens"), but the sample names can be anything you like.

### nofz

Sets the redshift distribution parameters. Currently we have only two: `zbins` is the number of (equal-number) redshift bins to use and `sigma` is the width of the z=0 redshift error distribution.

### 2pt

Sets the parameters used in the correlation function calculations. The binning (number, lower, upper) are set respectively by `tbins`, `tmin` and `tmax`. You'll also need to specify `ctype`, which gives one half of the two point correlatees (not sure if that's a word, but I'm using it anyway). When you're calling the code to do a two point measurement you'll need to give two config files (and so two samples, and two `ctypes`).

A few sensible pairings might be: `shear-position` (a.k.a. gamma_t), `position-position` (a.k.a. wtheta), `shear-shear` (a.k.a. xipm). You can also do the same with "ellipticity" in place of "shear".


## Commands

### Set up

Before you do anything you'll need to load the depedencies: `source setup_descsims`
Also make sure the repo root directory is in your python path.

### Redshift Distributions

The following generates some redshift distributions using a simple error model:

`python -m dc2_cosmopipe.scripts.nofz config/example_config.yaml`

This should load in one of the extragalactic catalogues (which one is defined in the config file; see below). Depending on the number of galaxies, this may take a little while.
Assuming the code runs to completion, you should end up with several new files in working directory (also set in the config file). Essentially what we're doing is first imposing some cuts to generate a semi-realistic sample; we divide the surviving galaxies into equal-number redshift bins by z_true, and stack Gaussian error distributions to give smooth mock photo-z distributions. Both sets of n(z) (the histograms of true z and the error-convolved versions) are saved as text files.

### 2pt Correlations

The following generates a particular set of two-point functions:

`python -m dc2_cosmopipe.scripts.2pt config/example_config1.yaml config/example_config2.yaml`

Each config file represents a sample, defined by its own set of cuts. The code loads the two catalogues defined by example_config1.yaml and example_config2.yaml and correlates the two. The binning and type of correlation (xipm, gammat, wtheta), and also which columns to use (e.g. true shear or ellipticity) are specified in the config files.
To generate a full 3x2pt datavector you need to call the above three times, once for each of the shear-shear, shear-position and position-position correlations.
The results are saved as text files in the working directory under names that (hopefully) make what they are self-evident.

This section of code is really just a wrapper for treecorr, which is what's used for the pair counting. The idea is (a) to make it easy to rerun with different simulations and/or samples with minimal hassle and (b) to nicely handle things like randoms and mean shear subtraction (if and when it's needed) in a consistent way.
