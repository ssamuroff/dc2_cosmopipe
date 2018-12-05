# DC2 cosmology

A place for DC2 related analysis and pipeline code.

Everything that follows assuming you're running on NERSC.

## Commands

### Set up

Before you do anything you'll need to load the depedencies: `source setup_descsims`

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
