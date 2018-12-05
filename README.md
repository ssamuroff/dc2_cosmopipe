# DC2 cosmology

A place for DC2 related analysis and pipeline code.

Everything that follows assuming you're running on NERSC.

## Commands

Before you do anything you'll need to load the depedencies: `source setup_descsims`

The following generates some redshift distributions using a simple error model. 

`python -m dc2_cosmopipe.scripts.nofz config/example_config.yaml`

This should load in one of the extragalactic catalogues (which one is defined in the config file; see below). Depending on the number of galaxies, this may take a little while.
Assuming the code runs to completion, you should end up with several new files in working directory (also set in the config file). Essentially what we're doing is first imposing some cuts to generate a semi-realistic sample; we divide the surviving galaxies into equal-number redshift bins by z_true, and stack Gaussian error distributions to give smooth mock photo-z distributions. Both sets of n(z) (the histograms of true z and the error-convolved versions) are saved as text files.


