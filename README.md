# barycorrpy (v0.4.4)

### Please join the google group for updates regarding bug reports, new versions etc:
To sign up for updates, please join the Google Group linked here -
https://groups.google.com/g/barycorrpy

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.1115856.svg)](https://doi.org/10.5281/zenodo.1115856)
[![Build Status](https://travis-ci.com/shbhuk/barycorrpy.svg?branch=master)](https://travis-ci.com/shbhuk/barycorrpy)
[![PyPI version](https://badge.fury.io/py/barycorrpy.svg)](https://badge.fury.io/py/barycorrpy)
[![ASCL:1808.001](https://img.shields.io/badge/ascl-1808.001-blue.svg?colorB=262255)](http://ascl.net/1808.001)
<!-- [![PyPI downloads](https://img.shields.io/pypi/dm/barycorrpy.svg)](https://pypistats.org/packages/barycorrpy) -->
<!-- [![Hits](https://hits.seeyoufarm.com/api/count/incr/badge.svg?url=https%3A%2F%2Fgithub.com%2Fshbhuk%2Fbarycorrpy&count_bg=%23C83D50&title_bg=%23555555&icon=&icon_color=%23E7E7E7&title=hits&edge_flat=false)](https://hits.seeyoufarm.com) -->

[Barycorrpy](http://iopscience.iop.org/article/10.3847/2515-5172/aaa4b7) is the Python version of Jason Eastman and Jason Wright's IDL code [BaryCorr](http://astroutils.astronomy.ohio-state.edu/exofast/pro/exofast/bary/zbarycorr.pro) based on [Wright and Eastman (2014)](https://arxiv.org/pdf/1409.4774.pdf). BCPy can be used to calculate the barycentric velocity correction for a star with an accuracy well below 1 cm/s .
To do this, it takes into consideration the following effects-

1. Revolution of the Earth to consider position and velocity of the geocenter with respect to the Solar System barycenter
2. Rotation of the Earth
3. Precession, nutation and polar motion of the Earth, along with the above to calculate the position and velocity of the observatory with respect to the geocenter
4. Gravitational time dilation due to objects of the Solar System
5. Leap second offset
6. Proper motion and systemic radial velocity of the star
7. Parallax
8. Shapiro delay



The installation instructions and the guide on how to run and use the code are explained in the [wiki](https://github.com/shbhuk/barycorrpy/wiki

There is also a video tutorial (courtesy of the [Exoplanet Modeling and Analysis Center](https://emac.gsfc.nasa.gov/) ) available [here](https://www.youtube.com/watch?v=5SqmL6TdJjs), describing how to get started with barycorrpy.

It is important to note that the code does not just output a barycentric velocity which must be subtracted from the raw Radial Velocity. It outputs the net radial velocity after correcting for barycentric correction. This is because the correction involves a cross term due to the relativistic addition. Therefore include the zmeas in the input parameters.



### JDUTC to BJDTDB converter
As explained in [Eastman et al. 2010](http://adsabs.harvard.edu/abs/2010PASP..122..935E), we also include a JDUTC to BJDTDB time converter.

We include the following corrections -

1. Clock Correction - To correct for difference between UTC and TDB time scales.
2. Geometric Correction - Light travel time from observatory to Solar System Barycenter.
3. Einstein Correction - Relativistic correction since Earth is not an inertial frame.

The output of our function - utc_tdb.JDUTC_to_BJDTDB() matches the [web applet converter](http://astroutils.astronomy.ohio-state.edu/time/utc2bjd.html) to about 0.2ms . Therefore for applications requiring higher precision it should not be used.
We advocate not using the standard Astropy converter this due to the same leap second staleness issue. Therefore for cases requiring such high precision, the leap second should be checked for and be updated as is done by this code.


### INSTALLATION INSTRUCTIONS

The instructions for installation and getting started for this package are detailed in the [wiki](https://github.com/shbhuk/barycorrpy/wiki).

## Citation
To cite using this code you can cite this paper - [RNAAS](http://iopscience.iop.org/article/10.3847/2515-5172/aaa4b7).

Guidelines to cite this package can be found [here](https://github.com/AASJournals/Tutorials/blob/master/Repositories/CitingRepositories.md).

To sign up for updates, please join the Google Group linked here -
https://groups.google.com/forum/#!forum/barycorrpy
