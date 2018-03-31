# barycorrpy
Version 0.2

[![DOI](https://zenodo.org/badge/106101961.svg)](https://zenodo.org/badge/latestdoi/106101961)


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

It is important to note that the code does not just output a barycentric velocity which must be subtracted from the raw Radial Velocity. It outputs the net radial velocity after correcting for barycentric correction. This is because the correction involves a cross term due to the relativistic addition. Therefore include the zmeas in the input parameters.


### Leap Second Management

When converting UTC to TDB ([different time standards explained](http://www.cv.nrao.edu/~rfisher/Ephemerides/times.html#TDB)), we need to inlcude for the leap second correction.   
We do not use Astropy for this correction due to the lack of an efficient mechanism to update the files when a new leap second is announced. Leap seconds are hard coded into Astropy's ERFA routines, and thus to update for a new leap second the user would have to update Astropy and re-compile it.   
In this routine `utc_tdb.py`, we incorporate this is a stand alone file which is checked every time the code is run. Details for this can be found in the [wiki](https://github.com/shbhuk/barycorrpy/wiki).

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
