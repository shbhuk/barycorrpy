# barycorrpy
Version 0.1
UNDER TESTING

[![DOI](https://zenodo.org/badge/106101961.svg)](https://zenodo.org/badge/latestdoi/106101961)


Barycorrpy is the Python version of Jason Eastman and Jason Wright's IDL code [BaryCorr](http://astroutils.astronomy.ohio-state.edu/exofast/pro/exofast/bary/zbarycorr.pro) based on [Wright and Eastman (2014)](https://arxiv.org/pdf/1409.4774.pdf) . BCPy can be used to calculate the barycentric velocity correction for a star with an accuracy of ~ 1 cm/s . 
To do this, it takes into consideration the following - 

1. Revolution of the Earth to consider position and velocity of the Geo-Center (GC) wrt Solar System Barycenter (SSB)
2. Rotation of the Earth
3. Precession, Nutation and Polar motion of the Earth, along with **2.** to calculate the position and velocity of the observatory wrt the Geo-Center.
4. Gravitational Time dilation due to objects of the Solar System
5. Leap Second offset between UTC and TAI time standards
6. Proper motion and systemic radial velocity of the star
7. Parallax correction
8. Shapiro delay



The installation instructions and the guide on how to run and use the code are explained in the [Wiki](https://github.com/shbhuk/barycorrpy/wiki).


### Leap Second Management

When converting UTC to TDB ([different time standards explained](http://www.cv.nrao.edu/~rfisher/Ephemerides/times.html#TDB)), we need to inlcude for the leap second correction.   
We do not use Astropy for this correction due to the lack of an efficient mechanism to update the files when a new leap second is announced. Leap seconds are hard coded into Astropy's ERFA routines, and thus to update for a new leap second the user would have to update Astropy and re-compile it.   
In this routine `utc_tdb.py`, we incorporate this is a stand alone file which is checked every time the code is run. 


### INSTALLATION INSTRUCTIONS

The instructions for installation and getting started for this package are detailed in the [Wiki](https://github.com/shbhuk/barycorrpy/wiki).
