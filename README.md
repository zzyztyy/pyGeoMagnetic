# pyGeoMagnetic
This is a Python program about geomagnetic field based on IGRF (International Geomagnetic Reference Field) and Apex model. It can transform coordinate between GeoGraphical and GeoMagnetic.  
The source code Apex model come from Fortran file. In this package, we reform them by using Python. Although Python can translate Fortran by `pyf` module, it needs Fortran compiler. This package can run __without Fortran compiler__.  
****
## How to use it?  
Download latest released version and unzip it. Then Install (requires NumPy and pyIGRF before installation):

    python setup.py install
Here is an example:
```python
import pyGeoMagApex as pgma
# GeoGraphic to GeoMagnetic
mlat, mlon = pgma.gd2qd(lat=39.3, lon=116.1)
# GeoMagnetic to GeoGraphic
glat, glon = pgma.qd2gd(mlat=23.1, mlon=-171.1)
# altitude and date can be set, default as 0km and 2005
mlat, mlon = pgma.gd2qd(lat=20.2, lon=103.1, alt=100, date=2006.)
```
****
## Apex Model and Modified Apex Model
This package is refer to [Magnetic-Apex and Quasi-Dipole coordinate][apex]. This paper give some Fortran file about two models called Apex model and modified Apex model which includes too many functions and complex math and physics. So in this package, we change the structure but keep the mathematic method but only for Apex model.
Modified Apex model is too difficult but fast. We will do it next setp *(Maybe)*.

[apex]: https://doi.org/10.1029/2010JA015326
