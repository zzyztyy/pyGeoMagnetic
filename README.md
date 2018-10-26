# pyGeoMagnetic

## What is pyGeoMagnetic?  
It is a Python program about geomagnetic field based on IGRF12 and Apex model. We can calculate magnetic field intensity and transform coordinate between GeoGraphical and GeoMagnetic.  
The source code of IGRF12 and Apex model come from Fortran file. In this package, we reform them to Python file. Although Python can translate Fortran by pyf module, it needs Fortran compiler. This package can run without Fortran compiler.  
****
## How to use it?  
main.py includes some functions which can be used directly.  
*We will give some example when this package is completed.*  
****
## What we did?
This model is called [Magnetic-Apex and Quasi-Dipole coordinate][apex]. It includes too many functions and complex math and physics. So in this package, we change the structure but keep the mathematic method.
The model has two algorithm. One is simple, basic but slow, the other one is complex but fast. We would make the simple one firstly and make another one later.
### The Apex model includes 4 steps
 - **Get Magnetic Intensity**
 - **Map Magnetic Line to Apex Position**
 - **Calculate Magnetic Coordinate**
 - **Smooth Base Vectors**  
 *more information needs to be added after complete this part.*
 
[apex]: https://doi.org/10.1029/2010JA015326

