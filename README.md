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
### Part 1  Re-code IGRF12 model from Fortran to Python3
###### The IGRF12 model includes 4 functions
  - **IGRF12** *Main function for printing and calculate other component*  
	We would rather replace this by different function than re-code it.
	Beacuse there are many lines for printing on screen or getting legal parameter.
  - **DMDDEC**    *From "Degree(int) Minute(int)"  to "Degree(float)"*  
	We do not need it.
  - **DDECDM** *From "Degree(float)" to "Degree(int) Minute(int)"*  
	We do not need it.
  - **igrf12syn** *Math part of the model*  
	We use a python package called goto to re-code in every single line. Then We remove goto and restruct it.   
	We would replace old parameter and logic by Python's for the model can be understanding more easily.  
	We separate some auxiliary functions from igrf12syn and make two function in main.py to use the previous method.  
### Part 2  Reform the Fortran file of Apex model
This model is called [Magnetic-Apex and Quasi-Dipole coordinate][apex]. It includes too many functions and complex math and physics. So in this package, we change the structure but keep the mathematic method.
The model has two algorithm. One is simple, basic but slow, the other one is complex but fast. We would make the simple one firstly and make another one later.
###### The Apex model includes 4 steps
 - **Get Magnetic Intensity**
 - **Map Magnetic Line to Apex Position**
 - **Calculate Magnetic Coordinate**
 - **Smooth Base Vectors**  
 *more information needs to be added after complete this part.*
 
[apex]: https://doi.org/10.1029/2010JA015326

