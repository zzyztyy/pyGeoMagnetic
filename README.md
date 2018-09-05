# IGRF12

Re-code IGRF12 model from Fortran to Python3
*********************
The IGRF12 model include 4 function:
  - **IGRF12** *Main function for printing and calculate other component*
	We would rather replace this by different function than recode it.
	Beacuse there are many lines for printing on screen or getting legal parameter.
  - **DMDDEC**    *From "Degree(int) Minute(int)"  to "Degree(float)"*
	We would replace it.
  - **DDECDM** *From "Degree(float)" to "Degree(int) Minute(int)"*
	We would raplace it.
  - **igrf12syn** *Math part of the model*
	We use a python package called goto to re-code in every single line.
	We Remove goto and restruct it.
	We would Replace old parameter and logic by Python's for the model can be understanding more easily.
*****************************
We would also add some auxiliary functions. For example:
 - Transform between geolat&geolon and maglat&maglon
 - Track the magnetic line
******************************
*PS: Maybe pyapex can be references.*

