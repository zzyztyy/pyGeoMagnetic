# IGRF12
Re-code IGRF12 model from Fortran to Python3

The IGRF12 model include 4 function:
  1. IGRF12: main function for printing and calculate other composize
  2. DMDDEC: checking
  3. DDECDM: checking
  4. igrf12syn: math part of the model

Now we trying to re-code by 3 steps:
  1. Use a python package called goto to re-code in every single line.
  2. Remove goto and restruct.
  3. Replace old parameter and logic by Python's for the model can be understanding more easily.
 
Done:
  2-2
  3-2
  4-2

Doing:
  4-3

Waitting list:
  Almost everything!
