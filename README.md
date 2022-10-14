# ms-sentry version 1.1.0

This script is meant to be deployed on Thermo Fisher Orbitrap LC-MS instrument PCs to automate some basic data quality control.

Dependencies:

Python 3.9.xx (Incompatible with Python 3.10.xx because of Pythonnet)

numPy (1.22.x)

pandas (1.4.x)

matplotlib (3.5.x)

openpyxl (3.0.x)

pythonnnet (2.5.x)

  Note: pythonnet is a dependency of fisher_py, which is the python wrapper of the ThermoRawFileReader DLLS that this program relies on. The pythonnet version for Python   3.10 is incompatible with fisher_py, however, pythonnet does not officially support Python 3.9. Therefore, please install the 3.9 pre-release wheel file directly (can   be found here: https://www.lfd.uci.edu/~gohlke/pythonlibs/#pythonnet). 
  
fisher_py (1.0.xx)
