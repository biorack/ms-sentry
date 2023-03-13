# MS-Sentry version 2.0.1

This script is meant to be deployed on Thermo Fisher Orbitrap LC-MS instrument PCs to automate basic data quality control. Still in development.

## Installation

Install required python dependencies and [ThermoRawFileParser](https://github.com/compomics/ThermoRawFileParser) 

## Usage:

Run ms-sentry.py in Windows command line and provide the experimental directory and optional parameters as input. 

```
usage: ms_sentry.py [-h] [-num NUM_FILES] [-min MIN_FILE_AGE] [-skip SKIP_BLANKS] [-export EXPORT_NUM] directory

MS-Sentry generates plots and tables useful for performing quality control on Thermo Orbitrap data.

positional arguments:
  directory             experimental directory containing Thermo raw files

optional arguments:
  -h, --help            show this help message and exit

  -num NUM_FILES, --num_files NUM_FILES
                        number of files to analyze in dataset. default is all files currently in directory
  -min MIN_FILE_AGE, --min_file_age MIN_FILE_AGE
                        minimum file age in minutes. default is 30
  -skip SKIP_BLANKS, --skip_blanks SKIP_BLANKS
                        skip blank injections. default is True
  -export EXPORT_NUM, --export_num EXPORT_NUM
                        export partial analysis every n number of files. default is None
```

## Dependencies:

Python 3.9.xx (Incompatible with Python 3.10.xx because of Pythonnet)

numPy (1.22.x)

pandas (1.4.x)

matplotlib (3.5.x)

openpyxl (3.0.x)

pythonnnet (2.5.x)

  Note: pythonnet is a dependency of fisher_py, which is the python wrapper of the ThermoRawFileReader DLLS that this program relies on. The pythonnet version for Python   3.10 is incompatible with fisher_py, however, pythonnet does not officially support Python 3.9. Therefore, please install the 3.9 pre-release wheel file directly (can   be found here: https://www.lfd.uci.edu/~gohlke/pythonlibs/#pythonnet). 
  
fisher_py (1.0.xx)

pymzml (2.5.x)
