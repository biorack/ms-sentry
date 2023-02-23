#convert raw file to mzml using ThermoRawFileParser 

import os
import subprocess

current_dir = os.path.dirname(__file__)
thermo_parser_path = os.path.join(current_dir, 'ThermoRawFileParser\ThermoRawFileParser.exe')
out_dir = os.path.join(current_dir, 'mzml_temp')

def _raw_to_mzml(file_path):

    command = '{converter} -i={input} -o={output} -L=1- -l=2 -f=2'.format(converter=thermo_parser_path, input=file_path, output=out_dir)
    conversion_log = subprocess.run(command, capture_output=True)
    
    return conversion_log
