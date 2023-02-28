#data model to store files and file data
import os
import re
import glob
from pathlib import Path
import extract_mzml_data as exdata
import metatlas.metatlas.tools.validate_filenames as validate

filename_categories_vocab = ['ISTD', 'QC', 'InjBl'] #Words that will be searched for in filename to determine file category. case insensitive
ignore_errors = ['(Filename and parent directory do not contain the same batch fields\.)', '(Parent directory contains .* fields but the minimum allowed is 9\.)']
extraction_control_vocab = 'ExCtrl'

class File():

    def __init__(self, path, name):
        self.path = path
        self.name = name
        self.polarity = self.get_file_polarity()
        self.chromatography = self.get_file_chromatography()
        self.category = self.get_file_category()
        self.run_num = self.get_file_run_number()
        self.ms_num = self.get_file_msnum()
        self.group_name = self.get_group_name()
        self.ms1_data = exdata._extract_ms1_data(self.path)
        self.ms1_tic = exdata._extract_ms1_tic(self.path)

    def get_file_polarity(self):
        return self.name.split('_')[9]

    def get_file_msnum(self):
        return self.name.split('_')[10]
        
    def get_file_category(self):
        """Take filename string as input, return the 'sample type' (S1, ISTD, QC, InjBl) as string"""

        group_field = self.name.split('_')[12]
        optional_field = self.name.split('_')[14]

        search_fields = group_field + optional_field
        file_category = [ele for ele in filename_categories_vocab if (ele.upper() in search_fields.upper())]

        if file_category == []:
            file_category = ['S1']

        if extraction_control_vocab in group_field:
            file_category = ['ExCtrl']

        return file_category[0]

    def get_file_chromatography(self):
        return self.name.split('_')[7]

    def get_file_run_number(self):
        run_num_str = self.name.split('_')[15]
        if not run_num_str.isnumeric():
            run_num = int(re.sub('\\D', '', run_num_str)) #for JGI naming scheme (ex. Run4)
        else:
            run_num = int(run_num_str) #for EGSB naming scheme (ex. '100')

        return run_num

    def get_group_name(self):
        return self.name.split('_')[12]

class RawDataset():

    def __init__(self, path):
        self.path = path
        self.analyzed_files = []
        self.files_to_analyze = []

        self.set_files_to_analyze()

        self.chromatography = self.get_chromatography()

    def set_files_to_analyze(self):
        files = sorted(glob.glob(self.path + '\\*.raw'), key=os.path.getmtime) 
        files = [file for file in files if file not in self.analyzed_files]
        self.files_to_analyze = files

    def add_analyzed_file(self, file):
        self.files_to_analyze.remove(file)
        self.analyzed_files.append(file)

    def get_chromatography(self):
        test_file = self.files_to_analyze[0]
        test_name = os.path.basename(test_file)
        if 'c18' in test_name.split('_')[7].lower():
            chromatography = 'c18'
        if 'hilic' in test_name.split('_')[7].lower():
            chromatography = 'hilic'
        return chromatography

    def validate_experiment_filenames(self):
        filename_warnings = []
        filename_errors = []
        for file in set(self.files_to_analyze + self.analyzed_files):

            file_name_path = Path(os.path.basename(file))

            warnings, errors = validate.get_validation_messages(file_name_path, minimal=True)
            filename_warnings.append({'file_name':file_name_path, 'warnings':warnings})
            filename_errors.append({'file_name':file_name_path, 'errors':errors})
            
        return filename_warnings, filename_errors


    


