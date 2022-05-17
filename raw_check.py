from ms_sentry_util import RawDataset
import os
import sys
import pandas as pd

phe_params = {'Name' : '13C,15N-Phenylalanine', 'Positive M/Z' : 176.1135, 'Negative M/Z' : 174.0989, 'Ideal RT' : 9.00, 'Ideal Intensity' : 5000000} 
ab_params = {'Name' : 'ABMBA', 'Positive M/Z' : 229.9811, 'Negative M/Z' : 227.9665, 'Ideal RT' : 5.00, 'Ideal Intensity' : 5000000}
trp_params = {'Name' : '13C,15N-Tryptophan', 'Positive M/Z' : 218.1281, 'Negative M/Z' : 216.1136, 'Ideal RT' : 10.00, 'Ideal Intensity' : 5000000}
cyt_params = {'Name' : '13C2,15N3-Cytosine', 'Positive M/Z' : 117.0484, 'Negative M/Z' : 115.0338, 'Ideal RT' : 5.00, 'Ideal Intensity' : 5000000}

compounds = {phe_params['Name'] : phe_params, ab_params['Name'] : ab_params, trp_params['Name'] : trp_params, cyt_params['Name'] : cyt_params}

print('-------------------------------------------------------------------------------------')
print('Please provide the path to the directory where the experiment .raw files are stored.')
print('Use Windows file-path format')

while True:
	raw_path = input('Path to Experiment Directory:')

	if not os.path.exists(raw_path):
		print('This path does not exist')
		continue 

	else:
		expdir = os.path.abspath(raw_path)
		print('Path to directory saved!')
		break

print('-------------------------------------------------------------------------------------')
print('Please enter the number of most recent files you would like to analyze.')
print('Any files modified less than 30 minutes ago have automatically been excluded')

while True:
	files_number = input('Number of most recent files to be analyzed:')

	if not files_number.isnumeric():
		print('Input not integer')
		continue

	else:
		break
		
exp1 = RawDataset()
exp1.set_file_paths(expdir, files_num = int(files_number))
exp1.set_file_names()

print('-------------------------------------------------------------------------------------')
print('Checking filenames and data collection mode (profile/centroid)...')

name_conformation_report = exp1.check_file_names()
file_centroid_report = exp1.check_centroid()

for name, result in name_conformation_report.items():
	if result != "":
		print('Non-conforming filenames have been detected')
		print('Reason(s): ' + str(result))
		summary_result_filenames = False
		break

	else:
		summary_result_filenames = True

if summary_result_filenames == True:
	print('Filenames are conforming')

if summary_result_filenames == False: 

	while True:
		user_export = input('Would you like to export failure report? (y/n):')

		if user_export == 'y':
			print('Exporting...')
			name_report = pd.DataFrame.from_dict(name_conformation_report, orient = 'index')

			if not os.path.isdir(expdir + '\\failure_report'):
				os.mkdir(expdir + '\\failure_report')

			name_report.to_csv(expdir + '\\failure_report\\filenames.csv')
			break

		if user_export == 'n':
			break

		else:
			print('Input not recognized')
			continue

	while True:
		user_continue = input('Would you like to try continuing with analysis? (y/n):')

		if user_continue == 'y':
			print('Continuing...')
			break

		if user_continue == 'n':
			sys.exit('Exiting...')
			break

		else:
			print('Input not recognized')
			continue

for name, result in file_centroid_report.items():
	if result != True:
		print('Warning: data collected in profile mode')
		summary_result_centroid = False
		break

	else:
		summary_result_centroid = True


if summary_result_centroid == True:
	print('Data collected in centroid mode')

if summary_result_centroid == False: 

	while True:
		user_continue = input('Would you like to continue with analysis? (y/n):')

		if user_continue == 'y':
			print('Continuing...')
			break

		if user_continue == 'n':
			sys.exit('Exiting...')
			break

		else:
			print('Input not recognized')
			continue

print('-------------------------------------------------------------------------------------')
print('Extracting data from .raw files. This can take a while...')
exp1.get_data(compounds)
exp1.make_dfs_from_data()

print('Exporting plots and tables')
exp1.make_qc_plots(expdir)
exp1.export_dfs(expdir)

print('-------------------------------------------------------------------------------------')
print('Done!')
print('-------------------------------------------------------------------------------------')