import re
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

from fisher_py_mod import Raw

from fisher_py.data.business import TraceType
from fisher_py.data import ToleranceUnits, Device
from fisher_py.raw_file_reader.raw_file_access import RawFileAccess

import os
import time
import glob

import pandas as pd
import openpyxl

polarities = ['FPS', 'POS', 'NEG'] #Acceptable polarity values present in the filename. 
underscore_num_set = 15 #Number of underscores that should be in each filename so that untargeted jobs are able to process them after upload


def ppm_diff(observed, theoretical):
	"""Take two numbers as arguments, return PPM difference."""
	return (theoretical - observed) / theoretical * 1000000


def closest(scn, t):
	"""Take a numpy 1d array and a number, return value closest to number contained in the array."""
	idx = (np.abs(scn - t)).argmin()
	return scn[idx]


def get_chromatogram_from_mass(file, mz_t, ppm_tolerance, ms_filter: str):
	"""Take Raw() object and parameters, return 1d arrays of retention time and intensity.
	
	Arguments:
	file -- .raw file from either RawFile() or Raw() object
	mz_t -- theoretical m/z which is used as the MassRange mass 
	ppm_tolerance -- PPM error range of theoretical m/z for MassRange
	ms_filter -- ms or ms2 
	"""
	rt, i = file.get_chromatogram(mz_t, ppm_tolerance, TraceType.MassRange, ToleranceUnits.ppm, ms_filter)
	return rt, i


def get_scan_data_from_chromatogram(file, rt, i, ms_filter: str):
	"""Take Raw() object and parameters, return scan data"""
	max_i_idx = np.argmax(i)
	max_i = i[max_i_idx]

	if ms_filter == 'ms':
		mz, i2, charges, real_rt = file.get_scan_ms1(rt[max_i_idx])
		return mz, i2, charges, real_rt, max_i

	if ms_filter == 'ms2':
		mz, i2, charges, real_rt = file.get_scan_ms2(rt[max_i_idx])
		return mz, i2, charges, real_rt, max_i


def get_compound_data(file, polarity, mz_t, ppm_tolerance, include_ms2_data: bool):
	rt, i = get_chromatogram_from_mass(file, mz_t, ppm_tolerance, 'ms')
	mz, i2, charges, real_rt, max_i = get_scan_data_from_chromatogram(file, rt, i, 'ms')

	mz_o = closest(mz, mz_t)
	ppmdif = ppm_diff(mz_o, mz_t)

	ms1_data = (polarity, mz_o, ppmdif, real_rt, max_i)
	ms2_data = ('NA')

	if include_ms2_data:
		rt_ms2, i_ms2 = get_chromatogram_from_mass(file, mz_t, ppm_tolerance, 'ms2')
		max_i_idx_precursor = np.argmax(i_ms2)
		max_i_precursor = i_ms2[max_i_idx_precursor]

		ms2_data = max_i_precursor

	return ms1_data, ms2_data


def make_figure(df, y: str, include_samples: bool):
	title = df.at[0, 'Compound Name']

	istd_df_pos = df.loc[df['File Category'].isin(['ISTD']) & df['Polarity'].isin(['POS'])]
	istd_df_neg = df.loc[df['File Category'].isin(['ISTD']) & df['Polarity'].isin(['NEG'])]
	fig, axs = plt.subplots(2, sharex=True)

	if include_samples:
		sample_df_pos = df.loc[df['File Category'].isin(['S1']) & df['Polarity'].isin(['POS'])]
		sample_df_neg = df.loc[df['File Category'].isin(['S1']) & df['Polarity'].isin(['NEG'])]

		axs[0].scatter(sample_df_pos['Run Number'], sample_df_pos[y], color='orange', label='Sample Inj.')
		axs[0].scatter(istd_df_pos['Run Number'], istd_df_pos[y], color='blue', label='ISTD Inj.')
		
		axs[1].scatter(sample_df_neg['Run Number'], sample_df_neg[y], color='orange', label='Sample Inj.')
		axs[1].scatter(istd_df_neg['Run Number'], istd_df_neg[y], color='blue', label='ISTD Inj.')

	else:
		axs[0].scatter(istd_df_pos['Run Number'], istd_df_pos[y], color='blue', label='ISTD Inj.')
		axs[1].scatter(istd_df_neg['Run Number'], istd_df_neg[y], color='blue', label='ISTD Inj.')

	axs[0].margins(y=2)
	axs[1].margins(y=2)
	axs[0].set_title('POS')
	axs[1].set_title('NEG')

	y_mean_pos = [np.mean(istd_df_pos[y]) for i in istd_df_pos[y]]
	y_mean_neg = [np.mean(istd_df_neg[y]) for i in istd_df_neg[y]]
	axs[0].plot(istd_df_pos['Run Number'], y_mean_pos, color='red', linestyle='--', label='ISTD Inj Mean')
	axs[1].plot(istd_df_neg['Run Number'], y_mean_neg, color='red', linestyle='--', label='ISTD Inj Mean')
	fig.text(0.06, 0.5, y, ha='center', va='center', rotation='vertical')
	axs[1].set(xlabel='Run Number')
	axs[0].legend(loc="upper right", prop={'size': 6})
	axs[1].legend(loc="upper right", prop={'size': 6})
	fig.suptitle(title)

	return fig


def get_file_age(path):
	modification_time = os.path.getmtime(path)

	current_time = time.time()

	mins_old = (current_time - modification_time) / 60

	return mins_old


class RawDataset(object):

	def __Init__(self, file_paths, file_names, file_data, file_dfs):
		self.file_paths = file_paths
		self.file_names = file_names
		self.file_data = file_data
		self.file_dfs = file_dfs

	def set_file_paths(self, path, files_num = 0):
		"""
		Set the list of raw file paths that will be analyzed. Remove any files from the list that are under 30 minutes old.
        
		Arguments:
		path -- path to experiment directory containing raw files
		files_num -- number of most recent files to analyze. default is 0, which includes all files in directory.
		"""
		file_paths = sorted(glob.iglob(path + '\\*.RAW'), key=os.path.getmtime)
        
		while get_file_age(file_paths[-1]) < 30:
			file_paths.remove(file_paths[-1])
            
		file_paths_len = len(file_paths)
            
		if files_num >= file_paths_len or files_num == 0:
			files_include = 0
		else:
			files_include = file_paths_len - files_num
        
		self.file_paths = file_paths[files_include:file_paths_len]

	def set_file_names(self):
		self.file_names = [os.path.basename(x) for x in self.file_paths]

	@staticmethod
	def _get_file_polarity(name):
		return name.split('_')[9]
	
	@staticmethod
	def _get_file_msnumber(name):
		return name.split('_')[10]

	@staticmethod
	def _get_file_category(name):
		return name.split('_')[14].split('-')[-1]

	@staticmethod
	def _get_file_run_number(name):
		run_num_str = name.split('_')[15]
		run_num = int(re.sub('\\D', '', run_num_str))

		return run_num

	def check_centroid(self):
		"""
		Checks if the first scan of each file is centroid. 
		
		Returns: Dictionary with file name as keys and True or False as values. True means the first scan of the file is centroid.
		"""
		centroid_report = {}
		for file in self.file_paths:
			raw_file = RawFileAccess(file)
			raw_file.select_instrument(Device.MS, 1)

			centroid = raw_file.is_centroid_scan_from_scan_number(1)
			
			if centroid:
				res = True
			else:
				res = False
			
			centroid_report[os.path.basename(file)] = res
			
			raw_file.dispose()
			del raw_file

		return centroid_report

	def check_file_names(self):
		"""
		Checks that each file name has the correct number of underscores
		and has the polarity descriptor in the correct location.
		
		Returns: True if all files are conforming and False if any are not.
		"""
		report = {}
		for name in self.file_names:
			
			failure_report = []
			underscore_num = name.count('_')
			polarity = self._get_file_polarity(name)

			if underscore_num != underscore_num_set:
				res = False
				underscore_fail = 'Incorrect Number of Underscores'
				failure_report.append(underscore_fail) 

			if not [True for ele in polarities if (ele in polarity)]:
				res = False
				polarity_fail = 'Polarity Descriptor Invalid'
				failure_report.append(polarity_fail)
			else:
				res = True
				
			if res == False:
				report[name] = failure_report
			else:
				report[name] = res
				
		return report

	def get_data(self, compounds, ppm_tolerance=10):

		def istd_data_update(run_num, file_name, file_category, compound_name, polarity,
							mz_observed, ppm, retention_time, ms1_intensity, ms2_precursor_intensity):

			run_numbers.append(run_num)
			file_names.append(file_name)
			file_categories.append(file_category)
			compound_names.append(compound_name)
			pols.append(polarity)
			mzs.append(mz_observed)
			ppms.append(ppm)
			rts.append(retention_time)
			ims1s.append(ms1_intensity)
			precursor_ims2s.append(ms2_precursor_intensity)

		self.file_data = {}
		istd_data = dict.fromkeys(['Run Number', 'File Name',
								'File Category', 'Compound Name', 'Polarity', 'Observed M/Z', 'PPM Error', 'RT',
								'MS1 Intensity', 'MS2 Precursor Intensity'])

		for c in compounds:

			istd_data = istd_data.copy()
			istd = compounds[c]
			compound_name, mz_t_pos, mz_t_neg, rt_t, i_t = istd.values()
			(run_numbers, file_names, file_categories, compound_names, pols, mzs, 
			ppms, rts, ims1s, precursor_ims2s) = [], [], [], [], [], [], [], [], [], []

			for f in self.file_paths:
				file = Raw(f)
				file_name = os.path.basename(f)
				file_polarity = self._get_file_polarity(file_name)
				file_category = self._get_file_category(file_name)
				file_msnumber = self._get_file_msnumber(file_name)
				run_num = self._get_file_run_number(file_name)
				include_ms2 = False
				
				if file_msnumber == 'MSMS':
					include_ms2 = True
					
				if file_polarity == 'FPS':
					pos_data_ms1, pos_data_ms2 = get_compound_data(file, 'POS', mz_t_pos, ppm_tolerance, include_ms2)
					neg_data_ms1, neg_data_ms2 = get_compound_data(file, 'NEG', mz_t_neg, ppm_tolerance, include_ms2)

					pos, pos_mz_o, pos_ppmdif, pos_real_rt, pos_max_ms1_i = pos_data_ms1
					pos_max_precursor_i = pos_data_ms2
					neg, neg_mz_o, neg_ppmdif, neg_real_rt, neg_max_ms1_i = neg_data_ms1
					neg_max_precursor_i = neg_data_ms2

					istd_data_update(run_num, file_name, file_category, c, pos, pos_mz_o, pos_ppmdif, pos_real_rt, 
									pos_max_ms1_i, pos_max_precursor_i)
					istd_data_update(run_num, file_name, file_category, c, neg, neg_mz_o, neg_ppmdif, neg_real_rt,
									neg_max_ms1_i, neg_max_precursor_i)

				if file_polarity == 'POS':
					pos_data_ms1, pos_data_ms2 = get_compound_data(file, 'POS', mz_t_pos, ppm_tolerance, include_ms2)

					pos, pos_mz_o, pos_ppmdif, pos_real_rt, pos_max_ms1_i = pos_data_ms1
					pos_max_precursor_i = pos_data_ms2

					istd_data_update(run_num, file_name, file_category, c, pos, pos_mz_o, pos_ppmdif, pos_real_rt,
									pos_max_ms1_i, pos_max_precursor_i)

				if file_polarity == 'NEG':
					neg_data_ms1, neg_data_ms2 = get_compound_data(file, 'NEG', mz_t_neg, ppm_tolerance, True)

					neg, neg_mz_o, neg_ppmdif, neg_real_rt, neg_max_ms1_i = neg_data_ms1
					neg_max_precursor_i = neg_data_ms2

					istd_data_update(run_num, file_name, file_category, c, neg, neg_mz_o, neg_ppmdif, neg_real_rt,
									neg_max_ms1_i, neg_max_precursor_i)

				istd_data['Run Number'] = run_numbers
				istd_data['File Name'] = file_names
				istd_data['File Category'] = file_categories
				istd_data['Compound Name'] = compound_names
				istd_data['Polarity'] = pols
				istd_data['Observed M/Z'] = mzs
				istd_data['PPM Error'] = ppms
				istd_data['RT'] = rts
				istd_data['MS1 Intensity'] = ims1s
				istd_data['MS2 Precursor Intensity'] = precursor_ims2s

				file.close()
				del file

			self.file_data[c] = istd_data

	def make_dfs_from_data(self):

		self.file_dfs = {}

		for key in self.file_data:
			to_df = self.file_data[key]
			df = pd.DataFrame.from_dict(to_df, orient='index').transpose()

			self.file_dfs[key] = df

	def make_qc_plots(self, path):

		if not os.path.isdir(path + '\\qc_output'):
			os.mkdir(path + '\\qc_output')

		output_dir = (path + '\\qc_output')
		
		pdf = PdfPages(output_dir + '/qc_plots.pdf')

		for compound in self.file_dfs:
			df = self.file_dfs[compound]

			fig1 = make_figure(df, 'MS1 Intensity', False)
			fig2 = make_figure(df, 'PPM Error', True)
			fig3 = make_figure(df, 'RT', True)

			pdf.savefig(fig1)
			pdf.savefig(fig2)
			pdf.savefig(fig3)

		pdf.close()

	def export_dfs(self, extension: str, path):

		if not os.path.isdir(path + '\\qc_output'):
			os.mkdir(path + '\\qc_output')

		output_dir = (path + '\\qc_output')

		if extension == 'CSV':

			for compound in self.file_dfs:
				self.file_dfs[compound].to_csv(output_dir + '/' + compound + '.csv', index=False)

		if extension == 'XLSX':

			writer = pd.ExcelWriter(output_dir + '/ISTDS.xlsx')

			for compound in self.file_dfs:
				self.file_dfs[compound].to_excel(writer, sheet_name=compound, index=False, engine=openpyxl)

			writer.save()
