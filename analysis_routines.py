#routines for analyzing the mzml data after extraction

import extract_mzml_data as exdata
import numpy as np
import pandas as pd
import os
from utils import ppm_diff, find_nearest, _is_blank, print_progress_bar
from raw_file_validation import get_raw_file_age, check_file_collection_method
from raw_to_mzml import _raw_to_mzml
import plots
from dataset import File, filename_categories_vocab

def initialize_data_storage(current_dir):

    ms1_peak_cols = ['file_name',
                     'run_num',
                     'file_category',
                     'polarity',
                     'compound_name',
                     'retention_time',
                     'theoretical_mz',
                     'observed_mz',
                     'ppm_error',
                     'observed_intensity']

    ms1_tic_cols = ['file_name', 
                    'run_num',
                    'file_category',
                    'polarity',
                    'group',
                    'ms1_tic']

    ms2_peak_cols = ['file_name', 
                    'run_num',
                    'file_category',
                    'polarity',
                    'theoretical_mz', 
                    'observed_mz',
                    'ppm_error',
                    'observed_intensity']

    pd.DataFrame(columns=ms1_peak_cols).to_csv(os.path.join(current_dir, 'data_store/ms1_peak_data.csv'), index=False)
    pd.DataFrame(columns=ms1_tic_cols).to_csv(os.path.join(current_dir, 'data_store/ms1_tic_data.csv'), index=False)
    pd.DataFrame(columns=ms2_peak_cols).to_csv(os.path.join(current_dir, 'data_store/ms2_peak_data.csv'), index=False)


def _store_data(experiment_ms1_peak_data, experiment_ms1_tic_data, experiment_ms2_peak_data, current_dir):

    if any(experiment_ms1_peak_data):

        ms1_df = pd.DataFrame(experiment_ms1_peak_data)
        ms1_df.to_csv(os.path.join(current_dir, 'data_store/ms1_peak_data.csv'), mode='a', index=False, header=False)

    if any(experiment_ms1_tic_data):

        ms1_tic_df = pd.DataFrame(experiment_ms1_tic_data)
        ms1_tic_df.to_csv(os.path.join(current_dir, 'data_store/ms1_tic_data.csv'), mode='a', index=False, header=False)

    if any(experiment_ms2_peak_data):

        ms2_df = pd.DataFrame(experiment_ms2_peak_data)
        ms2_df.to_csv(os.path.join(current_dir, 'data_store/ms2_peak_data.csv'), mode='a', index=False, header=False)

def read_data_storage(current_dir):

    ms1_df = pd.read_csv(os.path.join(current_dir, 'data_store/ms1_peak_data.csv'))
    ms1_tic_df = pd.read_csv(os.path.join(current_dir, 'data_store/ms1_tic_data.csv'))
    ms2_df = pd.read_csv(os.path.join(current_dir, 'data_store/ms2_peak_data.csv'))

    return ms1_df, ms1_tic_df, ms2_df

def close_data_storage(current_dir):
    os.remove(os.path.join(current_dir, 'data_store/ms1_peak_data.csv'))
    os.remove(os.path.join(current_dir, 'data_store/ms1_tic_data.csv'))
    os.remove(os.path.join(current_dir, 'data_store/ms2_peak_data.csv'))

def _collect_ms1_peak_data(file, atlas_df, polarity,use_rt_filtering=True, tolerance=0.0015):

    file_peak_data = []
    for idx, row in atlas_df.iterrows():

        theoretical_mz = row['{pol}_mz'.format(pol=polarity)]

        if use_rt_filtering:
            eic = exdata._get_ms1_eic(file.ms1_data, theoretical_mz, rt=row['ideal_rt'], tolerance=tolerance)
        else:
            eic = exdata._get_ms1_eic(file.ms1_data, theoretical_mz, tolerance=tolerance)
        
        rt, mz, i = exdata._get_peak_data(eic)

        peak_data = {'file_name':file.name,
                     'run_num':file.run_num,
                     'file_category':file.category,
                     'polarity':polarity,
                     'compound_name':row['compound_name'],
                     'retention_time':rt,
                     'theoretical_mz':theoretical_mz,
                     'observed_mz':mz,
                     'ppm_error':ppm_diff(mz, theoretical_mz),
                     'observed_intensity':i}

        file_peak_data.append(peak_data)

    return file_peak_data

def _collect_ms2_null_peaks(ms2_diagnostic, file):
    not_found_peaks = []
    for mz in ms2_diagnostic[file.polarity]['diagnostic_ions']:

        not_found_peaks.append({'file_name':file.name, 
                                'run_num':file.run_num,
                                'file_category':file.category,
                                'polarity':file.polarity,
                                'theoretical_mz':mz, 
                                'observed_mz':np.nan,
                                'ppm_error':np.nan,
                                'observed_intensity':np.nan})

    return not_found_peaks

def _collect_matching_ms2_peaks(ms2_spec, ms2_diagnostic, file, frag_tolerance=0.02):

    found_peaks = []
    for mz in ms2_diagnostic[file.polarity]['diagnostic_ions']:
        idx = np.isclose(ms2_spec[0], mz, atol=frag_tolerance)

        if idx.any():

            filtered_mzs = ms2_spec[0][idx]
            filtered_is = ms2_spec[1][idx]

            nearest_mz_idx = find_nearest(filtered_mzs, mz)
            nearest_mz = filtered_mzs[nearest_mz_idx]
            nearest_mz_intensity = filtered_is[nearest_mz_idx]

            found_peaks.append({'file_name':file.name, 
                                'run_num':file.run_num,
                                'file_category':file.category,
                                'polarity':file.polarity,
                                'theoretical_mz':mz, 
                                'observed_mz':nearest_mz,
                                'ppm_error':ppm_diff(nearest_mz, mz),
                                'observed_intensity':nearest_mz_intensity})
        else:
            found_peaks.append({'file_name':file.name,
                                'run_num':file.run_num, 
                                'file_category':file.category,
                                'polarity':file.polarity,
                                'theoretical_mz':mz, 
                                'observed_mz':np.nan,
                                'ppm_error':np.nan,
                                'observed_intensity':np.nan})

    return found_peaks

def analyze_file_ms1_tic(file):

    file_tic_data = {'file_name':file.name, 
                    'run_num':file.run_num,
                    'file_category':file.category,
                    'polarity':file.polarity,
                    'group':file.group_name,
                    'ms1_tic':file.ms1_tic}

    return file_tic_data

def analyze_file_ms1_data(file, atlas_df, tolerance=0.0015):

    if file.polarity == 'FPS':
        file_peak_data = _collect_ms1_peak_data(file, atlas_df, 'POS', tolerance=tolerance) + _collect_ms1_peak_data(file, atlas_df, 'NEG', tolerance=tolerance)
    else:
        file_peak_data = _collect_ms1_peak_data(file, atlas_df, file.polarity, tolerance=tolerance)

    return file_peak_data

def analyze_file_ms2_data(file, ms2_diagnostic, tolerance=0.0015, frag_tolerance=0.02):

    file_ms2_peaks = exdata._extract_ms2_data(file.path, ms2_diagnostic[file.polarity]['pmz'], tolerance=tolerance)

    if len(file_ms2_peaks) > 0:
        ms2_intensity_sums = [ms2[1].sum() for ms2 in file_ms2_peaks]
        ms2_i_max_idx = np.argmax(ms2_intensity_sums)

        max_ms2_spec = file_ms2_peaks[ms2_i_max_idx]
        file_matching_peaks =_collect_matching_ms2_peaks(max_ms2_spec, ms2_diagnostic, file, frag_tolerance=frag_tolerance)

    else:
        file_matching_peaks = _collect_ms2_null_peaks(ms2_diagnostic, file)

    return file_matching_peaks

def analyze_experiment(dataset, atlas_df, ms2_diagnostic, current_dir, args, exported_files, tolerance=0.0015, frag_tolerance=0.005):

    experiment_ms1_peak_data = []
    experiment_ms2_peak_data = []
    experiment_ms1_tic_data = []

    analaysis_interrupted = False

    for file_path in dataset.files_to_analyze:
    
        file_age = get_raw_file_age(file_path)
        if file_age < args.min_file_age:
            continue

        file_name = os.path.basename(file_path).split('.')[0]
        if args.skip_blanks:
            is_blank = _is_blank(file_name, filename_categories_vocab)
            if is_blank:
                dataset.add_analyzed_file(file_path)
                print_progress_bar(len(dataset.analyzed_files), args.num_files, prefix = 'Progress:', suffix = 'Complete', length = 50)

                if len(dataset.analyzed_files) >= args.num_files:
                    break
                elif (len(dataset.analyzed_files) - exported_files) == args.export_num:
                    break
                else:
                    continue
        try:
            centroid_check = check_file_collection_method(file_path)
            if not centroid_check:
                print('\nData collected in {filename} is not centroid, check method!'.format(filename=os.path.basename(file_path)))
                analaysis_interrupted = True
                break
        except:
            print('\nUnable to determine collection method, check that {filename} is not truncated or corrupt'.format(filename=os.path.basename(file_path)))
            analaysis_interrupted = True
            break

        conversion_log = _raw_to_mzml(file_path)
        mzml_dir = conversion_log.args.split('-o=')[1].split(' ')[0]
        mzml_file_path = os.path.join(mzml_dir, (file_name + '.mzML'))

        file = File(mzml_file_path, file_name)

        file_ms1_data = analyze_file_ms1_data(file, atlas_df, tolerance=tolerance)
        file_ms1_tic = analyze_file_ms1_tic(file)

        experiment_ms1_tic_data.append(file_ms1_tic)
        experiment_ms1_peak_data = experiment_ms1_peak_data + file_ms1_data

        if file.ms_num == 'MS2' or file.ms_num == 'MSMS':

            file_matching_peaks = analyze_file_ms2_data(file, ms2_diagnostic, tolerance=tolerance, frag_tolerance=frag_tolerance)
            experiment_ms2_peak_data = experiment_ms2_peak_data + file_matching_peaks

        os.remove(mzml_file_path)
        dataset.add_analyzed_file(file_path)
        print_progress_bar(len(dataset.analyzed_files), args.num_files, prefix = 'Progress:', suffix = 'Complete', length = 50)

        if len(dataset.analyzed_files) >= args.num_files:
            break
        elif (len(dataset.analyzed_files) - exported_files) == args.export_num:
            break

    _store_data(experiment_ms1_peak_data, experiment_ms1_tic_data, experiment_ms2_peak_data, current_dir)

    return analaysis_interrupted

def export_tables_plots(current_dir, args, file_name_warnings_df, atlas_df, ms2_diagnostic, exported_files):

    ms1_peak_data, ms1_tic_data, ms2_peak_data = read_data_storage(current_dir)
            
    qc_output_dir = os.path.join(args.directory + '\\qc_output_{exported_files}'.format(exported_files=exported_files))

    if not os.path.isdir(qc_output_dir):
        os.mkdir(qc_output_dir)

    file_name_warnings_df.to_csv(os.path.join(qc_output_dir, 'file_name_warnings_report.csv'))

    ms1_peak_data.to_csv(os.path.join(qc_output_dir, 'ms1_data_sheet.csv'))
    plots.make_ms1_qc_plots(ms1_peak_data, atlas_df, qc_output_dir)
    plots.make_ms1_qc_plots(ms1_peak_data, atlas_df, qc_output_dir, logy=True)

    plots.make_unlabeled_ms1_qc_plots(ms1_peak_data, atlas_df, qc_output_dir)
    plots.make_unlabeled_ms1_qc_plots(ms1_peak_data, atlas_df, qc_output_dir, logy=True)

    plots.make_ms1_tic_qc_plots(ms1_tic_data, qc_output_dir, logy=False)
    plots.make_ms1_tic_qc_plots(ms1_tic_data, qc_output_dir, logy=True)

    if not ms2_peak_data.empty:
        ms2_peak_data.to_csv(os.path.join(qc_output_dir, 'ms2_data_sheet.csv'))
        plots.make_ms2_qc_plots(ms2_peak_data, ms2_diagnostic, qc_output_dir)
        plots.make_ms2_qc_plots(ms2_peak_data, ms2_diagnostic, qc_output_dir, logy=True)
