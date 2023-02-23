# given the temp mzMl file_path, extract ms1 and ms2 data

import pymzml
import numpy as np 
from utils import find_nearest

def _extract_ms1_data(file_path):

    run = pymzml.run.Reader(file_path)
    ms1_peaks = []

    for spec in run:
        if spec.ms_level == 1:
            ms1_peaks.append([np.repeat(spec.scan_time_in_minutes(), spec.i.shape[0]), spec.mz,  spec.i])

    run.close()

    ms1_rts = np.concatenate([rtmzi[0] for rtmzi in ms1_peaks])
    ms1_mzs = np.concatenate([rtmzi[1] for rtmzi in ms1_peaks])
    ms1_is = np.concatenate([rtmzi[2] for rtmzi in ms1_peaks])

    ms1_data = (ms1_rts, ms1_mzs, ms1_is)

    return ms1_data

def _get_ms1_eic(ms1_data, mz, tolerance=0.02, rt=None, rt_window = 2):
    
    mz_filter = np.isclose(mz, ms1_data[1], atol=tolerance)
    if rt != None:
        rt_filter = np.isclose(rt, ms1_data[0], atol=rt_window)
        eic_filter = mz_filter * rt_filter
    else:
        eic_filter = mz_filter

    ms1_eic = (ms1_data[0][eic_filter], ms1_data[1][eic_filter], ms1_data[2][eic_filter])

    return ms1_eic

def _get_peak_data(ms1_eic):

    if not any(ms1_eic[2]):
        peak_data = (np.nan, np.nan, np.nan)

    else:
        max_i_idx = np.argmax(ms1_eic[2])

        peak_rt = ms1_eic[0][max_i_idx]
        peak_mz = ms1_eic[1][max_i_idx]
        peak_i = ms1_eic[2][max_i_idx]

        peak_data = (peak_rt, peak_mz, peak_i)

    return peak_data

def _extract_ms1_tic(file_path):

    run = pymzml.run.Reader(file_path)
    ms1_tic = []
    ms1_tic_rts = []

    for spec in run:
        if spec.ms_level == 1:
            ms1_tic.append(spec.TIC)
            ms1_tic_rts.append(spec.scan_time_in_minutes())

    run.close()
        
    ms1_tic_data = (ms1_tic, ms1_tic_rts)

    return ms1_tic_data

def _extract_ms2_data(file_path, pmz, tolerance=0.02):

    run = pymzml.run.Reader(file_path)
    ms2_data = []

    for spec in run:
        if spec.ms_level == 2:
            mz = spec.selected_precursors[0]['mz']
            if np.isclose(pmz, mz, atol=tolerance):
                ms2_data.append(np.array([spec.mz, spec.i]))

    run.close()
    
    return ms2_data

