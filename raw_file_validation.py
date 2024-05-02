from fisher_py.raw_file_reader import RawFileReaderAdapter
from fisher_py.data import Device
from random import randrange
import time
import os

def get_raw_file_age(file_path):

    current_time = time.time()
    modified_time = os.path.getmtime(file_path)
    age_minutes = (current_time - modified_time) / 60
    return age_minutes

def _centroid_data_collection(raw_file, scan_number):
    scan_statistics = raw_file.get_scan_stats_for_scan_number(scan_number)
    if scan_statistics.is_centroid_scan:
        return True
    else:
        return False

def check_file_collection_method(file_path, sampling_number=10):
    try: 
        raw_file = RawFileReaderAdapter.file_factory(file_path)
        raw_file.select_instrument(Device.MS, 1)
    except:
        message = "Corrupted file detected, {}".format(file_path)
        raise Exception(message)

    centroid_results = []
    max_scan = raw_file.run_header_ex.last_spectrum

    for i in range(sampling_number):
        random_scan = randrange(1, max_scan)
        centroid_results.append(_centroid_data_collection(raw_file, random_scan))

    if not all(centroid_results):
        return False
    else:
        return True
