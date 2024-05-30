
import numpy as np
import pandas as pd
import os
import pytest
from analysis_routines import initialize_data_storage, _store_data, read_data_storage, close_data_storage, _collect_ms1_peak_data
import sys


# Initialize data storage Tests

def test_file_creation(request):
    rootdir = request.config.rootdir
    initialize_data_storage(rootdir)

   # Check if the CSV files were created
    
    assert os.path.isfile(os.path.join(rootdir, 'data_store/ms1_peak_data.csv'))
    assert os.path.isfile(os.path.join(rootdir, 'data_store/ms1_tic_data.csv'))
    assert os.path.isfile(os.path.join(rootdir, 'data_store/ms2_peak_data.csv'))


def test_column_names(request):
    rootdir = request.config.rootdir
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

    
    # Read the CSV files and check column names
    ms1_peak_df = pd.read_csv(os.path.join(rootdir, 'data_store/ms1_peak_data.csv'))
    ms1_tic_df = pd.read_csv(os.path.join(rootdir, 'data_store/ms1_tic_data.csv'))
    ms2_peak_df = pd.read_csv(os.path.join(rootdir, 'data_store/ms2_peak_data.csv'))

    assert list(ms1_peak_df.columns) == ms1_peak_cols
    assert list(ms1_tic_df.columns) == ms1_tic_cols
    assert list(ms2_peak_df.columns) == ms2_peak_cols


# test Store data


def test_file_creation_and_appending(request):
    rootdir = request.config.rootdir
    # Sample data
    #higher_dir = os.path.dirname(rootdir)
    experiment_ms1_peak_data = [{'col1': 1, 'col2': 'A'}, {'col1': 2, 'col2': 'B'}]
    experiment_ms1_tic_data = [{'col1': 3, 'col2': 'C'}, {'col1': 4, 'col2': 'D'}]
    experiment_ms2_peak_data = [{'col1': 5, 'col2': 'E'}, {'col1': 6, 'col2': 'F'}]

 
    _store_data(experiment_ms1_peak_data, experiment_ms1_tic_data, experiment_ms2_peak_data, rootdir)

    # Check if the CSV files were created and contain the correct data
    ms1_peak_df = pd.read_csv(os.path.join(rootdir, 'data_store/ms1_peak_data.csv'))
    ms1_tic_df = pd.read_csv(os.path.join(rootdir, 'data_store/ms1_tic_data.csv'))
    ms2_df = pd.read_csv(os.path.join(rootdir, 'data_store/ms2_peak_data.csv'))

    assert not ms1_peak_df.empty
    assert not ms1_tic_df.empty
    assert not ms2_df.empty


def test_handling_empty_data(request):
   
    rootdir = request.config.rootdir
    _store_data([], [], [], rootdir)

    # Check if the CSV files were created 
    assert os.path.isfile(os.path.join(rootdir, 'data_store/ms1_peak_data.csv'))
    assert os.path.isfile(os.path.join(rootdir, 'data_store/ms1_tic_data.csv'))
    assert os.path.isfile(os.path.join(rootdir, 'data_store/ms2_peak_data.csv'))

# Test Read data storage

def test_read_data_storage(request):

    rootdir = request.config.rootdir

    ms1_df, ms1_tic_df, ms2_df = read_data_storage(rootdir)

    assert ms1_df.equals(pd.read_csv(os.path.join(rootdir, 'data_store/ms1_peak_data.csv')))
    assert ms1_tic_df.equals(pd.read_csv(os.path.join(rootdir, 'data_store/ms1_tic_data.csv')))
    assert ms2_df.equals(pd.read_csv(os.path.join(rootdir, 'data_store/ms2_peak_data.csv')))

#Test close data storage
    
def test_close_data_storage(request):

    rootdir = request.config.rootdir 
    
    close_data_storage(rootdir)

    assert not os.path.exists(os.path.join(rootdir, 'data_store/ms1_peak_data.csv'))
    assert not os.path.exists(os.path.join(rootdir, 'data_store/ms1_tic_data.csv'))
    assert not os.path.exists(os.path.join(rootdir, 'data_store/ms2_peak_data.csv'))

#Test _collect_ms1_peak_data work in progress
    
'''# Sample data for testing
sample_file = ...  # Define your sample file object
sample_atlas_df = ...  # Define your sample atlas DataFrame
def test__collect_ms1_peak_data()
    
    @pytest.fixture
    def _get_ms1_eic()
        return "tests/test_data/test_file_ms1_eic.npz"
    
    @pytest.fixture
    def _get_peak_data()
        test_peak = "tests/test_data/test_file_ms1_peak0.npz" + "tests/test_data/test_file_ms1_peak1.npz" + "tests/test_data/test_file_ms1_peak2.npz"
        return test_peak

# Test peak data extraction
def test_peak_data_extraction():
    peak_data = _collect_ms1_peak_data(sample_file, sample_atlas_df, polarity='positive')
    assert isinstance(peak_data, list)
    assert all(isinstance(entry, dict) for entry in peak_data)
    # Add more specific assertions about peak data extraction

# Test peak data attributes
def test_peak_data_attributes():
    peak_data = _collect_ms1_peak_data(sample_file, sample_atlas_df, polarity='positive')
    assert all(attr in peak_data[0] for attr in ['file_name', 'run_num', 'file_category', 'polarity',
                                                 'compound_name', 'retention_time', 'theoretical_mz',
                                                 'observed_mz', 'ppm_error', 'observed_intensity'])

# Test polarity handling
def test_polarity_handling():
    peak_data = _collect_ms1_peak_data(sample_file, sample_atlas_df, polarity='positive')
    assert all(entry['polarity'] == 'positive' for entry in peak_data)

# Test RT filtering
def test_rt_filtering():
    # Assuming you have a test file and atlas with known peak RTs for positive polarity
    # You can create a test case where you know which peaks should be filtered in or out
    # Then compare the extracted peak data with the expected results

# Test tolerance handling
#def test_tolerance_handling():
    # Similar to RT filtering, you can create test cases with known peaks within/outside tolerance
    # And compare the extracted peak data with the expected results

# Test error handling
# You can write test cases to check how the function behaves with invalid input data or edge cases

# Test performance
# You can create large datasets and measure the execution time of the function
'''