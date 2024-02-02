
import numpy as np
import pandas as pd
import os
import pytest
from analysis_routines import initialize_data_storage, _store_data


# Initialize data storage Tests


'''@pytest.fixture
def whatToName_dir(tmp_path):
    return tmp_path
    # Create the "data_store" subdirectory
'''
@pytest.fixture
def current_dir(tmp_path):
    # Create the "data_store" subdirectory within the temporary directory
    data_store_dir = tmp_path / "data_store"
    data_store_dir.mkdir()
    return data_store_dir

def test_file_creation(current_dir):
    print ("hello")
    #assert os.path.isdir(whatToName_dir)
    #fun = "trythis"
    higher_dir = os.path.dirname(current_dir)
    initialize_data_storage(higher_dir)

   # Check if the CSV files were created
    assert os.path.isfile(os.path.join(current_dir, 'ms1_peak_data.csv'))
    assert os.path.isfile(os.path.join(current_dir, 'ms1_tic_data.csv'))
    assert os.path.isfile(os.path.join(current_dir, 'ms2_peak_data.csv'))
'''
def test_column_names(current_dir):
    higher_dir = os.path.dirname(current_dir)
    ms1_peak_cols, ms1_tic_cols, ms2_peak_cols =  initialize_data_storage(higher_dir)

    # Read the CSV files and check column names
    ms1_peak_df = pd.read_csv(os.path.join(current_dir, 'data_store/ms1_peak_data.csv'))
    ms1_tic_df = pd.read_csv(os.path.join(current_dir, 'data_store/ms1_tic_data.csv'))
    ms2_peak_df = pd.read_csv(os.path.join(current_dir, 'data_store/ms2_peak_data.csv'))

    assert list(ms1_peak_df.columns) == ms1_peak_cols
    assert list(ms1_tic_df.columns) == ms1_tic_cols
    assert list(ms2_peak_df.columns) == ms2_peak_cols


# test Store data


def test_file_creation_and_appending(current_dir):
    # Sample data
    higher_dir = os.path.dirname(current_dir)
    experiment_ms1_peak_data = [{'col1': 1, 'col2': 'A'}, {'col1': 2, 'col2': 'B'}]
    experiment_ms1_tic_data = [{'col1': 3, 'col2': 'C'}, {'col1': 4, 'col2': 'D'}]
    experiment_ms2_peak_data = [{'col1': 5, 'col2': 'E'}, {'col1': 6, 'col2': 'F'}]

 
    _store_data(experiment_ms1_peak_data, experiment_ms1_tic_data, experiment_ms2_peak_data, higher_dir)

    # Check if the CSV files were created and contain the correct data
    ms1_peak_df = pd.read_csv(os.path.join(current_dir, 'data_store/ms1_peak_data.csv'))
    ms1_tic_df = pd.read_csv(os.path.join(current_dir, 'data_store/ms1_tic_data.csv'))
    ms2_df = pd.read_csv(os.path.join(current_dir, 'data_store/ms2_peak_data.csv'))

    assert not ms1_peak_df.empty
    assert not ms1_tic_df.empty
    assert not ms2_df.empty

    # Add more specific checks as needed

def test_handling_empty_data(current_dir):
    # empty data
    higher_dir = os.path.dirname(current_dir)
    _store_data([], [], [], higher_dir)

    # Check if the CSV files were created (they should be empty)
    assert os.path.isfile(os.path.join(current_dir, 'data_store/ms1_peak_data.csv'))
    assert os.path.isfile(os.path.join(current_dir, 'data_store/ms1_tic_data.csv'))
    assert os.path.isfile(os.path.join(current_dir, 'data_store/ms2_peak_data.csv'))

    # Test Read data storage

def test_read_data_storage(current_dir):

    read_data_storage(current_dir)
'''



