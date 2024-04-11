import pymzml
import numpy as np 
import os
from utils import find_nearest
from extract_mzml_data import _extract_ms1_data, _get_ms1_eic
import pytest
# test for _extract_ms1_data

def test_return_value_equals_real_value(request):

    rootdir = request.config.rootdir

    real_rts = np.load(os.path.join(rootdir, "tests/test_data/test_file_rts.npy"))
    real_is = np.load(os.path.join(rootdir, "tests/test_data/test_file_is.npy"))
    real_mzs = np.load(os.path.join(rootdir, "tests/test_data/test_file_mzs.npy"))

    tes_ms1_rts, tes_ms1_mzs, tes_ms1_is  = _extract_ms1_data("./tests/test_data/20210929_JGI-AK-TH_DS_503916_WRFungi_final_QE-HF_C18_USDAY63680_FPS_MS1_31-P_CsAWO-pellet_3_Rg80to1200-CE102040-Csubverm-QC-C18QU_Run248.mzML")
    
   

    assert np.allclose(tes_ms1_rts, real_rts, atol=1e-6)
    assert np.allclose(tes_ms1_is, real_is, atol=1e-6)
    assert np.allclose(tes_ms1_mzs, real_mzs, atol=1e-6)

def test_no_negatives(request):

    rootdir = request.config.rootdir
   
    ms1_data = _extract_ms1_data("./tests/test_data/20210929_JGI-AK-TH_DS_503916_WRFungi_final_QE-HF_C18_USDAY63680_FPS_MS1_31-P_CsAWO-pellet_3_Rg80to1200-CE102040-Csubverm-QC-C18QU_Run248.mzML")

    ms1_rts = np.load(os.path.join(rootdir, "tests/test_data/test_file_rts.npy"))
    ms1_is = np.load(os.path.join(rootdir, "tests/test_data/test_file_is.npy"))
    ms1_mzs = np.load(os.path.join(rootdir, "tests/test_data/test_file_mzs.npy"))

    assert (ms1_rts >= 0).all()
    assert (ms1_is >= 0).all()
    assert (ms1_mzs >= 0).all()

def test_same_size(request):
    rootdir = request.config.rootdir
   
    ms1_data = _extract_ms1_data("./tests/test_data/20210929_JGI-AK-TH_DS_503916_WRFungi_final_QE-HF_C18_USDAY63680_FPS_MS1_31-P_CsAWO-pellet_3_Rg80to1200-CE102040-Csubverm-QC-C18QU_Run248.mzML")

    ms1_rts = np.load(os.path.join(rootdir, "tests/test_data/test_file_rts.npy"))
    ms1_is = np.load(os.path.join(rootdir, "tests/test_data/test_file_is.npy"))
    ms1_mzs = np.load(os.path.join(rootdir, "tests/test_data/test_file_mzs.npy"))

    assert np.shape(ms1_rts) == np.shape(ms1_is) == np.shape(ms1_mzs)
  


    #test for _get_ms1_eic
    
def test_eic_return_value_equals_real_value(request):

    rootdir = request.config.rootdir

    real_vals1 = np.load(os.path.join(rootdir, "tests/test_data/test_file_ms1_eic1.npz"))["arr_0"]
    real_vals2 = np.load(os.path.join(rootdir, "tests/test_data/test_file_ms1_eic2.npz"))["arr_0"]
    real_vals3 = np.load(os.path.join(rootdir, "tests/test_data/test_file_ms1_eic3.npz"))["arr_0"]

    test_ms1_data = np.load(os.path.join(rootdir, "tests/test_data/test_ms1_data_file.npz"))["arr_0"]

    test_vals1, test_vals2, test_vals3 = _get_ms1_eic(test_ms1_data, 176.1135)
    
    assert np.allclose(test_vals1, real_vals1, atol=1e-6)
    assert np.allclose(test_vals2, real_vals2, atol=1e-6)
    assert np.allclose(test_vals3, real_vals3, atol=1e-6)

    '''
def test_eic_no_negatives(request):

    rootdir = request.config.rootdir
   
    ms1_data = _extract_ms1_data("./tests/test_data/20210929_JGI-AK-TH_DS_503916_WRFungi_final_QE-HF_C18_USDAY63680_FPS_MS1_31-P_CsAWO-pellet_3_Rg80to1200-CE102040-Csubverm-QC-C18QU_Run248.mzML")

    ms1_rts = np.load(os.path.join(rootdir, "tests/test_data/test_file_rts.npy"))
    ms1_is = np.load(os.path.join(rootdir, "tests/test_data/test_file_is.npy"))
    ms1_mzs = np.load(os.path.join(rootdir, "tests/test_data/test_file_mzs.npy"))

    assert (ms1_rts >= 0).all()
    assert (ms1_is >= 0).all()
    assert (ms1_mzs >= 0).all()

def test_eic_same_size(request):
    rootdir = request.config.rootdir
   
    ms1_data = _extract_ms1_data("./tests/test_data/20210929_JGI-AK-TH_DS_503916_WRFungi_final_QE-HF_C18_USDAY63680_FPS_MS1_31-P_CsAWO-pellet_3_Rg80to1200-CE102040-Csubverm-QC-C18QU_Run248.mzML")

    ms1_rts = np.load(os.path.join(rootdir, "tests/test_data/test_file_rts.npy"))
    ms1_is = np.load(os.path.join(rootdir, "tests/test_data/test_file_is.npy"))
    ms1_mzs = np.load(os.path.join(rootdir, "tests/test_data/test_file_mzs.npy"))

    tolerance = 1e-6
    #how do i adjust the tolerance
    assert np.all(np.abs(np.diff(ms1_rts)) <= tolerance)
    assert np.all(np.abs(np.diff(ms1_is)) <= tolerance)
    assert np.all(np.abs(np.diff(ms1_mzs)) <= tolerance)'''