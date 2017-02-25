# -*-coding:Utf-8 -*
# ====================================================================
# This routine tests the C++ routines and compare results with the
# python version of the code.
# ====================================================================
# Standard packages
# ====================================================================
import sys
import pickle
import numpy as np
# ====================================================================
# Non-standard packages
# ====================================================================
import libcpm
# ====================================================================
#   Functions
# ====================================================================
def run_test(**kwargs):

    # Load data from Python version, for comparison
    file = open("temp/train_lim4c", "r")
    train_lim = pickle.load(file)
    file.close()

    if train_lim == None: train_lim = [-1, -1]

    file = open("temp/predictor_flux_matrix4c", "r")
    predictor_flux_matrix = pickle.load(file)
    file.close()

    file = open("temp/target_flux4c", "r")
    target_flux = pickle.load(file)
    file.close()

    file = open("temp/covar_list4c", "r")
    covar_list = pickle.load(file)
    file.close()

    file = open("temp/time4c", "r")
    time = pickle.load(file)
    file.close()

    file = open("temp/result4c", "r")
    result = pickle.load(file)
    file.close()

    if covar_list == None: covar_list = [1.0]

    l2vector = [1000]

    predictor_flux_matrix_temp = np.reshape(predictor_flux_matrix, predictor_flux_matrix.shape[0]*predictor_flux_matrix.shape[1])

    # Run C/C++ code
    list = [train_lim, target_flux.tolist(), predictor_flux_matrix_temp.tolist(), covar_list, time.tolist(), l2vector]
    result_c = libcpm.linear_least_squares_wrap(list)

    # Compare the results
    result_c = np.array(result_c[0])
    result_c = np.atleast_2d(result_c).T
    np.testing.assert_almost_equal(result_c, result, decimal=7)


# ====================================================================
# Main function
# ====================================================================
if __name__=="__main__":
    run_test()










