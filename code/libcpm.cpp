// File written by Clement Ranc.
// This file is a small piece of the translation in C++ of the
// code K2-CPM (Wang et al., 2016).
//     https://github.com/jvc2688/K2-CPM
//==================================================================//

#include <iostream>
#include<fstream>
#include<iomanip>
#include <vector>
#include "matrix.h"
#include "libcpm.h"

using namespace std;

//==================================================================//
// Functions
//==================================================================//
void linear_least_squares(Table* a, Table* y, const Table* yvar, const Table* l2, Table* result){
/*
    Solver of linear systems using Cholesky's decomposition. Let's define
    a matrix A (dimension n1 x n2) and two vectors X (dimension n2) and Y
    (dimension n1). This function solves the equation AX = Y. The
    following steps are considered:

    1/ Add the observational uncertainties to the data.

    2/ Compute the square matrix A^T A and the vector A^T Y.

    3/ Find the square matrix L so that A^T A = L L^T, where L^T is the
    transposition of L.

    4/ Find the solution X of the system A^T A X = L L^T X = A^T Y.

    The solver assume that the matrix A^T A is a square matrix, symmetric
    and positive-definite.

    Inputs
    ------
    a -- Pointer to Table of dimension n_data x n_predictors
        The basis matrix.
    y -- Pointer to Table with dimension n_data
        The observations.
    yvar -- Pointer to Table with dimension n_data
        The observational variance of the points y.
    l2 -- Pointer to Table with dimension n_predictors
        The L2 regularization strength.
    result -- Pointer to Table with dimension n_predictors
        The solution will be written in this Table.
*/

    // Declarations and initialization
    // -------------------------------
    int i, j, k, dim1a, dim2a, x;
    double s;

    dim1a = a->get_size1();
    dim2a = a->get_size2();
    x = y->get_size1();
    assert(dim1a == x);
    x = y->get_size2();
    assert(x == 1);
    x = y->get_size3();
    assert(x == 1);

    x = yvar->get_size1();
    assert(dim1a == x);
    x = yvar->get_size2();
    assert(x == 1);
    x = yvar->get_size3();
    assert(x == 1);

    x = l2->get_size1();
    assert(dim2a == x);
    x = l2->get_size2();
    assert(x == 1);
    x = l2->get_size3();
    assert(x == 1);

    Matrix ata(dim2a);
    Table cia(dim1a, dim2a), at(dim2a, dim1a), ciy(dim1a), b(dim2a);

    // Incorporate observational uncertainties
    // ---------------------------------------
    for(i = 0; i < dim1a; ++i){
        for(j = 0; j < dim2a; ++j){
            cia.set(i, j) = (*a)(i, j) / (*yvar)(i);
            at.set(j, i) = (*a)(i, j);  // compute transpose of a
        }
    }
    ciy = (*y) / (*yvar);

    // Compute the pre-factor
    // ----------------------
    for(i = 0; i < dim2a; ++i){
        for(j = 0; j < dim2a; ++j){
            s = 0;
            for(k = 0; k < dim1a; ++k) s += (at)(i, k) * (*a)(k, j);
            ata.set(i, j) = s;
        }
    }

    for(i = 0; i < dim2a; ++i){
        s = 0;
        for(j = 0; j < dim1a; ++j) s += (at)(i, j) * ciy(j);
        b.set(i) = s;
    }

    // Incorporate any L2 regularization
    // ---------------------------------
    for(i = 0; i < dim2a; ++i) { ata.set(i, i) += (*l2)(i); }

    // Solve the equations overwriting the matrix and tables
    // -----------------------------------------------------
    ata.cholesky_factor();
    ata.cholesky_solve(b);
    *result = b;
}
//==================================================================//
void fit_target_no_train(Table& target_flux, Table& predictor_flux_matrix,
        Table& time, Table& covar_list, Table& l2_vector, double* train_lim,
        Table& result){
/*
    Fit the fluxes of the pixels.

    Input
    -----
    target_flux -- Table, dimension n_data
        The target flux.
    predictor_flux_matrix -- Table, dimension n_data x n_predictors
        The flux of nearby stars used in the fitting process.
    time -- Table, dimension n_data
        Date of the observations (unit ?).
    covar_list -- Table, dimension n_data
        List of the standart deviation for the predictors.
    l2_vector -- Table, dimension n_predictors
        Array of L2 regularization strength.
    train_lim -- array of double, dimension 2.
        The dates between train_lim[0] and train_lim[1] are excluded from
        the fit.
    result -- Table, dimension n_predictors
        Result of the fit will be written in this Table.
*/

    // Declarations and initializations
    int i, j, n_dates, n_trainmask, n_pred;

    n_dates = time.get_size1();
    n_pred = predictor_flux_matrix.get_size2();

    // Create a mask
    if ((train_lim[0]>0) && (train_lim[1]>train_lim[0])) {
        n_trainmask = 0;
        for(i=0; i<n_dates; ++i) {
            if ((time(i)<train_lim[0]) || (time(i)>train_lim[1])) ++n_trainmask;
        }
    }
    else n_trainmask = n_dates;

    Table predictor_flux_matrix_curr(n_trainmask, n_pred), target_flux_curr(n_trainmask), covar_list_curr(n_trainmask);
    if(n_dates == n_trainmask) {
        predictor_flux_matrix_curr = predictor_flux_matrix;
        target_flux_curr = target_flux;
        covar_list_curr = covar_list;
    }
    else {
        for(i=0; i<n_trainmask; ++i) {
            if ((time(i)<train_lim[0]) || (time(i)>train_lim[1])) {
                for(j=0; j<n_pred; ++j) {
                    predictor_flux_matrix_curr.set(i, j) = predictor_flux_matrix(i, j);
                }
                target_flux_curr.set(i) = target_flux(i);
                covar_list_curr.set(i) = covar_list(i);
            }
        }
    }

    covar_list_curr = pow(covar_list_curr, 2);
    linear_least_squares(&predictor_flux_matrix_curr, &target_flux_curr, &covar_list_curr, &l2_vector, &result);

}
//==================================================================//
vector< vector<double> > linear_least_squares_wrap(vector< vector<double> > input_matrix){
/*
    This function will be removed. It allows to test the routines from Python.

    The input comes from Python and includes real tables from Dun's routine.
    This function is only used to test the C/C++/Python results.
*/

    vector< vector<double> > result_global;

    // ----------------------------------------
    // To test the function fit_target_no_train
    // ----------------------------------------

    // Declaration and initialisations
    int i, j, n_dates, n_pred, n_covar, n_l2;
    double train_lim[2];
    vector<double> result_fit_target_no_train;

    train_lim[0] = input_matrix[0][0];
    train_lim[1] = input_matrix[0][1];

    n_dates = input_matrix[1].size();
    Table time(n_dates), target_flux(n_dates);
    for(i = 0; i < n_dates; ++i) {
        time.set(i) = input_matrix[4][i];
        target_flux.set(i) = input_matrix[1][i];
    }

    n_pred = input_matrix[2].size() / n_dates;
    Table predictor_flux_matrix(n_dates, n_pred);
    for(i=0; i<n_dates; ++i){
        for(j=0; j<n_pred; ++j){
            predictor_flux_matrix.set(i, j) = input_matrix[2][i * n_pred + j];
        }
    }

    n_covar = input_matrix[3].size();
    Table covar_list(n_dates);
    if (n_covar != n_dates) covar_list = input_matrix[3][0];
    else for(i = 0; i < n_dates; ++i) covar_list.set(i) = input_matrix[3][i];

    n_l2 = input_matrix[5].size();
    Table l2_vector(n_pred);
    if (n_l2 != n_pred) l2_vector = input_matrix[5][0];
    else for(i = 0; i < n_pred; ++i) l2_vector.set(i) = input_matrix[5][i];

    Table result(n_pred);
    result = 0.0;
    fit_target_no_train(target_flux, predictor_flux_matrix, time, covar_list, l2_vector, train_lim, result);

    for(i=0; i<n_pred; ++i) result_fit_target_no_train.push_back(result(i));

    result_global.push_back(result_fit_target_no_train);

    // ----------------------------------------
    // To test the solver using Choleski method
    // ----------------------------------------
    Matrix * m_test;
    m_test = new Matrix(3);
    m_test->set(0,0) = 25;
    m_test->set(0,1) = 15;
    m_test->set(0,2) = -5;
    m_test->set(1,0) = 15;
    m_test->set(1,1) = 18;
    m_test->set(1,2) = 0;
    m_test->set(2,0) = -5;
    m_test->set(2,1) = 0;
    m_test->set(2,2) = 11;

    Table * y_test;
    y_test = new Table(3);
    y_test->set(0) = -2;
    y_test->set(1) = 0.25;
    y_test->set(2) = 10.2;
//    cout << endl << "y_test --> " << *y_test;

//    cout << (*m_test);
//    m_test->cholesky_factor();
//    cout << (*m_test);
//    cout << (*y_test);
//    m_test->cholesky_solve(*y_test);
//    cout << (*y_test);

    delete m_test;
    delete y_test;

    // ----------------
    // Return to python
    // ----------------
    return result_global;
}
