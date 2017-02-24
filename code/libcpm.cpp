// File written by Clement Ranc.
// This file is a small piece of the translation in C++ of the
// code K2-CPM (Wang et al., 2016).
//     https://github.com/jvc2688/K2-CPM
//==================================================================//

#include <iostream>
#include <vector>
#include "matrix.h"
#include "libcpm.h"

using namespace std;

//==================================================================//
// Functions
//==================================================================//
void linear_least_squares(Matrix* a, Table* y, const Table* yvar, const Table* l2){
/*
    Solver of linear systems using Cholesky's decomposition. Let's define
    a matrix A (dimension n x n) and two vectors X and Y (each one of
    dimension n). The solver executes the following tasks:

    1/ Find the matrix L (dimension n x n) so that A = L L^T, where L^T is
    the transposition of L.

    2/ Find the solution X of the system AX = L L^T X = Y.

    The solver assume that the matrix A is a square matrix, symmetric and
    positive-definite.

    Inputs
    ------
    a -- Pointer to Matrix of dimension n x n
        The basis matrix, with n the number of data.
    y -- Pointer to Table with dimension n
        The observations.
    yvar -- Pointer to Table with dimension n
        The observational variance of the points y.
    l2 -- Pointer to Table with dimension n
        The L2 regularization strength.

    Return: N/A
    -------
        Overwrite the Table y by the solution X.
*/

    // Declarations and initialization
    // -------------------------------
    int i, j, dim1a, x;

    dim1a = a->get_size1();
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
    assert(dim1a == x);
    x = l2->get_size2();
    assert(x == 1);
    x = l2->get_size3();
    assert(x == 1);

    Matrix cia(dim1a), at(dim1a), ata(dim1a);
    Table ciy(dim1a), b(dim1a);

    // Incorporate observational uncertainties
    // ---------------------------------------
    for(i = 0; i < dim1a; ++i){
        for(j = 0; j < dim1a; ++j){
            cia.set(i, j) = (*a)(i, j) / (*yvar)(i);
            at.set(i, j) = (*a)(j, i);  // compute transpose of a
        }
    }
    ciy = (*y) / (*yvar);

    // Compute the pre-factor
    // ----------------------
    ata = at * cia;

    // Incorporate any L2 regularization
    // ---------------------------------
    for(i = 0; i < dim1a; ++i) { ata.set(i, i) += (*l2)(i); }

    // Solve the equations overwriting the matrix and tables
    // -----------------------------------------------------
    ata.cholesky_factor();
    b = at * (*y);
    ata.cholesky_solve(b);
    *y = b;
}
//==================================================================//
void linear_least_squares(Matrix* a, Table* y, const double yvar, const double l2){
/*
    Solver of linear systems using Cholesky's decomposition. Let's define
    a matrix A (dimension n x n) and two vectors X and Y (each one of
    dimension n). The solver executes the following tasks:

    1/ Find the matrix L (dimension n x n) so that A = L L^T, where L^T is
    the transposition of L.

    2/ Find the solution X of the system AX = L L^T X = Y.

    The solver assume that the matrix A is a square matrix, symmetric and
    positive-definite.

    Inputs
    ------
    a -- Pointer to Matrix of dimension n x n
        The basis matrix, with n the number of data.
    y -- Pointer to Table with dimension n
        The observations.
    yvar -- double
        The observational variance of the points y.
    l2 -- double
        The L2 regularization strength.

    Return: N/A
    -------
        Overwrite the Table y by the solution X.
*/

    // Declarations and initialization
    // -------------------------------
    int dim1a;
    dim1a = a->get_size1();
    Table tab(dim1a), tabl2(dim1a);

    tab = yvar;
    tabl2 = l2;

    // Solve the system
    // ----------------
    linear_least_squares(a, y, &tab, &tabl2);
}
//==================================================================//
void linear_least_squares(Matrix* a, Table* y, const double yvar, const Table* l2){
/*
    Solver of linear systems using Cholesky's decomposition. Let's define
    a matrix A (dimension n x n) and two vectors X and Y (each one of
    dimension n). The solver executes the following tasks:

    1/ Find the matrix L (dimension n x n) so that A = L L^T, where L^T is
    the transposition of L.

    2/ Find the solution X of the system AX = L L^T X = Y.

    The solver assume that the matrix A is a square matrix, symmetric and
    positive-definite.

    Inputs
    ------
    a -- Pointer to Matrix of dimension n x n
        The basis matrix, with n the number of data.
    y -- Pointer to Table with dimension n
        The observations.
    yvar -- double
        The observational variance of the points y.
    l2 -- Pointer to Table with dimension n
        The L2 regularization strength.

    Return: N/A
    -------
        Overwrite the Table y by the solution X.
*/

    // Declarations and initialization
    // -------------------------------
    int dim1a;
    dim1a = a->get_size1();
    Table tab(dim1a);

    tab = yvar;

    // Solve the system
    // ----------------
    linear_least_squares(a, y, &tab, l2);
}
//==================================================================//
void linear_least_squares(Matrix* a, Table* y, const Table* yvar, const double l2){
/*
    Solver of linear systems using Cholesky's decomposition. Let's define
    a matrix A (dimension n x n) and two vectors X and Y (each one of
    dimension n). The solver executes the following tasks:

    1/ Find the matrix L (dimension n x n) so that A = L L^T, where L^T is
    the transposition of L.

    2/ Find the solution X of the system AX = L L^T X = Y.

    The solver assume that the matrix A is a square matrix, symmetric and
    positive-definite.

    Inputs
    ------
    a -- Pointer to Matrix of dimension n x n
        The basis matrix, with n the number of data.
    y -- Pointer to Table with dimension n
        The observations.
    yvar -- Pointer to Table with dimension n
        The observational variance of the points y.
    l2 -- double
        The L2 regularization strength.

    Return: N/A
    -------
        Overwrite the Table y by the solution X.
*/

    // Declarations and initialization
    // -------------------------------
    int dim1a;
    dim1a = a->get_size1();
    Table tabl2(dim1a);

    tabl2 = l2;

    // Solve the system
    // ----------------
    linear_least_squares(a, y, yvar, &tabl2);
}
//==================================================================//
vector< vector<double> > linear_least_squares_wrap(vector< vector<double> > input_matrix){
/*
    This function will be removed. It allows to test the routines and
    gives an interface with python.

    For now, nothing from the input is used.
*/
    // Declarations
    // ------------
    Matrix * m_test;
    Table * y_test;
    vector< vector<double> > result;

    // Initializations
    // ---------------
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

    y_test = new Table(3);
    y_test->set(0) = -2;
    y_test->set(1) = 0.25;
    y_test->set(2) = 10.2;
    // cout << endl << "y_test --> " << *y_test;

    // Least square solver
    // -------------------
    linear_least_squares(m_test, y_test, 1.0, 1.0);

    // Release memory
    //---------------
    delete m_test;
    delete y_test;

    return result;
}











