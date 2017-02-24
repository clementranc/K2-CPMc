// File written by Clement Ranc.
// This file is a small piece of the translation in C++ of the
// code K2-CPM (Wang et al., 2016).
//     https://github.com/jvc2688/K2-CPM
//==================================================================//

#ifndef __MATRIX_H_
#define __MATRIX_H_
#include "table.h"

class Matrix: public Table {
    // Constructors
    // ------------
    public:
    explicit Matrix(int);
    Matrix (const Table&);
    Matrix (const Matrix&);

    virtual ~Matrix();

    // Definitions
    // -----------
    void operator=(const Matrix&);
    void operator=(const Table&);
    void operator=(double);

    virtual Table operator*(const Table&) const;

    // Useful functions
    // ----------------
    protected:
    virtual void print(ostream&) const;

    public:
    void cholesky_factor();  // Caution: overwrite the matrix by lower Choleski factor
    void cholesky_solve(Table&);  // Caution: overwrite argument


};

#endif
