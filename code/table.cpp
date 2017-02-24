// File written by Clement Ranc.
// This file is a small piece of the translation in C++ of the
// code K2-CPM (Wang et al., 2016).
//     https://github.com/jvc2688/K2-CPM
//==================================================================//

#include<fstream>
#include<iomanip>
#include<cmath>
#include "table.h"

//==================================================================//
Table::Table(int d1, int d2, int d3) : size1(d1),
        size2(d2), size3(d3), tab(0x0) {
    assert((d1>0) && (d2>0) && (d3>0));
    tab = new double [d1 * d2 * d3];
}
//==================================================================//
Table::Table(const Table& table_in) : size1(table_in.size1),
        size2(table_in.size2), size3(table_in.size3), tab(0x0) {

    assert((size1>0) && (size2>0) && (size3>0));
    int size = size1 * size2 * size3;
    tab = new double [size];
    assert(table_in.tab != 0x0);
    for (int i=0; i<size; i++) tab[i] = table_in.tab[i];
}
//==================================================================//
Table::~Table() {
  if (tab != 0x0) delete[] tab ;
}
//==================================================================//
void Table::operator=(const Table& table_in) {
    assert(size1 == table_in.size1);
    assert(size2 == table_in.size2);
    assert(size3 == table_in.size3);
    assert(table_in.tab != 0x0);
    assert(tab != 0x0);

    int size = size1 * size2 * size3;
    for (int i=0; i<size; i++) tab[i] = table_in.tab[i];
}
 //==================================================================//
void Table::operator=(double x) {
    assert(tab != 0x0);

    int size = size1 * size2 * size3;
    for (int i=0; i<size; i++) tab[i] = x;
}
//==================================================================//
void Table::print(ostream& ost) const {
    assert(tab != 0x0);

    int tabdim = 0;
    if (size1 > 1) tabdim++;
    if (size2 > 1) tabdim++;
    if (size3 > 1) tabdim++;

    ost << "Table of ";
    if (size1 > 1) {
        ost << size1;
        if (tabdim>1) {
            ost << "x";
            tabdim--;
        }
    }
    if (size2>1) {
        ost << size2;
        if (tabdim>1) ost << "x";
    }
    if (size3>1) ost << size3;
    ost << " elements" << endl;

    ost << setprecision(3);

    if (size3 == 1) {
        for (int i=0; i<size1; i++) {
            for (int j=0; j<size2; j++) {
                ost << (*this)(i, j) << '\t';
	        }
            if (size2 >1) ost << endl;
        }
        ost << endl;
    }
    else {
        for (int i=0; i<size1; i++) {
            ost << "i=" << i << '\n';
            for (int j=0; j<size2; j++) {
                for (int k=0; k<size3; k++) {
                    ost << (*this)(i, j, k) << '\t';
                }
                ost << endl;
            }
            ost << endl ;
        }
        ost << endl ;
    }
}
//==================================================================//
ostream& operator<<(ostream& ost, const Table& table_in ) {
    assert(table_in.tab != 0x0) ;
    table_in.print(ost) ;
    return ost;
}
//==================================================================//
Table operator-(const Table& table_in) {
    int s1 = table_in.get_size1();
    int s2 = table_in.get_size2();
    int s3 = table_in.get_size3();
    int size = s1 * s2 * s3;

    Table result(s1, s2, s3);
    for (int i=0; i<size; i++) result.tab[i] = -table_in.tab[i];
    return result;
}
//==================================================================//
Table operator+(const Table& t1, const Table& t2) {
    int s1 = t1.get_size1();
    int s2 = t1.get_size2();
    int s3 = t1.get_size3();
    assert ((t2.get_size1() == s1) && (t2.get_size2() == s2)
      && (t2.get_size3() == s3));
    int size = s1 * s2 * s3;

    Table result(s1, s2, s3);
    for (int i=0; i<size; i++)
        result.tab[i]  = t1.tab[i] + t2.tab[i];
    return result;
}
//==================================================================//
Table operator+(const Table& t1, double x) {
    int s1 = t1.get_size1();
    int s2 = t1.get_size2();
    int s3 = t1.get_size3();

    Table result(s1, s2, s3);
    result = x;
    return t1 + result;
}
//==================================================================//
Table operator+(double x, const Table& t1) {
    return t1 + x;
}
//==================================================================//
Table operator-(const Table& t1, const Table& t2) {
    return t1 + (-t2);
}
//==================================================================//
Table operator-(const Table& t1, double x) {
    return t1 + (-x);
}
//==================================================================//
Table operator-(double x, const Table& t1) {
    return x + (-t1);
}
//==================================================================//
Table operator*(const Table& t1, double x) {
    int s1 = t1.get_size1();
    int s2 = t1.get_size2();
    int s3 = t1.get_size3();
    int size = s1 * s2 * s3;

    Table result(s1, s2, s3);
    for (int i=0; i<size; i++) result.tab[i] = t1.tab[i] * x;
    return result;
}
//==================================================================//
Table operator*(double x, const Table& t1) {
    return t1 * x;
}
//==================================================================//
Table operator/(const Table& t1, const Table& t2) {
    int s1 = t1.get_size1();
    int s2 = t1.get_size2();
    int s3 = t1.get_size3();
    assert ((t2.get_size1() == s1) && (t2.get_size2() == s2)
      && (t2.get_size3() == s3));
    int size = s1 * s2 * s3;

    Table result(s1, s2, s3);
    for (int i=0; i<size; i++) result.tab[i]  = t1.tab[i] / t2.tab[i];
    return result;
}
//==================================================================//
Table operator/(const Table& t1, double x) {
    assert(x > 1e-16);
    return t1 * (1.0 / x);
}
//==================================================================//
Table operator/(double x, const Table& t1) {
    int s1 = t1.get_size1();
    int s2 = t1.get_size2();
    int s3 = t1.get_size3();

    Table result(s1, s2, s3) ;
    result = x;
    return result/t1;
}
//==================================================================//
Table sqrt(const Table& table_in) {
    int s1 = table_in.get_size1();
    int s2 = table_in.get_size2();
    int s3 = table_in.get_size3();
    int size = s1 * s2 * s3;

    Table result(s1, s2, s3);
    for (int i=0; i<size; i++) result.tab[i] = sqrt(table_in.tab[i]);
    return result;
}
//==================================================================//
Table pow(const Table& table_in, double x) {
  int s1 = table_in.get_size1();
  int s2 = table_in.get_size2();
  int s3 = table_in.get_size3();
  int size = s1 * s2 * s3;

  Table result(s1, s2, s3);
  for (int i=0; i<size; i++) result.tab[i] = pow(table_in.tab[i], x);
  return result;
}
