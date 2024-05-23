#ifndef QR_HPP
#define QR_HPP

#include <iostream>
#include <fstream>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/lu.hpp>
#include <boost/numeric/ublas/triangular.hpp>

namespace ublas = boost::numeric::ublas;

void qr(ublas::matrix<double> &A, ublas::matrix<double> &Q);
void small_qr(ublas::matrix<double> &A, ublas::matrix<double> &Q, int n);
void householder(ublas::matrix<double> &Q, ublas::vector<double> v, int n);
double cnum(ublas::matrix<double> &A);
inline double sign(double x);
ublas::matrix<double> cut(ublas::matrix<double> &A, int n);
void embed(ublas::matrix<double> &A, ublas::matrix<double> &M, int n);

#endif
