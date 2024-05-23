#include <iostream>
#include <fstream>
#include <cmath>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include "utility.hpp"

namespace ublas = boost::numeric::ublas;

// double R(double *a, double *x) {
  // double r = 0;
  // for (int i = 0; i < 2; ++i) {
    // r += (a[i]-x[i])*(a[i]-x[i]);
  // }
  // return sqrt(r);
// }

bool check(int K, double **old_X, double *old_M, double **new_X, double *new_M) {
  bool b = false;
  double O[2]{0};
  for (int k = 0; k < K; ++k) {
    if (R(O, new_X[k]) < 2) {
      if (abs(old_M[k]-new_M[k]) > 1.E-4) {
        b = true;
        break;
      }
      if (abs(old_X[k][0]-new_X[k][0]) > 1.E-4) {
        b = true;
        break;
      }
      if (abs(old_X[k][1]-new_X[k][1]) > 1.E-4) {
        b = true;
        break;
      }
    }
  }
  return b;
}

// void update(int K, double **X, double *M, ublas::vector<double> dXM) {
  // for (int k = 0; k < K; ++k) {
    // X[k][0] += dXM(3*k+0);
    // X[k][1] += dXM(3*k+1);
    // // X[4*k][2] += dXM[4*k+2];
    // M[k] += dXM(3*k+2);
  // }
// }

std::ostringstream XM(int K, double **X, double *M) {
  std::ostringstream oss;
  for (int k = 0; k < K; ++k) {
    oss << X[k][0] << " ";
    oss << X[k][1] << " ";
    oss << M[k] << '\n';
  }
  return oss;
}

double total_mass(int K, double *M) {
  double v = 0;
  for (int k = 0; k < K; ++k) {
    v += M[k];
  }
  return v;
}

double partial_mass(int K, double **X, double *M) {
  double v = 0;
  double O[2]{0};
  for (int k = 0; k < K; ++k) {
    if (R(O, X[k]) < 2) v += M[k];
  }
  return v;
}

// std::ostream &operator<<(std::ostream &os, position &XM) {
  // XM.print(&os);
  // return os;
// }
