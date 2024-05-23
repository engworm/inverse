#ifndef UTILITY_HPP
#define UTILITY_HPP

#include <iostream>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/multiprecision/cpp_int.hpp>
#include <boost/multiprecision/cpp_dec_float.hpp>
#include <boost/math/constants/constants.hpp>


namespace ublas = boost::numeric::ublas;

template <typename T>
T R(T *a, T *x) {
  T r = 0;
  for (int i = 0; i < 2; ++i) {
    r += (a[i]-x[i])*(a[i]-x[i]);
  }
  return sqrt(r);
}

template <typename T>
T ball_potential(T *a, T cx, T cy) {
  T c[2];
  c[0] = cx; c[1] = cy;

  T r = R(c, a);
  T v = 0.25/(M_PI*r);

  return v;
}

template <typename T>
T* ball_gravity(T *a, T cx, T cy) {
  T c[2];
  c[0] = cx; c[1] = cy;

  T r = R(c, a);
  T *g = new T[2];
  g[0] = 0.25*(a[0]-c[0])/(M_PI*r*r*r);
  g[1] = 0.25*(a[1]-c[1])/(M_PI*r*r*r);

  return g;

}


template <typename T>
bool check(int K, T **old_X, T *old_M, T **new_X, T *new_M) {
  bool b = false;
  T O[2]{0};
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

template <typename T>
void update(int K, T **X, T *M, ublas::vector<T> dXM) {
  for (int k = 0; k < K; ++k) {
    X[k][0] += dXM(3*k+0);
    X[k][1] += dXM(3*k+1);
    M[k] += dXM(3*k+2);
  }
}

template <typename T>
std::ostringstream XM(int K, T **X, T *M) {
  std::ostringstream oss;
  for (int k = 0; k < K; ++k) {
    oss << X[k][0] << " ";
    oss << X[k][1] << " ";
    oss << M[k] << '\n';
  }
  return oss;
}


template <typename T>
T total_mass(int K, T *M) {
  T v = 0;
  for (int k = 0; k < K; ++k) {
    v += M[k];
  }
  return v;
}

template <typename T>
T partial_mass(int K, T **X, T *M) {
  T v = 0;
  T O[2]{0};
  for (int k = 0; k < K; ++k) {
    if (R(O, X[k]) < 2) v += M[k];
  }
  return v;
}

#endif
