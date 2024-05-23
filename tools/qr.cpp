#include <iostream>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/lu.hpp>
#include <boost/numeric/ublas/triangular.hpp>
#include <boost/multiprecision/cpp_int.hpp>
#include <boost/multiprecision/cpp_dec_float.hpp>
#include "qr.hpp"

namespace ublas = boost::numeric::ublas;
namespace bmp = boost::multiprecision;

void qr(ublas::matrix<double> &A, ublas::matrix<double> &Q) {
  int N = A.size1();
  for (int i = 0; i < N-1; ++i) {
    small_qr(A, Q, i);
  }
}

void small_qr(ublas::matrix<double> &A, ublas::matrix<double> &Q, int n) {
  ublas::matrix<double> M = cut(A, n);

  int N = M.size1();
  ublas::vector<double> a = ublas::column(M, 0);

  ublas::vector<double> u = a;
  double r = sign(u(0))*ublas::norm_2(a);
  u(0) += r;

  householder(Q, u, n);

  for (int i = 1; i < N; ++i) {
    ublas::vector<double> b = ublas::column(M, i);
    ublas::column(M, i) -= (ublas::inner_prod(u, b)/ublas::inner_prod(u, a))*u;
  }

  M(0, 0) = -r;
  for (int i = 1; i < N; ++i) M(i, 0) = 0;

  // // AにMを埋め込む
  embed(A, M, n);

}

void householder(ublas::matrix<double> &Q, ublas::vector<double> v, int n) {
  int N = v.size();

  v /= ublas::norm_2(v);

  ublas::zero_matrix<double> Z(N, N);
  ublas::matrix<double> H = Z;
  for (int i = 0; i < N; ++i) {
    for (int j = 0; j < N; ++j) {
      H(i, j) = -2.*v(i)*v(j);
    }
  }

  // std::cout << "H" << std::endl;
  // std::cout << H << std::endl;

  int L = Q.size1();
  ublas::identity_matrix<double> I(L, L);
  ublas::zero_matrix<double> ZZ(L, L);
  ublas::matrix<double> cH = ZZ;

  embed(cH, H, n);
  cH += I;

  // std::cout << "cH" << std::endl;
  // std::cout << cH << std::endl;

  // std::cout << "Q" << std::endl;
  // std::cout << Q << std::endl;
  // std::cout << ublas::prod(Q, ublas::trans(Q)) << std::endl;

  Q = ublas::prod(cH, Q);

  // std::cout << "Q" << std::endl;
  // std::cout << Q << std::endl;

}

double cnum(ublas::matrix<double> &A) {
  int N = A.size1();
  ublas::identity_matrix<double> I(N, N);
  ublas::matrix<double> Q = I;
  for (int i = 0; i < 20; ++i) {
    Q = I;
    qr(A, Q);
    A = ublas::prod(A, Q);
  }
  return abs(A(0, 0)/A(N-1, N-1));
}

// ublas::matrix<double> convert_mptodp(ublas::matrix<bmp::cpp_dec_float_100> &A) {
  // int N = A.size1();
  // ublas::matrix<double> M(N, N);
  // for (int i = 0; i < N; ++i) {
    // for (int j = 0; j < N; ++j) {
      // M(i, j) = (double)A(i, j);
    // }
  // }
  // return M;
// }

// ublas::matrix<double> convert_mmptodp(ublas::matrix<mpf_dec_float_500> &A) {
  // int N = A.size1();
  // ublas::matrix<double> M(N, N);
  // for (int i = 0; i < N; ++i) {
    // for (int j = 0; j < N; ++j) {
      // M(i, j) = (double)A(i, j);
    // }
  // }
  // return M;
// }


inline double sign(double x) {
  if (x >= 0) return 1.;
  else return -1.;
}

ublas::matrix<double> cut(ublas::matrix<double> &A, int n) {
  int N = A.size1();
  ublas::matrix<double> M(N-n, N-n);
  for (int i = 0; i < N-n; ++i) {
    for (int j = 0; j < N-n; ++j) {
      M(i, j) = A(n+i, n+j);
    }
  }
  return M;
}

void embed(ublas::matrix<double> &A, ublas::matrix<double> &M, int n) {
  int N = A.size1();
  for (int i = 0; i < N-n; ++i) {
    for (int j = 0; j < N-n; ++j) {
      A(n+i, n+j) = M(i, j);
    }
  }
}
