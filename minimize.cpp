#include <iostream>
#include <cmath>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/lu.hpp>
#include <boost/numeric/ublas/triangular.hpp>
#include <boost/multiprecision/cpp_int.hpp>
#include <boost/multiprecision/cpp_dec_float.hpp>
#include <boost/math/constants/constants.hpp>
#include "minimize.hpp"
#include "tools/utility.hpp"
#include "tools/qr.hpp"

namespace ublas = boost::numeric::ublas;
namespace bmp = boost::multiprecision;

// /************************************************************
// * position
// *************************************************************/

// position::position() {return;}

/************************************************************
* traveller
*************************************************************/

traveller::traveller(int K) : K(K) {return;}

double traveller::cost(double **, double *) {return 0;}

double traveller::condition(double **X, double *M) {
  ublas::matrix<double> J = gradJ(X, M);
  ublas::matrix<double> H = 2*ublas::prod(ublas::trans(J), J);
  return cnum(H);
}

ublas::vector<double> traveller::newton(double **X, double *M, double &COND) {
  ublas::matrix<double> J = gradJ(X, M);
  ublas::vector<double> D = difD(X, M);

  ublas::matrix<double> H = 2*ublas::prod(ublas::trans(J), J);
  ublas::vector<double> Y = 2*ublas::prod(ublas::trans(J), D);


  ublas::lu_factorize(H);
  ublas::inplace_solve(H, Y, ublas::unit_lower_tag());
  ublas::inplace_solve(H, Y, ublas::upper_tag());

  ublas::matrix<double> cH = 2*ublas::prod(ublas::trans(J), J);
  COND = cnum(cH);   

  return Y;
}

ublas::vector<double> traveller::levenberg(double **X, double *M, double &COND, int rp) {
  ublas::matrix<double> J = gradJ(X, M);
  ublas::vector<double> D = difD(X, M);

  ublas::matrix<double> H = 2*ublas::prod(ublas::trans(J), J);
  ublas::vector<double> Y = 2*ublas::prod(ublas::trans(J), D);

  ublas::matrix<double> DH = diagonal(H);

  double **cX = new double*[K];
  cX[0] = new double[K*2];
  for (int k = 1; k < K; ++k) cX[k] = cX[k-1]+2;

  double *cM = new double[K];

  for (int k = 0; k < K; ++k) {
    cX[k][0] = X[k][0];
    cX[k][1] = X[k][1];
    cM[k] = M[k];
  }

  double ii = 0;
  ublas::matrix<double> cH;
  ublas::vector<double> cY;
  double r;
  ublas::vector<double> dXM(3*K);
  do {
    if (ii == 100) {
      std::cerr << "reguralized param is too large!" << std::endl;
      exit(1);
    }
    r = pow(double(10), -rp+ii);

    cH = H+r*DH;
    cY = Y;

    ublas::lu_factorize(cH);
    ublas::inplace_solve(cH, cY, ublas::unit_lower_tag());
    ublas::inplace_solve(cH, cY, ublas::upper_tag());

    update(K, cX, cM, cY);

    ii++;
  }while(cost(X, M) < cost(cX, cM));


  cH = H+r*DH;
  COND = cnum(cH);

  return cY;

}

ublas::matrix<double> traveller::gradJ(double **, double *) {
  ublas::matrix<double> J;
  return J;
}

ublas::vector<double> traveller::difD(double **, double *) {
  ublas::vector<double> D;
  return D;
}

ublas::matrix<double> traveller::diagonal(ublas::matrix<double> &A) {
  ublas::zero_matrix Z(3*K, 3*K);
  ublas::matrix<double> DA = Z;
  for (int i = 0; i < 3*K; ++i) {
    DA(i, i) = A(i, i);
  }
  return DA;
}

/************************************************************
* rookie
*************************************************************/

rookie::rookie(int N, int K, double **A, double *P) : traveller(K), N(N), K(K), A(A), P(P) {return;}

double rookie::cost(double **X, double *M) {
  double v = 0;
  for (int n = 0; n < N; ++n) {
    double p = pointmass_potential(A[n], X, M);
    v += (P[n]-p)*(P[n]-p);
  }
  v /= K;
  return v;
}

double rookie::pointmass_potential(double *a, double **X, double *M) {
  double p = 0;
  for (int k = 0; k < K; ++k) {
    double r = R(a, X[k]);
    p += M[k]/r;
  }
  p *= 0.25;
  p /= M_PI;
  return p;
}


ublas::matrix<double> rookie::gradJ(double **X, double *M) {
  ublas::matrix<double> J(N, 3*K);
  for (int n = 0; n < N; ++n) {
    for (int k = 0; k < K; ++k) {
      double r = R(A[n], X[k]);
      J(n, 3*k+0) = 0.25*M[k]*(A[n][0]-X[k][0])/(M_PI*r*r*r);
      J(n, 3*k+1) = 0.25*M[k]*(A[n][1]-X[k][1])/(M_PI*r*r*r);
      // J(n, 4*k+2) = 0.25*(A[n][2]-X[k][2])/(M_PI*r*r*r);
      J(n, 3*k+2) = 0.25/(M_PI*r);
    }
  }
  return J;
}

ublas::vector<double> rookie::difD(double **X, double *M) {
  ublas::vector<double> D(N);
  for (int n = 0; n < N; ++n) {
    double p = pointmass_potential(A[n], X, M);
    D(n) = P[n]-p;
  }
  return D;
}

/************************************************************
* veteran
*************************************************************/

veteran::veteran(int N, int K, double **A, double **G) : traveller(K), N(N), K(K), A(A), G(G) {
  return;
}

double veteran::cost(double **X, double *M) {
  double v = 0;
  for (int n = 0; n < N; ++n) {
    double *g = pointmass_gravity(A[n], X, M);
    for (int i = 0; i < 2; ++i) {
      v += (G[n][i]-g[i])*(G[n][i]-g[i]);
    }
  }
  v /= K;
  return v;
}

double* veteran::pointmass_gravity(double *a, double **X, double *M) {
  double *g = new double[2];
  g[0] = 0; g[1] = 0;

  for (int k = 0; k < K; ++k) {
    double r = R(a, X[k]);
    for (int i = 0; i < 2; ++i) {
      g[i] += 0.25*M[k]*(a[i]-X[k][i])/(M_PI*r*r*r);
    }
  }
  return g;
}

ublas::matrix<double> veteran::gradJ(double **X, double *M) {
  ublas::matrix<double> J(2*N, 3*K);
  for (int n = 0; n < N; ++n) {
    for (int k = 0; k < K; ++k) {
      double r = R(A[n], X[k]);
      J(2*n+0, 3*k+0) = 0.25*M[k]*(-r*r+3*(A[n][0]-X[k][0])*(A[n][0]-X[k][0]))/(M_PI*r*r*r*r*r);
      J(2*n+0, 3*k+1) = 0.25*M[k]*3*(A[n][0]-X[k][0])*(A[n][1]-X[k][1])/(M_PI*r*r*r*r*r);
      J(2*n+0, 3*k+2) = 0.25*(A[n][0]-X[k][0])/(M_PI*r*r*r);

      J(2*n+1, 3*k+0) = 0.25*M[k]*3*(A[n][0]-X[k][0])*(A[n][1]-X[k][1])/(M_PI*r*r*r*r*r);
      J(2*n+1, 3*k+1) = 0.25*M[k]*(-r*r+3*(A[n][1]-X[k][1])*(A[n][1]-X[k][1]))/(M_PI*r*r*r*r*r);
      J(2*n+1, 3*k+2) = 0.25*(A[n][1]-X[k][1])/(M_PI*r*r*r);
    }
  }
  return J;
}

ublas::vector<double> veteran::difD(double **X, double *M) {
  ublas::vector<double> D(2*N);
  for (int n = 0; n < N; ++n) {
    double *g = pointmass_gravity(A[n], X, M);
    D(2*n+0) = G[n][0]-g[0];
    D(2*n+1) = G[n][1]-g[1];
  }
  return D;
}
