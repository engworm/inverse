#ifndef MINIMIZE_CPP
#define MINIMIZE_CPP

#include <iostream>
#include <cmath>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/lu.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/multiprecision/cpp_int.hpp>
#include <boost/multiprecision/cpp_dec_float.hpp>
#include <boost/math/constants/constants.hpp>
#include "tools/utility.hpp"

namespace ublas = boost::numeric::ublas;
namespace bmp = boost::multiprecision;

namespace MINIMIZE {

  template<typename T>
  T pi = static_cast<T>(3.141592653589793238462643383279502884197169399375105820974944592307816406286208998628034825342117068);

  //------------------------------------------------------------
  // traveller
  //------------------------------------------------------------

  template <typename T>
  class traveller {
    public:
      traveller(int K) : K(K) {return;}

    public:
      virtual T cost(T **, T *) {return 0;}

      ublas::vector<T> newton(T **X, T *M) {
        ublas::matrix<T> J = gradJ(X, M);
        ublas::vector<T> D = difD(X, M);

        ublas::matrix<T> H = 2*ublas::prod(ublas::trans(J), J);
        ublas::vector<T> Y = 2*ublas::prod(ublas::trans(J), D);


        ublas::lu_factorize(H);
        ublas::inplace_solve(H, Y, ublas::unit_lower_tag());
        ublas::inplace_solve(H, Y, ublas::upper_tag());

        ublas::matrix<T> cH = 2*ublas::prod(ublas::trans(J), J);

        return Y;
      }

      ublas::vector<T> levenberg(T **X, T *M, int rp) {
        ublas::matrix<T> J = gradJ(X, M);
        ublas::vector<T> D = difD(X, M);

        ublas::matrix<T> H = 2*ublas::prod(ublas::trans(J), J);
        ublas::vector<T> Y = 2*ublas::prod(ublas::trans(J), D);

        ublas::matrix<T> DH = diagonal(H);

        T **cX = new T*[K];
        cX[0] = new T[K*2];
        for (int k = 1; k < K; ++k) cX[k] = cX[k-1]+2;

        T *cM = new T[K];

        for (int k = 0; k < K; ++k) {
          cX[k][0] = X[k][0];
          cX[k][1] = X[k][1];
          cM[k] = M[k];
        }

        T ii = 0;
        ublas::matrix<T> cH;
        ublas::vector<T> cY;
        T r;
        ublas::vector<T> dXM(3*K);
        do {
          if (ii == 100) {
            std::cerr << "reguralized param is too large!" << std::endl;
            exit(1);
          }
          r = pow(T(10), -rp+ii);

          cH = H+r*DH;
          cY = Y;

          ublas::lu_factorize(cH);
          ublas::inplace_solve(cH, cY, ublas::unit_lower_tag());
          ublas::inplace_solve(cH, cY, ublas::upper_tag());

          update(K, cX, cM, cY);

          ii++;
        }while(cost(X, M) < cost(cX, cM));


        cH = H+r*DH;

        return cY;

      }
  
    private:

      virtual ublas::matrix<T> gradJ(T **X, T *M) {
        ublas::matrix<T> J;
        return J;
      } 

      virtual ublas::vector<T> difD(T **X, T *M) {
        ublas::vector<T> D;
        return D;
      }

      ublas::matrix<T> diagonal(ublas::matrix<T> &A) {
        ublas::zero_matrix Z(3*K, 3*K);
        ublas::matrix<T> DA = Z;
        for (int i = 0; i < 3*K; ++i) {
          DA(i, i) = A(i, i);
        }
        return DA;
      }

    private:
      int K;
  };

  //------------------------------------------------------------
  // rookie
  //------------------------------------------------------------

  template <typename T>
  class rookie : public traveller<T> {
    public:
      rookie(int N, int K, T **A, T *P) : traveller<T>(K), N(N), K(K), A(A), P(P) {return;}
      // ~rookie();

    public:
      virtual T cost(T **X, T *M) override {
        T v = 0;
        for (int n = 0; n < N; ++n) {
          T p = pointmass_potential(A[n], X, M);
          v += (P[n]-p)*(P[n]-p);
        }
        v /= N;
        return v;
      }

  
    private:
      T pointmass_potential(T *a, T **X, T *M) {
        T p = 0;
        for (int k = 0; k < K; ++k) {
          T r = R(a, X[k]);
          // p += M[k]*log(r);
          p += M[k]/r;
        }
        // p *= -0.25/pi<T>; //
        p *= 0.25/pi<T>; //
        return p;
      } 

      virtual ublas::matrix<T> gradJ(T **X, T *M) override {
        ublas::matrix<T> J(N, 3*K);
        for (int n = 0; n < N; ++n) {
          for (int k = 0; k < K; ++k) {
            T r = R(A[n], X[k]);
            J(n, 3*k+0) = 0.25*M[k]*(A[n][0]-X[k][0])/(pi<T>*r*r*r);
            J(n, 3*k+1) = 0.25*M[k]*(A[n][1]-X[k][1])/(pi<T>*r*r*r);
            J(n, 3*k+2) = 0.25/(pi<T>*r);

            // J(n, 3*k+0) = -0.5*M[k]*(A[n][0]-X[k][0])/(pi<T>*r*r);
            // J(n, 3*k+1) = -0.5*M[k]*(A[n][1]-X[k][1])/(pi<T>*r*r);
            // J(n, 3*k+2) = -0.5*log(r)/pi<T>;
          }
        }
        return J;
      } 

      virtual ublas::vector<T> difD(T **X, T *M) override {
        ublas::vector<T> D(N);
        for (int n = 0; n < N; ++n) {
          T p = pointmass_potential(A[n], X, M);
          D(n) = P[n]-p;
        }
        return D;
      }

    private:
      int N;
      int K;
      T **A;
      T *P;
  };

  //------------------------------------------------------------
  // veteran
  //------------------------------------------------------------

  template <typename T>
  class veteran : public traveller<T> {

    public:
      veteran(int N, int K, T **A, T **G) : traveller<T>(K), N(N), K(K), A(A), G(G) {return;}

    public:
      virtual T cost(T **X, T *M) override {
        T v = 0;
        for (int n = 0; n < N; ++n) {
          T *g = pointmass_gravity(A[n], X, M);
          for (int i = 0; i < 2; ++i) {
            v += (G[n][i]-g[i])*(G[n][i]-g[i]);
          }
        }
        v /= N;
        return v;
      }

    private:

      T* pointmass_gravity(T *a, T **X, T *M) {
        T *g = new T[2];
        g[0] = 0; g[1] = 0;

        for (int k = 0; k < K; ++k) {
          T r = R(a, X[k]);
          for (int i = 0; i < 2; ++i) {
            // g[i] += -0.5*M[k]*(a[i]-X[k][i])/(pi<T>*r*r);
            g[i] += 0.25*M[k]*(a[i]-X[k][i])/(pi<T>*r*r*r);
          }
        }
        return g;
      }

      virtual ublas::matrix<T> gradJ(T **X, T *M) override {
        ublas::matrix<T> J(2*N, 3*K);
        for (int n = 0; n < N; ++n) {
          for (int k = 0; k < K; ++k) {
            T r = R(A[n], X[k]);
            J(2*n+0, 3*k+0) = 0.25*M[k]*(-r*r+3*(A[n][0]-X[k][0])*(A[n][0]-X[k][0]))/(pi<T>*r*r*r*r*r);
            J(2*n+0, 3*k+1) = 0.25*M[k]*3*(A[n][0]-X[k][0])*(A[n][1]-X[k][1])/(pi<T>*r*r*r*r*r);
            J(2*n+0, 3*k+2) = 0.25*(A[n][0]-X[k][0])/(pi<T>*r*r*r);

            J(2*n+1, 3*k+0) = 0.25*M[k]*3*(A[n][0]-X[k][0])*(A[n][1]-X[k][1])/(pi<T>*r*r*r*r*r);
            J(2*n+1, 3*k+1) = 0.25*M[k]*(-r*r+3*(A[n][1]-X[k][1])*(A[n][1]-X[k][1]))/(pi<T>*r*r*r*r*r);
            J(2*n+1, 3*k+2) = 0.25*(A[n][1]-X[k][1])/(pi<T>*r*r*r);

            // J(2*n+0, 3*k+0) = -0.5*M[k]*((A[n][0]-X[k][0])*(A[n][0]-X[k][0])-(A[n][1]-X[k][1])*(A[n][1]-X[k][1]))/(pi<T>*r*r*r*r);
            // J(2*n+0, 3*k+1) = -0.5*M[k]*(A[n][0]-X[k][0])*(A[n][1]-X[k][1])/(pi<T>*r*r*r*r);
            // J(2*n+0, 3*k+2) = -0.5*(A[n][0]-X[k][0])/(pi<T>*r*r);

            // J(2*n+1, 3*k+0) = -0.5*M[k]*(A[n][0]-X[k][0])*(A[n][1]-X[k][1])/(pi<T>*r*r*r*r);
            // J(2*n+1, 3*k+1) = -0.5*M[k]*((A[n][1]-X[k][1])*(A[n][1]-X[k][1])-(A[n][0]-X[k][0])*(A[n][0]-X[k][0]))/(pi<T>*r*r*r*r);
            // J(2*n+1, 3*k+2) = -0.5*(A[n][1]-X[k][1])/(pi<T>*r*r);
          }
        }
        return J;
      } 

      virtual ublas::vector<T> difD(T **X, T *M) override {
        ublas::vector<T> D(2*N);
        for (int n = 0; n < N; ++n) {
          T *g = pointmass_gravity(A[n], X, M);
          D(2*n+0) = G[n][0]-g[0];
          D(2*n+1) = G[n][1]-g[1];
        }
        return D;
      }

    private:
      int N;
      int K;
      T **A;
      T **G;
  };

}

#endif
