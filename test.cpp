#include <iostream>
#include <cmath>
#include <fstream>
#include <tuple>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/multiprecision/cpp_int.hpp>
#include <boost/multiprecision/cpp_dec_float.hpp>
#include <boost/math/constants/constants.hpp>
#include <boost/hana.hpp>
#include <boost/hana/ext/std/tuple.hpp>
#include "minimize.hpp"
#include "tools/utility.hpp"
#include "kansoku.hpp"

namespace bmp = boost::multiprecision;

template <typename T>
T pi = static_cast<T>(3.141592653589793238462643383279502884197169399375105820974944592307816406286208998628034825342117068);

template <typename T>
void minimize_procedure(char kansoku, char precision, char &method, int N, int K, T R, int tol) {
  //------------------------------------------------------------
  // 観測値
  //------------------------------------------------------------
  T **A = new T*[N];
  A[0] = new T[N*2];
  for (int n = 1; n < N; ++n) A[n] = A[n-1]+2;

  for (int n = 0; n < N; ++n) {
    T phi = 2*n*pi<T>/N;
    A[n][0] = R*cos(phi);
    A[n][1] = R*sin(phi);
  }

  T *P = new T[N];

  for (int n = 0; n < N; ++n) {
    // T p = ball_potential(A[n], 0, 0);
    double a[2] = {static_cast<double>(A[n][0]), static_cast<double>(A[n][1])};
    T p = static_cast<T>(KANSOKU::potential(a, KANSOKU::elliptic));
    P[n] = round(p*1.E+4)/1.E+4;
    // std::cerr << "row  : " << p << std::endl;
    // std::cerr << "round: " << P[n] << std::endl;
    // P[n] = p;
  }

  T **G = new T*[N];
  G[0] = new T[N*2];
  for (int n = 1; n < N; ++n) G[n] = G[n-1]+2;

  for (int n = 0; n < N; ++n) {
    // g = ball_gravity(A[n], 0, 0);
    double a[2] = {static_cast<double>(A[n][0]), static_cast<double>(A[n][1])};
    double *v = KANSOKU::gravity(a, KANSOKU::elliptic);
    T g[2] = {static_cast<T>(v[0]), static_cast<T>(v[1])};
    G[n][0] = round(g[0]*1.E+4)/1.E+4;
    G[n][1] = round(g[1]*1.E+4)/1.E+4;
    // std::cerr << "row  :" << g[0] << std::endl;
    // std::cerr << "round: " << G[n][0] << std::endl;
    // G[n][0] = g[0];
    // G[n][1] = g[1];
  }


  //------------------------------------------------------------
  // FILE & initialize traveller
  //------------------------------------------------------------
  MINIMIZE::traveller<T> *minimize;

  std::ofstream LOG;
  std::ofstream DATA;
  if (kansoku == 'p') {
    if (precision == 'd') {
      LOG.open("pote/dp/log/log.txt");
      DATA.open("pote/dp/log/data.txt");
      minimize = new MINIMIZE::rookie<T>(N, K, A, P);
    }
    else if (precision == 'm') {
      LOG.open("pote/mp/log/log.txt");
      DATA.open("pote/mp/log/data.txt");
      minimize = new MINIMIZE::rookie<T>(N, K, A, P);
    }
  }
  else if (kansoku == 'g') {
    if (precision == 'd') {
      LOG.open("grav/dp/log/log.txt");
      DATA.open("grav/dp/log/data.txt");
      minimize = new MINIMIZE::veteran<T>(N, K, A, G);
    }
    if (precision == 'm') {
      LOG.open("grav/mp/log/log.txt");
      DATA.open("grav/mp/log/data.txt");
      minimize = new MINIMIZE::veteran<T>(N, K, A, G);
    }
  }

  //------------------------------------------------------------
  // 初期値
  //------------------------------------------------------------

  T **X = new T*[K];
  X[0] = new T[K*2];
  for (int k = 1; k < K; ++k) X[k] = X[k-1]+2;

  T *M = new T[K];

  for (int k = 0; k < K; ++k) {
    T phi = 2*k*pi<T>/K;
    X[k][0] = 0.8*cos(phi);
    X[k][1] = 0.8*sin(phi);
    M[k] = 40*pi<T>/K;
  }


  //------------------------------------------------------------
  // minimize session
  //------------------------------------------------------------


  T cost = minimize->cost(X, M);
  LOG << "before ##############################" << std::endl;
  LOG << "cost: " << cost << std::endl;
  LOG << "total mass  : " << total_mass(K, M) << std::endl;
  LOG << "partial mass: " << partial_mass(K, X, M) << std::endl;
  LOG << "XM" << std::endl;
  LOG << XM(K, X, M).str();
  LOG << "#####################################" << std::endl;


  cost = minimize->cost(X, M);
  DATA << 0 << " " << cost << " " << NAN << std::endl;
  
  T **cX = new T*[K];
  cX[0] = new T[K*2];
  for (int k = 1; k < K; ++k) cX[k] = cX[k-1]+2;

  T *cM = new T[K];

  int roop = 1;
  T v1, v2;
  bool change;
  do {
    if (roop == 300) {
      std::cerr << "too long" << std::endl;
      break;
    }

    for (int k = 0; k < K; ++k) {
      cX[k][0] = X[k][0];
      cX[k][1] = X[k][1];
      cM[k] = M[k];
    }

    ublas::vector<T> dXM;
    if (method == 'n') {
      dXM = minimize->newton(cX, cM);
    }
    else if (method == 'l') {
      dXM = minimize->levenberg(cX, cM, tol);
    }
    else {
      std::cerr << "choose your method correctly" << std::endl;
      break;
    }

    // cX, cMは変更されているので初期化し直す
    for (int k = 0; k < K; ++k) {
      cX[k][0] = X[k][0];
      cX[k][1] = X[k][1];
      cM[k] = M[k];
    }

    update(K, cX, cM, dXM);

    // LOGFILE
    cost = minimize->cost(cX, cM);
    LOG << "-----------------------------------" << std::endl;
    LOG << "counter: " << roop << std::endl;
    LOG << "cost: " << cost << std::endl;
    LOG << "total mass  : " << total_mass(K, cM) << std::endl;
    LOG << "partial mass: " << partial_mass(K, cX, cM) << std::endl;
    // LOG << "condition num: " << cond << std::endl;
    LOG << "XM" << std::endl;
    LOG << XM(K, cX, cM).str();
    LOG << "-----------------------------------" << std::endl;

    DATA << roop << " " << cost << " " << NAN << std::endl;

    roop++;

    change = check(K, X, M, cX, cM);

    for (int k = 0; k < K; ++k) {
      X[k][0] = cX[k][0];
      X[k][1] = cX[k][1];
      M[k] = cM[k];
    }

  }while(roop <= 20);
  // }while(change || (roop < 100));


  std::cout << XM(K, X, M).str();

  LOG << "after ##############################" << std::endl;
  LOG << "cost: " << minimize->cost(X, M) << std::endl;
  LOG << "total mass  : " << total_mass(K, M) << std::endl;
  LOG << "partial mass: " << partial_mass(K, X, M) << std::endl;
  LOG << "XM" << std::endl;
  LOG << XM(K, X, M).str();
  LOG << "#####################################" << std::endl;

  delete[] A[0];
  delete[] A;
  delete[] P;
  delete[] G[0];
  delete[] G;

  delete[] cX[0];
  delete[] cX;
  delete[] cM;

  delete[] X[0];
  delete[] X;
  delete[] M;

  LOG.close();
  DATA.close();

}


int main(int argc, char *argv[]) {

  char *boundary = argv[1];
  char *precision = argv[2];
  char *method = argv[3];
  int N = strtol(argv[4], nullptr, 10);
  int K = strtol(argv[5], nullptr, 10);
  double R = strtod(argv[6], nullptr);
  int e = strtol(argv[7], nullptr, 10);

  if (boundary[0] == 'p') {
    if (precision[0] == 'd') {
      minimize_procedure('p', 'd', method[0], N, K, R, e);
    }
    else if (precision[0] == 'm') {
      minimize_procedure('p', 'm', method[0], N, K, (bmp::cpp_dec_float_50)R, e);
    }
  }
  else if (boundary[0] == 'g') {
    if (precision[0] == 'd') {
      minimize_procedure('g', 'd', method[0], N, K, R, e);
    }
    else if (precision[0] == 'm') {
      minimize_procedure('g', 'm', method[0], N, K, (bmp::cpp_dec_float_50)R, e);
    }
  }
  else {
    std::cerr << "check your command line arg 'boundary' again" << std::endl;
    exit(1);
  }

  return 0;
}

