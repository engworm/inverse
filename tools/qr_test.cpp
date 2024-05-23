#include <iostream>
#include <fstream>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/lu.hpp>
#include <boost/numeric/ublas/triangular.hpp>
#include "qr.hpp"

namespace ublas = boost::numeric::ublas;

int main() {
  int K = 2;
  ublas::matrix<double> H(K, K);
  std::ifstream HESSIAN("hessian.txt");

  for (int i = 0; i < K; ++i) {
    for (int j = 0; j < K; ++j) {
      HESSIAN >> H(i, j);
    }
  }

  std::cout << H << std::endl;

  ublas::identity_matrix<double> I(K, K);
  ublas::matrix<double> Q = I;

  std::cout << cnum(H) << std::endl;
  std::cout << H << std::endl;


  // std::cout << cnum(H) << std::endl;
}
