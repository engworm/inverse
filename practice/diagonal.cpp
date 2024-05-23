#include <iostream>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>

namespace ublas = boost::numeric::ublas;

int main() {
  ublas::matrix<double> A(2,2);
  std::cout << A << std::endl;
  return 0;
}
