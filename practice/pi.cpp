#include <iostream>
#include <iomanip>
#include <boost/multiprecision/cpp_dec_float.hpp>
#include <boost/math/constants/constants.hpp>

namespace bmp = boost::multiprecision;

int main() {

  bmp::cpp_dec_float_100 pi = boost::math::constants::pi<bmp::cpp_dec_float_100>();
  std::cout << std::setprecision(100);
  std::cout << pi << std::endl;

  return 0;

}