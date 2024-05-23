#include <iostream>
#include <boost/multiprecision/cpp_bin_float.hpp>

namespace bmp = boost::multiprecision;

int main() {
  bmp::cpp_dec_float_100 x = 1;
  std::cout << x << std::endl;
}