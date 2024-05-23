#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <random>
// #include <exception>
#include <stdexcept>
#include <boost/multiprecision/cpp_int.hpp>
#include <boost/multiprecision/cpp_dec_float.hpp>
#include <boost/math/constants/constants.hpp>
#include <boost/random.hpp>
#include <boost/random/random_device.hpp>
#include "param.hpp"

namespace bmp = boost::multiprecision;
using namespace MY_PARAM;


int main() {
  //------------------------------------------------------------
  // distribution
  //------------------------------------------------------------
  boost::random::uniform_real_distribution<bmp::cpp_dec_float_100> radi(0, BODY);
  boost::random::uniform_real_distribution<bmp::cpp_dec_float_100> theta(0., 0.5*M_PI);
  boost::random::uniform_real_distribution<bmp::cpp_dec_float_100> phi(0., 2*M_PI);
  boost::random::uniform_real_distribution<bmp::cpp_dec_float_100> mass(0, BODY*BODY*M_PI/K);   // 第二引数に注意


  //------------------------------------------------------------
  // initial value (random, out)
  //------------------------------------------------------------
  bmp::cpp_dec_float_100 **X = new bmp::cpp_dec_float_100*[K];
  X[0] = new bmp::cpp_dec_float_100[K*2];
  for (int k = 1; k < K; ++k) X[k] = X[k-1]+2;

  bmp::cpp_dec_float_100 *M = new bmp::cpp_dec_float_100[K];

  boost::random::random_device iseed;
  boost::random::independent_bits_engine<boost::random::mt19937, 2048, bmp::cpp_int> iengine(iseed);

  std::cout << std::setprecision(std::numeric_limits<bmp::cpp_dec_float_100>::digits10 + 1);
  for (int k = 0; k < K; ++k) {
    bmp::cpp_dec_float_100 r = radi(iengine);
    bmp::cpp_dec_float_100 p = phi(iengine);

    std::cout << r*cos(p) << " ";
    std::cout << r*sin(p) << " ";
    std::cout << mass(iengine) << std::endl;
    // std::cout << TOTALMASS/K << std::endl;

  }
  
  return 0;

}
