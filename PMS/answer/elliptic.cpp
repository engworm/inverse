#include <iostream>
#include <cmath>

int main() {
  double h = 1.E-2;
  double x = -2;
  double y = -2;
  for (int i = 0; i < 4/h; ++i) {
    for (int j = 0; j < 4/h; ++j) {
      double xi = x+j*h;
      double yj = y+i*h;
      int v = 0;
      if (0.5*xi*xi+yj*yj < 1) v = 10;
      std::cout << v << " ";
    }
    std::cout << std::endl;
  }
  return 0;
}
