#include <iostream>
#include <cmath>

int main() {
  // if (-0.5 < y[1] && sqrt((T)3)*y[0]+y[1] < 1 && -sqrt((T)3)*y[0]+y[1] < 1) v = 1;
  double h = 1.E-2;
  double x = -2;
  double y = -2;
  for (int i = 0; i < 4/h; ++i) {
    for (int j = 0; j < 4/h; ++j) {
      double xi = x+j*h;
      double yj = y+i*h;
      int v = 0;
      if (-0.5 < yj && sqrt(3)*xi+yj < 1 && -sqrt(3)*xi+yj < 1) v = 1;
      std::cout << v << " ";
    }
    std::cout << std::endl;
  }
  return 0;
}
