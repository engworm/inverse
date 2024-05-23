#include <iostream>

auto lambda = [](){
  double o[2]{0};
  std::cout << o[0] << " " << o[1] << std::endl;
};

int main() {
  lambda();
  return 0;
}
