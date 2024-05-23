#include <iostream>

template <typename T>
class base {
  public:
    base(T x); 
    void get();
  private:
    T x;
};

template <>
class base<double> {
  public:
    base(double x) : x(x) {std::cout << "mogumogu" << std::endl;}

    void get() {
      std::cout << 2*x << std::endl;
    }

  private:
    double x;
};

int main() {

  // base<int> test(3);
  // test.get();

  base<double> test2(3.2);
  test2.get();

}