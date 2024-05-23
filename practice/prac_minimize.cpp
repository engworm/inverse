#include <iostream>
#include "prac_minimize.hpp"

traveller::traveller() {return;}

double traveller::cost(double **, double *) {return 0;}

rookie::rookie(int N, int K, double **A, double *P) :N(N), K(K), A(A), P(P) {
  return;
}

double rookie::cost(double **X, double *M) {
  return 1.;
}
