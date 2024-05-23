#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <iomanip>

int main(int argc, char *argv[]) {
  //------------------------------------------------------------
  // param
  //------------------------------------------------------------
  double a = -3; double b = 3;
  double c = -3; double d = 3;

  double h = 1.E-2;
  int xN = (b-a)/h;
  int yN = (d-c)/h;

  int K = std::strtol(argv[1], nullptr, 10);  

  //------------------------------------------------------------
  // 格子点に再配分する
  //------------------------------------------------------------
  double **X = new double*[K];
  X[0] = new double[K*2];
  for (int k = 1; k < K; ++k) X[k] = X[k-1]+2;

  double *M = new double[K];

  for (int k = 0; k < K; ++k) {
    for (int j = 0; j < 2; ++j) {
      std::cin >> X[k][j];
    }
    std::cin >> M[k];
  }


  double **S = new double*[xN];
  S[0] = new double[xN*yN];
  for (int i = 1; i < xN; ++i) S[i] = S[i-1]+yN;


  for (int i = 0; i < xN; ++i) {
    for (int j = 0; j < yN; ++j) S[i][j] = 0.;
  }

  for (int k = 0; k < K; ++k) {
    int xn = X[k][0]/h+xN/2;
    int yn = X[k][1]/h+yN/2;
    if (0 <= xn && xn < xN && 0<= yn && yn < yN) {
      S[xn][yn] += M[k];
    }
  }

  for (int i = 0; i < xN; ++i) {
    for (int j = 0; j < yN; ++j) {
      std::cout << S[j][i]/(h*h) << " ";
    }
    std::cout << std::endl;
  }

  delete[] X[0];
  delete[] X;
  delete[] M;

  delete[] S[0];
  delete[] S;


}
