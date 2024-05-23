#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>

int main(int argc, char *argv[]) {
  //------------------------------------------------------------
  // param
  //------------------------------------------------------------

  double a = -3; double b = 3;
  double c = -3; double d = 3;

  int counter = std::strtol(argv[1], nullptr, 10);

  double ep = 1.E-5;

  double h =  1.E-2;
  int xN = (b-a)/h;
  int yN = (d-c)/h;

  //------------------------------------------------------------
  // scatter
  //------------------------------------------------------------

  double **S = new double*[xN];
  S[0] = new double[xN*yN];
  for (int i = 1; i < xN; ++i) S[i] = S[i-1]+yN;

  for (int i = 0; i < xN; ++i) {
    for (int j = 0; j < yN; ++j) {
      std::cin >> S[i][j];
      S[i][j] *= h*h;
    }
  }

  // 変化量を記録する matrix
  double **dS = new double*[xN];
  dS[0] = new double[xN*yN];
  for (int i = 1; i < xN; ++i) dS[i] = dS[i-1]+yN;

  for (int i = 0; i < xN; ++i) {
    for (int j = 0; j < yN; ++j) dS[i][j] = 0.;
  }

  // change?
  bool change = false;

  // scatter
  for (int i = 0; i < xN; ++i) {
    for (int j = 0; j < yN; ++j) {

      if (S[i][j] == 0) continue;

      double d = S[i][j]-10*h*h;
      if (d > ep) {
        change = true;

        dS[i][j] += -(d+ep);

        double share = (d+ep)*0.25;
        dS[i-1][j] += share; dS[i+1][j] += share;
        dS[i][j-1] += share; dS[i][j+1] += share;

      }
    }
  }

  if (!change) {
    std::cerr << counter << " step doesn't change matrix state" << std::endl;
    return 1;
  }

  // update
  for (int i = 0; i < xN; ++i) {
    for (int j = 0; j < yN; ++j) S[i][j] += dS[i][j];
  }

  for (int i = 0; i < xN; ++i) {
    for (int j = 0; j < yN; ++j) {
      std::cout << S[i][j]/(h*h) << " ";
    }
    std::cout << std::endl;
  }

  
  delete[] S[0];
  delete[] S;

  delete[] dS[0];
  delete[] dS;

  return 0;
}
