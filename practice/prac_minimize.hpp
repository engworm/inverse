#include <iostream>

#ifndef MINIMIZE_CPP
#define MINIMIZE_CPP

class traveller {
  public:
    traveller();

  public:
    virtual double cost(double **, double *);
};

class rookie : public traveller {
  public:
    rookie(int N, int K, double **A, double *P);
    ~rookie();

  public:
    double cost(double **X, double *M) override;
  
  private:
    int N;
    int K;
    double **A;
    double *P;
};

class veteran : public traveller {

};


#endif
