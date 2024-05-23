#ifndef KANSOKU_HPP
#define KANSOKU_HPP

#include <iostream>
#include <cmath>
#include <boost/multiprecision/cpp_int.hpp>
#include <boost/multiprecision/cpp_dec_float.hpp>
#include <boost/math/constants/constants.hpp>

namespace bmp = boost::multiprecision;

namespace KANSOKU {

  template<typename T>
  T pi = static_cast<T>(3.141592653589793238462643383279502884197169399375105820974944592307816406286208998628034825342117068);

  template<typename T>
  T distance(T *x, T *y) {
    T r = sqrt((x[0]-y[0])*(x[0]-y[0])+(x[1]-y[1])*(x[1]-y[1]));
    return r;
  }

  template<typename T>
  T E(T *x, T *y) {
    T r = distance(x, y);
    return 0.25/(pi<T>*r);
  }

  template<typename T>
  T* dE(T *x, T *y) {
    T r = distance(x, y);
    T *v = new T[2];
    v[0] = 0.25*(x[0]-y[0])/(pi<T>*r*r*r);
    v[1] = 0.25*(x[1]-y[1])/(pi<T>*r*r*r);
    return v;
  }
  template<typename T>
  int elliptic(T *y) {
    int v = 0;
    if (0.5*y[0]*y[0]+y[1]*y[1] < 1) v = 10;
    return v;
  }

  template<typename T>
  int ball(T *y) {
    int v = 0;
    if (sqrt(y[0]*y[0]+y[1]*y[1]) < sqrt(1/pi<T>)) v = 1;
    return v;
  }

  template<typename T>
  int square(T *y) {
    int v = 0;
    if (abs(y[0]) < 1 && abs(y[1]) < 1) v = 1;
    return v;
  }

  template<typename T>
  int triangle(T *y) {
    int v = 0;
    if (-0.5 < y[1] && sqrt((T)3)*y[0]+y[1] < 1 && -sqrt((T)3)*y[0]+y[1] < 1) v = 1;
    return v;
  }

  template<typename T>
  T rho(T *x, T *y, int (*domain)(T *)) {
    return E(x, y)*domain(y);
  }

  template<typename T>
  T *drho(T *x, T *y, int(*domain)(T *)) {
    T *dr = new T[2];
    for (int i = 0; i < 2; ++i) {
      dr[i] = dE(x, y)[i]*domain(y);
    }
    return dr;
  }

  // それぞれ地点xでのポテンシャル，重力を返す

  template<typename T>
  T potential(T *x, int (*domain)(T *)) {

    T a = -2; T b = 2;
    T c = -2; T d = 2;

    T dx = (b-a)/100.;
    T dy = (d-c)/100.;

    T area = 0.;

    for (int i = 0; i < 100; ++i) {
      for (int j = 0; j < 100; ++j) {
        T lx, rx, ly, ry;
        lx = a+i*dx; rx = a+(i+1)*dx;
        ly = c+j*dy; ry = c+(j+1)*dy;

        T smallarea = 0.;
        T *y = new T[2];
        y[0] = 0.5*(lx+rx); y[1] = 0.5*(ly+ry);

        smallarea += rho(x, y, domain);
        smallarea *= dx*dy;

        area += smallarea;

      }
    }

    return area;
  }


  template<typename T>
  T *gravity(T *x, int (*domain)(T *)) {

    T a = -2; T b = 2;
    T c = -2; T d = 2;

    T dx = (b-a)/100.;
    T dy = (d-c)/100.;

    T *g = new T[2];
    g[0] = 0; g[1] = 0;

    for (int i = 0; i < 100; ++i) {
      for (int j = 0; j < 100; ++j) {
        T lx, rx, ly, ry;
        lx = a+i*dx; rx = a+(i+1)*dx;
        ly = c+j*dy; ry = c+(j+1)*dy;

        T gx = 0.; T gy = 0.;
        T *y = new T[2];

        y[0] = 0.5*(lx+rx);
        y[1] = 0.5*(ly+ry);

        gx += drho(x, y, domain)[0];
        gy += drho(x, y, domain)[1];

        gx *= dx*dy;
        gy *= dx*dy;

        g[0] += gx;
        g[1] += gy;

      }
    }

    return g;
  }
}


#endif
