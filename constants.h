#ifndef CONSTANTS_H
#define CONSTANTS_H


#include <cmath>

struct constants
{
  static const double em;// = 1.8228883e03;//from DVR-3D
  static const double htocm;// = 2.19474316e05;//from DVR-3D
  static const double atob;// = 1.889726124626;
  static const double pi2; //= sqrt(M_PI);
  static const double pi4; //= sqrt(sqrt(M_PI));
  enum matr_type {PURE,TRANSFORMED,FOR_AAT};
};


#endif
