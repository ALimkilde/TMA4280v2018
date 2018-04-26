#ifndef TMA_QUADRATURE_H_
#define TMA_QUADRATURE_H_

#include <tma/types.h>
#include <math.h>

namespace tma
{

struct QR
{

  struct interval
  {
    uint num_from_deg(uint d){
      return ceil((d+1)/2.0);
    }

    void get_quadrature(uint d, double* x, double* w){
      uint n = num_from_deg(d);
      if (n == 1){
        x[0] = 0.5;
        w[0] = 1;
      } else if (n == 2){
        x[0] = 0.5 - 0.5/sqrt(3.0);
        x[1] = 0.5 + 0.5/sqrt(3.0);
        w[0] = 0.5;
        w[1] = 0.5;
      } else if (n == 5){
        x[0] = 0.5;
        x[1] = 0.5-0.5/3.0*sqrt(5.0-2.0*sqrt(10.0/7.0));
        x[2] = 0.5+0.5/3.0*sqrt(5-2*sqrt(10.0/7.0));
        x[3] = 0.5-0.5/3.0*sqrt(5+2*sqrt(10.0/7.0));
        x[4] = 0.5+0.5/3.0*sqrt(5+2*sqrt(10.0/7.0));

        w[0] = 0.5*128/225;
        w[1] = 0.5*(322+13*sqrt(70))/900;
        w[2] = 0.5*(322+13*sqrt(70))/900;
        w[3] = 0.5*(322-13*sqrt(70))/900;
        w[4] = 0.5*(322-13*sqrt(70))/900;
      }
    }
  };

  interval interval_;

  struct triangle
  {

  };
};
} /* namespace tma */

#endif /* TMA_QUADRATURE_H_ */
