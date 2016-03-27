#include <cmath>

struct vdc
{
  vdc(int base = 2) : base(base){};

  double operator()(int n){
    int i = 1, r = n, m;
    double res = 0.;
    while(r>0){
      m = r % base;
      res += pow(base,-i)*((double)m);
      r = (r-m)/base;
      i++;
    }
    return res;
  }

private:
  int base;
};
