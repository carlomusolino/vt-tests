#ifndef FVCC_1D_REC_HELPERS_HH
#define FVCC_1D_REC_HELPERS_HH

#include <cmath>

namespace advect_1d {

  static inline constexpr int sign(const double& x){return (0 < x) - (x < 0);};

  static inline constexpr double max(const double& x, const double& y){ return (x>y) ? x : y;};
  
  template<typename... Targs>
  inline constexpr double max(const double& x, const double & y, Targs... Fargs)
  {
    return max(x,max(y,Fargs...));
  }
  
  
  inline constexpr double min(const double& x, const double& y){ return (x>y) ? y : x;};
  
  template<typename... Targs>
  inline constexpr double min(const double& x, const double & y, Targs... Fargs)
  {
    return min(x,min(y,Fargs...));
  }
  
  
  static inline double minmod(const double& x, const double& y){
    return 0.5*(sign(x)+sign(y))*min(std::fabs(x),std::fabs(y)) ;
  };
  
  static inline double mc2(const double& x, const double& y){
    return 0.5*(sign(x)+sign(y))*min(2.*std::fabs(x),2.*std::fabs(y), 0.5*std::fabs(x+y));};
  
}

#endif
