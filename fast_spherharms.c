#include <stdio.h>
#include <gsl/gsl_sf_legendre.h>
#include <gsl/gsl_sf_trig.h>
#include <gsl/gsl_pow_int.h>
#include <math.h>


double sin_term(int l,int m,double phi,double theta){
 double out = gsl_pow_int(-1.0,m)*gsl_sf_legendre_sphPlm(l,m,gsl_sf_cos(theta))*gsl_sf_sin(m*phi);
 return out;
}

double cos_term(int l,int m, double phi, double theta){
  double out= gsl_pow_int(-1.0,m)*gsl_sf_legendre_sphPlm(l,m,gsl_sf_cos(theta))*gsl_sf_cos(m*phi);
  return out; 
}




float acos_fast(float x) {
  float negate = (float)x < 0;
  x =  x < 0 ? -x : x;
  float ret = -0.0187293;
  ret = ret * x;
  ret = ret + 0.0742610;
  ret = ret * x;
  ret = ret - 0.2121144;
  ret = ret * x;
  ret = ret + 1.5707288;
  ret = ret * sqrt(1.0-x);
  ret = ret - 2 * negate * ret;
  return negate * 3.14159265358979 + ret;
}


float fast_grf(int L, float x, float y,float z,float rands[]){
  float result = 0.0;
  float phi =atan2(y,x);
  float theta = acos_fast(z);
  float temp ;  
  register int l = 0;
  register int m = 1; 
  for (l=0; l<= L;l++){
    result += gsl_sf_legendre_sphPlm(l,0,gsl_sf_cos(theta))*rands[2*L*l];
    temp = 0.0;
    for(m=1;m<=l;m++){
      temp += (1 - ((m & 1) << 1))*gsl_sf_legendre_sphPlm(l,m,gsl_sf_cos(theta))*(rands[2*m+2*L*l+l]*gsl_sf_sin(m*phi)+rands[1+2*m+2*L*l]*gsl_sf_cos(m*phi));
    }
    result += 1.41421*temp; 
  }
  return result;
}

