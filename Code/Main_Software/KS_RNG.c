
#include "required_libs.h"

void KS_RNG(int *N, double *z, double *Xb, double *lambda);
int leftmost(double u, double l);
int rightmost(double u, double l);

void KS_RNG(int *N, double *z, double *Xb, double *lambda) {
  
  // Declarations
  int i, ok;
  double r, x, y, u1, u2;
  
  GetRNGstate();
  
  for( i=0; i<*N; i++ ) {
    
    r = fabs( z[i] - Xb[i] ) + 0.0000001; // Adding small positive constant to avoid dividing by zero
    ok = 0;
    
    // Rejection sampling
    while ( ok == 0 ) {
      
      y = -1.0;
      while( y < 0 ) { // Avoid the erroneous situation where y < 0
        x = pow( rnorm( 0.0, 1.0 ), 2.0 );
        y = 1.0 + ( x - sqrt(x*(4.0*r + x)) )/(2.0*r);
      }
      
      u1 = runif( 0.0, 1.0 );
      if( u1 <= 1.0/(1.0 + y) ) {
        lambda[i] = r/y;
      } else {
        lambda[i] = r*y;
      }
      
      u2 = runif( 0.0, 1.0 );
      if( lambda[i] > 4.0/3.0 ) {
        ok = rightmost(u2, lambda[i]);
      } else {
        ok = leftmost(u2, lambda[i]);
      }
      
    }
    
  }
  
  PutRNGstate();
  
  return;
    
}

// rightmost function for lambda acceptance
int rightmost(double u, double l) {
   
  double z = 1.0, x = exp( -l/2.0 ), j = 0;
  
  while ( 1 ) {
    
    j = j + 1.0;
    z = z - ( pow(j + 1.0, 2.0) * pow(x, pow(j + 1.0, 2.0) - 1.0) );
    if ( z > u ) return 1;
    
    j = j + 1.0;
    z = z + ( pow(j + 1.0, 2.0) * pow(x, pow(j + 1.0, 2.0) - 1.0) );
    if ( z < u ) return 0;
      
  }
  
}

// leftmost function for lambda acceptance
int leftmost(double u, double l) {
  
  double H = ( log(2.0)/2.0 ) + ( 2.5*log(M_PI) ) - ( 2.5*log(l) ) - ( pow(M_PI, 2.0)/(2.0*l) ) + ( l/2.0 ), 
    lu = log( u ), z = 1.0, x = exp( -pow(M_PI, 2.0) / (2.0*l) ), k = l / pow(M_PI, 2.0), j = 0.0;
  
  while ( 1 ) {
    
    j = j + 1.0;
    z = z - k*pow(x, ( pow(j, 2.0) - 1.0 ));
    if ( (H + log(z)) > lu ) return 1;
    
    j = j + 1.0;
    z = z + k*pow(x, ( pow(j, 2.0) - 1.0 ));
    if ( (H + log(z)) < lu ) return 0;
      
  }
  
}
