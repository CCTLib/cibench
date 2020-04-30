// $Id: generic_blas_vaxpy3.h,v 1.3 2009-09-15 20:48:42 bjoo Exp $

/*! @file
 *  @brief Generic Scalar VAXPY routine
 *
 */

#ifndef QDP_GENERIC_BLAS_VAXPY3
#define QDP_GENERIC_BLAS_VAXPY3

namespace QDP {
// (Vector) out = (Scalar) (*scalep) * (Vector) InScale + (Vector) Add
inline
void vaxpy3(REAL *Out,REAL  *scalep,REAL *InScale, REAL *Add,int n_4vec)
{

  REAL a=*scalep;

  int len=n_4vec*24;
  for(int i=0; i < len; i++) {
    Out[i]=a*InScale[i]+Add[i];
  }

#if 0
  cout << "generic vaxpy3" << endl;

   double a;
   double x0r;
   double x0i;
  
   double x1r;
   double x1i;
  
   double x2r;
   double x2i;
  
   double y0r;
   double y0i;
  
   double y1r;
   double y1i;
  
   double y2r;
   double y2i;
  
   double z0r;
   double z0i;
  
   double z1r;
   double z1i;
  
   double z2r;
   double z2i;
  
  a = *scalep;
  
   int index_x = 0;
   int index_y = 0;
   int index_z = 0;
  
   int counter;

  // OK. Each n_4vec a length 24 vector. Here we've unrolled
  // 6 deep, so we need to do 4* this 
  int len=4*n_4vec;

  for( counter = 0; counter < len; counter++) {
    x0r = (double)InScale[index_x++];
    y0r = (double)Add[index_y++];
    z0r = a*x0r + y0r;
    Out[index_z++] =(REAL) z0r;
    
    x0i = (double)InScale[index_x++];
    y0i = (double)Add[index_y++];
    z0i = a*x0i + y0i;
    Out[index_z++] =(REAL) z0i;
    
    x1r = (double)InScale[index_x++];
    y1r = (double)Add[index_y++];
    z1r = a*x1r + y1r;
    Out[index_z++] = (REAL)z1r;
    
    x1i = (double)InScale[index_x++];
    y1i = (double)Add[index_y++];
    z1i = a*x1i + y1i;
    Out[index_z++] = (REAL)z1i;
    
    x2r = (double)InScale[index_x++];     
    y2r = (double)Add[index_y++];
    z2r = a*x2r + y2r;
    Out[index_z++] = (REAL)z2r;
    
    x2i = (double)InScale[index_x++];
    y2i = (double)Add[index_y++];
    z2i = a*x2i + y2i;  
    Out[index_z++] = (REAL)z2i;
  }
#endif

}


} // namespace QDP;

#endif // guard
