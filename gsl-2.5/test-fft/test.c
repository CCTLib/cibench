#include <stdio.h>
#include <math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_fft_complex.h>

#define REAL(z,i) ((z)[2*(i)])
#define IMAG(z,i) ((z)[2*(i)+1])

#define N (0x10000)
int
main (void)
{
  int i; double data[4*N];

  for (i = 0; i < N; i++)
    {
       REAL(data,i) = 0.0; IMAG(data,i) = 0.0;
    }

  REAL(data,0) = 1.0;

  for (i = 1; i <= 10; i++)
    {
       REAL(data,i) = REAL(data,N-i) = 1.0;
    }

#if 0
  for (i = 0; i < 128; i++)
    {
      printf ("%d %e %e\n", i, 
              REAL(data,i), IMAG(data,i));
    }
  printf ("\n");
#endif

  for (i = 0; i < 15000; i++)
    gsl_fft_complex_radix2_transform (data, 1, N, -1);

#if 0
  for (i = 0; i < 128; i++)
    {
      printf ("%d %e %e\n", i, 
              REAL(data,i)/sqrt(128), 
              IMAG(data,i)/sqrt(128));
    }
#endif

  return 0;
}
