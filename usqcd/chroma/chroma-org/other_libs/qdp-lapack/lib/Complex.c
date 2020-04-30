/* $Id: Complex.c,v 1.4 2009-10-21 20:50:57 kostas Exp $ */
/* ================================================
 *
 * Based on slu_complex.c implemented by Sherry Li.
 *
 * This file defines common arithmetic operations for complex type.
 *
 * Note: _primme prefix is given to avoid conflicts with
 * some f2c and g2c implementations of z_* functions
 * A. Stathopoulos, Sept 15, 2006
 *
 * ================================================ */

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "qdp-lapack_Complex.h"

/*------------------------------------------------------------------------
 * Single Precision complex
 *-----------------------------------------------------------------------*/
/* Complex multiplication c= a*b */
void c_mul_primme(Complex_C *c, Complex_C *a, Complex_C *b)
{
	c->r = a->r*b->r - a->i*b->i ;
	c->i = a->r*b->i + a->i*b->r ;
}

/* Complex Division c = a/b */
void c_div_primme(Complex_C *c, Complex_C *a, Complex_C *b)
{
	float ratio, den;
	float abr, abi, cr, ci;

	if( (abr = b->r) < 0.)
		abr = - abr;
	if( (abi = b->i) < 0.)
		abi = - abi;
	if( abr <= abi ) {
		if (abi == 0) {
			fprintf(stderr, "c_div_primme.c: division by zero\n");
			exit(-1);
		}
		ratio = b->r / b->i ;
		den = b->i * (1 + ratio*ratio);
		cr = (a->r*ratio + a->i) / den;
		ci = (a->i*ratio - a->r) / den;
	} else {
		ratio = b->i / b->r ;
		den = b->r * (1 + ratio*ratio);
		cr = (a->r + a->i*ratio) / den;
		ci = (a->i - a->r*ratio) / den;
	}
	c->r = cr;
	c->i = ci;
}


/* Returns sqrt(c.r^2 + c.i^2) */
float c_abs_primme(Complex_C c)
{
	float temp;
	float real = c.r;
	float imag = c.i;

	if (real < 0) real = -real;
	if (imag < 0) imag = -imag;
	if (imag > real) {
		temp = real;
		real = imag;
		imag = temp;
	}
	if ((real+imag) == real) return(real);

	temp = imag/real;
	temp = real*sqrt(1.0 + temp*temp);  /*overflow!!*/
	return (temp);
}


/* Approximates the abs */
/* Returns abs(c.r) + abs(c.i) */
float c_abs1_primme(Complex_C c)
{
	float real = c.r;
	float imag = c.i;

	if (real < 0) real = -real;
	if (imag < 0) imag = -imag;

	return (real + imag);
}

/* Return the exponentiation */
void c_exp_primme(Complex_C *r, Complex_C *z)
{
	float expx;

	expx = exp(z->r);
	r->r = expx * cos(z->i);
	r->i = expx * sin(z->i);
}

/* Return the complex conjugate */
void s_cnjg_primme(Complex_C *r, Complex_C *z)
{
	r->r = z->r;
	r->i = -z->i;
}

/* Return the imaginary part */
float s_imag_primme(Complex_C *z)
{
	return (z->i);
}

/*------------------------------------------------------------------------
 * Double Precision complex
 *-----------------------------------------------------------------------*/
/* Complex multiplication c= a*b */
void z_mul_primme(Complex_Z *c, Complex_Z *a, Complex_Z *b)
{       
        c->r = a->r*b->r - a->i*b->i ;
        c->i = a->r*b->i + a->i*b->r ;
}                       

/* Complex Division c = a/b */
void z_div_primme(Complex_Z *c, Complex_Z *a, Complex_Z *b)
{
	double ratio, den;
	double abr, abi, cr, ci;

	if( (abr = b->r) < 0.)
		abr = - abr;
	if( (abi = b->i) < 0.)
		abi = - abi;
	if( abr <= abi ) {
		if (abi == 0) {
			fprintf(stderr, "z_div_primme.c: division by zero\n");
			exit(-1);
		}
		ratio = b->r / b->i ;
		den = b->i * (1 + ratio*ratio);
		cr = (a->r*ratio + a->i) / den;
		ci = (a->i*ratio - a->r) / den;
	} else {
		ratio = b->i / b->r ;
		den = b->r * (1 + ratio*ratio);
		cr = (a->r + a->i*ratio) / den;
		ci = (a->i - a->r*ratio) / den;
	}
	c->r = cr;
	c->i = ci;
}


/* Returns sqrt(z.r^2 + z.i^2) */
double z_abs_primme(Complex_Z z)
{
	double temp;
	double real = z.r;
	double imag = z.i;

	if (real < 0) real = -real;
	if (imag < 0) imag = -imag;
	if (imag > real) {
		temp = real;
		real = imag;
		imag = temp;
	}
	if ((real+imag) == real) return(real);

	temp = imag/real;
	temp = real*sqrt(1.0 + temp*temp);  /*overflow!!*/
	return (temp);
}


/* Approximates the abs */
/* Returns abs(z.r) + abs(z.i) */
double z_abs1_primme(Complex_Z z)
{
	double real = z.r;
	double imag = z.i;

	if (real < 0) real = -real;
	if (imag < 0) imag = -imag;

	return (real + imag);
}

/* Return the exponentiation */
void z_exp_primme(Complex_Z *r, Complex_Z *z)
{
	double expx;

	expx = exp(z->r);
	r->r = expx * cos(z->i);
	r->i = expx * sin(z->i);
}

/* Return the complex conjugate */
void d_cnjg_primme(Complex_Z *r, Complex_Z *z)
{
	r->r = z->r;
	r->i = -z->i;
}

/* Return the imaginary part */
double d_imag_primme(Complex_Z *z)
{
	return (z->i);
}
