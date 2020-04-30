#include "sse_mult_su3_an.h"
#include <xmmintrin.h>

#ifdef __cplusplus
extern "C" { 
#endif

void
intrin_sse_mult_su3_an(su3_matrixf* aa, su3_matrixf* bb, su3_matrixf* cc)
{

	 /* XMM Variables */
	 __m128 xmm2, xmm3, xmm0, xmm1, xmm6, xmm7, xmm4, xmm5;

	xmm0 = _mm_loadl_pi(xmm0, (__m64 *)&((bb)->e[0][0]) );
	xmm1 = _mm_loadl_pi(xmm1, (__m64 *)&((bb)->e[1][0]) );
	xmm2 = _mm_loadl_pi(xmm2, (__m64 *)&((bb)->e[2][0]) );
	xmm0 = _mm_loadh_pi(xmm0, (__m64 *)&((bb)->e[0][1]) );
	xmm1 = _mm_loadh_pi(xmm1, (__m64 *)&((bb)->e[1][1]) );
	xmm2 = _mm_loadh_pi(xmm2, (__m64 *)&((bb)->e[2][1]) );
	xmm3 = _mm_load_ss((float *)&((aa)->e[0][0].real) );
	xmm6 = _mm_load_ss((float *)&((aa)->e[1][0].real) );
	xmm4 = _mm_load_ss((float *)&((aa)->e[0][1].real) );
	xmm7 = _mm_load_ss((float *)&((aa)->e[2][1].real) );
	xmm5 = _mm_load_ss((float *)&((aa)->e[0][2].real) );
	xmm3 = _mm_shuffle_ps( xmm3, xmm3, 0x00 );
	xmm6 = _mm_shuffle_ps( xmm6, xmm6, 0x00 );
	xmm4 = _mm_shuffle_ps( xmm4, xmm4, 0x00 );
	xmm3 = _mm_mul_ps( xmm3, xmm0 );
	xmm7 = _mm_shuffle_ps( xmm7, xmm7, 0x00 );
	xmm6 = _mm_mul_ps( xmm6, xmm1 );
	xmm5 = _mm_shuffle_ps( xmm5, xmm5, 0x00 );
	xmm4 = _mm_mul_ps( xmm4, xmm0 );
	xmm3 = _mm_add_ps( xmm3, xmm6 );
	xmm7 = _mm_mul_ps( xmm7, xmm2 );
	xmm5 = _mm_mul_ps( xmm5, xmm0 );
	xmm4 = _mm_add_ps( xmm4, xmm7 );
	xmm6 = _mm_load_ss((float *)&((aa)->e[1][2].real) );
	xmm7 = _mm_load_ss((float *)&((aa)->e[2][0].real) );
	xmm6 = _mm_shuffle_ps( xmm6, xmm6, 0x00 );
	xmm7 = _mm_shuffle_ps( xmm7, xmm7, 0x00 );
	xmm6 = _mm_mul_ps( xmm6, xmm1 );
	xmm7 = _mm_mul_ps( xmm7, xmm2 );
	xmm5 = _mm_add_ps( xmm5, xmm6 );
	xmm3 = _mm_add_ps( xmm3, xmm7 );
	xmm6 = _mm_load_ss((float *)&((aa)->e[1][1].real) );
	xmm7 = _mm_load_ss((float *)&((aa)->e[2][2].real) );
	xmm6 = _mm_shuffle_ps( xmm6, xmm6, 0x00 );
	xmm7 = _mm_shuffle_ps( xmm7, xmm7, 0x00 );
	xmm6 = _mm_mul_ps( xmm6, xmm1 );
	xmm7 = _mm_mul_ps( xmm7, xmm2 );
	xmm4 = _mm_add_ps( xmm4, xmm6 );
	xmm5 = _mm_add_ps( xmm5, xmm7 );
	xmm6 = _mm_load_ss((float *)&((aa)->e[0][0].imag) );
	xmm7 = _mm_load_ss((float *)&((aa)->e[1][1].imag) );
	xmm0 = _mm_shuffle_ps( xmm0, xmm0, 0xb1 );
	xmm1 = _mm_shuffle_ps( xmm1, xmm1, 0xb1 );
	xmm2 = _mm_shuffle_ps( xmm2, xmm2, 0xb1 );
	xmm6 = _mm_shuffle_ps( xmm6, xmm6, 0x00 );
	xmm7 = _mm_shuffle_ps( xmm7, xmm7, 0x00 );
	 	 xmm0 = _mm_xor_ps( xmm0, _sse_sgn24.xmm );
	 	 xmm1 = _mm_xor_ps( xmm1, _sse_sgn24.xmm );
	 	 xmm2 = _mm_xor_ps( xmm2, _sse_sgn24.xmm );
	xmm6 = _mm_mul_ps( xmm6, xmm0 );
	xmm7 = _mm_mul_ps( xmm7, xmm1 );
	xmm3 = _mm_add_ps( xmm3, xmm6 );
	xmm4 = _mm_add_ps( xmm4, xmm7 );
	xmm6 = _mm_load_ss((float *)&((aa)->e[2][2].imag) );
	xmm7 = _mm_load_ss((float *)&((aa)->e[0][1].imag) );
	xmm6 = _mm_shuffle_ps( xmm6, xmm6, 0x00 );
	xmm7 = _mm_shuffle_ps( xmm7, xmm7, 0x00 );
	xmm6 = _mm_mul_ps( xmm6, xmm2 );
	xmm7 = _mm_mul_ps( xmm7, xmm0 );
	xmm5 = _mm_add_ps( xmm5, xmm6 );
	xmm4 = _mm_add_ps( xmm4, xmm7 );
	xmm6 = _mm_load_ss((float *)&((aa)->e[1][0].imag) );
	xmm7 = _mm_load_ss((float *)&((aa)->e[0][2].imag) );
	xmm6 = _mm_shuffle_ps( xmm6, xmm6, 0x00 );
	xmm7 = _mm_shuffle_ps( xmm7, xmm7, 0x00 );
	xmm6 = _mm_mul_ps( xmm6, xmm1 );
	xmm7 = _mm_mul_ps( xmm7, xmm0 );
	xmm3 = _mm_add_ps( xmm3, xmm6 );
	xmm5 = _mm_add_ps( xmm5, xmm7 );
	xmm0 = _mm_load_ss((float *)&((aa)->e[2][0].imag) );
	xmm6 = _mm_load_ss((float *)&((aa)->e[1][2].imag) );
	xmm7 = _mm_load_ss((float *)&((aa)->e[2][1].imag) );
	xmm0 = _mm_shuffle_ps( xmm0, xmm0, 0x00 );
	xmm6 = _mm_shuffle_ps( xmm6, xmm6, 0x00 );
	xmm7 = _mm_shuffle_ps( xmm7, xmm7, 0x00 );
	xmm0 = _mm_mul_ps( xmm0, xmm2 );
	xmm6 = _mm_mul_ps( xmm6, xmm1 );
	xmm7 = _mm_mul_ps( xmm7, xmm2 );
	xmm3 = _mm_add_ps( xmm3, xmm0 );
	xmm5 = _mm_add_ps( xmm5, xmm6 );
	xmm4 = _mm_add_ps( xmm4, xmm7 );
	_mm_storel_pi((__m64 *)&((cc)->e[0][0]), xmm3 );
	_mm_storel_pi((__m64 *)&((cc)->e[1][0]), xmm4 );
	_mm_storel_pi((__m64 *)&((cc)->e[2][0]), xmm5 );
	_mm_storeh_pi((__m64 *)&((cc)->e[0][1]), xmm3 );
	_mm_storeh_pi((__m64 *)&((cc)->e[1][1]), xmm4 );
	_mm_storeh_pi((__m64 *)&((cc)->e[2][1]), xmm5 );
	xmm0 = _mm_loadl_pi(xmm0, (__m64 *)&((bb)->e[0][2]) );
	xmm1 = _mm_loadl_pi(xmm1, (__m64 *)&((bb)->e[1][2]) );
	xmm2 = _mm_loadl_pi(xmm2, (__m64 *)&((bb)->e[2][2]) );
	xmm0 = _mm_shuffle_ps( xmm0, xmm0, 0x44 );
	xmm1 = _mm_shuffle_ps( xmm1, xmm1, 0x44 );
	xmm2 = _mm_shuffle_ps( xmm2, xmm2, 0x44 );
	xmm3 = _mm_load_ss((float *)&((aa)->e[0][0].real) );
	xmm7 = _mm_load_ss((float *)&((aa)->e[0][1].real) );
	xmm3 = _mm_shuffle_ps( xmm3, xmm7, 0x00 );
	xmm4 = _mm_load_ss((float *)&((aa)->e[1][0].real) );
	xmm7 = _mm_load_ss((float *)&((aa)->e[1][1].real) );
	xmm4 = _mm_shuffle_ps( xmm4, xmm7, 0x00 );
	xmm3 = _mm_mul_ps( xmm3, xmm0 );
	xmm4 = _mm_mul_ps( xmm4, xmm1 );
	xmm3 = _mm_add_ps( xmm3, xmm4 );
	xmm5 = _mm_load_ss((float *)&((aa)->e[2][0].real) );
	xmm7 = _mm_load_ss((float *)&((aa)->e[2][1].real) );
	xmm5 = _mm_shuffle_ps( xmm5, xmm7, 0x00 );
	xmm5 = _mm_mul_ps( xmm5, xmm2 );
	xmm3 = _mm_add_ps( xmm3, xmm5 );
	xmm1 = _mm_shuffle_ps( xmm1, xmm0, 0x44 );
	xmm7 = _mm_load_ss((float *)&((aa)->e[0][2].real) );
	xmm6 = _mm_load_ss((float *)&((aa)->e[1][2].real) );
	xmm6 = _mm_shuffle_ps( xmm6, xmm7, 0x00 );
	xmm6 = _mm_mul_ps( xmm6, xmm1 );
	xmm0 = _mm_shuffle_ps( xmm0, xmm0, 0xB1 );
	 	 xmm0 = _mm_xor_ps( xmm0, _sse_sgn24.xmm );
	xmm1 = _mm_shuffle_ps( xmm1, xmm1, 0x11 );
	 	 xmm1 = _mm_xor_ps( xmm1, _sse_sgn24.xmm );
	xmm2 = _mm_shuffle_ps( xmm2, xmm2, 0xB1 );
	 	 xmm2 = _mm_xor_ps( xmm2, _sse_sgn24.xmm );
	xmm4 = _mm_load_ss((float *)&((aa)->e[0][0].imag) );
	xmm7 = _mm_load_ss((float *)&((aa)->e[0][1].imag) );
	xmm4 = _mm_shuffle_ps( xmm4, xmm7, 0x00 );
	xmm4 = _mm_mul_ps( xmm4, xmm0 );
	xmm3 = _mm_add_ps( xmm3, xmm4 );
	xmm5 = _mm_load_ss((float *)&((aa)->e[1][0].imag) );
	xmm7 = _mm_load_ss((float *)&((aa)->e[1][1].imag) );
	xmm5 = _mm_shuffle_ps( xmm5, xmm7, 0x00 );
	xmm5 = _mm_mul_ps( xmm5, xmm1 );
	xmm3 = _mm_add_ps( xmm3, xmm5 );
	xmm5 = _mm_load_ss((float *)&((aa)->e[2][0].imag) );
	xmm7 = _mm_load_ss((float *)&((aa)->e[2][1].imag) );
	xmm5 = _mm_shuffle_ps( xmm5, xmm7, 0x00 );
	xmm5 = _mm_mul_ps( xmm5, xmm2 );
	xmm3 = _mm_add_ps( xmm3, xmm5 );
	_mm_storel_pi((__m64 *)&((cc)->e[0][2]), xmm3 );
	_mm_storeh_pi((__m64 *)&((cc)->e[1][2]), xmm3 );
	xmm1 = _mm_shuffle_ps( xmm1, xmm0, 0x44 );
	xmm7 = _mm_load_ss((float *)&((aa)->e[0][2].imag) );
	xmm5 = _mm_load_ss((float *)&((aa)->e[1][2].imag) );
	xmm5 = _mm_shuffle_ps( xmm5, xmm7, 0x00 );
	xmm5 = _mm_mul_ps( xmm5, xmm1 );
	xmm6 = _mm_add_ps( xmm6, xmm5 );
	xmm2 = _mm_shuffle_ps( xmm2, xmm2, 0xB4 );
	 	 xmm2 = _mm_xor_ps( xmm2, _sse_sgn3.xmm );
	xmm7 = _mm_loadl_pi(xmm7, (__m64 *)&((aa)->e[2][2]) );
	xmm7 = _mm_shuffle_ps( xmm7, xmm7, 0x05 );
	xmm7 = _mm_mul_ps( xmm7, xmm2 );
	xmm6 = _mm_add_ps( xmm6, xmm7 );
	xmm7 = xmm6 ; 
	xmm7 = _mm_shuffle_ps( xmm7, xmm7, 0xEE );
	xmm6 = _mm_add_ps( xmm6, xmm7 );
	_mm_storel_pi((__m64 *)&((cc)->e[2][2]), xmm6 );
}

#ifdef __cplusplus
}; 
#endif

