#include <stdio.h>
#include <stdint.h>
#include "params.h"
#include "poly.h"
#include "polyvec.h"

static void test_poly_mul(poly *r, poly *a, poly *b)
{
	poly_ntt(a);
	poly_ntt(b);

	poly_basemul_montgomery(r,a,b);

	poly_invntt_tomont(r);
	poly_reduce(r);
}

int main()
{
	poly a,b,r;

  	for(int i=0; i<KYBER_N; i++){
    	a.coeffs[i] = 0;
    	b.coeffs[i] = 0;
  	}

  	a.coeffs[1] = 2;
  	b.coeffs[1] = 3;
  	b.coeffs[3] = 2;

  	test_poly_mul(&r,&a,&b);

  	for(int i=0; i<KYBER_N; i++)
    	printf("'%d', ",r.coeffs[i]);

  	return 0;
}