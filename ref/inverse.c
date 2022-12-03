#include <stdint.h>
#include <stdlib.h>
#include "params.h"
#include "poly.h"
#include "polyvec.h"
#include "ntt.h"
#include "inverse.h"

static int16_t fqinv(int16_t x)
{
	int32_t res = 1, z = ((x%KYBER_Q)+KYBER_Q)%KYBER_Q, b = KYBER_Q-2;
	while(b>0){
		if(b&1) res = res*z % KYBER_Q;
		z = z*z % KYBER_Q;
		b = b>>1;
	}
	return (int16_t)res;
}

static void baseinv(int16_t r[2], const int16_t a[2], int16_t zeta)
{
	int16_t d = fqmul(a[1], a[1]);
	d = fqmul(d, zeta);
	d -= fqmul(a[0], a[0]);
	d = fqinv(d);
	r[0] = fqmul(-a[0], d);
	r[1] = fqmul( a[1], d);
}

static void poly_ntt_inv(poly *r, const poly *a)
{
	unsigned int i;
	for(i=0;i<KYBER_N/4;i++) {
	   baseinv(&r->coeffs[4*i], &a->coeffs[4*i], zetas[64+i]);
	   baseinv(&r->coeffs[4*i+2], &a->coeffs[4*i+2], -zetas[64+i]);
  }
}

void mat_inv_2(polyvec *b, polyvec *a)
{
	if(KYBER_K != 2){
		printf("ERROR: KYBER_K != 2.");
		return;
	}

	poly den,t1,t2;
	poly_basemul_montgomery(&t1,&a[0].vec[0],&a[1].vec[1]);
	poly_basemul_montgomery(&t2,&a[1].vec[0],&a[0].vec[1]);
	poly_sub(&t1,&t1,&t2);
	poly_reduce(&t1);
	poly_ntt_inv(&den,&t1);

	poly_basemul_montgomery(&b[0].vec[0],&a[1].vec[1],&den);
	poly_basemul_montgomery(&b[1].vec[1],&a[0].vec[0],&den);
	poly_basemul_montgomery(&b[1].vec[0],&a[1].vec[0],&den);
	poly_basemul_montgomery(&b[0].vec[1],&a[0].vec[1],&den);

	for(int i=0; i<256; i++){
		b[0].vec[1].coeffs[i] *= -1;
		b[1].vec[0].coeffs[i] *= -1;
	}

	poly_reduce(&b[0].vec[0]); poly_reduce(&b[0].vec[1]);
	poly_reduce(&b[1].vec[0]); poly_reduce(&b[1].vec[1]);
}


