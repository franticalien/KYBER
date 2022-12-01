#include <stdio.h>
#include <stdint.h>
#include "params.h"
#include "poly.h"
#include "polyvec.h"
#include "ntt.h"

static void test_poly_mul(poly *r, poly *a, poly *b)
{
	poly_ntt(a);
	poly_ntt(b);

	poly_basemul_montgomery(r,a,b);

	poly_invntt_tomont(r);
	poly_reduce(r);
}


static int16_t fqinv(int16_t x)
{
	int32_t res = 1, z = x, b = KYBER_Q-2;
	while(b>0){
		if(b&1) res = res*z % KYBER_Q;
		z = z*z % KYBER_Q;
		b = b>>1;
	}
	return (int16_t)res;
}

static void poly_div(poly *q, poly *r, poly *a, poly *b)
{
	int16_t d = 255;
	while(b->coeffs[d]==0 && d>=0) d--;

	if(d<0){
		printf("ERROR: Division by 0."); 
		return;
	}

	int32_t tmp = 0;
	int32_t z = fqinv(b->coeffs[d]);

	for(int i=0; i<256; i++){
		q->coeffs[i] = 0;
		r->coeffs[i] = a->coeffs[i];
	}

	for(int i=255; i>=d; i--){
		tmp =  z*(r->coeffs[i]) % KYBER_Q;
		q->coeffs[i-d] = tmp;
		for(int j=0; i-d+j<256; j++){
			tmp = (q->coeffs[i-d])*(b->coeffs[j]) % KYBER_Q;
			r->coeffs[i-d+j] -= tmp;
			r->coeffs[i-d+j] = (KYBER_Q + r->coeffs[i-d+j]) % KYBER_Q;
		}
	}
}

int main()
{
	poly a,b,q,r,c;

	for(int i=0; i<250; i++){
		a.coeffs[i] = (KYBER_Q + (a.coeffs[i] % KYBER_Q)) % KYBER_Q;
		b.coeffs[i] = (KYBER_Q + (b.coeffs[i] % KYBER_Q)) % KYBER_Q;
	}

	for(int i=100; i<256; i++)
		b.coeffs[i] = 0;

  	poly_div(&q,&r,&a,&b);

  	test_poly_mul(&c,&b,&q);

  	for(int i=0; i<256; i++)
  		c.coeffs[i] = (KYBER_Q + (c.coeffs[i] + r.coeffs[i]) % KYBER_Q) % KYBER_Q;

  	for(int i=0; i<256; i++)
  		if(c.coeffs[i] % KYBER_Q != a.coeffs[i] % KYBER_Q){
  			printf("FAIL");
  			return 0;
  		}
  	printf("SUCCESS.");

  	return 0;
}