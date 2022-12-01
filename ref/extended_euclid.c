#include <stdio.h>
#include <stdint.h>
#include "params.h"
#include "poly.h"
#include "polyvec.h"

static void poly_mul(poly *r, poly *a, poly *b)
{
	poly_ntt(a);
	poly_ntt(b);

	poly_basemul_montgomery(r,a,b);

	poly_invntt_tomont(r);
	poly_reduce(r);
}

/*
* Checks if poly  is identically zero or not
*/
uint8_t check_zero(poly *a)
{
    uint8_t res = 1;
    for(int i=0; i<KYBER_N; i++)
    {
        if(a->coeffs[i]!=0)
            res = 0;
    }
    return res;
}
/*
* Sets all coeffs equal to 0
*/
static void clean(poly *a)
{
    for(int i=0; i<KYBER_N; i++)
    {
        a->coeffs[i]=0;
    }
}


/*
* Sets a = b
*/
static void poly_setEqual(poly *a, poly *b)
{
    for(int i=0; i<KYBER_N; i++)
    {
        a->coeffs[i] = b->coeffs[i];
    }
}
/*
* Finds STD for which sa+tb=d
*/
static void poly_extended_gcd(poly *a, poly *b, poly *s, poly *t, poly *d)
{
    poly *s1;
    poly *s2;
    poly *t1;
    poly *t2;
    poly *q;
    poly *r;
    clean(t1);
    clean(t2);
    clean(s1);
    clean(s2);
    clean(q);
    clean(r);
    s2->coeffs[0] = 1;
    t1->coeffs[1] = 1;
    while(check_zero(a))
    {
        clean(q);
        clean(r);
        poly_divide(a,b,q,r);
        poly *temp;
        clean(temp);
        poly_mul(temp,q,s1);
        poly_sub(s,s2,temp);
        clean(temp);
        poly_mul(temp,q,t1);
        poly_sub(t,t2,temp);
        poly_setEqual(a,b);
        poly_setEqual(b,r);
        poly_setEqual(s2,s1);
        poly_setEqual(s1,s);
        poly_setEqual(t2,t1);
        poly_setEqual(t1,t);
    }
    poly_setEqual(d,a);
    poly_setEqual(s,s2);
    poly_setEqual(t,t2);
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

  	poly_mul(&r,&a,&b);

  	for(int i=0; i<KYBER_N; i++)
    	printf("'%d', ",r.coeffs[i]);

  	return 0;
}