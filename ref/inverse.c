#include <stdio.h>
#include <stdint.h>
#include "params.h"
#include "poly.h"
#include "polyvec.h"
#include "ntt.h"
#include <stdlib.h>

void poly_print(poly *a)
{
	for(int i=0 ; i<KYBER_N ;i++)
	{
		printf("%d ", a->coeffs[i]);
	}
	printf("\n");
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
void clean(poly *a)
{
    for(int i=0; i<KYBER_N; i++)
    {
        a->coeffs[i]=0;
    }
}


/*
* Sets a = b
*/
void poly_setEqual(poly *a, poly *b)
{
    for(int i=0; i<KYBER_N; i++)
    {
        a->coeffs[i] = b->coeffs[i];
    }
}

void poly_red(poly *a)
{
	for(int i=0; i<256; i++)
		a->coeffs[i] = (KYBER_Q + (a->coeffs[i] % KYBER_Q)) % KYBER_Q;
}

static void poly_mul(poly *r, poly *a, poly *b)
{
	poly x,y;
	poly_setEqual(&x,a);
	poly_setEqual(&y,b);
	poly_ntt(&x);
	poly_ntt(&y);

	poly_basemul_montgomery(r,&x,&y);

	poly_invntt_tomont(r);
	poly_red(r);
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

static void poly_div(poly *q, poly *r, poly *a, poly *b, uint8_t flag)
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
	if(flag)   //ensure d>=1
	{
		if(d==0)
		{
			clean(r);
			poly_setEqual(q,a);
			return;
		}
		q->coeffs[256-d] = z;
		for(int j=0; j<d; j++){
			tmp = (q->coeffs[256-d])*(b->coeffs[j]) % KYBER_Q;
			r->coeffs[256-d+j] -= tmp;
		}
		poly_red(r);
	}

	for(int i=255; i>=d; i--){
		tmp =  z*(r->coeffs[i]) % KYBER_Q;
		q->coeffs[i-d] = tmp;
		for(int j=0; i-d+j<256; j++){
			tmp = (q->coeffs[i-d])*(b->coeffs[j]) % KYBER_Q;
			r->coeffs[i-d+j] -= tmp;
		}
		poly_red(r);
	}
}


/*
* Finds STD for which sa+tb=d
*/
static void poly_extended_gcd(poly *a, poly *b, poly *s, poly *t, poly *d)
{
    poly *s1 = malloc(sizeof(poly));
    poly *s2 = malloc(sizeof(poly));
    poly *t1 = malloc(sizeof(poly));
    poly *t2 = malloc(sizeof(poly));
    poly *q = malloc(sizeof(poly));
    poly *r = malloc(sizeof(poly));
    clean(t1); clean(t2);
    clean(s1); clean(s2);
    clean(q); clean(r);
    s2->coeffs[0] = 1;
    t1->coeffs[0] = 1;

    while(!check_zero(b))
    {
		clean(q); clean(r);
        poly_div(q,r,a,b,0);

        poly *temp = malloc(sizeof(poly));
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

static void poly_gcd(poly *a, poly *b, poly *x, poly *y, poly *d)
{
	poly d0,d1,d2,x0,x1,x2,y0,y1,y2,q,r,tmp;
	poly_setEqual(&d0,a);
	poly_setEqual(&d1,b);
	clean(&x0); clean(&x1);
	clean(&y0); clean(&y1);
	x0.coeffs[0] = 1;
	y1.coeffs[0] = 1;

	while(!check_zero(&d1)){
		poly_div(&q,&r,&d0,&d1,0);

		poly_setEqual(&d2,&d1); poly_setEqual(&x2,&x1); poly_setEqual(&y2,&y1);

		poly_setEqual(&d1,&r);

		poly_mul(&tmp,&q,&x1);
		poly_sub(&x1,&x0,&tmp);

		poly_mul(&tmp,&q,&y1);
		poly_sub(&y1,&y0,&tmp);

		poly_red(&x1);
		poly_red(&y1);

		poly_setEqual(&d0,&d2); poly_setEqual(&x0,&x2); poly_setEqual(&y0,&y2);
	}

	poly_setEqual(d,&d0); poly_setEqual(x,&x0); poly_setEqual(y,&y0);
}

static void poly_inv(poly *b, poly *c)
{
	poly a;
	clean(&a);
	a.coeffs[0]=1;

	poly_red(b);

	poly q,r;
	clean(&q);
	clean(&r);
	poly_div(&q,&r,&a,b,1);

	poly s,t,d;
	clean(&s);
	clean(&t);
	clean(&d);
	poly_gcd(b,&r,&s,&t,&d);

	poly temp;
	clean(&temp);
	poly_mul(&temp,&t,&q);
	poly_sub(c,&s,&temp);
	poly_red(c);
}

int main()
{
	poly a,b,s,t,d;
	clean(&a);
	clean(&b);
	clean(&s);
	clean(&t);
	clean(&d);	

  	a.coeffs[4] = 1;
	a.coeffs[1] = 0;
	a.coeffs[0] = 0;
  	b.coeffs[0] = 1;
	b.coeffs[1] = 0; 
  	b.coeffs[2] = 0;

 //  	poly_extended_gcd(&a,&b,&s,&t,&d);

 //  	poly_print(&s);
	// poly_print(&t);
	// poly_print(&d);


	poly_inv(&a,&s);
	poly_print(&s);

  	return 0;
}


// int main()
// {
// 	poly a,b,q,r,c;

// 	for(int i=0; i<250; i++){
// 		a.coeffs[i] = (KYBER_Q + (a.coeffs[i] % KYBER_Q)) % KYBER_Q;
// 		b.coeffs[i] = (KYBER_Q + (b.coeffs[i] % KYBER_Q)) % KYBER_Q;
// 	}

// 	for(int i=100; i<256; i++)
// 		b.coeffs[i] = 0;

//   	poly_div(&q,&r,&a,&b);

//   	test_poly_mul(&c,&b,&q);

//   	for(int i=0; i<256; i++)
//   		c.coeffs[i] = (KYBER_Q + (c.coeffs[i] + r.coeffs[i]) % KYBER_Q) % KYBER_Q;

//   	for(int i=0; i<256; i++)
//   		if(c.coeffs[i] % KYBER_Q != a.coeffs[i] % KYBER_Q){
//   			printf("FAIL");
//   			return 0;
//   		}
//   	printf("SUCCESS.");

//   	return 0;
// }