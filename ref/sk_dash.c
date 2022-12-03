#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "kem.h"
#include "params.h"
#include "indcpa.h"
#include "polyvec.h"
#include "poly.h"
#include "ntt.h"
#include "symmetric.h"
#include "randombytes.h"
#include "inverse.h"

#define MAX(x, y) (((x) > (y)) ? (x) : (y))
#define MIN(x, y) (((x) < (y)) ? (x) : (y))

/*
Sets sk = sk' given pk; sk'=(A^-1)t
*/
static void indcpa_get_sk(uint8_t sk[KYBER_INDCPA_SECRETKEYBYTES], 
				   const uint8_t pk[KYBER_INDCPA_PUBLICKEYBYTES])
{
	uint8_t seed[KYBER_SYMBYTES];
	polyvec pkpv, a[KYBER_K], a_inv[KYBER_K], skpv;

	unpack_pk(&pkpv,seed,pk);
	
	gen_matrix(a, seed,0);
	mat_inv_2(a_inv,a);

	// matrix-vector multiplication
    for(int i=0;i<KYBER_K;i++){
      polyvec_basemul_acc_montgomery(&skpv.vec[i], &a_inv[i], &pkpv);
      poly_tomont(&skpv.vec[i]);
    }

    pack_sk(sk, &skpv);
}

/*
Sets sk = sk' given pk; sk'=(A^-1)t and adds random symbytes to sk.
*/
static int crypto_kem_get_sk(uint8_t *sk,
                       uint8_t *pk)
{
	size_t i;
	indcpa_get_sk(sk,pk);
	for(i=0;i<KYBER_INDCPA_PUBLICKEYBYTES;i++)
    	sk[i+KYBER_INDCPA_SECRETKEYBYTES] = pk[i];
  	hash_h(sk+KYBER_SECRETKEYBYTES-2*KYBER_SYMBYTES, pk, KYBER_PUBLICKEYBYTES);
  	/* Value z for pseudo-random output on reject */
  	randombytes(sk+KYBER_SECRETKEYBYTES-KYBER_SYMBYTES, KYBER_SYMBYTES);
	return 0;

}

static int16_t get_diff(poly *p, poly *q, poly *e)
{
	poly_sub(e,p,q);
	poly_reduce(e);

	int16_t diff = 0;
	for(int i=0; i<256; i++)
		diff = MAX(diff, abs(e->coeffs[i]));

	return diff;
}

static int16_t get_error(uint8_t *m, uint8_t *m_dash, poly *e, int flag) // uint8_t m[CRYPTO_BYTES]
{
	uint8_t pk[CRYPTO_PUBLICKEYBYTES];
    uint8_t sk[CRYPTO_SECRETKEYBYTES];
    uint8_t sk_dash[CRYPTO_SECRETKEYBYTES];
  	uint8_t ct[CRYPTO_CIPHERTEXTBYTES];
  	poly p,q;

  	//Alice generates a public key
  	crypto_kem_keypair_mod(pk, sk, flag);

  	//Bob derives a secret key and creates a response
  	crypto_kem_enc_mod(ct, m, pk, &p);

  	//We create another secret key from pk
  	crypto_kem_get_sk(sk_dash, pk);

  	//Alice uses Bobs response and new secret key to regenerate response
  	crypto_kem_dec_mod(m_dash, ct, sk_dash, &q);

  	return get_diff(&p,&q,e);
}

int main()
{
	printf("\nTEST_SK_DASH");

	poly e;
	uint8_t m[CRYPTO_BYTES],m_dash[CRYPTO_BYTES],flag;

	flag = 0; // pk: t = As + e.
	//flag = 1; // pk: t = As + 0.
	int16_t diff = get_error(m,m_dash,&e,(int)flag);
	
	printf("\n-----------------------------------------------------------\n");
	printf("m_orig = [ ");
	for(int i=0;i<CRYPTO_BYTES;i++)
		printf("%d ",m[i]);
	printf("]\n\n");
	printf("m_dash = [ ");
	for(int i=0;i<CRYPTO_BYTES;i++)
		printf("%d ",m_dash[i]);
	printf("]\n\n");
	printf("error  = [ ");
	for(int i=0; i<256; i++)
		printf("%d ",e.coeffs[i]);
	printf("]\n\n||error|| = %d; [q/4] = %d",diff,KYBER_Q/4);
	printf("\n-----------------------------------------------------------\n");
	

	return 0;
}
