#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include "kem.h"
#include "params.h"
#include "indcpa.h"
#include "polyvec.h"
#include "poly.h"
#include "ntt.h"
#include "symmetric.h"
#include "randombytes.h"
#include "inverse.h"

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

static int16_t calc_diff(uint8_t *m1, uint8_t *m2)
{
	int16_t diff = 0;
	for(int i=0; i<CRYPTO_BYTES; i++)
	{
		if(diff < m1[i] - m2[i])
			diff = m1[i] - m2[i];
		if(diff < m2[i] - m1[i])
			diff = m2[i] - m1[i];
	}
	return diff;
}

static int16_t get_m_dash(uint8_t *m, uint8_t *m_dash) // uint8_t m[CRYPTO_BYTES]
{
	uint8_t pk[CRYPTO_PUBLICKEYBYTES];
    uint8_t sk[CRYPTO_SECRETKEYBYTES];
    uint8_t sk_dash[CRYPTO_SECRETKEYBYTES];
  	uint8_t ct[CRYPTO_CIPHERTEXTBYTES];

	for(int i=0;i<CRYPTO_SECRETKEYBYTES;i++)
	{
		sk_dash[i]=0;
	}

  	//Alice generates a public key
  	crypto_kem_keypair(pk, sk);

  	//Bob derives a secret key and creates a response
  	crypto_kem_enc(ct, m, pk);

  	//We create another secret key from pk
  	crypto_kem_get_sk(sk_dash, pk);

	// printf("\nSK:\n ");
	// for(int i=0;i<KYBER_SECRETKEYBYTES;i++)
	// {
	// 	printf("%d ",sk[i]);
	// }
	// printf("\nSK_dash:\n ");
	// for(int i=0;i<KYBER_SECRETKEYBYTES;i++)
	// {
	// 	printf("%d ",sk_dash[i]);
	// }	
	// printf("\n");

	// uint16_t flag =1;
	// for(int i=0;i<KYBER_SECRETKEYBYTES;i++)
	// {
	// 	if(sk[i]!=sk_dash[i])
	// 	{
	// 		flag=0;
			
	// 		printf("\nNot equal at index: %d\n", i);
	// 	}

	// }		
	// if(flag)
	// {
	// 	printf("\n-------------SK=SK_DASH!!-----------\n");
	// }
	// else
	// {
	// 	printf("\n-------------SK!=SK_DASH!!-----------\n");
	// }

  	//Alice uses Bobs response and new secret key to regenerate response
  	crypto_kem_dec(m_dash, ct, sk_dash);

  	return calc_diff(m,m_dash);
}

int main()
{
	printf("TEST_SK_DASH\n");

	uint8_t m[CRYPTO_BYTES],m_dash[CRYPTO_BYTES];
	
	printf("\n-----------------------------------------------------------\n");

	get_m_dash(m,m_dash);

	printf("m_orig = [ ");
	for(int i=0;i<CRYPTO_BYTES;i++)
		printf("%d ",m[i]);

	printf("]\n\n");
	
	printf("m_dash = [ ");
	for(int i=0;i<CRYPTO_BYTES;i++)
		printf("%d ",m_dash[i]);

	printf("]\n-----------------------------------------------------------\n");
	

	return 0;
}
