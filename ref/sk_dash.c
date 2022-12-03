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
#include "inverse.c"

void indcpa_get_sk(uint8_t sk[KYBER_INDCPA_SECRETKEYBYTES], 
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

int crypto_kem_get_sk(uint8_t *sk,
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

int16_t calc_diff(uint8_t *m1, uint8_t *m2)
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

int16_t get_error(uint8_t *m) // uint8_t m[CRYPTO_BYTES]
{
	uint8_t pk[CRYPTO_PUBLICKEYBYTES];
    uint8_t sk[CRYPTO_SECRETKEYBYTES];
    uint8_t sk_dash[CRYPTO_SECRETKEYBYTES];
  	uint8_t ct[CRYPTO_CIPHERTEXTBYTES];
  	uint8_t m_dash[CRYPTO_BYTES];

  	//Alice generates a public key
  	crypto_kem_keypair(pk, sk);

  	//Bob derives a secret key and creates a response
  	crypto_kem_enc(ct, m, pk);

  	//We create another secret key from pk
  	crypto_kem_get_sk(sk_dash, pk);

  	//Alice uses Bobs response and new secret key to regenerate response
  	crypto_kem_dec(m_dash, ct, sk_dash);

  	return calc_diff(m,m_dash);
}

int main()
{
	printf("TEST_SK_DASH.");
	return 0;
}