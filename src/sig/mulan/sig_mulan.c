#include <stdlib.h>

#include <oqs/sig_mulan.h>

#ifdef OQS_ENABLE_SIG_MULAN

OQS_SIG *OQS_SIG_mulan_new() {

	OQS_SIG *sig = malloc(sizeof(OQS_SIG));
	if (sig == NULL) {
		return NULL;
	}
	sig->method_name = OQS_SIG_alg_mulan;
	sig->alg_version = "fzhang test";

	sig->claimed_nist_level = 2;
	sig->euf_cma = true;

	sig->length_public_key = OQS_SIG_mulan_length_public_key;
	sig->length_secret_key = OQS_SIG_mulan_length_secret_key;
	sig->length_signature = OQS_SIG_mulan_length_signature;

	sig->keypair = OQS_SIG_mulan_keypair;
	sig->sign = OQS_SIG_mulan_sign;
	sig->verify = OQS_SIG_mulan_verify;

	return sig;
}

/*int PQCLEAN_DILITHIUM2_CLEAN_crypto_sign_keypair(uint8_t *pk, uint8_t *sk);
int PQCLEAN_DILITHIUM2_CLEAN_crypto_sign_signature(uint8_t *sig, size_t *siglen, const uint8_t *m, size_t mlen, const uint8_t *sk);
int PQCLEAN_DILITHIUM2_CLEAN_crypto_sign_verify(const uint8_t *sig, size_t siglen, const uint8_t *m, size_t mlen, const uint8_t *pk);

OQS_API OQS_STATUS OQS_SIG_dilithium_2_keypair(uint8_t *public_key, uint8_t *secret_key) {
	return (OQS_STATUS) PQCLEAN_DILITHIUM2_CLEAN_crypto_sign_keypair(public_key, secret_key);
}
OQS_API OQS_STATUS OQS_SIG_dilithium_2_sign(uint8_t *signature, size_t *signature_len, const uint8_t *message, size_t message_len, const uint8_t *secret_key) {
	return (OQS_STATUS) PQCLEAN_DILITHIUM2_CLEAN_crypto_sign_signature(signature, signature_len, message, message_len, secret_key);
}
OQS_API OQS_STATUS OQS_SIG_dilithium_2_verify(const uint8_t *message, size_t message_len, const uint8_t *signature, size_t signature_len, const uint8_t *public_key) {
	return (OQS_STATUS) PQCLEAN_DILITHIUM2_CLEAN_crypto_sign_verify(signature, signature_len, message, message_len, public_key);
}*/
int mulan_sig_keygen(unsigned char * pk, unsigned long long * pk_bytes, unsigned char * sk, unsigned long long * sk_bytes); 	
int mulan_sig_sign(unsigned char * sk, unsigned long long sk_bytes, unsigned char * m, unsigned long long mlen, unsigned char * sm, unsigned long long * smlen);
int mulan_sig_verf(unsigned char * pk, unsigned long long pk_bytes, unsigned char * sm, unsigned long long smlen, unsigned char * m, unsigned long long mlen);

OQS_API OQS_STATUS OQS_SIG_mulan_keypair(uint8_t *public_key, uint8_t *secret_key) {
	unsigned long long pk_bytes, sk_bytes;
	int ret = mulan_sig_keygen((unsigned char*)public_key, &pk_bytes, (unsigned char*)secret_key, &sk_bytes);
	if(ret == 0)
		return OQS_SUCCESS;
	else
		return OQS_ERROR;
	
}
OQS_API OQS_STATUS OQS_SIG_mulan_sign(uint8_t *signature, size_t *signature_len, const uint8_t *message, size_t message_len, const uint8_t *secret_key) {
	int ret = mulan_sig_sign((unsigned char*)secret_key, 3056, (unsigned char *)message, message_len, (unsigned char*)signature, (unsigned long long*)signature_len);
	if(ret == 0)
		return OQS_SUCCESS;
	else
		return OQS_ERROR;
	
}
OQS_API OQS_STATUS OQS_SIG_mulan_verify(const uint8_t *message, size_t message_len, const uint8_t *signature, size_t signature_len, const uint8_t *public_key) {
	int ret = mulan_sig_verf((unsigned char*)public_key, 1312,  (unsigned char*)signature, signature_len, (unsigned char*)message, message_len);
	if(ret == 1)
		return OQS_SUCCESS;
	else
		return OQS_ERROR;
}
#endif
