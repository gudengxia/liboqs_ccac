#ifndef OQS_SIG_MULAN_H
#define OQS_SIG_MULAN_H

#include <oqs/oqs.h>

#ifdef OQS_ENABLE_SIG_MULAN
#define OQS_SIG_mulan_length_public_key 1312
#define OQS_SIG_mulan_length_secret_key 3056
#define OQS_SIG_mulan_length_signature 2573

OQS_SIG *OQS_SIG_mulan_new(void);
OQS_API OQS_STATUS OQS_SIG_mulan_keypair(uint8_t *public_key, uint8_t *secret_key);
OQS_API OQS_STATUS OQS_SIG_mulan_sign(uint8_t *signature, size_t *signature_len, const uint8_t *message, size_t message_len, const uint8_t *secret_key);
OQS_API OQS_STATUS OQS_SIG_mulan_verify(const uint8_t *message, size_t message_len, const uint8_t *signature, size_t signature_len, const uint8_t *public_key);
#endif

#endif
