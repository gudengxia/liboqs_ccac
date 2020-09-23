/**************************************************************************
 * Alg.cpp (Version 1.1.1) created on May 12, 2019.
 *************************************************************************/
#include "Alg.h"
#include "api.h"

/*puchar_byts_t aigis_sig_get_pk_byts()
{
	return AIGIS_SIG_PUBLICKEYBYTES;
}
puchar_byts_t aigis_sig_get_sk_byts()
{
	return AIGIS_SIG_SECRETKEYBYTES;
}
puchar_byts_t aigis_sig_get_sn_byts()
{
	return AIGIS_SIG_BYTES;
}*/

int aigis_sig_keygen(puchar_t pk, puchar_byts_t* pk_byts,
	puchar_t sk, puchar_byts_t* sk_byts)
{
	*pk_byts = AIGIS_SIG_PUBLICKEYBYTES;
	*sk_byts = AIGIS_SIG_SECRETKEYBYTES;
	return msig_keygen(pk,sk);
}

int aigis_sig_sign(puchar_t sk, puchar_byts_t sk_byts,
	puchar_t m, puchar_byts_t m_byts,
	puchar_t sn, puchar_byts_t* sn_byts)
{
	if(sk_byts < 3056)
		return 1;
	return msig_sign(sk, m, m_byts,sn, sn_byts);
}

int aigis_sig_verf(
	puchar_t pk, puchar_byts_t pk_byts,
	puchar_t sn, puchar_byts_t sn_byts,
	puchar_t m, puchar_byts_t m_byts)
{
	if(pk_byts < 1312)
		return 0;
	return msig_verf(pk, sn, sn_byts, m, m_byts);
}