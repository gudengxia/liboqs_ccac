//
//  rng.h
//
//  Created by Bassham, Lawrence E (Fed) on 8/29/17.
//  Copyright © 2017 Bassham, Lawrence E (Fed). All rights reserved.
//

#ifndef rng_h
#define rng_h

#include <stdio.h>
#include "Alg.h"

/**********************************************************
* 函数： aigis_rand_get_sd_byts
* 功能：获取随机种子字节长度
* 返回：随机种子（字节）长度
*********************************************************/
PKC_ALG_API puchar_byts_t aigis_rand_get_sd_byts(void);

/**********************************************************
* 函数： aigis_rand_init
* 功能：根据输入的随机种子初始化随机数发生器
* 输入： s：随机种子
* s_byts：随机种子字节长度
* 返回：成功执行返回 0，否则返回错误代码（负数）
* 备注：在使用 aigis_rand_byts 函数之前必须调用该函数。可在调用密
* 钥生成函数之前调用该函数，以免影响性能。
*********************************************************/
PKC_ALG_API puchar_byts_t aigis_rand_init(unsigned char * s, unsigned long long s_byts);
/**********************************************************
* 函数： aigis_rand_byts
* 功能： 根据输入随机数字节长度输出随机数字节序列
* 输入： r_byts：随机数字节长度
* 输出： r：随机数
* 返回：成功执行返回 0，否则返回错误代码（负数）
*********************************************************/
PKC_ALG_API puchar_byts_t aigis_rand_byts(unsigned long long r_byts, unsigned char * r);
#endif /* rng_h */
