
#include <immintrin.h>
#include "params.h"
#include "ntt.h"


static const uint64_t mont_invn = (((uint64_t)MONT*MONT % PARAM_Q) * (PARAM_Q-1) % PARAM_Q) * ((PARAM_Q-1) >> 8) % PARAM_Q;

#if PARAM_Q == 2021377
static const uint32_t zetas[PARAM_N] = {1562548,518470,697898,862629,1367459,1539276,1513857,1662806,929015,1757045,1879015,449873,75689,1125711,1680345,620849,769419,486664,1389778,658915,1319993,73499,1391732,1199964,291970,655587,966181,128755,288564,10420,1980158,1011904,1937906,838813,854780,1453936,1704819,1740984,86645,1360044,115556,1570480,1655800,272433,1245520,1190005,238406,1726139,1013693,1948648,1020399,1544116,1120075,656153,591869,1620799,275832,517427,1601944,1555925,1293833,1705829,1357642,142050,739420,1568070,1535360,740638,57925,1038012,65439,1844105,673379,1768997,924638,1986117,1394208,1277276,129269,1760277,1173604,1161770,1897168,807697,965038,1876057,1963820,1794916,924093,251419,168030,1073286,1902394,347156,1488477,511116,572755,1686880,268077,53223,1268228,579769,1043786,272581,1574784,1729984,568576,276296,1095755,282107,158374,915466,1569380,908136,972609,923797,466409,1762448,798650,436051,1275685,1122838,1862,1854194,1432015,1507507,1452715,1170924,137295,531590,556763,1442250,896280,320184,333460,1993546,622613,1352919,881664,1176558,1936677,2011958,1357750,534023,142791,40293,638104,1519860,1189220,1763667,792470,1813814,830483,1256948,1537350,64760,561409,823180,786453,1106713,1491299,1582163,822179,1663832,1269819,84100,780824,310495,1043416,763923,1440072,1308437,1369984,1027053,641681,932722,1248044,318540,1777818,702544,1566714,1301662,265980,696370,1576958,449193,1228202,1635455,1143957,1349609,120737,1115065,1815624,573533,10820,1911846,533321,1147868,1126927,145151,641139,275750,276830,1257214,988074,1857331,105366,1608247,1752751,817865,294374,1145376,1447053,647982,1517128,301974,233775,1669708,1146108,1913137,707228,1147423,349817,1972001,777351,1874015,964313,161863,1142539,1331457,1604014,1320129,1103939,1236477,447210,1613614,1666811,51306,383284,1573619,677023,994549,23785,210391,461525,1779756,430663,84620,1731642,1784991,147098,942182,1953450,1853187,1567373,1541031};

static const uint32_t zetas_inv[PARAM_N] = {480346,454004,168190,67927,1079195,1874279,236386,289735,1936757,1590714,241621,1559852,1810986,1997592,1026828,1344354,447758,1638093,1970071,354566,407763,1574167,784900,917438,701248,417363,689920,878838,1859514,1057064,147362,1244026,49376,1671560,873954,1314149,108240,875269,351669,1787602,1719403,504249,1373395,574324,876001,1727003,1203512,268626,413130,1916011,164046,1033303,764163,1744547,1745627,1380238,1876226,894450,873509,1488056,109531,2010557,1447844,205753,906312,1900640,671768,877420,385922,793175,1572184,444419,1325007,1755397,719715,454663,1318833,243559,1702837,773333,1088655,1379696,994324,651393,712940,581305,1257454,977961,1710882,1240553,1937277,751558,357545,1199198,439214,530078,914664,1234924,1198197,1459968,1956617,484027,764429,1190894,207563,1228907,257710,832157,501517,1383273,1981084,1878586,1487354,663627,9419,84700,844819,1139713,668458,1398764,27831,1687917,1701193,1125097,579127,1464614,1489787,1884082,850453,568662,513870,589362,167183,2019515,898539,745692,1585326,1222727,258929,1554968,1097580,1048768,1113241,451997,1105911,1863003,1739270,925622,1745081,1452801,291393,446593,1748796,977591,1441608,753149,1968154,1753300,334497,1448622,1510261,532900,1674221,118983,948091,1853347,1769958,1097284,226461,57557,145320,1056339,1213680,124209,859607,847773,261100,1892108,744101,627169,35260,1096739,252380,1347998,177272,1955938,983365,1963452,1280739,486017,453307,1281957,1879327,663735,315548,727544,465452,419433,1503950,1745545,400578,1429508,1365224,901302,477261,1000978,72729,1007684,295238,1782971,831372,775857,1748944,365577,450897,1905821,661333,1934732,280393,316558,567441,1166597,1182564,83471,1009473,41219,2010957,1732813,1892622,1055196,1365790,1729407,821413,629645,1947878,701384,1362462,631599,1534713,1251958,1400528,341032,895666,1945688,1571504,142362,264332,1092362,358571,507520,482101,653918,1158748,1323479,1331599};

#elif PARAM_Q == 3870721

static const uint32_t zetas[PARAM_N] = {2337707,2505409,267692,529914,420735,181988,2608440,3865338,3665767,288746,2524026,3008396,901579,70491,1821213,1437514,3375394,502705,3475623,3513653,1833017,3651222,947790,1966036,2704588,2850143,3030905,1622520,3210245,3127826,292206,3096784,3201921,3867412,1705316,2917474,2975359,2004421,2812268,890313,2511631,3623292,2803099,2903766,1596209,2040136,3468632,2156661,2913824,2560388,1214035,3468039,575792,2926910,3407464,2292204,2285761,2338667,63216,3835938,3204529,1818443,3786633,3241498,944328,616348,2927622,64038,1171534,1361903,2827360,3144828,2738981,1714811,3625146,89505,2787809,2363190,2513795,3306399,1418851,1206903,926563,211044,466372,3410093,1353383,3610570,934100,2471859,2037600,2996463,1698492,525418,1662944,1981925,1210222,1813802,314420,2466015,3516872,3320431,1355971,1500137,493991,36365,3235243,214827,2544017,1739057,945221,1038283,2889903,3364214,1674857,1434035,1665177,2651227,1575769,1155464,467835,1713031,2041544,408424,137443,2029527,2115209,2293884,2137416,3189891,2471629,2229785,2611740,2394735,2287191,2862622,300090,1004990,401830,143957,2910193,3787906,3628164,3171269,2239135,3038465,601725,2887353,2766912,1622354,2989501,1339396,1939160,2386893,103181,2793304,911193,3295333,3025653,2513246,314427,939239,57676,2293294,2833811,2842292,3139575,2705158,1290463,3780876,1462003,668827,1850975,2327221,2910099,2724881,418972,957090,321362,2898276,3523069,1463158,3818473,453440,1891547,1601731,529312,3301251,1117070,3520718,634170,1958581,929634,1133255,3807619,1159272,3292496,3530590,927442,3686531,2605292,384058,1415774,1040397,3663661,2332173,1131260,680774,1186917,3736575,1064994,2954460,3051663,1162037,2962553,2130376,1717870,3565361,2935922,2347272,1768863,3125776,1686747,3137894,2993356,1574419,1073008,1262182,183934,914847,3373156,3688758,2538361,614066,3211143,3565127,1322591,3426188,2951336,172348,3747492,3719872,2962113,778168,2880082,1051508,3741079,1816757,763621,328987,2831790,1276220,135870,3388537,3034187,2419032};

static const uint32_t zetas_inv[PARAM_N] = {1451689,836534,482184,3734851,2594501,1038931,3541734,3107100,2053964,129642,2819213,990639,3092553,908608,150849,123229,3698373,919385,444533,2548130,305594,659578,3256655,1332360,181963,497565,2955874,3686787,2608539,2797713,2296302,877365,732827,2183974,744945,2101858,1523449,934799,305360,2152851,1740345,908168,2708684,819058,916261,2805727,134146,2683804,3189947,2739461,1538548,207060,2830324,2454947,3486663,1265429,184190,2943279,340131,578225,2711449,63102,2737466,2941087,1912140,3236551,350003,2753651,569470,3341409,2268990,1979174,3417281,52248,2407563,347652,972445,3549359,2913631,3451749,1145840,960622,1543500,2019746,3201894,2408718,89845,2580258,1165563,731146,1028429,1036910,1577427,3813045,2931482,3556294,1357475,845068,575388,2959528,1077417,3767540,1483828,1931561,2531325,881220,2248367,1103809,983368,3268996,832256,1631586,699452,242557,82815,960528,3726764,3468891,2865731,3570631,1008099,1583530,1475986,1258981,1640936,1399092,680830,1733305,1576837,1755512,1841194,3733278,3462297,1829177,2157690,3402886,2715257,2294952,1219494,2205544,2436686,2195864,506507,980818,2832438,2925500,2131664,1326704,3655894,635478,3834356,3376730,2370584,2514750,550290,353849,1404706,3556301,2056919,2660499,1888796,2207777,3345303,2172229,874258,1833121,1398862,2936621,260151,2517338,460628,3404349,3659677,2944158,2663818,2451870,564322,1356926,1507531,1082912,3781216,245575,2155910,1131740,725893,1043361,2508818,2699187,3806683,943099,3254373,2926393,629223,84088,2052278,666192,34783,3807505,1532054,1584960,1578517,463257,943811,3294929,402682,2656686,1310333,956897,1714060,402089,1830585,2274512,966955,1067622,247429,1359090,2980408,1058453,1866300,895362,953247,2165405,3309,668800,773937,3578515,742895,660476,2248201,839816,1020578,1166133,1904685,2922931,219499,2037704,357068,395098,3368016,495327,2433207,2049508,3800230,2969142,862325,1346695,3581975,204954,5383,1262281,3688733,3449986,3340807,3603029,951197};
#endif

void aigis_ntt(uint32_t a[PARAM_N])
{

   int i,j;
   __m128i *p128a = (__m128i *) a;
   __m256i *p256a = (__m256i *) a;
   
#ifndef _WIN64
   __m256i vq4x = _mm256_set_epi32(0, PARAM_Q, 0, PARAM_Q, 0, PARAM_Q, 0, PARAM_Q);
   __m256i v2q4x = _mm256_set_epi32(0, 2 * PARAM_Q, 0, 2 * PARAM_Q, 0, 2 * PARAM_Q, 0, 2 * PARAM_Q);
   __m256i vqinv4x = _mm256_set_epi32(0, QINV, 0, QINV, 0, QINV, 0, QINV);
#else
   __m256i vq4x = _mm256_set1_epi64x(PARAM_Q);
   __m256i v2q4x = _mm256_set1_epi64x(2 * PARAM_Q);
   __m256i vqinv4x = _mm256_set1_epi64x(QINV);
#endif
   __m256i r[PARAM_N/4],t[12],pzeta4x[7];
   
#ifndef _WIN64
   pzeta4x[0] = _mm256_set_epi32(0, zetas[1], 0, zetas[1], 0, zetas[1], 0, zetas[1]);
   pzeta4x[1] = _mm256_set_epi32(0, zetas[2], 0, zetas[2], 0, zetas[2], 0, zetas[2]);
   pzeta4x[2] = _mm256_set_epi32(0, zetas[3], 0, zetas[3], 0, zetas[3], 0, zetas[3]);
   pzeta4x[3] = _mm256_set_epi32(0, zetas[4], 0, zetas[4], 0, zetas[4], 0, zetas[4]);
   pzeta4x[4] = _mm256_set_epi32(0, zetas[5], 0, zetas[5], 0, zetas[5], 0, zetas[5]);
   pzeta4x[5] = _mm256_set_epi32(0, zetas[6], 0, zetas[6], 0, zetas[6], 0, zetas[6]);
   pzeta4x[6] = _mm256_set_epi32(0, zetas[7], 0, zetas[7], 0, zetas[7], 0, zetas[7]);
#else
   pzeta4x[0] = _mm256_set1_epi64x(zetas[1]);
   pzeta4x[1] = _mm256_set1_epi64x(zetas[2]);
   pzeta4x[2] = _mm256_set1_epi64x(zetas[3]);
   pzeta4x[3] = _mm256_set1_epi64x(zetas[4]);
   pzeta4x[4] = _mm256_set1_epi64x(zetas[5]);
   pzeta4x[5] = _mm256_set1_epi64x(zetas[6]);
   pzeta4x[6] = _mm256_set1_epi64x(zetas[7]);
#endif

	for(i=0; i< 8; i++)
	{
	   //level 0 
	   t[0] = _mm256_cvtepu32_epi64(p128a[i]);
	   t[1] = _mm256_cvtepu32_epi64(p128a[i+8]);
	   t[2] = _mm256_cvtepu32_epi64(p128a[i+16]);
	   t[3] = _mm256_cvtepu32_epi64(p128a[i+24]);
	   t[4] = _mm256_cvtepu32_epi64(p128a[i+32]);
	   t[5] = _mm256_cvtepu32_epi64(p128a[i+40]);
	   t[6] = _mm256_cvtepu32_epi64(p128a[i+48]);
	   t[7] = _mm256_cvtepu32_epi64(p128a[i+56]);
	   
	   //mul 
	   t[4] = _mm256_mul_epu32(t[4],pzeta4x[0]);
	   t[5] = _mm256_mul_epu32(t[5],pzeta4x[0]);
	   t[6] = _mm256_mul_epu32(t[6],pzeta4x[0]);
	   t[7] = _mm256_mul_epu32(t[7],pzeta4x[0]);
	   
	   //reduce
	   t[8] = _mm256_mul_epu32(t[4],vqinv4x);
	   t[9] = _mm256_mul_epu32(t[5],vqinv4x);
	   t[10] = _mm256_mul_epu32(t[6],vqinv4x);
	   t[11] = _mm256_mul_epu32(t[7],vqinv4x);
	   t[8] = _mm256_mul_epu32(t[8],vq4x);
	   t[9] = _mm256_mul_epu32(t[9],vq4x);
	   t[10] = _mm256_mul_epu32(t[10],vq4x);
	   t[11] = _mm256_mul_epu32(t[11],vq4x);
	   
	   t[8] = _mm256_add_epi64(t[4],t[8]);
	   t[9] = _mm256_add_epi64(t[5],t[9]);
	   t[10] = _mm256_add_epi64(t[6],t[10]);
	   t[11] = _mm256_add_epi64(t[7],t[11]);
	   
	   t[8] = _mm256_srli_epi64(t[8],32);
	   t[9] = _mm256_srli_epi64(t[9],32);
	   t[10] = _mm256_srli_epi64(t[10],32);
	   t[11] = _mm256_srli_epi64(t[11],32);
	   
	   //butterfly
	   t[4] =  _mm256_add_epi64(t[0],v2q4x);
	   t[5] =  _mm256_add_epi64(t[1],v2q4x);
	   t[6] =  _mm256_add_epi64(t[2],v2q4x);
	   t[7] =  _mm256_add_epi64(t[3],v2q4x);
	   
	   t[0] =  _mm256_add_epi64(t[0],t[8]);
	   t[1] =  _mm256_add_epi64(t[1],t[9]);
	   t[2] =  _mm256_add_epi64(t[2],t[10]);
	   t[3] =  _mm256_add_epi64(t[3],t[11]);
	   t[4] =  _mm256_sub_epi64(t[4],t[8]);
	   t[5] =  _mm256_sub_epi64(t[5],t[9]);
	   t[6] =  _mm256_sub_epi64(t[6],t[10]);
	   t[7] =  _mm256_sub_epi64(t[7],t[11]);
	   
	   //level 1
	   //mul 
	   t[2] = _mm256_mul_epu32(t[2],pzeta4x[1]);
	   t[3] = _mm256_mul_epu32(t[3],pzeta4x[1]);
	   t[6] = _mm256_mul_epu32(t[6],pzeta4x[2]);
	   t[7] = _mm256_mul_epu32(t[7],pzeta4x[2]);
	   
	   //reduce
	   t[8] = _mm256_mul_epu32(t[2],vqinv4x);
	   t[9] = _mm256_mul_epu32(t[3],vqinv4x);
	   t[10] = _mm256_mul_epu32(t[6],vqinv4x);
	   t[11] = _mm256_mul_epu32(t[7],vqinv4x);
	   t[8] = _mm256_mul_epu32(t[8],vq4x);
	   t[9] = _mm256_mul_epu32(t[9],vq4x);
	   t[10] = _mm256_mul_epu32(t[10],vq4x);
	   t[11] = _mm256_mul_epu32(t[11],vq4x);
	   
	   t[8] = _mm256_add_epi64(t[2],t[8]);
	   t[9] = _mm256_add_epi64(t[3],t[9]);
	   t[10] = _mm256_add_epi64(t[6],t[10]);
	   t[11] = _mm256_add_epi64(t[7],t[11]);
	   
	   t[8] = _mm256_srli_epi64(t[8],32);
	   t[9] = _mm256_srli_epi64(t[9],32);
	   t[10] = _mm256_srli_epi64(t[10],32);
	   t[11] = _mm256_srli_epi64(t[11],32);
	   
	   //butterfly
	   t[2] =  _mm256_add_epi64(t[0],v2q4x);
	   t[3] =  _mm256_add_epi64(t[1],v2q4x);
	   t[6] =  _mm256_add_epi64(t[4],v2q4x);
	   t[7] =  _mm256_add_epi64(t[5],v2q4x);
	   
	   t[0] =  _mm256_add_epi64(t[0],t[8]);
	   t[1] =  _mm256_add_epi64(t[1],t[9]);
	   t[4] =  _mm256_add_epi64(t[4],t[10]);
	   t[5] =  _mm256_add_epi64(t[5],t[11]);
	   t[2] =  _mm256_sub_epi64(t[2],t[8]);
	   t[3] =  _mm256_sub_epi64(t[3],t[9]);
	   t[6] =  _mm256_sub_epi64(t[6],t[10]);
	   t[7] =  _mm256_sub_epi64(t[7],t[11]);
	   
	   
	   //level 2  
	   //mul 
	   t[1] = _mm256_mul_epu32(t[1],pzeta4x[3]);
	   t[3] = _mm256_mul_epu32(t[3],pzeta4x[4]);
	   t[5] = _mm256_mul_epu32(t[5],pzeta4x[5]);
	   t[7] = _mm256_mul_epu32(t[7],pzeta4x[6]);
	   
	   //reduce
	   t[8] = _mm256_mul_epu32(t[1],vqinv4x);
	   t[9] = _mm256_mul_epu32(t[3],vqinv4x);
	   t[10] = _mm256_mul_epu32(t[5],vqinv4x);
	   t[11] = _mm256_mul_epu32(t[7],vqinv4x);
	   t[8] = _mm256_mul_epu32(t[8],vq4x);
	   t[9] = _mm256_mul_epu32(t[9],vq4x);
	   t[10] = _mm256_mul_epu32(t[10],vq4x);
	   t[11] = _mm256_mul_epu32(t[11],vq4x);
	   
	   t[8] = _mm256_add_epi64(t[1],t[8]);
	   t[9] = _mm256_add_epi64(t[3],t[9]);
	   t[10] = _mm256_add_epi64(t[5],t[10]);
	   t[11] = _mm256_add_epi64(t[7],t[11]);
	   
	   t[8] = _mm256_srli_epi64(t[8],32);
	   t[9] = _mm256_srli_epi64(t[9],32);
	   t[10] = _mm256_srli_epi64(t[10],32);
	   t[11] = _mm256_srli_epi64(t[11],32);
	   
	   //butterfly
	   t[1] =  _mm256_add_epi64(t[0],v2q4x);
	   t[3] =  _mm256_add_epi64(t[2],v2q4x);
	   t[5] =  _mm256_add_epi64(t[4],v2q4x);
	   t[7] =  _mm256_add_epi64(t[6],v2q4x);
	   
	   r[i + 0] =  _mm256_add_epi64(t[0],t[8]);
	   r[i + 16] =  _mm256_add_epi64(t[2],t[9]);
	   r[i + 32] =  _mm256_add_epi64(t[4],t[10]);
	   r[i + 48] =  _mm256_add_epi64(t[6],t[11]);
	   r[i + 8] =  _mm256_sub_epi64(t[1],t[8]);
	   r[i + 24] =  _mm256_sub_epi64(t[3],t[9]);
	   r[i + 40] =  _mm256_sub_epi64(t[5],t[10]);
	   r[i + 56] =  _mm256_sub_epi64(t[7],t[11]);
	 }
	
	j=0;
	for(i=0; i< 64; i+= 8,j++)
	{
	//level 3 
	   t[0] = r[i];
	   t[1] = r[i+1];
	   t[2] = r[i+2];
	   t[3] = r[i+3];
	   
#ifndef _WIN64
	   pzeta4x[0] = _mm256_set_epi32(0, zetas[8 + j], 0, zetas[8 + j], 0, zetas[8 + j], 0, zetas[8 + j]);
#else
	   pzeta4x[0] = _mm256_set1_epi64x(zetas[8+j]); 
#endif
	   
	   //mul 
	   t[4] = _mm256_mul_epu32(r[i+4],pzeta4x[0]);
	   t[5] = _mm256_mul_epu32(r[i+5],pzeta4x[0]);
	   t[6] = _mm256_mul_epu32(r[i+6],pzeta4x[0]);
	   t[7] = _mm256_mul_epu32(r[i+7],pzeta4x[0]);
	   
	   //reduce
	   t[8] = _mm256_mul_epu32(t[4],vqinv4x);
	   t[9] = _mm256_mul_epu32(t[5],vqinv4x);
	   t[10] = _mm256_mul_epu32(t[6],vqinv4x);
	   t[11] = _mm256_mul_epu32(t[7],vqinv4x);
	   t[8] = _mm256_mul_epu32(t[8],vq4x);
	   t[9] = _mm256_mul_epu32(t[9],vq4x);
	   t[10] = _mm256_mul_epu32(t[10],vq4x);
	   t[11] = _mm256_mul_epu32(t[11],vq4x);
	   
	   t[8] = _mm256_add_epi64(t[4],t[8]);
	   t[9] = _mm256_add_epi64(t[5],t[9]);
	   t[10] = _mm256_add_epi64(t[6],t[10]);
	   t[11] = _mm256_add_epi64(t[7],t[11]);
	   
	   t[8] = _mm256_srli_epi64(t[8],32);
	   t[9] = _mm256_srli_epi64(t[9],32);
	   t[10] = _mm256_srli_epi64(t[10],32);
	   t[11] = _mm256_srli_epi64(t[11],32);
	   
	   //butterfly
	   t[4] =  _mm256_add_epi64(r[i+0],v2q4x);
	   t[5] =  _mm256_add_epi64(r[i+1],v2q4x);
	   t[6] =  _mm256_add_epi64(r[i+2],v2q4x);
	   t[7] =  _mm256_add_epi64(r[i+3],v2q4x);
	   
	   t[0] =  _mm256_add_epi64(r[i+0],t[8]);
	   t[1] =  _mm256_add_epi64(r[i+1],t[9]);
	   t[2] =  _mm256_add_epi64(r[i+2],t[10]);
	   t[3] =  _mm256_add_epi64(r[i+3],t[11]);
	   t[4] =  _mm256_sub_epi64(t[4],t[8]);
	   t[5] =  _mm256_sub_epi64(t[5],t[9]);
	   t[6] =  _mm256_sub_epi64(t[6],t[10]);
	   t[7] =  _mm256_sub_epi64(t[7],t[11]);
	   
	   //level 4
#ifndef _WIN64
	   pzeta4x[0] = _mm256_set_epi32(0, zetas[16 + 2 * j], 0, zetas[16 + 2 * j], 0, zetas[16 + 2 * j], 0, zetas[16 + 2 * j]);
	   pzeta4x[1] = _mm256_set_epi32(0, zetas[17 + 2 * j], 0, zetas[17 + 2 * j], 0, zetas[17 + 2 * j], 0, zetas[17 + 2 * j]);
#else
	   pzeta4x[0] = _mm256_set1_epi64x(zetas[16+2*j]);
	   pzeta4x[1] = _mm256_set1_epi64x(zetas[17 + 2 * j]);
#endif
	   //mul 
	   t[2] = _mm256_mul_epu32(t[2],pzeta4x[0]);
	   t[3] = _mm256_mul_epu32(t[3],pzeta4x[0]);
	   t[6] = _mm256_mul_epu32(t[6],pzeta4x[1]);
	   t[7] = _mm256_mul_epu32(t[7],pzeta4x[1]);
	   
	   //reduce
	   t[8] = _mm256_mul_epu32(t[2],vqinv4x);
	   t[9] = _mm256_mul_epu32(t[3],vqinv4x);
	   t[10] = _mm256_mul_epu32(t[6],vqinv4x);
	   t[11] = _mm256_mul_epu32(t[7],vqinv4x);
	   t[8] = _mm256_mul_epu32(t[8],vq4x);
	   t[9] = _mm256_mul_epu32(t[9],vq4x);
	   t[10] = _mm256_mul_epu32(t[10],vq4x);
	   t[11] = _mm256_mul_epu32(t[11],vq4x);
	   
	   t[8] = _mm256_add_epi64(t[2],t[8]);
	   t[9] = _mm256_add_epi64(t[3],t[9]);
	   t[10] = _mm256_add_epi64(t[6],t[10]);
	   t[11] = _mm256_add_epi64(t[7],t[11]);
	   
	   t[8] = _mm256_srli_epi64(t[8],32);
	   t[9] = _mm256_srli_epi64(t[9],32);
	   t[10] = _mm256_srli_epi64(t[10],32);
	   t[11] = _mm256_srli_epi64(t[11],32);
	   
	   //butterfly
	   t[2] =  _mm256_add_epi64(t[0],v2q4x);
	   t[3] =  _mm256_add_epi64(t[1],v2q4x);
	   t[6] =  _mm256_add_epi64(t[4],v2q4x);
	   t[7] =  _mm256_add_epi64(t[5],v2q4x);
	   
	   t[0] =  _mm256_add_epi64(t[0],t[8]);
	   t[1] =  _mm256_add_epi64(t[1],t[9]);
	   t[4] =  _mm256_add_epi64(t[4],t[10]);
	   t[5] =  _mm256_add_epi64(t[5],t[11]);
	   t[2] =  _mm256_sub_epi64(t[2],t[8]);
	   t[3] =  _mm256_sub_epi64(t[3],t[9]);
	   t[6] =  _mm256_sub_epi64(t[6],t[10]);
	   t[7] =  _mm256_sub_epi64(t[7],t[11]);
	   
	   
	   //level 5
#ifndef _WIN64
	   pzeta4x[0] = _mm256_set_epi32(0, zetas[32 + 4 * j], 0, zetas[32 + 4 * j], 0, zetas[32 + 4 * j], 0, zetas[32 + 4 * j]);
	   pzeta4x[1] = _mm256_set_epi32(0, zetas[33 + 4 * j], 0, zetas[33 + 4 * j], 0, zetas[33 + 4 * j], 0, zetas[33 + 4 * j]);
	   pzeta4x[2] = _mm256_set_epi32(0, zetas[34 + 4 * j], 0, zetas[34 + 4 * j], 0, zetas[34 + 4 * j], 0, zetas[34 + 4 * j]);
	   pzeta4x[3] = _mm256_set_epi32(0, zetas[35 + 4 * j], 0, zetas[35 + 4 * j], 0, zetas[35 + 4 * j], 0, zetas[35 + 4 * j]);
#else
	   pzeta4x[0] = _mm256_set1_epi64x(zetas[32 + 4 * j]);
	   pzeta4x[1] = _mm256_set1_epi64x(zetas[33 + 4 * j]);
	   pzeta4x[2] = _mm256_set1_epi64x(zetas[34 + 4 * j]);
	   pzeta4x[3] = _mm256_set1_epi64x(zetas[35 + 4 * j]);
#endif
	   
	   //mul 
	   t[1] = _mm256_mul_epu32(t[1],pzeta4x[0]);
	   t[3] = _mm256_mul_epu32(t[3],pzeta4x[1]);
	   t[5] = _mm256_mul_epu32(t[5],pzeta4x[2]);
	   t[7] = _mm256_mul_epu32(t[7],pzeta4x[3]);
	   
	   //reduce
	   t[8] = _mm256_mul_epu32(t[1],vqinv4x);
	   t[9] = _mm256_mul_epu32(t[3],vqinv4x);
	   t[10] = _mm256_mul_epu32(t[5],vqinv4x);
	   t[11] = _mm256_mul_epu32(t[7],vqinv4x);
	   t[8] = _mm256_mul_epu32(t[8],vq4x);
	   t[9] = _mm256_mul_epu32(t[9],vq4x);
	   t[10] = _mm256_mul_epu32(t[10],vq4x);
	   t[11] = _mm256_mul_epu32(t[11],vq4x);
	   
	   t[8] = _mm256_add_epi64(t[1],t[8]);
	   t[9] = _mm256_add_epi64(t[3],t[9]);
	   t[10] = _mm256_add_epi64(t[5],t[10]);
	   t[11] = _mm256_add_epi64(t[7],t[11]);
	   
	   t[8] = _mm256_srli_epi64(t[8],32);
	   t[9] = _mm256_srli_epi64(t[9],32);
	   t[10] = _mm256_srli_epi64(t[10],32);
	   t[11] = _mm256_srli_epi64(t[11],32);
	   
	   //butterfly
	   t[1] =  _mm256_add_epi64(t[0],v2q4x);
	   t[3] =  _mm256_add_epi64(t[2],v2q4x);
	   t[5] =  _mm256_add_epi64(t[4],v2q4x);
	   t[7] =  _mm256_add_epi64(t[6],v2q4x);

	   t[1] =  _mm256_sub_epi64(t[1],t[8]);
	   t[3] =  _mm256_sub_epi64(t[3],t[9]);
	   t[5] =  _mm256_sub_epi64(t[5],t[10]);
	   t[7] =  _mm256_sub_epi64(t[7],t[11]);
	   
	   t[8] =  _mm256_add_epi64(t[0],t[8]);
	   t[9] =  _mm256_add_epi64(t[2],t[9]);
	   t[10] =  _mm256_add_epi64(t[4],t[10]);
	   t[11] =  _mm256_add_epi64(t[6],t[11]);
	   
	   
	  //level 6
	  t[0] = _mm256_permute2x128_si256(t[8],t[1],0x20);
       t[1] = _mm256_permute2x128_si256(t[8],t[1],0x31); 
       t[2] = _mm256_permute2x128_si256(t[9],t[3],0x20);
       t[3] = _mm256_permute2x128_si256(t[9],t[3],0x31);
       t[4] = _mm256_permute2x128_si256(t[10],t[5],0x20);
       t[5] = _mm256_permute2x128_si256(t[10],t[5],0x31);  
       t[6] = _mm256_permute2x128_si256(t[11],t[7],0x20);
       t[7] = _mm256_permute2x128_si256(t[11],t[7],0x31);
       
#ifndef _WIN64
	   pzeta4x[0] = _mm256_set_epi32(0, zetas[65 + 8 * j], 0, zetas[65 + 8 * j], 0, zetas[64 + 8 * j], 0, zetas[64 + 8 * j]);
	   pzeta4x[1] = _mm256_set_epi32(0, zetas[67 + 8 * j], 0, zetas[67 + 8 * j], 0, zetas[66 + 8 * j], 0, zetas[66 + 8 * j]);
	   pzeta4x[2] = _mm256_set_epi32(0, zetas[69 + 8 * j], 0, zetas[69 + 8 * j], 0, zetas[68 + 8 * j], 0, zetas[68 + 8 * j]);
	   pzeta4x[3] = _mm256_set_epi32(0, zetas[71 + 8 * j], 0, zetas[71 + 8 * j], 0, zetas[70 + 8 * j], 0, zetas[70 + 8 * j]);
#else
	   pzeta4x[0] = _mm256_set_epi64x(zetas[65 + 8 * j], zetas[65 + 8 * j], zetas[64 + 8 * j], zetas[64 + 8 * j]);
	   pzeta4x[1] = _mm256_set_epi64x(zetas[67 + 8 * j], zetas[67 + 8 * j], zetas[66 + 8 * j], zetas[66 + 8 * j]);
	   pzeta4x[2] = _mm256_set_epi64x(zetas[69 + 8 * j], zetas[69 + 8 * j], zetas[68 + 8 * j], zetas[68 + 8 * j]);
	   pzeta4x[3] = _mm256_set_epi64x(zetas[71 + 8 * j], zetas[71 + 8 * j], zetas[70 + 8 * j], zetas[70 + 8 * j]);
#endif
       
       //mul
       t[1] = _mm256_mul_epu32(t[1],pzeta4x[0]);
       t[3] = _mm256_mul_epu32(t[3],pzeta4x[1]);
       t[5] = _mm256_mul_epu32(t[5],pzeta4x[2]);
       t[7] = _mm256_mul_epu32(t[7],pzeta4x[3]);
		  		   
	  //reduce
	   t[8] = _mm256_mul_epu32(t[1],vqinv4x);
	   t[9] = _mm256_mul_epu32(t[3],vqinv4x);
	   t[10] = _mm256_mul_epu32(t[5],vqinv4x);
	   t[11] = _mm256_mul_epu32(t[7],vqinv4x);
	   t[8] = _mm256_mul_epu32(t[8],vq4x);
	   t[9] = _mm256_mul_epu32(t[9],vq4x);
	   t[10] = _mm256_mul_epu32(t[10],vq4x);
	   t[11] = _mm256_mul_epu32(t[11],vq4x);
	   
	   t[8] = _mm256_add_epi64(t[1],t[8]);
	   t[9] = _mm256_add_epi64(t[3],t[9]);
	   t[10] = _mm256_add_epi64(t[5],t[10]);
	   t[11] = _mm256_add_epi64(t[7],t[11]);
	   
	   t[8] = _mm256_srli_epi64(t[8],32);
	   t[9] = _mm256_srli_epi64(t[9],32);
	   t[10] = _mm256_srli_epi64(t[10],32);
	   t[11] = _mm256_srli_epi64(t[11],32);
		  
       //butterfly
	   t[1] =  _mm256_add_epi64(t[0],v2q4x);
	   t[3] =  _mm256_add_epi64(t[2],v2q4x);
	   t[5] =  _mm256_add_epi64(t[4],v2q4x);
	   t[7] =  _mm256_add_epi64(t[6],v2q4x);
	   
	   t[0] =  _mm256_add_epi64(t[0],t[8]);
	   t[2] =  _mm256_add_epi64(t[2],t[9]);
	   t[4] =  _mm256_add_epi64(t[4],t[10]);
	   t[6] =  _mm256_add_epi64(t[6],t[11]);
	   t[1] =  _mm256_sub_epi64(t[1],t[8]);
	   t[3] =  _mm256_sub_epi64(t[3],t[9]);
	   t[5] =  _mm256_sub_epi64(t[5],t[10]);
	   t[7] =  _mm256_sub_epi64(t[7],t[11]);
	  
	  //level 7
	  t[8] = _mm256_permute4x64_epi64(t[1],0xb1);
	  t[9] = _mm256_permute4x64_epi64(t[3],0xb1);
	  t[10] = _mm256_permute4x64_epi64(t[5],0xb1);
	  t[11] = _mm256_permute4x64_epi64(t[7],0xb1);
	  t[1] = _mm256_blend_epi32(t[8],t[0],0xcc);
	  t[3] = _mm256_blend_epi32(t[9],t[2],0xcc);
	  t[5] = _mm256_blend_epi32(t[10],t[4],0xcc);
	  t[7] = _mm256_blend_epi32(t[11],t[6],0xcc);
	  
	  t[1] = _mm256_permute4x64_epi64(t[1],0xb1);
	  t[3] = _mm256_permute4x64_epi64(t[3],0xb1);
	  t[5] = _mm256_permute4x64_epi64(t[5],0xb1);
	  t[7] = _mm256_permute4x64_epi64(t[7],0xb1);
	  t[0] = _mm256_blend_epi32(t[0],t[8],0xcc);
	  t[2] = _mm256_blend_epi32(t[2],t[9],0xcc);
	  t[4] = _mm256_blend_epi32(t[4],t[10],0xcc);
	  t[6] = _mm256_blend_epi32(t[6],t[11],0xcc);
	  
#ifndef _WIN64	  
	  pzeta4x[0] = _mm256_set_epi32(0, zetas[131 + 16 * j], 0, zetas[130 + 16 * j], 0, zetas[129 + 16 * j], 0, zetas[128 + 16 * j]);
	  pzeta4x[1] = _mm256_set_epi32(0, zetas[135 + 16 * j], 0, zetas[134 + 16 * j], 0, zetas[133 + 16 * j], 0, zetas[132 + 16 * j]);
	  pzeta4x[2] = _mm256_set_epi32(0, zetas[139 + 16 * j], 0, zetas[138 + 16 * j], 0, zetas[137 + 16 * j], 0, zetas[136 + 16 * j]);
	  pzeta4x[3] = _mm256_set_epi32(0, zetas[143 + 16 * j], 0, zetas[142 + 16 * j], 0, zetas[141 + 16 * j], 0, zetas[140 + 16 * j]);
#else
	  pzeta4x[0] = _mm256_set_epi64x(zetas[131 + 16*j],zetas[130 + 16*j],zetas[129 + 16*j],zetas[128 + 16*j]);
	  pzeta4x[1] = _mm256_set_epi64x(zetas[135 + 16 * j], zetas[134 + 16 * j], zetas[133 + 16 * j], zetas[132 + 16 * j]);
	  pzeta4x[2] = _mm256_set_epi64x(zetas[139 + 16 * j], zetas[138 + 16 * j], zetas[137 + 16 * j], zetas[136 + 16 * j]);
	  pzeta4x[3] = _mm256_set_epi64x(zetas[143 + 16 * j], zetas[142 + 16 * j], zetas[141 + 16 * j], zetas[140 + 16 * j]);
#endif
	  //mul
       t[1] = _mm256_mul_epu32(t[1],pzeta4x[0]);
       t[3] = _mm256_mul_epu32(t[3],pzeta4x[1]);
       t[5] = _mm256_mul_epu32(t[5],pzeta4x[2]);
       t[7] = _mm256_mul_epu32(t[7],pzeta4x[3]);
		  		   
	  //reduce
	   t[8] = _mm256_mul_epu32(t[1],vqinv4x);
	   t[9] = _mm256_mul_epu32(t[3],vqinv4x);
	   t[10] = _mm256_mul_epu32(t[5],vqinv4x);
	   t[11] = _mm256_mul_epu32(t[7],vqinv4x);
	   t[8] = _mm256_mul_epu32(t[8],vq4x);
	   t[9] = _mm256_mul_epu32(t[9],vq4x);
	   t[10] = _mm256_mul_epu32(t[10],vq4x);
	   t[11] = _mm256_mul_epu32(t[11],vq4x);
	   
	   t[8] = _mm256_add_epi64(t[1],t[8]);
	   t[9] = _mm256_add_epi64(t[3],t[9]);
	   t[10] = _mm256_add_epi64(t[5],t[10]);
	   t[11] = _mm256_add_epi64(t[7],t[11]);
	   
	   t[8] = _mm256_srli_epi64(t[8],32);
	   t[9] = _mm256_srli_epi64(t[9],32);
	   t[10] = _mm256_srli_epi64(t[10],32);
	   t[11] = _mm256_srli_epi64(t[11],32);
		  
       //butterfly
	   t[1] =  _mm256_add_epi64(t[0],v2q4x);
	   t[3] =  _mm256_add_epi64(t[2],v2q4x);
	   t[5] =  _mm256_add_epi64(t[4],v2q4x);
	   t[7] =  _mm256_add_epi64(t[6],v2q4x);
	   
	   t[0] =  _mm256_add_epi64(t[0],t[8]);
	   t[2] =  _mm256_add_epi64(t[2],t[9]);
	   t[4] =  _mm256_add_epi64(t[4],t[10]);
	   t[6] =  _mm256_add_epi64(t[6],t[11]);
	   t[1] =  _mm256_sub_epi64(t[1],t[8]);
	   t[3] =  _mm256_sub_epi64(t[3],t[9]);
	   t[5] =  _mm256_sub_epi64(t[5],t[10]);
	   t[7] =  _mm256_sub_epi64(t[7],t[11]);
	   
	  
	  //store
	  t[1] = _mm256_slli_epi64(t[1],32);
	  t[3] = _mm256_slli_epi64(t[3],32); 
	  t[5] = _mm256_slli_epi64(t[5],32); 
	  t[7] = _mm256_slli_epi64(t[7],32);
	  
	  p256a[4*j] = _mm256_or_si256(t[0],t[1]);
	  p256a[4*j+1] = _mm256_or_si256(t[2],t[3]);
	  p256a[4*j+2] = _mm256_or_si256(t[4],t[5]);
	  p256a[4*j+3] = _mm256_or_si256(t[6],t[7]);
	  
      }
}
void aigis_invntt(uint32_t a[PARAM_N])
{

	int i,j;
   __m128i *p128a = (__m128i *) a;
   __m256i *p256a = (__m256i *) a;
   
#ifndef _WIN64
   __m256i vq4x = _mm256_set_epi32(0, PARAM_Q, 0, PARAM_Q, 0, PARAM_Q, 0, PARAM_Q);
   __m256i v512q4x = _mm256_set_epi32(0, 512 * PARAM_Q, 0, 512 * PARAM_Q, 0, 512 * PARAM_Q, 0, 512 * PARAM_Q);
   __m256i vqinv4x = _mm256_set_epi32(0, QINV, 0, QINV, 0, QINV, 0, QINV);
   __m256i vf4x = _mm256_set_epi32(0, mont_invn, 0, mont_invn, 0, mont_invn, 0, mont_invn);
#else
   __m256i vq4x = _mm256_set1_epi64x(PARAM_Q);
   __m256i v512q4x = _mm256_set1_epi64x(512 * PARAM_Q);
   __m256i vqinv4x = _mm256_set1_epi64x(QINV);
   __m256i vf4x = _mm256_set1_epi64x(mont_invn);
#endif
   __m256i vidx = _mm256_set_epi32(0,0,0,0,7,5,3,1);
   
   __m256i r[PARAM_N/4];
   __m256i t[12],pzeta4x[7];
  
   j=0; 
    for(i = 0; i < PARAM_N/4; i += 8,j++) 
   {   
   
   //level 0
#ifndef _WIN64
		pzeta4x[0] = _mm256_set_epi32(0, zetas_inv[16 * j + 3], 0, zetas_inv[16 * j + 2], 0, zetas_inv[16 * j + 1], 0, zetas_inv[16 * j]);
		pzeta4x[1] = _mm256_set_epi32(0, zetas_inv[16 * j + 7], 0, zetas_inv[16 * j + 6], 0, zetas_inv[16 * j + 5], 0, zetas_inv[16 * j + 4]);
		pzeta4x[2] = _mm256_set_epi32(0, zetas_inv[16 * j + 11], 0, zetas_inv[16 * j + 10], 0, zetas_inv[16 * j + 9], 0, zetas_inv[16 * j + 8]);
		pzeta4x[3] = _mm256_set_epi32(0, zetas_inv[16 * j + 15], 0, zetas_inv[16 * j + 14], 0, zetas_inv[16 * j + 13], 0, zetas_inv[16 * j + 12]);
#else
       pzeta4x[0] = _mm256_set_epi64x(zetas_inv[16*j+3],zetas_inv[16*j+2],zetas_inv[16*j+1],zetas_inv[16*j]);
       pzeta4x[1] = _mm256_set_epi64x(zetas_inv[16*j+7],zetas_inv[16*j+6],zetas_inv[16*j+5],zetas_inv[16*j+4]);
       pzeta4x[2] = _mm256_set_epi64x(zetas_inv[16*j+11],zetas_inv[16*j+10],zetas_inv[16*j+9],zetas_inv[16*j+8]);
       pzeta4x[3] = _mm256_set_epi64x(zetas_inv[16*j+15],zetas_inv[16*j+14],zetas_inv[16*j+13],zetas_inv[16*j+12]);
#endif 
       t[8] = p256a[4*j];
       t[9] = p256a[4*j+1];
       t[10] = p256a[4*j+2];
       t[11] = p256a[4*j+3];
       t[1] = _mm256_srli_epi64(t[8],32);
       t[3] = _mm256_srli_epi64(t[9],32);
       t[5] = _mm256_srli_epi64(t[10],32);
       t[7] = _mm256_srli_epi64(t[11],32);
        
       //butterfly
       
       t[0] = _mm256_add_epi64(t[8],t[1]);
       t[2] = _mm256_add_epi64(t[9],t[3]);
       t[4] = _mm256_add_epi64(t[10],t[5]);
       t[6] = _mm256_add_epi64(t[11],t[7]);
       
       t[8] = _mm256_add_epi64(t[8],v512q4x);
       t[9] = _mm256_add_epi64(t[9],v512q4x);
       t[10] = _mm256_add_epi64(t[10],v512q4x);
       t[11] = _mm256_add_epi64(t[11],v512q4x);
       
       t[1] = _mm256_sub_epi64(t[8],t[1]);
       t[3] = _mm256_sub_epi64(t[9],t[3]);
       t[5] = _mm256_sub_epi64(t[10],t[5]);
       t[7] = _mm256_sub_epi64(t[11],t[7]);
  
	  //mul 
	  t[1] = _mm256_mul_epu32(t[1],pzeta4x[0]);
	  t[3] = _mm256_mul_epu32(t[3],pzeta4x[1]);
	  t[5] = _mm256_mul_epu32(t[5],pzeta4x[2]);
	  t[7] = _mm256_mul_epu32(t[7],pzeta4x[3]);
	  
    		   
	  //reduce
	  t[8] = _mm256_mul_epu32(t[1],vqinv4x);
	  t[9] = _mm256_mul_epu32(t[3],vqinv4x);
	  t[10] = _mm256_mul_epu32(t[5],vqinv4x);
	  t[11] = _mm256_mul_epu32(t[7],vqinv4x);
	  
	  
	  t[8] = _mm256_mul_epu32(t[8],vq4x);
	  t[9] = _mm256_mul_epu32(t[9],vq4x);
	  t[10] = _mm256_mul_epu32(t[10],vq4x);
	  t[11] = _mm256_mul_epu32(t[11],vq4x);
	  
	  t[1] = _mm256_add_epi64(t[1],t[8]);
	  t[3] = _mm256_add_epi64(t[3],t[9]);
	  t[5] = _mm256_add_epi64(t[5],t[10]);
	  t[7] = _mm256_add_epi64(t[7],t[11]);
	  
	  t[1] = _mm256_srli_epi64(t[1],32);
	  t[3] = _mm256_srli_epi64(t[3],32);
	  t[5] = _mm256_srli_epi64(t[5],32);
	  t[7] = _mm256_srli_epi64(t[7],32);
     
       //level 1
#ifndef _WIN64
       pzeta4x[0] = _mm256_set_epi32(0,zetas_inv[129 + 8*j],0,zetas_inv[129+8*j],0,zetas_inv[128+8*j],0,zetas_inv[128+8*j]);
       pzeta4x[1] = _mm256_set_epi32(0,zetas_inv[131 + 8*j],0,zetas_inv[131+8*j],0,zetas_inv[130+8*j],0,zetas_inv[130+8*j]);
       pzeta4x[2] = _mm256_set_epi32(0,zetas_inv[133 + 8*j],0,zetas_inv[133+8*j],0,zetas_inv[132+8*j],0,zetas_inv[132+8*j]);
       pzeta4x[3] = _mm256_set_epi32(0,zetas_inv[135 + 8*j],0,zetas_inv[135+8*j],0,zetas_inv[134+8*j],0,zetas_inv[134+8*j]);
#else
	  pzeta4x[0] = _mm256_set_epi64x(zetas_inv[129 + 8*j],zetas_inv[129+8*j],zetas_inv[128+8*j],zetas_inv[128+8*j]);
	  pzeta4x[1] = _mm256_set_epi64x(zetas_inv[131 + 8 * j], zetas_inv[131 + 8 * j], zetas_inv[130 + 8 * j], zetas_inv[130 + 8 * j]);
	  pzeta4x[2] = _mm256_set_epi64x(zetas_inv[133 + 8 * j], zetas_inv[133 + 8 * j], zetas_inv[132 + 8 * j], zetas_inv[132 + 8 * j]);
	  pzeta4x[3] = _mm256_set_epi64x(zetas_inv[135 + 8 * j], zetas_inv[135 + 8 * j], zetas_inv[134 + 8 * j], zetas_inv[134 + 8 * j]);
#endif
       t[8] = _mm256_permute4x64_epi64(t[1],0xb1);
       t[9] = _mm256_permute4x64_epi64(t[3],0xb1);
       t[10] = _mm256_permute4x64_epi64(t[5],0xb1);
       t[11] = _mm256_permute4x64_epi64(t[7],0xb1);
       
       t[1] = _mm256_blend_epi32(t[8],t[0],0xcc);
       t[3] = _mm256_blend_epi32(t[9],t[2],0xcc);
       t[5] = _mm256_blend_epi32(t[10],t[4],0xcc);
       t[7] = _mm256_blend_epi32(t[11],t[6],0xcc);
       
       t[1] = _mm256_permute4x64_epi64(t[1],0xb1);
       t[3] = _mm256_permute4x64_epi64(t[3],0xb1);
       t[5] = _mm256_permute4x64_epi64(t[5],0xb1);
       t[7] = _mm256_permute4x64_epi64(t[7],0xb1);
       
       t[8] = _mm256_blend_epi32(t[0],t[8],0xcc);
       t[9] = _mm256_blend_epi32(t[2],t[9],0xcc);
       t[10] = _mm256_blend_epi32(t[4],t[10],0xcc);
       t[11] = _mm256_blend_epi32(t[6],t[11],0xcc);
  
       //butterfly
       t[0] = _mm256_add_epi64(t[8],t[1]);
       t[2] = _mm256_add_epi64(t[9],t[3]);
       t[4] = _mm256_add_epi64(t[10],t[5]);
       t[6] = _mm256_add_epi64(t[11],t[7]);
       
       t[8] =  _mm256_add_epi64(t[8],v512q4x);
       t[9] =  _mm256_add_epi64(t[9],v512q4x);
       t[10] =  _mm256_add_epi64(t[10],v512q4x);
       t[11] =  _mm256_add_epi64(t[11],v512q4x);
       
       
       t[1] =  _mm256_sub_epi64(t[8],t[1]);
       t[3] =  _mm256_sub_epi64(t[9],t[3]);
       t[5] =  _mm256_sub_epi64(t[10],t[5]);
       t[7] =  _mm256_sub_epi64(t[11],t[7]);
       
       
	  
	  //mul
	  t[1] = _mm256_mul_epu32(t[1],pzeta4x[0]);
	  t[3] = _mm256_mul_epu32(t[3],pzeta4x[1]); 
	  t[5] = _mm256_mul_epu32(t[5],pzeta4x[2]); 
	  t[7] = _mm256_mul_epu32(t[7],pzeta4x[3]); 
	   
	    		   
	  //reduce
	  t[8] = _mm256_mul_epu32(t[1],vqinv4x);
	  t[9] = _mm256_mul_epu32(t[3],vqinv4x);
	  t[10] = _mm256_mul_epu32(t[5],vqinv4x);
	  t[11] = _mm256_mul_epu32(t[7],vqinv4x);
	  
	  t[8] = _mm256_mul_epu32(t[8],vq4x);
	  t[9] = _mm256_mul_epu32(t[9],vq4x);
	  t[10] = _mm256_mul_epu32(t[10],vq4x);
	  t[11] = _mm256_mul_epu32(t[11],vq4x);
	  
	  t[1] = _mm256_add_epi64(t[1],t[8]);
	  t[3] = _mm256_add_epi64(t[3],t[9]);
	  t[5] = _mm256_add_epi64(t[5],t[10]);
	  t[7] = _mm256_add_epi64(t[7],t[11]);
	  
	  t[1] = _mm256_srli_epi64(t[1],32);
	  t[3] = _mm256_srli_epi64(t[3],32);
	  t[5] = _mm256_srli_epi64(t[5],32);
	  t[7] = _mm256_srli_epi64(t[7],32);
	  
	  
	  t[8] = _mm256_permute2x128_si256(t[0],t[1],0x20);
	  t[1] = _mm256_permute2x128_si256(t[0],t[1],0x31);
	  t[9] = _mm256_permute2x128_si256(t[2],t[3],0x20);
	  t[3] = _mm256_permute2x128_si256(t[2],t[3],0x31);
	  t[10] = _mm256_permute2x128_si256(t[4],t[5],0x20);
	  t[5] = _mm256_permute2x128_si256(t[4],t[5],0x31);
	  t[11] = _mm256_permute2x128_si256(t[6],t[7],0x20);
	  t[7] = _mm256_permute2x128_si256(t[6],t[7],0x31);
	  
	  
	  //level 2
#ifndef _WIN64
	  pzeta4x[0] = _mm256_set_epi32(0, zetas_inv[192 + 4 * j], 0, zetas_inv[192 + 4 * j], 0, zetas_inv[192 + 4 * j], 0, zetas_inv[192 + 4 * j]);
	  pzeta4x[1] = _mm256_set_epi32(0, zetas_inv[193 + 4 * j], 0, zetas_inv[193 + 4 * j], 0, zetas_inv[193 + 4 * j], 0, zetas_inv[193 + 4 * j]);
	  pzeta4x[2] = _mm256_set_epi32(0, zetas_inv[194 + 4 * j], 0, zetas_inv[194 + 4 * j], 0, zetas_inv[194 + 4 * j], 0, zetas_inv[194 + 4 * j]);
	  pzeta4x[3] = _mm256_set_epi32(0, zetas_inv[195 + 4 * j], 0, zetas_inv[195 + 4 * j], 0, zetas_inv[195 + 4 * j], 0, zetas_inv[195 + 4 * j]);
#else
	  pzeta4x[0] = _mm256_set1_epi64x(zetas_inv[192+4*j]);
	  pzeta4x[1] = _mm256_set1_epi64x(zetas_inv[193 + 4 * j]);
	  pzeta4x[2] = _mm256_set1_epi64x(zetas_inv[194 + 4 * j]);
	  pzeta4x[3] = _mm256_set1_epi64x(zetas_inv[195 + 4 * j]);
#endif
	  
	  //butterfly
       t[0] = _mm256_add_epi64(t[8],t[1]);
       t[2] = _mm256_add_epi64(t[9],t[3]);
       t[4] = _mm256_add_epi64(t[10],t[5]);
       t[6] = _mm256_add_epi64(t[11],t[7]);
       
       t[8] =  _mm256_add_epi64(t[8],v512q4x);
       t[9] =  _mm256_add_epi64(t[9],v512q4x);
       t[10] =  _mm256_add_epi64(t[10],v512q4x);
       t[11] =  _mm256_add_epi64(t[11],v512q4x);
       
       
       t[1] =  _mm256_sub_epi64(t[8],t[1]);
       t[3] =  _mm256_sub_epi64(t[9],t[3]);
       t[5] =  _mm256_sub_epi64(t[10],t[5]);
       t[7] =  _mm256_sub_epi64(t[11],t[7]);
       
       //mul
	  t[1] = _mm256_mul_epu32(t[1],pzeta4x[0]);
	  t[3] = _mm256_mul_epu32(t[3],pzeta4x[1]); 
	  t[5] = _mm256_mul_epu32(t[5],pzeta4x[2]); 
	  t[7] = _mm256_mul_epu32(t[7],pzeta4x[3]); 
	   
	    		   
	  //reduce
	  t[8] = _mm256_mul_epu32(t[1],vqinv4x);
	  t[9] = _mm256_mul_epu32(t[3],vqinv4x);
	  t[10] = _mm256_mul_epu32(t[5],vqinv4x);
	  t[11] = _mm256_mul_epu32(t[7],vqinv4x);
	  
	  t[8] = _mm256_mul_epu32(t[8],vq4x);
	  t[9] = _mm256_mul_epu32(t[9],vq4x);
	  t[10] = _mm256_mul_epu32(t[10],vq4x);
	  t[11] = _mm256_mul_epu32(t[11],vq4x);
	  
	  t[1] = _mm256_add_epi64(t[1],t[8]);
	  t[3] = _mm256_add_epi64(t[3],t[9]);
	  t[5] = _mm256_add_epi64(t[5],t[10]);
	  t[7] = _mm256_add_epi64(t[7],t[11]);
	  
	  t[1] = _mm256_srli_epi64(t[1],32);
	  t[3] = _mm256_srli_epi64(t[3],32);
	  t[5] = _mm256_srli_epi64(t[5],32);
	  t[7] = _mm256_srli_epi64(t[7],32);
	  
	  
	  //level 3
#ifndef _WIN64
	  pzeta4x[0] = _mm256_set_epi32(0, zetas_inv[224 + 2 * j], 0, zetas_inv[224 + 2 * j], 0, zetas_inv[224 + 2 * j], 0, zetas_inv[224 + 2 * j]);
	  pzeta4x[1] = _mm256_set_epi32(0, zetas_inv[225 + 2 * j], 0, zetas_inv[225 + 2 * j], 0, zetas_inv[225 + 2 * j], 0, zetas_inv[225 + 2 * j]);
#else
	  pzeta4x[0] = _mm256_set1_epi64x(zetas_inv[224+2*j]);
	  pzeta4x[1] = _mm256_set1_epi64x(zetas_inv[225 + 2 * j]);
#endif
       
       //butterfly
       t[8] =  _mm256_add_epi64(t[0],v512q4x);
       t[9] =  _mm256_add_epi64(t[1],v512q4x);
       t[10] =  _mm256_add_epi64(t[4],v512q4x);
       t[11] =  _mm256_add_epi64(t[5],v512q4x);
       
       t[0] = _mm256_add_epi64(t[0],t[2]);
       t[1] = _mm256_add_epi64(t[1],t[3]);
       t[4] = _mm256_add_epi64(t[4],t[6]);
       t[5] = _mm256_add_epi64(t[5],t[7]);
       
       
       t[2] =  _mm256_sub_epi64(t[8],t[2]);
       t[3] =  _mm256_sub_epi64(t[9],t[3]);
       t[6] =  _mm256_sub_epi64(t[10],t[6]);
       t[7] =  _mm256_sub_epi64(t[11],t[7]);
       
       //mul
	  t[2] = _mm256_mul_epu32(t[2],pzeta4x[0]);
	  t[3] = _mm256_mul_epu32(t[3],pzeta4x[0]); 
	  t[6] = _mm256_mul_epu32(t[6],pzeta4x[1]); 
	  t[7] = _mm256_mul_epu32(t[7],pzeta4x[1]); 
	   
	    		   
	  //reduce
	  t[8] = _mm256_mul_epu32(t[2],vqinv4x);
	  t[9] = _mm256_mul_epu32(t[3],vqinv4x);
	  t[10] = _mm256_mul_epu32(t[6],vqinv4x);
	  t[11] = _mm256_mul_epu32(t[7],vqinv4x);
	  
	  t[8] = _mm256_mul_epu32(t[8],vq4x);
	  t[9] = _mm256_mul_epu32(t[9],vq4x);
	  t[10] = _mm256_mul_epu32(t[10],vq4x);
	  t[11] = _mm256_mul_epu32(t[11],vq4x);
	  
	  t[2] = _mm256_add_epi64(t[2],t[8]);
	  t[3] = _mm256_add_epi64(t[3],t[9]);
	  t[6] = _mm256_add_epi64(t[6],t[10]);
	  t[7] = _mm256_add_epi64(t[7],t[11]);
	  
	  t[2] = _mm256_srli_epi64(t[2],32);
	  t[3] = _mm256_srli_epi64(t[3],32);
	  t[6] = _mm256_srli_epi64(t[6],32);
	  t[7] = _mm256_srli_epi64(t[7],32);
	  
	  
	  //level 4
#ifndef _WIN64
	  pzeta4x[0] = _mm256_set_epi32(0, zetas_inv[240 + j], 0, zetas_inv[240 + j], 0, zetas_inv[240 + j], 0, zetas_inv[240 + j]);
#else
	  pzeta4x[0] = _mm256_set1_epi64x(zetas_inv[240+j]);
#endif
       
       //butterfly
       t[8] =  _mm256_add_epi64(t[0],v512q4x);
       t[9] =  _mm256_add_epi64(t[1],v512q4x);
       t[10] =  _mm256_add_epi64(t[2],v512q4x);
       t[11] =  _mm256_add_epi64(t[3],v512q4x);
       
       r[i+0] = _mm256_add_epi64(t[0],t[4]);
       r[i+1] = _mm256_add_epi64(t[1],t[5]);
       r[i+2] = _mm256_add_epi64(t[2],t[6]);
       r[i+3] = _mm256_add_epi64(t[3],t[7]);
       
       
       t[4] =  _mm256_sub_epi64(t[8],t[4]);
       t[5] =  _mm256_sub_epi64(t[9],t[5]);
       t[6] =  _mm256_sub_epi64(t[10],t[6]);
       t[7] =  _mm256_sub_epi64(t[11],t[7]);
       
       //mul
	  t[4] = _mm256_mul_epu32(t[4],pzeta4x[0]);
	  t[5] = _mm256_mul_epu32(t[5],pzeta4x[0]); 
	  t[6] = _mm256_mul_epu32(t[6],pzeta4x[0]); 
	  t[7] = _mm256_mul_epu32(t[7],pzeta4x[0]); 
	   
	    		   
	  //reduce
	  t[8] = _mm256_mul_epu32(t[4],vqinv4x);
	  t[9] = _mm256_mul_epu32(t[5],vqinv4x);
	  t[10] = _mm256_mul_epu32(t[6],vqinv4x);
	  t[11] = _mm256_mul_epu32(t[7],vqinv4x);
	  
	  t[8] = _mm256_mul_epu32(t[8],vq4x);
	  t[9] = _mm256_mul_epu32(t[9],vq4x);
	  t[10] = _mm256_mul_epu32(t[10],vq4x);
	  t[11] = _mm256_mul_epu32(t[11],vq4x);
	  
	  t[4] = _mm256_add_epi64(t[4],t[8]);
	  t[5] = _mm256_add_epi64(t[5],t[9]);
	  t[6] = _mm256_add_epi64(t[6],t[10]);
	  t[7] = _mm256_add_epi64(t[7],t[11]);
	  
	  r[i+4] = _mm256_srli_epi64(t[4],32);
	  r[i+5] = _mm256_srli_epi64(t[5],32);
	  r[i+6] = _mm256_srli_epi64(t[6],32);
	  r[i+7] = _mm256_srli_epi64(t[7],32); 
      }
    
    
#ifndef _WIN64
	  pzeta4x[0] = _mm256_set_epi32(0, zetas_inv[248], 0, zetas_inv[248], 0, zetas_inv[248], 0, zetas_inv[248]);
	  pzeta4x[1] = _mm256_set_epi32(0, zetas_inv[249], 0, zetas_inv[249], 0, zetas_inv[249], 0, zetas_inv[249]);
	  pzeta4x[2] = _mm256_set_epi32(0, zetas_inv[250], 0, zetas_inv[250], 0, zetas_inv[250], 0, zetas_inv[250]);
	  pzeta4x[3] = _mm256_set_epi32(0, zetas_inv[251], 0, zetas_inv[251], 0, zetas_inv[251], 0, zetas_inv[251]);
	  pzeta4x[4] = _mm256_set_epi32(0, zetas_inv[252], 0, zetas_inv[252], 0, zetas_inv[252], 0, zetas_inv[252]);
	  pzeta4x[5] = _mm256_set_epi32(0, zetas_inv[253], 0, zetas_inv[253], 0, zetas_inv[253], 0, zetas_inv[253]);
	  pzeta4x[6] = _mm256_set_epi32(0, zetas_inv[254], 0, zetas_inv[254], 0, zetas_inv[254], 0, zetas_inv[254]);
#else
	  pzeta4x[0] = _mm256_set1_epi64x(zetas_inv[248]);
	  pzeta4x[1] = _mm256_set1_epi64x(zetas_inv[249]);
	  pzeta4x[2] = _mm256_set1_epi64x(zetas_inv[250]);
	  pzeta4x[3] = _mm256_set1_epi64x(zetas_inv[251]);
	  pzeta4x[4] = _mm256_set1_epi64x(zetas_inv[252]);
	  pzeta4x[5] = _mm256_set1_epi64x(zetas_inv[253]);
	  pzeta4x[6] = _mm256_set1_epi64x(zetas_inv[254]);
#endif

    for(i=0; i< 8; i++)
	{
	   //level 5   
	  //butterfly
       t[8] =  _mm256_add_epi64(r[i],v512q4x);
       t[9] =  _mm256_add_epi64(r[i+16],v512q4x);
       t[10] =  _mm256_add_epi64(r[i+32],v512q4x);
       t[11] =  _mm256_add_epi64(r[i+48],v512q4x);
       
       t[0] = _mm256_add_epi64(r[i],   r[i+8]);
       t[2] = _mm256_add_epi64(r[i+16],r[i+24]);
       t[4] = _mm256_add_epi64(r[i+32],r[i+40]);
       t[6] = _mm256_add_epi64(r[i+48],r[i+56]);
       
       
       t[1] =  _mm256_sub_epi64(t[8],r[i+8]);
       t[3] =  _mm256_sub_epi64(t[9],r[i+24]);
       t[5] =  _mm256_sub_epi64(t[10],r[i+40]);
       t[7] =  _mm256_sub_epi64(t[11],r[i+56]);
       
       //mul
	  t[1] = _mm256_mul_epu32(t[1],pzeta4x[0]);
	  t[3] = _mm256_mul_epu32(t[3],pzeta4x[1]); 
	  t[5] = _mm256_mul_epu32(t[5],pzeta4x[2]); 
	  t[7] = _mm256_mul_epu32(t[7],pzeta4x[3]); 
	   
	    		   
	  //reduce
	  t[8] = _mm256_mul_epu32(t[1],vqinv4x);
	  t[9] = _mm256_mul_epu32(t[3],vqinv4x);
	  t[10] = _mm256_mul_epu32(t[5],vqinv4x);
	  t[11] = _mm256_mul_epu32(t[7],vqinv4x);
	  
	  t[8] = _mm256_mul_epu32(t[8],vq4x);
	  t[9] = _mm256_mul_epu32(t[9],vq4x);
	  t[10] = _mm256_mul_epu32(t[10],vq4x);
	  t[11] = _mm256_mul_epu32(t[11],vq4x);
	  
	  t[1] = _mm256_add_epi64(t[1],t[8]);
	  t[3] = _mm256_add_epi64(t[3],t[9]);
	  t[5] = _mm256_add_epi64(t[5],t[10]);
	  t[7] = _mm256_add_epi64(t[7],t[11]);
	  
	  t[1] = _mm256_srli_epi64(t[1],32);
	  t[3] = _mm256_srli_epi64(t[3],32);
	  t[5] = _mm256_srli_epi64(t[5],32);
	  t[7] = _mm256_srli_epi64(t[7],32);
	  
	  
	  //level 6   
	  //butterfly
       t[8] =  _mm256_add_epi64(t[0],v512q4x);
       t[9] =  _mm256_add_epi64(t[1],v512q4x);
       t[10] =  _mm256_add_epi64(t[4],v512q4x);
       t[11] =  _mm256_add_epi64(t[5],v512q4x);
       
       t[0] = _mm256_add_epi64(t[0],t[2]);
       t[1] = _mm256_add_epi64(t[1],t[3]);
       t[4] = _mm256_add_epi64(t[4],t[6]);
       t[5] = _mm256_add_epi64(t[5],t[7]);
       
       
       t[2] =  _mm256_sub_epi64(t[8],t[2]);
       t[3] =  _mm256_sub_epi64(t[9],t[3]);
       t[6] =  _mm256_sub_epi64(t[10],t[6]);
       t[7] =  _mm256_sub_epi64(t[11],t[7]);
       
       //mul
	  t[2] = _mm256_mul_epu32(t[2],pzeta4x[4]);
	  t[3] = _mm256_mul_epu32(t[3],pzeta4x[4]); 
	  t[6] = _mm256_mul_epu32(t[6],pzeta4x[5]); 
	  t[7] = _mm256_mul_epu32(t[7],pzeta4x[5]); 
	   
	    		   
	  //reduce
	  t[8] = _mm256_mul_epu32(t[2],vqinv4x);
	  t[9] = _mm256_mul_epu32(t[3],vqinv4x);
	  t[10] = _mm256_mul_epu32(t[6],vqinv4x);
	  t[11] = _mm256_mul_epu32(t[7],vqinv4x);
	  
	  t[8] = _mm256_mul_epu32(t[8],vq4x);
	  t[9] = _mm256_mul_epu32(t[9],vq4x);
	  t[10] = _mm256_mul_epu32(t[10],vq4x);
	  t[11] = _mm256_mul_epu32(t[11],vq4x);
	  
	  t[2] = _mm256_add_epi64(t[2],t[8]);
	  t[3] = _mm256_add_epi64(t[3],t[9]);
	  t[6] = _mm256_add_epi64(t[6],t[10]);
	  t[7] = _mm256_add_epi64(t[7],t[11]);
	  
	  t[2] = _mm256_srli_epi64(t[2],32);
	  t[3] = _mm256_srli_epi64(t[3],32);
	  t[6] = _mm256_srli_epi64(t[6],32);
	  t[7] = _mm256_srli_epi64(t[7],32);
	  
	  
	  //level 7   
	  //butterfly
       t[8] =  _mm256_add_epi64(t[0],v512q4x);
       t[9] =  _mm256_add_epi64(t[1],v512q4x);
       t[10] =  _mm256_add_epi64(t[2],v512q4x);
       t[11] =  _mm256_add_epi64(t[3],v512q4x);
       
       t[0] = _mm256_add_epi64(t[0],t[4]);
       t[1] = _mm256_add_epi64(t[1],t[5]);
       t[2] = _mm256_add_epi64(t[2],t[6]);
       t[3] = _mm256_add_epi64(t[3],t[7]);
       
       
       t[4] =  _mm256_sub_epi64(t[8],t[4]);
       t[5] =  _mm256_sub_epi64(t[9],t[5]);
       t[6] =  _mm256_sub_epi64(t[10],t[6]);
       t[7] =  _mm256_sub_epi64(t[11],t[7]);
       
       //mul the first half by _f
       t[0] = _mm256_mul_epu32(t[0],vf4x);
	  t[1] = _mm256_mul_epu32(t[1],vf4x); 
	  t[2] = _mm256_mul_epu32(t[2],vf4x); 
	  t[3] = _mm256_mul_epu32(t[3],vf4x);
	   
	  //reduce
	  t[8] = _mm256_mul_epu32(t[0],vqinv4x);
	  t[9] = _mm256_mul_epu32(t[1],vqinv4x);
	  t[10] = _mm256_mul_epu32(t[2],vqinv4x);
	  t[11] = _mm256_mul_epu32(t[3],vqinv4x);
	  
	  t[8] = _mm256_mul_epu32(t[8],vq4x);
	  t[9] = _mm256_mul_epu32(t[9],vq4x);
	  t[10] = _mm256_mul_epu32(t[10],vq4x);
	  t[11] = _mm256_mul_epu32(t[11],vq4x);
	  
	  t[0] = _mm256_add_epi64(t[0],t[8]);
	  t[1] = _mm256_add_epi64(t[1],t[9]);
	  t[2] = _mm256_add_epi64(t[2],t[10]);
	  t[3] = _mm256_add_epi64(t[3],t[11]);
	  
	  /*t[4] = _mm256_srli_epi64(t[4],32);
	  t[5] = _mm256_srli_epi64(t[5],32);
	  t[6] = _mm256_srli_epi64(t[6],32);
	  t[7] = _mm256_srli_epi64(t[7],32); */
       
       
       //mul
	  t[4] = _mm256_mul_epu32(t[4],pzeta4x[6]);
	  t[5] = _mm256_mul_epu32(t[5],pzeta4x[6]); 
	  t[6] = _mm256_mul_epu32(t[6],pzeta4x[6]); 
	  t[7] = _mm256_mul_epu32(t[7],pzeta4x[6]); 
	   
	    		   
	  //reduce
	  t[8] = _mm256_mul_epu32(t[4],vqinv4x);
	  t[9] = _mm256_mul_epu32(t[5],vqinv4x);
	  t[10] = _mm256_mul_epu32(t[6],vqinv4x);
	  t[11] = _mm256_mul_epu32(t[7],vqinv4x);
	  
	  t[8] = _mm256_mul_epu32(t[8],vq4x);
	  t[9] = _mm256_mul_epu32(t[9],vq4x);
	  t[10] = _mm256_mul_epu32(t[10],vq4x);
	  t[11] = _mm256_mul_epu32(t[11],vq4x);
	  
	  t[4] = _mm256_add_epi64(t[4],t[8]);
	  t[5] = _mm256_add_epi64(t[5],t[9]);
	  t[6] = _mm256_add_epi64(t[6],t[10]);
	  t[7] = _mm256_add_epi64(t[7],t[11]);
	  
	  /*t[4] = _mm256_srli_epi64(t[4],32);
	  t[5] = _mm256_srli_epi64(t[5],32);
	  t[6] = _mm256_srli_epi64(t[6],32);
	  t[7] = _mm256_srli_epi64(t[7],32);*/
	  
	  
	  t[0] = _mm256_permutevar8x32_epi32(t[0],vidx);
	  t[1] = _mm256_permutevar8x32_epi32(t[1],vidx);
	  t[2] = _mm256_permutevar8x32_epi32(t[2],vidx);
	  t[3] = _mm256_permutevar8x32_epi32(t[3],vidx);
	  t[4] = _mm256_permutevar8x32_epi32(t[4],vidx);
	  t[5] = _mm256_permutevar8x32_epi32(t[5],vidx);
	  t[6] = _mm256_permutevar8x32_epi32(t[6],vidx);
	  t[7] = _mm256_permutevar8x32_epi32(t[7],vidx);
	  
	  p128a[i] = _mm256_extracti128_si256(t[0],0);
	  p128a[i+8] = _mm256_extracti128_si256(t[1],0);
	  p128a[i+16] = _mm256_extracti128_si256(t[2],0);
	  p128a[i+24] = _mm256_extracti128_si256(t[3],0);
	  p128a[i+32] = _mm256_extracti128_si256(t[4],0);
	  p128a[i+40] = _mm256_extracti128_si256(t[5],0);
	  p128a[i+48] = _mm256_extracti128_si256(t[6],0);
	  p128a[i+56] = _mm256_extracti128_si256(t[7],0);
	}
}
