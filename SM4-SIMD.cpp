#include <iostream>
#include <immintrin.h>
#include<stdio.h>
#include "SM4-SIMD.h"

#define u8 unsigned char
#define u32 unsigned int

u32 loopLeft(u32 a, int length) {
	return ((a << length) | a >> (32 - length));
}

u32 functionT(u32 b) {
	u8 a[4];
	short i;
	a[0] = b / 0x1000000;
	a[1] = b / 0x10000;
	a[2] = b / 0x100;
	a[3] = b;
	b = Sbox[a[0]] * 0x1000000 + Sbox[a[1]] * 0x10000 + Sbox[a[2]] * 0x100 + Sbox[a[3]];
	b = b ^ loopLeft(b, 13) ^ loopLeft(b, 23);
	return b;
}

//密钥扩展算法
void getRK(u32 MK[], u32 K[], u32 RK[]) {
	int i;
	for (i = 0; i < 4; i++) {
		K[i] = MK[i] ^ FK[i];
	}
	for (i = 0; i < 32; i++) {
		K[(i + 4) % 4] = K[i % 4] ^ functionT(K[(i + 1) % 4] ^ K[(i + 2) % 4] ^ K[(i + 3) % 4] ^ CK[i]);
		RK[i] = K[(i + 4) % 4];
	}
}

void SM4(u8* m, u8* c, u32* RK, int mode) {
    __m256i X[4], R[4];
    __m256i Mask;
    Mask = _mm256_set1_epi32(0xFF);
    R[0] = _mm256_loadu_si256((const __m256i*)m + 0);
    R[1] = _mm256_loadu_si256((const __m256i*)m + 1);
    R[2] = _mm256_loadu_si256((const __m256i*)m + 2);
    R[3] = _mm256_loadu_si256((const __m256i*)m + 3);
    X[0] = _mm256_unpacklo_epi64(_mm256_unpacklo_epi32(R[0], R[1]), _mm256_unpacklo_epi32(R[2], R[3]));
    X[1] = _mm256_unpackhi_epi64(_mm256_unpacklo_epi32(R[0], R[1]), _mm256_unpacklo_epi32(R[2], R[3]));
    X[2] = _mm256_unpacklo_epi64(_mm256_unpackhi_epi32(R[0], R[1]), _mm256_unpackhi_epi32(R[2], R[3]));
    X[3] = _mm256_unpackhi_epi64(_mm256_unpackhi_epi32(R[0], R[1]), _mm256_unpackhi_epi32(R[2], R[3]));
    __m256i vindex = _mm256_setr_epi8(3, 2, 1, 0, 7, 6, 5, 4, 11, 10, 9, 8, 15, 14, 13, 12, 3, 2, 1, 0, 7, 6, 5, 4, 11, 10, 9, 8, 15, 14, 13, 12);
    X[0] = _mm256_shuffle_epi8(X[0], vindex);
    X[1] = _mm256_shuffle_epi8(X[1], vindex);
    X[2] = _mm256_shuffle_epi8(X[2], vindex);
    X[3] = _mm256_shuffle_epi8(X[3], vindex);

    for (int i = 0; i < 32; i++) {
        __m256i k =_mm256_set1_epi32((mode == 0) ? RK[i] : RK[31 - i]);//加密or解密
        R[0] = _mm256_xor_si256(_mm256_xor_si256(X[1], X[2]),_mm256_xor_si256(X[3], k));
        R[1] = _mm256_xor_si256(X[0], _mm256_i32gather_epi32((const int*)BOX0,_mm256_and_si256(R[0], Mask), 4));
        R[0] = _mm256_srli_epi32(R[0], 8);
        R[1] = _mm256_xor_si256(R[1], _mm256_i32gather_epi32((const int*)BOX1, _mm256_and_si256(R[0], Mask), 4));
        R[0] = _mm256_srli_epi32(R[0], 8);
        R[1] = _mm256_xor_si256(R[1], _mm256_i32gather_epi32((const int*)BOX2, _mm256_and_si256(R[0], Mask), 4));
        R[0] = _mm256_srli_epi32(R[0], 8);
        R[1] = _mm256_xor_si256(R[1], _mm256_i32gather_epi32((const int*)BOX3, _mm256_and_si256(R[0], Mask), 4));
        X[0] = X[1];
        X[1] = X[2];
        X[2] = X[3];
        X[3] = R[1];
    }

    X[0] = _mm256_shuffle_epi8(X[0], vindex);
    X[1] = _mm256_shuffle_epi8(X[1], vindex);
    X[2] = _mm256_shuffle_epi8(X[2], vindex);
    X[3] = _mm256_shuffle_epi8(X[3], vindex);
    _mm256_storeu_si256((__m256i*)c + 0, _mm256_unpacklo_epi64(_mm256_unpacklo_epi32(X[3], X[2]), _mm256_unpacklo_epi32(X[1], X[0])));
    _mm256_storeu_si256((__m256i*)c + 1, _mm256_unpackhi_epi64(_mm256_unpacklo_epi32(X[3], X[2]), _mm256_unpacklo_epi32(X[1], X[0])));
    _mm256_storeu_si256((__m256i*)c + 2, _mm256_unpacklo_epi64(_mm256_unpackhi_epi32(X[3], X[2]), _mm256_unpackhi_epi32(X[1], X[0])));
    _mm256_storeu_si256((__m256i*)c + 3, _mm256_unpackhi_epi64(_mm256_unpackhi_epi32(X[3], X[2]), _mm256_unpackhi_epi32(X[1], X[0])));
}

int main()
{
	u32 X[4]; // 明文 
    u32 MK[4] = { 0x01234567,0x89abcdef,0xfedcba98,0x76543210 }; // 密钥 
	u32 RK[32]; // 轮密钥  
	u32 K[4]; // 中间数据 
	u32 Y[4]; // 密文 
	short i; // 临时变量 
    u8 in[16 * 8] = { 0x01, 0x23, 0x45, 0x67, 0x89, 0xab, 0xcd, 0xef, 0xfe, 0xdc, 0xba, 0x98, 0x76, 0x54, 0x32, 0x10 };
	printf("**************生成轮密钥*****************\n");
	getRK(MK, K, RK);
	for (i = 0; i < 32; i++) {
		printf("[%2d]：%08x    ", i, RK[i]);
		if (i % 4 == 3)	printf("\n");
	}
    SM4(in, in, RK, 0);
    printf("**************生成加密结果*****************\n");
    for (int j = 0; j < 8; j++) {
        for (int i = 0; i < 16; i++) {
            printf("%02x ", in[i + 16 * j]);
        }
        printf("\n");
    }
}

