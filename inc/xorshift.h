/*  Header file for Sebastiano Vigna's xorshift1024*φ
    I'm keeping things as close to the original as possible*/

#ifndef XORSHIFT_H
#define XORSHIFT_H

#include <stdint.h>

extern uint64_t s[16];
extern int p;

uint64_t next(void);
void jump(void);
double nextU01();
void initxorshift(int n);

#endif /* XORSHIFT_H */