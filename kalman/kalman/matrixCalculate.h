#pragma once
#include <cstdlib>
#include <cstdio>

/********¼ÆËã3¡Á3¾ØÕóµÄÄæ¾ØÕó*********/
void inv33(double *m);

/********¼ÆËãn¡Án¾ØÕóµÄÄæ¾ØÕó*********/
void dcinv(double a[], int n);

/**************¾ØÕó³Ë·¨º¯Êı(result=m1*m2')**************/
void mXtrm(int n1, int n2, int n3, double *m1, double *m2, double *result);

/**************¾ØÕó³Ë·¨º¯Êı(c=base+a*b)**************/
void mAddMult(int n, int m, int l, double *base, double *a, double *b, double *c);

/**************¾ØÕó³Ë·¨º¯Êı(result=m1*m2)**************/
void mXm(int n1, int n2, int n3, double *m1, double *m2, double *result);

void mAdd(int n1, int n2, int n3, double *m1, double *m2, double *result);

void white(double Qk[]);
void printMat(double A[], int m, int n);