#include "pch.h"
#include "matrixCalculate.h"
#include <iostream>

void mAdd(int n1, int n2, int n3, double *m1, double *m2, double *result)
{
	if (n1 != n2 && n2 != n3) {
		std::cout << "Matrix Add Error" << std::endl;
		return;
	}

	for (int i = 0; i < n1 * n2; ++i)
		*result = *(m1++) + *(m2++);
}

/********¼ÆËã3¡Á3¾ØÕóµÄÄæ¾ØÕó*********/
void inv33(double *m)
{
	int i;
	double rr;
	double qa[9];
	rr = m[0] * m[4] * m[8] + m[2] * m[3] * m[7] + m[1] * m[5] * m[6]
		- m[2] * m[4] * m[6] - m[1] * m[3] * m[8] - m[0] * m[5] * m[7];
	qa[0] = (m[4] * m[8] - m[5] * m[7]) / rr;
	qa[1] = -(m[1] * m[8] - m[2] * m[7]) / rr;
	qa[2] = (m[1] * m[5] - m[2] * m[4]) / rr;
	qa[3] = -(m[3] * m[8] - m[5] * m[6]) / rr;
	qa[4] = (m[0] * m[8] - m[2] * m[6]) / rr;
	qa[5] = -(m[0] * m[5] - m[2] * m[3]) / rr;
	qa[6] = (m[3] * m[7] - m[4] * m[6]) / rr;
	qa[7] = -(m[0] * m[7] - m[1] * m[6]) / rr;
	qa[8] = (m[0] * m[4] - m[1] * m[3]) / rr;
	for (i = 0; i < 9; i++)
		m[i] = qa[i];
}

/********¼ÆËãn¡Án¾ØÕóµÄÄæ¾ØÕó*********/
void dcinv(double a[], int n)
{
	int *is, *js, i, j, k, l, u, v;
	double d, p;

	is = (int *)malloc(n * sizeof(int));
	js = (int *)malloc(n * sizeof(int));

	for (k = 0; k <= n - 1; k++)
	{
		d = 0.0;
		for (i = k; i <= n - 1; i++)
			for (j = k; j <= n - 1; j++)
			{
				l = i * n + j;
				p = fabs(a[l]);
				if (p > d)
				{
					d = p;
					is[k] = i;
					js[k] = j;
				}
			}
		if (d + 1.0 == 1.0)
		{
			free(is);
			free(js);
			printf("err* *not inv\n");
		}
		if (is[k] != k)
			for (j = 0; j <= n - 1; j++)
			{
				u = k * n + j;
				v = is[k] * n + j;
				p = a[u];
				a[u] = a[v];
				a[v] = p;
			}
		if (js[k] != k)
			for (i = 0; i <= n - 1; i++)
			{
				u = i * n + k;
				v = i * n + js[k];
				p = a[u];
				a[u] = a[v];
				a[v] = p;
			}
		l = k * n + k;
		a[l] = 1.0 / a[l];
		for (j = 0; j <= n - 1; j++)
			if (j != k)
			{
				u = k * n + j;
				a[u] = a[u] * a[l];
			}
		for (i = 0; i <= n - 1; i++)
			if (i != k)
				for (j = 0; j <= n - 1; j++)
					if (j != k)
					{
						u = i * n + j;
						a[u] = a[u] - a[i*n + k] * a[k*n + j];
					}
		for (i = 0; i <= n - 1; i++)
			if (i != k)
			{
				u = i * n + k;
				a[u] = -a[u] * a[l];
			}
	}

	for (k = n - 1; k >= 0; k--)
	{
		if (js[k] != k)
			for (j = 0; j <= n - 1; j++)
			{
				u = k * n + j;
				v = js[k] * n + j;
				p = a[u];
				a[u] = a[v];
				a[v] = p;
			}
		if (is[k] != k)
			for (i = 0; i <= n - 1; i++)
			{
				u = i * n + k;
				v = i * n + is[k];
				p = a[u];
				a[u] = a[v];
				a[v] = p;
			}
	}

	free(is);
	free(js);
	
}

/**************¾ØÕó³Ë·¨º¯Êý(result=m1*m2')**************/
void mXtrm(int n1, int n2, int n3, double *m1, double *m2, double *result)
{
	int i, j, k;
	for (i = 0; i < n1; i++)
		for (j = 0; j < n3; j++)
		{
			result[i*n3 + j] = 0.0;
			for (k = 0; k < n2; k++)
				result[i*n3 + j] = result[i*n3 + j] + m1[i*n2 + k] * m2[j*n2 + k];
		}
}

/**************¾ØÕó³Ë·¨º¯Êý(c=base+a*b)**************/
void mAddMult(int n, int m, int l, double *base, double *a, double *b, double *c)
{
	int i, j, k;
	for (i = 0; i < n; i++)
		for (j = 0; j < l; j++)
		{
			c[i*l + j] = base[i*l + j];
			for (k = 0; k < m; k++)
				c[i*l + j] += a[i*m + k] * b[l*k + j];
		}
}

/**************¾ØÕó³Ë·¨º¯Êý(result=m1*m2)**************/
void mXm(int n1, int n2, int n3, double *m1, double *m2, double *result)
{
	int i, j, k;
	for (i = 0; i < n1; i++)
		for (j = 0; j < n3; j++)
		{
			result[i*n3 + j] = 0.0;
			for (k = 0; k < n2; k++)
				result[i*n3 + j] = result[i*n3 + j] + m1[i*n2 + k] * m2[k*n3 + j];
		}
}


/********°×ÔëÉù·¢ÉúÆ÷*********/
void white(double Qk[])
{
	static double s0 = 65536.0, w0 = 2053.0, v0 = 13849.0, r0 = 0.0, t0 = 0, EGx, m0;
	double G[] = { 10, 10, 15, 0.2, 0.2, 0.2, 0.1, 0.1, 0.1, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01 };


	for (int i = 0; i < 12; i++)
	{
		r0 = w0 * r0 + v0;
		m0 = r0 / s0;
		r0 = r0 - m0 * s0;
		t0 = t0 + r0 / s0;
	}
	EGx = t0 - 6.0;

	mXtrm(15, 1, 15, G, G, Qk);
	for (int i = 0; i < 225; ++i)
		Qk[i] *= EGx;
}

void printMat(double A[], int m, int n)
{
	for (int i = 0; i < m; ++i) {
		for (int j = 0; j < n; ++j)
			std::cout << A[i*n + j] << " ";
		std::cout << std::endl;
	}
}
