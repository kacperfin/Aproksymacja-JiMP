#include "makespl.h"
#include "piv_ge_solver.h"

#include <stdio.h>
#include <stdlib.h>
#include <float.h>

/* UWAGA: liczbę używanych f. bazowych można ustawić przez wartość
          zmiennej środowiskowej APPROX_BASE_SIZE
*/

/*
 * Funkcje bazowe: n - liczba funkcji a,b - granice przedzialu aproksymacji i
 * - numer funkcji x - wspolrzedna dla ktorej obliczana jest wartosc funkcji
 */

//Funkcja obl. wartosc wielomianu Czebyszewa
double fT(int k, double x)
{
	if(k==0) return 1;
	if(k==1) return x;
	
	return 2*x*fT(k-1, x) - fT(k-2, x);
}

//Pierwsza pochodna
double dfT1(int k, double x)
{
	if(k==0) return 0;
	if(k==1) return 1;

	return 2*x*dfT1(k-1, x) + 2*fT(k-1, x) - dfT1(k-2, x);
}

//Druga pochodna
double dfT2(int k, double x)
{
	if(k==0) return 0;
	if(k==1) return 0;

	return 4*dfT1(k-1, x) + 2*x*dfT2(k-1, x) - dfT2(k-2, x);
}

double dfT3(int k, double x)
{
	if(k==0) return 0;
	if(k==1) return 0;

	return 6*dfT2(k-1, x) + 2*x*dfT3(k-1, x) - dfT3(k-2, x);
}

void
make_spl(points_t * pts, spline_t * spl)
{

	matrix_t       *eqs= NULL;
	double         *x = pts->x;
	double         *y = pts->y;
	double		a = x[0];
	double		b = x[pts->n - 1];
	int		i, j, k;
	int		nb = pts->n - 3 > 10 ? 10 : pts->n - 3;
  char *nbEnv= getenv( "APPROX_BASE_SIZE" );

	if( nbEnv != NULL && atoi( nbEnv ) > 0 )
		nb = atoi( nbEnv );

	eqs = make_matrix(nb, nb + 1);

	for (j = 0; j < nb; j++) {
		for (i = 0; i < nb; i++)
			for (k = 0; k < pts->n; k++)
				add_to_entry_matrix(eqs, j, i, fT(i, x[k]) * fT(j, x[k]));

		for (k = 0; k < pts->n; k++)
			add_to_entry_matrix(eqs, j, nb, y[k] * fT(j, x[k]));
	}

	if (piv_ge_solver(eqs)) {
		spl->n = 0;
		return;
	}

	if (alloc_spl(spl, nb) == 0) {
		for (i = 0; i < spl->n; i++) {
			double xx = spl->x[i] = a + i*(b-a)/(spl->n-1);
			xx+= 10.0*DBL_EPSILON;  // zabezpieczenie przed ulokowaniem punktu w poprzednim przedziale
			spl->f[i] = 0;
			spl->f1[i] = 0;
			spl->f2[i] = 0;
			spl->f3[i] = 0;
			for (k = 0; k < nb; k++) {
				double		ck = get_entry_matrix(eqs, k, nb);
				spl->f[i]  += ck * fT  (k, xx);
				spl->f1[i] += ck * dfT1 (k, xx);
				spl->f2[i] += ck * dfT2(k, xx);
				spl->f3[i] += ck * dfT3(k, xx);
			}
		}
	}

}
