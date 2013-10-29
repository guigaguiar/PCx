double dnrm2(int *n, double *x, int *incx);
double ddot(int *n, double *dx, int *incx, double *dy, int *incy); // arquivo cblas
void   sortd(double *va, int *v, int n);
void   both(MMTtype *M, int *Acol, double *v1, double *v2, double *x, double *v, double *y, double *w);
void   cholesky(int aux6, int pp2, int p, double s1, double v, double valgaux[], double valg[], double valg1[], double wwaux[], double l1[], double rr1[], double rr2[], double rr4[], double *rr3, double u[], double dd12[]);
void   resolucao_sistema(int aux6, int p, double s1, double sti[], double ssf[], double ssf2[], double dd12[], double valg1[], double rr4[]);
