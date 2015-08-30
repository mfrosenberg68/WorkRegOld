/* "Rauschfunktionen" */
extern double aint(double x);
extern double randlc(double& x, const double a);
extern double Xi1double(double x, int N = 100);
extern double Xi2double(double x, int N = 100);
extern double Xi3double(double x, int N = 100);
extern mp_real Xi1(mp_real x, int N = 100);
extern mp_real Xi2(mp_real x, int N = 100);
extern mp_real Xi3(mp_real x, int N = 100);

/* Funktionen der Modelldaten */
extern mp_real mpe2(mp_real x, mp_real xi);
extern mp_real mpe2exakt(mp_real x);

/* reine regul. Lösungen; keine Datenfehler */
extern mp_real mpRegexakt(mp_real x_, mp_real b_);
