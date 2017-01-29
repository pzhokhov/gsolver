
#include "extmath.h"


void expandtwice(f_complex* A, int stride, int dist, int Nin, int Nout, int centerpoint)
{
	for (int n=0; n<Nout; n++)
	{
	 f_complex* B = A+n*dist;
	 if (centerpoint == 0) 
	 {
	  for (int i=0; i<Nin/2; i++)    B[i*stride] = B[2*i*stride]; 
	  for (int i=Nin/2; i<Nin; i++)  B[i*stride] = 0; 
	 }
	 if (centerpoint == 1) 
	 {
	  for (int i=0; i<Nin/4; i++)
	  {
		  B[i*stride] = B[2*i*stride];
		  B[(Nin-i)*stride] = B[(Nin-2*i)*stride]; 
	  }
	  for (int i=Nin/4; i<3*Nin/4; i++)  A[i*stride + n*dist] = 0;
	 }
	 if (centerpoint == 2) 
	 {
	  for (int i=0; i<Nin/4; i++)
	  {
		  B[(Nin/2-i)*stride] = B[(Nin/2-2*i)*stride];
		  B[(Nin/2+i)*stride] = B[(Nin/2+2*i)*stride];   
	  }
	  for (int i=0; i<=Nin/4; i++)
	  {
		  B[i*stride] = 0;
		  B[(Nin-i)*stride] = 0;
	  }
	 }
	}

	if (centerpoint*(centerpoint-2) > 0) throw "expandtwice(): Unknown value of centerpoint (should be 0, 1 or 2)";
}

void contracttwice(f_complex* A, int stride, int dist, int Nin, int Nout, int centerpoint)
{
	for (int n=0; n<Nout; n++)
	{
	 f_complex* B = A+n*dist;
	 if (centerpoint == 0) 
	 {
	  for (int i=Nin/2-1; i>=0; i--)    B[2*i*stride]     = B[i*stride]; 
	  for (int i=0; i<Nin/2-1; i++)      B[(2*i+1)*stride] = (B[2*i*stride] + B[(2*i+2)*stride])/(float_type)2.0;
	  B[Nin-1]=0;
	 }

	 if (centerpoint == 1) 
	 {
	  for (int i=Nin/4-1; i>=0; i--)
	  {
		  B[2*i*stride] = B[i*stride];
		  B[(Nin-1-2*i)*stride] = B[(Nin-1-i)*stride]; 
	  }
	  for (int i=0; i<Nin/2-1; i++)      B[(2*i+1)*stride] = (B[2*i*stride] + B[(2*i+2)*stride])/(float_type)2.0;

	  for (int i=Nin/4; i<3*Nin/4; i++)  A[i*stride + n*dist] = 0;
	 }
	 if (centerpoint == 2) 
	 {
	  for (int i=Nin/4-1; i>=0; i--)
	  {
		 B[(Nin/2-2*i)*stride] = B[(Nin/2-i)*stride];
		 B[(Nin/2+2*i)*stride] = B[(Nin/2+i)*stride]; 
	  }
	  for (int i=0; i<Nin/2-1; i++)      B[(2*i+1)*stride] = (B[(2*i)*stride] + B[(2*i+2)*stride])/(float_type)2.0;
	 }
	}

	if (centerpoint*(centerpoint-2) > 0) throw "expandtwice(): Unknown value of centerpoint (should be 0, 1 or 2)";
}



void least_squares(float_type* x, float_type* y, float_type* out, int N, int nN, int xstride, int xdist, int ystride, int ydist, int ostride, int odist)
{
	for (int n=0; n<nN; n++)
	{
		float_type Sx=0, Sxx=0;
		float_type Sy=0, Sxy=0;
		for (int i=0; i<N; i++)
		{
			float_type xc =  x[n*xdist + i*xstride];
			float_type yc =  y[n*ydist + i*ystride];
			Sx += xc;
			Sxx += xc*xc;
			Sxy += xc*yc;
			Sy += yc;
		}
		float_type a = (Sxy*N - Sx*Sy)/(Sxx*N - Sx*Sx);
		float_type b = (Sy - a*Sx)/N;
		out[n*odist] = a;
		out[n*odist+ostride]=b;
	}
}

#ifdef _WIN32
 inline double exp2(double x) {return exp(log(2.0)*x);}
#endif
float_type oddpow(float_type x, int pw)
{
   if (x >= 0) return pow(x,pw);    
   if (x < 0 && pw%2==0) return NAN;
   return (-pow(-x, pw)); 
}

float_type oddroot(float_type x, int pw)
{
   if (x >= 0) return pow(x,((float_type)(1.0/pw)));    
   if (x < 0 && pw%2==0) return NAN;
   return (float_type)(-pow(-x, ((float_type)(1.0/pw)))); 
}


double ellipk(double z)
{
	if (z < 0 || z > 1) throw "ellipk(): Complete elliptical integral of first kind is defined for 0 < z < 1";
	double ap = 1, bp = sqrt(1-z), cp = sqrt(z);
	double an = 0, bn =0,          cn = 0;
	do
	{
		an = 0.5*(ap+bp);
		bn = sqrt(ap*bp);
		cn = 0.5*(ap-bp); 
		ap=an; bp=bn; cp=cn;
	} while (cn > EM_ELLIP_TOLERANCE);
	return M_PI/2.0/an;
}

double ellipe(double z)
{
	if (z < 0 || z > 1) throw "ellipe(): Complete elliptical intgral of second kind is defined for 0 < z < 1";
	double ap = 1, bp = sqrt(1-z), cp = sqrt(z);
	double an = 0, bn = 0,         cn = 0;
	int i = 0;
	do
	{
		an = 0.5*(ap+bp);
		bn = sqrt(ap*bp);
		cn = 0.5*(ap-bp);
		i++;
		z += exp2((double)i)*cn*cn; 
		ap=an; bp=bn; cp=cn;
	} while (cn > EM_ELLIP_TOLERANCE);
	return M_PI/2.0/an*(1-z/2);
}

void ellipke(double z, double* k, double* e)
{
	if (z < 0 || z > 1) throw "ellipke(): Complete elliptical intgrals of both kinds are defined for 0 < z < 1";
	double ap = 1, bp = sqrt(1-z), cp = sqrt(z);
	double an = 0, bn = 0,         cn = 0;
	int i = 0;
	do
	{
		an = 0.5*(ap+bp);
		bn = sqrt(ap*bp);
		cn = 0.5*(ap-bp);
		i++;
		z += exp2((double)i)*cn*cn; 
		ap=an; bp=bn; cp=cn;
	} while (cn > EM_ELLIP_TOLERANCE);
	(*k) =  M_PI/2.0/an;
	(*e) =  (*k)*(1-z/2);
}


double dawson(double x)
{
    /* Initialized data */

    static double zero = 0.;
    static double p1[10] = { .9999999999999999444888487,
            -.1388680862539319882387189,.04701390228872047265251676,
            -.002843881214410084990427962,4.07205792429155824171949e-4,
            -1.238777833290491215723104e-5,9.282648725834448502710848e-7,
            -1.348483044559394160161606e-8,4.185720653743377217952454e-10,
            -2.690203987887047893118563e-12 };
    static double q1[10] = { 1.,.5277985804127346830538769,
            .1322129558972101326386194,.020742277464144764226317,
            .00226061077235076702922345,1.78910965284246248671986e-4,
            1.038676337674144218263657e-5,4.322878276786317666016833e-7,
            1.192668463722972539511844e-8,1.712571708546905554147296e-10 };
    static double p2[9] = { -1.662798629229032210119498,
            -107.998245924983567789468,96.92308277747642719646142,
            4.703418187014092088915617,-14.65360740701534125740579,
            5.313652262936985781749399,6.760560926522734659371849,
            5.149051989461839173856105,.4997537232238673104989246 };
    static double q2[8] = { .4658884381436620841787643,
            10556.53012109847077226732,-2.576680879849772232148552,
            209.4725618926938466302091,57.15518351555917320183654,
            285.0942952341033560514912,-171.8459791160867737858099,
            .2825051295956025332323235 };
    static double p3[10] = { -4.551695032550948205596341,
            -18.66471233384938699373378,-7.363156691268305298336827,
            -66.84072403376967486110516,48.45072650814914538841549,
            26.97905867354676345826192,-33.50441498205924517606035,
            7.509644598389196179510918,-1.484323418233439717539567,
            .4999998109248588312736316 };
    static double q3[9] = { 44.78209080259717467242805,
            99.86071980394520863910665,14.02383731261493848840871,
            3488.177588222863562350534,-9.188713852932158809849738,
            1240.185000099171645615572,-68.80249525045122638289279,
            -2.312515753851451361100544,.2500414923699223884723252 };
    static double p4[9] = { -70.78530823653408887707884,
            23.73886730846445658471564,-8.366401463009488281841186,
            -27.31844677062325033034539,-5.744126756544925882508323,
            -6.42407786948386783087983,-4.499897677567821219213328,
            -2.499999993894251870685252,.4999999999999950039963891 };
    static double q4[8] = { 2319.719562454406798224224,
            -7.723796528921189175065364,202.8581650918903136471269,
            -22.39780498932484675833619,-16.00301283836633459145559,
            -7.009377334838185147347644,-2.500001663946111785108428,
            .7500000000081665924023255 };
    static double one = 1.;
    static double half = .5;
    static double twent5 = 25.;
    static double twel25 = 12.25;
    static double six25 = 6.25;
    static double zcon = .5000000000000000277555756;
    static double sys071 = 1.34078079e155;
    static double sys073 = 2.2227587494851e-162;

    /* System generated locals */
    double ret_val;

    /* Local variables */
    static double frac, sump, sumq;
    static int i__;
    static double y, w2;

    if (fabs(x) > sys071) {
        goto L40;
    }
    if (fabs(x) < sys073) {
        goto L45;
    }
    y = x * x;
    if (y >= six25) {
        goto L10;
    }
/* L5: */
    sump = ((((((((p1[9] * y + p1[8]) * y + p1[7]) * y + p1[6]) * y + p1[5]) *
             y + p1[4]) * y + p1[3]) * y + p1[2]) * y + p1[1]) * y + p1[0];
    sumq = ((((((((q1[9] * y + q1[8]) * y + q1[7]) * y + q1[6]) * y + q1[5]) *
             y + q1[4]) * y + q1[3]) * y + q1[2]) * y + q1[1]) * y + q1[0];
    ret_val = x * sump / sumq;
    goto L50;
L10:
    if (y >= twel25) {
        goto L20;
    }
    frac = zero;
    for (i__ = 1; i__ <= 8; ++i__) {
        frac = q2[i__ - 1] / (p2[i__ - 1] + y + frac);
/* L15: */
    }
    ret_val = (p2[8] + frac) / x;
    goto L50;
L20:
    if (y >= twent5) {
        goto L30;
    }
    frac = zero;
    for (i__ = 1; i__ <= 9; ++i__) {
        frac = q3[i__ - 1] / (p3[i__ - 1] + y + frac);
/* L25: */
    }
    ret_val = (p3[9] + frac) / x;
    goto L50;
L30:
    w2 = one / x / x;
    frac = zero;
    for (i__ = 1; i__ <= 8; ++i__) {
        frac = q4[i__ - 1] / (p4[i__ - 1] + y + frac);
/* L35: */
    }
    frac = p4[8] + frac;
    ret_val = (zcon + half * w2 * frac) / x;
    goto L50;
L40:
    ret_val = half / x;
    goto L50;
L45:
    ret_val = x;
L50:
    return ret_val;
} /* sf12d_c */


