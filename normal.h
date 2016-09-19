#ifndef NORMAL_H
#define NORMAL_H

#include <stddef.h>
#include <float.h>
#include <math.h>

// interface
static double cdf(double x);
static double logcdf(double x);


#define GSL_ROOT6_DBL_EPSILON  2.4607833005759251e-03
#define M_SQRTPI   1.77245385090551602729816748334      /* sqrt(pi) */

#ifndef M_PI
#define M_PI       3.14159265358979323846264338328      /* pi */
#endif

#define GAUSS_EPSILON  (DBL_EPSILON / 2)
#define GAUSS_XUPPER (8.572)
#define GAUSS_XLOWER (-37.519)

#define GAUSS_SCALE (16.0)

static double
get_del (double x, double rational)
{
  double xsq = 0.0;
  double del = 0.0;
  double result = 0.0;

  xsq = floor (x * GAUSS_SCALE) / GAUSS_SCALE;
  del = (x - xsq) * (x + xsq);
  del *= 0.5;

  result = exp (-0.5 * xsq * xsq) * exp (-1.0 * del) * rational;

  return result;
}

#ifndef M_1_SQRT2PI
#define M_1_SQRT2PI (M_2_SQRTPI * M_SQRT1_2 / 2.0)
#endif

#define SQRT32 (4.0 * M_SQRT2)

// inline static double cdf(double x)
// {
//     return (1 + erf(x/sqrt(2))) / 2;
// }

/*
 * Normal cdf for fabs(x) < 0.66291
 */
static double
gauss_small (const double x)
{
  unsigned int i;
  double result = 0.0;
  double xsq;
  double xnum;
  double xden;

  const double a[5] = {
    2.2352520354606839287,
    161.02823106855587881,
    1067.6894854603709582,
    18154.981253343561249,
    0.065682337918207449113
  };
  const double b[4] = {
    47.20258190468824187,
    976.09855173777669322,
    10260.932208618978205,
    45507.789335026729956
  };

  xsq = x * x;
  xnum = a[4] * xsq;
  xden = xsq;

  for (i = 0; i < 3; i++)
    {
      xnum = (xnum + a[i]) * xsq;
      xden = (xden + b[i]) * xsq;
    }

  result = x * (xnum + a[3]) / (xden + b[3]);

  return result;
}

/*
 * Normal cdf for 0.66291 < fabs(x) < sqrt(32).
 */
static double
gauss_medium (const double x)
{
  unsigned int i;
  double temp = 0.0;
  double result = 0.0;
  double xnum;
  double xden;
  double absx;

  const double c[9] = {
    0.39894151208813466764,
    8.8831497943883759412,
    93.506656132177855979,
    597.27027639480026226,
    2494.5375852903726711,
    6848.1904505362823326,
    11602.651437647350124,
    9842.7148383839780218,
    1.0765576773720192317e-8
  };
  const double d[8] = {
    22.266688044328115691,
    235.38790178262499861,
    1519.377599407554805,
    6485.558298266760755,
    18615.571640885098091,
    34900.952721145977266,
    38912.003286093271411,
    19685.429676859990727
  };

  absx = fabs (x);

  xnum = c[8] * absx;
  xden = absx;

  for (i = 0; i < 7; i++)
    {
      xnum = (xnum + c[i]) * absx;
      xden = (xden + d[i]) * absx;
    }

  temp = (xnum + c[7]) / (xden + d[7]);

  result = get_del (x, temp);

  return result;
}

/*
 * Normal cdf for
 * {sqrt(32) < x < GAUSS_XUPPER} union { GAUSS_XLOWER < x < -sqrt(32) }.
 */
static double
gauss_large (const double x)
{
  int i;
  double result;
  double xsq;
  double temp;
  double xnum;
  double xden;
  double absx;

  const double p[6] = {
    0.21589853405795699,
    0.1274011611602473639,
    0.022235277870649807,
    0.001421619193227893466,
    2.9112874951168792e-5,
    0.02307344176494017303
  };
  const double q[5] = {
    1.28426009614491121,
    0.468238212480865118,
    0.0659881378689285515,
    0.00378239633202758244,
    7.29751555083966205e-5
  };

  absx = fabs (x);
  xsq = 1.0 / (x * x);
  xnum = p[5] * xsq;
  xden = xsq;

  for (i = 0; i < 4; i++)
    {
      xnum = (xnum + p[i]) * xsq;
      xden = (xden + q[i]) * xsq;
    }

  temp = xsq * (xnum + p[4]) / (xden + q[4]);
  temp = (M_1_SQRT2PI - temp) / absx;

  result = get_del (x, temp);

  return result;
}

double
cdf (const double x)
{
  double result;
  double absx = fabs (x);

  if (absx < GAUSS_EPSILON)
    {
      result = 0.5;
      return result;
    }
  else if (absx < 0.66291)
    {
      result = 0.5 + gauss_small (x);
      return result;
    }
  else if (absx < SQRT32)
    {
      result = gauss_medium (x);

      if (x > 0.0)
        {
          result = 1.0 - result;
        }

      return result;
    }
  else if (x > GAUSS_XUPPER)
    {
      result = 1.0;
      return result;
    }
  else if (x < GAUSS_XLOWER)
    {
      result = 0.0;
      return result;
    }
  else
    {
      result = gauss_large (x);

      if (x > 0.0)
        {
          result = 1.0 - result;
        }
    }

  return result;
}

static double erfc8_sum(double x)
{
  /* estimates erfc(x) valid for 8 < x < 100 */
  /* This is based on index 5725 in Hart et al */

  static double P[] = {
      2.97886562639399288862,
      7.409740605964741794425,
      6.1602098531096305440906,
      5.019049726784267463450058,
      1.275366644729965952479585264,
      0.5641895835477550741253201704
  };
  static double Q[] = {
      3.3690752069827527677,
      9.608965327192787870698,
      17.08144074746600431571095,
      12.0489519278551290360340491,
      9.396034016235054150430579648,
      2.260528520767326969591866945,
      1.0
  };
  double num=0.0, den=0.0;
  int i;

  num = P[5];
  for (i=4; i>=0; --i) {
      num = x*num + P[i];
  }
  den = Q[6];
  for (i=5; i>=0; --i) {
      den = x*den + Q[i];
  }

  return num/den;
}

inline static double log_erfc8(double x)
{
  double e;
  e = erfc8_sum(x);
  e = log(e) - x*x;
  return e;
}

inline static double erfc8(double x)
{
  double e;
  e = erfc8_sum(x);
  e *= exp(-x*x);
  return e;
}

/* data for a Chebyshev series over a given interval */

typedef struct {

  double * c;   /* coefficients                */
  size_t order; /* order of expansion          */
  double a;     /* lower interval point        */
  double b;     /* upper interval point        */

  /* The following exists (mostly) for the benefit
   * of the implementation. It is an effective single
   * precision order, for use in single precision
   * evaluation. Users can use it if they like, but
   * only they know how to calculate it, since it is
   * specific to the approximated function. By default,
   * order_sp = order.
   * It is used explicitly only by the gsl_cheb_eval_mode
   * functions, which are not meant for casual use.
   */
  size_t order_sp;

  /* Additional elements not used by specfunc */

  double * f;   /* function evaluated at chebyschev points  */
} cheb_series;

inline static double cheb_eval(const cheb_series * cs, const double x)
{
  size_t i;
  double d1 = 0.0;
  double d2 = 0.0;

  double y = (2.0 * x - cs->a - cs->b) / (cs->b - cs->a);
  double y2 = 2.0 * y;

  for (i = cs->order; i >= 1; i--)
    {
      double temp = d1;
      d1 = y2 * d1 - d2 + cs->c[i];
      d2 = temp;
    }

  return y * d1 - d2 + 0.5 * cs->c[0];
}

/* Chebyshev fit for erfc((t+1)/2), -1 < t < 1
 */
static double erfc_xlt1_data[20] = {
  1.06073416421769980345174155056,
 -0.42582445804381043569204735291,
  0.04955262679620434040357683080,
  0.00449293488768382749558001242,
 -0.00129194104658496953494224761,
 -0.00001836389292149396270416979,
  0.00002211114704099526291538556,
 -5.23337485234257134673693179020e-7,
 -2.78184788833537885382530989578e-7,
  1.41158092748813114560316684249e-8,
  2.72571296330561699984539141865e-9,
 -2.06343904872070629406401492476e-10,
 -2.14273991996785367924201401812e-11,
  2.22990255539358204580285098119e-12,
  1.36250074650698280575807934155e-13,
 -1.95144010922293091898995913038e-14,
 -6.85627169231704599442806370690e-16,
  1.44506492869699938239521607493e-16,
  2.45935306460536488037576200030e-18,
 -9.29599561220523396007359328540e-19
};

static cheb_series erfc_xlt1_cs = {
  erfc_xlt1_data,
  19,
  -1, 1,
  12
};

/* Chebyshev fit for erfc(x) exp(x^2), 1 < x < 5, x = 2t + 3, -1 < t < 1
 */
static double erfc_x15_data[25] = {
  0.44045832024338111077637466616,
 -0.143958836762168335790826895326,
  0.044786499817939267247056666937,
 -0.013343124200271211203618353102,
  0.003824682739750469767692372556,
 -0.001058699227195126547306482530,
  0.000283859419210073742736310108,
 -0.000073906170662206760483959432,
  0.000018725312521489179015872934,
 -4.62530981164919445131297264430e-6,
  1.11558657244432857487884006422e-6,
 -2.63098662650834130067808832725e-7,
  6.07462122724551777372119408710e-8,
 -1.37460865539865444777251011793e-8,
  3.05157051905475145520096717210e-9,
 -6.65174789720310713757307724790e-10,
  1.42483346273207784489792999706e-10,
 -3.00141127395323902092018744545e-11,
  6.22171792645348091472914001250e-12,
 -1.26994639225668496876152836555e-12,
  2.55385883033257575402681845385e-13,
 -5.06258237507038698392265499770e-14,
  9.89705409478327321641264227110e-15,
 -1.90685978789192181051961024995e-15,
  3.50826648032737849245113757340e-16
};
static cheb_series erfc_x15_cs = {
  erfc_x15_data,
  24,
  -1, 1,
  16
};

/* Chebyshev fit for erfc(x) x exp(x^2), 5 < x < 10, x = (5t + 15)/2, -1 < t < 1
 */
static double erfc_x510_data[20] = {
  1.11684990123545698684297865808,
  0.003736240359381998520654927536,
 -0.000916623948045470238763619870,
  0.000199094325044940833965078819,
 -0.000040276384918650072591781859,
  7.76515264697061049477127605790e-6,
 -1.44464794206689070402099225301e-6,
  2.61311930343463958393485241947e-7,
 -4.61833026634844152345304095560e-8,
  8.00253111512943601598732144340e-9,
 -1.36291114862793031395712122089e-9,
  2.28570483090160869607683087722e-10,
 -3.78022521563251805044056974560e-11,
  6.17253683874528285729910462130e-12,
 -9.96019290955316888445830597430e-13,
  1.58953143706980770269506726000e-13,
 -2.51045971047162509999527428316e-14,
  3.92607828989125810013581287560e-15,
 -6.07970619384160374392535453420e-16,
  9.12600607264794717315507477670e-17
};
static cheb_series erfc_x510_cs = {
  erfc_x510_data,
  19,
  -1, 1,
  12
};

inline static double gsl_sf_erfc(double x)
{
  const double ax = fabs(x);
  double e_val;

  /* CHECK_POINTER(result) */

  if(ax <= 1.0) {
    double t = 2.0*ax - 1.0;
    e_val = cheb_eval(&erfc_xlt1_cs, t);
  }
  else if(ax <= 5.0) {
    double ex2 = exp(-x*x);
    double t = 0.5*(ax-3.0);
    e_val = ex2 * cheb_eval(&erfc_x15_cs, t);
  }
  else if(ax < 10.0) {
    double exterm = exp(-x*x) / ax;
    double t = (2.0*ax - 15.0)/5.0;
    e_val = exterm * cheb_eval(&erfc_x510_cs, t);
  }
  else {
    e_val = erfc8(ax);
  }

  if(x < 0.0)
    return 2.0 - e_val;
  else
    return e_val;
}

inline static double logerfc(double x)
{
  /* CHECK_POINTER(result) */

  if(x*x < 10.0*GSL_ROOT6_DBL_EPSILON) {
    const double y = x / M_SQRTPI;
    /* series for -1/2 Log[Erfc[Sqrt[Pi] y]] */
    const double c3 = (4.0 - M_PI)/3.0;
    const double c4 = 2.0*(1.0 - M_PI/3.0);
    const double c5 = -0.001829764677455021;  /* (96.0 - 40.0*M_PI + 3.0*M_PI*M_PI)/30.0  */
    const double c6 =  0.02629651521057465;   /* 2.0*(120.0 - 60.0*M_PI + 7.0*M_PI*M_PI)/45.0 */
    const double c7 = -0.01621575378835404;
    const double c8 =  0.00125993961762116;
    const double c9 =  0.00556964649138;
    const double c10 = -0.0045563339802;
    const double c11 =  0.0009461589032;
    const double c12 =  0.0013200243174;
    const double c13 = -0.00142906;
    const double c14 =  0.00048204;
    double series = c8 + y*(c9 + y*(c10 + y*(c11 + y*(c12 + y*(c13 + c14*y)))));
    series = y*(1.0 + y*(1.0 + y*(c3 + y*(c4 + y*(c5 + y*(c6 + y*(c7 + y*series)))))));
    return -2.0 * series;
  }
  else if(x > 8.0) {
    return log_erfc8(x);
  }
  else {
    return log(gsl_sf_erfc(x));
  }
}

/* Implements log(e^x + e^y).
 */
inline static double logaddexp(double x, double y)
{
    double tmp = x - y;

    if (x == y)
        return x + M_LN2;

    if (tmp > 0)
        return x + log1p(exp(-tmp));
    else if (tmp <= 0)
        return y + log1p(exp(tmp));

    return tmp;
}

/* Implements log(sx * e^x + sy * e^y).
 *
 * It assumes that sx * e^x + sy * e^y > 0.
 */
inline static double logaddexps(double x, double y, double sx, double sy)
{
    double tmp = x - y;

    double sxx = log(fabs(sx)) + x;
    double syy = log(fabs(sy)) + y;

    if (sxx == syy)
    {
        if (sx * sy > 0)
            return sxx + M_LN2;
        return -DBL_MAX;
    }

    if (sx > 0 && sy > 0)
    {
        if (tmp > 0)
            return sxx + log1p((sy/sx) * exp(-tmp));
        else if (tmp <= 0)
            return syy + log1p((sx/sy) * exp(tmp));
    }
    else if (sx > 0)
        return sxx + log1p((sy/sx) * exp(-tmp));
    else
        return syy + log1p((sx/sy) * exp(tmp));
    return tmp;
}

inline static double logerf(double x)
{
    return logaddexps(0, logerfc(x), 1, -1);
}

// static double logcdf(double x)
// {
//     if (x > 0)
//         return logaddexp(0, logerf(x/sqrt(2))) - log(2);
//     return logaddexps(0, logerfc(x/sqrt(2)), 1, -0.5);
// }

static double logcdf(double a)
{

    double log_LHS,		/* we compute the left hand side of the approx (LHS) in one shot */
     last_total = 0,		/* variable used to check for convergence */
	right_hand_side = 1,	/* includes first term from the RHS summation */
	numerator = 1,		/* numerator for RHS summand */
	denom_factor = 1,	/* use reciprocal for denominator to avoid division */
	denom_cons = 1.0 / (a * a);	/* the precomputed division we use to adjust the denominator */
    long sign = 1, i = 0;

    if (a > 6) {
	return -cdf(-a);     /* log(1+x) \approx x */
    }
    if (a > -20) {
	return log(cdf(a));
    }
    log_LHS = -0.5 * a * a - log(-a) - 0.5 * log(2 * M_PI);

    while (fabs(last_total - right_hand_side) > DBL_EPSILON) {
	i += 1;
	last_total = right_hand_side;
	sign = -sign;
	denom_factor *= denom_cons;
	numerator *= 2 * i - 1;
	right_hand_side += sign * numerator * denom_factor;

    }
    return log_LHS + log(right_hand_side);
}

#endif
