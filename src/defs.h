/* optionmatrix:                                                            */
/*                                                                          */
/* Options & Futures Matrix Modeler                                         */
/* View and Control Theoretical Option Chains                               */
/*                                                                          */
/* File: defs.h of optionmatrix                                             */
/*                                                                          */
/* Copyright (c) Anthony Bradford. 2012.                                    */
/* http://opensourcefinancialmodels.com                                     */
/* info@opensourcefinancialmodels.com                                       */
/*                                                                          */
/* optionmatrix may be freely redistributed.                                */
/* See file COPYING included with this distribution for license information */

/* 
   optionmatrix is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 3 of the License, or
   (at your option) any later version.

   optionmatrix program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "../config.h"

#include <stdio.h>
#include <ctype.h>
#include <time.h>
#include <pthread.h>

#include <map>
#include <string>
#include <vector>

#ifdef HAVE_LOCALE_H
# include <locale.h>
#else
# error Sorry, this code requires <locale.h>
#endif

#ifdef HAVE_STDLIB_H
# include <stdlib.h>
#else
# error Sorry, this code requires <stdlib.h>
#endif

#ifdef HAVE_LIBM
# include <math.h>
#else
# error Sorry, this code requires the 'm' library (-lm)
#endif

#ifdef HAVE_FLOAT_H
# include <float.h>
#else
# error Sorry, this code requires <float.h>
#endif

#ifdef HAVE_STDDEF_H
# include <stddef.h>
#else
# error Sorry, this code requires <stddef.h>
#endif

#ifdef HAVE_UNISTD_H
# include <unistd.h>
#else
# error Sorry, this code requires <unistd.h>
#endif

#ifdef HAVE_STRING_H
# include <string.h>
#endif

#ifdef HAVE_STRINGS_H
# include <strings.h>
#endif

#ifdef HAVE_LIBGSL
# include "gsl/gsl_errno.h"
#endif

#ifndef HAVE_LIBPTHREAD
# error Sorry, missing pthread library which is needed by this code
#endif

#ifndef HAVE_MEMSET
# error Sorry, missing memset which is needed by this code
#endif

#ifndef HAVE_MODF
# error Sorry, missing modf which is needed by this code
#endif

#ifndef HAVE_POW
# error Sorry, missing pow which is needed by this code
#endif

#ifndef HAVE_SQRT
# error Sorry, missing sqrt which is needed by this code
#endif

#ifndef HAVE_LOCALECONV
# error Sorry, missing localeconv which is needed by this code
#endif

#include "gtk_include.h"

/* models to use. comment out to exclude */
#define ABRADFORD
#define SPINSKY
#define METAOPTIONS
#define FINRECIPES
//#define DUMMYTESTMODELS     // For testing

using namespace std;

enum { OPTION_CLASS = 0, FUTURES_CLASS, BOND_CLASS, TERMSTRUCTURE_CLASS };

// date securities via decimal/date entry dating or date to calendar engine ( e.g. 3rd friday + 1 )..
enum { DECIMALS, CALENDAR };

enum { LEG1, LEG2 };

// screen style formats
enum { 
       // begin screens for options
       DECIMAL_GREEKS = 0, DECIMAL_OPTIONS, 
       CALENDAR_OPTIONS1,  CALENDAR_OPTIONS2,
       CALENDAR_OPTIONS3,  CALENDAR_OPTIONS4,
       CALENDAR_OPTIONS5,  CALENDAR_OPTIONS6,
       CALENDAR_CUSTOM,
       OPTIONS_FORMAT_END,

       // begin screens for futures
       DECIMAL_FUTURE, CALENDAR_FUTURE,
       FUTURES_FORMAT_END,

       // next are properties, keep separate from formats
       // by FUTURES_FORMAT_END 

       DEMO_OPTIONS, DEMO_FUTURES, 
       OPTIONS_CALENDAR, FUTURES_CALENDAR,
       PROPERTIES,

       // always keep QUIT as last enum
       QUIT,

     };

enum { NORMAL_DISTRIBUTION = 0 };

enum { POLYNOMIAL_APPROX6, POLYNOMIAL_APPROX4, RATIONAL_APPROX7, RATIONAL_APPROX4, HART,
       ROMBERGS_METHOD, SIMPSONS_RULE, TRAPEZOID_RULE, ERF };
       //SIMPSONS_RULE, TRAPEZOID_RULE, ERF };

struct _data {

  // input data
  int modeltype;
  int term_model;

  double price;
  double strike;
  double rate;
  double volatility;
  double t[3];
  double dividend;

  pthread_mutex_t mutexCashflow;

  vector<double> amounts;
  vector<double> times;
  vector<double> times_adjusted;

  vector<double> coupon;
  vector<double> coupon_times;
  vector<double> coupon_times_adjusted;

  vector<double> principal;
  vector<double> principal_times;
  vector<double> principal_times_adjusted;

  vector<double> *generic_amounts;
  vector<double> *generic_times;
  vector<double> *generic_times_adjusted;

#ifdef FINRECIPES
  class term_structure_class *term;
#endif

  double UseZ;
  double UseB;
  double UseJ;
  double UseP;
  double UseQ;
  double UseR;
  double UseS;
  double UseT;
  int UsePound;  // '#'

  int steps;

  // computed data

  double te;  // time elapsed in decimal since user input to minus from t
  double te2; // time elapsed in decimal since user input to minus from t2
  double te3; // time elapsed in decimal since user input to minus from t3
  double te4;

  double call;
  double put;

  double future;

  double calldelta;
  double putdelta;

  double gamma;
  double vega;

  double calltheta;
  double puttheta;

  double callrho;
  double putrho;

  // Term Stuctures
  double discount_t1;
  double discount_t2;
  double spot_t1;
  double spot_t2;
  double forward;

  // Bonds
  double bond_price;
  double pv_discrete;
  double pv_continous;
  double irr;
  double irr_discrete;
  bool   uirr;
  double YTMDiscrete;
  double YTMContinous;
  double durationContinous;
  double durationDiscrete;

  double durationMacaulayDiscrete;
  double durationModifiedDiscrete;

  double convexityContinous;
  double convexityDiscrete;

  bool debug;
  bool isnan;

};

typedef std::map<std::string, int> treeToIndex;

struct _properties {
   
    struct _data data;

    int decimalorcalendar;

    int modeltype;
    int occurence_day;
    int occurence_in_month;
    int occurence_plus_offset;
    int expiration_hour,expiration_minute,expiration_second;
    unsigned int expiration_time;
    double adjusttime;
  
    // Leg 1 or single option or future
    int *time_to_expiration,*days_to_expiration,*expiration_month,*expiration_year,*start_time_to_expiration,*start_days_to_expiration,*start_expiration_month,*start_expiration_year;

    // Leg 2
    int *time_to_expiration2,*days_to_expiration2,*expiration_month2,*expiration_year2,*start_time_to_expiration2,*start_days_to_expiration2,*start_expiration_month2,*start_expiration_year2;

    // All Months 0, Jan cycle 1, Feb cycle 2, Mar cycle 3
    int optionscycle;
    int skipmonth;   // Leg 1
    int skipmonth2;  // Leg 2
    int day_offset_counter;
    double day_offset;

    time_t starttime;
    time_t starttime2;
    time_t starttime3;
    time_t starttime4;

    time_t updatetime;

    int calendarSelected;

    bool realTimeBleeding;
    int updatedelayseconds;
    bool highlightInTheMoney;
    bool spreads;

    int precision;
    int strikestoogle;
    double customstrike;
    double strike_offset;
    double strike_offset2;
    double discard;
    bool verticalSpread;

    struct _int_to_name_nonconst *termstructure_name_list; 

    int format;

    int distribution_type;
    double distribution_mean;
    double distribution_variance;
    int integration_type;

    bool textExport;
    struct lconv * lc;

    // curses specific
    int start;
    int mode;               // incrementor/decrementor 1 = small, 2 = medium, 3 = large
    int demotype;

    // GTK+ specific
    struct gtkinfo GtkInfo;
    struct elementListWithGroup *listModelsForGroups;
    treeToIndex TreeToIndex;
};

struct integratable_distributions {

  double (*constantvalue) (double);
  double (*integrationpart) (double,double,double);
  char   des[30];

};

// list of numerical integration methods
struct numerical_integration_method {

  double (*method) (const double a, const double b, int n, double (*fun) (double,double,double),const double parm2, const double parm3);

  char           des[30];
  mutable int    resolution;
  bool           allowOnlyEven;
  int            UpperLimit;

};

/* list of all financial models */
enum { 

#ifdef ABRADFORD

       BLACKSCHOLES,           // by Anthony Bradford
       MERTON73,               // by Anthony Bradford
       BLACK76,                // by Anthony Bradford

#endif

#ifdef SPINSKY

       AMERICAN_CRR,           // by Seth Pinsky
       AMERICAN_EQUIPROB,      // by Seth Pinsky
       EURO_CRR,               // by Seth Pinsky
       EURO_EQUIPROB,          // by Seth Pinsky

#endif

#ifdef METAOPTIONS
  
       BLACKSCHOLES2,              // by metaoptions
       GBS,                        // by metaoptions
       MERTON73B,                  // by metaoptions
       GFRENCH,                    // by metaoptions
       GCARRY,                     // by metaoptions
       BAWAMERICAN,                // by metaoptions
       BSAMERICAN,                 // by metaoptions
       BLACK76B,                   // by metaoptions
       ASSETORNOTHING,             // by metaoptions
       EXECUTIVE,                  // by metaoptions
       GARMANKOHLHAGEN,            // by metaoptions
       CASHORNOTHING,              // by metaoptions
       GAPOPTION,                  // by metaoptions
       SUPERSHARE,                 // by metaoptions
       SUPERSHARE2,                // by metaoptions
       VASICEKBONDPRICE,           // by metaoptions
       VASICEKBONDOPTION,          // by metaoptions
       TIMESWITCHOPTION,           // by metaoptions
       FOREQUOPTINDOMCUR,          // by metaoptions
       QUANTO,                     // by metaoptions
       EQUITYLINKEDFXO,            // by metaoptions
       SPREADAPPROXIMATION,        // by metaoptions
       JUMPDIFFUSION,              // by metaoptions
       BISECTION,                  // by metaoptions - computes volatility
       NEWTONRAPHSON,              // by metaoptions - computes volatility
       BAWBISECTION,               // by metaoptions - computes volatility
       BSBISECTION,                // by metaoptions - computes volatility
       AMERICANEXCHANGEOPTION,     // by metaoptions
       EUROPEANEXCHANGEOPTION,     // by metaoptions
       MILTERSENSWARTZ,            // by metaoptions
       PARTIALTIMETWOASSETBARRIER, // by metaoptions
       TAKEOVERFXOPTION,           // by metaoptions
       TWOASSETBARRIER,            // by metaoptions
       TWOASSETCASHORNOTHING,      // by metaoptions
       TWOASSETCORRELATION,        // by metaoptions
       FLOATINGSTRIKELOOKBACK,     // by metaoptions
       FLOATINGSTRIKELOOKBACK2,    // by metaoptions
       FIXEDSTRIKELOOKBACK,        // by metaoptions
       OPTIONSONTHEMAXMIN,         // by metaoptions
       STANDARDBARRIER,            // by metaoptions
       DOUBLEBARRIER,              // by metaoptions
       SOFTBARRIER,                // by metaoptions
       BINARYBARRIER,              // by metaoptions has 28 int states
       DISCRETEADJUSTEDBARRIER,    // by metaoptions

       // Crashes - Segmentation fault
       //BARRIERBINOMINAL,           // BarrierBinominal.c

       CONVERTIBLEBOND,            // ConvertibleBond.c
       CRRBINOMINAL,               // CRRBinominal.c

       // Missing Transpose() function...
       //IMPLIEDTRINOMINALTREE,      // ImpliedTrinominalTree.c  - computes volatility

       THREEDIMENSIONALBINOMINAL,  // ThreeDimensionalBinominal.c
       TRINOMINALTREE,             // TrinominalTree.c
       LOOKBARRIER,                // by metaoptions
       PARTIALTIMEBARRIER,         // has multiple text states, some call abort()?
       ROLLGESKEWHALEY,            // by metaoptions
       EXTREMESPREADOPTION,        // by metaoptions
       EXTREMESPREADOPTION2,       // by metaoptions
       PARTIALFIXEDLB,             // by metaoptions
       PARTIALFLOATLB,             // by metaoptions
       PARTIALFLOATLB2,            // by metaoptions
       EXTENDIBLEWRITER,           // by metaoptions
       CALLSONOPTIONS,             // by metaoptions OptionsOnOptions.c
       PUTSONOPTIONS,              // by metaoptions OptionsOnOptions.c
       LEVYASIAN,                  // by metaoptions
       GEOMETRICAVERAGERATEOPTION, // by metaoptions
       FORWARDSTARTOPTION,         // by metaoptions
       SWAPOPTION,                 // by metaoptions
       TURNBULLWAKEMANASIAN,       // by metaoptions

       // Crashes with positive dividend. Dividend disabled
       EXCHANGEEXCHANGEOPTION,     // by metaoptions
       SIMPLECHOOSER,              // by metaoptions
       COMPLEXCHOOSER,             // by metaoptions uses 3 bleeding time vars

#endif

#ifdef FINRECIPES

       BLACKSCHOLES3,          // option_price_call_black_scholes()
       PERPETUAL,              // option_price_american_perpetual_call()
       FUTOPTEURBLACK,         // futures_option_price_call_european_black()
       FUTOPTAMBINOMIAL,       // futures_option_price_call_american_binomial()

       // prototyped in fin_recipes.h but not defined in the library
       //AMERBJERKSUNDSTENSLAND, // option_price_american_call_approximated_bjerksund_stensland()

       AMERBINOMIAL,           // option_price_call_american_binomial() signature 1 without dividend
       AMERBINOMIALDIV,        // option_price_call_american_binomial() signature 2 with dividend
       AMERICANTRINOMIAL,      // option_price_call_american_trinomial()
       AMERBAW,                // option_price_american_call_approximated_baw()
       AMERPUTAPPROXJOHNSON,   // option_price_american_put_approximated_johnson()       

#ifdef HAVE_LIBGSL

       AMERPUTAPPROXGESKEJOHNSON, // option_price_american_put_approximated_geske_johnson()
       HESTON,                    // heston_call_option_price()

#endif

       BONDZEROBLACK,          // bond_option_price_call_zero_black_scholes()
       BONDAMERBINOMIAL,       // bond_option_price_call_american_binomial()
       BOND_ZERO_AM_RENDLEMAN_BARTTER, // bond_option_price_call_zero_american_rendleman_bartter()
       BSPAYOUT,               // option_price_european_call_payout()
       EUROBIONMIAL,           // option_price_call_european_binomial()
       ASIANGEOMETRICAVG,      // option_price_asian_geometric_average_price_call()
       EUROLOOKBACK,           // option_price_european_lookback_call()
       EUROLOOKBACK2,          // option_price_european_lookback_call()
       MERTONJUMPDIFF,         // option_price_call_merton_jump_diffusion()
       CURRAMBINOMIAL,         // currency_option_price_call_american_binomial()
       CURREURO,               // currency_option_price_call_european()
       ROLLGESKEWHALEY2,       // option_price_american_call_one_dividend()
       EUROBINOMIAL1P,         // option_price_call_european_binomial_single_period()
       EUROBINOMIALMP,         // option_price_call_european_binomial_multi_period_given_ud()
       AMBINOMIAL,             // option_price_generic_binomial()
       EUROSIM,                // option_price_call_european_simulated()
       SIMEUROGENERIC,         // derivative_price_simulate_european_option_generic()
       SIMEUROGENERICCV,       // derivative_price_simulate_european_option_generic_with_control_variate()
       SIMEUROGENERICAV,       // derivative_price_simulate_european_option_generic_with_antithetic_variate()
       SIMPRICEPATH,           // derivative_price_simulate_european_option_generic()
       SIMPRICEPATHCONTROLVARIATE,  // derivative_price_simulate_european_option_generic_with_control_variate()
       DISTLOGRAND,             // simulate_lognormal_random_variable()
       AMFINITEDIFFEXP,        // option_price_call_american_finite_diff_explicit()
       EUROFINITEDIFFEXP,      // option_price_call_european_finite_diff_explicit()

#ifdef HAVE_NEWMAT_NEWMAT_H

       AMFINDIFFIMP,           // findiff_imp_am_call_newmat.cc findiff_imp_am_put_newmat.cc  
       EURFINDDIFFIMP,         // findiff_imp_eur_call_newmat.cc findiff_imp_eur_put_newmat.cc

#endif

#ifdef HAVE_ITPP_ITBASE_H

       AMFINDIFFIMPPUT,        // option_price_put_american_finite_diff_implicit_itpp()

#endif

       // puts prototyped in fin_recipes.h but not defined in the library
       IMPLIEDNEWTON,          // option_price_implied_volatility_call_black_scholes_newton()

       // puts prototyped in fin_recipes.h but not defined in the library
       IMPLIEDBISECTIONS,      //option_price_implied_volatility_call_black_scholes_bisections()

       BONDZEROVASICEK,        // bond_option_price_call_zero_vasicek()
       EURODIVIDENDS,          // option_price_european_call_dividends()

       // Throws std::bad_alloc and freezes...
       //AMDISDIVSBINOMIAL,      // option_price_call_american_discrete_dividends_binomial()

       AMPROPORTDIVSBINOMIAL,  // option_price_call_american_proportional_dividends_binomial()
       BERMUDIANBINOMIAL,      // option_price_call_bermudan_binomial()
       BSCOUPONBOND,           // bond_option_price_call_coupon_bond_black_scholes()
       WARRANT_NO_DIV,         // warrant_price_adjusted_black_scholes()
       WARRANT_DIV,            // warrant_price_adjusted_black_scholes()

       // Segmentation fault, To implement un-comment EURBOND_HO_LEE reference in defaults.c
       //EURBOND_HO_LEE,

       TERMFLAT,                // term_structure_class_flat.cc
       TERMCIR,                 // termstru_discfact_cir.cc
       TERMVASICEK,             // termstru_discfact_vasicek.cc
       TERMNELSONSIEGEL,        // termstru_yield_nelson_siegel.cc
       TERMSVENSSON,            // termstru_yield_svensson.cc
       TERMCUBICSPLINE,         // termstru_discfact_cubic_spline.cc
       TERMINTERPOLATED,        // termstru_yield_interpolated.cc

       // prototyped in fin_recipes_extra.h but not defined in the library
       //TERMDISESTCIR,          // term_structure_discount_factor_estimated_cir()

       // prototyped in fin_recipes_extra.h but not defined in the library
       //TERMYIELDBLISS          // term_structure_yield_bliss()

       /*
         Not implemented...

         termstru_transforms.cc
         
         simulate_european_options_generic_routine_price_sequence.cc:
                            vector<double>prices = 
                            simulate_lognormally_distributed_sequence(
                                    S,r,sigma,time,no_steps);

         vector< vector<double> > interest_rate_trees_gbm_build(
                            const double& r0,
                            const double& u,
                            const double& d,
                            const int& n);

         double interest_rate_trees_gbm_value_of_cashflows(
                            const vector<double>& cflow,
                            const vector< vector<double> >& r_tree,
                            const double& q);

         double interest_rate_trees_gbm_value_of_callable_bond(
                            const vector<double>& cflows,
                            const vector< vector<double> >& r_tree,
                            const double& q,
                            const int& first_call_time,
                            const double& call_price);

         vector< vector<term_structure_class_ho_lee> > 
                            term_structure_ho_lee_build_term_structure_tree(
                            term_structure_class* initial,
                            const int& no_steps,
                            const double& delta,
                            const double& pi);
        */


       FUTURES,                // futures_price()

#endif

#ifdef ABRADFORD

       FUTURES2,              // by Anthony Bradford futures_price()
       BACHELIER,
       BACHELIERMODIFIED,
       SPRENKLE,
       BONESS,
       SAMUELSON,

#endif

#ifdef FINRECIPES

       BONDS,
       BONDSTERM,
       BONDSPRINCIPAL,
       IRR,

#endif

#ifdef DUMMYTESTMODELS

       TESTOPTION1,
       TESTOPTION2,
       TESTOPTION3,
       TESTOPTION4,
       TESTOPTIONONEDIVIDEND1,
       TESTOPTIONONEDIVIDEND2,
       TESTFUTURES1,
       TESTFUTURES2,

#endif

     };

/* _FUTURES is used for internal calendar control. Don't remove */
enum { _FUTURES = -1 };

struct _int_to_name          { const char string[30]; };
struct _int_to_name_nonconst { char string[30]; };
struct _int_to_function { double (*fun) (const double&, const double&); };

struct option_algorithm {

  int modeltype;
  char des[40];
  char source[80];
  char curses_des[240]; // 240 = 3 lines of text

  // source path should have a separator ';' to include multiple files in 1 reference...
  // there are instances where more than 2 source files need to be displayed...
  char sourceCode[200];
  char sourceCode2[200];
  char category[200];

  int ReservedNotUsed;
  bool supportCND;
  bool supportSteps;
  mutable int steps;
  int assetClass;
  bool perpetual;

  bool supportPrice;
  bool supportRate;
  bool supportVolatility;
  bool supportStrikes;
  bool allowNegativeOptions;
  bool supportCalls;
  bool supportPuts;
  bool produceCallDelta; 
  bool producePutDelta; 
  bool produceGamma;
  bool produceVega;
  bool produceCallTheta;
  bool producePutTheta;
  bool produceCallRho;
  bool producePutRho;
  bool failsOnMeanVarianceChanges;

  /*
    supportDividend:
        0 - no dividend support
        1 - dividend supported
        2 - force initial default value of dividend
        3 - vector<double>& times, vector<double>& amounts
        4 - vector<double>& potential_exercise_times
        5 - vector<double>& coupon_times, vector<double>& coupon_amounts
        6 - vector<double>& coupon_times, vector<double>& coupon_amounts, 
            vector<double>& principal_times, vector<double>& principal_amounts
  */
  int supportDividend;
  double defaultDividend;
  char supportTime1des[50];
  /* 
     supportTime2
        0 - t2 not used
        1 - t2 used, model supports decimal/date entry dating and not calendar dating
        2 - t2 used, model supports decimal/date entry and calendar dating...
  */
  int supportTime2;
  char supportTime2des[50];

  int supportTime3;
  char supportTime3des[50];

  char price[50];
  char dividend[50];
  char call[50];
  char put[50];
  char strike[50];
  char volatility[50];

  /*
    0 = not used
    1 = used as double
    2 = used as int
  */
  int iUseZ;
  bool bZallow0Negative;
  char UseZdes[50];
  double Zdefault;
  double Zmax;

  /*
    0 = not used
    1 = used as double
    2 = used as int
  */
  int iUseB;
  bool bBallow0Negative;
  char UseBdes[50];
  double Bdefault;
  double Bmax;

  /*
    0 = not used
    1 = used as double
    2 = used as int
  */
  int iUseJ;
  bool bJallow0Negative;
  char UseJdes[50];
  double Jdefault;
  double Jmax;

  /*
    0 = not used
    1 = used as double
    2 = used as int
  */
  int iUseP;
  bool bPallow0Negative;
  char UsePdes[50];
  double Pdefault;
  double Pmax;  

  /*
    0 = not used
    1 = used as double
    2 = used as int
  */
  int iUseQ;
  bool bQallow0Negative;
  char UseQdes[50];
  double Qdefault;
  double Qmax;  

  /*
    0 = not used
    1 = used as double
    2 = used as int
  */
  int iUseR;
  bool bRallow0Negative;
  char UseRdes[50];
  double Rdefault;
  double Rmax;  

  /*
    0 = not used
    1 = used as double
    2 = used as int
  */
  int iUseS;
  bool bSallow0Negative;
  char UseSdes[50];
  double Sdefault;
  double Smax;  

  /*
    0 = not used
    1 = used as double
    2 = used as int
  */
  int iUseT;
  bool bTallow0Negative;
  char UseTdes[50];
  double Tdefault;
  double Tmax;  

  /* bUsePound '#' is for option states */
  int bUsePound;
  bool bPoundallow0Negative;
  char UsePounddes[50];
  int Pounddefault;
  int Poundmax;

  bool bUseStateNames;
  mutable const struct _int_to_name *StateNames;

};

struct _strike_control {

  double xcontrol;
  double incrementor;
  int strikes5or1;
  double retdiscard;
  char des[5];
  int precision;
  double sliderScale;

};

struct term_structure {
  int modeltype;
};

struct elementList
{
  char elementName[200];
};

struct elementListWithGroup
{
  char groupName[200];
  char elementName[200];
  int index;
};
