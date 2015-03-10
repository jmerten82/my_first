#ifndef   	LITTLEASTRO_H_
# define   	LITTLEASTRO_H_

#include <iostream>
#include <fstream>
#include <gsl/gsl_math.h>
#include <gsl/gsl_integration.h>

/**
  * Class implementing few elementary routines for a computing the
  * angular-diameter distance in a Friedmann model.
  */

class cosModel
{
  protected:
    static double const epsilon, tiny, tolerance;
    double OmegaMatter, OmegaLambda, wDarkEnergy;
    double expansionFunction (double a);
    double gsl_integrateQag
      (double (*function) (double, void*), double, double);
  public:
    /**
      * Constructor: Initialise the cosmological model with present-day
      * density parameter and cosmological constant. Contributions by
      * relativistic matter are neglected. A constant equation-of-state
      * parameter can be given.
      */
    cosModel (double Omega0 = 0.3, double Lambda0 = 0.7, double wDE = -1.0);

    /**
      * For other than the simplest cases, the angular-diameter distance
      * needs to be numerically integrated. This function returns the
      * integral kernel.
      */
    double angularDistanceKernel (double x);

    /**
      * Returns the angular-diameter distance between the two redshifts
      * zLower and zUpper > zLower. It is given in units of the Hubble
      * radius c/H_0.
      */
    double angularDistance (double zLower, double zUpper);
};

/**
  * Not a class member. Needed as an external integration kernel
  * to be used by the GSL integration routine.
  */

double gsl_angularDistanceKernel (double x, void *p);


#endif 	    /* !LITTLEASTRO_H_ */


