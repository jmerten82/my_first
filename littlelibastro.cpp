#include <littlelibastro.h>

using namespace std;

/**
  * Implementation of the class members.
  */

cosModel::cosModel (double Omega0, double Lambda0, double wDE)
: OmegaMatter (Omega0), OmegaLambda (Lambda0), wDarkEnergy (wDE)
{}

double cosModel::gsl_integrateQag
  (double (*function) (double, void*), double a, double b )
{
  gsl_function gFunction;
  gFunction.function = function;
  gFunction.params = this;
  double result, error;
  const size_t nPoints = 64;
  gsl_integration_workspace * work = gsl_integration_workspace_alloc (nPoints);
  gsl_integration_qag
    (&gFunction, a, b, tiny, tolerance, nPoints, 1, work, &result, &error);
  gsl_integration_workspace_free (work);
  return result;
}

double cosModel::angularDistance (double zLower, double zUpper)
{
  double curvature = OmegaMatter + OmegaLambda -1.0;
  int curvatureFlag = 0;
  if (curvature < 0.0)
    curvatureFlag = -1;
  if (curvature > 0.0)
    curvatureFlag = 1;
  if (fabs (curvature) <= epsilon)
    curvatureFlag = 0;
  curvature = sqrt (fabs (curvature));
  double distance = 0.0;

  if (OmegaLambda <= epsilon)
  {
    if (fabs (OmegaMatter-1.0) <= epsilon)
      distance =
        2.0/(1.0+zUpper)*(1.0/sqrt(1.0+zLower)-1.0/sqrt(1.0+zUpper));
    else
      distance =
        2.0/gsl_pow_2 (OmegaMatter)/(1.0+zLower)/gsl_pow_2 (1.0+zUpper)*
          (sqrt(1.0+OmegaMatter*zLower)*(2.0-OmegaMatter+OmegaMatter*zUpper)-
           sqrt(1.0+OmegaMatter*zUpper)*(2.0-OmegaMatter+OmegaMatter*zLower));
  }
  else
  {
    distance = gsl_integrateQag (&gsl_angularDistanceKernel, 1.0+zLower, 1.0+zUpper);
    if (curvature*distance > epsilon && curvatureFlag == -1)
      distance = sinh (curvature*distance)/curvature;
    if (curvature*distance > epsilon && curvatureFlag == 1)
      distance = sin(curvature*distance)/curvature;
    distance /= (1.0+zUpper);
  }
  return distance;
}

double cosModel::expansionFunction (double a)
{ 
  double root =
    OmegaMatter/gsl_pow_3 (a) +
    OmegaLambda/pow (a, 3.0*(1.0+wDarkEnergy)) +
    (1.0-OmegaMatter-OmegaLambda)/gsl_pow_2 (a);
  return sqrt (root);
}

double cosModel::angularDistanceKernel (double x)
{ return 1.0/expansionFunction (1.0/x); }


double gsl_angularDistanceKernel (double x, void *p)
{
  cosModel * cosmology = reinterpret_cast <cosModel *> (p);
  return cosmology->angularDistanceKernel (x);
}

/**
  * Static constants controlling the internal precision.
  */

double const cosModel::epsilon = 0.01;
double const cosModel::tolerance = 1.0e-4;
double const cosModel::tiny = 1.0e-6;
