#include "tposrestr.h"

//---------------------------------------------------------------------------------------------------
double TPosRestr::dLength_SQR() const
{
  return dDistance_SQR(patAtom1, dX0, dY0, dZ0);
}

//---------------------------------------------------------------------------------------------------
double TPosRestr::dLength() const
{
  return sqrt(dLength_SQR());
}

//---------------------------------------------------------------------------------------------------
double TPosRestr::dCalcEnergy() const
{
  double dEfficientDistance = dLength();
  return dKr*dEfficientDistance*dEfficientDistance;
}

//---------------------------------------------------------------------------------------------------
void TPosRestr::CalcGradients() const
{
  patAtom1->vecEnergyGradient -= 2.*dKr*vecMakeVector(patAtom1, dX0, dY0, dZ0);
}

//---------------------------------------------------------------------------------------------------
