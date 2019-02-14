#include "thbond.h"

//---------------------------------------------------------------------------------------------------
const double THBond::dREq       = 1.55,
	     THBond::dK         = 1.,
	     THBond::dXiR       = 2.,
	     THBond::dXiTheta	= 2.,
	     THBond::dThreshold = -1.;
//---------------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------
double THBond::dCosTheta() const
{// Returns the cosine of angle between parent, donor and acceptor: P---D:::A.
  TVector vLeft  = vecMakeVector(patDonor, patParent), 
	  vRight = vecMakeVector(patDonor, patAcceptor);
  vLeft.Normalize();
  vRight.Normalize();
  
  double dCosAngle = vLeft*vRight;
  if     (dCosAngle >  1.) dCosAngle =  1.;
  else if(dCosAngle < -1.) dCosAngle = -1.;
  
  return dCosAngle;
}
//---------------------------------------------------------------------------------------------------
double THBond::dLength_SQR() const
{
  return dDistance_SQR(patDonor, patAcceptor);
}

//---------------------------------------------------------------------------------------------------
double THBond::dLength() const
{
  return sqrt(dLength_SQR());
}

//---------------------------------------------------------------------------------------------------
double THBond::dCalcEnergy() const
{
  return -dK*exp(-dXiR*(dLength()-dREq)*(dLength()-dREq) - dXiTheta*dCosTheta());
}

//---------------------------------------------------------------------------------------------------
void THBond::CalcGradients() const
{
  TVector vecNorm21 = vecMakeVector(patDonor, patParent),
	  vecNorm23 = vecMakeVector(patDonor, patAcceptor);
  double dLen21 = vecNorm21.dLength(),
	 dLen23 = vecNorm23.dLength();
	 
  vecNorm21 = vecNorm21/dLen21;	// Normalize by hand.
  vecNorm23 = vecNorm23/dLen23;
  
  // Whoooah:
  double dEnergyDepth = -dK*exp(-dXiR*(dLen23-dREq)*(dLen23-dREq) - dXiTheta*dCosTheta());
  TVector vecGrad1 = -dXiTheta*((vecNorm21*vecNorm23)*vecNorm21 - vecNorm23)/dLen21,
	  vecGrad2 =  dXiTheta*((vecNorm21*vecNorm23)*vecNorm21 - vecNorm23)/dLen21
		     +dXiTheta*((vecNorm21*vecNorm23)*vecNorm23 - vecNorm21)/dLen23
		     +2.*dXiR*(dLen23 - dREq)*(vecNorm21*vecNorm23-1.)*vecNorm23,
	  vecGrad3 = -2.*dXiR*(dLen23 - dREq)*(vecNorm21*vecNorm23-1.)*vecNorm23
		     -dXiTheta*((vecNorm21*vecNorm23)*vecNorm23 - vecNorm21)/dLen23;

  if(TConsts::bIsNan((vecGrad1 + vecGrad2 + vecGrad3).dLength())) return;	// Protection against numerical instabilities.

  patParent  ->vecEnergyGradient -= dEnergyDepth*vecGrad1;
  patDonor   ->vecEnergyGradient -= dEnergyDepth*vecGrad2;
  patAcceptor->vecEnergyGradient -= dEnergyDepth*vecGrad3;
}

//---------------------------------------------------------------------------------------------------
