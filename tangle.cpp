#include "tangle.h"

//---------------------------------------------------------------------------------------------------
double TAngle::dValue() const
{
  TVector vLeft  = vecMakeVector(*patAtom2, *patAtom1), 
	  vRight = vecMakeVector(*patAtom2, *patAtom3);
  vLeft.Normalize();
  vRight.Normalize();
  
  double  dProduct = vLeft*vRight;
  if(     dProduct >  1.) dProduct =  1.;
  else if(dProduct < -1.) dProduct = -1.;
  
  return acos(dProduct);
}

//---------------------------------------------------------------------------------------------------
double TAngle::dCalcEnergy() const
{
  double dThetaEfficient = dValue() - apAngleInfo.dThetaEq*TConsts::PI/180.;
  return apAngleInfo.dKTheta*dThetaEfficient*dThetaEfficient;
}

//---------------------------------------------------------------------------------------------------
bool TAngle::bContains(const TBond& _boBond) const
{
  bool bRetVal = false;
  if((patAtom1 == _boBond.patAtom1 && patAtom2 == _boBond.patAtom2) ||
     (patAtom1 == _boBond.patAtom2 && patAtom2 == _boBond.patAtom1) ||
     (patAtom2 == _boBond.patAtom1 && patAtom3 == _boBond.patAtom2) ||
     (patAtom2 == _boBond.patAtom2 && patAtom3 == _boBond.patAtom1) ) bRetVal = true;
    
  return bRetVal;
}

//---------------------------------------------------------------------------------------------------
bool TAngle::bDihedralForming(const TAngle& _anAngle) const
{
  bool bRetVal = false;
  
  if(patAtom1 == _anAngle.patAtom1 && patAtom2 == _anAngle.patAtom2 && patAtom3 == _anAngle.patAtom3);
  else if((patAtom1 == _anAngle.patAtom2 && patAtom2 == _anAngle.patAtom1) ||
	  (patAtom1 == _anAngle.patAtom2 && patAtom2 == _anAngle.patAtom3) ||
	  (patAtom2 == _anAngle.patAtom1 && patAtom3 == _anAngle.patAtom2) ||
	  (patAtom2 == _anAngle.patAtom3 && patAtom3 == _anAngle.patAtom2)  ) bRetVal = true;
  
  return bRetVal;
}

//---------------------------------------------------------------------------------------------------
void TAngle::CalcGradients()
{
  TVector vecDirect21 = vecMakeVector(patAtom2, patAtom1),
	  vecDirect23 = vecMakeVector(patAtom2, patAtom3);
	  
  double dLen21 = vecDirect21.dLength(),
	 dLen23 = vecDirect23.dLength();
	 
  vecDirect21 = vecDirect21/dLen21;
  vecDirect23 = vecDirect23/dLen23;

  double  dCosTheta = vecDirect21*vecDirect23;
  if(     dCosTheta > 1. ) dCosTheta = 1.;
  else if(dCosTheta < -1.) dCosTheta = -1.;
  
  double dEnergyDepth = 2.*apAngleInfo.dKTheta*(acos(dCosTheta) - apAngleInfo.dThetaEq*TConsts::PI/180.),
	 dSinTheta = (vecDirect21.vecCrossProd(vecDirect23)).dLength();
	 
  if(fabs(dSinTheta) < 0.00001 || TConsts::bIsNan(dSinTheta)) return;	// Protection against numerical instabilities.
	 
  TVector vecGrad1 = ((vecDirect21*vecDirect23)*vecDirect21 - vecDirect23)/dLen21,
	  vecGrad3 = ((vecDirect21*vecDirect23)*vecDirect23 - vecDirect21)/dLen23;
	  
	  vecGrad1 *= dEnergyDepth/dSinTheta;
	  vecGrad3 *= dEnergyDepth/dSinTheta;

  patAtom1->vecEnergyGradient += vecGrad1;
  patAtom3->vecEnergyGradient += vecGrad3;
  patAtom2->vecEnergyGradient -= vecGrad1 + vecGrad3;		// Newton's 3rd law ;)
}

//---------------------------------------------------------------------------------------------------
