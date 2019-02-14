#include "telectr.h"

//---------------------------------------------------------------------------------------------------
void TElectr::Scale14()
{
  dQQ /= 1.2;
}

//---------------------------------------------------------------------------------------------------
double TElectr::dCalcEnergy() const
{
// Some polarizable medium. Does influence energy, but not gradients:
  if(TConsts::USEBORN)
  {
    return dQQ/(sqrt(dLength_SQR())*TConsts::EpsilonRNA) +	
	      -dQQ*(1./TConsts::EpsilonRNA - 1./TConsts::EpsilonWater)/dEffectiveBornDistance();
  }
  else return dQQ/(sqrt(dLength_SQR())*TConsts::EpsilonRNA);	// Simplified approach.
}

//---------------------------------------------------------------------------------------------------
void TElectr::CalcGradients() const
{
  TVector vec12 = vecMakeVector(patAtom1, patAtom2);
  double dLen_SQR = vec12.dLength_SQR();
	 
  vec12 *= dQQ/(dLen_SQR*sqrt(dLen_SQR)*TConsts::EpsilonRNA);
  
  patAtom1->vecEnergyGradient += vec12;
  patAtom2->vecEnergyGradient -= vec12;
}

//---------------------------------------------------------------------------------------------------
double TElectr::dEffectiveBornDistance() const
{
  double dLen_SQR = dLength_SQR(),
	 R1R2	  = (patAtom1->dBornRadius)*(patAtom2->dBornRadius);
  return sqrt(dLen_SQR + R1R2*exp(-dLen_SQR/(4.*R1R2)));
}

//---------------------------------------------------------------------------------------------------
TElectr elMakeElectrPair(TAtom& _atAtomA, TAtom& _atAtomB)
{
  /**************************************************************************************************************
  * It would be convenient to use some simple form of electrostatic energy, like const/r, so we define dQQ
  * as the product of two interacting charges times Coulomb constant.
  **************************************************************************************************************/  
  return TElectr(_atAtomA, _atAtomB, TConsts::COULOMB*_atAtomA.dCharge*_atAtomB.dCharge);
}
    
//---------------------------------------------------------------------------------------------------

