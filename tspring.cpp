#include "tspring.h"

//---------------------------------------------------------------------------------------------------
double TSpring::dCalcEnergy() const
{
  double dEfficientDistance = dLength() - bpBondInfo.dREq;
  return bpBondInfo.dKr*dEfficientDistance*dEfficientDistance;
}

//---------------------------------------------------------------------------------------------------
void TSpring::CalcGradients() const
{
  TVector vecDirect12;
  vecDirect12 = vecMakeVector(patAtom1, patAtom2);
  double dBondLength = vecDirect12.dLength();
  vecDirect12.Normalize();
  
  patAtom1->vecEnergyGradient -= (2.*(bpBondInfo.dKr)*(dBondLength - bpBondInfo.dREq)) * vecDirect12;
  patAtom2->vecEnergyGradient += (2.*(bpBondInfo.dKr)*(dBondLength - bpBondInfo.dREq)) * vecDirect12;
}

//---------------------------------------------------------------------------------------------------
