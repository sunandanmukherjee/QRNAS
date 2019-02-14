#include "tbond.h"

//---------------------------------------------------------------------------------------------------
double TBond::dLength_SQR() const
{
  return dDistance_SQR(patAtom1, patAtom2);
}

//---------------------------------------------------------------------------------------------------
double TBond::dLength() const
{
  return sqrt(dLength_SQR());
}

//---------------------------------------------------------------------------------------------------
double TBond::dCalcEnergy() const
{
  double dEfficientDistance = dLength() - bpBondInfo.dREq;
  return bpBondInfo.dKr*dEfficientDistance*dEfficientDistance;
}

//---------------------------------------------------------------------------------------------------
void TBond::CalcGradients()
{
  TVector vecDirect12 = vecMakeVector(patAtom1, patAtom2);
  double dBondLength = vecDirect12.dLength();
  vecDirect12.Normalize();
  
  patAtom1->vecEnergyGradient -= (2.*(bpBondInfo.dKr)*(dBondLength - bpBondInfo.dREq)) * vecDirect12;
  patAtom2->vecEnergyGradient += (2.*(bpBondInfo.dKr)*(dBondLength - bpBondInfo.dREq)) * vecDirect12;
}

//---------------------------------------------------------------------------------------------------
