#include "tlennard.h"

//---------------------------------------------------------------------------------------------------
void TLennard::Scale14()
{
  LJInfo.dLJ6  /= 2.;
  LJInfo.dLJ12 /= 2.;
}

//---------------------------------------------------------------------------------------------------
double TLennard::dCalcEnergy() const
{
  double dLen_SQR = dLength_SQR();
  
  if(dLen_SQR > TConsts::CUTOFF_SQR) return 0.;
  
/**************************************************************************
 * Naiive cutoff for Lennard-Jones follows. Needs further improvement.
 **************************************************************************/
   return LJInfo.dLJ6 /(dLen_SQR*dLen_SQR*dLen_SQR) + 
          LJInfo.dLJ12/(dLen_SQR*dLen_SQR*dLen_SQR*dLen_SQR*dLen_SQR*dLen_SQR);
}

//---------------------------------------------------------------------------------------------------
void TLennard::CalcGradients() const
{
  TVector vec12 = vecMakeVector(patAtom1, patAtom2);
  double dLen_SQR = vec12.dLength_SQR();

  if(dLen_SQR > TConsts::CUTOFF_SQR) return;
  
  vec12 *= (6. *LJInfo.dLJ6 /(dLen_SQR*dLen_SQR*dLen_SQR*dLen_SQR) + 
   	    12.*LJInfo.dLJ12/(dLen_SQR*dLen_SQR*dLen_SQR*dLen_SQR*dLen_SQR*dLen_SQR*dLen_SQR));
  patAtom1->vecEnergyGradient += vec12;
  patAtom2->vecEnergyGradient -= vec12;
}

//---------------------------------------------------------------------------------------------------
TLennard ljMakeLJPair(TAtom& _atAtomA, TAtom& _atAtomB)
{
 /**************************************************************************************************************
  * Object of type vdWInfo contains the properties of single *atom*. Therefore, to get LJ6 and LJ12 parameters
  * we must do *pairwise* Lorenz-Berthelot calculation.
  **************************************************************************************************************/  
  double dRStar   = (_atAtomA.vdWInfo.dR0 + _atAtomB.vdWInfo.dR0),		 // Lorenz-Berthelot
	 dEpsilonAB = sqrt(_atAtomA.vdWInfo.dEpsilon*_atAtomB.vdWInfo.dEpsilon); // I hope its correct
  double dLJ6  = -2.*dEpsilonAB*pow(dRStar,  6),
	 dLJ12 =     dEpsilonAB*pow(dRStar, 12);
  return TLennard(_atAtomA, _atAtomB, LJParam(dLJ6, dLJ12));
}
//---------------------------------------------------------------------------------------------------
