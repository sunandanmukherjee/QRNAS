#include "tdihed.h"

//---------------------------------------------------------------------------------------------------
double TDihedral::dValue() const
{
  TVector vec12 = vecMakeVector(patAtom1, patAtom2), 	
	  vec23 = vecMakeVector(patAtom2, patAtom3), 
	  vec34 = vecMakeVector(patAtom3, patAtom4);

  TVector vecNormal123 = vec12.vecCrossProd(vec23),
	  vecNormal234 = vec23.vecCrossProd(vec34);	  
  vecNormal123.Normalize();
  vecNormal234.Normalize();
  
  double  dCosPhi = vecNormal123*vecNormal234;
  if(     dCosPhi >  1.) dCosPhi =  1.;
  else if(dCosPhi < -1.) dCosPhi = -1.;
  
  return acos(dCosPhi) * (vec23*(vecNormal123.vecCrossProd(vecNormal234)) >= 0. ? 1.:-1.);
}

//---------------------------------------------------------------------------------------------------
double TDihedral::dCalcEnergy() const
{
  return dpDihedralInfo.dKI/double(dpDihedralInfo.iDiv)*(1.+cos(dpDihedralInfo.dPN*dValue()
	    -dpDihedralInfo.dPhase*TConsts::PI/180.));
}

//---------------------------------------------------------------------------------------------------
void TDihedral::CalcGradients()
{
  TVector vec12 = vecMakeVector(patAtom1, patAtom2), 	
	  vec23 = vecMakeVector(patAtom2, patAtom3), 
	  vec34 = vecMakeVector(patAtom3, patAtom4);

  TVector vecNormal123 = vec12.vecCrossProd(vec23),
	  vecNormal234 = vec23.vecCrossProd(vec34);	   
	  
  double dLen123 = vecNormal123.dLength(),
	 dLen234 = vecNormal234.dLength();
	 
  vecNormal123 /= dLen123;	// Normalized, without an explicit call of TVector::Normalize()
  vecNormal234 /= dLen234;	// for questionable speed gain.
  
  double dPhi = acos(vecNormal123*vecNormal234) * (vec23*(vecNormal123.vecCrossProd(vecNormal234)) >= 0. ? 1.:-1.);
  double dSinPhi = sin(dPhi);
  
  if(fabs(dSinPhi) < 0.00001 || TConsts::bIsNan(dSinPhi)) return;
    
 double dEnergyDepth = dpDihedralInfo.dKI*dpDihedralInfo.dPN*sin(dPhi*dpDihedralInfo.dPN-dpDihedralInfo.dPhase*TConsts::PI/180.)/
			  (double(dpDihedralInfo.iDiv)*dSinPhi);

/***************************************************************************************************
 * Now follows the preliminary version of formulae. Due to (theoretical) instability at dPhi=0
 * another approach was employed below (but unused since the first one is fine).
 **************************************************************************************************/
  TVector a = (-vecNormal234 + vecNormal123*cos(dPhi))/dLen123,
	  b = (-vecNormal123 + vecNormal234*cos(dPhi))/dLen234;

  patAtom1->vecEnergyGradient += dEnergyDepth*(vec23.vecCrossProd(a));
  patAtom2->vecEnergyGradient += dEnergyDepth*(a.vecCrossProd(vec12 + vec23) + vec34.vecCrossProd(b));
  patAtom3->vecEnergyGradient += dEnergyDepth*(b.vecCrossProd(vec23 + vec34) + vec12.vecCrossProd(a));
  patAtom4->vecEnergyGradient += dEnergyDepth*(vec23.vecCrossProd(b));
  return;

/*
  double dLen23  = vec23.dLength();
  patAtom1->vecEnergyGradient -= dEnergyDepth*(dLen23/dLen123*vecNormal123);
  patAtom2->vecEnergyGradient += dEnergyDepth*((vec12*vec23/(dLen23*dLen23)+1.)*dLen23*vecNormal123/dLen123 + vec23*vec34/(dLen23*dLen234)*vecNormal234);
  patAtom3->vecEnergyGradient -= dEnergyDepth*((vec23*vec34/(dLen23*dLen23)+1.)*dLen23*vecNormal234/dLen234 + vec12*vec23/(dLen23*dLen123)*vecNormal123);
  patAtom4->vecEnergyGradient += dEnergyDepth*(dLen23/dLen234*vecNormal234);
*/
}
//---------------------------------------------------------------------------------------------------
    