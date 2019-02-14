#include "timpro.h"

//---------------------------------------------------------------------------------------------------
double TImproper::dValue() const
{
  TVector vec12 = vecMakeVector(patAtom3, patAtom1), 	
	  vec23 = vecMakeVector(patAtom1, patAtom2), 
	  vec34 = vecMakeVector(patAtom2, patAtom4);

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
double TImproper::dCalcEnergy() const
{
  return ipImproperInfo.dKI*(1.+cos(ipImproperInfo.dPN*dValue()
	    -ipImproperInfo.dPhase*TConsts::PI/180.));
}

//---------------------------------------------------------------------------------------------------
void TImproper::CalcGradients()
{
  TVector vec12 = vecMakeVector(patAtom3, patAtom1), 	
	  vec23 = vecMakeVector(patAtom1, patAtom2), 
	  vec34 = vecMakeVector(patAtom2, patAtom4);

  TVector vecNormal123 = vec12.vecCrossProd(vec23),
	  vecNormal234 = vec23.vecCrossProd(vec34);	   
	  
  double dLen123 = vecNormal123.dLength(),	//consider +SMALLVAL here
	 dLen234 = vecNormal234.dLength();
	 
  vecNormal123 /= dLen123;	// Normalized, without an explicit call of TVector::Normalize()
  vecNormal234 /= dLen234;	// for questionable speed gain
  
  double dPhi = acos(vecNormal123*vecNormal234) * (vec23*(vecNormal123.vecCrossProd(vecNormal234)) >= 0. ? 1.:-1.);
  double dSinPhi = sin(dPhi);
  
  if(fabs(dSinPhi) < 0.00001 || TConsts::bIsNan(dSinPhi)) return;
    
  double dEnergyDepth = ipImproperInfo.dKI*ipImproperInfo.dPN*sin(dPhi*ipImproperInfo.dPN-
			  ipImproperInfo.dPhase*TConsts::PI/180.)/dSinPhi;

/***************************************************************************************************
 * The 'worse' method for calcn of gradients.
 **************************************************************************************************/

  TVector a = (-vecNormal234 + vecNormal123*cos(dPhi))/dLen123,
	  b = (-vecNormal123 + vecNormal234*cos(dPhi))/dLen234;

  patAtom3->vecEnergyGradient += dEnergyDepth*(vec23.vecCrossProd(a));
  patAtom2->vecEnergyGradient += dEnergyDepth*(a.vecCrossProd(vec12 + vec23) + vec34.vecCrossProd(b));
  patAtom1->vecEnergyGradient += dEnergyDepth*(b.vecCrossProd(vec23 + vec34) + vec12.vecCrossProd(a));
  patAtom4->vecEnergyGradient += dEnergyDepth*(vec23.vecCrossProd(b));
}
    
//---------------------------------------------------------------------------------------------------
    