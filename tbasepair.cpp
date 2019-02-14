#include "tbasepair.h"

//---------------------------------------------------------------------------------------------------
const double TBasePair::dK = 5.;		//[kcal/mol/Ang^2] Force constant for off-plane distances
const double TBasePair::dThreshold = 2.;	//[kcal/mol] Energy threshold - when a base pair is a base pair?
//---------------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------
TBasePair::TBasePair(THBond* _phbHBond1, THBond* _phbHBond2, THBond* _phbHBond3)
{
  if(_phbHBond1 != _phbHBond2 && _phbHBond2 != _phbHBond3 && _phbHBond1 != _phbHBond3)
  {
    if(_phbHBond1 != NULL)
    {
      AtomPtrs.push_back(_phbHBond1->patParent);
      AtomPtrs.push_back(_phbHBond1->patDonor);
      AtomPtrs.push_back(_phbHBond1->patAcceptor);      
    }
    if(_phbHBond2 != NULL)
    {
      AtomPtrs.push_back(_phbHBond2->patParent);
      AtomPtrs.push_back(_phbHBond2->patDonor);
      AtomPtrs.push_back(_phbHBond2->patAcceptor);      
    }
    if(_phbHBond3 != NULL)
    {
      AtomPtrs.push_back(_phbHBond3->patParent);
      AtomPtrs.push_back(_phbHBond3->patDonor);
      AtomPtrs.push_back(_phbHBond3->patAcceptor);      
    }
    
    if(!(n = AtomPtrs.size())) 
      std::cerr << "Internal Error: Cannot build the pair of bases." << std::endl;
    
    Fit3DPlane();
  }
  else std::cerr << "Internal Warning: Attempt to build a trivial pair of bases." << std::endl;
}

//---------------------------------------------------------------------------------------------------
TBasePair::TBasePair(const std::vector<TAtom*> &_AtomPtrs)
{
  if(!_AtomPtrs.empty())
  {
    AtomPtrs = _AtomPtrs;
    n = AtomPtrs.size();
    Fit3DPlane();    

    for(unsigned int i = 0; i < n; ++i)
      for(unsigned int j = i; j < n; ++j)
      {
	if(AtomPtrs[i]->sGetResIdent() != AtomPtrs[j]->sGetResIdent() &&
	    Textutils::bIsWCHydrogenBond(AtomPtrs[i]->sGetName(), AtomPtrs[i]->sGetResName(), 
					 AtomPtrs[j]->sGetName(), AtomPtrs[j]->sGetResName()) )
	{
	  BaseSprings.push_back(TSpring(AtomPtrs[i], AtomPtrs[j], BondParam(10., 1.65)));
	}
      }
  }
  else std::cerr << "Internal Warning: Attempt to build an empty pair of bases." << std::endl;
}
//---------------------------------------------------------------------------------------------------
void TBasePair::CalcCenterOfMass()
{// It's not precisely the center of mass (because atoms have no masses) but rather geometrical center.
  Xm = 0., Ym = 0., Zm = 0.;
  for(unsigned int i = 0; i < n; ++i)
  {
    double Xk, Yk, Zk;
    AtomPtrs[i]->GetCoords(Xk, Yk, Zk);
    Xm += Xk;
    Ym += Yk;
    Zm += Zk;
  }
  Xm /= double(n);	// Mean coordinates
  Ym /= double(n);
  Zm /= double(n);  
}

//---------------------------------------------------------------------------------------------------
double TBasePair::Fit3DPlane()
{/******************************************************************************************************
 * Returns the same as TotalSQRDistance(), but for optimal {a, b} (also calculates the optimal {a, b}
 * and stores them as this->a and this->b). It's a solution to 3rd degree polynomial.
 * Save your time and don't stare at this for too long.
 ******************************************************************************************************/
  CalcCenterOfMass();
  
  double Sxx = 0., Syy = 0., Szz = 0., Sxy = 0., Sxz = 0., Syz = 0.;
  const unsigned int n = AtomPtrs.size();  
  for(unsigned int i = 0; i < n; ++i)
  {
    double xk, yk, zk;
    AtomPtrs[i]->GetCoords(xk, yk, zk);
    xk -= Xm;		// Reduced coords
    yk -= Ym;
    zk -= Zm;
    
    Sxx += xk*xk;
    Syy += yk*yk;
    Szz += zk*zk;
    Sxy += xk*yk;
    Sxz += xk*zk;
    Syz += yk*zk;
  }
  
  const double c0 = Syz*(Sxy*Sxy - Sxz*Sxz) + Sxy*Sxz*(Szz - Syy),
	       c1 = Sxy*Sxy*Sxy + Sxy*(Sxz*Sxz - 2.*Syz*Syz - Szz*Szz) + Sxy*(Sxx*Szz + Syy*Szz - Sxx*Syy) + 
		    Sxz*Syz*(Syy + Szz - 2.*Sxx),
	       c2 = Syz*Syz*Syz + Syz*(Sxz*Sxz - 2.*Sxy*Sxy - Sxx*Sxx) + Syz*(Sxx*Szz + Sxx*Syy - Syy*Szz) +
		    Sxy*Sxz*(Sxx + Syy - 2.*Szz),
	       c3 = Sxy*(Syz*Syz - Sxz*Sxz) + Sxz*Syz*(Sxx - Syy);

  const double r = c2/c3,
	       s = c1/c3,
	       t = c0/c3;
  
  const double p = s - r*r/3.,
	       q = 2.*r*r*r/27. - r*s/3. + t,
	       R = q*q/4. + p*p*p/27.;

  double SumDelta;	// Sum over all AtomPtrs of squared distances to the plane
  
  if(R >= 0.)
  {
    a = -r/3. + pow(-q/2. + sqrt(R), 0.333333333) + pow(-q/2. - sqrt(R), 0.333333333);
    b = (Sxy*Syz*a*a + (Syz*Syz - Sxy*Sxy)*a - Sxy*Syz)/
	((Syz*(Sxx - Syy) - Sxy*Sxz)*a + Sxy*(Syy - Szz) + Sxz*Syz);
    SumDelta = TotalSQRDistance(a, b);
  }
  else
  {
    const double rho = sqrt(-p*p*p/27.),
		 phi = acos(-q/(2.*rho));
    
    const double a1 = -r/3. + 2.*pow(rho, 0.333333333)*cos(phi/3.),
		 a2 = -r/3. + 2.*pow(rho, 0.333333333)*cos((phi + 2.*TConsts::PI)/3.),
		 a3 = -r/3. + 2.*pow(rho, 0.333333333)*cos((phi + 4.*TConsts::PI)/3.);
		 
    const double b1 = (Sxy*Syz*a1*a1 + (Syz*Syz - Sxy*Sxy)*a1 - Sxy*Syz)/
		      ((Syz*(Sxx - Syy) - Sxy*Sxz)*a1 + Sxy*(Syy - Szz) + Sxz*Syz),
		 b2 = (Sxy*Syz*a2*a2 + (Syz*Syz - Sxy*Sxy)*a2 - Sxy*Syz)/
		      ((Syz*(Sxx - Syy) - Sxy*Sxz)*a2 + Sxy*(Syy - Szz) + Sxz*Syz),
		 b3 = (Sxy*Syz*a3*a3 + (Syz*Syz - Sxy*Sxy)*a3 - Sxy*Syz)/
		      ((Syz*(Sxx - Syy) - Sxy*Sxz)*a3 + Sxy*(Syy - Szz) + Sxz*Syz);
		 
    const double SumDelta1 = TotalSQRDistance(a1, b1),
		 SumDelta2 = TotalSQRDistance(a2, b2),
		 SumDelta3 = TotalSQRDistance(a3, b3);
		 
    a = a1, b = b1; SumDelta = SumDelta1;
    if(SumDelta2 < SumDelta1) {a = a2; b = b2; SumDelta = SumDelta2;}
    if(SumDelta3 < SumDelta1 && SumDelta3 < SumDelta2) {a = a3; b = b3; SumDelta = SumDelta3;}
  } 

/*// If you need the explicit equation for the plane: A*x + B*y + C*z + 1 = 0
  double lambda = 0.;
  for(unsigned int i = 0; i < n; ++i)
  {
    double Xk, Yk, Zk;
    AtomPtrs[i]->GetCoords(Xk, Yk, Zk);
    lambda += a*Xk + b*Yk + Zk;
  }
  lambda /= double(n);
  
  double A = -a/lambda,
	 B = -b/lambda,
	 C = -1./lambda;*/
  
  return SumDelta;  
}

//---------------------------------------------------------------------------------------------------
double TBasePair::TotalSQRDistance(double _a, double _b) const
{// Returns sum of distances of all atoms from the plane given by numbers {a, b}, and containing the center of mass
  double SumDelta = 0.;
  for(unsigned int i = 0; i < n; ++i)
  {
    double xk, yk, zk;
    AtomPtrs[i]->GetCoords(xk, yk, zk);
    xk -= Xm;		// Reduced coords
    yk -= Ym;
    zk -= Zm;
    SumDelta += (_a*xk + _b*yk + zk)*(_a*xk + _b*yk + zk);
  }
  return SumDelta/(_a*_a + _b*_b + 1.);
}

//---------------------------------------------------------------------------------------------------
double TBasePair::TotalSQRDistance() const
{// Same thing as above, but overloaded without arguments (uses this->a and this->b instead of _a and _b)
  double SumDelta = 0.;
  for(unsigned int i = 0; i < n; ++i)
  {
    double xk, yk, zk;
    AtomPtrs[i]->GetCoords(xk, yk, zk);
    xk -= Xm;		// Reduced coords
    yk -= Ym;
    zk -= Zm;
    SumDelta += (a*xk + b*yk + zk)*(a*xk + b*yk + zk);
  }
  return SumDelta/(a*a + b*b + 1.);
}

//---------------------------------------------------------------------------------------------------
double TBasePair::dCalcEnergy() const
{
  double dEnergy = dK*TotalSQRDistance();
  for(std::list<TSpring>::const_iterator I = BaseSprings.begin(), end = BaseSprings.end(); I != end; ++I)
    dEnergy += I->dCalcEnergy();  
  
  return (TConsts::bIsNan(dEnergy)) ? 0. : dEnergy;
}

//---------------------------------------------------------------------------------------------------
void TBasePair::CalcGradients()
{
  Fit3DPlane();
  
  if(!TConsts::bIsNan(a) && !TConsts::bIsNan(b))
    for(unsigned int i = 0; i < AtomPtrs.size(); ++i)
    {
      double xk, yk, zk;
      AtomPtrs[i]->GetCoords(xk, yk, zk);
      xk -= Xm;		// Reduced coords
      yk -= Ym;
      zk -= Zm;
      AtomPtrs[i]->vecEnergyGradient += (2.*dK*(a*xk + b*yk + zk)/(a*a + b*b + 1.))*TVector(a, b, 1.);
    }
  
  for(std::list<TSpring>::const_iterator I = BaseSprings.begin(), end = BaseSprings.end(); I != end; ++I)
    I->CalcGradients();
}

//---------------------------------------------------------------------------------------------------
