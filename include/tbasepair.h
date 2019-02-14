#pragma once
#ifndef  BASEPAIR_H
#define  BASEPAIR_H

#include <list>
#include <vector>

#include <cassert>
#include <cmath>

#include "tatom.h"
#include "tconsts.h"
#include "thbond.h"
#include "textutils.h"
#include "tspring.h"
#include "tvector.h"
#include "params.h"

class TBasePair
{
  public:
    static const double dThreshold;	// - When possible base pair is a true base pair? 
				        // - When its energy is < dThreshold!
 
    TBasePair(THBond* _phbHBond1, THBond* _phbHBond2, THBond* _phbHBond3 = NULL);
    TBasePair(const std::vector<TAtom*>& _AtomPtrs);
    TBasePair(): AtomPtrs(), BaseSprings() {}

    double dCalcEnergy() const;
    void   CalcGradients();
    
  private:
    std::vector<TAtom*> AtomPtrs;
    std::list<TSpring> BaseSprings;	// It cannot be a vector, since operator= would involve assignment of const pointers
    double Xm, Ym, Zm;			// Coords of center of mass
    double a, b;
    unsigned int n;			// A constant, equal to AtomPtrs.size()
    
    double Fit3DPlane();
    double TotalSQRDistance(double _a, double _b) const;
    double TotalSQRDistance() const;
    void CalcCenterOfMass();
    
    static const double dK;
};

#endif /*BASEPAIR_H*/