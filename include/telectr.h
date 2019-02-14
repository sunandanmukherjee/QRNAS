#pragma once
#ifndef  TELECTR_H
#define  TELECTR_H

#include <cmath>

#include "tnonbond.h"
#include "tvector.h"
#include "tconsts.h"

class TElectr: public TNonbond
{
  private:
    double dQQ;					// Value such that: Energy = dQQ/r
    double dEffectiveBornDistance() const;	// Calculate Born distance (the 'f_AB' parameter)

  public:
    TElectr(TAtom& _atAtom1, TAtom& _atAtom2, double _dQQ =-1.): 
      TNonbond(&_atAtom1, &_atAtom2), dQQ(_dQQ) {}
    TElectr(TAtom* _patAtom1, TAtom* _patAtom2, double _dQQ =-1.):
      TNonbond(_patAtom1, _patAtom2), dQQ(_dQQ) {}
    TElectr(): TNonbond(), dQQ(-1.) {}
    
    void Scale14();		// Scale interaction down when between atoms separated by 3 bonds
    
    double dCalcEnergy() const;
    void   CalcGradients() const; // Adds its contribution to TAtoms' vecEnergyGradient-s
};

#endif /*TELECTR_H*/