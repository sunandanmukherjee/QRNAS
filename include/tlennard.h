#pragma once
#ifndef  TLENNARD_H
#define  TLENNARD_H

#include "tnonbond.h"
#include "tvector.h"
#include "params.h"

class TLennard: public TNonbond
{
  private:
    LJParam LJInfo;
  
  public:
    TLennard(TAtom& _atAtom1, TAtom& _atAtom2, LJParam _LJInfo = LJParam()): 
      TNonbond(&_atAtom1, &_atAtom2), LJInfo(_LJInfo) {}
    TLennard(TAtom* _patAtom1, TAtom* _patAtom2, LJParam _LJInfo = LJParam()): 
      TNonbond(_patAtom1, _patAtom2), LJInfo(_LJInfo) {}
    TLennard(): TNonbond(), LJInfo() {}

    void Scale14();		// Scale interaction down when between atoms separated by 3 bonds

    double dCalcEnergy() const;
    void   CalcGradients() const;
};

#endif /*TLENNARD_H*/