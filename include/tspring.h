#pragma once
#ifndef  TSPRING_H
#define  TSPRING_H

#include "tnonbond.h"
#include "params.h"
#include "tvector.h"

class TSpring: public TNonbond
{
  private:
    BondParam bpBondInfo;

  public:
    TSpring(TAtom& _atAtom1, TAtom& _atAtom2, BondParam _bpBondInfo): 
      TNonbond(&_atAtom1, &_atAtom2), bpBondInfo(_bpBondInfo){}
    TSpring(TAtom* _patAtom1, TAtom* _patAtom2, BondParam _bpBondInfo): 
      TNonbond(_patAtom1, _patAtom2), bpBondInfo(_bpBondInfo) {}
    TSpring(): TNonbond(), bpBondInfo() {}

    double dCalcEnergy() const;
    void CalcGradients() const;
};

#endif /*TSPRING_H*/