#pragma once
#ifndef  TBOND_H
#define  TBOND_H

#include <cmath>

#include "tatom.h"
#include "tconsts.h"
#include "textutils.h"
#include "params.h"

class TBond
{
  public:
    TAtom * const patAtom1;
    TAtom * const patAtom2;

    TBond(TAtom& _atAtom1, TAtom& _atAtom2, BondParam _bpBondInfo): patAtom1(&_atAtom1), patAtom2(&_atAtom2), 
      bpBondInfo(_bpBondInfo) {}
    TBond(TAtom* _patAtom1, TAtom* _patAtom2, BondParam _bpBondInfo): patAtom1(_patAtom1), patAtom2(_patAtom2), 
      bpBondInfo(_bpBondInfo) {}
    TBond(const TBond& _boBond): patAtom1(_boBond.patAtom1), patAtom2(_boBond.patAtom2), 
      bpBondInfo(_boBond.bpBondInfo) {}
    TBond(): patAtom1(NULL), patAtom2(NULL), bpBondInfo() {}

    double dLength_SQR() const;
    double dLength() const;
    double dCalcEnergy() const;
    void CalcGradients();

  private:
    TBond& operator=(const TBond&);
    BondParam bpBondInfo;
				
};

#endif /*TBOND_H*/