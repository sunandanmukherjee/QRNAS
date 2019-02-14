#pragma once
#ifndef  TANGLE_H
#define  TANGLE_H

#include <cmath>

#include "tatom.h"
#include "tbond.h"
#include "tconsts.h"
#include "textutils.h"
#include "tvector.h"
#include "params.h"

class TAngle
{
  public:
    TAtom * const patAtom1;
    TAtom * const patAtom2;
    TAtom * const patAtom3;

    TAngle(TAtom& _atAtom1, TAtom& _atAtom2, TAtom& _atAtom3, AngleParam _apAngleInfo): 
      patAtom1(&_atAtom1), patAtom2(&_atAtom2), patAtom3(&_atAtom3), apAngleInfo(_apAngleInfo) {}
    TAngle(TAtom* _patAtom1, TAtom* _patAtom2, TAtom* _patAtom3, AngleParam _apAngleInfo):
      patAtom1(_patAtom1), patAtom2(_patAtom2), patAtom3(_patAtom3), apAngleInfo(_apAngleInfo) {}
    TAngle(const TAngle& _anAngle): patAtom1(_anAngle.patAtom1), 
      patAtom2(_anAngle.patAtom2), patAtom3(_anAngle.patAtom3), apAngleInfo(_anAngle.apAngleInfo) {}
    TAngle(): patAtom1(NULL), patAtom2(NULL), patAtom3(NULL), apAngleInfo() {}

    double dValue() const;
    double dCalcEnergy() const;
    bool bContains(const TBond& _boBond) const;
    bool bDihedralForming(const TAngle& _anAngle) const;	//returns false if angles are isomorphic; otherwise tells whether angle forms dihedral with the given one
    void CalcGradients();

  private:
    TAngle& operator=(const TAngle&);
    AngleParam apAngleInfo;				
};

#endif /*TANGLE_H*/