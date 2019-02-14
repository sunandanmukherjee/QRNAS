#pragma once
#ifndef  TNONBOND_H
#define  TNONBOND_H

#include <cmath>

#include "tatom.h"

class TNonbond
{
  public:
    TAtom * const patAtom1;
    TAtom * const patAtom2;

    TNonbond(TAtom* _patAtom1, TAtom* _patAtom2): patAtom1(_patAtom1), patAtom2(_patAtom2) {}
    TNonbond(TAtom& _atAtom1, TAtom& _atAtom2): patAtom1(&_atAtom1), patAtom2(&_atAtom2) {}
    TNonbond(): patAtom1(NULL), patAtom2(NULL) {}
    virtual double dCalcEnergy()   const =0;
    virtual void   CalcGradients() const =0;

    double dLength()     const;
    double dLength_SQR() const;
};

#endif /*TNONBOND_H*/