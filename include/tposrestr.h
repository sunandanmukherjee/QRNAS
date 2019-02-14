#pragma once
#ifndef  TPOSRESTR_H
#define  TPOSRESTR_H

#include <cmath>

#include "tatom.h"

class TPosRestr
{
  public:
    TAtom * const patAtom1;
    TPosRestr(TAtom* _patAtom1, double _dX0, double _dY0, double _dZ0, double _dKr):
      patAtom1(_patAtom1), dKr(_dKr), dX0(_dX0), dY0(_dY0), dZ0(_dZ0) {}
    TPosRestr(TAtom& _atAtom1, double _dX0, double _dY0, double _dZ0, double _dKr):
      patAtom1(&_atAtom1), dKr(_dKr), dX0(_dX0), dY0(_dY0), dZ0(_dZ0) {}    

    double dLength()     const;
    double dLength_SQR() const;

    double dCalcEnergy()   const;
    void   CalcGradients() const;
    
  private:
    const double dKr;
    const double dX0, dY0, dZ0;
};

#endif /*TPOSRESTR_H*/