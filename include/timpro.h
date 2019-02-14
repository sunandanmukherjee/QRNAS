#pragma once
#ifndef  TIMPROPER_H
#define  TIMPROPER_H

#include "tatom.h"
#include "params.h"

class TImproper
{
  public:
    TAtom * const patAtom1;
    TAtom * const patAtom2;
    TAtom * const patAtom3;
    TAtom * const patAtom4;

    TImproper(TAtom& _atAtom1, TAtom& _atAtom2, TAtom& _atAtom3, TAtom& _atAtom4, ImproperParam _ipImproperInfo):
      patAtom1(&_atAtom1), patAtom2(&_atAtom2), patAtom3(&_atAtom3), patAtom4(&_atAtom4), ipImproperInfo(_ipImproperInfo) {}    
    TImproper(TAtom* _patAtom1, TAtom* _patAtom2, TAtom* _patAtom3, TAtom* _patAtom4, ImproperParam _ipImproperInfo):
      patAtom1(_patAtom1), patAtom2(_patAtom2), patAtom3(_patAtom3), patAtom4(_patAtom4), ipImproperInfo(_ipImproperInfo) {}
    TImproper(const TImproper& _imImproper): patAtom1(_imImproper.patAtom1), patAtom2(_imImproper.patAtom2),
      patAtom3(_imImproper.patAtom3), patAtom4(_imImproper.patAtom4), ipImproperInfo(_imImproper.ipImproperInfo) {}
    TImproper(): patAtom1(NULL), patAtom2(NULL), patAtom3(NULL), patAtom4(NULL), ipImproperInfo() {}

    double dValue() const;
    double dCalcEnergy() const;
    void CalcGradients();

  private:
    TImproper& operator=(const TImproper&);
    ImproperParam ipImproperInfo;				
};

#endif /*TIMPROPER_H*/