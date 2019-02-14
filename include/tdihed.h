#pragma once
#ifndef  TDIHEDRAL_H
#define  TDIHEDRAL_H

#include <cmath>

#include "tatom.h"
#include "tvector.h"
#include "params.h"

class TDihedral
{
  public:
    TAtom * const patAtom1;
    TAtom * const patAtom2;
    TAtom * const patAtom3;
    TAtom * const patAtom4;

    TDihedral(TAtom& _atAtom1, TAtom& _atAtom2, TAtom& _atAtom3, TAtom& _atAtom4, DihedralParam _dpDihedralInfo):
      patAtom1(&_atAtom1), patAtom2(&_atAtom2), patAtom3(&_atAtom3), patAtom4(&_atAtom4), dpDihedralInfo(_dpDihedralInfo) {}
    TDihedral(TAtom* _patAtom1, TAtom* _patAtom2, TAtom* _patAtom3, TAtom* _patAtom4, DihedralParam _dpDihedralInfo):
      patAtom1(_patAtom1), patAtom2(_patAtom2), patAtom3(_patAtom3), patAtom4(_patAtom4), dpDihedralInfo(_dpDihedralInfo) {}
    TDihedral(const TDihedral& _diDihedral): patAtom1(_diDihedral.patAtom1), patAtom2(_diDihedral.patAtom2), 
      patAtom3(_diDihedral.patAtom3), patAtom4(_diDihedral.patAtom4), dpDihedralInfo(_diDihedral.dpDihedralInfo) {}
    TDihedral():
      patAtom1(NULL), patAtom2(NULL), patAtom3(NULL), patAtom4(NULL), dpDihedralInfo() {}
    
    double dValue() const;
    double dCalcEnergy() const;
    void    CalcGradients();    

  private:
    TDihedral& operator=(const TDihedral&);
    DihedralParam dpDihedralInfo;				
};

#endif /*TDIHEDRAL_H*/