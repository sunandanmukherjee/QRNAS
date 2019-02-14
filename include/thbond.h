#pragma once
#ifndef  HBOND_H
#define  HBOND_H

#include "tatom.h"
#include "tvector.h"

class THBond
{
  friend class TBasePair;
  public:
    TAtom * const patParent;		// Parent atom (to whom hydrogen binds)
    TAtom * const patDonor;		// Hydrogen
    TAtom * const patAcceptor;		// Nitrogen or oxygen

    static const double dThreshold;
    
    THBond(TAtom& _patParent, TAtom& _patDonor, TAtom& _patAcceptor): 
      patParent(&_patParent), patDonor(&_patDonor), patAcceptor(&_patAcceptor) {}
    THBond(TAtom* _patParent, TAtom* _patDonor, TAtom* _patAcceptor): 
      patParent(_patParent), patDonor(_patDonor), patAcceptor(_patAcceptor) {}
    THBond(): patParent(NULL), patDonor(NULL), patAcceptor(NULL) {}

    double dCalcEnergy() const;
    void   CalcGradients() const;
    
    double dCosTheta()   const;
    double dLength()     const;
    double dLength_SQR() const;
  private:
    static const double dREq, dK, dXiR, dXiTheta;

//    THBond& operator=(const THBond&);
};

#endif /*HBOND_H*/