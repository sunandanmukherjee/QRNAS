#pragma once
#ifndef  TVECTOR_H
#define  TVECTOR_H

#include <cmath>

#include "tconsts.h"

class TAtom;

class TVector
{
  protected:
    double dX, dY, dZ;
  public:
    TVector(const double _dX=0., const double _dY=0., const double _dZ=0.): dX(_dX), dY(_dY), dZ(_dZ) {}

    friend double operator*(const TVector& _tvVec1, const TVector& _tvVec2);	//scalar product
    friend TVector operator*(const double _dMul, const TVector& _tvVec);
    friend TAtom& operator+=(TAtom& _atAtom, const TVector& _vecVector);
    friend TAtom& operator-=(TAtom& _atAtom, const TVector& _vecVector);
    
    TVector vecCrossProd(const TVector& _tvVec) const;	//v.vecCrossProd(u) returns (v x u)
    TVector operator-() const;			//unary minus
    TVector operator+(const TVector& _tvVec) const;
    TVector operator-(const TVector& _tvVec) const;	//binary minus
    TVector operator*(const double& _dMul) const;
    TVector operator/(const double& _dDiv) const;
    TVector& operator+=(const TVector& _tvVec);
    TVector& operator-=(const TVector& _tvVec);
    TVector& operator*=(const double& _dMul);
    TVector& operator/=(const double& _dDiv);
    
    double dLength() const;
    double dLength_SQR() const;	//such that: length()*length() = length_SQR()
    
    void Normalize();
    
};

#endif /*TVECTOR_H*/