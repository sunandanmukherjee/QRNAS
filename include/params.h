#pragma once
#ifndef  PARAMS_H
#define  PARAMS_H

struct vdWParam
{
  double dR0, dEpsilon;
  vdWParam(double _dR0 =-1., double _dEpsilon =-1.): dR0(_dR0), dEpsilon(_dEpsilon) {}
};

struct BondParam
{
  double dKr, dREq;
  BondParam(const double _dKr =-1., const double _dREq =-1.): dKr(_dKr), dREq(_dREq) {}
};

struct AngleParam
{
  double dKTheta, dThetaEq;
  AngleParam(const double _dKTheta =-1., const double _dThetaEq =-1.): dKTheta(_dKTheta), dThetaEq(_dThetaEq) {}
};

struct ImproperParam
{
  double dKI, dPhase, dPN;
  ImproperParam(const double _dKI =-1., const double _dPhase =-1., const double& _dPN =-1.): 
      dKI(_dKI), dPhase(_dPhase), dPN(_dPN) {}
};

struct DihedralParam
{
  int iDiv;
  double dKI, dPhase, dPN;
  DihedralParam(const int _iDiv =-1, const double _dKI =-1., const double _dPhase =-1., const double _dPN =-1.): 
      iDiv(_iDiv), dKI(_dKI), dPhase(_dPhase), dPN(_dPN) {}
};

struct LJParam
{
  double dLJ6, dLJ12;	// Values such that: Energy = dLJ6/r^6 + dLJ12/r^12
  LJParam(double _dLJ6 =-1., double _dLJ12 =-1.):  dLJ6(_dLJ6), dLJ12(_dLJ12) {}
};

#endif /*PARAMS_H*/