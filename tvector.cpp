#include "tvector.h"

//---------------------------------------------------------------------------------------------------
TVector TVector::vecCrossProd(const TVector& _tvVec) const
{
  return TVector(dY*_tvVec.dZ - dZ*_tvVec.dY, dZ*_tvVec.dX - dX*_tvVec.dZ, dX*_tvVec.dY - dY*_tvVec.dX);
}
//---------------------------------------------------------------------------------------------------

double operator*(const TVector& _tvVec1, const TVector& _tvVec2)
{
  return _tvVec1.dX*_tvVec2.dX + _tvVec1.dY*_tvVec2.dY + _tvVec1.dZ*_tvVec2.dZ;
}

//---------------------------------------------------------------------------------------------------
TVector operator*(const double _dMul, const TVector& _tvVec)
{
  return _tvVec*_dMul;
}

//---------------------------------------------------------------------------------------------------
TVector TVector::operator-() const			//unary minus
{
  return TVector(-dX, -dY, -dZ);
}

//---------------------------------------------------------------------------------------------------
TVector TVector::operator+(const TVector& _tvVec) const
{
  return TVector(dX + _tvVec.dX, dY + _tvVec.dY, dZ + _tvVec.dZ);
}

//---------------------------------------------------------------------------------------------------
TVector TVector::operator-(const TVector& _tvVec) const	//binary minus
{
  return TVector(*this) + (-_tvVec);
}

//---------------------------------------------------------------------------------------------------
TVector TVector::operator*(const double& _dMul) const
{
  return TVector(dX*_dMul, dY*_dMul, dZ*_dMul);
}

//---------------------------------------------------------------------------------------------------
TVector TVector::operator/(const double& _dDiv) const
{
  return TVector(dX/_dDiv, dY/_dDiv, dZ/_dDiv);  
}

//---------------------------------------------------------------------------------------------------
TVector& TVector::operator+=(const TVector& _tvVec)
{
  dX += _tvVec.dX;
  dY += _tvVec.dY;
  dZ += _tvVec.dZ;
  return *this;
}

//---------------------------------------------------------------------------------------------------
TVector& TVector::operator-=(const TVector& _tvVec)
{
  dX -= _tvVec.dX;
  dY -= _tvVec.dY;
  dZ -= _tvVec.dZ;
  return *this;
}

//---------------------------------------------------------------------------------------------------
TVector& TVector::operator*=(const double& _dMul)
{
  dX *= _dMul;
  dY *= _dMul;
  dZ *= _dMul;  
  return *this;
}

//---------------------------------------------------------------------------------------------------
TVector& TVector::operator/=(const double& _dDiv)
{
  dX /= _dDiv;
  dY /= _dDiv;
  dZ /= _dDiv;    
  return *this;
}

//---------------------------------------------------------------------------------------------------
double TVector::dLength() const
{
  return sqrt(dLength_SQR());
}

//---------------------------------------------------------------------------------------------------
double TVector::dLength_SQR() const
{
  return (*this)*(*this);
}

//---------------------------------------------------------------------------------------------------
void TVector::Normalize()
{
  (*this) /= dLength();
}
//---------------------------------------------------------------------------------------------------
