#include "tatom.h"

//---------------------------------------------------------------------------------------------------
bool TAtom::bRecognizeCoords()
{
  dX = Textutils::dToDouble(sPDB_tag.substr(30,8));
  dY = Textutils::dToDouble(sPDB_tag.substr(38,8));
  dZ = Textutils::dToDouble(sPDB_tag.substr(46,8));
  
  return true;
}

//---------------------------------------------------------------------------------------------------
bool TAtom::bRecognizeNumber()
{
  iNum = Textutils::iToInt(sPDB_tag.substr(6,5));
  return iNum > 0;
}

//---------------------------------------------------------------------------------------------------
bool TAtom::bRecognizeResidue()
{
  bool bRetVal = false;
  std::string sResTemp = sPDB_tag.substr(17,3); 
  Textutils::bCropSpaces(sResTemp);

  Textutils::bCorrectResName(sResTemp);
  
  if(sResTemp == "")
  {
    std::cerr << "Error: no residue name (characters 18-20) in line: " << std::endl
	      << sPDB_tag << std::endl
	      << "                 ^^^ residue name expected here." << std::endl;    
    return false;
  }
  
/******************************************************************************************************
*  Now: check that all files describing given residue are available to read.
*  The files are RES.bnd, RES.angl and RES.dihd by default, where RES is the residue name.
*  Impropers (RES.impr) and RES.pdb are optional, therefore not checked here.
******************************************************************************************************/
  if(Textutils::bFileExists(TConsts::PREPDIR+sResTemp +".dscr") &&
     Textutils::bFileExists(TConsts::BONDSDIR+sResTemp+".bnd" ) &&
     Textutils::bFileExists(TConsts::ANGLEDIR+sResTemp+".angl") &&
     Textutils::bFileExists(TConsts::DIHEDDIR+sResTemp+".dihd")  )
  {
    sResidue = sResTemp;
    bRetVal = true;
  }
  else
  {
    std::string sResTempUpper = Textutils::sToUpper(sResTemp);
    if(Textutils::bFileExists(TConsts::PREPDIR+sResTempUpper +".dscr") &&
       Textutils::bFileExists(TConsts::BONDSDIR+sResTempUpper+".bnd" ) &&
       Textutils::bFileExists(TConsts::ANGLEDIR+sResTempUpper+".angl") &&
       Textutils::bFileExists(TConsts::DIHEDDIR+sResTempUpper+".dihd")  )
    {
      std::cerr << "Error: unknown residue: " << sResTemp << (". Did you mean: "+sResTempUpper+"?") << std::endl;
    }
    else std::cerr << "Error: unknown residue: " << sResTemp << std::endl;
  }
  
/******************************************************************************************************
*  Now: recognize 'residue  identity', i.e. its number with chain identifier.
*  Reserved field (column 21) and insertion code (col 27) are omitted
******************************************************************************************************/
  sResIdent = sPDB_tag.substr(20, 7);
  if(Textutils::sCropSpaces(sResIdent) == "")
  {
    std::cerr << "Error: no sequence identifier (characters 21-27) in line: " << std::endl
	      << sPDB_tag << std::endl
	      << "                     ^^^^^ chain letter code + residue number expected here." << std::endl;
    bRetVal = false;
  }
  else
  {
    iResidueNumber = Textutils::iToInt(sPDB_tag.substr(22, 4));
    cChainCode     = sPDB_tag[21];
  }
  
  return bRetVal;
}

//---------------------------------------------------------------------------------------------------
bool TAtom::bRecognizeName()
{
  sName = sPDB_tag.substr(12,4); 
  Textutils::bCropSpaces(sName);  
  Textutils::bCorrectAtomName(sName);		// Change the name, according to our internal standards.
  return true;
}

//---------------------------------------------------------------------------------------------------
bool TAtom::bRecognizeType()
{
/******************************************************************************************************
 *  Now: read the partial atomic charges from appropriate *.dscr file.
 *  The path of the file must be like: PROGRAMPATH/forcefield/PREPI/XYZ.dscr where XYZ is the residue name
 ******************************************************************************************************/
  bool bRetVal = false;
  std::string sLine, sAtomType; 
  std::ifstream ifsIn((TConsts::PREPDIR+sResidue+".dscr").c_str());
  if(ifsIn != NULL)
  {
    while(std::getline(ifsIn, sLine))
    {
      if(sLine.length() > 0)
      {
	std::istringstream issIn(sLine);
	std::string sAtomNameTemp;
	issIn >> sAtomNameTemp;
	Textutils::bCorrectAtomName(sAtomNameTemp);

	if(sAtomNameTemp == sName)
	{
	  issIn >> sType >> dCharge >> vdWInfo.dR0 >> vdWInfo.dEpsilon;
	  bRetVal = true;
	  
	  break;
	}
      }
    }
    ifsIn.close();
  }
  
  return bRetVal;
}

//---------------------------------------------------------------------------------------------------
bool TAtom::bRecognizeRestrain()
{
  bool bRetVal = true;
  if(TConsts::POSRESTRAINTS)
  {
    double dOccupancy = Textutils::dToDouble(sPDB_tag.substr(54,6));
    if(dOccupancy > -1e-4)
    {
      if(dOccupancy < 1e-4)
      {
	dPosRestr = -1.; // Completely fixed
	if(TConsts::TALKATIVE > 1) std::clog << "Atom with position fixed: " << sPDB_tag.substr(0,30) << std::endl;
      }
      else if(dOccupancy > 0.9999) dPosRestr = 0.; // Completely free
      else 
      {
	double bFactor = Textutils::dToDouble(sPDB_tag.substr(60,6)); // Spring attached
	if(fabs(bFactor) > 1e-4) dPosRestr = 1./bFactor;
	else dPosRestr = 0.;
	if(TConsts::TALKATIVE > 1) std::clog << "Atom restrained (" << dPosRestr << " kcal/mol/Ang^2): " 
	                                     << sPDB_tag.substr(0,30) << std::endl;
      }
      if(TConsts::bIsNan(dPosRestr)) dPosRestr = 0.;
      bRetVal = true;
    }
    else bRetVal = false;
  }
  return bRetVal;  
}

//---------------------------------------------------------------------------------------------------
TAtom& operator+=(TAtom& _atAtom, const TVector& _vecVector)
{
  _atAtom.dX += _vecVector.dX;
  _atAtom.dY += _vecVector.dY;
  _atAtom.dZ += _vecVector.dZ;
  return _atAtom;
}

//---------------------------------------------------------------------------------------------------
TAtom& operator-=(TAtom& _atAtom, const TVector& _vecVector)
{
  _atAtom.dX -= _vecVector.dX;
  _atAtom.dY -= _vecVector.dY;
  _atAtom.dZ -= _vecVector.dZ;
  return _atAtom;
}

//---------------------------------------------------------------------------------------------------
TVector vecMakeVector(const TAtom& _atAtom1, const TAtom& _atAtom2)
{
  double dX1, dX2, dY1, dY2, dZ1, dZ2;
  _atAtom1.GetCoords(dX1, dY1, dZ1);
  _atAtom2.GetCoords(dX2, dY2, dZ2);
  return TVector(dX2 - dX1, dY2 - dY1, dZ2 - dZ1);
}

//---------------------------------------------------------------------------------------------------
TVector vecMakeVector(const TAtom* _patAtom1, const TAtom* _patAtom2)
{
  double dX1, dX2, dY1, dY2, dZ1, dZ2;
  _patAtom1->GetCoords(dX1, dY1, dZ1);
  _patAtom2->GetCoords(dX2, dY2, dZ2);
  return TVector(dX2 - dX1, dY2 - dY1, dZ2 - dZ1);
}

//---------------------------------------------------------------------------------------------------
TVector vecMakeVector(const TAtom* _patAtom1, const double _dX, const double _dY, const double _dZ)
{
  double dX1, dY1, dZ1;
  _patAtom1->GetCoords(dX1, dY1, dZ1);
  return TVector(_dX - dX1, _dY - dY1, _dZ - dZ1);
}

//---------------------------------------------------------------------------------------------------
double dDistance_SQR(const TAtom& _atAtom1, const TAtom& _atAtom2)
{
  double dX(_atAtom1.dX-_atAtom2.dX), dY(_atAtom1.dY-_atAtom2.dY), dZ(_atAtom1.dZ-_atAtom2.dZ);
  return dX*dX + dY*dY + dZ*dZ;
}

//---------------------------------------------------------------------------------------------------
double dDistance_SQR(const TAtom& _atAtom1, const double _dX, const double _dY, const double _dZ)
{
  double dX(_atAtom1.dX-_dX), dY(_atAtom1.dY-_dY), dZ(_atAtom1.dZ-_dZ);
  return dX*dX + dY*dY + dZ*dZ;
}

//---------------------------------------------------------------------------------------------------
double dDistance_SQR(const TAtom* _patAtom1, const TAtom* _patAtom2)
{
  double dX(_patAtom1->dX-_patAtom2->dX), dY(_patAtom1->dY-_patAtom2->dY), dZ(_patAtom1->dZ-_patAtom2->dZ);
  return dX*dX + dY*dY + dZ*dZ;
}

//---------------------------------------------------------------------------------------------------
double dDistance_SQR(const TAtom* _patAtom1, const double _dX, const double _dY, const double _dZ)
{
  double dX(_patAtom1->dX-_dX), dY(_patAtom1->dY-_dY), dZ(_patAtom1->dZ-_dZ);
  return dX*dX + dY*dY + dZ*dZ;  
}

//---------------------------------------------------------------------------------------------------
bool TAtom::bAssignTag(std::string _sLine)
{
  bool bRetVal = false;  
  if(_sLine.length() >= 54 && (_sLine.substr(0,4) == "ATOM" || _sLine.substr(0,6) == "HETATM"))
  {
    Textutils::DosToUnix(_sLine);
    sPDB_tag = _sLine + (_sLine.length() < 80 ? Textutils::sEmptySpaces(80 - _sLine.length()) : "");
    bRetVal = true;
  }
  return bRetVal;
}

//---------------------------------------------------------------------------------------------------
bool TAtom::bWrite(std::ostream& _osOut_) const
{
  bool bHetatm = sResidue != "A" && sResidue != "RA3" && sResidue != "RA5" && sResidue != "RAN" &&
		 sResidue != "C" && sResidue != "RC3" && sResidue != "RC5" && sResidue != "RCN" &&
		 sResidue != "G" && sResidue != "RG3" && sResidue != "RG5" && sResidue != "RGN" &&
		 sResidue != "U" && sResidue != "RU3" && sResidue != "RU5" && sResidue != "RUN" ;

  return (_osOut_ << (bHetatm ? "HETATM" : "ATOM  ")
		  << std::setfill(' ') << std::right 
		  << std::setw(5) << iNum%100000 << ' '
		  << std::left  << std::setw(4) << ((sName.length() < 4 ? " " : "") + sName) << ' '
		  << std::right << std::setw(3) << sResidue 
		  << sResIdent << "   "
		  << std::right << std::fixed << std::setprecision(3) 
		  << std::setw(8) << dX 
		  << std::setw(8) << dY 
		  << std::setw(8) << dZ 
		  << (sPDB_tag.length() >= 80 ? sPDB_tag.substr(54, 26) : "") << std::endl);
}

//---------------------------------------------------------------------------------------------------
bool TAtom::bRecognize()
{
/*****************************************************************************************************
 * Now follows recognition of various atomic properties, based on data from sPDB_tag.
 * Sequence of actions IS important since recognition of atom type requires knowledge of residue name;
 * recognition of van der Waals params requires atom type, etc, but C++ is known for 
 * short-circuited evaluation of such expressions so the code below is OK. (C++03 standard, 6.5.13).
 *****************************************************************************************************/

 return sPDB_tag.length() >= 80 && (sPDB_tag.substr(0,4) == "ATOM" || sPDB_tag.substr(0,6) == "HETATM") &&
	bRecognizeCoords() && bRecognizeNumber() && bRecognizeResidue() && bRecognizeName() && 
	bRecognizeName()   && bRecognizeType()   && bRecognizeRestrain();
}

//---------------------------------------------------------------------------------------------------
std::string TAtom::sGetResIdent() const
{
  return sResIdent;
}

//---------------------------------------------------------------------------------------------------
bool TAtom::bSetResIdent(const std::string& _sNewResIdent)
{
  if(_sNewResIdent.length() > 7)
  {
    std::cerr << "Internal Warning: residue identifier exceeds 7 characters." << std::endl;
    sPDB_tag.replace(20, 7, _sNewResIdent.substr(0, 7));
  }
  else sPDB_tag.replace(20, 7, Textutils::sEmptySpaces(7-_sNewResIdent.length())+_sNewResIdent);  
  return bRecognizeResidue();
}

//---------------------------------------------------------------------------------------------------
std::string TAtom::sGetResName() const
{
  return sResidue;
}

//---------------------------------------------------------------------------------------------------
bool TAtom::bSetResName(std::string _sNewResName)
{
  if(_sNewResName.length() > 3)
  {
    std::cerr << "Internal Warning: residue name exceeds 3 characters." << std::endl;
    sPDB_tag.replace(17, 3, _sNewResName.substr(0, 3));
  }
  else sPDB_tag.replace(17, 3, Textutils::sEmptySpaces(3-_sNewResName.length())+_sNewResName);
  return bRecognizeResidue();
}

//---------------------------------------------------------------------------------------------------
std::string TAtom::sGetName() const
{
  return sName;
}

//---------------------------------------------------------------------------------------------------
void TAtom::SetName(std::string _sName)
{
  std::ostringstream ossOut;
  ossOut << std::setfill(' ') << std::left << std::setw(4) << ((_sName.length() < 4 ? " " : "") + _sName);
  
  sPDB_tag.replace(12, 4, ossOut.str());
  bRecognizeName();
}

//---------------------------------------------------------------------------------------------------
std::string TAtom::sGetType() const
{
  return sType;
}

//---------------------------------------------------------------------------------------------------
void TAtom::GetCoords(double& _dX_, double& _dY_, double& _dZ_) const
{
  _dX_ = dX;
  _dY_ = dY;
  _dZ_ = dZ;
}

//---------------------------------------------------------------------------------------------------
int  TAtom::iGetResidueNumber() const
{
  return iResidueNumber;
}

//---------------------------------------------------------------------------------------------------
char TAtom::cGetChainCode() const
{
  return cChainCode;
}

//---------------------------------------------------------------------------------------------------
bool TAtom::bTransform(TAtom _at1,      TAtom _at2,      TAtom _at3,
		       TAtom _at1Prime, TAtom _at2Prime, TAtom _at3Prime)
{
 /**************************************************************************************************************
  * Transforms *this atom to a 'primed' frame of reference. Atoms _at1, _at2, _at3... are some examples 
  * from non-primed frame and _at1Prime, _at2Prime and _atPrime3 are the same atoms in 'primed' frame, given to
  * find the orientation between both frames. Approximately, of course.
  **************************************************************************************************************/
  bool bRetVal = true;

  TVector vec1ToThis   = vecMakeVector(_at1, *this),	// Vector pointing from _at1 to *this
	  vec1To1Prime = vecMakeVector(_at1, _at1Prime),// Vector pointing from _at1 to _at1Prime
  /* {e1,e2,e3} - some orthonormal basis in nonprimed frame follows. */
	  e1	          = vecMakeVector(_at1, _at2),
	  e2              = vecMakeVector(_at1, _at3),
	  e3,
  /* {e1Prime,e2Prime,e3Prime} - corresponding basis in primed frame. */
	  e1Prime         = vecMakeVector(_at1Prime, _at2Prime),
	  e2Prime         = vecMakeVector(_at1Prime, _at3Prime),
	  e3Prime;
  
  e1.Normalize();	//
  e2 -= (e2*e1)*e1;	// Gram-Schmidt orthogonalization
  e2.Normalize();	// Now normalization...
  e3 = e1.vecCrossProd(e2);
  e3.Normalize();	// Just to get sure.
  
  e1Prime.Normalize();	// Analogous procedure in primed frame.
  e2Prime -= (e2Prime*e1Prime)*e1Prime;
  e2Prime.Normalize();
  e3Prime = e1Prime.vecCrossProd(e2Prime);
  e3Prime.Normalize();
	  
  double c1 = e1*vec1ToThis,	//
	 c2 = e2*vec1ToThis,	// Coefficients of expansion of vec1ToThis in the basis of {ei}
	 c3 = e3*vec1ToThis;	//
  
  //Now move your nonprimed ass:
  (*this) += -vec1ToThis	// First, move it to the position of _at1
	     +vec1To1Prime	// ...then move to the position of _at1Prime
	     + c1*e1Prime + c2*e2Prime + c3*e3Prime;	// ...then rebuild your position with primed components
  
  return bRetVal;
}

//---------------------------------------------------------------------------------------------------
std::string TAtom::sGetPDB_tag() const
{
  return sPDB_tag;
}

//---------------------------------------------------------------------------------------------------
double TAtom::dGetvdWRadius() const
{
  return vdWInfo.dR0;
}

//---------------------------------------------------------------------------------------------------
