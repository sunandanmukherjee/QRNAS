#include "textutils.h"

namespace Textutils
{
  const char *pcBar1 = "-----------------------------------------------------------------------------\n",
	     *pcBar2 = "=============================================================================\n";
  
//---------------------------------------------------------------------------------------------------
  std::string sEmptySpaces(const int _iN)	//what are we living for?
  {
    if(_iN <= 0) return "";
    else return std::string(_iN, ' ');
  }

//---------------------------------------------------------------------------------------------------
  template <class T> std::string sToString(const T& _tArg_)
  {
    std::ostringstream ossOut;
    ossOut << _tArg_;
    return (ossOut.str());
  }

//---------------------------------------------------------------------------------------------------
  double dToDouble(const std::string& _sStr)
  {
    std::istringstream issIn(_sStr);
    double dRetVal;
    if(!(issIn >> dRetVal)) return -1.;
    else return dRetVal;
  }
  
//---------------------------------------------------------------------------------------------------
  int iToInt(const std::string& _sStr)
  {
    std::istringstream issIn(_sStr);
    int iRetVal;
    if(!(issIn >> iRetVal)) return -1;
    else return iRetVal;    
  }

//---------------------------------------------------------------------------------------------------
  bool bCropSpaces(std::string& _sStr_)
  {// Trims leading and trailing spaces.
    bool bRetVal = false;
  
    unsigned int uStartPos = _sStr_.find_first_not_of(" \t");
    unsigned int uEndPos = _sStr_.find_last_not_of(" \t");

    if((std::string::npos == uStartPos ) || (std::string::npos == uEndPos)) _sStr_ = "";
    else
    {
      _sStr_ = _sStr_.substr(uStartPos, uEndPos - uStartPos + 1);
      bRetVal = true;
    }
    
    return bRetVal;
  }

  std::string sCropSpaces(std::string _sStr_)
  {// Clone of bCropSpaces, overloaded to allow for piping.
    bCropSpaces(_sStr_);
    return _sStr_;
  }
  
//---------------------------------------------------------------------------------------------------
  std::string sMakeOutfileName(std::string _sStr_)
  {// Given something like "../structures/infile.pdb" make "../structures/infile_out.pdb"
    if(_sStr_.length() > 4 && sToUpper(_sStr_.substr(_sStr_.length()-4, 4)) == ".PDB")
    {
      _sStr_.insert(_sStr_.length()-4, "_out");
    }
    else _sStr_ = _sStr_ + "_out";
      
    return _sStr_;
  }

//---------------------------------------------------------------------------------------------------
  std::string sMakeOutfileName(std::string _sStr_, int _iNum)
  {// Given something like "../structures/file.pdb", 7 make "../structures/file_7.pdb"
    if(_sStr_.length() > 4 && sToUpper(_sStr_.substr(_sStr_.length()-4, 4)) == ".PDB")
    {
      _sStr_.insert(_sStr_.length()-4, "_" + sToString(_iNum));
    }
    else _sStr_ = _sStr_ + "_" + sToString(_iNum);
      
    return _sStr_;
  }

//---------------------------------------------------------------------------------------------------
  bool bExtractPath(char* pcArgv1)
  {// Given sth like "/usr/bin/QRNA/executable.exe" make "/usr/bin/QRNA"
    bool bRetVal = true;    
    TConsts::PROGRAMPATH = std::string(pcArgv1);
    if(TConsts::PROGRAMPATH == "") bRetVal = false;
    else
    {
      unsigned int uEndPos = TConsts::PROGRAMPATH.find_last_of("/\\");
      if(uEndPos == std::string::npos) bRetVal = false;
      else TConsts::PROGRAMPATH = TConsts::PROGRAMPATH.substr(0, uEndPos);
    }
    return bRetVal;
  }

//---------------------------------------------------------------------------------------------------
  void ToUpper(std::string& _sStr_)
  {// ::toupper can be a macro. To avoid all its drawbacks we write it in parentheses, to ensure
   // that a function under same name is involved instead.
    std::transform(_sStr_.begin(), _sStr_.end(), _sStr_.begin(), (::toupper));
  }
  
  std::string sToUpper(std::string _sStr)
  {// Same as above, but can be piped.
    std::transform(_sStr.begin(), _sStr.end(), _sStr.begin(), (::toupper));
    return _sStr;
  }

//---------------------------------------------------------------------------------------------------
  void DosToUnix(std::string& _sPDBTag_)
  {
    if(_sPDBTag_.length() > 1 && _sPDBTag_[_sPDBTag_.length()-2] == '\r') // is this OK? to be checked
      _sPDBTag_ = _sPDBTag_.substr(0, _sPDBTag_.length()-2)+'\n';
  }

//---------------------------------------------------------------------------------------------------
  bool bFileExists(const std::string& _sFileName)
  {
    return (std::ifstream(_sFileName.c_str()) != NULL);
  }
  
//---------------------------------------------------------------------------------------------------
  bool bCorrectAtomName(std::string& _sAtomName_)
  {//corrects atom names; 'false' means: 'nothing changed'
    bool bRetVal = true;
    
    if(_sAtomName_.empty()) return false;
    
    ToUpper(_sAtomName_);
    
    int iStar = _sAtomName_.find('*');
    if(iStar != (int)std::string::npos) _sAtomName_[iStar] = '\'';	//replace stars with primes
    
         if(_sAtomName_ == "O1P")  _sAtomName_ = "OP1";
    else if(_sAtomName_ == "O2P")  _sAtomName_ = "OP2";
    else if(_sAtomName_ == "O3P")  _sAtomName_ = "OP3";
    else if(_sAtomName_ == "PA")   _sAtomName_ = "P";
    else if(_sAtomName_ == "1H2")  _sAtomName_ = "H21";
    else if(_sAtomName_ == "2H2")  _sAtomName_ = "H22";
    else if(_sAtomName_ == "1H4")  _sAtomName_ = "H41";
    else if(_sAtomName_ == "2H4")  _sAtomName_ = "H42";
    else if(_sAtomName_ == "1H6")  _sAtomName_ = "H61";
    else if(_sAtomName_ == "2H6")  _sAtomName_ = "H62";
    else if(_sAtomName_ == "1HM'") _sAtomName_ = "HM'1";
    else if(_sAtomName_ == "2HM'") _sAtomName_ = "HM'2";
    else if(_sAtomName_ == "3HM'") _sAtomName_ = "HM'3";
    else if(_sAtomName_ == "1H5'") _sAtomName_ = "H5'1";
    else if(_sAtomName_ == "2H5'") _sAtomName_ = "H5'2";
    else if(_sAtomName_ == "1H5M") _sAtomName_ = "H5M1";
    else if(_sAtomName_ == "2H5M") _sAtomName_ = "H5M2";
    else if(_sAtomName_ == "3H5M") _sAtomName_ = "H5M3";
    else if(_sAtomName_ == "H5'")  _sAtomName_ = "H5'1";
    else if(_sAtomName_ == "H5''") _sAtomName_ = "H5'2";
    else if(_sAtomName_ == "H2'1") _sAtomName_ = "H2'";
    else if(_sAtomName_ == "HO2'") _sAtomName_ = "HO'2";
    else if(_sAtomName_ == "HO5'") _sAtomName_ = "H5T";
    else if(_sAtomName_ == "HO3'") _sAtomName_ = "H3T";
    else if(_sAtomName_ == "CM'")  _sAtomName_ = "CM2";
    else if(_sAtomName_ == "C5M")  _sAtomName_ = "C51";
    else bRetVal = false;
    
    return bRetVal;
  }

//---------------------------------------------------------------------------------------------------
  bool bCorrectResName(std::string& _sResName_)
  {
    bool bRetVal = true;
     
 	 if(_sResName_ == "RA" ) _sResName_ = "A";
    else if(_sResName_ == "RC" ) _sResName_ = "C";
    else if(_sResName_ == "RG" ) _sResName_ = "G";
    else if(_sResName_ == "RU" ) _sResName_ = "U";
    else if(_sResName_ == "2RA") _sResName_ = "RIA";
    else if(_sResName_ == "G6A") _sResName_ = "6GA";
    else if(_sResName_ == "DMA") _sResName_ = "MA6";
    else if(_sResName_ == "DMG") _sResName_ = "M2G";
    else if(_sResName_ == "6TA") _sResName_ = "T6A";
    else if(_sResName_ == "5MC") _sResName_ = "M5C";
    else if(_sResName_ == "DAG") _sResName_ = "PQ1";
    else if(_sResName_ == "DHU") _sResName_ = "H2U";
    else if(_sResName_ == "MMU") _sResName_ = "2MU";
    else if(_sResName_ == "Q" || _sResName_ == "QUG") _sResName_ = "QUO";
    else bRetVal = false;
    
    return bRetVal;
  }

//---------------------------------------------------------------------------------------------------
  bool bIsHydrogenDonor(const std::string& _sAtomType)
  {
    bool bRetVal = false;
    if(_sAtomType == "H"  ||
       _sAtomType == "HO" ) bRetVal = true;
    
    return bRetVal;
  }

//---------------------------------------------------------------------------------------------------
  bool bIsHydrogenAcceptor(const std::string& _sAtomType)
  {
    bool bRetVal = false;
    if(_sAtomType == "O"  ||
       _sAtomType == "NC" ||
       _sAtomType == "N"  ) bRetVal = true;
    
    return bRetVal;
  }

//---------------------------------------------------------------------------------------------------
  bool bIsWCHydrogenBond(std::string _sAtomName1, std::string _sResName1, 
			 std::string _sAtomName2, std::string _sResName2)
  {// Are the two given atoms bonded by a hydrogen bond when forming a canonical base pair?
    if((_sResName2 == "A" || _sResName2 == "RA3" || _sResName2 == "RA5" || _sResName2 == "RAN" ||
        _sResName2 == "DA"|| _sResName2 == "DA3" || _sResName2 == "DA5" || _sResName2 == "DAN") &&
       (_sResName1 == "U" || _sResName1 == "RU3" || _sResName1 == "RU5" || _sResName1 == "RUN" ||
        _sResName1 == "T" || _sResName1 == "DT"  || _sResName1 == "DT3" || _sResName1 == "DT5" || _sResName1 == "DTN"))
    {
      std::swap(_sResName1, _sResName2);
      std::swap(_sAtomName1, _sAtomName2);
    }
    // So the kosher order is {A,U} (and {G,C}), not {U,A} (nor {C,G}).

    if((_sResName1 == "A" || _sResName1 == "RA3" || _sResName1 == "RA5" || _sResName1 == "RAN" ||
        _sResName1 == "DA"|| _sResName1 == "DA3" || _sResName1 == "DA5" || _sResName1 == "DAN") &&
       (_sResName2 == "U" || _sResName2 == "RU3" || _sResName2 == "RU5" || _sResName2 == "RUN" ||
        _sResName2 == "T" || _sResName2 == "DT"  || _sResName2 == "DT3" || _sResName2 == "DT5" || _sResName2 == "DTN"))
    {
      bool bRetVal = false;
      if((_sAtomName1 == "H61" && _sAtomName2 == "O4") ||
	 (_sAtomName1 == "N1"  && _sAtomName2 == "H3") ) bRetVal = true;
      return bRetVal;
    } 

    if((_sResName1 == "C" || _sResName1 == "RC3" || _sResName1 == "RC5" || _sResName1 == "RCN" ||
        _sResName1 == "DC"|| _sResName1 == "DC3" || _sResName1 == "DC5" || _sResName1 == "DCN") &&
       (_sResName2 == "G" || _sResName2 == "RG3" || _sResName2 == "RG5" || _sResName2 == "RGN" ||
        _sResName2 == "DG"|| _sResName2 == "DG3" || _sResName2 == "DG5" || _sResName2 == "DGN"))
    {
      std::swap(_sResName1, _sResName2);
      std::swap(_sAtomName1, _sAtomName2);      
    }
    if((_sResName2 == "C" || _sResName2 == "RC3" || _sResName2 == "RC5" || _sResName2 == "RCN" ||
        _sResName2 == "DC"|| _sResName2 == "DC3" || _sResName2 == "DC5" || _sResName2 == "DCN") &&
       (_sResName1 == "G" || _sResName1 == "RG3" || _sResName1 == "RG5" || _sResName1 == "RGN" ||
        _sResName1 == "DG"|| _sResName1 == "DG3" || _sResName1 == "DG5" || _sResName1 == "DGN"))
    {
      bool bRetVal = false;
      if((_sAtomName1 == "H21" && _sAtomName2 == "O2" ) ||
	 (_sAtomName1 == "H1"  && _sAtomName2 == "N3" ) ||
	 (_sAtomName1 == "O6"  && _sAtomName2 == "H41") ) bRetVal = true;
      return bRetVal;
    } 
    else return false;
  }
//---------------------------------------------------------------------------------------------------
  bool bIsInBase(const std::string& _sAtomName, const std::string& _sResName)
  {// Is the atom in a basic part of nucleotide?
    bool bRetVal = false;
    if(_sResName == "A" || _sResName == "RA3" || _sResName == "RA5" || _sResName == "RAN" ||
       _sResName == "DA"|| _sResName == "DA3" || _sResName == "DA5" || _sResName == "DAN")
    {
      if(_sAtomName == "N9" ||
	 _sAtomName == "C8" ||
	 _sAtomName == "N7" ||
	 _sAtomName == "C5" ||
	 _sAtomName == "C6" ||
	 _sAtomName == "N6" ||
	 _sAtomName == "N1" ||
	 _sAtomName == "C2" ||
	 _sAtomName == "N3" ||
	 _sAtomName == "C4" ||
	 _sAtomName == "H8" ||
	 _sAtomName == "H61"||
	 _sAtomName == "H62"||
	 _sAtomName == "H2" ) bRetVal = true;
    }
    else if(_sResName == "C" || _sResName == "RC3" || _sResName == "RC5" || _sResName == "RCN" ||
            _sResName == "DC"|| _sResName == "DC3" || _sResName == "DC5" || _sResName == "DCN")
    {
      if(_sAtomName == "N1" ||
	 _sAtomName == "C2" ||
	 _sAtomName == "O2" ||
	 _sAtomName == "N3" ||
	 _sAtomName == "C4" ||
	 _sAtomName == "N4" ||
	 _sAtomName == "C5" ||
	 _sAtomName == "C6" ||
	 _sAtomName == "H41"||
	 _sAtomName == "H42"||
	 _sAtomName == "H6" ||
	 _sAtomName == "H5" ) bRetVal = true;      
    }
    else if(_sResName == "G" || _sResName == "RG3" || _sResName == "RG5" || _sResName == "RGN" ||
            _sResName == "DG"|| _sResName == "DG3" || _sResName == "DG5" || _sResName == "DGN")
    {
      if(_sAtomName == "N9" ||
	 _sAtomName == "C8" ||
	 _sAtomName == "N7" ||
	 _sAtomName == "C5" ||
	 _sAtomName == "C6" ||
	 _sAtomName == "O6" ||
	 _sAtomName == "N1" ||
	 _sAtomName == "C2" ||
	 _sAtomName == "N2" ||
	 _sAtomName == "N3" ||
	 _sAtomName == "C4" ||
	 _sAtomName == "H8" ||
	 _sAtomName == "H21"||
	 _sAtomName == "H22"||
	 _sAtomName == "H1" ) bRetVal = true;      
    }
    else if(_sResName == "U" || _sResName == "RU3" || _sResName == "RU5" || _sResName == "RUN" ||
            _sResName == "T" || _sResName == "DT"  || _sResName == "DT3" || _sResName == "DT5" || _sResName == "DTN")
    {
      if(_sAtomName == "N1" ||
	 _sAtomName == "C2" ||
	 _sAtomName == "O2" ||
	 _sAtomName == "N3" ||
	 _sAtomName == "C4" ||
	 _sAtomName == "O4" ||
	 _sAtomName == "C5" ||
	 _sAtomName == "C6" ||
	 _sAtomName == "C7" || // Specific to thymine.
	 _sAtomName == "H71"|| // Specific to thymine.
	 _sAtomName == "H72"|| // Specific to thymine.
	 _sAtomName == "H73"|| // Specific to thymine.
	 _sAtomName == "H5" ||
	 _sAtomName == "H6" ||
	 _sAtomName == "H3" ) bRetVal = true;      
    }
    
    return bRetVal;
  }
//---------------------------------------------------------------------------------------------------
  bool bValidateVienna(const std::string& _sVienna)
  {
    bool bRetVal = true;
    
    int iBrackets = 0;
    for(unsigned int i = 0; i < _sVienna.length(); ++i)
    {
      if(     _sVienna[i] == '(') iBrackets++;
      else if(_sVienna[i] == ')') iBrackets--; 
      else if(_sVienna[i] == '.') {}
      else
      {
	bRetVal = false;
	break;
      }
      if(iBrackets < 0)
      {
	bRetVal = false;
	break;
      }
    }
    if(iBrackets) bRetVal = false;
    
    return bRetVal;
  }

//---------------------------------------------------------------------------------------------------
}
