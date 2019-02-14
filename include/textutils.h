#pragma once

#ifndef  TXTUTLS_H
#define  TXTUTLS_H

#include <sstream>
#include <fstream>
#include <iomanip>

#include <algorithm>

#include <cstring>

#include "tconsts.h"

namespace Textutils
{
  extern const char *pcBar1, *pcBar2;

  std::string sEmptySpaces(const int _iN);		//what are we living for?
							//returns a std::string with just _n blank spaces
  
  template <class T> std::string sToString(const T& _tArg_);	//universal converter to std::string
  double dToDouble(const std::string& _sStr);
  int iToInt(const std::string& _sStr);
  
  bool bCropSpaces(std::string& _sStr_);
  std::string sCropSpaces(std::string _sStr_);
  
  std::string sMakeOutfileName(std::string _sStr_);
  std::string sMakeOutfileName(std::string _sStr_, int _iNum);
  bool bExtractPath(char* pcArgv1);
  
  void ToUpper(std::string& _sStr_);
  std::string sToUpper(std::string _sStr);

  void DosToUnix(std::string& _sPDBTag_);
  bool bFileExists(const std::string& _sFileName);
  bool bCorrectAtomName(std::string& _sAtomName_);
  bool bCorrectResName(std::string& _sResName_);
  
  bool bIsHydrogenDonor(const std::string& _sAtomType);
  bool bIsHydrogenAcceptor(const std::string& _sAtomType);
  bool bIsWCHydrogenBond(std::string _sAtomName1, std::string _sResName1, 
			 std::string _sAtomName2, std::string _sResName2);
  bool bIsInBase(const std::string& _sAtomName, const std::string& _sResName);
  
  bool bValidateVienna(const std::string& _sVienna);
}

#endif /*TXTUTLS_H*/