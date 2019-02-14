#include <sstream>

#include "tatom.h"
#include "bprestr.h"

//---------------------------------------------------------------------------------------------------
TBasePairRestrDescr::TBasePairRestrDescr(std::string _sLine): bIsUsed(false), sLine(_sLine)
{
  for(unsigned int u = 0; u < _sLine.length(); ++u)
  {
    if(_sLine[u] == '/') _sLine[u] = ' ';
    if(_sLine[u] == ',') _sLine[u] = '.';    
  }
  std::istringstream issIn(_sLine);
  std::string sLabel, sChainCode1, sChainCode2;
  if(issIn >> sLabel >> sChainCode1 >> iResNumber1 >> sChainCode2 >> iResNumber2)
  {
    if(sLabel == "BASEPAIR" && 
       !sChainCode1.empty() && !sChainCode2.empty() && iResNumber1 > 0 && iResNumber2 > 0)
    {
      cChainCode1 = sChainCode1[0];
      cChainCode2 = sChainCode2[0];
      isOK = true;
    }
    else isOK = false;
  }
  else isOK = false;
}

//---------------------------------------------------------------------------------------------------
bool TBasePairRestrDescr::bAtomFits(const TAtom* pat) const
{
  return pat != NULL && (
    (pat->iGetResidueNumber() == iResNumber1 && pat->cGetChainCode() == cChainCode1) ||
    (pat->iGetResidueNumber() == iResNumber2 && pat->cGetChainCode() == cChainCode2) );
}

//---------------------------------------------------------------------------------------------------
bool TBasePairRestrDescr::bIsOK() const
{
  return isOK;
}

//---------------------------------------------------------------------------------------------------
const std::string& TBasePairRestrDescr::sGetLine() const
{
  return sLine;
}

//---------------------------------------------------------------------------------------------------
