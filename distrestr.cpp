#include <sstream>

#include "tatom.h"
#include "distrestr.h"

//---------------------------------------------------------------------------------------------------
TDistRestrDescr::TDistRestrDescr(std::string _sLine): bIsUsed(false), sLine(_sLine)
{
  for(unsigned int u = 0; u < _sLine.length(); ++u)
  {
    if(_sLine[u] == '/'     ) _sLine[u] = ' ';
    else if(_sLine[u] == ',') _sLine[u] = '.';    
  }
  std::istringstream issIn(_sLine);
  std::string sLabel, sChainCode1, sChainCode2;
  double dMinDist, dMaxDist;
  if(issIn >> sLabel >> sChainCode1 >> iResNumber1 >> sAtomName1
                     >> sChainCode2 >> iResNumber2 >> sAtomName2
                     >> dMinDist >> dMaxDist >> bp.dKr)
  {
    if(sLabel == "DISTANCE" && 
       !sChainCode1.empty() && !sChainCode2.empty() &&
       iResNumber1 > 0 && iResNumber2 > 0 && sAtomName1 != "" && sAtomName2 != "" &&
       dMinDist <= dMaxDist && dMinDist >= 0. && bp.dKr > 0.)
    {
      cChainCode1 = sChainCode1[0];
      cChainCode2 = sChainCode2[0];
      bp.dREq = (dMinDist + dMaxDist)/2.;
      Textutils::bCorrectAtomName(sAtomName1);
      Textutils::bCorrectAtomName(sAtomName2);
      isOK = true;
    }
    else isOK = false;
  }
  else isOK = false;
}

//---------------------------------------------------------------------------------------------------
bool TDistRestrDescr::bAtomFits(const TAtom* pat) const
{
  return pat != NULL && (
    (pat->sGetName() == sAtomName1 && pat->iGetResidueNumber() == iResNumber1 && pat->cGetChainCode() == cChainCode1) ||
    (pat->sGetName() == sAtomName2 && pat->iGetResidueNumber() == iResNumber2 && pat->cGetChainCode() == cChainCode2) );
}

//---------------------------------------------------------------------------------------------------
bool TDistRestrDescr::bIsOK() const
{
  return isOK;
}

//---------------------------------------------------------------------------------------------------
const std::string& TDistRestrDescr::sGetLine() const
{
  return sLine;
}

//---------------------------------------------------------------------------------------------------
