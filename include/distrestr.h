#pragma once
#ifndef  DISTRESTR_H
#define  DISTRESTR_H

#include <string>

#include "params.h"

class TAtom;

class TDistRestrDescr
{
  public:
    explicit TDistRestrDescr(std::string _sLine); // Copying the argument _does_ make sense here.
    bool bAtomFits(const TAtom* pat) const;
    bool bIsOK() const;
    const std::string& sGetLine() const;
    BondParam bp;
    bool bIsUsed;
  private:
    bool isOK;
    std::string sAtomName1, sAtomName2;
    char cChainCode1, cChainCode2;
    int iResNumber1, iResNumber2;
    std::string sLine;
};

#endif /*DISTRESTR_H*/
