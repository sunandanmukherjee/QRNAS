#pragma once
#ifndef  BPRESTR_H
#define  BPRESTR_H

#include <string>

class TAtom;

class TBasePairRestrDescr
{
  public:
    explicit TBasePairRestrDescr(std::string _sLine); // Copying the argument _does_ make sense here.
    bool bAtomFits(const TAtom* pat) const;
    bool bIsOK() const;
    const std::string& sGetLine() const;
    bool bIsUsed;
  private:
    bool isOK;
    char cChainCode1, cChainCode2;
    int iResNumber1, iResNumber2;
    std::string sLine;
};

#endif /*BPRESTR_H*/
