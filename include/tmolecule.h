#pragma once
#ifndef  TMOLECULE_H
#define  TMOLECULE_H

#ifndef  SEQUENTIAL
  #include <pthread.h>
#endif /*SEQUENTIAL*/

#include <iostream>	
#include <fstream>
#include <istream>
#include <ostream>
#include <vector>
#include <iterator>
#include <map>
#include <stack>
#include <cassert>
#include <list>	

#include "tatom.h"
#include "tvector.h"

#include "tbond.h"
#include "tangle.h"
#include "tdihed.h"
#include "timpro.h"

#include "tnonbond.h"
#include "tposrestr.h"
#include "telectr.h"
#include "tlennard.h"
#include "thbond.h"
#include "tspring.h"
#include "tbasepair.h"

#include "params.h"

class TMolecule
{
  private:
    std::vector<TAtom>   Atoms;		// Table that contains descr. of atoms in molecule
    std::list<TBond>     Bonds;		// All bonds in molecule
    std::list<TAngle>    Angles;	// All planar angles in molecule
    std::list<TDihedral> Dihedrals;	// Dihedral angles in molecule
    std::list<TImproper> Impropers;	// Improper dihedral angles in molecule
    std::list<TPosRestr> PosRestrs;     // Positional restraints
    
    mutable std::vector<std::vector<char> > vvContacts;	// Atoms' contacts within given molecule
    static std::vector<
             std::pair<
               std::list<TElectr>::const_iterator, std::list<TElectr>::const_iterator> > IteratorsPARALLEL;
    
#ifndef  SEQUENTIAL
    static pthread_mutex_t mutexPARALLEL;
    static double dElectrEnergyPARALLEL;
#endif /*SEQUENTIAL*/
    
/********************************************************************************
 * Static member TMolecule::Nonbond is a std::vector where pairs of atoms 
 * are stored (ones that interact in non-covalent way). It is static, so
 * common for all molecules.
 ********************************************************************************/
    static std::list<TElectr>   Electros;
    static std::list<TLennard>  Lennards;
    static std::list<THBond>    HBonds;
    static std::list<TSpring>   Springs;
    static std::list<TBasePair> BasePairs;
    
    bool bMakeIntraBonds(std::map<std::string, int>& AtomContentMap,
			 const std::vector<std::pair<std::string,std::string> >& ResSequence);
    bool bMakeInterBonds(std::map<std::string, int>& AtomContentMap,
			 const std::vector<std::pair<std::string,std::string> >& ResSequence);
    bool bMakeIntraAngles(std::map<std::string, int>& AtomContentMap,
			  const std::vector<std::pair<std::string,std::string> >& ResSequence);
    bool bMakeInterAngles(std::map<std::string, int>& AtomContentMap,
			  const std::vector<std::pair<std::string,std::string> >& ResSequence);
    bool bAddAngle(const TBond& _boBond1, const TBond& _boBond2);
    bool bMakeIntraDihedrals(std::map<std::string, int>& AtomContentMap,
			     const std::vector<std::pair<std::string,std::string> >& ResSequence);
    bool bMakeInterDihedrals(std::map<std::string, int>& AtomContentMap,
			     const std::vector<std::pair<std::string,std::string> >& ResSequence);
    bool bAddDihedral(const TAngle& _anAngle1, const TAngle& _anAngle2);
    bool bMakeIntraImpropers(std::map<std::string, int>& AtomContentMap,
			     const std::vector<std::pair<std::string,std::string> >& ResSequence);
    bool bMakePosRestraints();

    bool bMakeElectrosAndLennardsWithin();
    bool bMakeHBondsWithin();
    
    static bool bMakeElectrosAndLennardsExtern(TMolecule& _molMol1, TMolecule& _molMol2);
    static bool bMakeHBondsExtern(             TMolecule& _molMol1, TMolecule& _molMol2);
    static bool bMakeSprings(                  TMolecule& _molMol1, TMolecule& _molMol2);
    static bool bMakeBPRestraints(             TMolecule& _molMol1, TMolecule& _molMol2);

    bool bBuildContacts(std::map<std::string, int>& AtomContentMap) const;
    void ClearContacts() const;
    
    // Static private - a hack to ensure no acces to data of TMolecule object.
    // Although it could be defined outside the class..:
    static int iFillMap(std::map<std::string, TAtom>& _ResPDBMap_, const std::string& _sResName);
    
    static bool bCompleteResidue(std::vector<TAtom>& AtomsNew, 
				 std::map<std::string, int>&   CurrentResContent,
				 std::map<std::string, TAtom>& CurrentResPDB,
			   const std::string& sCurrentResName,
			   const std::string& sCurrentResIdent);
			  
    bool bAddMissingAtoms();
    bool bAddOP3();
    
    bool bIsPhosphorusAt5Term() const;
    bool bIsOP3At5Term() const;
    
    bool bCorrectTerminalResidueNames();
    void CorrectAtomNumbers();

/********************************************************************************
 * A common hack follows: declaration of operator= as a private member.
 * To avoid being used and discourage the compiler from defining one.
 * No definition (in *.cpp) is needed for this operator since it is never used.
 ********************************************************************************/
    TMolecule& operator=(const TMolecule& _moMol);	// Avoid copying! External pointers to Atoms[] elements don't like it!
    int  iReadAtoms(std::istream& _isIn);
    bool bRecognizeResidues();
    bool bRecognizeNames();
    bool bRecognizeAtoms();
		 
  public:
    TMolecule(): Atoms(), Bonds(), Angles(), Dihedrals(), Impropers(), vvContacts() {}
    TMolecule(const TMolecule& _moMol);		// Dummy copy constructor, never use it! See body for details.
    
    void CalcBornRadii();
    bool bRead( std::istream& _isIn);
    bool bWrite(std::ostream& _osOut_) const;
    bool bEmpty() const;
    unsigned int uGetNumberOfAtoms() const;    
    static unsigned int uGetNumberOfElectros();
    static unsigned int uGetNumberOfLennards();
    static unsigned int uGetNumberOfHBonds();
    static unsigned int uGetNumberOfSprings();
    unsigned int uGetNumberOfPosRestrs();
    static unsigned int uGetNumberOfBasePairs();

    // A wrap around bMakeElectrosAndLennardsExtern(..), bMakeHBondsExtern(..), 
    // bMakeElectrosAndLennardsWithin() and bMakeHBondsWithin():
    friend bool bMakeNonbondsBetween(TMolecule& _molMol1, TMolecule& _molMol2);
    friend bool bAddSpecificBasePair(TMolecule& _molMol1, const std::string& _sResIdent1, 
				     TMolecule& _molMol2, const std::string& _sResIdent2);

    bool bImposeVienna(const std::string& _Vienna);
				     
    static bool bMakeBasePairs();

    void   CalcCovalentGradients();
    bool   bClearGradients();
    double dCalcCovalentEnergy(double& _dMolEnergy_) const;
    double dCalcBondsEnergy() const;
    double dCalcAnglesEnergy() const;
    double dCalcDihedralsEnergy() const;
    double dCalcImpropersEnergy() const;
    double dCalcPosRestrEnergy() const;
   
    void   Shift(double _dFactor);
    double dAllVectorsNorm();

    static bool   bCalcNonbondGradients();
#ifdef SEQUENTIAL
    static double dCalcNonbondEnergy(double& _dMolEnergy_);
#endif /*SEQUENTIAL*/
    
#ifndef SEQUENTIAL
    static double dCalcNonbondEnergyPARALLEL(double& _dMolEnergy_);
    static bool bMakeIteratorsPARALLEL();
    static void* pvCalcIntervalEnergyPARALLEL(void* arg);
#endif /*SEQUENTIAL*/
};

#endif /*TMOLECULE_H*/