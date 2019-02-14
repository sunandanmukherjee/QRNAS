#include "tmolecule.h"

std::list<TElectr>   TMolecule::Electros;
std::list<TLennard>  TMolecule::Lennards;
std::list<THBond>    TMolecule::HBonds;
std::list<TSpring>   TMolecule::Springs;
std::list<TBasePair> TMolecule::BasePairs;

std::vector<std::pair<
  std::list<TElectr>::const_iterator, std::list<TElectr>::const_iterator> > TMolecule::IteratorsPARALLEL;

#ifndef SEQUENTIAL
  pthread_mutex_t TMolecule::mutexPARALLEL = PTHREAD_MUTEX_INITIALIZER;
  double TMolecule::dElectrEnergyPARALLEL;
#endif /*SEQUENTIAL*/
//---------------------------------------------------------------------------------------------------
TMolecule::TMolecule(const TMolecule& _moMol)
{
  /********************************************************************************
  * This copy constructor must be defined due to std::list<> expectations
  * but SHOULD NEVER BE USED unless an empty molecule is about to be copied.
  * (it would cause invalidation of all Bonds, Angles, Dihedrals... since
  * they all contain pointers to Atoms[] elements)
  ********************************************************************************/
  if(!(_moMol.bEmpty()))
  {
    std::cerr << "Fatal Internal Error: forbidden copying of nonempty molecule!" << std::endl;
    assert(0);
  }
}

//---------------------------------------------------------------------------------------------------
int TMolecule::iReadAtoms(std::istream& _isIn)
{
  int iRetVal = 0;
  std::string sLine;
  
  while(std::getline(_isIn, sLine))	// Read single line from input PDB
  {
    TAtom atTemp;
    if(atTemp.bAssignTag(sLine))
    {					// It fits as a description of atom
      Atoms.push_back(atTemp);
      ++iRetVal;
    }
    else if(sLine.length() >= 3 && sLine.substr(0,3) == "TER") break;
  }
  
  return iRetVal;
}

//---------------------------------------------------------------------------------------------------
bool TMolecule::bMakeIntraBonds(std::map<std::string, int>& AtomContentMap,
			   const std::vector<std::pair<std::string,std::string> >& ResSequence)
{// A piece of code not to be fucked with unless necessary.
  bool bRetVal = true;
  if(TConsts::TALKATIVE > 0) std::clog << "Building intraresidual bonds..." << std::endl;
  
  for(unsigned int j = 0; j < ResSequence.size(); ++j)
  {
    std::ifstream ifsIn((TConsts::BONDSDIR+ResSequence[j].second+".bnd").c_str());
    
    if(ifsIn != NULL)
    {
      std::string sBondLine;
      while(std::getline(ifsIn, sBondLine))
      {
	std::istringstream issIn(sBondLine);
	std::string sAtom1Temp, sAtom2Temp;
	double dKrTemp, dREqTemp;
	issIn >> sAtom1Temp >> sAtom2Temp >> dKrTemp >> dREqTemp;

	Textutils::bCorrectAtomName(sAtom1Temp);
	Textutils::bCorrectAtomName(sAtom2Temp);
	
	if(AtomContentMap.find(ResSequence[j].first+sAtom1Temp) != AtomContentMap.end() &&
	    AtomContentMap.find(ResSequence[j].first+sAtom2Temp) != AtomContentMap.end()) 
	{  
	  Bonds.push_back(TBond(Atoms[AtomContentMap[ResSequence[j].first+sAtom1Temp]], 
				Atoms[AtomContentMap[ResSequence[j].first+sAtom2Temp]], 
				BondParam(dKrTemp, dREqTemp) ));
	  if(TConsts::TALKATIVE > 2) 
	    std::clog << "Bond added: " << sAtom1Temp << " " << sAtom2Temp 
	    << " in residue: " << ResSequence[j].first << std::endl;
	}  
	else
	{
	  std::cerr << "Warning: Bond expected, but some of the atoms are missing: " << sAtom1Temp << " " << sAtom2Temp 
		    << ", residue No. " << ResSequence[j].first << std::endl;
	  bRetVal = false;		// Unknown atom name in *.bnd file.
	}
      }
    }
    else
    {
      std::cerr << "Error: Cannot read file: " << (TConsts::BONDSDIR+ResSequence[j].second+".bnd") << std::endl;
      bRetVal = false;
      break;		// Unknown residue is a good reason to break.
    }    
  }  
  
  return bRetVal;
}

//---------------------------------------------------------------------------------------------------
bool TMolecule::bMakeInterBonds(std::map<std::string, int>& AtomContentMap,
			   const std::vector<std::pair<std::string,std::string> >& ResSequence)
{
  bool bRetVal = true;
  if(TConsts::TALKATIVE > 0) std::clog << "Building interresidual bonds..." << std::endl;

  for(unsigned int j = 0; j < ResSequence.size()-1; ++j)
  {  
    if(AtomContentMap.find(ResSequence[j].first+"O3'") !=  AtomContentMap.end() && 
	AtomContentMap.find(ResSequence[j+1].first+"P") !=  AtomContentMap.end())
    {
      Bonds.push_back(TBond(Atoms[AtomContentMap[ResSequence[j].first+"O3'"]], 
			    Atoms[AtomContentMap[ResSequence[j+1].first+"P"]], 
			    TConsts::GeneralBondsData[TConsts::AtomTypeNumbers[Atoms[AtomContentMap[ResSequence[j].first+"O3'"]].sGetType()]]\
						     [TConsts::AtomTypeNumbers[Atoms[AtomContentMap[ResSequence[j+1].first+"P"]].sGetType()]]));
      if(TConsts::TALKATIVE > 2) 
	std::clog << "Interresidual bond added: " 
		  << Atoms[AtomContentMap[ResSequence[j].first+"O3'"]].sGetPDB_tag().substr(0,27) << " --- "
		  << Atoms[AtomContentMap[ResSequence[j+1].first+"P"]].sGetPDB_tag().substr(0,27)
		  << std::endl;
    }
    else 
    {
      std::cerr << "Error: Cannot link residues: " << ResSequence[j].first << " and "
		    << ResSequence[j+1].first << std::endl;
      bRetVal = false;
    }
  }
  return bRetVal;
}

//---------------------------------------------------------------------------------------------------
bool TMolecule::bMakeIntraAngles(std::map<std::string, int>& AtomContentMap,
			    const std::vector<std::pair<std::string,std::string> >& ResSequence)
{
  bool bRetVal = true;
  if(TConsts::TALKATIVE > 0) std::clog << "Building intraresidual angles..." << std::endl;
  
  for(unsigned int j = 0; j < ResSequence.size(); ++j)
  {
    std::ifstream in((TConsts::ANGLEDIR+ResSequence[j].second+".angl").c_str());
    
    if(in != NULL)
    {
      std::string sAnglLine;
      while(std::getline(in, sAnglLine))
      {
	std::istringstream issIn(sAnglLine);
	std::string sAtom1Temp, sAtom2Temp, sAtom3Temp;
	double dKTheta, dThetaEq;
	issIn >> sAtom1Temp >> sAtom2Temp >> sAtom3Temp >> dKTheta >> dThetaEq;

	Textutils::bCorrectAtomName(sAtom1Temp);
	Textutils::bCorrectAtomName(sAtom2Temp);
	Textutils::bCorrectAtomName(sAtom3Temp);
	
	if(AtomContentMap.find(ResSequence[j].first+sAtom1Temp) != AtomContentMap.end() &&
	   AtomContentMap.find(ResSequence[j].first+sAtom2Temp) != AtomContentMap.end() &&
	   AtomContentMap.find(ResSequence[j].first+sAtom3Temp) != AtomContentMap.end() ) 
	{  
	  Angles.push_back(TAngle(Atoms[AtomContentMap[ResSequence[j].first+sAtom1Temp]], 
				  Atoms[AtomContentMap[ResSequence[j].first+sAtom2Temp]],
				  Atoms[AtomContentMap[ResSequence[j].first+sAtom3Temp]],
				  AngleParam(dKTheta, dThetaEq) ));
	  if(TConsts::TALKATIVE > 2) 
	    std::clog << "Angle added: " << sAtom1Temp << " " << sAtom2Temp << " " << sAtom3Temp << " in res: " 
	              << ResSequence[j].first << std::endl;
	}  
	else
	{
	  std::cerr << "Warning: Angle expected, but some of the atoms are missing: " << sAtom1Temp << " " << sAtom2Temp << " " << sAtom3Temp 
		    << ", residue No. " << ResSequence[j].first << std::endl;
	  bRetVal = false;		// Unknown atom name in *.angl file.
	}
      }
    }
    else
    {
      std::cerr << "Error: Cannot read file: " << (TConsts::ANGLEDIR+ResSequence[j].second+".angl") << std::endl;
      bRetVal = false;
      break;		// Unknown residue is a good reason to break.
    }
  }  

  return bRetVal;
}

//---------------------------------------------------------------------------------------------------
bool TMolecule::bMakeInterAngles(std::map<std::string, int>& AtomContentMap,
			    const std::vector<std::pair<std::string,std::string> >& ResSequence)
{
  bool bRetVal = true;
  if(TConsts::TALKATIVE > 0) std::clog << "Building interresidual angles..." << std::endl;
  TConsts::iAnglesInit();

  for(std::list<TBond>::iterator I = Bonds.begin(); I != Bonds.end(); ++I)
  {
    if(I->patAtom1->sGetResIdent() != I->patAtom2->sGetResIdent())  // Interresidual bond - a good source of interresidual angles.
    for(std::list<TBond>::iterator J = Bonds.begin(); J != Bonds.end(); ++J)
    {
      bAddAngle(*I, *J);
    }
  }

  return bRetVal;
}

//---------------------------------------------------------------------------------------------------
bool TMolecule::bAddAngle(const TBond& _boBond1, const TBond& _boBond2)
{
  bool bRetVal = true;
  
  if((_boBond1.patAtom1 == _boBond2.patAtom1 && _boBond1.patAtom2 == _boBond2.patAtom2) ||
     (_boBond1.patAtom2 == _boBond2.patAtom1 && _boBond1.patAtom1 == _boBond2.patAtom2)) bRetVal = false; //case:same bonds
  else 
  {
    TAtom *patAng1 = NULL, *patAng2 = NULL, *patAng3 = NULL;
    if(_boBond1.patAtom1 == _boBond2.patAtom1)
    {
      patAng1 = _boBond1.patAtom2;
      patAng2 = _boBond1.patAtom1;
      patAng3 = _boBond2.patAtom2;
    }
    else if(_boBond1.patAtom1 == _boBond2.patAtom2)
    {
      patAng1 = _boBond1.patAtom2;
      patAng2 = _boBond1.patAtom1;
      patAng3 = _boBond2.patAtom1;
    }
    else if(_boBond1.patAtom2 == _boBond2.patAtom1)
    {
      patAng1 = _boBond1.patAtom1;
      patAng2 = _boBond1.patAtom2;
      patAng3 = _boBond2.patAtom2;
    }
    else if(_boBond1.patAtom2 == _boBond2.patAtom2)
    {
      patAng1 = _boBond1.patAtom1;
      patAng2 = _boBond1.patAtom2;
      patAng3 = _boBond2.patAtom1;
    }
    
    if(patAng1 != NULL && patAng2 != NULL && patAng3 != NULL)
    {
      bRetVal = true;
      AngleParam apTemp = TConsts::GeneralAnglesData[TConsts::AtomTypeNumbers[patAng1->sGetType()]]\
						    [TConsts::AtomTypeNumbers[patAng2->sGetType()]]\
						    [TConsts::AtomTypeNumbers[patAng3->sGetType()]];
      Angles.push_back(TAngle(patAng1, patAng2, patAng3, apTemp));
      if(TConsts::TALKATIVE > 2)
      {
	std::clog << "Interresidual angle added: " 
		  << patAng1->sGetName() << " (res " << patAng1->sGetResIdent() << ") --- "
		  << patAng2->sGetName() << " (res " << patAng2->sGetResIdent() << ") --- "
		  << patAng3->sGetName() << " (res " << patAng3->sGetResIdent() << ")" << std::endl;
      }
    }
    else bRetVal = false;
  }
  return bRetVal;
}

//---------------------------------------------------------------------------------------------------
bool TMolecule::bMakeIntraDihedrals(std::map<std::string, int>& AtomContentMap,
			       const std::vector<std::pair<std::string,std::string> >& ResSequence)
{
  bool bRetVal = true;
  if(TConsts::TALKATIVE > 0) std::clog << "Building intraresidual dihedrals..." << std::endl;
  
  /*************************************************************************************
   * Now follows the formation of 'intraresidual' dihedrals.
   ************************************************************************************/
  
  for(unsigned int j = 0; j < ResSequence.size(); ++j)
  {
    std::ifstream ifsIn((TConsts::DIHEDDIR+ResSequence[j].second+".dihd").c_str());
    
    if(ifsIn != NULL)
    {
      std::string sDihLine;
      while(std::getline(ifsIn, sDihLine))
      {
	std::istringstream issIn(sDihLine);
	std::string sAtom1Temp, sAtom2Temp, sAtom3Temp, sAtom4Temp;
	int iDiv; double dKI, dPhase, dPN;
	issIn >> sAtom1Temp >> sAtom2Temp >> sAtom3Temp >> sAtom4Temp >> iDiv >> dKI >> dPhase >> dPN;

	Textutils::bCorrectAtomName(sAtom1Temp);
	Textutils::bCorrectAtomName(sAtom2Temp);
	Textutils::bCorrectAtomName(sAtom3Temp);
	Textutils::bCorrectAtomName(sAtom4Temp);
	
	if(AtomContentMap.find(ResSequence[j].first+sAtom1Temp) != AtomContentMap.end() &&
	   AtomContentMap.find(ResSequence[j].first+sAtom2Temp) != AtomContentMap.end() &&
	   AtomContentMap.find(ResSequence[j].first+sAtom3Temp) != AtomContentMap.end() && 
	   AtomContentMap.find(ResSequence[j].first+sAtom4Temp) != AtomContentMap.end() ) 
	{  
	  Dihedrals.push_back(TDihedral(Atoms[AtomContentMap[ResSequence[j].first+sAtom1Temp]], 
					Atoms[AtomContentMap[ResSequence[j].first+sAtom2Temp]],
					Atoms[AtomContentMap[ResSequence[j].first+sAtom3Temp]],
					Atoms[AtomContentMap[ResSequence[j].first+sAtom4Temp]],
					DihedralParam(iDiv, dKI, dPhase, dPN) ));
        if(TConsts::TALKATIVE > 2) 
	  std::clog << "Dihedral added: " << sAtom1Temp << " " << sAtom2Temp << " " << sAtom3Temp  << " "
		    << sAtom4Temp << " in res No. " << ResSequence[j].first << std::endl;
	}  
	else
	{
	  std::cerr << "Warning: Dihedral expected, but some of the atoms are missing: " << sAtom1Temp << " " << sAtom2Temp << " " << sAtom3Temp 
		     << " " << sAtom4Temp << ", residue No. " << ResSequence[j].first << std::endl;
	  bRetVal = false;		// Unknown atom name in *.dihd file.
	}
      }
    }
    else
    {
      std::cerr << "Error: cannot read file: " << (TConsts::DIHEDDIR+ResSequence[j].second+".dihd") << std::endl;
      bRetVal = false;
      break;		// Unknown residue is a good reason to break.
    }
  }
  
  return bRetVal;
}

//---------------------------------------------------------------------------------------------------
bool TMolecule::bMakeInterDihedrals(std::map<std::string, int>& AtomContentMap,
			       const std::vector<std::pair<std::string,std::string> >& ResSequence)
{
  bool bRetVal = true;
  if(TConsts::TALKATIVE > 0) std::clog << "Building interresidual dihedrals..." << std::endl;
  
  TConsts::iDihedralsInit();
  
  for(std::list<TBond>::iterator K = Bonds.begin(); K != Bonds.end(); ++K)
  {
    if(K->patAtom1->sGetResIdent() != K->patAtom2->sGetResIdent())  // Interresidual bond - a good source of interresidual dihedrals.
    {
      for(std::list<TAngle>::iterator L = Angles.begin(); L != Angles.end(); ++L)
      {
	if(L->bContains(*K))
	{
	  for(std::list<TAngle>::iterator M = Angles.begin(); M != Angles.end(); ++M)
	    if(!(M->bContains(*K))) bAddDihedral(*L, *M);
	  for(std::list<TAngle>::iterator M = Angles.begin(); M != L; ++M)
	    if(M->bContains(*K)) bAddDihedral(*L, *M);
	}
      }
    }
  }
  
  return bRetVal;
}

//---------------------------------------------------------------------------------------------------
bool TMolecule::bAddDihedral(const TAngle& _anAngle1, const TAngle& _anAngle2)
{
  bool bRetVal = true;
  if((_anAngle1.patAtom1 == _anAngle2.patAtom1 && _anAngle1.patAtom2 == _anAngle2.patAtom2 && _anAngle1.patAtom3 == _anAngle2.patAtom3) || 
     (_anAngle1.patAtom1 == _anAngle2.patAtom3 && _anAngle1.patAtom2 == _anAngle2.patAtom2 && _anAngle1.patAtom3 == _anAngle2.patAtom1)) 
      bRetVal = false;	//the same angle
  else 
  {
    TAtom *patDih1 = NULL, *patDih2 = NULL, 
	  *patDih3 = NULL, *patDih4 = NULL;
    
    if(     _anAngle1.patAtom1 == _anAngle2.patAtom2 && _anAngle1.patAtom2 == _anAngle2.patAtom1)
    {
      patDih1 = _anAngle2.patAtom3,
      patDih2 = _anAngle2.patAtom2,
      patDih3 = _anAngle2.patAtom1,
      patDih4 = _anAngle1.patAtom3;
    }
    else if(_anAngle1.patAtom1 == _anAngle2.patAtom2 && _anAngle1.patAtom2 == _anAngle2.patAtom3)
    {
      patDih1 = _anAngle2.patAtom1,
      patDih2 = _anAngle2.patAtom2,
      patDih3 = _anAngle2.patAtom3,
      patDih4 = _anAngle1.patAtom3;
    }
    else if(_anAngle1.patAtom2 == _anAngle2.patAtom1 && _anAngle1.patAtom3 == _anAngle2.patAtom2)
    {
      patDih1 = _anAngle1.patAtom1,
      patDih2 = _anAngle2.patAtom1,
      patDih3 = _anAngle2.patAtom2,
      patDih4 = _anAngle2.patAtom3;
    }
    else if(_anAngle1.patAtom2 == _anAngle2.patAtom3 && _anAngle1.patAtom3 == _anAngle2.patAtom2)
    {
      patDih1 = _anAngle1.patAtom1,
      patDih2 = _anAngle2.patAtom3,
      patDih3 = _anAngle2.patAtom2,
      patDih4 = _anAngle2.patAtom1;
    }
    
    if(patDih1 != NULL && patDih2 != NULL && patDih3 != NULL && patDih4 != NULL)
    {
      if(TConsts::TALKATIVE > 2)
      std::clog << "Interresidual dihedral added: " 
	        << patDih1->sGetName() << " (res " << patDih1->sGetResIdent() << ") --- "
		<< patDih2->sGetName() << " (res " << patDih2->sGetResIdent() << ") --- "
		<< patDih3->sGetName() << " (res " << patDih3->sGetResIdent() << ") --- "
		<< patDih4->sGetName() << " (res " << patDih4->sGetResIdent() << ")" << std::endl;

      for(std::list<DihedralParam>::iterator I = TConsts::GeneralDihedralsData[TConsts::AtomTypeNumbers[patDih1->sGetType()]]\
									      [TConsts::AtomTypeNumbers[patDih2->sGetType()]]\
									      [TConsts::AtomTypeNumbers[patDih3->sGetType()]]\
									      [TConsts::AtomTypeNumbers[patDih4->sGetType()]].begin();
					    I != TConsts::GeneralDihedralsData[TConsts::AtomTypeNumbers[patDih1->sGetType()]]\
									      [TConsts::AtomTypeNumbers[patDih2->sGetType()]]\
									      [TConsts::AtomTypeNumbers[patDih3->sGetType()]]\
									      [TConsts::AtomTypeNumbers[patDih4->sGetType()]].end();
					    ++I )
      {
	I->dPN = fabs(I->dPN);	// If multiple fourier coeffs appear some of them are negative.
	Dihedrals.push_back(TDihedral(patDih1, patDih2, patDih3, patDih4, *I));
      }
      
    }
    else bRetVal = false;
  }
  
  return bRetVal;
}

//---------------------------------------------------------------------------------------------------
bool TMolecule::bMakeIntraImpropers(std::map<std::string, int>& AtomContentMap,
			       const std::vector<std::pair<std::string,std::string> >& ResSequence)
{
  bool bRetVal = true;
  if(TConsts::TALKATIVE > 0) std::clog << "Building intraresidual impropers..." << std::endl;
  /*************************************************************************************
   * Now follows the formation of 'intraresidual' impropers.
   ************************************************************************************/
  
  for(unsigned int j = 0; j < ResSequence.size(); ++j)
  {
    std::ifstream ifsIn((TConsts::IMPRODIR+ResSequence[j].second+".impr").c_str());
    
    if(ifsIn != NULL)
    {
      std::string sImprLine;
      while(std::getline(ifsIn, sImprLine))
      {
	std::istringstream issIn(sImprLine);
	std::string sAtom1Temp, sAtom2Temp, sAtom3Temp, sAtom4Temp;
	double dKI, dPhase, dPN;
	issIn >> sAtom1Temp >> sAtom2Temp >> sAtom3Temp >> sAtom4Temp >> dKI >> dPhase >> dPN;

	Textutils::bCorrectAtomName(sAtom1Temp);
	Textutils::bCorrectAtomName(sAtom2Temp);
	Textutils::bCorrectAtomName(sAtom3Temp);
	Textutils::bCorrectAtomName(sAtom4Temp);
	
	if(AtomContentMap.find(ResSequence[j].first+sAtom1Temp) != AtomContentMap.end() &&
	   AtomContentMap.find(ResSequence[j].first+sAtom2Temp) != AtomContentMap.end() &&
	   AtomContentMap.find(ResSequence[j].first+sAtom3Temp) != AtomContentMap.end() && 
	   AtomContentMap.find(ResSequence[j].first+sAtom4Temp) != AtomContentMap.end() ) 
	{  
	  Impropers.push_back(TImproper(Atoms[AtomContentMap[ResSequence[j].first+sAtom1Temp]], 
					Atoms[AtomContentMap[ResSequence[j].first+sAtom2Temp]],
					Atoms[AtomContentMap[ResSequence[j].first+sAtom3Temp]],
					Atoms[AtomContentMap[ResSequence[j].first+sAtom4Temp]],
					ImproperParam(dKI, dPhase, dPN) ));
	  if(TConsts::TALKATIVE > 2)
	    std::clog << "Improper added: " << sAtom1Temp << " " << sAtom2Temp << " " << sAtom3Temp  << " "
		      << sAtom4Temp << " in res No. " << ResSequence[j].first << std::endl;
	}  
	else
	{
	  std::cerr << "Warning: Improper expected, but some of the atoms are missing: " << sAtom1Temp << " " << sAtom2Temp << " " << sAtom3Temp 
		     << " " << sAtom4Temp << ", residue No. " << ResSequence[j].first << std::endl;
	  bRetVal = false;		// Unknown atom name in *.impr file.
	}
      }
    }
    else
    {
      std::cerr << "Warning: no impropers for residue " << ResSequence[j].second << " have been predefined." << std::endl;
      // Unknown residue is a NOT a reason to break, since impropers are NOT mandatory.
    }
  }
  
  return bRetVal;
}

//---------------------------------------------------------------------------------------------------
bool TMolecule::bMakePosRestraints()
{
  for(unsigned int i = 0; i < Atoms.size(); ++i)
  {
    if(Atoms[i].dPosRestr > 0.)
    {
      double dX0, dY0, dZ0;
      Atoms[i].GetCoords(dX0, dY0, dZ0);
      PosRestrs.push_back(TPosRestr(Atoms[i], dX0, dY0, dZ0, Atoms[i].dPosRestr));
    }
  }
  return true;
}
    
//---------------------------------------------------------------------------------------------------
bool TMolecule::bMakeElectrosAndLennardsWithin()
{
  bool bRetVal = true;

  for(unsigned int i = 0; i < vvContacts.size(); ++i)
    for(unsigned int j = i+1; j < vvContacts[i].size(); ++j)
    {
      if(vvContacts[i][j] >= 3)
      {
	TElectr elToAdd  = elMakeElectrPair(Atoms[i], Atoms[j]);
	TLennard ljToAdd = ljMakeLJPair(    Atoms[i], Atoms[j]); 
	
	if(vvContacts[i][j] == 3)
	{
	  elToAdd.Scale14();
	  ljToAdd.Scale14();
	}
	
	if(TConsts::ELECTR)
	  Electros.push_back(elToAdd);
	
	// Now follows magic ass-factor of 5: (the same one appears in bMakeElectrosAndLennardsExtern()).
	// It gives us some margin when seeking Lennard-Jones-bound atom pairs.
	if(TConsts::VANDERWAALS)
	  if(ljToAdd.dLength_SQR() < TConsts::CUTOFF_SQR + 5.)
	    Lennards.push_back(ljToAdd);	  
      }
    }
  
  return bRetVal;
}
//---------------------------------------------------------------------------------------------------
bool TMolecule::bMakeHBondsWithin()
{
  bool bRetVal = true;

  for(unsigned int i = 0; i < vvContacts.size(); ++i)		// Possible hydrogen donors.
    for(unsigned int j = 0; j < vvContacts[i].size(); ++j)	// Possible hydrogen acceptors.
    {
      if(Textutils::bIsHydrogenDonor(Atoms[i].sGetType())    && 
	 Textutils::bIsHydrogenAcceptor(Atoms[j].sGetType()) && 
	 vvContacts[i][j] == 4)
      {
	for(unsigned int k = 0; k < vvContacts[i].size(); ++k)	// Possible hydrogen 'parents', i.e. P-D...A
	{
	  if(vvContacts[i][k] == 1)
	  {
	    THBond hbTemp(Atoms[k], Atoms[i], Atoms[j]);
	    if(hbTemp.dCalcEnergy() < THBond::dThreshold)	// Energy threshold for a hydrogen bond.
	    {
	      HBonds.push_back(hbTemp);
	      
	      if(TConsts::TALKATIVE > 1)
	      {
		std::clog << "Intrachain hydrogen bond detected: "
			  << Atoms[i].sGetPDB_tag().substr(0,27) << " - "
			  << Atoms[j].sGetPDB_tag().substr(0,27) << std::endl;
	      }
	      
	      break;
	    }
	  }
	}
      }
    }
    
  return bRetVal;
}

//---------------------------------------------------------------------------------------------------
bool TMolecule::bMakeElectrosAndLennardsExtern(TMolecule& _molMol1, TMolecule& _molMol2)
{
  bool bRetVal = true;
  for(unsigned int i = 0; i < _molMol1.Atoms.size(); ++i)
    for(unsigned int j = 0; j < _molMol2.Atoms.size(); ++j)
    {
      if(TConsts::ELECTR) 
	TMolecule::Electros.push_back(elMakeElectrPair(_molMol1.Atoms[i], _molMol2.Atoms[j]));    
					      
      // Now magic constant of 10 follows: (its another occurrence is in bMakeElectrosAndLennardsWithin()).
      // It gives us some margin when defining Lennard-Jones-bound atom pairs.
      if(TConsts::VANDERWAALS) 
	if(dDistance_SQR(_molMol1.Atoms[i], _molMol2.Atoms[j]) < TConsts::CUTOFF_SQR + 10.)
	  TMolecule::Lennards.push_back(ljMakeLJPair(_molMol1.Atoms[i], _molMol2.Atoms[j]));	  
    }
  
  return bRetVal;
}

//---------------------------------------------------------------------------------------------------
bool TMolecule::bMakeHBondsExtern(TMolecule& _molMol1, TMolecule& _molMol2)
{
  bool bRetVal = true;
  for(unsigned int i = 0; i < _molMol1.Atoms.size(); ++i)
    for(unsigned int j = 0; j < _molMol2.Atoms.size(); ++j)
    {
      TAtom *patDonor    = NULL,
	    *patAcceptor = NULL;
      TMolecule* pmolDonor = NULL;
      
      if(     Textutils::bIsHydrogenDonor(_molMol1.Atoms[i].sGetType()))
      {
	patDonor = &(_molMol1.Atoms[i]);
	pmolDonor = &_molMol1;
	patAcceptor = &(_molMol2.Atoms[j]);
      }
      else if(Textutils::bIsHydrogenDonor(_molMol2.Atoms[j].sGetType()))
      {
	patDonor = &(_molMol2.Atoms[j]);
	pmolDonor = &_molMol2;
	patAcceptor = &(_molMol1.Atoms[i]);
      }
      if(patAcceptor != NULL)
	  if(!(Textutils::bIsHydrogenAcceptor(patAcceptor->sGetType()))) patAcceptor = NULL;
      
      // "People always told me: be careful what you do":
      if(patDonor != NULL && patAcceptor != NULL && pmolDonor != NULL && patDonor != patAcceptor)
      {
	// Now find the hydrogen's ``parent'' atom:
	for(std::list<TBond>::iterator K = pmolDonor->Bonds.begin(); K != pmolDonor->Bonds.end(); ++K)
	{
	  TAtom *patParent = NULL;
	  if(     K->patAtom1 == patDonor) patParent = K->patAtom2;
	  else if(K->patAtom2 == patDonor) patParent = K->patAtom1;
	    
	  if(patParent != NULL)
	  {
	    THBond hbTemp(patParent, patDonor, patAcceptor);
	    if(hbTemp.dCalcEnergy() < THBond::dThreshold)	// Criterion for an H-bond to exist.
	    {
	      TMolecule::HBonds.push_back(hbTemp);

	      if(TConsts::TALKATIVE > 2)
	      {
		std::clog << "Interchain hydrogen bond detected: "
			  << patDonor   ->sGetPDB_tag().substr(0,27) << " - "
			  << patAcceptor->sGetPDB_tag().substr(0,27) << std::endl;
	      }
	      
	      break;
	    }
	  }
	}
      }
    }
  
  return bRetVal;
}

//---------------------------------------------------------------------------------------------------
bool TMolecule::bMakeSprings(TMolecule& _molMol1, TMolecule& _molMol2)
{
  bool bRetVal = true;
  for(unsigned int i = 0; i < _molMol1.Atoms.size(); ++i)
    for(unsigned int j = 0; j < _molMol2.Atoms.size(); ++j)
    {
      TAtom *pat1 = &(_molMol1.Atoms[i]),
	    *pat2 = &(_molMol2.Atoms[j]);
      if(pat1 == pat2) continue;
      for(unsigned int k = 0; k < TConsts::DistanceRestraintsDescriptors.size(); ++k)
      {
	if(TConsts::DistanceRestraintsDescriptors[k].bAtomFits(pat1) && 
	   TConsts::DistanceRestraintsDescriptors[k].bAtomFits(pat2) &&
	  !TConsts::DistanceRestraintsDescriptors[k].bIsUsed)
	{
	  Springs.push_back(TSpring(pat1, pat2, TConsts::DistanceRestraintsDescriptors[k].bp));
	  TConsts::DistanceRestraintsDescriptors[k].bIsUsed = true;
	  break;
	}
      }
    }
  
  return bRetVal;
}

//---------------------------------------------------------------------------------------------------
bool TMolecule::bMakeBPRestraints(TMolecule& _molMol1, TMolecule& _molMol2)
{
  for(unsigned int k = 0u; k < TConsts::BPRestraintsDescriptors.size(); ++k)
  {
    if(TConsts::BPRestraintsDescriptors[k].bIsUsed) continue;
    std::vector<TAtom*> AtomPtrs;
    
    bool bSthAdded1 = false;
    for(unsigned int i = 0u; i < _molMol1.Atoms.size(); ++i)
    {
      
      if(TConsts::BPRestraintsDescriptors[k].bAtomFits(&_molMol1.Atoms[i]) && 
         Textutils::bIsInBase(_molMol1.Atoms[i].sGetName(), _molMol1.Atoms[i].sGetResName()))
      {
        AtomPtrs.push_back(&_molMol1.Atoms[i]);
        bSthAdded1 = true;
      }
    }
    
    if(!bSthAdded1) continue;
    
    if(&_molMol1 != &_molMol2) // Nonidentical molecules
    {
      for(unsigned int i = 0u; i < _molMol2.Atoms.size(); ++i)
      {
        if(TConsts::BPRestraintsDescriptors[k].bAtomFits(&_molMol2.Atoms[i]) && 
           Textutils::bIsInBase(_molMol2.Atoms[i].sGetName(), _molMol2.Atoms[i].sGetResName()) &&
           AtomPtrs.front()->sGetResIdent() != _molMol2.Atoms[i].sGetResIdent())
        {
          AtomPtrs.push_back(&_molMol2.Atoms[i]);
        }
      }
    }
    
    bool bTwoDifferent = false;
    for(unsigned int i = 0u; i < AtomPtrs.size(); ++i)
    {
      if(AtomPtrs.front()->sGetResIdent() != AtomPtrs[i]->sGetResIdent())
      {
        bTwoDifferent = true;
        break;
      }
    }
   
    if(bTwoDifferent)
    {
      TMolecule::BasePairs.push_back(TBasePair(AtomPtrs));
      TConsts::BPRestraintsDescriptors[k].bIsUsed = true;
    }
  }
    
  return true;
}

//---------------------------------------------------------------------------------------------------
bool TMolecule::bMakeBasePairs()
{
  std::vector<std::vector<THBond*> > PossibleSS;	// Future candidates for secondary struct.
  
  for(std::list<THBond>::iterator I = HBonds.begin(); I != HBonds.end(); ++I)
  {
    bool bAdded = false;
    for(unsigned int j = 0; j < PossibleSS.size(); ++j)
    {
      std::string sRes1Ident     = PossibleSS[j][0]->patDonor->sGetResIdent(),
		  sRes2Ident     = PossibleSS[j][0]->patAcceptor->sGetResIdent(),
		  sCurrRes1Ident = I->patDonor->sGetResIdent(),
		  sCurrRes2Ident = I->patAcceptor->sGetResIdent();
      
      if((sRes1Ident == sCurrRes1Ident && sRes2Ident == sCurrRes2Ident) ||
	 (sRes1Ident == sCurrRes2Ident && sRes2Ident == sCurrRes1Ident))
      {
	if(sCurrRes1Ident != sCurrRes2Ident)
	{
	  PossibleSS[j].push_back(&(*I));
	  bAdded = true;
	}
	else std::cerr << "Internal Warning: attempt to build a trivial base pair." << std::endl;
      }
    }
    if(!bAdded) PossibleSS.push_back(std::vector<THBond*>(1, &(*I)));
  }
  
  for(unsigned int j = 0; j < PossibleSS.size(); ++j)
  {
    if(PossibleSS[j].size() > 1)
    {
      TBasePair bpToAdd(PossibleSS[j][0], PossibleSS[j][1], (PossibleSS[j].size() > 2 ? PossibleSS[j][2] : NULL));
      if(bpToAdd.dCalcEnergy() < TBasePair::dThreshold)
      {
	BasePairs.push_back(bpToAdd);      
	if(TConsts::TALKATIVE > 1) 
	  std::clog << "H-bonded base pair detected: " 
		<< (PossibleSS[j][0]->patDonor->sGetResName() + PossibleSS[j][0]->patDonor->sGetResIdent()) 
		<< "--- "
		<< (PossibleSS[j][0]->patAcceptor->sGetResName() + PossibleSS[j][0]->patAcceptor->sGetResIdent())
		<< std::endl;
      }
    }
  }
  return true;
}

//---------------------------------------------------------------------------------------------------
bool TMolecule::bBuildContacts(std::map<std::string, int>& AtomContentMap) const
{
  bool bRetVal = true;
  
/******************************************************************************************
 * Matrix vvContacts describes relation between pair of atoms (i,j):
 * 0:  i == j  
 * 1:  i and j are separated by one bond 
 * 2:  i and j separated by two bonds;
 * 3:  i and j separated by 3 bonds (so called '1-4 interaction');
 * 4:  i and j are separated by >= 4 bonds, interacting normally.
 *****************************************************************************************
 * If either two cases apply (as 1 and 3 for vicinal carbons in cyclobutane) 
 * then the smaller value is taken (1 in the example of cyclobutane).
 *****************************************************************************************/
  vvContacts = std::vector<std::vector<char> >(Atoms.size(), std::vector<char>(Atoms.size(), 4));

  for(std::list<TDihedral>::const_iterator I = Dihedrals.begin(); I != Dihedrals.end(); ++I)
  {
    if(AtomContentMap.find(I->patAtom1->sGetResIdent()+I->patAtom1->sGetName()) == AtomContentMap.end() ||
       AtomContentMap.find(I->patAtom4->sGetResIdent()+I->patAtom4->sGetName()) == AtomContentMap.end() )
    {  
      bRetVal = false;
      std::cerr << "Internal error: dihedral " 
		<< I->patAtom1->sGetName() << "(" << I->patAtom1->sGetResIdent() << ") -- " 
		<< I->patAtom2->sGetName() << "(" << I->patAtom2->sGetResIdent() << ") -- " 
		<< I->patAtom3->sGetName() << "(" << I->patAtom3->sGetResIdent() << ") -- "
		<< I->patAtom4->sGetName() << "(" << I->patAtom4->sGetResIdent() << ") contains unknown atom(s)." << std::endl;
    }
    else
    {
      vvContacts[AtomContentMap[I->patAtom1->sGetResIdent()+I->patAtom1->sGetName()]]\
		[AtomContentMap[I->patAtom4->sGetResIdent()+I->patAtom4->sGetName()]] = 
      vvContacts[AtomContentMap[I->patAtom4->sGetResIdent()+I->patAtom4->sGetName()]]\
		[AtomContentMap[I->patAtom1->sGetResIdent()+I->patAtom1->sGetName()]] = 3;    
    }
  }

  for(std::list<TAngle>::const_iterator I = Angles.begin(); I != Angles.end(); ++I)
  {
    if(AtomContentMap.find(I->patAtom1->sGetResIdent()+I->patAtom1->sGetName()) == AtomContentMap.end() ||
       AtomContentMap.find(I->patAtom3->sGetResIdent()+I->patAtom3->sGetName()) == AtomContentMap.end() )
    {  
      bRetVal = false;
      std::cerr << "Internal error: angle " 
		<< I->patAtom1->sGetName() << "(" << I->patAtom1->sGetResIdent() << ") -- " 
		<< I->patAtom2->sGetName() << "(" << I->patAtom2->sGetResIdent() << ") -- " 
		<< I->patAtom3->sGetName() << "(" << I->patAtom3->sGetResIdent() << ") contains unknown atom(s)." << std::endl;
    }
    else
    {
      vvContacts[AtomContentMap[I->patAtom1->sGetResIdent()+I->patAtom1->sGetName()]]\
		[AtomContentMap[I->patAtom3->sGetResIdent()+I->patAtom3->sGetName()]] = 
      vvContacts[AtomContentMap[I->patAtom3->sGetResIdent()+I->patAtom3->sGetName()]]\
		[AtomContentMap[I->patAtom1->sGetResIdent()+I->patAtom1->sGetName()]] = 2;
    }
  }

  for(std::list<TBond>::const_iterator I = Bonds.begin(); I != Bonds.end(); ++I)
  {
    if(AtomContentMap.find(I->patAtom1->sGetResIdent()+I->patAtom1->sGetName()) == AtomContentMap.end() ||
       AtomContentMap.find(I->patAtom2->sGetResIdent()+I->patAtom2->sGetName()) == AtomContentMap.end() )
    {  
      bRetVal = false;
      std::cerr << "Internal error: bond " 
		<< I->patAtom1->sGetName() << "(" << I->patAtom1->sGetResIdent() << ") -- " 
		<< I->patAtom2->sGetName() << "(" << I->patAtom2->sGetResIdent() << ") contains unknown atom(s)." << std::endl;
    }
    else
    {
      vvContacts[AtomContentMap[I->patAtom1->sGetResIdent()+I->patAtom1->sGetName()]]\
		[AtomContentMap[I->patAtom2->sGetResIdent()+I->patAtom2->sGetName()]] = 
      vvContacts[AtomContentMap[I->patAtom2->sGetResIdent()+I->patAtom2->sGetName()]]\
		[AtomContentMap[I->patAtom1->sGetResIdent()+I->patAtom1->sGetName()]] = 1;
    }
  }
  
  for(unsigned int i = 0; i < vvContacts.size(); ++i)
  vvContacts[i][i] = 0;
  
  return bRetVal;
}

//---------------------------------------------------------------------------------------------------
void TMolecule::ClearContacts() const
{// Release some memory.
  vvContacts.clear();
}

//---------------------------------------------------------------------------------------------------
bool TMolecule::bRecognizeResidues()
{
  bool bRetVal = true;
  for(unsigned int i = 0; i < Atoms.size(); ++i)
    if (!(Atoms[i].bRecognizeResidue()))
      bRetVal = false;
  
  return bRetVal;
}

//---------------------------------------------------------------------------------------------------
bool TMolecule::bRecognizeNames()
{// Recognize all atoms' names given the sPDB_tag.
  for(unsigned int i = 0; i < Atoms.size(); ++i)
    Atoms[i].bRecognizeName();
  return true;
}

//---------------------------------------------------------------------------------------------------
bool TMolecule::bRecognizeAtoms()
{
  std::vector<TAtom> OnlyGoodAtoms;
  bool bRetVal = true;
  for(unsigned int i = 0; i < Atoms.size(); ++i)
    if (Atoms[i].bRecognize())
    {
      OnlyGoodAtoms.push_back(Atoms[i]);
    }
    else 
    {
      bRetVal = false;
      std::cerr << "Error: recognition of atom " << Atoms[i].sGetPDB_tag().substr(0,27) << "... failed. Atom erased." << std::endl;
    }
    
  Atoms = OnlyGoodAtoms;
  return bRetVal;  
}

//---------------------------------------------------------------------------------------------------
void TMolecule::CalcBornRadii()
{// At development stage - use at your own risk
  for(unsigned int i = 0; i < Atoms.size(); ++i)
    Atoms[i].dBornRadius = 1./Atoms[i].dGetvdWRadius();
  
  for(unsigned int i = 0; i < Atoms.size(); ++i)
    for(unsigned int j = i; j < Atoms.size(); ++j)
    {
       double dAj       = Atoms[j].dGetvdWRadius(),		// van der Waals radius of j-th atom.
	      dLen_SQR  = dDistance_SQR(Atoms[i], Atoms[j]),	// Square of distance between atoms i and j.
	      dLen      = sqrt(dLen_SQR),			// Distance itsef from i to j.
	      dVal;						// The correction being introduced.
	      
       if(Atoms[i].dGetvdWRadius() + dAj > dLen) // Overlap of spheres - some correction is needed (factor of 1/2).
	 dVal = dAj/(4.*(dLen_SQR-dAj*dAj)) - 1./(8.*dLen)*log((dLen-dAj)/(dLen+dAj));
       else					// No overlap - ordinary formulae apply.
	 dVal = dAj/(2.*(dLen_SQR-dAj*dAj)) - 1./(4.*dLen)*log((dLen-dAj)/(dLen+dAj));

       Atoms[i].dBornRadius -= dVal;
       Atoms[j].dBornRadius -= dVal;
    }

  for(unsigned int i = 0; i < Atoms.size(); ++i)
    if(TConsts::bIsNan(Atoms[i].dBornRadius) || Atoms[i].dBornRadius == 0.) 
      Atoms[i].dBornRadius = Atoms[i].dGetvdWRadius();
    else Atoms[i].dBornRadius = 1./Atoms[i].dBornRadius;
}

//---------------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------
bool TMolecule::bRead(std::istream& _isIn)
{
  if(!iReadAtoms(_isIn)) 
    return false;
  
  if(!bRecognizeResidues())
  {
    std::cerr << "Error: unknown residue(s). Molecule discarded." << std::endl;
    Atoms.clear();
    return false;
  }
  
  bRecognizeNames();
  bCorrectTerminalResidueNames();
  bRecognizeAtoms();
  bAddMissingAtoms();
  CorrectAtomNumbers();

/******************************************************************************************
 * Now the AtomContentMap is created to simplify the search for pointers to given atoms.
 * The index of this map is a std::string of the form: <residuenumber>+<atomname>.
 * <residuenumber> is given by TAtom's method sGetResIdent().
 * The insertion code is included in <residuenumber> if applicable.
 *****************************************************************************************/
  std::map<std::string, int> AtomContentMap;
  
/******************************************************************************************
 * The following vector of pairs stores the residue identity (i.e. chain letter + 
 * res. number) and the residue name. Will be of use when creating bonds, angles,
 * impopers and dihedrals.
 *
 * Both references to ResSequence and AtomContentMap are passed to appropriate functions
 * (bMakeBonds, bMakeAngles, etc.)
 *****************************************************************************************/
  std::vector<std::pair<std::string,std::string> > ResSequence;
  if(!(Atoms.empty())) ResSequence.push_back(std::make_pair(Atoms[0].sGetResIdent(),Atoms[0].sGetResName()));
  for(unsigned int i = 0; i < Atoms.size(); ++i)
  {
    if(AtomContentMap.find(Atoms[i].sGetResIdent()+Atoms[i].sGetName()) == AtomContentMap.end()) 
      AtomContentMap[Atoms[i].sGetResIdent()+Atoms[i].sGetName()] = i;
    else std::cerr << "Error: coincidence of atom names: in residue " << Atoms[i].sGetResIdent() 
		   << " name " << Atoms[i].sGetName() << " is ambiguous." << std::endl;
    
    if(Atoms[i].sGetResIdent() != ResSequence.back().first)
      ResSequence.push_back(std::make_pair(Atoms[i].sGetResIdent(), Atoms[i].sGetResName()));
  }
  std::clog << "Number of atoms:    " << Atoms.size() << std::endl
	    << "Number of residues: " << ResSequence.size() << std::endl;
  
  bool bRetIntraBonds     = bMakeIntraBonds(    AtomContentMap, ResSequence),
       bRetIntraAngles    = bMakeIntraAngles(   AtomContentMap, ResSequence),
       bRetIntraDihedrals = bMakeIntraDihedrals(AtomContentMap, ResSequence),
       bRetIntraImpropers = bMakeIntraImpropers(AtomContentMap, ResSequence),
       bRetInterBonds     = bMakeInterBonds(    AtomContentMap, ResSequence),
       bRetInterAngles    = bMakeInterAngles(   AtomContentMap, ResSequence),
       bRetInterDihedrals = bMakeInterDihedrals(AtomContentMap, ResSequence),
       bRetIntraPosRestr  = bMakePosRestraints();
  
  bBuildContacts(AtomContentMap);	// Will be deleted during creation of nonbonds 'within'.

  if(TConsts::TALKATIVE > 0)
  std::clog << "Number of bonds built:     " << Bonds.size()     << std::endl
	    << "Number of angles built:    " << Angles.size()    << std::endl
	    << "Number of dihedrals built: " << Dihedrals.size() << std::endl 
	    << "Number of impropers built: " << Impropers.size() << std::endl  
            << Textutils::pcBar1 << "End of molecule (chain)" << std::endl << Textutils::pcBar1;

	    
  return bRetIntraBonds && bRetIntraAngles && bRetIntraDihedrals && bRetIntraImpropers &&
	 bRetInterBonds && bRetInterAngles && bRetInterDihedrals && bRetIntraPosRestr;
}

//---------------------------------------------------------------------------------------------------
bool TMolecule::bWrite(std::ostream& _osOut_) const
{
  bool bRetVal = true;
  for(unsigned int i = 0; i < Atoms.size(); ++i)
    if(!Atoms[i].bWrite(_osOut_)) bRetVal = false;
  _osOut_ << "TER" << std::endl;
  return bRetVal;
}

//---------------------------------------------------------------------------------------------------
bool TMolecule::bEmpty() const
{
  return Atoms.empty();
}

//---------------------------------------------------------------------------------------------------
unsigned int TMolecule::uGetNumberOfAtoms() const
{
  return Atoms.size();
}
//---------------------------------------------------------------------------------------------------
unsigned int TMolecule::uGetNumberOfElectros()
{
  return Electros.size();
}

//---------------------------------------------------------------------------------------------------
unsigned int TMolecule::uGetNumberOfLennards()
{
  return Lennards.size();
}

//---------------------------------------------------------------------------------------------------
unsigned int TMolecule::uGetNumberOfHBonds()
{
  return HBonds.size();
}

//---------------------------------------------------------------------------------------------------
unsigned int TMolecule::uGetNumberOfSprings()
{
  return Springs.size();
}

//---------------------------------------------------------------------------------------------------
unsigned int TMolecule::uGetNumberOfPosRestrs()
{
  unsigned int uConstraints = 0u;
  for(unsigned int i = 0u; i < Atoms.size(); ++i)
    if(Atoms[i].dPosRestr < 0.) ++uConstraints;
  return PosRestrs.size() + uConstraints;
}

//---------------------------------------------------------------------------------------------------
unsigned int TMolecule::uGetNumberOfBasePairs()
{
  return BasePairs.size();
}

//---------------------------------------------------------------------------------------------------
int TMolecule::iFillMap(std::map<std::string, TAtom>& _ResPDBMap_, const std::string& _sResName)
{
  int iRetVal = 0;
  std::string sLine;
  
  std::ifstream ifsIn((TConsts::PDBDIR+_sResName+".pdb").c_str());
  if(ifsIn != NULL)
  {
    while(std::getline(ifsIn, sLine))	//read single line from PDB
    {
      TAtom atTemp;
      if(atTemp.bAssignTag(sLine))
      {					//it fits as a description of atom
	if(atTemp.bRecognize())
	{
	  if(_ResPDBMap_.find(atTemp.sGetName()) == _ResPDBMap_.end())  //Occupancy handled herein.
	  {
	    _ResPDBMap_[atTemp.sGetName()]=atTemp;
	    ++iRetVal;
	  }
	  else std::cerr << "Internal Warning: ambiguous atom description: " << sLine << (" in "+TConsts::PDBDIR+_sResName+".pdb") << std::endl;
	}
	else 
	  std::cerr << "Internal Warning: Unresolved atom description: " << sLine << (" in "+TConsts::PDBDIR+_sResName+".pdb") << std::endl;
      }
      else if(sLine.length() >= 3 && sLine.substr(0,3) == "TER") break;
    }
  }
  else std::cerr << ("Internal Warning: File "+TConsts::PDBDIR+_sResName+".pdb not found.") << std::endl;
  
  return iRetVal;
}

//---------------------------------------------------------------------------------------------------
bool TMolecule::bCompleteResidue(std::vector<TAtom>& AtomsNew, 
				 std::map<std::string, int>& CurrentResContent,
				 std::map<std::string, TAtom>& CurrentResPDB,
				 const std::string& sCurrentResName,
				 const std::string& sCurrentResIdent)
{
  bool bRetVal = true;

  std::ifstream ifsIn((TConsts::DIHEDDIR+sCurrentResName+".dihd").c_str());
  if(ifsIn != NULL)
  {
    std::string sDihLine;
    while(std::getline(ifsIn, sDihLine))
    {
      std::istringstream issIn(sDihLine);
      std::string sAtom1Name, sAtom2Name, sAtom3Name, sAtom4Name;
      issIn >> sAtom1Name >> sAtom2Name >> sAtom3Name >> sAtom4Name;

      Textutils::bCorrectAtomName(sAtom1Name);
      Textutils::bCorrectAtomName(sAtom2Name);
      Textutils::bCorrectAtomName(sAtom3Name);
      Textutils::bCorrectAtomName(sAtom4Name);
      
      if(CurrentResPDB.find(sAtom1Name) == CurrentResPDB.end() ||
	 CurrentResPDB.find(sAtom2Name) == CurrentResPDB.end() ||
	 CurrentResPDB.find(sAtom3Name) == CurrentResPDB.end() ||
	 CurrentResPDB.find(sAtom4Name) == CurrentResPDB.end() )
      {
	std::cerr << "Internal Warning: cannot find atom(s): "
		  << (sAtom1Name+", "+sAtom2Name+", "+sAtom3Name+" or "+sAtom4Name)
		  << (" in internal "+TConsts::PDBDIR+sCurrentResName+".pdb file") << std::endl;
	continue;
      }
      
      unsigned int uMatch = 0;	//number of atoms in given dihedral that can be found in current residue
      bool bMatch1 = false, bMatch2 = false, bMatch3 = false, bMatch4 = false;
      if(CurrentResContent.find(sAtom1Name) != CurrentResContent.end())
      {
	bMatch1 = true;
	++uMatch;
      }
      if(CurrentResContent.find(sAtom2Name) != CurrentResContent.end())
      {
	bMatch2 = true;
	++uMatch;
      }
      if(CurrentResContent.find(sAtom3Name) != CurrentResContent.end())
      {
	bMatch3 = true;
	++uMatch;
      }
      if(CurrentResContent.find(sAtom4Name) != CurrentResContent.end())
      {
	bMatch4 = true;
	++uMatch;
      }
      
      if(uMatch == 3)
      {//Exactly 3 atoms match and one does not
	if(!bMatch2) sAtom2Name.swap(sAtom1Name);	//
	if(!bMatch3) sAtom3Name.swap(sAtom1Name);	// We want the first (sAtom1Temp) atom to be the one to add
	if(!bMatch4) sAtom4Name.swap(sAtom1Name);	//
	
	TAtom atAtomToAdd(CurrentResPDB[sAtom1Name]);
	atAtomToAdd.bTransform(CurrentResPDB[sAtom2Name], CurrentResPDB[sAtom3Name], CurrentResPDB[sAtom4Name],
			       AtomsNew[CurrentResContent[sAtom2Name]], AtomsNew[CurrentResContent[sAtom3Name]], AtomsNew[CurrentResContent[sAtom4Name]]);
	if(!(atAtomToAdd.bSetResIdent(sCurrentResIdent)))
	  std::cerr << "Error: Something is wrong with residue numbering." << std::endl;
			       
	CurrentResContent[sAtom1Name] = AtomsNew.size();
	AtomsNew.push_back(atAtomToAdd);
	if(TConsts::TALKATIVE > 0) 
	{
	  std::clog << "Missing atom added: ";
	  atAtomToAdd.bWrite(std::clog);
	}
      }
    }
  }
  else
  {
    std::cerr << "Internal Warning: cannot read file: " << (TConsts::DIHEDDIR+sCurrentResName+".dihd") << std::endl;
    bRetVal = false;
  }
    
  return bRetVal;
}

//---------------------------------------------------------------------------------------------------
bool TMolecule::bAddMissingAtoms()
{
  bool bRetVal = true;
  std::vector<TAtom> AtomsNew;
  std::map<std::string, int> CurrentResContent;
  std::string sCurrentResIdent = Atoms.front().sGetResIdent(), 
	      sCurrentResName  = Atoms.front().sGetResName();
  for(unsigned int i = 0; i < Atoms.size(); ++i)
  {
    if(sCurrentResIdent == Atoms[i].sGetResIdent())
    {
      if(CurrentResContent.find(Atoms[i].sGetName()) == CurrentResContent.end())
      {// By the way we got a nice handling of occupancy: only first ocurrence of atom counts.
	CurrentResContent[Atoms[i].sGetName()] = AtomsNew.size();
	AtomsNew.push_back(Atoms[i]);
      }
    }
    if(sCurrentResIdent != Atoms[i].sGetResIdent() || i+1 == Atoms.size())
    {// New residue began. First, complete the previous one...
      std::map<std::string, TAtom> CurrentResPDB;
      if(iFillMap(CurrentResPDB, sCurrentResName))
      {
	bCompleteResidue(AtomsNew, CurrentResContent, CurrentResPDB, sCurrentResName, sCurrentResIdent);
      }
      else std::cerr << "Warning: prospective completion of residue " << sCurrentResName << " not possible." << std::endl;
      CurrentResContent.clear();
      
      //...now we can start the new one.
      sCurrentResIdent = Atoms[i].sGetResIdent();
      sCurrentResName  = Atoms[i].sGetResName();
      if(i+1 < Atoms.size()) i--;
    }
  }
  
  Atoms = AtomsNew;

  if(bIsPhosphorusAt5Term() && !bIsOP3At5Term()) bAddOP3(); // Common error - missing OP3 - needs special treatment

  return bRetVal;
}

//---------------------------------------------------------------------------------------------------
bool TMolecule::bAddOP3()
{
  const std::string s5TermResIdent = Atoms.front().sGetResIdent();

  Atoms.insert(Atoms.begin(), TAtom()); // From now on Atoms[0] will be the OP3
  TAtom *patP       = NULL, 
	*patOP1     = NULL, 
	*patOP2     = NULL, 
	*patO5Prime = NULL, 
	*patC5Prime = NULL, 
	*patOP3     = &Atoms[0];
	
  for(unsigned int i = 1; i < Atoms.size(); ++i)
    if(s5TermResIdent == Atoms[i].sGetResIdent())
    {
      if(     Atoms[i].sGetName() == "P"  ) patP       = &Atoms[i];
      else if(Atoms[i].sGetName() == "OP1") patOP1     = &Atoms[i];
      else if(Atoms[i].sGetName() == "OP2") patOP2     = &Atoms[i];
      else if(Atoms[i].sGetName() == "O5'") patO5Prime = &Atoms[i];
      else if(Atoms[i].sGetName() == "C5'") patC5Prime = &Atoms[i];
    }
    else break;
  
  if(patP == NULL || patOP1 == NULL || patOP2 == NULL || patO5Prime == NULL || patC5Prime == NULL)
  {
    Atoms.erase(Atoms.begin());		// Remove what have ye just added
    return false;
  }
  
  TVector v1 = vecMakeVector(patP, patOP1),
	  v2 = vecMakeVector(patP, patOP2),
	  v3 = vecMakeVector(patP, patO5Prime);
	  
  *patOP3 = *patOP2;
  patOP3->SetName("OP3");  
  *patOP3 -= v2 + (v1 + v2 + v3);

  Bonds.push_back(TBond(patOP3, patP, 
			    TConsts::GeneralBondsData[TConsts::AtomTypeNumbers["O2"]]\
						     [TConsts::AtomTypeNumbers["P" ]]));
  Angles.push_back(TAngle(patOP3, patP, patOP1, 
			    TConsts::GeneralAnglesData[TConsts::AtomTypeNumbers["O2"]]\
						      [TConsts::AtomTypeNumbers["P" ]]\
						      [TConsts::AtomTypeNumbers["O2"]]));
  Angles.push_back(TAngle(patOP3, patP, patOP2, 
			    TConsts::GeneralAnglesData[TConsts::AtomTypeNumbers["O2"]]\
						      [TConsts::AtomTypeNumbers["P" ]]\
						      [TConsts::AtomTypeNumbers["O2"]]));
  Angles.push_back(TAngle(patOP3, patP, patO5Prime, 
			    TConsts::GeneralAnglesData[TConsts::AtomTypeNumbers["O2"]]\
						      [TConsts::AtomTypeNumbers["P" ]]\
						      [TConsts::AtomTypeNumbers["OS"]]));
						      
  Dihedrals.push_back(TDihedral(patOP3, patP, patO5Prime, patC5Prime, 
		      *TConsts::GeneralDihedralsData[TConsts::AtomTypeNumbers["O2"]]\
						    [TConsts::AtomTypeNumbers["P" ]]\
						    [TConsts::AtomTypeNumbers["OS"]]\
						    [TConsts::AtomTypeNumbers["CT"]].begin()));
   // Hopefully O2-P-OS-CT dihedral has just one Fourier component
   
  if(TConsts::TALKATIVE > 0)
  {
    std::clog << "Missing atom added: ";
    patOP3->bWrite(std::clog);
  }
  return true;
}

//---------------------------------------------------------------------------------------------------
bool TMolecule::bIsPhosphorusAt5Term() const
{//is there a phosphate group at 5'-terminus?
  bool bRetVal = false;
  const std::string s5TermResIdent = Atoms.front().sGetResIdent();
  for(unsigned int i = 0; i < Atoms.size(); ++i)
    if(s5TermResIdent == Atoms[i].sGetResIdent())
    {
      if(Atoms[i].sGetName()=="P")
      {
	bRetVal = true;
	break;
      }
    }
    else break;
  
  return bRetVal;
}

//---------------------------------------------------------------------------------------------------
bool TMolecule::bIsOP3At5Term() const
{// Self explanatory. Note that it does not apply if 5'-end is capped!
  bool bRetVal = false;
//  if(Atoms.front().sGetResName() == )
  const std::string s5TermResIdent = Atoms.front().sGetResIdent();
  for(unsigned int i = 0; i < Atoms.size(); ++i)
    if(s5TermResIdent == Atoms[i].sGetResIdent())
    {
      if(Atoms[i].sGetName()=="OP3")
      {
	bRetVal = true;
	break;
      }
    }
    else break;
  
  return bRetVal;  
}
//---------------------------------------------------------------------------------------------------
bool TMolecule::bCorrectTerminalResidueNames()
{// Method that corrects the exterior residues' names (i.e. at 5' and 3').
  bool bRetVal = true;

  const std::string s5TermResIdent = Atoms.front().sGetResIdent();
  const std::string s3TermResIdent = Atoms.back().sGetResIdent();  
  
  if(s5TermResIdent == s3TermResIdent && !bIsPhosphorusAt5Term())
  {// Only one residue - change its name to RXN. (UNTESTED)
    for(unsigned int i = 0; i < Atoms.size(); ++i)
    {
      bool bCorrected = false;
      if(s5TermResIdent != Atoms[i].sGetResIdent())
      {
	std::cerr << "Warning: Single nucleotide expected in residue: " << s5TermResIdent << std::endl;
	break;
      }
      
      if(Atoms[i].sGetResName() == "A")
      {
	Atoms[i].bSetResName("RAN");
	bCorrected = true;
      }
      else if(Atoms[i].sGetResName() == "C")
      {
	Atoms[i].bSetResName("RCN");
	bCorrected = true;
      }
      else if(Atoms[i].sGetResName() == "G")
      {
	Atoms[i].bSetResName("RGN");
	bCorrected = true;
      }
      else if(Atoms[i].sGetResName() == "U")
      {
	Atoms[i].bSetResName("RUN");
	bCorrected = true;
      }
      
      if(bCorrected && TConsts::TALKATIVE > 0)
	  std::clog << "Atom's "<< i+1 << " " << Atoms[i].sGetName() << " residue renamed to " << Atoms[i].sGetResName() << std::endl;
    }    
  }
  else
  {
    if(bIsPhosphorusAt5Term())
    {// Case: phosphate at the 5'-terminus. 
      for(unsigned int i = 0; i < Atoms.size(); ++i)
      {
	bool bCorrected5 = false;
	if(s5TermResIdent != Atoms[i].sGetResIdent()) break;
	
	if(Atoms[i].sGetResName() == "RA5")
	{
	  Atoms[i].bSetResName("A");
	  bCorrected5 = true;
	}
	else if(Atoms[i].sGetResName() == "RC5")
	{
	  Atoms[i].bSetResName("C");
	  bCorrected5 = true;
	}
	else if(Atoms[i].sGetResName() == "RG5")
	{
	  Atoms[i].bSetResName("G");
	  bCorrected5 = true;
	}
	else if(Atoms[i].sGetResName() == "RU5")
	{
	  Atoms[i].bSetResName("U");
	  bCorrected5 = true;
	}
	
	if(bCorrected5 && TConsts::TALKATIVE > 0)
	    std::clog << "Atom's "<< i+1 << " " << Atoms[i].sGetName() << " residue renamed to " << Atoms[i].sGetResName() << std::endl;
      }
    }
    else
    {// Case: HO- at the 5'-terminus (H5T)
      for(unsigned int i = 0; i < Atoms.size(); ++i)
      {
	bool bCorrected5 = false;
	if(s5TermResIdent != Atoms[i].sGetResIdent()) break;
	
	if(Atoms[i].sGetResName() == "A")
	{
	  Atoms[i].bSetResName("RA5");
	  bCorrected5 = true;
	}
	else if(Atoms[i].sGetResName() == "C")
	{
	  Atoms[i].bSetResName("RC5");
	  bCorrected5 = true;
	}
	else if(Atoms[i].sGetResName() == "G")
	{
	  Atoms[i].bSetResName("RG5");
	  bCorrected5 = true;
	}
	else if(Atoms[i].sGetResName() == "U")
	{
	  Atoms[i].bSetResName("RU5");
	  bCorrected5 = true;
	}
	
	if(bCorrected5 && TConsts::TALKATIVE > 0)
	    std::clog << "Atom's "<< i+1 << " " << Atoms[i].sGetName() << " residue renamed to " << Atoms[i].sGetResName() << std::endl;
      }
    }

    for(int i = Atoms.size()-1; i >= 0; --i)
    {
	bool bCorrected3 = false;
	if(s3TermResIdent != Atoms[i].sGetResIdent()) break;
	
	if(Atoms[i].sGetResName() == "A" || Atoms[i].sGetResName() == "RA5")
	{
	  Atoms[i].bSetResName("RA3");
	  bCorrected3 = true;
	}
	else if(Atoms[i].sGetResName() == "C" || Atoms[i].sGetResName() == "RC5")
	{
	  Atoms[i].bSetResName("RC3");
	  bCorrected3 = true;
	}
	else if(Atoms[i].sGetResName() == "G" || Atoms[i].sGetResName() == "RG5")
	{
	  Atoms[i].bSetResName("RG3");
	  bCorrected3 = true;
	}
	else if(Atoms[i].sGetResName() == "U" || Atoms[i].sGetResName() == "RU5")
	{
	  Atoms[i].bSetResName("RU3");
	  bCorrected3 = true;
	}
	
	if(bCorrected3 && TConsts::TALKATIVE > 0)
	    std::clog << "Atom's "<< i+1 << " " << Atoms[i].sGetName() << " residue renamed to " << Atoms[i].sGetResName() << std::endl;
    }
  }
  return bRetVal;
}

//---------------------------------------------------------------------------------------------------
void TMolecule::CorrectAtomNumbers()
{// To be used only once during the program execution! Otherwise static variable must be reset.
  static int i = 0;
  for(unsigned int j = 0; j < Atoms.size(); ++j, ++i)
    Atoms[j].iNum = i+1;
}

//---------------------------------------------------------------------------------------------------
bool bMakeNonbondsBetween(TMolecule& _molMol1, TMolecule& _molMol2)
{
  bool bRetVal;
  if(&_molMol1 != &_molMol2)
  {
    bool bRetElectrosAndLennardsExtern = TMolecule::bMakeElectrosAndLennardsExtern(_molMol1, _molMol2),
	 bRetHBondsExtern 	       = (TConsts::HBONDS ? TMolecule::bMakeHBondsExtern(_molMol1, _molMol2) : true);
    bRetVal = bRetElectrosAndLennardsExtern && bRetHBondsExtern;
  }
  else 
  {
    bool bRetElectrosAndLennardsWithin = _molMol1.bMakeElectrosAndLennardsWithin(),
	 bRetHBondsWithin              = (TConsts::HBONDS ? _molMol1.bMakeHBondsWithin() : true);
    _molMol1.ClearContacts();
    bRetVal = bRetElectrosAndLennardsWithin && bRetHBondsWithin;
  }
  
  bool bRetSprings  = TMolecule::bMakeSprings(_molMol1, _molMol2),
       bRetBPRestrs = TMolecule::bMakeBPRestraints(_molMol1, _molMol2);
  return bRetVal && bRetSprings && bRetBPRestrs;
}

//---------------------------------------------------------------------------------------------------
bool bAddSpecificBasePair(TMolecule& _molMol1, const std::string& _sResIdent1, 
			  TMolecule& _molMol2, const std::string& _sResIdent2)
{
  std::vector<TAtom*> AtomPtrs;
  
  for(unsigned int i = 0; i < _molMol1.Atoms.size(); ++i)
  {
    if(_molMol1.Atoms[i].sGetResIdent() == _sResIdent1 && 
        Textutils::bIsInBase(_molMol1.Atoms[i].sGetName(), _molMol1.Atoms[i].sGetResName()))
      AtomPtrs.push_back(&_molMol1.Atoms[i]);
  }
  for(unsigned int i = 0; i < _molMol2.Atoms.size(); ++i)
  {
    if(_molMol2.Atoms[i].sGetResIdent() == _sResIdent2 && 
        Textutils::bIsInBase(_molMol2.Atoms[i].sGetName(), _molMol2.Atoms[i].sGetResName()))
      AtomPtrs.push_back(&_molMol2.Atoms[i]);    
  }
    
  TMolecule::BasePairs.push_back(TBasePair(AtomPtrs));
  return !AtomPtrs.empty();
}

//---------------------------------------------------------------------------------------------------
bool TMolecule::bImposeVienna(const std::string& _Vienna)
{
  bool bRetVal = true;
  std::stack<std::string> Stack;
  std::vector<std::string> ResSequence;
  
  std::string sCurrentResIdent = Atoms.front().sGetResIdent();
  ResSequence.push_back(sCurrentResIdent);
  for(unsigned int i = 0; i < Atoms.size(); ++i)
  {
    if(Atoms[i].sGetResIdent() != sCurrentResIdent)
    {
      sCurrentResIdent = Atoms[i].sGetResIdent();
      ResSequence.push_back(sCurrentResIdent);
    }
  }
  
  if((bRetVal = (_Vienna.length() == ResSequence.size())))
  {
    for(unsigned int i = 0; i < _Vienna.length(); ++i)
    {
	if(_Vienna[i] == '(') Stack.push(ResSequence[i]);
	if(_Vienna[i] == ')')
	{
	  bAddSpecificBasePair(*this, ResSequence[i], *this, Stack.top());
	  Stack.pop();
	}
    }
  }
  
  return bRetVal;
}

//---------------------------------------------------------------------------------------------------
void TMolecule::CalcCovalentGradients()
{
  for(std::list<TBond>::iterator I = Bonds.begin(); I != Bonds.end(); ++I)
    I->CalcGradients();

  for(std::list<TAngle>::iterator I = Angles.begin(); I != Angles.end(); ++I)
    I->CalcGradients();

  for(std::list<TDihedral>::iterator I = Dihedrals.begin(); I != Dihedrals.end(); ++I)
    I->CalcGradients();
  
  for(std::list<TImproper>::iterator I = Impropers.begin(); I != Impropers.end(); ++I)
    I->CalcGradients();

  for(std::list<TPosRestr>::const_iterator I = PosRestrs.begin(), end = PosRestrs.end(); I != end; ++I)
    I->CalcGradients();
}

//---------------------------------------------------------------------------------------------------
bool TMolecule::bClearGradients()
{
  for(unsigned int i = 0; i < Atoms.size(); ++i)
    Atoms[i].vecEnergyGradient = TVector(0., 0., 0.);  
  return true;
}

//---------------------------------------------------------------------------------------------------
double TMolecule::dCalcCovalentEnergy(double& _dMolEnergy_) const
{
  _dMolEnergy_ = dCalcBondsEnergy() + dCalcAnglesEnergy() + dCalcDihedralsEnergy() + dCalcImpropersEnergy();
  return _dMolEnergy_ + dCalcPosRestrEnergy();
}

//---------------------------------------------------------------------------------------------------
double TMolecule::dCalcBondsEnergy() const
{
  double dBondsEnergy  = 0.;
  for(std::list<TBond>::const_iterator I = Bonds.begin(); I != Bonds.end(); ++I)
    dBondsEnergy += I->dCalcEnergy();
  return dBondsEnergy;
}

//---------------------------------------------------------------------------------------------------
double TMolecule::dCalcAnglesEnergy() const
{  
  double dAnglesEnergy = 0.;
  for(std::list<TAngle>::const_iterator I = Angles.begin(); I != Angles.end(); ++I)
    dAnglesEnergy  += I->dCalcEnergy();
  return dAnglesEnergy;
}

//---------------------------------------------------------------------------------------------------
double TMolecule::dCalcDihedralsEnergy() const
{
  double dDihedsEnergy = 0.;
  for(std::list<TDihedral>::const_iterator I = Dihedrals.begin(); I != Dihedrals.end(); ++I)
  {
    double dSingleDihedEnergy = I->dCalcEnergy();
    dDihedsEnergy += (TConsts::bIsNan(dSingleDihedEnergy) ? 0. : dSingleDihedEnergy);	// Naiive protection against numerical instabilities.
  }
  return dDihedsEnergy;
}

//---------------------------------------------------------------------------------------------------
double TMolecule::dCalcImpropersEnergy() const
{
  double dImprosEnergy = 0.;
  for(std::list<TImproper>::const_iterator I = Impropers.begin(); I != Impropers.end(); ++I)
  {
    double dSingleImproEnergy = I->dCalcEnergy();
    dImprosEnergy += (TConsts::bIsNan(dSingleImproEnergy) ? 0. : dSingleImproEnergy);
  }
  return dImprosEnergy;
}

//---------------------------------------------------------------------------------------------------
double TMolecule::dCalcPosRestrEnergy() const
{
  double dPosRestrEnergy = 0.;
  for(std::list<TPosRestr>::const_iterator I = PosRestrs.begin(), end = PosRestrs.end(); I != end; ++I)
  {
    double dSinglePREnergy = I->dCalcEnergy();
    dPosRestrEnergy += (TConsts::bIsNan(dSinglePREnergy) ? 0. : dSinglePREnergy);
  }
  return dPosRestrEnergy;
}

//---------------------------------------------------------------------------------------------------
void TMolecule::Shift(double _dFactor)
{
  for(unsigned int i = 0; i < Atoms.size(); ++i)
    if(Atoms[i].dPosRestr > -1e-4) Atoms[i] += _dFactor*Atoms[i].vecEnergyGradient;
}

//---------------------------------------------------------------------------------------------------
double TMolecule::dAllVectorsNorm()
{
  double dTotalNorm = 0.;	//Norm of the overall 3N-dimensional shift vector for given molecule
  for(unsigned int i = 0; i < Atoms.size(); ++i)
    dTotalNorm += Atoms[i].vecEnergyGradient.dLength_SQR();
  return dTotalNorm;
}

//---------------------------------------------------------------------------------------------------
bool TMolecule::bCalcNonbondGradients()
{
  for(std::list<TElectr>::const_iterator I = Electros.begin(), end = Electros.end(); I != end; ++I)
    I->CalcGradients();

  for(std::list<TLennard>::const_iterator I = Lennards.begin(), end = Lennards.end(); I != end; ++I)
    I->CalcGradients();

  for(std::list<THBond>::const_iterator I = HBonds.begin(), end = HBonds.end(); I != end; ++I)
    I->CalcGradients();  
  
  for(std::list<TSpring>::const_iterator I = Springs.begin(), end = Springs.end(); I != end; ++I)
    I->CalcGradients();  

  for(std::list<TBasePair>::iterator I = BasePairs.begin(), end = BasePairs.end(); I != end; ++I)
    I->CalcGradients();  // This one cannot be const-iterated.

  return true;	// Unsupported yet.
}

//---------------------------------------------------------------------------------------------------
#ifdef SEQUENTIAL
  double TMolecule::dCalcNonbondEnergy(double& _dMolEnergy_)
  {
    double dMolEnergy = 0., dRestrEnergy = 0.;
    
    for(std::list<TElectr>::const_iterator I = Electros.begin(), end = Electros.end(); I != end; ++I)
      dMolEnergy += I->dCalcEnergy();

    for(std::list<TLennard>::const_iterator I = Lennards.begin(), end = Lennards.end(); I != end; ++I)
      dMolEnergy += I->dCalcEnergy();
    
    for(std::list<THBond>::const_iterator I = HBonds.begin(), end = HBonds.end(); I != end; ++I)
      dMolEnergy += I->dCalcEnergy();

    for(std::list<TSpring>::const_iterator I = Springs.begin(), end = Springs.end(); I != end; ++I)
      dRestrEnergy += I->dCalcEnergy();

    for(std::list<TBasePair>::const_iterator I = BasePairs.begin(), end = BasePairs.end(); I != end; ++I)
      dRestrEnergy += I->dCalcEnergy();

    _dMolEnergy_ = dMolEnergy;
    return dMolEnergy + dRestrEnergy;
  }
#endif /*SEQUENTIAL*/

//---------------------------------------------------------------------------------------------------
#ifndef SEQUENTIAL
  double TMolecule::dCalcNonbondEnergyPARALLEL(double& _dMolEnergy_)
  {
    double dMolEnergy = 0., dRestrEnergy = 0.;
    dElectrEnergyPARALLEL = 0.;
    std::vector<pthread_t> Threads;
    
    if(TConsts::ELECTR && !Electros.empty() && TConsts::NUMTHREADS > 1)
    {
      Threads.resize(TConsts::NUMTHREADS-1);
      bMakeIteratorsPARALLEL();
      for(unsigned int j = 0u; j < Threads.size(); ++j)
        pthread_create(&Threads[j], NULL, pvCalcIntervalEnergyPARALLEL, (void*) &IteratorsPARALLEL[j]);
    }
    else
    {
      for(std::list<TElectr>::const_iterator I = Electros.begin(), end = Electros.end(); I != end; ++I)
        dElectrEnergyPARALLEL += I->dCalcEnergy();
    }
    
    for(std::list<TLennard>::const_iterator I = Lennards.begin(), end = Lennards.end(); I != end; ++I)
      dMolEnergy += I->dCalcEnergy();
    
    for(std::list<THBond>::const_iterator I = HBonds.begin(), end = HBonds.end(); I != end; ++I)
      dMolEnergy += I->dCalcEnergy();

    for(std::list<TSpring>::const_iterator I = Springs.begin(), end = Springs.end(); I != end; ++I)
      dRestrEnergy += I->dCalcEnergy();

    for(std::list<TBasePair>::const_iterator I = BasePairs.begin(), end = BasePairs.end(); I != end; ++I)
      dRestrEnergy += I->dCalcEnergy();

    
    for(unsigned int j = 0u; j < Threads.size(); ++j)
      pthread_join(Threads[j], NULL);
    
     dMolEnergy += dElectrEnergyPARALLEL;
    _dMolEnergy_ = dMolEnergy;
    
    return dMolEnergy + dRestrEnergy;
  }

  //---------------------------------------------------------------------------------------------------
  bool TMolecule::bMakeIteratorsPARALLEL()
  {
    static bool bDone = false;	// To avoid multiple calls.
    if(bDone) return false;
    
    const unsigned int N = Electros.size();
    const unsigned int iPerThread = double(N)/double(TConsts::NUMTHREADS-1);
    
    unsigned int j = 1u;
    std::list<TElectr>::const_iterator itBegin = Electros.begin(), 
				      itEnd   = Electros.begin();
    do
    {
      ++j;
      ++itEnd;
      if(!(j%iPerThread))
      {
	IteratorsPARALLEL.push_back(std::make_pair(itBegin, itEnd));
	itBegin = itEnd;
      }
    } while(j < N && itEnd != Electros.end());
    
    return bDone = true;
  }

  //---------------------------------------------------------------------------------------------------
  void* TMolecule::pvCalcIntervalEnergyPARALLEL(void* arg)
  {// arg is an argument of type std::pair<std::list<TElectr>::const_iterator, std::list<TElectr>::const_iterator>*
  // that describes interval of Electros to be processed
    std::list<TElectr>::const_iterator 
      itBegin = ((std::pair<std::list<TElectr>::const_iterator, std::list<TElectr>::const_iterator>*)arg)->first, 
      itEnd   = ((std::pair<std::list<TElectr>::const_iterator, std::list<TElectr>::const_iterator>*)arg)->second;
    
    double dIntervalEnergy = 0.;
    for(std::list<TElectr>::const_iterator I = itBegin; I != itEnd; ++I)
      dIntervalEnergy += I->dCalcEnergy();
    
    pthread_mutex_lock(&mutexPARALLEL);
    dElectrEnergyPARALLEL += dIntervalEnergy;
    pthread_mutex_unlock(&mutexPARALLEL);
    
    return NULL;
  }
#endif /*SEQUENTIAL*/
//---------------------------------------------------------------------------------------------------
