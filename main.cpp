/***************************************************************************************************
* QRNA - Quick Refinement of Nucleic Acids by Juliusz Stasiewicz (jstasiewicz@genesilico.pl)
* main.cpp - this is the main program file.
***************************************************************************************************/

#include "main.h"

//---------------------------------------------------------------------------------------------------
int main(int iArgC, char* ppcArgV[])
{
  if(bParseArgs(iArgC, ppcArgV))
  {
    if(TConsts::INPUTPDB == "") std::cerr << "Error: no input PDB file specified." << std::endl;
    else 
    {
      TConsts::CheckAndDisplayParams();
      
      try{bReadPDB(TConsts::INPUTPDB);}
      catch(const std::bad_alloc& BA) 
      {
	std::cerr << "Fatal Error: Out of memory. Terminating." << std::endl;
	_exit(EXIT_FAILURE); 	// Achtung: nonstandard function used. Ordinary exit() can lead to exception-loops.
      }
      
      std::clog << "Number of molecules (chains) read: " << Molecules.size() << std::endl 
		<< Textutils::pcBar1;

      if(!Molecules.empty()) bMinimize();
      else std::cerr << "Warning: empty system, cannot proceed." << std::endl;
    }
  }
    
  return EXIT_SUCCESS;
}

//---------------------------------------------------------------------------------------------------
bool bParseArgs(int iArgC, char* ppcArgV[])
{// Returns true if it read anything that made sense.
  bool bRetVal = false;
  
  Textutils::bExtractPath(ppcArgV[0]);
  TConsts::bSetPaths();
  
  if(iArgC == 1) 
  {
    DisplayInfo();
    return false;
  }
  
  std::string sInput, sOutput, sDistRestraints;
  for(int i = 1; i < iArgC; ++i)
  {
    if(!std::strcmp(ppcArgV[i], "-c") || !std::strcmp(ppcArgV[i], "-C"))
    {
      if(i + 1 < iArgC)
      {
	TConsts::PARAMSFILE = ppcArgV[++i];
	bRetVal = true;
      }
      else std::cerr << "Error: expected config file name after '-c'." << std::endl;
    }    
    else if(!std::strcmp(ppcArgV[i], "-i") || !std::strcmp(ppcArgV[i], "-I"))
    {
      if(i + 1 < iArgC)
      {
	sInput = ppcArgV[++i];
	bRetVal = true;
      }
      else std::cerr << "Error: expected input PDB file name after '-i'." << std::endl;
    }
    else if(!std::strcmp(ppcArgV[i], "-o") || !std::strcmp(ppcArgV[i], "-O"))
    {
      if(i + 1 < iArgC)	
      {
	sOutput = ppcArgV[++i];
	bRetVal = true;
      }
      else std::cerr << "Warning: expected file name after '-o'." << std::endl;
    }
    else if(!std::strcmp(ppcArgV[i], "-m") || !std::strcmp(ppcArgV[i], "-M"))
    {
      if(i + 1 < iArgC)	
      {
	sDistRestraints = ppcArgV[++i];
	bRetVal = true;
      }
      else std::cerr << "Warning: expected file name after '-m'." << std::endl;
    }
    else if(!std::strcmp(ppcArgV[i], "-P"))
    {
      TConsts::POSRESTRAINTS = true;
      bRetVal = true;
    }
    else
    {
      std::cerr << "Warning: unknown option: " << ppcArgV[i] << std::endl;
      DisplayInfo();
      return false;
    }
  }
  
  if((!TConsts::bUserInit()) && (TConsts::PARAMSFILE != "")) 
    std::cerr << "Warning: cannot read from config file: " << TConsts::PARAMSFILE << std::endl;

  if(sDistRestraints != "")
  {
    if(TConsts::RESTRFILE != "" && sDistRestraints != TConsts::RESTRFILE)
      std::cerr << "Warning: conflicting restraints descriptions (in config file and parameters line). The latter is used." << std::endl;
    TConsts::RESTRFILE = sDistRestraints;
  }
  if((!TConsts::bRestrInit()) && (TConsts::RESTRFILE != "")) 
    std::cerr << "Warning: cannot read restraints from file: " << TConsts::RESTRFILE << std::endl;
  
  if(sInput != "")
  {
    if(TConsts::INPUTPDB != "" && TConsts::INPUTPDB != sInput)
      std::cerr << "Warning: input PDB file names (from config file and parameters line) disagree. The latter is used." << std::endl;
    TConsts::INPUTPDB = sInput;		// Override the filenames specified in PARAMSFILE
  }
  if(sOutput != "") TConsts::OUTPUTPDB = sOutput;	// Override the filenames specified in PARAMSFILE
  else if(TConsts::OUTPUTPDB == "") TConsts::OUTPUTPDB = Textutils::sMakeOutfileName(TConsts::INPUTPDB);
    
  return bRetVal;
}

//---------------------------------------------------------------------------------------------------
void DisplayInfo()
{
  std::clog << "QRNA - Quick Refinement of Nucleic Acids (version 0.3)\n" 
	    << "     by Juliusz Stasiewicz (jstasiewicz@genesilico.pl)\n\n"
	    << "To use type:\n"
	    << "  QRNA -i <input PDBfile> [-o <output PDBfile>] [-c <configfile>] [-P] [-m <restraintsfile>]" << std::endl
	    << "OR specify <input PDBfile>, <output PDBfile> and <restraintsfile> in <configfile> and type just:" << std::endl
	    << "  QRNA -c <configfile>" << std::endl;
}

//---------------------------------------------------------------------------------------------------
bool bReadPDB(const std::string& _sFileName)
{
  std::ifstream ifsIn(_sFileName.c_str());
  bool bRetVal = false;
  
  if(ifsIn != NULL)
  {
    TConsts::bDefaultDataInit();	

    std::clog << Textutils::pcBar2 << "Reading PDB file: " << _sFileName << std::endl << Textutils::pcBar2;
    Molecules.push_back(TMolecule());
/*************************************************************************************************
 * Molecules cannot be copied! It would invalidate all bonds, angles, dihedrals,
 * etc. inside. The tricky way around is to store them in a list, so that reallocation is avoided.
 *************************************************************************************************/
    while(Molecules.back().bRead(ifsIn) || (!TConsts::PEDANTICMOL && !Molecules.back().bEmpty()) || ifsIn.good())
    {
      if(!(Molecules.back().bEmpty()))
      {
        Molecules.push_back(TMolecule());
        bRetVal = true;
      }
    }
    Molecules.pop_back();
    
    TConsts::bDefaultDataCleanUp();    
    MakeAllExternNonbonds();
  }
  else std::cerr << Textutils::pcBar2 << "Error: cannot read: " << _sFileName << std::endl << Textutils::pcBar2;  
  
  return bRetVal;
}

//---------------------------------------------------------------------------------------------------
bool bWritePDB(const std::string& _sFileName)
{ // Return value is not supported here.
  bool bRetVal = true;
  std::clog << "Writing PDB file: " << _sFileName << " ...";
  std::ofstream ofsOut(_sFileName.c_str());
  for(std::list<TMolecule>::iterator J = Molecules.begin(); J != Molecules.end(); ++J)
    if(!(J->bWrite(ofsOut))) bRetVal = false;
    
  if(bRetVal) std::clog << "Done." << std::endl;
  else std::cerr << "FAILED!" << std::endl;
  
  return bRetVal;
}

//---------------------------------------------------------------------------------------------------
void MakeAllExternNonbonds()
{
  if(TConsts::TALKATIVE > 0) 
    std::clog << Textutils::pcBar1 << "Building interresidual nonbonded pairs & H-bonds (it may take a while)..." << std::endl;
  
  for(std::list<TMolecule>::iterator I = Molecules.begin(); I != Molecules.end(); ++I)
    for(std::list<TMolecule>::iterator J = I; J != Molecules.end(); ++J)
      bMakeNonbondsBetween(*I, *J);
  
  for(unsigned int k = 0; k < TConsts::DistanceRestraintsDescriptors.size(); ++k)
  {
    if(!TConsts::DistanceRestraintsDescriptors[k].bIsUsed)
      std::cerr << "Warning: bad restrain description: " 
                << TConsts::DistanceRestraintsDescriptors[k].sGetLine() << std::endl;
  }
  for(unsigned int k = 0; k < TConsts::BPRestraintsDescriptors.size(); ++k)
  {
    if(!TConsts::BPRestraintsDescriptors[k].bIsUsed) 
      std::cerr << "Warning: bad base pair description: " 
                << TConsts::BPRestraintsDescriptors[k].sGetLine() << std::endl;
  }

  if(TConsts::SSCONSTR)
  {
    if(TConsts::SSDETECT) TMolecule::bMakeBasePairs();
    else if(!TConsts::SECSTRUCT.empty() && !Molecules.empty()) 
      Molecules.front().bImposeVienna(TConsts::SECSTRUCT); // SS can be provided for only one chain.
  }

  unsigned int uPosRestr = 0u;
  for(std::list<TMolecule>::iterator I = Molecules.begin(); I != Molecules.end(); ++I)
    uPosRestr += I->uGetNumberOfPosRestrs();
  
  std::clog << "Number of electrostatic pairs:   " << TMolecule::uGetNumberOfElectros()  << std::endl  
	    << "Number of van der Waals pairs:   " << TMolecule::uGetNumberOfLennards()  << std::endl  
	    << "Number of H-bonds built:         " << TMolecule::uGetNumberOfHBonds()    << std::endl
	    << "Number of spring restraints:     " << TMolecule::uGetNumberOfSprings()   << std::endl  
	    << "Number of positional restraints: " << uPosRestr << std::endl;
	    
  std::clog << "Number of base pairs built:      " << TMolecule::uGetNumberOfBasePairs() << std::endl;  
}

//---------------------------------------------------------------------------------------------------
double dTotalEnergy(double& _dMolEnergy_)
{// Returns total energy of the system.
  double dMolNonbondEnergy, dMolCovalentEnergy = 0.;
  
  double dNonbondEnergy 
  #ifdef SEQUENTIAL
                          = TMolecule::dCalcNonbondEnergy(dMolNonbondEnergy);
  #else
                          = TMolecule::dCalcNonbondEnergyPARALLEL(dMolNonbondEnergy);
  #endif
  
  double dCovalentEnergy = 0.;
  
  for(std::list<TMolecule>::iterator J = Molecules.begin(); J != Molecules.end(); ++J)
  {
    double dSingleMolCovalentEnergy;
    dCovalentEnergy += J->dCalcCovalentEnergy(dSingleMolCovalentEnergy);
    dMolCovalentEnergy += dSingleMolCovalentEnergy;
  }
  
  _dMolEnergy_ = dMolNonbondEnergy + dMolCovalentEnergy ;
  return dCovalentEnergy + dNonbondEnergy;
}

//---------------------------------------------------------------------------------------------------
bool bMinimize()
{
  for(unsigned int i = 1; i <= TConsts::NSTEPS; ++i)
  {
    for(std::list<TMolecule>::iterator J = Molecules.begin(); J != Molecules.end(); ++J)
    {
      J->bClearGradients();
      J->CalcCovalentGradients();
      if(TConsts::USEBORN && TConsts::ELECTR) J->CalcBornRadii();
    }
    TMolecule::bCalcNonbondGradients();
    if(TConsts::BLOWGUARD) CalcAverageNorm();

    if(TConsts::TALKATIVE > 0) std::clog << "Performing minimization step: " << i << ". ";
    double dMolEnergy;
    double dTotalEnergy = dGoldenMinimStep(dMolEnergy);
    if(TConsts::TALKATIVE > 0) 
      std::clog << "Total energy = " << std::setprecision(4) << dTotalEnergy << " kcal/mol " 
	        << "(" << dMolEnergy << " without restraints)" << std::endl;
    
    if(!(i%TConsts::WRITEFREQ))
    {
      if(TConsts::TRAJECTORY) bWritePDB(Textutils::sMakeOutfileName(TConsts::OUTPUTPDB, i));
      else bWritePDB(TConsts::OUTPUTPDB);	//each TConsts::WRITEFREQ steps write some output
    }
  }
  
  return true;
}

//---------------------------------------------------------------------------------------------------
double dGoldenMinimStep(double& _dMolEnergy_)
{
  static const double dGolden = .5*(sqrt(5.)-1.);
  double dStep = TConsts::MAXSTEP,
	    xL = dStep*(1.-dGolden),
	    xR = dStep*dGolden;
  double dMolEnergy;

  double dEnergyStart = dTotalEnergy(dMolEnergy);
  
  bShiftAllMolecules(-xR);
  double dEnergyR = dTotalEnergy(dMolEnergy);
  
  bShiftAllMolecules(xR-xL);
  double dEnergyL = dTotalEnergy(dMolEnergy);
  
  int iStepNo = 0;
  while(fabs(dStep) > TConsts::MINSTEP)
  {
    if(TConsts::TALKATIVE > 1) std::clog << std::endl << "Golden section step #" << ++iStepNo;
    if(dEnergyL < dEnergyR || (dEnergyL > dEnergyStart && dEnergyR > dEnergyStart))
    {
      dStep = xR;
      xR = xL;
      dEnergyR = dEnergyL;
      xL = dStep*(1.-dGolden);
      bShiftAllMolecules(xR-xL);
      dEnergyL = dTotalEnergy(dMolEnergy);	
    }
    else
    {
      dStep = dStep - xL;
      xL = xR;
      dEnergyL = dEnergyR;
      xR = dStep*dGolden;
      bShiftAllMolecules(-xR);
      dEnergyR = dTotalEnergy(dMolEnergy);
      bShiftAllMolecules(xR-xL);
    }
  }
  bShiftAllMolecules(xL-dStep/2.);		// The middle of interval.

  if(TConsts::TALKATIVE > 1) std::clog << ". Golden section search reached desired accuracy." << std::endl;

  _dMolEnergy_ = dMolEnergy;
  return (dEnergyL+dEnergyR)/2.;	// Return approximate energy of the system.
}
//---------------------------------------------------------------------------------------------------
void CalcAverageNorm()
{// Average norm of vecEnergyGradient-s. Useful if too big steps must be avoided.
  double dTotalNorm = 0;
  unsigned int N = 0;
  
  for(std::list<TMolecule>::iterator J = Molecules.begin(); J != Molecules.end(); ++J)
  {
    dTotalNorm += J->dAllVectorsNorm();
    N += J->uGetNumberOfAtoms();
  }
    
  dAverageNorm = sqrt(dTotalNorm/double(N));	//The same norm, but now per single atom. (Average shift)
}
//---------------------------------------------------------------------------------------------------
bool bShiftAllMolecules(double _dFactor)
{
  if(TConsts::bIsNan(_dFactor))
  {
    std::cerr << "Internal Warning: NaNs encountered. Possible disruption of system." << std::endl;
    return false;
  }

  if(TConsts::BLOWGUARD) _dFactor /= dAverageNorm;	//Rescaling - to avoid too big steps; (avg) max 1. angstr per atom

  for(std::list<TMolecule>::iterator J = Molecules.begin(); J != Molecules.end(); ++J)
    J->Shift(_dFactor);
  
  return true;
}
