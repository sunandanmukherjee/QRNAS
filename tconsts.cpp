#include "textutils.h"
#include "tconsts.h"

#include <cstdlib>


namespace TConsts
{
  std::map<std::string, int> AtomTypeNumbers;	  //associative table that translates atom names (C1, N*, HW, etc) 
						  //of type std::string into internal numbers (int)
  std::vector<vdWParam> GeneralvanderWaalsData;
  std::vector<std::vector<BondParam> > GeneralBondsData;
  std::vector<std::vector<std::vector<AngleParam> > > GeneralAnglesData;
  std::vector<std::vector<std::vector<std::vector<ImproperParam> > > > GeneralImpropersData;
  std::vector<std::vector<std::vector<std::vector<std::list<DihedralParam> > > > > GeneralDihedralsData;

  std::vector<TDistRestrDescr> DistanceRestraintsDescriptors;
  std::vector<TBasePairRestrDescr> BPRestraintsDescriptors;
  
//---------------------------------------------------------------------------------------------------
  const long double    PI = 3.1415926535897932384626433832795028841968;
  const long double SQRT2 = 1.41421356237309504880168872420969807857;
  
  const double COULOMB = 332.063711;		// Coulomb's electrostatic const "k" [Ang*kcal/mol/e^2]
  const double kB      = 0.0019872065;		// Boltzmann const [kcal * K^-1 * mol^-1]
  const double PVUnit  = 0.00001457;		// Constant such that PV*[atm]*[Ang^3]=[kcal/mol]
  const double EpsilonRNA = 20;			// Dielectric const. for average RNA (after Z-J Tan & S-J Chen)
  const double EpsilonWater = 81.;		// Dielectric const. for water

  std::string ATOMSFILE, VDWFILE, BONDSFILE, ANGLEFILE, IMPROFILE, DIHEDFILE;
  std::string ATOMSDIR, VDWDIR, BONDSDIR, ANGLEDIR, IMPRODIR, DIHEDDIR, PDBDIR, PREPDIR;
  std::string INPUTPDB, OUTPUTPDB, PARAMSFILE, RESTRFILE, QRNAS_FF_DIR, PROGRAMPATH;
  
  
//---------------------------------------------------------------------------------------------------
/**************************************************************************
 * Default values for user-specifiable parameters follow.
 **************************************************************************/  

  unsigned int NUMTHREADS   
  #ifdef SEQUENTIAL
                            = 1u;
  #else
                            = 4u;
  #endif /*SEQUENTIAL*/

  unsigned int NSTEPS       = 100000u;		// Max number of MM steps.
  unsigned int TALKATIVE    = 1u;
  unsigned int WRITEFREQ    = 100u;
  unsigned int PEDANTICMOL  = 0u;
  
  bool   TRAJECTORY         = false;
  bool   USEBORN	    = false;
  bool   SSCONSTR	    = true;
  bool   SSDETECT	    = true;
  bool   HBONDS		    = true;
  bool   ELECTR		    = true;
  bool   VANDERWAALS	    = true;
  bool   BLOWGUARD          = true;

  std::string SECSTRUCT;

  bool   POSRESTRAINTS      = false;
  
  double CUTOFF             = 12.;		// Cutoff radius [Angstr\"om].
  double CUTOFF_SQR         = CUTOFF*CUTOFF;
  
  double MAXSTEP	    = 1.e-4;
  double MINSTEP	    = 1.e-8;
  
//---------------------------------------------------------------------------------------------------
  bool bSetPaths()
  {
    if (std::getenv("QRNAS_FF_DIR"))
    {
        const char* env_ff = std::getenv("QRNAS_FF_DIR");
        QRNAS_FF_DIR = std::string(env_ff);
    }
    if (QRNAS_FF_DIR != "")
    {
        ATOMSDIR  = QRNAS_FF_DIR + "/ATOMS/";
        VDWDIR    = QRNAS_FF_DIR + "/VDW/";
        BONDSDIR  = QRNAS_FF_DIR + "/BONDS/";
        ANGLEDIR  = QRNAS_FF_DIR + "/ANGLE/";
        IMPRODIR  = QRNAS_FF_DIR + "/IMPRO/";
        DIHEDDIR  = QRNAS_FF_DIR + "/DIHED/";
        PDBDIR    = QRNAS_FF_DIR + "/PDB/";
        PREPDIR   = QRNAS_FF_DIR + "/PREPI/";
    }
    else if (PROGRAMPATH != "")
    {
        ATOMSDIR  = PROGRAMPATH + "/forcefield/ATOMS/";
        VDWDIR    = PROGRAMPATH + "/forcefield/VDW/";
        BONDSDIR  = PROGRAMPATH + "/forcefield/BONDS/";
        ANGLEDIR  = PROGRAMPATH + "/forcefield/ANGLE/";
        IMPRODIR  = PROGRAMPATH + "/forcefield/IMPRO/";
        DIHEDDIR  = PROGRAMPATH + "/forcefield/DIHED/";
        PDBDIR    = PROGRAMPATH + "/forcefield/PDB/";
        PREPDIR   = PROGRAMPATH + "/forcefield/PREPI/";
    }

    ATOMSFILE = ATOMSDIR + "gen_atoms.dat";
    VDWFILE   = VDWDIR   + "gen_vdw.dat";
    BONDSFILE = BONDSDIR + "gen_bonds.dat";
    ANGLEFILE = ANGLEDIR + "gen_angles.dat";
    IMPROFILE = IMPRODIR + "gen_impro.dat";
    DIHEDFILE = DIHEDDIR + "gen_dihed.dat";
     
    if (!std::ifstream(ATOMSFILE.c_str()).good())
    {
        std::cerr << std::endl<<"The path to forcefield is directory is not configured"<<std::endl;
        std::cerr << "It can be configured either by setting QRNAS_FF_DIR in bashrc"<<std::endl;
        std::cerr << "Or by copying or setting symbolic link to forcefield in current path"<<std::endl<<std::endl;
        _exit(EXIT_FAILURE);
    }
    
    return true;
  }

//---------------------------------------------------------------------------------------------------
  bool bUserInit(const std::string& _sFileName)	
  {// Read user-specified params.
    bool bRetVal = false;

    std::ifstream ifsIn(_sFileName.c_str());
    if(ifsIn != NULL)
    {
      if(TALKATIVE > 0)
        std::clog << Textutils::pcBar2 << "Reading config file: " 
                  << _sFileName << std::endl << Textutils::pcBar2;
      
      std::string sLine;
      while(std::getline(ifsIn, sLine))
      {
	if(sLine.length() >= 9 && Textutils::sToUpper(sLine.substr(0, 9)) != "OUTPUTPDB" 
			       && Textutils::sToUpper(sLine.substr(0, 8)) != "INPUTPDB" 
			       && Textutils::sToUpper(sLine.substr(0, 9)) != "RESTRFILE")
	{
	  // File names cannot be modified as below:
          int iComma = sLine.find(',');
	  if(iComma != (int)std::string::npos) sLine[iComma] = '.';	// Replace ',' with '.' if some dumbass uses commas as decimal points
	  int iEqual = sLine.find_first_of('=');
	  if(iEqual != (int)std::string::npos) sLine[iEqual] = ' ';	// Replace '=' with ' ' to allow for syntax like a=7.6
	}
	int iComment = sLine.find_first_of("%#");
	if(iComment != (int)std::string::npos)
	  sLine = sLine.substr(0, iComment);
	
	std::istringstream issIn(sLine);
	std::string sFlag;
	issIn >> sFlag;
	Textutils::ToUpper(sFlag);
	
	if(sFlag == "NSTEPS")
	{
	  int iNSteps;
	  if(issIn >> iNSteps && iNSteps > 0)
	    NSTEPS = iNSteps;
	  else
	    std::cerr << "Error: invalid value " << iNSteps << " of " << sFlag << ". Using default " << NSTEPS << std::endl;
	}
	else if(sFlag == "VERBOSE") 
	{
	  int iTalkative;
	  if(issIn >> iTalkative && iTalkative >= 0)
	    TALKATIVE = iTalkative;
	  else
	    std::cerr << "Error: invalid value " << iTalkative << " of " << sFlag << ". Using default " << TALKATIVE << std::endl;
	}
	else if(sFlag == "WRITEFREQ") 
	{
	  int iWritefreq;
	  if(issIn >> iWritefreq && iWritefreq > 0)
	    WRITEFREQ = iWritefreq;
	  else
	    std::cerr << "Error: invalid value " << iWritefreq << " of " << sFlag << ". Using default " << WRITEFREQ << std::endl;
	}
	else if(sFlag == "TRAJECTORY") 
	{
	  int iTrajectory;
	  if(issIn >> iTrajectory && iTrajectory >= 0)
	    TRAJECTORY = iTrajectory;
	  else
	    std::cerr << "Error: invalid value " << iTrajectory << " of " << sFlag << ". Using default " << TRAJECTORY << std::endl;
	}
	else if(sFlag == "INPUTPDB") 
	{
	  std::string sInputpdb;
	  if(issIn >> sInputpdb && sInputpdb != "")
	    INPUTPDB = sInputpdb;
	  else
	    std::cerr << "Error: invalid value " << sInputpdb << " of " << sFlag << ". Using " << INPUTPDB << std::endl;
	}
	else if(sFlag == "OUTPUTPDB") 
	{
	  std::string sOutputpdb;
	  if(issIn >> sOutputpdb && sOutputpdb != "")
	    OUTPUTPDB = sOutputpdb;
	  else
	    std::cerr << "Error: invalid value " << sOutputpdb << " of " << sFlag << ". Using default " << OUTPUTPDB << std::endl;
	}
	else if(sFlag == "RESTRFILE") 
	{
	  std::string sRestrFile;
	  if(issIn >> sRestrFile && sRestrFile != "")
	    RESTRFILE = sRestrFile;
	  else
	    std::cerr << "Error: invalid value " << sRestrFile << " of " << sFlag << ". Assuming no pairwise restraints." << std::endl;
	}
	else if(sFlag == "CUTOFF") 
	{
	  double dCutoff;
	  if(issIn >> dCutoff && dCutoff >= 0.)
	  {
	    CUTOFF = dCutoff;
	    CUTOFF_SQR = CUTOFF*CUTOFF;	    
	  }
	  else
	    std::cerr << "Error: invalid value " << dCutoff << " of " << sFlag << ". Using default " << CUTOFF << std::endl;
	}
	else if(sFlag == "MAXSTEP") 
	{
	  double dMaxstep;
	  if(issIn >> dMaxstep && dMaxstep >= 0.)
	    MAXSTEP = dMaxstep;
	  else
	    std::cerr << "Error: invalid value " << dMaxstep << " of " << sFlag << ". Using default " << MAXSTEP << std::endl;
	}
	else if(sFlag == "MINSTEP") 
	{
	  double dMinstep;
	  if(issIn >> dMinstep && dMinstep >= 0. && dMinstep < MAXSTEP)
	    MINSTEP = dMinstep;
	  else
	    std::cerr << "Error: invalid value " << dMinstep << " of " << sFlag << ". Using default " << MINSTEP << std::endl;
	}
	else if(sFlag == "PEDANTICMOL") 
	{
	  int iPedanticmol;
	  if(issIn >> iPedanticmol && iPedanticmol >= 0)
	    PEDANTICMOL = iPedanticmol;
	  else
	    std::cerr << "Error: invalid value " << iPedanticmol << " of " << sFlag << ". Using default " << PEDANTICMOL << std::endl;
	}
	else if(sFlag == "USEBORN") 
	{
	  int iUseBorn;
	  if(issIn >> iUseBorn && iUseBorn >= 0)
	    USEBORN = iUseBorn;
	  else
	    std::cerr << "Error: invalid value " << iUseBorn << " of " << sFlag << ". Using default " << USEBORN << std::endl;
	}
	else if(sFlag == "SSCONSTR") 
	{
	  int iSSConstr;
	  if(issIn >> iSSConstr && iSSConstr >= 0)
	    SSCONSTR = iSSConstr;
	  else
	    std::cerr << "Error: invalid value " << iSSConstr << " of " << sFlag << ". Using default " << SSCONSTR << std::endl;
	}
	else if(sFlag == "SSDETECT") 
	{
	  int iSSDetect;
	  if(issIn >> iSSDetect && iSSDetect >= 0)
	    SSDETECT = iSSDetect;
	  else
	    std::cerr << "Error: invalid value " << iSSDetect << " of " << sFlag << ". Using default " << SSDETECT << std::endl;
	}
	else if(sFlag == "SECSTRUCT") 
	{
	  std::string sSecStruct;
	  if(issIn >> sSecStruct && Textutils::bValidateVienna(sSecStruct))
	    SECSTRUCT = sSecStruct;
	  else
	    std::cerr << "Error: invalid value " << sSecStruct << " of " << sFlag << ". None assumed." << std::endl;
	}
	else if(sFlag == "HBONDS") 
	{
	  int iHBonds;
	  if(issIn >> iHBonds && iHBonds >= 0)
	    HBONDS = iHBonds;
	  else
	    std::cerr << "Error: invalid value " << iHBonds << " of " << sFlag << ". Using default " << HBONDS << std::endl;
	}
	else if(sFlag == "ELECTR") 
	{
	  int iElectr;
	  if(issIn >> iElectr && iElectr >= 0)
	    ELECTR = iElectr;
	  else
	    std::cerr << "Error: invalid value " << iElectr << " of " << sFlag << ". Using default " << ELECTR << std::endl;
	}
	else if(sFlag == "VANDERWAALS" || sFlag == "VDW") 
	{
	  int iVdw;
	  if(issIn >> iVdw && iVdw >= 0)
	    VANDERWAALS = iVdw;
	  else
	    std::cerr << "Error: invalid value " << iVdw << " of " << sFlag << ". Using default " << VANDERWAALS << std::endl;
	}
	else if(sFlag == "NUMTHREADS") 
	{
	  int iNumThreads;
	  if(issIn >> iNumThreads && iNumThreads > 0) 
	    NUMTHREADS = iNumThreads;
	  else
	    std::cerr << "Error: invalid value " << iNumThreads << " of " << sFlag << ". Using default " << NUMTHREADS << std::endl;
	}
	else if(sFlag == "BLOWGUARD") 
	{
	  int iBlowGuard;
	  if(issIn >> iBlowGuard && iBlowGuard >= 0) 
	    BLOWGUARD = iBlowGuard;
	  else
	    std::cerr << "Error: invalid value " << iBlowGuard << " of " << sFlag << ". Using default " << BLOWGUARD << std::endl;
	}
	else if(sFlag == "POSRESTRAINTS") 
	{
	  int iPosRes;
	  if(issIn >> iPosRes && iPosRes >= 0) 
	    POSRESTRAINTS = iPosRes;
	  else
	    std::cerr << "Error: invalid value " << iPosRes << " of " << sFlag << ". Using default " << POSRESTRAINTS << std::endl;
	}
	else if(sFlag != "")
	{
	  std::cerr << "Error: unknown flag: " << sFlag << " in parameters file." << std::endl;
	}
      }
      bRetVal = true;
    }
    return bRetVal;
  }

//---------------------------------------------------------------------------------------------------
  bool bRestrInit(const std::string& _sFileName)
  {
    bool bRetVal = false;
    std::ifstream ifsIn(_sFileName.c_str());
    if(ifsIn != NULL)
    {
      if(TALKATIVE > 0)
        std::clog << Textutils::pcBar2 << "Reading pairwise restraints from file: " 
                  << _sFileName << std::endl << Textutils::pcBar2;
      
      std::string sLine;
      while(std::getline(ifsIn, sLine))
      {
	Textutils::bCropSpaces(sLine);
	Textutils::ToUpper(sLine);
	bool bOKLine = false;
	
	int iComment = sLine.find_first_of("%#");
	if(iComment != (int)std::string::npos)
	  sLine = sLine.substr(0, iComment);
	if(sLine.length() > 8)
	{
	  if(sLine.substr(0,8) == "DISTANCE")
	  {
	    TDistRestrDescr drdDistRestrDescrTemp(sLine);
	    if(drdDistRestrDescrTemp.bIsOK())
	    {
	      DistanceRestraintsDescriptors.push_back(drdDistRestrDescrTemp);
	      if(TALKATIVE > 1) std::clog << "Pairwise atom restraint read: " << sLine << std::endl;
	      bOKLine = true;
	    }
	  }
	  if(sLine.substr(0,8) == "BASEPAIR")
	  {
	    TBasePairRestrDescr bprdRestrDescrTemp(sLine);
	    if(bprdRestrDescrTemp.bIsOK())
	    {
	      BPRestraintsDescriptors.push_back(bprdRestrDescrTemp);
	      if(TALKATIVE > 1) std::clog << "Basepair restraint read: " << sLine << std::endl;
	      bOKLine = true;
	    }
	  }
	}
	
	if(!bOKLine && !sLine.empty()) std::cerr << "Warning: misformatted restraint data: " << sLine << std::endl;
      }
      bRetVal = true;
    }
    return bRetVal;
  }
//---------------------------------------------------------------------------------------------------

  void CheckAndDisplayParams()
  {
    if(TALKATIVE > 0) std::clog << Textutils::pcBar1;

    if(WRITEFREQ == 0)
    {
      WRITEFREQ = NSTEPS;
      std::cerr << "Warning: WRITEFREQ cannot be equal 0. Adjusted to WRITEFREQ = NSTEPS." << std::endl;
    }
    if(WRITEFREQ > NSTEPS)
    {
      WRITEFREQ = NSTEPS;
      std::cerr << "Warning: WRITEFREQ is greater than NSTEPS. Adjusted to WRITEFREQ = NSTEPS." << std::endl;
    }
    if(MAXSTEP <= MINSTEP)
    {
      MINSTEP = MAXSTEP/10.;
      std::cerr << "Warning: MAXSTEP must be greater than MINSTEP. Adjusted to MINSTEP = MAXSTEP/10." << std::endl;
    }
    if(USEBORN && !ELECTR)
    {
      USEBORN = false;
      std::cerr << "Warning: cannot use Born solvent when electrostatics is off. Adjusted to USEBORN = 0." << std::endl;
    }
    if(SSDETECT && !SSCONSTR)
    {
      SSDETECT = false;
      std::cerr << "Warning: Pointless to detect sec. structure when SSCONSTR is off. Adjusted to SSDETECT = 0." << std::endl;
    }
    if(!SECSTRUCT.empty() && !SSCONSTR)
    {
      SSCONSTR = true;
      std::cerr << "Warning: Secondary struct. specified, but SSCONSTR is off. Adjusted to SSCONSTR = 1." << std::endl;
    }
    if(!BPRestraintsDescriptors.empty() && !SSCONSTR)
    {
      SSCONSTR = true;
      std::cerr << "Warning: Basepair restraints specified, but SSCONSTR is off. Adjusted to SSCONSTR = 1." << std::endl;
    }
    if(SSDETECT && !SECSTRUCT.empty())
    {
      SSDETECT = false;
      std::cerr << "Warning: Specification of SECSTRUCT overrides SSDETECT. Adjusted to SSDETECT = 0." << std::endl;
    }
    if(SSDETECT && !BPRestraintsDescriptors.empty())
    {
      SSDETECT = false;
      std::cerr << "Warning: Specification of basepair restraints overrides SSDETECT. Adjusted to SSDETECT = 0." << std::endl;
    }
    if(!SECSTRUCT.empty() && !BPRestraintsDescriptors.empty())
    {
      SECSTRUCT.clear();
      std::cerr << "Warning: Base pair restraints override SECSTRUCT." << std::endl;
    }
    if(SSCONSTR && !HBONDS && SSDETECT)
    {
      SSCONSTR = false;
      std::cerr << "Warning: Cannot detect secondary structure with HBONDS = 0. Adjusted to SSCONSTR = 0." << std::endl;
    }
    if(NUMTHREADS > 1 && !ELECTR)
    {
      NUMTHREADS = 1u;
      std::cerr << "Warning: Currently only electrostatics can be parallelized (and you set ELECTR = 0). Adjusted to NUMTHREADS = 1." << std::endl;
    }
    #ifdef SEQUENTIAL
    if(NUMTHREADS > 1)
    {
      NUMTHREADS = 1u;
      std::cerr << "Warning: Non-parallel version of QRNA cannot spawn threads." << std::endl;
    }
    #endif /*SEQUENTIAL*/
    if(NUMTHREADS >= uGetNumCores())
    {
      std::cerr << "Warning: More threads than cores available." << std::endl;
    }
    
    if(TALKATIVE > 0)
    {
      std::clog << "Parameters to be used:" << std::endl
	        << "Number of minimization steps (NSTEPS): " << NSTEPS << std::endl
	        << "Golden section step sizes: maximal (MAXSTEP): " << MAXSTEP << std::endl
	        << "                           minimal (MINSTEP): " << MINSTEP << "" << std::endl
	        << "Verbosity (VERBOSE): " << TALKATIVE << std::endl
	        << "Write output every: " << WRITEFREQ << " steps (WRITEFREQ)" << std::endl
	        << "Record trajectory (TRAJECTORY): " << (TRAJECTORY ? "enabled" : "disabled") << std::endl
	        << "Pedantic treatment of molecules (PEDANTICMOL): " << (PEDANTICMOL ? "enabled" : "disabled") << std::endl
  	        << "Van der Waals cutoff radius (CUTOFF): " << CUTOFF << std::endl
	        << "Born generalized solvent (USEBORN): " << (USEBORN ? "enabled" : "disabled") << std::endl
	        << "Constraints on secondary struct. (SSCONSTR): " << (SSCONSTR ? "enabled" : "disabled") << std::endl
	        << "Auto-detection of secondary struct. (SSDETECT): " << (SSDETECT ? "enabled" : "disabled") << std::endl
	        << "Secondary structure (SECSTRUCT): " << (SECSTRUCT != "" ? SECSTRUCT : "none specified" ) << std::endl
	        << "Hydrogen bonds (HBONDS): " << (HBONDS ? "enabled" : " disabled") << std::endl
	        << "Electrostatics (ELECTR): " << (ELECTR ? "enabled" : " disabled") << std::endl
	        << "van der Waals interactions (VANDERWAALS): " << (VANDERWAALS ? "enabled" : "disabled") << std::endl
	        << "Anti-blow-up protection (BLOWGUARD): " << (BLOWGUARD ? "enabled" : "disabled") << std::endl
	        << "Occupancy as positional restraint (POSRESTRAINTS): " << (POSRESTRAINTS ? "enabled" : "disabled") << std::endl
	        << "Restraints description (RESTRFILE): " << (RESTRFILE == "" ? "none" : RESTRFILE) << std::endl
	        << "Number of parallel threads: " << NUMTHREADS << std::endl;
		      
      std::clog << std::endl << Textutils::pcBar1;
    }
  }
  
//---------------------------------------------------------------------------------------------------
  bool bDefaultDataCleanUp()			//release some memory
  {
    AtomTypeNumbers.clear();
    GeneralvanderWaalsData.clear();
    GeneralBondsData.clear();
    GeneralAnglesData.clear();
    GeneralImpropersData.clear();
    GeneralDihedralsData.clear();
    
    return true;
  }

//---------------------------------------------------------------------------------------------------
  int iAtomsInit(const std::string& _sFileName)
  {
    int iRetVal = 0;
    std::string sLine;
    
    std::ifstream ifsIn(_sFileName.c_str());
    if(ifsIn != NULL)
    {
      while(std::getline(ifsIn, sLine))
      {
	if(sLine == "END") break;
	else if(sLine.length() >= 2)
	{
	  std::string sAtomCode = sLine.substr(0,2);		//numbers found by eyeballing parm99 file, should be checked!
	  Textutils::bCropSpaces(sAtomCode);
		  
	  AtomTypeNumbers[sAtomCode] = iRetVal++;
	}
      }
      ifsIn.close();
    }
    else std::cerr << "Internal Error: cannot read file: " << _sFileName << std::endl;
    
    return iRetVal;
  }

//---------------------------------------------------------------------------------------------------
  int ivdWInit(const std::string& _sFileName, std::vector<vdWParam>& _vanderWaalsData_)
  {
    int iRetVal = 0;
    std::string sLine;
    
    std::ifstream ifsIn(_sFileName.c_str());
    if(ifsIn != NULL)
    {
      _vanderWaalsData_.resize(AtomTypeNumbers.size());		//as many atoms in vdWData as many atom types
      while(std::getline(ifsIn, sLine))
      {
	if(sLine == "END") break;
	else if(sLine.length() >= 30)
	{
	  std::string sAtomCode = sLine.substr(2,2);
	  Textutils::bCropSpaces(sAtomCode);
	  
	  if(AtomTypeNumbers.find(sAtomCode) != AtomTypeNumbers.end())
	  {
	    double _dR0 = Textutils::dToDouble(sLine.substr(10,10)),
		    _dEpsilon = Textutils::dToDouble(sLine.substr(20,10));
		  
	    _vanderWaalsData_[AtomTypeNumbers[sAtomCode]] = vdWParam(_dR0, _dEpsilon);
	    ++iRetVal;
	  }
	  else if(TALKATIVE > 0) 
	    std::cerr << "Internal Warning: unknown atom type : " << sAtomCode << std::endl;
	}
      }
    }
    else    std::cerr << "Internal Error: cannot read file: " << _sFileName << std::endl;
    
    return iRetVal;
  }
  
//---------------------------------------------------------------------------------------------------
  int iBondsInit(const std::string& _sFileName, std::vector<std::vector<BondParam> >& _BondsData_)
  {
    int iRetVal = 0;
    std::string sLine;
    
    std::ifstream ifsIn(_sFileName.c_str());
    if(ifsIn != NULL)
    {
      _BondsData_.resize(AtomTypeNumbers.size(), std::vector<BondParam>(AtomTypeNumbers.size()));	//resize BondsData to a square matrix of size N(AtomTypes)^2
      
      while(std::getline(ifsIn, sLine))
      {
	if(sLine == "END") break;
	else if(sLine.length() >= 25)
	{
	  std::string sAtomCode1 = sLine.substr(0,2), sAtomCode2 = sLine.substr(3,2);
	  Textutils::bCropSpaces(sAtomCode1);
	  Textutils::bCropSpaces(sAtomCode2);
	  
	  if(AtomTypeNumbers.find(sAtomCode1) != AtomTypeNumbers.end() && 
	     AtomTypeNumbers.find(sAtomCode2) != AtomTypeNumbers.end() )
	  {
	    double _dKr = Textutils::dToDouble(sLine.substr(5,10)),
		    _dREq = Textutils::dToDouble(sLine.substr(15,10));
		  
	    _BondsData_[AtomTypeNumbers[sAtomCode1]][AtomTypeNumbers[sAtomCode2]] = BondParam(_dKr, _dREq);
	    _BondsData_[AtomTypeNumbers[sAtomCode2]][AtomTypeNumbers[sAtomCode1]] = BondParam(_dKr, _dREq);
	    ++iRetVal;
	  }
	  else if(TALKATIVE > 0) 
	    std::cerr << "Internal Warning: unknown atom type(s) in bond: " << sAtomCode1 << " " << sAtomCode2 << std::endl;
	}
      }
    }
    else    std::cerr << "Internal Error: cannot read file: " << _sFileName << std::endl;    
    
    return iRetVal;
  }
  
//---------------------------------------------------------------------------------------------------
  int iAnglesInit(const std::string& _sFileName, std::vector<std::vector<std::vector<AngleParam> > >& _AnglesData_)
  {
    int iRetVal = 0;
    std::string sLine;
    
    std::ifstream ifsIn(_sFileName.c_str());
    if(ifsIn != NULL)
    {
      _AnglesData_.resize(AtomTypeNumbers.size());
      for(unsigned int i = 0; i < _AnglesData_.size(); ++i)
	_AnglesData_[i].resize(AtomTypeNumbers.size(), std::vector<AngleParam>(AtomTypeNumbers.size()));	//resize AnglesData to a cubic matrix of size N(AtomTypes)^3
      
      while(std::getline(ifsIn, sLine))
      {
	if(sLine == "END") break;
	else if(sLine.length() >= 28)
	{
	  std::string sAtomCode1 = sLine.substr(0,2), 
		      sAtomCode2 = sLine.substr(3,2),
		      sAtomCode3 = sLine.substr(6,2);
	  Textutils::bCropSpaces(sAtomCode1);
	  Textutils::bCropSpaces(sAtomCode2);
	  Textutils::bCropSpaces(sAtomCode3);
	  
	  if(AtomTypeNumbers.find(sAtomCode1) != AtomTypeNumbers.end() && 
	      AtomTypeNumbers.find(sAtomCode2) != AtomTypeNumbers.end() && 
	      AtomTypeNumbers.find(sAtomCode3) != AtomTypeNumbers.end())
	  {
	    double _dKTheta = Textutils::dToDouble(sLine.substr(8,10)),
		    _dThetaEq = Textutils::dToDouble(sLine.substr(18,10));
		  
	    _AnglesData_[AtomTypeNumbers[sAtomCode1]][AtomTypeNumbers[sAtomCode2]][AtomTypeNumbers[sAtomCode3]] = AngleParam(_dKTheta, _dThetaEq);
	    _AnglesData_[AtomTypeNumbers[sAtomCode3]][AtomTypeNumbers[sAtomCode2]][AtomTypeNumbers[sAtomCode1]] = AngleParam(_dKTheta, _dThetaEq);

	    ++iRetVal;
	  }
	  else if(TALKATIVE > 0) 
	    std::cerr << "Internal Warning: unknown atoms in angle: " << sAtomCode1 << " " << sAtomCode2 
			  << " " << sAtomCode3 << std::endl;
	}
      }
    }
    else    std::cerr << "Internal Error: cannot read file: " << _sFileName << std::endl;   

    return iRetVal;      
  }

//---------------------------------------------------------------------------------------------------
  int iImpropersInit(const std::string& _sFileName, std::vector<std::vector<std::vector<std::vector<ImproperParam> > > >& _ImpropersData_)
  {
    int iRetVal = 0;
    std::string sLine;
    
    std::ifstream ifsIn(_sFileName.c_str());
    if(ifsIn != NULL)
    {
      _ImpropersData_.resize(AtomTypeNumbers.size());
      for(unsigned int i = 0; i < _ImpropersData_.size(); ++i)
	_ImpropersData_[i].resize(AtomTypeNumbers.size());	// Resize ImpropersData to a hypercubic matrix of size N(AtomTypes)^4
      for(unsigned int i = 0; i < _ImpropersData_.size(); ++i)
	for(unsigned int j = 0; j < _ImpropersData_[i].size(); ++j)
	  _ImpropersData_[i][j].resize(AtomTypeNumbers.size(), std::vector<ImproperParam>(AtomTypeNumbers.size()));
      
      while(std::getline(ifsIn, sLine))
      {
	if(sLine == "END") break;
	else if(sLine.length() >= 60)
	{
	  std::string sAtomCode1 = sLine.substr(0,2), 
		      sAtomCode2 = sLine.substr(3,2),
		      sAtomCode3 = sLine.substr(6,2),
		      sAtomCode4 = sLine.substr(9,2);
	  Textutils::bCropSpaces(sAtomCode1);
	  Textutils::bCropSpaces(sAtomCode2);
	  Textutils::bCropSpaces(sAtomCode3);
	  Textutils::bCropSpaces(sAtomCode4);
	  
	  if((AtomTypeNumbers.find(sAtomCode1) != AtomTypeNumbers.end() || sAtomCode1 == "X") && 
	      (AtomTypeNumbers.find(sAtomCode2) != AtomTypeNumbers.end() || sAtomCode2 == "X") && 
	      (AtomTypeNumbers.find(sAtomCode3) != AtomTypeNumbers.end() || sAtomCode3 == "X") &&
	      (AtomTypeNumbers.find(sAtomCode4) != AtomTypeNumbers.end() || sAtomCode4 == "X"))
	  {
	    double _dKI    = Textutils::dToDouble(sLine.substr(15,15)),
		    _dPhase = Textutils::dToDouble(sLine.substr(30,15)),
		    _dPN    = Textutils::dToDouble(sLine.substr(45,15));
	      
	    unsigned int  ibeg = 0, iend = _ImpropersData_.size(),
			  jbeg = 0, jend = _ImpropersData_[0].size(),
			  kbeg = 0, kend = _ImpropersData_[0][0].size(),
			  lbeg = 0, lend = _ImpropersData_[0][0][0].size();
	    
	    if(sAtomCode1 != "X")
	    {
	      ibeg = AtomTypeNumbers[sAtomCode1];
	      iend = ibeg + 1;
	    }
	    if(sAtomCode2 != "X")
	    {
	      jbeg = AtomTypeNumbers[sAtomCode2];
	      jend = jbeg + 1;
	    }
	    if(sAtomCode3 != "X")
	    {
	      kbeg = AtomTypeNumbers[sAtomCode3];
	      kend = kbeg + 1;
	    }
	    if(sAtomCode4 != "X")
	    {
	      lbeg = AtomTypeNumbers[sAtomCode4];
	      lend = lbeg + 1;
	    }
	    for(unsigned int i = ibeg; i < iend; ++i)
	      for(unsigned int j = jbeg; j < jend; ++j)
		for(unsigned int k = kbeg; k < kend; ++k)
		  for(unsigned int l = lbeg; l < lend; ++l)
		  {	      
		    _ImpropersData_[i][j][k][l] = ImproperParam(_dKI, _dPhase, _dPN);
		    _ImpropersData_[l][k][j][i] = ImproperParam(_dKI, _dPhase, _dPN);
		  }
	    ++iRetVal;
	  }
	  else if(TALKATIVE > 0) 
	    std::cerr << "Internal Warning: unknown atom(s) in improper: " << sAtomCode1 << " " << sAtomCode2 
			  << " " << sAtomCode3 << " " << sAtomCode4 << std::endl;
	}
      }
    }
    else    std::cerr << "Internal Error: cannot read file: " << _sFileName << std::endl; 
    
    return iRetVal;            
  }
  
//---------------------------------------------------------------------------------------------------
  int iDihedralsInit(const std::string& _sFileName, std::vector<std::vector<std::vector<std::vector<std::list<DihedralParam> > > > >& _DihedralsData_)
  {// Genius-made piece of code; don't fuck with it unless necessary!
    int iRetVal = 0;
    std::string sLine;
    
    std::ifstream ifsIn(_sFileName.c_str());
    if(ifsIn != NULL)
    {
      _DihedralsData_.resize(AtomTypeNumbers.size());
      
      for(unsigned int i = 0; i < _DihedralsData_.size(); ++i)
	_DihedralsData_[i].resize(AtomTypeNumbers.size());	// Resize DihedralsData to a hypercubic matrix of size N(AtomTypes)^4
      for(unsigned int i = 0; i < _DihedralsData_.size(); ++i)
	for(unsigned int j = 0; j < _DihedralsData_[i].size(); ++j)
	  _DihedralsData_[i][j].resize(AtomTypeNumbers.size(), std::vector<std::list<DihedralParam> >(AtomTypeNumbers.size()));
      
      while(std::getline(ifsIn, sLine))
      {
	if(sLine == "END") break;
	else if(sLine.length() >= 60)
	{
	  std::string sAtomCode1 = sLine.substr(0,2), 
		      sAtomCode2 = sLine.substr(3,2),
		      sAtomCode3 = sLine.substr(6,2),
		      sAtomCode4 = sLine.substr(9,2);
	  Textutils::bCropSpaces(sAtomCode1);
	  Textutils::bCropSpaces(sAtomCode2);
	  Textutils::bCropSpaces(sAtomCode3);
	  Textutils::bCropSpaces(sAtomCode4);
	  
	  if( (AtomTypeNumbers.find(sAtomCode1) != AtomTypeNumbers.end() || sAtomCode1 == "X") && 
	      (AtomTypeNumbers.find(sAtomCode2) != AtomTypeNumbers.end() || sAtomCode2 == "X") && 
	      (AtomTypeNumbers.find(sAtomCode3) != AtomTypeNumbers.end() || sAtomCode3 == "X") &&
	      (AtomTypeNumbers.find(sAtomCode4) != AtomTypeNumbers.end() || sAtomCode4 == "X"))
	  {
	    int    _iDiv   = Textutils::iToInt(sLine.substr(11,4));
	    double _dKI    = Textutils::dToDouble(sLine.substr(15,15)),
		    _dPhase = Textutils::dToDouble(sLine.substr(30,15)),
		    _dPN    = Textutils::dToDouble(sLine.substr(45,15));
	      
	    unsigned int  ibeg = 0, iend = _DihedralsData_.size(),
			  jbeg = 0, jend = _DihedralsData_[0].size(),
			  kbeg = 0, kend = _DihedralsData_[0][0].size(),
			  lbeg = 0, lend = _DihedralsData_[0][0][0].size();
	    
	    if(sAtomCode1 != "X")
	    {
	      ibeg = AtomTypeNumbers[sAtomCode1];
	      iend = ibeg + 1;
	    }
	    if(sAtomCode2 != "X")
	    {
	      jbeg = AtomTypeNumbers[sAtomCode2];
	      jend = jbeg + 1;
	    }
	    if(sAtomCode3 != "X")
	    {
	      kbeg = AtomTypeNumbers[sAtomCode3];
	      kend = kbeg + 1;
	    }
	    if(sAtomCode4 != "X")
	    {
	      lbeg = AtomTypeNumbers[sAtomCode4];
	      lend = lbeg + 1;
	    }
	    for(unsigned int i = ibeg; i < iend; ++i)
	      for(unsigned int j = jbeg; j < jend; ++j)
		for(unsigned int k = kbeg; k < kend; ++k)
		  for(unsigned int l = lbeg; l < lend; ++l)
		  {	      
		    if(!_DihedralsData_[i][j][k][l].empty() && _DihedralsData_[i][j][k][l].back().dPN >= 0.)
		    {						// The case is that some general dihedral description (e.g X-X-C-C)
		      _DihedralsData_[i][j][k][l].clear();	// occured before a better, more specific one (e.g C-N*-CT-O)
		      _DihedralsData_[l][k][j][i].clear();	// so must be discarded.
		    }
		    _DihedralsData_[i][j][k][l].push_back(DihedralParam(_iDiv, _dKI, _dPhase, _dPN));
		    if(l != i || k != j) _DihedralsData_[l][k][j][i].push_back(DihedralParam(_iDiv, _dKI, _dPhase, _dPN));		    
		  }
	    ++iRetVal;
	  }
	  else if(TALKATIVE > 0) 
	    std::cerr << "Internal Warning: unknown atom(s) in dihedral: " << sAtomCode1 << " " << sAtomCode2 
			  << " " << sAtomCode3 << " " << sAtomCode4 << std::endl;
	}
      }
    }
    else    std::cerr << "Internal Error: cannot read file: " << _sFileName << std::endl;    
    
    return iRetVal;            
  }

//---------------------------------------------------------------------------------------------------
  bool bDefaultDataInit()
  {
    int iAtoms  = iAtomsInit(),
        iBonds  = iBondsInit(),
        iAngles = iAnglesInit(),
        iDihs   = iDihedralsInit(),
        iImpros = iImpropersInit(),
        ivdWs   = ivdWInit();
      
    if(TALKATIVE > 0)
    {
      std::clog << "Initializing..."		 		     << std::endl;
      std::clog << "Generic atom descriptions read:     " << iAtoms  << std::endl;
      std::clog << "Generic bond descriptions read:     " << iBonds  << std::endl;
      std::clog << "Generic angle descriptions read:    " << iAngles << std::endl;
      std::clog << "Generic dihedral descriptions read: " << iDihs   << std::endl;
      std::clog << "Generic improper descriptions read: " << iImpros << std::endl;
      std::clog << "Generic vdW descriptions read:      " << ivdWs   << std::endl;
      std::clog << Textutils::pcBar1;
    }
    return !(AtomTypeNumbers.empty());
  }
//---------------------------------------------------------------------------------------------------
  unsigned int uGetNumCores() 
  {
    #ifdef WIN32
      SYSTEM_INFO Sysinfo;
      GetSystemInfo(&Sysinfo);
      return Sysinfo.dwNumberOfProcessors;
    #elif MACOS
      int nm[2];
      size_t Len = 4;
      uint32_t Count;

      nm[0] = CTL_HW; 
      nm[1] = HW_AVAILCPU;
      sysctl(nm, 2, &Count, &Len, NULL, 0);

      if(Count < 1) 
      {
	nm[1] = HW_NCPU;
	sysctl(nm, 2, &Count, &Len, NULL, 0);
	if(Count < 1) Count = 1;
      }
      return Count;
    #else
      return sysconf(_SC_NPROCESSORS_ONLN);
    #endif
  }
//---------------------------------------------------------------------------------------------------
}
