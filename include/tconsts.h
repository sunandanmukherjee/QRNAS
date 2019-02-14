#pragma once
#ifndef  TCONSTS_H
#define  TCONSTS_H

#include <fstream>
#include <iostream>
#include <sstream>
#include <cstring>
#include <map>
#include <vector>
#include <list>

#include <cmath>

#ifdef _WIN32			// For multithreading. Anyway, only pthreads are supported now.
#include <windows.h>
#elif MACOS
#include <sys/param.h>
#include <sys/sysctl.h>
#else
#include <unistd.h>
#endif

#include "params.h"
#include "distrestr.h"
#include "bprestr.h"

/**************************************************************************
 * A namespace containing all the constants
 * (mathematical, physical, as well as user-specified parameters).
 * It was designed to be an abstract class (that's the reason for 'T' 
 * in its name), but namespace turned out to be equally good.
 **************************************************************************/

namespace TConsts
{
  extern std::map<std::string, int> AtomTypeNumbers;	// Associative table that translates atom names (C1, N*, HW, etc) 
							// of type std::string into internal numbers (int).
  extern std::vector<vdWParam> GeneralvanderWaalsData;
  extern std::vector<std::vector<BondParam> > GeneralBondsData;
  extern std::vector<std::vector<std::vector<AngleParam> > > GeneralAnglesData;
  extern std::vector<std::vector<std::vector<std::vector<ImproperParam> > > > GeneralImpropersData;
  extern std::vector<std::vector<std::vector<std::vector<std::list<DihedralParam> > > > > GeneralDihedralsData;

  extern std::vector<TDistRestrDescr> DistanceRestraintsDescriptors;
  extern std::vector<TBasePairRestrDescr> BPRestraintsDescriptors;
  
  extern const long double PI;		// Probably 'long' keyword does nothing.
  extern const long double SQRT2;
  extern const double COULOMB;
  extern const double kB;
  extern const double PVUnit;
  extern const double EpsilonRNA;
  extern const double EpsilonWater;
  
  extern std::string ATOMSFILE, VDWFILE, BONDSFILE, ANGLEFILE, IMPROFILE, DIHEDFILE;
  extern std::string ATOMSDIR, VDWDIR, BONDSDIR, ANGLEDIR, IMPRODIR, DIHEDDIR, PDBDIR, PREPDIR;
  extern std::string INPUTPDB, OUTPUTPDB, PARAMSFILE, RESTRFILE, PROGRAMPATH;
  
  extern unsigned int NUMTHREADS;

  extern unsigned int NSTEPS;
  
/**************************************************************************
 * TALKATIVE={0,1,2} - whether to report all (2), most (1)
 * or nearly none (0) activities to std::clog
 **************************************************************************/
  extern unsigned int TALKATIVE;
/**************************************************************************
 * WRITEFREQ describes how often (in steps) a PDB file should be written.
 **************************************************************************/
  extern unsigned int WRITEFREQ;
/**************************************************************************
 * PEDANTICMOL={0,1} defines whether bReadMolecule==false is 
 * a reason to stop further processing
 **************************************************************************/
  extern unsigned int PEDANTICMOL;
/**************************************************************************
 * Cut-off distance for dispersive interactions.
 **************************************************************************/
  extern double CUTOFF;
  extern double CUTOFF_SQR;
  extern bool   TRAJECTORY;
  extern bool   USEBORN;	// Born generalized solvent.
  extern bool   SSCONSTR;	// Constraints on secondary structure?
  extern bool   SSDETECT;
  extern bool   HBONDS;
  extern bool   ELECTR;
  extern bool   VANDERWAALS;
  extern bool   BLOWGUARD;
  
  extern std::string SECSTRUCT;

  extern bool   POSRESTRAINTS;
  
  extern double MAXSTEP;
  extern double MINSTEP;

  bool bSetPaths();
  bool bUserInit(const std::string& _sFileName=PARAMSFILE);
  bool bRestrInit(const std::string& _sFileName=RESTRFILE);
  void CheckAndDisplayParams();
  
  int iAtomsInit(const std::string& _sFileName=ATOMSFILE);
  
  int ivdWInit(const std::string& _sFileName=VDWFILE, 
	       std::vector<vdWParam>& _vanderWaalsData_ = GeneralvanderWaalsData);
	       
  int iBondsInit(const std::string& _sFileName=BONDSFILE, 
		 std::vector<std::vector<BondParam> >& _BondsData_ = GeneralBondsData);
		 
  int iAnglesInit(const std::string& _sFileName=ANGLEFILE, 
		  std::vector<std::vector<std::vector<AngleParam> > >& _AnglesData_ = GeneralAnglesData);
		  
  int iImpropersInit(const std::string& _sFileName=IMPROFILE, 
		     std::vector<std::vector<std::vector<std::vector<ImproperParam> > > >& _ImpropersData_
			= GeneralImpropersData);
			
  int iDihedralsInit(const std::string& _sFileName=DIHEDFILE, 
		     std::vector<std::vector<std::vector<std::vector<std::list<DihedralParam> > > > >& _DihedralsData_
			= GeneralDihedralsData);
			
  bool bDefaultDataInit();
  bool bDefaultDataCleanUp();
  inline bool bIsNan(const double _dVal){volatile double _dVolatVal = _dVal; return _dVolatVal != _dVolatVal;}
  
  unsigned int uGetNumCores();
}

#endif /*TCONSTS_H*/
