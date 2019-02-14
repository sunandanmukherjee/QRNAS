#pragma once
#ifndef  TATOM_H
#define  TATOM_H

#include <cmath>		// For sqrt().

#include <istream>		// For atom's self-read-in.
#include <ostream>		// For atom's self-print-out.
#include <sstream>		// For reading info from strings by >> operators.
#include <fstream>		// For atom's self-recognition.
#include <iostream>		// For general output.	
#include <iomanip>		// For formatting of output.

#include <cstring>
#include <map>
#include <algorithm>
#include <vector>

#include <cstdlib>		// For exit on error!

#include "tvector.h"
#include "textutils.h"
#include "tconsts.h"

#include "params.h"

class TElectr;
class TLennard;

class TAtom
{
  private:
    double dX, dY, dZ;
    std::string sPDB_tag,	// All the crap from PDB.
		sResidue,	// Residue name, with spaces cropped.
		sResIdent,	// Residue 'identity', i.e. number, insertion code, etc.
		sName;	        // Atom's name within its residue.

    std::string sType;		// Atom's type, same as used by AMBER in parm99.dat
    
    int iResidueNumber;
    char cChainCode;

    double dCharge;		// Atomic partial charge, in elementary charge units.
    vdWParam vdWInfo;		// Atom's van der Waals data.

    bool bRecognizeCoords();	// Reads coordinates from sPDB_tag.
    bool bRecognizeNumber();	// Reads iNum from sPDB_tag.
    bool bRecognizeType();	// Finds atomic dCharge, van-der-Waals params and sType.
    bool bRecognizeRestrain();

  public:
    double dPosRestr;		// How much can we move it?
    int iNum;			// Atom's number, as given in PDB
    double dBornRadius;		// Atom's Born radius (for treatment of solvent)
    TVector vecEnergyGradient;	// Gradient of total energy with respect to this atom's coords (dX, dY, dZ).
    
    friend TAtom& operator+=(TAtom& _atAtom, const TVector& _vecVector);	// Shifts the atom by vector _vecVector.
    friend TAtom& operator-=(TAtom& _atAtom, const TVector& _vecVector);	//
    
    friend TVector vecMakeVector(const TAtom& _atAtom1, const TAtom& _atAtom2);	   // Useful, and (pity) it cannot be
    friend TVector vecMakeVector(const TAtom* _patAtom1, const TAtom* _patAtom2);  // defined in any more elegant way.
    friend TVector vecMakeVector(const TAtom* _patAtom1, const double _dX, const double _dY, const double _dZ);
    friend double dDistance_SQR(const TAtom& _atAtom1, const TAtom& _atAtom2);	   // Also useful and needs no calculation of
    friend double dDistance_SQR(const TAtom* _patAtom1, const TAtom* _patAtom2);   // vector between atoms (speedup).
    friend double dDistance_SQR(const TAtom& _atAtom1, const double _dX, const double _dY, const double _dZ);
    friend double dDistance_SQR(const TAtom* _patAtom1, const double _dX, const double _dY, const double _dZ);

    TAtom(std::string _sPDB_tag = ""): dX(-1.), dY(-1.), dZ(-1.), 		// Default Ctor
	    sPDB_tag(_sPDB_tag), dPosRestr(0.), dBornRadius(0.) {}

    bool bAssignTag(std::string _sLine);	// Setter for sPDB_tag
    bool bWrite(std::ostream& _osOut) const;	// Writes atom's data to given stream (in PDB format)
    bool bRecognizeResidue();			// Read sResidue and sResIdent from sPDB_tag
    bool bRecognizeName();	// Reads sName from sPDB_tag.
    bool bRecognize();		// Self-assignment of all data, based on PDB_tag + some info from TConsts
				// It just calls all other RecognizeXXXX methods.

    std::string sGetResIdent() const;				// Getter for residue number and chain letter code, etc.
    bool 	bSetResIdent(const std::string& _sNewResIdent);	// Setter for residue number and chain letter code, etc.
    std::string sGetResName() const;		      		// Getter for residue name; needs an atom to be recognized!
    bool        bSetResName(std::string _sNewResName);   	// Setter for atom's residue name.
    std::string sGetName() const;				// Obvious
    void        SetName(std::string _sName);			// Obvious
    std::string sGetType() const;				// ...
    void 	 GetCoords(double& _dX_, double& _dY_, double& _dZ_) const;	// Writes atom's coords in the given references.
    double      dGetvdWRadius() const;
    std::string sGetPDB_tag() const;		      		// Getter for sPDB_tag.
    
    int  iGetResidueNumber() const;
    char cGetChainCode() const;
    
    bool bTransform(TAtom _at1,      TAtom _at2,      TAtom _at3,
		    TAtom _at1Prime, TAtom _at2Prime, TAtom _at3Prime); // A rarely-used method. Translates and rotates *this in a way that
									// _at1 becomes _at1Prime, _at2 --- _at2Prime, _at3 --- _at3Prime    

    friend TElectr  elMakeElectrPair(TAtom& _atAtomA, TAtom& _atAtomB);  // extern function - find it in telectr.cpp
    friend TLennard ljMakeLJPair(    TAtom& _atAtomA, TAtom& _atAtomB);  // extern function - find it in tlennard.cpp
};

#endif /*TATOM_H*/