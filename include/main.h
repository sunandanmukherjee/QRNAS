#pragma once
#ifndef  MAIN_H
#define  MAIN_H

#include <iostream>
#include <fstream>
#include <cstring>

#include <list>

#include "tconsts.h"
#include "textutils.h"
#include "tmolecule.h"

std::list<TMolecule> Molecules;
double dAverageNorm;

bool   bParseArgs(int iArgC, char* ppcArgV[]);
void   DisplayInfo();

bool   bReadPDB(const std::string& _sFileName);
bool   bWritePDB(const std::string& _sFileName);

void   MakeAllExternNonbonds();
double dTotalEnergy(double& _dMolEnergy_);
bool   bMinimize();
double dGoldenMinimStep(double& _dMolEnergy_);
void   CalcAverageNorm();
bool   bShiftAllMolecules(double _dFactor);

#endif /*MAIN_H*/
