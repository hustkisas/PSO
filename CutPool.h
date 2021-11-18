#pragma once
#ifndef CUTPOOL_H_HEADER_INCLUDED_BECEF599
#define CUTPOOL_H_HEADER_INCLUDED_BECEF599

#include "Global.h"

class CutPool
{
public:
	CutPool(void);
	~CutPool(void);

	string strCutIndentifier;
	map<string, CUT> cuts;
	CutPool(string str):strCutIndentifier(str){};
	map<string, PAIR_DOUBLE_INT> cutTightness;
	void setCutIdentifier(string str);
	int updateCutPool(CPXENVptr& env, CPXLPptr& lp);
	void displayCutTightness();
	int getNumCuts();	

	int removeCutFromCPLEX(CPXENVptr& env, CPXLPptr& lp, string& vecNames);
	int removeInactiveCuts(CPXENVptr& env, CPXLPptr& lp, int numInactiveRounds);
	int removeInactiveCutsRandomly(CPXENVptr& env, CPXLPptr& lp, int numInactiveRounds);
	int removeInactiveCutsFromPool(int numInactiveRounds);
	int extractCutsFromModel(CPXENVptr& env, CPXLPptr& lp);




};

#endif
