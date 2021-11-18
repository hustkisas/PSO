#pragma once
#ifndef VACCINATIONINSTANCEGENERATOR_HEADER_INCLUDED_BECEF599
#define VACCINATIONINSTANCEGENERATOR_HEADER_INCLUDED_BECEF599
#include "Global.h"
#include "stocc.h"
#include "OptModel.h"
using namespace std;

typedef  pair<  pair<string, string>  , double>   SCENARIO_TRIPLE;
typedef vector<SCENARIO_TRIPLE>   SCENARIO_ROW;

struct VacPolicy{
public:
	int policyID;
	int nVaccinatedChildren;
	int nVaccinatedAdults;
	int nVaccinatedElderly;
	int nChildren;
	int nAdults;
	int nElderly;
	int idx;
	VacPolicy(int pID, int nVacChi, int nVacAd, int nVacEld, int nChi, int nAd, int nEld):policyID(pID),nVaccinatedChildren(nVacChi),
		nVaccinatedAdults(nVacAd),nVaccinatedElderly(nVacEld),nChildren(nChi),nAdults(nAd),nElderly(nEld),idx(0){};

};

struct  HouseHold{
public:
	int hhID;
	int size;
	int nChildren;
	int nAdults;
	int nElderly;
	float proportion;
	vector <VacPolicy> vecPolicies;
	HouseHold(int _hhID, int _size, int _nChi, int _nAd, int _nEld, float _prop):hhID(_hhID),size(_size), nChildren(_nChi),nAdults(_nAd),nElderly(_nEld),proportion(_prop){};
};



class VaccinationInstanceGenerator
{
public:
	VaccinationInstanceGenerator(int _s, int _nr, int _n, int _m, int _k, double _epsilon = 0);
	VaccinationInstanceGenerator(){};
	~VaccinationInstanceGenerator(void);

	CPXENVptr cpxenv;
	CPXLPptr  cpxlp;


	int seed;
	int n;
	int nr;
	int m;
	int k;
	double epsilon;

	int bigM;;

	StochasticLib1 randGen_eps;
	StochasticLib1 randGen_m;
	StochasticLib1 randGen_b;
	StochasticLib1 randGen_s;
	StochasticLib1 randGen_t;




	map<int, HouseHold>  hhs; 
	map<int, VacPolicy> vacPolicies;
	double avg_hh_size;

	// SMPS
	vector< pair<double, SCENARIO_ROW>  > scenarios;


	int initilizeCPXlp(string modelName);
	int initilizeCPXenv();

	float get_u();
	float get_s();
	float get_a(HouseHold hh, VacPolicy p);


	bool readFamilyTypes();
	bool createEquaConstraints();
	bool createObjective(OptModel::MODEL_TYPE modelType);
	bool  createSenarioConstraints();
	bool  createCardConstraints();
	bool generateModel(string& modelName, OptModel::MODEL_TYPE modelType);
	bool complementModel();
	double getOriginalObjValueForComplementedModel(CPXENVptr& cpx, CPXLPptr& lp);
	
	
	// SMPS
	void ReadSMPS_StoFile(string& fileName);
	void ReadSMPS_CorFile(string& fileName, OptModel::MODEL_TYPE modelType);
	void BuildModelFromSMPSFiles(string& fileName, OptModel::MODEL_TYPE modelType);
};



#endif

