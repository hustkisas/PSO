#ifndef WOLSEYSTRENGTHENING_H_HEADER_INCLUDED_BECEF599
#define WOLSEYSTRENGTHENING_H_HEADER_INCLUDED_BECEF599
#include "OptModel.h"
class WolseyStrengthening
{
public:
	WolseyStrengthening(OptModel* _model);
	WolseyStrengthening(){};
public:
	~WolseyStrengthening(void);
	static int iteration_cnt;
	static vector<double> beta;
	int timeLimit;
	int timeElapsed;
	void strengthenModel(int method);
	double calculateBetaByFormula(int& nzcnt,  int* rmatind, double *rmatval , double & lp_optimal);
	double calculateBetaBySolvingMIP(int& nzcnt,  int* rmatind, double *rmatval , double & lp_optimal,  OptModel::MODEL_TYPE modelType = OptModel::MODEL_TYPE_LP);
	void removeRedundantConsts(vector<string>& rowNums);
	void removeObsoleteVars(vector<string>& varNames);


	OptModel* model_p;
	int runTest(string& instanceName, ofstream& log, int method =0);
	int setTimeLimit(int timLmt);

	//July 07 2012, calculate beta using the k-th largest values
	double calculateBetaBySolvingMIP_rev1(int& nzcnt,  int* rmatind, double *rmatval , double & lp_optimal, OptModel::MODEL_TYPE modelType = OptModel::MODEL_TYPE_LP);
	double calculateMinForScenarioByLP( int enforced_scenario_num);

	//July 13 2012, adding 2-constraint set inequalities
	double calculateBetaBySolvingMIP_2ConstSet(int& nzcnt,  int* rmatind, double *rmatval , double & lp_optimal, OptModel::MODEL_TYPE modelType= OptModel::MODEL_TYPE_LP);
	double calculateScenarioLBBy2Const(int& nzcnt,  int* rmatind, double *rmatval , int enforced_scenario_num);
	map<int, OptModel*> LBmodels;
	int runTestImproved(string& instanceName, ofstream& log, int method);

	// July 16 2012, adding 2-constraint set inequalities after strengthening
	int runTestImprovedNew(string& instanceName, ofstream& log, int method);

	// July 19 2012, 
	void strengthenModelWith2ConstSet(int method);
	OptModel* min_model;

	// July 20 2012, add resource cuts
	void strengthenModelWithResourceCut(int method);
	int runTestResourceCut(string& instanceName, ofstream& log, int method);

	// Aug 2 2012, combine resource cuts with 2-const cuts
	void strengthenModelWithResourceCut2SetCut(int method);
	int runTestResourceCut2ConstCut(string& instanceName, ofstream& log, int method);
	int runTestResourceCut2ConstCutSeparation(string& instanceName, ofstream& log, int method);
	void strengthenModelWithResourceCut2SetCutSeparation(int method);
	double calculateBetaBySolvingMIP_rev2_addingViolated2ConstCuts(int& nzcnt,  int* rmatind, double *rmatval , double & lp_optimal, OptModel::MODEL_TYPE modelType= OptModel::MODEL_TYPE_LP);
	double calculateMinForScenarioByLP_withViolated2ConstCuts( int enforced_scenario_num);
	void strengthenModelNew( int method);

	// Sep 4 2012, paper first review, strengthen model with theoretical bound
	void strengthenModelWiTheoreticalBound();
	int runTestWiTheoreticalBound(string& instanceName, ofstream& log, int method);









};

#endif

