#ifndef MIXINGSETDECOPOSITION_H_HEADER_INCLUDED_BECEF599
#define MIXINGSETDECOPOSITION_H_HEADER_INCLUDED_BECEF599
#include "Global.h"
#include "OptModel.h"

struct BRANCHCALLBACKINFO;
class MixingSetDecomposition
{
public:
	MixingSetDecomposition(OptModel* _model);
	MixingSetDecomposition();
public:
	~MixingSetDecomposition(void);
	OptModel* model_p;
	OptModel model_mininization;

	static int Prev_Node_Num; // 2013, used in branch-and-cut
	static double B_B_Cut_Gen_Time; //April 2013, count cut generation time during branch and cut

	map< int, multimap<double,int> > mapViolatedScenario_Mins;

	int addSeperateCuts(string& cstType);
	vector<CUT> genertaeMixingSetCuts(int ind_violated, CUT_COEFICIENTS &alpha, double* nodex, BRANCHCALLBACKINFO*nodeInfo=NULL);
	double calculateMinForScenarioWithResConsByFormula(CUT_COEFICIENTS& c, int scenario_index);
	double calculateMinForScenarioByFormula(CUT_COEFICIENTS& costVec, int enforced_scenario_num);
	double calculateMinForScenario(CUT_COEFICIENTS& costVec, int scenario_number, BRANCHCALLBACKINFO* p_nodeInfo);
	double calculateMinForScenarioByLP(CUT_COEFICIENTS& costVec, int enforced_scenario_num, bool flag_new_cost_vec = false);
	CUT makeType9Cut(vector<int>& indexSet, CUT_COEFICIENTS& vec_mins, CUT_COEFICIENTS& vec_alpha, int k);
	CUT seperateMostViolatedMVPCut(CUT_COEFICIENTS& alpha, CUT_COEFICIENTS& h, double* nodex, int k);
	int runTest(string& instanceName, ofstream& log);

	// Sep 1 2012, rewrite separation algorithm
	CUT seperateMVPCut(CUT_COEFICIENTS& alpha, CUT_COEFICIENTS& h, double* nodex, int k);
	// Sep 10 2012, calculate the time spent on separation
	double TIME_CALCULATING_BETA;
	double TIME_CUT_GENERATION;

	// April 1, 2013, rewrite the mixing set, use formula to calculate beta. And need to do branch-and-cut
	vector<CUT> generateMixingSetCuts(int ind_violated, CUT_COEFICIENTS &alpha, double* nodex, BRANCHCALLBACKINFO*nodeInfo=NULL, int cutType =0);
	int runTest2013(string& instanceName, ofstream& log); 
	int callbackSeparateCuts(double* &nodex, int &nzcnt, double &rhs, int &sense, int * &cutind, double * &cutval);

};

#endif

