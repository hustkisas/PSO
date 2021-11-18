#ifndef OPTMODEL_H_HEADER_INCLUDED_BECEF599
#define OPTMODEL_H_HEADER_INCLUDED_BECEF599

#include "Global.h"
#include "DataModel.h"
#include "CutPool.h"

struct BRANCHCALLBACKINFO;

class OptModel{
public:
	string modelName;
	CPXENVptr cpxenv;
	CPXLPptr  cpxlp;
	CPXLPptr  cpxlp_backup;
	DataModel* dm_p;
	bool is_a_copy;

	vector<int> idx_X;
	vector<int> idx_Z; 
	vector<int> idx_S;
	map<int, int> map_scenario_rownum; 

	map<string, int> var_name_idx_table; 
	map<int, string> idx_var_name_table; 

	map<int, int> map_z_idx_row_num;
	map<int, int> map_row_num_z_idx;

	// solution statistics
	double MIP_optimal_value;
	int MIP_solving_time;
	int MIP_BB_nodes;
	double MIP_GAP;
	double MIP_Root_Obj_Val;
	double MIP_INCUMBENT_OBJ_VAL;
	int LP_solving_time;
	double LP_optimal_value;

	// solution improvement
	double Prev_MIP_optimal_value;
	double Prev_LP_optimal_value;

	//callback info
	int num_one_at_curr_node;
	int num_rootnode_cut;
	int num_IC_rootnode;
	int num_MIX_rootnode;

	// cut pool
	CutPool ICCutPool;
	CutPool MixCutPool;
	CutPool MIRCutPool;
	CutPool CHCutPool;
	CutPool GeneralCutPool;
	int updateCutPool(string cutType = string("ALL"));
	int storeCuts(string cutType);
	int addCutPoolToModel(string poolName, int asConst);
	int removeInactiveCuts(string cutType, int numRoundsInactive);
	void setGeneralCutPoolIdentifier(string identifier);




	enum MODEL_TYPE {MODEL_TYPE_MIP, MODEL_TYPE_LP};
	enum CUT_TYPE{CUT_AS_FORMULATION, CUT_AS_USERCUT, CUT_AS_LAZYCONSTRAIN};

	MODEL_TYPE modelType;


	OptModel(){};
	OptModel(string _modelName, int CPX_SCRIND = 0);
	OptModel(CPXENVptr cpxenv, CPXLPptr  cpxlp);
	OptModel(CPXENVptr _cpxenv, string _modelName=string("copy"));
	OptModel(OptModel& _model);



	void ClearModel(string& _modelName);

	~OptModel();

	int initilizeCPXenv(int CPX_SCRIND = 0);
	int initilizeCPXlp(string modelName);
	int populateModelWithoutSlackVars( DataModel& dm, OptModel::MODEL_TYPE modeltype, char cardinalitySign, bool resConst, bool negativeBound);

	void writeModel2File(string modelName);
	void setDataModel(DataModel* dm);



	CPXENVptr& getCPXEnv(void);
	CPXLPptr& getCPXLP(void);
	vector<int>& getIndex_X(void);
	vector<int>& getIndex_Z(void);
	vector<int>& getIndex_S(void);
	void writeModel2File();

	int solve(void);
	int setCPXParameter(int whichparam, int newvalue);
	int setCPXDoubleParameter(int whichparam, double newvalue);
	vector<PAIR_INDEX_VALUE> getSolutionValues(int varType=1);
	double getDoubleSolQuality(int paraName);

	int mapVarIndexName(void);
	void mapScenarioRownum();
	void mapZIdxRowNum();

	void printSolution(int varType=0/*0:all 1:x 2:z*/);
	void printObjectiveValue(void);
	int changeLP2MIP(int changeType = 0); //0: only z; 1: both x and z
	int printTableau(string tableauFileName = "tableau.csv", CPXCENVptr env = NULL, CPXLPptr  lp=NULL);
	int getNumRows(void);
	int getNumCols(void);
	bool isIntegerFeasible(int varRange = 0 /*0: z; 1: x and z*/);

	int changeMIP2LP(void);
	int changeProbTypeMIP2LP();
	int populateModelFromFile(string lpFileName);
	int getNumAllCuts();
	void mapVarNameIndex();

	int populateModelWithReformulation( DataModel& dm, OptModel::MODEL_TYPE modeltype);
	int addReformAsCut( DataModel& dm, OptModel::CUT_TYPE cutType);



	//callback related
	int prepareCallbackCutCoefsReform(CPXCENVptr env,  int branch_index, int &nzcnt, double &rhs, int &sense, int * &cutind, double * &cutval);
	int prepareCallbackCutCoefsReformNew(CPXCENVptr env,  int branch_index, int &nzcnt, double &rhs, int &sense, int * &cutind, double * &cutval);


	int setCutCallBackFunction(int (CPXPUBLIC *cutcallback) (CALLBACK_CUT_ARGS), void *cbhandle);
	int setNodeCallBackFunction(int(CPXPUBLIC *nodecallback) (CALLBACK_NODE_ARGS), void *cbhandle);
	int setBranchCallBackFunction(int(CPXPUBLIC *branchcallback) (CALLBACK_BRANCH_ARGS), void *cbhandle);
	int setSolveCallBackFunction(int(CPXPUBLIC *solvecallback) (CALLBACK_SOLVE_ARGS), void *cbhandle);
	int setIncumbentCallBackFunction(int (CPXPUBLIC *incumbentcallback) (CALLBACK_INCUMBENT_ARGS), void *cbhandle);

	int getNumCutsInCutPool(string cutType);
	int findRowNumByRowName(string& rowName);


	multimap<double, int> getAllFractionalZ();
	void printSolution2File(string fileName, int varType=0/*all solutions*/, CPXCENVptr env=NULL, CPXLPptr  lp=NULL);
// 	int prepareCallbackCutCoefsIntersection(CPXENVptr& env,  CPXLPptr & nodelp_p , int &nzcnt, double &rhs, int &sense, int * &cutind, double * &cutval);

	CUT substituteSlackVarInCut(CPXENVptr &env,  CPXLPptr & lp, CUT_COEFICIENTS& cutCoefs);
	char ** getColumnNames(CPXCENVptr &env , CPXCLPptr &lp);


	//callback static functions
	static bool CallbackIsHeuristicSolution(CPXENVptr &env, CPXLPptr& lp, void * cbdata, int wherefrom, OptModel* model);
	bool callbackAddACut(void * cbdata, int wherefrom, CUT& cut);
	int callbackAddReverseCHCuts(void * cbdata, int wherefrom,int branch_index);



	// Decomposition related functions
	int initilizeMasterModel( DataModel& dm, OptModel::MODEL_TYPE modeltype, char cardinalitySign);
	double IsXValidForScenarioZ(int ind_scenario, double* x, double z);
	double IsXValidForScenarioZ(int ind_scenario, double* x);
	bool IsXFeasibleSolution(double* x);
	bool addCallbackSeperateCuts(CPXCENVptr &env, void * cbdata, int wherefrom, OptModel* model);
	bool addCutToMasterModel(CUT& cut, CPXCENVptr &env, void * cbdata, int wherefrom);
	double testCutValidity(CUT& cut, double* nodex);
	vector<int> checkViolatedSenarioByX(double* nodex);


	bool addAConstToCPXModel(CUT& cut, string& constName, int cstType  = 0);
	int removeConstraintFromCPLEX(string& constraintName);
	int removeVariableFromCPLEX(string& varName);

	// hard instance
	int initilizeMasterModelWithResConst( DataModel& dm, OptModel::MODEL_TYPE modeltype, char cardinalitySign, bool negativeBound);

	static int CUT_SELECTION_METHOD;
	static int NUM_CUT_PER_ITER ;
	static int NUM_ITERATIONS;

	// MIR cuts
	int addAllMIRCuts();

	// Convex Hull 
	int AddConvexHullToModel(int cstType = 0);

	// mixing set with more cost vectors
	vector< pair<string, int> > readRowNamesFromModel();
	map< string, int> readRowNamesMapFromModel();

	int addSeperateCuts(string cstType= "MIR");
	vector<int> checkViolatedConstrainsByX(string& cstType, double* nodex);
	vector<int> checkSortedViolatedConstrainsByX(string& cstType, double* nodex);

	double IsXValidForConstrain(int rowNum, double* x);
	CUT_COEFICIENTS readAlphaFromConstrain(int rowNum);
	
	// mixing set on child nodes
	int callbackAddMixingSet(void * cbdata, int wherefrom, BRANCHCALLBACKINFO*nodeInfo);
	CUT readScenarioFromModel(int scenarioNum);
	int findSceNumByRowNum(int rowNum);


	// network flow model
	vector<int> idx_T;
	int populateModelNetworkFow( DataModel& dm, OptModel::MODEL_TYPE modeltype, char cardinalitySign);

	// wacky formulation 
	map<string, int> idx_map_X;
	map<string, int> idx_map_Z;
	map<string, int> idx_map_W;
	int populateWeckyModel( DataModel& dm, OptModel::MODEL_TYPE modeltype, char cardinalitySign, bool resConst);


	// record integer solutions that have been found
	vector<PAIR_VALUE_INDEX> feasibleSols;

	// aggregation heuristics
	vector<VEC_INT> aggGroups;
	bool Agg_generateAggregationModel();
	double calculateDistance(CUT &scenario1, CUT &scenario2, int method=1);
	double Agg_populateAggregationModel(DataModel& dm, OptModel::MODEL_TYPE modeltype, char cardinalitySign, bool resConst, int method);
	int Agg_buildInitialAggGroups(DataModel& dm, int groupSize, int normType);
	int Agg_buildInitialAggGroupsNew(DataModel& dm, int groupSize, int normType);
	int Agg_getViolatedConsts(multimap<double, pair<int,int> >& scenarios_picked, int method);
	int Agg_updateGroups(multimap<double, pair<int,int> > violated_indicies, int num_to_pick);
	int Agg_updateGroupsFromSmallestViolation(multimap<double, pair<int,int> > violated_indicies, int num_to_pick);
	int Agg_populatePartialFixedModel( DataModel& dm, OptModel::MODEL_TYPE modeltype, char cardinalitySign, bool resConst, multimap<double, pair<int,int> > violated_indicies);
	int AddObjectiveCut( double objVal, char sign);
	double IsXValidForScenarioZNormalized(int ind_scenario, double* x);
	bool IsCurrentSolutionFeasible();
	void Agg_recordGroupInfo(multimap<double, pair<int, int> >& violated_indicies, int num_to_pick);


	double calculateLBbyEuclideanDistance();

public:
	double startingtime;
	double startTimer(void);
	double getElapsedTime(void);
	CUT_COEFICIENTS getObjective();

	string getVarNameByIndex(int index);

	// Sep 7 2012, CKVLP first review
	static double MIN_CUT_VIOLATION;

	// April 7 2013, CKVLP second review
	int removeConstraintsFromCPLEXbyKeyWord(string& key);
};

#endif

