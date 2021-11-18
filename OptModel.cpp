#include "OptModel.h"
#include <algorithm>
#include <ilcplex/cplex.h>
#include "callbacks.h"


int OptModel::CUT_SELECTION_METHOD = 1;
int OptModel::NUM_CUT_PER_ITER =1;
int OptModel::NUM_ITERATIONS = 10;
extern string CAT_NAME;
double OptModel::MIN_CUT_VIOLATION = 0.01;

OptModel::OptModel(string _modelName, int CPX_SCRIND){
	modelName = _modelName;
	num_one_at_curr_node = 0;
	num_rootnode_cut = 0;
	num_MIX_rootnode = 0;
	num_IC_rootnode = 0;
	is_a_copy = false;

	MIP_optimal_value = 0;
	MIP_GAP = -1;
	MIP_Root_Obj_Val = 0;
	Prev_MIP_optimal_value = 0;
	LP_optimal_value =0 ;
	Prev_LP_optimal_value = 0;

	ICCutPool.setCutIdentifier("IC");
	MixCutPool.setCutIdentifier("MIX");
	MIRCutPool.setCutIdentifier("MIR");
	CHCutPool.setCutIdentifier("CH");


	initilizeCPXenv(CPX_SCRIND);
	initilizeCPXlp(modelName);
}

OptModel::OptModel(OptModel& _model)
{
	int status;
	initilizeCPXenv();
	idx_X = _model.idx_X;
	idx_Z = _model.idx_Z; 
	idx_S = _model.idx_X;
	map_scenario_rownum = _model.map_scenario_rownum; 
	cpxlp = CPXcloneprob( cpxenv, _model.getCPXLP(), &status);

}



int OptModel::initilizeCPXenv( int CPX_SCRIND)
{
	cpxenv = NULL;
	int           status = 0;

	/* Initialize the CPLEX environment */
	cpxenv = CPXopenCPLEX (&status);

	if ( cpxenv == NULL ) {
		char  errmsg[1024];
		fprintf (stderr, "Could not open CPLEX environment.\n");
		CPXgeterrorstring (cpxenv, status, errmsg);
		fprintf (stderr, "%s", errmsg);
		return 0;
	}

	/* Turn on output to the screen */
	status = CPXsetintparam (cpxenv, CPX_PARAM_SCRIND, CPX_SCRIND);
	if ( status ) {
		fprintf (stderr, 
			"Failure to turn on screen indicator, error %d.\n", status);
		return 0;
	}

	status = CPXsetintparam(cpxenv, CPX_PARAM_THREADS, 1);
	status = CPXsetintparam(cpxenv, CPX_PARAM_CLOCKTYPE, 1);

	return 1;
}


int OptModel::initilizeCPXlp(string modelName)
{
	int status = 0;
	if (cpxenv == NULL)
	{
		cout<< "CPXenv hasn't been initialized yet when initializing CPXlp. Now Initilize CPXenv"<<endl;
		initilizeCPXenv();
	}

	cpxlp = NULL;

	/* Create the problem. */
	cpxlp = CPXcreateprob (cpxenv, &status, modelName.c_str());
	if ( cpxlp == NULL ) {
		fprintf (stderr, "Failed to create LP.\n");
		return 0;
	}

	return 1;

}


void OptModel::ClearModel(string& _modelName){
	modelName = _modelName;
	num_one_at_curr_node = 0;
	num_rootnode_cut = 0;
	num_MIX_rootnode = 0;
	num_IC_rootnode = 0;

	MIP_optimal_value = 0;
	Prev_MIP_optimal_value = 0;
	LP_optimal_value =0 ;
	Prev_LP_optimal_value = 0;

	ICCutPool.setCutIdentifier("IC");
	MixCutPool.setCutIdentifier("MIX");
	MIRCutPool.setCutIdentifier("MIR");
	CHCutPool.setCutIdentifier("CH");

	CPXfreeprob(cpxenv, &cpxlp);
//	initilizeCPXenv();
	initilizeCPXlp(_modelName);
}


int OptModel::populateModelWithoutSlackVars( DataModel& dm, OptModel::MODEL_TYPE modeltype, char cardinalitySign, bool resConst, bool negativeBound)
{
	setDataModel(&dm);
	idx_X.clear();
	idx_Z.clear();
	idx_S.clear();
	double *obj     = NULL;
	double *lb      = NULL;
	double *ub      = NULL;
	char   *ctype   = NULL;
	int    *rmatind = NULL;
	double *rmatval = NULL;
	char** colname = NULL;

	int    status = 0; 
	int num_cols = DataModel::n + DataModel::m;

	CPXchgobjsen (cpxenv, cpxlp, 1); /* Minimization problem */

	/* Allocate colcnt-sized arrays */

	obj     = (double *) malloc (num_cols* sizeof(double));
	lb      = (double *) malloc (num_cols * sizeof(double));
	ub      = (double *) malloc (num_cols * sizeof(double));
	ctype   = (char *)   malloc (num_cols * sizeof(char));
	colname = (char**) malloc(num_cols*sizeof(char*));

	if ( obj     == NULL ||
		lb      == NULL ||
		ub      == NULL ||
		ctype   == NULL ) {
			fprintf (stderr, "Could not allocate colcnt arrays\n");
			status = CPXERR_NO_MEMORY;
			return 0;
	}

	string* strColNames = new string[num_cols];
	for (int i =0; i < num_cols; i ++)
	{
		if (i < DataModel::n )
		{
			obj[i] = dm.costVec[i].first;                      // fill in obj coefficients
			idx_X.push_back(i);                        // store index for x variables
			strColNames[i] = "x"+ToString(i);
			colname[i] = new char[8];
			strcpy(colname[i], strColNames[i].c_str()); // make name for x
			ctype[i] = 'C';
		}
		else{
			obj[i] = 0;
			idx_Z.push_back(i);
			strColNames[i] = "z" + ToString(i-DataModel::n);
			colname[i] = new char[8];
			strcpy(colname[i], strColNames[i].c_str()); // make name for z
			if (modeltype == MODEL_TYPE_LP)
				ctype[i] = 'C';
			else if( modeltype == MODEL_TYPE_MIP)
				ctype[i] = 'B';
		}
	}
	delete []strColNames;

	for (int i =0; i < num_cols; i ++)       // fill in upper bound
	{
		if (i<DataModel::n)
			ub[i] = CPX_INFBOUND;
		else 
			ub[i] = 1;
	}

	for (int i =0; i < num_cols; i ++)       // fill in upper bound
	{
		if (i<DataModel::n)
			lb[i] = 0;
		else
			lb[i] = - CPX_INFBOUND;
	}

	if (modeltype == MODEL_TYPE_LP)
	{
		if (negativeBound)
			status = CPXnewcols (cpxenv, cpxlp, num_cols, obj, lb, ub, NULL, colname);
		else 
			status = CPXnewcols (cpxenv, cpxlp, num_cols, obj, NULL, ub, NULL, colname);
	}
	else if (modeltype == MODEL_TYPE_MIP)
	{
		status = CPXnewcols (cpxenv, cpxlp, num_cols, obj, NULL, ub, ctype, colname);
	}
	if ( status ) {
		fprintf (stderr, "Could not add new columns.\n");
		return 0;
	}

	modelType = modeltype;

	/* adding rows to model */
	for (int i = 0; i < DataModel::m; i++) 
	{

		int totaNumCols = DataModel::n +1 ;
		rmatind = (int * )   malloc (totaNumCols * sizeof(int));
		rmatval = (double *) malloc (totaNumCols * sizeof(double));

		VEC_MTX_ROW rowEntries = dm.Arows[i];  // fill in the entries for a row
		double rhs[1];
		char   sense[1];

		switch (dm.senses[i])
		{
		case 'L':
			for (int j =0; j <DataModel::n; j ++)
			{
				PAIR_INDEX_VALUE entryPair = rowEntries[j];
				rmatind[j]= entryPair.first;
				rmatval[j] = - entryPair.second; 
			}
			rmatind[DataModel::n] = DataModel::n + i;
			rmatval[DataModel::n] = dm.rhs[i];
			rhs[0] = -dm.rhs[i];
			sense[0] = 'G';

			break;
		case 'G':
			for (int j =0; j <DataModel::n; j ++)
			{
				PAIR_INDEX_VALUE entryPair = rowEntries[j];
				rmatind[j] = entryPair.first;
				rmatval[j]= entryPair.second; 
			}
			rmatind[DataModel::n] = DataModel::n + i;
			rmatval[DataModel::n] = dm.rhs[i];
			rhs[0] = dm.rhs[i];
			sense[0] = 'G';

			break;
		case 'E':
			for (int j =0; j <DataModel::n; j ++)
			{
				PAIR_INDEX_VALUE entryPair = rowEntries[j];
				rmatind[j] = entryPair.first;
				rmatval[j] = entryPair.second; 
			}
			rmatind[DataModel::n] = DataModel::n + i;
			rmatval[DataModel::n] = dm.rhs[i];
			rhs[0] = dm.rhs[i];
			sense[0] = 'E';

			break;
		}
		int    rmatbeg[1];
		rmatbeg[0] = 0;
		char** rowName = new char*; 
		string cstName =  /* i == dm.Arows.size()-1 ?  string("RES")  : */string("SN") +ToString(i);
		char* pStr = new char[20];
		strcpy(pStr, cstName.c_str());
		rowName[0] = pStr;
		status = CPXaddrows (cpxenv, cpxlp, 0, 1, totaNumCols , rhs,
			sense, rmatbeg, rmatind, rmatval,
			NULL, rowName);
		if ( status ) {
			fprintf (stderr, "Could not add new rows.\n");
			return 0;
		}
		delete rowName;
		free_and_null((char **)&rmatind); 
		free_and_null((char **)&rmatval); 

	}

	/* add resource constraint*/
	if (resConst)
	{
		rmatind = (int * )   malloc (DataModel::n * sizeof(int));
		rmatval = (double *) malloc (DataModel::n * sizeof(double));
		for (int i = 0; i < DataModel::n ; i++)
		{
			rmatind[i] = i;
			rmatval[i] = 1;
		}
		double rhs[1];
		rhs[0] = 1;
		char   sense[1];
		sense[0] = 'E';
		int    rmatbeg[1];
		rmatbeg[0] = 0;
		char ** rowName = new char*;
		rowName[0] = "RES";

		status = CPXaddrows (cpxenv, cpxlp, 0, 1, DataModel::n , rhs,
			sense, rmatbeg, rmatind, rmatval, NULL, rowName);
	}
	/* adding cardinality constraint */
	lb[0] = 0;
	ub[0] = 0; // exactly k violated
	int totaNumCols = DataModel::m;

	/* prepare row data*/

	rmatind = (int * )   malloc (totaNumCols * sizeof(int));
	rmatval = (double *) malloc (totaNumCols * sizeof(double));

	for (int j =0; j <DataModel::m; j ++)
	{
		rmatind[j] = idx_Z[j];
		rmatval[j] = 1;
	}
	int    rmatbeg[1];
	double rhs[1];
	char   sense[1];
	rhs[0] = DataModel::k;
	sense[0] = cardinalitySign;
	rmatbeg[0] = 0;
	char ** rowName = new char*;
	rowName[0] = "CARD";
	status = CPXaddrows (cpxenv, cpxlp, 0, 1, DataModel::m, rhs,
		sense, rmatbeg, rmatind, rmatval,
		NULL, rowName);
	delete rowName;
	if ( status ) {
		fprintf (stderr, "Could not add new rows.\n");
		return 0;
	}

	string fileName ;
	if (modeltype == MODEL_TYPE_LP)
		fileName = modelName + "R.lp" ;
	else if (modeltype == MODEL_TYPE_MIP)
		fileName = modelName + "MIP.lp";
	writeModel2File(fileName.c_str());
	//CPXwriteprob(cpxenv, lp, "initialLPRelaxation.lp", NULL);

	free_and_null((char **)&obj);
	free_and_null((char **)&lb);
	free_and_null((char **)&ub);
	free_and_null((char **)&ctype);
	free_and_null((char **)&colname); 
	free_and_null((char **)&rmatind); 
	free_and_null((char **)&rmatval); 

	return 0;
}



int OptModel::populateWeckyModel( DataModel& dm, OptModel::MODEL_TYPE modeltype, char cardinalitySign, bool resConst)
{
	setDataModel(&dm);
	idx_X.clear();
	idx_Z.clear();
	idx_S.clear();
	double *obj     = NULL;
	double *lb      = NULL;
	double *ub      = NULL;
	char   *ctype   = NULL;
	int    *rmatind = NULL;
	double *rmatval = NULL;
	char** colname = NULL;

	int    status = 0;
	int n = DataModel::n; int K = DataModel::k; int m = DataModel::m;
	int num_x = n; int num_w = (1+K)*n*(m-K); int num_z = m-K + K*(m-K);
	int num_col = num_x + num_w + num_z ;

	CPXchgobjsen (cpxenv, cpxlp, 1); /* Minimization problem */

	/* Allocate colcnt-sized arrays */

	obj     = (double *) malloc (num_col * sizeof(double));
	lb      = (double *) malloc (num_col * sizeof(double));
	ub      = (double *) malloc (num_col * sizeof(double));
	ctype   = (char *)   malloc (num_col * sizeof(char));
	colname = (char**) malloc(num_col*sizeof(char*));

	if ( obj     == NULL ||
		lb      == NULL ||
		ub      == NULL ||
		ctype   == NULL ) {
			fprintf (stderr, "Could not allocate colcnt arrays\n");
			status = CPXERR_NO_MEMORY;
			return 0;
	}

	string* strColNames = new string[num_col];
	for ( int i = 0; i < num_x; i++)
	{
			obj[i] = dm.costVec[i].first;                      // fill in obj coefficients
			idx_X.push_back(i);                        // store index for x variables
			strColNames[i] = "x"+ToString(i);
			idx_map_X.insert(make_pair(strColNames[i], i));
			colname[i] = new char[8];
			strcpy(colname[i], strColNames[i].c_str()); // make name for x
			ctype[i] = 'C';
	}
	int idx_cnt = num_x;
	for (int i = 0; i < m-K; i++) // for the first m-K scenarios. each scenario can be assigned onto only one candidate constraint
	{
		for (int j =0; j<n; j++)
		{
			obj[idx_cnt] = 0;
			strColNames[idx_cnt] = "wi" +ToString(i) + "k"+ToString(i) + "j"+ToString(j);
			idx_map_W.insert(make_pair(strColNames[idx_cnt],idx_cnt));
			colname[idx_cnt] = new char[8];
			strcpy(colname[idx_cnt], strColNames[idx_cnt].c_str()); 
			ctype[idx_cnt] = 'C';
			idx_cnt++;
		}
	}
	for (int i = m-K; i < m ; i++)
	{
		for (int k = 0; k < m-K; k++)
		{
			for (int j = 0; j < n; j++)
			{
				obj[idx_cnt] = 0;
				strColNames[idx_cnt] = "wi"+ToString(i) + "k"+ ToString(k) + "j"+ToString(j);
				idx_map_W.insert(make_pair(strColNames[idx_cnt],idx_cnt));
				colname[idx_cnt] = new char[8];
				strcpy(colname[idx_cnt], strColNames[idx_cnt].c_str());
				ctype[i] = 'C';
				idx_cnt++;
			}
		}
	}

	for (int i = 0; i < m-K; i++)
	{
		obj[idx_cnt] = 0;
		idx_Z.push_back(idx_cnt);
		strColNames[idx_cnt] = "zi"+ToString(i) + "k"+ ToString(i) ;
		idx_map_Z.insert(make_pair(strColNames[idx_cnt],idx_cnt));
		colname[idx_cnt] = new char[8];
		strcpy(colname[idx_cnt], strColNames[idx_cnt].c_str());
		if (modeltype == MODEL_TYPE_LP)
			ctype[idx_cnt] = 'C';
		else if( modeltype == MODEL_TYPE_MIP)
			ctype[idx_cnt] = 'B';
		idx_cnt++;
	}
	for (int i = m-K; i < m ; i++)
	{
		for (int k =0; k< m-K; k++)
		{
			obj[idx_cnt] = 0;
			idx_Z.push_back(idx_cnt);
			strColNames[idx_cnt] = "zi"+ToString(i) + "k"+ ToString(k) ;
			idx_map_Z.insert(make_pair(strColNames[idx_cnt],idx_cnt));
			colname[idx_cnt] = new char[8];
			strcpy(colname[idx_cnt], strColNames[idx_cnt].c_str());
			if (modeltype == MODEL_TYPE_LP)
				ctype[idx_cnt] = 'C';
			else if( modeltype == MODEL_TYPE_MIP)
				ctype[idx_cnt] = 'B';
			idx_cnt++;
		}
	}

	delete []strColNames;

	for (int i =0; i < num_col; i ++)       // fill in upper bound
	{
		if (i<num_x + num_w)
			ub[i] = CPX_INFBOUND;
		else
			ub[i] = 1;
	}

	if (modeltype == MODEL_TYPE_LP)
	{
		status = CPXnewcols (cpxenv, cpxlp, num_col, obj, NULL, ub, NULL, colname);
	}
	else if (modeltype == MODEL_TYPE_MIP)
	{
		status = CPXnewcols (cpxenv, cpxlp, num_col, obj, NULL, ub, ctype, colname);
	}
	if ( status ) {
		fprintf (stderr, "Could not add new columns.\n");
		return 0;
	}

	modelType = modeltype;

	/* adding rows to model */
	for (int k = 0; k < m-K; k++) 
	{
		int totaNumCols = n*(1 + K );
		rmatind = (int * )   malloc (totaNumCols * sizeof(int));
		rmatval = (double *) malloc (totaNumCols * sizeof(double));
		int col_idx_cnt = 0;
		double rhs[1];
		char   sense[1];

		for(int i =0; i < m; i++)
		{
			if (i < m-K && i != k)
				continue;

			VEC_MTX_ROW rowEntries = dm.Arows[i];  

			for (int j = 0; j < n; j++)
			{
				PAIR_VALUE_INDEX entryPair = rowEntries[j];
				string varName = "wi"+ToString(i) + "k" + ToString(k) + "j" + ToString(j);
				map<string, int>::iterator itr  = idx_map_W.find(varName);
				rmatind[col_idx_cnt] = itr->second;
				rmatval[col_idx_cnt] = (float)((int)(entryPair.first/dm.rhs[i]*100+0.5))/100.0;
				col_idx_cnt++;
			}
		}
		rhs[0] = 1;
		sense[0] = 'G';
		int    rmatbeg[1];
		rmatbeg[0] = 0;
		char** rowName = new char*; 
		string cstName =  /* i == dm.Arows.size()-1 ?  string("RES")  : */string("Const") +ToString(k);
		char* pStr = new char[20];
		strcpy(pStr, cstName.c_str());
		rowName[0] = pStr;
		status = CPXaddrows (cpxenv, cpxlp, 0, 1, totaNumCols , rhs,
			sense, rmatbeg, rmatind, rmatval,
			NULL, rowName);
		if ( status ) {
			fprintf (stderr, "Could not add new rows.\n");
			return 0;
		}
		delete rowName;
		free_and_null((char **)&rmatind); 
		free_and_null((char **)&rmatval); 

	}

	/* add unitary constraint*/ // every constraint can be chosen from only one scenario
	for (int k = 0; k < m-K; k++)
	{
		int totaNumCols = 1+ K ;
		rmatind = (int * )   malloc (totaNumCols * sizeof(int));
		rmatval = (double *) malloc (totaNumCols * sizeof(double));
		int col_idx_cnt = 0;
		double rhs[1];
		char   sense[1];

		for (int i = 0; i < m ; i ++)
		{
			if ( i < m -K && i != k)
			continue;
			string varName = "zi" + ToString(i) + "k" + ToString(k);
			map<string, int>::iterator itr = idx_map_Z.find(varName);
			rmatind[col_idx_cnt] = itr->second;
			rmatval[col_idx_cnt] = 1;
			col_idx_cnt++;
		}
		rhs[0] = 1;
		sense[0] = 'E';
		int    rmatbeg[1];
		rmatbeg[0] = 0;
		char** rowName = new char*; 
		string cstName =  /* i == dm.Arows.size()-1 ?  string("RES")  : */string("UnitaryK") +ToString(k);
		char* pStr = new char[20];
		strcpy(pStr, cstName.c_str());
		rowName[0] = pStr;
		status = CPXaddrows (cpxenv, cpxlp, 0, 1, totaNumCols , rhs,
			sense, rmatbeg, rmatind, rmatval,
			NULL, rowName);
		if ( status ) {
			fprintf (stderr, "Could not add new rows.\n");
			return 0;

		delete rowName;
		free_and_null((char **)&rmatind); 
		free_and_null((char **)&rmatval); 
		}
	}

	// add unitary constraint for scenario
	for (int i = m - K; i < m ; i++)
	{
		int totaNumCols = m - K ;
		rmatind = (int * )   malloc (totaNumCols * sizeof(int));
		rmatval = (double *) malloc (totaNumCols * sizeof(double));
		int col_idx_cnt = 0;
		double rhs[1];
		char   sense[1];
		for (int k = 0; k < m -K ; k ++)
		{
			string varName = "zi" + ToString(i) + "k" + ToString(k);
			map<string, int>::iterator itr = idx_map_Z.find(varName);
			rmatind[col_idx_cnt] = itr->second;
			rmatval[col_idx_cnt] = 1;
			col_idx_cnt++;		
		}
		rhs[0] = 1;
		sense[0] = 'L';
		int    rmatbeg[1];
		rmatbeg[0] = 0;
		char** rowName = new char*; 
		string cstName =  /* i == dm.Arows.size()-1 ?  string("RES")  : */string("Unitaryi") +ToString(i);
		char* pStr = new char[20];
		strcpy(pStr, cstName.c_str());
		rowName[0] = pStr;
		status = CPXaddrows (cpxenv, cpxlp, 0, 1, totaNumCols , rhs,
			sense, rmatbeg, rmatind, rmatval,
			NULL, rowName);
		if ( status ) {
			fprintf (stderr, "Could not add new rows.\n");
			return 0;
		}
		delete rowName;
		free_and_null((char **)&rmatind); 
		free_and_null((char **)&rmatval); 

	}

	/* adding bound constraint by z*/
	for (int i = 0; i < m; i++)
	{
		for(int k = 0; k < m-K; k++)
		{
			string varName = "zi" + ToString(i) + "k" + ToString(k);
			map<string, int>::iterator itrZ = idx_map_Z.find(varName);
			if (itrZ == idx_map_Z.end())
				continue;
			int ind_z = itrZ->second;
			for (int j = 0 ; j < n ; j ++)
			{
				string varName = "wi"+ToString(i) + "k" + ToString(k) + "j" + ToString(j);
				map<string, int>::iterator itrW  = idx_map_W.find(varName);
				if (itrW == idx_map_W.end())
					continue;
				int ind_w = itrW->second;


				int totaNumCols = 2 ;
				rmatind = (int * )   malloc (totaNumCols * sizeof(int));
				rmatval = (double *) malloc (totaNumCols * sizeof(double));
				double rhs[1];
				char   sense[1];
				rmatind[0] = ind_w;
				rmatval[0] = 1;
				rmatind[1] = ind_z;
				rmatval[1] = -1;
				rhs[0] = 0;
				sense[0] = 'L';
				int    rmatbeg[1];
				rmatbeg[0] = 0;
				char** rowName = new char*; 
				string cstName =  string("BdZ") +"w" + "i"+ToString(i) + "k" + ToString(k) + "j" + ToString(j);
				char* pStr = new char[20];
				strcpy(pStr, cstName.c_str());
				rowName[0] = pStr;
				status = CPXaddrows (cpxenv, cpxlp, 0, 1, totaNumCols , rhs,
					sense, rmatbeg, rmatind, rmatval,
					NULL, rowName);
				if ( status ) {
					fprintf (stderr, "Could not add new rows.\n");
					return 0;
				}
				delete rowName;
				free_and_null((char **)&rmatind); 
				free_and_null((char **)&rmatval); 
			}
		}
	}

	/* adding bound constraint by x*/
	for (int i = 0; i < m; i++)
	{
		for(int k = 0; k < m-K; k++)
		{
			for (int j = 0 ; j < n ; j ++)
			{
				string varName = "wi"+ToString(i) + "k" + ToString(k) + "j" + ToString(j);
				map<string, int>::iterator itrW  = idx_map_W.find(varName);
				if (itrW == idx_map_W.end())
					continue;
				int ind_w = itrW->second;
				string varName1 = "x" + ToString(j) ;
				map<string, int>::iterator itrX  = idx_map_X.find(varName1);
				int ind_x = itrX->second;

				int totaNumCols = 2 ;
				rmatind = (int * )   malloc (totaNumCols * sizeof(int));
				rmatval = (double *) malloc (totaNumCols * sizeof(double));
				double rhs[1];
				char   sense[1];
				rmatind[0] = ind_w;
				rmatval[0] = 1;
				rmatind[1] = ind_x;
				rmatval[1] = -1;
				rhs[0] = 0;
				sense[0] = 'L';
				int    rmatbeg[1];
				rmatbeg[0] = 0;
				char** rowName = new char*; 
				string cstName =  string("BdX") +"w" + "i"+ToString(i) + "k" + ToString(k) + "j" + ToString(j);
				char* pStr = new char[20];
				strcpy(pStr, cstName.c_str());
				rowName[0] = pStr;
				status = CPXaddrows (cpxenv, cpxlp, 0, 1, totaNumCols , rhs,
					sense, rmatbeg, rmatind, rmatval,
					NULL, rowName);
				if ( status ) {
					fprintf (stderr, "Could not add new rows.\n");
					return 0;
				}
				delete rowName;
				free_and_null((char **)&rmatind); 
				free_and_null((char **)&rmatval); 
			}
		}
	}

	/* adding bound constraint by x and z*/
	for (int i = 0; i < m; i++)
	{
		for(int k = 0; k < m-K; k++)
		{
			for (int j = 0 ; j < n ; j ++)
			{
				string varNameW = "wi"+ToString(i) + "k" + ToString(k) + "j" + ToString(j);
				map<string, int>::iterator itrW  = idx_map_W.find(varNameW);
				if (itrW == idx_map_W.end())
					continue;
				int ind_w = itrW->second;
				string varNameX = "x" + ToString(j) ;
				map<string, int>::iterator itrX  = idx_map_X.find(varNameX);
				int ind_x = itrX->second;
				string varNameZ = "zi" + ToString(i) + "k" + ToString(k);
				map<string, int>::iterator itrZ = idx_map_Z.find(varNameZ);
				int ind_z = itrZ->second;

				int totaNumCols = 3 ;
				rmatind = (int * )   malloc (totaNumCols * sizeof(int));
				rmatval = (double *) malloc (totaNumCols * sizeof(double));
				double rhs[1];
				char   sense[1];
				rmatind[0] = ind_w;
				rmatval[0] = -1;
				rmatind[1] = ind_x;
				rmatval[1] = 1;
				rmatind[2] = ind_z;
				rmatval[2] = 1;
				rhs[0] = 1;
				sense[0] = 'L';
				int    rmatbeg[1];
				rmatbeg[0] = 0;
				char** rowName = new char*; 
				string cstName =  string("BdXZ") +"w" + "i"+ToString(i) + "k" + ToString(k) + "j" + ToString(j);
				char* pStr = new char[20];
				strcpy(pStr, cstName.c_str());
				rowName[0] = pStr;
				status = CPXaddrows (cpxenv, cpxlp, 0, 1, totaNumCols , rhs,
					sense, rmatbeg, rmatind, rmatval,
					NULL, rowName);
				if ( status ) {
					fprintf (stderr, "Could not add new rows.\n");
					return 0;
				}
				delete rowName;
				free_and_null((char **)&rmatind); 
				free_and_null((char **)&rmatval); 
			}
		}
	}

	string fileName ;
	if (modeltype == MODEL_TYPE_LP)
		fileName = modelName + "R.lp" ;
	else if (modeltype == MODEL_TYPE_MIP)
		fileName = modelName + "MIP.lp";
	writeModel2File(fileName.c_str());


	//CPXwriteprob(cpxenv, lp, "initialLPRelaxation.lp", NULL);

	free_and_null((char **)&obj);
	free_and_null((char **)&lb);
	free_and_null((char **)&ub);
	free_and_null((char **)&ctype);
	free_and_null((char **)&colname); 
	free_and_null((char **)&rmatind); 
	free_and_null((char **)&rmatval); 

	return 0;
}





int OptModel::populateModelNetworkFow( DataModel& dm, OptModel::MODEL_TYPE modeltype, char cardinalitySign)
{
	setDataModel(&dm);
	idx_X.clear();
	idx_Z.clear();
	idx_S.clear();
	idx_T.clear();

	double *obj     = NULL;
	double *lb      = NULL;
	double *ub      = NULL;
	char   *ctype   = NULL;
	int    *rmatind = NULL;
	double *rmatval = NULL;
	char** colname = NULL;

	int    status = 0;
	int num_cols = DataModel::n + 3*DataModel::m;

	CPXchgobjsen (cpxenv, cpxlp, 1); /* Minimization problem */

	/* Allocate colcnt-sized arrays */

	obj     = (double *) malloc (num_cols * sizeof(double));
	lb      = (double *) malloc (num_cols* sizeof(double));
	ub      = (double *) malloc (num_cols * sizeof(double));
	ctype   = (char *)   malloc (num_cols * sizeof(char));
	colname = (char**) malloc(num_cols*sizeof(char*));

	if ( obj     == NULL ||
		lb      == NULL ||
		ub      == NULL ||
		ctype   == NULL ) {
			fprintf (stderr, "Could not allocate colcnt arrays\n");
			status = CPXERR_NO_MEMORY;
			return 0;
	}

	string* strColNames = new string[num_cols];
	for (int i =0; i < num_cols; i ++)
	{
		if (i < DataModel::n )
		{
			obj[i] = dm.costVec[i].first;                      // fill in obj coefficients
			idx_X.push_back(i);                        // store index for x variables
			strColNames[i] = "x"+ToString(i);
			colname[i] = new char[8];
			strcpy(colname[i], strColNames[i].c_str()); // make name for x
			ctype[i] = 'C';
		}
		else if ( i >= DataModel::n && i < DataModel::n + DataModel::m)
		{
			obj[i] = 0;
			idx_S.push_back(i);
			strColNames[i] = "s" + ToString(i-DataModel::n);
			colname[i] = new char[8];
			strcpy(colname[i], strColNames[i].c_str()); // make name for z
			ctype[i] = 'C';
		}
		else if (i >= DataModel::n + DataModel::m && i < DataModel::n + 2*DataModel::m)
		{
			obj[i] = 0;
			idx_T.push_back(i);
			strColNames[i] = "t" + ToString(i-DataModel::n - DataModel::m);
			colname[i] = new char[8];
			strcpy(colname[i], strColNames[i].c_str()); // make name for z
			ctype[i] = 'C';
		}
		else
		{
			obj[i] = 0;
			idx_Z.push_back(i);
			strColNames[i] = "z" + ToString(i-(DataModel::n + 2*DataModel::m));
			colname[i] = new char[8];
			strcpy(colname[i], strColNames[i].c_str()); // make name for z
			if (modeltype == MODEL_TYPE_LP)
				ctype[i] = 'C';
			else if( modeltype == MODEL_TYPE_MIP)
				ctype[i] = 'B';
		}
	}
	delete []strColNames;

	for (int i =0; i < num_cols; i ++)       // fill in upper bound
	{
		if (i<DataModel::n + 2*DataModel::m)
			ub[i] = CPX_INFBOUND;
		else
			ub[i] = 1;
	}

	if (modeltype == MODEL_TYPE_LP)
	{
		status = CPXnewcols (cpxenv, cpxlp, num_cols, obj, NULL, ub, NULL, colname);
	}
	else if (modeltype == MODEL_TYPE_MIP)
	{
		status = CPXnewcols (cpxenv, cpxlp, num_cols, obj, NULL, ub, ctype, colname);
	}
	if ( status ) {
		fprintf (stderr, "Could not add new columns.\n");
		return 0;
	}

	modelType = modeltype;

	/* adding rows to model */
	for (int i = 0; i < DataModel::m; i++) 
	{

		int totaNumCols = DataModel::n + 2;
		rmatind = (int * )   malloc (totaNumCols * sizeof(int));
		rmatval = (double *) malloc (totaNumCols * sizeof(double));

		VEC_MTX_ROW rowEntries = dm.Arows[i];  // fill in the entries for a row
		double rhs[1];
		char   sense[1];

		for (int j =0; j <DataModel::n; j ++)
		{
			PAIR_VALUE_INDEX entryPair = rowEntries[j];
			rmatind[j] = entryPair.second;
			rmatval[j] = entryPair.first;
		}
		rmatind[DataModel::n] = DataModel::n + i;
		rmatval[DataModel::n] = 1;
		rmatind[DataModel::n + 1] = DataModel::n + DataModel::m + i;
		rmatval[DataModel::n + 1] = -1;
		rhs[0] = dm.rhs[i];
		sense[0] = 'E';

		int    rmatbeg[1];
		rmatbeg[0] = 0;
		char** rowName = new char*; 
		string cstName = string("SN") +ToString(i);
		char* pStr = new char[20];
		strcpy(pStr, cstName.c_str());
		rowName[0] = pStr;
		status = CPXaddrows (cpxenv, cpxlp, 0, 1, totaNumCols , rhs,
			sense, rmatbeg, rmatind, rmatval,
			NULL, rowName);
		if ( status ) {
			fprintf (stderr, "Could not add new rows.\n");
			return 0;
		}
		delete rowName;
		free_and_null((char **)&rmatind); 
		free_and_null((char **)&rmatval); 

	}

	// add fixed charge constraints

	for (int i = 0; i < DataModel::m; i++) 
	{

		int totaNumCols = 2;
		rmatind = (int * )   malloc (totaNumCols * sizeof(int));
		rmatval = (double *) malloc (totaNumCols * sizeof(double));

		double rhs[1];
		char   sense[1];

		rmatind[0] = DataModel::n + i;
		rmatval[0] = 1;
		rmatind[1] = DataModel::n + 2*DataModel::m + i;
		rmatval[1] = -dm.rhs[i];
		rhs[0] = 0;
		sense[0] = 'L';

		int    rmatbeg[1];
		rmatbeg[0] = 0;
		char** rowName = new char*; 
		string cstName = string("UB") +ToString(i);
		char* pStr = new char[20];
		strcpy(pStr, cstName.c_str());
		rowName[0] = pStr;
		status = CPXaddrows (cpxenv, cpxlp, 0, 1, totaNumCols , rhs,
			sense, rmatbeg, rmatind, rmatval,
			NULL, rowName);
		if ( status ) {
			fprintf (stderr, "Could not add new rows.\n");
			return 0;
		}
		delete rowName;
		free_and_null((char **)&rmatind); 
		free_and_null((char **)&rmatval); 
	}


	/* adding cardinality constraint */
	int totaNumCols = DataModel::m ;
	rmatind = (int * )   malloc (totaNumCols * sizeof(int));
	rmatval = (double *) malloc (totaNumCols * sizeof(double));

	for (int j =0; j <DataModel::m; j ++)
	{
		rmatind[j] = DataModel::n + 2*DataModel::m + j;
		rmatval[j] = 1;
	}
	int    rmatbeg[1];
	double rhs[1];
	char   sense[1];
	rhs[0] = DataModel::k;
	sense[0] = cardinalitySign;
	rmatbeg[0] = 0;
	char ** rowName = new char*;
	rowName[0] = "CARD";
	status = CPXaddrows (cpxenv, cpxlp, 0, 1, DataModel::m, rhs,
		sense, rmatbeg, rmatind, rmatval,
		NULL, rowName);
	delete rowName;
	if ( status ) {
		fprintf (stderr, "Could not add new rows.\n");
		return 0;
	}

	string fileName ;
	if (modeltype == MODEL_TYPE_LP)
		fileName = modelName + "R.lp" ;
	else if (modeltype == MODEL_TYPE_MIP)
		fileName = modelName + "MIP.lp";
	writeModel2File(fileName.c_str());
	//CPXwriteprob(cpxenv, lp, "initialLPRelaxation.lp", NULL);

	free_and_null((char **)&obj);
	free_and_null((char **)&lb);
	free_and_null((char **)&ub);
	free_and_null((char **)&ctype);
	free_and_null((char **)&colname); 
	free_and_null((char **)&rmatind); 
	free_and_null((char **)&rmatval); 

	return 1;

}


int OptModel::initilizeMasterModel( DataModel& dm, OptModel::MODEL_TYPE modeltype, char cardinalitySign)
{
	idx_X.clear();
	idx_Z.clear();
	idx_S.clear();
	double *obj     = NULL;
	double *lb      = NULL;
	double *ub      = NULL;
	char   *ctype   = NULL;
	int    *rmatind = NULL;
	double *rmatval = NULL;
	char** colname = NULL;

	int    status = 0;

	CPXchgobjsen (cpxenv, cpxlp, 1); /* Minimization problem */

	/* Allocate colcnt-sized arrays */

	obj     = (double *) malloc (DataModel::NUM_COLS * sizeof(double));
	lb      = (double *) malloc (DataModel::NUM_COLS * sizeof(double));
	ub      = (double *) malloc (DataModel::NUM_COLS * sizeof(double));
	ctype   = (char *)   malloc (DataModel::NUM_COLS * sizeof(char));
	colname = (char**) malloc(DataModel::NUM_COLS*sizeof(char*));

	if ( obj     == NULL ||
		lb      == NULL ||
		ub      == NULL ||
		ctype   == NULL ) {
			fprintf (stderr, "Could not allocate colcnt arrays\n");
			status = CPXERR_NO_MEMORY;
			return 0;
	}

	string* strColNames = new string[DataModel::NUM_COLS];
	for (int i =0; i < DataModel::NUM_COLS; i ++)
	{
		if (i < DataModel::n )
		{
			obj[i] = dm.costVec[i].first;                      // fill in obj coefficients
			idx_X.push_back(i);                        // store index for x variables
			strColNames[i] = "x"+ToString(i);
			colname[i] = new char[8];
			strcpy(colname[i], strColNames[i].c_str()); // make name for x
			ctype[i] = 'C';
		}
		else{
			obj[i] = 0;
			idx_Z.push_back(i);
			strColNames[i] = "z" + ToString(i-DataModel::n);
			colname[i] = new char[8];
			strcpy(colname[i], strColNames[i].c_str()); // make name for z
			if (modeltype == MODEL_TYPE_LP)
				ctype[i] = 'C';
			else if( modeltype == MODEL_TYPE_MIP)
				ctype[i] = 'I';
		}
	}
	delete []strColNames;

	for (int i =0; i < DataModel::NUM_COLS; i ++)       // fill in upper bound
	{
		if (i<DataModel::n)
			ub[i] = CPX_INFBOUND;
		else
			ub[i] = 1 ; //CPX_INFBOUND;  // master model is MIP 
	}

	if (modeltype == MODEL_TYPE_LP)
	{
		status = CPXnewcols (cpxenv, cpxlp, DataModel::NUM_COLS, obj, NULL, ub, NULL, colname);
	}
	else if (modeltype == MODEL_TYPE_MIP)
	{
		status = CPXnewcols (cpxenv, cpxlp, DataModel::NUM_COLS, obj, NULL, ub, ctype, colname);
	}
	if ( status ) {
		fprintf (stderr, "Could not add new columns.\n");
		return 0;
	}

	modelType = modeltype;


	/* adding cardinality constraint */

	rmatind = (int * )   malloc (DataModel::m * sizeof(int));
	rmatval = (double *) malloc (DataModel::m  * sizeof(double));

	for (int j =0; j <DataModel::m; j ++)
	{
		rmatind[j] = idx_Z[j];
		rmatval[j] = 1;
	}
	int    rmatbeg[1];
	double rhs[1];
	char   sense[1];
	rhs[0] = DataModel::k;
	sense[0] = 'L';
	rmatbeg[0] = 0;
	status = CPXaddrows (cpxenv, cpxlp, 0, 1, DataModel::m , rhs,
		sense, rmatbeg, rmatind, rmatval,
		NULL, NULL);
	if ( status ) {
		fprintf (stderr, "Could not add new rows.\n");
		return 0;
	}

	free_and_null((char **)&rmatind); 
	free_and_null((char **)&rmatval); 


	string fileName ;
	if (modeltype == MODEL_TYPE_LP)
		fileName = modelName + "R.lp" ;
	else if (modeltype == MODEL_TYPE_MIP)
		fileName = modelName + "MIP.lp";
	writeModel2File(fileName.c_str());
	//CPXwriteprob(cpxenv, lp, "initialLPRelaxation.lp", NULL);

	free_and_null((char **)&obj);
	free_and_null((char **)&lb);
	free_and_null((char **)&ub);
	free_and_null((char **)&ctype);
	free_and_null((char **)&colname); 
	free_and_null((char **)&rmatind); 
	free_and_null((char **)&rmatval); 

	return 1;

}

int OptModel::initilizeMasterModelWithResConst( DataModel& dm, OptModel::MODEL_TYPE modeltype, char cardinalitySign, bool negativeBound)
{
	dm_p = &dm;
	idx_X.clear();
	idx_Z.clear();
	idx_S.clear();
	double *obj     = NULL;
	double *lb      = NULL;
	double *ub      = NULL;
	char   *ctype   = NULL;
	int    *rmatind = NULL;
	double *rmatval = NULL;
	char** colname = NULL;

	int    status = 0;

	CPXchgobjsen (cpxenv, cpxlp, 1); /* Minimization problem */

	/* Allocate colcnt-sized arrays */

	obj     = (double *) malloc (DataModel::NUM_COLS * sizeof(double));
	lb      = (double *) malloc (DataModel::NUM_COLS * sizeof(double));
	ub      = (double *) malloc (DataModel::NUM_COLS * sizeof(double));
	ctype   = (char *)   malloc (DataModel::NUM_COLS * sizeof(char));
	colname = (char**) malloc(DataModel::NUM_COLS*sizeof(char*));

	if ( obj     == NULL ||
		lb      == NULL ||
		ub      == NULL ||
		ctype   == NULL ) {
			fprintf (stderr, "Could not allocate colcnt arrays\n");
			status = CPXERR_NO_MEMORY;
			return 0;
	}

	string* strColNames = new string[DataModel::NUM_COLS];
	for (int i =0; i < DataModel::NUM_COLS; i ++)
	{
		if (i < DataModel::n )
		{
			obj[i] = dm.costVec[i].first;                      // fill in obj coefficients
			idx_X.push_back(i);                        // store index for x variables
			strColNames[i] = "x"+ToString(i);
			colname[i] = new char[8];
			strcpy(colname[i], strColNames[i].c_str()); // make name for x
			ctype[i] = 'C';
		}
		else{
			obj[i] = 0;
			idx_Z.push_back(i);
			strColNames[i] = "z" + ToString(i-DataModel::n);
			colname[i] = new char[8];
			strcpy(colname[i], strColNames[i].c_str()); // make name for z
			if (modeltype == MODEL_TYPE_LP)
				ctype[i] = 'C';
			else if( modeltype == MODEL_TYPE_MIP)
				ctype[i] = 'I';
		}
	}
	delete []strColNames;

	for (int i =0; i < DataModel::NUM_COLS; i ++)       // fill in upper bound
	{
		if (i<DataModel::n)
			ub[i] = CPX_INFBOUND;
		else
			ub[i] = 1 ; //CPX_INFBOUND;  // master model is MIP 
	}

	for (int i =0; i < DataModel::NUM_COLS; i ++)       // fill in upper bound
	{
		if (i<DataModel::n)
			lb[i] = 0;
		else
			lb[i] = - CPX_INFBOUND;
	}

	if (modeltype == MODEL_TYPE_LP)
	{
		if (negativeBound)
			status = CPXnewcols (cpxenv, cpxlp, DataModel::NUM_COLS, obj, lb, ub, NULL, colname);
		else 
			status = CPXnewcols (cpxenv, cpxlp, DataModel::NUM_COLS, obj, NULL, ub, NULL, colname);
	}
	else if (modeltype == MODEL_TYPE_MIP)
	{
		status = CPXnewcols (cpxenv, cpxlp, DataModel::NUM_COLS, obj, NULL, ub, ctype, colname);
	}
	if ( status ) {
		fprintf (stderr, "Could not add new columns.\n");
		return 0;
	}

	modelType = modeltype;


	/* adding cardinality constraint */

	rmatind = (int * )   malloc (DataModel::m * sizeof(int));
	rmatval = (double *) malloc (DataModel::m  * sizeof(double));

	for (int j =0; j <DataModel::m; j ++)
	{
		rmatind[j] = idx_Z[j];
		rmatval[j] = 1;
	}
	int    rmatbeg[1];
	double rhs[1];
	char   sense[1];
	rhs[0] = DataModel::k;
	sense[0] = 'L';
	rmatbeg[0] = 0;
	char*  rowName[1];
	rowName[0] = "CARD";
	status = CPXaddrows (cpxenv, cpxlp, 0, 1, DataModel::m , rhs,
		sense, rmatbeg, rmatind, rmatval,
		NULL, rowName);
	if ( status ) {
		fprintf (stderr, "Could not add new rows.\n");
		return 0;
	}
	free_and_null((char **)&rmatind); 
	free_and_null((char **)&rmatval); 


	/* adding resource constraint */

	rmatind = (int * )   malloc (DataModel::n * sizeof(int));
	rmatval = (double *) malloc (DataModel::n  * sizeof(double));

	for (int j =0; j <DataModel::n; j ++)
	{
		rmatind[j] = idx_X[j];
		rmatval[j] = 1;
	}
	rhs[0] = 1;
	sense[0] = DataModel::res_const_sign == 0 ? 'E' : 'L';
	rmatbeg[0] = 0;
	rowName[0] = "RES";
	status = CPXaddrows (cpxenv, cpxlp, 0, 1, DataModel::n , rhs,
		sense, rmatbeg, rmatind, rmatval,
		NULL, rowName);
	if ( status ) {
		fprintf (stderr, "Could not add new rows.\n");
		return 0;
	}

	free_and_null((char **)&rmatind); 
	free_and_null((char **)&rmatval); 


	string fileName ;
	if (modeltype == MODEL_TYPE_LP)
		fileName = modelName + "R.lp" ;
	else if (modeltype == MODEL_TYPE_MIP)
		fileName = modelName + "MIP.lp";
	writeModel2File(fileName.c_str());
	//CPXwriteprob(cpxenv, lp, "initialLPRelaxation.lp", NULL);

	free_and_null((char **)&obj);
	free_and_null((char **)&lb);
	free_and_null((char **)&ub);
	free_and_null((char **)&ctype);
	free_and_null((char **)&colname); 
	free_and_null((char **)&rmatind); 
	free_and_null((char **)&rmatval); 

	return 1;

}



void OptModel::writeModel2File()
{

	string fileName = modelName + ".lp";
	CPXwriteprob(cpxenv, cpxlp, fileName.c_str(), NULL);

}

void OptModel::writeModel2File(string fileName)
{

	CPXwriteprob(cpxenv, cpxlp, fileName.c_str(), NULL);

}
CPXENVptr& OptModel::getCPXEnv(void)
{
	return cpxenv;
}

CPXLPptr& OptModel::getCPXLP(void)
{
	return cpxlp;
}

vector<int>& OptModel::getIndex_X(void)
{
	return idx_X;
}

vector<int>& OptModel::getIndex_Z(void)
{
	return idx_Z;
}


vector<int>& OptModel::getIndex_S(void)
{
	return idx_S;
}

int OptModel::solve()
{

	int probType =CPXgetprobtype(cpxenv, cpxlp);
	int status=0;
	int solve_status = 0;
	double timestamp;
	CPXgettime(cpxenv, &timestamp);
	double starttime =timestamp; 
	if (probType == CPXPROB_MILP)
	{
//		setBranchCallBackFunction(myBranchCallbackRecordCPXRootLPOptVal, this);
		solve_status = CPXmipopt (cpxenv, cpxlp);
		CPXgettime(cpxenv, &timestamp);
		MIP_solving_time = timestamp - starttime;
		Prev_MIP_optimal_value = MIP_optimal_value;
		CPXgetbestobjval(cpxenv, cpxlp, &MIP_optimal_value);
		CPXgetsolnpoolobjval(cpxenv, cpxlp, -1, &MIP_INCUMBENT_OBJ_VAL);
		CPXgetmiprelgap(cpxenv, cpxlp, &MIP_GAP);
		MIP_BB_nodes = CPXgetnodecnt(cpxenv, cpxlp) ;
	}
	else if (probType == CPXPROB_LP)
	{
		status = CPXlpopt(cpxenv, cpxlp);
		CPXgettime(cpxenv, &timestamp);
		LP_solving_time = timestamp - starttime;
//		solve_status = CPXgetstat(cpxenv, cpxlp); 
		Prev_LP_optimal_value = LP_optimal_value;
		status = CPXgetobjval(cpxenv, cpxlp, &LP_optimal_value);
	}
	return solve_status;


}

int OptModel::setCPXParameter(int whichparam, int newvalue)
{
	return CPXsetintparam(cpxenv, whichparam, newvalue); 
}

vector<PAIR_INDEX_VALUE> OptModel::getSolutionValues(int varType)
{
	int num_cols  = CPXgetnumcols(cpxenv,  cpxlp) ; 
	double *val_vars = new double[num_cols];
	CPXgetx(cpxenv, cpxlp, val_vars, 0, num_cols-1);

	vector<PAIR_INDEX_VALUE> sol;
	if (varType == 1)
	{
		for (int j = 0; j < idx_X.size(); j++)
		{
			sol.push_back(PAIR_INDEX_VALUE(idx_X[j], val_vars[idx_X[j]]));
		}
	}
	else if (varType ==2)
	{
		for (int j = 0; j < idx_Z.size(); j++)
		{
			sol.push_back(PAIR_INDEX_VALUE(idx_Z[j], val_vars[idx_Z[j]]));
		}
	}
	else if (varType == 3)
	{

		for (int j = 0; j < idx_S.size(); j++)
		{
			sol.push_back(PAIR_INDEX_VALUE(idx_S[j], val_vars[idx_S[j]]));
		}
	}
	delete []val_vars;
	return sol;
}

int OptModel::mapVarIndexName(void) // variable subscript -- column index
{
	idx_Z.clear();
	idx_X.clear();
	idx_S.clear();
	int num_rows = CPXgetnumrows(cpxenv,  cpxlp);
	int num_cols  = CPXgetnumcols(cpxenv,  cpxlp) ; // need to know why
	int surplus;

	int status = CPXgetcolname (cpxenv, cpxlp, NULL, NULL, 0, &surplus, 0,num_cols-1);
	int cur_colnamespace = - surplus;
	if ( cur_colnamespace > 0 ) {
		char          **cur_colname      = (char **) malloc (sizeof(char *)*(num_cols));
		char          *cur_colnamestore = (char *)  malloc (cur_colnamespace);
		if ( cur_colname      == NULL ||
			cur_colnamestore == NULL   ) {
				fprintf (stderr, "Failed to get memory for column names.\n");
				status = -1;
				return 0;
		}
		status = CPXgetcolname (cpxenv, cpxlp, cur_colname, cur_colnamestore, 
			cur_colnamespace, &surplus, 0, num_cols-1);
		// map name with column index
		for (int i = 0; i <num_cols ; i ++)
		{
			if (cur_colname[i][0] == 'x')
			{
				idx_X.push_back(i);
			}
			else if (cur_colname[i][0] == 'z')
			{
				idx_Z.push_back(i);
			}
			else if (cur_colname[i][0] == 's')
			{
				idx_S.push_back(i);
			}
			string varName = cur_colname[i];
			var_name_idx_table.insert(make_pair(varName, i));
			idx_var_name_table.insert(make_pair(i, varName));
		}
		free_and_null( (char **)&cur_colname);
		free_and_null((char **)&cur_colnamestore);
	}
	DataModel::n = idx_X.size();
	DataModel::m = idx_Z.size();
	return 1;
}

void OptModel::printSolution(int varType)
{
	int num_cols  = CPXgetnumcols(cpxenv,  cpxlp) ; 
	double *val_vars = new double[num_cols];
	CPXgetx(cpxenv, cpxlp, val_vars, 0, num_cols-1);
	cout<<endl;
	if (varType == 0)
	{
		for (int i = 0; i < idx_X.size(); i++)
		{
			cout << "x"<<i << " = " << val_vars[idx_X[i]]<<endl;
		}
		for (int i = 0; i < idx_Z.size(); i++)
		{
			cout << "z"<<i << " = " << val_vars[idx_Z[i]]<<endl;
		}
	}
	else if (varType == 1)
	{
		for (int i = 0; i < idx_X.size(); i++)
		{
			cout << "x"<<i << " = " << val_vars[idx_X[i]]<<endl;
		}

	}
	else if (varType == 2)
	{
		for (int i = 0; i < idx_Z.size(); i++)
		{
			cout << "z"<<i << " = " << val_vars[idx_Z[i]]<<endl;
		}

	}
	delete[] val_vars;
}


void OptModel::printSolution2File(string fileName, int varType, CPXCENVptr env, CPXLPptr  lp)
{
	if (env == NULL)
		env = cpxenv;
	if (lp == NULL)
		lp = cpxlp;

	ofstream log(fileName.c_str(),ios_base::out|ios_base::app);

	vector < pair<string, int> > cNames;
	int num_cols = CPXgetnumcols(cpxenv,  cpxlp);
	int surplus;

	int status = CPXgetcolname (cpxenv, cpxlp, NULL, NULL, 0, &surplus, 0,num_cols-1);
	int cur_colnamespace = - surplus;
	if ( cur_colnamespace > 0 ) 
	{
		char          **cur_colname = (char **) malloc (sizeof(char *)*(num_cols));
		char          *cur_colnamestore = (char *)  malloc (cur_colnamespace);
		if ( cur_colname      == NULL ||
			cur_colnamestore == NULL   ) {
				fprintf (stderr, "Failed to get memory for column names.\n");
				status = -1;
				cout<< endl <<" error in printTableau"<<endl;
				return ;
		}
		status = CPXgetcolname (cpxenv, cpxlp, cur_colname, cur_colnamestore, 
			cur_colnamespace, &surplus, 0, num_cols-1);
		// map name with column index
		for (int c = 0; c <num_cols ; c ++)
		{
			string rowName = cur_colname[c];
			cNames.push_back(make_pair(rowName,c));
		}
		free_and_null( (char **)&cur_colname);
		free_and_null((char **)&cur_colnamestore);
	}


	double *val_vars = new double[num_cols];
	CPXgetx(env, lp, val_vars, 0, num_cols-1);
	cout<<endl;
	if (varType == 0)
	{
		for (int i = 0; i < num_cols; i++)
		{
			string varName = cNames[i].first;
			log << varName<< " = "<< val_vars[i]<< endl;
		}
	}
	else if (varType == 1)
	{
		for (int i = 0; i < num_cols; i++)
		{
			string varName = cNames[i].first;
			if (varName.substr(0,1) != "x")
				continue;
			log << varName<< " = " << val_vars[i]<< endl;
		}
	}
	else if (varType == 2)
	{
		for (int i = 0; i < num_cols; i++)
		{
			string varName = cNames[i].first;
			if (varName.substr(0,1) != "z")
				continue;
			log << varName<< " = " << val_vars[i]<< endl;
		}
	}
	delete[] val_vars;
	log.close();
}



void OptModel::printObjectiveValue(void)
{
	double objval_p;
	int status = CPXgetobjval(cpxenv, cpxlp, &objval_p);
	cout << endl<<"optimal objective value: " << objval_p <<endl;

}

int OptModel::changeLP2MIP(int changeType)
{
	int status = 0;
	if (changeType == 0 )
	{
		int cnt = idx_Z.size();
		int *indices = new int[idx_Z.size()];
		for (int i = 0; i < idx_Z.size(); i++)
		{ 
			indices[i] = idx_Z[i];
		}
		char* xctype= new char[cnt];
		for (int i = 0; i < idx_Z.size(); i ++)
		{
			xctype[i] = CPX_INTEGER;
		}
		status = CPXchgctype( cpxenv,  cpxlp,  cnt,   indices,   xctype);
		if ( status ) {
			fprintf (stderr, "Failed to change variable type.\n");
		}
		delete[] indices;
		modelType= MODEL_TYPE_MIP;
	}
	else if (changeType == 1)
	{
		int cnt = idx_Z.size() + idx_X.size();
		int *indices = new int[cnt];
		for (int i = 0; i < idx_X.size(); i++)
		{ 
			indices[i] = idx_X[i];
		}
		for (int i = 0; i < idx_Z.size(); i++)
		{ 
			indices[idx_X.size() +i] =   idx_Z[i];
		}
		char* xctype= new char[cnt];
		for (int i = 0; i <cnt; i ++)
		{
			xctype[i] = CPX_INTEGER;
		}
		status = CPXchgctype( cpxenv,  cpxlp,  cnt,   indices,   xctype);
		if ( status ) {
			fprintf (stderr, "Failed to change variable type.\n");
		}
		delete[] indices;
		modelType= MODEL_TYPE_MIP;

	}

	return status;
	
}

OptModel::~OptModel(){
	/* Free up the problem as allocated by CPXcreateprob, if necessary */
	int status;
	if ( cpxlp != NULL ) {
		status = CPXfreeprob (cpxenv, &cpxlp);
		if ( status ) {
			fprintf (stderr, "CPXfreeprob failed, error code %d.\n", status);
		}
	}

	/* Free up the CPLEX environment, if necessary */

	if ( cpxenv != NULL  && is_a_copy == false ) {
		status = CPXcloseCPLEX (&cpxenv);

		/* Note that CPXcloseCPLEX produces no output,
		so the only way to see the cause of the error is to use
		CPXgeterrorstring.  For other CPLEX routines, the errors will
		be seen if the CPX_PARAM_SCRIND indicator is set to CPX_ON. */

		if ( status > 0 ) {
			char  errmsg[1024];
			fprintf (stderr, "Could not close CPLEX environment.\n");
// 			CPXgeterrorstring (cpxenv, status, errmsg);
// 			fprintf (stderr, "%s", errmsg);
		}
	}

}


int OptModel::printTableau(string tableauFileName, CPXCENVptr env, CPXLPptr  lp)
{
		if (env == NULL)
			env = cpxenv;
		if (lp == NULL)
			lp = cpxlp;

		int num_rows = CPXgetnumrows(env,  lp);
		int num_cols  = CPXgetnumcols(env,  lp) ; 
		int surplus;
		vector<string> varNames;

		int status = CPXgetcolname (env, lp, NULL, NULL, 0, &surplus, 0,num_cols-1);
		int cur_colnamespace = - surplus;
		if ( cur_colnamespace > 0 ) {
			char          **cur_colname      = (char **) malloc (sizeof(char *)*(num_cols));
			char          *cur_colnamestore = (char *)  malloc (cur_colnamespace);
			if ( cur_colname      == NULL ||
				cur_colnamestore == NULL   ) {
					fprintf (stderr, "Failed to get memory for column names.\n");
					status = -1;
					cout<< endl <<" error in printTableau"<<endl;
					return 0;
			}
			status = CPXgetcolname (env, lp, cur_colname, cur_colnamestore, 
				cur_colnamespace, &surplus, 0, num_cols-1);
			// map name with column index
			for (int i = 0; i <num_cols ; i ++)
			{
				string varName = cur_colname[i];
				varNames.push_back(varName);
			}
			free_and_null( (char **)&cur_colname);
			free_and_null((char **)&cur_colnamestore);

		}

		/***************get head ****************************/
		int* head = new int[num_rows]; // head is the 0-th column of the tableau
		double *head_val = new double[num_rows];
		status = CPXgetbhead (env, lp, head, head_val);

		ofstream log(tableauFileName.c_str() ,ios_base::app);
		streambuf * old_buf = std::cout.rdbuf(log.rdbuf()); 

		double *objval_p = new double;
		status = CPXgetobjval(env, lp,  objval_p);
		cout<<"objective value, " << *objval_p<<endl;

		double* tabealu;
		tabealu = (double*)malloc(num_rows*(num_cols+num_rows)*sizeof(double));
		for (int i = 0; i < num_cols; i++)
		{
			CPXbinvacol(env, lp, i, tabealu+i*num_rows);
		}
		char* sense = new char[num_rows];
		CPXgetsense(env, lp, sense, 0, num_rows -1);

		for (int i = num_cols; i < num_cols + num_rows; i++)
		{
			CPXbinvcol(env, lp, i - num_cols, tabealu + i*num_rows);
			if (sense[i - num_cols] == 'G')
			{
				for (int j = 0; j < num_rows; j ++)
				{
					tabealu[i*num_rows + j] = - tabealu[i*num_rows + j];
				}
			}
		}

		cout << endl;
		cout << "Tableau,   x,   " ;
		for ( int j = 0; j < num_cols ; j ++)
		{
			cout << varNames[j] << ",  ";
		}
		cout << endl;
		for (int i = 0; i < num_rows; i++)
		{
			if (head[i] < 0)
			{
				cout << head[i] <<",  " << head_val[i]<< ",  ";
			}
			else{
				cout << varNames[head[i]] <<",  " << head_val[i]<< ",  ";
			}
			for (int j = 0; j < num_cols + num_rows; j ++)
			{
				cout << tabealu[j*num_rows +i] << ", " ;
			}
			cout<<endl;
		}
		cout << endl;
		cout.rdbuf( old_buf);

		delete []head;
		delete []head_val;
		delete objval_p;
		delete []sense;
		free_and_null((char**) &tabealu);
		return 1;

}

int OptModel::getNumRows(void)
{
	return CPXgetnumrows(cpxenv, cpxlp);
}

int OptModel::getNumCols(void)
{
	return CPXgetnumcols(cpxenv, cpxlp);
}

bool OptModel::isIntegerFeasible(int varRange)
{
	int num_cols  = CPXgetnumcols(cpxenv,  cpxlp) ; 
	double *val_vars = new double[num_cols];
	CPXgetx(cpxenv, cpxlp, val_vars, 0, num_cols-1);
	
	if (varRange == 0)
	{
		for (int i = 0; i < idx_Z.size(); i++)
		{
			int idx = idx_Z[i];
			if (fabs(fmod(val_vars[idx],1)) > SMALL_DECIMAL)
			{
				delete []val_vars;
				return false;
			}
		}
	}
	else if (varRange == 1)
	{
		for (int i = 0; i < num_cols; i++)
		{
			if (fabs(fmod(val_vars[i],1)) > SMALL_DECIMAL)
			{
				delete []val_vars;
				return false;
			}
		}
	}
	delete []val_vars;
	return true;
}


int OptModel::getNumAllCuts()
{
	int numCuts = 0;
	int cntCut = 0;

	CPXgetnumcuts(cpxenv, cpxlp, CPX_CUT_CLIQUE, &cntCut);
	numCuts += cntCut;
	CPXgetnumcuts(cpxenv, cpxlp, CPX_CUT_COVER, &cntCut);
	numCuts += cntCut;
	CPXgetnumcuts(cpxenv, cpxlp, CPX_CUT_DISJ, &cntCut);
	numCuts += cntCut;
	CPXgetnumcuts(cpxenv, cpxlp, CPX_CUT_FLOWCOVER, &cntCut);
	numCuts += cntCut;
	CPXgetnumcuts(cpxenv, cpxlp, CPX_CUT_FLOWPATH, &cntCut);
	numCuts += cntCut;
	CPXgetnumcuts(cpxenv, cpxlp, CPX_CUT_FRAC, &cntCut);
	numCuts += cntCut;
	CPXgetnumcuts(cpxenv, cpxlp, CPX_CUT_GUBCOVER, &cntCut);
	numCuts += cntCut;
	CPXgetnumcuts(cpxenv, cpxlp, CPX_CUT_IMPLBD, &cntCut);
	numCuts += cntCut;
	CPXgetnumcuts(cpxenv, cpxlp, CPX_CUT_LOCALCOVER, &cntCut);
	numCuts += cntCut;
	CPXgetnumcuts(cpxenv, cpxlp, CPX_CUT_MIR, &cntCut);
	numCuts += cntCut;
	CPXgetnumcuts(cpxenv, cpxlp, CPX_CUT_OBJDISJ, &cntCut);
	numCuts += cntCut;
	CPXgetnumcuts(cpxenv, cpxlp, CPX_CUT_SOLNPOOL, &cntCut);
	numCuts += cntCut;
	CPXgetnumcuts(cpxenv, cpxlp, CPX_CUT_TABLE, &cntCut);
	numCuts += cntCut;
	CPXgetnumcuts(cpxenv, cpxlp, CPX_CUT_TIGHTEN, &cntCut);
	numCuts += cntCut;
	CPXgetnumcuts(cpxenv, cpxlp, CPX_CUT_USER, &cntCut);
	numCuts += cntCut;
	CPXgetnumcuts(cpxenv, cpxlp, CPX_CUT_ZEROHALF, &cntCut);
	numCuts += cntCut;

	return numCuts;

}



int OptModel::changeMIP2LP(void)
{
	int status = CPXchgprobtype(cpxenv, cpxlp, CPXPROB_LP);
	return status;
	/*
	int cnt = idx_Z.size();
	int *indices = new int[idx_Z.size()];
	for (int i = 0; i < idx_Z.size(); i++)
	{ 
		indices[i] = idx_Z[i];
	}
	char* xctype= new char[cnt];
	for (int i = 0; i < idx_Z.size(); i ++)
	{
		xctype[i] = CPX_CONTINUOUS;
	}
	int status = 0;
	status = CPXchgctype( cpxenv,  cpxlp,  cnt,   indices,   xctype);
	if ( status ) {
		fprintf (stderr, "Failed to change variable type.\n");
	}
	delete[] indices;
	modelType= MODEL_TYPE_LP;
	return status;
	*/
}

int OptModel::changeProbTypeMIP2LP(void)
{
	int status = CPXchgprobtype(cpxenv,cpxlp, CPXPROB_LP) ;
	if ( status ) {
		fprintf (stderr, "Failed to change variable type.\n");
	}
	modelType= MODEL_TYPE_LP;
	return status;
}


int OptModel::populateModelFromFile(string lpFileName)
{
	if (cpxenv == NULL)
		initilizeCPXenv();
	initilizeCPXlp("modelonfile");
	int status = CPXreadcopyprob (cpxenv, cpxlp, lpFileName.c_str(), NULL);
	if(status !=0)
		return -1;
//	mapVarNameIndex();
	mapVarIndexName();
	mapZIdxRowNum();
	mapScenarioRownum();
	char* buf_str = new char[50];
	int surplus_p;
	CPXgetprobname(cpxenv, cpxlp, buf_str, 50, &surplus_p);

	string probName(buf_str);

	DataModel::r = getInstanceFloatParameterValue(probName, "r");
	if (DataModel::k == 0 )
	{
		DataModel::k = getInstanceIntParameterValue(probName, "k");
	}

/*
	// find value of k
	size_t found = lpFileName.find('k');
	char str[40];
	lpFileName.copy(str, 100 , found+1);

	char * pEnd;
	long int subscript = strtol(str , &pEnd, 10);
	DataModel::k = subscript;
	*/


	// n
	DataModel::n = idx_X.size();
	// m
	DataModel::m = idx_Z.size();
	// NUM_COLS
	DataModel::NUM_COLS = DataModel::m + DataModel::n; 

	if (CPXgetprobtype(cpxenv, cpxlp) == CPXPROB_LP)
	{
		modelType = MODEL_TYPE_LP;
	}
	else
		modelType = MODEL_TYPE_MIP;

	return 0;

}

void  OptModel::mapVarNameIndex()
{
	idx_Z.clear();
	idx_X.clear();
	idx_S.clear();
	int num_rows = CPXgetnumrows(cpxenv,  cpxlp);
	int num_cols  = CPXgetnumcols(cpxenv,  cpxlp) ; // need to know why
	int surplus;

	int status = CPXgetcolname (cpxenv, cpxlp, NULL, NULL, 0, &surplus, 0,num_cols-1);
	int cur_colnamespace = - surplus;
	if ( cur_colnamespace > 0 ) {
		char          **cur_colname      = (char **) malloc (sizeof(char *)*(num_cols));
		char          *cur_colnamestore = (char *)  malloc (cur_colnamespace);
		if ( cur_colname      == NULL ||
			cur_colnamestore == NULL   ) {
				fprintf (stderr, "Failed to get memory for column names.\n");
				status = -1;
				cout<< endl <<" error in printTableau"<<endl;
				return;
		}
		status = CPXgetcolname (cpxenv, cpxlp, cur_colname, cur_colnamestore, 
			cur_colnamespace, &surplus, 0, num_cols-1);
		// map name with column index
		map<int,int> mapIndex_X;
		map<int,int> mapIndex_Z;
		map<int,int> mapIndex_S;
		for (int i = 0; i <num_cols ; i ++)
		{
			if (cur_colname[i][0] == 'x')
			{
				char * pEnd;
				long int subscript = strtol (cur_colname[i]+1,&pEnd,10);
				mapIndex_X.insert(make_pair(subscript,i));
			}
			else if (cur_colname[i][0] == 'z')
			{
				char * pEnd;
				int subscript = strtol (cur_colname[i]+1,&pEnd,10);
				mapIndex_Z.insert(make_pair(subscript,i));
			}
			else if (cur_colname[i][0] == 's')
			{
				char * pEnd;
				int subscript = strtol (cur_colname[i]+1,&pEnd,10);
				mapIndex_S.insert(make_pair(subscript,i));
			}
			string varName = cur_colname[i];

			var_name_idx_table.insert(make_pair(varName, i));
		}

		for (map<int,int>::iterator itr  = mapIndex_X.begin(); itr != mapIndex_X.end(); ++itr)
		{
			idx_X.push_back(itr->second);
		}
		for (map<int,int>::iterator itr  = mapIndex_Z.begin(); itr != mapIndex_Z.end(); ++itr)
		{
			idx_Z.push_back(itr->second);
		}
		for (map<int,int>::iterator itr  = mapIndex_S.begin(); itr != mapIndex_S.end(); ++itr)
		{
			idx_S.push_back(itr->second);
		}
		free_and_null( (char **)&cur_colname);
		free_and_null((char **)&cur_colnamestore);
	}

}


int OptModel::setCutCallBackFunction(int (CPXPUBLIC *cutcallback) (CALLBACK_CUT_ARGS), void *cbhandle)
{
    int status = 0;//CPXsetcutcallbackfunc(cpxenv, cutcallback, cbhandle);
	return status;
}


int OptModel::setIncumbentCallBackFunction(int (CPXPUBLIC *incumbentcallback) (CALLBACK_INCUMBENT_ARGS), void *cbhandle)
{
	int status = CPXsetincumbentcallbackfunc(cpxenv, incumbentcallback, cbhandle);
	return status;
}

int OptModel::setNodeCallBackFunction(int(CPXPUBLIC *nodecallback) (CALLBACK_NODE_ARGS), void *cbhandle)
{
	int status = CPXsetnodecallbackfunc(cpxenv, nodecallback, cbhandle);
	return status;
}


int OptModel::setBranchCallBackFunction(int(CPXPUBLIC *branchcallback) (CALLBACK_BRANCH_ARGS), void *cbhandle)
{
	int status = CPXsetbranchcallbackfunc(cpxenv, branchcallback, cbhandle);
	return status;
}

int OptModel::setSolveCallBackFunction(int(CPXPUBLIC *solvecallback) (CALLBACK_SOLVE_ARGS), void *cbhandle)
{
	int status = CPXsetsolvecallbackfunc(cpxenv, solvecallback, cbhandle);
	return status;
}

int OptModel::populateModelWithReformulation( DataModel& dm, OptModel::MODEL_TYPE modeltype)
{
	idx_X.clear();
	idx_Z.clear();
	idx_S.clear();
	double *obj     = NULL;
	double *lb      = NULL;
	double *ub      = NULL;
	char   *ctype   = NULL;
	int    *rmatind = NULL;
	double *rmatval = NULL;
	char** colname = NULL;

	int    status = 0;

	CPXchgobjsen (cpxenv, cpxlp, 1); /* Minimization problem */

	/* Allocate colcnt-sized arrays */

	obj     = (double *) malloc (DataModel::NUM_COLS * sizeof(double));
	lb      = (double *) malloc (DataModel::NUM_COLS * sizeof(double));
	ub      = (double *) malloc (DataModel::NUM_COLS * sizeof(double));
	ctype   = (char *)   malloc (DataModel::NUM_COLS * sizeof(char));
	colname = (char**) malloc(DataModel::NUM_COLS*sizeof(char*));

	if ( obj     == NULL ||
		lb      == NULL ||
		ub      == NULL ||
		ctype   == NULL ) {
			fprintf (stderr, "Could not allocate colcnt arrays\n");
			status = CPXERR_NO_MEMORY;
			return 0;
	}

	string* strColNames = new string[DataModel::NUM_COLS];
	for (int i =0; i < DataModel::NUM_COLS; i ++)
	{
		if (i < DataModel::n )
		{
			obj[i] = dm.costVec[i].first;                      // fill in obj coefficients
			idx_X.push_back(i);                        // store index for x variables
			strColNames[i] = "x"+ToString(i);
			colname[i] = new char[8];
			strcpy(colname[i], strColNames[i].c_str()); // make name for x
			ctype[i] = 'C';
		}
		else{
			obj[i] = 0;
			idx_Z.push_back(i);
			strColNames[i] = "z" + ToString(i-DataModel::n);
			colname[i] = new char[8];
			strcpy(colname[i], strColNames[i].c_str()); // make name for z
			if (modeltype == MODEL_TYPE_LP)
				ctype[i] = 'C';
			else if( modeltype == MODEL_TYPE_MIP)
				ctype[i] = 'I';
		}
	}
	delete []strColNames;

	for (int i =0; i < DataModel::NUM_COLS; i ++)       // fill in upper bound
	{
		if (i<DataModel::n)
			ub[i] = CPX_INFBOUND;
		else
			ub[i] = CPX_INFBOUND;//1;
	}

	if (modeltype == MODEL_TYPE_LP)
	{
		status = CPXnewcols (cpxenv, cpxlp, DataModel::NUM_COLS, obj, NULL, NULL, NULL, colname);
	}
	else if (modeltype == MODEL_TYPE_MIP)
	{
		status = CPXnewcols (cpxenv, cpxlp, DataModel::NUM_COLS, obj, NULL, NULL, ctype, colname);
	}
	if ( status ) {
		fprintf (stderr, "Could not add new columns.\n");
		return 0;
	}

	modelType = modeltype;

	/* adding rows to model */
	for (int i = 0; i < dm.Arows.size(); i++)  // the last row is the resource constrain
	{
		string s_name = "s"+ToString(i);
		colname[0] = new char[8];
		strcpy(colname[0], s_name.c_str());

		/* add new column for slack var*/
		status = CPXnewcols (cpxenv, cpxlp, 1, NULL, NULL, NULL, NULL, colname);
		if ( status ) {
			fprintf (stderr, "Could not add new columns.\n");
			return 0;
		}

		int totaNumCols = CPXgetnumcols (cpxenv, cpxlp) ;
		idx_S.push_back(totaNumCols-1);  // store s index

		/* prepare row data*/

		rmatind = (int * )   malloc (totaNumCols * sizeof(int));
		rmatval = (double *) malloc (totaNumCols * sizeof(double));

		VEC_MTX_ROW rowEntries = dm.Arows[i];  // fill in the entries for a row
		double rhs[1];
		char   sense[1];

		switch (dm.senses[i])
		{
		case 'L':
			for (int j =0; j <rowEntries.size(); j ++)
			{
				PAIR_VALUE_INDEX entryPair = rowEntries[j];
				rmatind[j] = entryPair.second;
				rmatval[j] = - entryPair.first;
			}
			rmatind[rowEntries.size()] = totaNumCols-1; // add slack var
			rmatval[rowEntries.size()] = -1;                     // add slack var
			rhs[0] = -dm.rhs[i];
			sense[0] = 'E';

			break;
		case 'G':
			for (int j =0; j <rowEntries.size(); j ++)
			{
				PAIR_VALUE_INDEX entryPair = rowEntries[j];
				rmatind[j] = entryPair.second;
				rmatval[j] = entryPair.first;
			}
			rmatind[rowEntries.size()] = totaNumCols-1; // add slack var
			rmatval[rowEntries.size()] = -1;                     // add slack var
			rhs[0] = dm.rhs[i];
			sense[0] = 'E';

			break;
		case 'E':
			for (int j =0; j <rowEntries.size(); j ++)
			{
				PAIR_VALUE_INDEX entryPair = rowEntries[j];
				rmatind[j] = entryPair.second;
				rmatval[j] = entryPair.first;
			}
			rmatind[rowEntries.size()] = totaNumCols-1; // add slack var
			rmatval[rowEntries.size()] = 0;                     // add slack var
			rhs[0] = dm.rhs[i];
			sense[0] = 'E';

			break;
		}
		int    rmatbeg[1];
		rmatbeg[0] = 0;
		status = CPXaddrows (cpxenv, cpxlp, 0, 1, rowEntries.size() + 1, rhs,
			sense, rmatbeg, rmatind, rmatval,
			NULL, NULL);
		if ( status ) {
			fprintf (stderr, "Could not add new rows.\n");
			return 0;
		}

		free_and_null((char **)&rmatind); 
		free_and_null((char **)&rmatval); 

		// add the reformulation row

// 		if (i == dm.Arows.size() -1 )// the last row is resource constrain, so skip it
// 		{
// 			break;
// 		}

		s_name = "s"+ToString(DataModel::m+i);
		colname[0] = new char[8];
		strcpy(colname[0], s_name.c_str());

		/* add new column for slack var*/
		status = CPXnewcols (cpxenv, cpxlp, 1, NULL, NULL, NULL, NULL, colname);
		if ( status ) {
			fprintf (stderr, "Could not add new columns.\n");
			return 0;
		}

		totaNumCols = CPXgetnumcols (cpxenv, cpxlp) ;
		idx_S.push_back(totaNumCols-1);  // store s index

		/* prepare row data*/

		rmatind = (int * )   malloc (totaNumCols * sizeof(int));
		rmatval = (double *) malloc (totaNumCols * sizeof(double));

		rowEntries = dm.Arows[i];  // fill in the entries for a row
// 		double rhs[1];
// 		char   sense[1];

		double max_a_i = 0; 
		for (int k = 0; k < rowEntries.size()-1; k++)
		{
			PAIR_VALUE_INDEX entryPair = rowEntries[k];
			if (entryPair.first > max_a_i)
				max_a_i = entryPair.first;
		}

		for (int j =0; j <rowEntries.size(); j ++)
		{
			PAIR_VALUE_INDEX entryPair = rowEntries[j];
			if(j == DataModel::n + i)
			{
				rmatind[j] = entryPair.second;
				rmatval[j] = 1/SMALL_DECIMAL;
				//break;
			}
			else{
				rmatind[j] = entryPair.second;
				rmatval[j] = entryPair.first;
			}
		}
		rmatind[rowEntries.size()] = totaNumCols-1; // add slack var
		rmatval[rowEntries.size()] = 1;                     // add slack var
		rhs[0] = 1/SMALL_DECIMAL + dm.rhs[i];
		sense[0] = 'E';

// 		int    rmatbeg[1];
		rmatbeg[0] = 0;
		status = CPXaddrows (cpxenv, cpxlp, 0, 1, rowEntries.size() + 1, rhs,
			sense, rmatbeg, rmatind, rmatval,
			NULL, NULL);
		if ( status ) {
			fprintf (stderr, "Could not add new rows.\n");
			return 0;
		}

		free_and_null((char **)&rmatind); 
		free_and_null((char **)&rmatval); 

	}

	/* adding cardinality constraint */
	/* add new column for slack var*/
	lb[0] = 0;
	ub[0] = 0; // exactly k violated
	string s_name = "s"+ToString(idx_S.size());
	colname[0] = new char[8];
	strcpy(colname[0], s_name.c_str());
	status = CPXnewcols (cpxenv, cpxlp, 1, NULL, lb, ub, NULL, colname);
	if ( status ) {
		fprintf (stderr, "Could not add new columns.\n");
		return 0;
	}
	int totaNumCols = CPXgetnumcols (cpxenv, cpxlp) ;
	idx_S.push_back(totaNumCols-1);  // store s index

	/* prepare row data*/

	rmatind = (int * )   malloc (totaNumCols * sizeof(int));
	rmatval = (double *) malloc (totaNumCols * sizeof(double));

	for (int j =0; j <DataModel::m; j ++)
	{
		rmatind[j] = idx_Z[j];
		rmatval[j] = 1;
	}
	rmatind[DataModel::m] = totaNumCols-1; // add slack var
	rmatval[DataModel::m] = 1;                     // add slack var
	int    rmatbeg[1];
	double rhs[1];
	char   sense[1];
	rhs[0] = DataModel::k;
	sense[0] = 'E';
	rmatbeg[0] = 0;
	status = CPXaddrows (cpxenv, cpxlp, 0, 1, DataModel::m + 1, rhs,
		sense, rmatbeg, rmatind, rmatval,
		NULL, NULL);
	if ( status ) {
		fprintf (stderr, "Could not add new rows.\n");
		return 0;
	}

	free_and_null((char **)&rmatind); 
	free_and_null((char **)&rmatval); 


	// add z <= 1 
	rmatind = (int * )   malloc (2 * sizeof(int));
	rmatval = (double *) malloc (2 * sizeof(double));

	for (int i =0 ; i < DataModel::m; i++)
	{
		string s_name = "s"+ToString(idx_S.size());
		colname[0] = new char[8];
		strcpy(colname[0], s_name.c_str());

		/* add new column for slack var*/
		status = CPXnewcols (cpxenv, cpxlp, 1, NULL, NULL, NULL, NULL, colname);
		if ( status ) {
			fprintf (stderr, "Could not add new columns.\n");
			return 0;
		}

		int totaNumCols = CPXgetnumcols (cpxenv, cpxlp) ;
		idx_S.push_back(totaNumCols-1);  // store s index

		/* prepare row data*/


		rmatind[0] = idx_Z[i];
		rmatval[0] = 1;

		rmatind[1] = totaNumCols-1; // add slack var
		rmatval[1] = 1;

		int    rmatbeg[1];
		double rhs[1];
		char   sense[1];
		rhs[0] = 1;
		sense[0] = 'E';
		rmatbeg[0] = 0;
		status = CPXaddrows (cpxenv, cpxlp, 0, 1, 2, rhs,
			sense, rmatbeg, rmatind, rmatval,
			NULL, NULL);
		if ( status ) {
			fprintf (stderr, "Could not add new rows.\n");
			return 0;
		}

	}

	string fileName ;
	if (modeltype == MODEL_TYPE_LP)
		fileName = modelName + "R.lp" ;
	else if (modeltype == MODEL_TYPE_MIP)
		fileName = modelName + "MIP.lp";
	writeModel2File(fileName.c_str());
	//CPXwriteprob(cpxenv, lp, "initialLPRelaxation.lp", NULL);

	free_and_null((char **)&obj);
	free_and_null((char **)&lb);
	free_and_null((char **)&ub);
	free_and_null((char **)&ctype);
	free_and_null((char **)&colname); 
	free_and_null((char **)&rmatind); 
	free_and_null((char **)&rmatval); 

	return 1;

}



int OptModel::addReformAsCut( DataModel& dm, OptModel::CUT_TYPE cutType)
{
	int    *rmatind = NULL;
	double *rmatval = NULL;

	int    status = 0;
	int totaNumCols = DataModel::n + DataModel::m;

	/* adding rows to model */
	for (int i = 0; i < dm.Arows.size(); i++) 
	{
		/* prepare row data*/

		rmatind = (int * )   malloc (totaNumCols * sizeof(int));
		rmatval = (double *) malloc (totaNumCols * sizeof(double));

		VEC_MTX_ROW rowEntries = dm.Arows[i];  // fill in the entries for a row
		double rhs[1];
		char   sense[1];

		double max_a_i = 0; 
		for (int k = 0; k < rowEntries.size(); k++)
		{
			PAIR_VALUE_INDEX entryPair = rowEntries[k];
			if (entryPair.first > max_a_i)
				max_a_i = entryPair.first;
		}

		for (int j =0; j <rowEntries.size(); j ++)
		{
			PAIR_VALUE_INDEX entryPair = rowEntries[j];
			if(j == DataModel::n + i)
			{
				rmatind[j] = entryPair.second;
				rmatval[j] = max_a_i;
			}
			else{
				rmatind[j] = entryPair.second;
				rmatval[j] = entryPair.first;
			}
		}
		rhs[0] = max_a_i + dm.rhs[i];
		sense[0] = 'L';

		int    rmatbeg[1];
		rmatbeg[0] = 0;
		if (cutType == OptModel::CUT_AS_USERCUT)
			status = CPXaddusercuts (cpxenv, cpxlp, 1, rowEntries.size() , rhs,	sense, rmatbeg, rmatind, rmatval, NULL);
		// 			status = CPXaddusercuts(cpxenv, mip, 1, ind, rhs, sense, rmatbeg, rmatind, rmatval, row_name);
		else  if(cutType == OptModel::CUT_AS_FORMULATION)
			status = CPXaddrows (cpxenv,  cpxlp, 0, 1, rowEntries.size() , rhs,	sense, rmatbeg, rmatind, rmatval, NULL, NULL);
		else if(cutType == OptModel::CUT_AS_LAZYCONSTRAIN)
			status = CPXaddlazyconstraints(cpxenv, cpxlp, 1 , rowEntries.size() ,  rhs,  sense,  rmatbeg,  rmatind, rmatval, NULL);

		if ( status ) {
			fprintf (stderr, "Could not add new rows.\n");
			return 0;
		}

		free_and_null((char **)&rmatind); 
		free_and_null((char **)&rmatval); 

	}

	return 1;

}

OptModel::OptModel(CPXENVptr _cpxenv, CPXLPptr  _cpxlp)
{
	cpxenv = _cpxenv;
	cpxlp = _cpxlp;
	num_one_at_curr_node = 0;
	num_rootnode_cut = 0;
	num_IC_rootnode = 0;
	num_MIX_rootnode = 0;
	MIP_optimal_value = 0;
	Prev_MIP_optimal_value = 0;
	LP_optimal_value =0 ;
	Prev_LP_optimal_value = 0;
	ICCutPool.setCutIdentifier("IC");
	MixCutPool.setCutIdentifier("MIX");

}


OptModel::OptModel(CPXENVptr _cpxenv, string _modelName)
{
	cpxenv = _cpxenv;
	is_a_copy = true;
	initilizeCPXlp(_modelName);
	num_one_at_curr_node = 0;
	num_rootnode_cut = 0;
	num_IC_rootnode = 0;
	num_MIX_rootnode = 0;
	MIP_optimal_value = 0;
	Prev_MIP_optimal_value = 0;
	LP_optimal_value =0 ;
	Prev_LP_optimal_value = 0;
	ICCutPool.setCutIdentifier("IC");
	MixCutPool.setCutIdentifier("MIX");

}

void OptModel::setDataModel(DataModel* dm)
{
	dm_p = dm;
}

// int OptModel::prepareCallbackCutCoefsReform(CPXCENVptr env,  int branch_index, int &nzcnt, double &rhs, int &sense, int * &cutind, double * &cutval)
// {
// 	int rowNum = 0;
// 	for (int i = 0 ; i < idx_Z.size(); i++)
// 	{
// 		if (idx_Z[i]==branch_index)
// 		{
// 			rowNum = i;
// 			break;
// 		}
// 	}
// 
// 	int    status = 0;
// 	int totaNumCols = DataModel::n ;
// 
// 	/* prepare row data*/
// 
// 	cutind = (int * )   malloc (totaNumCols * sizeof(int));
// 	cutval = (double *) malloc (totaNumCols * sizeof(double));
// 
// 	VEC_MTX_ROW rowEntries = dm_p->Arows[rowNum];  // fill in the entries for a row
// 
// 	nzcnt = DataModel::n;
// 
// 	for (int j =0; j <DataModel::n; j ++)
// 	{
// 		PAIR_VALUE_INDEX entryPair = rowEntries[j];
// 		cutind[j] = entryPair.second;
// 		cutval[j] = entryPair.first;
// 	}
// 	rhs = dm_p->rhs[rowNum];
// 	sense = 'L';
// 
// 
// 	return 1;
// }



int OptModel::prepareCallbackCutCoefsReform(CPXCENVptr env,  int branch_index, int &nzcnt, double &rhs, int &sense, int * &cutind, double * &cutval)
{
	int rowNum = 0;

	for (int i = 0 ; i < idx_Z.size(); i++)
	{
		if (idx_Z[i]==branch_index)
		{
			rowNum = i;
			break;
		}
	}
	string scenarioName = "SN" + ToString(rowNum);
	int    status = 0;
	int surplus;
	int totaNumCols = DataModel::n ;

	/* prepare row data*/
	int num_rows = CPXgetnumrows(cpxenv, cpxlp);
	int num_cols = CPXgetnumcols(cpxenv, cpxlp);

	cutind = (int * )   malloc (num_cols * sizeof(int));
	cutval = (double *) malloc (num_cols * sizeof(double));


	status = CPXgetrowname (cpxenv, cpxlp, NULL, NULL, 0, &surplus, 0,num_rows-1);
	int cur_rownamespace = - surplus;
	if ( cur_rownamespace > 0 ) 
	{
		char          **cur_rowname = (char **) malloc (sizeof(char *)*(num_rows));
		char          *cur_rownamestore = (char *)  malloc (cur_rownamespace);
		if ( cur_rowname      == NULL ||
			cur_rownamestore == NULL   ) {
				fprintf (stderr, "Failed to get memory for column names.\n");
				status = -1;
				return 1;
		}
		status = CPXgetrowname (cpxenv, cpxlp, cur_rowname, cur_rownamestore, 
			cur_rownamespace, &surplus, 0, num_rows-1);
		// map name with column index
		for (int r = 0; r <num_rows ; r ++)
		{
			string rowName = cur_rowname[r];
			if (rowName != scenarioName )// only look at scenario constrains
				continue;

			int t_nzcnt;
			int * t_rmatbeg = new int(0);
			int * t_rmatind = new int[num_cols]();
			double * t_rmatval = new double[num_cols]();
			int t_rmatspace = num_cols;
			int t_surplus; 
			double t_rhs;

			CPXgetrows(cpxenv, cpxlp, &t_nzcnt, t_rmatbeg, t_rmatind, t_rmatval, t_rmatspace, &t_surplus, r, r);// get row coefficients from original formulation
			CPXgetrhs(cpxenv, cpxlp, &t_rhs, r, r);
			int counter = 0;
			for (int i = 0; i < t_nzcnt; i ++)
			{
				if ( t_rmatind[i] >= DataModel::n) // only read x coefficients
					continue; 

				cutind[counter] = t_rmatind[i];
				cutval[counter] = t_rmatval[i];
				counter ++;
			}
			nzcnt = counter; 
			sense = 'L';
			rhs = t_rhs;

			break;
		}

		free_and_null((char**) &cur_rowname);
		free_and_null((char**) &cur_rownamestore);

		return 0;
	}
	return 1;
}



int OptModel::prepareCallbackCutCoefsReformNew(CPXCENVptr env,  int branch_index, int &nzcnt, double &rhs, int &sense, int * &cutind, double * &cutval)
{
	int rowNum = 0;
	int num_cols = getNumCols();

	map<int,int>::iterator itr = map_z_idx_row_num.find(branch_index); // note that map_z_idx_row_num might have not been initiated yet.
	if (itr != map_z_idx_row_num.end())
	{
		rowNum = itr->second;
	}
	else
		return 0;

	int    status = 0;
	int totaNumCols = DataModel::n ;


	int t_nzcnt;
	int * t_rmatbeg = new int(0);
	int * t_rmatind = new int[num_cols]();
	double * t_rmatval = new double[num_cols]();
	cutind = new int[totaNumCols];
	cutval = new double[totaNumCols];
	int t_rmatspace = num_cols;
	int t_surplus; 
	double t_rhs;

	CPXgetrows(cpxenv, cpxlp, &t_nzcnt, t_rmatbeg, t_rmatind, t_rmatval, t_rmatspace, &t_surplus, rowNum, rowNum);// get row coefficients from original formulation
	CPXgetrhs(cpxenv, cpxlp, &t_rhs, rowNum, rowNum);
	int counter = 0;
	for (int i = 0; i < t_nzcnt; i ++)
	{
		if ( t_rmatind[i] >= DataModel::n) // only read x coefficients
			continue; 

		cutind[counter] = t_rmatind[i];
		cutval[counter] = t_rmatval[i];
		counter ++;
	}
	nzcnt = counter; 
	sense = 'L';
	rhs = t_rhs;
	delete t_rmatbeg ;
	delete[] t_rmatind;
	delete[] t_rmatval;
	return 1;
}




int OptModel::callbackAddMixingSet(void * cbdata, int wherefrom, BRANCHCALLBACKINFO*nodeInfo)
{
	static int counter = 0;
	int num_cols = CPXgetnumcols(cpxenv, cpxlp);
	double* nodex= new double[num_cols];
	int status = CPXgetcallbacknodex(cpxenv, cbdata, wherefrom, nodex, 0, num_cols-1);

	string constraintType("MIR");
	vector<int> inds_violated = checkViolatedConstrainsByX(constraintType, nodex);

	if (inds_violated.empty() )
	{
		delete []nodex;
		return false;
	}

	for (int i = 0; i <=  inds_violated.size(); i++)// -1 means cost vector
	{
		int ind_violated =  i == inds_violated.size() ? -1 : inds_violated[i];
		if (ind_violated != -1)
		{
			int sceNum = findSceNumByRowNum(ind_violated);
			int z_index = idx_Z[sceNum];
			vector<int>::iterator itr1  = find(nodeInfo->vec_index_ones.begin(), nodeInfo->vec_index_ones.end(), z_index);
			vector<int>::iterator itr0  = find(nodeInfo->vec_index_zeros.begin(), nodeInfo->vec_index_zeros.end(), z_index);
			if ( itr0 != nodeInfo->vec_index_zeros.end() || itr1 != nodeInfo->vec_index_ones.end() )
				continue;
		}

		CUT_COEFICIENTS costAlpha = readAlphaFromConstrain(ind_violated);
		vector<CUT> cuts /*= genertaeCMIXcuts(ind_violated, costAlpha, nodex, nodeInfo)*/;
		if (cuts.empty())
		{
			continue;
		}
		else
		{
			// add cuts to master model
			for (int i = 0; i < cuts.size(); i++)
			{
				string str = "MIX"+ToString(num_MIX_rootnode);
				callbackAddACut(cbdata, wherefrom, cuts[i]);
				num_MIX_rootnode++;
				// 				writeModel2File("masterModel.lp");
			}

			delete []nodex;
			return cuts.size();	

		}
	}


	return -1;


}

int OptModel::callbackAddReverseCHCuts(void * cbdata, int wherefrom,int branch_index)
{
	int rowNum = 0;
	for (int i = 0 ; i < idx_Z.size(); i++)
	{
		if (idx_Z[i]==branch_index)
		{
			rowNum = i;
			break;
		}
	}
	string scenarioName = "MIR" + ToString(rowNum);
	int    status = 0;
	int surplus;
	int totaNumCols = DataModel::n ;

	/* prepare row data*/
	int num_rows = CPXgetnumrows(cpxenv, cpxlp);
	int num_cols = CPXgetnumcols(cpxenv, cpxlp);

	int *cutind = new int[ DataModel::n ];
	double* cutval = new double [DataModel::n];
	int nzcnt ;
	double rhs;

	status = CPXgetrowname (cpxenv, cpxlp, NULL, NULL, 0, &surplus, 0,num_rows-1);
	int cur_rownamespace = - surplus;
	char          **cur_rowname = (char **) malloc (sizeof(char *)*(num_rows));
	char          *cur_rownamestore = (char *)  malloc (cur_rownamespace);
	if ( cur_rownamespace > 0 ) 
	{
		if ( cur_rowname      == NULL ||
			cur_rownamestore == NULL   ) {
				fprintf (stderr, "Failed to get memory for column names.\n");
				status = -1;
				return 1;
		}
		status = CPXgetrowname (cpxenv, cpxlp, cur_rowname, cur_rownamestore, 
			cur_rownamespace, &surplus, 0, num_rows-1);
		// map name with column index
		for (int r = 0; r <num_rows ; r ++)
		{
			string rowName = cur_rowname[r];
			if (rowName != scenarioName )// only look at scenario constrains
				continue;

			int t_nzcnt;
			int * t_rmatbeg = new int(0);
			int * t_rmatind = new int[num_cols]();
			double * t_rmatval = new double[num_cols]();
			int t_rmatspace = num_cols;
			int t_surplus; 
			double t_rhs;

			CPXgetrows(cpxenv, cpxlp, &t_nzcnt, t_rmatbeg, t_rmatind, t_rmatval, t_rmatspace, &t_surplus, r, r);// get row coefficients from original formulation
			CPXgetrhs(cpxenv, cpxlp, &t_rhs, r, r);
			int counter = 0;
			for (int i = 0; i < t_nzcnt; i ++)
			{
				if ( t_rmatind[i] >= DataModel::n) // only read x coefficients
					continue; 

				cutind[counter] = t_rmatind[i];
				cutval[counter] = t_rmatval[i];
				counter ++;
			}
 			nzcnt = counter; 
			char sense = 'L';
			rhs = t_rhs;
			CPXcutcallbackaddlocal(cpxenv, cbdata, wherefrom, nzcnt, rhs, sense, cutind, cutval);
			delete [] t_rmatbeg;
			delete [] t_rmatind;
			break;
		}

	}


	// add ConvexHull cuts
	double *nodex = new double[num_cols];
	status = CPXgetcallbacknodex(cpxenv, cbdata, wherefrom, nodex, 0, num_cols -1);

	double *b_j = new double[DataModel::n];
	memset( b_j , 0 , sizeof(double)*DataModel::n);

	// fill in b_j
	for (int i = 0 ; i < nzcnt; i++)
	{
		int index = cutind[i];
		b_j[index] = cutval[i];
	}
	delete []cutval;
	delete []cutind;

	double b_0 = rhs;
	for (int r = 0 ; r < num_rows; r ++)
	{
		string rowName = cur_rowname[r];
		if (rowName.find("MIR") == string::npos)// only look at scenario constrains
			continue; 

		int t_nzcnt;
		int * t_rmatbeg = new int(0);
		int * t_rmatind = new int[num_cols]();
		double * t_rmatval = new double[num_cols]();
		int t_rmatspace = num_cols;
		int t_surplus; 
		double t_rhs;

		CPXgetrows(cpxenv, cpxlp, &t_nzcnt, t_rmatbeg, t_rmatind, t_rmatval, t_rmatspace, &t_surplus, r, r);// get row coefficients from original formulation
		CPXgetrhs(cpxenv, cpxlp, &t_rhs, r, r);

		double * a_j = new double[DataModel::n];
		int * index_a_j = new int[DataModel::n];
		memset(a_j, 0 , sizeof(double)*DataModel::n);
		double a_0 = t_rhs;

		int index_z ;
		for (int i = 0; i < t_nzcnt; i++) // fill coefs for a_j
		{
			int index = t_rmatind[i];
			if (index >= DataModel::n)
			{
				index_z = index;
				continue;
			}
			a_j[index] = t_rmatval[i];
			index_a_j[index] = i;
		}
		
		// scale a_j
		multimap<double, int > sorted_t_j;
		for(int i = 0; i < DataModel::n; i ++)
		{
			if (a_j[i] < SMALL_DECIMAL || b_j[i] < SMALL_DECIMAL)
				continue;

			sorted_t_j.insert(make_pair(b_0*a_j[i]/b_j[i], i));
		}

		//for the first k a_j that is strictly smaller than rhs
		for (multimap<double, int>::iterator itr = sorted_t_j.begin(); itr != sorted_t_j.end(); ++itr)// for every a_i < a_0, generate one facet
		{
			if (itr->first >= a_0) 
				continue;
			// prepare cut coefficients :  add coef for i such that a_i and b_i are both non-zero.
			CUT_COEFICIENTS coefs;
			int index_i = itr->second;
			for (multimap<double, int>::iterator itrBigger = ++itr ; itrBigger != sorted_t_j.end() ; ++itrBigger) // for every entry that is bigger than a_i
			{
				int index_j = itrBigger->second;
				double value_coef = (a_j[index_j]/b_j[index_j] - a_j[index_i]/b_j[index_i] )*b_0  ; 
				coefs.push_back(pair<int, double>(index_j, value_coef));
			}


			// prepare cut coefficients : add coef for i such that only one of a_i and b_i is zero
			for (int i = 0; i < DataModel::n; i++)
			{
				if (b_j[i] == 0 && a_j[i] > SMALL_DECIMAL)          // if only one of a_i and b_i is zero
				{
					coefs.push_back(pair<int, double>( i, a_j[i]));
				}
			}
			
			// prepare cut coefficients : add coef for z
			double newRhs = a_0 - b_0*a_j[index_i]/b_j[index_i];

			coefs.push_back(pair<int, double>(index_z, newRhs));
			CUT cut(coefs, newRhs);
			if (testCutValidity(cut, nodex) < -0.5)
			{
				callbackAddACut( cbdata, wherefrom, cut);
			}
		}

		delete  []t_rmatind;
		delete  []t_rmatval;
		delete t_rmatbeg;

	}
	free_and_null((char**) &cur_rowname);
	free_and_null((char**) &cur_rownamestore);
	delete []nodex;

	delete b_j;
	return 1;
}


// int OptModel::prepareCallbackCutCoefsIntersection(CPXENVptr& env,  CPXLPptr & nodelp_p , int &nzcnt, double &rhs, int &sense, int * &cutind, double * &cutval)
// {
// 
// 
// 	multimap<double, int> map_fracs = getAllFractionalZ(); // map_fracs stores the row index! not column index, ranked by closeness to 0.5
// 	if (map_fracs.size() <DataModel::k+1 - num_one_at_curr_node)
// 	{
// 		cout <<endl<< "num fracs: "<<map_fracs.size()<<"  num fixed 1: "<< num_one_at_curr_node <<"  :*********   Not enough fractionals  **********"<<endl;
// 		return 1;
// 	}
// 
// 	int cnt_frac = 0;
// 	vector<int> idx_frac_z ;
// 	for(multimap<double, int>::iterator itr = map_fracs.begin(); itr != map_fracs.end(); ++itr)
// 	{
// 		if (itr->first <= 0.5- SMALL_DECIMAL /*&& cnt_frac <= K+5*/)
// 			//			if (itr->first <= 1- INT_TOLERANCE  && itr->first >=  INT_TOLERANCE /*&& cnt_frac <= K+5*/)
// 		{
// 			idx_frac_z.push_back(itr->second);
// 			cnt_frac ++;
// 		}
// 		else if (cnt_frac <= DataModel::k+1)
// 		{
// 			idx_frac_z.push_back(itr->second);
// 		}
// 	}
// 
// 	vector<CUT> cuts;
//     cuts = selectCutsFromDiffFracSet(env, nodelp_p, idx_frac_z,  2);
// 
// 	for (int i = 0; i < cuts.size(); i ++)
// 	{
// 		CUT_COEFICIENTS cut = cuts[i].first;
// 		int rowNum = 0;
// 		int    status = 0;
// 		int totaNumCols = cut.size() ;
// 
// 		/* prepare row data*/
// 
// 		cutind = (int * )   malloc (totaNumCols * sizeof(int));
// 		cutval = (double *) malloc (totaNumCols * sizeof(double));
// 
// 		nzcnt = totaNumCols;
// 
// 		for (int j =0; j <totaNumCols; j ++)
// 		{
// 			PAIR_INDEX_VALUE entryPair = cut[j];
// 			cutind[j] = entryPair.first;
// 			cutval[j] = entryPair.second;
// 		}
// 		rhs = cuts[i].second;
// 		sense = 'G';
// 
// 		free_and_null((char**) &cutind);
// 		free_and_null((char**) &cutval);
// 	}
// 	return 0;
// }



int OptModel::setCPXDoubleParameter(int whichparam, double newvalue)
{
	return CPXsetdblparam(cpxenv, whichparam, newvalue); 
}

multimap<double, int> OptModel::getAllFractionalZ()
{

// 	CPXwriteprob(env, lp, "rootnodelp.lp",NULL);
// 	printTableau("testTableau.csv", env, lp);
	int surplus; 
	int num_cols = CPXgetnumcols(cpxenv, cpxlp);
	int status = CPXgetcolname (cpxenv, cpxlp, NULL, NULL, 0, &surplus, 0,num_cols-1);

	multimap<double, int> map_fracs;
	/***************get head ****************************/
	int num_rows = CPXgetnumrows(cpxenv,  cpxlp);

	int* head = new int[num_rows]; // head is the 0-th column of the tableau
	double *head_val = new double[num_rows];
	status = CPXgetbhead (cpxenv, cpxlp, head, head_val);
	for (int i = 0; i < num_rows; i++)
	{
		int index = head[i];
		double val = head_val[i];
		vector<int>::iterator itr = find(idx_Z.begin(), idx_Z.end(), index);
		if (itr == idx_Z.end())
			continue;
		if( val >= -SMALL_DECIMAL && val <= SMALL_DECIMAL )// degenerate
			continue;
		map_fracs.insert(make_pair(fabs(val-0.5), i));
// 		map_fracs.insert(make_pair(fabs(1-val), i));
	}
	delete []head;
	delete []head_val;
	return map_fracs;
}



CUT OptModel::substituteSlackVarInCut(CPXENVptr &env,  CPXLPptr & lp, CUT_COEFICIENTS& cutCoefs)
{
	CUT_COEFICIENTS slackCoefs;
	for (int i = 0; i < cutCoefs.size(); i++)// pick up slack vars
	{
		if (cutCoefs[i].first < 0 )
		{
			slackCoefs.push_back(cutCoefs[i]);
		}
	}

	if (slackCoefs.size() >0)
	{
		int num_cols = CPXgetnumcols (env, lp) ;
		int num_rows = CPXgetnumrows(env, lp);
		double* original_var_coef = new double[num_cols];
		memset(original_var_coef, 0 , sizeof(double)*num_cols);
		for (int i = 0; i < cutCoefs.size(); i++)// put coefficients of original vars into an array
		{
			if (cutCoefs[i].first >= 0)
			{
				original_var_coef[cutCoefs[i].first] = cutCoefs[i].second;
			}
		}

		// the following variables are for CPXgetrows CPXgetrhs
		int nzcnt;
		int * rmatbeg = new int(0);
		int * rmatind = new int[num_cols]();
		double * rmatval = new double[num_cols]();
		int rmatspace = num_cols;
		int surplus; 
		double rhs;

		double newCutRHS = 1;
		for (int i = 0; i < slackCoefs.size(); i ++) // for every slack vars in the cut
		{
			int row_index =  - slackCoefs[i].first -1;
			double slackCoef = slackCoefs[i].second;
			CPXgetrows(env, lp, &nzcnt, rmatbeg, rmatind, rmatval, rmatspace, &surplus, row_index, row_index);			// get the row corresponding to the slack var
			CPXgetrhs(env, lp, &rhs, row_index, row_index);
			char sense ;
			CPXgetsense(env, lp, &sense, row_index, row_index);

			if (sense == 'L' || sense == 'E')
			{
				for (int j = 0; j < nzcnt; j++)
				{
					original_var_coef[rmatind[j]] += - slackCoef*rmatval[j];
				}
				newCutRHS += - slackCoef*rhs;
			}
			else
			{
				for (int j = 0; j < nzcnt; j++)
				{
					original_var_coef[rmatind[j]] += slackCoef*rmatval[j];
				}
				newCutRHS += slackCoef*rhs;
			}
		}

		// put substituted coefficients into new cut
		CUT_COEFICIENTS newCutCoefs;
		for (int  i = 0; i < num_cols; i++)
		{
			if ( fabs( original_var_coef[i] ) < SMALL_DECIMAL )// skip zero entries
				continue;

			newCutCoefs.push_back(make_pair(i, original_var_coef[i]));
		}

		delete []original_var_coef;
		delete rmatbeg;
		delete []rmatind;
		delete []rmatval;
		return CUT(newCutCoefs, newCutRHS);
	}
	else
		return CUT(cutCoefs,1);

}

char ** OptModel::getColumnNames(CPXCENVptr &env , CPXCLPptr & lp)
{
	if (env == NULL)
		env = cpxenv;
	if (lp == NULL)
		lp = cpxlp;

	int num_cols  = CPXgetnumcols(env,  lp) ; // need to know why
	int surplus;
	char          **cur_colname;
	int status = CPXgetcolname (env, lp, NULL, NULL, 0, &surplus, 0,num_cols-1);
	int cur_colnamespace = - surplus;
	if ( cur_colnamespace > 0 ) {
		cur_colname      = (char **) malloc (sizeof(char *)*(num_cols));
		char          *cur_colnamestore = (char *)  malloc (cur_colnamespace);
		if ( cur_colname      == NULL ||
			cur_colnamestore == NULL   ) {
				fprintf (stderr, "Failed to get memory for column names.\n");
				status = -1;
				cout<< endl <<" error in printTableau"<<endl;
				return 0;
		}
		status = CPXgetcolname (env, lp, cur_colname, cur_colnamestore, 
			cur_colnamespace, &surplus, 0, num_cols-1);
 		free_and_null((char **)&cur_colname);
		free_and_null((char **)&cur_colnamestore);

	}
	return cur_colname;
}


bool OptModel::CallbackIsHeuristicSolution(CPXENVptr &env, CPXLPptr& lp, void * cbdata, int wherefrom, OptModel* model)
{
	int num_cols = CPXgetnumcols(env, lp);
	double *lb = new double[num_cols];
	double *ub = new double[num_cols];

	CPXgetcallbacknodelb(env, cbdata, wherefrom, lb, 0, num_cols-1);
	CPXgetcallbacknodeub(env, cbdata, wherefrom, ub, 0, num_cols-1);

	int z_size = model->idx_Z.size();
	for (int i = 0; i < z_size; i ++)
	{
		int index  = model->idx_Z[i];
		if (lb[index] != ub[index])
		{
			delete []lb;
			delete []ub;
			return false;
		}
	}
	delete []lb;
	delete []ub;
	return true;

}


bool OptModel::addCallbackSeperateCuts(CPXCENVptr &env, void * cbdata, int wherefrom, OptModel* model)
{
	CPXCLPptr origlp;
	int status = CPXgetcallbacklp (env, cbdata, wherefrom, &origlp);

	int cols = CPXgetnumcols(env, origlp);
	double* nodex= new double[cols];
	status = CPXgetcallbacknodex (env, cbdata, wherefrom, nodex, 0, cols-1);

	vector<int> inds_violated = checkViolatedSenarioByX(nodex);

	if (inds_violated.empty() )
	{
		delete []nodex;
		return false;
	}

	for (int i = 0; i < inds_violated.size(); i++)
	{
		int ind_violated = inds_violated[i];
		vector<CUT> cuts;// = genertaeCMIXcuts(ind_violated, nodex);
		if (cuts.empty())
		{
			continue;
		}
		else
		{
			// add cuts to master model
			for (int i = 0; i < cuts.size(); i++)
			{
				addCutToMasterModel(cuts[i], env, cbdata, wherefrom);
			}

			delete []nodex;

			return true;	

		}
	}

	return false;

}
				

// int OptModel::addSeperateCuts()
// {
// 	int cols = CPXgetnumcols(cpxenv, cpxlp);
// 	double* nodex= new double[cols];
// 	int status = CPXgetx (cpxenv, cpxlp, nodex, 0, cols-1);
// 
// 	vector<int> inds_violated = checkViolatedSenarioByX(nodex);
// 
// 	if (inds_violated.empty() )
// 	{
// 		delete []nodex;
// 		return false;
// 	}
// 
// 	for (int i = 0; i <=  inds_violated.size(); i++)// -1 means cost vector
// 	{
// 		int ind_violated =  i == inds_violated.size() ? -1 : inds_violated[i];
// 		vector<CUT> cuts = genertaeCMIXcuts(ind_violated, nodex);
// 		if (cuts.empty())
// 		{
// 			continue;
// 		}
// 		else
// 		{
// 			// add cuts to master model
// 			for (int i = 0; i < cuts.size(); i++)
// 			{
// 				string str = "MIX"+ToString(num_MIX_rootnode);
// 				addAConstToCPXModel(cuts[i],str);
// 				num_MIX_rootnode++;
// // 				writeModel2File("masterModel.lp");
// 			}
// 
// 			delete []nodex;
// 			return cuts.size();	
// 
// 		}
// 	}
// 
// 	return -1;
// 
// }




int OptModel::findSceNumByRowNum(int rowNum)
{
	map< string, int> rNames;
	int surplus;
	string rowName;

	int status = CPXgetrowname (cpxenv, cpxlp, NULL, NULL, 0, &surplus, rowNum,rowNum);
	int cur_rownamespace = - surplus;
	if ( cur_rownamespace > 0 ) 
	{
		char          **cur_rowname = (char **) malloc (sizeof(char *)*(1));
		char          *cur_rownamestore = (char *)  malloc (cur_rownamespace);
		if ( cur_rowname      == NULL ||
			cur_rownamestore == NULL   ) {
				fprintf (stderr, "Failed to get memory for column names.\n");
				status = -1;
				cout<< endl <<" error in printTableau"<<endl;
				return -1;
		}
		status =	 (cpxenv, cpxlp, cur_rowname, cur_rownamestore, 
			cur_rownamespace, &surplus, rowNum, rowNum);
		// map name with column index
		rowName = cur_rowname[0];
		free_and_null( (char **)&cur_rowname);
		free_and_null((char **)&cur_rownamestore);

	}
	size_t t =rowName.find_first_of("0123456789");
	string strSceNum = rowName.substr(t);
	return atoi(strSceNum.c_str());
}


int OptModel::findRowNumByRowName(string& rowName)
{
	int rowNum;
	int num_rows = getNumRows();
	int surplus;

	int status = CPXgetrowname (cpxenv, cpxlp, NULL, NULL, 0, &surplus, 0,num_rows-1);
	int cur_rownamespace = - surplus;
	if ( cur_rownamespace > 0 ) 
	{
		char          **cur_rowname = (char **) malloc (sizeof(char *)*(num_rows));
		char          *cur_rownamestore = (char *)  malloc (cur_rownamespace);
		if ( cur_rowname      == NULL ||
			cur_rownamestore == NULL   ) {
				fprintf (stderr, "Failed to get memory for column names.\n");
				status = -1;
				cout<< endl <<" error in reading column names"<<endl;
				return -1;
		}
		status = CPXgetrowname (cpxenv, cpxlp, cur_rowname, cur_rownamestore, 
			cur_rownamespace, &surplus, 0, num_rows-1);
		// map name with column index

		for (int r = 0; r <num_rows ; r ++)
		{
			string name = cur_rowname[r];
			if ( name == rowName)
			{
				rowNum = r;
				break;
			}
		}
		free_and_null( (char **)&cur_rowname);
		free_and_null((char **)&cur_rownamestore);
	}
	return rowNum;
}



void OptModel::mapScenarioRownum()
{

	int num_rows = getNumRows();
	int surplus;
	int num_cols = CPXgetnumcols(cpxenv, cpxlp);

	int status = CPXgetrowname (cpxenv, cpxlp, NULL, NULL, 0, &surplus, 0,num_rows-1);
	int cur_rownamespace = - surplus;
	if ( cur_rownamespace > 0 ) 
	{
		char          **cur_rowname = (char **) malloc (sizeof(char *)*(num_rows));
		char          *cur_rownamestore = (char *)  malloc (cur_rownamespace);
		if ( cur_rowname      == NULL ||
			cur_rownamestore == NULL   ) {
				fprintf (stderr, "Failed to get memory for column names.\n");
				status = -1;
				cout<< endl <<" error in reading column names"<<endl;
				return;
		}
		status = CPXgetrowname (cpxenv, cpxlp, cur_rowname, cur_rownamestore, 
			cur_rownamespace, &surplus, 0, num_rows-1);
		// map name with column index

		for (int r = 0; r <num_rows ; r ++)
		{
			string rowName = cur_rowname[r];
			string strSceNum;
			size_t t =rowName.find_first_of("0123456789");
			if (t != string::npos)
			{
				size_t t2 = rowName.find_first_not_of("0123456789");
				strSceNum = rowName.substr(t, t2-t);
				rowName.erase(t);
			}
			if (rowName != "SN")
				continue;
			
			int sceNum = atoi(strSceNum.c_str());
			map_scenario_rownum.insert(make_pair(sceNum, r));
		}
		free_and_null( (char **)&cur_rowname);
		free_and_null((char **)&cur_rownamestore);
	}
}



void OptModel::mapZIdxRowNum()
{

	int num_rows = getNumRows();
	int surplus;
	int num_cols = CPXgetnumcols(cpxenv, cpxlp);

	int status = CPXgetrowname (cpxenv, cpxlp, NULL, NULL, 0, &surplus, 0,num_rows-1);
	int cur_rownamespace = - surplus;
	if ( cur_rownamespace > 0 ) 
	{
		char          **cur_rowname = (char **) malloc (sizeof(char *)*(num_rows));
		char          *cur_rownamestore = (char *)  malloc (cur_rownamespace);
		if ( cur_rowname      == NULL ||
			cur_rownamestore == NULL   ) {
				fprintf (stderr, "Failed to get memory for column names.\n");
				status = -1;
				cout<< endl <<" error in reading column names"<<endl;
				return;
		}
		status = CPXgetrowname (cpxenv, cpxlp, cur_rowname, cur_rownamestore, 
			cur_rownamespace, &surplus, 0, num_rows-1);
		// map name with column index

		for (int r = 0; r <num_rows ; r ++)
		{
			string rowName = cur_rowname[r];
			string cpyRowName = rowName;
			size_t t =rowName.find_first_of("0123456789");
			if(t== string::npos)
				continue;
			rowName.erase(t);
			if (rowName != "SN")
				continue;

			int t_nzcnt;
			int * t_rmatbeg = new int(0);
			int * t_rmatind = new int[num_cols]();
			double * t_rmatval = new double[num_cols]();
			int t_rmatspace = num_cols;
			int t_surplus; 

			CPXgetrows(cpxenv, cpxlp, &t_nzcnt, t_rmatbeg, t_rmatind, t_rmatval, t_rmatspace, &t_surplus, r, r);  // get row coefficients from original formulation
			int counter = 0;
			for (int i = 0; i < t_nzcnt; i ++)
			{
				if ( t_rmatind[i] < DataModel::n) 
					continue; 
//				cpyRowName
				map_z_idx_row_num.insert( make_pair (t_rmatind[i], r ) );
				map_row_num_z_idx.insert( make_pair ( r, t_rmatind[i] ) );
				break;
			}
			delete []t_rmatbeg;
			delete []t_rmatind;
			delete []t_rmatval;

		}
		free_and_null( (char **)&cur_rowname);
		free_and_null((char **)&cur_rownamestore);
	}
}

bool OptModel::addCutToMasterModel(CUT& cut, CPXCENVptr &env, void * cbdata, int wherefrom)
{
	CUT_COEFICIENTS vec_coefs = cut.first;
	int nzcnt = vec_coefs.size();
	double * coefs_p = new double[nzcnt];
	int * idx_p = new int[nzcnt];
	for (int i = 0; i < nzcnt; i ++)
	{
		coefs_p[i] = vec_coefs[i].second;
		idx_p[i] = vec_coefs[i].first;
	}
	char sense = 'G';
	double rhs = cut.second;

	return CPXcutcallbackadd(env, cbdata, wherefrom, nzcnt, rhs, sense, idx_p, coefs_p, false);

}

bool OptModel::addAConstToCPXModel(CUT& cut, string& constName, int cstType)
{
	CUT_COEFICIENTS vec_coefs = cut.first;
	int nzcnt = vec_coefs.size();
	int status =0;
	/* prepare row data*/

	int* rmatind = (int * )   malloc (nzcnt * sizeof(int));
	double* rmatval = (double *) malloc (nzcnt * sizeof(double));

	// fill in the entries for a row
	int ind = 0;
	for (int j =0; j <vec_coefs.size(); j ++)
	{
		pair<int,double> entryPair = vec_coefs[j];
		if(entryPair.second == 0) // skip zero elements
			continue;
		rmatind[ind] = entryPair.first;
		rmatval[ind] = entryPair.second;
		ind ++;
	}

	int    rmatbeg[1];
	double rhs[1];
	char   sense[1];
	rhs[0] = cut.second;
	sense[0] = 'G';
	rmatbeg[0] = 0;

 	char *row_name[1];
	char str[20];
	strcpy(str, constName.c_str());
 	row_name[0]= str;

	if (cstType == 0)// add as constraint
	{
		status = CPXaddrows (cpxenv, cpxlp, 0, 1, ind , rhs,	sense, rmatbeg, rmatind, rmatval, NULL, row_name);
	}
	else if (cstType == 1) // add as user cut
	{
		status = CPXaddusercuts (cpxenv, cpxlp, 1, ind, rhs,	sense, rmatbeg, rmatind, rmatval, row_name);
	}

	if ( status ) {
		fprintf (stderr, "Could not add new rows.\n");
		return 0;
	}

	// 	CPXwriteprob(cpxenv, lp, "cutModel.lp", NULL);

	free_and_null((char**) &rmatind);
	free_and_null((char**) &rmatval);
	return 1;
}


bool OptModel::callbackAddACut(void * cbdata, int wherefrom, CUT& cut)
{
	CUT_COEFICIENTS vec_coefs = cut.first;
	int nzcnt = vec_coefs.size();
	int status =0;
	/* prepare row data*/

	int* cutind = (int * )   malloc (nzcnt * sizeof(int));
	double* cutval = (double *) malloc (nzcnt * sizeof(double));

	// fill in the entries for a row
	int ind = 0;
	for (int j =0; j <vec_coefs.size(); j ++)
	{
		pair<int,double> entryPair = vec_coefs[j];
		if(entryPair.second == 0) // skip zero elements
			continue;
		cutind[ind] = entryPair.first;
		cutval[ind] = entryPair.second;
		ind ++;
	}

	double rhs = cut.second;
	int   sense = 'G';

	status = CPXcutcallbackaddlocal(cpxenv, cbdata, wherefrom, ind, rhs, sense, cutind, cutval);

	if ( status ) {
		fprintf (stderr, "Could not add new rows.\n");
		return 0;
	}

	// 	CPXwriteprob(cpxenv, lp, "cutModel.lp", NULL);

	free_and_null((char**) &cutind);
	free_and_null((char**) &cutval);
	return 1;
}



double OptModel::testCutValidity(CUT& cut, double* nodex)
{
	CUT_COEFICIENTS& coefs = cut.first;
	double rhs = cut.second;
	double lhs = 0;
	for (int i = 0; i < coefs.size(); i ++)
	{
		double val = coefs[i].second;
		int index = coefs[i].first;
		lhs += val* nodex[index];
	}
	return lhs - rhs;
}

vector<int> OptModel::checkViolatedSenarioByX(double* nodex)
{
	// sort z values in solution and start from the smallest z, possible 0
	multimap<double, int> zValues;
	for (int i = 0; i < idx_Z.size(); i++)
	{
		int ind = idx_Z[i];
		zValues.insert(make_pair(nodex[ind], i));
	}

	// check if scenario z<1 is violated
	vector<int> violatedScenarios;
	int ind_violated = -1;
	for (multimap<double, int>::iterator itr = zValues.begin(); itr !=  zValues.end(); ++itr)
	{
		double val = itr->first;
		if (val >= 1- SMALL_DECIMAL) // only look at z that is smaller than 1
			break;

		// test if scenario z is violated
		int ind = itr->second;
		double violation = IsXValidForScenarioZ(ind, nodex, val);

		if (violation >= -SMALL_DECIMAL)
			continue;
		
		ind_violated =ind;
		violatedScenarios.push_back(ind_violated);

	}

	return violatedScenarios;
}

// to do , should return the most violated ones
vector<int> OptModel::checkViolatedConstrainsByX(string& cstType, double* nodex)
{

	vector< pair<string, int> > rowNameNum = readRowNamesFromModel();
	// check if scenario z<1 is violated
	multimap<double,int> sortedViolation;
	vector<int> violatedScenarios;
	int ind_violated = -1;
	for (int i = 0 ; i <  rowNameNum.size(); i++)
	{
		string rowName = rowNameNum[i].first;
		size_t t =rowName.find_first_of("0123456789");
		if (t != string::npos)
		{
			rowName.erase(t);
		}

		if ((rowName!= "SN" && rowName!= "MIR" && rowName != cstType &&  cstType != "ALL") || rowName == "MIX" ||  rowName == "CARD" || rowName == "c" )
			continue;

		int rowNum =rowNameNum[i].second;

		// test if this constrain is violated
		double violation = IsXValidForConstrain(rowNum, nodex);

		if (violation >= -SMALL_DECIMAL)
			continue;

		ind_violated =rowNum;
		violatedScenarios.push_back(ind_violated);

	}

	return violatedScenarios;
}


vector<int> OptModel::checkSortedViolatedConstrainsByX(string& cstType, double* nodex)
{

	vector< pair<string, int> > rowNameNum = readRowNamesFromModel();
	// check if scenario z<1 is violated
	multimap<double,int> sortedViolation;
	int ind_violated = -1;
	for (int i = 0 ; i <  rowNameNum.size(); i++)
	{
		string rowName = rowNameNum[i].first;
		size_t t =rowName.find_first_of("0123456789");
		if (t != string::npos)
		{
			rowName.erase(t);
		}

		if ((rowName!= "SN" && rowName!= "MIR" && rowName != cstType && cstType != "ALL") || rowName == "MIX" ||  rowName == "CARD" )
			continue;

		int rowNum =rowNameNum[i].second;

		// test if this constrain is violated
		double violation = IsXValidForConstrain(rowNum, nodex);

		if (violation >= -SMALL_DECIMAL)
			continue;

		ind_violated =rowNum;
		sortedViolation.insert(make_pair(violation, ind_violated));
	}
	vector<int> violatedScenarios;
	for (multimap<double, int>::iterator itr = sortedViolation.begin(); itr!= sortedViolation.end(); ++itr)
	{
		violatedScenarios.push_back(itr->second);
	}
	return violatedScenarios;
}

double OptModel::IsXValidForConstrain( int rowNum, double* x)
{
	int num_cols = CPXgetnumcols(cpxenv, cpxlp);
	int nzcnt;
	int * rmatbeg = new int(0);
	int * rmatind = new int[num_cols]();
	double * rmatval = new double[num_cols]();
	int rmatspace = num_cols;
	int surplus; 
	double rhs;

	CPXgetrows(cpxenv, cpxlp, &nzcnt, rmatbeg, rmatind, rmatval, rmatspace, &surplus, rowNum, rowNum);

	double lhs = 0; 
	for (int i = 0; i < nzcnt; i++)
	{
		int index = rmatind[i];
		if (index >= DataModel::n)
			continue;
		lhs += x[index]*rmatval[i];
	}

	CPXgetrhs(cpxenv, cpxlp, &rhs, rowNum, rowNum);
	return lhs - rhs;
}



double OptModel::IsXValidForScenarioZ(int ind_scenario, double* x, double z)
{
	 VEC_MTX_ROW row =  dm_p->Arows[ind_scenario];
	
	 double LHS =0; 
	 for (int i = 0; i < idx_X.size(); i++)
	 {
		 int ind = idx_X[i];
		 LHS += row[ind].second*x[ind];
	 }
	//LHS += z*dm_p->rhs[ind_scenario]; // now we only check if ax>=b, instead of ax+zb >= b

	return LHS - dm_p->rhs[ind_scenario];
}


double OptModel::IsXValidForScenarioZ(int ind_scenario, double* x)
{
	VEC_MTX_ROW row =  dm_p->Arows[ind_scenario];

	double LHS =0; 
	for (int i = 0; i < idx_X.size(); i++)
	{
		int ind = idx_X[i];
		LHS += row[ind].second*x[ind];
	}
	//LHS += z*dm_p->rhs[ind_scenario]; // now we only check if ax>=b, instead of ax+zb >= b

	return LHS - dm_p->rhs[ind_scenario];
}

double OptModel::IsXValidForScenarioZNormalized(int ind_scenario, double* x)
{
	VEC_MTX_ROW row =  dm_p->Arows[ind_scenario];

	double LHS =0; 
	for (int i = 0; i < idx_X.size(); i++)
	{
		int ind = idx_X[i];
		double coeff = (float)((int)(row[ind].second/dm_p->rhs[ind_scenario]*100+0.5))/100.0;
		LHS += x[ind]*coeff;//row[ind].second*x[ind]/dm_p->rhs[ind_scenario];
	}
	//LHS += z*dm_p->rhs[ind_scenario]; // now we only check if ax>=b, instead of ax+zb >= b

	return LHS - 1;
}


bool OptModel::IsXFeasibleSolution(double* x)
{
	int neg_cnt = 0;
	for (int i = 0; i < DataModel::m; i++)
	{
		if (IsXValidForScenarioZ(i, x) < -SMALL_DECIMAL )
		{
			neg_cnt ++;
		}
	}
	return neg_cnt <= DataModel::k;
}

bool OptModel::IsCurrentSolutionFeasible()
{
	int num_cols = getNumCols();
	double *x = new double[num_cols];
	CPXgetx(cpxenv, cpxlp, x, 0, num_cols -1);
	int neg_cnt = 0;

	for (int i = 0; i < DataModel::m; i++)
	{
		if (IsXValidForScenarioZNormalized(i, x) < -SMALL_DECIMAL )
		{
			neg_cnt ++;
		}
	}
	return neg_cnt <= DataModel::k;
}

double OptModel::getDoubleSolQuality(int paraName)
{
	double quality_p;
	CPXgetdblquality(cpxenv, cpxlp, &quality_p, paraName);
	return quality_p;
}



int OptModel::updateCutPool(string cutType)
{
	if (cutType == string("IC"))
	{
		ICCutPool.updateCutPool(cpxenv, cpxlp);
	}
	else if (cutType == string("MIX"))
	{
		MixCutPool.updateCutPool(cpxenv, cpxlp);
	}
	else if( cutType == string("MIR"))
		MIRCutPool.updateCutPool(cpxenv, cpxlp);
	else if ( cutType == string("CH"))
	{
		CHCutPool.updateCutPool(cpxenv, cpxlp);
	}
	else if (cutType == string("ALL"))
	{
		ICCutPool.updateCutPool(cpxenv, cpxlp);
		MixCutPool.updateCutPool(cpxenv, cpxlp);
		MIRCutPool.updateCutPool(cpxenv, cpxlp);
		CHCutPool.updateCutPool(cpxenv, cpxlp);
	}
	else
	{
		GeneralCutPool.updateCutPool(cpxenv, cpxlp);
	}

	return 0;
}


int OptModel::storeCuts(string cutType)
{
	if (cutType == string("IC"))
	{
		ICCutPool.extractCutsFromModel(cpxenv, cpxlp);
	}
	else if (cutType == string("MIX"))
	{
		MixCutPool.extractCutsFromModel(cpxenv, cpxlp);
	}
	else if( cutType == string("MIR"))
		MIRCutPool.extractCutsFromModel(cpxenv, cpxlp);
	else if ( cutType == string("CH"))
	{
		CHCutPool.extractCutsFromModel(cpxenv, cpxlp);
	}
	else if (cutType == string("ALL"))
	{
		ICCutPool.extractCutsFromModel(cpxenv, cpxlp);
		MixCutPool.extractCutsFromModel(cpxenv, cpxlp);
// 		MIRCutPool.extractCutsFromModel(cpxenv, cpxlp);
		CHCutPool.extractCutsFromModel(cpxenv, cpxlp);
	}
	else
	{
		GeneralCutPool.setCutIdentifier(cutType);
		GeneralCutPool.extractCutsFromModel(cpxenv,cpxlp);
	}

	return 0;
}


int OptModel::addAllMIRCuts()
{
	for (int i = 0; i < DataModel::m ; i++)
	{
		VEC_MTX_ROW row = dm_p->Arows[i];
		multimap<double, PAIR_VALUE_INDEX> sortedAi; 
		for (int j = 0; j < DataModel::n; j++) // find the smallest aij
		{
			sortedAi.insert(make_pair( row[j].second, row[j] ));
		}
		PAIR_VALUE_INDEX min_aij = sortedAi.begin()->second;
		
		int index_j = min_aij.second;
		double a_j = min_aij.first;

		CUT_COEFICIENTS coefs;
		CUT cut;
		int index_zi = idx_Z[i];
		double t = dm_p->rhs[i] - a_j;
		if (t <=0) // then this constrain is always satisfied, set z_i = 0
		{	
			coefs.push_back(pair<int, double>(index_zi, -1));
			cut = CUT(coefs , 0);
		}
		else
		{
			for (int k = 0 ; k < DataModel::n; k++) // replace xi with 1-sum{x}
			{
				if (k == index_j)
					continue;

				double a_i = row[k].second;
				int index_i = idx_X[row[k].first];
				coefs.push_back(pair<int, double>(index_i, a_i - a_j));
			}
			coefs.push_back(pair<int, double>(index_zi, t));
			cut = CUT(coefs, t);
		}
		string str = string("MIR")+ToString(i);
		addAConstToCPXModel(cut, str);
	}
	return 0;
}


int OptModel::AddConvexHullToModel(int cstType)
{
	int num_rows = CPXgetnumrows(cpxenv,  cpxlp);
	int num_cols = CPXgetnumcols(cpxenv, cpxlp);
	int surplus;
	vector<string> rowNames;

	int status = CPXgetrowname (cpxenv, cpxlp, NULL, NULL, 0, &surplus, 0,num_rows-1);
	int cur_rownamespace = - surplus;
	if ( cur_rownamespace > 0 ) 
	{
		char          **cur_rowname = (char **) malloc (sizeof(char *)*(num_rows));
		char          *cur_rownamestore = (char *)  malloc (cur_rownamespace);
		if ( cur_rowname      == NULL ||
			cur_rownamestore == NULL   ) {
				fprintf (stderr, "Failed to get memory for column names.\n");
				status = -1;
				cout<< endl <<" error in printTableau"<<endl;
				return 0;
		}
		status = CPXgetrowname (cpxenv, cpxlp, cur_rowname, cur_rownamestore, 
			cur_rownamespace, &surplus, 0, num_rows-1);
		// map name with column index
		for (int r = 0; r <num_rows ; r ++)
		{
			string rowName = cur_rowname[r];
			if (rowName.find("MIR") == string::npos)// only look at scenario constrains
				continue; 

			// add Convex hull
			int nzcnt;
			int * rmatbeg = new int(0);
			int * rmatind = new int[num_cols]();
			double * rmatval = new double[num_cols]();
			int rmatspace = num_cols;
			int surplus; 
			double rhs;

			CPXgetrows(cpxenv, cpxlp, &nzcnt, rmatbeg, rmatind, rmatval, rmatspace, &surplus, r, r);// for every row, add the convex hull
			CPXgetrhs(cpxenv, cpxlp, &rhs, r, r);
			if (rhs <= SMALL_DECIMAL)
				continue;

			multimap<double, int > sorted_a_j;
			for(int j = 0; j < nzcnt; j ++)
			{
				sorted_a_j.insert(make_pair(rmatval[j], j));
			}

			//for the first k a_j that is strictly smaller than rhs
			for (multimap<double, int>::iterator itr = sorted_a_j.begin(); itr != sorted_a_j.end(); ++itr)// for every a_i < a_0, generate one facet
			{
				if (itr->first >= rhs) 
					continue;
				// prepare cut coefficients
				CUT_COEFICIENTS coefs;
				double a_i = itr->first;
				int index_a_i = itr->second;
				for (int j = 0 ; j <  nzcnt; j++) // for every entry in the row from cplex model
				{
					if ( rmatval[j] <= a_i || j == index_a_i)// don't consider 
						continue;
					double value_coef = rmatval[j] - a_i; 
					int index_a_j = rmatind[j];
					coefs.push_back(pair<int, double>(index_a_j, value_coef));
				}
				double newRhs = rhs - a_i;
				CUT cut(coefs, newRhs);
				string cstName = string("CH") +string("r") + ToString(r) + string("i") +ToString(index_a_i);
				addAConstToCPXModel(cut, cstName, cstType);
			}

			delete  []rmatind;
			delete  []rmatval;
		}
		free_and_null( (char **)&cur_rowname);
		free_and_null( (char **)&cur_rownamestore);

	}
	return 0;
};



vector< pair<string, int> > OptModel::readRowNamesFromModel()
{
	vector < pair<string, int> > rNames;
	int num_rows = CPXgetnumrows(cpxenv,  cpxlp);
	int surplus;

	int status = CPXgetrowname (cpxenv, cpxlp, NULL, NULL, 0, &surplus, 0,num_rows-1);
	int cur_rownamespace = - surplus;
	if ( cur_rownamespace > 0 ) 
	{
		char          **cur_rowname = (char **) malloc (sizeof(char *)*(num_rows));
		char          *cur_rownamestore = (char *)  malloc (cur_rownamespace);
		if ( cur_rowname      == NULL ||
			cur_rownamestore == NULL   ) {
				fprintf (stderr, "Failed to get memory for column names.\n");
				status = -1;
				cout<< endl <<" error in printTableau"<<endl;
				return rNames;
		}
		status = CPXgetrowname (cpxenv, cpxlp, cur_rowname, cur_rownamestore, 
			cur_rownamespace, &surplus, 0, num_rows-1);
		// map name with column index
		for (int r = 0; r <num_rows ; r ++)
		{
			string rowName = cur_rowname[r];
			rNames.push_back(make_pair(rowName,r));
		}
		free_and_null( (char **)&cur_rowname);
 		free_and_null((char **)&cur_rownamestore);

	}
	return rNames;
};


map< string, int> OptModel::readRowNamesMapFromModel()
{
	map< string, int> rNames;
	int num_rows = CPXgetnumrows(cpxenv,  cpxlp);
	int surplus;

	int status = CPXgetrowname (cpxenv, cpxlp, NULL, NULL, 0, &surplus, 0,num_rows-1);
	int cur_rownamespace = - surplus;
	if ( cur_rownamespace > 0 ) 
	{
		char          **cur_rowname = (char **) malloc (sizeof(char *)*(num_rows));
		char          *cur_rownamestore = (char *)  malloc (cur_rownamespace);
		if ( cur_rowname      == NULL ||
			cur_rownamestore == NULL   ) {
				fprintf (stderr, "Failed to get memory for column names.\n");
				status = -1;
				cout<< endl <<" error in printTableau"<<endl;
				return rNames;
		}
		status = CPXgetrowname (cpxenv, cpxlp, cur_rowname, cur_rownamestore, 
			cur_rownamespace, &surplus, 0, num_rows-1);
		// map name with column index
		for (int r = 0; r <num_rows ; r ++)
		{
			string rowName = cur_rowname[r];
			rNames.insert(make_pair(rowName,r));
		}
		free_and_null( (char **)&cur_rowname);
		free_and_null((char **)&cur_rownamestore);

	}
	return rNames;
};


CUT_COEFICIENTS OptModel::readAlphaFromConstrain(int rowNum)// returned vector includes 0 coefficients
{

	if(rowNum == -1 )
	{
		CUT_COEFICIENTS vec_alpha;
		for (int i = 0 ; i < DataModel::n ; i++)
		{
			vec_alpha.push_back(make_pair(dm_p->costVec[i].second, dm_p->costVec[i].first));
		}
		return vec_alpha;
	}
	else
	{
		int num_cols = CPXgetnumcols(cpxenv, cpxlp);
		int nzcnt;
		int * rmatbeg = new int(0);
		int * rmatind = new int[num_cols]();
		double * rmatval = new double[num_cols]();
		int rmatspace = num_cols;
		int surplus; 
		double rhs;

		CPXgetrows(cpxenv, cpxlp, &nzcnt, rmatbeg, rmatind, rmatval, rmatspace, &surplus, rowNum, rowNum);
		double *p_coefs = new double[DataModel::n]();
		memset(p_coefs,0 , sizeof(double)*DataModel::n);
		CUT_COEFICIENTS coefs;
		for (int i = 0; i < nzcnt; i++)
		{
			int index = rmatind[i];
			if ( index >= DataModel::n)
				continue;
			p_coefs[index] = rmatval[i];
		}
		for(int i = 0; i < DataModel::n; i++)
		{
			coefs.push_back(pair<int,double>(i,p_coefs[i]));
		}
		delete []p_coefs;
		return coefs;
	}
}

CUT OptModel::readScenarioFromModel(int rowNum)
{

	if(rowNum == -1 )
	{
		CUT_COEFICIENTS vec_alpha;
		for (int i = 0 ; i < DataModel::n ; i++)
		{
			vec_alpha.push_back(make_pair(dm_p->costVec[i].second, dm_p->costVec[i].first));
		}
		return CUT(vec_alpha,0);
	}
	else
	{
		int num_cols = CPXgetnumcols(cpxenv, cpxlp);
		int nzcnt;
		int * rmatbeg = new int(0);
		int * rmatind = new int[num_cols]();
		double * rmatval = new double[num_cols]();
		int rmatspace = num_cols;
		int surplus; 
		double rhs;
		
		CPXgetrows(cpxenv, cpxlp, &nzcnt, rmatbeg, rmatind, rmatval, rmatspace, &surplus, rowNum, rowNum);
		CPXgetrhs(cpxenv, cpxlp, &rhs, rowNum, rowNum);
		double *p_coefs = new double[DataModel::n]();
		memset(p_coefs,0 , sizeof(double)*DataModel::n);
		CUT_COEFICIENTS coefs;
		for (int i = 0; i < nzcnt; i++)
		{
			int index = rmatind[i];
			if ( index >= DataModel::n)
				continue;
			p_coefs[index] = rmatval[i];
		}
		for(int i = 0; i < DataModel::n; i++)
		{
			coefs.push_back(pair<int,double>(i,p_coefs[i]));
		}
		delete []p_coefs;
		return CUT(coefs,rhs);
	}


}



int OptModel::addCutPoolToModel(string poolName, int asConst)
{
	if (poolName == string("IC"))
	{
		for(map<string, CUT>::iterator itr = ICCutPool.cuts.begin(); itr != ICCutPool.cuts.end(); ++itr)
		{
			CUT cut = itr->second;
			string cstName = itr->first;
			addAConstToCPXModel(cut, cstName, asConst);
		}	
	}
	else if (poolName == string("MIX"))
	{
		for(map<string, CUT>::iterator itr = MixCutPool.cuts.begin(); itr != MixCutPool.cuts.end();++itr)
		{
			CUT cut = itr->second;
			string cstName = itr->first;
			addAConstToCPXModel(cut, cstName, asConst);
		}
	}
	else if( poolName == string("MIR"))
	{
		for(map<string, CUT>::iterator itr = MIRCutPool.cuts.begin(); itr != MIRCutPool.cuts.end();++itr)
		{
			CUT cut = itr->second;
			string cstName = itr->first;
			addAConstToCPXModel(cut, cstName, asConst);
		}	}
	else if ( poolName == string("CH"))
	{
		for(map<string, CUT>::iterator itr = CHCutPool.cuts.begin(); itr != CHCutPool.cuts.end();++itr)
		{
			CUT cut = itr->second;
			string cstName = itr->first;
			addAConstToCPXModel(cut, cstName, asConst);
		}	}
	else if (poolName == string("ALL"))
	{
		for(map<string, CUT>::iterator itr = ICCutPool.cuts.begin(); itr != ICCutPool.cuts.end();++itr)
		{
			CUT cut = itr->second;
			string cstName = itr->first;
			addAConstToCPXModel(cut, cstName, asConst);
		}
		for(map<string, CUT>::iterator itr = MixCutPool.cuts.begin(); itr != MixCutPool.cuts.end();++itr)
		{
			CUT cut = itr->second;
			string cstName = itr->first;
			addAConstToCPXModel(cut, cstName, asConst);
		}
		for(map<string, CUT>::iterator itr = CHCutPool.cuts.begin(); itr != CHCutPool.cuts.end(); ++itr)
		{
			CUT cut = itr->second;
			string cstName = itr->first;
			addAConstToCPXModel(cut, cstName, asConst);
		}
	}

	return 0;


}


int OptModel::removeInactiveCuts(string cutType, int numRoundsInactive)
{
	if (cutType == string("IC"))
	{
		ICCutPool.removeInactiveCuts(cpxenv, cpxlp, numRoundsInactive);
	}
	else if (cutType == string("MIX"))
	{
		MixCutPool.removeInactiveCuts(cpxenv, cpxlp, numRoundsInactive);
	}
	else if( cutType == string("MIR"))
		MIRCutPool.removeInactiveCuts(cpxenv, cpxlp, numRoundsInactive);
	else if ( cutType == string("CH"))
	{
		CHCutPool.removeInactiveCuts(cpxenv, cpxlp, numRoundsInactive);
	}
	else if (cutType == string("ALL"))
	{
		ICCutPool.removeInactiveCuts(cpxenv, cpxlp, numRoundsInactive);
		MixCutPool.removeInactiveCuts(cpxenv, cpxlp, numRoundsInactive);
// 		MIRCutPool.removeInactiveCuts(cpxenv, cpxlp, numRoundsInactive);
		CHCutPool.removeInactiveCuts(cpxenv, cpxlp, numRoundsInactive);
	}
	else
	{
		GeneralCutPool.removeInactiveCuts(cpxenv, cpxlp, numRoundsInactive);
	}

	return 0;
}



bool OptModel::Agg_generateAggregationModel()
{
	// calculate norm
	double** p = new double*[DataModel::m];
	for (int i = 0; i < DataModel::m; i++)
	{
		p[i] = new double[DataModel::m];
		memset(p[i], 0, sizeof(double)*DataModel::m);
	}

	for (int r = 0; r  < DataModel::m; r++)
	{
		CUT scenario1 = readScenarioFromModel(r);
		for (int c = 0; c< DataModel::m; c++)
		{
			if (c <= r )
			{
				p[r][c] = p[c][r];
				continue;
			}
			CUT scenario2 = readScenarioFromModel(c);
			p[r][c] = calculateDistance(scenario1,scenario2,1);
		}
	}
	

	int group_size = 2; int num_group = ceil(((double)DataModel::m)/group_size);
	vector<VEC_INT> groups;
	bool * taken_p = new bool[DataModel::m];
	memset(taken_p, false, sizeof(bool)*DataModel::m);
	for (int i = 0; i < DataModel::m; i++)
	{
		if(taken_p[i])
			continue;
		map<double, int> sortedDist;
		for (int j = i+1; j < DataModel::m; j++)
		{
			if (taken_p[j])
				continue;
			sortedDist.insert(make_pair(p[i][j], j));
		}
		int idx = sortedDist.begin()->second;
		taken_p[i] = true; taken_p[idx]= true;
		vector<int> group;
		group.push_back(i);
		group.push_back(idx);
		groups.push_back(group);
	}
	delete []taken_p;
	OptModel aggregationModel("aggregationModel");
	aggregationModel.Agg_populateAggregationModel(*dm_p, OptModel::MODEL_TYPE_MIP, 'L', false,	1);
	aggregationModel.solve();
	


	for (int i = 0; i < DataModel::m; i++)
	{
		delete p[i];
	}
	delete []p;

	return true;
}


double OptModel::calculateDistance(CUT &scenario1, CUT &scenario2, int method)
{
	double norm = 0;
	if (method == 2)
	{
		double sum_sq = 0;
		for (int i = 0; i < DataModel::n; i++)
		{
			double diff = ((scenario1.first)[i].second)/(scenario1.second) - ((scenario2.first)[i].second)/(scenario2.second);
			sum_sq += diff*diff;
		}
		return sqrt(sum_sq);
	}	
	else if (method == 1)
	{
		double sum = 0;
		for (int i = 0; i < DataModel::n; i++)
		{
			double diff = ((scenario1.first)[i].second)/(scenario1.second) - ((scenario2.first)[i].second)/(scenario2.second);
			sum += fabs(diff);
		}
		return sum;
	}

    return 0;
}



double OptModel::Agg_populateAggregationModel(DataModel& dm, OptModel::MODEL_TYPE modeltype, char cardinalitySign, bool resConst,	int method)
{
	setDataModel(&dm);
	idx_X.clear();
	idx_Z.clear();
	idx_S.clear();
	double *obj     = NULL;
	double *lb      = NULL;
	double *ub      = NULL;
	char   *ctype   = NULL;
	int    *rmatind = NULL;
	double *rmatval = NULL;
	char** colname = NULL;

	int    status = 0;
	int num_cols = DataModel::n + aggGroups.size();

	CPXchgobjsen (cpxenv, cpxlp, 1); /* Minimization problem */

	/* Allocate colcnt-sized arrays */

	obj     = (double *) malloc (num_cols * sizeof(double));
	lb      = (double *) malloc (num_cols* sizeof(double));
	ub      = (double *) malloc (num_cols * sizeof(double));
	ctype   = (char *)   malloc (num_cols * sizeof(char));
	colname = (char**) malloc(num_cols*sizeof(char*));

	if ( obj     == NULL ||
		lb      == NULL ||
		ub      == NULL ||
		ctype   == NULL ) {
			fprintf (stderr, "Could not allocate colcnt arrays\n");
			status = CPXERR_NO_MEMORY;
			return 0;
	}

	string* strColNames = new string[num_cols];
	for (int i =0; i < num_cols; i ++)
	{
		if (i < DataModel::n )
		{
			obj[i] = dm.costVec[i].first;                      // fill in obj coefficients
			idx_X.push_back(i);                        // store index for x variables
			strColNames[i] = "x"+ToString(i);
			colname[i] = new char[8];
			strcpy(colname[i], strColNames[i].c_str()); // make name for x
			ctype[i] = 'C';
		}
		else{
			obj[i] = 0;
			idx_Z.push_back(i);
			strColNames[i] = "z" + ToString(i-DataModel::n);
			colname[i] = new char[8];
			strcpy(colname[i], strColNames[i].c_str()); // make name for z
			if (modeltype == MODEL_TYPE_LP)
				ctype[i] = 'C';
			else if( method == 2 &&modeltype == MODEL_TYPE_MIP)
				ctype[i] = 'B';
			else if( method == 1 && modeltype == MODEL_TYPE_MIP)
				ctype[i] = 'I';
		}
	}
	delete []strColNames;

	for (int i =0; i <num_cols; i ++)       // fill in upper bound
	{
		if (i<DataModel::n)
			ub[i] = CPX_INFBOUND;
		else 
		{
			if (method ==1)
				ub[i] = aggGroups[ i - DataModel::n].size();
			else if (method == 2)
				ub[i] = 1;

		}
	}

	if (modeltype == MODEL_TYPE_LP)
	{
		status = CPXnewcols (cpxenv, cpxlp, num_cols, obj, NULL, ub, NULL, colname);
	}
	else if (modeltype == MODEL_TYPE_MIP)
	{
		status = CPXnewcols (cpxenv, cpxlp, num_cols, obj, NULL, ub, ctype, colname);
	}
	if ( status ) {
		fprintf (stderr, "Could not add new columns.\n");
		return 0;
	}

	modelType = modeltype;

	/* adding rows to model */
	for (int i = 0; i < aggGroups.size(); i++) 
	{
		VEC_INT group = aggGroups[i];
		int totaNumCols = DataModel::n +1 ;
		rmatind = (int * )   malloc (totaNumCols * sizeof(int));
		rmatval = (double *) malloc (totaNumCols * sizeof(double));
		memset(rmatval, 0, sizeof(double)*totaNumCols);

		for (int c = 0; c < DataModel::n; c++)
		{
			double a_j = 0;
			int idx_a_j;
			for (int r = 0; r < group.size(); r++)
			{
				int idx_row = group[r];
				VEC_MTX_ROW rowEntries = dm.Arows[idx_row];
				PAIR_INDEX_VALUE entryPair = rowEntries[c];
				if (method ==1 ) // aggregation
				{
					rmatval[c] += (float)((int)(entryPair.second/dm.rhs[i]*100+0.5))/100.0;
				}
				else if (method ==2) // maximum
				{
					if ((float)((int)(entryPair.second/dm.rhs[i]*100+0.5))/100.0 > rmatval[c])
					{
						rmatval[c] = (float)((int)(entryPair.second/dm.rhs[i]*100+0.5))/100.0;
					}
				}
				rmatind[c] = c;
			}
		}

		rmatind[DataModel::n] = DataModel::n +i ;
		rmatval[DataModel::n] = 1;
		double rhs[1];
		rhs[0] = method == 1 ? group.size() : 1;
		char   sense[1];
		sense[0] = 'G';

		int    rmatbeg[1];
		rmatbeg[0] = 0;
		char** rowName = new char*; 
		string cstName = string("SN") +ToString(i);
		char* pStr = new char[20];
		strcpy(pStr, cstName.c_str());
		rowName[0] = pStr;
		status = CPXaddrows (cpxenv, cpxlp, 0, 1, totaNumCols , rhs,
			sense, rmatbeg, rmatind, rmatval,
			NULL, rowName);
		if ( status ) {
			fprintf (stderr, "Could not add new rows.\n");
			return 0;
		}
		delete rowName;
		free_and_null((char **)&rmatind); 
		free_and_null((char **)&rmatval); 

	}

	/* add resource constraint*/
	if (resConst)
	{
		rmatind = (int * )   malloc (DataModel::n * sizeof(int));
		rmatval = (double *) malloc (DataModel::n * sizeof(double));
		for (int i = 0; i < DataModel::n ; i++)
		{
			rmatind[i] = i;
			rmatval[i] = 1;
		}
		double rhs[1];
		rhs[0] = 1;
		char   sense[1];
		sense[0] = 'E';
		int    rmatbeg[1];
		rmatbeg[0] = 0;
		char ** rowName = new char*;
		rowName[0] = "RES";

		status = CPXaddrows (cpxenv, cpxlp, 0, 1, DataModel::n , rhs,
			sense, rmatbeg, rmatind, rmatval, NULL, rowName);
	}



	/* adding cardinality constraint */
	/* add new column for slack var*/
	lb[0] = 0;
	ub[0] = 0; // exactly k violated
	int totaNumCols = aggGroups.size() ;

	/* prepare row data*/

	rmatind = (int * )   malloc (totaNumCols * sizeof(int));
	rmatval = (double *) malloc (totaNumCols * sizeof(double));

	for (int j =0; j <totaNumCols; j ++)
	{
		rmatind[j] = idx_Z[j];
		rmatval[j] = method == 1 ? 1 : aggGroups[j].size();
	}
	int    rmatbeg[1];
	double rhs[1];
	char   sense[1];
	rhs[0] = DataModel::k;
	sense[0] = cardinalitySign;
	rmatbeg[0] = 0;
	char ** rowName = new char*;
	rowName[0] = "CARD";
	status = CPXaddrows (cpxenv, cpxlp, 0, 1, totaNumCols, rhs,
		sense, rmatbeg, rmatind, rmatval,
		NULL, rowName);
	delete rowName;
	if ( status ) {
		fprintf (stderr, "Could not add new rows.\n");
		return 0;
	}

	string fileName; string suffix;
	if (method == 1)
		fileName = "aggregation";
	else
		fileName = "aggMax";
	if (modeltype == MODEL_TYPE_LP)
		suffix = "R.lp" ;
	else if (modeltype == MODEL_TYPE_MIP)
		suffix = "MIP.lp";
	fileName = fileName + suffix;
	writeModel2File(fileName.c_str());

	free_and_null((char **)&obj);
	free_and_null((char **)&lb);
	free_and_null((char **)&ub);
	free_and_null((char **)&ctype);
	free_and_null((char **)&colname); 
	free_and_null((char **)&rmatind); 
	free_and_null((char **)&rmatval); 
	return 0;
}



int OptModel::Agg_buildInitialAggGroups(DataModel& dm, int groupSize, int normType)
{
	dm_p = &dm;
	double** p = new double*[DataModel::m];
	for (int i = 0; i < DataModel::m; i++)
	{
		p[i] = new double[DataModel::m];
		memset(p[i], 0, sizeof(double)*DataModel::m);
	}

	multimap<double, pair<int, int> >  sortedDistance;
	for (int r = 0; r  < DataModel::m; r++)
	{
		LHS lhs = dm_p->Arows[r];
		RHS rhs = dm_p->rhs[r];
		INEQUALITY scenario1(make_pair(lhs,rhs));
		for (int c = 0; c< DataModel::m; c++)
		{
			if (c <= r )
			{
				p[r][c] = p[c][r];
				continue;
			}
			 INEQUALITY scenario2(make_pair(dm_p->Arows[c],dm_p->rhs[c]));
			p[r][c] = calculateDistance(scenario1,scenario2,normType);
			sortedDistance.insert( make_pair( p[r][c],  pair<int,int>(r,c)));
		}
	}


	int num_group = ceil(((double)DataModel::m)/groupSize);
	int * taken_p = new int[DataModel::m]; // every entry tells which group this element goes
	memset(taken_p, -1, sizeof(int)*DataModel::m);

	for (multimap<double, pair<int, int> >::iterator itr = sortedDistance.begin(); itr != sortedDistance.end(); ++itr)
	{
		int idx_c1 = itr->second.first;
		int idx_c2 = itr->second.second;
		if (taken_p[idx_c1] == -1 && taken_p[idx_c2] == -1 )// none of the scenario has been used yet
		{
			vector<int> group;
			group.push_back(idx_c1);
			group.push_back(idx_c2);
			if (aggGroups.size() < num_group) // if not all groups have been generated, put the new group in a new slot.
			{
				aggGroups.push_back(group);
				taken_p[idx_c1] = aggGroups.size() -1;
				taken_p[idx_c2] = aggGroups.size() -1;
			}
			else // if there are already num_group there, just find a group that has space for two and put them there.
			{
				for (int i = 0 ; i < num_group ; i++)
				{
					if(aggGroups[i].size() > groupSize -2)
						continue;
					else
					{
						aggGroups[i].insert( aggGroups[i].end(), group.begin(), group.end());
						taken_p[idx_c1] = i;
						taken_p[idx_c2] = i;
						break;
					}

				}

			}
		}
		else 
		{
			if (taken_p[idx_c1] != -1 && taken_p[idx_c2] != -1) // if both scenarios have been used
			{
				continue;
			}
			else  // if exactly one of the two scenarios is used
			{
				int idx_group = taken_p[idx_c1] == -1 ?  taken_p[idx_c2] :  taken_p[idx_c1];
				int idx_scenario_not_used  = taken_p[idx_c1] == -1 ? idx_c1 : idx_c2;
				if (idx_group != num_group -1 &&  aggGroups[idx_group].size()< groupSize)
				{
					aggGroups[idx_group].push_back(idx_scenario_not_used);
					taken_p[idx_scenario_not_used] = idx_group;
				}
				else if (idx_group == num_group -1 &&  aggGroups[idx_group].size()< groupSize)
				{
					aggGroups[idx_group].push_back(idx_scenario_not_used);
					taken_p[idx_scenario_not_used] = idx_group;
				}
				else
				{
					continue;
				}
			}
		}
	}
	delete []taken_p;
	return 0;
}


int OptModel::Agg_buildInitialAggGroupsNew(DataModel& dm, int groupSize, int normType)
{
	dm_p = &dm;
	multimap<double, pair<int, int> >  sortedDistance;
	for (int r = 0; r  < DataModel::m; r++)
	{
		LHS lhs = dm_p->Arows[r];
		RHS rhs = dm_p->rhs[r];
		INEQUALITY scenario1(make_pair(lhs,rhs));
		for (int c = 0; c< DataModel::m; c++)
		{
			if (c <= r )
			{
				continue;
			}
			INEQUALITY scenario2(make_pair(dm_p->Arows[c],dm_p->rhs[c]));
			double distance = calculateDistance(scenario1,scenario2,normType);
			sortedDistance.insert( make_pair( distance,  pair<int,int>(r,c)));
		}
	}


	int num_group = DataModel::m/2;    //ceil(((double)DataModel::m)/groupSize);
	int * taken_p = new int[DataModel::m]; // every entry tells which group this element goes
	memset(taken_p, -1, sizeof(int)*DataModel::m);

	for (multimap<double, pair<int, int> >::iterator itr = sortedDistance.begin(); itr != sortedDistance.end(); ++itr)
	{
		int idx_c1 = itr->second.first;
		int idx_c2 = itr->second.second;
		if (taken_p[idx_c1] == -1 && taken_p[idx_c2] == -1 )// none of the scenario has been used yet
		{
			vector<int> group;
			group.push_back(idx_c1);
			group.push_back(idx_c2);
// 			if (aggGroups.size() < num_group) // if not all groups have been generated, put the new group in a new slot.
// 			{
				aggGroups.push_back(group);
				taken_p[idx_c1] = aggGroups.size() -1;
				taken_p[idx_c2] = aggGroups.size() -1;
// 			}
// 			else // if there are already num_group there, just find a group that has space for two and put them there.
// 			{
// 				for (int i = 0 ; i < num_group ; i++)
// 				{
// 					if(aggGroups[i].size() > groupSize -2)
// 						continue;
// 					else
// 					{
// 						aggGroups[i].insert( aggGroups[i].end(), group.begin(), group.end());
// 						taken_p[idx_c1] = i;
// 						taken_p[idx_c2] = i;
// 						break;
// 					}
// 
// 				}
// 
// 			}
		}
		else 
		{
			if (taken_p[idx_c1] != -1 && taken_p[idx_c2] != -1) // if both scenarios have been used
			{
				continue;
			}
			else  // if exactly one of the two scenarios is used
			{
				int idx_group = taken_p[idx_c1] == -1 ?  taken_p[idx_c2] :  taken_p[idx_c1];
				int idx_scenario_not_used  = taken_p[idx_c1] == -1 ? idx_c1 : idx_c2;
				aggGroups[idx_group].push_back(idx_scenario_not_used);
				taken_p[idx_scenario_not_used] = idx_group;
			}
		}
	}
	delete []taken_p;
	return 0;
}

int OptModel::Agg_getViolatedConsts(multimap<double, pair<int,int> >& scenarios_picked, int method)// aggregation when method = 1; taking max when method =2
{
	int num_cols = getNumCols();
	double *x = new double[num_cols];
	CPXgetx(cpxenv, cpxlp, x, 0, num_cols -1);
	vector<int> violatedZ;
	vector<int> satisfiedZ;
	for (int i = DataModel::n; i < DataModel::n + aggGroups.size(); i++)
	{
		if (x[i]>= 1-SMALL_DECIMAL) // find z =1
		{
			violatedZ.push_back(i-DataModel::n);
		}
		else
			satisfiedZ.push_back( i-DataModel::n);
	}

	int sum = 0;
	for (int i = 0; i < violatedZ.size(); i++)
	{
		int idx_group = violatedZ[i];
		if (method ==1)
			sum += x[idx_group + DataModel::n];
		else
			sum += aggGroups[idx_group].size();
	}
	int num_to_pick = DataModel::k - sum +1;
	if (num_to_pick <1)
	{
		num_to_pick =1;
	}

	int cnt = 0;
	for (int i = 0; i < satisfiedZ.size(); i++)
	{
		int idx_group = satisfiedZ[i];
		vector<int> group = aggGroups[idx_group];
		if (group.size() == 1)
			continue;
		for (int j = 0; j < group.size(); j++)
		{
			int idx_scenario = group[j];
			double violation = IsXValidForScenarioZNormalized(idx_scenario, x);
			if (violation <= -SMALL_DECIMAL)//violated
			{
				scenarios_picked.insert( make_pair( violation, pair<int,int>(idx_scenario, idx_group) ) );
				cnt++;
			}
		}
	}

	if (method == 1 /*&& scenarios_picked.empty()*/)
	{
		cnt = 0 ;
		for (int i = 0; i < violatedZ.size(); i++)
		{
			int idx_group = violatedZ[i];
			vector<int> group = aggGroups[idx_group];
			if (group.size() == 1)
				continue;
			for (int j = 0; j < group.size(); j++)
			{
				int idx_scenario = group[j];
				double violation = IsXValidForScenarioZNormalized(idx_scenario, x);
				if (violation <= -SMALL_DECIMAL)//violated
				{
					scenarios_picked.insert( make_pair( violation, pair<int,int>(idx_scenario, idx_group) ) );
					cnt++;
				}
			}
		}

	}

	delete []x;
	return num_to_pick;
}


int OptModel::Agg_updateGroups(multimap<double, pair<int,int> > violated_indicies, int num_to_pick)
{
	int cnt =0;
	for (multimap<double, pair<int,int> >::iterator itr = violated_indicies.begin(); itr != violated_indicies.end(); ++itr)
	{
		if (cnt == num_to_pick)
			break;
		int idx_scenario = itr->second.first;
		int idx_group = itr->second.second;
		// remove the scenario from the group
		vector<int> group = aggGroups[idx_group];
		for (vector<int>::iterator itr = group.begin(); itr != group.end(); ++itr)
		{
			if ( *itr == idx_scenario)
			{
				group.erase(itr);
				if (group.empty())
				{
					int temp = 0;
					vector<VEC_INT>::iterator itrG ;
					for ( itrG = aggGroups.begin(); itrG != aggGroups.end(); ++itrG)
					{
						if (temp == idx_group)
						{
							break;
						}
						temp++;
					}
					aggGroups.erase(itrG);
				}
				else
					aggGroups[idx_group] = group;
				cnt++;
				break;
			}
		}
	}

			// add the scenario to the groups
	cnt = 0;
	for (multimap<double, pair<int,int> >::iterator itr = violated_indicies.begin(); itr != violated_indicies.end(); ++itr)
	{
		if (cnt == num_to_pick)
			break;
		int idx_scenario = itr->second.first;
		vector<int> group;
		group.push_back(idx_scenario);
		aggGroups.push_back(group);
		cnt++;
	}

	return 0;
}


int OptModel::Agg_updateGroupsFromSmallestViolation(multimap<double, pair<int,int> > violated_indicies, int num_to_pick)
{
	int cnt =0;
	for (multimap<double, pair<int,int> >::reverse_iterator itr = violated_indicies.rbegin(); itr != violated_indicies.rend(); ++itr)
	{
		if (cnt == num_to_pick)
			break;
		int idx_scenario = itr->second.first;
		int idx_group = itr->second.second;
		// remove the scenario from the group
		vector<int> group = aggGroups[idx_group];
		for (vector<int>::iterator itr = group.begin(); itr != group.end(); ++itr)
		{
			if ( *itr == idx_scenario)
			{
				group.erase(itr);
				if (group.empty())
				{
					int temp = 0;
					vector<VEC_INT>::iterator itrG ;
					for ( itrG = aggGroups.begin(); itrG != aggGroups.end(); ++itrG)
					{
						if (temp == idx_group)
						{
							break;
						}
						temp++;
					}
					aggGroups.erase(itrG);
				}
				else
					aggGroups[idx_group] = group;
				cnt++;
				break;
			}
		}
	}

	// add the scenario to the groups
	cnt = 0;
	for (multimap<double, pair<int,int> >::reverse_iterator itr = violated_indicies.rbegin(); itr != violated_indicies.rend(); ++itr)
	{
		if (cnt == num_to_pick)
			break;
		int idx_scenario = itr->second.first;
		vector<int> group;
		group.push_back(idx_scenario);
		aggGroups.push_back(group);
		cnt++;
	}

	return 0;
}


int OptModel::Agg_populatePartialFixedModel( DataModel& dm, OptModel::MODEL_TYPE modeltype, char cardinalitySign, bool resConst, multimap<double, pair<int,int> > violated_indicies)
{
	setDataModel(&dm);
	idx_X.clear();
	idx_Z.clear();
	idx_S.clear();
	double *obj     = NULL;
	double *lb      = NULL;
	double *ub      = NULL;
	char   *ctype   = NULL;
	int    *rmatind = NULL;
	double *rmatval = NULL;
	char** colname = NULL;

	int    status = 0; 
	int num_cols = DataModel::n + DataModel::m;

	CPXchgobjsen (cpxenv, cpxlp, 1); /* Minimization problem */

	/* Allocate colcnt-sized arrays */

	obj     = (double *) malloc (num_cols* sizeof(double));
	lb      = (double *) malloc (num_cols * sizeof(double));
	ub      = (double *) malloc (num_cols * sizeof(double));
	ctype   = (char *)   malloc (num_cols * sizeof(char));
	colname = (char**) malloc(num_cols*sizeof(char*));

	if ( obj     == NULL ||
		lb      == NULL ||
		ub      == NULL ||
		ctype   == NULL ) {
			fprintf (stderr, "Could not allocate colcnt arrays\n");
			status = CPXERR_NO_MEMORY;
			return 0;
	}

	string* strColNames = new string[num_cols];
	for (int i =0; i < num_cols; i ++)
	{
		if (i < DataModel::n )
		{
			obj[i] = dm.costVec[i].first;                      // fill in obj coefficients
			idx_X.push_back(i);                        // store index for x variables
			strColNames[i] = "x"+ToString(i);
			colname[i] = new char[8];
			strcpy(colname[i], strColNames[i].c_str()); // make name for x
			ctype[i] = 'C';
		}
		else{
			obj[i] = 0;
			idx_Z.push_back(i);
			strColNames[i] = "z" + ToString(i-DataModel::n);
			colname[i] = new char[8];
			strcpy(colname[i], strColNames[i].c_str()); // make name for z
			if (modeltype == MODEL_TYPE_LP)
				ctype[i] = 'C';
			else if( modeltype == MODEL_TYPE_MIP)
				ctype[i] = 'B';
		}
	}
	delete []strColNames;

	for (int i =0; i < num_cols; i ++)       // fill in upper bound
	{
		if (i<DataModel::n)
			ub[i] = CPX_INFBOUND;
		else 
			ub[i] = 0;
	}

	for (multimap<double, pair<int,int> >::iterator itr =  violated_indicies.begin(); itr != violated_indicies.end(); ++itr)
	{
		pair<int, int> entry = itr->second;
		int idx = DataModel::n + entry.first;
		ub[idx] = 1;
	}
	


	if (modeltype == MODEL_TYPE_LP)
	{
		status = CPXnewcols (cpxenv, cpxlp, num_cols, obj, NULL, ub, NULL, colname);
	}
	else if (modeltype == MODEL_TYPE_MIP)
	{
		status = CPXnewcols (cpxenv, cpxlp, num_cols, obj, NULL, ub, ctype, colname);
	}
	if ( status ) {
		fprintf (stderr, "Could not add new columns.\n");
		return 0;
	}

	modelType = modeltype;

	/* adding rows to model */
	for (int i = 0; i < DataModel::m; i++) 
	{

		int totaNumCols = DataModel::n +1 ;
		rmatind = (int * )   malloc (totaNumCols * sizeof(int));
		rmatval = (double *) malloc (totaNumCols * sizeof(double));

		VEC_MTX_ROW rowEntries = dm.Arows[i];  // fill in the entries for a row
		double rhs[1];
		char   sense[1];

		switch (dm.senses[i])
		{
		case 'L':
			for (int j =0; j <DataModel::n; j ++)
			{
				PAIR_INDEX_VALUE entryPair = rowEntries[j];
				rmatind[j]= entryPair.first;
				rmatval[j] = - (float)((int)(entryPair.second/dm.rhs[i]*100+0.5))/100.0; 
			}
			rmatind[DataModel::n] = DataModel::n + i;
			rmatval[DataModel::n] = 1;
			rhs[0] = -1;
			sense[0] = 'G';

			break;
		case 'G':
			for (int j =0; j <DataModel::n; j ++)
			{
				PAIR_INDEX_VALUE entryPair = rowEntries[j];
				rmatind[j] = entryPair.first;
				rmatval[j]= (float)((int)(entryPair.second/dm.rhs[i]*100+0.5))/100.0; 
			}
			rmatind[DataModel::n] = DataModel::n + i;
			rmatval[DataModel::n] = 1;
			rhs[0] = 1;
			sense[0] = 'G';

			break;
		case 'E':
			for (int j =0; j <DataModel::n; j ++)
			{
				PAIR_INDEX_VALUE entryPair = rowEntries[j];
				rmatind[j] = entryPair.first;
				rmatval[j] = (float)((int)(entryPair.second/dm.rhs[i]*100+0.5))/100.0; 
			}
			rmatind[DataModel::n] = DataModel::n + i;
			rmatval[DataModel::n] = 1;
			rhs[0] = 1;
			sense[0] = 'E';

			break;
		}
		int    rmatbeg[1];
		rmatbeg[0] = 0;
		char** rowName = new char*; 
		string cstName =  /* i == dm.Arows.size()-1 ?  string("RES")  : */string("SN") +ToString(i);
		char* pStr = new char[20];
		strcpy(pStr, cstName.c_str());
		rowName[0] = pStr;
		status = CPXaddrows (cpxenv, cpxlp, 0, 1, totaNumCols , rhs,
			sense, rmatbeg, rmatind, rmatval,
			NULL, rowName);
		if ( status ) {
			fprintf (stderr, "Could not add new rows.\n");
			return 0;
		}
		delete rowName;
		free_and_null((char **)&rmatind); 
		free_and_null((char **)&rmatval); 

	}

	/* add resource constraint*/
	if (resConst)
	{
		rmatind = (int * )   malloc (DataModel::n * sizeof(int));
		rmatval = (double *) malloc (DataModel::n * sizeof(double));
		for (int i = 0; i < DataModel::n ; i++)
		{
			rmatind[i] = i;
			rmatval[i] = 1;
		}
		double rhs[1];
		rhs[0] = 1;
		char   sense[1];
		sense[0] = 'E';
		int    rmatbeg[1];
		rmatbeg[0] = 0;
		char ** rowName = new char*;
		rowName[0] = "RES";

		status = CPXaddrows (cpxenv, cpxlp, 0, 1, DataModel::n , rhs,
			sense, rmatbeg, rmatind, rmatval, NULL, rowName);
	}
	/* adding cardinality constraint */
	lb[0] = 0;
	ub[0] = 0; // exactly k violated
	int totaNumCols = DataModel::m;

	/* prepare row data*/

	rmatind = (int * )   malloc (totaNumCols * sizeof(int));
	rmatval = (double *) malloc (totaNumCols * sizeof(double));

	for (int j =0; j <DataModel::m; j ++)
	{
		rmatind[j] = idx_Z[j];
		rmatval[j] = 1;
	}
	int    rmatbeg[1];
	double rhs[1];
	char   sense[1];
	rhs[0] = DataModel::k;
	sense[0] = cardinalitySign;
	rmatbeg[0] = 0;
	char ** rowName = new char*;
	rowName[0] = "CARD";
	status = CPXaddrows (cpxenv, cpxlp, 0, 1, DataModel::m, rhs,
		sense, rmatbeg, rmatind, rmatval,
		NULL, rowName);
	delete rowName;
	if ( status ) {
		fprintf (stderr, "Could not add new rows.\n");
		return 0;
	}

	string fileName ;
	if (modeltype == MODEL_TYPE_LP)
		fileName = modelName + "R.lp" ;
	else if (modeltype == MODEL_TYPE_MIP)
		fileName = modelName + "MIP.lp";
	writeModel2File(fileName.c_str());
	//CPXwriteprob(cpxenv, lp, "initialLPRelaxation.lp", NULL);

	free_and_null((char **)&obj);
	free_and_null((char **)&lb);
	free_and_null((char **)&ub);
	free_and_null((char **)&ctype);
	free_and_null((char **)&colname); 
	free_and_null((char **)&rmatind); 
	free_and_null((char **)&rmatval); 

	return 0;
}



int OptModel::AddObjectiveCut( double objVal, char sign)
{
		int totaNumCols = DataModel::n  ;
		int* rmatind = (int * )   malloc (totaNumCols * sizeof(int));
		double* rmatval = (double *) malloc (totaNumCols * sizeof(double));

		for (int j =0; j <DataModel::n; j ++)
		{
			rmatind[j]= j;
			rmatval[j] = 1; 
		}
		double rhs[1];
		char   sense[1];
		rhs[0] = objVal;
		sense[0] = sign;
		int    rmatbeg[1];
		rmatbeg[0] = 0;
		char** rowName = new char*; 
		string cstName = "ObjCut";
		char* pStr = new char[20];
		strcpy(pStr, cstName.c_str());
		rowName[0] = pStr;
		int status = CPXaddrows (cpxenv, cpxlp, 0, 1, totaNumCols , rhs,
			sense, rmatbeg, rmatind, rmatval,
			NULL, rowName);
		if ( status ) {
			fprintf (stderr, "Could not add new rows.\n");
			return 0;
		}
		delete rowName;
		free_and_null((char **)&rmatind); 
		free_and_null((char **)&rmatval); 

		return 0;
}

void OptModel::Agg_recordGroupInfo(multimap<double, pair<int, int> >& violated_indicies, int num_to_pick)
{
	string logName = CAT_NAME+ "groupinfo.csv";
	ofstream log(logName.c_str(),ios_base::out|ios_base::app);
	for (int i = 0; i < aggGroups.size(); i++)
	{
		log << "GroupNum," << i << ",[" ;
		VEC_INT group = aggGroups[i];
		for (int j = 0; j < group.size(); j++)
			log << group[j]<<":";
		log << "],";
	}
    log<<endl;
	int num_cols = getNumCols();
	double *x = new double[num_cols];
	CPXgetx(cpxenv, cpxlp, x, 0, num_cols -1);
	for(int i = 0; i < aggGroups.size(); i++)
	{
		log << x[DataModel::n + i] <<",";
		log<< IsXValidForConstrain(i, x) <<",";
		log << "[";
		VEC_INT group = aggGroups[i];
		for (int j = 0; j < group.size(); j++)
		{
			int idx_scenario = group[j];
			double violation =IsXValidForScenarioZNormalized(idx_scenario,x);
			log << (violation <= SMALL_DECIMAL ? 1 : 0 );
			log << ":";
		}
		log << "],"; 
	}

	int cnt = 0; 
	log<< endl ;
	for(multimap<double, pair<int, int> >::iterator itr = violated_indicies.begin(); itr != violated_indicies.end(); ++itr)
	{ 
		if (cnt == num_to_pick)
			break;
		log << "Scenario_ID, " << itr->second.first<< ",Group_ID," <<itr->second.second<<",";
		cnt ++;
	}
	log <<endl;

	log.close();

}

double OptModel::calculateLBbyEuclideanDistance()
{
	multimap<double, int> sortedDistance;
	for(int i = 0; i < DataModel::m; i++)
	{
		double sqSum = 0;
		VEC_MTX_ROW coefs= dm_p->Arows[i];
		for (int j = 0; j < DataModel::n; j++)
		{
			sqSum += coefs[j].second*coefs[j].second;
		}
		sortedDistance.insert(make_pair(1/sqrt(sqSum), i));
	}
	
	int cnt = 0;
	double k_minimum ; int index;
	for (multimap<double, int>::reverse_iterator ritr = sortedDistance.rbegin(); ritr != sortedDistance.rend(); ++ritr, cnt++)
	{
		if (cnt == DataModel::k)
		{
			k_minimum = ritr->first;
			index = ritr->second;
			break;
		}
	}
	
	double min_c = 1000;
	for (int i = 0; i <dm_p->costVec.size(); i++ )
	{
		if (dm_p->costVec[i].second < min_c)
		{
			min_c = dm_p->costVec[i].first;
		}
	}
	
	return min_c*k_minimum;

}

void OptModel::setGeneralCutPoolIdentifier(string identifier)
{
	GeneralCutPool.setCutIdentifier(identifier);
}


int OptModel::getNumCutsInCutPool(string cutType)
{
	int num_cuts = 0;
	if (cutType == string("IC"))
	{
		num_cuts = ICCutPool.getNumCuts();
	}
	else if (cutType == string("MIX"))
	{
		num_cuts = MixCutPool.getNumCuts();
	}
	else if( cutType == string("MIR"))
		num_cuts = MIRCutPool.getNumCuts();
	else if ( cutType == string("CH"))
	{
		num_cuts = CHCutPool.getNumCuts();
	}
	else if (cutType == string("ALL"))
	{
		num_cuts = ICCutPool.getNumCuts();
		num_cuts += MixCutPool.getNumCuts();
		// 		MIRCutPool.extractCutsFromModel(cpxenv, cpxlp);
		num_cuts += CHCutPool.getNumCuts();
	}
	else
	{
		//		num_cuts = GeneralCutPool.getNumCuts();
	}

	return num_cuts;
}


double OptModel::startTimer(void)
{
	double timestamp;
	CPXgettime(cpxenv,&timestamp);
	startingtime = timestamp;
	return timestamp;
}

double OptModel::getElapsedTime(void)
{
	double timestamp;
	CPXgettime(cpxenv,&timestamp);
	 return timestamp - startingtime;
}

CUT_COEFICIENTS OptModel::getObjective()
{
	double *obj = new double[DataModel::n];
	CPXgetobj(cpxenv, cpxlp,obj, 0, DataModel::n-1);
	CUT_COEFICIENTS coefs;
	for (int i = 0; i < DataModel::n; i++)
	{
		coefs.push_back(make_pair(i, obj[i]));
	}
	delete []obj;
	return coefs;
}

int OptModel::removeConstraintFromCPLEX(string& constraintName)
{
	int num_rows = CPXgetnumrows(cpxenv,  cpxlp);
	int surplus;
	vector<string> rowNames;

	int status = CPXgetrowname (cpxenv, cpxlp, NULL, NULL, 0, &surplus, 0,num_rows-1);
	int cur_rownamespace = - surplus;
	if ( cur_rownamespace > 0 ) 
	{
		char          **cur_rowname = (char **) malloc (sizeof(char *)*(num_rows));
		char          *cur_rownamestore = (char *)  malloc (cur_rownamespace);
		if ( cur_rowname      == NULL ||
			cur_rownamestore == NULL   ) {
				fprintf (stderr, "Failed to get memory for column names.\n");
				status = -1;
				cout<< endl <<" error in printTableau"<<endl;
				return 0;
		}
		status = CPXgetrowname (cpxenv, cpxlp, cur_rowname, cur_rownamestore, 
			cur_rownamespace, &surplus, 0, num_rows-1);
		// map name with column index
		for (int i = 0; i <num_rows ; i ++)
		{
			string rowName = cur_rowname[i];
			if(rowName != constraintName)
				continue;

			// remove from CPLEX model
			int rlist = i; int rlen = 1;
			CPXpivotin(cpxenv, cpxlp, &rlist, rlen);

			int status = CPXdelrows (cpxenv, cpxlp, i, i);
			if (status ==0)
			{
//				cout << rowName << " was successfully removed"<< endl;
			}
			else
				cout << rowName << " was NOT successfully removed"<< endl;

			break;
		}
		free_and_null( (char **)&cur_rowname);
		free_and_null( (char **)&cur_rownamestore);

	}
	return 0;
}


int OptModel::removeConstraintsFromCPLEXbyKeyWord(string& key)
{
	bool found = true;
	while(found) {
		found = false;
		int num_rows = CPXgetnumrows(cpxenv,  cpxlp);
		int surplus;
		vector<string> rowNames;

		int status = CPXgetrowname (cpxenv, cpxlp, NULL, NULL, 0, &surplus, 0,num_rows-1);
		int cur_rownamespace = - surplus;
		if ( cur_rownamespace > 0 ) 
		{
			char          **cur_rowname = (char **) malloc (sizeof(char *)*(num_rows));
			char          *cur_rownamestore = (char *)  malloc (cur_rownamespace);
			if ( cur_rowname      == NULL ||
				cur_rownamestore == NULL   ) {
					fprintf (stderr, "Failed to get memory for column names.\n");
					status = -1;
					cout<< endl <<" error in printTableau"<<endl;
					return 0;
			}
			status = CPXgetrowname (cpxenv, cpxlp, cur_rowname, cur_rownamestore, 
				cur_rownamespace, &surplus, 0, num_rows-1);
			// map name with column index
			for (int i = 0; i <num_rows ; i ++)
			{
				string rowName = cur_rowname[i];
				if(rowName.find(key) == std::string::npos)
					continue;

				// remove from CPLEX model
				int rlist = i; int rlen = 1;
				CPXpivotin(cpxenv, cpxlp, &rlist, rlen);

				int status = CPXdelrows (cpxenv, cpxlp, i, i);
				if (status ==0)
				{
					cout << rowName << " was successfully removed"<< endl;
					found = true;
					break;
				}
				else
					cout << rowName << " was NOT successfully removed"<< endl;

			}
			free_and_null( (char **)&cur_rowname);
			free_and_null( (char **)&cur_rownamestore);
		}
	}
	return 0;
}


int OptModel::removeVariableFromCPLEX(string& varName)
{
	int num_cols = CPXgetnumcols(cpxenv,  cpxlp);
	int surplus;
	vector<string> rowNames;

	int status = CPXgetcolname (cpxenv, cpxlp, NULL, NULL, 0, &surplus, 0,num_cols-1);
	int cur_colnamespace = - surplus;
	if ( cur_colnamespace > 0 ) 
	{
		char          **cur_colname = (char **) malloc (sizeof(char *)*(num_cols));
		char          *cur_colnamestore = (char *)  malloc (cur_colnamespace);
		if ( cur_colname      == NULL ||
			cur_colnamestore == NULL   ) {
				fprintf (stderr, "Failed to get memory for column names.\n");
				status = -1;
				cout<< endl <<" error in printTableau"<<endl;
				return 0;
		}
		status = CPXgetcolname (cpxenv, cpxlp, cur_colname, cur_colnamestore, 
			cur_colnamespace, &surplus, 0, num_cols-1);
		// map name with column index
		for (int i = 0; i <num_cols ; i ++)
		{
			string colName = cur_colname[i];
			if(colName != varName)
				continue;

			int status = CPXdelcols (cpxenv, cpxlp, i, i);
			if (status ==0)
			{
				//				cout << rowName << " was successfully removed"<< endl;
			}
			else
				cout << colName << " was NOT successfully removed"<< endl;

			break;
		}
		free_and_null( (char **)&cur_colname);
		free_and_null( (char **)&cur_colnamestore);

	}
	return 0;
}


string OptModel::getVarNameByIndex(int index)
{
	int num_cols  = CPXgetnumcols(cpxenv,  cpxlp) ; // need to know why
	int surplus;
	int status = CPXgetcolname (cpxenv, cpxlp, NULL, NULL, 0, &surplus, 0,num_cols-1);
	int cur_colnamespace = - surplus;
	string varName;
	if ( cur_colnamespace > 0 ) {
		char          **cur_colname      = (char **) malloc (sizeof(char *)*(num_cols));
		char          *cur_colnamestore = (char *)  malloc (cur_colnamespace);
		if ( cur_colname      == NULL ||
			cur_colnamestore == NULL   ) {
				fprintf (stderr, "Failed to get memory for column names.\n");
				status = -1;
				return 0;
		}
		status = CPXgetcolname (cpxenv, cpxlp, cur_colname, cur_colnamestore, 
			cur_colnamespace, &surplus, 0, num_cols-1);

		varName = cur_colname[index];
	}

	return varName;
}
