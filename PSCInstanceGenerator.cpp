#include "PSCInstanceGenerator.h"

PSCInstanceGenerator::PSCInstanceGenerator(void)
{
}

PSCInstanceGenerator::~PSCInstanceGenerator(void)
{
}
int PSCInstanceGenerator::initilizeCPXenv()
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
	status = CPXsetintparam (cpxenv, CPX_PARAM_SCRIND, CPX_ON);
	if ( status ) {
		fprintf (stderr, 
			"Failure to turn on screen indicator, error %d.\n", status);
		return 0;
	}

	status = CPXsetintparam(cpxenv, CPX_PARAM_THREADS, 1);
	status = CPXsetintparam(cpxenv, CPX_PARAM_CLOCKTYPE, 1);

	return 1;
}


int PSCInstanceGenerator::initilizeCPXlp(string modelName)
{
	int status = 0;
	if (cpxenv == NULL)
	{
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




bool PSCInstanceGenerator::generateModel(string& modelName, OptModel::MODEL_TYPE modelType, bool flag_random_cost)
{
	initilizeCPXenv();
	initilizeCPXlp(modelName);
	createObjective(modelType, num_levels, flag_random_cost);
	createSetCoverConstraints(1, c);
//	createRandomSetCoverConstraints(1, c);
	if (num_levels == 2)
	{
		createSetCoverConstraints(2);
	}
	createCardConstraints(1, coverage);
	if (num_levels == 2)
	{
		createCardConstraints(2, 1);
	}
 	string modelFileName = modelName + ".sav";
 	CPXwriteprob(cpxenv, cpxlp, modelFileName.c_str(), NULL);
	return true;
}


bool PSCInstanceGenerator::createObjective(OptModel::MODEL_TYPE modelType, int num_level, bool flag_random){

	double *obj     = NULL;
	double *lb      = NULL;
	double *ub      = NULL;
	char   *ctype   = NULL;
	int    *rmatind = NULL;
	double *rmatval = NULL;
	char** colname = NULL;

	int    status = 0; 
	int num_cols = n + n*num_level ;

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

	CRandomMersenne randGen(15);
	int cnt_idx = 0;
	for (int i=0  ; i < n; i++)
	{
		int r = i/grid_num_col;
		int c = i%grid_num_col;
		obj[cnt_idx] = flag_random ?  randGen.IRandom(1,100): 1;                      // fill in obj coefficients
		string varname = "x["+ ToString(r) + "][" + ToString(c)+"]";
		colname[cnt_idx] = new char[12];
		strcpy(colname[cnt_idx], varname.c_str()); // make name for x
		if (modelType == OptModel::MODEL_TYPE_LP)
			ctype[i] = 'C';
		else if( modelType == OptModel::MODEL_TYPE_MIP)
			ctype[i] = 'B';
		cnt_idx++;

	}
	for (int i = n; i < num_cols ; i++)
	{
		obj[i] = 0;
		string varName = "z" + ToString((i-n)/n+1) +"["  + ToString( (i%n)/grid_num_col)+"]["  + ToString( (i%n)%grid_num_col ) + "]";
		colname[i] = new char[8];
		strcpy(colname[i], varName.c_str()); // make name for z
		if (modelType == OptModel::MODEL_TYPE_LP)
			ctype[i] = 'C';
		else if( modelType == OptModel::MODEL_TYPE_MIP)
			ctype[i] = 'B';
	}


	for (int i =0; i < num_cols; i ++)       // fill in upper bound
	{
		ub[i] = 1;
	}

	if (modelType == OptModel::MODEL_TYPE_LP)
	{
		status = CPXnewcols (cpxenv, cpxlp, num_cols, obj, NULL, ub, NULL, colname);
	}
	else if (modelType == OptModel::MODEL_TYPE_MIP)
	{
		status = CPXnewcols (cpxenv, cpxlp, num_cols, obj, NULL, ub, ctype, colname);
	}
	if ( status ) {
		fprintf (stderr, "Could not add new columns.\n");
		return 0;
	}
    return true;
}


bool PSCInstanceGenerator::createSetCoverConstraints(int signal_level, int multi_cover_cnt)// a_ij is determined by signal level and the distance from the signal source
{
	int    *rmatind = NULL;
	double *rmatval = NULL;
	/* adding rows to model */

	int totaNumCols = n;
	for (int i = 0; i < n; i++)
	{
		rmatind = (int * )   malloc (totaNumCols * sizeof(int));
		rmatval = (double *) malloc (totaNumCols * sizeof(double));
		double rhs[1];
		char   sense[1];

		int cnt_nz = 0;
		for (int j = 0; j < n; j++)
		{
			if(i == j )
				continue;
			rmatind[cnt_nz] =  j;
			rmatval[cnt_nz] = get_a_ij(signal_level, i,j); 
			cnt_nz ++;
		}
		rmatind[cnt_nz] = (signal_level)*n + i;
		rmatval[cnt_nz] = multi_cover_cnt;
		rhs[0] = multi_cover_cnt;
		sense[0] = 'G';
		int    rmatbeg[1];
		rmatbeg[0] = 0;
		string cstName = "SN"+ToString(signal_level)+ "_" + ToString(i);
		char* pStr = new char[20];
		strcpy(pStr, cstName.c_str());
		int status = CPXaddrows (cpxenv, cpxlp, 0, 1, totaNumCols , rhs,
			sense, rmatbeg, rmatind, rmatval,
			NULL, &pStr);
		if ( status ) {
			fprintf (stderr, "Could not add new rows.\n");
			return 0;
		}
		free_and_null((char **)&rmatind); 
		free_and_null((char **)&rmatval); 
	}
	return true;
}

bool PSCInstanceGenerator::createRandomSetCoverConstraints(int signal_level, int multi_cover_cnt)// a_ij is generated randomly
{
	int    *rmatind = NULL;
	double *rmatval = NULL;
	/* adding rows to model */

	int totaNumCols = n;
	CRandomMersenne randGen(17); //*******
	for (int i = 0; i < n; i++)
	{
		rmatind = (int * )   malloc (totaNumCols * sizeof(int));
		rmatval = (double *) malloc (totaNumCols * sizeof(double));
		double rhs[1];
		char   sense[1];

		int cnt_nz = 0;
		map<int,int> ones; //*********
		for (int l = 0; l < 10; l++)
		{
			ones.insert(make_pair(randGen.IRandom(0, n-1), 1));
		}

		for (int j = 0; j < n; j++)
		{
			if(i == j )
				continue;
			rmatind[cnt_nz] =  j;
			rmatval[cnt_nz] = randGen.IRandom(0,1);//ones.find(j) == ones.end()? 0 : 1; //get_a_ij(signal_level, i,j); 
			cnt_nz ++;
		}
		rmatind[cnt_nz] = (signal_level)*n + i;
		rmatval[cnt_nz] = multi_cover_cnt;
		rhs[0] = multi_cover_cnt;
		sense[0] = 'G';
		int    rmatbeg[1];
		rmatbeg[0] = 0;
		string cstName = "SN"+ToString(signal_level)+ "_" + ToString(i);
		char* pStr = new char[20];
		strcpy(pStr, cstName.c_str());
		int status = CPXaddrows (cpxenv, cpxlp, 0, 1, totaNumCols , rhs,
			sense, rmatbeg, rmatind, rmatval,
			NULL, &pStr);
		if ( status ) {
			fprintf (stderr, "Could not add new rows.\n");
			return 0;
		}
		free_and_null((char **)&rmatind); 
		free_and_null((char **)&rmatval); 
	}
	return true;

}


bool PSCInstanceGenerator::get_a_ij(int level, int idx_signal, int idx_sensor)
{
	int x_signal = idx_signal/grid_num_col;
	int y_signal= idx_signal%grid_num_col;
	int x_sensor = idx_sensor/grid_num_col;
	int y_sensor= idx_sensor%grid_num_col;

	double distance = grid_edge_length*sqrt ((double)abs(x_signal-x_sensor)*abs(x_signal-x_sensor) +  (double)abs(y_signal-y_sensor)*abs(y_signal-y_sensor) );
	if (level == 1 )
	{
		return distance <= distance_level_1;
	}
	else if (level == 2)
	{
		return distance <= distance_level_2;
	}
    return true;
}



bool  PSCInstanceGenerator::createCardConstraints(int signal_level, double cr){

	int    *rmatind = NULL;
	double *rmatval = NULL;
	/* adding rows to model */

	int totaNumCols = n;
	rmatind = (int * )   malloc (totaNumCols * sizeof(int));
	rmatval = (double *) malloc (totaNumCols * sizeof(double));
	double rhs[1];
	char   sense[1];

	for (int i = 0; i < totaNumCols; i++)
	{
		rmatind[i] = ( signal_level )*n + i;
		rmatval[i] = 1; 
	}
	rhs[0] = n*(1-cr);
	sense[0] = 'L';
	int    rmatbeg[1];
	rmatbeg[0] = 0;
	string cstName = "CARD"+ToString(signal_level);
	char* pStr = new char[20];
	strcpy(pStr, cstName.c_str());
	int status = CPXaddrows (cpxenv, cpxlp, 0, 1, totaNumCols , rhs,
		sense, rmatbeg, rmatind, rmatval,
		NULL, &pStr);
	if ( status ) {
		fprintf (stderr, "Could not add new rows.\n");
		return 0;
	}
	free_and_null((char **)&rmatind); 
	free_and_null((char **)&rmatval); 
	return true;
}
