#include "PortfolioInstanceGenerator.h"
#include "stocc.h"

double PortfolioInstanceGenerator::A_MAX_RAND = 1.5;
double PortfolioInstanceGenerator::A_MIN_RAND = 0.8;
vector<double> PortfolioInstanceGenerator::ro;

PortfolioInstanceGenerator::PortfolioInstanceGenerator(void)
{
}

PortfolioInstanceGenerator::~PortfolioInstanceGenerator(void)
{
}


int PortfolioInstanceGenerator::initilizeCPXenv()
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


int PortfolioInstanceGenerator::initilizeCPXlp(string modelName)
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

bool PortfolioInstanceGenerator::createObjective(OptModel::MODEL_TYPE modelType){
	double *obj     = NULL;
	double *lb      = NULL;
	double *ub      = NULL;
	char   *ctype   = NULL;
	int    *rmatind = NULL;
	double *rmatval = NULL;
	char** colname = NULL;

	int    status = 0; 
	int num_cols = n + m;

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
			obj[i] = 1;                      // fill in obj coefficients
			strColNames[i] = "x"+ToString(i);
			colname[i] = new char[8];
			strcpy(colname[i], strColNames[i].c_str()); // make name for x
			ctype[i] = 'C';
		}
		else{
			obj[i] = 0;
			strColNames[i] = "z" + ToString(i-DataModel::n);
			colname[i] = new char[8];
			strcpy(colname[i], strColNames[i].c_str()); // make name for z
			if (modelType == OptModel::MODEL_TYPE_LP)
				ctype[i] = 'C';
			else if( modelType == OptModel::MODEL_TYPE_MIP)
				ctype[i] = 'B';
		}
	}
	delete []strColNames;

	for (int i =0; i < num_cols; i ++)       // fill in upper bound
	{
		if (i<n)
			ub[i] = CPX_INFBOUND;
		else 
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

	free_and_null((char **)&obj);
	free_and_null((char **)&lb);
	free_and_null((char **)&ub);
	free_and_null((char **)&ctype);
	free_and_null((char **)&colname); 
	free_and_null((char **)&rmatind); 
	free_and_null((char **)&rmatval); 
    return true;
}


bool PortfolioInstanceGenerator::createRandomObjective(OptModel::MODEL_TYPE modelType){  // not using unit vector as the objective constraint
	double *obj     = NULL;
	double *lb      = NULL;
	double *ub      = NULL;
	char   *ctype   = NULL;
	int    *rmatind = NULL;
	double *rmatval = NULL;
	char** colname = NULL;

	int    status = 0; 
	int num_cols = n + m;

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

	CRandomMersenne randGen(seed +1);
	string* strColNames = new string[num_cols];
	for (int i =0; i < num_cols; i ++)
	{
		if (i < DataModel::n )
		{
			obj[i] = randGen.IRandom(1,100); //1;                      // fill in obj coefficients
			strColNames[i] = "x"+ToString(i);
			colname[i] = new char[8];
			strcpy(colname[i], strColNames[i].c_str()); // make name for x
			ctype[i] = 'C';
		}
		else{
			obj[i] = 0;
			strColNames[i] = "z" + ToString(i-DataModel::n);
			colname[i] = new char[8];
			strcpy(colname[i], strColNames[i].c_str()); // make name for z
			if (modelType == OptModel::MODEL_TYPE_LP)
				ctype[i] = 'C';
			else if( modelType == OptModel::MODEL_TYPE_MIP)
				ctype[i] = 'B';
		}
	}
	delete []strColNames;

	for (int i =0; i < num_cols; i ++)       // fill in upper bound
	{
		if (i<n)
			ub[i] = CPX_INFBOUND;
		else 
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

	free_and_null((char **)&obj);
	free_and_null((char **)&lb);
	free_and_null((char **)&ub);
	free_and_null((char **)&ctype);
	free_and_null((char **)&colname); 
	free_and_null((char **)&rmatind); 
	free_and_null((char **)&rmatval); 
    return true;
}


bool PortfolioInstanceGenerator::createObjectiveNormalDistribution(OptModel::MODEL_TYPE modelType){  // not using unit vector as the objective constraint
	double *obj     = NULL;
	double *lb      = NULL;
	double *ub      = NULL;
	char   *ctype   = NULL;
	int    *rmatind = NULL;
	double *rmatval = NULL;
	char** colname = NULL;

	int    status = 0; 
	int num_cols = n + m;

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

	CRandomMersenne randPara(seed+2);
	vector<double> mu;
	for (int i =0; i < n; i++)
	{
		double val = A_MIN_RAND + ((A_MAX_RAND - A_MIN_RAND)*randPara.Random());
		mu.push_back(val);
	}
	string* strColNames = new string[num_cols];
	for (int i =0; i < num_cols; i ++)
	{
		if (i < DataModel::n )
		{
			obj[i] = 2-mu[i];//randGen.IRandom(1,100); //1;                      // fill in obj coefficients
			strColNames[i] = "x"+ToString(i);
			colname[i] = new char[8];
			strcpy(colname[i], strColNames[i].c_str()); // make name for x
			ctype[i] = 'C';
		}
		else{
			obj[i] = 0;
			strColNames[i] = "z" + ToString(i-DataModel::n);
			colname[i] = new char[8];
			strcpy(colname[i], strColNames[i].c_str()); // make name for z
			if (modelType == OptModel::MODEL_TYPE_LP)
				ctype[i] = 'C';
			else if( modelType == OptModel::MODEL_TYPE_MIP)
				ctype[i] = 'B';
		}
	}
	delete []strColNames;

	for (int i =0; i < num_cols; i ++)       // fill in upper bound
	{
		if (i<n)
			ub[i] = CPX_INFBOUND;
		else 
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

	free_and_null((char **)&obj);
	free_and_null((char **)&lb);
	free_and_null((char **)&ub);
	free_and_null((char **)&ctype);
	free_and_null((char **)&colname); 
	free_and_null((char **)&rmatind); 
	free_and_null((char **)&rmatval); 
    return true;
}



void PortfolioInstanceGenerator::addScenarioConstraints()
{
	double rhs_MIN = r;
	double A_MAX_RAND = 1.5;
	double A_MIN_RAND = 0.8;


	// generate random generators
	CRandomMersenne** pRandGen = new CRandomMersenne*[n] ;
	for (int i = 0; i < n; i++)
	{
		pRandGen[i]= new CRandomMersenne(seed + i*3);
	}

	int    *rmatind = NULL;
	double *rmatval = NULL;

	// generate matrix A
	ro = vector<double>(n,0);
	for (int i = 0; i<m; i++)
	{
		for (int t = 0 ; t < nr; t++)
		{
			double rhsV=r;
			rmatind = new int[n+1];
			rmatval = new double[n+1];
			for (int j=0; j<n; j++)
			{
				double coef =0;
				double temp = A_MIN_RAND + ((A_MAX_RAND - A_MIN_RAND)*pRandGen[j]->Random());
				coef = temp;//(float)((int)(temp/rhsV*10000+0.5))/10000.0;//
				rmatind[j] = j;
				rmatval[j] = coef;
				ro[j] +=  coef;
			}
			rmatind[n] = n + i;
			rmatval[n] = rhsV;

			double rhs[1];
			char   sense[1];
			rhs[0] = rhsV;
			sense[0] = 'G';

			int    rmatbeg[1];
			rmatbeg[0] = 0;

			string cstName = "SN"+ToString(i) + "R" + ToString(t);
			char* pStr = new char[20];
			strcpy(pStr, cstName.c_str());
			int status = CPXaddrows (cpxenv, cpxlp, 0, 1, n+1, rhs,
				sense, rmatbeg, rmatind, rmatval,
				NULL, &pStr);
			if ( status ) {
				fprintf (stderr, "Could not add new rows.\n");
				return;
			}
			delete []pStr;
			delete []rmatval;
			delete []rmatind;
		}
	}
	for(int j = 0; j<n; j++)
		ro[j]=ro[j]/m;
}

void PortfolioInstanceGenerator::addScenarioConstraintsNormalDistribution()
{
	double rhs_MIN = r;


	// generate random generators
	StochasticLib1** pRandGen = new StochasticLib1*[n] ;
	for (int i = 0; i < n; i++)
	{
		pRandGen[i]= new StochasticLib1(seed + i*3);
	}

	CRandomMersenne randPara(seed+2);
	vector<double> mu;
	vector<double> delta;
	for (int i =0; i < n; i++)
	{
		double val = A_MIN_RAND + ((A_MAX_RAND - A_MIN_RAND)*randPara.Random());
		mu.push_back(val);
		double val_delta = 0.1*val + ((0.8*val - 0.1*val)*randPara.Random());
		delta.push_back(val_delta);
	}

	int    *rmatind = NULL;
	double *rmatval = NULL;

	// generate matrix A
	for (int i = 0; i<m; i++)
	{
		for (int t = 0 ; t < nr; t++)
		{
			double rhsV=r;
			rmatind = new int[n+1];
			rmatval = new double[n+1];
			for (int j=0; j<n; j++)
			{
				double coef =0;
				//coef = A_MIN_RAND + ((A_MAX_RAND - A_MIN_RAND)*pRandGen[j]->Random());
				coef = pRandGen[j]->NormalTrunc(mu[j],delta[j],delta[j]);
				rmatind[j] = j;
				rmatval[j] = coef;
			}
			rmatind[n] = n + i;
			rmatval[n] = rhsV;

			double rhs[1];
			char   sense[1];
			rhs[0] = rhsV;
			sense[0] = 'G';

			int    rmatbeg[1];
			rmatbeg[0] = 0;

			string cstName = "SN"+ToString(i) + "R" + ToString(t);
			char* pStr = new char[20];
			strcpy(pStr, cstName.c_str());
			int status = CPXaddrows (cpxenv, cpxlp, 0, 1, n+1, rhs,
				sense, rmatbeg, rmatind, rmatval,
				NULL, &pStr);
			if ( status ) {
				fprintf (stderr, "Could not add new rows.\n");
				return;
			}
			delete []pStr;
			delete []rmatval;
			delete []rmatind;
		}
	}
}


void PortfolioInstanceGenerator::addScenarioConstraintsRandomRHS()
{
	double rhs_MIN = r;
	double A_MAX_RAND = 1.5;
	double A_MIN_RAND = 0.8;


	// generate random generators
	CRandomMersenne** pRandGen = new CRandomMersenne*[n] ;
	for (int i = 0; i < n; i++)
	{
		pRandGen[i]= new CRandomMersenne(seed + i*3);
	}
	CRandomMersenne randRHS(seed +9);

	int    *rmatind = NULL;
	double *rmatval = NULL;
	
	vector< vector<double> > matrixA;
	for (int i = 0; i < nr; i++)
	{
		matrixA.push_back(vector<double>(n,0));
	}
	// generate matrix A
	for (int i = 0; i<m; i++)
	{
		for (int t = 0 ; t < nr; t++)
		{
			double rhsV=r;
			rmatind = new int[n+1];
			rmatval = new double[n+1];
			for (int j=0; j<n; j++)
			{
				double coef =0;
				if (i == 0)
				{
					coef = A_MIN_RAND + ((A_MAX_RAND - A_MIN_RAND)*pRandGen[j]->Random());
					matrixA[t][j] = coef;
				}
				else
				{
					coef = matrixA[t][j];
				}
				rmatind[j] = j;
				rmatval[j] = coef;
			}
			rmatind[n] = n + i;
			rmatval[n] = 1 + ((A_MAX_RAND - 1)*randRHS.Random());  //rhsV;

			double rhs[1];
			char   sense[1];
			rhs[0] = rmatval[n];  //rhsV;
			sense[0] = 'G';

			int    rmatbeg[1];
			rmatbeg[0] = 0;

			string cstName = "SN"+ToString(i) + "R" + ToString(t);
			char* pStr = new char[20];
			strcpy(pStr, cstName.c_str());
			int status = CPXaddrows (cpxenv, cpxlp, 0, 1, n+1, rhs,
				sense, rmatbeg, rmatind, rmatval,
				NULL, &pStr);
			if ( status ) {
				fprintf (stderr, "Could not add new rows.\n");
				return;
			}
			delete []pStr;
			delete []rmatval;
			delete []rmatind;
		}
	}
}

void PortfolioInstanceGenerator::addDominantScenarioConstraintsRandomRHS()
{
	double rhs_MIN = r;
	double A_MAX_RAND = 1.5;
	double A_MIN_RAND = 0.8;


	// generate random generators
	CRandomMersenne** pRandGen = new CRandomMersenne*[n] ;
	for (int i = 0; i < n; i++)
	{
		pRandGen[i]= new CRandomMersenne(seed + i*3);
	}
	CRandomMersenne randRHS(seed +9);

	int    *rmatind = NULL;
	double *rmatval = NULL;

	// generate matrix A
	vector< vector<double> > matrixA;
	for (int i = 0; i < nr; i++)
	{
		vector<double> row;
		for (int j = 0 ; j < n; j++)
		{
			double coef = A_MIN_RAND + ((A_MAX_RAND - A_MIN_RAND)*pRandGen[j]->Random());
			row.push_back(coef);
		}
		matrixA.push_back(row);
	}

	// generate random rhs
	vector< vector<double> > matrixRHS;
	for (int i = 0; i < m; i++)
	{
		vector<double> row;
		for (int j = 0 ; j < nr; j++)
		{
			double rhs = 1 + ((A_MAX_RAND - 1)*randRHS.Random());
			row.push_back(rhs);
		}
		matrixRHS.push_back(row);
	}

	// generate matrix A
	for (int t = 0 ; t < nr; t++)
	{
		multimap<double, int > sortedRHS;
		for (int i = 0; i<m; i++)
		{
			 double rhs = matrixRHS[i][t];
			sortedRHS.insert(make_pair(rhs, i));
		}

		int cnt = 0;
		double h_k = 0;
		for (multimap<double, int>::reverse_iterator itr = sortedRHS.rbegin(); itr != sortedRHS.rend(); ++itr)
		{
			cnt ++;
			if (cnt == k +1)
			{
				h_k = itr->first;
				break;
			}
		}

		cnt = 0;
		for (multimap<double, int>::reverse_iterator itr = sortedRHS.rbegin(); itr != sortedRHS.rend(); ++itr)
		{
			if (cnt == k )
				break;
			rmatind = new int[n+1];
			rmatval = new double[n+1];
			for (int c =0; c<n; c++)
			{
				rmatind[c] = c;
				rmatval[c] = matrixA[t][c];
			}
			int i = itr->second;
			double h_i  = itr->first;
			rmatind[n] = n + i;
			rmatval[n] = h_i - h_k;

			double rhs[1];
			char   sense[1];
			rhs[0] =  h_i;
			sense[0] = 'G';

			int    rmatbeg[1];
			rmatbeg[0] = 0;

			string cstName = "SN"+ToString(i) + "R" + ToString(t);
			char* pStr = new char[20];
			strcpy(pStr, cstName.c_str());
			int status = CPXaddrows (cpxenv, cpxlp, 0, 1, n+1, rhs,
				sense, rmatbeg, rmatind, rmatval,
				NULL, &pStr);
			if ( status ) {
				fprintf (stderr, "Could not add new rows.\n");
				return;
			}
			delete []pStr;
			delete []rmatval;
			delete []rmatind;
			cnt++;
		}

	}
}

void PortfolioInstanceGenerator::addCardinalityConstraint()
{
	int    *rmatind = new int[m];
	double *rmatval = new double[m];
	
	// generate matrix A
	for (int i = 0; i<m; i++)
	{
		rmatind[i] = n + i;
		rmatval[i] = 1;

	}
	double rhs[1];
	char   sense[1];
	rhs[0] = k;
	sense[0] = 'L';

	int    rmatbeg[1];
	rmatbeg[0] = 0;

	string cstName ="CARD";
	char* pStr = new char[20];
	strcpy(pStr, cstName.c_str());
	int status = CPXaddrows (cpxenv, cpxlp, 0, 1, m, rhs,
		sense, rmatbeg, rmatind, rmatval,
		NULL, &pStr);
	if ( status ) {
		fprintf (stderr, "Could not add new rows.\n");
		return;
	}

	delete []rmatval;
	delete []rmatind;
}



void PortfolioInstanceGenerator::addResourceConstraint()
{
	int    *rmatind = new int[n];
	double *rmatval = new double[n];

	// generate matrix A
	for (int i = 0; i<n; i++)
	{
		rmatind[i] = i;
		rmatval[i] = 1;

	}
	double rhs[1];
	char   sense[1];
	rhs[0] = 1;
	sense[0] = 'E';

	int    rmatbeg[1];
	rmatbeg[0] = 0;

	string cstName ="RES";
	char* pStr = new char[20];
	strcpy(pStr, cstName.c_str());
	int status = CPXaddrows (cpxenv, cpxlp, 0, 1, n, rhs,
		sense, rmatbeg, rmatind, rmatval,
		NULL, &pStr);
	if ( status ) {
		fprintf (stderr, "Could not add new rows.\n");
		return;
	}

	delete []rmatval;
	delete []rmatind;
}




void PortfolioInstanceGenerator::buildRandomLHSModel(string modelName, OptModel::MODEL_TYPE modelType)
{
	initilizeCPXenv();
	initilizeCPXlp(modelName);

	createObjective(modelType);
	addScenarioConstraints();
	addCardinalityConstraint();
	string modelFileName = modelName + ".sav";
	CPXwriteprob(cpxenv, cpxlp, modelFileName.c_str(), NULL);
}


void PortfolioInstanceGenerator::buildRandomLHSModelWRes(string modelName, OptModel::MODEL_TYPE modelType)
{
	initilizeCPXenv();
	initilizeCPXlp(modelName);

	createRandomObjective(modelType);
	addScenarioConstraints();
	addResourceConstraint();
	addCardinalityConstraint();
	string modelFileName = modelName + ".sav";
	CPXwriteprob(cpxenv, cpxlp, modelFileName.c_str(), NULL);
}

void PortfolioInstanceGenerator::buildRandomLHSModelWResNormalDistribution(string modelName, OptModel::MODEL_TYPE modelType) //Jim's suggestion
{
	initilizeCPXenv();
	initilizeCPXlp(modelName);

	createObjectiveNormalDistribution(modelType);
	addScenarioConstraintsNormalDistribution();
	addResourceConstraint();
	addCardinalityConstraint();
	string modelFileName = modelName + ".sav";
	CPXwriteprob(cpxenv, cpxlp, modelFileName.c_str(), NULL);
}


void PortfolioInstanceGenerator::buildRandomRHSModel(string modelName, OptModel::MODEL_TYPE modelType)
{
	initilizeCPXenv();
	initilizeCPXlp(modelName);

	createObjective(modelType);
	addScenarioConstraintsRandomRHS();
	addCardinalityConstraint();
	string modelFileName = modelName + ".sav";
	CPXwriteprob(cpxenv, cpxlp, modelFileName.c_str(), NULL);
}

void PortfolioInstanceGenerator::buildTighterRandomRHSModel(string modelName, OptModel::MODEL_TYPE modelType)
{
	initilizeCPXenv();
	initilizeCPXlp(modelName);

	createObjective(modelType);
	addDominantScenarioConstraintsRandomRHS();
	addCardinalityConstraint();
	string modelFileName = modelName + ".sav";
	CPXwriteprob(cpxenv, cpxlp, modelFileName.c_str(), NULL);
}


