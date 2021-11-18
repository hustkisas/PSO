#include "Global.h"
#include "VaccinationInstanceGenerator.h"
#include "stocc.h"
#include "randomc.h"
#include "OptModel.h"

VaccinationInstanceGenerator::VaccinationInstanceGenerator(int _s, int _nr, int _n, int _m, int _k, double _epsilon)
:seed(_s),nr(_nr),n(_n),m(_m),k(_k), epsilon(_epsilon), avg_hh_size(0)
{
	bigM = 100000;
	randGen_eps = StochasticLib1(seed + 1);
	randGen_m = StochasticLib1(seed + 2);
	randGen_b = StochasticLib1(seed + 3);
	randGen_s = StochasticLib1(seed + 4);
	randGen_t = StochasticLib1(seed + 5);
	if (epsilon > SMALL_DECIMAL)
	{
		k = floor(m*epsilon);
		DataModel::k = k ;
	}
}

VaccinationInstanceGenerator::~VaccinationInstanceGenerator(void)
{
}

int VaccinationInstanceGenerator::initilizeCPXenv()
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


int VaccinationInstanceGenerator::initilizeCPXlp(string modelName)
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


bool VaccinationInstanceGenerator::readFamilyTypes(){
	fstream file;
	//read household data
	file.open("hh.txt", ios::in);
	if (!file){
		cerr << "Can not open hh.txt" <<  endl;
		exit(0);
	}

	string strLine;
	getline(file,strLine);
	if (strLine.empty())
	{
		cout<<endl<<"Wrong Data"<<endl;
		exit(0);
	}


	int cnt_f = 0; 
	int cnt_policy =0;
	avg_hh_size = 0;
	while (!strLine.empty())
	{
		int hhsize=0; int nChildren=0; int nAdults =0; int nElderly = 0;float prop = 0;
		istringstream lineStream(strLine) ;
		lineStream>>hhsize>>nChildren>>nAdults>>nElderly>>prop;
		HouseHold hh(cnt_f,hhsize,nChildren,nAdults,nElderly,prop);
		double test = hhsize*prop;
		avg_hh_size += test;
		for (int i =0; i <=nChildren; i++)
		{
			for (int j =0; j <=nAdults; j ++ )
			{
				for (int k =0; k<= nElderly; k++)
				{
					VacPolicy  p(cnt_policy,i,j,k,nChildren,nAdults,nElderly);
					vacPolicies.insert(make_pair(cnt_policy,p));
					hh.vecPolicies.push_back(p);
					cnt_policy++;
				}
			}
		}
		hhs.insert(make_pair(cnt_f,hh));
		cnt_f ++;
		getline(file,strLine);
	}
	n = vacPolicies.size();
	DataModel::n = n;
	return true;
}



bool VaccinationInstanceGenerator::createObjective(OptModel::MODEL_TYPE modelType){

	double *obj     = NULL;
	double *lb      = NULL;
	double *ub      = NULL;
	char   *ctype   = NULL;
	int    *rmatind = NULL;
	double *rmatval = NULL;
	char** colname = NULL;

	int    status = 0; 
	int num_cols = n + m ;

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

	int cnt_idx = 0;
	map<int, HouseHold>::iterator itr=hhs.begin();
	for ( ; itr != hhs.end(); ++itr)
	{
		HouseHold& hh = itr->second;
		float propotion_hh = hh.proportion;
		vector<VacPolicy>& vecPolicies = hh.vecPolicies;
		int cnt_p = 0;
		for (vector<VacPolicy>::iterator itrPolicy = vecPolicies.begin(); itrPolicy != vecPolicies.end(); itrPolicy++)
		{
			VacPolicy p = *itrPolicy;
			int pID = p.policyID;
			int nVacChildre = p.nVaccinatedChildren;
			int nVacAdults = p.nVaccinatedAdults;
			int nVacEldly = p.nVaccinatedElderly;

			obj[cnt_idx] = (nVacChildre + nVacAdults+ nVacEldly)*propotion_hh;                      // fill in obj coefficients
			string varname = "x["+ToString(hh.hhID)+"][" + ToString(cnt_p)+"]";
			colname[cnt_idx] = new char[12];
			strcpy(colname[cnt_idx], varname.c_str()); // make name for x
			ctype[cnt_idx] = 'C';
			(*itrPolicy).idx = cnt_idx;
			cnt_p++;
			cnt_idx++;

		}
	}
	for (int i = n; i < n+m ; i++)
	{
		obj[i] = 0;
		string varName = "z" + ToString( i-n);
		colname[i] = new char[8];
		strcpy(colname[i], varName.c_str()); // make name for z
		if (modelType == OptModel::MODEL_TYPE_LP)
			ctype[i] = 'C';
		else if( modelType == OptModel::MODEL_TYPE_MIP)
			ctype[i] = 'B';
	}


	for (int i =0; i < num_cols; i ++)       // fill in upper bound
	{
		if (i< n)
			ub[i] = 1;//CPX_INFBOUND;
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
    return true;
}


bool VaccinationInstanceGenerator::createEquaConstraints(){

	int    *rmatind = NULL;
	double *rmatval = NULL;
	/* adding rows to model */

	map<int,HouseHold>::iterator itr;
	for (itr = hhs.begin(); itr != hhs.end(); ++itr)
	{
		int totaNumCols = n;
		rmatind = (int * )   malloc (totaNumCols * sizeof(int));
		rmatval = (double *) malloc (totaNumCols * sizeof(double));
		double rhs[1];
		char   sense[1];

		HouseHold &hh = itr->second;
		vector<VacPolicy>& policies = hh.vecPolicies;
		vector<VacPolicy>::iterator itrP;
		int cnt_nz = 0;
		for (itrP = policies.begin(); itrP!= policies.end(); ++itrP)
		{
			VacPolicy& p = *itrP;
			rmatind[cnt_nz] = p.idx;
			rmatval[cnt_nz] = 1; 
			cnt_nz ++;
		}
		rhs[0] = 1;
		sense[0] = 'E';
		int    rmatbeg[1];
		rmatbeg[0] = 0;
		string cstName = "EQHH"+ToString(hh.hhID);
		char* pStr = new char[20];
		strcpy(pStr, cstName.c_str());
		int status = CPXaddrows (cpxenv, cpxlp, 0, 1, cnt_nz, rhs,
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


float VaccinationInstanceGenerator::get_u(){
	float t  = randGen_s.IRandom(0,10000)<=5000 ? 0.7 : 1.3;
	return t;
}

float VaccinationInstanceGenerator::get_s(){
	float t  = randGen_t.IRandom(0,10000)<=5000 ? 0.7 : 1.3;
	return t;
}

float VaccinationInstanceGenerator::get_a(HouseHold hh, VacPolicy p){
	double m = randGen_m.NormalTruncMy(1, 0.5, 0, 1000000000000000);
	double eps = randGen_eps.NormalTruncMy(0.85, 0.32, 0, 1);
	double b = randGen_b.NormalTruncMy(0.6, 0.32, 0, 1);
	double hf = hh.proportion;
	double u1 = get_u();
	double u2 = get_u();
	double u3 = get_u();
	double s1 = get_s();
	double s2 = get_s();
	double s3 = get_s();
	int f1 = p.nChildren;
	int f2 = p.nAdults;
	int f3 = p.nElderly;
	int v1 = p.nVaccinatedChildren;
	int v2 = p.nVaccinatedAdults;
	int v3 = p.nVaccinatedElderly;

	double term1 = 0;
	term1 = u1*s1*((1-b)*(f1 - v1*eps) + b*v1*eps*(1-eps))
		+ u2*s2*((1-b)*(f2 - v2*eps) + b*v2*eps*(1-eps))
		+u2*s2*((1-b)*(f3 - v3*eps) + b*v3*eps*(1-eps));
	double term2 = 0;
	term2 = b*
		(
		u1*s1*(f1 - v1*eps)*(f1 - v1*eps)
		+u2*s1*(f1 - v1*eps)*(f2 - v2*eps)
		+u3*s1*(f1 - v1*eps)*(f3 - v3*eps)
		+u1*s2*(f2 - v2*eps)*(f1 - v1*eps)
		+u2*s2*(f2 - v2*eps)*(f2 - v2*eps)
		+u3*s2*(f2 - v2*eps)*(f3 - v3*eps)
		+u1*s3*(f3 - v3*eps)*(f1 - v1*eps)
		+u2*s3*(f3 - v3*eps)*(f2 - v2*eps)
		+u3*s3*(f3 - v3*eps)*(f3 - v3*eps)
		);
	double product = m*hf*(term1 + term2)/avg_hh_size;
	return product;
}


bool  VaccinationInstanceGenerator::createSenarioConstraints(){
	
	int    *rmatind = NULL;
	double *rmatval = NULL;
	/* adding rows to model */
	
	for (int i = 0; i < m; i++)
	{
		int totaNumCols = n + 1;
		rmatind = (int * )   malloc (totaNumCols * sizeof(int));
		rmatval = (double *) malloc (totaNumCols * sizeof(double));
		double rhs[1];
		char   sense[1];

		map<int, HouseHold>::iterator itr;
		int cnt_nz = 0;
		for (itr = hhs.begin(); itr != hhs.end(); ++itr)
		{
			HouseHold& hh = itr->second;
			vector<VacPolicy>& policies = hh.vecPolicies;
			vector<VacPolicy>::iterator itrP;
			for (itrP=policies.begin(); itrP != policies.end(); ++itrP)
			{
				VacPolicy p = *itrP;
				int pID = p.policyID;
				float a = get_a(hh, p);
				rmatind[cnt_nz] = p.idx;
				rmatval[cnt_nz] = a; 
				cnt_nz ++;
			}
		}
		rmatind[cnt_nz] = n + i;
		rmatval[cnt_nz] = -bigM;
		rhs[0] = 1;
		sense[0] = 'L';
		int    rmatbeg[1];
		rmatbeg[0] = 0;
		string cstName = "SN"+ToString(i);
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


bool  VaccinationInstanceGenerator::createCardConstraints(){

	int    *rmatind = NULL;
	double *rmatval = NULL;
	/* adding rows to model */

	int totaNumCols = m;
	rmatind = (int * )   malloc (totaNumCols * sizeof(int));
	rmatval = (double *) malloc (totaNumCols * sizeof(double));
	double rhs[1];
	char   sense[1];

	for (int i = 0; i < m; i++)
	{
		rmatind[i] = n + i;
		rmatval[i] = 1; 
	}
	rhs[0] = k;
	sense[0] = 'L';
	int    rmatbeg[1];
	rmatbeg[0] = 0;
	string cstName = "CARD";
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

bool VaccinationInstanceGenerator::generateModel(string& modelName, OptModel::MODEL_TYPE modelType)
{
	initilizeCPXenv();
	initilizeCPXlp(modelName);
	readFamilyTypes();
	createObjective(modelType);
	createSenarioConstraints();
	createCardConstraints();
	createEquaConstraints();
	string modelFileName = modelName + ".sav";
	CPXwriteprob(cpxenv, cpxlp, modelFileName.c_str(), NULL);

// 	CPXmipopt(cpxenv, cpxlp);
// 	double obj_val;
// 	CPXgetobjval(cpxenv,cpxlp, &obj_val);

	complementModel();
	modelFileName = modelName+"c.sav";
	CPXwriteprob(cpxenv, cpxlp, modelFileName.c_str(), NULL);
// 	CPXmipopt(cpxenv, cpxlp);
// 	double b = getOriginalObjValueForComplementedModel(cpxenv, cpxlp);
	return true;
}

bool VaccinationInstanceGenerator::complementModel()
{
	// modify the objective function
// 	CPXchgobjsen(cpxenv, cpxlp, CPX_MAX);
	double* obj = new double[n];
	int* indicies = new int[n];
	CPXgetobj(cpxenv, cpxlp,obj, 0, n-1);
	for (int i = 0; i < n; i++)
	{
		indicies[i] = i;
		obj[i] = - obj[i];
	}
	CPXchgobj(cpxenv, cpxlp, n,indicies,obj );


	//modify the equation and inequality constraint
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
				return false;
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

	for (int i = 0; i < rNames.size(); i++)
	{
		int num_cols = CPXgetnumcols(cpxenv, cpxlp);
		int nzcnt;
		int * rmatbeg = new int(0);
		int * rmatind = new int[num_cols]();
		double * rmatval = new double[num_cols]();
		int rmatspace = num_cols;
		int surplus; 
		CPXgetrows(cpxenv, cpxlp, &nzcnt, rmatbeg, rmatind, rmatval, rmatspace, &surplus, i, i);

		string rowName = rNames[i].first;
		if (rowName.substr(0,2) == string("SN"))  // modify scenario constraint
		{
			double sum_a = 0;
			int idx_z = 0;
			for (int j = 0; j < nzcnt; j++)
			{
				if(rmatind[j] < n)             // a_fv
					sum_a += rmatval[j];
				else                                 // coefficient for z
					idx_z = rmatind[j];
			}
			double new_coef_z = sum_a -1;
			CPXchgcoef(cpxenv, cpxlp, i, idx_z, new_coef_z);
			CPXchgcoef(cpxenv, cpxlp, i, -1, new_coef_z);
			char sense = 'G';
			CPXchgsense(cpxenv, cpxlp, 1, &i, &sense);
		}
		else if (rowName.substr(0,2) == string("CA"))
		{
			continue;
		}
		else
		{
			CPXchgcoef(cpxenv, cpxlp, i, -1, nzcnt-1);
		}
	}
	return true;
}


double VaccinationInstanceGenerator::getOriginalObjValueForComplementedModel(CPXENVptr& env, CPXLPptr& lp)
{
	double* x = new double[n];
	double* c = new double[n];
	CPXgetx(env, lp, x,0, n-1);
	CPXgetobj(env, lp, c, 0, n-1);
	double obj=0;
	for (int i = 0; i < n; i++)
	{
		obj += -(1-x[i])*c[i];
	}
	delete[] x;
	delete[] c;
	return obj;
}


void VaccinationInstanceGenerator::ReadSMPS_StoFile(string& fileName)
{
	fstream file;
	file.open((fileName+string(".sto")).c_str(), ios::in);
	if (!file){
		cerr << "Can not open sto file, " <<  endl;
		return;
	}
	string strLine;
	double prob=0;
	vector<double> coefs;
	while (getline(file,strLine))
	{
		string str;
		istringstream buff(strLine);
		buff >> str;
		if(str == "SC")
		{
			buff >> str;
			buff >>str;
			buff >> prob;
			break;
		}
	}

	SCENARIO_ROW row;
	while (getline(file, strLine))
	{
		string str ;
		istringstream buff(strLine);
		buff >>str;
		if (str == "ENDATA")
		{
			scenarios.push_back(make_pair(prob, row) );
			row.clear();
			break;
		}
		if (str == "SC")
		{
			buff >> str;
			buff >>str;
			scenarios.push_back(make_pair(prob, row) );

			buff >> prob;
			row.clear();
			continue;
		}

		string str_col_idx = str; 
		string str_row_idx;
		double val;
		buff >>str_row_idx >> val;
// 		int idx_col = atoi(str_col_idx..c_str());
// 		int idx_row = atoi(str_row_idx.substr(1).c_str());
// 		pair<int, int> index(idx_col, idx_row);
		pair<string, string> str_index(str_col_idx, str_row_idx);
		row.push_back(make_pair(str_index, val));
	}

	DataModel::m = scenarios.size();
	DataModel::k = DataModel::m*epsilon;
}


void VaccinationInstanceGenerator::ReadSMPS_CorFile(string& fileName, OptModel::MODEL_TYPE modelType)
{
	OptModel model(cpxenv);
	int result= rename( (fileName + string(".cor")).c_str(), (fileName + string(".mps")).c_str());
	model.populateModelFromFile((fileName + string(".mps")).c_str());
	n = model.idx_X.size();
	m = scenarios.size();
	k = floor(m*epsilon);
	nr = 1;

	// adding columns
	double *obj     = NULL;
	double *lb      = NULL;
	double *ub      = NULL;
	char   *ctype   = NULL;
	int    *rmatind = NULL;
	double *rmatval = NULL;
	char** colname = NULL;

	int    status = 0; 
	int num_cols = scenarios.size() ;

	obj     = (double *) malloc (num_cols* sizeof(double));
	lb      = (double *) malloc (num_cols * sizeof(double));
	ub      = (double *) malloc (num_cols * sizeof(double));
	ctype   = (char *)   malloc (num_cols * sizeof(char));
	colname = (char**) malloc(num_cols*sizeof(char*));

	for (int i = 0; i< num_cols; i++)
	{
		obj[i] = 0;
		string strColName = "z" + ToString(i);
		colname[i] = new char[8];
		strcpy(colname[i], strColName.c_str()); // make name for z
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
		status = CPXnewcols (model.getCPXEnv(), model.getCPXLP(), num_cols, obj, NULL, ub, NULL, colname);
	}
	else if (modelType == OptModel::MODEL_TYPE_MIP)
	{
		status = CPXnewcols (model.getCPXEnv(), model.getCPXLP(), num_cols, obj, NULL, ub, ctype, colname);
	}
	if ( status ) {
		fprintf (stderr, "Could not add new columns.\n");
		return ;
	}
	
	free_and_null((char **)&obj);
	free_and_null((char **)&lb);
	free_and_null((char **)&ub);
	free_and_null((char **)&ctype);
	free_and_null((char **)&colname); 
	free_and_null((char **)&rmatind); 
	free_and_null((char **)&rmatval); 
	
	// add scenarios
	for (int i = 0; i < scenarios.size(); i++)
	{
		double rhs[1];
		char   sense[1];

		SCENARIO_ROW coefs = scenarios[i].second;
		int    *rmatind = new int[coefs.size() + 1 ];
		double *rmatval = new double[coefs.size() + 1 ];

		for (int j = 0; j < coefs.size(); j ++)
		{
			string col_name = coefs[j].first.first;
			int idx_col = model.var_name_idx_table.find(col_name)->second;
			string row_name = coefs[j].first.second;
			double val = coefs[j].second;
			model.removeConstraintFromCPLEX(row_name);
			rmatind[ j ] = idx_col;
			rmatval[ j ] = val;
		}
		rmatind[coefs.size()] = DataModel::n + i;
		rmatval[coefs.size()] = -bigM;
		rhs[0] = 1;
		sense[0] = 'L';
		int    rmatbeg[1];
		rmatbeg[0] = 0;
		string cstName = "SN"+ToString(i);
		char* pStr = new char[20];
		strcpy(pStr, cstName.c_str());
		int status = CPXaddrows (model.getCPXEnv(), model.getCPXLP(), 0, 1, coefs.size() + 1 , rhs,
			sense, rmatbeg, rmatind, rmatval,
			NULL, &pStr);
		if ( status ) {
			fprintf (stderr, "Could not add new rows.\n");
			return ;
		}
		delete rmatind; 
		delete rmatval; 
	}
	cpxlp = CPXcloneprob( cpxenv, model.getCPXLP(), &status);
}


void VaccinationInstanceGenerator::BuildModelFromSMPSFiles(string& fileName, OptModel::MODEL_TYPE modelType)
{
	initilizeCPXenv();
	initilizeCPXlp(fileName);
	ReadSMPS_StoFile(fileName);
	ReadSMPS_CorFile(fileName, modelType);
	createCardConstraints();
 	string modelFileName = fileName + ".sav";
// 	CPXwriteprob(cpxenv, cpxlp, modelFileName.c_str(), NULL);

//    	CPXlpopt(cpxenv, cpxlp);
//   	double obj_val;
// 	CPXgetobjval(cpxenv,cpxlp, &obj_val);
	complementModel();
// 	modelFileName = fileName+"c.lp";
	CPXwriteprob(cpxenv, cpxlp, modelFileName.c_str(), NULL);
//  	CPXlpopt(cpxenv, cpxlp);
//  	CPXgetobjval(cpxenv,cpxlp, &obj_val);
// 	double b = getOriginalObjValueForComplementedModel(cpxenv, cpxlp);

}


