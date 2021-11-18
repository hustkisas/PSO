#include "TransportationInstanceGenerator.h"
#include "GlobalFunctions.h"

TransportationInstanceGenerator::TransportationInstanceGenerator(void)
{
}

TransportationInstanceGenerator::~TransportationInstanceGenerator(void)
{
}

void TransportationInstanceGenerator::readDataFile(string& fileName){

	fstream file;
	file.open(fileName.c_str(), ios::in);

	if (!file){
		cerr << "Can not open parameter file, use the default value" <<  endl;
	}
	else{
		cout<<"file opened successfully!"<<endl;
		string strLine;
		getline(file,strLine);
		if (strLine.empty())
		{
			cout<<endl<<"Wrong Data"<<endl;
			exit(0);
		}

		getline(file,strLine);
		istringstream is_epsilon(strLine) ;
		is_epsilon >>  epsilon;

		getline(file,strLine);
		istringstream is_num_supplier(strLine) ;
		is_num_supplier >> num_supplier;

		getline(file,strLine);
		istringstream is_num_demand(strLine) ;
		is_num_demand  >> num_demand;

		getline(file,strLine);
		istringstream is_num_scenario(strLine) ;
		is_num_scenario  >> num_scenario;

		string front_delim = string("[");
		string end_delim = string("]");
		vector<string> sections = ReadSectionsFromFile(fileName,front_delim ,end_delim );

		if (fileName.substr(0,1) == string("g"))  // if the instance is general probability
		{
			string t = sections[0];
			istringstream is_probability(sections[0]); 
			while (getline(is_probability,strLine,','))
			{
				double p = 0;
				istringstream buf(strLine);
				buf>> p;
				vecScenarioProbability.push_back(p);
			}
			sections.erase(sections.begin());
		}

		// read capacity
		string t = sections[0];
		istringstream is_capactiy(t);
		while (getline(is_capactiy,strLine,','))
		{
			int p = 0;
			istringstream buf(strLine);
			buf>> p;
			vecSupplierCapacity.push_back(p);
		}

		// read shipping cost 
		vector<string> strCostsSuppliers = ReadSectionsFromString(sections[1],front_delim, end_delim);
		for (int i = 0; i < strCostsSuppliers.size(); i++)
		{
			string strCostPerSupplier = strCostsSuppliers[i];
			istringstream is_costPerSupplier(strCostPerSupplier);
			vector<int> vecCosts;
			while (getline(is_costPerSupplier, strLine, ','))
			{
				int c = 0;
				istringstream buff(strLine);
				buff >> c ;
				vecCosts.push_back(c);
			}
			vecShippingCosts.push_back(vecCosts);
	  }

	// read demand
		vector<string> strDemands = ReadSectionsFromString(sections[2],front_delim, end_delim);
		for (int i = 0; i < strDemands.size(); i++)
		{
			string strDemand = strDemands[i];
			istringstream is_demand(strDemand);
			vector<int> vecDemands;
			while (getline(is_demand, strLine, ','))
			{
				int d = 0;
				istringstream buff(strLine);
				buff >> d ;
				vecDemands.push_back(d);
			}
			vecDemandScenarios.push_back(vecDemands);
		}
	}
}


int TransportationInstanceGenerator::initilizeCPXenv()
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


int TransportationInstanceGenerator::initilizeCPXlp(string modelName)
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



void TransportationInstanceGenerator::buildModel(string& modelName, OptModel::MODEL_TYPE modelType)
{
	initilizeCPXenv();
	initilizeCPXlp(modelName);
	createObjective(modelType);
	addScenarioConstraints();
	addCapacityConstraints();
	addCardinalityConstraint();
//	CPXmipopt(cpxenv, cpxlp);
	string modelFileName = modelName + ".sav";
	CPXwriteprob(cpxenv, cpxlp, modelFileName.c_str(), NULL);
}



void TransportationInstanceGenerator::buildTighterModel(string& modelName, OptModel::MODEL_TYPE modelType)
{
	initilizeCPXenv();
	initilizeCPXlp(modelName);
	createObjective(modelType);
	addScenarioDominatingConstraints();
	addCapacityConstraints();
	addCardinalityConstraint();
	//	CPXmipopt(cpxenv, cpxlp);
	string modelFileName = modelName + ".sav";
	CPXwriteprob(cpxenv, cpxlp, modelFileName.c_str(), NULL);
}


bool TransportationInstanceGenerator::createObjective(OptModel::MODEL_TYPE modelType){

	double *obj     = NULL;
	double *lb      = NULL;
	double *ub      = NULL;
	char   *ctype   = NULL;
	int    *rmatind = NULL;
	double *rmatval = NULL;
	char** colname = NULL;

	int    status = 0; 
	int num_cols = num_supplier*num_demand + num_scenario ;

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

	for (int s =0; s < num_supplier; s++)
	{
		vector<int> vecCost = vecShippingCosts[s];
		for (int d = 0; d < vecCost.size(); d++)
		{
			double cost = vecCost[d];
			obj[s*num_demand + d] = cost;
			string strColName = "x["+ ToString(s)+"][" +ToString(d) +"]";
			colname[s*num_demand + d] = new char[8];
			strcpy(colname[s*num_demand + d], strColName.c_str()); // make name for x
			ctype[s*num_demand + d] = 'C';
		}
	}

	for (int i = 0; i< num_scenario; i++)
	{
		obj[num_supplier*num_demand+i] = 0;
		string strColName = "z" + ToString(i);
		colname[num_supplier*num_demand+i] = new char[8];
		strcpy(colname[num_supplier*num_demand+i], strColName.c_str()); // make name for z
		if (modelType == OptModel::MODEL_TYPE_LP)
			ctype[num_supplier*num_demand+i] = 'C';
		else if( modelType == OptModel::MODEL_TYPE_MIP)
			ctype[num_supplier*num_demand+i] = 'B';
	}

	for (int i =0; i < num_cols; i ++)       // fill in upper bound
	{
		if (i<num_supplier*num_demand )
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

void TransportationInstanceGenerator::addScenarioConstraints()
{
	int    *rmatind = NULL;
	double *rmatval = NULL;

	for (int i = 0; i<num_scenario; i++)
	{
		for (int d = 0 ; d < num_demand; d++)
		{
			rmatind = new int[num_supplier+1];
			rmatval = new double[num_supplier+1];
			for (int s=0; s<num_supplier; s++)
			{
				rmatind[s] = num_demand*s + d;
				rmatval[s] = 1;
			}
			rmatind[num_supplier] = num_demand*num_supplier + i;
			rmatval[num_supplier] = vecDemandScenarios[i][d];

			double rhs[1];
			char   sense[1];
			rhs[0] =  vecDemandScenarios[i][d];
			sense[0] = 'G';

			int    rmatbeg[1];
			rmatbeg[0] = 0;

			string cstName = "SN"+ToString(i) + "R" + ToString(d);
			char* pStr = new char[20];
			strcpy(pStr, cstName.c_str());
			int status = CPXaddrows (cpxenv, cpxlp, 0, 1, num_supplier+1, rhs,
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


void TransportationInstanceGenerator::addScenarioDominatingConstraints()
{
	int    *rmatind = NULL;
	double *rmatval = NULL;

	for (int d = 0 ; d < num_demand; d++)
	{
		multimap<double, int > sortedRHS;
		for (int i = 0; i<num_scenario; i++)
		{
			sortedRHS.insert(make_pair(vecDemandScenarios[i][d], i));
		}

		int cnt = 0;
		double h_k = 0;
		int k = floor(epsilon*num_scenario);
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
			rmatind = new int[num_supplier+1];
			rmatval = new double[num_supplier+1];
			for (int s=0; s<num_supplier; s++)
			{
				rmatind[s] = num_demand*s + d;
				rmatval[s] = 1;
			}
			int i = itr->second;
			double h_i  = itr->first;
			rmatind[num_supplier] = num_demand*num_supplier + i;
			rmatval[num_supplier] = h_i - h_k;

			double rhs[1];
			char   sense[1];
			rhs[0] =  h_i;
			sense[0] = 'G';

			int    rmatbeg[1];
			rmatbeg[0] = 0;

			string cstName = "SN"+ToString(i) + "R" + ToString(d);
			char* pStr = new char[20];
			strcpy(pStr, cstName.c_str());
			int status = CPXaddrows (cpxenv, cpxlp, 0, 1, num_supplier+1, rhs,
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


void TransportationInstanceGenerator::addCapacityConstraints()
{
	int    *rmatind = NULL;
	double *rmatval = NULL;

	for (int i = 0; i<num_supplier; i++)
	{
		rmatind = new int[num_demand];
		rmatval = new double[num_demand];
		for (int d = 0 ; d < num_demand; d++)
		{
			rmatind[d] = num_demand*i + d;
			rmatval[d] = 1;
		}
		double rhs[1];
		char   sense[1];
		rhs[0] =  vecSupplierCapacity[i];
		sense[0] = 'L';

		int    rmatbeg[1];
		rmatbeg[0] = 0;

		string cstName = "CAP"+ToString(i);
		char* pStr = new char[20];
		strcpy(pStr, cstName.c_str());
		int status = CPXaddrows (cpxenv, cpxlp, 0, 1, num_supplier+1, rhs,
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



void TransportationInstanceGenerator::addCardinalityConstraint()
{
	int    *rmatind = new int[num_scenario];
	double *rmatval = new double[num_scenario];

	// generate matrix A
	for (int i = 0; i<num_scenario; i++)
	{
		rmatind[i] = num_supplier*num_demand + i;
		rmatval[i] = 1;

	}
	double rhs[1];
	char   sense[1];
	rhs[0] = floor(epsilon*num_scenario);
	sense[0] = 'L';

	int    rmatbeg[1];
	rmatbeg[0] = 0;

	string cstName ="CARD";
	char* pStr = new char[20];
	strcpy(pStr, cstName.c_str());
	int status = CPXaddrows (cpxenv, cpxlp, 0, 1, num_scenario, rhs,
		sense, rmatbeg, rmatind, rmatval,
		NULL, &pStr);
	if ( status ) {
		fprintf (stderr, "Could not add new rows.\n");
		return;
	}

	delete []rmatval;
	delete []rmatind;
}

