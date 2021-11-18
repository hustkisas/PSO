#include "CutPool.h"

CutPool::CutPool(void)
{
}

CutPool::~CutPool(void)
{
}

void CutPool::setCutIdentifier(string str)
{
	strCutIndentifier = str;
}


int CutPool::updateCutPool(CPXENVptr& env, CPXLPptr& lp)
{
	int num_rows = CPXgetnumrows(env,  lp);
	int surplus;
	vector<string> rowNames;

	int status = CPXgetrowname (env, lp, NULL, NULL, 0, &surplus, 0,num_rows-1);
	int cur_rownamespace = - surplus;
	if ( cur_rownamespace > 0 ) {
		char          **cur_rowname = (char **) malloc (sizeof(char *)*(num_rows));
		char          *cur_rownamestore = (char *)  malloc (cur_rownamespace);
		if ( cur_rowname      == NULL ||
			cur_rownamestore == NULL   ) {
				fprintf (stderr, "Failed to get memory for column names.\n");
				status = -1;
				cout<< endl <<" error in printTableau"<<endl;
				return 0;
		}
		status = CPXgetrowname (env, lp, cur_rowname, cur_rownamestore, 
			cur_rownamespace, &surplus, 0, num_rows-1);
		// map name with column index
		for (int i = 0; i <num_rows ; i ++)
		{
			string rowName = cur_rowname[i];
			rowNames.push_back(rowName);
		}
		free_and_null( (char **)&cur_rowname);
 		free_and_null( (char **)&cur_rownamestore);

	}

	double *slacks = new double[num_rows];
	CPXgetslack(env, lp, slacks, 0, num_rows-1);
	
	for (int i =0; i < num_rows; i++)
	{
		string nameRowi = rowNames[i];
		if (nameRowi.find(strCutIndentifier) == string::npos)
			continue;
		
		double slack = slacks[i];
		map<string, PAIR_DOUBLE_INT>::iterator itr = cutTightness.find(nameRowi) ;
		if (itr == cutTightness.end())
		{
			int cnt = (fabs(slack) >= SMALL_DECIMAL); // =1 if not tight
			PAIR_DOUBLE_INT tightnessPair(slacks[i],cnt);
			cutTightness.insert(make_pair(nameRowi, tightnessPair));
		}
		else
		{
			itr->second.first = slack;
			itr->second.second = fabs(slack) <= SMALL_DECIMAL ? 0 :  itr->second.second+1;
		}
	}

	return 0;
}

void CutPool::displayCutTightness()
{
	cout << endl << "Cut Tightness: " <<endl;

	map<string, PAIR_DOUBLE_INT>::iterator itr; 
	for ( itr = cutTightness.begin() ; itr != cutTightness.end(); ++itr)
	{
		cout << itr->first << "\t current slack: " << itr->second.first << "\t # iters not tight: " << itr->second.second<<endl;
	}
	cout << endl;
}

int CutPool::removeCutFromCPLEX(CPXENVptr& env, CPXLPptr& lp, string& cutName)
{
	int num_rows = CPXgetnumrows(env,  lp);
	int surplus;
	vector<string> rowNames;

	int status = CPXgetrowname (env, lp, NULL, NULL, 0, &surplus, 0,num_rows-1);
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
		status = CPXgetrowname (env, lp, cur_rowname, cur_rownamestore, 
			cur_rownamespace, &surplus, 0, num_rows-1);
		// map name with column index
		for (int i = 0; i <num_rows ; i ++)
		{
			string rowName = cur_rowname[i];
			if(rowName != cutName)
				continue;

			// remove from CPLEX model
			int status = CPXdelrows (env, lp, i, i);
			if (status ==0)
			{
// 				cout << rowName << " was successfully removed"<< endl;
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

int CutPool::removeInactiveCuts(CPXENVptr& env, CPXLPptr& lp, int numInactiveRounds)
{
	vector<string> inactiveCutNames;
	for (map<string, PAIR_DOUBLE_INT>::iterator itr  = cutTightness.begin(); itr != cutTightness.end(); ++itr)
	{
		if (itr->second.second >= numInactiveRounds )
		{
			inactiveCutNames.push_back(itr->first);
		}
	}

	for (int i = 0; i < inactiveCutNames.size(); i++)
	{
		string name = inactiveCutNames[i];
		removeCutFromCPLEX(env, lp, name);
		cutTightness.erase(name);
		cuts.erase(name);
	}

	return inactiveCutNames.size();
}


int CutPool::removeInactiveCutsFromPool(int numInactiveRounds)
{
	vector<string> inactiveCutNames;
	for (map<string, PAIR_DOUBLE_INT>::iterator itr  = cutTightness.begin(); itr != cutTightness.end(); ++itr)
	{
		if (itr->second.second >= numInactiveRounds )
		{
			inactiveCutNames.push_back(itr->first);
		}
	}

	for (int i = 0; i < inactiveCutNames.size(); i++)
	{
		string name = inactiveCutNames[i];
// 		cutTightness.erase(name);
		cuts.erase(name);
	}

	return inactiveCutNames.size();
}



int CutPool::removeInactiveCutsRandomly(CPXENVptr& env, CPXLPptr& lp, int numInactiveRounds)
{
	vector<string> inactiveCutNames;
	for (map<string, PAIR_DOUBLE_INT>::iterator itr  = cutTightness.begin(); itr != cutTightness.end(); ++itr)
	{
		if (itr->second.second >= numInactiveRounds )
		{
			inactiveCutNames.push_back(itr->first);
		}
	}

	srand(cutTightness.size());
	for (int i = 0; i < inactiveCutNames.size(); i++)
	{
		if(rand()%20 ==1)
			continue;
		string name = inactiveCutNames[i];
		removeCutFromCPLEX(env, lp, name);
		cutTightness.erase(name);
	}

	return inactiveCutNames.size();
}

int CutPool::extractCutsFromModel(CPXENVptr& env, CPXLPptr& lp)
{
	int num_rows = CPXgetnumrows(env,  lp);
	int num_cols = CPXgetnumcols(env, lp);

	int surplus;
	vector<string> rowNames;

	int status = CPXgetrowname (env, lp, NULL, NULL, 0, &surplus, 0,num_rows-1);
	int cur_rownamespace = - surplus;
	if ( cur_rownamespace > 0 ) {
		char          **cur_rowname = (char **) malloc (sizeof(char *)*(num_rows));
		char          *cur_rownamestore = (char *)  malloc (cur_rownamespace);
		if ( cur_rowname      == NULL ||
			cur_rownamestore == NULL   ) {
				fprintf (stderr, "Failed to get memory for column names.\n");
				status = -1;
				cout<< endl <<" error in printTableau"<<endl;
				return 0;
		}
		status = CPXgetrowname (env, lp, cur_rowname, cur_rownamestore, 
			cur_rownamespace, &surplus, 0, num_rows-1);
		// map name with column index
		for (int i = 0; i <num_rows ; i ++)
		{
			string rowName = cur_rowname[i];
			rowNames.push_back(rowName);
		}
		free_and_null( (char **)&cur_rowname);
		/*free_and_null( (char **)&cur_rownamespace);*/
	}

	int nzcnt;
	int * rmatbeg = new int(0);
	int * rmatind = new int[num_cols]();
	double * rmatval = new double[num_cols]();
	for (int i =0; i < num_rows; i++)
	{
		string nameRowi = rowNames[i];
		if (nameRowi.find(strCutIndentifier) == string::npos)
			continue;
		
		int rmatspace = num_cols;
		int surplus; 
		double rhs;

		CPXgetrows(env, lp, &nzcnt, rmatbeg, rmatind, rmatval, rmatspace, &surplus, i, i);
		CPXgetrhs(env, lp, &rhs, i, i) ;
		CUT_COEFICIENTS coefs;
		for (int j = 0; j < nzcnt; j++)
		{
			coefs.push_back(pair<int, double>(rmatind[j], rmatval[j]));
		}
		CUT cut(coefs, rhs);
		cuts.insert(make_pair(nameRowi,cut));
	}
	delete rmatbeg;
	delete []rmatind;
	delete []rmatval;

	return cuts.size();
}




int CutPool::getNumCuts()
{
	return cuts.size();
}

