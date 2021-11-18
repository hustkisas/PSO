#include "WolseyStrengthening.h"
#include "Global.h"
#include "OptModel.h"
#include "TwoConstrainSetCutGenerator.h"
#include "ResourceCutGenerator.h"
#include "PortfolioInstanceGenerator.h"

extern string CAT_NAME;
int WolseyStrengthening::iteration_cnt = 0;
vector<double> WolseyStrengthening::beta= vector<double>(2*DataModel::m, 0);

WolseyStrengthening::WolseyStrengthening(OptModel* _model)
{
	model_p = _model;
	WolseyStrengthening::beta= vector<double>(DataModel::m, 0);
	timeLimit = 1/SMALL_DECIMAL;
	timeElapsed = 0;

}

void WolseyStrengthening::strengthenModel(int method)
{
	min_model = new OptModel(*model_p);
	if (method == 3)
	{
		TwoConstrainSetCutGenerator cutGen(min_model);
		cutGen.addAllPairFacets(OptModel::CUT_AS_FORMULATION);
	}
	LHS lhs = model_p->getObjective();
	RHS rhs = model_p->LP_optimal_value;
	CUT cut(lhs,rhs);
	string cutName = "objcut";
	min_model->addAConstToCPXModel(cut, cutName, 0);

	int num_rows = model_p->getNumRows();
	int surplus;
	CPXENVptr& cpxenv =model_p->getCPXEnv();
	CPXLPptr&  cpxlp = model_p->getCPXLP();
	double lp_optimal = model_p->LP_optimal_value;
	int num_cols = CPXgetnumcols(cpxenv, cpxlp);
	vector< double > gs;
	vector< int > list_i;
	vector< int > list_j;
	vector<double> list_rhs;

	int status = CPXgetrowname (cpxenv, cpxlp, NULL, NULL, 0, &surplus, 0,num_rows-1);
	vector<string> rowsToDelete;
	vector<string> varsToDelete;
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
			if (t != string::npos)
				rowName.erase(t);
			if (rowName != "SN")
				continue;

			int nzcnt;
			int * rmatbeg = new int(0);
			int * rmatind = new int[num_cols]();
			double * rmatval = new double[num_cols]();
			int rmatspace = num_cols;
			int surplus; 
			double rhs;

			CPXgetrows(cpxenv, cpxlp, &nzcnt, rmatbeg, rmatind, rmatval, rmatspace, &surplus, r, r);
			CPXgetrhs(cpxenv, cpxlp, &rhs, r, r);
			int index_z = 0;
			for (int j = 0 ; j < nzcnt; j++)
			{
				if (rmatind[j] >=  DataModel::n)
				{
					index_z = rmatind[j];
					break;
				}
			}
			double g_i ;
			if (method == 0)
				g_i = calculateBetaBySolvingMIP(nzcnt,  rmatind, rmatval, lp_optimal); // 
			else if(method == 1 || method == 3 )
				g_i = calculateBetaBySolvingMIP_rev1(nzcnt,  rmatind, rmatval, lp_optimal);
			else if(method == 2)
				g_i = calculateBetaBySolvingMIP_2ConstSet(nzcnt,  rmatind, rmatval, lp_optimal);
			gs.push_back(g_i);
			list_i.push_back(r);
			list_j.push_back(index_z);
			list_rhs.push_back(rhs);
			if(rhs - g_i < 0)
			{
 				rowsToDelete.push_back(cpyRowName);
				map<int, string>::iterator itr  = model_p->idx_var_name_table.find(index_z);
				if (itr!= model_p->idx_var_name_table.end())
				{
					varsToDelete.push_back(itr->second);
				}
			}
			delete rmatbeg;
			delete []rmatind;
			delete []rmatval;
		}
		free_and_null( (char **)&cur_rowname);
		free_and_null((char **)&cur_rownamestore);
	}

	// change the coeffients
	int * rowlist = new int[list_i.size()];
	int * collist = new int[list_j.size()];
	double * vallist = new double[gs.size()];
	for (int k = 0 ; k < gs.size(); k++)
	{
		rowlist[k] = list_i[k];
		collist[k] = list_j[k];
		vallist[k] = (list_rhs[k] - gs[k])<=0 ? list_rhs[k] : list_rhs[k] - gs[k];
	}
	CPXchgcoeflist(cpxenv, cpxlp, gs.size(), rowlist, collist, vallist);
	delete [] rowlist;
	delete []collist;
	delete []vallist;

	for (map<int, OptModel*>::iterator itr = LBmodels.begin(); itr != LBmodels.end(); ++itr)
	{
		OptModel* p = itr->second;
		delete p;
	}
	LBmodels.clear();

}


void WolseyStrengthening::strengthenModelWiTheoreticalBound()
{
	int num_rows = model_p->getNumRows();
	int surplus;
	CPXENVptr& cpxenv =model_p->getCPXEnv();
	CPXLPptr&  cpxlp = model_p->getCPXLP();
	double lp_optimal = model_p->LP_optimal_value;
	int num_cols = CPXgetnumcols(cpxenv, cpxlp);
	vector< double > gs;
	vector< int > list_i;
	vector< int > list_j;
	vector<double> list_rhs;


	int status = CPXgetrowname (cpxenv, cpxlp, NULL, NULL, 0, &surplus, 0,num_rows-1);
	vector<string> rowsToDelete;
	vector<string> varsToDelete;
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

		for (int r = 0; r <num_rows ; r ++)
		{
			string rowName = cur_rowname[r];
			string cpyRowName = rowName;
			size_t t =rowName.find_first_of("0123456789");
			if (t != string::npos)
				rowName.erase(t);
			if (rowName != "SN")
				continue;

			int nzcnt;
			int * rmatbeg = new int(0);
			int * rmatind = new int[num_cols]();
			double * rmatval = new double[num_cols]();
			int rmatspace = num_cols;
			int surplus; 
			double rhs;

			map<double,int> sortedRo;
			CPXgetrows(cpxenv, cpxlp, &nzcnt, rmatbeg, rmatind, rmatval, rmatspace, &surplus, r, r);
			CPXgetrhs(cpxenv, cpxlp, &rhs, r, r);
			int index_z = 0;
			for (int j = 0 ; j < nzcnt; j++)
			{
				if (rmatind[j] >=  DataModel::n)
				{
					index_z = rmatind[j];
					break;
				}
				sortedRo.insert(make_pair(rmatval[j]/PortfolioInstanceGenerator::ro[j],0));
			}
			double g_i = sortedRo.begin()->first;
			g_i = g_i*((double)(DataModel::m-DataModel::k))/(DataModel::m-g_i*DataModel::k);
			gs.push_back(g_i);
			list_i.push_back(r);
			list_j.push_back(index_z);
			list_rhs.push_back(rhs);
			if(rhs - g_i < 0)
			{
 				rowsToDelete.push_back(cpyRowName);
				map<int, string>::iterator itr  = model_p->idx_var_name_table.find(index_z);
				if (itr!= model_p->idx_var_name_table.end())
				{
					varsToDelete.push_back(itr->second);
				}
			}
			delete rmatbeg;
			delete []rmatind;
			delete []rmatval;
		}
		free_and_null( (char **)&cur_rowname);
		free_and_null((char **)&cur_rownamestore);
	}

	// change the coeffients
	int * rowlist = new int[list_i.size()];
	int * collist = new int[list_j.size()];
	double * vallist = new double[gs.size()];
	for (int k = 0 ; k < gs.size(); k++)
	{
		rowlist[k] = list_i[k];
		collist[k] = list_j[k];
		vallist[k] = (list_rhs[k] - gs[k])<=0 ? list_rhs[k] : list_rhs[k] - gs[k];
	}
	CPXchgcoeflist(cpxenv, cpxlp, gs.size(), rowlist, collist, vallist);
	delete [] rowlist;
	delete []collist;
	delete []vallist;

}



void WolseyStrengthening::strengthenModelWith2ConstSet(int method)
{
	// build basic model with all 2-const set inequalities
	min_model = new OptModel(*model_p);
	TwoConstrainSetCutGenerator cutGen(min_model);
	cutGen.addAllPairFacets(OptModel::CUT_AS_FORMULATION);
	LHS lhs = model_p->getObjective();
	RHS rhs = model_p->LP_optimal_value;
	CUT cut(lhs,rhs);
	string cutName = "objcut";
	min_model->addAConstToCPXModel(cut, cutName, 0);

	int num_rows = model_p->getNumRows();
	int surplus;
	CPXENVptr& cpxenv =model_p->getCPXEnv();
	CPXLPptr&  cpxlp = model_p->getCPXLP();
	double lp_optimal = model_p->LP_optimal_value;
	int num_cols = CPXgetnumcols(cpxenv, cpxlp);
	vector< double > gs;
	vector< int > list_i;
	vector< int > list_j;
	vector<double> list_rhs;

	int status = CPXgetrowname (cpxenv, cpxlp, NULL, NULL, 0, &surplus, 0,num_rows-1);
	vector<string> rowsToDelete;
	vector<string> varsToDelete;
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
			if (t != string::npos)
				rowName.erase(t);
			if (rowName != "SN")
				continue;

			int nzcnt;
			int * rmatbeg = new int(0);
			int * rmatind = new int[num_cols]();
			double * rmatval = new double[num_cols]();
			int rmatspace = num_cols;
			int surplus; 
			double rhs;

			CPXgetrows(cpxenv, cpxlp, &nzcnt, rmatbeg, rmatind, rmatval, rmatspace, &surplus, r, r);
			CPXgetrhs(cpxenv, cpxlp, &rhs, r, r);
			int index_z = 0;
			for (int j = 0 ; j < nzcnt; j++)
			{
				if (rmatind[j] >=  DataModel::n)
				{
					index_z = rmatind[j];
					break;
				}
			}
			double g_i ;
			if (method == 0)
				g_i = calculateBetaBySolvingMIP(nzcnt,  rmatind, rmatval, lp_optimal); // 
			else if(method == 1)
				g_i = calculateBetaBySolvingMIP_rev1(nzcnt,  rmatind, rmatval, lp_optimal);
			else if(method == 2)
				g_i = calculateBetaBySolvingMIP_2ConstSet(nzcnt,  rmatind, rmatval, lp_optimal);
			gs.push_back(g_i);
			list_i.push_back(r);
			list_j.push_back(index_z);
			list_rhs.push_back(rhs);
			if(rhs - g_i < 0)
			{
				rowsToDelete.push_back(cpyRowName);
				map<int, string>::iterator itr  = model_p->idx_var_name_table.find(index_z);
				if (itr!= model_p->idx_var_name_table.end())
				{
					varsToDelete.push_back(itr->second);
				}
			}
			delete rmatbeg;
			delete []rmatind;
			delete []rmatval;
		}
		free_and_null( (char **)&cur_rowname);
		free_and_null((char **)&cur_rownamestore);
	}

	// change the coeffients
	int * rowlist = new int[list_i.size()];
	int * collist = new int[list_j.size()];
	double * vallist = new double[gs.size()];
	for (int k = 0 ; k < gs.size(); k++)
	{
		rowlist[k] = list_i[k];
		collist[k] = list_j[k];
		vallist[k] = (list_rhs[k] - gs[k])<=0 ? list_rhs[k] : list_rhs[k] - gs[k];
	}
	CPXchgcoeflist(cpxenv, cpxlp, gs.size(), rowlist, collist, vallist);
	delete [] rowlist;
	delete []collist;
	delete []vallist;

	for (map<int, OptModel*>::iterator itr = LBmodels.begin(); itr != LBmodels.end(); ++itr)
	{
		OptModel* p = itr->second;
		delete p;
	}
	
	delete min_model;
}


void WolseyStrengthening::strengthenModelWithResourceCut(int method)
{
	// build basic model with all 2-const set inequalities
	min_model = new OptModel(*model_p);
	ResourceCutGenerator cutGen(min_model);
	cutGen.addAllScenarioCuts();
	LHS lhs = model_p->getObjective();
	RHS rhs = model_p->LP_optimal_value;
	CUT cut(lhs,rhs);
	string cutName = "objcut";
	min_model->addAConstToCPXModel(cut, cutName, 0);

	int num_rows = model_p->getNumRows();
	int surplus;
	CPXENVptr& cpxenv =model_p->getCPXEnv();
	CPXLPptr&  cpxlp = model_p->getCPXLP();
	double lp_optimal = model_p->LP_optimal_value;
	int num_cols = CPXgetnumcols(cpxenv, cpxlp);
	vector< double > gs;
	vector< int > list_i;
	vector< int > list_j;
	vector<double> list_rhs;

	int status = CPXgetrowname (cpxenv, cpxlp, NULL, NULL, 0, &surplus, 0,num_rows-1);
	vector<string> rowsToDelete;
	vector<string> varsToDelete;
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
			if (t != string::npos)
				rowName.erase(t);
			if (rowName != "SN")
				continue;

			int nzcnt;
			int * rmatbeg = new int(0);
			int * rmatind = new int[num_cols]();
			double * rmatval = new double[num_cols]();
			int rmatspace = num_cols;
			int surplus; 
			double rhs;

			CPXgetrows(cpxenv, cpxlp, &nzcnt, rmatbeg, rmatind, rmatval, rmatspace, &surplus, r, r);
			CPXgetrhs(cpxenv, cpxlp, &rhs, r, r);
			int index_z = 0;
			for (int j = 0 ; j < nzcnt; j++)
			{
				if (rmatind[j] >=  DataModel::n)
				{
					index_z = rmatind[j];
					break;
				}
			}
			double g_i ;
			if (method == 0)
				g_i = calculateBetaBySolvingMIP(nzcnt,  rmatind, rmatval, lp_optimal); // 
			else if(method == 1)
				g_i = calculateBetaBySolvingMIP_rev1(nzcnt,  rmatind, rmatval, lp_optimal);
// 			else if(method == 2)
// 				g_i = calculateBetaBySolvingMIP_2ConstSet(nzcnt,  rmatind, rmatval, lp_optimal);
			gs.push_back(g_i);
			list_i.push_back(r);
			list_j.push_back(index_z);
			list_rhs.push_back(rhs);
			if(rhs - g_i < 0)
			{
				rowsToDelete.push_back(cpyRowName);
				map<int, string>::iterator itr  = model_p->idx_var_name_table.find(index_z);
				if (itr!= model_p->idx_var_name_table.end())
				{
					varsToDelete.push_back(itr->second);
				}
			}
			delete rmatbeg;
			delete []rmatind;
			delete []rmatval;
		}
		free_and_null( (char **)&cur_rowname);
		free_and_null((char **)&cur_rownamestore);
	}

	// change the coeffients
	int * rowlist = new int[list_i.size()];
	int * collist = new int[list_j.size()];
	double * vallist = new double[gs.size()];
	for (int k = 0 ; k < gs.size(); k++)
	{
		rowlist[k] = list_i[k];
		collist[k] = list_j[k];
		vallist[k] = (list_rhs[k] - gs[k])<=0 ? list_rhs[k] : list_rhs[k] - gs[k];
	}
	CPXchgcoeflist(cpxenv, cpxlp, gs.size(), rowlist, collist, vallist);
	delete [] rowlist;
	delete []collist;
	delete []vallist;

	for (map<int, OptModel*>::iterator itr = LBmodels.begin(); itr != LBmodels.end(); ++itr)
	{
		OptModel* p = itr->second;
		delete p;
	}

	delete min_model;
}


void WolseyStrengthening::strengthenModelWithResourceCut2SetCut(int method)
{
	// build basic model with all 2-const set inequalities
	min_model = new OptModel(*model_p);
	ResourceCutGenerator cutGen(min_model);
	cutGen.addAllScenarioCuts();
	TwoConstrainSetCutGenerator twoConstSetGen(min_model);
	twoConstSetGen.addAllPairFacets(OptModel::CUT_AS_FORMULATION);

	LHS lhs = model_p->getObjective();
	RHS rhs = model_p->LP_optimal_value;
	CUT cut(lhs,rhs);
	string cutName = "objcut";
	min_model->addAConstToCPXModel(cut, cutName, 0);

	int num_rows = model_p->getNumRows();
	int surplus;
	CPXENVptr& cpxenv =model_p->getCPXEnv();
	CPXLPptr&  cpxlp = model_p->getCPXLP();
	double lp_optimal = model_p->LP_optimal_value;
	int num_cols = CPXgetnumcols(cpxenv, cpxlp);
	vector< double > gs;
	vector< int > list_i;
	vector< int > list_j;
	vector<double> list_rhs;

	int status = CPXgetrowname (cpxenv, cpxlp, NULL, NULL, 0, &surplus, 0,num_rows-1);
	vector<string> rowsToDelete;
	vector<string> varsToDelete;
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
			if (t != string::npos)
				rowName.erase(t);
			if (rowName != "SN")
				continue;

			int nzcnt;
			int * rmatbeg = new int(0);
			int * rmatind = new int[num_cols]();
			double * rmatval = new double[num_cols]();
			int rmatspace = num_cols;
			int surplus; 
			double rhs;

			CPXgetrows(cpxenv, cpxlp, &nzcnt, rmatbeg, rmatind, rmatval, rmatspace, &surplus, r, r);
			CPXgetrhs(cpxenv, cpxlp, &rhs, r, r);
			int index_z = 0;
			for (int j = 0 ; j < nzcnt; j++)
			{
				if (rmatind[j] >=  DataModel::n)
				{
					index_z = rmatind[j];
					break;
				}
			}
			double g_i ;
			if (method == 0)
				g_i = calculateBetaBySolvingMIP(nzcnt,  rmatind, rmatval, lp_optimal); // 
			else if(method == 1)
				g_i = calculateBetaBySolvingMIP_rev1(nzcnt,  rmatind, rmatval, lp_optimal);
			// 			else if(method == 2)
			// 				g_i = calculateBetaBySolvingMIP_2ConstSet(nzcnt,  rmatind, rmatval, lp_optimal);
			gs.push_back(g_i);
			list_i.push_back(r);
			list_j.push_back(index_z);
			list_rhs.push_back(rhs);
			if(rhs - g_i < 0)
			{
				rowsToDelete.push_back(cpyRowName);
				map<int, string>::iterator itr  = model_p->idx_var_name_table.find(index_z);
				if (itr!= model_p->idx_var_name_table.end())
				{
					varsToDelete.push_back(itr->second);
				}
			}
			delete rmatbeg;
			delete []rmatind;
			delete []rmatval;
		}
		free_and_null( (char **)&cur_rowname);
		free_and_null((char **)&cur_rownamestore);
	}

	// change the coeffients
	int * rowlist = new int[list_i.size()];
	int * collist = new int[list_j.size()];
	double * vallist = new double[gs.size()];
	for (int k = 0 ; k < gs.size(); k++)
	{
		rowlist[k] = list_i[k];
		collist[k] = list_j[k];
		vallist[k] = (list_rhs[k] - gs[k])<=0 ? list_rhs[k] : list_rhs[k] - gs[k];
	}
	CPXchgcoeflist(cpxenv, cpxlp, gs.size(), rowlist, collist, vallist);
	delete [] rowlist;
	delete []collist;
	delete []vallist;

	for (map<int, OptModel*>::iterator itr = LBmodels.begin(); itr != LBmodels.end(); ++itr)
	{
		OptModel* p = itr->second;
		delete p;
	}

	delete min_model;
}


void WolseyStrengthening::strengthenModelNew( int method)
{
	// build basic model with all 2-const set inequalities
	min_model = new OptModel(*model_p);
	ResourceCutGenerator cutGen(min_model);
	cutGen.addAllScenarioCuts();
	
	if (method == 1)
	{
		TwoConstrainSetCutGenerator twoConstSetGen(min_model);
		for (int i =0; i < 10; i++)
		{
			if(twoConstSetGen.addViolatedFacets())
				break;
		}
	}

 	TwoConstrainSetCutGenerator twoConstSetGen(min_model);
 	twoConstSetGen.addAllPairFacets(OptModel::CUT_AS_FORMULATION);

	LHS lhs = model_p->getObjective();
	RHS rhs = model_p->LP_optimal_value;
	CUT cut(lhs,rhs);
	string cutName = "objcut";
	min_model->addAConstToCPXModel(cut, cutName, 0);

	int num_rows = model_p->getNumRows();
	int surplus;
	CPXENVptr& cpxenv =model_p->getCPXEnv();
	CPXLPptr&  cpxlp = model_p->getCPXLP();
	double lp_optimal = model_p->LP_optimal_value;
	int num_cols = CPXgetnumcols(cpxenv, cpxlp);
	vector< double > gs;
	vector< int > list_i;
	vector< int > list_j;
	vector<double> list_rhs;

	int status = CPXgetrowname (cpxenv, cpxlp, NULL, NULL, 0, &surplus, 0,num_rows-1);
	vector<string> rowsToDelete;
	vector<string> varsToDelete;
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
			if (t != string::npos)
				rowName.erase(t);
			if (rowName != "SN")
				continue;

			int nzcnt;
			int * rmatbeg = new int(0);
			int * rmatind = new int[num_cols]();
			double * rmatval = new double[num_cols]();
			int rmatspace = num_cols;
			int surplus; 
			double rhs;

			CPXgetrows(cpxenv, cpxlp, &nzcnt, rmatbeg, rmatind, rmatval, rmatspace, &surplus, r, r);
			CPXgetrhs(cpxenv, cpxlp, &rhs, r, r);
			int index_z = 0;
			for (int j = 0 ; j < nzcnt; j++)
			{
				if (rmatind[j] >=  DataModel::n)
				{
					index_z = rmatind[j];
					break;
				}
			}
			double g_i ;
			if (method == 0)
				g_i = calculateBetaBySolvingMIP(nzcnt,  rmatind, rmatval, lp_optimal); // 
			else if(method == 1)
				g_i = calculateBetaBySolvingMIP_rev1(nzcnt,  rmatind, rmatval, lp_optimal);
			// 			else if(method == 2)
			// 				g_i = calculateBetaBySolvingMIP_2ConstSet(nzcnt,  rmatind, rmatval, lp_optimal);
			gs.push_back(g_i);
			list_i.push_back(r);
			list_j.push_back(index_z);
			list_rhs.push_back(rhs);
			if(rhs - g_i < 0)
			{
				rowsToDelete.push_back(cpyRowName);
				map<int, string>::iterator itr  = model_p->idx_var_name_table.find(index_z);
				if (itr!= model_p->idx_var_name_table.end())
				{
					varsToDelete.push_back(itr->second);
				}
			}
			delete rmatbeg;
			delete []rmatind;
			delete []rmatval;
		}
		free_and_null( (char **)&cur_rowname);
		free_and_null((char **)&cur_rownamestore);
	}

	// change the coeffients
	int * rowlist = new int[list_i.size()];
	int * collist = new int[list_j.size()];
	double * vallist = new double[gs.size()];
	for (int k = 0 ; k < gs.size(); k++)
	{
		rowlist[k] = list_i[k];
		collist[k] = list_j[k];
		vallist[k] = (list_rhs[k] - gs[k])<=0 ? list_rhs[k] : list_rhs[k] - gs[k];
	}
	CPXchgcoeflist(cpxenv, cpxlp, gs.size(), rowlist, collist, vallist);
	delete [] rowlist;
	delete []collist;
	delete []vallist;

	for (map<int, OptModel*>::iterator itr = LBmodels.begin(); itr != LBmodels.end(); ++itr)
	{
		OptModel* p = itr->second;
		delete p;
	}

	delete min_model;
}



void WolseyStrengthening::strengthenModelWithResourceCut2SetCutSeparation(int method)
{
	// first add all resource cuts, then solve each LP with violated 2-const cuts
	min_model = new OptModel(*model_p);
	ResourceCutGenerator cutGen(min_model);
	cutGen.addAllScenarioCuts();
	TwoConstrainSetCutGenerator twoConstSetGen(min_model);
	twoConstSetGen.addAllPairFacets(OptModel::CUT_AS_FORMULATION);

	LHS lhs = model_p->getObjective();
	RHS rhs = model_p->LP_optimal_value;
	CUT cut(lhs,rhs);
	string cutName = "objcut";
	min_model->addAConstToCPXModel(cut, cutName, 0);

	int num_rows = model_p->getNumRows();
	int surplus;
	CPXENVptr& cpxenv =model_p->getCPXEnv();
	CPXLPptr&  cpxlp = model_p->getCPXLP();
	double lp_optimal = model_p->LP_optimal_value;
	int num_cols = CPXgetnumcols(cpxenv, cpxlp);
	vector< double > gs;
	vector< int > list_i;
	vector< int > list_j;
	vector<double> list_rhs;

	int status = CPXgetrowname (cpxenv, cpxlp, NULL, NULL, 0, &surplus, 0,num_rows-1);
	vector<string> rowsToDelete;
	vector<string> varsToDelete;
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
			if (t != string::npos)
				rowName.erase(t);
			if (rowName != "SN")
				continue;

			int nzcnt;
			int * rmatbeg = new int(0);
			int * rmatind = new int[num_cols]();
			double * rmatval = new double[num_cols]();
			int rmatspace = num_cols;
			int surplus; 
			double rhs;

			CPXgetrows(cpxenv, cpxlp, &nzcnt, rmatbeg, rmatind, rmatval, rmatspace, &surplus, r, r);
			CPXgetrhs(cpxenv, cpxlp, &rhs, r, r);
			int index_z = 0;
			for (int j = 0 ; j < nzcnt; j++)
			{
				if (rmatind[j] >=  DataModel::n)
				{
					index_z = rmatind[j];
					break;
				}
			}
			double g_i ;
			if (method == 0)
				g_i = calculateBetaBySolvingMIP(nzcnt,  rmatind, rmatval, lp_optimal); // 
			else if(method == 1)
				g_i = calculateBetaBySolvingMIP_rev2_addingViolated2ConstCuts(nzcnt,  rmatind, rmatval, lp_optimal);
			// 			else if(method == 2)
			// 				g_i = calculateBetaBySolvingMIP_2ConstSet(nzcnt,  rmatind, rmatval, lp_optimal);
			gs.push_back(g_i);
			list_i.push_back(r);
			list_j.push_back(index_z);
			list_rhs.push_back(rhs);
			if(rhs - g_i < 0)
			{
				rowsToDelete.push_back(cpyRowName);
				map<int, string>::iterator itr  = model_p->idx_var_name_table.find(index_z);
				if (itr!= model_p->idx_var_name_table.end())
				{
					varsToDelete.push_back(itr->second);
				}
			}
			delete rmatbeg;
			delete []rmatind;
			delete []rmatval;
		}
		free_and_null( (char **)&cur_rowname);
		free_and_null((char **)&cur_rownamestore);
	}

	// change the coeffients
	int * rowlist = new int[list_i.size()];
	int * collist = new int[list_j.size()];
	double * vallist = new double[gs.size()];
	for (int k = 0 ; k < gs.size(); k++)
	{
		rowlist[k] = list_i[k];
		collist[k] = list_j[k];
		vallist[k] = (list_rhs[k] - gs[k])<=0 ? list_rhs[k] : list_rhs[k] - gs[k];
	}
	CPXchgcoeflist(cpxenv, cpxlp, gs.size(), rowlist, collist, vallist);
	delete [] rowlist;
	delete []collist;
	delete []vallist;

	for (map<int, OptModel*>::iterator itr = LBmodels.begin(); itr != LBmodels.end(); ++itr)
	{
		OptModel* p = itr->second;
		delete p;
	}

	delete min_model;
}

void WolseyStrengthening::removeRedundantConsts(vector<string>& constNames)
{
	for (int i =0; i < constNames.size(); i++)
	{
		string rowName = constNames[i];
		model_p->removeConstraintFromCPLEX(rowName);
	}
}

void WolseyStrengthening::removeObsoleteVars(vector<string>& varNames)
{
	for (int i =0; i < varNames.size(); i++)
	{
		string varName = varNames[i];
		model_p->removeVariableFromCPLEX(varName);
	}
	model_p->mapVarIndexName();
}


double WolseyStrengthening::calculateBetaByFormula(int& nzcnt,  int* rmatind, double *rmatval , double & lp_optimal)
{
	map<double, int> sorted_a_j;
	for (int j = 0 ; j < nzcnt; j++)
	{
		if (rmatind[j] < DataModel::n)
			sorted_a_j.insert(make_pair(rmatval[j], 0));
	}
	return  lp_optimal*sorted_a_j.begin()->first;
}	

double WolseyStrengthening::calculateBetaBySolvingMIP(int& nzcnt,  int* rmatind, double *rmatval , double & lp_optimal, OptModel::MODEL_TYPE modelType)
{
	int status;
	OptModel subProblem("subproblem");
	subProblem.cpxlp = CPXcloneprob( subProblem.getCPXEnv(), model_p->getCPXLP(), &status);

	//add objective cut
	LHS lhs = subProblem.getObjective();
	RHS rhs = lp_optimal;
	CUT cut(lhs,rhs);
	string cutName = "objcut";
	subProblem.addAConstToCPXModel(cut, cutName, 0);
	subProblem.writeModel2File("subProblem.lp");


	double* values = new double[DataModel::n]();
	int* indices = new int[DataModel::n]();
	for (int j = 0; j < DataModel::n; j++)
	{
		indices[j] = j;
	}
	for (int j = 0 ; j < nzcnt; j++)
	{
		if (rmatind[j] >= DataModel::n)
			break;
// 		indices[j] = rmatind[j];
// 		values[j] = rmatval[j];

		values[  rmatind[j] ] = rmatval[j];

	}
	CPXchgobj(subProblem.getCPXEnv(), subProblem.getCPXLP(), DataModel::n, indices, values);
	delete []values;
	delete []indices;
	
	if (modelType == OptModel::MODEL_TYPE_MIP)
	{
		//subProblem.mapVarNameIndex();
		subProblem.mapVarIndexName();
		subProblem.changeLP2MIP();
	}
 	subProblem.writeModel2File("subproblem.lp");
	subProblem.solve();

	double beta =0;
	if (modelType == OptModel::MODEL_TYPE_MIP)
		beta = subProblem.MIP_optimal_value;
	else
		beta = subProblem.LP_optimal_value;
	iteration_cnt++;
	return beta;
}


double WolseyStrengthening::calculateBetaBySolvingMIP_rev1(int& nzcnt,  int* rmatind, double *rmatval , double & lp_optimal, OptModel::MODEL_TYPE modelType)
{

	int status;
// 	OptModel subProblem("subproblem");
// 	subProblem.cpxlp = CPXcloneprob( subProblem.getCPXEnv(), model_p->getCPXLP(), &status);
 	min_model->writeModel2File("subProblem.lp");
// change objective function
	double* values = new double[DataModel::n]();
	int* indices = new int[DataModel::n]();
	for (int j = 0; j < DataModel::n; j++)
	{
		indices[j] = j;
	}
	for (int j = 0 ; j < nzcnt; j++)
	{
		if (rmatind[j] >= DataModel::n)
			break;
		// 		indices[j] = rmatind[j];
		// 		values[j] = rmatval[j];

		values[  rmatind[j] ] = rmatval[j];
	}
	CPXchgobj(min_model->getCPXEnv(), min_model->getCPXLP(), DataModel::n, indices, values);
	delete []values;
	delete []indices;


		// calculate the minimal for every scenarios with cost vector set to row(the violated one)
	multimap<double,int> map_mins;
	for (int i = 0; i < DataModel::m; i ++)
	{
		double minVal; 
		minVal = calculateMinForScenarioByLP( i);
		map_mins.insert(make_pair(minVal,i));
	}

	int i=0;double min_k =0;
	for(multimap<double,int>::reverse_iterator ritr = map_mins.rbegin(); ritr!= map_mins.rend(); ++ritr, i++)
	{
		if(i == DataModel::k )
		{
			min_k = ritr->first;
			break;
		}
	}
	return min_k;
}


double WolseyStrengthening::calculateBetaBySolvingMIP_rev2_addingViolated2ConstCuts(int& nzcnt,  int* rmatind, double *rmatval , double & lp_optimal, OptModel::MODEL_TYPE modelType)
{

	int status;
	// 	OptModel subProblem("subproblem");
	// 	subProblem.cpxlp = CPXcloneprob( subProblem.getCPXEnv(), model_p->getCPXLP(), &status);
	min_model->writeModel2File("subProblem.lp");
	// change objective function
	double* values = new double[DataModel::n]();
	int* indices = new int[DataModel::n]();
	for (int j = 0; j < DataModel::n; j++)
	{
		indices[j] = j;
	}
	for (int j = 0 ; j < nzcnt; j++)
	{
		if (rmatind[j] >= DataModel::n)
			break;
		// 		indices[j] = rmatind[j];
		// 		values[j] = rmatval[j];

		values[  rmatind[j] ] = rmatval[j];
	}
	CPXchgobj(min_model->getCPXEnv(), min_model->getCPXLP(), DataModel::n, indices, values);
	delete []values;
	delete []indices;


	// calculate the minimal for every scenarios with cost vector set to row(the violated one)
	multimap<double,int> map_mins;
	for (int i = 0; i < DataModel::m; i ++)
	{
		double minVal; 
		minVal = calculateMinForScenarioByLP_withViolated2ConstCuts( i);
		map_mins.insert(make_pair(minVal,i));
	}

	int i=0;double min_k =0;
	for(multimap<double,int>::reverse_iterator ritr = map_mins.rbegin(); ritr!= map_mins.rend(); ++ritr, i++)
	{
		if(i == DataModel::k )
		{
			min_k = ritr->first;
			break;
		}
	}
	return min_k;
}


double WolseyStrengthening::calculateMinForScenarioByLP_withViolated2ConstCuts( int enforced_scenario_num)
{
	/*
	// change objective function
	double* values = new double[DataModel::n]();
	int* indices = new int[DataModel::n]();
	for (int j = 0 ; j < costVec.size(); j++)
	{
	if (costVec[j].first >= DataModel::n)
	break;
	indices[j] = costVec[j].first;
	values[j] = costVec[j].second;
	}
	CPXchgobj(model_mininization.getCPXEnv(), model_mininization.getCPXLP(), DataModel::n, indices, values);
	delete []values;
	delete []indices;
	*/
	// add scenario constraint
	int rowNum = model_p->map_scenario_rownum[enforced_scenario_num];
	CUT scenario = model_p->readScenarioFromModel(rowNum);// index here is row number
	string cstName = "MIXMIN";
	min_model->addAConstToCPXModel(scenario,  cstName, 0);
	//	model_mininization.writeModel2File("minmization.lp"); 
	double gap = 1000;
	while (true)
	{
		min_model->solve();
		gap = min_model->LP_optimal_value -  min_model->Prev_LP_optimal_value;
		if (gap <= 0.1)
			break;
		TwoConstrainSetCutGenerator cutgen(min_model);
		cutgen.addViolatedFacets();
	}
	min_model->removeConstraintFromCPLEX(cstName);
	return min_model->LP_optimal_value;
}


double WolseyStrengthening::calculateMinForScenarioByLP( int enforced_scenario_num)
{
/*
	// change objective function
	double* values = new double[DataModel::n]();
	int* indices = new int[DataModel::n]();
	for (int j = 0 ; j < costVec.size(); j++)
	{
		if (costVec[j].first >= DataModel::n)
			break;
		indices[j] = costVec[j].first;
		values[j] = costVec[j].second;
	}
	CPXchgobj(model_mininization.getCPXEnv(), model_mininization.getCPXLP(), DataModel::n, indices, values);
	delete []values;
	delete []indices;
*/
	// add scenario constraint
	int rowNum = model_p->map_scenario_rownum[enforced_scenario_num];
	CUT scenario = model_p->readScenarioFromModel(rowNum);// index here is row number
	string cstName = "MIXMIN";
	min_model->addAConstToCPXModel(scenario,  cstName, 0);
//	model_mininization.writeModel2File("minmization.lp");
	min_model->solve();
	min_model->removeConstraintFromCPLEX(cstName);
	return min_model->LP_optimal_value;
}




WolseyStrengthening::~WolseyStrengthening(void)
{
}

int WolseyStrengthening::runTestImproved(string& instanceName, ofstream& log, int method)
{
	cout << endl<<"****** Improved WolseyStrengthening started ******"<<endl;

	string iterLogName  = CAT_NAME+"iters.txt";
	ofstream iter_log(iterLogName.c_str(),ios_base::out|ios_base::app);

	OptModel WSmodel("WSModel", 1);
	WSmodel.populateModelFromFile(instanceName+".sav");
	WSmodel.changeMIP2LP();
	model_p = &WSmodel;
	WSmodel.startTimer();
	int cnt = 0;
	while (true)
	{
		if (WSmodel.getElapsedTime() >= 1000000)
		{
			break;
		}
		WSmodel.writeModel2File("strengthen_iter.lp");
		WSmodel.solve();
		iter_log << WSmodel.LP_optimal_value << ", time elapsed, "<< WSmodel.getElapsedTime()<<endl;
		if (fabs(WSmodel.LP_optimal_value - WSmodel.Prev_LP_optimal_value)/(fabs(WSmodel.LP_optimal_value) + SMALL_DECIMAL) < 0.0001)
			break;
		strengthenModel(1);
		cout << "Strengthening Iteration " << cnt  << "   Time elapsed: "<< WSmodel.getElapsedTime()<<endl;
		cnt++;
	}
	iter_log.close();
	log << " WS_Time, "<< WSmodel.getElapsedTime()<<",";
	log<< WSmodel.LP_optimal_value << ",";

	// start improved strengthening
	WSmodel.startTimer();
	WSmodel.Prev_LP_optimal_value = 0;
	WSmodel.LP_optimal_value = 0;
	cnt = 0;
	while (true)
	{
		if (WSmodel.getElapsedTime() >= 1000000)
		{
			break;
		}
		WSmodel.writeModel2File("strengthen_iter.lp");
		WSmodel.solve();
		iter_log << WSmodel.LP_optimal_value << ", time elapsed, "<< WSmodel.getElapsedTime()<<endl;
		if (fabs(WSmodel.LP_optimal_value - WSmodel.Prev_LP_optimal_value)/(fabs(WSmodel.LP_optimal_value) + SMALL_DECIMAL) < 0.0001)
			break;
		strengthenModel(2);
		cout << "Strengthening Iteration " << cnt  << "   Time elapsed: "<< WSmodel.getElapsedTime()<<endl;
		cnt++;
	}
	iter_log.close();
	log << " ImprovedWS_Time, "<< WSmodel.getElapsedTime()<<",ImprovedLB,";
	log<< WSmodel.LP_optimal_value << ",";

	
	
	WSmodel.changeLP2MIP(0);
	WSmodel.writeModel2File(CAT_NAME + "stren1.lp");
	WSmodel.solve();
	log <<"Root_Val,"<< WSmodel.MIP_Root_Obj_Val<<",MIP_Opt,"<< WSmodel.MIP_optimal_value<<","<<WSmodel.MIP_INCUMBENT_OBJ_VAL<<","<<WSmodel.MIP_GAP<<","<<WSmodel.MIP_solving_time<<","<<WSmodel.MIP_BB_nodes<<",";
	return 0;
}


int WolseyStrengthening::runTestImprovedNew(string& instanceName, ofstream& log, int method)
{
	cout << endl<<"****** Improved WolseyStrengthening started ******"<<endl;

	string iterLogName  = CAT_NAME+"iters.txt";
	ofstream iter_log(iterLogName.c_str(),ios_base::out|ios_base::app);

	OptModel WSmodel("WSModel", 1);
	WSmodel.populateModelFromFile(instanceName+".sav");
	WSmodel.changeMIP2LP();
	model_p = &WSmodel;
	WSmodel.startTimer();
	int cnt = 0;
	while (true)
	{
		if (WSmodel.getElapsedTime() >= 1000000)
		{
			break;
		}
		WSmodel.writeModel2File("strengthen_iter.lp");
		WSmodel.solve();
		iter_log << WSmodel.LP_optimal_value << ", time elapsed, "<< WSmodel.getElapsedTime()<<endl;
		if (fabs(WSmodel.LP_optimal_value - WSmodel.Prev_LP_optimal_value)/(fabs(WSmodel.LP_optimal_value) + SMALL_DECIMAL) < 0.0001)
			break;
		strengthenModel(1);
		cout << "Strengthening Iteration " << cnt  << "   Time elapsed: "<< WSmodel.getElapsedTime()<<endl;
		cnt++;
	}
	iter_log.close();
	log << " WS_Time, "<< WSmodel.getElapsedTime()<<",";
	log<< WSmodel.LP_optimal_value << ",";

	// start improved strengthening
	TwoConstrainSetCutGenerator cutGen(model_p);
	WSmodel.startTimer();
	WSmodel.Prev_LP_optimal_value = 0;
	WSmodel.LP_optimal_value = 0;
	cnt = 0;
	while (true)
	{
		if (WSmodel.getElapsedTime() >= 1000000)
		{
			break;
		}
		WSmodel.solve();
		iter_log << WSmodel.LP_optimal_value << ", time elapsed, "<< WSmodel.getElapsedTime()<<endl;
		if (fabs(WSmodel.LP_optimal_value - WSmodel.Prev_LP_optimal_value)/(fabs(WSmodel.LP_optimal_value) + SMALL_DECIMAL) < 0.0001)
			break;
		cutGen.addViolatedFacet();
		cout << "Strengthening Iteration " << cnt  << "   Time elapsed: "<< WSmodel.getElapsedTime()<<endl;
		cnt++;
	}
	iter_log.close();
	log << " ImprovedWS_Time, "<< WSmodel.getElapsedTime()<<",ImprovedLB,";
	log<< WSmodel.LP_optimal_value << ",";



	WSmodel.changeLP2MIP(0);
	WSmodel.writeModel2File(CAT_NAME + "stren1.lp");
	WSmodel.solve();
	log <<"Root_Val,"<< WSmodel.MIP_Root_Obj_Val<<",MIP_Opt,"<< WSmodel.MIP_optimal_value<<","<<WSmodel.MIP_INCUMBENT_OBJ_VAL<<","<<WSmodel.MIP_GAP<<","<<WSmodel.MIP_solving_time<<","<<WSmodel.MIP_BB_nodes<<",";
	return 0;
}




int WolseyStrengthening::runTest(string& instanceName, ofstream& log, int method)
{
	cout << endl<<"****** WolseyStrengthening started ******"<<endl;

	string iterLogName  = CAT_NAME+"iters.txt";
	ofstream iter_log(iterLogName.c_str(),ios_base::out|ios_base::app);

	OptModel WSmodel("WSModel", 1);
	WSmodel.populateModelFromFile(instanceName+".sav");
	WSmodel.changeMIP2LP();
	model_p = &WSmodel;
	WSmodel.startTimer();
	int cnt = 0;
	while (true)
	{
		if (WSmodel.getElapsedTime() >= 10000)
		{
			break;
		}
		WSmodel.writeModel2File("strengthen_iter.lp");
		WSmodel.solve();
		iter_log << WSmodel.LP_optimal_value << ", time elapsed, "<< WSmodel.getElapsedTime()<<endl;
		if (fabs(WSmodel.LP_optimal_value - WSmodel.Prev_LP_optimal_value)/(fabs(WSmodel.LP_optimal_value) + SMALL_DECIMAL) < 0.0001)
			break;
		strengthenModel(method);
		cout << "Strengthening Iteration " << cnt  << "   Time elapsed: "<< WSmodel.getElapsedTime()<<endl;
		cnt++;
	}
	iter_log.close();
	log << " WS_Time, "<< WSmodel.getElapsedTime()<<",";
	log<< WSmodel.LP_optimal_value << ",";
	WSmodel.changeLP2MIP(0);
	WSmodel.writeModel2File(CAT_NAME + "WS.lp");
	WSmodel.solve();
	log <<"Root_Val,"<< WSmodel.MIP_Root_Obj_Val<<",MIP_Opt,"<< WSmodel.MIP_optimal_value<<","<<WSmodel.MIP_INCUMBENT_OBJ_VAL<<","<<WSmodel.MIP_GAP<<","<<WSmodel.MIP_solving_time<<","<<WSmodel.MIP_BB_nodes<<",";
	return 0;
}

int WolseyStrengthening::runTestWiTheoreticalBound(string& instanceName, ofstream& log, int method)
{
	cout << endl<<"****** Strengthening With Theoretical Bound Started ******"<<endl;
	OptModel WSmodel("WSModel", 1);
	WSmodel.populateModelFromFile(instanceName+".sav");
	WSmodel.changeMIP2LP();
	model_p = &WSmodel;
	WSmodel.startTimer();
	strengthenModelWiTheoreticalBound();
	WSmodel.solve();
	log<<"LP after strengthening, "<< WSmodel.LP_optimal_value << ",";
	WSmodel.changeLP2MIP(0);
//	WSmodel.setCPXParameter(CPX_PARAM_NODELIM, 0);
	WSmodel.writeModel2File(CAT_NAME + "WS.lp");
	WSmodel.solve();
	log <<"Root_Val,"<< WSmodel.MIP_Root_Obj_Val<<",MIP_Opt,"<< WSmodel.MIP_optimal_value<<","<<WSmodel.MIP_INCUMBENT_OBJ_VAL<<","<<WSmodel.MIP_GAP<<","<<WSmodel.MIP_solving_time<<","<<WSmodel.MIP_BB_nodes<<",";
	return 0;
}




int WolseyStrengthening::runTestResourceCut(string& instanceName, ofstream& log, int method)
{
	cout << endl<<"****** WolseyStrengthening started ******"<<endl;

	string iterLogName  = CAT_NAME+"iters.txt";
	ofstream iter_log(iterLogName.c_str(),ios_base::out|ios_base::app);

	OptModel WSmodel("WSModel", 1);
	WSmodel.populateModelFromFile(instanceName+".sav");
	WSmodel.changeMIP2LP();
	model_p = &WSmodel;
	WSmodel.startTimer();
	int cnt = 0;
	while (true)
	{
		if (WSmodel.getElapsedTime() >= 10000)
		{
			break;
		}
		WSmodel.writeModel2File("strengthen_iter.lp");
		WSmodel.solve();
		iter_log << WSmodel.LP_optimal_value << ", time elapsed, "<< WSmodel.getElapsedTime()<<endl;
		if (fabs(WSmodel.LP_optimal_value - WSmodel.Prev_LP_optimal_value)/(fabs(WSmodel.LP_optimal_value) + SMALL_DECIMAL) < 0.0001)
			break;
		strengthenModelWithResourceCut(method);
		cout << "Strengthening Iteration " << cnt  << "   Time elapsed: "<< WSmodel.getElapsedTime()<<endl;
		cnt++;
	}
	iter_log.close();
	log << " WS_Time, "<< WSmodel.getElapsedTime()<<",";
	log<< WSmodel.LP_optimal_value << ",";
	WSmodel.changeLP2MIP(0);
	WSmodel.writeModel2File(CAT_NAME + "WS.lp");
	WSmodel.solve();
	log <<"Root_Val,"<< WSmodel.MIP_Root_Obj_Val<<",MIP_Opt,"<< WSmodel.MIP_optimal_value<<","<<WSmodel.MIP_INCUMBENT_OBJ_VAL<<","<<WSmodel.MIP_GAP<<","<<WSmodel.MIP_solving_time<<","<<WSmodel.MIP_BB_nodes<<",";
	return 0;
}



int WolseyStrengthening::runTestResourceCut2ConstCut(string& instanceName, ofstream& log, int method)
{
	cout << endl<<"****** WolseyStrengthening started ******"<<endl;

	string iterLogName  = CAT_NAME+"iters.txt";
	ofstream iter_log(iterLogName.c_str(),ios_base::out|ios_base::app);

	OptModel WSmodel("WSModel", 1);
	WSmodel.populateModelFromFile(instanceName+".sav");
	WSmodel.changeMIP2LP();
	model_p = &WSmodel;
	WSmodel.startTimer();
	int cnt = 0;
	while (true)
	{
		if (WSmodel.getElapsedTime() >= 10000)
		{
			break;
		}
		WSmodel.writeModel2File("strengthen_iter.lp");
		WSmodel.solve();
		iter_log << WSmodel.LP_optimal_value << ", time elapsed, "<< WSmodel.getElapsedTime()<<endl;
		if (fabs(WSmodel.LP_optimal_value - WSmodel.Prev_LP_optimal_value)/(fabs(WSmodel.LP_optimal_value) + SMALL_DECIMAL) < 0.01)
			break;
		strengthenModelWithResourceCut2SetCut(method);
		cout << "Strengthening Iteration " << cnt  << "   Time elapsed: "<< WSmodel.getElapsedTime()<<endl;
		cnt++;
	}
	iter_log.close();
	log << " WS_Time, "<< WSmodel.getElapsedTime()<<",";
	log<< WSmodel.LP_optimal_value << ",";
	WSmodel.changeLP2MIP(0);
	WSmodel.writeModel2File(CAT_NAME + "Res2SetCut.lp");
	WSmodel.solve();
	log <<"Root_Val,"<< WSmodel.MIP_Root_Obj_Val<<",MIP_Opt,"<< WSmodel.MIP_optimal_value<<","<<WSmodel.MIP_INCUMBENT_OBJ_VAL<<","<<WSmodel.MIP_GAP<<","<<WSmodel.MIP_solving_time<<","<<WSmodel.MIP_BB_nodes<<",";
	return 0;
}


int WolseyStrengthening::runTestResourceCut2ConstCutSeparation(string& instanceName, ofstream& log, int method)
{
	cout << endl<<"****** WolseyStrengthening started ******"<<endl;

	string iterLogName  = CAT_NAME+"iters.txt";
	ofstream iter_log(iterLogName.c_str(),ios_base::out|ios_base::app);

	OptModel WSmodel("WSModel", 1);
	WSmodel.populateModelFromFile(instanceName+".sav");
	WSmodel.changeMIP2LP();
	model_p = &WSmodel;
	WSmodel.startTimer();
	int cnt = 0;
	while (true)
	{
		if (WSmodel.getElapsedTime() >= 10000)
		{
			break;
		}
		WSmodel.writeModel2File("strengthen_iter.lp");
		WSmodel.solve();
		iter_log << WSmodel.LP_optimal_value << ", time elapsed, "<< WSmodel.getElapsedTime()<<endl;
		if (fabs(WSmodel.LP_optimal_value - WSmodel.Prev_LP_optimal_value)/(fabs(WSmodel.LP_optimal_value) + SMALL_DECIMAL) < 0.01)
			break;
		strengthenModelNew(method);
		cout << "Strengthening Iteration " << cnt  << "   Time elapsed: "<< WSmodel.getElapsedTime()<<endl;
		cnt++;
	}
	iter_log.close();
	log << " WS_Time, "<< WSmodel.getElapsedTime()<<",";
	log<< WSmodel.LP_optimal_value << ",";
	WSmodel.changeLP2MIP(0);
	WSmodel.writeModel2File(CAT_NAME + "Res2SetCut.lp");
	WSmodel.solve();
	log <<"Root_Val,"<< WSmodel.MIP_Root_Obj_Val<<",MIP_Opt,"<< WSmodel.MIP_optimal_value<<","<<WSmodel.MIP_INCUMBENT_OBJ_VAL<<","<<WSmodel.MIP_GAP<<","<<WSmodel.MIP_solving_time<<","<<WSmodel.MIP_BB_nodes<<",";
	return 0;
}


int WolseyStrengthening::setTimeLimit(int timLmt)
{
	timeLimit = timLmt;
	return 0;
}


double WolseyStrengthening::calculateBetaBySolvingMIP_2ConstSet(int& nzcnt,  int* rmatind, double *rmatval , double & lp_optimal, OptModel::MODEL_TYPE modelType)
{
	// calculate the minimal for every scenarios with cost vector set to row(the violated one)
	multimap<double,int> map_mins;
	for (int i = 0; i < DataModel::m; i ++)
	{
		double minVal; 
		minVal = calculateScenarioLBBy2Const(nzcnt, rmatind, rmatval, i);
		map_mins.insert(make_pair(minVal,i));
	}

	int i=0;double min_k =0;
	for(multimap<double,int>::reverse_iterator ritr = map_mins.rbegin(); ritr!= map_mins.rend(); ++ritr, i++)
	{
		if(i == DataModel::k )
		{
			min_k = ritr->first;
			break;
		}
	}
	return min_k;
}


double WolseyStrengthening::calculateScenarioLBBy2Const(int& nzcnt,  int* rmatind, double *rmatval , int enforced_scenario_num)
{
	OptModel* min_model_p = NULL;
	if(LBmodels.find(enforced_scenario_num) == LBmodels.end() )
	{
		min_model_p = new OptModel(*model_p);
		LBmodels.insert(make_pair(enforced_scenario_num, min_model_p));
		int rowNum =model_p->map_scenario_rownum[enforced_scenario_num];
		CUT scenario = min_model_p->readScenarioFromModel(rowNum);// index here is row number
		string cstName = "MIXMIN";
		min_model_p->addAConstToCPXModel(scenario,  cstName, 0);
		TwoConstrainSetCutGenerator cutGen(min_model_p);
		// add cuts without lifted coefficients
		for (map<int,int>::iterator itr = model_p->map_scenario_rownum.begin(); itr!= model_p->map_scenario_rownum.end(); ++itr)
		{
			int i = itr->second; 
			int j = enforced_scenario_num;
			if (i == j )
				continue;
			vector<TwoConstrainSetCutGenerator::Facet> cuts;
			double* a=NULL; double* b=NULL; double* r=NULL; 
			cutGen.get2Constrains(i , j , a , b ,  r);
			cuts = cutGen.generateFacetsForOneChoiceNoZ2(a, b, r, i, j );
			cutGen.addFacetsToModel(cuts, OptModel::CUT_AS_FORMULATION);
			delete []a;
			delete []b;
			delete r;
		}
		// add cuts with lifted coeffcients
	}
	else
		min_model_p = LBmodels[enforced_scenario_num];

// 
// 	int status;
// 	OptModel subProblem("subproblem");
// 	subProblem.cpxlp = CPXcloneprob( subProblem.getCPXEnv(), model_p->getCPXLP(), &status);
// 	LHS lhs = subProblem.getObjective();
// 	RHS rhs = lp_optimal;
// 	CUT cut(lhs,rhs);
// 	string cutName = "objcut";
// 	subProblem.addAConstToCPXModel(cut, cutName, 0);
//	subProblem.writeModel2File("subProblem.lp");
	// change objective function
	double* values = new double[DataModel::n]();
	int* indices = new int[DataModel::n]();
	for (int j = 0; j < DataModel::n; j++)
	{
		indices[j] = j;
	}
	for (int j = 0 ; j < nzcnt; j++)
	{
		if (rmatind[j] >= DataModel::n)
			break;
		// 		indices[j] = rmatind[j];
		// 		values[j] = rmatval[j];

		values[  rmatind[j] ] = rmatval[j];
	}
	CPXchgobj( min_model_p->getCPXEnv(), min_model_p->getCPXLP(), DataModel::n, indices, values);
	delete []values;
	delete []indices;



	//	min_model_p->writeModel2File("minmization.lp");
	min_model_p->solve();
	return min_model_p->LP_optimal_value;
}

