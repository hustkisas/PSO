#include "ResourceCutGenerator.h"
#include "OptModel.h"

ResourceCutGenerator::ResourceCutGenerator(void)
{
}

ResourceCutGenerator::~ResourceCutGenerator(void)
{
}


void ResourceCutGenerator::addAllScenarioCuts()
{
	vector<CUT> cuts;
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
			vector<CUT> t_cuts= genCutForOneScenario(rmatind, rmatval, index_z, rhs);
			cuts.insert(cuts.end(), t_cuts.begin(), t_cuts.end());
			delete rmatbeg;
			delete []rmatind;
			delete []rmatval;
		}
		free_and_null( (char **)&cur_rowname);
		free_and_null((char **)&cur_rownamestore);
	}

	// add the cuts
	string cstName("ResCut");
	for (int i = 0; i < cuts.size(); i++)
	{
		CUT cut = cuts[i];
		model_p->addAConstToCPXModel(cut, cstName);
	}
	
}

vector<CUT>  ResourceCutGenerator::genCutForOneScenario(int* index, double* vals, int idx_z, double rhs )
{
	multimap<double, int> map_val_idx;
	for (int i = 0; i < DataModel::n; i++)
	{
		map_val_idx.insert(make_pair( vals[i], index[i] ));
	}

	vector<CUT> cuts;
	for (multimap<double, int>::iterator itr  = map_val_idx.begin(); itr != map_val_idx.end(); ++itr )
	{
		CUT_COEFICIENTS coefs;
		double subtrahend = itr->first;
		multimap<double, int>::iterator itr_sub = itr;
		itr_sub++;
		for (; itr_sub!= map_val_idx.end(); ++itr_sub)
		{
			double coef = itr_sub->first - subtrahend;
			int idx = itr_sub->second;
			coefs.push_back(make_pair(idx, coef));
		}
		coefs.push_back(make_pair(idx_z, rhs - subtrahend));
		CUT cut(coefs, rhs-subtrahend);
		cuts.push_back(cut);
	}

	return cuts;
}
