#include "IntersectionCutGenerator.h"
#include "OptModel.h"

IntersectionCutGenerator::	IntersectionCutGenerator(OptModel* _optModel)
:model_p(_optModel){
}

IntersectionCutGenerator::~IntersectionCutGenerator(void)
{
}

int IntersectionCutGenerator::addOnRoundIC2Root(int numCutPerIteration)
{
	multimap<double, int> map_fracs = model_p->getAllFractionalZ(); // map_fracs stores the row index! not column index, ranked by closeness to 0.5
	if (map_fracs.size() <DataModel::k+1 - model_p->num_one_at_curr_node)
	{
		cout <<endl<< "num fracs: "<<map_fracs.size()<<"  num fixed 1: "<< model_p->num_one_at_curr_node <<"  :*********   Not enough fractionals  **********"<<endl;
		return -1;
	}

	int cnt_frac = 0;
	vector<int> idx_frac_z ;
	for(multimap<double, int>::iterator itr = map_fracs.begin(); itr != map_fracs.end(); ++itr)
	{
		if (itr->first <= 0.5- SMALL_DECIMAL /*&& cnt_frac <= K+5*/)
			//			if (itr->first <= 1- INT_TOLERANCE  && itr->first >=  INT_TOLERANCE /*&& cnt_frac <= K+5*/)
		{
			idx_frac_z.push_back(itr->second);
			cnt_frac ++;
		}
		else if (cnt_frac <= DataModel::k+1)
		{
			idx_frac_z.push_back(itr->second);
		}
	}

	vector<CUT> cuts;
	cuts = generateICFromDiffFracSet(idx_frac_z,  numCutPerIteration);
	if (cuts.empty())
	{
		return 0;
	}

	for (int i = 0; i < cuts.size(); i ++)
	{
		string str_row_name = "IC"+ToString(model_p->num_IC_rootnode);
		model_p->addAConstToCPXModel(cuts[i], str_row_name, 0);
		model_p->num_IC_rootnode++;
	}
	return cuts.size();
}


vector<CUT> IntersectionCutGenerator::generateICFromDiffFracSet(vector<int>& idx_fractionals, int num_cuts)
{
	list<VEC_INT>list_frac_sets;
	vector<CUT_COEFICIENTS> allCutCoefs;
	list_frac_sets  = partitionFractionalSet(idx_fractionals, DataModel::k+1 - model_p->num_one_at_curr_node, num_cuts);
	int cnt_cut = 0;
	// 	cout << "FractSet: "<<endl;
	for (list<VEC_INT>::iterator itr = list_frac_sets.begin(); itr != list_frac_sets.end(); ++itr)
	{
		if (cnt_cut >= num_cuts)
			break;

		vector<int> vec_frac_set = *itr;
		// 		map<int,int> sortedVec;
		// 		cout<<"------"<<endl;
		// 		for (int p = 0; p < vec_frac_set.size(); p++)
		// 		{
		// 			sortedVec.insert(make_pair(vec_frac_set[p],1));
		// 		}
		// 		for (map<int,int>::iterator itr = sortedVec.begin(); itr != sortedVec.end(); ++itr)
		// 		{
		// 			cout<< itr->first << "\t";
		// 		}
		// 		cout <<endl;

		CUT_COEFICIENTS cutCoefs = makeIntersectionCutsFromTableau(vec_frac_set);
		allCutCoefs.push_back(cutCoefs);
		// 		CUT cut = substituteSlackVarInCut(env, lp, cutCoefs);
		// 		allCuts.push_back(cut);
		cnt_cut++;
	}

	multimap<double, int > rankedListProduct;
	for(int i = 0; i < allCutCoefs.size(); i ++)
	{
		CUT_COEFICIENTS cut = allCutCoefs[i];
		double  squa_sum=0;
		for (int j = 0; j < cut.size(); j ++)
		{
			if (cut[j].second >= SMALL_DECIMAL)
			{
				squa_sum += (cut[j].second)*(cut[j].second);
			}
		}
		rankedListProduct.insert(make_pair(sqrt(squa_sum), i));
	}

	vector<CUT> strongCuts;
	cnt_cut =0;
	map<double, int> mapCoefSum; int prevSize = 0;
	for (multimap<double, int >::iterator  itr  =  rankedListProduct.begin(); itr  != rankedListProduct.end(); ++itr)
	{
		if (cnt_cut >= num_cuts)
			break;
		int index = itr->second;

		prevSize = mapCoefSum.size();
		CUT_COEFICIENTS cut = allCutCoefs[index];
		double coef_key = 0;
		for (int i = 0; i < cut.size(); i ++)
			coef_key += (cut[i].first+1)*cut[i].second;
		mapCoefSum.insert(make_pair(coef_key,0));
		if (mapCoefSum.size() > prevSize)
		{
			strongCuts.push_back(model_p->substituteSlackVarInCut(model_p->getCPXEnv(), model_p->getCPXLP(), cut));
			cnt_cut++;
		}
		else
			cout << endl<<"duplicate cut found  " <<endl;
	}
	return strongCuts;
}


list<VEC_INT> IntersectionCutGenerator::partitionFractionalSet(vector<int>& idx_fractionals, int set_size, int num_set)
{
	int frac_size= idx_fractionals.size();
	int step_size = 0;
	if ( num_set != 1)
	{
		step_size = 	(frac_size-set_size)/(num_set-1) == 0? 1 : (frac_size-set_size)/(num_set-1);
	}
	int num_shift = step_size ==0 ? 0 : (frac_size - set_size)/step_size;

	list<VEC_INT> list_sets;
	for (int i = 0; i <= num_shift; i++)
	{
		VEC_INT frac_set;
		frac_set.insert(frac_set.begin(), idx_fractionals.begin()+ i*step_size, idx_fractionals.begin()+ i*step_size+set_size);
		list_sets.push_back(frac_set);
	}
	return list_sets;
}



vector<pair<int, double> > IntersectionCutGenerator::makeIntersectionCutsFromTableau( vector<int>& idx_disjunctions)
{
	CPXENVptr &env = model_p->getCPXEnv();  
	CPXLPptr & lp=model_p->getCPXLP();
	//model_p->printTableau("tab.csv",env,lp);
	int num_rows = CPXgetnumrows(env,  lp);
	int num_cols  = CPXgetnumcols(env,  lp) ; 

	int *col_basis_status = new int[num_cols];
	int *row_bais_status = new int[num_rows];
	double *val_vars = new double[num_cols];
	CPXgetbase(env,  lp, col_basis_status,  row_bais_status);
	CPXgetx(env, lp, val_vars, 0, num_cols-1);

	// get head of tableau
	int* head = new int[num_rows]; // head is the 0-th column of the tableau
	double *head_val = new double[num_rows];
	int status = CPXgetbhead (env, lp, head, head_val);

	// 	vector<double> f;
	// 	for (int i = 0; i < idx_disjunctions.size(); i++)
	// 	{
	// 		f.push_back(head_val[idx_disjunctions[i]]);
	// 	}

	/* make gauge function */
	vector<pair<int, double> > gauge;
	for( int i =0; i < idx_disjunctions.size(); i++ )
	{
		int row_index = idx_disjunctions[i];
		double coef = -1/head_val[row_index];
		gauge.push_back(pair<int,double>(row_index,coef));
	}

	double* tableauCol;
	vector<pair<int, double> > cutCoefs;
	tableauCol = (double*)malloc(num_rows*sizeof(double));

	// process columns of original variables
	for (int j  =0; j < num_cols; j ++)
	{
		if(col_basis_status[j] == CPX_BASIC) // only look at non-basic row
			continue;
		CPXbinvacol(env, lp, j , tableauCol);

		double max_pi = - CPX_INFBOUND;
		for (int i = 0; i < gauge.size(); i ++)
		{
			pair<int,double> g = gauge[i];
			double pi = g.second*( - tableauCol[g.first]);

			if (pi > max_pi)
				max_pi = pi;
		}
		/* if all pi is negative, that means that, to hit any boundary of max-lattice-free set, 
		it has to go the direction opposite to the ray. so the ray direction is unbounded. */
		if (max_pi < 0) // ignore the ray that goes infinitely. 
		{
			//			max_pi =0;
			cout << endl<<"Note: there is max_pi negative"<<endl;
		}
		cutCoefs.push_back(pair<int,double>(j,max_pi));
	}

	// processing columns of slack variables
	char* sense = new char[num_rows];
	CPXgetsense(env, lp, sense, 0, num_rows -1);
	for (int j  =0; j < num_rows; j ++)
	{
		if (sense[j] == 'E')
			continue;
		if(row_bais_status[j] == CPX_BASIC) // only look at non-basic row
			continue;
		CPXbinvcol(env, lp, j , tableauCol);

		double max_pi = - CPX_INFBOUND;
		for (int i = 0; i < gauge.size(); i ++)
		{
			pair<int,double> g = gauge[i];
			double pi = 0 ;
			if (sense[j] == 'G'){
				pi= g.second*tableauCol[g.first];
			}
			else {
				pi= g.second*( - tableauCol[g.first]);
			}

			if (pi > max_pi)
				max_pi = pi;
		}
		/* if all pi is negative, that means that, to hit any boundary of max-lattice-free set, 
		it has to go the direction opposite to the ray. so the ray direction is unbounded. */
		if (max_pi < 0) // ignore the ray that goes infinitely. 
		{
			// 			max_pi =0;
			cout << endl<<"Note: there is max_pi negative"<<endl;
		}
		cutCoefs.push_back(pair<int,double>(-j-1,max_pi));// slack var uses negative index
	}
	delete []col_basis_status;
	delete []row_bais_status;
	delete []val_vars;
	delete []head;
	delete []head_val;
	delete []sense;
	free_and_null((char **)&tableauCol);

	return cutCoefs;
}

