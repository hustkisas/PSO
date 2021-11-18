#include "MixingSetDecomposition.h"
#include "OptModel.h"
#include "callbacks.h"
#include "WolseyStrengthening.h"


int MixingSetDecomposition::Prev_Node_Num = -1;
double MixingSetDecomposition::B_B_Cut_Gen_Time = 0;

MixingSetDecomposition::MixingSetDecomposition(OptModel* _model)
:model_p(_model), model_mininization("MixingSubproblem"){
	int status;
	model_mininization.cpxlp = CPXcloneprob( model_mininization.getCPXEnv(), model_p->getCPXLP(), &status);
	TIME_CALCULATING_BETA = 0;
	TIME_CUT_GENERATION = 0;

}

MixingSetDecomposition::MixingSetDecomposition()
{
	OptModel t("MixingSubproblem");
	model_mininization = t;
	int status;
	TIME_CALCULATING_BETA = 0;
	TIME_CUT_GENERATION = 0;

};



MixingSetDecomposition::~MixingSetDecomposition(void)
{
}


int MixingSetDecomposition::addSeperateCuts(string& cstType)
{
	static int counter = 0;
	int cols = model_p->getNumCols();
	double* nodex= new double[cols];
	int status = CPXgetx (model_p->getCPXEnv(), model_p->getCPXLP(), nodex, 0, cols-1);

	vector<int> inds_violated = model_p->checkViolatedConstrainsByX(cstType, nodex);

	if (inds_violated.empty() )
	{
		delete []nodex;
		return false;
	}


	int itr_cnt = 0;
	multimap<double,CUT> sortedCuts;
	for (int i = 0; i <  inds_violated.size(); i++)// -1 means cost vector
	{
		if (model_p->getElapsedTime() >10000 || itr_cnt >=3 )
		{
			break;
		}
		int ind_violated =  i == inds_violated.size() ? -1 : inds_violated[i];
		CUT_COEFICIENTS costAlpha = model_p->readAlphaFromConstrain(ind_violated);
//		vector<CUT> cuts = genertaeMixingSetCuts(ind_violated, costAlpha, nodex);  
		vector<CUT> cuts = generateMixingSetCuts(ind_violated, costAlpha, nodex);  
		if (cuts.empty())
		{
			continue;
		}
		else
		{
			for (int i = 0; i < cuts.size(); i++)
			{
				CUT cut = cuts[i];
				double violation = model_p->testCutValidity(cut, nodex) ;
				sortedCuts.insert(make_pair(violation, cut));
			}
			itr_cnt ++;
		}
	}

				// add cuts to master model
	int num_cut =0;
	for (multimap<double,CUT>::iterator itr = sortedCuts.begin(); itr != sortedCuts.end(); ++itr,num_cut++)
	{
		if(num_cut<= 5)// add 10 most violated cuts
		{
			string str = "MIX"+ToString(model_p->num_MIX_rootnode);
			model_p->addAConstToCPXModel(itr->second,str);
			model_p->num_MIX_rootnode++;
		// 				writeModel2File("masterModel.lp");
		}
	}

	delete []nodex;
	return sortedCuts.size()>= 10? 10: sortedCuts.size();	
}


//April 1, 2013, this function is called in cut callback
int MixingSetDecomposition::callbackSeparateCuts(double* &nodex, int &nzcnt, double &rhs, int &sense, int * &cutind, double * &cutval)
{
//	vector<int> inds_violated = model_p->checkViolatedConstrainsByX(string("SN"), nodex); //2013 commented out
	vector<int> inds_violated;
	for(int i = 30; i < 30+DataModel::m; i++)
		inds_violated.push_back(i);

	if (inds_violated.empty() )
	{
		delete []nodex;
		return false;
	}


	multimap<double,CUT> sortedCuts;
	int itr_cnt = 0;
	for (int i = 0; i <  inds_violated.size(); i++)// -1 means cost vector
	{
		if (model_p->getElapsedTime() >10000 || itr_cnt >=3)
		{
			break;
		}
		int ind_violated =  i == inds_violated.size() ? -1 : inds_violated[i];
		CUT_COEFICIENTS costAlpha = model_p->readAlphaFromConstrain(ind_violated);
//		vector<CUT> cuts = genertaeMixingSetCuts(ind_violated, costAlpha, nodex);  
		vector<CUT> cuts = generateMixingSetCuts(ind_violated, costAlpha, nodex, NULL, 2);  
		if (cuts.empty())
		{
			continue;
		}
		else
		{
			itr_cnt ++;
			for (int i = 0; i < cuts.size(); i++)
			{
				CUT cut = cuts[i];
				double violation = model_p->testCutValidity(cut, nodex) ;
				sortedCuts.insert(make_pair(violation, cut));
			}
		}
	}

	if(sortedCuts.size()==0)
		return 0;
	CUT cut = sortedCuts.begin()->second;
	CUT_COEFICIENTS vec_coefs = cut.first;
	nzcnt = vec_coefs.size();
	int status =0;
	/* prepare row data*/

	cutind = (int * )   malloc (nzcnt * sizeof(int));
	cutval = (double *) malloc (nzcnt * sizeof(double));

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

	rhs = cut.second;
	sense = 'G';

	return 1;
}



vector<CUT> MixingSetDecomposition::genertaeMixingSetCuts(int ind_violated, CUT_COEFICIENTS &alpha, double* nodex, BRANCHCALLBACKINFO*nodeInfo)
{
	// calculate the minimal for every scenarios with cost vector set to row(the violated one)
	map< int, multimap<double,int> >::iterator itrFound = mapViolatedScenario_Mins.find(ind_violated);
	multimap<double, int> map_mins;
	if (itrFound == mapViolatedScenario_Mins.end())
	{

		for (int i = 0; i < DataModel::m; i ++)
		{
			double minVal; 
			if (nodeInfo == NULL)
			{
				minVal = calculateMinForScenarioByLP(alpha, i, i == 0 );
// 				minVal = calculateMinForScenarioByFormula(alpha,i);
// 				int a = 0;
			}
			else
				minVal = calculateMinForScenario(alpha, i, nodeInfo);
			map_mins.insert(make_pair(minVal,i));
		}
		if (nodeInfo == NULL)
		{
			mapViolatedScenario_Mins.insert(make_pair( ind_violated, map_mins));
		}
	}
	else
		map_mins = itrFound->second;
	// be careful  of k here, k-i
	// make CMIX cuts 
	CUT_COEFICIENTS vec_mins; int counter =0;
	for (multimap<double, int>::reverse_iterator ritr = map_mins.rbegin(); ritr != map_mins.rend(); ++ritr, counter++)
	{
		if (ritr->first < SMALL_DECIMAL)
		{
			continue;
		}
// 		if ( counter == DataModel::k)
// 		{
// 			double beta_i = WolseyStrengthening::beta[ind_violated] ;
// 			double h_i =  ritr->first ;
// 			if (WolseyStrengthening::beta[ind_violated] > ritr->first )
// 				vec_mins.push_back(PAIR_INDEX_VALUE(ritr->second, WolseyStrengthening::beta[ind_violated]));
// 			else
// 				vec_mins.push_back(PAIR_INDEX_VALUE(ritr->second, ritr->first));
// 		}
// 		else
 			vec_mins.push_back(PAIR_INDEX_VALUE(ritr->second, ritr->first));
	}

	// figure out the "k" here if working on child nodes
	int k;
	if (nodeInfo == NULL)
		k = DataModel::k;
	else
		k = DataModel::k - nodeInfo->vec_index_ones.size();

	// make type (10) inequalities
	vector<CUT> cuts10;
	if( k > vec_mins.size()) // if there are less mins than k
		return cuts10; 

	for (int i = 1 ; i  < k ; i++)
	{
		vector<int> scenarioIndexSet;
		scenarioIndexSet.push_back(0);
		scenarioIndexSet.push_back(i);
		CUT cut = makeType9Cut(scenarioIndexSet, vec_mins, alpha, k); 
		double violation = model_p->testCutValidity(cut, nodex) ;
		if ( violation < -OptModel::MIN_CUT_VIOLATION )
		{
 			cuts10.push_back(cut);
		}

	}


	CUT cut = seperateMVPCut(alpha, vec_mins, nodex, k); //makeType9Cut(scenarioIndexSet, vec_mins, alpha, k); 
	double violation = model_p->testCutValidity(cut, nodex) ;
	if ( violation <  -OptModel::MIN_CUT_VIOLATION)
	{
 	//	cuts10.push_back(cut);
	}
	return cuts10;

}

//Apr 2, 2013. for second submission. 
vector<CUT> MixingSetDecomposition::generateMixingSetCuts(int ind_violated, CUT_COEFICIENTS &alpha, double* nodex, BRANCHCALLBACKINFO*nodeInfo, int cutType)
{
	// calculate the minimal for every scenarios with cost vector set to row(the violated one)
	map< int, multimap<double,int> >::iterator itrFound = mapViolatedScenario_Mins.find(ind_violated);
	multimap<double, int> map_mins;
	if (itrFound == mapViolatedScenario_Mins.end())
	{

		for (int i = 0; i < DataModel::m; i ++)
		{
			double minVal; 
			if (nodeInfo == NULL)
			{
				minVal = calculateMinForScenarioByLP(alpha, i, i == 0 );
// 				minVal = calculateMinForScenarioByFormula(alpha,i);
//				minVal = calculateMinForScenarioWithResConsByFormula(alpha, i);
// 				int a = 0;
			}
			else
				minVal = calculateMinForScenario(alpha, i, nodeInfo);
			map_mins.insert(make_pair(minVal,i));
		}
		if (nodeInfo == NULL)
		{
			mapViolatedScenario_Mins.insert(make_pair( ind_violated, map_mins));
		}
	}
	else
		map_mins = itrFound->second;
	// be careful  of k here, k-i
	// make CMIX cuts 
	CUT_COEFICIENTS vec_mins; int counter =0;
	for (multimap<double, int>::reverse_iterator ritr = map_mins.rbegin(); ritr != map_mins.rend(); ++ritr, counter++)
	{
		if (ritr->first < SMALL_DECIMAL)
		{
			continue;
		}
		vec_mins.push_back(PAIR_INDEX_VALUE(ritr->second, ritr->first));
	}

	// figure out the "k" here if working on child nodes
	int k;
	if (nodeInfo == NULL)
		k = DataModel::k;
	else
		k = DataModel::k - nodeInfo->vec_index_ones.size();

	// make type (10) inequalities
	vector<CUT> cuts;
	if( k > vec_mins.size()) // if there are less mins than k
		return cuts; 

	multimap<double,CUT> sortedType9Cuts;
	if(cutType == 0 || cutType == 1)
	{
		for (int i = 1 ; i  < k ; i++)
		{
			vector<int> scenarioIndexSet;
			scenarioIndexSet.push_back(0);
			scenarioIndexSet.push_back(i);
			CUT cut = makeType9Cut(scenarioIndexSet, vec_mins, alpha, k); 
			double violation = model_p->testCutValidity(cut, nodex) ;
			if ( violation < -OptModel::MIN_CUT_VIOLATION )
			{
				sortedType9Cuts.insert(make_pair(violation, cut));
			}
		}
		int num_cut = 0;
		for(multimap<double,CUT>::iterator itr = sortedType9Cuts.begin(); itr != sortedType9Cuts.end(); ++itr)
		{
			if(num_cut < 2 )
			{
				cuts.push_back(itr->second);
				num_cut ++;
			}
		}
	}

	if(cutType == 0 || cutType == 2)
	{
		CUT cut = seperateMVPCut(alpha, vec_mins, nodex, k); //makeType9Cut(scenarioIndexSet, vec_mins, alpha, k); 
		double violation = model_p->testCutValidity(cut, nodex) ;
		if ( violation <  -OptModel::MIN_CUT_VIOLATION)
		{
 			cuts.push_back(cut);
		}
	}
	return cuts;

}


double MixingSetDecomposition::calculateMinForScenarioWithResConsByFormula(CUT_COEFICIENTS& c, int scenario_index)
{
	int rowNum =model_p->map_scenario_rownum[scenario_index];
	CUT scenario = model_p->readScenarioFromModel(rowNum);// index here is row number

	VEC_MTX_ROW a = scenario.first;
	double rhs = scenario.second;
	vector<int> index_greater;
	vector<int> index_less;
	for (int j = 0; j < DataModel::n; j ++ )
	{
		if (a[j].second<rhs)
			index_less.push_back(j);
		else
			index_greater.push_back(j);
	}

	// make a sorted list for all the cost at extreme points
	map<double, int> sorted_cost;
	// extreme points by sum x = 1 with axles 
	for (int j = 0; j < index_greater.size(); j++)
	{
		double c_j  = c[index_greater[j]].second;
		sorted_cost.insert(make_pair(c_j, 1));
	}

	// extreme points by sum x= 1 and  ax >= 1

	for (int i = 0; i < index_greater.size(); i ++)
	{
		int index_i = index_greater[i];
		for (int j = 0; j < index_less.size(); j++)
		{
			int index_j = index_less[j];
			double cost_ij =  c[index_i].second*( rhs - a[index_j].second)/( a[index_i].second - a[index_j].second ) + c[index_j].second*( a[index_i].second - rhs)/(a[index_i].second-a[index_j].second);
			sorted_cost.insert(make_pair(cost_ij, 1));
		}
	}

	// extreme points by ax>=1 if sum x <=1
	if (DataModel::res_const_sign == 1 ) // if sum x <=1
	{
		for (int i = 0; i < index_greater.size(); i++)
		{
			int index_i = index_greater[i];
			double cost = rhs*c[index_i].second/a[index_i].second;
			sorted_cost.insert(make_pair(cost, 1));
		}
	}

	return sorted_cost.begin()->first;
}




double MixingSetDecomposition::calculateMinForScenarioByFormula(CUT_COEFICIENTS& costVec, int enforced_scenario_num)
{
	int rowNum =model_p->map_scenario_rownum[enforced_scenario_num];
	CUT scenario = model_p->readScenarioFromModel(rowNum);// index here is row number

	VEC_MTX_ROW row = scenario.first;
	map<double, int> map_ratio;
	for (int i = 0; i < costVec.size(); ++i)
	{
		map_ratio.insert(make_pair(costVec[i].second/row[i].second, i));
	}
	double rhs = scenario.second;

	if (map_ratio.empty())
	{
		return -1;
	}
	else
		return rhs*map_ratio.begin()->first;
}



double MixingSetDecomposition::calculateMinForScenario(CUT_COEFICIENTS& costVec, int scenario_number, BRANCHCALLBACKINFO* p_nodeInfo)
{
	// figure out which variables are fixed to 0 or 1
	map<string, int> map_const_name;
	for (int i = 0; i < p_nodeInfo->vec_index_zeros.size(); i++)
	{
		int ind_Z = p_nodeInfo->vec_index_zeros[i];
		int rowNum = 0;
		for (int i = 0 ; i < model_p->idx_Z.size(); i++)
		{
			if (model_p->idx_Z[i]==ind_Z)
			{
				rowNum = i;
				if (rowNum == scenario_number)
				{
					return -1;
				}
				break;
			}
		}
		string scenarioName = "MIR" + ToString(rowNum);
		map_const_name.insert(make_pair(scenarioName, 1));
	}

	for (int i = 0; i < p_nodeInfo->vec_index_ones.size(); i++)
	{
		int ind_Z = p_nodeInfo->vec_index_ones[i];
		int rowNum = 0;
		for (int i = 0 ; i < model_p->idx_Z.size(); i++)
		{
			if (model_p->idx_Z[i]==ind_Z)
			{
				rowNum = i;
				if (rowNum == scenario_number)
				{
					return -1;
				}
				break;
			}
		}
		string scenarioName = "MIR" + ToString(rowNum);
		map_const_name.insert(make_pair(scenarioName, 0));
	}


	// constrain name for scenario_number
	string scenarioName = "MIR" + ToString(scenario_number);
	map_const_name.insert(make_pair(scenarioName, 0));
	//constrain name for resource constrain
	map_const_name.insert(make_pair(string("RES"),2));

	IloEnv SepCutEnv;
	IloModel SepCutModel(SepCutEnv);
	IloNumVarArray x(SepCutEnv, DataModel::n, 0, IloInfinity, ILOFLOAT);
	// create objective 
	IloNumArray c(SepCutEnv);
	for (int i = 0; i < costVec.size(); ++i)
	{
		c.add(costVec[i].second);
	}
	SepCutModel.add(IloMinimize(SepCutEnv, IloScalProd(c,x)));

	// create constrain
	// scenario constrain
	map<string, int>   mapRowNameIndex = model_p->readRowNamesMapFromModel();
	for (map<string, int>::iterator itr = map_const_name.begin(); itr!= map_const_name.end(); ++itr)
	{
		string constName =itr->first;
		int rowIndexForScenario = mapRowNameIndex.find(constName)->second;
		CUT t_a = model_p->readScenarioFromModel(rowIndexForScenario);
		CUT_COEFICIENTS lhs = t_a.first;
		double rhs = t_a.second;
		IloNumArray a(SepCutEnv);
		for (int i = 0; i < lhs.size(); i++)
			a.add(lhs[i].second);
		if (itr->second == 0)
			SepCutModel.add(IloScalProd(a, x) >=rhs ) ;
		else if (itr->second == 1)
			SepCutModel.add(IloScalProd(a, x) <=rhs ) ;
		else
			SepCutModel.add(IloScalProd(a, x) ==rhs ) ;
	}


	IloCplex SepCutSolver(SepCutModel);
	if (SepCutSolver.solve())
	{
		return SepCutSolver.getObjValue();
	}
	else
		return -1;
}


double MixingSetDecomposition::calculateMinForScenarioByLP(CUT_COEFICIENTS& costVec, int enforced_scenario_num, bool flag_new_cost_vec)
{
	IloEnv env;
	IloTimer timer(env);
    timer.start();

	// change objective function
	if (flag_new_cost_vec)
	{
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

	}

	// add scenario constraint
	int rowNum =model_p->map_scenario_rownum[enforced_scenario_num];
	CUT scenario = model_p->readScenarioFromModel(rowNum);// index here is row number
	string cstName = "MIXMIN";
	model_mininization.addAConstToCPXModel(scenario,  cstName, 0);
	model_mininization.writeModel2File("minmization.lp");
	int status = model_mininization.solve();
	if(status != 0)
	{
		cout << "Infeasibility when calculating beta! enforced scenario: " << enforced_scenario_num << endl;
		exit(0);
	}
	model_mininization.removeConstraintFromCPLEX(cstName);
	model_mininization.writeModel2File("mintest.lp");
	TIME_CALCULATING_BETA += timer.getTime();

	return model_mininization.LP_optimal_value;
}



CUT MixingSetDecomposition::makeType9Cut(vector<int>& indexSet, CUT_COEFICIENTS& vec_mins, CUT_COEFICIENTS& vec_alpha, int k)
{
	CUT_COEFICIENTS coefs=vec_alpha;

	for (int i = 0; i < indexSet.size(); i ++)
	{
		int scenario_idx_t = indexSet[i];
		int scenario_idx_t1 = i < indexSet.size()-1 ? indexSet[i+1] : k;
		double coef = vec_mins[scenario_idx_t].second - vec_mins[scenario_idx_t1].second;
		int index_z = model_p->idx_Z[vec_mins[scenario_idx_t].first];
		coefs.push_back(PAIR_INDEX_VALUE(index_z,coef));
	}

	double h1 = vec_mins[indexSet[0]].second;

	return make_pair(coefs, h1);

}

// all index in this function start from 1, not 0 !!!!!
CUT MixingSetDecomposition::seperateMostViolatedMVPCut(CUT_COEFICIENTS& alpha, CUT_COEFICIENTS& h, double* nodex, int k)
{
	const int p = k;
	double *PiA = new double[p+2]; // PiA[0] is not used here, to keep subscripts consistent with the subscripts in separation algorithm
	int *Array_Prev_Index = new int[p+2]; // to store the previous node on the longest path 
	PiA[p+1] = 0;
	Array_Prev_Index[p+1] = -1;

	for (int j = p; j >=1; j--)  
	{
		map<double, int> max_PiAi;
		for (int i = j+1; i <= p+1 ; i++)
		{
			int Zind_j = model_p->idx_Z[h[j-1].first]; 
			double t =  (h[j-1].second -h[i-1].second)*(1- nodex[Zind_j]) + PiA[i]; // use h[j-1] and h[i-1], because h starts from h[0]
			max_PiAi.insert(make_pair(t, i));
		}
		PiA[j] = max_PiAi.rbegin()->first;
		Array_Prev_Index[j] = max_PiAi.rbegin()->second;
	}


	// calculate alpha*x^
	double LHS = h[0].second + PiA[1];
	for (int i = 0; i < alpha.size(); i++)
	{
		int Xind = model_p->idx_X[alpha[i].first];
		double product = alpha[i].second*nodex[Xind];
		LHS = LHS - product;
	}

	double RHS = h[(1) - 1].second - h[(p+1) - 1].second;

	CUT_COEFICIENTS cut_coefs;
	int current_node = 1;
	if(LHS - RHS > 0.00001)// violated cut is identified
	{
		cut_coefs = alpha; // add alpha*x to the cut
		//construct cut by tracking down the longest path.
		while(current_node <= p ){
			int Zind = model_p->idx_Z[ h[current_node - 1].first ];
			int prev_node = Array_Prev_Index[current_node];
			double coef = h[current_node - 1].second - h[prev_node-1].second;
			cut_coefs.push_back(make_pair(Zind, coef));
			current_node = prev_node;
		}
	}

	CUT cut(cut_coefs, h[(1) - 1].second); 

	delete []PiA;
	delete []Array_Prev_Index;
	return cut;
}

CUT MixingSetDecomposition::seperateMVPCut(CUT_COEFICIENTS& alpha, CUT_COEFICIENTS& h, double* nodex, int k)
{
	IloEnv env;
	IloTimer timer(env);
	timer.start();

	// first convert data to MVP format
    double u = h[0].second	- h[k].second;
	// vec_a is the vector for MVP
	vector< pair<double,int> >vec_a;
	vec_a.push_back(make_pair(0,-1));
	for (int i = 0; i < k; i++)
	{
		vec_a.push_back(make_pair(h[k-1-i].second-h[k].second,h[k-1-i].first));
	}

	vector<double> pi(k+1,0);
	vector<int> prevIndex(k+1,0);
	// calculate longest path
	for(int i =1; i< k+1; i++)
	{
		map<double,int> sorted;
		int Zindex = model_p->idx_Z[vec_a[i].second];
		double zVal = 1-nodex[Zindex];
		for (int j=0; j<i; j++)
		{
			double val = pi[j]+(vec_a[i].first -vec_a[j].first)*zVal;
			int prev_ind = j;
			sorted.insert(make_pair(val,prev_ind));
		}
		pi[i] = sorted.rbegin()->first;
		prevIndex[i] = sorted.rbegin()->second;
	}

	// calculate alpha*x^
	double LHS = h[0].second + pi[k];
	for (int i = 0; i < alpha.size(); i++)
	{
		int Xind = model_p->idx_X[alpha[i].first];
		double product = alpha[i].second*nodex[Xind];
		LHS = LHS - product;
	}

	// check violation
	CUT cut;
	CUT_COEFICIENTS cut_coefs;
	double rhs = 0;
	if (LHS - u > SMALL_DECIMAL)
	{
		cut_coefs = alpha; // add alpha*x to the cut
		int idx = k; // the last layer is always picked
		while (idx != 0)
		{
			int index = vec_a[idx].second;
			int z_Idx = model_p->idx_Z[index];
			int prev_idx = prevIndex[idx];
			double val = vec_a[idx].first - (prev_idx==0? 0 : vec_a[prev_idx].first);
			cut_coefs.push_back(make_pair(z_Idx,val));
			if (prev_idx == 0)
			{
				rhs = h[0].second;
			}
			idx = prev_idx;
		}
		cut = CUT(cut_coefs, rhs);
	}

	TIME_CUT_GENERATION += timer.getTime();

	return cut;
}


int MixingSetDecomposition::runTest(string& instanceName, ofstream& log)
{
	cout << endl<<"****** MixingSetDecomposition started ******"<<endl;

	OptModel model("mixSet");
	model.populateModelFromFile(instanceName+".sav");
	int status;
	model_p = &model;
 	model_p->setCPXParameter(CPX_PARAM_PRELINEAR, CPX_OFF);
 	model_p->setCPXParameter(CPX_PARAM_MIPCBREDLP, CPX_OFF);
 	model.setCPXParameter(CPX_PARAM_PREIND, CPX_OFF);

	model_mininization.cpxlp = CPXcloneprob( model_mininization.getCPXEnv(), model.getCPXLP(), &status);
 	model_mininization.setCPXParameter(CPX_PARAM_PRELINEAR, CPX_OFF);
 	model_mininization.setCPXParameter(CPX_PARAM_MIPCBREDLP, CPX_OFF);
 	model_mininization.setCPXParameter(CPX_PARAM_PREIND, CPX_OFF);
	model_mininization.changeMIP2LP();
// 	model.changeMIP2LP();

	model.startTimer();
	for (int i =0; i < 1000 ; i ++)
	{
		if(model.getElapsedTime() > 10000) 
			break;
		model.solve();
		cout << "Strengthening Iteration " << i  << " LP value "<<model.LP_optimal_value << "   Time elapsed: "<< model.getElapsedTime()<<endl;
		cout << model.LP_optimal_value<<endl;
		string cutType("ALL");
		int flag = addSeperateCuts(cutType);
		model.storeCuts("MIX");
		if (flag<=0)
		{
			model.updateCutPool("ALL");
			break;
		}
	}
	log << " time_on_mixing, "<< model.getElapsedTime()<<",";
	log << model.getNumCutsInCutPool("MIX")<< "," << model.LP_optimal_value << ",";
	log << "time on beta," << TIME_CALCULATING_BETA<<",time on separation, " << TIME_CUT_GENERATION<<",";
//	return 0;
 	model.removeInactiveCuts("ALL", 1);
// 	model.addCutPoolToModel("MIX",1);
	model.changeLP2MIP();
	model.setCPXParameter(CPX_PARAM_PRELINEAR, CPX_ON);
	model.setCPXParameter(CPX_PARAM_MIPCBREDLP, CPX_ON);
	model.setCPXParameter(CPX_PARAM_PREIND, CPX_ON);
 	string mixModelName = "MIP.lp";
 	model.writeModel2File(mixModelName);

	// 	model.setBranchCallBackFunction(myBranchCallback, NULL);
	// 	model.setCutCallBackFunction(myCutCallbackReverseSign, &model);
	model.solve();
	log <<"CPX Root,"<<model.MIP_Root_Obj_Val <<",Opt,"<<model.MIP_optimal_value << ",incumb,"<<model.MIP_INCUMBENT_OBJ_VAL<<",gap,"<<model.MIP_GAP<<",time," <<model.MIP_solving_time << ",nodes,"<< model.MIP_BB_nodes<<",";
	return 0;
}


int MixingSetDecomposition::runTest2013(string& instanceName, ofstream& log)
{
	cout << endl<<"****** MixingSetDecomposition started ******"<<endl;

	OptModel model("mixSet");
	model.populateModelFromFile(instanceName+".sav");
	int status;
	model_p = &model;
 	model_p->setCPXParameter(CPX_PARAM_PRELINEAR, CPX_OFF);
 	model_p->setCPXParameter(CPX_PARAM_MIPCBREDLP, CPX_OFF);
 	model.setCPXParameter(CPX_PARAM_PREIND, CPX_OFF);

	model_mininization.cpxlp = CPXcloneprob( model_mininization.getCPXEnv(), model.getCPXLP(), &status);
	string key("SN");
	model_mininization.removeConstraintsFromCPLEXbyKeyWord(key);
	string card("CARD");
	model_mininization.removeConstraintFromCPLEX(card);
	model_mininization.writeModel2File("newtest.lp");

	model_mininization.writeModel2File("new.lp");
 	model_mininization.setCPXParameter(CPX_PARAM_PRELINEAR, CPX_OFF);
 	model_mininization.setCPXParameter(CPX_PARAM_MIPCBREDLP, CPX_OFF);
// 	model_mininization.setCPXParameter(CPX_PARAM_PREIND, CPX_OFF);
	model_mininization.changeMIP2LP();
// 	model.changeMIP2LP();


	model.startTimer();
	for (int i =0; i < 1000 ; i ++)
	{
		if(model.getElapsedTime() > 10000) 
			break;
		model.solve();
		cout << "Strengthening Iteration " << i  << " LP value "<<model.LP_optimal_value << "   Time elapsed: "<< model.getElapsedTime()<<endl;
		cout << model.LP_optimal_value<<endl;
		string cutType("ALL");
		int flag = addSeperateCuts(cutType);
		model.storeCuts("MIX");
		if (flag<=0)
		{
			model.updateCutPool("ALL");
			break;
		}
	}
	log << " time_on_mixing, "<< model.getElapsedTime()<<",";
	log << "num_root_node_cuts_added," << model.getNumCutsInCutPool("MIX")<< ", root_lp" << model.LP_optimal_value << ",";
	log <<"time on separation, " << TIME_CUT_GENERATION<<",";
//	return 0;
 	model.removeInactiveCuts("ALL", 1);
	log << "num_root_node_cuts_left," << model.getNumCutsInCutPool("MIX")<<",";
// 	model.addCutPoolToModel("MIX",1);
	model.changeLP2MIP();
	model.setCutCallBackFunction(myCutCallbackMixingSet,(void *)this);
//	model.setCPXParameter(CPX_PARAM_PRELINEAR, CPX_ON);
//	model.setCPXParameter(CPX_PARAM_MIPCBREDLP, CPX_ON);
	model.setCPXParameter(CPX_PARAM_PREIND, CPX_ON);
 	string mixModelName = "MIP.lp";
 	model.writeModel2File(mixModelName);

	// 	model.setBranchCallBackFunction(myBranchCallback, NULL);
	// 	model.setCutCallBackFunction(myCutCallbackReverseSign, &model);
	model.solve();
	log <<"CPX Root,"<<model.MIP_Root_Obj_Val <<",Opt,"<<model.MIP_optimal_value << ",incumb,"<<model.MIP_INCUMBENT_OBJ_VAL<<",gap,"<<model.MIP_GAP<<",time," <<model.MIP_solving_time <<", Cut_Gen_Time," <<MixingSetDecomposition::B_B_Cut_Gen_Time<< ",nodes,"<< model.MIP_BB_nodes<<",";
	return 0;
}
