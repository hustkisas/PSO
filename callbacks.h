#ifndef CALLBACKS_H_HEADER_INCLUDED_BECEF599
#define CALLBACKS_H_HEADER_INCLUDED_BECEF599


#include "Global.h"
#include "DataModel.h"
#include "OptModel.h"
#include "MixingSetDecomposition.h"
#include "CPXTimer.h"

struct CPXRootNodeSimplexInfo{
	int cntItr;
	double mipBound;
	double solvingTime;
	int numCuts;
	double intOpt;
	CPXRootNodeSimplexInfo( int cnt, double bd, double stime, int numc, double opt): cntItr(cnt), mipBound(bd), solvingTime(stime), numCuts(numc), intOpt(opt){ };
} ;

struct CUTCALLBACKINFO{
	vector<int> idx_X;
	vector<int> idx_Z;
	vector<int> idx_S;
	CUTCALLBACKINFO(vector<int> x, vector<int> z, vector<int> s):idx_X(x),idx_Z(z), idx_S(s){};
};


struct BRANCHCALLBACKINFO{
	int branch_index;
	char lu;
	int bd;
	int num_cut_added;
	int num_fixed_ones;
	int num_fixed_zeros;
	bool flag_fixed_at_1;
	bool flag_CH_added;
	bool flag_MIX_added;
	vector<int> vec_index_ones;
	vector<int> vec_index_zeros;

	bool flag_mix_added_on_path;
 	
	BRANCHCALLBACKINFO(int ind, char c, int bound,  int num_0, int num_1, bool flag_direction )
		:branch_index(ind), lu(c), bd(bound), num_fixed_zeros(num_0), num_fixed_ones(num_1),flag_CH_added(false){
			num_cut_added = 0;
			flag_fixed_at_1 = flag_direction;
			flag_mix_added_on_path = false;
			flag_MIX_added = false;
	}
	BRANCHCALLBACKINFO(){};

}   ;


// cut callback global vars
static bool FLAG_ROOT_IC_CUT_ADDED = false;



static int CPXPUBLIC 
myBranchCallback (CPXCENVptr env,
				  void       *cbdata,
				  int        wherefrom,
				  void       *cbhandle,
				  int        type,
				  int        sos,
				  int        nodecnt,
				  int        bdcnt,
				  const double     *nodeest,
				  const int        *nodebeg,
				  const int        *indices,
				  const char       *lu,
				  const int        *bd,
				  int        *useraction_p)

{
	// 	cout <<endl << "Begin of mybranchcallback" <<endl;
	// 	cout << "node count: " << nodecnt<<endl;

	// to do : 1, branch: if there are k zs are fixed to 1 then, only generate one node // done
	//to do : 2, use sum z <= k -t
	// to do : 3, check if the intersection cut is correct

	if (indices[0] < DataModel::n)
	{
		return 0;
	}

	void *p; 
	CPXgetcallbacknodeinfo(env, cbdata, wherefrom, 0, CPX_CALLBACK_INFO_NODE_USERHANDLE, &p);
	BRANCHCALLBACKINFO * t = (BRANCHCALLBACKINFO*)p;

	int node_seq;
	CPXgetcallbacknodeinfo( env,  cbdata, wherefrom, 0 /*nodeindex*/, CPX_CALLBACK_INFO_NODE_SEQNUM,  (void*)(&node_seq));
	// 	cout << "CPX_CALLBACK_INFO_NODE_SEQNUM: "<< node_seq << endl;
	int parent_num_ones = 0;
	int parent_num_zeros = 0;
	vector<int> vec_ind_ones;
	vector<int> vec_ind_zeros;
	if (node_seq != 0)//if not root node
	{
		parent_num_ones = t->num_fixed_ones;
		parent_num_zeros = t->num_fixed_zeros;
		vec_ind_ones = t->vec_index_ones;
		vec_ind_zeros = t->vec_index_zeros;
	}
	for (int i = 0; i < nodecnt ; i++)
	{
		int seqnum1; 
		bool flag_fixed_direction = false;
		int num_fixed_var_0 = parent_num_zeros;
		int num_fixed_var_1 = parent_num_ones;
		if (lu[i] == 'L' && bd[i] == 1)// if branch up
		{
			num_fixed_var_1 = parent_num_ones +1;
			if (num_fixed_var_1 > DataModel::k) // if already k zs are fixed to one, then don't create the node
			{
// 				break;
			}
			flag_fixed_direction = true;
			vec_ind_ones.push_back(indices[i]);
		}
		else // if branch down
		{
			num_fixed_var_0 = parent_num_zeros +1;
			flag_fixed_direction = false;
			vec_ind_zeros.push_back(indices[i]);
		}
		BRANCHCALLBACKINFO 	*p_branch_info = new BRANCHCALLBACKINFO(indices[i], lu[i], bd[i], num_fixed_var_0, num_fixed_var_1, flag_fixed_direction);
		p_branch_info->vec_index_zeros = vec_ind_zeros;
		p_branch_info->vec_index_ones = vec_ind_ones;
		if (node_seq != 0)
		{
			if(t->flag_CH_added || t->flag_mix_added_on_path)
				p_branch_info->flag_mix_added_on_path = true;
		}
        int status = 0; //CPXbranchcallbackbranchbds (env, cbdata, wherefrom,
			//nodeest[i], 1, indices+i, lu+i, bd+i,
			//p_branch_info, &seqnum1);
// 		cout << "sequence node created: " << seqnum1<<endl;
// 		cout<< "z index " <<indices[i]<<" is fixed to" <<flag_fixed_direction<<endl;
	}


	// 	CPXgetcallbacknodeinfo( env,  cbdata, wherefrom, 0 /*nodeindex*/, CPX_CALLBACK_INFO_NODE_VAR,  (void*)result_p);
	// 	cout << "CPX_CALLBACK_INFO_NODE_VAR: "<< *result_p<<endl;
	// 	int branch_index = *result_p;
	// 	CPXgetcallbacknodeinfo( env,  cbdata, wherefrom, 0 /*nodeindex*/, CPX_CALLBACK_INFO_NODE_SEQNUM,  (void*)result_p);
	// 	cout << "CPX_CALLBACK_INFO_NODE_SEQNUM: "<< *result_p << endl;
	// 	delete result_p;
	// 
	// 	CPXLPptr    nodelp_p ;
	// 	CPXgetcallbacknodelp(env, cbdata, wherefrom,  &nodelp_p);
	// 
	// 	double *lb = new double;
	// 	CPXgetlb(env, nodelp_p, lb, branch_index, branch_index);
	// 	cout << "lower bound: " << *lb << endl;
	// 
	// 	double *ub = new double;
	// 	CPXgetub(env, nodelp_p, ub, branch_index, branch_index);
	// 	cout << "upper bound: " << *ub << endl;

	*useraction_p = CPX_CALLBACK_SET;

	// 	cout <<endl << "End of  mybranchcallback" <<endl;
	return 0;
}




static int CPXPUBLIC 
myCutCallbackReverseSign (CPXCENVptr env,
						  void       *cbdata,
						  int        wherefrom,
						  void       *cbhandle,
						  int        *useraction_p)

{
	cout << "inside myCutCallback" <<endl;
	OptModel* model = (OptModel*)cbhandle;

	int* result_p = new int;
	CPXgetcallbacknodeinfo( env,  cbdata, wherefrom, 0 /*nodeindex*/, CPX_CALLBACK_INFO_NODE_SEQNUM,  (void*)result_p);
// 	cout << "Solving Node: "<< *result_p << endl;
	*useraction_p = CPX_CALLBACK_DEFAULT;
	if (*result_p == 0)//root node
	{
		return 0;
	}

	CPXgetcallbacknodeinfo( env,  cbdata, wherefrom, 0 /*nodeindex*/, CPX_CALLBACK_INFO_NODE_DEPTH,  (void*)result_p);
	if (*result_p >= DataModel::k)
	{
// 		return 0;
	}


	CPXgetcallbacknodeinfo( env,  cbdata, wherefrom, 0 /*nodeindex*/, CPX_CALLBACK_INFO_NODE_VAR,  (void*)result_p);
	int branch_index = *result_p;
	
	cout << "CPX_CALLBACK_INFO_NODE_VAR: "<< *result_p<<endl;

	if (*result_p == -1 || branch_index < DataModel::n) // if not branching on variable or not branching on z
	{
		return 0;
	}


	void *p; 
	CPXgetcallbacknodeinfo(env, cbdata, wherefrom, 0, CPX_CALLBACK_INFO_NODE_USERHANDLE, &p);
	BRANCHCALLBACKINFO * t = (BRANCHCALLBACKINFO*)p;
	if (t->lu =! 'L' || t->bd != 1 || t->num_cut_added > 0) // only add reformulation cut to z = 1 node
	{
		return 0;
	}

	int nzcnt; double rhs; int sense; int * cutind; double * cutval;
	cout << endl<< "Cut callback on " << model->getVarNameByIndex(branch_index) << endl;
//	model->prepareCallbackCutCoefsReform(env, branch_index,nzcnt, rhs, sense, cutind, cutval);
	model->prepareCallbackCutCoefsReformNew(env, branch_index,nzcnt, rhs, sense, cutind, cutval);
	CPXcutcallbackaddlocal(env, cbdata, wherefrom, nzcnt, rhs, sense, cutind, cutval);
	/*
	for ( int i = 0; i < nzcnt ; i ++)
	{
		cout << model->getVarNameByIndex( cutind[i]) << "\t"<< cutval[i]<<endl;
	}
	t->num_cut_added ++ ;
*/
/************************************** write down node LP to check correctness **************************/
// 	CPXLPptr nodeLP;
// 	CPXgetcallbacknodelp( env,  cbdata,  wherefrom,  &nodeLP);
// 	string modelName = model->getVarNameByIndex(branch_index) + ".lp";
// 	CPXwriteprob(env, nodeLP, modelName.c_str(), NULL);
	




	// 	CPXLPptr    nodelp_p ;
	// 	CPXgetcallbacknodelp(env, cbdata, wherefrom, &nodelp_p);
	// 	CPXwriteprob(env, nodelp_p, "MIRnodelp.lp", NULL);
	*useraction_p = CPX_CALLBACK_SET;

	delete result_p;
	delete []cutval;
	delete []cutind;

	return 0;

}

static int CPXPUBLIC 
myCutCallbackMixingSet (CPXCENVptr env,
						  void       *cbdata,
						  int        wherefrom,
						  void       *cbhandle,
						  int        *useraction_p)

{
	MixingSetDecomposition* ms = (MixingSetDecomposition*)cbhandle;

	int* result_p = new int;
	CPXgetcallbacknodeinfo( env,  cbdata, wherefrom, 0 /*nodeindex*/, CPX_CALLBACK_INFO_NODE_SEQNUM,  (void*)result_p);
	*useraction_p = CPX_CALLBACK_DEFAULT;
	if (*result_p == 0)//root node
	{
		return 0;
	}

	CPXgetcallbacknodeinfo( env,  cbdata, wherefrom, 0 /*nodeindex*/, CPX_CALLBACK_INFO_NODE_DEPTH,  (void*)result_p);
	if (*result_p >= 5)//DataModel::k)
	{
 		return 0;
	}
//	cout << "inside myCutCallback" <<endl;

	int node_num; 
	CPXgetcallbacknodeinfo(env, cbdata, wherefrom, 0, CPX_CALLBACK_INFO_NODE_SEQNUM, &node_num);
	if(MixingSetDecomposition::Prev_Node_Num == node_num)
		return 0;
	else
		MixingSetDecomposition::Prev_Node_Num = node_num;

	void *p; 
	MixingSetDecomposition* mixingset_p = (MixingSetDecomposition*)cbhandle;

	OptModel* model_p = mixingset_p->model_p;
	int cols = model_p->getNumCols();
	double* nodex= new double[cols];
	CPXgetcallbacknodex(env, cbdata, wherefrom, nodex, 0, cols-1);

	//then run the same process as the root node.
	int nzcnt; double rhs; int sense; int * cutind; double * cutval;
	CPXTimer timer;
	timer.startTimer();
	if(mixingset_p->callbackSeparateCuts(nodex, nzcnt, rhs, sense, cutind, cutval))
	{
		CPXcutcallbackaddlocal(env, cbdata, wherefrom, nzcnt, rhs, sense, cutind, cutval);
//		cout << "Cuts added in branch-and-bound"<<endl;
		delete result_p;
		delete []cutval;
		delete []cutind;
	}
	*useraction_p = CPX_CALLBACK_SET;
	MixingSetDecomposition::B_B_Cut_Gen_Time += timer.getElapsedTime();

	return 0;

}




static int CPXPUBLIC 
myCutCallbackReverseCH (CPXCENVptr env,
						  void       *cbdata,
						  int        wherefrom,
						  void       *cbhandle,
						  int        *useraction_p)

{
	// 	cout << "inside myCutCallback" <<endl;
	OptModel* model = (OptModel*)cbhandle;

	int* result_p = new int;
	CPXgetcallbacknodeinfo( env,  cbdata, wherefrom, 0 /*nodeindex*/, CPX_CALLBACK_INFO_NODE_SEQNUM,  (void*)result_p);
	// 	cout << "CPX_CALLBACK_INFO_NODE_SEQNUM: "<< *result_p << endl;
	*useraction_p = CPX_CALLBACK_DEFAULT;
	if (*result_p == 0)//root node
	{
		return 0;
	}

	CPXgetcallbacknodeinfo( env,  cbdata, wherefrom, 0 /*nodeindex*/, CPX_CALLBACK_INFO_NODE_DEPTH,  (void*)result_p);
	if (*result_p >= 5/*DataModel::k*/)
	{
		return 0;
	}


	CPXgetcallbacknodeinfo( env,  cbdata, wherefrom, 0 /*nodeindex*/, CPX_CALLBACK_INFO_NODE_VAR,  (void*)result_p);
	int branch_index = *result_p;
	// 	cout << "CPX_CALLBACK_INFO_NODE_VAR: "<< *result_p<<endl;

	if (*result_p == -1) // if not branching on variable
	{
		return 0;
	}


	void *p; 
	CPXgetcallbacknodeinfo(env, cbdata, wherefrom, 0, CPX_CALLBACK_INFO_NODE_USERHANDLE, &p);
	BRANCHCALLBACKINFO * t = (BRANCHCALLBACKINFO*)p;
	if (t->lu =! 'L' || t->bd != 1 || t->num_cut_added > 0) // only add reformulation cut to z = 1 node
	{
		return 0;
	}

	if (t->flag_CH_added)
		return 0;
	if(t->num_fixed_ones >=5)
		return 0;



	model->callbackAddReverseCHCuts(cbdata, wherefrom, branch_index);
	t->flag_CH_added = true;
	t->num_cut_added ++ ;

	// 	CPXLPptr    nodelp_p ;
	// 	CPXgetcallbacknodelp(env, cbdata, wherefrom, &nodelp_p);
	// 	CPXwriteprob(env, nodelp_p, "MIRnodelp.lp", NULL);
	*useraction_p = CPX_CALLBACK_SET;

	delete result_p;

	return 0;

}

static int CPXPUBLIC 
myBranchCallbackRecordCPXRootLPOptVal (CPXCENVptr env,
									   void       *cbdata,
									   int        wherefrom,
									   void       *cbhandle,
									   int        type,
									   int        sos,
									   int        nodecnt,
									   int        bdcnt,
									   const double     *nodeest,
									   const int        *nodebeg,
									   const int        *indices,
									   const char       *lu,
									   const int        *bd,
									   int        *useraction_p)
{
	int node_seq;
	CPXgetcallbacknodeinfo( env,  cbdata, wherefrom, 0 /*nodeindex*/, CPX_CALLBACK_INFO_NODE_SEQNUM,  &node_seq);
//	cout << "CPX_CALLBACK_INFO_NODE_SEQNUM: "<< node_seq << endl;

	OptModel* model  = (OptModel*)cbhandle;

	if (node_seq == 0)//root node
	{
		double val =0; 
		CPXgetcallbacknodeobjval(env, cbdata, wherefrom, &val);
		//CPXgetcallbacknodeinfo(env, cbdata, wherefrom, 0, CPX_CALLBACK_INFO_NODE_OBJVAL, &val);
		model->MIP_Root_Obj_Val = val;
	}
	else{// non-root node
	}

	*useraction_p = CPX_CALLBACK_DEFAULT;
	return 0;


}


static int CPXPUBLIC 
myCutCallbackRecordCPLEXRootLP (CPXCENVptr env,
								void       *cbdata,
								int        wherefrom,
								void       *cbhandle,
								int        *useraction_p)

{
	int node_seq;
	CPXgetcallbacknodeinfo( env,  cbdata, wherefrom, 0 /*nodeindex*/, CPX_CALLBACK_INFO_NODE_SEQNUM,  &node_seq);
	cout << "CPX_CALLBACK_INFO_NODE_SEQNUM: "<< node_seq << endl;

	OptModel* model  = (OptModel*)cbhandle;
	string modelName("CPXRoot.lp"); 

	if (node_seq == 0)//root node
	{
		CPXLPptr    nodelp_p ;

		CPXgetcallbacknodelp(env, cbdata, wherefrom, &nodelp_p);
		CPXwriteprob(env, nodelp_p, modelName.c_str(), NULL);
	}
	else{// non-root node
	}

	*useraction_p = CPX_CALLBACK_SET;
	return 0;


}

#endif

