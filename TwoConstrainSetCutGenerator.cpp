#include "TwoConstrainSetCutGenerator.h"
#include "randomc.h"
#include "OptModel.h"
extern string CAT_NAME;

TwoConstrainSetCutGenerator::~TwoConstrainSetCutGenerator(void)
{
}


int TwoConstrainSetCutGenerator::get2Constrains( int ind_1, int ind_2, double* &a, double* & b, double* & r){
	CPXENVptr cpxenv = model->getCPXEnv();
	CPXLPptr  lp = model->getCPXLP();

	int num_cols  = CPXgetnumcols(cpxenv,  lp) ; 

	int nzcnt_p = 0 ;
	int* rmatind = new int[num_cols]();
	double* rmatval = new double[num_cols]();
	int rmatspace = num_cols;
	int surplus_p=0;
	int r_begin = ind_1;
	int r_end = ind_1;

	int* rmatbeg = new int(0);

	CPXgetrows     (cpxenv,  lp, &nzcnt_p, 
		rmatbeg, rmatind, rmatval,
		rmatspace, &surplus_p,
		r_begin, r_end);

	double rhs = 0;
	CPXgetrhs(cpxenv,  lp, &rhs, r_begin, r_end);
	r = new double;
	*r= rhs;

	//put coefficients in a 
	a = new double[DataModel::n+1]();
	memset(a,0,sizeof(a)); 
	for (int i = 0; i< nzcnt_p; i++)
	{
		if (rmatind[i] >= DataModel::n && rmatval[i]>0.000001) // if it is the coefficients for z
		{
			a[DataModel::n] = rmatval[i];
		}
		else if(rmatind[i] < DataModel::n)
			a[rmatind[i]] = rmatval[i];
	}
	memset(rmatind,0,sizeof(rmatind)); 
	memset(rmatval,0,sizeof(rmatval)); 
	//fetch coefficients for b
	r_begin = ind_2;
	r_end = ind_2;
	CPXgetrows     (cpxenv,  lp, &nzcnt_p, 
		rmatbeg, rmatind, rmatval,
		rmatspace, &surplus_p,
		r_begin, r_end);

	//put coefficients in b
	b = new double[DataModel::n+1]();
	memset(b,0,sizeof(b)); 
	for (int i = 0; i< nzcnt_p; i++)
	{
		if (rmatind[i] >= DataModel::n && rmatval[i]>0.000001) // if it is the coefficients for z
		{
			b[DataModel::n] = rmatval[i];
		}
		else if(rmatind[i] < DataModel::n)
			b[rmatind[i]] = rmatval[i];
	}

	delete []rmatind;
	delete []rmatval;
	return 0;

}


vector<TwoConstrainSetCutGenerator::Facet> TwoConstrainSetCutGenerator::generateFacetsForOneChoice(double* a, double* b, double *r , int ind_1, int ind_2, int n)
{
	map<int,int> indSet1, indSet2;//indSet2 is bi> ai 
	for (int i = 0; i < n; i++)
	{
		if (b[i]>a[i]+0.000001)
		{
			indSet2.insert(make_pair(i,i));
		}
		else
			indSet1.insert(make_pair(i,i));
	}
	
	vector<Facet> cuts;
	for (map<int,int>::iterator itr = indSet2.begin(); itr != indSet2.end(); ++itr)
	{
		int ind = itr->first;
		double lambda1 = b[ind];
		double lambda2 = a[ind];

		map<int,int> geqSet, lessSet;
		for (int i = 0; i < n; i++)
		{
			if (a[i]/b[i] >= a[ind]/b[ind])
			{
				geqSet.insert(make_pair(i,i));
			}
			else
				lessSet.insert(make_pair(i,i));
		}


		Facet cut;
		for (map<int,int>::iterator itr1 = geqSet.begin(); itr1!= geqSet.end(); ++itr1)
		{
			int index = itr1->first;
			cut.pi.insert(make_pair(index, a[index]*lambda1));
		}

		for (map<int,int>::iterator itr2 = lessSet.begin(); itr2!= lessSet.end(); ++itr2)
		{
			int index = itr2->first;
			cut.pi.insert(make_pair(index, b[index]*lambda2));
		}
		
		cut.p=pair<int,double>( ind_1 , (*r)*(lambda1-lambda2));
		cut.pi_0 = lambda1*(*r);

		//lift z2
		map<double, int> sortedList;
		for (map<int, double>::iterator itrPi= cut.pi.begin(); itrPi != cut.pi.end(); ++itrPi)
		{
			double pi = itrPi->second;
			int index = itrPi->first;
			sortedList.insert(make_pair(pi*(*r)/a[index], 0));
		}
		sortedList.insert(make_pair(cut.p.second, 0));

		cut.q = pair<int, double>( ind_2, cut.pi_0 - sortedList.begin()->first);

		cuts.push_back(cut);
	}
	
	return cuts;
}

vector<TwoConstrainSetCutGenerator::Facet> TwoConstrainSetCutGenerator::generateFacetsForOneChoiceNewLifting(double* a, double* b, double *r , int ind_1, int ind_2, int n)
{

	map<int,int> indSet1, indSet2;//indSet2 is bi> ai 
	for (int i = 0; i < n; i++)
	{
		if (b[i]>a[i]+0.000001)
		{
			indSet2.insert(make_pair(i,i));
		}
		else
			indSet1.insert(make_pair(i,i));
	}

	vector<Facet> cuts;
	for (map<int,int>::iterator itr = indSet2.begin(); itr != indSet2.end(); ++itr)
	{
		int ind = itr->first;
		double lambda1 = b[ind];
		double lambda2 = a[ind];

		map<int,int> geqSet, lessSet;
		for (int i = 0; i < n; i++)
		{
			if (a[i]/b[i] >= a[ind]/b[ind])
			{
				geqSet.insert(make_pair(i,i));
			}
			else
				lessSet.insert(make_pair(i,i));
		}


		Facet cut;
		for (map<int,int>::iterator itr1 = geqSet.begin(); itr1!= geqSet.end(); ++itr1)
		{
			int index = itr1->first;
			cut.pi.insert(make_pair(index, a[index]*lambda1));
		}

		for (map<int,int>::iterator itr2 = lessSet.begin(); itr2!= lessSet.end(); ++itr2)
		{
			int index = itr2->first;
			cut.pi.insert(make_pair(index, b[index]*lambda2));
		}

		cut.p=pair<int,double>( ind_1 , (*r)*(lambda1-lambda2));
		cut.pi_0 = lambda1*(*r);


		//lift z2 
		int idx_z1 = model->getIndex_Z()[ind_1];
		int idx_z2 = model->getIndex_Z()[ind_2];
		double coefZ1, coefZ2;
		int status;
		status = CPXgetcoef (model->getCPXEnv(), model->getCPXLP(), ind_1, idx_z1, &coefZ1);
		status = CPXgetcoef (model->getCPXEnv(), model->getCPXLP(), ind_2, idx_z2, &coefZ2);

		// case: z2=1, z1=0
		vector<double> vec_a(a, a+ n);
		vector<double> vec_b(b, b+ n);
		double t1 = *r;
		double t2 = *r - coefZ2;
		// make the rhs to be equal
		for (int i =0 ; i < vec_a.size(); i++)
			vec_a[i] = t2 <= SMALL_DECIMAL ? vec_a[i]*1 : vec_a[i]*t2;
		for (int i =0 ; i < vec_b.size(); i++)
			vec_b[i] = vec_b[i]*t1;
		double temp = t2;
		t2 = t1*t2;
		t1 = t2==0 ? t1 : t1*temp;

		map<double, int> sortedList;// collect basic feasible solutions
		for (int i = 0; i < vec_a.size(); i++)
		{
			double x_val = t1/vec_a[i] > t2/vec_b[i] ? t1/vec_a[i] : t2/vec_b[i];
			double obj_val = x_val*cut.pi[i];
			sortedList.insert(make_pair(obj_val, 0));
		}
		if (fabs(*r - coefZ2) >= 0.00001)
		{
			for (int i = 0; i < vec_a.size()-1; i++)
			{
				for (int j = i+1; j < vec_a.size(); j++)
				{
					bool A = vec_a[i] >= vec_b[i];
					bool B = vec_a[j] < vec_b[j];
					if (A == B)
					{ 
						double x_i = (t1*vec_b[j] - t2*vec_a[j])/(vec_a[i]*vec_b[j] - vec_a[j]*vec_b[i]);
						double x_j = (t2*vec_a[i] - t1*vec_b[i])/(vec_a[i]*vec_b[j] - vec_a[j]*vec_b[i]);
						double obj_val = cut.pi[i]*x_i + cut.pi[j]*x_j;
						sortedList.insert(make_pair(obj_val,0));
					}
				}
			}
		}

		// case: z2=1, z1=1;
		vec_a = vector<double>(a, a+ n);
		vec_b = vector<double>(b, b+ n);
		t1 = *r- coefZ1;
		t2 = *r - coefZ2;
		// make the rhs to be equal
		for (int i =0 ; i < vec_a.size(); i++)
			vec_a[i] = t2==0 ? vec_a[i] : vec_a[i]*t2;
		for (int i =0 ; i < vec_b.size(); i++)
			vec_b[i] = t1 == 0 ? vec_b[i] : vec_b[i]*t1;
		temp = t2;
		t2 = t1==0? t2 : t1*t2;
		t1 = temp == 0 ? t1 : t1*temp;

		for (int i = 0; i < vec_a.size(); i++)
		{
			double x_val = t1/vec_a[i] > t2/vec_b[i] ? t1/vec_a[i] : t2/vec_b[i];
			double obj_val = x_val*cut.pi[i] + cut.p.second;
			sortedList.insert(make_pair(obj_val, 0));
		}
		if ((fabs(*r - coefZ1) >= 0.00001) && (fabs(*r - coefZ2) >= 0.00001))
		{
			for (int i = 0; i < vec_a.size()-1; i++)
			{
				for (int j = i+1; j < vec_a.size(); j++)
				{
					bool A = vec_a[i] >= vec_b[i];
					bool B = vec_a[j] < vec_b[j];
					if (A == B)
					{ 
						double x_i = (t1*vec_b[j] - t2*vec_a[j])/(vec_a[i]*vec_b[j] - vec_a[j]*vec_b[i]);
						double x_j = (t2*vec_a[i] - t1*vec_b[i])/(vec_a[i]*vec_b[j] - vec_a[j]*vec_b[i]);
						double obj_val = cut.pi[i]*x_i + cut.pi[j]*x_j + cut.p.second;
						sortedList.insert(make_pair(obj_val,0));
					}
				}
			}
		}

		cut.q = pair<int, double>( ind_2, cut.pi_0 - sortedList.begin()->first);

		cuts.push_back(cut);
	}

	return cuts;
}

/*
vector<TwoConstrainSetCutGenerator::Facet> TwoConstrainSetCutGenerator::generateFacetsForOneChoiceDownLifting(double* a, double* b, double *r , int ind_1, int ind_2, int n)
{

	// test -----
	a[0]= 7; a[1]=6;a[2]=1;
	b[0]=2; b[1] =7;b[2]=9;
	n=3;
	// test

	int idx_z1 = model->getIndex_Z()[ind_1];
	int idx_z2 = model->getIndex_Z()[ind_2];
	double coefZ1, coefZ2;
	int status;
	status = CPXgetcoef (model->getCPXEnv(), model->getCPXLP(), ind_1, idx_z1, &coefZ1);
	status = CPXgetcoef (model->getCPXEnv(), model->getCPXLP(), ind_2, idx_z2, &coefZ2);
	//test
	coefZ1 = 

	vector<double> vec_a(a, a+ n);
	vector<double> vec_b(b, b+ n);
	double multiplier = fabs(*r-coefZ2)<SMALL_DECIMAL ? 1 : *r/((*r) - coefZ2);
	for (int i = 0; i < vec_b.size(); i++) // normalize the rhs of bx>=(r-c2) to be r
		vec_b[i] = multiplier*vec_b[i];

	map<int,int> indSet1, indSet2;//indSet2 is bi> ai 
	for (int i = 0; i < n; i++)
	{
		if (vec_b[i]>vec_a[i]+0.000001)
		{
			indSet2.insert(make_pair(i,i));
		}
		else
			indSet1.insert(make_pair(i,i));
	}

	vector<Facet> cuts;
	for (map<int,int>::iterator itr = indSet2.begin(); itr != indSet2.end(); ++itr)
	{
		int ind = itr->first;
		double lambda1 = vec_b[ind];
		double lambda2 = vec_a[ind];

		map<int,int> geqSet, lessSet;
		for (int i = 0; i < n; i++)
		{
			if (vec_a[i]/vec_b[i] >= vec_a[ind]/vec_b[ind])
			{
				geqSet.insert(make_pair(i,i));
			}
			else
				lessSet.insert(make_pair(i,i));
		}


		Facet cut;
		for (map<int,int>::iterator itr1 = geqSet.begin(); itr1!= geqSet.end(); ++itr1)
		{
			int index = itr1->first;
			cut.pi.insert(make_pair(index, vec_a[index]*lambda1));
		}

		for (map<int,int>::iterator itr2 = lessSet.begin(); itr2!= lessSet.end(); ++itr2)
		{
			int index = itr2->first;
			cut.pi.insert(make_pair(index, vec_b[index]*lambda2));
		}

		cut.p=pair<int,double>( ind_1 , (*r)*(lambda1-lambda2));
		cut.pi_0 = lambda1*(*r);


		//down lift z2 

		// case: z2=0, z1=0
		vec_a = vector<double>(a, a+ n);
		vec_b = vector<double>(b, b+ n);
		double t1 = *r;
		double t2 = *r;
		
		map<double, int> sortedList;// collect basic feasible solutions
		for (int i = 0; i < vec_a.size(); i++)
		{
			double x_val = t1/vec_a[i] > t2/vec_b[i] ? t1/vec_a[i] : t2/vec_b[i];
			double obj_val = x_val*cut.pi[i];
			sortedList.insert(make_pair(obj_val, 0));
		}
		for (int i = 0; i < vec_a.size()-1; i++)
		{
			for (int j = i+1; j < vec_a.size(); j++)
			{
				bool A = vec_a[i] >= vec_b[i];
				bool B = vec_a[j] < vec_b[j];
				if (A == B)
				{ 
					double x_i = (t1*vec_b[j] - t2*vec_a[j])/(vec_a[i]*vec_b[j] - vec_a[j]*vec_b[i]);
					double x_j = (t2*vec_a[i] - t1*vec_b[i])/(vec_a[i]*vec_b[j] - vec_a[j]*vec_b[i]);
					double obj_val = cut.pi[i]*x_i + cut.pi[j]*x_j;
					sortedList.insert(make_pair(obj_val,0));
				}
			}
		}

		// case: z2=0, z1=1;
		vec_a = vector<double>(a, a+ n);
		vec_b = vector<double>(b, b+ n);
		t1 = *r- coefZ1;
		t2 = *r;
		// make the rhs to be equal
		for (int i =0 ; i < vec_a.size(); i++)
			vec_a[i] = t2==0 ? vec_a[i] : vec_a[i]*t2;
		for (int i =0 ; i < vec_b.size(); i++)
			vec_b[i] = t1 == 0 ? vec_b[i] : vec_b[i]*t1;
		double temp = t2;
		t2 = t1==0? t2 : t1*t2;
		t1 = temp == 0 ? t1 : t1*temp;

		for (int i = 0; i < vec_a.size(); i++)
		{
			double x_val = t1/vec_a[i] > t2/vec_b[i] ? t1/vec_a[i] : t2/vec_b[i];
			double obj_val = x_val*cut.pi[i] + cut.p.second;
			sortedList.insert(make_pair(obj_val, 0));
		}
		if ((fabs(*r - coefZ1) >= 0.00001) )
		{
			for (int i = 0; i < vec_a.size()-1; i++)
			{
				for (int j = i+1; j < vec_a.size(); j++)
				{
					bool A = vec_a[i] >= vec_b[i];
					bool B = vec_a[j] < vec_b[j];
					if (A == B)
					{ 
						double x_i = (t1*vec_b[j] - t2*vec_a[j])/(vec_a[i]*vec_b[j] - vec_a[j]*vec_b[i]);
						double x_j = (t2*vec_a[i] - t1*vec_b[i])/(vec_a[i]*vec_b[j] - vec_a[j]*vec_b[i]);
						double obj_val = cut.pi[i]*x_i + cut.pi[j]*x_j + cut.p.second;
						sortedList.insert(make_pair(obj_val,0));
					}
				}
			}
		}
		double zeta = sortedList.begin()->first;
		double q = cut.pi_0 - zeta;
		cut.q = pair<int, double>( ind_2,  -q );
		cut.pi_0 = cut.pi_0 - q;
		cuts.push_back(cut);
	}

	return cuts;
}
*/


vector<TwoConstrainSetCutGenerator::Facet> TwoConstrainSetCutGenerator::generateFacetsForOneChoiceNoZ2(double* a, double* b, double *r , int ind_1, int ind_2, int n)
{
	map<int,int> indSet1, indSet2;//indSet2 is bi> ai 
	for (int i = 0; i < n; i++)
	{
		if (b[i]>a[i]+0.000001)
		{
			indSet2.insert(make_pair(i,i));
		}
		else
			indSet1.insert(make_pair(i,i));
	}

	vector<Facet> cuts;
	for (map<int,int>::iterator itr = indSet2.begin(); itr != indSet2.end(); ++itr)
	{
		int ind = itr->first;
		double lambda1 = b[ind];
		double lambda2 = a[ind];

		map<int,int> geqSet, lessSet;
		for (int i = 0; i < n; i++)
		{
			if (a[i]/b[i] >= a[ind]/b[ind])
			{
				geqSet.insert(make_pair(i,i));
			}
			else
				lessSet.insert(make_pair(i,i));
		}


		Facet cut;
		for (map<int,int>::iterator itr1 = geqSet.begin(); itr1!= geqSet.end(); ++itr1)
		{
			int index = itr1->first;
			cut.pi.insert(make_pair(index, a[index]*lambda1));
		}

		for (map<int,int>::iterator itr2 = lessSet.begin(); itr2!= lessSet.end(); ++itr2)
		{
			int index = itr2->first;
			cut.pi.insert(make_pair(index, b[index]*lambda2));
		}

		cut.p=pair<int,double>( ind_1 , (*r)*(lambda1-lambda2));
		cut.pi_0 = lambda1*(*r);

		//lift z2
// 		map<double, int> sortedList;
// 		for (map<int, double>::iterator itrPi= cut.pi.begin(); itrPi != cut.pi.end(); ++itrPi)
// 		{
// 			double pi = itrPi->second;
// 			int index = itrPi->first;
// 			sortedList.insert(make_pair(pi*(*r)/a[index], 0));
// 		}
// 		sortedList.insert(make_pair(cut.p.second, 0));

		cut.q = pair<int, double>( ind_2, 0 /*cut.pi_0 - sortedList.begin()->first*/);

		cuts.push_back(cut);
	}

	return cuts;
}


int TwoConstrainSetCutGenerator::addFacetsToModel(vector<Facet>& cuts, OptModel::CUT_TYPE cutType){

	/* prepare row data*/

	int* rmatind = (int * )   malloc ((DataModel::n+2) * sizeof(int));
	double* rmatval = (double *) malloc ((DataModel::n+2)  * sizeof(double));

	for (int i = 0; i < cuts.size(); i++)
	{
		Facet cut = cuts[i];
		// fill in the entries for x
		int ind = 0;
		for (map<int,double>::iterator itr = cut.pi.begin(); itr != cut.pi.end(); ++itr)
		{
			if(itr->second == 0) // skip zero elements
				continue;
			rmatind[ind] = itr->first;
			rmatval[ind] = itr->second;
			ind ++;
		}
		
		// fill in the entries for z
		if (cut.p.first<cut.q.first)
		{
			rmatind[ind] = model->getIndex_Z()[cut.p.first];
			rmatval[ind] = cut.p.second;
		}
		else{
			rmatind[ind] = model->getIndex_Z()[cut.q.first];
			rmatval[ind] = cut.q.second;
		}
		ind++;
		if (cut.p.first<cut.q.first)
		{
			rmatind[ind] = model->getIndex_Z()[cut.q.first];
			rmatval[ind] = cut.q.second;
		}
		else{
			rmatind[ind] = model->getIndex_Z()[cut.p.first];
			rmatval[ind] = cut.p.second;
		}

		int    rmatbeg[1];
		double rhs[1];
		char   sense[1];
		rhs[0] = cut.pi_0;
		sense[0] = 'G';
		rmatbeg[0] = 0;

		string str_row_name = "2SET"+ToString(cut.p.first)+"_"+ToString(cut.q.first)+"_"+ToString(i);
		char **row_name = new (char*);
		row_name[0]= new char[8];
		strcpy(row_name[0], str_row_name.c_str());

		int status = 0;
//		status = CPXaddrows (model->getCPXEnv(), model->getCPXLP(), 0, 1, ind+1 , rhs,	sense, rmatbeg, rmatind, rmatval, NULL, row_name);
		if (cutType == OptModel::CUT_AS_USERCUT)
		    status = CPXaddusercuts (model->getCPXEnv(), model->getCPXLP(), 1, ind+1 , rhs,	sense, rmatbeg, rmatind, rmatval, row_name);
// 			status = CPXaddusercuts(cpxenv, mip, 1, ind, rhs, sense, rmatbeg, rmatind, rmatval, row_name);
		else  if(cutType == OptModel::CUT_AS_FORMULATION)
 			status = CPXaddrows (model->getCPXEnv(), model->getCPXLP(), 0, 1, ind+1 , rhs,	sense, rmatbeg, rmatind, rmatval, NULL, row_name);
		else if(cutType == OptModel::CUT_AS_LAZYCONSTRAIN)
			status = CPXaddlazyconstraints(model->getCPXEnv(), model->getCPXLP(), 1 , ind+1 ,  rhs,  sense,  rmatbeg,  rmatind, rmatval, row_name);

	delete []row_name;

		if ( status ) {
			fprintf (stderr, "Could not add new rows.\n");
			return 0;
		}

	}

	free_and_null((char**) &rmatind);
	free_and_null((char**) &rmatval);

	return 1;

}


int TwoConstrainSetCutGenerator::genAddFacets(int numChoices){
	int m = DataModel::m;
	CRandomMersenne randGen(1);
	
	map<string, int> uniqueList;
	int k=0;
	int totalNumber = nCr(DataModel::m, 2);//DataModel::m;
	while( k < numChoices && k < totalNumber)
	{
		int i = randGen.IRandom(0, m-2);
		int j = randGen.IRandom( i+1 ,m-1 );
		string key = ToString(i)+"-"+ToString(j);
		uniqueList.insert(make_pair(key,0));

		if(k < uniqueList.size()){//
			double* a=NULL; double* b=NULL; double* r=NULL; 
			get2Constrains(i ,j , a , b ,  r);
			vector<Facet> cuts1 = generateFacetsForOneChoice(a, b, r, i, j);
			vector<Facet> cuts2 = generateFacetsForOneChoice(b, a, r, j, i);
			vector<Facet> cuts;
			cuts.insert(cuts.end(), cuts1.begin(), cuts1.end());
			cuts.insert(cuts.end(), cuts2.begin(), cuts2.end());
			addFacetsToModel(cuts, OptModel::CUT_AS_USERCUT);
			k++;
		}
	}
	return 1;
	

}

/*  currently the most violated facet is added*/
int TwoConstrainSetCutGenerator::addViolatedFacet(void)
{
	multimap<double,PAIR_INT_INT> fracZsSorted = selectFracZs(10);
	vector<PAIR_INT_INT> zS;
	for (multimap<double,PAIR_INT_INT>::iterator itr = fracZsSorted.begin(); itr != fracZsSorted.end(); ++itr)
	{
		zS.push_back(itr->second);
	}

	vector<Facet> cuts;
	for (int i = 0; i < zS.size()-1; i++)
	{
		for (int j = i+1; j < zS.size(); j++)
		{
			double* a=NULL; double* b=NULL; double* r=NULL; 
			get2Constrains(zS[i].first ,zS[j].first , a , b ,  r);
			vector<Facet> cuts1 = generateFacetsForOneChoiceNewLifting(a, b, r, zS[i].first, zS[j].first );
			vector<Facet> cuts2 = generateFacetsForOneChoiceNewLifting(b, a, r, zS[j].first , zS[i].first);
			cuts.insert(cuts.end(), cuts1.begin(), cuts1.end());
			cuts.insert(cuts.end(), cuts2.begin(), cuts2.end());


		}
	}
	pair<double,Facet> maxVioCut =  selectMostViolatedFacet(cuts);
	if (maxVioCut.first <=0)
	{
		return 1;
	}
	vector<Facet> vecCut;
	vecCut.push_back(maxVioCut.second);
	addFacetsToModel(vecCut, OptModel::CUT_AS_FORMULATION);
	return 0;
}

int TwoConstrainSetCutGenerator::addViolatedFacets(void)
{
	model->solve();
	model->writeModel2File("testtt.lp");
	multimap<double,PAIR_INT_INT> fracZsSorted = selectFracZs(50);
	vector<PAIR_INT_INT> zS;
	for (multimap<double,PAIR_INT_INT>::iterator itr = fracZsSorted.begin(); itr != fracZsSorted.end(); ++itr)
	{
		zS.push_back(itr->second);
	}

	map<double,Facet> mapViolatedFacets;
	for (int i = 0; i < zS.size()-1; i++)
	{
		vector<Facet> temp_cuts;
		for (int j = i+1; j < zS.size(); j++)
		{
			double* a=NULL; double* b=NULL; double* r=NULL; 
			get2Constrains(zS[i].first ,zS[j].first , a , b ,  r);
			vector<Facet> cuts1 = generateFacetsForOneChoiceNewLifting(a, b, r, zS[i].first, zS[j].first );
			vector<Facet> cuts2 = generateFacetsForOneChoiceNewLifting(b, a, r, zS[j].first , zS[i].first);
			temp_cuts.insert(temp_cuts.end(), cuts1.begin(), cuts1.end());
			temp_cuts.insert(temp_cuts.end(), cuts2.begin(), cuts2.end());
		}
		pair<double,Facet> maxVioCut =  selectMostViolatedFacet(temp_cuts);
		if (maxVioCut.first > 0.05 )
		{
			mapViolatedFacets.insert(maxVioCut);
		}
	}

	if (mapViolatedFacets.size() ==0)
	{
		return 0;
	}
	vector<Facet> vecCut;
	for (map<double, Facet>::iterator itr = mapViolatedFacets.begin(); itr != mapViolatedFacets.end(); ++itr)
	{
		vecCut.push_back(itr->second);
	}
	addFacetsToModel(vecCut, OptModel::CUT_AS_FORMULATION);
	return 1;
}



multimap<double,PAIR_INT_INT> TwoConstrainSetCutGenerator::selectFracZs(int num_fracZ)
{
	vector<PAIR_INDEX_VALUE> zS = model->getSolutionValues(2);
	multimap<double,PAIR_INT_INT> fracZsSorted;
	for (int i = 0; i < zS.size(); i++)
	{
		PAIR_INDEX_VALUE p= zS[i];
		if(p.second>=SMALL_DECIMAL && p.second <= 1-SMALL_DECIMAL){
			fracZsSorted.insert(make_pair(fabs(0.5-p.second), PAIR_INT_INT(i,p.first)));
		}
	}
	multimap<double,PAIR_INT_INT>::iterator itr = fracZsSorted.begin();
	for(int i = 0; (i < num_fracZ && itr != fracZsSorted.end() ); i++)
		++itr;
	fracZsSorted.erase(itr, fracZsSorted.end());
	return fracZsSorted;
}

pair<double,TwoConstrainSetCutGenerator::Facet> TwoConstrainSetCutGenerator::selectMostViolatedFacet(vector<TwoConstrainSetCutGenerator::Facet>& facets)
{
	vector<PAIR_INDEX_VALUE> x = model->getSolutionValues(1);
	vector<PAIR_INDEX_VALUE> z = model->getSolutionValues(2);
	multimap<double, TwoConstrainSetCutGenerator::Facet> sortedVioList;
	for (int i = 0; i < facets.size(); i++)
	{
		Facet cut = facets[i];
		double violation = calcViolationForCut(cut, x, z);
		sortedVioList.insert(make_pair(violation, cut));
	}
	return pair<double,Facet>(sortedVioList.rbegin()->first, sortedVioList.rbegin()->second);
}

double TwoConstrainSetCutGenerator::calcViolationForCut(Facet& cut, vector<PAIR_INDEX_VALUE>& x, vector<PAIR_INDEX_VALUE>& z)
{
	double summand = 0;
	for (int i = 0; i < x.size(); i++)
	{
		int x_index = x[i].first;
		map<int, double>::iterator itrCoef = cut.pi.find(x_index);
		if ( itrCoef != cut.pi.end())
		{
			summand += (x[i].second)*(itrCoef->second);
		}
	}

	int z1_index = cut.p.first;
	summand += (z[z1_index].second)*(cut.p.second);
	int z2_index = cut.q.first;
	summand += (z[z2_index].second)*(cut.q.second);
	return cut.pi_0 - summand ;
}


int TwoConstrainSetCutGenerator::addAllPairFacets(OptModel::CUT_TYPE cutType)
{
	vector<int> zS= model->getIndex_Z();
	vector<Facet> cuts;
	for (int i = 0; i < zS.size()-1; i++)
	{
		for (int j = i+1; j < zS.size(); j++)
		{
			double* a=NULL; double* b=NULL; double* r=NULL; 
			get2Constrains(i , j , a , b ,  r);
//  			vector<Facet> cuts1 = generateFacetsForOneChoice(a, b, r, i, j );
//  			vector<Facet> cuts2 = generateFacetsForOneChoice(b, a, r, j, i);
 			vector<Facet> cuts1 = generateFacetsForOneChoiceNewLifting(a, b, r, i, j );
 			vector<Facet> cuts2 = generateFacetsForOneChoiceNewLifting(b, a, r, j, i);
// 			vector<Facet> cuts1 = generateFacetsForOneChoiceDownLifting(a, b, r, i, j );
// 			vector<Facet> cuts2 = generateFacetsForOneChoiceDownLifting(b, a, r, j, i);

			cuts.insert(cuts.end(), cuts1.begin(), cuts1.end());
			cuts.insert(cuts.end(), cuts2.begin(), cuts2.end());
			delete []a;
			delete []b;
			delete r;
		}
	}

	addFacetsToModel(cuts, cutType);
	return 0;
}


int TwoConstrainSetCutGenerator::runTestAddingAllPairs(string& instanceName, ofstream& log, int method)
{
	cout << endl<<"****** Two-Constraint Set started ******"<<endl;

	OptModel TCSmodel("TCSmodel", 1);
	TCSmodel.populateModelFromFile(instanceName+".sav");
	TCSmodel.changeMIP2LP();
	model = &TCSmodel;
	TCSmodel.startTimer();
	int cnt = 0;
	
 	addAllPairFacets(OptModel::CUT_AS_FORMULATION);
	TCSmodel.changeLP2MIP(0);
	TCSmodel.writeModel2File(CAT_NAME + "stren1.lp");
	TCSmodel.solve();
	log <<"Root_Val,"<< TCSmodel.MIP_Root_Obj_Val<<",MIP_Opt,"<< TCSmodel.MIP_optimal_value<<","<<TCSmodel.MIP_INCUMBENT_OBJ_VAL<<","<<TCSmodel.MIP_GAP<<","<<TCSmodel.MIP_solving_time<<","<<TCSmodel.MIP_BB_nodes<<",";
	return 0;
}

int TwoConstrainSetCutGenerator::runTestAddingViolatedCuts(string& instanceName, ofstream& log, int method)
{
	cout << endl<<"****** Two-Constraint Set started ******"<<endl;

	OptModel TCSmodel("TCSmodel", 1);
	TCSmodel.populateModelFromFile(instanceName+".sav");
	TCSmodel.changeMIP2LP();
	model = &TCSmodel;
	TCSmodel.startTimer();
	int cnt = 0;
	while(true)
	{
		TCSmodel.solve();
		if (fabs(TCSmodel.LP_optimal_value - TCSmodel.Prev_LP_optimal_value)/(fabs(TCSmodel.LP_optimal_value) + SMALL_DECIMAL) < 0.0001)
			break;
		addViolatedFacet();
	}

	TCSmodel.changeLP2MIP(0);
	TCSmodel.writeModel2File(CAT_NAME + "stren1.lp");
	TCSmodel.solve();
	log <<"Root_Val,"<< TCSmodel.MIP_Root_Obj_Val<<",MIP_Opt,"<< TCSmodel.MIP_optimal_value<<","<<TCSmodel.MIP_INCUMBENT_OBJ_VAL<<","<<TCSmodel.MIP_GAP<<","<<TCSmodel.MIP_solving_time<<","<<TCSmodel.MIP_BB_nodes<<",";
	return 0;
}

