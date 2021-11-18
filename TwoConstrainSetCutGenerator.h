#ifndef TWOCONSTRAINSETCUTGENERATOR_H_HEADER_INCLUDED_BECEF599
#define TWOCONSTRAINSETCUTGENERATOR_H_HEADER_INCLUDED_BECEF599
#include "OptModel.h"


class TwoConstrainSetCutGenerator
{
public:
	struct Facet{
		map<int, double> pi;
		pair<int,double> p;
		pair<int,double> q;
		double pi_0;
	};
	
	OptModel* model;
	TwoConstrainSetCutGenerator(OptModel* md):model(md){};
	TwoConstrainSetCutGenerator(){};
	~TwoConstrainSetCutGenerator(void);
	int get2Constrains(int ind_1, int ind_2, double* &a, double* &b, double* & r);
	vector<Facet> generateFacetsForOneChoice(double* a, double* b, double *r, int ind_1, int ind_2, int n = DataModel::n);
	int addFacetsToModel(vector<Facet>& cuts, OptModel::CUT_TYPE cutType);
	int genAddFacets(int numChoices);

	int addViolatedFacet(void);
	multimap<double,PAIR_INT_INT> selectFracZs(int num_fracZ);
	pair<double,Facet> selectMostViolatedFacet(vector<Facet>& facets);
	double calcViolationForCut(Facet& cut, vector<PAIR_INDEX_VALUE>& x, vector<PAIR_INDEX_VALUE>& z);
	int addAllPairFacets(OptModel::CUT_TYPE cutType);

	// 2012 July 13, for improved strengthening
	vector<TwoConstrainSetCutGenerator::Facet> generateFacetsForOneChoiceNoZ2(double* a, double* b, double *r , int ind_1, int ind_2, int n = DataModel::n);

	// 2012 July 17, add 2-const cuts after strengthening. new lifting since the coefficients of z and rhs are not equal now
	vector<TwoConstrainSetCutGenerator::Facet> generateFacetsForOneChoiceNewLifting(double* a, double* b, double *r , int ind_1, int ind_2, int n = DataModel::n);
	int runTestAddingAllPairs(string& instanceName, ofstream& log, int method);
	int runTestAddingViolatedCuts(string& instanceName, ofstream& log, int method);
	vector<TwoConstrainSetCutGenerator::Facet> generateFacetsForOneChoiceDownLifting(double* a, double* b, double *r , int ind_1, int ind_2, int n = DataModel::n);

	// 2012 Aug 2, combining this with resource cuts.
	int addViolatedFacets(void);



};

#endif
