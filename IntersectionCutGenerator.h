#ifndef INTERSECTIONCUTGENERATOR_H_HEADER_INCLUDED_BECEF599
#define INTERSECTIONCUTGENERATOR_H_HEADER_INCLUDED_BECEF599
#include "Global.h"
class OptModel;
class IntersectionCutGenerator
{
public:
	IntersectionCutGenerator(OptModel* _optModel);
public:
	~IntersectionCutGenerator(void);

	OptModel* model_p;
	int addOnRoundIC2Root(int numCutPerIteration);
	vector<CUT> generateICFromDiffFracSet(vector<int>& idx_fractionals, int num_cuts);
	list<VEC_INT> partitionFractionalSet(vector<int>& idx_fractionals, int set_size, int num_set);
	vector<pair<int, double> > makeIntersectionCutsFromTableau( vector<int>& idx_disjunctions);


};

#endif
