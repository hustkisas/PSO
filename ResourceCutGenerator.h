#ifndef RESOURCECUTGENERATOR_HEADER_INCLUDED_BECEF599
#define RESOURCECUTGENERATOR_HEADER_INCLUDED_BECEF599
#include "Global.h"
#include "CutGenerator.h"

class ResourceCutGenerator:
	public CutGenerator

{
public:
	ResourceCutGenerator(void);
	ResourceCutGenerator(OptModel* p):CutGenerator(p){};
	~ResourceCutGenerator(void);
	void addAllScenarioCuts();
	vector<CUT> genCutForOneScenario(int* index, double* vals, int idx_z, double rhs );
};

#endif
