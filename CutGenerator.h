#ifndef CUTGENERATOR_H_HEADER_INCLUDED_BECEF599
#define CUTGENERATOR_H_HEADER_INCLUDED_BECEF599
#include "Global.h"
class OptModel;

class CutGenerator
{
public:
	CutGenerator(void);
	~CutGenerator(void);

public:
	struct Facet{
		map<int, double> pi;
		pair<int,double> p;
		pair<int,double> q;
		double pi_0;
	};
	OptModel* model_p;
	CutGenerator(OptModel* md):model_p(md){};

};

#endif