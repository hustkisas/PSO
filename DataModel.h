#ifndef DATAMODEL_H_HEADER_INCLUDED_BECEF599
#define DATAMODEL_H_HEADER_INCLUDED_BECEF599

#include "Global.h"

class DataModel
{
public:
	static int n;
	static int m;
	static double epsilon;
	static int k;
	static int nr;
	static int num_level;
	static int c; // used in multi-cover
	static double r; 
	static int seed;
	static int NUM_COLS;
	static int res_const_sign;
	VEC_MTX_ROW costVec; 
	vector<VEC_MTX_ROW>  Arows;
	vector<char> senses;
	vector<double> rhs;

	int readDataFile(string fileName);

};

/****************** PSO **************************/




#endif 

