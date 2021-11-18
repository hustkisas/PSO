#ifndef DATAGENERATOR_H_HEADER_INCLUDED_BECEF599
#define DATAGENERATOR_H_HEADER_INCLUDED_BECEF599
#include "Global.h"
// #define A_MIN_RAND 1//r/10+1
// #define A_MAX_RAND 2*r//1.5*r

static double a_MIN_RAND = 8000;
static double a_MAX_RAN =  15000;
static double b_MIN_RAND = 1.1;
static double b_MAX_RAND =1.1;
static double c_MAX_RAND = 2000;
static double c_MIN_RAND = 1000;


class DataGenerator
{
public:
	DataGenerator(void);
	~DataGenerator(void);

	DataGenerator(int _s, int _nr, int _n, int _m, int _k, double _r, int sign):seed(_s),nr(_nr),n(_n),m(_m),k(_k),r(_r),res_const_sign(sign){};
	int seed;
	int nr; // number of rows in one scenario
	int n;
	int m;
	int k;
	double r;
	int res_const_sign;
	void generateProbDataWithoutResConst(std::string fileName="data.csv");
	void generateProbDataWithoutResConstMultiRows(std::string fileName="data.csv");
	void generateProbDataWithResConstSingleRow(string fileName);
	void generateProbDataSingleRowTest(string fileName);
	void generateProbDataSingleRowWithUnitCostVector(string fileName);
	void generateIntegerDataSingleRowWithUnitCostVector(string fileName);


};

#endif

