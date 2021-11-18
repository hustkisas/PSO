#pragma once
#ifndef PORTFOLIOINSTANCE_HEADER_INCLUDED_BECEF599
#define PORTFOLIOINSTANCE_HEADER_INCLUDED_BECEF599
#include "OptModel.h"
#include "randomc.h"

class PortfolioInstanceGenerator
{
public:
	PortfolioInstanceGenerator(void);
	~PortfolioInstanceGenerator(void);
	PortfolioInstanceGenerator(int _s, int _nr, int _n, int _m, int _k, double _r, int sign):seed(_s),nr(_nr),n(_n),m(_m),k(_k),r(_r),res_const_sign(sign){};
	int seed;
	int nr; // number of rows in one scenario
	int n;
	int m;
	int k;
	double r;
	int res_const_sign;

	static double A_MAX_RAND ;
	static double A_MIN_RAND;


	CPXENVptr cpxenv;
	CPXLPptr  cpxlp;
	int initilizeCPXlp(string modelName);
	int initilizeCPXenv();

	bool createObjective(OptModel::MODEL_TYPE modelType);
	bool createRandomObjective(OptModel::MODEL_TYPE modelType);
	void addScenarioConstraints();
	void addScenarioConstraintsRandomRHS();
	void addCardinalityConstraint();
	void buildRandomRHSModel(string modelName, OptModel::MODEL_TYPE modelType);
	void addResourceConstraint();


	void addDominantScenarioConstraintsRandomRHS();
	void buildTighterRandomRHSModel(string modelName, OptModel::MODEL_TYPE modelType);

	void buildRandomLHSModel(string modelName, OptModel::MODEL_TYPE modelType);

	void buildRandomLHSModelWRes(string modelName, OptModel::MODEL_TYPE modelType);

	// Sep 3, 2012. New sampling technique, Jim's suggestion
	void buildRandomLHSModelWResNormalDistribution(string modelName, OptModel::MODEL_TYPE modelType); //Jim's suggestion
	void addScenarioConstraintsNormalDistribution();
	bool createObjectiveNormalDistribution(OptModel::MODEL_TYPE modelType);
	
	// Sep 4, 2012. Strengthen model with theoretial bound. This part record the average of each column
	static vector<double> ro;


};

#endif

