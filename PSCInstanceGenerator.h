#pragma once
#ifndef PSCINSTANCEGENERATOR_HEADER_INCLUDED_BECEF599
#define PSCINSTANCEGENERATOR_HEADER_INCLUDED_BECEF599
#include "Global.h"
#include "stocc.h"
#include "OptModel.h"
using namespace std;


class PSCInstanceGenerator
{
public:
	PSCInstanceGenerator(void);
	PSCInstanceGenerator(int num_row, int num_col, int _c/*multi-cover*/, int numLevels, double cr): grid_num_row(num_row), grid_num_col(num_col), c(_c), num_levels(numLevels), coverage(cr)
		{ 
			n = grid_num_col*grid_num_row; 
			m = num_levels*grid_num_col*grid_num_row; 
			grid_edge_length = 20;
			distance_level_2 = 73;
			distance_level_1 = 30; 
	};
	~PSCInstanceGenerator(void);
	int m; // 
	int n; // 
	int c; //multi-cover
	int num_levels;
	double coverage;
	int grid_num_row;
	int grid_num_col;
	double grid_edge_length;
	double distance_level_1;
	double distance_level_2;


	CPXENVptr cpxenv;
	CPXLPptr  cpxlp;

	int initilizeCPXlp(string modelName);
	int initilizeCPXenv();
	bool generateModel(string& modelName, OptModel::MODEL_TYPE modelType, bool flag_random_cost);
	bool createObjective(OptModel::MODEL_TYPE modelType, int num_level, bool flag_random = false);
	bool createSetCoverConstraints(int signal_level, int multi_cover_cnt = 3);
	bool get_a_ij(int level, int idx_signal, int idx_sensor);
	bool  createCardConstraints(int signal_level , double cr);
	bool createRandomSetCoverConstraints(int signal_level, int multi_cover_cnt);


};

#endif