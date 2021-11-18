#pragma once
#ifndef TRANSPORTATIONINSTANCEGENERATOR_H
#define TRANSPORTATIONINSTANCEGENERATOR_H
#include "Global.h"
#include "OptModel.h"

class TransportationInstanceGenerator
{
public:
	TransportationInstanceGenerator(void);
public:
	~TransportationInstanceGenerator(void);

	double epsilon ;
	int num_supplier;
	int num_demand;
	int num_scenario;

	vector<double> vecScenarioProbability;
	vector<int> vecSupplierCapacity;
	vector< vector<int> > vecDemandScenarios;
	vector< vector<int> > vecShippingCosts;

	void readDataFile(string& fileName);

	CPXENVptr cpxenv;
	CPXLPptr  cpxlp;

	int initilizeCPXlp(string modelName);
	int initilizeCPXenv();


	bool createObjective(OptModel::MODEL_TYPE modelType);
	void buildModel(string& modelName, OptModel::MODEL_TYPE modelType);
	void addScenarioConstraints();
	void addCapacityConstraints();
	void addCardinalityConstraint();
	void addScenarioDominatingConstraints();
	void buildTighterModel(string& modelName, OptModel::MODEL_TYPE modelType);



};

#endif


