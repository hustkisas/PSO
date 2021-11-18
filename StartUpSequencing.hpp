//
//  StartUpSequencing.hpp
//  Power System Optimization
//
//  Created by QQ on 9/3/17.
//  Copyright © 2017 QQ. All rights reserved.
//

#ifndef StartUpSequencing_hpp
#define StartUpSequencing_hpp

#include "Global.h"
#include "PowerSystemDataModel.hpp"
#endif /* StartUpSequencing_hpp */


class StartUpSequencing
{
public:
    PowerSystemDataModel& dataModel;
    IloEnv cplexEnv;
    IloModel optModel;
    IloCplex mipSolver;
    NumVarMatrix x; // x[i,t]: whether generator i is started in time t
    NumVarMatrix y; // y[i,t]: whether laod i is picke up at time t
    IloNumVar varTotalTime; // total restoration time
    
    
    int timeHorizen;
    int numTimePeriods;
    
    
    StartUpSequencing(PowerSystemDataModel& dm);
    ~StartUpSequencing();
    
    
    bool initilizeModel();
    bool addCapacityCurveConst();
    bool addObjMinTotalTime();
    bool buildModel();
    bool solveModel();
    bool writeCplexModel(std::string fileName = std::string("cplexModel.lp"));
    bool printSolution(std::string fileName = std::string(""));
    bool generateSequence();
    
};
