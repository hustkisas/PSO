//
//  StartUpSequencing.cpp
//  Power System Optimization
//
//  Created by QQ on 9/3/17.
//  Copyright © 2017 QQ. All rights reserved.
//

#include "StartUpSequencing.hpp"

StartUpSequencing::StartUpSequencing(PowerSystemDataModel& dm):dataModel(dm)
{
    //estimate time horizen
    int total_time= 0;
    for(MAP_GEN::iterator itr=dataModel.generators.begin(); itr!=dataModel.generators.end();itr++)
    {
        total_time += itr->second->crank_time + itr->second->ramp_time;
    }
    
    timeHorizen = total_time*1.3;
    numTimePeriods = total_time/TIME_BUCKET;
};

StartUpSequencing::~StartUpSequencing(){};



bool StartUpSequencing::initilizeModel()
{
    try{
        cplexEnv = IloEnv();
        optModel = IloModel(cplexEnv);
        x = NumVarMatrix(cplexEnv,dataModel.generators.size());
        y = NumVarMatrix(cplexEnv,dataModel.criticalLoads.size());
        varTotalTime = IloNumVar(cplexEnv,0);
        
        int cnt = 0;
        for (MAP_GEN::iterator itr = dataModel.generators.begin(); itr!= dataModel.generators.end(); itr++) {
            x[cnt] = IloNumVarArray(cplexEnv, numTimePeriods, 0, 1, ILOINT);
            for (int i=0; i<numTimePeriods; i++) {
                x[cnt][i].setObject(itr->second);
                string varName = std::string("x[")+ ToString(cnt+1) + std::string("][") + ToString(i+1)+std::string("]");
                const char *cstr = varName.c_str();
                x[cnt][i].setName(cstr);
            }
            cnt++;
        }
        cnt = 0;
        for (MAP_BUS::iterator itr = dataModel.criticalLoads.begin(); itr!= dataModel.criticalLoads.end(); itr++) {
            y[cnt] = IloNumVarArray(cplexEnv, numTimePeriods, 0, 1, ILOINT);
            for (int i=0; i<numTimePeriods; i++) {
                y[cnt][i].setObject(itr->second);
                string varName = std::string("y[")+ ToString(cnt+1) + std::string("][") + ToString(i+1)+std::string("]");
                const char *cstr = varName.c_str();
                y[cnt][i].setName(cstr);

            }
            cnt++;
        }

    }
    catch(IloException& e) {
        cerr  << " ERROR: " << e << endl;
        throw;
    }
    catch(...) {
        cerr  << " ERROR" << endl;
        throw;
    }
  
    return true;
};

bool StartUpSequencing::addCapacityCurveConst()
{
    try{
        for (int t=0; t<=numTimePeriods-1; t++) {
            IloExpr lhs(cplexEnv);
            for (int tt=0; tt<=t ; tt++) {
                for (MAP_GEN::iterator itr=dataModel.generators.begin(); itr!=dataModel.generators.end(); itr++) {// generators
                    double p_tt =0;
                    struct Generator* pGen = itr->second;
					if (TIME_BUCKET*(t - tt + 1) > pGen->getCrankTime() + pGen->getRampTime()) // after ramping
						p_tt = pGen->p_max;
					else if (TIME_BUCKET*(t - tt + 1) > pGen->getCrankTime() && TIME_BUCKET*(t - tt + 1) <= pGen->getCrankTime() + pGen->getRampTime())// ramping
					{
						if (TIME_BUCKET*(t - tt + 1) - pGen->getCrankTime() - pGen->getRampTime() >= 0)// last time bucket
						{
							double t0 = int(pGen->getRampTime()) % TIME_BUCKET;
							p_tt = pGen->ramp_rate*TIME_BUCKET*t0*(1 - t0 / (TIME_BUCKET * 2));
						}
						else if (TIME_BUCKET*(t - tt + 1) - pGen->getCrankTime() < TIME_BUCKET) // first time bucket
						{
							double t0 = TIME_BUCKET*(t - tt + 1) - pGen->getCrankTime();
							p_tt = (t0 / TIME_BUCKET)*(t0 / TIME_BUCKET)*pGen->ramp_rate*TIME_BUCKET / 2;
						}
						else // middle
							p_tt = (TIME_BUCKET*(2 * t - 2 * tt + 1) - 2 * pGen->getCrankTime())*pGen->ramp_rate / 2;
					}
					else// cranking
						p_tt = -pGen->crank_p;
                    lhs += x[itr->first-1][tt]*p_tt;
                    cplexEnv.out()<< "t="<<t <<"\ttt="<<tt<<endl;
                }
                for (MAP_BUS::iterator itr=dataModel.criticalLoads.begin(); itr!=dataModel.criticalLoads.end(); itr++) {// generators
                    double d_tt =0;
                    struct Bus* pBus = itr->second;
                    d_tt = -pBus->critical_load;
                    lhs += y[itr->first-1][tt]*d_tt;
                }
                
            }
            string constName = std::string("Cap[")+ ToString(t)+std::string("]");
            IloRange capConstraint(lhs>=0);
            capConstraint.setName(constName.c_str());
            optModel.add(capConstraint);
        }
        //add cardinality constraints
        for (int i = 0; i<x.getSize(); i++) { // each generator is started in only one time period
            optModel.add(IloSum(x[i])==1);
        }
        for (int i = 0; i<y.getSize(); i++) { // each load is started being picked up in only one time period
            optModel.add(IloSum(y[i])==1);
        }
    }
    catch(IloException& e) {
        cerr  << " ERROR: " << e << endl;
        throw;
    }
    catch(...) {
        cerr  << " ERROR" << endl;
        throw;
    }
    
    return true;

};


bool StartUpSequencing::addObjMinTotalTime()
{
    try{
        IloIntArray coef(cplexEnv,numTimePeriods);
        for (int i=0; i<numTimePeriods; i++) {
            coef[i]=(i+1)*TIME_BUCKET;
        }
        for (int i=0; i<dataModel.generators.size(); i++) { // x<=xx
            optModel.add(IloScalProd(coef, x[i])<= varTotalTime);
        }
        for (int i=0; i<dataModel.criticalLoads.size(); i++) { // y<=xx
            optModel.add(IloScalProd(coef, y[i])<= varTotalTime);
        }
        optModel.add(IloMinimize(cplexEnv,varTotalTime));
    }
    catch(IloException& e) {
        cerr  << " ERROR: " << e << endl;
        throw;
    }
    catch(...) {
        cerr  << " ERROR" << endl;
        throw;
    }
    
    return true;

};

bool StartUpSequencing::buildModel()
{
    try{
        addCapacityCurveConst();
        addObjMinTotalTime();
    }
    catch(IloException& e) {
        cerr  << " ERROR: " << e << endl;
        throw;
    }
    catch(...) {
        cerr  << " ERROR" << endl;
        throw;
    }
   
    return true;
};


bool StartUpSequencing::solveModel()
{
    try{
        mipSolver = IloCplex(cplexEnv);
        mipSolver.extract(optModel);
        mipSolver.solve();
        writeCplexModel();
    }
    catch(IloException& e) {
        cerr  << " ERROR: " << e << endl;
        throw;
    }
    catch(...) {
        cerr  << " ERROR" << endl;
        throw;
    }
    
    return true;
};


bool StartUpSequencing::printSolution(std::string fileName)
{
    try{
		if (fileName.empty())
		{
			fileName = "sequence.txt";
		}
 		ofstream log(fileName.c_str(), ios_base::out | ios_base::trunc);
       // out x
		map<int, IloNumVar> sequences;
		map<int, IloNumVar> gen_sequences;
		map<int, IloNumVar> load_sequences;
        for (int i = 0; i < x.getSize(); i++)
        {
            for(int t=0; t <x[i].getSize();t++)
            {
                if (mipSolver.getValue(x[i][t])>=1) {
					sequences.insert(make_pair(t, x[i][t]));
					gen_sequences.insert(make_pair(t, x[i][t]));
//                    cout << x[i][t].getName()<<"=" << mipSolver.getValue(x[i][t])<<endl;
                }
            }
        }
        // output y
        for (int i = 0; i < y.getSize(); i++)
        {
            for(int t=0; t <y[i].getSize();t++)
            {
                if (mipSolver.getValue(y[i][t])>=1) {
					sequences.insert(make_pair(t, y[i][t]));
					load_sequences.insert(make_pair(t, y[i][t]));
//                    cout << y[i][t].getName()<<"=" << mipSolver.getValue(y[i][t])<<endl;
                }
            }
        }
		
		for (map<int, IloNumVar>::iterator itr = sequences.begin(); itr != sequences.end(); itr++)
		{
			int t = itr->first; 
			string var_name = itr->second.getName();
			void *p = itr->second.getObject();
			if (var_name.find('x')!= string::npos)//if generator
			{
				Generator* pGen = (Generator*)p;
				log << t <<'\t'<< pGen->gen_num<<"\tG"<<endl ;
			}
			else
			{
				Bus* pBus = (Bus*)p;
				log << t << '\t' << pBus->bus_num << "\tL" << endl;
			}
		}
        log.close();        
    }
    catch(IloException& e) {
        cerr  << " ERROR: " << e << endl;
        throw;
    }
    catch(...) {
        cerr  << " ERROR" << endl;
        throw;
    }
    
    return true;
};


bool StartUpSequencing::generateSequence()
{
    try{
        initilizeModel();
        buildModel();
        solveModel();
    }
    catch(IloException& e) {
        cerr  << " ERROR: " << e << endl;
        throw;
    }
    catch(...) {
        cerr  << " ERROR" << endl;
        throw;
    }
    
    return true;
};


bool StartUpSequencing::writeCplexModel(std::string fileName)
{
    try{
        mipSolver.exportModel(fileName.c_str());
    }
    catch(IloException& e) {
        cerr  << " ERROR: " << e << endl;
        throw;
    }
    catch(...) {
        cerr  << " ERROR" << endl;
        throw;
    }
    
    return true;


};




