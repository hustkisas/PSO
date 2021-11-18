//
//  PowerSystemDataModel.cpp
//  Power System Optimization
//
//  Created by QQ on 8/30/17.
//  Copyright © 2017 QQ. All rights reserved.
//

#include "PowerSystemDataModel.hpp"

PowerSystemDataModel::PowerSystemDataModel(void){};
PowerSystemDataModel::~PowerSystemDataModel(void)
{
    for(MAP_BUS::iterator itr = buses.begin(); itr!= buses.end(); itr++)
        delete itr->second;
    for(MAP_BRANCH::iterator itr = branches.begin(); itr!= branches.end(); itr++)
        delete itr->second;
    for(MAP_GEN::iterator itr = generators.begin(); itr!= generators.end(); itr++)
        delete itr->second;

};


bool PowerSystemDataModel::readMatPowerFile(std::string fileName)
{
    fstream file;
    file.open(fileName.c_str(), ios::in);
    
    if (!file){
        cerr << "Can not open parameter file, use the default value" <<  endl;
		return false;
    }
    else{
        cout<<"file opened successfully!"<<endl;
        string strLine;
        getline(file,strLine);
        if (strLine.empty())
        {
            cout<<endl<<"Wrong Data"<<endl;
            exit(0);
        }
        
        /************************ read bus info ********************/
        while (getline(file,strLine))// read every line of matrix
        {
            if(strLine.find("%")!=std::string::npos)
                break;
            istringstream stream_line(strLine);
            int bus_i, type; double Pd,Qd,Gs,Bs,area,Vm, Va, baseKV,zone, Vmax, Vmin;
            stream_line >> bus_i>>type>>Pd>>Qd>>Gs>>Bs>>area>>Vm>>Va>>baseKV>>zone>> Vmax>>Vmin;
            struct Bus* pBus = new Bus(bus_i,type,Pd,Qd,Gs,Bs,area,Vm,Va,baseKV,zone,Vmax,Vmin);
            buses.insert(std::make_pair(bus_i,pBus));
        }

        /************************ read generator info ********************/
        int cnt = 1;
        while (getline(file,strLine))// read every line of matrix
        {
            if(strLine.find("%")!=std::string::npos)
                break;
            
            istringstream stream_line(strLine);
            int bus_num; double Pg, Qg,Qmax, Qmin, Vg, mBase, status,Pmax, Pmin, Pc1, Pc2, Qc1min, Qc1max, Qc2min, Qc2max, ramp_agc, ramp_10, ramp_30, ramp_q,apf;
            stream_line >>  bus_num>> Pg>> Qg>>Qmax>> Qmin>> Vg>> mBase>> status>>Pmax>> Pmin>> Pc1>> Pc2>> Qc1min>> Qc1max>> Qc2min>> Qc2max>> ramp_agc>> ramp_10>> ramp_30>> ramp_q>>apf;
            struct Generator *pGen = new struct Generator( cnt, bus_num, Pg, Qg,Qmax, Qmin, Vg, mBase, status,Pmax, Pmin, Pc1, Pc2, Qc1min, Qc1max, Qc2min, Qc2max, ramp_agc, ramp_10, ramp_30, ramp_q,apf);
            generators.insert(std::make_pair(cnt,pGen));
            cnt++;
        }

        /************************ read branch info ********************/
        cnt=1;
        while (getline(file,strLine))// read every line of matrix
        {
            if(strLine.find("%")!=std::string::npos)
                break;
            
            istringstream stream_line(strLine);
            int fbus,tbus; double r, x, b, rateA, rateB, rateC, ratio, angle, status, angmin, angmax;
            
            
            stream_line >> fbus>>tbus>> r>> x>> b>> rateA>> rateB>> rateC>> ratio>> angle>> status>> angmin>> angmax;
            struct Branch *pBranch = new struct Branch( cnt, fbus,tbus, r, x, b, rateA, rateB, rateC, ratio, angle, status, angmin, angmax);
            branches.insert(std::make_pair(cnt,pBranch));
            cnt++;
        }
   
    }
    file.close();
    //     remove(fileName.c_str());
    return true;

};

bool PowerSystemDataModel::readRestorationDataFile(std::string filename)
{
    fstream file;
    file.open(filename.c_str(), ios::in);
    
    if (!file){
        cerr << "Can not open parameter file, use the default value" <<  endl;
		return false;
    }
    else{
        cout<<"file opened successfully!"<<endl;
        string strLine;
        getline(file,strLine);
        if (strLine.empty())
        {
            cout<<endl<<"Wrong Data"<<endl;
            exit(0);
        }
        
        /************************ read generation restart curve data ********************/
        while (getline(file,strLine))// read every line of matrix
        {
            if(strLine.find("%")!=std::string::npos || strLine.find("#")!=std::string::npos)
                break;

            istringstream stream_line(strLine);
            int GenID,UnitBus,UnitType; double CrankPower,CrankTime,RampingRate,RampingTime,Pmax,EarliestTime,LatestTime, Priority;
            stream_line >> GenID>>UnitBus>>UnitType>> CrankPower>>CrankTime>>RampingRate>>RampingTime>>Pmax>>EarliestTime>>LatestTime>>Priority;
            MAP_GEN::iterator itr = generators.find(GenID);
            if (itr == generators.end())
            {
                cout << "Generator "<< GenID << " is not in matpower file"<<endl;
                return false;
            }
            struct Generator *p = itr->second;
            p->is_black= (UnitType>=1);
            p->crank_p = CrankPower;
            p->crank_time = CrankTime;
            p->ramp_rate = RampingRate;
            p->ramp_time = RampingTime;
            p->p_max = Pmax;
            p->earliest_time = EarliestTime;
            p->latest_time = LatestTime;
            p->restoration_priority = Priority;
            if(p->is_black)
                bs_gens.insert(make_pair(GenID, p));
            else
                nbs_gens.insert(make_pair(GenID, p));
        }
        /************************ read critical load data ********************/
        while (getline(file,strLine))// read every line of matrix
        {
            if(strLine.find("%")!=std::string::npos || strLine.find("#")!=std::string::npos)
                break;

            istringstream stream_line(strLine);
            int busID; double critical_load;
            stream_line >> busID>>critical_load;
            MAP_BUS::iterator itr = buses.find(busID);
            if (itr == buses.end())
            {
                cout << "Bus "<< busID << " is not in matpower file"<<endl;
                return false;
            }
            else if(critical_load<=0)
                continue;
            struct Bus *p = itr->second;
            p->critical_load = critical_load;
            criticalLoads.insert(make_pair(criticalLoads.size()+1, p));
        }
        
    }
    file.close();
    //     remove(fileName.c_str());
    return true;

}

bool PowerSystemDataModel::initilizeNetwork()
{
    for(MAP_BRANCH::iterator itr = branches.begin(); itr!=branches.end(); ++itr)
    {
        struct Branch* pBranch = itr->second;
        int from_bus = pBranch->from_bus;
        int to_bus = pBranch->to_bus;
        buses[from_bus]->out_arcs.insert(pBranch->branch_num);
        buses[to_bus]->in_arcs.insert(pBranch->branch_num);
    }
    
    return true;
};
    
    
    
    
    
    
    
    
    
