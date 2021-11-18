#include "Global.h"
#include <string>
//#include <ilcplex/ilocplex.h>
#include <ilcplex/cplex.h>
#include <stdio.h>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <list>
#include <cmath>
#include "callbacks.h"
#include "randomc.h"
#include "CPXTimer.h"
#include "PowerSystemDataModel.hpp"
#include "StartUpSequencing.hpp"

const int MODEL_TYPE_LP = 0;
const int  MODEL_TYPE_MIP = 1;


typedef   vector<int>   VEC_INT;  
typedef vector<pair<int, double> > CUT_COEFICIENTS;
static int CNT_CUT = 0; 


int BASE_SEED =0;
int RAND_CUT_LMT = 300;
int CUT_SELECTION_METHOD= 1;
int NUM_CUT_PER_ITER = 10;
int NUM_ITERATIONS = 5;
string CAT_NAME = "test";



int main(int argc,char *argv[])
{

	int r =1;
	string instanceName("test_");
    
    PowerSystemDataModel dataModel;
    if(!dataModel.readMatPowerFile(instanceName+"matpower.dat"))
        exit(1);
    if(!dataModel.readRestorationDataFile(instanceName+"restoration.dat"))
        exit(1);
    StartUpSequencing sequencer(dataModel);
    sequencer.generateSequence();
    sequencer.printSolution();
    

    
    
    /*
	if (argc == 2)
	{
		modelFileName = argv[1];
		CAT_NAME = modelFileName ;
	}
	else if (argc == 3)
	{
		modelFileName = argv[1];
		DataModel::epsilon = atof(argv[2]);
		CAT_NAME = modelFileName ;
	}
	else if (argc > 3)
	{

	}*/
    
//	ofstream log(logName.c_str(),ios_base::out|ios_base::app);
//	log<< DataModel::n<<","<<DataModel::m<<","<<DataModel::c <<"," << DataModel::r<<",";

	CPXTimer timer;
	timer.startTimer();

	return 0;



}

