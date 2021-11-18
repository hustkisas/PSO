#include "DataModel.h"


int DataModel::n =0;
int DataModel::m =0;
double DataModel::epsilon = 0;
int DataModel::k = 0;
int DataModel::nr = 1;
double DataModel::r = -1;
int DataModel::seed = 1;
int DataModel::NUM_COLS = 0;
int DataModel::res_const_sign = 0;
int DataModel::num_level = 1;
int DataModel::c = 1;

int DataModel::readDataFile(string fileName){

	DataModel::r = getInstanceFloatParameterValue(fileName, "r");

	fstream file;
	file.open(fileName.c_str(), ios::in);

	if (!file){
		cerr << "Can not open parameter file, use the default value" <<  endl;
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
        bool is_first_comment_line = true;
		while (getline(file,strLine))// read every line of matrix
		{
            if (strLine.find("%") && is_first_comment_line)
            {
                is_first_comment_line = false;
                continue;
            }
            else if(strLine.find("%") && !is_first_comment_line)
                break;
            
			istringstream stream_line(strLine);
            int bus_i, type; double Pd,Qd,Gs,Bs,area,Vm, Va, baseKV,zone, Vmax, Vmin;
            stream_line >> bus_i>>type>>Pd>>Qd>>Gs>>Bs>>area>>Vm>>Va>>baseKV>>zone>> Vmax>>Vmin;
            
		}

	}
	file.close();
// 	remove(fileName.c_str());
    return 0;
}

