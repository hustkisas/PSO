#ifndef GLOBAL_H_HEADER_INCLUDED_BECEF599
#define GLOBAL_H_HEADER_INCLUDED_BECEF599

#include <ilconcert/iloenv.h>
#include <ilcplex/ilocplex.h>
#include <stdio.h>
#include <vector>
#include <map>
#include <utility>
#include <list>
#include <set>
#pragma warning(disable: 4996) 

typedef std::pair<int,double> PAIR_INDEX_VALUE;
typedef std::pair<int,int> PAIR_INT_INT;
typedef std::pair<double, int> PAIR_VALUE_INDEX;
typedef std::pair<double, int> PAIR_DOUBLE_INT;
typedef std::vector<PAIR_INDEX_VALUE> VEC_MTX_ROW;
typedef std::vector<double> VEC_DOUBLE;
typedef std::vector<int> VEC_INT;
typedef std::set<int> SET_BUS;
typedef std::set<int> SET_LINE;
typedef std::set<int> SET_INT;
typedef std::map<int,struct Bus*> MAP_BUS;
typedef std::map<int,struct Branch*> MAP_BRANCH;
typedef std::map<int,struct Generator*> MAP_GEN;

typedef std::vector<std::pair<int, double> > CUT_COEFICIENTS;
typedef double CUT_RHS;
typedef std::pair< CUT_COEFICIENTS, CUT_RHS> CUT;

typedef std::vector<std::pair<int, double> > LHS;
typedef double RHS;
typedef std::pair< LHS, RHS> INEQUALITY;

typedef IloArray<IloIntArray> IloIntMatrix; 
typedef IloArray<IloNumArray> IloNumMatrix;
typedef IloArray<IloNumVarArray> NumVarMatrix;
typedef IloArray<IloIntVarArray> IntVarMatrix;

const int TIME_BUCKET = 3;

struct Generator{
    int gen_num;
    int host_bus;
    double Pg;
    double Qg;
    double q_max;
    double q_min;
    double Vg; //Vg, voltage magnitude setpoint (p.u.)
    double mBase; // mBase, total MVA base of this machine, defaults to baseMVA
    int status; // >0 in-srervice
    double p_max;
    double p_min;
    double Pc1; //lower real power output of PQ capability curve (MW)
    double Pc2; //upper real power output of PQ capability curve (MW)
    double Qc1min; //minimum reactive power output at Pc1 (MVAr)
    double Qc1max; //maximum reactive power output at Pc1 (MVAr)
    double Qc2min; //minimum reactive power output at Pc2 (MVAr)
    double Qc2max; //maximum reactive power output at Pc2 (MVAr)
    double ramp_agc; // ramp rate for load following/AGC (MW/min)
    double ramp_10; // ramp rate for 10 minute reserves (MW)
    double ramp_30; // ramp rate for 30 minute reserves (MW)
    double ramp_q; //ramp rate for reactive power (2 sec timescale) (MVAr/min)
    double apf; //APF, area participation factor

    double ramp_rate;// MW/min
    double ramp_time;//in minutes
    bool is_black; //true=is black start
    double crank_time;// in minutes
    double crank_p; // MW
    double earliest_time;
    double latest_time;
    double restoration_priority;
    
    Generator(
        int _gen_num,
        int _host_bus,
        double _Pg,
        double _Qg,
        double _q_max,
        double _q_min,
        double _Vg,
        double _mBase,
        int _status,
        double _p_max,
        double _p_min,
        double _Pc1,
        double _Pc2,
        double _Qc1min,
        double _Qc1max,
        double _Qc2min,
        double _Qc2max,
        double _ramp_agc,
        double _ramp_10,
        double _ramp_30,
        double _ramp_q,
        double _apf)
       :gen_num(_gen_num),host_bus(_host_bus),Pg(_Pg), Qg(_Qg), q_max(_q_max), q_min(_q_min), Vg(_q_min), mBase(_mBase),    status(_status), p_max(_p_max), p_min(_p_min), Pc1(_Pc1), Pc2(_Pc2), Qc1min(_Qc1min), Qc1max(_Qc1max), Qc2min(_Qc2min),
    Qc2max(_Qc2max), ramp_agc(_ramp_agc), ramp_10(_ramp_10), ramp_30(_ramp_30), ramp_q(_ramp_q),apf(_apf){};
    
	int getRampTime() { return ramp_time;/* / TIME_BUCKET;*/ };
	int getCrankTime() { return crank_time;/* / TIME_BUCKET;*/ };
    
    
    
};

struct Bus{
    int bus_num;
    int bus_type;
    double Pd;
    double Qd;
    double Gs; //shunt conductance
    double Bs; // shunt susceptance
    int area_id;
    double Vm; // voltage magnitute
    double Va; // voltage angle
    double baseKV; // base voltage (kV)
    int zone; //loss zone (positive integer)
    double Vmax; //maximum voltage magnitude (p.u.)
    double Vmin; //minimum voltage magnitude (p.u.)
    //restoration data
    double critical_load;
    SET_LINE out_arcs;
    SET_LINE in_arcs;
/*    Bus (
         int _bus_num,int _bus_type, double _Pd, double _Qd, double _Gs, double _Bs, int _area_id, double _Vm, double _Va, double _baseKV,
         int _zone, double _Vmax, double _Vmin);*/
    
    Bus (
         int _bus_num,int _bus_type, double _Pd, double _Qd, double _Gs, double _Bs, int _area_id, double _Vm, double _Va, double _baseKV,int _zone, double _Vmax, double _Vmin):
     bus_num(_bus_num), bus_type(_bus_type), Pd (_Pd), Qd(_Qd), Gs(_Gs), Bs(_Bs), area_id(_area_id), Vm(_Vm), Va(_Va), baseKV(_baseKV),zone(_zone), Vmax(_Vmax), Vmin(_Vmin)
    {};
};

struct Branch{
    int branch_num, from_bus,to_bus;
    double r, x, b, rateA, rateB, rateC, ratio, angle, status, angle_min, angle_max;
    Branch(
           int _line_num, int _fbus,int _tbus,
           double _r, double _x, double _b,double _rateA,double _rateB,double _rateC,double _ratio,double _angle,double _status,double _angmin, double _angmax
           ): branch_num(_line_num), from_bus(_fbus),to_bus(_tbus), r(_r), x(_x), b(_b), rateA(_rateA), rateB(_rateB), rateC(_rateC), ratio(_ratio), angle(_angle), status(_status), angle_min(_angmin), angle_max(_angmax){};
    
    
};



#define nCr(n, r) (factorial(n) / factorial(n-r) / factorial(r))
#define nPr(n, r) (factorial(n) / factorial(n-r))
#define SMALL_DECIMAL 0.000001


static std::string ToString(int n)
{
	std::ostringstream oss;
	oss << n;
	std::string s = oss.str();
	return s;
}

static std::string FloatToString(float n)
{
	std::ostringstream oss;
	oss << n;
	std::string s = oss.str();
	return s;
}

static int factorial(int x) {
	int fac = 1;
	for (int i=2; i<=x; i++) fac *= i;
	return fac;
}


ILOSTLBEGIN



static void free_and_null (char **ptr)
{
	if ( *ptr != NULL ) {
		free (*ptr);
		*ptr = NULL;
	}
}

static int getInstanceIntParameterValue(std::string str, std::string paraName)
{
	std::size_t found1 = str.find_first_of(paraName.c_str());
	std::size_t found2=str.find_first_of("rnmkwcs",found1+1);
	std::string paraValue = str.substr(found1+1, found2- found1-1);
	return atoi(paraValue.c_str());

}
static double getInstanceFloatParameterValue(std::string str, std::string paraName)
{
	std::size_t found1 = str.find_first_of(paraName.c_str());
	std::size_t found2=str.find_first_of("rnmkwcs",found1+1);
	std::string paraValue = str.substr(found1+1, found2- found1-1);
	return atof(paraValue.c_str());

}

#endif 
