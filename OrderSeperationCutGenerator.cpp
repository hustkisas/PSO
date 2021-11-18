#include "Global.h"
#include "OrderSeperationCutGenerator.h"
#include "OptModel.h"
#include "DataModel.h"

OrderSeperationCutGenerator::OrderSeperationCutGenerator(OptModel* _model)
{
	model = _model;
}

OrderSeperationCutGenerator::~OrderSeperationCutGenerator(void)
{
}

bool OrderSeperationCutGenerator::seperateCurrentSolution()
{

	/* Initilize model and create variables */
	IloEnv env = IloEnv();
	IloModel SeperationModel(env, "Seperation");
	IloCplex SeperationSolver(SeperationModel);
	int dimension = DataModel::n;
	IloNumVarArray alpha(env, dimension, 0, IloInfinity, ILOFLOAT);
	for (int i = 0; i < alpha.getSize(); i++)
	{
		string name = "a" + ToString(i);
		alpha[i].setName(name.c_str());
	}
	IloNumVar beta(env, 0, IloInfinity);
	beta.setName("b");

	// set Speration objective
	double* x = new double[DataModel::n];
	CPXgetx(model->cpxenv, model->cpxlp, x, 0, DataModel::n -1);
	IloExpr expr(env);
	for (int i = 0; i < dimension; i++)
	{
		expr += x[i]*alpha[i];
	}
	SeperationModel.add(IloMinimize(env, expr - beta));
	delete []x;

	IloNumVarArray betaArray(env, DataModel::m, 0, IloInfinity, ILOFLOAT);
	IloNumVarArray z(env, DataModel::m,0, 1, IloNumVar::Bool);
	for (int i = 0; i < DataModel::m; i++)
	{
		string strBetai= "b" + ToString(i);
		betaArray[i].setName(strBetai.c_str());
		string strZi = "z" + ToString(i);
		z[i].setName(strZi.c_str());
	}
	
	/**************************   build model  ******************************/

	DataModel*dm = model->dm_p;
	for (int i = 0 ; i < DataModel::m ; i++)
	{
		double b_i = dm->rhs[i];
		for (int j = 0; j < DataModel::n; j++)
		{
			double a_ij = dm->Arows[i][j].second;
			SeperationModel.add(betaArray[i]<=(b_i/a_ij)*alpha[j]);
		}
	}

	int BigM = 1/SMALL_DECIMAL;
	for (int i = 0; i < DataModel::m; i ++)
	{
		SeperationModel.add(beta <= betaArray[i] + BigM*z[i]);
	}

	SeperationModel.add(IloSum(z)==DataModel::m - DataModel::k -1);
	SeperationModel.add(IloSum(alpha) == 1);
	SeperationSolver.exportModel("seperation.lp");
	SeperationSolver.solve();
	if(!SeperationSolver.getCplexStatus() == IloCplex::Optimal)
		cout << "Seperation is not solved"<<endl;
	if (SeperationSolver.getObjValue() >= - SMALL_DECIMAL)
	{
		return false;
	}
	IloNumArray val_alpha(env);
	IloNum val_beta;
	IloNumArray val_betaArray(env);
	IloNumArray val_z(env);
	SeperationSolver.getValues(val_z, z);
	SeperationSolver.getValues(val_betaArray, betaArray);
	SeperationSolver.getValues(val_alpha, alpha);
	val_beta = SeperationSolver.getValue(beta);
	cout << "beta: " << val_beta<<endl;
	cout << "beta_i" << val_betaArray<<endl;
	cout << "z: "  << val_z <<endl;
	cout << "alpha" << val_alpha<<endl;
	CUT_COEFICIENTS coefs; CUT_RHS rhs = val_beta;
	for (int i = 0; i < val_alpha.getSize(); i++)
	{
		coefs.push_back(make_pair(i, val_alpha[i]));
	}
	CUT cut(coefs, rhs);
//	model->addAConstToCPXModel(cut, string("test"), 0);
	return true;
}
