#include "DataGenerator.h"
#include <fstream>
#include "randomc.h"

using namespace std;

DataGenerator::DataGenerator(void)
{
}

DataGenerator::~DataGenerator(void)
{
}


void DataGenerator::generateProbDataWithoutResConst(string fileName)
{

	int A_MAX_RAND = m;
	int A_MIN_RAND = 1;
	fstream file;
	file.open(fileName.c_str(), ios::out);
	file<<n<<endl;
	file<<m<<endl;
	file<<k<<endl;

	vector<double> rhs;
	vector<VEC_DOUBLE> matrixA;
	CRandomMersenne** pRandGen = new CRandomMersenne*[n] ;
	for (int i = 0; i < n; i++)
	{
		pRandGen[i]= new CRandomMersenne(seed + i*3);
	}


	for (int i = 0; i<m; i++)
	{
		for (int t = 0 ; t < nr; t++)
		{
			double rhsV=0;
			while (rhsV >= -0.001 && rhsV <= 0.001)
			{		
				rhsV = A_MIN_RAND + (int)((A_MAX_RAND-A_MIN_RAND) * pRandGen[0]->Random());
			}
			rhs.push_back(rhsV);
			VEC_DOUBLE row_i;
			for (int j=0; j<n; j++)
			{
				double cof =0;
				if (j == 0 )
				{
					cof = A_MIN_RAND + (int)((A_MAX_RAND - A_MIN_RAND)*pRandGen[j]->Random());//randGen.IRandom(1,2);//randGen.IRandom(MIN_RAND_NUM,MAX_RAND_NUM); //1;
				}
				else
				{
					cof = A_MIN_RAND + (int)((A_MAX_RAND - A_MIN_RAND)*pRandGen[j]->Random());
				}
				row_i.push_back(cof);
			}
			matrixA.push_back(row_i);
		}
	}

	//generate cost vector
	CRandomMersenne randGenCost(seed+20);
	vector<double> costVec;
	for (int i = 0; i < n; i++)
	{
		costVec.push_back(double(randGenCost.IRandom(A_MIN_RAND , A_MAX_RAND)));
	}


	// write cost vector
	for (int i  = 0 ; i < n; i++)
	{
		file<<costVec[i]<<"\t";
	}
	file<<endl;

	for (int i =0; i < matrixA.size() ; i ++)// write sampling matrix
	{
		//write x
		for (int j = 0; j < n; j ++)
		{
			file<<matrixA[i][j]<<"\t";
		}
		// write sign
		file<<"G"<< "\t";
		// write rhs
		file<<rhs[i]<<endl;
	}



	file.close();

	cout<< "file generated " <<endl;

	for (int i = 0; i < n ; i++)
	{
		delete pRandGen[i];
	}

	delete []pRandGen;
}

void DataGenerator::generateProbDataWithoutResConstMultiRows(string fileName)
{

	nr = r;
	int A_MAX_RAND = 200;
	int A_MIN_RAND = 20;
	fstream file;
	file.open(fileName.c_str(), ios::out);
	file<<n<<endl;
	file<<m<<endl;
	file<<k<<endl;

	vector<double> rhs;
	vector<VEC_DOUBLE> matrixA;
	CRandomMersenne** pRandGen = new CRandomMersenne*[n] ;
	for (int i = 0; i < n; i++)
	{
		pRandGen[i]= new CRandomMersenne(seed + i*3);
	}


	for (int i = 0; i<m; i++)
	{
		for (int t = 0 ; t < nr; t++)
		{
			double rhsV=0;
			while (rhsV >= -0.001 && rhsV <= 0.001)
			{		
				rhsV = A_MIN_RAND + (int)((A_MAX_RAND-A_MIN_RAND) * pRandGen[0]->Random());
			}
			rhs.push_back(rhsV);
			VEC_DOUBLE row_i;
			for (int j=0; j<n; j++)
			{
				double cof =0;
				if (j == 0 )
				{
					cof = A_MIN_RAND + (int)((A_MAX_RAND - A_MIN_RAND)*pRandGen[j]->Random());//randGen.IRandom(1,2);//randGen.IRandom(MIN_RAND_NUM,MAX_RAND_NUM); //1;
				}
				else
				{
					cof = A_MIN_RAND + (int)((A_MAX_RAND - A_MIN_RAND)*pRandGen[j]->Random());
				}
				row_i.push_back(cof);
			}
			matrixA.push_back(row_i);
		}
	}

	//generate cost vector
	CRandomMersenne randGenCost(seed+20);
	vector<double> costVec;
	for (int i = 0; i < n; i++)
	{
		costVec.push_back(double(randGenCost.IRandom(A_MIN_RAND , A_MAX_RAND)));
	}


	// write cost vector
	for (int i  = 0 ; i < n; i++)
	{
		file<<costVec[i]<<"\t";
	}
	file<<endl;

	for (int i =0; i < matrixA.size() ; i ++)// write sampling matrix
	{
		//write x
		for (int j = 0; j < n; j ++)
		{
			file<<matrixA[i][j]<<"\t";
		}
		// write sign
		file<<"G"<< "\t";
		// write rhs
		file<<rhs[i]<<endl;
	}

	file.close();

	cout<< "file generated " <<endl;

	for (int i = 0; i < n ; i++)
	{
		delete pRandGen[i];
	}

	delete []pRandGen;



}

void DataGenerator::generateProbDataWithResConstSingleRow(string fileName)
{
	double rhs_MIN = r;
	double A_MAX_RAND = 1.5;
	double A_MIN_RAND = 0.8;
	fstream file;
	file.open(fileName.c_str(), ios::out);
	file<<n<<endl;
	file<<m<<endl;
	file<<k<<endl;

	vector<double> rhs;
	vector<VEC_DOUBLE> matrixA;
	CRandomMersenne** pRandGen = new CRandomMersenne*[n] ;
	for (int i = 0; i < n; i++)
	{
		pRandGen[i]= new CRandomMersenne(seed + i*3);
	}


	for (int i = 0; i<m; i++)
	{
		for (int t = 0 ; t < nr; t++)
		{
			double rhsV=0;
			while (rhsV >= -0.001 && rhsV <= 0.001)
			{		
				rhsV =rhs_MIN; /*A_MIN_RAND + (int)((A_MAX_RAND-A_MIN_RAND) * pRandGen[0]->Random())*/
			}
			rhs.push_back(rhsV);
			VEC_DOUBLE row_i;
			for (int j=0; j<n; j++)
			{
				double cof =0;
				if (j == 0 )
				{
					cof = A_MIN_RAND + ((A_MAX_RAND - A_MIN_RAND)*pRandGen[j]->Random());//randGen.IRandom(1,2);//randGen.IRandom(MIN_RAND_NUM,MAX_RAND_NUM); //1;
					cof = (float)((int)(cof*10000+0.5))/10000.0;
				}
				else
				{
					cof = A_MIN_RAND + ((A_MAX_RAND - A_MIN_RAND)*pRandGen[j]->Random());
					cof = (float)((int)(cof*10000+0.5))/10000.0;
				}
				row_i.push_back(cof);
			}
			matrixA.push_back(row_i);
		}
	}

	//generate cost vector
	CRandomMersenne randGenCost(seed+20);
	vector<double> costVec;
	for (int i = 0; i < n; i++)
	{
		costVec.push_back(double(randGenCost.IRandom(1 , 50)));
	}


	// write cost vector
	for (int i  = 0 ; i < n; i++)
	{
		file<<costVec[i]<<"\t";
	}
	file<<endl;

	for (int i =0; i < matrixA.size() ; i ++)// write sampling matrix
	{
		//write x
		for (int j = 0; j < n; j ++)
		{
			file<<matrixA[i][j]<<"\t";
		}
		// write sign
		file<<"G"<< "\t";
		// write rhs
		file<<rhs[i]<<endl;
	}

	file<<endl;

	file.close();

	cout<< "file generated " <<endl;

	for (int i = 0; i < n ; i++)
	{
		delete pRandGen[i];
	}

	delete []pRandGen;
}


void DataGenerator::generateProbDataSingleRowTest(string fileName)
{
	fstream file;
	file.open(fileName.c_str(), ios::out);
	file<<n<<endl;
	file<<m<<endl;
	file<<k<<endl;

	vector<double> rhs;
	vector<VEC_DOUBLE> matrixA;
	CRandomMersenne** pRandGen = new CRandomMersenne*[n] ;
	for (int i = 0; i < n; i++)
	{
		pRandGen[i]= new CRandomMersenne(seed + i*3);
	}


	for (int i = 0; i<m; i++)
	{
		for (int t = 0 ; t < nr; t++)
		{
			double rhsV=0;
			while (rhsV >= -0.001 && rhsV <= 0.001)
			{		
				rhsV =2; 
			}
			rhs.push_back(rhsV);
			VEC_DOUBLE row_i;
			for (int j=0; j<n; j++)
			{
				double cof =0;
				cof = pRandGen[j]->IRandom(1,2) == 1 ? 1:3   ;
				row_i.push_back(cof);
			}
			matrixA.push_back(row_i);
		}
	}

	//generate cost vector
	CRandomMersenne randGenCost(seed+20);
	vector<double> costVec;
	for (int i = 0; i < n; i++)
	{
		costVec.push_back(1);
	}


	// write cost vector
	for (int i  = 0 ; i < n; i++)
	{
		file<<costVec[i]<<"\t";
	}
	file<<endl;

	for (int i =0; i < matrixA.size() ; i ++)// write sampling matrix
	{
		//write x
		for (int j = 0; j < n; j ++)
		{
			file<<matrixA[i][j]<<"\t";
		}
		// write sign
		file<<"G"<< "\t";
		// write rhs
		file<<rhs[i]<<endl;
	}

	file.close();
	cout<< "file generated " <<endl;

	for (int i = 0; i < n ; i++)
	{
		delete pRandGen[i];
	}
	delete []pRandGen;
}



void DataGenerator::generateIntegerDataSingleRowWithUnitCostVector(string fileName)
{
	double rhs_MIN = r;
	double A_MAX_RAND = 2;
	double A_MIN_RAND = 0.00;
	fstream file;
	file.open(fileName.c_str(), ios::out);
	file<<n<<endl;
	file<<m<<endl;
	file<<k<<endl;

	vector<double> rhs;
	vector<VEC_DOUBLE> matrixA;
	// genearte random generators
	CRandomMersenne** pRandGen = new CRandomMersenne*[n] ;
	for (int i = 0; i < n; i++)
	{
		pRandGen[i]= new CRandomMersenne(seed + i*3);
	}


	// generate matrix A
	for (int i = 0; i<m; i++)
	{
		for (int t = 0 ; t < nr; t++)
		{
			double rhsV=r;
			rhs.push_back(rhsV);
			VEC_DOUBLE row_i;
			for (int j=0; j<n; j++)
			{
				double cof =0;
				cof = pRandGen[j]->Random() <= 0.5 ? A_MIN_RAND : A_MAX_RAND; //A_MIN_RAND + ((A_MAX_RAND - A_MIN_RAND)*pRandGen[j]->Random());
				row_i.push_back(cof);
			}
			matrixA.push_back(row_i);
		}
	}

	//generate cost vector
	CRandomMersenne randGenCost(seed+20);
	vector<double> costVec;
	for (int i = 0; i < n; i++)
	{
		costVec.push_back(1/*double(randGenCost.IRandom(1 , 50))*/);
	}


	// write cost vector
	for (int i  = 0 ; i < n; i++)
	{
		file<<costVec[i]<<"\t";
	}
	file<<endl;

	// write matrix A 
	for (int i =0; i < matrixA.size() ; i ++)
	{
		//write x
		for (int j = 0; j < n; j ++)
		{
			file<<matrixA[i][j]<<"\t";
		}
		// write sign
		file<<"G"<< "\t";
		// write rhs
		file<<rhs[i]<<endl;
	}
}


void DataGenerator::generateProbDataSingleRowWithUnitCostVector(string fileName)
{
	double rhs_MIN = r;
	double A_MAX_RAND = 1.5;
	double A_MIN_RAND = 0.8;
	fstream file;
	file.open(fileName.c_str(), ios::out);
	file<<n<<endl;
	file<<m<<endl;
	file<<k<<endl;

	vector<double> rhs;
	vector<VEC_DOUBLE> matrixA;
	// genearte random generators
	CRandomMersenne** pRandGen = new CRandomMersenne*[n] ;
	for (int i = 0; i < n; i++)
	{
		pRandGen[i]= new CRandomMersenne(seed + i*3);
	}


	// generate matrix A
	for (int i = 0; i<m; i++)
	{
		for (int t = 0 ; t < nr; t++)
		{
			double rhsV=r;
			rhs.push_back(rhsV);
			VEC_DOUBLE row_i;
			for (int j=0; j<n; j++)
			{
				double cof =0;
				cof = A_MIN_RAND + ((A_MAX_RAND - A_MIN_RAND)*pRandGen[j]->Random());
				row_i.push_back(cof);
			}
			matrixA.push_back(row_i);
		}
	}

	//generate cost vector
	CRandomMersenne randGenCost(seed+20);
	vector<double> costVec;
	for (int i = 0; i < n; i++)
	{
		costVec.push_back(1/*double(randGenCost.IRandom(1 , 50))*/);
	}


	// write cost vector
	for (int i  = 0 ; i < n; i++)
	{
		file<<costVec[i]<<"\t";
	}
	file<<endl;

	// write matrix A 
	for (int i =0; i < matrixA.size() ; i ++)
	{
		//write x
		for (int j = 0; j < n; j ++)
		{
			file<<(float)((int)(matrixA[i][j]/rhs[i]*10000+0.5))/10000.0 /*matrixA[i][j]*/<<"\t";
		}
		// write sign
		file<<"G"<< "\t";
		// write rhs
		file<<rhs[i]<<endl;
	}
}

