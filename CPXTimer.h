#ifndef CPXTIMER_H_HEADER_INCLUDED_BECEF599
#define CPXTIMER_H_HEADER_INCLUDED_BECEF599
#include <ilcplex/cplex.h>
#include <ilconcert/iloenv.h>

class CPXTimer
{
public:
	CPXTimer(void);
	~CPXTimer(void);
	CPXENVptr cpxenv;
	int initilizeCPXenv(int CPX_SCRIND = 0);

	double startingtime;
	double startTimer(void);
	double getElapsedTime(void);
};

#endif



