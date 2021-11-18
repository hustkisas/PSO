#include "CPXTimer.h"


CPXTimer::CPXTimer(void)
{
	initilizeCPXenv();
}

CPXTimer::~CPXTimer(void)
{
}

double CPXTimer::startTimer(void)
{
	double timestamp;
	CPXgettime(cpxenv,&timestamp);
	startingtime = timestamp;
	return timestamp;
}

double CPXTimer::getElapsedTime(void)
{
	double timestamp;
	CPXgettime(cpxenv,&timestamp);
	return timestamp - startingtime;
}


int CPXTimer::initilizeCPXenv( int CPX_SCRIND)
{
	cpxenv = NULL;
	int           status = 0;

	/* Initialize the CPLEX environment */
	cpxenv = CPXopenCPLEX (&status);

	if ( cpxenv == NULL ) {
		char  errmsg[1024];
		fprintf (stderr, "Could not open CPLEX environment.\n");
		CPXgeterrorstring (cpxenv, status, errmsg);
		fprintf (stderr, "%s", errmsg);
		return 0;
	}

	/* Turn on output to the screen */
	status = CPXsetintparam (cpxenv, CPX_PARAM_SCRIND, CPX_SCRIND);
	if ( status ) {
		fprintf (stderr, 
			"Failure to turn on screen indicator, error %d.\n", status);
		return 0;
	}

	status = CPXsetintparam(cpxenv, CPX_PARAM_THREADS, 1);
	status = CPXsetintparam(cpxenv, CPX_PARAM_CLOCKTYPE, 1);

	return 1;
}


