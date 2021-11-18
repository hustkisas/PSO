#ifndef INTERSECTIONCUT_H
#define INTERSECTIONCUT_H

#include "Global.h"


class IntersectionCut{

public:
	string slackVarName;
	double rhs;
	map<int,double> ind_coefs;
	map<string, double> name_coefs;
	map<string, double> origin_var_coefs;
	IntersectionCut(string sName, double _rhs):slackVarName(sName),  rhs(_rhs) {}
	void addName_coef(string varName, double coef){
		name_coefs.insert(make_pair(varName, coef));
	}
	void addInd_ceof(int index, double coef){
		ind_coefs.insert(make_pair(index, coef));
	}

	bool convertindex2VarName(CPXCENVptr& cpxenv, CPXLPptr& lp){
		int num_cols  = CPXgetnumcols(cpxenv,  lp) ; // need to know why
		int surplus;

		int status = CPXgetcolname (cpxenv, lp, NULL, NULL, 0, &surplus, 0,num_cols-1);
		int cur_colnamespace = - surplus;
		if ( cur_colnamespace > 0 ) {
			char          **cur_colname      = (char **) malloc (sizeof(char *)*(num_cols));
			char          *cur_colnamestore = (char *)  malloc (cur_colnamespace);
			if ( cur_colname      == NULL ||
				cur_colnamestore == NULL   ) {
					fprintf (stderr, "Failed to get memory for column names.\n");
					status = -1;
					return false;
			}
			status = CPXgetcolname (cpxenv, lp, cur_colname, cur_colnamestore, 
				cur_colnamespace, &surplus, 0, num_cols-1);
			// map column index to var name
			for (map<int,double>::iterator itr = ind_coefs.begin(); itr != ind_coefs.end(); ++itr)
			{
				int index = itr->first;
				double v_coef = itr->second;
				string varName = cur_colname[index];
				addName_coef(varName, v_coef);
			}
			free_and_null( (char **)&cur_colname);
			free_and_null((char **)&cur_colnamestore);
		}
		return true;
	}

	
};



#endif 
