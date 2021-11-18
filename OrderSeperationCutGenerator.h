#ifndef ORDERSEPERATIONCUTGENERATOR_HEADER_INCLUDED_BECEF599
#define ORDERSEPERATIONCUTGENERATOR_HEADER_INCLUDED_BECEF599
class OptModel;
class OrderSeperationCutGenerator
{
public:
	OrderSeperationCutGenerator(OptModel* _model);
public:
	~OrderSeperationCutGenerator(void);

	bool seperateCurrentSolution();


	OptModel* model ;
};


#endif

