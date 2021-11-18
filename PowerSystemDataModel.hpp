//
//  PowerSystemDataModel.hpp
//  Power System Optimization
//
//  Created by QQ on 8/30/17.
//  Copyright © 2017 QQ. All rights reserved.
//

#ifndef PowerSystemDataModel_hpp
#define PowerSystemDataModel_hpp


#include <global.h>



class PowerSystemDataModel
{
public:
    PowerSystemDataModel(void);
    ~PowerSystemDataModel(void);
    
    MAP_BUS buses;
    MAP_BRANCH branches;
    MAP_GEN generators;
    MAP_GEN bs_gens;
    MAP_GEN nbs_gens;
    MAP_BUS criticalLoads;
    
    
    bool readMatPowerFile(std::string fileName);
    bool readRestorationDataFile(std::string filename);
    bool initilizeNetwork(void);
};

#endif /* PowerSystemDataModel_hpp */
