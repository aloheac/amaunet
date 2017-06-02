//
// Created by loheac on 2/28/17.
//

#ifndef AMAUNETC_EXPRESSIONSERIALIZATION_H
#define AMAUNETC_EXPRESSIONSERIALIZATION_H

#include <iostream>
#include <fstream>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/string.hpp>
#include <boost/serialization/shared_ptr.hpp>
#include "PTSymbolicObjects.h"

int saveSumToFile( Sum &expr, std::string filename );

Sum loadSumFromFile( std::string filename );

int splitSumToFiles( Sum &expr, int blockSize, std::string saveDir );

Sum loadAndEvaluateSumFromFiles( std::string saveDir, int numberOfFiles, int EXPANSION_ORDER_IN_A, int POOL_SIZE );

#endif //AMAUNETC_EXPRESSIONSERIALIZATION_H
