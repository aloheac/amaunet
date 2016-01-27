#include <iostream>
#include "PTSymbolicObjects.h"

using namespace std;


int main( int argc, char** argv ) {
	std::string dat = pt_util::loadExpressionFile( "/home/loheac/cuda-workspace/AmaunetEval/src/DATA.dat" );
	ExpressionInterpreter expint = ExpressionInterpreter();
	Sum expr = expint.parseExpression( dat );
	cout << expr << endl;
}




