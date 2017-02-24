//
// Created by loheac on 2/22/17.
//

#ifndef AMAUNETC_FEYNMANDIAGRAM_H
#define AMAUNETC_FEYNMANDIAGRAM_H

#include <vector>
#include <map>
#include "PathIntegration.h"

class Vertex {

    friend class FeynmanDiagram;

public:

    Vertex( unsigned int id );

    void connectTo( unsigned int id );

    unsigned int getID() const;

    void setID( unsigned int id );

    std::string to_string();

private:

    unsigned int vertexID;

    std::vector<unsigned int> connectedVertices;
};

class FeynmanDiagram {

public:

    FeynmanDiagram();

    void addVertex( Vertex v );

    void addInfinityLoop();

    int getInfinityLoopCount();

    bool connect( unsigned int vertexA, unsigned int vertexB );

    void transformIndices( std::map<unsigned int, unsigned int> transformationMap );

    bool isIdenticalTo( FeynmanDiagram otherDiagram );

    bool isSimilarTo( FeynmanDiagram otherDiagram );

    std::string to_string();

private:

    int infinityLoopCount;

    std::vector<Vertex> vertices;

};

FeynmanDiagram constructDiagram( DeltaContractionSet indexSet );

bool compareContractionSetsViaDiagrams( DeltaContractionSet setA, DeltaContractionSet setB );
        
#endif //AMAUNETC_FEYNMANDIAGRAM_H
