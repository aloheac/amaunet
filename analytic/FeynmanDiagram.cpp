//
// Created by loheac on 2/22/17.
//

#include <sstream>
#include <vector>
#include <set>
#include <algorithm>
#include <iostream>
#include "PathIntegration.h"
#include "FeynmanDiagram.h"

using namespace std;

Vertex::Vertex( unsigned int id ) : vertexID( id ) { }

void Vertex::connectTo( unsigned int id ) {
    connectedVertices.push_back( id );
    sort( connectedVertices.begin(), connectedVertices.end() );
}

unsigned int Vertex::getID() const {
    return vertexID;
}

void Vertex::setID( unsigned int id ) {
    vertexID = id;
}

string Vertex::to_string() {
    stringstream ss;
    ss << vertexID << " --> {";

    for ( vector<unsigned int>::iterator iter = connectedVertices.begin(); iter != connectedVertices.end(); ++iter ) {
        ss << " " << *iter << " ";
    }

    ss << "}";
    return ss.str();
}

FeynmanDiagram::FeynmanDiagram() : infinityLoopCount( 0 ) {
    vertices = vector<Vertex>();
}

bool cmp_vertex( const Vertex &a, const Vertex &b ) {
    return a.getID() < b.getID();
}

void FeynmanDiagram::addVertex( Vertex v ) {
    for ( vector<Vertex>::iterator iter = vertices.begin(); iter != vertices.end(); ++iter ) {
        if ( iter->getID() == v.getID() ) {
            cout << "***ERROR: A vertex with ID " << v.getID() << " already exists in Feynman diagram." << endl;
        }
    }

    vertices.push_back( v );

    sort( vertices.begin(), vertices.end(), cmp_vertex );
}

void FeynmanDiagram::addInfinityLoop() {
    infinityLoopCount++;
}

int FeynmanDiagram::getInfinityLoopCount() {
    return infinityLoopCount;
}

bool FeynmanDiagram::connect( unsigned int vertexA, unsigned int vertexB ) {
    int verticesConnected = 0;

    for ( vector<Vertex>::iterator iter = vertices.begin(); iter != vertices.end(); ++iter ) {
        if ( iter->getID() == vertexA ) {
            iter->connectTo( vertexB );
            verticesConnected++;
        }

        if ( iter->getID() == vertexB ) {
            iter->connectTo( vertexA );
            verticesConnected++;
        }
    }

    if ( verticesConnected == 2 ) {
        return true;
    } else {
        return false;
    }
}

void FeynmanDiagram::transformIndices( std::map<unsigned int, unsigned int> transformationMap ) {
    for ( vector<Vertex>::iterator iter = vertices.begin(); iter != vertices.end(); ++iter ) {
        iter->setID( transformationMap[ iter->getID() ] );

        for ( int i = 0; i < iter->connectedVertices.size(); i++ ) {
            iter->connectedVertices[ i ] = transformationMap[ iter->connectedVertices[ i ] ];
        }

        sort( iter->connectedVertices.begin(), iter->connectedVertices.end() );
    }

    sort( vertices.begin(), vertices.end(), cmp_vertex );
}

bool FeynmanDiagram::isIdenticalTo( FeynmanDiagram otherDiagram ) {
    if ( vertices.size() != otherDiagram.vertices.size() ) return false;
    if ( infinityLoopCount != otherDiagram.infinityLoopCount ) return false;

    for ( int i = 0; i < vertices.size(); i++ ) {
        if ( vertices[i].getID() != otherDiagram.vertices[i].getID() ) return false;
        if ( vertices[i].connectedVertices != otherDiagram.vertices[i].connectedVertices ) return false;
    }

    return true;
}

bool FeynmanDiagram::isSimilarTo( FeynmanDiagram otherDiagram ) {
    if ( vertices.size() != otherDiagram.vertices.size() ) return false;

    // Collect the all node indices that are present across all diagrams.
    set<int> validNodeIndicesSet;

    // Loop over elements of the set of nodes in this diagram.
    // Insert them into the set of indices for this diagram.
    for ( vector<Vertex>::iterator iter = vertices.begin(); iter != vertices.end(); ++iter ) {
        validNodeIndicesSet.insert( iter->getID() );
    }

    // Loop over elements of the set of nodes in the other diagram.
    // Insert them into the set of indices for the other diagram.
    for ( vector<Vertex>::iterator iter = otherDiagram.vertices.begin(); iter != otherDiagram.vertices.end(); ++iter ) {
        validNodeIndicesSet.insert( iter->getID() );
    }

    // Copy elements of one of the sets of indices to a vector; we require an ordered set in order to generate
    // permutations from the data structure.
    vector<int> validNodeIndices;
    for ( set<int>::iterator iter = validNodeIndicesSet.begin(); iter != validNodeIndicesSet.end(); ++iter ) {
        validNodeIndices.push_back( *iter );
    }

    vector<int> originalIndexPermutation( validNodeIndices );  // Store the original permutation of indices; we need
                                                               // this to generate the mapping from the original index
                                                               // to the permutation of the index.

    FeynmanDiagram originalDiagram( otherDiagram );

    do {
        // Generate the index combination map.
        map<unsigned int, unsigned int> indexMap;
        for ( unsigned int i = 0; i < validNodeIndices.size(); i++ ) {
            indexMap.insert( pair<unsigned int, unsigned int>( originalIndexPermutation[ i ], validNodeIndices[ i ] ) );
        }

        // Rename node indices of the other diagram according to the generated map.
        otherDiagram.transformIndices( indexMap );

        if ( isIdenticalTo( otherDiagram ) ) return true;

        otherDiagram = originalDiagram;

    } while ( next_permutation( validNodeIndices.begin(), validNodeIndices.end() ) );

    return false;
}

string FeynmanDiagram::to_string() {
    stringstream ss;
    ss << "FeynmanDiagram[";

    for ( vector<Vertex>::iterator iter = vertices.begin(); iter != vertices.end(); ++iter ) {
        ss << " " << iter->to_string() << " ";
    }

    ss << "]";
    return ss.str();
}

FeynmanDiagram constructDiagram( DeltaContractionSet indexSet ) {
    // Iterate over the set of contractions to gather what indices are present.
    set<int> presentIndicesSet;
    for ( vector<IndexContraction>::iterator iter = indexSet.getIteratorBegin(); iter != indexSet.getIteratorEnd(); ++iter ) {
        if ( iter->i != iter->j ) {
            presentIndicesSet.insert( iter->i );
            presentIndicesSet.insert( iter->j );
        }
    }

    // Load the set of indices into a vector such that they're ordered, and sort the indices.
    vector<int> presentIndices;
    for ( set<int>::iterator iter = presentIndicesSet.begin(); iter != presentIndicesSet.end(); ++iter ) {
        presentIndices.push_back( *iter );
    }

    sort( presentIndices.begin(), presentIndices.end() );

    // Create vertices and insert them into the diagram.
    FeynmanDiagram diagram;
    for ( vector<int>::iterator iter = presentIndices.begin(); iter!= presentIndices.end(); ++iter ) {
        diagram.addVertex( Vertex( *iter ) );
    }

    // Connect nodes in the Feynman diagram, and note where "infinity loops" are present, as appropriate.
    for ( vector<IndexContraction>::iterator iter = indexSet.getIteratorBegin(); iter != indexSet.getIteratorEnd(); ++iter ) {
        if ( iter->i == iter->j ) {
            diagram.addInfinityLoop();
        } else {
            diagram.connect( iter->i, iter->j );
        }
    }

    return diagram;
}

bool compareContractionSetsViaDiagrams( DeltaContractionSet setA, DeltaContractionSet setB ) {
    FeynmanDiagram diagramA = constructDiagram( setA );
    FeynmanDiagram diagramB = constructDiagram( setB );

    return diagramA.isSimilarTo( diagramB );
}