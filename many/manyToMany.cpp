/* Copyright (C) 2005, 2006, 2007, 2008 
 * Robert Geisberger, Dominik Schultes, Peter Sanders,
 * Universitaet Karlsruhe (TH)
 *
 * This file is part of Contraction Hierarchies.
 *
 * Contraction Hierarchies is free software; you can redistribute it
 * and/or modify it under the terms of the GNU Affero General Public License
 * as published by the Free Software Foundation; either version 3 of
 * the License, or (at your option) any later version.
 *
 * Contraction Hierarchies is distributed in the hope that it will be
 * useful, but WITHOUT ANY WARRANTY; without even the implied warranty
 * of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License
 * along with Contraction Hierarchies; see the file COPYING; if not,
 * see <http://www.gnu.org/licenses/>.
 */

#include <iostream>
#include <iomanip>
#include <fstream>
using namespace std;

/**
 * Performing the bucket scans is a crucial part of the many-to-many computation.
 * It can be switched off to allow time measurements of the forward search
 * without accounting for the bucket scans.
 */
const bool performBucketScans = true;


#include "../config.h"
#include "../stats/utils.h"
Counter counter;
//#include "../datastr/general.h"
#include "../datastr/graph/graph.h"
#include "../datastr/graph/SearchGraph.h"

typedef datastr::graph::SearchGraph TransitGraph;

//#include "../stats/log.h"
#include "../processing/DijkstraCH.h"
#include "manyToMany.h"



typedef datastr::graph::SearchGraph MyGraph;
typedef DijkstraCHManyToManyFW DijkstraManyToManyFW;
typedef DijkstraCHManyToManyBW DijkstraManyToManyBW;
inline NodeID mapNodeID(const MyGraph *const g, const NodeID u) {
    // map the actual node ID to the node ID that is used internally by highway-node routing
    return g->mapExtToIntNodeID(u);
}


/** Returns a random node identifier between 0 and n-1. */
inline NodeID randomNodeID(NodeID n) {
    NodeID x = (NodeID)(rand() / (double)(RAND_MAX+1.0) * n);
    cerr << "generated rand node_id:" << x << endl;
    return x;
}


/** The main program. */
int main(int argc, char *argv[])
{
    // read arguments
    if (argc < 2) {
        cerr << "Too few arguments!" << endl
             << "Usage: ./manyToMany <filenameGraph> [<noOfSources>] [<noOfTargets>] [<flags>]" << endl
             << "       noOfSources: default = 1000" << endl
             << "       noOfTargets: default = 1000" << endl
             << "       flags: 1 = write source node IDs to cerr;" << endl
             << "              2 = write matrix to cerr;" << endl
             << "              4 = validate result;" << endl;
        exit(-1);
    }
    srand(time(NULL));

    const string filenameGraph = argv[1];

    NodeID noOfSources = 1000;
    if (argc >= 3) noOfSources = atoi(argv[2]);

    NodeID noOfTargets = 1000;
    if (argc >= 4) noOfTargets = atoi(argv[3]);

    LevelID earlyStopLevel = 10;

    bool writeSourceNodes = false;
    bool writeMatrix = false;
    bool validateResult = false;
    if (argc >= 5) {
        const int flags = atoi(argv[4]);
        writeSourceNodes = ((flags & 1) == 1);
        writeMatrix      = ((flags & 2) == 2);
        validateResult   = ((flags & 4) == 4);
    }

    VERBOSE( cerr << "read graph from '" << filenameGraph << "'" << endl
                  << "generate " << noOfSources << " random source nodes and " << noOfTargets << " random target nodes" << endl );
    VERBOSE( if (writeSourceNodes) cerr << "write source nodes to cerr" << endl );
    VERBOSE( if (writeMatrix)      cerr << "write matrix to cerr" << endl );
    VERBOSE( if (validateResult)   cerr << "validate result" << endl );
    VERBOSE( cerr << endl );


    // read graph
    ifstream inGraph( filenameGraph.c_str(), ios::binary );
    if (! inGraph) {
        cerr << "Input file '" << filenameGraph << "' not found." << endl;
        exit(-1);
    }
    MyGraph *const graph = new MyGraph(inGraph);
    inGraph.close();
    VERBOSE( cerr << "Graph read." << endl );


    // create many-to-many object
    ManyToMany<MyGraph, DijkstraManyToManyFW, DijkstraManyToManyBW, performBucketScans> mtm(graph, earlyStopLevel);

    // prepare input
    const NodeID noOfNodes = graph->noOfNodes();
    
    cerr << "Generating sources" << endl;
    vector<NodeID> sources;
    for (NodeID i = 0; i < noOfSources; i++) sources.push_back(mapNodeID(graph, randomNodeID(noOfNodes)));

    vector<NodeID> targets;
    cerr << "Generating targets" << endl;
    for (NodeID i = 0; i < noOfTargets; i++) targets.push_back(mapNodeID(graph, randomNodeID(noOfNodes)));

    if (writeSourceNodes) {
        for (NodeID i = 0; i < noOfSources; i++) {
            if (i > 0) cerr << endl;
            cerr << sources[i];
        }
    }

    // compute matrix
    Matrix<EdgeWeight> matrix(noOfSources, noOfTargets);
    mtm.computeMatrix( sources, targets, matrix );

    if (validateResult) {
        // compute reference solution
        Matrix<EdgeWeight> matrixRef(noOfSources, noOfTargets);
        mtm.computeMatrixNaive( sources, targets, matrixRef );

        // check
        if (matrix == matrixRef) {
            VERBOSE( cerr << "Solution validated." << endl );
        }
        else {
            cerr << "Wrong solution!" << endl;
        }
    }

    if (writeMatrix) {
        VERBOSE( cerr << "write matrix" << endl );
        cerr << matrix;
    }

    delete graph;
}

// doesn't look nice, but required by the compiler (gcc 4)
const EdgeWeight Weight::MAX_VALUE;
