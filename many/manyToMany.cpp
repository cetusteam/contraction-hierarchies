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
#include <sys/types.h>
#include <fcntl.h>
#include <google/protobuf/io/coded_stream.h>
#include <google/protobuf/io/zero_copy_stream.h>
#include <google/protobuf/io/zero_copy_stream_impl.h>

using namespace std;
using namespace google::protobuf::io;

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
#include "protos/ch/many2many.pb.h"

typedef datastr::graph::SearchGraph MyGraph;
typedef DijkstraCHManyToManyFW DijkstraManyToManyFW;
typedef DijkstraCHManyToManyBW DijkstraManyToManyBW;

static void readNodes(const MyGraph *const g, const google::protobuf::RepeatedField<uint32_t>& in, vector<NodeID>& vs) {
    vs.resize(in.size());
    for (auto i=0; i < in.size(); i++) {
        // map the actual node ID to the node ID that is used internally by
        // highway-node routing
        vs[i] = g->mapExtToIntNodeID(in.Get(i));
    }
}

static void readMessage(const string& filename, ::google::protobuf::Message& message) {
    int fd = open(filename.c_str(), O_RDONLY);
    FileInputStream raw_input(fd);
    raw_input.SetCloseOnDelete(true);
    message.ParseFromZeroCopyStream(&raw_input);
}

static void writeMessage(const string& filename, ::google::protobuf::Message& message) {
    int fd = open(filename.c_str(), O_WRONLY);
    FileOutputStream raw_output(fd);
    raw_output.SetCloseOnDelete(true);
    message.SerializeToZeroCopyStream(&raw_output);
}

MyGraph *const readGraph(const string& inGraphFn) {
    ifstream inGraph(inGraphFn, ios::binary);
    if (!inGraph) {
        cerr << "Input file '" << inGraphFn << "' couldn't be opened." << endl;
        exit(-1);
    }
    return new MyGraph(inGraph);
}


/** The main program. */
int main(int argc, char *argv[])
{
    // read arguments
    if (argc < 2) {
        cerr << "Too few arguments!" << endl
             << "Usage: ./manyToMany <inSgrGraph> <inNodes> <outMatrix>" << endl;
        exit(-1);
    }

    const string inGraphFn   = argv[1];
    const string inNodesFn   = argv[2];
    const string outMatrixFn = argv[3];

    const bool writeSourceNodes = false;
    const bool writeMatrix = false;
    const bool validateResult = false;

    VERBOSE( cerr << "read graph from '" << inGraphFn << "'" << endl);
    VERBOSE( cerr << "read nodes from '" << inNodesFn << "'" << endl);
    VERBOSE( cerr << "write matrix to '" << outMatrixFn << "'" << endl);
    VERBOSE( if (writeSourceNodes) cerr << "write source nodes to cerr" << endl );
    VERBOSE( if (writeMatrix)      cerr << "write matrix to cerr" << endl );
    VERBOSE( if (validateResult)   cerr << "validate result" << endl );

    uint32_t sz_sources, sz_targets;
    cin >> sz_sources >> sz_targets;

    // read graph
    MyGraph *const graph = readGraph(inGraphFn);

    VERBOSE( cerr << "Graph read." << endl );

    // read the input
    ch::protos::Nodes nodes;
    readMessage(inNodesFn, nodes);

    vector<NodeID> sources;
    readNodes(graph, nodes.sources(), sources);

    vector<NodeID> targets;
    readNodes(graph, nodes.targets(), targets);

    // create many-to-many object
    LevelID earlyStopLevel = 10;
    ManyToMany<MyGraph, DijkstraManyToManyFW, DijkstraManyToManyBW, performBucketScans> mtm(graph, earlyStopLevel);

    // compute matrix
    Matrix<EdgeWeight> matrix(sources.size(), targets.size());
    mtm.computeMatrix( sources, targets, matrix );

    // generate and write output
    ch::protos::Matrix outMatrix;
    for (int i=0; i < sources.size(); i++) {
        ch::protos::Row* row = outMatrix.add_rows();
        for (int j=0; j < targets.size(); j++) {
            row->add_data(matrix.value(i, j));
        }
    }
    writeMessage(outMatrixFn, outMatrix);

    if (validateResult) {
        // compute reference solution
        Matrix<EdgeWeight> matrixRef(sources.size(), targets.size());
        mtm.computeMatrixNaive( sources, targets, matrixRef);

        // check
        if (matrix == matrixRef) {
            VERBOSE( cerr << "Solution validated." << endl );
        } else {
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
