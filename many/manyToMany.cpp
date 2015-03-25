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
#include <unistd.h>
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
#define SCALE 1


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
#include "options_parser.h"

typedef datastr::graph::SearchGraph MyGraph;
typedef DijkstraCHManyToManyFW DijkstraManyToManyFW;
typedef DijkstraCHManyToManyBW DijkstraManyToManyBW;

static void readNodes(const MyGraph *const g, const google::protobuf::RepeatedField<uint32_t>& in, vector<NodeID>* vs) {
    vs->resize(in.size());
    for (auto i=0; i < in.size(); i++) {
        // map the actual node ID to the node ID that is used internally by
        // highway-node routing
        (*vs)[i] = g->mapExtToIntNodeID(in.Get(i));
    }
}

static void readMessage(const string& filename, ::google::protobuf::Message* message) {
    int fd = open(filename.c_str(), O_RDONLY);
    message->ParseFromFileDescriptor(fd);
    close(fd);
}

static void writeMessage(const string& filename, ::google::protobuf::Message* message) {
    int fd = open(filename.c_str(), O_WRONLY);
    message->SerializeToFileDescriptor(fd);
    close(fd);
}

MyGraph *const readGraph(const string& inGraphFn) {
    ifstream inGraph(inGraphFn, ios::binary);
    if (!inGraph) {
        cerr << "Input file '" << inGraphFn << "' couldn't be opened." << endl;
        exit(-1);
    }
    return new MyGraph(inGraph);
}

class many2many_opts : public options_parser_t {
public:

    string input_fn;
    string nodes_fn;
    string matrix_fn;

    many2many_opts() : options_parser_t("i:n:o:") {}
    void parse() override {
        // input params
        set_option('i', &input_fn);
        set_option('n', &nodes_fn);

        // output params
        set_option('o', &matrix_fn);
    }
    void
    run_checks() override {
        check_not_empty(input_fn, "input file not provided (-i)");
        check_not_empty(nodes_fn, "input nodes not provided (-n)");
        check_not_empty(matrix_fn, "output matrix file name not provided (-o)");
    }

    void
    usage(const char *p) override {
        cerr << "usage: " << p << " -i graph.sgr -n many2many-nodes.pb -o matrix.pb" << endl;
    }
};


/** The main program. */
int main(int argc, char *argv[])
{
    GOOGLE_PROTOBUF_VERIFY_VERSION;

    many2many_opts opts;
    opts.parse_options(argc, argv);


    const bool writeSourceNodes = false;
    const bool writeMatrix = false;
    const bool validateResult = false;

    VERBOSE( cerr << "read graph from '" << opts.input_fn << "'" << endl);
    VERBOSE( cerr << "read nodes from '" << opts.nodes_fn << "'" << endl);
    VERBOSE( cerr << "write matrix to '" << opts.matrix_fn << "'" << endl);
    VERBOSE( if (writeSourceNodes) cerr << "write source nodes to cerr" << endl );
    VERBOSE( if (writeMatrix)      cerr << "write matrix to cerr" << endl );
    VERBOSE( if (validateResult)   cerr << "validate result" << endl );

    // read graph
    MyGraph *const graph = readGraph(opts.input_fn);

    VERBOSE( cerr << "Graph read." << endl );

    // read the input
    protos::ch::Nodes nodes;
    readMessage(opts.nodes_fn, &nodes);

    vector<NodeID> sources;
    readNodes(graph, nodes.sources(), &sources);

    vector<NodeID> targets;
    readNodes(graph, nodes.targets(), &targets);

    // create many-to-many object
    LevelID earlyStopLevel = 10;
    ManyToMany<MyGraph, DijkstraManyToManyFW, DijkstraManyToManyBW, performBucketScans> mtm(graph, earlyStopLevel);

    // compute matrix
    Matrix<EdgeWeight> matrix(sources.size(), targets.size());
    mtm.computeMatrix( sources, targets, matrix );

    // generate and write output
    protos::ch::Matrix outMatrix;
    outMatrix.set_scale(SCALE);
    for (int i=0; i < sources.size(); i++) {
        protos::ch::Row* row = outMatrix.add_rows();
        for (int j=0; j < targets.size(); j++) {
            row->add_data(SCALE * matrix.value(i, j));
        }
    }
    writeMessage(opts.matrix_fn, &outMatrix);

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
