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

#ifndef _COMMAND_CONSTRUCT_H
#define _COMMAND_CONSTRUCT_H

#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>

#include "../Command.h"
#include "../datastr/graph/UpdateableGraph.h"
#include "../datastr/graph/SearchGraph.h"
#include "../io/createGraph.h"
#include "../io/output.h"
#include "../processing/DijkstraCH.h"
#include "../processing/ConstructCH.h"

const static bool COMPRESSED_SGR_FILE = false;

/**
 * Command-class for hierarchy construction with a given node order.
 * The main-procedure of this class is executed in ../main.cpp.
 */
namespace command {
    class Construct : public Command {
        public:
        typedef processing::ConstructCH<datastr::graph::UpdateableGraph> ProcessingConstruct;

        /**
         * The main-method is called from ../main.cpp if the first command-line
         * argument is "-c".
         */
        int main(int argc, char *argv[]) {

            // Parameter variables, specifiy input/output files
            // and contraction parameters.
            string ddsgFile("");      // input graph file
            string hcnFile("");      // input level file
            string importFile("");
            string exportFile("");

            // Encapsulates the contraction parameters, e.g. hop limits.
            ProcessingConstruct::ContractParameters contractParams;

            bool testAllNodes = false;
            string retrieve("");
            string sgrFile("");
            string edgesGraphFile("");
            string demandsFile("");
            string chFile("");
            datastr::graph::SearchGraphNodeOrder nodeOrder = datastr::graph::SGNO_LOGLEVEL_ORIGINAL;
            int i;

            int opt2 = -1;
            while ((opt2 = getopt(argc, argv, "f:h:m:k:ts:Lz:O:G:d:c:C:")) != -1)
            {
                switch(opt2)
                {
                    // Input graph in ddsg-format
                    case 'f':
                        ddsgFile = string(optarg);
                        break;
                    // Input node order as hcn-file. This file is usually created
                    // during the node ordering. Instead of a hcn-file, a already
                    // constructed search graph can be specified with -z
                    case 'h':
                        hcnFile = string(optarg);
                        break;
                    // Search graph in sgr-format. Note, this binary format
                    // depends on switches in config.h since it just serializes the
                    // graph representation.
                    // It is an input parameter if -h is not specified, otherwise an
                    // output parameter.
                    case 'z':
                        sgrFile = string(optarg);
                        break;

                    // limit of settled nodes in local searches during contraction
                    case 'm':
                        contractParams.maxSettledElim = atoi(optarg);
                        break;

                    // hop-/degree-limits specified as comma-separated list,
                    // e.g. "1,3.3,2,10,3,10,5" means hop-limit 1 until degree 3.3,
                    //                          then hop-limit 2 until degree 10,
                    //                          then hop-limit 3 until degree 10,
                    //                          then hop-limit 5
                    case 'k':
                        contractParams.maxHops.clear();
                        Command::createVector(string(optarg),contractParams.maxHops,(double)0);
                        break;

                    // perform forward- and backward searches from all nodes without and
                    // abort-on-success criterion. The statistics concerning settled nodes,...
                    // are stored to file to calculate worst case upper bounds.
                    case 't':
                        testAllNodes = true;
                        break;

                    // Retrieve witnesses and shortcuts from file. Those files contain The
                    // stored witnesses and shortcuts during node ordering. Note that
                    case 's':
                        retrieve = string(optarg);
                        break;

                    // omit local edge reduction
                    case 'L':
                        contractParams.localReduceEdges = !contractParams.localReduceEdges;
                        break;

                    // specifiy node order in the search graph
                    //   0 = node order from original graph
                    //   1 = node order according to level (0..(n-1))
                    //   2 = hybrid node order, partition nodes into partitons
                    //       of size n/2,n/4,n/8,... and use the original node order
                    //       within the partiton, see ../datastr/graph/SearchGraph.h
                    case 'O':
                        i = atoi(optarg);
                        if ( i == 0 ) nodeOrder = datastr::graph::SGNO_ORIGINAL;
                        else if ( i == 1 ) nodeOrder = datastr::graph::SGNO_LEVEL;
                        else if ( i == 2 ) nodeOrder = datastr::graph::SGNO_LOGLEVEL_ORIGINAL;
                        break;

                    // export edges of search graph to text file, requested by D.Delling
                    case 'G':
                        edgesGraphFile = string(optarg);
                        break;
                    // demands file, each line of the text-file specifies a s-t-query
                    // the runtime of each query is expored. Can be used to perform
                    // local queries.
                    case 'd':
                        demandsFile = string(optarg);
                        break;

                    // Core size, the contraction aborts if the remaining graph
                    // only contains this number of nodes. Used for assymetric
                    // many-to-many calculations.
                    case 'c':
                        contractParams.coreSize = atoi(optarg);
                        break;

                    // Export contraction hierarchy including node order and shortcuts
                    // to file
                    case 'C':
                        chFile = string(optarg);
                        break;
                }
            }


            // CH-export needs middle nodes of shortcuts, thus USE_CH_EXPAND
            #ifndef USE_CH_EXPAND
            if ( chFile != "" )
            {
                cerr << "#define USE_CH_EXPAND in config.h to export CH-file (-C)." << endl;
                exit(1);
            }
            #endif


            datastr::graph::SearchGraph* searchGraph = NULL;

            // Output messages are written to this stringstream. It will be written to
            // stdout and a logfile after the contraction and tests.
            stringstream ss;

            // number of s-t-pairs calculated for average query time and path expansion
            static const unsigned int noOfTestCases = 100000;
            // number of s-t-pairs whose shortest paths distance in the CH is compared to the
            // shortest paths distance in the original graph, see TEST below.
            static const unsigned int noOfTestCasesVerify = 100;
            stringstream ssName;

            ss << "noOfTestCases " << noOfTestCases << endl;
            ss << "noOfTestCasesVerify " << noOfTestCasesVerify << endl;
            ss << "sizeof Edge " << sizeof(Edge) << endl;

            // currently no on the fly edge reduction with witnesses
            if ( retrieve != "" ) contractParams.localReduceEdges = false;

            stPairs runs;       // s-t-pairs for tests

            // TEST
            // Computate some shortest paths distances in the original graph without using CHs.
            // Those distances can be compared to the results of the CH. Since the computation
            // in the original graph is time consuming, the distances are stored. The source and
            // target nodes are specified by the initalization of the random number generator
            // using srand(1);.
            // Perform test runs at the beginning because otherwise there is not enough free space.

            // The shortest paths distances are stored in this file
            string testLengthsFile(ddsgFile + ".test-lengths");
            ifstream inTest(testLengthsFile.c_str());
            vector<EdgeWeight> _test;
            if ( inTest.is_open() )
            {
                VERBOSE( cout << "Read test lengths from " << testLengthsFile << " ..." << endl; )
                while ( !inTest.eof() )
                {
                    EdgeWeight w;
                    inTest >> w >> ws;
                    _test.push_back(w);
                }
                inTest.close();
            }

            if ( _test.size() < noOfTestCasesVerify )
            {
                VERBOSE( cout << "TEST" << endl; )
                ifstream in(ddsgFile.c_str());
                if (!in.is_open()) { cerr << "Cannot open " << ddsgFile << endl; exit(1); }
                datastr::graph::UpdateableGraph* tGraph = importGraphListOfEdgesUpdateable(in, false, false, "");
                in.close();

                // Initalize random number generator to get same source and target nodes as
                // for the runs with the CH.
                srand(1);
                for (NodeID x = 0; x < noOfTestCases; x++) {
                    NodeID s = randomNodeID(tGraph->noOfNodes());
                    NodeID t = randomNodeID(tGraph->noOfNodes());
                    runs.push_back( stPair(s, t) );
                }


                VERBOSE( Percent percent(runs.size()); )
                processing::DijkstraCH<datastr::graph::UpdateableGraph, NormalPQueue, 2, false> dijkstraTest(tGraph);

                _test.resize(noOfTestCasesVerify);
                for (NodeID x = 0; x < noOfTestCasesVerify; x++) {
                    dijkstraTest.bidirSearch(runs[x].first, runs[x].second);
                    Path path;
                    dijkstraTest.pathTo(path, runs[x].second, -1);
                    _test[x] = path.length();
                    cout << "[" << x << "] " << runs[x].first << " -> " << runs[x].second << " length " << path.length() << endl;
                    //cout << path << endl;
                    dijkstraTest.clear();
                    VERBOSE( percent.printStatus(x); )
                }
                delete tGraph;

                VERBOSE( cout << "Write test lengths to " << testLengthsFile << " ..." << endl; )
                ofstream outTest(testLengthsFile.c_str());
                if (!outTest.is_open()) { cerr << "Cannot write to " << testLengthsFile << endl; exit(1); }
                for (NodeID x = 0; x < noOfTestCasesVerify; x++) {
                    outTest << _test[x] << endl;
                }
                outTest.close();
            }
            // TEST END


            // Construct hierarchy with given node order from a hcn-file.
            if ( hcnFile != "" )
            {
                // create a file name depending on the parameters
                ssName << hcnFile << "-" << contractParams.maxSettledElim << "-";
                copy(contractParams.maxHops.begin(), contractParams.maxHops.end(), ostream_iterator<double>(ssName,"-"));
                ssName << contractParams.localReduceEdges << "-" << nodeOrder << "-" << contractParams.coreSize;
                datastr::graph::UpdateableGraph* updGraph;

                // input original graph
                ifstream in(ddsgFile.c_str());
                if (!in.is_open()) { cerr << "Cannot open " << ddsgFile << endl; exit(1); }
                updGraph = importGraphListOfEdgesUpdateable(in, false, false, "");
                NodeID noOfEdgesOriginal = updGraph->noOfExistingEdges();
                in.close();

                // count unidirectional edges in original graph, needed for
                // the calculation of the space overhead in B/node.
                NodeID noOfUnidirectionalEdges = 0;
                for ( NodeID u = 0; u < updGraph->noOfNodes(); u++ )
                {
                    for ( EdgeID e = updGraph->firstLevelEdge(u); e < updGraph->lastEdge(u); e++ )
                    {
                        const Edge& edge = updGraph->edge(e);
                        if ( edge.isDirected(0) )
                        {
                            noOfUnidirectionalEdges++;
                        }
                    }
                }

                double time;
                // Construct with witness and shortcut infos. This is only and
                // preliminary implementation for testing purpose.
                if ( retrieve != "" )
                {
                    VERBOSE( cout << "Construct from " << hcnFile << " ..." << endl; )
                    ProcessingConstruct construct(updGraph);
                    in.open(hcnFile.c_str());
                    if (!in.is_open()) { cerr << "Cannot open " << hcnFile << endl; exit(1); }
                    construct.readLevels(in);
                    in.close();

                    ifstream inShortcuts((retrieve+".shortcuts").c_str());
                    if (!inShortcuts.is_open()) { cerr << "Cannot open " << (retrieve+".shortcuts") << endl; exit(1); }
                    VERBOSE( cout << "Retrieve shortcuts from " << (retrieve+".shortcuts") << endl; )
                    ifstream inWitnesses((retrieve+".witnesses").c_str());
                    if (!inWitnesses.is_open()) { cerr << "Cannot open " << (retrieve+".witnesses") << endl; exit(1); }
                    VERBOSE( cout << "Retrieve witnesses from " << (retrieve+".witnesses") << endl; )

                    time = timestamp();
                    construct.constructHierarchyWithWitnesses(contractParams,inShortcuts, inWitnesses);
                    time = timestamp() - time;
                    construct.clear();
                }

                // Construct CH only using the node order given by the hcn-file.
                else
                {
                    VERBOSE( cout << "Construct from " << hcnFile << " ..." << endl; )

                    // Construction class, contains code for construction.
                    ProcessingConstruct construct(updGraph);

                    // read in node order from hcn-file
                    in.open(hcnFile.c_str());
                    if (!in.is_open()) { cerr << "Cannot open " << hcnFile << endl; exit(1); }
                    construct.readLevels(in);
                    in.close();

                    // *** Perform construction ***
                    time = timestamp();
                    construct.constructHierarchy(contractParams);
                    time = timestamp() - time;

                    construct.clear();
                }

                // Convert UpdateableGraph to SearchGraph, a search-graph allows
                // faster queries and lower space consumption. It also changes the node
                // order that additionaly improves the avgerage query time.
                VERBOSE( cout << "Construct search-graph..." << endl; )
                searchGraph = new datastr::graph::SearchGraph(updGraph, nodeOrder);

                // Save search-graph to sgr-file. This is a binary format that depends
                // on the graph representation in the main memory. So it depends
                // on certain switches in ../config.h that specify the edge representation.
                if ( sgrFile != "" )
                {
                    if (sgrFile == "x")
                    {
                        sgrFile = ssName.str()+".sgr";
                        if ( COMPRESSED_SGR_FILE ) sgrFile.append(".gz");
                    }
                    VERBOSE( cout << "Serialize search-graph to " << sgrFile << " ..." << endl; )
                    ofstream out(sgrFile.c_str());
                    if (!out.is_open()) { cerr << "Cannot write to " << sgrFile << endl; exit(1); }
                    boost::iostreams::filtering_ostream zOut;
                    if ( COMPRESSED_SGR_FILE) zOut.push(boost::iostreams::gzip_compressor());
                    zOut.push(out);
                    searchGraph->serialize(zOut);
                }

                // save edges of search-graph to text file, requested by D.Delling
                if ( edgesGraphFile != "")
                {
                    if ( edgesGraphFile == "x" )
                    {
                        edgesGraphFile = ssName.str() + ".edges";
                    }
                    VERBOSE( cout << "Write edges of search-graph to " << edgesGraphFile << " ..." << endl; )
                    bool b = false;
                    CH_EXPAND( b = true; )
                    if ( !b ) cout << "Warning: no shortcut information available!" << endl;

                    vector<NodeID> mapIntToExtNodeID(searchGraph->noOfNodes());
                    for ( NodeID u = 0; u < searchGraph->noOfNodes(); u++ )
                    {
                        mapIntToExtNodeID[searchGraph->mapExtToIntNodeID(u)] = u;
                    }
                    ofstream out(edgesGraphFile.c_str());
                    if (!out.is_open()) { cerr << "Cannot write to " << edgesGraphFile << endl; exit(1); }
                    for ( NodeID u = 0; u < searchGraph->noOfNodes(); u++ )
                    {
                        EdgeID lastEdge = searchGraph->lastEdge(searchGraph->mapExtToIntNodeID(u));
                        for ( EdgeID e = searchGraph->firstLevelEdge(searchGraph->mapExtToIntNodeID(u)); e < lastEdge; e++ )
                        {
                            const Edge& edge = searchGraph->edge(e);
                            int dir = 3;
                            if (edge.isDirected(0)) dir -= 2;
                            if (edge.isDirected(1)) dir -= 1;
                            out << u << " " << mapIntToExtNodeID[edge.target()];
                            out << " " << edge.weight();
                            out << " " << dir;
                            if (edge.isShortcut()) out << " " << mapIntToExtNodeID[edge.shortcutMiddle()];
                            out << endl;
                        }
                    }
                    out.close();
                    out.open((edgesGraphFile+".noweight").c_str());
                    for ( NodeID u = 0; u < searchGraph->noOfNodes(); u++ )
                    {
                        EdgeID lastEdge = searchGraph->lastEdge(searchGraph->mapExtToIntNodeID(u));
                        for ( EdgeID e = searchGraph->firstLevelEdge(searchGraph->mapExtToIntNodeID(u)); e < lastEdge; e++ )
                        {
                            const Edge& edge = searchGraph->edge(e);
                            int dir = 3;
                            if (edge.isDirected(0)) dir -= 2;
                            if (edge.isDirected(1)) dir -= 1;
                            out << u << " " << mapIntToExtNodeID[edge.target()];
                            out << " " << dir;
                            if (edge.isShortcut()) out << " " << mapIntToExtNodeID[edge.shortcutMiddle()];
                            out << endl;
                        }
                    }
                    out.close();

                }

                // save contraction hierarchy to file
                // a special file format (.ch) is used, see ../docu/chDocu.html
                if ( chFile != "" )
                {
                    if (chFile == "x")
                    {
                        chFile = ssName.str()+".ch";
                    }
                    VERBOSE( cout << "Serialize contraction hierarchy to " << chFile << " ..." << endl; )
                    ofstream out(chFile.c_str());
                    if (!out.is_open()) { cerr << "Cannot write to " << chFile << endl; exit(1); }
                    exportContractionHierarchy( out, updGraph );
                }


                // calculate statistics: number of shortcut edges
                // this information is saved in the log file
                NodeID noOfShortcutEdges = 0;
                for ( NodeID u = 0; u < searchGraph->noOfNodes(); u++ )
                {
                    for ( EdgeID e = searchGraph->firstLevelEdge(u); e < searchGraph->lastEdge(u); e++ )
                    {
                        const Edge& edge = searchGraph->edge(e);
                        if ( edge.isShortcut() )
                        {
                            noOfShortcutEdges++;
                        }
                    }
                }


                // Write parameters and statistics to string-stream, it will be written
                // to a logfile and to stdout. This avois outputting the same information twice.
                ss << "maxSettledElim " << contractParams.maxSettledElim << endl;
                ss << "maxHops "; copy(contractParams.maxHops.begin(), contractParams.maxHops.end(), ostream_iterator<double>(ss,"-")); ss << endl;
                ss << "localReduceEdges " << contractParams.localReduceEdges << endl;
                ss << "construction time: " <<  time << " s" << endl;
                ss << "nodeOrder: " << nodeOrder << endl;
                ss << "#edges: " << searchGraph->noOfEdges() << " / " << updGraph->noOfExistingEdges() << endl;
                ss << "#edges original unidir./bidir.: " << noOfUnidirectionalEdges << " / " << noOfEdgesOriginal << endl;
                ss << "space overhead [B/node] unidir./bidir.: " <<
                        ((((double)sizeof(Edge)*searchGraph->noOfEdges())
                        -((double)8*noOfUnidirectionalEdges))/searchGraph->noOfNodes()) << " / " <<
                    ((((double)sizeof(Edge)*searchGraph->noOfEdges())
                        -((double)8*noOfEdgesOriginal))/searchGraph->noOfNodes()) << endl;
                ss << "#shortcut edges: " << noOfShortcutEdges << endl;

                delete updGraph;
            }

            // if no node order is given by a hcn-file, the no hierarchy construction is performed
            // and the search-graph is deserialized from a binary file
            // Note: only sgr-files created by the same binary should be used because
            //       the edge datastructure may differ
            else
            {
                VERBOSE( cout << "Deserialize search-graph from " << sgrFile << " ..." << endl; )

                ifstream in(sgrFile.c_str());
                if (!in.is_open()) { cerr << "Cannot open " << sgrFile << endl; exit(1); }
                boost::iostreams::filtering_istream zIn;
                if ( COMPRESSED_SGR_FILE) zIn.push(boost::iostreams::gzip_decompressor());
                zIn.push(in);
                searchGraph = new datastr::graph::SearchGraph(zIn);

                ssName << sgrFile;
                ss << "#edges: " << searchGraph->noOfEdges() << endl;
            }
            assert( searchGraph != NULL );

            // Prepare source and target pairs for statistic runs. They are specified by the
            // initalization of the random number generator: srand(1)
            if ( runs.size() == 0 )
            {
                srand(1);
                for (NodeID x = 0; x < noOfTestCases; x++) {
                    NodeID s = randomNodeID(searchGraph->noOfNodes());
                    NodeID t = randomNodeID(searchGraph->noOfNodes());
                    runs.push_back( stPair(s, t) );
                }
            }

            // Prepare source and target pairs for cache warmup runs. They are specified by the
            // initalization of the random number generator: srand(12)
            // Different pairs as for the statistic runs are used to measure query times
            // on a server with permanently arriving, different requests
            srand(12);
            stPairs warmup;
            for (NodeID x = 0; x < noOfTestCases; x++) {
                NodeID s = randomNodeID(searchGraph->noOfNodes());
                NodeID t = randomNodeID(searchGraph->noOfNodes());
                warmup.push_back( stPair(s, t) );
            }

            // node mapping (external to internal IDs)
            // note: in UpdateableGraph the original IDs are used,
            //       mapping only necessary for SearchGraph
            // the node mapping is necessary to get the same source-target pairs used during
            // the comparison test runs, see above (TEST)
            for (NodeID x = 0; x < runs.size(); x++) {
                assert( runs[x].first < searchGraph->noOfNodes() );
                assert( runs[x].second < searchGraph->noOfNodes() );
                runs[x].first = searchGraph->mapExtToIntNodeID(runs[x].first);
                runs[x].second = searchGraph->mapExtToIntNodeID(runs[x].second);


                assert( warmup[x].first < searchGraph->noOfNodes() );
                assert( warmup[x].second < searchGraph->noOfNodes() );
                warmup[x].first = searchGraph->mapExtToIntNodeID(warmup[x].first);
                warmup[x].second = searchGraph->mapExtToIntNodeID(warmup[x].second);
            }

            DijkstraSearchCH dijkstra(searchGraph);

            double singleTime;

            // Cache warmup, the required run time does not count
            Checksum checkSum1 = 0;
            VERBOSE( cout << "Cache warmup" << flush; )
            singleTime = timestamp();
            for ( NodeID x = 0; x < warmup.size(); x++ )
            {
                EdgeWeight dist = dijkstra.bidirSearch(warmup[x].first, warmup[x].second);
                if (dist == Weight::MAX_VALUE) dist = 0;
                checkSum1 += dist;
                dijkstra.clear();
            }
            singleTime = timestamp() - singleTime;
            VERBOSE( cout << " done, took " << singleTime << " seconds, checksum " << checkSum1 << "." << endl; )


            // Test runs to measure the average runtime, only the sum of all shortest paths
            // distances is calculated as a checksum. Additional statistics as settled nodes,...
            // are calculated in a second run
            VERBOSE( cout << "Test run without statistics" << flush; )
            double time = timestamp();
            Checksum checkSum2 = 0;
            for (NodeID x = 0; x < runs.size(); x++) {
                EdgeWeight dist = dijkstra.bidirSearch(runs[x].first, runs[x].second);
                if (dist == Weight::MAX_VALUE) dist = 0;
                checkSum2 += dist;
                dijkstra.clear();
            }
            time = timestamp() - time;
            VERBOSE( cout << " done, took " << time << " seconds, checksum " << checkSum2 << "." << endl; )

            ss << "avg query time: " <<  ((time/(runs.size()))*1000) << " ms" << endl;


            // If the program is compiled with edges supporting path expansion
            // (switch USE_CH_EXPAND in ../config.h) we perform test runs with
            // path expansion. The same source and target pairs as for the average query time
            // are used.
            bool expandPaths = false;
            CH_EXPAND( expandPaths = true; )

            double timeExpand = 0;
            Checksum checkSum4 = 0;
            if ( expandPaths )
            {
                timeExpand = timestamp();
                VERBOSE( cout << "Test run with path expansion " << flush; )
                for (NodeID x = 0; x < runs.size(); x++) {
                    dijkstra.bidirSearch(runs[x].first, runs[x].second);
                    Path path;
                    dijkstra.pathTo(path, runs[x].second, -1, true, true /* expand */);
                    EdgeWeight dist = path.length();
                    if (dist == Weight::MAX_VALUE) dist = 0;
                    checkSum4 += dist;
                    dijkstra.clear();
                }
                timeExpand = timestamp() - timeExpand;
                VERBOSE( cout << " done, took " << timeExpand << " seconds, checksum " << checkSum4 << "." << endl; )

                // now subtract time without expansion from the measured time with expansion
                // to get sole expansion time
                timeExpand -= time;
            }


            // Now we perform test runs with the same soure and target pairs as in the two
            // loops above. Now, the runtime is not relevant, we now gather statistics like
            // the number of settled nodes. In case of path expansion, we also check the
            // expanded paths. And we compare the first noOfTestCasesVerify shortest paths
            // distances to the distances computated without CHs.
            VERBOSE( cout << "Test run with statistics" << flush; )
            double time2 = timestamp();

            Checksum checkSum3 = 0;
            unsigned int settledNodesMin = UINT_MAX;
            unsigned int settledNodesMax = 0;
            unsigned int settledNodesSum = 0;
            unsigned int settledMinusStalledNodesMin = UINT_MAX;
            unsigned int settledMinusStalledNodesMax = 0;
            unsigned int settledMinusStalledNodesSum = 0;

            unsigned int pathNoOfEdgesMin = UINT_MAX;
            unsigned int pathNoOfEdgesMax = 0;
            unsigned int pathNoOfEdgesSum = 0;
            unsigned int pathExpandNoOfEdgesMin = UINT_MAX;
            unsigned int pathExpandNoOfEdgesMax = 0;
            unsigned int pathExpandNoOfEdgesSum = 0;

            COUNTING( counter.reset(); )
            for (NodeID x = 0; x < runs.size(); x++) {
                EdgeWeight dist = dijkstra.bidirSearch(runs[x].first, runs[x].second);
                if (dist == Weight::MAX_VALUE) dist = 0;
                checkSum3 += dist;

                // TEST
                // verify shortest paths distances
                if ( x < noOfTestCasesVerify )
                {
                    Path path;
                    dijkstra.pathTo(path, runs[x].second, -1);
                    if ( path.length() != _test[x] )
                    {
                        cout << endl;
                        cout << "x = " << x << endl;
                        cout << path.length() << " != " << _test[x] << endl;
                        cout << path << endl;
                        exit(1);
                    }
                }
                // TEST END


                unsigned int settledNodes = dijkstra.noOfSettledNodes();
                settledNodesSum += settledNodes;
                if (settledNodes > settledNodesMax) settledNodesMax = settledNodes;
                if (settledNodes < settledNodesMin) settledNodesMin = settledNodes;
                unsigned int settledMinusStalledNodes = dijkstra.noOfSettledMinusStalledNodes();
                settledMinusStalledNodesSum += settledMinusStalledNodes;
                if (settledMinusStalledNodes > settledMinusStalledNodesMax) settledMinusStalledNodesMax = settledMinusStalledNodes;
                if (settledMinusStalledNodes < settledMinusStalledNodesMin) settledMinusStalledNodesMin = settledMinusStalledNodes;


                // verify expanded paths
                if ( expandPaths )
                {

                    Path path;
                    //singleTime = timestamp();
                    dijkstra.pathTo(path, runs[x].second, -1, true, true /* expand */);
                    //timeExpand += timestamp() - singleTime;

                    pathExpandNoOfEdgesSum += path.noOfEdges();
                    if (path.noOfEdges() < pathExpandNoOfEdgesMin) pathExpandNoOfEdgesMin = path.noOfEdges();
                    if (path.noOfEdges() > pathExpandNoOfEdgesMax) pathExpandNoOfEdgesMax = path.noOfEdges();


                    // check expanded path to not expanded path
                    Path path2;
                    dijkstra.pathTo(path2, runs[x].second, -1, true, false /* do not expand */);

                    EdgeWeight w = 0;
                    for ( EdgeID index = 0; index < path.noOfEdges(); index++ )
                    {
                        const Edge& edge = searchGraph->edge(path.edge(index));
                        w += edge.weight();
                        assert( !edge.isShortcut() );
                    }
                    assert( path.isNotConnected() || w == path.length() );
                    assert( path2.isNotConnected() || path2.length() == path.length() );
                    assert( path2.isNotConnected() || path2.firstNode() == path.firstNode() );
                    assert( path2.isNotConnected() || path2.lastNode() == path.lastNode() );

                    // some stupid stuff so "w" is not removed by the optimizer
                    if ( !(path.isNotConnected() || w == path.length()) || path.length() != path2.length() )
                    {
                        cout << "err" << endl;
                        exit(1);
                    }
                    pathNoOfEdgesSum += path2.noOfEdges();
                    if (path2.noOfEdges() < pathNoOfEdgesMin) pathNoOfEdgesMin = path2.noOfEdges();
                    if (path2.noOfEdges() > pathNoOfEdgesMax) pathNoOfEdgesMax = path2.noOfEdges();
                }
                dijkstra.clear();
            }
            time2 = timestamp() - time2;
            VERBOSE( cout << " done, took " << time2 << " seconds, checksum " << checkSum3 << "." << endl; )
            ss << "#settled nodes min/avg/max: " << settledNodesMin << " / " << (((double)settledNodesSum/runs.size())) << " / " << settledNodesMax << endl;
            ss << "#settled minus stalled nodes min/avg/max: " << settledMinusStalledNodesMin << " / " << (((double)settledMinusStalledNodesSum/runs.size())) << " / " << settledMinusStalledNodesMax << endl;
            COUNTING( ss << "#relaxed edges avg: " << (counter.count(COUNT_RELAXED_EDGES,COUNT_AKKU)/runs.size()) << " / #successfully relaxed edges avg: " <<  (counter.count(COUNT_RELAXED_EDGES_SUCCESS,COUNT_AKKU)/runs.size()) << endl; )
            if ( expandPaths )
            {
                ss << "avg path expansion time: " << ((timeExpand/(runs.size()))*1000) << " ms" << endl;
                ss << "#edges in expanded path min/avg/max: " << pathExpandNoOfEdgesMin << " / " << ((double)pathExpandNoOfEdgesSum/(runs.size())) << " / " << pathExpandNoOfEdgesMax << endl;
                ss << "#edges in searched path min/avg/max: " << pathNoOfEdgesMin << " / " << ((double)pathNoOfEdgesSum/(runs.size())) << " / " << pathNoOfEdgesMax << endl;
            }

            // Option: -t
            // to calculate upper bounds on the search space size for arbitrary queries
            // we perform separated forward and backward searches from each node and
            // store the search space sizes. This is possible since the query times
            // for CHs are very small. The results are stored in separate files.
            // An external R-script (maxsettled*.r) can use them to create a plot.
            if (testAllNodes)
            {
                cout << "Test search spaces of all nodes" << endl;
                double timeTestAll = 0;
                vector<NodeID> settledNodes[2];
                vector<NodeID> settledMinusStalledNodes[2];
                COUNTING( vector<NodeID> relaxedEdges[2]; )
                COUNTING( vector<NodeID> successfullyRelaxedEdges[2]; )
                VERBOSE( Percent percent( searchGraph->noOfNodes() ) );
                for ( NodeID u = 0; u < searchGraph->noOfNodes(); u++ )
                {
                    for ( NodeID j = 0; j < 2; j++ )
                    {
                        COUNTING( counter.reset(); )
                        singleTime = timestamp();
                        dijkstra.searchWithoutTarget(u, j);
                        timeTestAll += timestamp() - singleTime;

                        NodeID i = dijkstra.noOfSettledNodes();
                        if (settledNodes[j].size() <= i) settledNodes[j].resize(i+1,0);
                        settledNodes[j][i]++;

                        i = dijkstra.noOfSettledMinusStalledNodes();
                        if (settledMinusStalledNodes[j].size() <= i) settledMinusStalledNodes[j].resize(i+1,0);
                        settledMinusStalledNodes[j][i]++;

                        COUNTING( i  = (NodeID)counter.count(COUNT_RELAXED_EDGES,COUNT_AKKU); )
                        COUNTING( if (relaxedEdges[j].size() <= i) relaxedEdges[j].resize(i+1,0); )
                        COUNTING( relaxedEdges[j][i]++; )

                        COUNTING( i  = (NodeID)counter.count(COUNT_RELAXED_EDGES_SUCCESS,COUNT_AKKU); )
                        COUNTING( if (successfullyRelaxedEdges[j].size() <= i) successfullyRelaxedEdges[j].resize(i+1,0); )
                        COUNTING( successfullyRelaxedEdges[j][i]++; )

                        singleTime = timestamp();
                        dijkstra.clear();
                        timeTestAll += timestamp() - singleTime;

                    }
                    VERBOSE( percent.printStatus(u) );
                }

                ss << "time to test all: " <<  timeTestAll << " s" << endl;

                // settled nodes
                ofstream out((ssName.str()+".stats.settled-nodes-forward").c_str());
                if (!out.is_open()) { cerr << "Cannot write to " << (ssName.str()+".stats.settled-nodes-forward") << endl; exit(1); }
                for ( NodeID i = 0; i < settledNodes[0].size(); i++ )
                {
                    out << settledNodes[0][i] << endl;
                }
                out.close();
                out.open((ssName.str()+".stats.settled-nodes-backward").c_str());
                if (!out.is_open()) { cerr << "Cannot write to " << (ssName.str()+".stats.settled-nodes-backward") << endl; exit(1); }
                for ( NodeID i = 0; i < settledNodes[1].size(); i++ )
                {
                    out << settledNodes[1][i] << endl;
                }
                out.close();

                // settled minus stalled nodes
                out.open((ssName.str()+".stats.settled-minus-stalled-nodes-forward").c_str());
                if (!out.is_open()) { cerr << "Cannot write to " << (ssName.str()+".stats.settled-minus-stalled-nodes-forward") << endl; exit(1); }
                for ( NodeID i = 0; i < settledMinusStalledNodes[0].size(); i++ )
                {
                    out << settledMinusStalledNodes[0][i] << endl;
                }
                out.close();
                out.open((ssName.str()+".stats.settled-minus-stalled-nodes-backward").c_str());
                if (!out.is_open()) { cerr << "Cannot write to " << (ssName.str()+".settled-minus-stalled-nodes-backward") << endl; exit(1); }
                for ( NodeID i = 0; i < settledMinusStalledNodes[1].size(); i++ )
                {
                    out << settledMinusStalledNodes[1][i] << endl;
                }
                out.close();

                // relaxed edges
                COUNTING( out.open((ssName.str()+".stats.relaxed-edges-forward").c_str());  )
                COUNTING( if (!out.is_open()) { cerr << "Cannot write to " << (ssName.str()+".stats.relaxed-edges-forward") << endl; exit(1); } )
                COUNTING( for ( NodeID i = 0; i < relaxedEdges[0].size(); i++ )             )
                COUNTING( {                                                                 )
                COUNTING(     out << relaxedEdges[0][i] << endl;                            )
                COUNTING( }                                                                 )
                COUNTING( out.close();                                                      )
                COUNTING( out.open((ssName.str()+".stats.relaxed-edges-backward").c_str()); )
                COUNTING( if (!out.is_open()) { cerr << "Cannot write to " << (ssName.str()+".stats.relaxed-edges-backward") << endl; exit(1); } )
                COUNTING( for ( NodeID i = 0; i < relaxedEdges[1].size(); i++ )             )
                COUNTING( {                                                                 )
                COUNTING(     out << relaxedEdges[1][i] << endl;                            )
                COUNTING( }                                                                 )
                COUNTING( out.close();                                                      )

                // successfully relaxed edges
                COUNTING( out.open((ssName.str()+".stats.successfully-relaxed-edges-forward").c_str());  )
                COUNTING( if (!out.is_open()) { cerr << "Cannot write to " << (ssName.str()+".stats.successfully-relaxed-edges-forward") << endl; exit(1); } )
                COUNTING( for ( NodeID i = 0; i < successfullyRelaxedEdges[0].size(); i++ )              )
                COUNTING( {                                                                              )
                COUNTING(     out << successfullyRelaxedEdges[0][i] << endl;                             )
                COUNTING( }                                                                              )
                COUNTING( out.close();                                                                   )
                COUNTING( out.open((ssName.str()+".stats.successfully-relaxed-edges-backward").c_str()); )
                COUNTING( if (!out.is_open()) { cerr << "Cannot write to " << (ssName.str()+".stats.successfully-relaxed-edges-backward") << endl; exit(1); } )
                COUNTING( for ( NodeID i = 0; i < successfullyRelaxedEdges[1].size(); i++ )              )
                COUNTING( {                                                                              )
                COUNTING(     out << successfullyRelaxedEdges[1][i] << endl;                             )
                COUNTING( }                                                                              )
                COUNTING( out.close();                                                                   )



                ss << "#total max settled nodes forward/backward: " << (settledNodes[0].size()-1) << " / " << (settledNodes[1].size()-1) << endl;
                ss << "#total max settled minus stalled nodes forward/backward: " << (settledMinusStalledNodes[0].size()-1) << " / " << (settledMinusStalledNodes[1].size()-1) << endl;
                COUNTING( ss << "#total max relaxed edges forward/backward: " << (relaxedEdges[0].size()-1) << " / " << (relaxedEdges[1].size()-1) << endl; )
                COUNTING( ss << "#total max successfully relaxed edges forward/backward: " << (successfullyRelaxedEdges[0].size()-1) << " / " << (successfullyRelaxedEdges[1].size()-1) << endl; )
            }

            // Option: -d
            // Perform queries for source and target pairs specified in a demands file.
            // This allows for local queries that better reflect real scenarios. The demands file
            // for the local queries need to be precomputed by a different program.
            // Note: It is wise to perform three runs for each demands file and to calculate
            //       the median-of-three since the time is measured for each run individually
            //       and almost certainly some outliers occur.
            if ( demandsFile != "" )
            {
                VERBOSE( cout << "Read demands from " << demandsFile << endl; )
                VERBOSE( cout << "Write demands results to " << (ssName.str()+"."+demandsFile+".results") << endl; )
                ifstream in ( demandsFile.c_str() );
                if (!in.is_open()) { cerr << "Cannot open " << demandsFile << endl; exit(1); }
                ofstream out ( (ssName.str()+"."+demandsFile+".results").c_str() );
                if (!out.is_open()) { cerr << "Cannot write to " << (ssName.str()+"."+demandsFile+".results") << endl; exit(1); }
                out << fixed << setprecision(7);
                while ( !in.eof() )
                {
                    NodeID i, s, t, category;
                    EdgeWeight dist;
                    in >> i >> s >> t >> dist >> category >> ws;
                    NodeID s2 = searchGraph->mapExtToIntNodeID(s);
                    NodeID t2 = searchGraph->mapExtToIntNodeID(t);

                    double time = timestamp();
                    EdgeWeight dist2 = dijkstra.bidirSearch(s2,t2);
                    dijkstra.clear();
                    time = timestamp() - time;

                    if ( dist != dist2 )
                    {
                        cerr << s << " -> " << t << ": expected distance " << dist << " got " << dist2 << endl;
                        exit(2);
                    }
                    out << i << " " << time << " " << category << endl;
                }
                in.close();
                out.close();
            }


            // output of the checksums (sum of shortest paths distances)
            // of the different runs
            ss << "checksums " << checkSum1 << " " << checkSum2 << " " << checkSum3;
            if ( expandPaths ) ss << " " << checkSum4;
            ss << endl;

            // write the log messages previously written to the stringstream
            // to stdout and to a log file
            VERBOSE( cout << ss.str(); )
            ofstream log((ssName.str()+".log").c_str());
            if (!log.is_open()) { cerr << "Cannot write to " << (ssName.str()+".log") << endl; exit(1); }
            log << ss.str();
            log.close();

            delete searchGraph;

            return 0;

        }

    };
}
#endif // _COMMAND_CONSTRUCT_H
