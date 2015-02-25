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

#ifndef _COMMAND_NODEORDER_H
#define _COMMAND_NODEORDER_H

#include <fstream>
#include <sstream>

//#include "../datastr/general.h"
#include "../io/createGraph.h"
#include "../io/output.h"

#include "../datastr/graph/SearchGraph.h"
#include "../Command.h"
#include "../processing/ConstructCH.h"

/**
 * Command-class for node order creation.
 * The main-procedure of this class is executed in ../main.cpp.
 */
namespace command {
    class NodeOrder : public Command {
        public:
        typedef processing::ConstructCH<datastr::graph::UpdateableGraph> ProcessingNodeOrder;

        /**
         * The main-method is called from ../main.cpp if the first command-line
         * argument is "-s".
         */
        int main(int argc, char *argv[]) {
            // The second command-line argument specifies the next choice,
            // possible are,
            //   -p: node ordering using a prority queue
            //   -l: partition a hcn-file with n levels into one with only a
            //       few levels
            //   -t: extract the level of a given node or the number of nodes
            //       in each existing level out of a hcn-file
            int opt1 = getopt(argc, argv, "plt");

            /* ************** *
             * priority queue *
             * ************** */
            // perform node ordering and construction using a priority queue
            // this subprogram is reached if the first two command-line arguments
            // are "-s -p". Subsequent arguments specifiy further parameters.
            if (opt1 == 'p') {

                // The variables below contain the parameters for node ordering
                // and input/output. They are specified using command-line
                // arguments. A short documentation of them can be found in
                // ../docu/node_order.html
                string file("graph.ddsg");      // input graph file
                string output("");   // output file (counter of nodes)
                string betweennessFile("");
                string reachFile("");
                string statsFile("");
                string logFile("");
                string exportSearchGraphFile("");

                // stores most of the relevant parameters for node ordering
                // especially the coefficients for the linear combination of
                // the priority terms and the limits to the local searches
                ProcessingNodeOrder::WeightCalculation calc;

                bool test = false;
                string sgrFile("");
                string chFile("");
                int opt2 = -1;
                while ((opt2 = getopt(argc, argv, "f:o:tb:r:s:i:l:x:y:z:w:j:v:n:q:p:m:k:V:S:TKLE:Z:C:e:")) != -1)
                {
                    switch(opt2)
                    {
                        // input graph in ddsg-format
                        case 'f':
                            file = string(optarg);
                            break;
                        // output of node order in hcn-format
                        case 'o':
                            output = string(optarg);
                            break;
                        // perform (time-consuming) tests during node ordering
                        case 't':
                            test = true;
                            break;
                        // betweenness centrality for each node, text format, one line per node
                        case 'b':
                            betweennessFile = string(optarg);
                            break;
                        // reach centrality for each node, text format, one line per node
                        case 'r':
                            reachFile = string(optarg);
                            break;
                        // prefix for statistical files, note that this support needs to be enabled
                        // in the template parameters of processing::ConstructCH
                        case 's':
                            statsFile = string(optarg);
                            break;
                        // log file
                        case 'l':
                            logFile = string(optarg);
                            break;

                        // *** The following coefficients are for the linear ***
                        // *** combination of the priority (weight).         ***
                        // coefficient for edge difference
                        case 'x':
                            calc.edgeDiffMult = atoi(optarg);
                            break;
                        // coefficient for search space size of local searches for contraction
                        case 'y':
                            calc.searchSpaceMult = atof(optarg);
                            break;
                        // coefficient for relative betweenness
                        case 'z':
                            calc.betweennessAdd = atoi(optarg);
                            break;
                        // coefficient for relative reach
                        case 'q':
                            calc.reachAdd = atoi(optarg);
                            break;
                        // coefficient for deleted neighbors
                        case 'w':
                            calc.delNeighbMult = atoi(optarg);
                            break;
                        // coefficient for number of new edges
                        case 'v':
                            calc.newEdgesMult = atoi(optarg);
                            break;
                        // coefficient for Voronoi region size
                        case 'V':
                            calc.voronoiMult = atoi(optarg);
                            break;
                        // coefficient for upper bound on search paths hops
                        case 'S':
                            calc.searchPathHopBorderMult = atoi(optarg);
                            break;
                        // in combination with "-S" does calculate an upper bound on the costs
                        // for a time-dependend search, requested by Veit Batz
                        case 'T':
                            calc.searchPathHopBorderOriginalEdges = true;
                            break;
                        // coefficient for sum of original edges new shortcuts represent
                        case 'e':
                            calc.shortcutOriginalEdgeSumMult = atoi(optarg);
                            break;
                        // limit of settled nodes in local searches during weight calculation
                        case 'n':
                            calc.maxSettledApprox = atoi(optarg);
                            break;
                        // limit of settled nodes in local searches during contraction
                        case 'm':
                            calc.maxSettledElim = atoi(optarg);
                            break;
                        // lazy updates, parameter also specifies the check intervals
                        // for complete recalculation of priority queue
                        case 'p':
                            calc.lazyUpdateRecalcLimit = atoi(optarg);
                            break;
                        // hop-/degree-limits specified as comma-separated list,
                        // e.g. "1,3.3,2,10,3,10,5" means hop-limit 1 until degree 3.3,
                        //                          then hop-limit 2 until degree 10,
                        //                          then hop-limit 3 until degree 10,
                        //                          then hop-limit 5
                        case 'k':
                            calc.maxHops.clear();
                            Command::createVector(string(optarg),calc.maxHops,(double)0);
                            break;
                        // perform update on all affected nodes after contraction
                        // this only works with hop-limits and is really slow
                        // should only be used for testing
                        case 'K':
                            calc.updateHops = true;
                            break;
                        // omit local edge reduction
                        case 'L':
                            calc.localReduceEdges = !calc.localReduceEdges;
                            break;
                        // export search graph in a text format
                        case 'E':
                            exportSearchGraphFile = string(optarg);
                            break;
                        // export search graph in binary format (no standardized format,
                        // can depend on switches in config.h)
                        case 'Z':
                            sgrFile = string(optarg);
                            break;
                        // export whole contraction hierarchy including node ordering
                        // and shortcuts to standardized format (../docu/chDocu.html)
                        case 'C':
                            chFile = string(optarg);
                            break;
                    }
                }

                // Output parameters to stdout
                VERBOSE( cout << "Order nodes with priority-queue..." << endl );
                VERBOSE( cout << "test " << test << endl );
                VERBOSE( cout << "edgeDiffMult " << calc.edgeDiffMult << endl );
                VERBOSE( cout << "searchSpaceMult " << calc.searchSpaceMult << endl );
                VERBOSE( cout << "betweennessAdd " << calc.betweennessAdd << endl );
                VERBOSE( cout << "reachAdd " << calc.reachAdd << endl );
                VERBOSE( cout << "delNeighbMult " << calc.delNeighbMult << endl );
                VERBOSE( cout << "newEdgesMult " << calc.newEdgesMult << endl );
                VERBOSE( cout << "maxSettledApprox " << calc.maxSettledApprox << endl );
                VERBOSE( cout << "lazyUpdateRecalcLimit " << calc.lazyUpdateRecalcLimit << endl );
                VERBOSE( cout << "maxSettledElim " << calc.maxSettledElim << endl );
                VERBOSE( cout << "maxHops "; copy(calc.maxHops.begin(), calc.maxHops.end(), ostream_iterator<double>(cout,"-")); cout << endl );
                VERBOSE( cout << "voronoiMult " << calc.voronoiMult << endl );
                VERBOSE( cout << "searchPathHopBorderMult " << calc.searchPathHopBorderMult << endl );
                VERBOSE( cout << "updateHops " << calc.updateHops << endl );
                VERBOSE( cout << "localReduceEdges " << calc.localReduceEdges << endl );
                VERBOSE( cout << "shortcutOriginalEdgeSumMult " << calc.shortcutOriginalEdgeSumMult << endl );
                VERBOSE( cout << "searchPathHopBorderOriginalEdges " << calc.searchPathHopBorderOriginalEdges << endl );


                // Warning messages if specified parameters do not work well with compilation parameters in config.h
                #ifdef COUNT_SHORTCUT_ORIGINAL_EDGES
                if (calc.shortcutOriginalEdgeSumMult == 0)
                {
                    cout << "WARNING: #define COUNT_SHORTCUT_ORIGINAL_EDGES in config.h set but shortcutOriginalEdgeSumMult (-e) == 0" << endl;
                    cout << "         this leads to decreased performance and increased memory consumption" << endl;
                }
                #else
                if (calc.shortcutOriginalEdgeSumMult != 0 ||  calc.searchPathHopBorderOriginalEdges )
                {
                    cerr << "#define COUNT_SHORTCUT_ORIGINAL_EDGES in config.h for shortcutOriginalEdgeSumMult (-e) != 0" << endl;
                    exit(1);
                }
                #endif

                // CH-export needs middle nodes of shortcuts, thus USE_CH_EXPAND
                #ifndef USE_CH_EXPAND
                if ( chFile != "" )
                {
                    cerr << "#define USE_CH_EXPAND in config.h to export CH-file (-C)." << endl;
                    exit(1);
                }
                #endif

                // The stringstream ss consists of a "-" separated list of all parameters and can be used as prefix
                // for filenames if "x" is specified as parameter value.
                stringstream ss;
                ss << "exact-" << calc.edgeDiffMult << "-" << calc.searchSpaceMult << "-" << calc.betweennessAdd << "-" << calc.delNeighbMult;
                ss << "-" << calc.newEdgesMult << "-" << calc.maxSettledApprox << "-" << calc.reachAdd << "-" << calc.lazyUpdateRecalcLimit;
                ss << "-" << calc.maxSettledElim << "-";
                copy(calc.maxHops.begin(), calc.maxHops.end(), ostream_iterator<double>(ss,"-"));
                ss << calc.voronoiMult << "-" << calc.searchPathHopBorderMult << "-" << calc.updateHops << "-" << calc.localReduceEdges;
                ss << "-" << calc.shortcutOriginalEdgeSumMult << "-" << calc.searchPathHopBorderOriginalEdges;
                if (output == "x")
                {
                    output = (ss.str() + ".hcn");
                }
                if (statsFile == "x")
                {
                    statsFile = (ss.str() + ".stats");
                }
                if (logFile == "x")
                {
                    logFile = (ss.str() + ".log");
                }
                if (exportSearchGraphFile == "x")
                {
                    exportSearchGraphFile = (ss.str() + ".ddsg");
                }
                if (sgrFile == "x")
                {
                    sgrFile = (ss.str() + ".sgr");
                }
                if (chFile == "x")
                {
                    chFile = (ss.str() + ".ch");
                }


                // Write parameters to log-file if requested.
                bool log = false;
                ofstream logOut;
                if (logFile != "")
                {
                    logOut.open(logFile.c_str());
                    if (!logOut.is_open()) { cerr << "Cannot write to " << logFile << endl; exit(1); }
                    log = true;

                    logOut << "test " << test << endl;
                    logOut << "edgeDiffMult " << calc.edgeDiffMult << endl;
                    logOut << "searchSpaceMult " << calc.searchSpaceMult << endl;
                    logOut << "betweennessAdd " << calc.betweennessAdd << endl;
                    logOut << "reachAdd " << calc.reachAdd << endl;
                    logOut << "delNeighbMult " << calc.delNeighbMult << endl;
                    logOut << "maxSettledApprox " << calc.maxSettledApprox << endl;
                    logOut << "lazyUpdateRecalcLimit " << calc.lazyUpdateRecalcLimit << endl;
                    logOut << "maxSettledElim " << calc.maxSettledElim << endl;
                    logOut << "maxHops "; copy(calc.maxHops.begin(), calc.maxHops.end(), ostream_iterator<double>(logOut,"-")); logOut << endl;
                    logOut << "voronoiMult " << calc.voronoiMult << endl;
                    logOut << "searchPathHopBorderMult " << calc.searchPathHopBorderMult << endl;
                    logOut << "updateHops " << calc.updateHops << endl;
                    logOut << "localReduceEdges " << calc.localReduceEdges << endl;
                    logOut << "shortcutOriginalEdgeSumMult " << calc.shortcutOriginalEdgeSumMult << endl;
                    logOut << "searchPathHopBorderOriginalEdges " << calc.searchPathHopBorderOriginalEdges << endl;
                    logOut << "statsFile " << statsFile << endl;
                }
                VERBOSE( cout << "output " << output << endl );


                /* ************* *
                 * Read in files *
                 * ************* */

                // Read original graph. Each node will be in level n. During the node ordering
                // each contratacted node will be put in its own level.
                VERBOSE( std::cout << "Open file " << file << "..." << std::endl; )
                std::ifstream in(file.c_str());
                if (!in.is_open()) { cerr << "Cannot open " << file << endl; exit(1); }
                datastr::graph::UpdateableGraph* graph = importGraphListOfEdgesUpdateable(in, false, false, "");
                in.close();


                // Read betwenness centrality if specified.
                vector<BetweennessValue>* betweenness = NULL;
                if (betweennessFile != "")
                {
                    betweenness = new vector<BetweennessValue>(graph->noOfNodes(),0)                    ;
                    VERBOSE( std::cout << "Import Betweenness centrality values from " << betweennessFile << "..." << std::endl; )
                    ifstream in2(betweennessFile.c_str());
                    if (!in2.is_open()) { cerr << "Cannot open " << betweennessFile << endl; exit(1); }
                    for (vector<BetweennessValue>::iterator iter = betweenness->begin(); iter < betweenness->end(); iter++) {
                        assert ( !in2.eof() );
                        in2 >> *iter >> std::ws;
                    }
                    in2.close();
                }

                // Read reach centrality if specified.
                vector<ReachValue>* reach = NULL;
                if (reachFile != "")
                {
                    reach = new vector<ReachValue>(graph->noOfNodes(),0);
                    VERBOSE( std::cout << "Import Reeach values from " << reachFile << "..." << std::endl; )
                    ifstream in2(reachFile.c_str());
                    if (!in2.is_open()) { cerr << "Cannot open " << reachFile << endl; exit(1); }
                    unsigned int buffer;
                    in2.read((char*)&buffer, sizeof(unsigned int)/sizeof(char));
                    assert( buffer == graph->noOfNodes() );
                    for (NodeID u = 0; u < graph->noOfNodes() ; u++) {
                        assert ( !in2.eof() );
                        in2.read((char*)&buffer, sizeof(unsigned int)/sizeof(char));
                        (*reach)[u] = buffer;
                    }
                    in2.close();
                }

                // Initialize pseudorandom number generator allow reproduction of results.
                // Currently, no random data is used during node ordering but it does not harm
                // either to specifiy the random seed.
                srand(22);

                // Object for node ordering.
                ProcessingNodeOrder c(graph);

                // Prepare output files for witness and shortcut information that can
                // be used in a later hierarchy construction step.
                // To export the wheter witness information, a template parameter
                // in the ProcessingNodeOrder class needs to be set. There is
                // currently node switch in ../config.h
                ofstream outShortcuts;
                ofstream outWitnesses;
                if ( c.savesShortcutsWitnesses() )
                {
                    outShortcuts.open((ss.str()+".shortcuts").c_str());
                    if (!outShortcuts.is_open()) { cerr << "Cannot write to " << (ss.str()+".shortcuts") << endl; exit(1); }
                    c.storeShortcuts(outShortcuts);
                    outWitnesses.open((ss.str()+".witnesses").c_str());
                    if (!outWitnesses.is_open()) { cerr << "Cannot write to " << (ss.str()+".witnesses") << endl; exit(1); }
                    c.storeWitnesses(outWitnesses);
                }

                /* ***************************************** *
                 * Node ordering and hierarchy construction. *
                 * ***************************************** */
                // This is the main step, a node order is created including a contraction hierarchy.
                VERBOSE( double timeCalc = timestamp(); )
                VERBOSE( cout << "Create node order and hierarchy using a priority queue..." << endl; )
                c.createHierarchy(calc, betweenness, reach, statsFile, test);
                VERBOSE( timeCalc = timestamp() - timeCalc; )

                // Now write the node order to a hcn-file.
                VERBOSE( cout << "Write node order (levels) to " << output << " ..." << endl; )
                ofstream out(output.c_str());
                if (!out.is_open()) { cerr << "Cannot write to " << output << endl; exit(1); }
                c.writeLevels(out);
                out.close();
                if ( outShortcuts.is_open() ) outShortcuts.close();
                if ( outWitnesses.is_open() ) outWitnesses.close();

                // additional output: search-graph to text format
                if ( exportSearchGraphFile != "" )
                {
                    VERBOSE( cout << "Export search-graph to " << exportSearchGraphFile << " ..." << endl; )
                    ofstream outSearchGraph(exportSearchGraphFile.c_str());
                    if (!outSearchGraph.is_open()) { cerr << "Cannot write to " << exportSearchGraphFile << endl; exit(1); }
                    exportSearchGraph(outSearchGraph, graph);
                }

                // additional output: search-graph to binary format (depends on switches in config.h)
                if ( sgrFile != "" )
                {
                    VERBOSE( cout << "Export search-graph to " << sgrFile << " ..." << endl; )
                    datastr::graph::SearchGraph* searchGraph = new datastr::graph::SearchGraph(graph, datastr::graph::SGNO_ORIGINAL);
                    ofstream out(sgrFile.c_str());
                    if (!out.is_open()) { cerr << "Cannot write to " << sgrFile << endl; exit(1); }
                    searchGraph->serialize(out);

                }

                // additional output: contraction hierarchy to binary format (../docu/chDocu.html)
                if ( chFile != "" )
                {
                    VERBOSE( cout << "Serialize contraction hierarchy to " << chFile << " ..." << endl; )
                    ofstream out(chFile.c_str());
                    if (!out.is_open()) { cerr << "Cannot write to " << chFile << endl; exit(1); }
                    exportContractionHierarchy( out, graph );
                }

                // statistic: count the total number of expanded edges
                // Only works correct for a n-level hierarchy otherwise edges may be counted twice.
                long totalExpandedEdges = 0;
                #ifdef COUNT_SHORTCUT_ORIGINAL_EDGES
                for ( NodeID u = 0; u < graph->noOfNodes(); u++ )
                {
                    for ( EdgeID e = graph->firstLevelEdge(u); e < graph->lastEdge(u); e++ )
                    {
                        const Edge& edge = graph->edge(e);
                        if ( edge.isShortcut() )
                        {
                            totalExpandedEdges += edge.shortcutOriginalEdgeCount();
                        }
                    }
                }
                #endif


                // Output runtimes and some statistics
                VERBOSE( cout.unsetf( ios_base::fixed ); )
                VERBOSE( cout << "node ordering time: " << timeCalc << " s" << endl );
                VERBOSE( cout << "max avg degree: " << c.lastMaxAvgDegree() << endl );
                VERBOSE( cout << "#edges: " << graph->noOfExistingEdges() << endl );
                VERBOSE( cout << "#total expanded edges: " << totalExpandedEdges << endl );

                if (log)
                {
                    logOut << "node ordering time: " << timeCalc << " s" << endl;
                    logOut << "max avg degree: " << c.lastMaxAvgDegree() << endl;
                    logOut << "#edges: " << graph->noOfExistingEdges() << endl;
                    logOut << "#total expanded shortcut edges: " << totalExpandedEdges << endl;
                    logOut.close();
                }

                // Cleanup memory.
                delete graph;
                if (betweenness != NULL) delete betweenness;
                if (reach != NULL) delete reach;
            }

            /* ******************** *
           * level classification *
           * ******************** */
            // Performs a transformation of a hcn-file with a single level per node
            // (== node order) into a hcn-file with only a few levels (option -l)
            // using that node order.
            // This is currently used to create hcn-files that can be used
            // with the implementation of HNR and TNR of Dominik Schultes.
            else if (opt1 == 'l') {
                string exact("in.hcn");     // input graph file
                string output("out.hcn");   // output file (counter of nodes)
                string levels("n//2x16");
                unsigned int a = 0;
                int opt2 = -1;
                while ((opt2 = getopt(argc, argv, "h:o:l:a:")) != -1)
                {
                    switch(opt2)
                    {
                        case 'h':
                            exact = string(optarg);
                            break;
                        case 'o':
                            output = string(optarg);
                            break;
                        case 'l':
                            levels = string(optarg);
                            break;
                        case 'a':
                            a = atoi(optarg);
                            break;
                    }
                }

                VERBOSE( cout << "levels " << levels << endl; )
                VERBOSE( cout << "Read from " << exact << " ..." << endl; )
                ifstream in(exact.c_str());
                if (!in.is_open()) { cerr << "Cannot open " << exact << endl; exit(1); }
                VERBOSE( cout << "Write to " << output << " ..." << endl; )
                ofstream out(output.c_str());
                if (!out.is_open()) { cerr << "Cannot write to " << output << endl; exit(1); }
                ProcessingNodeOrder::partitionExactInLevels(in,out,levels,a);
                in.close();
                out.close();
            }

            // Test a hcn-file. Either the level of a given node
            // is printed to stdout (-f and -n),
            // or the number of nodes in each existing level is printed
            // to stdout (-h).
            // Usefull for debugging purposes.
            if (opt1 == 't') {
                string file("");      // input graph file
                NodeID node = SPECIAL_NODEID;
                int opt2 = -1;
                string hcnFile("");
                while ((opt2 = getopt(argc, argv, "f:n:h:")) != -1)
                {
                    switch(opt2)
                    {
                        case 'f':
                            file = string(optarg);
                            break;
                        case 'n':
                            node = atoi(optarg);
                            break;
                        case 'h':
                            hcnFile = string(optarg);
                            break;
                    }
                }

                if (file != "")
                {
                    VERBOSE( cout << "node " << node << endl; )
                    VERBOSE( std::cout << "Open file " << file << "..." << std::endl; )
                    std::ifstream in(file.c_str());
                    if (!in.is_open()) { cerr << "Cannot open " << file << endl; exit(1); }

                    unsigned int buffer;
                    in.read((char*)&buffer, sizeof(unsigned int)/sizeof(char));
                    assert( node < buffer );
                    for ( unsigned int i = 0; i <= node; i++ )
                    {
                        in.read((char*)&buffer, sizeof(unsigned int)/sizeof(char));
                    }
                    cout << buffer << endl;
                }
                else if (hcnFile != "")
                {
                    vector<unsigned int> counter;
                    VERBOSE( std::cout << "Open file " << hcnFile << "..." << std::endl; )
                    std::ifstream in(hcnFile.c_str());
                    if (!in.is_open()) { cerr << "Cannot open " << hcnFile << endl; exit(1); }

                    unsigned int buffer;
                    in.read((char*)&buffer, sizeof(unsigned int)/sizeof(char));
                    unsigned int n = buffer;
                    VERBOSE( std::cout << "n = " << n << std::endl; )
                    for ( unsigned int i = 0; i < n; i++ )
                    {
                        in.read((char*)&buffer, sizeof(unsigned int)/sizeof(char));
                        while (counter.size() <= buffer) counter.push_back(0);
                        counter[buffer]++;
                    }
                    for ( unsigned int i = 0; i < counter.size(); i++ )
                    {
                        cout << "level " << i << ": " << counter[i] << endl;
                    }
                }
            }

            return 0;
        }
    };
}
#endif // _COMMAND_NODEORDER_H
