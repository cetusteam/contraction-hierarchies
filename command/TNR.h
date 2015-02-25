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

#ifndef _COMMAND_TNR_H
#define _COMMAND_TNR_H

#include "../Command.h"
#include "../datastr/graph/UpdateableGraph.h"
#include "../datastr/graph/SearchGraph.h"
#include "../io/createGraph.h"
#include "../io/output.h"
#include "../processing/DijkstraCH.h"

// Tests for applying CH to Transit Node Routing, these are
// not relevant for preprocess TNR, they only calculate
// interesting statistics, e.g. number of access nodes.
namespace command {
    class TNR : public Command {
        public:

        int main(int argc, char *argv[]) {

            string sgrFile("");
            
            int opt2 = -1;
            while ((opt2 = getopt(argc, argv, "z:")) != -1)
            {
                switch(opt2)
                {
                    case 'z':
                        sgrFile = string(optarg);
                        break;
                        
                }
            }

            stringstream ssName; 
            ssName << sgrFile;

            datastr::graph::SearchGraph* searchGraph;
            VERBOSE( cout << "Deserialize search-graph from " << sgrFile << " ..." << endl; )
            VERBOSE( cout << "Note: Node-ID must correspond to Level-ID!" << endl; )

            ifstream in(sgrFile.c_str());
            if (!in.is_open()) { cerr << "Cannot open " << sgrFile << endl; exit(1); }
            searchGraph = new datastr::graph::SearchGraph(in);

            processing::DijkstraCH<datastr::graph::SearchGraph, DynQueryPQueue, 2, true /*stall-on-demand*/, false/*deep stall-on-demand*/> dijkstra(searchGraph);

            cout << "Test |T|..." << endl;
            vector<bool> flag(searchGraph->noOfNodes(), false);
            vector<unsigned long long> sum(searchGraph->noOfNodes(), 0);
            NodeID T = searchGraph->noOfNodes();
            T >>= 1;
            unsigned int x = 0;
            while ( T > 1 )
            {
                x++;
                T >>= 1;
            }
            vector<NodeID> accessNodesMax(x, 0);
            vector<NodeID> accessNodesSum(x, 0);
            vector<NodeID> accessNodesCurrent(x, 0);
            vector<NodeID> settledNodesMax(x, 0);
            vector<NodeID> settledNodesSum(x, 0);
            vector<NodeID> settledNodesCurrent(x, 0);

            VERBOSE( Percent percent( searchGraph->noOfNodes() ) );
            double time = timestamp();
            
            for ( NodeID u = 0; u < searchGraph->noOfNodes(); u++ )
            {
                dijkstra.bidirSearch(u, SPECIAL_NODEID);
                const vector<NodeID>& settledNodes = dijkstra.settledNodes( 0 );
                
                for ( NodeID i = settledNodes.size()-1; ; i-- )
                {
                    NodeID v = settledNodes[i];
                    NodeID index = dijkstra.isReached( 0, v );
                    assert( index > 0 );
                    if (!flag[v] && !dijkstra.pqData( 0, index ).stalled() )
                    {
                        sum[v]++;
                    }
                    NodeID w = dijkstra.parentOf( v, 0 );
                    if (!flag[w]) flag[w] = true;
                    if (i == 0) break;
                }
                for ( NodeID i = 0; i < settledNodes.size(); i++ )
                {
                    flag[settledNodes[i]] = false;
                }
                
                fill( accessNodesCurrent.begin(), accessNodesCurrent.end(), 0);
                fill( settledNodesCurrent.begin(), settledNodesCurrent.end(), 0);
                
                for ( NodeID i = 0; i < settledNodes.size(); i++ )
                {
                    NodeID v = settledNodes[i];
                    NodeID index = dijkstra.isReached( 0, v );
                    assert( index > 0 );
                    if ( dijkstra.pqData( 0, index ).stalled() ) continue;
                    NodeID w = dijkstra.parentOf( v, 0 );
                    
                    T = searchGraph->noOfNodes();
                    T >>= 1;
                    x = 0;
                    while ( T > 1 )
                    {
                        NodeID border = searchGraph->noOfNodes()-T;
                        // if v is an access node to upper plain with T nodes
                        if ( v >= border && w < border )
                        {
                            accessNodesCurrent[x]++;
                        }
                        else if ( v < border )
                        {
                            settledNodesCurrent[x] = i;   
                        }
                        x++;
                        T >>= 1;
                    }
                    
                }

                for ( x = 0; x < accessNodesCurrent.size(); x++ )
                {
                    if ( accessNodesCurrent[x] > accessNodesMax[x] ) accessNodesMax[x] = accessNodesCurrent[x];
                    accessNodesSum[x] += accessNodesCurrent[x];
                    
                    if ( settledNodesCurrent[x] > settledNodesMax[x] ) settledNodesMax[x] = settledNodesCurrent[x];
                    settledNodesSum[x] += settledNodesCurrent[x];
                }
            
                dijkstra.clear();
                VERBOSE( percent.printStatus(u) );
            }
            time = timestamp() - time;
            VERBOSE( cout << "Test took " << time << " seconds." << endl; )

            for ( unsigned int i = sum.size()-2; ; i-- )
            {
                sum[i] += sum[i+1];
                if (i == 0) break;
            }
                
            ofstream log((ssName.str()+".tnr.access-nodes-1").c_str());
            if (!log.is_open()) { cerr << "Cannot write to " << (ssName.str()+".tnr.access-nodes-1") << endl; exit(1); }
            for ( NodeID u = 0; u < sum.size(); u++ )
            {
                log << sum[u] << endl;
            }
            log.close();
            
            log.open((ssName.str()+".tnr.access-nodes-sum").c_str());
            if (!log.is_open()) { cerr << "Cannot write to " << (ssName.str()+".tnr.access-nodes-sum") << endl; exit(1); }
            for ( x = 0; x < accessNodesCurrent.size(); x++ )
            {
                log << accessNodesSum[x] << endl;            
            }
            log.close();

            log.open((ssName.str()+".tnr.access-nodes-max").c_str());
            if (!log.is_open()) { cerr << "Cannot write to " << (ssName.str()+".tnr.access-nodes-max") << endl; exit(1); }
            for ( x = 0; x < accessNodesCurrent.size(); x++ )
            {
                log << accessNodesMax[x] << endl;            
            }
            log.close();
            

            log.open((ssName.str()+".tnr.settled-nodes-sum").c_str());
            if (!log.is_open()) { cerr << "Cannot write to " << (ssName.str()+".tnr.settled-nodes-sum") << endl; exit(1); }
            for ( x = 0; x < settledNodesCurrent.size(); x++ )
            {
                log << settledNodesSum[x] << endl;            
            }
            log.close();

            log.open((ssName.str()+".tnr.settled-nodes-max").c_str());
            if (!log.is_open()) { cerr << "Cannot write to " << (ssName.str()+".tnr.settled-nodes-max") << endl; exit(1); }
            for ( x = 0; x < settledNodesCurrent.size(); x++ )
            {
                log << settledNodesMax[x] << endl;            
            }
            log.close();

            delete searchGraph;

            return 0;

        }
        
    };
}
#endif // _COMMAND_TNR_H
