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

#ifndef _PROCESSING_DIJKSTRACH_H
#define _PROCESSING_DIJKSTRACH_H

#include <queue>
#include <stack>

#include "../datastr/graph/graph.h"
#include "../datastr/graph/pqNode.h"
#include "../datastr/searchSpaces.h"
#include "../datastr/graph/UpdateableGraph.h"
#include "../datastr/graph/SearchGraph.h"

namespace processing {

// The various Many-To-Many Modes. Used as template parameter for Dijkstra.
    const int MTMM_NONE = 0;
    const int MTMM_FW = 1;
    const int MTMM_BW = 2;


    /**
     * A versatile implementation of Dijkstra's algorithm.
     * This class originates from the DijkstraTemplate class in dijkstra.h
     * by Dominik Schultes and is modified for contraction hierarchies.
     * Can be used for
     * - normal version of Dijkstra's algorithm
     * - bidirectional version of Dijkstra's algorithm
     * - local searches during the construction of contraction Hierarchies
     * - distribution of Voronoi regions (based on shortest paths) among its neighboring regions
     * - multilevel query ("hwy search")
     * @param Graph the graph type (mostly SearchGraph)
     * @param PQueue the pqueue type (Normal, Hwy, or Constr)
     * @param searchDirections the number of search directions
     *                         (1 for normal Dijkstra, construction;
     *                          2 for bidir Dijkstra, hwy search)
     * @param stallOnDemand enable stall-on-demand technique
     * @param deepStallOnDemand also use parent pointers for stall-on-demand
     * @param manyToManyMode various modes related to many-to-many computations
     * @param localSearch local search, stop criterions like settled nodes limit or reached targets can be used.
     * @param countHops required for local search with hop-limits
     */
    template < typename Graph,
    typename PQueue,
    int searchDirections,
    bool stallOnDemand = false,
    bool deepStallOnDemand = false,
    int manyToManyMode = MTMM_NONE,
    bool localSearch = false,
    bool countHops = false>
    class DijkstraCH
    {
        /** Outputs data about the current state of the search for debugging purposes. */
        friend ostream& operator<<( ostream& os, const DijkstraCH &dijkstra ) {
            for (NodeID i=0; i<dijkstra._graph->noOfNodes(); i++) {
                os << "Node " << i << ": ";
                for (int j=0; j<searchDirections; j++) {
                    NodeID index = dijkstra.isReached(j, i);
                    if (index) {
                        if (! dijkstra.isSettled(j, i) ) os << "*";
                        os << dijkstra.pqKey(j, index) << "(p "
                                << dijkstra.pqData(j, index).parentNode() << " / "
                                << dijkstra.pqData(j, index).parentEdge() <<  ") | ";
                    }
                    else {
                        os << "X | ";
                    }
                }
                os << endl;
            }
            return os;
        }


public:
    typedef typename PQueue::PQElement PQElement;
    typedef typename PQueue::PQData PQData;
    typedef typename Graph::MyNode MyNode;

    /**
     * Constructor.
     * @param graph a pointer to the graph
     */
    DijkstraCH(Graph* graph)
    : _graph(graph),
      _upperBound(Weight::MAX_VALUE)
    {
        assert( (searchDirections >= 1) && (searchDirections <= 2) );
    }



    /**
    * Unidirectional search from a given node s.
    * Used by the construction and the normal version of Dijkstra's algorithm.
    * @param s source node
    * @param searchID search directions

    * The following parameters are only for local searches, template parameter localSearch enabled:
    * @param ignore ignore this node during the search (do not relax edges to this node)
    * @param maxKey stop the search if the smallest tentative distance in the pqueue is larger
    *               than this limit. (only for localSearch)
    * @param noOfTargets stop the search if this number of nodes with target-flag (in the graph)
    *               have been settled. (only for localSearch)
    * @param maxSettled stop the search if this number of nodes is settled. 0 = disabled (only for localSearch)
    * @param maxHops do only find paths with this number of hops. Template parameter countHops
    *               needs to be enabled. 0 = disabled (only for localSearch)
    */
    void searchWithoutTarget(const NodeID s, const int searchID = 0, const NodeID ignore = SPECIAL_NODEID,
                                const EdgeWeight maxKey = 0, const NodeID noOfTargets = 0,
                                const NodeID maxSettled = 0, const NodeID maxHops = 0)
    {
        assert( checkClean() );
        assert( (searchID >= 0) && (searchID <= 1) );

        insertStartNode(searchID, s);
        CurrentNode u;

        NodeID targetsRemaining = noOfTargets;
        NodeID remaining = maxSettled;
        while (! pqueue(searchID).empty() )
        {
            // *** deleteMin ***
            deleteMin(searchID, u);

            // a local search is limited
            if (localSearch)
            {
                // abort local search if distance too far
                if (u.dist > maxKey) break;

                // search until all targets reached
                if (pqData(searchID, u.index).isTarget())
                {
                    targetsRemaining--;
                    if (targetsRemaining == 0) break;
                }

                // only consider limited settled nodes, then stop search
                remaining--;
                if (remaining == 0) break;

                if (countHops)
                {
                    assert( maxHops == 0 || u.hops <= maxHops );
                    assert( u.nodeID == s || u.hops > 0 );
                    // do not relax edges after reaching maxHops count
                    if (maxHops > 0 && u.hops >= maxHops) continue;
                }
            }

            // *** relaxEdges ***
            relaxEdges(searchID, u, ignore);

        }
    }

    /**
    * Bidirectional search from a given node s to a given node t (and vice versa).
    * Used by the bidirectional version of Dijkstra's algorithm
    * and the multilevel query ("hwy search").
    * @param s the source node; SPECIAL_NODEID means: ignore forward direction
    * @param t the target node; SPECIAL_NODEID means: ignore backward direction
    * @return distance from s to t or "infinity" if s and t aren't connected
    */
    EdgeWeight bidirSearch(NodeID s, NodeID t) {
        assert( searchDirections == 2 );
        assert( checkClean() );


        // initialization of the hwy search
        _finished[0] = false;
        _finished[1] = false;

        _concludeHwySearch = false;

        if (s != SPECIAL_NODEID) insertStartNode(0,s);
        if (t != SPECIAL_NODEID) insertStartNode(1,t);

        CurrentNode u;
        int searchID = 0;
        bool pqEmpty[2] = {pqueue(0).empty(), pqueue(1).empty()};
        while (! (pqEmpty[0] && pqEmpty[1])) {

            // different strategies: "Which pqueue should be preferred ?"
            PQ_BIDIR_MIN( if (pqueue(1).min() < pqueue(0).min()) searchID = 1; else searchID = 0 );
            PQ_BIDIR_SIZE(
                    if (((pqueue(1).size() < pqueue(0).size()) && (! pqEmpty[1])) ||
                    pqEmpty[0]) searchID = 1; else searchID = 0 );
            PQ_BIDIR_ALTERNATE( if (! pqEmpty[1-searchID]) searchID = 1 - searchID );

            // If one search direction is finished, always choose the other one.
            if (_finished[0]) searchID = 1;
            if (_finished[1]) searchID = 0;

            // If the chosen pqueue is empty, there is nothing to do.
            if (pqEmpty[searchID]) break;

            // *** deleteMin ***
            bool abort = deleteMin(searchID, u);

            // *** relaxEdges ***
            relaxEdges( searchID, u );

            // Mark that we "want" to abort now and we will only go on
            // until we are sure that we haven't missed anything.
            if (abort) _concludeHwySearch = true;

            // Check a simple abort-on-success criterion.
            if ( pqueue(searchID).min() > _upperBound) {
                _finished[searchID] = true;
                VERBOSE2( cerr << searchID << " FINISHED ("
                        << pqueue(searchID).min() << " " << _upperBound << ")" << endl );
            }

            // Iff both search directions are finished, we may abort.
            if (_finished[0] && _finished[1]) break;

            pqEmpty[searchID] = pqueue(searchID).empty();
        }

        processUnsettledNodes(0);
        processUnsettledNodes(1);

        return _upperBound;
    }

    /**
    * Returns the distance from the source node s of the search that has been performed last
    * to the given node t.
    * @param t the target node
    * @param searchID the search direction that should be considered; special value -1 is used
    *                 to indicate that the result of a bidirectional search should be returned
    * @return distance from s to t or "infinity" if s and t aren't connected
    */
    EdgeWeight distanceTo(NodeID t, int searchID = 0) const {
        // special case: return result of the bidirectional search
        if (searchID == -1) return _upperBound;

        NodeID index = isReached(searchID, t);
        if (index) {
            return pqKey(searchID, index);
        }
        else {
            return Weight::MAX_VALUE; // "infinity" if not connected
        }
    }

    /** Returns the node ID of the parent of the given node t. */
    NodeID parentOf(const NodeID t, const int searchID = 0) const {
        const NodeID index = isReached(searchID, t);
        assert( index );
        return pqData(searchID, index).parentNode();
    }

    /**
    * Determines the shortest path from the source node s of the search
    * that has been performed last to the given node t.
    * @param path a reference to the path object that accepts the result
    * @param t the target node (not relevant if searchID = -1)
    * @param searchID the search direction that should be considered; special value -1 is used
    *                 to indicate that the result of a bidirectional search should be considered
    * @param forward true/false = return the path from s to t / from t to s
    * @param expand flag: expand the contracted path (i.e., return the expanded path)
    *                 This is a recursive unpacking routine that replaces each shortcut
    *                 by its two originating edges.
    * @return a reference to the given path object that accepts the result
     */
    Path& pathTo(Path& path, NodeID t, int searchID = 0, bool forward = true, bool expand = false) {
        path.clear();
        if (searchID == -1) { // special case: return result of the bidirectional search
            if (_upperBound == Weight::MAX_VALUE) {
                path.setNotConnected();
                return path;
            }

            if ( expand )
            {
                // The query algorithm is bidirectional. We have two
                // partial paths, one of the forward search and one of the
                // backward search.
                // Get forward path reversed, not expanded because
                // this is the natural way using the parent pointers.
                Path pathForwardReverse;
                pathTo(pathForwardReverse, _via, 0, false, false);
                // Expand forward path, this writes the expanded
                // path in a second object and implicitely reverses
                // the path.
                expandPath<true,false,true>(pathForwardReverse, path);
                assert( pathForwardReverse.length() == path.length() );

                // The backward path will be processed in code below
                // that is used for the unidirectional case.
                // We just need to set these variables accordingly.
                searchID = 1;
                t = _via;
                forward = !forward;

            }
            else
            {

                // first part of the path from s to the intermediate node 'via'
                pathTo(path, _via, 0, true, expand);

                // second part of the path from the intermediate node 'via' to t
                Path path2;
                pathTo(path2, _via, 1, false, expand);

                // combine both parts
                path.add(path2);

                // We have the path forwards.
                if (! forward) path.reverse(); // reverse, if required

                return path;
            }

        }
        else
        {
            // add first node to path, all following nodes are added
            // with an edge
            path.addFirstNode(t);
        }

        if ( expand )
        {
            // simple expand routine, because shortcut
            // from parent->child is original parent->middle->child

            // If we found a shortcut, we replace it by the two originating
            // edges. One edge is put on a stack and the second one is processed
            // in the next pass of the loop. _expandPathStack is this stack.
            assert( _expandPathStack.empty() );
            NodeID parent = SPECIAL_NODEID;
            EdgeID e = SPECIAL_NODEID;

            // If we found a shortcut, we replace it by the two originating
            // edges. One edge is put on a stack and the second one is processed
            // in the next pass of the loop. processSecondHalf indicates this.
            bool processSecondHalf = false;

            // Main loop, will run until all edges are replaced by shorcut edges.
            while (true)
            {
                // This variable determines which of the two edges of a shortcut is
                // accessed first in the edge array. Intuitively, the edge that
                // is processed first should be accessed first. But tests indicate
                // otherwise.
                //bool shortcutParent1Child2 = !forward;

                // To process an edge, we need the edge id of this edge
                // and the parent in the Dijkstra search graph to ensure
                // a correct order of the edges in the unpacked path.
                // The current position of the processing is marked by node t.
                // Either in the last pass of the loop an edge was specified
                // (processSecondHalf == true)...
                if ( !processSecondHalf )
                {
                    // ... or we take the next edge from the path found by the query algorithm if
                    // the stack is empty ...
                    if ( _expandPathStack.empty() )
                    {
                        NodeID index = isReached(searchID, t);
                        if ( ! index ) {
                            path.setNotConnected();
                            return path;
                        }
                        // determine parent
                        parent = pqData(searchID, index).parentNode();
                        if (parent == t) break;
                        e = pqData(searchID, index).parentEdge();
                    }
                    // ... otherwise we take the edge on top of the stack.
                    else
                    {
                        parent = _expandPathStack.top().first;
                        e = _expandPathStack.top().second;
                        _expandPathStack.pop();
                        //shortcutParent1Child2 = !shortcutParent1Child2;
                    }
                }
                // process edge
                expandEdge( path, t, parent, e, processSecondHalf, forward );
            }
            assert( _expandPathStack.empty() );
        }

        // No path expansion requested, only return path found by query algorith.
        else
        {
            // determine path backwards, starting with t
            while (true) {
                NodeID index = isReached(searchID, t);
                if ( ! index ) {
                    path.setNotConnected();
                    return path;
                }
                // determine parent
                NodeID parent = pqData(searchID, index).parentNode();
                if (parent == t) break;
                t = parent;
                EdgeID e = pqData(searchID, index).parentEdge();

                path.add( t, _graph->edge(e).weight(), e, _graph->edge(e).isShortcut() );
            }
        }

        // We have the path backwards.
        // reverse, if required
        // (Note that the expander cannot deal with a backward path since the
        //  included edge IDs would not correspond to this direction.)
        if (forward) path.reverse();
        return path;
    }

    /**
     * Expands a given path. This is almost the same path expansion routine as in pathTo() but
     * reads the path that should be expanded from a Path object instead of the Dijkstra search tree.
     * This is used for the forward search since an extracted path would be read from the Dijkstra
     * search three reversed, meaning from the target (or via node) to the source node. We do not want
     * to reverse the expanded path (100 times longer than the path found by the query algorithm).
     * So we store the reversed path found by the query algorithm, process it backwards to get an expanded
     * path forwards.
     * @param reverse should be expanded path be reversed to the input path
     * @param forward is the input path forward from source node to target node (or via node)
     * @param newPath should the pathExpanded return path be cleared and newly initalized
     * @param path a reference to the original path object
     * @param pathExpanded a reference to the path object that accepts the expanded result
     */
    template<bool reverse, bool forward, bool newPath>
    void expandPath(const Path& path, Path& pathExpanded)
    {
        if (newPath) pathExpanded.clear();
        // simple expand routine, because shortcut from parent->child is original parent->middle->child
        // assertion: edge-array is sorted by target
        if (path.isNotConnected())
        {
            pathExpanded.setNotConnected();
            return;
        }

        assert( !path.empty() );

        // The current position of the processing is marked by node t.
        // The current index in the input path is variable index.
        NodeID t;
        NodeID index;
        if ( !reverse )
        {
            t = path.firstNode();
            index = 0;
        }
        {
            t = path.lastNode();
            index = path.noOfEdges() - 1;
        }
        if (newPath) pathExpanded.addFirstNode(t);


        // If we found a shortcut, we replace it by the two originating
        // edges. One edge is put on a stack and the second one is processed
        // in the next pass of the loop. _expandPathStack is this stack.
        assert( _expandPathStack.empty() );
        NodeID parent = SPECIAL_NODEID;
        EdgeID e = SPECIAL_NODEID;



        // If we found a shortcut, we replace it by the two originating
        // edges. One edge is put on a stack and the second one is processed
        // in the next pass of the loop. processSecondHalf indicates this.
        bool processSecondHalf = false;

        // Main loop, will run until all edges are processed.
        // These can be on the stack or on the input path.
        while (!(_expandPathStack.empty() &&
                (( reverse && (index == SPECIAL_NODEID) )
                || ( !reverse && index >= path.noOfEdges() ))))
        {
            // This variable determines which of the two edges of a shortcut is
            // accessed first in the edge array. Intuitively, the edge that
            // is processed first should be accessed first. But tests indicate
            // otherwise.
            //bool shortcutParent1Child2 = !forward;

            // To process an edge, we need the edge id of this edge
            // and the parent in the Dijkstra search graph to ensure
            // a correct order of the edges in the unpacked path.
            // The current position of the processing is marked by node t.
            // Either in the last pass of the loop an edge was specified
            // (processSecondHalf == true)...
            if ( !processSecondHalf )
            {
                // ... or we take the next edge from the path found by the query algorithm if
                // the stack is empty ...
                if ( _expandPathStack.empty() )
                {
                    assert( index != SPECIAL_NODEID );
                    if ( !reverse )
                    {
                        e = path.edge(index);
                        index++;
                        parent = path.node(index);
                    }
                    else
                    {
                        e = path.edge(index);
                        parent = path.node(index);
                        if ( index > 0 ) index--;
                        else index = SPECIAL_NODEID;
                        //shortcutParent1Child2 = !shortcutParent1Child2;
                    }
                }
                // ... otherwise we take the edge on top of the stack.
                else
                {
                    parent = _expandPathStack.top().first;
                    e = _expandPathStack.top().second;
                    _expandPathStack.pop();
                    //shortcutParent1Child2 = !shortcutParent1Child2;
                }
            }

            // process edge
            expandEdge( pathExpanded, t, parent, e, processSecondHalf, forward );
        }
        assert( _expandPathStack.empty() );
    }

    /**
     * This subroutine is common to pathTo() and expandPath(), give an edge by its edge id e,
     * it is either a shortcut and gets expanded or not and will be appended to path.
     * A shortcut is expanded into two edges, the first is put on the stack, the second
     * is prepared to be processed next (processSecondHalf flag).
     * The variables have the same name as the ones used in pathTo() and expandPath().
     */
    void expandEdge(Path& path, NodeID& t, NodeID& parent, EdgeID& e, bool& processSecondHalf, const bool forward)
    {
        assert( parent != SPECIAL_NODEID );
        assert( e != SPECIAL_NODEID );
        const Edge& edge = _graph->edge(e);

        // A shortcut parent->child is split into two edges parent->middle->child
        if ( edge.isShortcut() )
        {
            NodeID middle = edge.shortcutMiddle();
            EdgeID firstEdge = _graph->firstLevelEdge(middle);

            // edge.shortcutEdge1() and edge.shortcutEdge2() are a relative index
            // into the adjacency array of the middle node. They are valid if
            // they are != edge.shortcutEdgeLimit()
            EdgeID eParent, eChild;
            if ( true /*|| shortcutParent1Child2*/ )
            {
                eParent = edge.shortcutEdge1();
                eChild = edge.shortcutEdge2();
            }
            else
            {
                eParent = edge.shortcutEdge2();
                eChild = edge.shortcutEdge1();
            }

            // eParent und eChild are valid relative indices
            if ( eParent != edge.shortcutEdgeLimit() && eChild != edge.shortcutEdgeLimit() )
            {
                eParent += firstEdge;
                eChild  += firstEdge;
                // Switch child <-> parent if required, the edges in the
                // expanded path need to be in the right order.
                if ( _graph->edge(eChild).target() == parent )
                {
                    EdgeID dummy = eParent;
                    eParent = eChild;
                    eChild = dummy;
                }
            }

            // at least on of eParent or eChild is invalid, need to scan through the EdgeID
            // array of the middle node to find the required edges. This can happen if there
            // are more than shortcutEdgeLimit() edges incident to a node.
            else
            {
                if ( eParent == edge.shortcutEdgeLimit() )
                {
                    eParent = SPECIAL_NODEID;
                    if ( eChild == edge.shortcutEdgeLimit() ) eChild = SPECIAL_NODEID;
                    else
                    {
                        eChild += firstEdge;
                        // Switch child <-> parent if required, the edges in the
                        // expanded path need to be in the right order.
                        if ( _graph->edge(eChild).target() == parent )
                        {
                            eParent = eChild;
                            eChild = SPECIAL_NODEID;
                        }
                        // Relative index to the child node was wrong.
                        else if ( _graph->edge(eChild).target() != t )
                        {
                            eChild = SPECIAL_NODEID;
                        }
                    }

                }
                else
                {
                    eParent += firstEdge;
                    // Switch child <-> parent if required, the edges in the
                    // expanded path need to be in the right order.
                    if ( _graph->edge(eParent).target() == t )
                    {
                        eChild = eParent;
                        eParent = SPECIAL_NODEID;
                    }
                    // Relative index to the parent node was wrong.
                    else if ( _graph->edge(eParent).target() != parent )
                    {
                        eParent = SPECIAL_NODEID;
                    }
                }

                NodeID directionParent = forward ? 1 : 0;
                NodeID directionChild  = 1-directionParent;
                EdgeID lastEdge = _graph->lastEdge(middle);
                // scan through all edges, starting at shortcutEdgeLimit since otherwise
                // the relative index would be valid.
                for ( EdgeID eMiddle = firstEdge + edge.shortcutEdgeLimit(); eMiddle < lastEdge; eMiddle++ )
                {
                    const Edge& edgeMiddle = _graph->edge(eMiddle);
                    if ( edgeMiddle.isDirected( directionParent ) && edgeMiddle.target() == parent )
                    {
                        eParent = eMiddle;
                        if ( eChild != SPECIAL_NODEID ) break;
                    }
                    if ( edgeMiddle.isDirected( directionChild ) && edgeMiddle.target() == t )
                    {
                        eChild = eMiddle;
                        if ( eParent != SPECIAL_NODEID ) break;
                    }
                }
            }
            assert( eParent != SPECIAL_NODEID );
            assert( eChild  != SPECIAL_NODEID );
            assert( eParent >= _graph->firstLevelEdge(middle) );
            assert( eParent <  _graph->lastEdge(middle) );
            assert( eChild  >= _graph->firstLevelEdge(middle) );
            assert( eChild  <  _graph->lastEdge(middle) );
            assert( _graph->edge(eParent).target() == parent );
            assert( _graph->edge(eChild).target() == t );
            assert( _graph->edge(eParent).weight() + _graph->edge(eChild).weight() == edge.weight() );

            // The first of the two edges is put on the stack ...
            _expandPathStack.push( make_pair( parent, eParent ) );
            // .. the other edge one is process in the next pass of the loop.
            processSecondHalf = true;
            parent = middle;
            e = eChild;
        }

        // Edge was no shortcut, can be added to the expanded path
        else
        {
            processSecondHalf = false;
            t = parent;
            path.add( t, edge.weight(), e, edge.isShortcut() );
        }
    }



    /**
    * Returns true iff s is the parent of u.
    * In the positive case, sets the given EdgeID e to the
    * ID of the edge (s,u).
    */
    bool directEdgeTo(const NodeID s, const NodeID u, EdgeID& e) const {
        const int searchID = 0;
        const NodeID index = isReached(searchID, u);
        assert( index != 0 );
        const PQData& data = pqData(searchID, index);
        if (data.parentNode() == s) {
            e = data.parentEdge();
            return true;
        }
        return false;
    }

    /** Returns the current upper bound. */
    EdgeWeight upperBound() const {return _upperBound;}

    /**
    * Initialises the upper bound to the given value.
    * Used when the main phase of the landmark query adopts the upper bound
    * from the initial phase.
    */
    void initUpperBound(EdgeWeight ub) {_upperBound = ub;}

    /** Clears the search in order to be able to start a new search. */
    void clear() {
        // clear both search directions
        clear(0);
        if (searchDirections == 2) clear(1);

        _upperBound = Weight::MAX_VALUE;
    }


    /**
    * Returns for a given search direction a reference to the vector
    * that contains the ids of all settled nodes.
    */
    vector<NodeID>& settledNodes(int searchID) { return _settledNodes[searchID]; }

    /**
     * Adds all edges of the search space of the given direction to the given vector.
     * Does not expand the edges.
     */
    void obtainSearchSpace(int searchID, vector<SearchSpaceEdge>& searchSpace) const {
        for (NodeID i = 0; i < _settledNodes[searchID].size(); i++) {
            const NodeID v = _settledNodes[searchID][i];
            const NodeID index = _graph->node(v).pqElement();
            const EdgeWeight k = pqKey(searchID, index);
            const PQData& data = pqData(searchID, index);
            if (data.isStartNode()) continue;

            searchSpace.push_back(
                SearchSpaceEdge(data.parentNode(), v, data.parentEdge(), k) );
        }
    }

    /**
    * Used for many-to-many computations.
    * Writes the search space of the last forward/backward search
    * to the given object.
    * In case of highway hierarchies, the complete search space is
    * considered as relevant.
    * In case of highway-node routing, only non-stalled nodes are relevant.
    */
    void obtainRelevantSearchSpace(SearchSpaces& ss) const {
        assert( (manyToManyMode == MTMM_FW) || (manyToManyMode == MTMM_BW) );
        const int searchID = (manyToManyMode == MTMM_FW) ? 0 : 1;
        for (NodeID i = 0; i < _settledNodes[searchID].size(); i++) {
            const NodeID u = _settledNodes[searchID][i];
            const NodeID index = _graph->node(u).pqElement();
            if (pqData(searchID, index).stalled()) continue;
            const EdgeWeight dist = pqKey(searchID, index);
            if (manyToManyMode == MTMM_FW) ss.addToSearchSpaceFW(u, dist);
            if (manyToManyMode == MTMM_BW) ss.addToSearchSpaceBW<false, false>(u, dist);
        }
    }


    /**
    * Returns the index of the corresponding pq element
    * if the specified node is reached; returns 0, otherwise.
    */
    NodeID isReached(int searchID, NodeID vID) const {
        return isReached(searchID, _graph->node(vID));
    }

    /**
    * Returns the index of the corresponding pq element
    * if the specified node is reached; returns 0, otherwise.
    */
    NodeID isReached(int searchID, const MyNode& v) const {
        NodeID index = v.pqElement();
        if (index == 0) return 0;
        if (pqueue(searchID).isDummy(index)) return 0;
        return index;
    }

    /**
    * Returns the index of the corresponding pq element
    * if the specified node is settled; returns 0, otherwise.
    */
    NodeID isSettled(int searchID, NodeID vID) const {
        NodeID index = _graph->node(vID).pqElement();
        if (index == 0) return 0;
        if (pqueue(searchID).isDummy(index)) return 0;

        // node IS reached ! Is it settled ?
        if (pqueue(searchID).elements()[index].hasBeenDeleted()) return index;
        return 0;
    }

    /** Returns the key of the specified pq element. */
    EdgeWeight pqKey(int searchID, NodeID index) const {
        assert( (searchID >= 0) && (searchID < searchDirections) );
        assert( (index > 0) && (index < pqueue(searchID).elements().size()) );
        assert( ! pqueue(searchID).isDummy(index) );

        return pqueue(searchID).elements()[index].key();
    }

    /** Returns a reference to the data object of the specified pq element. */
    PQData& pqData(int searchID, NodeID index) {
        assert( (searchID >= 0) && (searchID < searchDirections) );
        assert( (index > 0) && (index < pqueue(searchID).elements().size()) );
        assert( ! pqueue(searchID).isDummy(index) );

        return pqueue(searchID).elements()[index].data();
    }

    /** Returns a reference to the data object of the specified pq element. */
    const PQData& pqData(int searchID, NodeID index) const {
        assert( (searchID >= 0) && (searchID < searchDirections) );
        assert( (index > 0) && (index < pqueue(searchID).elements().size()) );
        assert( ! pqueue(searchID).isDummy(index) );

        return pqueue(searchID).elements()[index].data();
    }

    /**
    * Returns the number of elements in the priority queue.
    * Note that if two pqueues are used, both have the same number of elements.
    */

    NodeID noOfPQElements() const {
        return pqueue(0).elements().size();
    }

    /**
    * Returns the number of settled nodes in all search directions.
    */
    NodeID noOfSettledNodes() const {
        if ( searchDirections == 2 )
        {
            return _settledNodes[0].size() + _settledNodes[1].size();
        }
        else
        {
            return _settledNodes[0].size();
        }
    }

    /**
    * Returns the number of settled nodes that are not settled in all search directions.
    */
    NodeID noOfSettledMinusStalledNodes() const {
        NodeID result = 0;
        for ( int searchID = 0; searchID < searchDirections; searchID ++ )
        {
            for ( vector<NodeID>::const_iterator iter =  _settledNodes[searchID].begin();
                  iter != _settledNodes[searchID].end(); iter++ )
            {
                if ( !pqData( searchID, isReached( searchID, *iter ) ).stalled() ) result++;
            }
        }
        return result;
    }

    /**
    * Initialises the pqueue(s) so that they contain the given number of dummy elements.
    * This is useful when a search is performed after another search has not cleaned up.
    * In particular, this is used when the main phase of the landmark query is executed
    * so that the shortest path tree of the initial phase can be kept in order to be able
    * to compose the complete path at the end.
    */
    void initDummies(NodeID number) {
        VERBOSE2( cerr << "init dummies " << number << endl );
        pqueue(0).initDummies(number);
        if (searchDirections == 2) pqueue(1).initDummies(number);
    }


    /**
    * Insert a node with given distance and parent into the priority queue.
    * Used to update shortest paths based Voronoi regions.
    */
    void insertNode( const NodeID node, const EdgeWeight dist, const NodeID parent )
    {
        NodeID index = insert( 0, dist, node );
        pqData(0, index).updateParent(parent, SPECIAL_NODEID);
    }

    /**
     * Insert or update a node with given distance and parent into the priority queue.
     * Used to update shortest paths based Voronoi regions.
     */
    void updateNode( const NodeID node, const EdgeWeight dist, const NodeID parent )
    {
        NodeID e = isReached(0, node);
        if (! e) {
            // v has been unreached -> INSERT it
            e = insert( 0, dist, node );
        }
        else {
            // v has been reached -> try to DECREASE its KEY
            e = decreaseKey( 0, dist, node );
        }
        if (e) { // if the pqueue operation/the relaxation of the edge was successful
            COUNTING( counter.incDouble(COUNT_RELAXED_EDGES_SUCCESS) );

            PQData& data = pqData(0, e);
            // update the parent of the node v
            data.updateParent( parent, SPECIAL_NODEID );
        }
    }

    /**
     * Perform a single step of Dijkstra's algorithm. The topmost node in the priority
     * queue is deleted and returned (call-by-reference).
     * Used to update shortest paths based Voronoi regions.
     */
    bool searchNext(NodeID& node, EdgeWeight& dist, NodeID& parent ) {
        CurrentNode u;
        if (! pqueue(0).empty() ) {

            // *** deleteMin ***
            deleteMin(0, u);

            node = u.nodeID;
            dist = u.dist;
            parent = pqData(0, u.index).parentNode();

            return true;
        }
        return false;
    }


private:
    /**
    * Encapsulates data about the current node, i.e.,
    * the node that has just been deleted from the pqueue.
    */
    class CurrentNode
    {
        public:
            CurrentNode() {}
            CurrentNode(NodeID nID, EdgeWeight w) : nodeID(nID), dist(w) {}

            /** The id of this node. */
            NodeID nodeID;
            /** The id of the edge that is being relaxed. */
            EdgeID edgeID;
            /** The distance from the source node of the applicable search to this node. */
            EdgeWeight dist;
            /** The index of the element in the pqueue that represents this node. */
            NodeID index;
            /** The number of hops from the source node. */
            NodeID hops;
    };


    /** The graph. */
    Graph *const _graph;

    /** One priority queue for each search direction. */
    PQueue _pq[searchDirections];

    /** A list of the settled nodes (for each search direction). */
    vector<NodeID> _settledNodes[2];

    /**
    * An upper bound on the shortest path length.
    * Used during a bidirectional search.
    * Might be used to prune the search.
    * If the search is successful, equals to the actual shortest path distance.
    * If the search fails (not connected), equals to infinity (Weight::MAX_VALUE).
    */
    EdgeWeight _upperBound;

    /**
    * The intermediate node where the (tentative) shortest path
    * leaves one search scope and enters the other.
    * Used during a bidirectional search.
    */
    NodeID _via;

    /**
    * Indicates that both search scopes have met.
    * Now, we have to wait until the advanced abort-on-success criterion is fulfilled.
    */
    bool _concludeHwySearch;

    /** For each search direction, indicates that the search is finished. */
    bool _finished[2];


    /** Used for the stalling BFS (DM_QUERY). */
    queue< pair<NodeID, EdgeWeight> > _stallQueue;

    /** Used for path expansion. */
    stack< pair<NodeID,EdgeID> > _expandPathStack;

    /** Returns a reference to the specified pqueue. */
    PQueue& pqueue(int searchID) {
        assert( (searchID >= 0) && (searchID < searchDirections) );
        return _pq[searchID];
    }

    /** Returns a reference to the specified pqueue. */
    const PQueue& pqueue(int searchID) const {
        assert( (searchID >= 0) && (searchID < searchDirections) );
        return _pq[searchID];
    }

    /** Inserts the given node as start node of the specified search direction. */
    void insertStartNode(int searchID, NodeID nodeID) {
        NodeID index = insert(searchID, 0, nodeID);
        PQData& data = pqData(searchID, index);
        data.setStartNode();
    }

    /**
    * Inserts the given node with the given distance from the source node
    * into the pqueue of the specified search direction.
    * @return the index of the pq element that represents the inserted node
    */
    NodeID insert(int searchID, EdgeWeight dist, NodeID nodeID) {
        NodeID index = _graph->node(nodeID).pqElement();
        if (index == 0) {
            index = pqueue(searchID).insert(dist);
            _graph->node(nodeID).pqElement(index);
            if (searchDirections == 2) pqueue(1-searchID).insertDummy();
        }
        else {
            pqueue(searchID).insert(dist, index);
        }
        PQData& data = pqData(searchID, index);
        data.init(nodeID);
        if ( localSearch )
        {
            data.setTarget(_graph->node(nodeID).isTarget());
        }
        return index;
    }

    /**
    * Deletes the minimum element from the pqueue of the specified search direction.
    * @param searchID the search direction (forwards/backwards)
    * @param u a reference to the CurrentNode object that accepts the minimum element:
    *          sets index, nodeID, and dist
    * @return true iff 'successImpliesAbort' is turned on and both search scopes have met
    */
    bool deleteMin(int searchID, CurrentNode& u) {
        u.index = pqueue(searchID).deleteMin();
        return deleteElement( searchID, u );
    }

    /**
    * Deletes an arbitrary element from the specified pqueue.
    * @see deleteElement(...)
    * @return true iff 'successImpliesAbort' is turned on and both search scopes have met
    */
    bool deleteArbitrary(int searchID) {
        CurrentNode u;
        u.index = pqueue(searchID).deleteArbitrary();
        bool result = deleteElement( searchID, u );

        return result;
    }

    /**
    * Used as a subroutine by deleteMin and deleteArbitrary.
    * Performs some necessary steps directly after a node has been
    * deleted from a pqueue: the node is added to the list of settled nodes
    * and (if applicable) it is checked whether both search scopes have met.
    * @see deleteMin(...)
    * @see deleteArbitrary(...)
    * @param searchID the search direction (forwards/backwards)
    * @param u a reference to the CurrentNode object that contains the minimum element:
    *          expects index; sets nodeID and dist
    * @return true iff 'successImpliesAbort' is turned on and both search scopes have met
    */
    bool deleteElement(int searchID, CurrentNode& u) {
        // set nodeID and dist of the current node
        const PQData& data = pqData(searchID, u.index);
        u.nodeID = data.nodeID();
        u.dist = pqKey(searchID, u.index);

        // used for local search with hop-limit
        if (countHops)
        {
            u.hops = data.hops();
        }

        // add current node to the list of settled nodes
        _settledNodes[searchID].push_back(u.nodeID);

        if (searchDirections == 2) {
            // bidirectional search
            int opposite = 1 - searchID;
            // check if both search scopes have met
            NodeID index = isSettled(opposite, u.nodeID);
            if (index) {
                // The current node has been settled from the other direction as well.
                // A new path has been found.
                EdgeWeight newDist = u.dist + pqKey(opposite, index);
                // Check if the new path is shorter than the previously best one.
                if (newDist < _upperBound) {

                    VERBOSE2( cerr << searchID << " " << (_settledNodes[0].size() + _settledNodes[1].size())
                            << ": path found: ub = " << newDist << endl );
                    _upperBound = newDist;
                    _via = u.nodeID;
                }
                return true;
            }
        }
        return false;
    }

    /**
    * Decreases the key of a given node if applicable.
    * @param searchID the search direction (forwards/backwards)
    * @param newDist the new distance (key)
    * @param v the node whose key should be decreased
    * @return the index of the element whose key has been decreased
    *         or 0 if the key is not decreased.
    */
    NodeID decreaseKey(int searchID, EdgeWeight newDist, NodeID v)
    {
        // retrieve the element that represents the given node v
        const NodeID index = _graph->node(v).pqElement();
        const EdgeWeight key = pqKey(searchID, index);

        // no improvement, no decreaseKey
        if (newDist >= key) return 0;

        pqueue(searchID).decreaseKey(index, newDist);
        return index; // indicate that decreaseKey operation has been preformed
    }

    /**
    * Relaxes all (relevant) edges of the given node for the given search direction.
    */
    void relaxEdges(int searchID, CurrentNode& parent, NodeID ignore = SPECIAL_NODEID) {
        // No reference ! The memory of the underlying vector can be reallocated !
        const PQData parentData = pqData(searchID, parent.index);

        EdgeWeight parentDist = parent.dist;

        // stall-on-demand: Do not relax edges of a stalled node
        // since it is reached on a suboptimal (not shortest) path
        if ( stallOnDemand )
        {
            if ( parentData.stalled() ) return;
        }

        // Asymmetric many-to-many search. Edges into the core (topmost level)
        // are not relaxed on backward search (MTMM_BW) because the forward search
        // will relax into the core. This is correct since the core is an overlay
        // graph.
        if ( manyToManyMode == MTMM_BW )
        {
            // do not relax edges of a core node
            if ( _graph->node(parent.nodeID).isInCore() ) return;
        }

        // first and last edge of the current node in the current level
        EdgeID firstEdge = _graph->firstLevelEdge(parent.nodeID);
        EdgeID lastEdge = _graph->lastEdge(parent.nodeID);

        for (parent.edgeID = firstEdge; parent.edgeID < lastEdge; parent.edgeID++) {
            // note: parent.edgeID contains the id of the edge that should be relaxed
            // *** relaxEdge(searchID, parent, edgeLevel); ***

            // reference to the edge that should be relaxed
            Edge& edge = _graph->edge(parent.edgeID);

            NodeID vID = edge.target(); // the id of the target of the edge
            MyNode& v = _graph->node(vID); // reference to the target v of the edge
            NodeID index = 0;

            if (stallOnDemand)
            {
                // try to wake up a node v that can start a stalling process
                index = isReached(searchID, v);
                // Node v has to be the endpoint of an edge (u,v) that
                // points into the opposite direction
                // (because then the edge (v,u) which is used in the stalling process points into
                //  the right direction).
                // Furthermore, v has to be reached.
                if (edge.isDirected(1-searchID) && index) {
                    const PQData& data = pqData(searchID, index);
                    const bool stalled = data.stalled();
                    const EdgeWeight vKey = stalled ? data.stallKey() : pqKey(searchID, index);
                    const EdgeWeight newKey = vKey + edge.weight();
                    if (newKey < parentDist) { // check whether this node is 'promising'
                        stall(searchID, parent.nodeID, newKey);
                        assert( pqData(searchID, parent.index).stalled() );
                        return;
                    }
                }
            }

            // local search: if the current edge points to the node to be ignored, ignore it
            if ( localSearch )
            {
                if (vID == ignore) continue;
            }

            // if the current edge points to the wrong direction, ignore it
            if (! edge.isDirected(searchID) ) continue;

            // the new distance from the source of the search to the target of the edge
            EdgeWeight newDist = parentDist + edge.weight();

            COUNTING( counter.incDouble(COUNT_RELAXED_EDGES) );

            // perform the pqueue operation (insert or decreaseKey)
            NodeID e;
            index = isReached(searchID, v);
            if (! index) {
                // v has been unreached -> INSERT it
                e = insert( searchID, newDist, vID );
            }
            else {
                // v has been reached -> try to DECREASE its KEY
                e = decreaseKey( searchID, newDist, vID );
            }
            if (e) { // if the pqueue operation/the relaxation of the edge was successful
                COUNTING( counter.incDouble(COUNT_RELAXED_EDGES_SUCCESS) );
                PQData& data = pqData(searchID, e);
                // update the parent of the node v
                data.updateParent( parent.nodeID, parent.edgeID );

                // reached from a new parent -> stalling of this node
                // (if it has occured) is no longer valid!
                if (stallOnDemand)
                {
                    data.unstall();
                }

                // count hops of the current paths from the source node
                // used for local search with hop-limits
                if (countHops)
                {
                    data.setHops(parent.hops+1);
                }
            }
        }
    }

    /**
    * Processes all reached, but unsettled nodes for the given search direction.
    * Invokes deleteArbitrary.
    * @see deleteArbitrary(...)
    */
    void processUnsettledNodes(int searchID) {
        while (! pqueue(searchID).empty()) {
            deleteArbitrary(searchID);
        }
    }

    /**
    * Clears a given node vector.
    * Deletes all corresponding pqueue elements.
    */
    void clearNodeVector(int searchID, vector<NodeID>& nodeVector) {
        for (NodeID i = 0; i < nodeVector.size(); i++) {
            _graph->node( nodeVector[i] ).pqElement(0);
        }
        nodeVector.clear();
    }

    /**
     * Clears the given search direction. Especially pqElement entry
     * in the node data structure of the graph is cleared so the next
     * Dijkstra search can use it and decide wheter a node isReached().
     */
    void clear(int searchID) {
        assert( (searchID >= 0) && (searchID < searchDirections) );

        processUnsettledNodes(searchID);

        clearNodeVector(searchID, _settledNodes[searchID]);

        pqueue(searchID).clear();
    }


    /** Returns true iff the pqElement entries of all nodes have been reset. */
    bool checkClean() const {
        return true; // deactivate this check
        // This takes a lot of time.
        for (NodeID u = 0; u < _graph->noOfNodes(); u++) {
            if (_graph->node(u).pqElement() != 0) return false;
        }
        return true;
    }


    /**
    * Performs a stalling process, implemented by a BFS.
    * @param searchID the search direction (the direction of the stalling process
    *                is equal to the search direction)
    * @param u the 'stalling node' that starts the stalling process
    * @param key the key of the stalling node
    */
    void stall(const int searchID, NodeID u, EdgeWeight key) {
        COUNTING( counter.incDouble(COUNT_STALL_OPS) );
        NodeID index = isReached(searchID, u);
        assert( index > 0 );
        PQData& data = pqData(searchID, index);
        data.stallKey(key);
        _stallQueue.push(pair<NodeID, EdgeWeight>(u, key));
        if (deepStallOnDemand) stallParent(searchID, u, data, key);

        while (! _stallQueue.empty()) {
            COUNTING( counter.incDouble(COUNT_STALL_STEPS) );

            u = _stallQueue.front().first;
            key = _stallQueue.front().second;
            _stallQueue.pop();

            const EdgeID lastEdge = _graph->lastEdge(u);
            for (EdgeID e = _graph->firstLevelEdge(u); e < lastEdge; e++) {
                const Edge& edge = _graph->edge(e);
                const NodeID v = edge.target();

                // if the current edge points to the wrong direction, ignore it
                if (! edge.isDirected(searchID)) continue;

                NodeID index = isReached(searchID, v);
                if (index) {
                    const EdgeWeight newKey = key + edge.weight();
                    if (newKey < pqKey(searchID, index)) { // shorter path found?
                        PQData& data = pqData(searchID, index);

                        if (! data.stalled()) {
                            data.stallKey(newKey);
                            _stallQueue.push(pair<NodeID, EdgeWeight>(v, newKey));
                            if (deepStallOnDemand) stallParent(searchID, v, data, newKey);
                        }
                    }
                }
            }
        }
    }

    /**
    * Continue stalling process (see stall()) using the parent pointers.
    * This can increase the number of stalled nodes if search graphs are used that
    * store an edge only at the incident node with the smaller level.
    * This is called "deep stall-on-demand".
    */
    void stallParent(const int searchID, const NodeID u, const PQData& data, const EdgeWeight key)
    {
        NodeID index = isReached(searchID, data.parentNode());
        assert( index > 0 );
        PQData& parentData = pqData(searchID, index);
        if ( parentData.stalled() ) return;
        const EdgeID lastEdge = _graph->lastEdge(data.parentNode());
        for (EdgeID e = _graph->firstLevelEdge(data.parentNode()); e < lastEdge; e++) {
            const Edge& edge = _graph->edge(e);

            // if the current edge points to the wrong direction, ignore it
            // we need a edge from u to the parent => other search direction (from the view of the parent)
            if ( ! ( edge.isDirected(1-searchID)  && edge.target() == u) ) continue;

            const EdgeWeight newKey = key + edge.weight();
            if (newKey < pqKey(searchID, index)) { // shorter path found?
                parentData.stallKey(newKey);
                _stallQueue.push(pair<NodeID, EdgeWeight>(data.parentNode(), newKey));

            }
        }
    }

    };

} // namespace

/** Used for the many-to-many forward search. */
typedef processing::DijkstraCH<datastr::graph::SearchGraph, DynQueryPQueue, 2, true, true,
processing::MTMM_FW, false, false> DijkstraCHManyToManyFW;

/** Used for the many-to-many backward search. */
typedef processing::DijkstraCH<datastr::graph::SearchGraph, DynQueryPQueue, 2, true, true,
processing::MTMM_BW, false, false> DijkstraCHManyToManyBW;

/** Used for local searches during the contraction. */
typedef processing::DijkstraCH<datastr::graph::UpdateableGraph, LocalSearchPQueue,
    1, false, false, processing::MTMM_NONE, true, true > LocalDijkstraContract;

/** Used for upates of Voronoi regions (priority term) during node ordering */
typedef processing::DijkstraCH<datastr::graph::UpdateableGraph, NormalPQueue, 1, false, false,
    processing::MTMM_NONE, false, false > DijkstraUpdateVoronoi;

/** Used for bidirectional fast query in contraction hierarchies. */
typedef processing::DijkstraCH<datastr::graph::SearchGraph, DynQueryPQueue,
    2, true /*stall-on-demand*/, false/*deep stall-on-demand*/> DijkstraSearchCH;

/** Normal bidirectional query. Used for testing and debugging. */
typedef processing::DijkstraCH<datastr::graph::UpdateableGraph, NormalPQueue, 2, false> DijkstraSearchBidir;

#endif // _PROCESSING_DIJKSTRACH_H
