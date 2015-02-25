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


#ifndef _DATASTR_GRAPH_UPDATEABLEGRAPH_H
#define _DATASTR_GRAPH_UPDATEABLEGRAPH_H

/**
 * Represents a graph that allows changes (e.g., edge additions/removals)
 * at any point in time. Its interface is closly related to the Static/DynamicGraph.
 * Internally, an adjacency array is kept. The size of each edge group is always
 * a power of two. If there are less edges leaving a particular node, the data
 * structure contains holes. If an edge is added such that the existing capacity
 * is exceeded, the entire edge group is moved to the end of the edge array (which
 * is resized) and the capacity is doubled. Note that the capacity is not stored
 * explicitly (in order to save memory), but usually the smallest power of two 2^k
 * >= 'the number of edges' is assumed as capacity. We can mark that a node has an
 * 'increased capacity' of 2^[k+1}. This is useful when an edge is removed such that
 * the number of edges is reduced from 2^k + 1 to 2^k. Then, the old capacity of 2^[k+1}
 * is still available and not only 2^k. See also 'removeEdge'.
 * This class originated from updGraph.h and was modified for contraction hierarchies (CH).
 *
 * The adjacency array for each node is organized as follows:
 *      firstEdge
 *          ..
 *        (edges having their other incident node in a lower level)
 *          ..
 *      firstLevelEdge
 *          ..
 *        (edges having their other incident node in the same or higher level)
 *          ..
 *      lastEdge
 *
 * Note: This organization is the main difference to the graph datastructure
 *       used by the HH/HNR/TNR code of Dominik Schultes.
 *       The disadvantage is that direct access to the edges leading to a specific
 *       level are impossible. The advantage is that the organization is much simpler,
 *       allows a distinct level for each node and is sufficient for constructing
 *       contraction hierarchies (CH).
 */
namespace datastr { namespace graph {
class UpdateableGraph
{
public:
    /** Represents a Node in the UpdateableGraph. */
    class UpdNode
    {
    public:
        /** Default Constructor. */
        UpdNode() : _capacity(0), _target(false),
            _level(0), _firstEdge(0), _firstLevelEdgeOffset(0), _lastEdgeOffset(0), _pqElement(0) {}

        /**
       * Constructor.
       * @param fE index of the first edge leaving this node
       * @param k level of this node
       */
        UpdNode(const EdgeID fE, const LevelID k) {
            assert( (k >= 0) );
            _firstEdge = fE;
            _level = k;
            _firstLevelEdgeOffset = 0;
            _lastEdgeOffset = 0;
            _pqElement = 0; // null pointer means "not in pqueue"
        }

        /** Returns the index of the first edge leaving this node. */
        EdgeID firstEdge() const { return _firstEdge; }
        /** Sets the index of the first edge leaving this node. */
        void setFirstEdge(const EdgeID fE) {
            _firstEdge = fE;
        }

        /**
         * Returns the index of the first edge leaving this node having
         * the other incident node in the same or higher level than the
         * current node.
         */
        EdgeID firstLevelEdge() const { return _firstEdge + _firstLevelEdgeOffset ; }
        /**
         * Sets the index of the first edge leaving this node having
         * the other incident node in the same or higher level than the
         * current node.
         */
        void setFirstLevelEdge(const EdgeID fLE) {
            assert( fLE >= _firstEdge );
            _firstLevelEdgeOffset = fLE - _firstEdge;
        }

        /** Returns the index of the last edge leaving this node + 1. */
        EdgeID lastEdge() const { return _firstEdge + _lastEdgeOffset; }
        /** Sets the index of the last edge leaving this node + 1. */
        void setLastEdge(const EdgeID lE) {
            assert( lE >= _firstEdge + _firstLevelEdgeOffset );
            _lastEdgeOffset = lE - _firstEdge;
        }

        /** Returns the level of this nodes. */
        LevelID level() const { return _level; }
        /** Set the level of the node. A graph->changeNodeLevel(node) may be necessary. */
        void setLevel(const LevelID k) {
            _level = k;
        }

        /**
         * Returns true if the capacity of the node is not the smallest
         * power of two 2^k >= the number of edges, but 2^{k+1}.
         */
        bool increasedCapacity() const {return _capacity;}

        /** Marks that this node has an increased capacity. */
        void setIncreasedCapacity() {
            assert( ! increasedCapacity() );
            _capacity++;
        }

        /** Marks that this node does not have an increased capacity. */
        void unsetIncreasedCapacity() {
            assert( increasedCapacity() );
            _capacity--;
        }

        /** Returns wheter this node is a target in a local search. */
        bool isTarget() const { return _target; }
        /** Set the target flag for local searches. */
        void setTarget(bool v) { _target = v; }

        /** Returns the index of the corresponding elements in the pqueues. */
        NodeID pqElement() const {return _pqElement;}

        /** Sets the index of the corresponding elements in the pqueues. */
        void pqElement(const NodeID pqElement) {_pqElement = pqElement;}

        /** wheter this node is in the core, only relevant for SearchGraph */
        bool isInCore() const { assert(false); return false; }
        void setInCore(bool inCore) { assert(false); }

        // Serialization currently not supported.
        void serialize(ostream& out) const {
            assert(false);
        }

        // Serialization currently not supported.
        void deserialize(istream& in) {
            assert(false);
        }

    private:
        /**
         * Stores the 'increased capacity'-flag (bit 1 (= the LSB)),
         * the target-flag (bit 2) and the node level (bits 3-32),
         * and the first edge index (bits 6-32).
         */
        unsigned int _capacity:1;
        bool _target:1;
        LevelID _level:30;

        EdgeID _firstEdge;
        unsigned int _firstLevelEdgeOffset;
        unsigned int _lastEdgeOffset;
        NodeID _pqElement;
    };

    /** The node type used in this graph. */
    typedef UpdNode MyNode;


public:
    /** Constructor. Builds the graph. */
    UpdateableGraph(vector<CompleteEdge>& edges, const vector<LevelID>& nodeLevels) {
        construct(edges, nodeLevels);
    }

    /** Constructor. Deserializes the graph. */
    /*
    UpdateableGraph(istream& in) {
        deserialize(in);
    }
    */

    NodeID noOfNodes() const {return _nodes.size();}

    /** Returns the number of edges that memory has been allocated for. */
    EdgeID noOfEdges() const {return _edges.size();}

    /** Returns the index of the first edge of u. */
    EdgeID firstEdge(const NodeID u) const {
        return node(u).firstEdge();
    }

    /**
     * Returns the index of the first edge of u that leads to the
     * same or higher level.
     */
    EdgeID firstLevelEdge(const NodeID u) const {
        return node(u).firstLevelEdge();
    }

    /** Returns the index+1 of the last edge of u. */
    EdgeID lastEdge(const NodeID u) const {
        return node(u).lastEdge();
    }

    /**
    * Returns the number of edges stored at node u
    * irrespective of the direction.
    */
    EdgeID degree(const NodeID u) const {
        return lastEdge(u)-firstEdge(u);
    }

    /**
    * Returns the number of existing edges, i.e.,
    * noOfEdges() minus the number of non-used holes in
    * the edge array.
    * Note: this is a not a simple return of a variable
    */
    EdgeID noOfExistingEdges() const {
        EdgeID sum = 0;
        for (NodeID u = 0; u < noOfNodes(); u++) sum += degree(u);
        return sum;
    }

    const UpdNode& node(const NodeID u) const {
        assert( u < _nodes.size() );
        return _nodes[u];
    }

    UpdNode& node(const NodeID u) {
        assert( u < _nodes.size() );
        return _nodes[u];
    }

    const Edge& edge(const EdgeID e) const {
        assert( e < _edges.size() );
        return _edges[e];
    }

    Edge& edge(const EdgeID e) {
        assert( e < _edges.size() );
        return _edges[e];
    }

    /** Sort edges. */
    void sortEdges(EdgeID first, EdgeID last)
    {
        assert( first < last );
        assert( last <= _edges.size() );
        sort( _edges.begin() + first, _edges.begin() + last );
    }

    /**
    * The nodes are internally ordered by level.
    * Returns for a given original NodeID, the
    * internally used NodeID.
    */
    NodeID mapExtToIntNodeID(const NodeID ext) const {
        assert( ext < _mapExtToIntNodeIDs.size() );
        return _mapExtToIntNodeIDs[ext];
    }

    /**
     * The adjacency array of a node is partitioned into the edges leading
     * to lower levels and the edges leading to the same or higher levels.
     * This method ensures this partition if only the specified edge
     * may violate it.
     */
    void changeEdgeLevel(const NodeID u, const EdgeID e) {
        assert( (firstEdge(u) <= e) && (e < lastEdge(u)) );

        // edge e leads to a lower level and is in the wrong partiton
        if ( node(u).level() > node(edge(e).target()).level() && e >= firstLevelEdge(u) )
        {
            if (e > firstLevelEdge(u)) swap(_edges[e], _edges[firstLevelEdge(u)]);
            node(u).setFirstLevelEdge(firstLevelEdge(u)+1);
        }

        // edge e leads to the same or higher level and is in the wrong partiton
        else if ( node(u).level() <= node(edge(e).target()).level() && e < firstLevelEdge(u) )
        {
            if (e <  firstLevelEdge(u)-1) swap(_edges[e], _edges[firstLevelEdge(u)-1]);
            node(u).setFirstLevelEdge(firstLevelEdge(u)-1);
        }
    }

    /**
     * The adjacency array of a node is partitioned into the edges leading
     * to lower levels and the edges leading to the same or higher levels.
     * This method is called after the level of node u changed. It
     * ensures the correct partition of the edges if only edges incident
     * to node u may violate it.
     */
    void changeNodeLevel(const NodeID u)
    {
        // e is the index into the edge array leading to lower levels
        // f is the index into the edge array leading to higher levels
        // e starts at firstEdge and f at lastEdge-1 both move to the middle
        // and if a misplaced edge in both partitions is found, they are
        // swapped.
        EdgeID e = firstEdge(u);
        EdgeID f = lastEdge(u);

        if (e < f)
        {
            f--;
            // Work until index e and index f meet.
            while (e <= f)
            {
                // Find a misplaced edge in firstEdge .. firstLevelEdge-1
                while ( e <= f && node(u).level() > node(edge(e).target()).level() )
                {
                    // For each edge, check the reverse edge.
                    changeEdgeLevel(edge(e).target(),reverseEdge(u, e));
                    e++;
                }
                // Find a misplaced edge in firstLevelEdge .. lastEdge-1
                while ( e <= f && node(u).level() <= node(edge(f).target()).level() )
                {
                    // For each edge, check the reverse edge
                    changeEdgeLevel(edge(f).target(),reverseEdge(u, f));
                    if ( f == 0 ) break;
                    f--;
                }
                if ( f == 0 ) break;
                // If a misplaced edge in both partitions is found, swap them.
                if (e < f)
                {
                    swap(_edges[e], _edges[f]);
                    changeEdgeLevel(edge(e).target(),reverseEdge(u, e));
                    changeEdgeLevel(edge(f).target(),reverseEdge(u, f));
                    e++;
                    f--;
                }
            }

            // Update firstLevelEdge, the partition border.
            if (firstLevelEdge(u) != e) node(u).setFirstLevelEdge(e);
        }
        assert( firstEdge(u) <= firstLevelEdge(u) );
        assert( firstLevelEdge(u) <= lastEdge(u) );


        // Check the result of the above algorithm.
        #ifndef NDEBUG
        for ( EdgeID e = firstEdge(u); e < firstLevelEdge(u); e++ )
        {
            NodeID v = edge(e).target();
            assert( node(u).level() > node(v).level() );
            for ( EdgeID f = firstEdge(v); f < firstLevelEdge(v); f++ )
            {
                assert( node(v).level() > node(edge(f).target()).level() );
            }
            for ( EdgeID f = firstLevelEdge(v); f < lastEdge(v); f++ )
            {
                assert( node(v).level() <= node(edge(f).target()).level() );
            }
        }
        for ( EdgeID e = firstLevelEdge(u); e < lastEdge(u); e++ )
        {
            NodeID v = edge(e).target();
            assert( node(u).level() <= node(v).level() );
            for ( EdgeID f = firstEdge(v); f < firstLevelEdge(v); f++ )
            {
                assert( node(v).level() > node(edge(f).target()).level() );
            }
            for ( EdgeID f = firstLevelEdge(v); f < lastEdge(v); f++ )
            {
                assert( node(v).level() <= node(edge(f).target()).level() );
            }
        }
        #endif
    }

    /**
     * The adjacency array of a node is partitioned into the edges leading
     * to lower levels and the edges leading to the same or higher levels.
     * This method is called after the level of node u changed. It
     * ensures the correct partition of the edges if only the edges leading
     * to more important nodes incident to node u may violate it and the
     * partition is correct for the adjacency array of node u. So only the
     * reverse edges between firstLevelEdge .. lastEdge-1 need to be checked.
     * This special case occurs during node contraction. Not until after the
     * contraction of a node, its level in  the graph datastructure is updated.
     * Then, inductive, its own adjacency array is partitioned correct but the
     * adjacency array of its remaining neighbors need an update because from their
     * point of view the edge now leads to a lower level.
     */
    void changeNodeLevelOnlyReverseEdges(const NodeID u)
    {
        EdgeID lastEdge = this->lastEdge(u);
        for ( EdgeID e = firstLevelEdge(u); e < lastEdge; e++ )
        {
            changeEdgeLevel(edge(e).target(),reverseLevelEdge(u, e));
        }


        // Check the result of the above algorithm.
        #ifndef NDEBUG
        for ( EdgeID e = firstEdge(u); e < firstLevelEdge(u); e++ )
        {
            NodeID v = edge(e).target();
            assert( node(u).level() > node(v).level() );
            for ( EdgeID f = firstEdge(v); f < firstLevelEdge(v); f++ )
            {
                assert( node(v).level() > node(edge(f).target()).level() );
            }
            for ( EdgeID f = firstLevelEdge(v); f < this->lastEdge(v); f++ )
            {
                assert( node(v).level() <= node(edge(f).target()).level() );
            }
        }
        for ( EdgeID e = firstLevelEdge(u); e < this->lastEdge(u); e++ )
        {
            NodeID v = edge(e).target();
            assert( node(u).level() <= node(v).level() );
            for ( EdgeID f = firstEdge(v); f < firstLevelEdge(v); f++ )
            {
                assert( node(v).level() > node(edge(f).target()).level() );
            }
            for ( EdgeID f = firstLevelEdge(v); f < this->lastEdge(v); f++ )
            {
                assert( node(v).level() <= node(edge(f).target()).level() );
            }
        }
        #endif
    }

    /**
    * Adds the given edge leaving node u to the graph.
    * @return the number of places the edge group of u is moved by this operation;
    *   this is important when we have an id of an edge (u,v) before this operation
    *   is performed and want to use the correct edge id after this operation has
    *   been performed.
    */
    EdgeID addEdge(const NodeID u, const Edge& edge ) {
        const EdgeID size = degree(u);
        EdgeID capacity = edgeCapacity(u, size);
        const EdgeID oldFirstEdge = node(u).firstEdge();
        EdgeID newFirstEdge = oldFirstEdge;
        assert( size <= capacity );

        // ensure that there is room for the new edge
        if (size == capacity) {
            // Case 1: edge group is full!
            // -> move group to the end of the array, double the capacity
            assert( ! node(u).increasedCapacity() );

            // try to reuse edge space
            capacity <<= 1;

            newFirstEdge = _edges.size();
            _edges.insert(_edges.end(),
                          _edges.begin() + oldFirstEdge,
                          _edges.begin() + oldFirstEdge + size);
            appendClosedEdges(size, capacity); // fill the holes by clearly marked dummy edges
            node(u).setFirstEdge(newFirstEdge);
        }
        else {
            // Case 2: there is still room for more edges!
            if (node(u).increasedCapacity()) {
                assert( size <= (capacity >> 1) );
                // Let 2^k denote the smallest power of two >= size.
                // Since the node has an 'increased capacity', the capacity is 2^{k+1}.
                // If we add another edge and the new size exceeds 2^k, the node longer has
                // an 'increased capacity' (but a normal capacity).
                if (size == (capacity >> 1)) node(u).unsetIncreasedCapacity();
            }
        }

        // add the new edge
        assert( size < capacity );
        EdgeID e = newFirstEdge + size;
        assert( e < _edges.size() );
        _edges[e] = edge; // insert the new edge at the end of the edge group...
        node(u).setLastEdge(e+1);
        assert( firstEdge(u) <= e && e < lastEdge(u) );
        changeEdgeLevel(u, e); // ...and move it to the right level

        return (newFirstEdge - oldFirstEdge);
    }

    /**
    * Removes the specified edge (u,v) with ID e.
    * Leaves a new hole in the data structure. If the capacity
    * was 2^[k+1} and the size is now reduced from 2^k + 1 to 2^k,
    * we mark that the node u has an increased capacity of 2^[k+1}.
    * WARNING: Using only one bit to mark an increased capacity,
    * we cannot mark the case that the size is reduced to 2^{k-1}
    * (or even less), although there is space for 2^[k+1} edges.
    */
    void removeEdge(const NodeID u, const EdgeID e) {
        edge(e).makeClosed(); // mark as 'hole'

        if ( e < firstLevelEdge(u) )
        {
            if ( e < firstLevelEdge(u) - 1) swap( _edges[e], _edges[firstLevelEdge(u)-1] );
            if (firstLevelEdge(u) < lastEdge(u)) swap( _edges[firstLevelEdge(u)-1], _edges[lastEdge(u)-1] );
            node(u).setFirstLevelEdge(firstLevelEdge(u)-1);
        }
        else
        {
            if ( e < lastEdge(u)-1 ) swap( _edges[e], _edges[lastEdge(u)-1] );
        }
        node(u).setLastEdge(lastEdge(u)-1);

        // if applicable, set 'increased capacity'
        EdgeID size = degree(u);
        if (size == edgeCapacity(u, size)) {
            node(u).setIncreasedCapacity();
        }
    }


    /**
     * Add a shortcut edge at source and target and keeps invariant
     * invariant: for each pair of nodes (u,v) there is at most one edge (u,v)
     * Always the shortest edge is kept.
     * This method is used during the construction of a contraction hierarchy
     * to add shortcut edges that represent existing paths in the graph.
     * @return edge difference
     */
    int addShortcutEdge(NodeID u, const Edge& newEdge)
    {
        return doAddShortcutEdge<false>(u, newEdge);
    }
    /**
     * Same as addShortcutEdge, but only simulate the operation to calculate
     * the edge difference. This is used during node ordering as priority.
     */
    int addShortcutEdgeSimulate(NodeID u, const Edge& newEdge)
    {
        return doAddShortcutEdge<true>(u, newEdge);
    }


    const static int NOTHING = 0;
    const static int ONE_EDGE = 1;
    const static int TWO_EDGES = 2;
    const static int LOOK_FOR_SECOND_EDGE_FORWARD = -1;
    const static int LOOK_FOR_SECOND_EDGE_BACKWARD = -2;

    /**
     * Add a shortcut edge at source and target and keeps invariant
     * invariant: for each pair of nodes (u,v) there is at most one edge (u,v)
     * This method is used during the construction of a contraction hierarchy
     * to add shortcut edges that represent existing paths in the graph.
     * Always the shortest edge is kept.
     * Let newEdge = (u,v), the algorithm scans through the adjacency
     * array of u and looks for an edge (u,v) (and (v,u) with bidirectional
     * flag set).
     * For performance reasons, we only allow shortcut edges having their
     * source and target node in the same level.
     *
     * @param simulateOnly only simulated, return edge difference
     * @return edge difference, can be
     *        -2 two unidirectional edges replaced by a bidirectional
     *         0 could use the position of an existing edge
     *         2 new edge necessary
     *
     */
    template < bool simulateOnly >
    int doAddShortcutEdge(NodeID u, const Edge& newEdge)
    {
        // Flags indicating the directions for that we need a new edge.
        // If we found an existing edge with <= weight or we can used
        // the position of an existing edge, we unset these variables
        // accordingly.
        bool forward = true;
        bool backward = newEdge.isBidirected();

        // Current state of the search for a parallel edge.
        // NOTHING = no edge found
        // ONE_EDGE = one unidirectional edge found
        // TWO_EDGES = two unidirectional or one bidirectional edge found
        // LOOK_FOR_SECOND_EDGE_FORWARD = old edge is unidirectional backward, new edge is
        //             bidirectional. Need to check wheter a second unidirectional forward
        //             edge exists.
        // LOOK_FOR_SECOND_EDGE_BACKWARD = old edge is unidirectional forward, new edge is
        //             bidirectional. Need to check wheter a second unidirectional backward
        //             edge exists.
        int state = NOTHING;

        // If we replaced an existing unidirectional edge by a bidirectional edge
        // (state LOOK_FOR_SECOND_EDGE_FORWARD or LOOK_FOR_SECOND_EDGE_BACKWARD)
        // we need to look for a possible second edge. If we found such an edge,
        // and this edge is shorter than our new edge, we need to make our
        // new bidirectional edge one-way. These to variables store the indices
        // to these edges to perform makeOneWay().
        EdgeID changedEdgeID = SPECIAL_NODEID;
        EdgeID changedReverseEdgeID = SPECIAL_NODEID;

        // Only allow shortcut edges between nodes in the same level
        assert( node(u).level() == node(newEdge.target()).level() );
        EdgeID edgeID = firstLevelEdge(u);
        EdgeID lastEdge = this->lastEdge(u);

        int result = 0;

        // Scan through the adjacency array of
        for ( ; edgeID < lastEdge && (forward || backward) && state < TWO_EDGES; edgeID++)
        {
            Edge& edge = this->edge(edgeID);
            assert( !edge.isClosed() );
            if (edge.target() == newEdge.target())
            {
                // Found a matching edge.
                state++;
                if (edge.isBidirected()) state++;
                assert( state <= TWO_EDGES );

                // Existing edge with <= weight can be used, shortcut not necessary.
                if (edge.weight() <= newEdge.weight())
                {
                    forward = forward && !edge.isDirected(0);
                    backward = backward && !edge.isDirected(1);
                }

                // new edge has smaller weight
                else {
                    if (edge.isBidirected())
                    {
                        EdgeID reverseEdgeID = reverseLevelEdge(u, edgeID);
                        Edge& reverseEdge = this->edge(reverseEdgeID);
                        if (forward && backward)
                        {
                            // old and new edge bidirectional
                            if ( !simulateOnly )
                            {
                                edge.setWeight(newEdge.weight());
                                reverseEdge.setWeight(newEdge.weight());
                                if ( edge.type() != newEdge.type() )
                                {
                                    edge.setType(newEdge.type());
                                    reverseEdge.setType(newEdge.type());
                                }
                                edge.copyShortcutInfo(newEdge, false);
                                reverseEdge.copyShortcutInfo( newEdge, true);
                            }

                            // replace old edge by new edge
                            forward = false;
                            backward = false;
                        }
                        else if (forward)
                        {
                            // old edge bidirectional, new edge only forward
                            // => old edge becomes backward
                            if ( !simulateOnly )
                            {
                                edge.makeOneWay(1);
                                reverseEdge.makeOneWay(0);
                            }
                        }
                        else if (backward)
                        {
                            // old edge bidirectional, new edge only backward
                            // => old edge becomes forward
                            if ( !simulateOnly )
                            {
                                edge.makeOneWay(0);
                                reverseEdge.makeOneWay(1);
                            }
                        }
                    }

                    // old edge unidirectional
                    else
                    {
                        if ((edge.isDirected(0) && forward)
                            || (edge.isDirected(1) && backward))
                        {
                            // new edge covers all directions that the old edge covers
                            EdgeID reverseEdgeID = reverseLevelEdge(u, edgeID);
                            Edge& reverseEdge = this->edge(reverseEdgeID);

                            if ( !simulateOnly )
                            {
                                edge.setWeight(newEdge.weight());
                                reverseEdge.setWeight(newEdge.weight());
                                if ( edge.type() != newEdge.type() )
                                {
                                    edge.setType(newEdge.type());
                                    reverseEdge.setType(newEdge.type());
                                }
                                edge.copyShortcutInfo(newEdge, false);
                                reverseEdge.copyShortcutInfo( newEdge, true);
                            }

                            if (forward && backward)
                            {
                                // old edge is not bidirectional, new edge is bidirectional
                                if (state == ONE_EDGE)
                                {
                                    // There could be another old edge in the opposite
                                    // direction. We need to look for it.
                                    state = edge.isDirected(0)
                                        ? LOOK_FOR_SECOND_EDGE_BACKWARD
                                        : LOOK_FOR_SECOND_EDGE_FORWARD;
                                    changedEdgeID = edgeID;
                                    changedReverseEdgeID = reverseEdgeID;
                                }
                                if ( !simulateOnly )
                                {
                                    edge.makeTwoWay();
                                    reverseEdge.makeTwoWay();
                                }
                            }

                            // replace old edge by new edge, no need to add a new edge
                            forward = false;
                            backward = false;
                        }
                    }
                }
            }
        }

        // If a unidirectional edge was replaced by a bidirectional one
        // check for unidirectional edge in other direction.
        if ( state == LOOK_FOR_SECOND_EDGE_FORWARD || state == LOOK_FOR_SECOND_EDGE_BACKWARD )
        {
            for ( ; edgeID < lastEdge; edgeID++)
            {
                Edge& edge = this->edge(edgeID);
                assert( !edge.isClosed() );
                if (edge.target() == newEdge.target())
                {
                    assert( state != LOOK_FOR_SECOND_EDGE_FORWARD || ( edge.isDirected(0) && !edge.isDirected(1) ) );
                    assert( state != LOOK_FOR_SECOND_EDGE_BACKWARD || ( edge.isDirected(1) && !edge.isDirected(0) ) );

                    // this one is shorter than our new edge, reduce our edge to one-way (unidirectional)
                    if (edge.weight() < newEdge.weight())
                    {
                        if ( !simulateOnly )
                        {
                            Edge& changedEdge = this->edge(changedEdgeID);
                            Edge& changedReverseEdge = this->edge(changedReverseEdgeID);
                            if ( state == LOOK_FOR_SECOND_EDGE_FORWARD )
                            {
                                // changed edge only backwards
                                changedEdge.makeOneWay(1);
                                changedReverseEdge.makeOneWay(0);
                            }
                            else
                            {
                                // changed edge only forwards
                                changedEdge.makeOneWay(0);
                                changedReverseEdge.makeOneWay(1);
                            }
                        }
                    }
                    // our new edge covers this one, remove it
                    else
                    {
                        if ( !simulateOnly )
                        {
                            EdgeID reverseEdgeID = reverseEdge( u, edge );
                            removeEdge( u, edgeID );
                            removeEdge( newEdge.target(), reverseEdgeID );
                        }
                        result -= 2;
                    }
                    break;
                }
            }
        }

        // A new edge is necessary because there were no usable existing edges with
        // same source and target.
        if (forward || backward)
        {
            if ( !simulateOnly )
            {
                const Edge edgeFw(newEdge.target(), newEdge.weight(), newEdge.type(), forward, backward
                    ,newEdge.shortcutMiddle(), newEdge.shortcutEdge1(), newEdge.shortcutEdge2(),
                    newEdge.shortcutOriginalEdgeCount() );
                addEdge( u, edgeFw );
                addEdge( newEdge.target(), Edge(u, edgeFw) );
            }
            result += 2;
        }

        return result;
    }


    /**
    * Returns the ID of the edge (v,u) that corresponds to the
    * specified edge (u,v). If there are several reverse edges,
    * returns the one in the highest level.
    */
    EdgeID reverseEdge(const NodeID u, const EdgeID e) const {
        assert( ! edge(e).isClosed() );
        NodeID v = edge(e).target();
        // if there are several reverse edges,
        // take the one in the highest level
        // therefore: for-loop runs backwards
        for (EdgeID f = lastEdge(v)-1; f >= firstEdge(v); f--) {
            if (edge(f).isReverse(u, edge(e))) return f;
            if (f == 0) break;
        }

        // ERROR if no reverse edge is found
        cerr << "DEBUG: UpdateableGraph::reverseEdge" << endl
             << "       from " << u << " " << edge(e) << endl;
        printAllEdges(cerr, v);
        exit(-1);
    }


    /**
    * Returns the ID of the edge (v,u) that corresponds to the
    * specified edge (u,v). If there are several reverse edges,
    * returns the one in the highest level.
    */
    EdgeID reverseLevelEdge(const NodeID u, const EdgeID e) const {
        assert( ! edge(e).isClosed() );
        NodeID v = edge(e).target();
        // if there are several reverse edges,
        // take the one in the highest level
        // therefore: for-loop runs backwards
        for (EdgeID f = lastEdge(v)-1; f >= firstLevelEdge(v); f--) {
            if (edge(f).isReverse(u, edge(e))) return f;
            if (f == 0) break;
        }

        // ERROR if no reverse edge is found
        cerr << "DEBUG: UpdateableGraph::reverseLevelEdge" << endl
             << "       from " << u << " " << edge(e) << " [" << edge(e).type()  << "]" << endl;
        printAllLevelEdges(cerr, v);
        printAllEdges(cerr, v);
        exit(-1);
    }

    /**
    * Returns the ID of the edge (v,u) that corresponds to the
    * given edge (u,v). Considers only edges in level k.
    */
    EdgeID reverseEdge(const NodeID u, const Edge& e, const LevelID k = 0) const {
        const EdgeID last = lastEdge(e.target());
        for (EdgeID f = firstEdge(e.target()); f < last; f++) {
            if (edge(f).isReverse(u, e)) return f;
        }

        // ERROR if no reverse edge is found
        cerr << "DEBUG: UpdateableGraph::reverseEdge" << endl
             << "       from " << u << " " << e << " [" << e.type() << "]" << endl;
        printAllEdges(cerr, e.target());
        printAllEdges(cerr, u);
        exit(-1);
    }

    // Dynamization currently not supported.
    /** Returns a reference to the 'reliable levels'. */
    /*
    vector<char>& reliableLevels() {return _reliableLevels;}
    */


    // Dynamization currently not supported.
    /** Returns the reliable level of the specified node. */
    /*
    LevelID reliableLevel(const NodeID u) const {
        assert( u < _reliableLevels.size() );
        return _reliableLevels[u];
    }
    */

    // Serialization currently not supported.
    /** Serializes the graph to the given stream. */
    /*
    void serialize(ostream& out) {
        VERBOSE( cout << "datastr::graph::UpdateableGraph::serialize " << noOfNodes()
                      << " " << noOfExistingEdges() << " " << noOfEdges() << endl );
        writePrimitive(out, noOfNodes());
        VectorSerializer< UpdNode, NodeID, ComplexSerializer<UpdNode> >::serialize(out, _nodes);

        // serialize edges:
        // we have three different situations:
        // 1. now: edge groups might not be ordered by nodes, the reserved capacity might be
        //         greater than the assumed capacity (see also 'removeEdge')
        // 2. on disk: edge groups are stored in the right order, without holes
        // 3. after deserialization: edge groups are in the right order, with holes up to the next power of two
        // here, we determine the size of the edge array in situation 3 ...
        EdgeID futureEdgeArraySize = 0;
        for (NodeID u = 0; u < noOfNodes(); u++) {
            EdgeID capacity = degree(u);
            capacity = edgeCapacity(capacity);
            futureEdgeArraySize += capacity;
        }

        // ... and write it.
        writePrimitive(out, futureEdgeArraySize);

        // then, we write all edges (as in situation 2)
        for (NodeID u = 0; u < noOfNodes(); u++) {
            EdgeID last = lastEdge(u);
            for (EdgeID e = firstEdge(u); e < last; e++) edge(e).serialize(out);
        }
        // end of edge serialization

        VectorSerializer< NodeID, NodeID >::serialize(out, _mapExtToIntNodeIDs);
        VERBOSE( cout << "done." << endl );
        VERBOSE( printMemoryUsage(cout) );
    }
    */

    // Serialization currently not supported.
    /** Deserializes the graph from the given stream. */
    /*
    void deserialize(istream& in) {
        VERBOSE( cout << "datastr::graph::UpdateableGraph::deserialize " << flush );
        NodeID n;
        readPrimitive(in, n);
        _nodes.resize(n);

        VectorSerializer< UpdNode, NodeID, ComplexSerializer<UpdNode> >::deserialize(in, _nodes);

        EdgeID m;
        readPrimitive(in, m);
        // reserve some overhead for additional edges (but not too much, a factor of two would waste too much memory)
        _edges.reserve((EdgeID)(1.1 * m));
        _edges.resize(m);

        // deserialize edges, compute the first edge indices, and initialize the nodes
        EdgeID currentFirstEdge = 0;
        for (NodeID u = 0; u < noOfNodes(); u++) {
            node(u).setFirstEdge(currentFirstEdge);
            EdgeID last = lastEdge(u);
            for (EdgeID e = firstEdge(u); e < last; e++)
            {
                edge(e).deserialize(in);
            }
            assert( ! node(u).increasedCapacity() );
            EdgeID size = degree(u);
            EdgeID capacity = size;
            if (node(u).level() > 0) capacity = edgeCapacity(capacity);
            currentFirstEdge += capacity;
        }

        VectorSerializer< NodeID, NodeID >::deserialize(in, _mapExtToIntNodeIDs);
        //_affectedNodes.deserialize(in);
        VERBOSE( cout << noOfExistingEdges() << " " << noOfEdges() << endl; )
        VERBOSE( cout << "done." << endl );
    }
    */

    /** Not used. Provides same interface as the Static/DynamicGraph. */
    EdgeWeight dH(const int dir, const NodeID id, const LevelID level) const {
        assert( false );
        return 0;
    }

    /** Not used. Provides same interface as the Static/DynamicGraph. */
    LevelID coreLevel(const NodeID id) const {
        assert( false );
        return 0;
    }

    /** Not used. Provides same interface as the Static/DynamicGraph. */
    LevelID componentLevel(const NodeID id) const {
        assert( false );
        return 0;
    }

    /**
    * Returns true iff for any edge (u,v), the graph also contains a
    * corresponding edge (v,u).
    */
    bool checkReverseGraphExists() const {
        for (NodeID u = 0; u < noOfNodes(); u++) {
            for (EdgeID e = firstEdge(u); e < lastEdge(u); e++) {
                // if this method fails, the program terminates with an error message
                reverseEdge(u, edge(e));
            }
        }
        return true;
    }


private:
    vector<UpdNode> _nodes;

    vector<Edge> _edges;

    /** Maps original ('external') node IDs to the IDs used internally. */
    vector<NodeID> _mapExtToIntNodeIDs;

    /** For each node, its 'reliable level'. */
    //vector<char> _reliableLevels;

    /**
    * Returns the smallest power of two >= the given size,
    * i.e., returns the 'usual' edge capacity of an edge group
    * of the given size.
    */
    EdgeID edgeCapacity(const EdgeID size) const {
        EdgeID capacity = 1;
        while (capacity < size) capacity <<= 1;
        return capacity;
    }

    /**
    * Returns the actual edge capacity of node u whose edge
    * group has the given size. The result is equal to the
    * 'usual' edge capacity 2^k (see above) or 2^{k+1} iff
    * the node u has an 'increased capacity'.
    */
    EdgeID edgeCapacity(const NodeID u, const EdgeID size) const {
        EdgeID capacity = edgeCapacity(size);
        if (node(u).increasedCapacity()) capacity <<= 1;
        return capacity;
    }

    /**
    * To the last edge group in the edge array with size s and
    * capacity c, we add (c - s) 'holes'.
    */
    void appendClosedEdges(const EdgeID size, EdgeID capacity) {
        for (; capacity > size; capacity--) {
            // append a 'hole'
            _edges.push_back(Edge(SPECIAL_NODEID, Weight::MAX_VALUE, false, false, false));
        }
    }

    /**
     * Used during the initial construction of the graph. Sort the list of edges
     * after source,target and finally level.
     */
    struct CompareEdgesBySourceTargetLevel : public binary_function<CompleteEdge, CompleteEdge, bool>
    {
        UpdateableGraph* graph;
        bool operator()(const CompleteEdge& a, const CompleteEdge& b)
        {
            return (a.source() < b.source()
                || (a.source() == b.source() && graph->node(a.target()).level() < graph->node(b.target()).level()));
        }
    };


    /**
    * Initializes this UpdateableGraph from a given list of edges and node levels.
    * Note that the overlay graph hierarchy is constructed later.
    * Also the preprocessing of contraction hierarchies assumes that
    * the remaining nodes (at the beginning all nodes) are in the same highest
    * level.
    */
    void construct(vector<CompleteEdge>& edges, const vector<LevelID>& nodeLevels) {
        assert( ! nodeLevels.empty() );

        // reserve RAM in advance
        // necessary because memory space will get fragmented
        // and large graphs like ORTEC need more than 2GB
        //_edges.reserve(5 << 26);

        // currently no node mapping is used
        _mapExtToIntNodeIDs.resize(nodeLevels.size());
        _nodes.resize(nodeLevels.size());
        for (NodeID u = 0; u < nodeLevels.size(); u++) {
            _mapExtToIntNodeIDs[u] = u;
            _nodes[u].setLevel(nodeLevels[u]);
        }

        // re-arrange edges (sort)
        CompareEdgesBySourceTargetLevel comp;
        comp.graph = this;
        sort(edges.begin(), edges.end(), comp);

        // build adjacency array and indices
        // for each node into the adjacency array
        EdgeID e = 0;
        for (NodeID u = 0; u < nodeLevels.size(); u++) {
            UpdNode& nodeU = _nodes[u];
            nodeU.setFirstEdge(_edges.size());
            while ((e < edges.size()) && (edges[e].source() == u) && (node(edges[e].target()).level() < node(u).level())) {
                assert( edges[e].weight() > 0 );
                _edges.push_back(edges[e]);
                e++;
            }
            nodeU.setFirstLevelEdge(_edges.size());
            while ((e < edges.size()) && (edges[e].source() == u)) {
                assert( edges[e].weight() > 0 );
                _edges.push_back(edges[e]);
                e++;
            }
            nodeU.setLastEdge(_edges.size());
            EdgeID offset = _edges.size() - firstEdge(u);
            node(u).setIncreasedCapacity();
            appendClosedEdges(offset, edgeCapacity(u, offset));
        }

        assert( checkReverseGraphExists() );

        VERBOSE( printMemoryUsage(cout) );
    }

    /** Used for debugging purposes. Prints all edges (u,v). */
    void printAllEdges(ostream& out, const NodeID u) const {
        out << "all edges from " << u << ": ";
        bool first = true;
        EdgeID lastE = lastEdge(u);
        for (EdgeID e = firstEdge(u); e < lastE; e++) {
            if (! first) out << " | ";
            first = false;

            out << edge(e) << " [" << edge(e).type()  << "]";
        }
        out << endl;
    }

    /** Used for debugging purposes. Prints all edges (u,v). */
    void printAllLevelEdges(ostream& out, const NodeID u) const {
        out << "all level edges from " << u << ": ";
        bool first = true;
        EdgeID lastE = lastEdge(u);
        for (EdgeID e = firstLevelEdge(u); e < lastE; e++) {
            if (! first) out << " | ";
            first = false;

            out << edge(e) << " [" << edge(e).type()  << "]";
        }
        out << endl;
    }

    /** Prints the total memory usage of this UpdateableGraph. */
    void printMemoryUsage(ostream& out) const {
        out << "** memory usage on hard disk [MB (bytes per node)] **" << endl;

        unsigned int memoryTotal = 0;
        memoryTotal += printMemoryUsage(out, "edges", noOfExistingEdges() * sizeof(Edge));
        memoryTotal += printMemoryUsage(out, "node data", _nodes.size() * sizeof(UpdNode));
        memoryTotal += printMemoryUsage(out, "node ID mapping", _mapExtToIntNodeIDs.size() * sizeof(NodeID));

        printMemoryUsage(out, "TOTAL", memoryTotal);
        out << endl;
    }

    /** Prints the memory usage of one particular data structure of this UpdateableGraph. */
    unsigned int printMemoryUsage(ostream& out, const string descr, const unsigned int mem) const {
        unsigned int megaBytes = (unsigned int)((mem / (double)(1024*1024)) + 0.5);
        unsigned int bytesPerNode = (unsigned int)((mem / (double)noOfNodes()) + 0.5);
        out << "   " << descr << ": " << megaBytes << " (" << bytesPerNode << ")" << endl;
        return mem;
    }

};

} } // namespace

#endif // _DATASTR_GRAPH_UPDATEABLEGRAPH_H
