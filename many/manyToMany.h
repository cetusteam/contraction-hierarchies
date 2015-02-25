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


#ifndef MANYTOMANY_H
#define MANYTOMANY_H


/**
 * Represents a matrix that can store a result of
 * a many-to-many computation.
 */
template <typename value_type, typename size_type = unsigned int>
class Matrix
{
    friend ostream& operator<<( ostream& os, const Matrix& matrix ) {
        for (size_type i = 0; i < matrix.noOfEntries(); i++) {
            if (i > 0) {
                if (i % matrix.noOfCols() == 0) os << endl; else os << " ";
            }
            os << matrix._data[i];
        }
        return os;
    }
    
public:
    bool operator== (const Matrix<value_type, size_type>& matrix) const {
        if (noOfRows() != matrix.noOfRows()) return false;
        if (noOfCols() != matrix.noOfCols()) return false;
        return (_data == matrix._data);
    }

    /** Constructor. */
    Matrix(const size_type r, const size_type c) : _rows(r), _cols(c), _data(r*c) {}

    /** Inits all entries to a given value. */
    void init(const value_type& v) {
        for (size_type i = 0; i < noOfEntries(); i++) _data[i] = v;
    }

    size_type noOfRows() const {return _rows;}
    size_type noOfCols() const {return _cols;}
    size_type noOfEntries() const {return _data.size();}

    value_type value(const size_type r, const size_type c) const {return _data[index(r, c)];}

    /**
    * Sets a specified entry to the given value.
    * @param r the row index
    * @param c the column index
    * @param v the value
    */
    void set(const size_type r, const size_type c, const value_type& v) {
        _data[index(r, c)] = v;
    }

    /**
    * Sets a specified entry to the given value
    * if the new value is less than the old value.
    * @param r the row index
    * @param c the column index
    * @param v the value
    */
    void improve(const size_type r, const size_type c, const value_type& v) {
        improve(index(r, c), v);
    }

    /**
    * Sets a specified entry to the given value
    * if the new value is less than the old value.
    * @param i the entry index
    * @param v the value
    */
    void improve(const size_type i, const value_type& v) {
        assert( i < noOfEntries() );
        if (v < _data[i]) _data[i] = v;
    }

    /**
     * Returns the index of a specified entry.
     * @param r the row index
     * @param c the column index
     * @return the entry index
     */
    size_type index(const size_type r, const size_type c) const {
        assert( r < noOfRows() );
        assert( c < noOfCols() );
        const size_type i = r * noOfCols() + c;
        assert( i < noOfEntries() );
        return i;
    }
    
private:
    /** The number of rows. */
    const size_type _rows;

    /** The number of columns. */
    const size_type _cols;

    /**
     * All entries.
     * The unidimensional array is indexed as a two-dimensional array.
     * @see index
     */
    vector<value_type> _data;
};


/**
 * Provides methods to perform many-to-many computations.
 * @param Graph the graph type
 * @param DijkstraFW the type of Dijkstra used for the forward search
 * @param DijkstraBW the type of Dijkstra used for the backward search
 * @param performBucketScans performing the bucket scans is a crucial part of the many-to-many computation.
 *            It can be switched off to allow time measurements of the forward search
 *            without accounting for the bucket scans.
 */
template <typename Graph, typename DijkstraFW, typename DijkstraBW, bool performBucketScans = true>
class ManyToMany : public SearchSpaces
{
public:
    /** Constructor. */
    ManyToMany(Graph *const g, const LevelID earlyStopLevel) : _g(g), _dFW(g), _dBW(g), _earlyStopLevel(earlyStopLevel) {
        // early stop level currently not supported
        //_dFW.setEarlyStopLevel(earlyStopLevel);
        //_dBW.setEarlyStopLevel(earlyStopLevel);
    }

    /**
     * Reference implementation.
     * Performs point-to-point queries in order to fill the distance table.
     */
    void computeMatrixNaive(const vector<NodeID>& sources, const vector<NodeID>& targets, Matrix<EdgeWeight>& result) {
        assert( sources.size() == result.noOfRows() );
        assert( targets.size() == result.noOfCols() );

        VERBOSE( cout << "computing reference solution ..." << endl );

        VERBOSE( Percent progress(sources.size()) );
        for (NodeID u = 0; u < sources.size(); u++) {
            VERBOSE( progress.printStatus(u) );
            for (NodeID v = 0; v < targets.size(); v++) {
                const NodeID s = sources[u];
                const NodeID t = targets[v];
                const EdgeWeight w = _dFW.bidirSearch(s, t);
                _dFW.clear();
                result.set(u, v, w);
            }
        }

        VERBOSE( cout << "done." << endl );
    }

    /**
     * Efficient implementation.
     * Applies many-to-many algorithm in order to fill the distance table.
     */
    void computeMatrix(const vector<NodeID>& sources, const vector<NodeID>& targets, Matrix<EdgeWeight>& result) {
        assert( sources.size() == result.noOfRows() );
        assert( targets.size() == result.noOfCols() );

        VERBOSE( cout << "computing " << sources.size() << " x " << targets.size() << " table ..." << endl );
        COUNTING( counter.reset() );
        const double start = timestamp();

        // init table
        result.init(Weight::MAX_VALUE);

        // backward search
        VERBOSE( cout << "backward search" << endl );
        VERBOSE( Percent progress(targets.size()) );
        for (NodeID v = 0; v < targets.size(); v++) {
            VERBOSE( progress.printStatus(v) );
            const NodeID t = targets[v];
            setCurrentNode(v);
            _dBW.bidirSearch(SPECIAL_NODEID, t);
            _dBW.obtainRelevantSearchSpace(*this);
            _dBW.clear();
        }
        double elapsedTime = timestamp() - start;
        LOG_TIME( if (! performBucketScans) cerr << elapsedTime << " " );
        VERBOSE( cout << elapsedTime << " s" << endl );
        COUNTING( cout << "backward search space: " << (long long)counter.count(COUNT_DEL_MIN, COUNT_TEMP) << endl );
        COUNTING( counter.reset() );

        VERBOSE( cout << "  search space sizes: int " << _searchSpacesBwDynInt.size()
                      << " + short " << _searchSpacesBwDynShort.size() << endl );

        // only the 'int'-search space is used
        assert( _searchSpacesBwDynShort.size() == 0 );
        
        // sort and copy
        VERBOSE( cout << "sort and copy" << endl );
        _searchSpacesBwDynInt.sort(_g->noOfNodes());
        _searchSpacesBwInt = _searchSpacesBwDynInt;
        elapsedTime = timestamp() - start;
        LOG_TIME( if (! performBucketScans) cerr << elapsedTime << " " );
        VERBOSE( cout << elapsedTime << " s" << endl );

        // forward search
        VERBOSE( cout << "forward search" << endl );
        VERBOSE( progress.reinit(sources.size()) );
        COUNTING( unsigned long long bucketScans = 0 );
	    //COUNTING( unsigned long long bucketScansTop = 0 );
        for (NodeID u = 0; u < sources.size(); u++) {
            VERBOSE( progress.printStatus(u) );
            _dFW.bidirSearch(sources[u], SPECIAL_NODEID);
            _dFW.obtainRelevantSearchSpace(*this);
            _dFW.clear();

            if (performBucketScans) {
                const NodeID matrixIndexOffset = result.index(u, 0);
                for (NodeID i = 0; i < _searchSpaceFW.size(); i++) {
                    const NodeID via = _searchSpaceFW[i].target();
                    const EdgeWeight distFW = _searchSpaceFW[i].weight();
                
                    const ISSInt::const_iterator endInt = _searchSpacesBwInt.end(via);
                    for (ISSInt::const_iterator it = _searchSpacesBwInt.begin(via); it != endInt; it++) {
                        result.improve(matrixIndexOffset + it->origin(), distFW + it->dist());
                        COUNTING( bucketScans++ );
			            //COUNTING( if (_g->node(via).level() >= _earlyStopLevel) bucketScansTop++ );
                    }
                }
            }
            
            _searchSpaceFW.clear();
        }
        elapsedTime = timestamp() - start;
        LOG_TIME( cerr << elapsedTime << " " );
        LOG_TIME( if (performBucketScans) cerr << endl );
        VERBOSE( cout << elapsedTime << " s" << endl );
        COUNTING( cout << "forward search space: " << (long long)counter.count(COUNT_DEL_MIN, COUNT_TEMP) << endl );
        COUNTING( cout << "bucket scans (all): " << bucketScans << endl );
	    //COUNTING( cout << "bucket scans (top): " << bucketScansTop << endl );
    }
    
private:
    Graph *const _g;
    DijkstraFW _dFW;
    DijkstraBW _dBW;
    const LevelID _earlyStopLevel;
};

#endif // MANYTOMANY_H
