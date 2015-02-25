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

#ifndef COUNTER_H
#define COUNTER_H

#include <vector>


static const int COUNT_DEL_MIN = 0;
static const int COUNT_DECR_KEY = 1;
static const int COUNT_SKIPPED_EDGES = 2;
static const int COUNT_RELAXED_EDGES = 3;
static const int COUNT_RELAXED_EDGES_SUCCESS = 4;
static const int COUNT_DOUBLE_LOG_END = 5; // end of the counters that are written to the log file

static const int COUNT_STALL_OPS = 5;
static const int COUNT_STALL_STEPS = 6;
static const int COUNT_INCR_KEY = 7;
static const int COUNT_DOUBLE_END = 8;


static const int COUNT_TEMP = 0;
static const int COUNT_AKKU = 1;

static const int COUNT_IMPROVE_AFTER_MEETING = 0;
static const int COUNT_DEL_MIN_AFTER_MEETING = 1;
static const int COUNT_DEL_MIN_AFTER_IMPROVE = 2;
static const int COUNT_SINGLE_END = 3;

static const int COUNT_ENTRANCE_POINTS = 1;

static const unsigned int COUNT_LEVEL_END = 50;


/** 
 * Provides counters for one level.
 * The current implementation only supports a delMin counter
 * and an entrance point counter.
 */
class LevelCounter
{
public:
    LevelCounter() : _delMins(0), _entrancePoints(0) {}
    
    void inc(int i) {
        assert( (i == COUNT_DEL_MIN) || (i == COUNT_ENTRANCE_POINTS) );
        if (i == COUNT_DEL_MIN) _delMins++;
        if (i == COUNT_ENTRANCE_POINTS) _entrancePoints++;
    }

    void print(ostream& os, unsigned int level) {
        if ((_delMins == 0) && (_entrancePoints == 0)) return; // skip an unused counter
        os << "Level " << level << ": delMins = " << _delMins
           << "; entrance points = " << _entrancePoints << endl;
    }
    
private:
    double _delMins;
    double _entrancePoints;
};


/** Encapsulates several LevelCounter objects. */
class LevelCounters
{
    /** Prints all values of the level counters. */
    friend ostream& operator<<( ostream& os, LevelCounters &counter ) {
        for (unsigned int level = 0; level < COUNT_LEVEL_END; level++) {
            counter._levelCounter[level].print(os, level);
        }
        return os;
    }
    
public:
    /** Increments the specified counter of the given level. */
    void inc(unsigned int level, int i) {
        assert( level < COUNT_LEVEL_END );
        if (i == COUNT_ENTRANCE_POINTS) {
            for (unsigned int x = 0; x < _entrancePointLevels.size(); x++)
                if (_entrancePointLevels[x] == level) return;
            _entrancePointLevels.push_back(level);
            return;
        }
        _levelCounter[level].inc(i);
    }

    /** Dispatches the cached values. */
    void dispatch(int i) {
        assert( i == COUNT_ENTRANCE_POINTS );
        for (unsigned int x = 0; x < _entrancePointLevels.size(); x++)
            _levelCounter[_entrancePointLevels[x]].inc(i);
        _entrancePointLevels.clear();
    }
    
private:
    /** The level counters. */
    LevelCounter _levelCounter[COUNT_LEVEL_END];

    /** Caches all levels that the current node is an entrance point of. */
    vector<unsigned int> _entrancePointLevels;
};


/** Provides counters for various statistics. */
class Counter
{
    /** Prints all values of the single and level counters. */
    friend ostream& operator<<( ostream& os, Counter &counter ) {
        os << setprecision(20);
        os << "Counter" << endl;
        for (int i = 0; i < COUNT_DOUBLE_END; i++) {
            os << counter.descrDouble(i) << ": " << counter.count(i, COUNT_AKKU) << endl;
        }
        for (int i = 0; i < COUNT_SINGLE_END; i++) {
            os << counter.descrSingle(i) << ": " << counter.count(i) << endl;
        }
        os << counter._levelCounters;
        return os;
    }
    
public:
    Counter() {
        reset();
        _descrSingle[0] = "ImprovementsAfterMeeting";
        _descrSingle[1] = "DeleteMinAfterMeeting";
        _descrSingle[2] = "DeleteMinAfterImprovement";
        
        _descrDouble[0] = "DeleteMin";
        _descrDouble[1] = "DecreaseKey";
        _descrDouble[2] = "SkippedEdges";
        _descrDouble[3] = "RelaxedEdges";
        _descrDouble[4] = "SuccessfullyRelaxedEdges";
        _descrDouble[5] = "StallOperations";
        _descrDouble[6] = "StallSteps";
        _descrDouble[7] = "IncreaseKey";
    }

    /** Resets all counters. */
    void reset() {
        resetSingle();
        resetDouble();
    }

    /** Resets all single counters. */
    void resetSingle() {
        for (int i = 0; i < COUNT_SINGLE_END; i++) reset(i);
    }

    /**
     * Resets all double counters.
     * If "COUNT_TEMP" is given as argument, only the "temp" counters are reset
     * (but not the accumulated ones).
     */
    void resetDouble(int j = COUNT_AKKU) {
        for (int i = 0; i < COUNT_DOUBLE_END; i++) reset(i, j);
    }
    
    /** Resets the specified single counter. */
    void reset(int i) {
        assert( i < COUNT_SINGLE_END );
        _countsSingle[i] = 0;
    }

    /**
     * Resets the specified double counter.
     * If the second argument is COUNT_TEMP, only the "temp" counter is reset,
     * otherwise (COUNT_AKKU) both counters are reset.
     * @param i the counter ID
     * @param j COUNT_TEMP or COUNT_AKKU (see above)
     */
    void reset(int i, int j) {
        assert( i < COUNT_DOUBLE_END );
        assert( j < 2 );
        _countsDouble[COUNT_TEMP][i] = 0;
        if (j == COUNT_AKKU) _countsDouble[COUNT_AKKU][i] = 0;
    }

    /** Increments the specified single counter. */
    void incSingle(int i) {
        assert( i < COUNT_SINGLE_END );
        _countsSingle[i]++;
    }

    /** Adds the given value to the specified single counter. */
    void addSingle(int i, double value) {
        assert( i < COUNT_SINGLE_END );
        _countsSingle[i] += value;
    }

    /** Increments both specified double counters. */
    void incDouble(int i) {
        assert( i < COUNT_DOUBLE_END );
        _countsDouble[COUNT_TEMP][i]++;
        _countsDouble[COUNT_AKKU][i]++;
    }

    /** Returns a reference to the level counters. */
    LevelCounters& levels() {return _levelCounters;}

    /** Returns the value of the specified single counter. */
    double count(int i) const {
        assert( i < COUNT_SINGLE_END );
        return _countsSingle[i];
    }

    /**
     * Returns the value of the specified double counter.
     * @param i the counter ID
     * @param j COUNT_TEMP or COUNT_AKKU
     */
    double count(int i, int j) const {
        assert( i < COUNT_DOUBLE_END );
        return _countsDouble[j][i];
    }
    
    /** Returns a description of the meaning of the specified single counter. */
    string descrSingle(int i) const { return _descrSingle[i]; }

    /** Returns a description of the meaning of the specified double counter. */
    string descrDouble(int i) const { return _descrDouble[i]; }
    
private:
    /** The values of the single counters. */
    double _countsSingle[COUNT_SINGLE_END];
    
    /** The values of the double counters. */
    double _countsDouble[2][COUNT_DOUBLE_END];
    
    /** The descriptions of the meaning of the single counters. */
    string _descrSingle[COUNT_SINGLE_END];

    /** The descriptions of the meaning of the double counters. */
    string _descrDouble[COUNT_DOUBLE_END];
    
    /** The level counters. */
    LevelCounters _levelCounters;
};

#endif // COUNTER_H
