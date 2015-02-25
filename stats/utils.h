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

#ifndef UTILS_H
#define UTILS_H

#include <iostream>
#include <sys/time.h>


/** Returns a timestamp ('now') in seconds (incl. a fractional part). */
inline double timestamp() {
    struct timeval tp;
    gettimeofday(&tp, NULL);
    return double(tp.tv_sec) + tp.tv_usec / 1000000.;
}


/**
 * Provides methods that can be used to display the
 * progress of a procedure.
 */
class Percent
{
 private:
    typedef int percentInt;
    typedef unsigned int valueInt;
    
 public:
    /**
     * Constructor.
     * @param maxValue the value that corresponds to 100%
     *                 (or 0% if reverse counting is activated)
     * @param reverse instead of raising from 0 to maxValue,
     *                the value falls from maxValue to 0
     * @param step the progress is shown in steps of 'step' percent
     */
    Percent(valueInt maxValue, bool reverse = false, percentInt step = 2) {
        reinit(maxValue, reverse, step);
    }

    /** Reinitializes this object. */
    void reinit(valueInt maxValue, bool reverse = false, percentInt step = 2) {
        _maxValue = maxValue;
        _intervalPercent = _maxValue / 100;
        _nextThreshold = _intervalPercent;
        _reverse = reverse;
        _lastPercent = 0;
        _step = step;
        
        if (reverse) _nextThreshold = _maxValue - _intervalPercent;
    }

    /** If there has been significant progress, display it. */
    void printStatus(valueInt currentValue) {
        if ( ! _reverse) {
            if (currentValue >= _nextThreshold) {
                _nextThreshold += _intervalPercent;
                printPercent( currentValue / (double)_maxValue * 100 );
            }
            if (currentValue + 1 == _maxValue) finish();            
        }
        else {
            if (currentValue <= _nextThreshold) {
                _nextThreshold -= _intervalPercent;
                printPercent( 100 - (currentValue / (double)_maxValue * 100) );
            }
        }
    }
        
 private:
    valueInt _maxValue;
    valueInt _intervalPercent;
    valueInt _nextThreshold;
    bool _reverse;
    percentInt _lastPercent;
    percentInt _step;    

    /** Displays the new progress. */
    void printPercent(double percent) {
        while (percent >= _lastPercent+_step) {
            _lastPercent+=_step;
            if (_lastPercent % 10 == 0) {
                std::cout << " " << _lastPercent << "% ";
            }
            else {
                std::cout << ".";
            }
            std::cout.flush();
        }
    }

    void finish() {
        std::cout << " 100%" << std::endl;
    }
};


#endif // UTILS_H
