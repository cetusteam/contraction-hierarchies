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

#ifndef ELIMINATIONWEIGHT_H
#define ELIMINATIONWEIGHT_H

#include <limits>

/**
 * Encapsulates several special values of the elimination weights or 
 * priority terms for the node ordering. This class is required by the
 * BinaryHeap priority queue.
 * Here: elimination weight type = double
 */
class EliminationWeight
{
public:
    typedef double Type;
    
    /**
     * max value of an edge weight (due to the chosen type)
     * (used to represent 'infinity')
     */
    static const Type MAX_VALUE = __DBL_MAX__;
    static const Type MIN_VALUE = -__DBL_MAX__;
    
};

#endif // ELIMINATIONWEIGHT_H
