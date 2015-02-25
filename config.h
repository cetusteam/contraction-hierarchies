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

#ifndef CONFIG_H
#define CONFIG_H

// switch asserts on/off
#define NDEBUG

#include <cassert>

// select type of edge weights
//#include "floatEdgeWeight.h" // actually double !
#include "uintEdgeWeight.h"

// activate counting of several statistics
#define COUNTING(x)
#include "stats/counter.h"
#define LOG_TIME(x)

// activate verbose output
#define VERBOSE(x) x
#define VERBOSE2(x)
#define VERBOSE_CONTRACT(x)
#define VERBOSE_CONTRACT_VORONOI(x)

// on, store expand information, off do not
// required if paths found by query algorithm
// should be expanded to pahts in the original graph
// Note: affects the output of a sgr-file
#define USE_CH_EXPAND

#ifdef USE_CH_EXPAND
#define CH_EXPAND(x) x
#else
#define CH_EXPAND(x)
#endif

// count the number of original edges a shortcut represents
// required for node ordering with parameter (-e)
// This adds another 4 bytes to one edge, should be disabled for
// queries.
// Note: affects the output of a sgr-file
#define COUNT_SHORTCUT_ORIGINAL_EDGES

// activate counting of cache misses
//#define PAPI
//#include <papi.h>

// select strategy: during bidirectional search, we have to
// decide in which order the forward and the backward search are processed

// prefer smaller element at the heads of the priority queues
#define PQ_BIDIR_MIN(x)
// prefer smaller priority queue (with less elements)
#define PQ_BIDIR_SIZE(x)
// alternate forward and backward search
#define PQ_BIDIR_ALTERNATE(x) x

#endif // CONFIG_H
