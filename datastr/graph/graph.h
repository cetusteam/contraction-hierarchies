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


#ifndef GRAPH_H
#define GRAPH_H


#include <iostream>
#include <iomanip>
#include <vector>


using namespace std;

/** 
 * Type of a node ID.
 * Also used for everything that is bounded by the number of nodes,
 * e.g. the number of components.
 */
typedef unsigned int NodeID;

/** Type of an edge ID. */
typedef unsigned int EdgeID;

/** Type of a source/target node pair. */
typedef pair<NodeID, NodeID> stPair;
typedef vector<stPair> stPairs;

/** 
 * Special node ID which is normally not used
 * so that it can be used as flag.
 */
static const NodeID SPECIAL_NODEID = __INT_MAX__ * 2U + 1;

static const EdgeID MAX_EDGEID = __INT_MAX__ * 2U;

/** Type of a level ID. */
typedef int LevelID;

#include "edge.h"
#include "path.h"

#endif // GRAPH_H
