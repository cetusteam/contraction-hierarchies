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

using namespace std;

#include <iostream>
#include <iomanip>
#include <fstream>

#include "config.h"

#include "stats/counter.h"
#include "stats/utils.h"

Counter counter;

#include "command/NodeOrder.h"
#include "command/Construct.h"
#include "command/TNR.h"

/** The main program. */
int main(int argc, char *argv[])
{
    // The program is controlled by command-line arguments. The order of those
    // arguments is important. The first argument specififies the the Command-
    // class that is used.
    int opt = getopt(argc, argv, "sct");
    Command* m = NULL;
    int result = 0;
    switch (opt) {

        /* ************* *
         * node ordering *
         * ************* */
        case 's':
            m = new command::NodeOrder();
            break;

        /* ********************** *
         * hierarchy construction *
         * ********************** */
        case 'c':
            m = new command::Construct();
            break;
            
        /* ******************** *
         * transit node routing *
         * ******************** */
        // this does not actually perform TNR, only some preliminary tests
        // are currently performed
        case 't':
            m = new command::TNR();
            break;
            
    }
    if (m != NULL) {
        result = m->main(argc, argv);
        delete m;
    } else {
        cout << "For help see docu/index.html" << endl;
    }
    return result;
}


// doesn't look nice, but required by the compiler (gcc 4)
const EdgeWeight Weight::MAX_VALUE;
const EliminationWeight::Type EliminationWeight::MAX_VALUE;
const EliminationWeight::Type EliminationWeight::MIN_VALUE;
const EdgeWeight Weight::MAX_INTEGER;
const int datastr::graph::UpdateableGraph::LOOK_FOR_SECOND_EDGE_BACKWARD;
const int datastr::graph::UpdateableGraph::LOOK_FOR_SECOND_EDGE_FORWARD;
