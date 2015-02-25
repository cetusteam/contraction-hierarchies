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

#ifndef RANDOM_H
#define RANDOM_H

#include <utility>
#include <cassert>

/**
 * Generate a pseudo-random number that lies in the range [a,b).
 */
double random(double a, double b)
{
	assert(a < b);
	
	double res;

	do {
		// get a random value that lies in [0,b-a)
		res = ( (double)rand() / (double)(RAND_MAX + 1.0) ) * (b - a);
		// lift it such that it lies in [a,b)
		res += a;
	} while ( ! (a <= res && res < b) ); // if it does not lie in between [a,b) try again
	
	assert( a <= res && res < b );
	return res;
}

#endif /* RANDOM_H */
