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

#ifndef _COMMAND_H
#define _COMMAND_H

#include <boost/regex.hpp>
#include <math.h>

class Command {
    public:
        virtual ~Command() {};
        virtual int main(int argc, char *argv[]) = 0;

    /** create int vector from string, various formats supported
      * - comma-separated list
      * - a^b:c = a^b,a^b+1,...,a^c (ld(n) for c possible)
      * - axb = a,a,a,...,a (b times, n,n-1 for b possible)
      */
    template< typename T >
    static void createVector(const string& str, vector<T>& result, T n = 0) 
    {
        boost::regex pattern1("^((n|\\d+)\\*)?(\\d+)\\^(-?\\d+):(-?\\d+|(-?ld\\(n\\))\\+?(-?\\d+)?)(\\+(-?\\d+))?$");
        boost::regex pattern2("^(\\d+)x(-?\\d+|n|n-1)$");
        boost::regex pattern3("^n//(\\d+)x(\\d+)$");
        boost::smatch what;
        if (str == "") return;
        if (boost::regex_match(str, what, pattern1)) {
            double base = atof(what[3].str().data());
            double start = atof(what[4].str().data());
            double stop = 0;
            double step = 1;
            double factor = 1;
            if (what[2] == "n")
            {
                factor = n;
            }
            else if (what[1] != "")
            {
                factor = atoi(what[2].str().data());
            }
            if (what[6] == "ld(n)")
            {
                stop = log(n)/log(2);
            }
            else if (what[6] == "-ld(n)")
            {
                stop = -log(n)/log(2);
            }
            else
            {
                stop = atoi(what[5].str().data());
            }
            if (what[7] != "")
            {
                stop += atoi(what[7].str().data());
            }
            if (what[8] != "")
            {
                step = atoi(what[9].str().data());
            } 
            if (step > 0)
            {
                for (double i  = start; i <= stop; i += step)
                {
                    result.push_back((T)(factor*pow(base,i)));
                }
            }
            else
            {
                for (double i  = start; i >= stop; i += step)
                {
                    result.push_back((T)(factor*pow(base,i)));
                }
            }
        }
        else if (boost::regex_match(str, what, pattern2)) {
            T a = (T)atof(what[1].str().data());
            T times;
            if (what[2] == "n")
            {
                times = n;
            }
            else if (what[2] == "n-1")
            {
                times = n-1 ;
            }
            else
            {
                times = (T)atof(what[2].str().data());
            }
            for (T i = 0; i < times; i++)
            {
                result.push_back(a);
            }
        }
        else if (boost::regex_match(str, what, pattern3)) 
        {
            int div = atoi(what[1].str().data());
            int times = atoi(what[2].str().data());
            T step = n;
            for ( int i = 0; i < times; i++ ) step /= div;
            for ( int i = 1; i < times; i++ )
            {
                result.push_back(step);
                step *= div;
            }
        }
        else if (str == "n")
        {
            result.push_back(n);        
        }
        else 
        {
            char ch = ',';
            string::size_type pos = 0;
            string::size_type i;
            while ((i = str.find(ch,pos)) != string::npos) {
                string s = str.substr(pos, i-pos);
                if (s == "n")
                {
                    result.push_back(n);
                }
                else
                {
                    result.push_back((T)atof(s.data()));
                }
                pos = i+1;
            }
            result.push_back((T)atof(str.substr(pos).data()));
        }
    }
    
    };

#endif // _COMMAND_H
