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
#ifndef MISC_H_
#define MISC_H_


#include "../io/createGraph.h"
#include "../datastr/UpdateableGraph.h"


template <typename R, typename S, typename T>
class Triple
{
public:
	R first;
	S second;
	T third;

	Triple(R fst, S scnd, T thd)
	{
		first = fst;
		second = scnd;
		third = thd; 
	}
};

template <typename R, typename S, typename T, typename U>
class Quadruple
{
public:
	R first;
	S second;
	T third;
	U fourth;

	Quadruple(R fst, S scnd, T thd, U fth)
	{
		first = fst;
		second = scnd;
		third = thd; 
		fourth = fth;
	}
};


Edge *edgeFromTo(datastr::graph::UpdateableGraph *gr, NodeID src, NodeID tgt)
{
	EdgeID first_edge_id = gr->firstEdge(src);
	EdgeID last_edge_id = gr->lastEdge(src);
	for (EdgeID edge_id = first_edge_id; edge_id < last_edge_id; edge_id++)
	{
		Edge& edge = gr->edge(edge_id);
		if (edge.isDirected(FORWARD) && edge.target() == tgt) return &edge;
	}
	return NULL;
}

/** Create constant edge duration functions from the static edge weights. */
void createConstTimeFns(datastr::graph::UpdateableGraph *graph)
{
	for (NodeID n = 0 ; n < graph->noOfNodes() ; n++)
	{
		// iterate over all outgoing edges (u,v) of node u
		EdgeID first_edge_id = graph->firstEdge(n);
		EdgeID last_edge_id = graph->lastEdge(n);
		for (EdgeID edge_id = first_edge_id; edge_id < last_edge_id; edge_id++)
		{
			Edge& edge = graph->edge(edge_id);
			if (edge.getTimeFn() == NULL)
			{
				double weight = (double)edge.weight();
				assert (weight >= 0.0);

				vector< pair<double,double> > timeSample;
				timeSample.push_back( pair<double,double>(0.0, weight) );
				edge.setTimeFn(new TimeFn(timeSample));
			}
		}
	}
}

/**
 * Read an UpdateableGraph from a .ddsg file and a .hnc file.
 * @param gr_file	The full name of the .ddsg file.
 * @param ch_file	The full name of the .hnc file.
 */
datastr::graph::UpdateableGraph *readGraph(const char *gr_fln, const char *ch_fln)
{
	cout << "Reading graph data from '" << fln << "'...\n\n";

	ifstream inFile(fln);
	if ( ! inFile.is_open() ) {
		cerr << "Unable to open file.\n";
		exit(-1);
	}
	// read an updateable graph with static edge weights
	const string ch_file(ch_fln);
	UpdateableGraph *res = importGraphListOfEdgesUpdateable(inFile, false, false, ch_file);
	inFile.close();
	return res;
}



/**
 * Read time functions for a given graph from a .profile file.
 * @param gr	The graph.
 * @param fln	The full name of the .profile file.
 * @param unit	Value of a unit in seconds.
 */
void readTimeFns(datastr::graph::UpdateableGraph *gr, const char *fln, double unit)
{
	cout << "Reading time dependent edge weights from '" << fln << "'...\n";
	
	int n_removed_points = 0, n_const_fns = 0;

	ifstream inFile(fln);
	if ( ! inFile.is_open() )
	{
		cerr << "Unable to open file.\n";
		exit(-1);
	}

	// read duration functions from file and associate
	// them to the appropriate edges
	while ( ! inFile.eof() )
	{
		NodeID src, tgt;
		int emptyTime;
		vector< pair<double,double> > timeSample;


		inFile >> src;
		inFile >> tgt;

		Edge *edge = edgeFromTo(gr, src, tgt);
		assert (edge);
		
		inFile >> emptyTime;
		assert (edge->weight() == emptyTime);
		
		int x, x_before, x_next, x_first;

		inFile >> x;
		if ( !inFile ) break;
		pair<double,double> point(0, (double)x * ((1.0/TIME_UNIT) * unit) );
		timeSample.push_back(point);

		x_first = x;

		inFile >> x_next;
		for (int i=1; i<23; i++)
		{
			x_before = x;
			x = x_next;
			inFile >> x_next;

			if (x != x_before || x != x_next)
			{
				pair<double,double> point((double)(i*3600) / TIME_UNIT, (double)x * ((1.0/TIME_UNIT) * unit) );
				timeSample.push_back(point);
			}
			else n_removed_points++;
		}
		
		if (x_next != x || x_next != x_first) 
		{
			pair<double,double> point((double)(23*3600) / TIME_UNIT, (double)x_next * ((1.0/TIME_UNIT) * unit) );
			timeSample.push_back(point);
		}
		TimeFn *f = new TimeFn(timeSample);
		assert (f->isFifo());
		edge->setTimeFn(f);
	}
	inFile.close();

	int n_points = TimeFn::getTotalNPoints();
	int mem_usage = (int)((double)n_points / (double)(1024 * 1024)) * sizeof(double);

	int alloc_n_points = TimeFn::getTotalAllocNPoints();
	int alloc_mem_usage = (int)((double)alloc_n_points / (double)(1024 * 1024)) * sizeof(double);
	
	cout << "  total number of points is " << n_points << " requiring " << mem_usage << " MByte RAM.\n";
	cout << "  allocated total number of points is " << alloc_n_points << " requiring " << alloc_mem_usage << " MByte RAM.\n";
	cout << "  Note that " << n_removed_points << " redundant points have been removed from the input data.\n\n";
}





#endif /*MISC_H_*/
