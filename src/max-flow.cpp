#include <R.h>
#include "max-flow.hpp"
#include <boost/graph/copy.hpp>
using namespace boost;

inline double Graph_mf::calcD() {
  using namespace boost;
  double sum =0.0;
  graph_traits < NetGraph >::edge_iterator ei, eend;
  for(tie(ei,eend) = edges(gdual); ei != eend; ei++) {
    int s = source(*ei,gdual);
    int t = target(*ei,gdual);
    std::pair<Edge, bool> e = edge(vertex(s,gdual),
				   vertex(t,gdual),gdual);
    double l = get(edge_weight,gdual,e.first);
    double c = get(edge_capacity,gdual,e.first);
    sum += l*c;
  }
  return sum;
}


void Graph_mf::  max_concurrent_flow(std::vector<mf_demand> &demands,
				     double e) {

  int num_dem = demands.size();
  boost::copy_graph((NetGraph)*this,gdual);
  for(int i=0 ; i<num_dem; i++) {
    demands[i].flow = 0;
  }
  int N = num_vertices(*this);
  int m = num_edges(*this);
  number_flows = 0;
  double delta = pow(double(m) / (1.0 - e),-1.0/e);

  graph_traits < NetGraph >::edge_iterator ei, eend;

  for(tie(ei,eend) = edges(gdual); ei != eend; ei++) {
    int s = source(*ei,gdual);
    int t = target(*ei,gdual);
    std::pair<Edge, bool> e = edge(vertex(s,gdual),
				   vertex(t,gdual),gdual);
    double c = get(edge_capacity,gdual,e.first);
    put(edge_weight,gdual,e.first,delta/c);
  }
	

  double D;
  D=calcD();
  int doubreq =  2*int(ceil(1.0/e * log(m/(1-e))/log(1+e)));
	
  // assuming doubreq is about the maximum number of phases
  // then we want to only update the progress bar every 1%
  int updatepb = int(ceil(doubreq / 100.0));
  int phases =0;
  totalphases =0;
  
  std::vector<Vertex> penult(N);
  std::vector<double> dist(N);
  std::vector<mf_demand>::iterator vi,ve;


  std::vector<int> demand_index(num_dem);
  for(int i=0; i<num_dem; i++) {
    demand_index[i]=i;
  }

  // phases
  while(D < 1.0) {
    if(phases > doubreq) {
      for(vi=demands.begin(); vi < demands.end(); vi++) {
	vi->demand = vi->demand * 2;
      }
      phases = 0;
    }
    // steps
    random_shuffle(demand_index.begin(),demand_index.end());
	  
    for(int j=0; j<num_dem;j++) {
      int i= demand_index[j];
      mf_demand demand=demands[i];
      Vertex source = demand.source;
      Vertex sink = demand.sink;
      Vertex f,p;

      //iterations
      while( D < 1.0 && demand.demand > 0) {
	dijkstra_shortest_paths(gdual, source,
				predecessor_map(&penult[0]).distance_map(&dist[0]));

	// go through the path (backwards) and find minimum capacity
	f = sink;
	p = penult[f];
	std::pair<Edge, bool> ed = edge(p,f,gdual);

	double mincap=get(edge_capacity,gdual,ed.first);
	f = p;
	p = penult[p];
	double w =get(edge_weight,gdual,ed.first);
		  
	while(f != source) {
	  ed = edge(p,f,gdual);
	  double cap =get(edge_capacity,gdual,ed.first);
	  double w =get(edge_weight,gdual,ed.first);
	  mincap = cap < mincap? cap : mincap;
	  f = p;
	  p = penult[p];
	}

	// now we have the maximum flow we can push through this
	// step, and update demand (will add flow later)
	mincap = demand.demand < mincap ? demand.demand : mincap;
	demand.demand = demand.demand - mincap;

	// update each edge length = length (1 + (e*mincap) / capacity_e)
	// again go though the path backwards
	f = sink;
	p = penult[f];
	ed = edge(p,f,gdual);
		  
	double c=get(edge_capacity,gdual,ed.first);
	w =get(edge_weight,gdual,ed.first);
	w = w * (1 + (e * mincap) / c);
	put(edge_weight,gdual,ed.first,w);

	// note this records the path backwards! (slightly faster
	// when creating the path vector so no real reason to change
	std::list<Vertex> path;
	path.push_front(f);
	path.push_front(p);

	f = p;
	p = penult[p];

				 
	while(f != source) {
	  ed = edge(p,f,gdual);
	  c =get(edge_capacity,gdual,ed.first);
	  w =get(edge_weight,gdual,ed.first);
	  w = w * (1 + (e * mincap) / c);
	  put(edge_weight,gdual,ed.first,w);
	  path.push_front(p);
	  f = p;
	  p = penult[p];
	}
	
	if (demands[i].path_flow_map.find(path) == demands[i].path_flow_map.end()) {
	  demands[i].path_flow_map[path] = mincap;
	  number_flows++;
	} else {
	  demands[i].path_flow_map[path] += mincap;
	}


	demands[i].flow += mincap;

	D=calcD();
      }
    }
    phases++;
    totalphases++;
  }
  
}
