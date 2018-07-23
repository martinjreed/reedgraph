/*
    Copyright (C) 2012 Martin J Reed              martin@reednet.org.uk
    University of Essex, Colchester, UK

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along
    with this program; if not, write to the Free Software Foundation, Inc.,
    51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 */

#include "max-flow.hpp"
#include <float.h>
#include <boost/graph/copy.hpp>
#include <boost/lexical_cast.hpp>
using namespace boost;
#include <R.h>
#include <Rcpp.h>
#include <Rinternals.h>
#include <string>
#include <sstream>
#include <iostream>

template <typename T>
std::string to_string(T value)
{
  //create an output string stream
  std::ostringstream os ;
  
  //throw the value into the string stream
  os << value ;
  
  //convert the string stream into a string and return
  return os.str() ;
}

inline double Graph_mf::calcD() {
  using namespace boost;
  double sum =0.0;
  graph_traits < NetGraph >::edge_iterator ei, eend;
  for(tie(ei,eend) = edges(gdual); ei != eend; ei++) {
    double l = get(edge_weight,gdual,*ei);
    double c = get(edge_capacity,gdual,*ei);
    //int s = source(*ei,gdual);
    //int t = target(*ei,gdual);
    //Rprintf("%ld -> %ld, %lg, %lg\n",s,t,l,c);
    sum += l*c;
  }
  return sum;
}

double Graph_mf::assign_gflow(std::vector<mf_demand> &demands) {
  graph_traits < NetGraph >::edge_iterator ei, eend;
  for(tie(ei,eend) = edges(gflow); ei != eend; ei++) {
    //int s = source(*ei,gflow);
    //int t = target(*ei,gflow);
    //Rprintf("%ld -> %ld\n",s,t);
    put(edge_weight,gflow,*ei,0.0);
  }
  std::vector<mf_demand>::iterator di;

  for(di=demands.begin(); di != demands.end(); di++) {
    std::map<const std::list<Vertex>,double>::iterator mi;
    //Rprintf("(*di).path_flow_map.size()=%ld\n",(*di).path_flow_map.size());
    for(mi=(*di).path_flow_map.begin() ; mi != (*di).path_flow_map.end(); mi++) {
      std::list<Vertex>::const_iterator lh=(mi->first).begin();
      std::list<Vertex>::const_iterator lt=(mi->first).begin();
      lh++;
      for(;lh != mi->first.end(); lh++, lt++) {
	//Rprintf("%ld -> %ld\n",*lt,*lh);
	std::pair<Edge, bool> ed = edge(*lt,*lh,gflow);
	double w =get(edge_weight,gflow,ed.first);
	put(edge_weight,gflow,ed.first,w+mi->second);
      }
    }
  }
  gamma = DBL_MAX;
  for(tie(ei,eend) = edges(gflow); ei != eend; ei++) {
    double w = get(edge_weight,gflow,*ei);
    double c = get(edge_capacity,gflow,*ei);
    gamma = (c-w)/c < gamma ? (c-w)/c : gamma;
  }
  return(gamma);
}

inline double Graph_mf::calcLambda(std::vector<mf_demand> &demands) {
  lambda = DBL_MAX;
  gamma = 0;
  std::vector<mf_demand>::iterator di;
  for(di=demands.begin(); di != demands.end(); di++) {
    lambda = di->flow / di->demand < lambda ? di->flow / di->demand : lambda;
  }
  gamma= 1.0 - 1.0/lambda;
  return(lambda);
}
void Graph_mf:: sp_concurrent_flow(std::vector<mf_demand> &demands) {

  std::vector<mf_demand>::iterator di;
  for(di=demands.begin(); di != demands.end(); di++) {
    di->flow = 0;
    di->path_flow_map.erase(di->path_flow_map.begin(),di->path_flow_map.end());
  }

  int num_dem = demands.size();
  gflow.clear();
  boost::copy_graph((NetGraph)*this,gflow);
  gdual.clear();
  boost::copy_graph((NetGraph)*this,gdual);
  for(int i=0 ; i<num_dem; i++) {
    demands[i].flow = 0;
  }
  int N = num_vertices(*this);
  number_flows = 0;

  graph_traits < NetGraph >::edge_iterator ei, eend;
  // assign unit weigths to all edges
  for(tie(ei,eend) = edges(gdual); ei != eend; ei++) {
    int s = source(*ei,gdual);
    int t = target(*ei,gdual);
    std::pair<Edge, bool> e = edge(vertex(s,gdual),
				   vertex(t,gdual),gdual);
    put(edge_weight,gdual,e.first,1.0);
    e = edge(vertex(s,gflow),
	     vertex(t,gflow),gflow);
    put(edge_weight,gflow,e.first,0.0);
  }
  // for each demand
  for(int i=0; i<num_dem;i++) {
    mf_demand demand = demands[i];

    std::vector<Vertex> penult(N);
    std::vector<double> dist(N);
    std::vector<mf_demand>::iterator vi,ve;
    Vertex source = demand.source;
    Vertex sink = demand.sink;
    Vertex f,p;

    // calculated shortest path
    dijkstra_shortest_paths(gdual, source,
			    predecessor_map(&penult[0]).distance_map(&dist[0]));
    std::list<Vertex> path;

    f = sink;
    p = penult[f];

    std::pair<Edge, bool> ed = edge(p,f,gdual);

    path.push_front(f);
    path.push_front(p);
    f = p;
    p = penult[p];
    while(f != source) {
      ed = edge(p,f,gflow);
      double w =get(edge_weight,gflow,ed.first);
      put(edge_weight,gflow,ed.first,w+demand.demand);
      path.push_front(p);
      f = p;
      p = penult[p];
    }
    
    demands[i].path_flow_map[path] = demand.demand;
    demands[i].flow = demand.demand;
  }

  // calculate lambda
  lambda = DBL_MAX;
  gamma = 0;
  for(tie(ei,eend) = edges(gflow); ei != eend; ei++) {
    double w = get(edge_weight,gflow,*ei);
    double c = get(edge_capacity,gflow,*ei);
    lambda = c/w < lambda ? c/w : lambda;
  }
  gamma= 1.0 - 1.0/lambda;
}

double Graph_mf::calcBeta(std::vector<mf_demand> &demands) {
  double Alpha=0;
  int N = num_vertices(*this);

  std::vector<mf_demand>::iterator di;
  std::vector<Vertex> penult(N);
  std::vector<double> dist(N);

  for(di=demands.begin(); di != demands.end(); di++) {
    dijkstra_shortest_paths(gdual, di->source,
			    predecessor_map(&penult[0]).distance_map(&dist[0]));
    Alpha += (*di).demand * dist[di->sink];
  }
  double D= calcD();
  return(D/Alpha);
}


void Graph_mf::rescale_demands(std::vector<mf_demand> &demands,double scale) {
  std::vector<mf_demand>::iterator di;
  std::map<const std::list<Vertex>,double>::iterator mi;
  for(di=demands.begin(); di != demands.end(); di++) {
    (*di).demand *= scale;
  }
}

void Graph_mf::rescale_demands_flows(std::vector<mf_demand> &demands,double scale) {

  std::vector<mf_demand>::iterator di;
  std::map<const std::list<Vertex>,double>::iterator mi;
  for(di=demands.begin(); di != demands.end(); di++) {
    (*di).flow *= scale;
    for(mi=(*di).path_flow_map.begin() ; mi != (*di).path_flow_map.end(); mi++) {
      (*mi).second *= scale;
    }
  }
}

void Graph_mf:: max_concurrent_flow_prescaled(std::vector<mf_demand> &demands,
					      double e) {

  std::vector<mf_demand> save_demands = demands;
  long num_dem=demands.size();
  sp_concurrent_flow(demands);
  rescale_demands(demands,lambda);
  // Beta might be very high so obtain 2-approximate solution
  // (1+w) = 2 = (1-e)^-3
  double e2= 1- pow(2.0,-1.0/3.0);
  max_concurrent_flow(demands,e2);
  rescale_demands(demands,beta/2);

  max_concurrent_flow(demands,e);

  for(int i=0; i<num_dem; i++) {
    demands[i].demand = save_demands[i].demand;
  }
  lambda = calcLambda(demands);
  beta = calcBeta(demands);
  assign_gflow(demands);
}

void Graph_mf:: min_congestion_flow(std::vector<mf_demand> &demands,
				     double e) {
  max_concurrent_flow_prescaled(demands,e);
  rescale_demands_flows(demands,1/lambda);
  assign_gflow(demands);
}

void Graph_mf::  max_concurrent_flow(std::vector<mf_demand> &demands,
				     double e) {

  std::vector<mf_demand>::iterator di;
  for(di=demands.begin(); di != demands.end(); di++) {
    di->flow = 0;
    di->path_flow_map.erase(di->path_flow_map.begin(),di->path_flow_map.end());
  }

  int num_dem = demands.size();
  gflow.clear();
  boost::copy_graph((NetGraph)*this,gflow);
  gdual.clear();
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
	
  // assuming doubreq is about the maximum number of phases
  // then we want to only update the progress bar every 1%
  //int doubreq =  2*int(ceil(1.0/e * log(m/(1-e))/log(1+e)));
  //int updatepb = int(ceil(doubreq / 100.0));
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
    //Rprintf("D %lg \n",D);
    // this doubling could be included but it needs "undoubling" for
    // demands and flows at the end if prescaling is done this should
    // not be necessary
    /*
    if(phases > doubreq) {
      for(vi=demands.begin(); vi < demands.end(); vi++) {
	vi->demand = vi->demand * 2;
      }
      phases = 0;
      }*/
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
  
  double scalef = 1.0 / (log(1.0/delta) / log(1+e) );
  rescale_demands_flows(demands,scalef);
  lambda = calcLambda(demands);
  beta = calcBeta(demands);
  assign_gflow(demands);
  
}

void Graph_mf::max_concurrent_flow_int_prescaled(std::vector<mf_demand> &demands,
						 std::vector<mf_demand> &best_demands,
						 double e) {
  std::vector<mf_demand> unscaledDemands = demands;
  long num_dem=demands.size();
  sp_concurrent_flow(demands);
  rescale_demands(demands,lambda);
  // Beta might be very high so obtain 2-approximate solution
  // (1+w) = 2 = (1-e)^-3
  double e2= 1- pow(2.0,-1.0/3.0);
  max_concurrent_flow(demands,e2);
  //rescale_demands(demands,beta/2);
  rescale_demands(demands,beta*2);
  max_concurrent_flow_int(unscaledDemands,demands,best_demands,e);
  for(int i=0; i<num_dem; i++) {
    demands[i].demand = unscaledDemands[i].demand;
  }
  // calcLambda calculates gamma, but this is incorrect if prescaling, so
  // instead keep the value we had found in the algorithm.
  double saveGamma = gamma;
  lambda = calcLambda(demands);
  gamma = saveGamma;
  //Rprintf("lambda=%lg\n",lambda);
  beta = calcBeta(demands);
  //Rprintf("beta=%lg\n",beta);
  //assign_gflow(demands);

}

void Graph_mf::max_concurrent_flow_int(std::vector<mf_demand> &unscaledDemands,
					 std::vector<mf_demand> &demands,
					 std::vector<mf_demand> &best_demands,
					 double e) {
  using namespace std;
  gamma = -DBL_MAX;
  double best_mean_gamma = -DBL_MAX;
  // std::vector<mf_demand>::iterator di;
  // for(di=demands.begin(); di != demands.end(); di++) {
  //   di->flow = 0;
  //   di->path_flow_map.erase(di->path_flow_map.begin(),di->path_flow_map.end());
  // }
  assigned = 0;
  int num_dem = demands.size();

  gflow.clear();
  boost::copy_graph((NetGraph)*this,gflow);
  gdual.clear();
  boost::copy_graph((NetGraph)*this,gdual);
  for(int i=0 ; i<num_dem; i++) {
    demands[i].flow = 0;
    unscaledDemands[i].flow=0;
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
	
  // assuming doubreq is about the maximum number of phases
  // then we want to only update the progress bar every 1%
  //int doubreq =  2*int(ceil(1.0/e * log(m/(1-e))/log(1+e)));
  //int updatepb = int(ceil(doubreq / 100.0));
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
    //Rprintf("D=%lg\n",D);
    // this doubling could be included but it needs "undoubling" for
    // demands and flows at the end if prescaling is done this should
    // not be necessary
    /*
    if(phases > doubreq) {
      for(vi=demands.begin(); vi < demands.end(); vi++) {
	vi->demand = vi->demand * 2;
      }
      phases = 0;
      }*/
    // steps
    random_shuffle(demand_index.begin(),demand_index.end());
    NetGraph glimit;
    boost::copy_graph(gdual,glimit);
    for(tie(ei,eend) = edges(glimit); ei != eend; ei++) {
      put(edge_weight,glimit,*ei,0.0);	    
    }
    double local_gamma = DBL_MAX;

    long count=0;
    std::vector<mf_demand> int_demands(num_dem);
    for(int j=0; j<num_dem;j++) {
      int i= demand_index[j];
      int_demands[i].path_flow_map.clear();
      int_demands[i].flow = 0;
      mf_demand demand=demands[i];
      mf_demand unscaledDemand=unscaledDemands[i];
      Vertex vsource = demand.source;
      Vertex vsink = demand.sink;
      Vertex f,p;
      //Rprintf("Doing demand %d, %ld -> %ld\n",i,vsource,vsink);

      //iterations
      // we do not push through demand until it is full
      // rather push just once, hence change here
      //while( D < 1.0 && demand.demand > 0) {
	//Rprintf("demand.demand=%lg\n",demand.demand);
	NetGraph gtmp;
	boost::copy_graph(gdual,gtmp);
	for(tie(ei,eend) = edges(gtmp); ei != eend; ei++) {
	  int s = source(*ei,gtmp);
	  int t = target(*ei,gtmp);
	  std::pair<Edge, bool> ed = edge(vertex(s,gtmp),
					 vertex(t,gtmp),gtmp);
	  double c = get(edge_capacity,gtmp,ed.first);
	  std::pair<Edge, bool> et = edge(vertex(s,glimit),
					  vertex(t,glimit),glimit);
	  double w = get(edge_weight,glimit,et.first);
	  //Rprintf("%ld->%ld edge weight=%lg, for cap %lg, dem %lg\n",s,t,w,c,demand.demand);
	  if(c-w < unscaledDemand.demand) {
	    //Rprintf("%ld->%ld set to inf\n",s,t);
	    put(edge_weight,gtmp,ed.first,DBL_MAX);	    
	  }

	}
	
	dijkstra_shortest_paths(gtmp, vsource,
				predecessor_map(&penult[0]).distance_map(&dist[0]));
	//Rprintf("Penult:\n");
	for(int k=0;k<penult.size();k++) {
	  //Rprintf("%ld,",penult[k]);
	}
	//Rprintf("\nDist:\n");
	for(int k=0;k<dist.size();k++) {
	  //Rprintf("%lg,",dist[k]);
	}
	//Rprintf("distance to sink=%lg\n",dist[vsink]);
	if(dist[vsink] == DBL_MAX) {
	  //Rprintf("infinite\n");
	    
	  continue;
	  //demand.demand = 0;
	}
	count++;
	// go through the path (backwards) and find minimum capacity
	f = vsink;
	p = penult[f];
	std::pair<Edge, bool> ed = edge(p,f,gtmp);
	std::pair<Edge, bool> el = edge(p,f,glimit);

	double mincap=get(edge_capacity,gdual,ed.first);
	f = p;
	p = penult[p];
	double w =get(edge_weight,gdual,ed.first);
		  
	while(f != vsource) {
	  ed = edge(p,f,gdual);
	  double cap =get(edge_capacity,gdual,ed.first);
	  mincap = cap < mincap? cap : mincap;
	  f = p;
	  p = penult[p];
	}

	// now we have the maximum flow we can push through this
	// step, and update demand (will add flow later)
	// due to earlier check mincap = demand.demand always!
	// if this is not true, need to adapt later demands update!
	mincap = demand.demand < mincap ? demand.demand : mincap;
	demand.demand = demand.demand - mincap;

	// update each edge length = length (1 + (e*mincap) / capacity_e)
	// again go though the path backwards
	f = vsink;
	p = penult[f];
	ed = edge(p,f,gdual);
	el = edge(p,f,glimit);

	double c=get(edge_capacity,gdual,ed.first);
	w =get(edge_weight,gdual,ed.first);
	w = w * (1 + (e * mincap) / c);
	put(edge_weight,gdual,ed.first,w);
	double flow = get(edge_weight,glimit,el.first);
	put(edge_weight,glimit,el.first,flow+unscaledDemand.demand);
	// This bit will capture the gamma as we go through the edges
	double tmpGamma = 1.0 - (flow+unscaledDemand.demand)/c;
	if ( local_gamma > tmpGamma )
	  local_gamma = tmpGamma;
	//Rprintf("edge %ld|%ld flow:%lg+%lg\n",p,f,flow,mincap);
	std::list<Vertex> path;
	path.push_front(f);
	path.push_front(p);

	f = p;
	p = penult[p];

				 
	while(f != vsource) {
	  ed = edge(p,f,gdual);	  el = edge(p,f,glimit);
	  c =get(edge_capacity,gdual,ed.first);
	  w =get(edge_weight,gdual,ed.first);
	  w = w * (1 + (e * mincap) / c);
	  put(edge_weight,gdual,ed.first,w);
	  flow = get(edge_weight,glimit,el.first);
	  put(edge_weight,glimit,el.first,flow+unscaledDemand.demand);
	  // This bit will capture the gamma as we go through the edges
	  double tmpGamma = 1.0 - (flow+unscaledDemand.demand)/c;
	  if ( local_gamma > tmpGamma )
	    local_gamma = tmpGamma;
	  //Rprintf("edge %ld|%ld flow:%lg+%lg\n",p,f,flow,mincap);
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
	int_demands[i].path_flow_map.clear();
	int_demands[i].path_flow_map[path] = unscaledDemand.demand;
	int_demands[i].flow = mincap;

	demands[i].flow += mincap;
	D=calcD();

    }
    //D=1.0;
    //Rprintf("Assigned (old %ld) %ld/%ld, gamma=(old:%lg)%lg\n",assigned,count,num_dem,gamma,local_gamma);
    
    if(count > assigned) {
      assigned = count;
      gamma=local_gamma;
      for(int i=0 ; i<num_dem; i++) {
	best_demands[i].path_flow_map.clear();
	best_demands[i].flow=0;
	std::map<const std::list<Vertex>, double>::iterator pmi, pmend;
	for(pmi=int_demands[i].path_flow_map.begin();
	    pmi != int_demands[i].path_flow_map.end(); pmi++) {
	  std::list<Vertex> path = pmi->first;
	  double val = pmi->second;
	  best_demands[i].flow=val;
	  best_demands[i].path_flow_map[path]=val;
	}
      }
    } else if (count == assigned && gamma < local_gamma) {
      //Rprintf("found better gamma\n");
      gamma = local_gamma;
      for(int i=0 ; i<num_dem; i++) {
	best_demands[i].path_flow_map.clear();
	best_demands[i].flow=0;
	std::map<const std::list<Vertex>, double>::iterator pmi, pmend;
	for(pmi=int_demands[i].path_flow_map.begin();
	    pmi != int_demands[i].path_flow_map.end(); pmi++) {
	  std::list<Vertex> path = pmi->first;
	  double val = pmi->second;
	  best_demands[i].flow=val;
	  best_demands[i].path_flow_map[path]=val;
	}
      }
      
    } else if (count == assigned && gamma == local_gamma) {
      //Rprintf("gamma == localgamma\n");
      // we have the same gamma, let us try to minimize mean gamma
      double mean_gamma = 0;

      for(tie(ei,eend) = edges(glimit); ei != eend; ei++) {
	double flow = get(edge_weight,glimit,*ei);
	double c=get(edge_capacity,glimit,*ei);
	mean_gamma += (1.0 - flow/c);
      }
      mean_gamma = mean_gamma/m;
      if (best_mean_gamma < mean_gamma ) {
	//Rprintf("Found better mean gamma %lg %lg\n",best_mean_gamma,mean_gamma);
	best_mean_gamma = mean_gamma;
	for(int i=0 ; i<num_dem; i++) {
	  best_demands[i].path_flow_map.clear();
	  best_demands[i].flow=0;
	  std::map<const std::list<Vertex>, double>::iterator pmi, pmend;
	  for(pmi=int_demands[i].path_flow_map.begin();
	      pmi != int_demands[i].path_flow_map.end(); pmi++) {
	    std::list<Vertex> path = pmi->first;
	    double val = pmi->second;
	    best_demands[i].flow=val;
	    best_demands[i].path_flow_map[path]=val;
	  }
	}
	
      }
    }
    
    phases++;
    totalphases++;
  }
  //Rprintf("Assigned %ld/%ld\n",assigned,num_dem);
  // for(int i=0 ; i<num_dem; i++) {
  //   Rprintf("Demand %ld flow=%lg\n",i+1,best_demands[i].flow);
  // }
  double scalef = 1.0 / (log(1.0/delta) / log(1+e) );
  rescale_demands_flows(demands,scalef);
  lambda = calcLambda(demands);
  beta = calcBeta(demands);
  assign_gflow(best_demands);
  
}
