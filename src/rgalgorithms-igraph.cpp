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

#include <R.h>
#include <igraph.h>
#include <vector>
#include <map>
#include <algorithm>
// This cannot be included as igraph definitions clash with Rcpp and boost
// hence why this file has been split off from rgalgorithms.cpp
//#include "rgalgorithms.h"
#undef __BEGIN_DECLS
#undef __END_DECLS
#ifdef __cplusplus
# define __BEGIN_DECLS extern "C" {
# define __END_DECLS }
#else
# define __BEGIN_DECLS /* empty */
# define __END_DECLS /* empty */
#endif




class rg_demandi {
public:
  double flow;
  long source;
  long sink;
  double demand;
  std::map< std::vector< long > , double > path_flow_map;
};

double calcDi(igraph_t* graph, igraph_vector_t* weights, std::vector<double>* capacities) {
  double sum =0.0;
  igraph_eit_t ieit;
  igraph_integer_t from, to;

  igraph_eit_create(graph,igraph_ess_all(IGRAPH_EDGEORDER_ID),&ieit);
  IGRAPH_EIT_RESET(ieit);
  while(!IGRAPH_EIT_END(ieit)) {
    igraph_integer_t edgeid = IGRAPH_EIT_GET(ieit);
    double l = VECTOR(*weights)[(long)edgeid];
    double c = (*capacities)[(long)edgeid];
    sum += l*c;
    IGRAPH_EIT_NEXT(ieit);
  }
  return sum;
}

int updateExplicitFlowi(rg_demandi& demand,
			igraph_vector_t* vec,
			double flow) {
  std::vector<long> path;
  int newpath=0;
  for(int k=0;k<igraph_vector_size(vec);k++) {
    path.push_back((long)VECTOR(*vec)[k]);
  }
	
  if (demand.path_flow_map.find(path) == demand.path_flow_map.end()) {
    demand.path_flow_map[path] = flow;
    newpath=1;
  } else {
    demand.path_flow_map[path] += flow;
  }
  return newpath;
}

__BEGIN_DECLS



int igraph_get_shortest_paths_dijkstra_mod(const igraph_t *graph,
                                       igraph_vector_ptr_t *vertices,
				       igraph_vector_ptr_t *edges,
				       igraph_integer_t from,
				       igraph_vs_t to,
				       const igraph_vector_t *weights,
					   igraph_neimode_t mode);


void igraph_rg_fleischer_max_concurrent_flow(std::vector<long> vedges,
					     std::vector<double> capacities,
					     std::vector<rg_demandi> &demandsi,
					     std::vector<double> &lengths,
					     double e,
					     long N,
					     long &totalphases) {
  igraph_t graph;
  igraph_vector_t ivedges;
  int num_dem = demandsi.size();
  int m = vedges.size()/2;
  double delta = pow(double(m) / (1.0 - e),-1.0/e);
  lengths.resize(m);
  igraph_vector_init(&ivedges,vedges.size());
  for(int i=0;i<vedges.size();i++) {
    VECTOR(ivedges)[i]=vedges[i];
  }

  igraph_vector_t weights;
  igraph_vector_init(&weights,m);
  for(int i=0;i<m;i++) {
    VECTOR(weights)[i]=1.0;
  }
  
  
  igraph_empty(&graph,N,true);
  igraph_add_edges(&graph,&ivedges,0);
  igraph_vector_destroy(&ivedges);

  igraph_eit_t ieit;

  igraph_eit_create(&graph,igraph_ess_all(IGRAPH_EDGEORDER_ID),&ieit);
  while(!IGRAPH_EIT_END(ieit)) {
    igraph_integer_t edgeid = IGRAPH_EIT_GET(ieit);
    VECTOR(weights)[(long)edgeid] = delta / capacities[(long)edgeid];
    IGRAPH_EIT_NEXT(ieit);
  }
  igraph_eit_destroy(&ieit);

  double Di;
  Di=calcDi(&graph,&weights,&capacities);

  int doubreq =  2*int(ceil(1.0/e * log(m/(1-e))/log(1+e)));
	
  // assuming doubreq is about the maximum number of phases
  // then we want to only update the progress bar every 1%
  int updatepb = int(ceil(doubreq / 100.0));
  int phases =0;
  totalphases =0;
  std::vector<int> demand_index(num_dem);

  for(int i=0; i<num_dem; i++) {
    demand_index[i]=i;
  }

  //random_shuffle(demand_index.begin(),demand_index.end());
  if(delta == 0) {
    Rprintf("Error delta=0\n");
    // fix this just to "bomb out"
    Di=1.0;
  }

  std::vector< std::pair<int, double> > costpair(num_dem);
  std::vector<rg_demandi>::iterator vi;
  igraph_vector_ptr_t vecs;

  // phases
  while(Di < 1.0) {
    //Rprintf("\n new phase\n Doing demand:");
    if(phases > doubreq) {
      //Rprintf("DOubling %d %d\n",doubreq,totalphases);
      for(vi=demandsi.begin(); vi < demandsi.end(); vi++) {
	vi->demand = vi->demand * 2;
      }
      phases = 0;
    }

    random_shuffle(demand_index.begin(),demand_index.end());

    for(int j=0; j<num_dem;j++) {
      int i;
      i= demand_index[j];

      rg_demandi demandi=demandsi[i];

      long sourcei = demandi.source;
      long sinki = demandi.sink;

      //iterations
      while( Di < 1.0 && demandi.demand > 0) {
	//Rprintf("%lg\n",Di);
	igraph_vector_ptr_init(&vecs, 1);
	for (int ik=0; ik<igraph_vector_ptr_size(&vecs); ik++) {
	  VECTOR(vecs)[ik] = calloc(1, sizeof(igraph_vector_t));
	  igraph_vector_init((igraph_vector_t*)VECTOR(vecs)[ik], 0);
	}
		  
	igraph_get_shortest_paths_dijkstra_mod(&graph, &vecs, NULL,
					   (igraph_integer_t)sourcei, 
					   igraph_vss_1((igraph_integer_t)sinki),
					   &weights,IGRAPH_OUT);
	igraph_vector_t* vec = (igraph_vector_t*)VECTOR(vecs)[0];
	int k=0;
	igraph_es_t es;
	igraph_es_pairs_small(&es,IGRAPH_DIRECTED,
			      (long int) VECTOR(*vec)[k],(long int)VECTOR(*vec)[k+1],
			      -1);
	igraph_eit_create(&graph,es,&ieit);
	double mincapi = capacities[(long)IGRAPH_EIT_GET(ieit)];
	igraph_es_destroy(&es);
	igraph_eit_destroy(&ieit);

	for(k++;k<igraph_vector_size(vec)-1;k++) {
	  igraph_es_pairs_small(&es,IGRAPH_DIRECTED,
				(long int) VECTOR(*vec)[k],(long int)VECTOR(*vec)[k+1],
				-1);
	  igraph_eit_create(&graph,es,&ieit);
	  double cap = capacities[(long)IGRAPH_EIT_GET(ieit)];
	  mincapi = cap < mincapi? cap : mincapi;
	  igraph_es_destroy(&es);
	  igraph_eit_destroy(&ieit);
	}
	// now we have the maximum flow we can push through this
	// step, and update demand (will add flow later)
	mincapi = demandi.demand < mincapi ? demandi.demand : mincapi;
	demandi.demand = demandi.demand - mincapi;

	// update each edge length = length (1 + (e*mincap) / capacity_e)
	k=0;
	igraph_es_pairs_small(&es,IGRAPH_DIRECTED,
			      (long int) VECTOR(*vec)[k],(long int)VECTOR(*vec)[k+1],
			      -1);
	igraph_eit_create(&graph,es,&ieit);
	long c = capacities[(long)IGRAPH_EIT_GET(ieit)];
	double w = (double)VECTOR(weights)[(long)IGRAPH_EIT_GET(ieit)];
	w = w * (1 + (e * mincapi) / c);
	VECTOR(weights)[(long)IGRAPH_EIT_GET(ieit)]=w;
	igraph_es_destroy(&es);
	igraph_eit_destroy(&ieit);

	for(k++;k<igraph_vector_size(vec)-1;k++) {
	  igraph_es_pairs_small(&es,IGRAPH_DIRECTED,
				(long int) VECTOR(*vec)[k],(long int)VECTOR(*vec)[k+1],
				-1);
	  igraph_eit_create(&graph,es,&ieit);
	  c = capacities[(long)IGRAPH_EIT_GET(ieit)];
	  w = VECTOR(weights)[(long)IGRAPH_EIT_GET(ieit)];
	  w = w * (1 + (e * mincapi) / c);
	  VECTOR(weights)[(long)IGRAPH_EIT_GET(ieit)]=w;
	  igraph_es_destroy(&es);
	  igraph_eit_destroy(&ieit);

	}

	updateExplicitFlowi(demandsi[i],vec,mincapi);
		  
	demandsi[i].flow += mincapi;

	Di=calcDi(&graph,&weights,&capacities);
	for (int ik=0; ik<igraph_vector_ptr_size(&vecs); ik++) {
	  igraph_vector_destroy((igraph_vector_t*)VECTOR(vecs)[ik]);
	  free(VECTOR(vecs)[ik]);
	}
	igraph_vector_ptr_destroy(&vecs);

      }
    }
    phases++;
    totalphases++;
  }
  for(int i=0;i<m;i++) {
    lengths[i]=VECTOR(weights)[i];
  }
  igraph_vector_destroy(&weights);
  igraph_destroy(&graph);

}





__END_DECLS
