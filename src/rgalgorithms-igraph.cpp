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


#define IGRAPH_THREAD_LOCAL 

extern IGRAPH_THREAD_LOCAL igraph_interruption_handler_t 
  *igraph_i_interruption_handler;

#define IGRAPH_ALLOW_INTERRUPTION() \
       do { \
       if (igraph_i_interruption_handler) { if (igraph_allow_interruption(NULL) != IGRAPH_SUCCESS) return IGRAPH_INTERRUPTED; \
       } } while (0)


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

typedef struct igraph_2wheap_t {
  long int size;
  igraph_vector_t data;
  igraph_vector_long_t index;
  igraph_vector_long_t index2;
} igraph_2wheap_t;

int igraph_2wheap_init(igraph_2wheap_t *h, long int size);
void igraph_2wheap_destroy(igraph_2wheap_t *h);
int igraph_2wheap_clear(igraph_2wheap_t *h);
int igraph_2wheap_push_with_index(igraph_2wheap_t *h, 
				  long int idx, igraph_real_t elem);
igraph_bool_t igraph_2wheap_empty(const igraph_2wheap_t *h);
long int igraph_2wheap_size(const igraph_2wheap_t *h);
long int igraph_2wheap_max_size(const igraph_2wheap_t *h);
igraph_real_t igraph_2wheap_max(const igraph_2wheap_t *h);
long int igraph_2wheap_max_index(const igraph_2wheap_t *h);
igraph_real_t igraph_2wheap_deactivate_max(igraph_2wheap_t *h);
igraph_bool_t igraph_2wheap_has_elem(const igraph_2wheap_t *h, long int idx);
igraph_bool_t igraph_2wheap_has_active(const igraph_2wheap_t *h, long int idx);
igraph_real_t igraph_2wheap_get(const igraph_2wheap_t *h, long int idx);
igraph_real_t igraph_2wheap_delete_max(igraph_2wheap_t *h);
igraph_real_t igraph_2wheap_delete_max_index(igraph_2wheap_t *h, long int *idx);
int igraph_2wheap_modify(igraph_2wheap_t *h, long int idx, igraph_real_t elem);
int igraph_2wheap_check(igraph_2wheap_t *h);

#define IGRAPH_ALLOW_INTERRUPTION() \
       do { \
       if (igraph_i_interruption_handler) { if (igraph_allow_interruption(NULL) != IGRAPH_SUCCESS) return IGRAPH_INTERRUPTED; \
       } } while (0)


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

int igraph_get_shortest_paths_dijkstra_mod(const igraph_t *graph,
                                       igraph_vector_ptr_t *vertices,
				       igraph_vector_ptr_t *edges,
				       igraph_integer_t from,
				       igraph_vs_t to,
				       const igraph_vector_t *weights,
				       igraph_neimode_t mode) {
  /* Implementation details. This is the basic Dijkstra algorithm, 
     with a binary heap. The heap is indexed, i.e. it stores not only
     the distances, but also which vertex they belong to. The other
     mapping, i.e. getting the distance for a vertex is not in the
     heap (that would by the double-indexed heap), but in the result
     matrix.

     Dirty tricks:
     - the opposite of the distance is stored in the heap, as it is a
       maximum heap and we need a minimum heap.
     - we don't use IGRAPH_INFINITY in the distance vector during the
       computation, as IGRAPH_FINITE() might involve a function call 
       and we want to spare that. So we store distance+1.0 instead of 
       distance, and zero denotes infinity.
     - `parents' assigns the predecessors of all vertices in the
       shortest path tree to the vertices. In this implementation, the
       vertex ID + 1 is stored, zero means unreachable vertices.
  */
  
  long int no_of_nodes=igraph_vcount(graph);
  long int no_of_edges=igraph_ecount(graph);
  igraph_vit_t vit;
  igraph_2wheap_t Q;
  igraph_lazy_inclist_t inclist;
  igraph_vector_t dists;
  long int *parents;
  igraph_bool_t *is_target;
  long int i,to_reach;

  if (!weights) {
    return igraph_get_shortest_paths(graph, vertices, edges, from, to, mode);
  }
  
  if (igraph_vector_size(weights) != no_of_edges) {
    IGRAPH_ERROR("Weight vector length does not match", IGRAPH_EINVAL);
  }
  if (igraph_vector_min(weights) < 0) {
    IGRAPH_ERROR("Weight vector must be non-negative", IGRAPH_EINVAL);
  }

  IGRAPH_CHECK(igraph_vit_create(graph, to, &vit));
  IGRAPH_FINALLY(igraph_vit_destroy, &vit);

  if (vertices && IGRAPH_VIT_SIZE(vit) != igraph_vector_ptr_size(vertices)) {
    IGRAPH_ERROR("Size of `vertices' and `to' should match", IGRAPH_EINVAL);
  }
  if (edges && IGRAPH_VIT_SIZE(vit) != igraph_vector_ptr_size(edges)) {
    IGRAPH_ERROR("Size of `edges' and `to' should match", IGRAPH_EINVAL);
  }

  IGRAPH_CHECK(igraph_2wheap_init(&Q, no_of_nodes));
  IGRAPH_FINALLY(igraph_2wheap_destroy, &Q);
  IGRAPH_CHECK(igraph_lazy_inclist_init(graph, &inclist, mode));
  IGRAPH_FINALLY(igraph_lazy_inclist_destroy, &inclist);

  IGRAPH_VECTOR_INIT_FINALLY(&dists, no_of_nodes);
  //mjreed added
  for(i=0;i<no_of_nodes;i++) VECTOR(dists)[i]=-1.0;
  // end mjreed added

  parents = igraph_Calloc(no_of_nodes, long int);
  if (parents == 0) IGRAPH_ERROR("Can't calculate shortest paths", IGRAPH_ENOMEM);
  IGRAPH_FINALLY(igraph_free, parents);
  is_target = igraph_Calloc(no_of_nodes, igraph_bool_t);
  if (is_target == 0) IGRAPH_ERROR("Can't calculate shortest paths", IGRAPH_ENOMEM);
  IGRAPH_FINALLY(igraph_free, is_target);

  /* Mark the vertices we need to reach */
  to_reach=IGRAPH_VIT_SIZE(vit);
  for (IGRAPH_VIT_RESET(vit); !IGRAPH_VIT_END(vit); IGRAPH_VIT_NEXT(vit)) {
    if (!is_target[ (long int) IGRAPH_VIT_GET(vit) ]) {
      is_target[ (long int) IGRAPH_VIT_GET(vit) ] = 1;
    } else {
      to_reach--;		/* this node was given multiple times */
    }
  }

  //VECTOR(dists)[(long int)from] = 1.0;	/* zero distance */
  parents[(long int)from] = 0;
  VECTOR(dists)[(long int)from] = 0.0;	/* zero distance */
  igraph_2wheap_push_with_index(&Q, from, 0);
    
  while (!igraph_2wheap_empty(&Q) && to_reach > 0) {
    long int nlen, minnei=igraph_2wheap_max_index(&Q);
    igraph_real_t mindist=-igraph_2wheap_delete_max(&Q);
    igraph_vector_t *neis;

    IGRAPH_ALLOW_INTERRUPTION();

    if (is_target[minnei]) {
      is_target[minnei] = 0;
	  to_reach--;
	}

    /* Now check all neighbors of 'minnei' for a shorter path */
    neis=igraph_lazy_inclist_get(&inclist, minnei);
    nlen=igraph_vector_size(neis);
    for (i=0; i<nlen; i++) {
      long int edge=VECTOR(*neis)[i];
      long int tto=IGRAPH_OTHER(graph, edge, minnei);
      igraph_real_t altdist=mindist + VECTOR(*weights)[edge];
      igraph_real_t curdist=VECTOR(dists)[tto];
      //if (curdist==0) {
      if (curdist<0) {
        /* This is the first non-infinite distance */
        //VECTOR(dists)[tto] = altdist+1.0;
	VECTOR(dists)[tto] = altdist;
        parents[tto] = edge+1;
        IGRAPH_CHECK(igraph_2wheap_push_with_index(&Q, tto, -altdist));
	//} else if (altdist < curdist-1) {
      } else if (altdist < curdist) {
	      /* This is a shorter path */
        //VECTOR(dists)[tto] = altdist+1.0;
	VECTOR(dists)[tto] = altdist;
	parents[tto] = edge+1;
        IGRAPH_CHECK(igraph_2wheap_modify(&Q, tto, -altdist));
      }
    }
  } /* !igraph_2wheap_empty(&Q) */

  if (to_reach > 0) IGRAPH_WARNING("Couldn't reach some vertices");

  /* Reconstruct the shortest paths based on vertex and/or edge IDs */
  if (vertices || edges) {
    for (IGRAPH_VIT_RESET(vit), i=0; !IGRAPH_VIT_END(vit); IGRAPH_VIT_NEXT(vit), i++) {
      long int node=IGRAPH_VIT_GET(vit);
      igraph_vector_t *vvec=0, *evec=0;
      if (vertices) {
	vvec=(igraph_vector_t*)VECTOR(*vertices)[i];
	igraph_vector_clear(vvec);
      }
      if (edges) {
	evec=(igraph_vector_t*)VECTOR(*edges)[i];
	igraph_vector_clear(evec);
      }
      
      IGRAPH_ALLOW_INTERRUPTION();
      
      if (parents[node]>0) {
	long int size=0;
	long int act=node;
	long int edge;
	while (parents[act]) {
	  size++;
	  edge=parents[act]-1;
	  act=IGRAPH_OTHER(graph, edge, act);
	}
	if (vvec) { 
	  IGRAPH_CHECK(igraph_vector_resize(vvec, size+1)); 
	  VECTOR(*vvec)[size]=node;
	}
	if (evec) {
	  IGRAPH_CHECK(igraph_vector_resize(evec, size));
	}
	act=node;
	while (parents[act]) {
	  edge=parents[act]-1;
	  act=IGRAPH_OTHER(graph, edge, act);
	  size--;
	  if (vvec) { VECTOR(*vvec)[size]=act; }
	  if (evec) { VECTOR(*evec)[size]=edge; }
	}
      }
    }
  }
  
  igraph_lazy_inclist_destroy(&inclist);
  igraph_2wheap_destroy(&Q);
  igraph_vector_destroy(&dists);
  igraph_Free(is_target);
  igraph_Free(parents);
  igraph_vit_destroy(&vit);
  IGRAPH_FINALLY_CLEAN(6);
  
  return 0;
}



__END_DECLS
