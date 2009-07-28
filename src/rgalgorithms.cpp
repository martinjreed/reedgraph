#include "RBGL.hpp"
#include <boost/graph/dijkstra_shortest_paths.hpp>

// MATRIX(ptr,i,j,nc) {ptr[i * nc +j]}

extern "C" {
  typedef boost::graph_traits < Graph_dd >::edge_descriptor Edge;
  typedef boost::graph_traits < Graph_dd >::vertex_descriptor Vertex;

  // extract path from dijkstra_shortest_paths predecessor map
  // s source node, f destination node
  // p is predecessor map
  // d is distance map
  // It does not check for an infinite distance!
  void rg_cint_extractPath(Vertex s, 
						   Vertex f, 
						   std::vector<Vertex>& p, 
						   std::vector<double>& d,
						   std::vector<Vertex>& path) {
	using namespace boost;
	if( s == f ) {
	  path.resize(2);
	  path[0]=p[s];
	  path[1]=p[f];
	} else {
	  path.resize(1);
	  path[0]=f;
	  while(path[0] != s) {
		path.insert(path.begin(),p[f]);
		f=p[f];
	  }
	}
  }

  // set the weight along a path
  // g is graph
  // path is a vector of vertices (along path from source to dest)
  // val is the value to set along the path
  void rg_cint_setWeightOnPath(Graph_dd& g,std::vector<Vertex>& path,double val){
	using namespace boost;
	for(int i=0; i< path.size()-1; i++) {
	  
	  std::pair<Edge, bool> e = edge(vertex(path[i],g),vertex(path[i+1],g),g);
	  put(edge_weight,g,e.first,val);
	}
  }

  // calculate the number of edge disjoint paths in a graph
  // this is a primary function called by an R stub
  SEXP rg_num_edge_disjoint_paths_c(SEXP num_verts_in,
									SEXP num_edges_in,
									SEXP R_edges_in,
									SEXP R_weights_in) {
	using namespace boost;

	Graph_dd g(num_verts_in, num_edges_in, R_edges_in, R_weights_in);
	int N = num_vertices(g);
	std::vector<Vertex> p(N);
	std::vector<double> d(N);

	// set up a return matrix and set it to zero
	SEXP count;
	PROTECT(count = allocMatrix(INTSXP,N,N));
	int* countp;
	countp = INTEGER(count);
	for(int i=0;i<N;i++)
	  for(int j=0;j<N;j++)
		countp[i + j*N] = 0;

	graph_traits < Graph_dd >::vertex_iterator vi, viend;
	graph_traits < Graph_dd >::vertex_iterator vj, vjend;

	
	for (tie(vi, viend) = vertices(g); vi != viend; ++vi) {
	  dijkstra_shortest_paths(g, vertex(*vi, g),
							  predecessor_map(&p[0]).distance_map(&d[0]));
	  for (tie(vj,vjend) = vertices(g); vj != vjend; ++vj) {
		if (*vi != *vj ) {
		  if ( d[*vi] == std::numeric_limits<double>::max() ) {
			Rprintf("infinite for %d %d\n",*vi,*vj);
			countp[(*vi) + (*vj)*N]=0;
		  } else {
			Graph_dd gtmp(num_verts_in, num_edges_in, R_edges_in, R_weights_in);
			std::vector<Vertex> ptmp(p);
			std::vector<double> dtmp(d);
			while( dtmp[*vj] != std::numeric_limits<double>::max() ) {
			  countp[(*vi) + (*vj)*N]++;
			  std::vector<Vertex> path;
			  rg_cint_extractPath(*vi,*vj,ptmp,dtmp,path);
			  rg_cint_setWeightOnPath(gtmp,path,std::numeric_limits<double>::max());
			  dijkstra_shortest_paths(gtmp, vertex(*vi, gtmp),
									  predecessor_map(&ptmp[0]).distance_map(&dtmp[0]));
			  
			}
		  }
		  
		}
	  }
	}
	UNPROTECT(1);
	return(count);
  }
}
