#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/dijkstra_shortest_paths.hpp>
#include <Rmath.h>

#include <R.h>

#include "RBGL.hpp"
#include "Rcpp.h"


class Graph_rg: public boost::adjacency_list<boost::vecS, boost::vecS, boost::directedS, boost::property<boost::vertex_color_t,   boost::default_color_type>,  boost::property<boost::edge_weight_t, double,  boost::property<boost::edge_capacity_t, double> > >
{
  typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::directedS,
								 boost::property<boost::vertex_color_t, 
												 boost::default_color_type>,
								 boost::property<boost::edge_weight_t, double,
												 boost::property<boost::edge_capacity_t, double > > > Base;
  
public:
  typedef Base::graph_property_type graph_property_type;
  typedef Base::vertices_size_type vertices_size_type;
  typedef Base::edges_size_type edges_size_type;
  
  inline Graph_rg()
	: Base() { }
  inline Graph_rg(const graph_property_type& p)
	: Base(p) { }
  inline Graph_rg(const Base& x)
	: Base(x) { }
  inline Graph_rg(vertices_size_type num_vertices)
	: Base(num_vertices) { }
  inline Graph_rg(vertices_size_type num_vertices,
						  const graph_property_type& p)
	: Base(num_vertices, p) { }
  
  inline Graph_rg(SEXP num_verts_in,
				  SEXP num_edges_in,
				  SEXP R_edges_in,
				  SEXP R_weights_in,
				  SEXP capacities)
    : Base(Rcpp::as<int>(num_verts_in))
    {
        if (!Rf_isReal(R_weights_in)) error("R_weights_in should be Real");
        if (!Rf_isInteger(R_edges_in)) error("R_edges_in should be integer");
        if (!Rf_isReal(capacities)) error("capacities should be Real");
        int NE = Rf_asInteger(num_edges_in);
        int* edges_in = INTEGER(R_edges_in);
        if (Rf_isReal(R_weights_in)) {
		  double* weights_in = REAL(R_weights_in);
		  double* cap_in =REAL(capacities);
		  for (int i = 0;
			   i < NE ; 
			   i++, edges_in += 2, weights_in++, cap_in ++) {
			std::pair<edge_descriptor, bool> e =
			  boost::add_edge(*edges_in, *(edges_in+1),
							  *weights_in, *this);
			boost::put(boost::edge_capacity,*this,e.first,*cap_in);
		  }
        } 
    }
  
    inline Graph_rg(SEXP num_verts_in,
					SEXP num_edges_in,
					SEXP R_edges_in)
      : Base(Rcpp::as<int>(num_verts_in))
    {
	  if (!Rf_isInteger(R_edges_in)) error("R_edges_in should be integer");
	  int NE = Rf_asInteger(num_edges_in);
	  int* edges_in = INTEGER(R_edges_in);
	  for (int i = 0; i < NE ; i++, edges_in += 2) {
		boost::add_edge(*edges_in, *(edges_in+1), 1, *this);
	  }
    }
};

template <class T>
inline std::string to_string (const T& t)
{
  std::stringstream ss;
  ss << t;
  return ss.str();
}

inline std::string to_string(const int& t);


typedef boost::graph_traits < Graph_dd >::edge_descriptor Edge;
typedef boost::graph_traits < Graph_dd >::vertex_descriptor Vertex;

typedef boost::graph_traits < Graph_rg >::edge_descriptor rgEdge;
typedef boost::graph_traits < Graph_rg >::vertex_descriptor rgVertex;

  // these two probably obsolete, but keep for future sort of vector<pair>
struct less_cost {
public:
  bool operator()(const std::pair<int, double> a,
				  const std::pair<int, double> b) {
	return (a.second < b.second);
  }
};
struct great_cost {
public:
  bool operator()(const std::pair<int, double> a,
				  const std::pair<int, double> b) {
	return (a.second > b.second);
  }
};

class rg_demand {
 public:
  double flow;
  rgVertex source;
  rgVertex sink;
  double demand;
  std::map<const std::vector<rgVertex>,double> path_flow_map;
};



int updateExplicitFlow(rg_demand& demand,
					   std::vector<rgVertex>& penult,
					   double flow);
double calcD(Graph_rg& gdual);

