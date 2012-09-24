#include <vector>
#include <list>
#include <utility>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/dijkstra_shortest_paths.hpp>

class mf_demand;

typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::directedS,
			      boost::property<boost::vertex_color_t, 
					      boost::default_color_type>,
			      boost::property<boost::edge_weight_t, double,
					      boost::property<boost::edge_capacity_t, 
							      double > > > NetGraph;

class Graph_mf: public NetGraph {

public:
  typedef NetGraph::graph_property_type graph_property_type;
  typedef NetGraph::vertices_size_type vertices_size_type;
  typedef NetGraph::edges_size_type edges_size_type;
  
  long number_flows;
  inline Graph_mf()
    : NetGraph() { }
  inline Graph_mf(const graph_property_type& p)
  : NetGraph(p) { }
  inline Graph_mf(const NetGraph& x)
  : NetGraph(x) { }
  inline Graph_mf(vertices_size_type num_vertices)
  : NetGraph(num_vertices) { }
  inline Graph_mf(vertices_size_type num_vertices,
		  const graph_property_type& p)
    : NetGraph(num_vertices, p) { }

  inline Graph_mf(int num_vertices,
		  const std::vector< int > &edges,
		  const std::vector<double> &capacities,
		  const std::vector<int> demand_sources,
		  const std::vector<int> demand_sinks)
    : NetGraph(num_vertices)
  {
    int NE = edges.size()/2;
    for (int i = 0, j=0; i < NE ;  i++, j+=2 ) {
      std::pair<edge_descriptor, bool> e =
	boost::add_edge(edges[j], edges[j+1],
			0.0, *this);
      boost::put(boost::edge_capacity,*this,e.first,capacities[i]);
    }
  } 
  void max_concurrent_flow(std::vector<mf_demand> &demands,
		      double e=0.1);
  NetGraph gdual;
  long totalphases;
private:
  double calcD();
};

typedef boost::graph_traits < NetGraph >::edge_descriptor Edge;
typedef boost::graph_traits < NetGraph >::vertex_descriptor Vertex;

class mf_demand {
 public:
  double flow;
  Vertex source;
  Vertex sink;
  double demand;
  std::map<const std::list<Vertex>,double> path_flow_map;
};
