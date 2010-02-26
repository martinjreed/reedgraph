#include "RBGL.hpp"
#include <boost/graph/dijkstra_shortest_paths.hpp>
#include <math.h>
#include <string>
#include <R.h>
#include <R_ext/Parse.h>
#include <Rinternals.h>

template <class T>
inline std::string to_string (const T& t)
{
  std::stringstream ss;
  ss << t;
  return ss.str();
}

class Graph_rg
  : public boost::adjacency_list<boost::vecS, boost::vecS, boost::directedS,
								 boost::property<boost::vertex_color_t, 
												 boost::default_color_type>,
								 boost::property<boost::edge_weight_t, double,
												 boost::property<boost::edge_capacity_t, double> > >

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
	: Base(asInteger(num_verts_in))
    {
        if (!isReal(R_weights_in)) error("R_weights_in should be Real");
        if (!isInteger(R_edges_in)) error("R_edges_in should be integer");
        if (!isReal(capacities)) error("capacities should be Real");
        int NE = asInteger(num_edges_in);
        int* edges_in = INTEGER(R_edges_in);
        if (isReal(R_weights_in)) {
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
	  : Base(asInteger(num_verts_in))
    {
	  if (!isInteger(R_edges_in)) error("R_edges_in should be integer");
	  int NE = asInteger(num_edges_in);
	  int* edges_in = INTEGER(R_edges_in);
	  for (int i = 0; i < NE ; i++, edges_in += 2) {
		boost::add_edge(*edges_in, *(edges_in+1), 1, *this);
	  }
    }
};



extern "C" {

  void rg_cversion() {
	Rprintf("Version 1.1");
  }

  typedef boost::graph_traits < Graph_dd >::edge_descriptor Edge;
  typedef boost::graph_traits < Graph_dd >::vertex_descriptor Vertex;

  typedef boost::graph_traits < Graph_rg >::edge_descriptor rgEdge;
  typedef boost::graph_traits < Graph_rg >::vertex_descriptor rgVertex;

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


  double calcD(Graph_rg& gdual) {
	using namespace boost;
	double sum =0.0;
	graph_traits < Graph_rg >::edge_iterator ei, eend;
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


  class rg_demand {
	public:
	double flow;
	rgVertex source;
	rgVertex sink;
	double demand;
	std::map<const std::pair<rgVertex,rgVertex>,double> edge_flow_map;
	std::map<const std::vector<rgVertex>,double> path_flow_map;
  };

  std::pair<int,int> updateExplicitFlow(rg_demand& demand,
										std::vector<rgVertex>& penult,
										double flow) {
	rgVertex f,p;
	rgVertex source = demand.source;
	rgVertex sink = demand.sink;
	int newflows=0;
	int newpath=0;
	
	// note this records the path backwards! (slightly faster
	// when creating the path vector so no real reason to change
	std::vector<rgVertex> path;

	f = sink;
	p = penult[f];

	std::pair<rgVertex,rgVertex> key;
	key = std::pair<rgVertex,rgVertex>(p,f);
	path.push_back(f);
	path.push_back(p);
	// if key for this edge does not exist creat a new one
	if ( demand.edge_flow_map.find(key) == demand.edge_flow_map.end() ) {
	  demand.edge_flow_map[key] = flow;
	  newflows++;
	} else {
	  demand.edge_flow_map[key] += flow;
	}
	
	f = p;
	p = penult[p];

	while(f != source) {
	  key = std::pair<rgVertex,rgVertex>(p,f);
	  path.push_back(p);
	  if ( demand.edge_flow_map.find(key) == demand.edge_flow_map.end() ) {
		demand.edge_flow_map[key] = flow;
		newflows++;
	  } else {
		demand.edge_flow_map[key] += flow;
	  }
	  f = p;
	  p = penult[p];
	}
	
	if (demand.path_flow_map.find(path) == demand.path_flow_map.end()) {
	  demand.path_flow_map[path] = flow;
	  newpath=1;
	} else {
	  demand.path_flow_map[path] += flow;
	}

	return std::pair<int,int>(newflows,newpath);
  }
  double calcBeta(std::vector<rg_demand>& demands,
				  Graph_rg& graph) {
	using namespace boost;
	double alpha =0;
	double beta;
	int N = num_vertices(graph);
	std::vector<Vertex> p(N);
	std::vector<double> d(N);

	std::vector<rg_demand>::iterator vi,ve;
	for(vi=demands.begin(); vi < demands.end(); vi++) {
	  dijkstra_shortest_paths(graph, vi->source,
							  predecessor_map(&p[0]).distance_map(&d[0]));
	  alpha += vi->demand * d[vi->source];
	}
	beta = calcD(graph) /alpha;
	return beta;
  }
  SEXP rg_fleischer_max_concurrent_flow_c(SEXP num_verts_in,
										  SEXP num_edges_in,
										  SEXP R_edges_in,
										  SEXP R_weights_in,
										  SEXP capacities_in,
										  SEXP num_demands_in,
										  SEXP demands_sources_in,
										  SEXP demands_sinks_in,
										  SEXP demands_in,
										  SEXP Re,
										  SEXP Rupdateflow,
										  SEXP pb,
										  SEXP env,
										  SEXP Rprogress
										  ) {
	using namespace boost;

	bool progress=false;
	bool updateflow=asLogical(Rupdateflow);
	std::string title("");
	if ( isLogical(Rprogress) ) {
	  progress=asLogical(Rprogress);
	} else if ( isString(Rprogress) ) {
	  title = std::string(CHAR(STRING_ELT(Rprogress,0)));
	  progress=true;
	}
	Graph_rg gdual(num_verts_in, num_edges_in, R_edges_in, R_weights_in,capacities_in);
	
	double* dem_in = REAL(demands_in);
	int* dem_sources_in = INTEGER(demands_sources_in);
	int* dem_sinks_in = INTEGER(demands_sinks_in);
	int num_dem = asInteger(num_demands_in);
	std::vector<rg_demand> demands(num_dem);
	for(int i=0 ; i<num_dem; i++) {
	  demands[i].demand = dem_in[i];
	  demands[i].flow = 0;
	  demands[i].source = vertex(dem_sources_in[i],gdual);
	  demands[i].sink = vertex(dem_sinks_in[i],gdual);
	}
	int N = num_vertices(gdual);
	int m = num_edges(gdual);
	double e = asReal(Re);
	int num_dem_ed_flows =0;
	int num_dem_p_flows =0;
	double delta = pow(double(m) / (1.0 - e),-1.0/e);

	graph_traits < Graph_rg >::edge_iterator ei, eend;

	for(tie(ei,eend) = edges(gdual); ei != eend; ei++) {
	  int s = source(*ei,gdual);
	  int t = target(*ei,gdual);
	  std::pair<Edge, bool> e = edge(vertex(s,gdual),
									 vertex(t,gdual),gdual);
	  double c = get(edge_capacity,gdual,e.first);
	  put(edge_weight,gdual,e.first,delta/c);
	}
	

	double D;
	D=calcD(gdual);
	int doubreq =  int(ceil(2.0/e * log(m/(1-e))/log(1+e)));
	
	// assuming doubreq is about the maximum number of phases
	// then we want to only update the progress bar every 1%
	int updatepb = int(ceil(doubreq / 100.0));
	int phases =0;
	int totalphases =0;
  
	std::vector<Vertex> penult(N);
	std::vector<double> dist(N);
	std::vector<rg_demand>::iterator vi,ve;

	SEXP cmdsxp, cmdexpr, ansxp = R_NilValue;
	ParseStatus status;

	PROTECT(cmdsxp = allocVector(STRSXP,1));
	
	char cmd[256];
	//ansxp = eval(cmdexpr,env);

	double a = pow(2.0,-130);
	double b =1;
	double c = a +b;

	// phases
	while(D < 1.0) {
	  if(phases > doubreq) {
		//Rprintf("DOubling %d %d\n",doubreq,totalphases);
		for(vi=demands.begin(); vi < demands.end(); vi++) {
		  vi->demand = vi->demand * 2;
		}
		
		phases = 0;
	  }

	  // if progress only update about every 1%
	  if(progress && totalphases % updatepb ==0 ) {
		// calls R to update progress bar
		sprintf(cmd,"setTkProgressBar(pb,%d,label=\"%s\")",totalphases,title.c_str());
		SET_STRING_ELT(cmdsxp,0,mkChar(cmd));
		cmdexpr = PROTECT(R_ParseVector(cmdsxp, -1, &status, R_NilValue));
		for(int i = 0; i < length(cmdexpr); i++)
		  ansxp = eval(VECTOR_ELT(cmdexpr, i), env);
		//Rprintf("D=%lg\n",D);
		UNPROTECT(1);
	  }	  
		
	  // steps
	  for(vi=demands.begin(); vi < demands.end(); vi++) {
		rg_demand demand=*vi;
		rgVertex source = demand.source;
		rgVertex sink = demand.sink;
		rgVertex f,p;

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

		  f = p;
		  p = penult[p];
				 
		  while(f != source) {
			ed = edge(p,f,gdual);
			c =get(edge_capacity,gdual,ed.first);
			w =get(edge_weight,gdual,ed.first);
			w = w * (1 + (e * mincap) / c);
			put(edge_weight,gdual,ed.first,w);
			f = p;
			p = penult[p];
		  }
		  if(updateflow) {
			std::pair<int,int> upflowret = updateExplicitFlow(*vi,penult,mincap);
			num_dem_ed_flows += upflowret.first;
			num_dem_p_flows += upflowret.second;
		  }
		  vi->flow += mincap;

		  D=calcD(gdual);
		}
	  }
	  phases++;
	  totalphases++;
	}
	

	// need to return:
	// dual lengths
	// demand[].flow
	// demand[].[edgemap] key and values
    // demand[].[pathmap] key and values

	// set up return expressions
	SEXP duallenxp, flowsxp, lenxp, edmapkeyxp, edmapvalxp, pmapkeyxp, pmapvalxp, retlistxp;
	PROTECT(duallenxp = allocVector(REALSXP,m));
	PROTECT(flowsxp = allocVector(REALSXP,num_dem));
	PROTECT(lenxp = allocVector(INTSXP,3));
	PROTECT(edmapkeyxp = allocVector(INTSXP,num_dem_ed_flows * 3));
	PROTECT(edmapvalxp = allocVector(REALSXP,num_dem_ed_flows));
	PROTECT(pmapkeyxp = allocVector(VECSXP,num_dem_p_flows * 2));
	PROTECT(pmapvalxp = allocVector(REALSXP,num_dem_p_flows));
	PROTECT(retlistxp = allocVector(VECSXP,7));

	INTEGER(lenxp)[0] = num_dem_ed_flows;
	INTEGER(lenxp)[1] = num_dem_p_flows;
	INTEGER(lenxp)[2] = totalphases;

	int* edges_in = INTEGER(R_edges_in);
	for (int i = 0;
		 i < m ; 
		 i++, edges_in += 2) {
	  std::pair<rgEdge, bool> ed = edge(vertex(*edges_in,gdual),vertex(*(edges_in+1),gdual), gdual);
	  REAL(duallenxp)[i] = get(edge_weight,gdual,ed.first);
	}

	int i=0;
	int j=0;
	int k=0;
	int l=0;
	int n=0;
	for(vi=demands.begin(); vi < demands.end(); vi++, i++) {
	  REAL(flowsxp)[i]= vi->flow;
	  std::map<const std::pair<rgVertex,rgVertex>,double>::iterator mi,mend;
	  for(mi=vi->edge_flow_map.begin(); mi != vi->edge_flow_map.end() ; mi++) {
		INTEGER(edmapkeyxp)[j++]=i +1;
		INTEGER(edmapkeyxp)[j++]=mi->first.first +1;
		INTEGER(edmapkeyxp)[j++]=mi->first.second +1;
		REAL(edmapvalxp)[k++]=mi->second;
	  }
	  
	  std::map<const std::vector<rgVertex>, double>::iterator pmi, pmend;
	  for(pmi=vi->path_flow_map.begin(); pmi != vi->path_flow_map.end(); pmi++) {
		std::vector<rgVertex> path = pmi->first;
		double val = pmi->second;
		std::vector<rgVertex>::reverse_iterator vi, vend;
		std::string spath;
		std::string sdem = std::string(to_string<int>(i+1));
		vi=path.rbegin();
		spath.append(to_string<int>(*vi + 1));
		vi++;
		for(; vi != path.rend(); vi++) {
		  spath.append("|");
		  spath.append(to_string<int>(*vi + 1));
		}
		SEXP pathsxp, demsxp;
		// create the demand string identifier
		PROTECT(demsxp = allocVector(STRSXP,1));
		SET_STRING_ELT(demsxp,0,mkChar(sdem.c_str()));
		// create the path string
		PROTECT(pathsxp = allocVector(STRSXP,1));
		SET_STRING_ELT(pathsxp,0,mkChar(spath.c_str()));
		SET_VECTOR_ELT(pmapkeyxp,l++, demsxp);
		SET_VECTOR_ELT(pmapkeyxp,l++, pathsxp);
		REAL(pmapvalxp)[n++] = val;
		UNPROTECT(2);
	  }
		
	}
	
	SET_VECTOR_ELT(retlistxp,0,duallenxp);
	SET_VECTOR_ELT(retlistxp,1,flowsxp);
	SET_VECTOR_ELT(retlistxp,2,lenxp);
	SET_VECTOR_ELT(retlistxp,3,edmapkeyxp);
	SET_VECTOR_ELT(retlistxp,4,edmapvalxp);
	SET_VECTOR_ELT(retlistxp,5,pmapkeyxp);
	SET_VECTOR_ELT(retlistxp,6,pmapvalxp);
	UNPROTECT(9);
	return(retlistxp);
  
  }
}
