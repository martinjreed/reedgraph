// mjreed

#include "rgalgorithms.h"

extern "C" {

SEXP rg_test_c(SEXP v) {
  //  using namespace Rcpp;
  //using namespace std;

  RcppResultSet rs;
  rs.add("sum",1);
  return rs.getReturnList();
}

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
	std::map<const std::vector<rgVertex>,double> path_flow_map;
  };

  int updateExplicitFlow(rg_demand& demand,
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
	
	f = p;
	p = penult[p];

	while(f != source) {
	  path.push_back(p);
	  f = p;
	  p = penult[p];
	}
	
	if (demand.path_flow_map.find(path) == demand.path_flow_map.end()) {
	  demand.path_flow_map[path] = flow;
	  newpath=1;
	} else {
	  demand.path_flow_map[path] += flow;
	}

	return newpath;
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
  SEXP rg_fleischer_max_concurrent_flow_c_keep(SEXP num_verts_in,
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

	std::vector<int> demand_index(num_dem);
	for(int i=0; i<num_dem; i++) {
	  demand_index[i]=i;
	}

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
		sprintf(cmd,"setTxtProgressBar(pb,%d)",totalphases);
		SET_STRING_ELT(cmdsxp,0,mkChar(cmd));
		cmdexpr = PROTECT(R_ParseVector(cmdsxp, -1, &status, R_NilValue));
		for(int i = 0; i < length(cmdexpr); i++)
		  ansxp = eval(VECTOR_ELT(cmdexpr, i), env);
		//Rprintf("D=%lg\n",D);
		UNPROTECT(1);
	  }	  
		
	  // steps
	  random_shuffle(demand_index.begin(),demand_index.end());
	  
	  for(int j=0; j<num_dem;j++) {
		int i= demand_index[j];
		rg_demand demand=demands[i];
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
			num_dem_p_flows += updateExplicitFlow(demands[i],penult,mincap);
		  }
		  demands[i].flow += mincap;

		  D=calcD(gdual);
		}
	  }
	  phases++;
	  totalphases++;
	}
	

	// need to return:
	// dual lengths
	// demand[].flow
    // demand[].[pathmap] key and values

	// set up return expressions
	SEXP duallenxp, flowsxp, lenxp, pmapkeyxp, pmapvalxp, retlistxp;
	PROTECT(duallenxp = allocVector(REALSXP,m));
	PROTECT(flowsxp = allocVector(REALSXP,num_dem));
	PROTECT(lenxp = allocVector(INTSXP,3));
	PROTECT(pmapkeyxp = allocVector(VECSXP,num_dem_p_flows * 2));
	PROTECT(pmapvalxp = allocVector(REALSXP,num_dem_p_flows));
	PROTECT(retlistxp = allocVector(VECSXP,5));

	INTEGER(lenxp)[0] = num_dem_p_flows;
	INTEGER(lenxp)[1] = totalphases;

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
	SET_VECTOR_ELT(retlistxp,3,pmapkeyxp);
	SET_VECTOR_ELT(retlistxp,4,pmapvalxp);
	UNPROTECT(7);
	return(retlistxp);
  
  }

  std::pair<int,double> findshortestpathcost(std::vector< std::vector<int> >& paths,
					   std::vector<double>& vlengths) {
	int numpaths = paths.size();
	std::vector<int>::iterator vi;
	double sum=0.0;
	int minp=-1;
	double min=DBL_MAX;
	
	std::vector<long> trypath(numpaths);
	for(int i=0;i<numpaths;i++) {
	  trypath[i]=i;
	}
	random_shuffle(trypath.begin(),trypath.end());

	std::vector<long>::iterator tit,tend;
	
	tend=trypath.end();
	for(tit=trypath.begin();tit != tend;tit++) {
	  sum=0;
	  for(vi=paths[*tit].begin(); vi != paths[*tit].end(); vi++) {
		sum+=vlengths[*vi];
	  }
	  if( sum < min) {
		minp=*tit;
		min = sum;
	  }
	}
	return (std::pair<int,double>(minp,min));
  }

  std::pair<int,double> 
  findshortestpathcostopt(
						  std::vector< std::vector<int> >& paths,
						  std::vector<double>& vlengths,
						  std::vector<double>& vcapacity,
						  std::vector<double>& vweights,
						  double demand) {
	int numpaths = paths.size();
	std::vector<int>::iterator vi;
	double sum=0.0;
	int minp=-1;
	double min=DBL_MAX;
	
	std::vector<long> trypath(numpaths);
	for(int i=0;i<numpaths;i++) {
	  trypath[i]=i;
	}
	//random_shuffle(trypath.begin(),trypath.end());

	std::vector< std::pair<int,double> > best;
	std::vector< std::pair<int,double> >::iterator bit,bend;

	best.push_back(std::pair<int,double>(-1,DBL_MAX));

	std::vector<long>::iterator tit,tend;
	
	tend=trypath.end();
	for(tit=trypath.begin();tit != tend;tit++) {
	  sum=0;
	  for(vi=paths[*tit].begin(); vi != paths[*tit].end(); vi++) {
		sum+=vlengths[*vi];
	  }
	  if( sum < min) {
		minp=*tit;
		min = sum;
	  }
	  //	  if( sum < best.front().second * 0.9999999999999 ) {
	  if( sum < best.front().second) {
		best.clear();
		best.push_back(std::pair<int,double>(*tit,sum));
		//	  } else if(sum < best.front().second * 1.0000000000001) {
	  } else if(sum == best.front().second) {
		best.push_back(std::pair<int,double>(*tit,sum));
	  }
	}
	return (std::pair<int,double>(minp,min));
	
	bend=best.end();
	double bestcap=-DBL_MAX;
	for(bit=best.begin(); bit != bend; bit++){
	  double mincap = DBL_MAX;
	  for(vi=paths[(*bit).first].begin(); vi != paths[(*bit).first].end(); vi++) {
		if(vcapacity[*vi]-vweights[*vi]-demand < mincap)
		  mincap = vcapacity[*vi]-vweights[*vi]-demand;
	  }
	  if(mincap > bestcap){
		bestcap = mincap;
		minp=(*bit).first;
		min=(*bit).second;
	  }

	}
	return (std::pair<int,double>(minp,min));
  }


  int findshortestpath(std::vector< std::vector<int> >& paths,
					   std::vector<double>& vlengths) {
	int numpaths = paths.size();
	std::vector<int>::iterator vi;
	int minp=0;
	double min=DBL_MAX;

	std::vector<long> trypath(numpaths);
	for(int i=0;i<numpaths;i++) {
	  trypath[i]=i;
	}
	//random_shuffle(trypath.begin(),trypath.end());

	for(int j=0;j<numpaths;j++) {
	  int i=trypath[j];
	  double sum=0;
	  for(vi=paths[i].begin(); vi != paths[i].end(); vi++) {
		sum+=vlengths[*vi];
	  }
	  if( sum < min) {
		minp=i;
		min = sum;
	  }
	}
	return minp;
  }

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




  SEXP rg_max_concurrent_flow_capacity_restricted_c
  (
   SEXP RRdemandpaths,
   SEXP Rdemands,
   SEXP Rcapacity,
   SEXP Re,
   SEXP Rprogress,
   SEXP pb,
   SEXP env
   ) {

	bool progress=Rcpp::as<bool>(Rprogress);
	double e = Rcpp::as<double>(Re);
	std::vector<double> vcapacity = RcppVector<double>(Rcapacity).stlVector();
	std::vector<double> vdemands = RcppVector<double>(Rdemands).stlVector();
	int numD = vdemands.size();
	std::vector<double>::iterator vecdit;
	int M = vcapacity.size();
	//decode Rdemandpaths into paths
	std::vector<int>::iterator vi,ve;
	std::vector< std::vector< std::vector<int> > > paths;
	std::vector< std::vector< std::vector<int> > >::iterator dit;
	std::vector< std::vector<int> >::iterator pit;
	std::vector<int>::iterator eit;

	Rcpp::List demands(RRdemandpaths);
	int sz = demands.size();
	paths.resize(sz);
	for(int i=0;i<sz;i++) {
	  Rcpp::List demandpaths((SEXPREC*)demands[i]);
	  int sz = demandpaths.size();
	  paths[i].resize(sz);
	  for(int j=0;j<sz;j++) {
		paths[i][j]= Rcpp::as< std::vector<int> >(demandpaths[j]);
	  }
	}

	std::vector< std::vector<double> > pathflows(paths.size());
	for(int i=0;i<numD;i++) {
	  pathflows[i] = std::vector<double>(paths[i].size(),0.0);
	}

	int doubleCount=0;
	int phases=0;
	int totalphases=0;

	double delta = pow(double(M) / (1.0 - e),-1.0/e);
	std::vector<double> vlengths(M);
	for(int i=0;i<M;i++) {
	  vlengths[i]=delta / vcapacity[i];
	}

	int doubreq =  int(ceil(2.0/e * log(M/(1-e))/log(1+e)));
	int updatepb = int(ceil(doubreq / 100.0));

	double D=0.0;
	for(int i=0;i<M;i++) {
	  D+=vcapacity[i]*vlengths[i];
	}
	char cmd[256];
	SEXP cmdsxp, cmdexpr, ansxp = R_NilValue;
	ParseStatus status;

	PROTECT(cmdsxp = allocVector(STRSXP,1));
	int i=0;
	while( D < 1.0) {
	  if(phases > doubreq) {
		for(int i=0;i<numD;i++) {
		  vdemands[i] = vdemands[i] *2;
		}
		
		phases = 0;
	  }

	  std::vector<double> vcaptmp=vcapacity;

	  if(progress != false && totalphases % updatepb == 0) {
		sprintf(cmd,"setTxtProgressBar(pb,%d)",totalphases);
		SET_STRING_ELT(cmdsxp,0,mkChar(cmd));
		cmdexpr = PROTECT(R_ParseVector(cmdsxp, -1, &status, R_NilValue));
		for(int i = 0; i < length(cmdexpr); i++)
		  ansxp = eval(VECTOR_ELT(cmdexpr, i), env);
		//Rprintf("D=%lg\n",D);
		UNPROTECT(1);
	  }
	  
	  // note fix for rotation of demands
	  for(int j=0; j<numD;j++,i++) {
		i = i >= numD ? i - numD: i;
		double demand = vdemands[i];
		while( D < 1.0 && demand > 0.0) {
		  int p = findshortestpath(paths[i],vlengths);
		  //find minium of demand or capacity on path
		  double mincap=demand;
		  ve=paths[i][p].end();
		  for(vi=paths[i][p].begin(); vi != ve; vi++) {
			if(mincap > vcaptmp[*vi])
			  mincap = vcaptmp[*vi];
		  }

		  if(mincap <=0)
			break;
		  demand -= mincap;
		  int sz = paths[i][p].size();
		  for(int k=0; k < sz; k++) {
			double length = vlengths[paths[i][p][k]];
			length = length * (1.0 + (e * mincap) /
							   vcapacity[paths[i][p][k]]);
			vlengths[paths[i][p][k]] = length;
		  }
		  pathflows[i][p] += mincap;
		  D=0.0;
		  for(vi=paths[i][p].begin(); vi != ve; vi++) {
			vcaptmp[*vi] -= mincap;
		  }

		  for(int k=0;k<M;k++) {
			D+=vcapacity[k]*vlengths[k];
		  }
		}
	  }
	  i++;
	  phases++;
	  totalphases++;
	  
	}

	UNPROTECT(1);
	
	Rcpp::List retlist;
	retlist.push_back(totalphases,"totalphases");
	retlist.push_back(pathflows,"pathflows");
	retlist.push_back(vlengths,"vlengths");
	return Rcpp::wrap(retlist);
  }

  SEXP rg_fleischer_max_concurrent_flow_restricted_c
  (
   SEXP RRdemandpaths,
   SEXP Rdemands,
   SEXP Rcapacity,
   SEXP Re,
   SEXP Rprogress,
   SEXP pb,
   SEXP Rintgamma,
   SEXP env
   ) {

	bool progress=Rcpp::as<bool>(Rprogress);
	double e = Rcpp::as<double>(Re);
	double intgamma = Rcpp::as<double>(Rintgamma);
	bool calcgamma=false;
	if(intgamma >=0) 
	  calcgamma=true;
	std::vector<double> vcapacity = RcppVector<double>(Rcapacity).stlVector();
	std::vector<double> vdemands = RcppVector<double>(Rdemands).stlVector();
	std::vector<double> origdemands(vdemands);	

	int numD = vdemands.size();
	std::vector<double>::iterator vecdit;
	int M = vcapacity.size();
	//decode Rdemandpaths into paths
	std::vector<int>::iterator vi,ve;
	std::vector< std::vector< std::vector<int> > > paths;
	std::vector< std::vector< std::vector<int> > >::iterator dit;
	std::vector< std::vector<int> >::iterator pit;
	std::vector<int>::iterator eit;

	Rcpp::List demands(RRdemandpaths);
	int sz = demands.size();
	paths.resize(sz);
	for(int i=0;i<sz;i++) {
	  Rcpp::List demandpaths((SEXPREC*)demands[i]);
	  int sz = demandpaths.size();
	  paths[i].resize(sz);
	  for(int j=0;j<sz;j++) {
		paths[i][j]= Rcpp::as< std::vector<int> >(demandpaths[j]);
	  }
	}

	std::vector< std::vector<double> > pathflows(paths.size());
	for(int i=0;i<numD;i++) {
	  pathflows[i] = std::vector<double>(paths[i].size(),0.0);
	}

	int doubleCount=0;
	int phases=0;
	int totalphases=0;

	double delta = pow(double(M) / (1.0 - e),-1.0/e);
	std::vector<double> vlengths(M);
	for(int i=0;i<M;i++) {
	  vlengths[i]=delta / vcapacity[i];
	}

	int doubreq =  int(ceil(2.0/e * log(M/(1-e))/log(1+e)));
	int updatepb = int(ceil(doubreq / 100.0));

	double D=0.0;
	for(int i=0;i<M;i++) {
	  D+=vcapacity[i]*vlengths[i];
	}
	char cmd[256];
	SEXP cmdsxp, cmdexpr, ansxp = R_NilValue;
	ParseStatus status;

	PROTECT(cmdsxp = allocVector(STRSXP,1));
	int i=0;
	int countgamma=0;
	double bestgamma= -DBL_MAX;
	std::vector<int> demand_index(numD);
	for(int i=0; i<numD; i++) {
	  demand_index[i]=i;
	}
	std::vector<int> bestpaths;


	while( D < 1.0) {
	  if(phases > doubreq) {
		Rprintf("doubling!!\n");
		for(int i=0;i<numD;i++) {
		  vdemands[i] = vdemands[i] *2;
		}
		
		phases = 0;
	  }


	  if(progress != false && totalphases % updatepb == 0) {
		sprintf(cmd,"setTxtProgressBar(pb,%d)",totalphases);
		SET_STRING_ELT(cmdsxp,0,mkChar(cmd));
		cmdexpr = PROTECT(R_ParseVector(cmdsxp, -1, &status, R_NilValue));
		for(int i = 0; i < length(cmdexpr); i++)
		  ansxp = eval(VECTOR_ELT(cmdexpr, i), env);
		//Rprintf("D=%lg\n",D);
		UNPROTECT(1);
	  }
	  
	  std::vector<double> weights(M,0.0);
	  bool underdemand=true;

	  // note fix for random rotation of demand order
	  random_shuffle(demand_index.begin(),demand_index.end());

	  std::vector<int> pathrecord(numD);
	  for(int j=0; j<numD;j++,i++) {
		i= demand_index[j];
		double demand = vdemands[i];
		if(D >= 1.0)
		  underdemand=false;
		while( D < 1.0 && demand > 0.0) {
		  int p = findshortestpath(paths[i],vlengths);
		  pathrecord[i]=p;
		  //find minium of demand or capacity on path
		  double mincap=demand;
		  double minfree=DBL_MAX;
		  double smallestcap=0;
		  ve=paths[i][p].end();
		  for(vi=paths[i][p].begin(); vi != ve; vi++) {
			if(mincap > vcapacity[*vi])
			  mincap = vcapacity[*vi];
			if(minfree > vcapacity[*vi] - demand) {
			  minfree = vcapacity[*vi] - demand;
			  smallestcap = vcapacity[*vi];
			}
			
			if(calcgamma)
			  weights[*vi] += origdemands[i];
		  }
		  if(mincap < demand) {
			underdemand=false;
			//Rprintf("underdemand\n");
		  }

		  demand -= mincap;
		  int sz = paths[i][p].size();
		  for(int k=0; k < sz; k++) {
			double length = vlengths[paths[i][p][k]];
			length = length * (1.0 + (e * mincap) /
							   vcapacity[paths[i][p][k]]);
			vlengths[paths[i][p][k]] = length;
		  }
		  pathflows[i][p] += mincap;
		  D=0.0;
		  for(int k=0;k<M;k++) {
			D+=vcapacity[k]*vlengths[k];
		  }
		}
	  }
	  phases++;
	  totalphases++;
	  if(calcgamma && underdemand) {
		double gamma = DBL_MAX;
		for(int j=0; j< M; j++) {
		  double tmp = 1.0 - weights[j]/vcapacity[j];
		  if( gamma > tmp ) {
			gamma = tmp;
		  }
		}
		if(gamma >= intgamma * 0.9999999999) {
		  countgamma++;
		  //Rprintf("gamma %lg, matches int gamma of %lg\n",gamma, intgamma);
		}
		if(gamma > bestgamma) {
		  bestgamma = gamma;
		  bestpaths = pathrecord;
		  //Rprintf("bestgamma so far %lg\n",bestgamma);
		  //Rprintf("best paths: ");
		  //for(int j=0; j<numD; j++)
		  //Rprintf("%d,",bestpaths[j]);
		  //Rprintf("\n");
		}
	  }
	  
	}

	UNPROTECT(1);
	
	Rcpp::List retlist;
	retlist.push_back(totalphases,"totalphases");
	retlist.push_back(countgamma,"countgamma");
	retlist.push_back(bestgamma,"bestgamma");
	retlist.push_back(bestpaths,"bestpaths");
	retlist.push_back(pathflows,"pathflows");
	retlist.push_back(vlengths,"vlengths");
	return Rcpp::wrap(retlist);
  }



  SEXP rg_fleischer_max_concurrent_flow_c_old(SEXP num_verts_in,
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
		sprintf(cmd,"setTxtProgressBar(pb,%d)",totalphases);
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
			num_dem_p_flows += updateExplicitFlow(*vi,penult,mincap);
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
    // demand[].[pathmap] key and values

	// set up return expressions
	SEXP duallenxp, flowsxp, lenxp, pmapkeyxp, pmapvalxp, retlistxp;
	PROTECT(duallenxp = allocVector(REALSXP,m));
	PROTECT(flowsxp = allocVector(REALSXP,num_dem));
	PROTECT(lenxp = allocVector(INTSXP,3));
	PROTECT(pmapkeyxp = allocVector(VECSXP,num_dem_p_flows * 2));
	PROTECT(pmapvalxp = allocVector(REALSXP,num_dem_p_flows));
	PROTECT(retlistxp = allocVector(VECSXP,5));

	INTEGER(lenxp)[0] = num_dem_p_flows;
	INTEGER(lenxp)[1] = totalphases;

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
	SET_VECTOR_ELT(retlistxp,3,pmapkeyxp);
	SET_VECTOR_ELT(retlistxp,4,pmapvalxp);
	UNPROTECT(7);
	return(retlistxp);
  
  }


  double cost_on_path(Graph_rg gdual,int path[],int len,bool print,double rel) {

	using namespace boost;
	std::pair<Edge, bool> ed = edge(path[0]-1,path[1]-1,gdual);
	double cost =get(edge_weight,gdual,ed.first);
	Rprintf("|%ld|%ld|",path[0],path[1]);
	
	for(int i=2;i<len; i++) {
	  ed = edge(path[i-1]-1,path[i]-1,gdual);
	  cost +=get(edge_weight,gdual,ed.first);
	  Rprintf("%ld|",path[i]);
	}
	Rprintf(" = %lg\n",cost/rel);
	return cost;
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
										  SEXP Rprogress,
										  SEXP Rpermutation
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

	std::vector<int> permutation = RcppVector<int>(Rpermutation).stlVector();
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

	// to be removed
	std::vector< std::vector< std::vector<int> > > paths;

	std::vector< std::vector< int > > pathCount;
	int sz = demands.size();
	paths.resize(sz);
	pathCount.resize(sz);
	//end


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

	std::vector<int> demand_index(num_dem);
	if(permutation.size()==num_dem) {
	  demand_index = permutation;
	} else {
	  for(int i=0; i<num_dem; i++) {
		demand_index[i]=i;
	  }
	}

	//random_shuffle(demand_index.begin(),demand_index.end());
	if(delta == 0) {
	  Rprintf("Error delta=0\n");
	  // fix this just to "bomb out"
	  D=1.0;
	}

	std::vector< std::pair<int, double> > costpair(num_dem);

	if(true) {
	  
	  for(int x=0;x<num_dem;x++) {
		rg_demand demand=demands[x];
		rgVertex source = demand.source;
		rgVertex sink = demand.sink;
		
		dijkstra_shortest_paths(gdual, source,
								predecessor_map(&penult[0]).distance_map(&dist[0]));
		costpair[x]=std::pair<int,double>(x,dist[sink]);
		//Rprintf("distance for %ld is %lg\n",x,dist[sink]);
	  }
	  std::sort(costpair.begin(),costpair.end(),less_cost());
	}
	// phases
	while(D < 1.0) {
	  //Rprintf("\n new phase\n Doing demand:");
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
		sprintf(cmd,"setTxtProgressBar(pb,%d)",totalphases);
		SET_STRING_ELT(cmdsxp,0,mkChar(cmd));
		cmdexpr = PROTECT(R_ParseVector(cmdsxp, -1, &status, R_NilValue));
		for(int i = 0; i < length(cmdexpr); i++)
		  ansxp = eval(VECTOR_ELT(cmdexpr, i), env);
		//Rprintf("D=%lg\n",D);
		UNPROTECT(1);
	  }	  
	  if(true) {
		for(int x=0;x<num_dem;x++) {
		  rg_demand demand=demands[x];
		  rgVertex source = demand.source;
		  rgVertex sink = demand.sink;
		  
		  dijkstra_shortest_paths(gdual, source,
								  predecessor_map(&penult[0]).distance_map(&dist[0]));
		  costpair[x]=std::pair<int,double>(x,dist[sink]);
		  //Rprintf("distance for %ld is %lg\n",x,dist[sink]);
		}
		std::sort(costpair.begin(),costpair.end(),less_cost());
	  }
	  // steps
	  if(permutation.size()!=num_dem) {
		//random_shuffle(demand_index.begin()+2,demand_index.end()-6);
		//random_shuffle(demand_index.begin(),demand_index.end());
	  }
	  //else
		//random_shuffle(demand_index.begin(),demand_index.begin()+2);
		//random_shuffle(demand_index.begin(),demand_index.end());
	  
	  //	  for(int j=num_dem-1; j>=0;j--) {
	  for(int j=0; j<num_dem;j++) {
		//int i= demand_index[j];
		int i = costpair[j].first;
		//Rprintf("%ld, ",i);
		rg_demand demand=demands[i];
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

		  //remove
		  std::vector<int> tmppath(0);
		  //Rprintf("%ld,",(long)f);
		  tmppath.push_back(f);
		  //end
		  
		  f = p;
		  p = penult[p];
		  double w =get(edge_weight,gdual,ed.first);
		  
		  while(f != source) {
			//remove
			//Rprintf("%ld,",(long)f);
			tmppath.push_back(f);
			//end
			ed = edge(p,f,gdual);
			double cap =get(edge_capacity,gdual,ed.first);
			double w =get(edge_weight,gdual,ed.first);
			mincap = cap < mincap? cap : mincap;
			f = p;
			p = penult[p];
		  }
		  //remove
		  //Rprintf("%ld\n",(long)f);
		  tmppath.push_back(f);
		  //end
		  
		  if(false) {
			Rprintf("demand %ld\n",i);
			
			std::vector<int> mypath(tmppath.size());
			for(int x=0,y=tmppath.size()-1;x<tmppath.size();x++,y--) {
			  mypath[x]=tmppath[y]+1;
			}
			double lowest=cost_on_path(gdual,&mypath[0],tmppath.size(),true,1.0);
			Rprintf("%lg\n",lowest);
			int path1[3]= {1,3,2};
			cost_on_path(gdual,path1,3,true,lowest);
			
			int path2[3]= {1,9,2};
			cost_on_path(gdual,path2,3,true,lowest);
			
			int path3[4]= {1,4,3,2};
			cost_on_path(gdual,path3,4,true,lowest);
			
			int path4[5]= {1,7,6,9,2};
			cost_on_path(gdual,path4,5,true,lowest);
		  }


		  bool samepath=false;
		  for(int k=0;k<paths[i].size();k++) {
			samepath=true;
			for(int pi=0;pi<paths[i][k].size() && pi<tmppath.size();pi++){
			  if(paths[i][k][pi] != tmppath[pi])
				samepath=false;
			}
			if(samepath) {
			  pathCount[i][k]++;
			  break;
			}
		  }
		  if(!samepath) {
			pathCount[i].push_back(1);
			paths[i].push_back(tmppath);
		  }
		  //end
		  // now we have the maximum flow we can push through this
		  // step, and update demand (will add flow later)
		  mincap = demand.demand < mincap ? demand.demand : mincap;
		  demand.demand = demand.demand - mincap;
		  
		  /*if(demand.demand > 0) {
			Rprintf("demand over capacity\n");
			}*/
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
			num_dem_p_flows += updateExplicitFlow(demands[i],penult,mincap);
		  }
		  demands[i].flow += mincap;

		  D=calcD(gdual);
		}
	  }
	  phases++;
	  totalphases++;
	}

	/*
	for(int i=0;i<paths.size();i++) {
	  Rprintf("Demand %ld\n",i+1);
	  for(int j=0;j<paths[i].size();j++) {
		for(int k=paths[i][j].size()-1;k>=0;k--) {
		  Rprintf("%ld,",paths[i][j][k]);
		}
		Rprintf(": %ld\n",pathCount[i][j]);
	  }
	  }*/

	// need to return:
	// dual lengths
	// demand[].flow
    // demand[].[pathmap] key and values

	// set up return expressions
	SEXP duallenxp, flowsxp, lenxp, pmapkeyxp, pmapvalxp, retlistxp;
	PROTECT(duallenxp = allocVector(REALSXP,m));
	PROTECT(flowsxp = allocVector(REALSXP,num_dem));
	PROTECT(lenxp = allocVector(INTSXP,3));
	PROTECT(pmapkeyxp = allocVector(VECSXP,num_dem_p_flows * 2));
	PROTECT(pmapvalxp = allocVector(REALSXP,num_dem_p_flows));
	PROTECT(retlistxp = allocVector(VECSXP,5));

	INTEGER(lenxp)[0] = num_dem_p_flows;
	INTEGER(lenxp)[1] = totalphases;

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
	SET_VECTOR_ELT(retlistxp,3,pmapkeyxp);
	SET_VECTOR_ELT(retlistxp,4,pmapvalxp);
	UNPROTECT(7);
	return(retlistxp);
  
  }

  SEXP rg_fleischer_max_concurrent_flow_restricted_c_test
  (
   SEXP RRdemandpaths,
   SEXP Rdemands,
   SEXP Rcapacity,
   SEXP Re,
   SEXP Rprogress,
   SEXP pb,
   SEXP Rintgamma,
   SEXP Rbestpaths,
   SEXP env
   ) {

	bool progress=Rcpp::as<bool>(Rprogress);
	double e = Rcpp::as<double>(Re);
	double intgamma = Rcpp::as<double>(Rintgamma);
	bool calcgamma=false;
	if(intgamma >=0) 
	  calcgamma=true;
	std::vector<double> vcapacity = RcppVector<double>(Rcapacity).stlVector();
	std::vector<double> vdemands = RcppVector<double>(Rdemands).stlVector();
	std::vector<double> origdemands(vdemands);	

	int numD = vdemands.size();
	std::vector<double>::iterator vecdit;
	int M = vcapacity.size();
	//decode Rdemandpaths into paths
	std::vector<int>::iterator vi,ve;
	std::vector< std::vector< std::vector<int> > > paths;
	std::vector< std::vector< std::vector<int> > >::iterator dit;
	std::vector< std::vector<int> >::iterator pit;
	std::vector<int>::iterator eit;

	std::vector< std::vector< int > > pathCount;

	Rcpp::List demands(RRdemandpaths);
	int sz = demands.size();
	paths.resize(sz);
	pathCount.resize(sz);
	for(int i=0;i<sz;i++) {
	  Rcpp::List demandpaths((SEXPREC*)demands[i]);
	  int sz = demandpaths.size();
	  paths[i].resize(sz);
	  pathCount[i].resize(sz);
	  for(int j=0;j<sz;j++) {
		paths[i][j]= Rcpp::as< std::vector<int> >(demandpaths[j]);
		pathCount[i][j]=0;
	  }
	}

	std::vector< std::vector<double> > pathflows(paths.size());
	for(int i=0;i<numD;i++) {
	  pathflows[i] = std::vector<double>(paths[i].size(),0.0);
	}

	int doubleCount=0;
	int phases=0;
	int totalphases=0;

	double delta = pow(double(M) / (1.0 - e),-1.0/e);
	std::vector<double> vlengths(M);
	for(int i=0;i<M;i++) {
	  vlengths[i]=delta / vcapacity[i];
	}

	int doubreq =  int(ceil(2.0/e * log(M/(1-e))/log(1+e)));
	int updatepb = int(ceil(doubreq / 100.0));

	double D=0.0;
	for(int i=0;i<M;i++) {
	  D+=vcapacity[i]*vlengths[i];
	}
	char cmd[256];
	SEXP cmdsxp, cmdexpr, ansxp = R_NilValue;
	ParseStatus status;

	PROTECT(cmdsxp = allocVector(STRSXP,1));
	int i=0;
	int countgamma=0;
	double bestgamma= -DBL_MAX;
	std::vector<int> demand_index(numD);
	for(int i=0; i<numD; i++) {
	  demand_index[i]=i;
	}
	std::vector<int> bestpaths;


	while( D < 1.0) {
	  if(phases > doubreq) {
		Rprintf("doubling!!\n");
		for(int i=0;i<numD;i++) {
		  vdemands[i] = vdemands[i] *2;
		}
		
		phases = 0;
	  }


	  if(progress != false && totalphases % updatepb == 0) {
		sprintf(cmd,"setTxtProgressBar(pb,%d)",totalphases);
		SET_STRING_ELT(cmdsxp,0,mkChar(cmd));
		cmdexpr = PROTECT(R_ParseVector(cmdsxp, -1, &status, R_NilValue));
		for(int i = 0; i < length(cmdexpr); i++)
		  ansxp = eval(VECTOR_ELT(cmdexpr, i), env);
		Rprintf("D=%lg\n",D);
		UNPROTECT(1);
	  }
	  
	  std::vector<double> gweights(M,0.0);
	  std::vector<double> weights(M,0.0);
	  bool underdemand=true;

	  // note fix for random rotation of demand order
	  random_shuffle(demand_index.begin(),demand_index.end());

	  std::vector<int> pathrecord(numD);
	  for(int j=0; j<numD;j++,i++) {
		i= demand_index[j];
		double demand = vdemands[i];
		if(D >= 1.0)
		  underdemand=false;
		while( D < 1.0 && demand > 0.0) {
		  int p = findshortestpath(paths[i],vlengths);
		  /*std::pair<long, double> ppair = 
			findshortestpathcostopt(paths[i],
									vlengths,
									vcapacity,
									weights,
									vdemands[i]);
									int p = ppair.first;*/
		  if(underdemand)
			pathrecord[i]=p;
		  pathCount[i][p]++;
		  //find minium of demand or capacity on path
		  double mincap=demand;
		  double minfree=DBL_MAX;
		  double smallestcap=0;
		  ve=paths[i][p].end();
		  for(vi=paths[i][p].begin(); vi != ve; vi++) {
			if(mincap > vcapacity[*vi])
			  mincap = vcapacity[*vi];
			if(minfree > vcapacity[*vi] - demand) {
			  minfree = vcapacity[*vi] - demand;
			  smallestcap = vcapacity[*vi];
			}
			
			/*if(calcgamma)
			  weights[*vi] += origdemands[i];*/
		  }
		  if(mincap < demand) {
			underdemand=false;
			//Rprintf("underdemand\n");
		  }
		  for(vi=paths[i][p].begin(); vi != ve; vi++) {
			weights[*vi] += mincap;
		  }
		  demand -= mincap;
		  int sz = paths[i][p].size();
		  for(int k=0; k < sz; k++) {
			double length = vlengths[paths[i][p][k]];
			length = length * (1.0 + (e * mincap) /
							   vcapacity[paths[i][p][k]]);
			vlengths[paths[i][p][k]] = length;
		  }
		  pathflows[i][p] += mincap;
		  D=0.0;
		  for(int k=0;k<M;k++) {
			D+=vcapacity[k]*vlengths[k];
		  }
		}
	  }
	  phases++;
	  totalphases++;
	  if(calcgamma) {
		double gamma = DBL_MAX;
		for(int j=0; j<numD;j++) {
		  i= demand_index[j];
		  ve=paths[i][pathrecord[i]].end();
		  for(vi=paths[i][pathrecord[i]].begin(); vi != ve; vi++) {
			/*
			  if(minfree > vcapacity[*vi] - demand) {
			  minfree = vcapacity[*vi] - demand;
			  smallestcap = vcapacity[*vi];
			  }*/
			gweights[*vi] += origdemands[i];
		  }
		}		
		for(int j=0; j< M; j++) {
		  double tmp = 1.0 - gweights[j]/vcapacity[j];
		  if( gamma > tmp ) {
			gamma = tmp;
		  }
		}
		if(gamma >= intgamma * 0.9999999999) {
		  countgamma++;
		  //Rprintf("gamma %lg, matches int gamma of %lg\n",gamma, intgamma);
		}
		if(gamma > bestgamma) {
		  bestgamma = gamma;
		  bestpaths = pathrecord;
		  //Rprintf("bestgamma so far %lg\n",bestgamma);
		  //Rprintf("best paths: ");
		  //for(int j=0; j<numD; j++)
		  //Rprintf("%d,",bestpaths[j]);
		  //Rprintf("\n");
		}
	  }
	  
	}

	UNPROTECT(1);
	
	Rcpp::List retlist;
	retlist.push_back(pathCount,"pathcount");
	retlist.push_back(totalphases,"totalphases");
	retlist.push_back(countgamma,"countgamma");
	retlist.push_back(bestgamma,"bestgamma");
	retlist.push_back(bestpaths,"bestpaths");
	retlist.push_back(pathflows,"pathflows");
	retlist.push_back(vlengths,"vlengths");
	return Rcpp::wrap(retlist);
  }


  SEXP rg_max_concurrent_flow_int_c
  (
   SEXP RRdemandpaths,
   SEXP Rdemands,
   SEXP Rcapacity,
   SEXP Re,
   SEXP Rprogress,
   SEXP pb,
   SEXP Rintgamma,
   SEXP Rbestpaths,
   SEXP env
   ) {

	bool progress=Rcpp::as<bool>(Rprogress);
	double e = Rcpp::as<double>(Re);
	double intgamma = Rcpp::as<double>(Rintgamma);
	bool calcgamma=false;
	if(intgamma >=0) 
	  calcgamma=true;
	std::vector<double> vcapacity = RcppVector<double>(Rcapacity).stlVector();
	std::vector<double> vdemands = RcppVector<double>(Rdemands).stlVector();
	std::vector<double> origdemands(vdemands);	

	int numD = vdemands.size();
	std::vector<double>::iterator vecdit;
	int M = vcapacity.size();
	//decode Rdemandpaths into paths
	std::vector<int>::iterator vi,ve;
	std::vector< std::vector< std::vector<int> > > paths;
	std::vector< std::vector< int > > pathCount;
	std::vector< std::vector< std::vector<int> > >::iterator dit;
	std::vector< std::vector<int> >::iterator pit;
	std::vector<int>::iterator eit;

	Rcpp::List demands(RRdemandpaths);
	int sz = demands.size();
	paths.resize(sz);
	pathCount.resize(sz);
	for(int i=0;i<sz;i++) {
	  Rcpp::List demandpaths((SEXPREC*)demands[i]);
	  int sz = demandpaths.size();
	  paths[i].resize(sz);
	  pathCount[i].resize(sz);
	  for(int j=0;j<sz;j++) {
		paths[i][j]= Rcpp::as< std::vector<int> >(demandpaths[j]);
		pathCount[i][j]=0;
	  }
	}

	std::vector< std::vector<double> > pathflows(paths.size());
	for(int i=0;i<numD;i++) {
	  pathflows[i] = std::vector<double>(paths[i].size(),0.0);
	}

	int doubleCount=0;
	int phases=0;
	int phases2=0;
	int totalphases=0;

	double delta = pow(double(M) / (1.0 - e),-1.0/e);
	std::vector< std::vector<double> > dlengths(numD);
	for(int j=0; j<numD; j++) {
	  dlengths[j]=std::vector<double>(M);
	  for(int i=0;i<M;i++) {
		dlengths[j][i]=delta / vcapacity[j];
	  }
	}
	std::vector< std::vector<double> > dlengths2(numD);
	for(int j=0; j<numD; j++) {
	  dlengths2[j]=std::vector<double>(M);
	  for(int i=0;i<M;i++) {
		dlengths2[j][i]=delta / vcapacity[j];
	  }
	}

	std::vector< std::vector<double> > dlengths3(numD);
	for(int j=0; j<numD; j++) {
	  dlengths3[j]=std::vector<double>(M);
	  for(int i=0;i<M;i++) {
		dlengths3[j][i]=delta;
	  }
	}

	std::vector<double> vlengths(M);
	for(int i=0;i<M;i++) {
	  vlengths[i]=delta / vcapacity[i];
	}
	std::vector<double> vlengths2(M);
	for(int i=0;i<M;i++) {
	  vlengths2[i]=delta / vcapacity[i];
	}

	int doubreq =  int(ceil(2.0/e * log(M/(1-e))/log(1+e)));
	int updatepb = int(ceil(doubreq / 100.0));

	double D=0.0;
	for(int i=0;i<M;i++) {
	  D+=vcapacity[i]*vlengths[i];
	}
	char cmd[256];
	SEXP cmdsxp, cmdexpr, ansxp = R_NilValue;
	ParseStatus status;

	PROTECT(cmdsxp = allocVector(STRSXP,1));
	int i=0;
	int countgamma=0;
	double bestgamma= -DBL_MAX;
	double bestgammanormal= -DBL_MAX;
	double bestgammaindividual= -DBL_MAX;
	double bestgammajustgamma= -DBL_MAX;

	std::vector<int> demand_index(numD);
	for(int i=0; i<numD; i++) {
	  demand_index[i]=i;
	}
	std::vector<int> bestpaths;
	bool gammanormal=true;
	bool gammaindividual=false;//true;
	bool gammajustgamma=false;//true;
	// does not appear to work very well. (less_cost better than great_cost)?
	bool sort_order=true;
	std::vector<int> countunderv(1);

	int countnormal=0;
	int countindividual=0;
	int countjustgamma=0;

	countunderv[0]=0;
	while( D < 1.0) {
	  if(phases > doubreq) {
		Rprintf("doubling!!\n");
		for(int i=0;i<numD;i++) {
		  vdemands[i] = vdemands[i] *2;
		}
		
		phases = 0;
	  }
	  std::vector< std::pair<int,double> > orderv(numD);
	  if(sort_order) {
		for(int m=0;m<numD;m++) {
		  orderv[m].first=m;
		}
		random_shuffle(orderv.begin(),orderv.end());
	  }
	  

	  if(progress != false && totalphases % updatepb == 0) {
		sprintf(cmd,"setTxtProgressBar(pb,%d)",totalphases);
		SET_STRING_ELT(cmdsxp,0,mkChar(cmd));
		cmdexpr = PROTECT(R_ParseVector(cmdsxp, -1, &status, R_NilValue));
		for(int i = 0; i < length(cmdexpr); i++)
		  ansxp = eval(VECTOR_ELT(cmdexpr, i), env);
		//Rprintf("D=%lg\n",D);
		UNPROTECT(1);
	  }
	  
	  bool underdemand=true;

	  // note fix for random rotation of demand order

	  
	  std::vector<int> pathrecord(numD);
	  std::vector<int> pathrecord2(numD);
	  std::vector<int> pathrecord3(numD);
	  std::vector<double> weights(M,0.0);
	  std::vector<double> weights2(M,0.0);
	  std::vector<double> weights3(M,0.0);
	  std::vector<int> recdemands(0);


	  // OK for moment, but can be optimised later (could be
	  // done before start and only one weight updated each time)
	  if(gammajustgamma) {
		for(int j=0; j<numD;j++,i++) {
		  i= demand_index[j];
		  int pd = findshortestpath(paths[i],vlengths2);
		  pathrecord3[i]=pd;
		  ve=paths[i][pd].end();
		  for(vi=paths[i][pd].begin(); vi != ve; vi++) {
			weights3[*vi] += origdemands[i];
		  }
		}
	  }

	  
	  if(sort_order) {
		std::vector< std::pair<int,double> >::iterator pit,pend;
		pend=orderv.end();
		for(pit=orderv.begin(); pit != orderv.end();pit++) {
		  i=(*pit).first;
		  std::pair<long, double> p = findshortestpathcost(paths[i],vlengths);
		  (*pit).second = p.second;
		  //Rprintf("Cost for %ld=%lg\n",j,cost);
		}
		
		std::sort(orderv.begin(),orderv.end(),less_cost());
		for(int m=0;m<numD;m++) {
		  demand_index[m]=orderv[m].first;
		}
	  } else {
		//random_shuffle(demand_index.begin(),demand_index.end());
	  }
	  
	  for(int j=0; j<numD;j++,i++) {

		if(sort_order && false) {
		  std::vector< std::pair<int,double> >::iterator pit,pend;
		  pend=orderv.end();
		  for(pit=orderv.begin(); pit != pend;pit++) {
			long it=(*pit).first;
			std::pair<long, double> p = findshortestpathcost(paths[it],vlengths);
			(*pit).second = p.second;
		  }
		  
		  std::sort(orderv.begin(),orderv.end(),less_cost());
		  i = orderv[0].first;
		  orderv.erase(orderv.begin());
		} else {
		  i= demand_index[j];
		}
		recdemands.push_back(i);
		double demand = vdemands[i];
		double demand2 = vdemands[i];
		if(D >= 1.0)
		  underdemand=false;
		bool onlyfirst=true;
		int countunder=0;
		while( D < 1.0 && demand > 0.0) {
		  /*
		  std::pair<long, double> ppair = 
			findshortestpathcostopt(paths[i],
									vlengths,
									vcapacity,
									weights2,
									vdemands[i]);*/
		  std::pair<long, double> ppair2 = 
			findshortestpathcost(paths[i],
								 vlengths);
		  /*
		  Rprintf("%ld, %lg, %ld, %lg, %lg\n",
				  tmppair.first,
				  tmppair.second,
				  ppair.first,
				  ppair.second,tmppair.second - ppair.second);*/

		  int p = ppair2.first;
		  if(i == -1) {
			int testpath = 7;
			double path6cost=0;
			ve=paths[i][testpath].end();
			for(vi=paths[i][testpath].begin(); vi != ve; vi++) {
			  path6cost += vlengths[*vi];
			}
			Rprintf("%ld,\t %lg,\t %lg,\t %lg\n",
					p,path6cost,ppair2.second,path6cost/ppair2.second);
		  }
			//int p = findshortestpath(paths[i],vlengths);
		  pathCount[i][p]++;
		  //find minimum of demand or capacity on path
		  double mincap=demand;
		  ve=paths[i][p].end();
		  for(vi=paths[i][p].begin(); vi != ve; vi++) {
			if(mincap > vcapacity[*vi])
			  mincap = vcapacity[*vi];
		  }

		  demand -= mincap;
		  
		  if(demand <=0) {
			pathrecord2[i]=p;
			for(vi=paths[i][p].begin(); vi != ve; vi++) {
				weights2[*vi] += origdemands[i];
			}
		  }		  
		  int sz = paths[i][p].size();
		  for(int k=0; k < sz; k++) {
			double length = vlengths[paths[i][p][k]];
			length *= (1.0 + (e * mincap) /
					   vcapacity[paths[i][p][k]]);
			vlengths[paths[i][p][k]] = length;
		  }

		  // New
		  if(demand > 0.0) {
			countunder++;
		  } else {
			
			if(gammaindividual) {
			  int pd = findshortestpath(paths[i],dlengths[i]);
			
			  pathrecord[i]=pd;
			  double mincap=origdemands[i];
			  double minfree=DBL_MAX;
			  double smallestcap=0;
			  double penalty=0;
			  ve=paths[i][pd].end();
			  for(vi=paths[i][pd].begin(); vi != ve; vi++) {
				if(mincap > vcapacity[*vi])
				  mincap = vcapacity[*vi];
				
				if(minfree > vcapacity[*vi] - weights[*vi] - origdemands[i]) {
				  minfree = vcapacity[*vi] - weights[*vi] - origdemands[i];
				}
			  }
			  //if(sort_order) {
			  //orderv[i]=std::pair<int,double>(i,minfree);
			  //}
			  if(minfree < 0) {
				//penalty=-minfree;
				//penalty=mincap;
				penalty=origdemands[i];
				minfree=0;
			  }	
			  
			  ve=paths[i][pd].end();
			  for(vi=paths[i][pd].begin(); vi != ve; vi++) {
				weights[*vi] += origdemands[i];//minfree;
			  }
			// First update all lengths as normal G&K algorithm
			  sz = paths[i][pd].size();
			  
			  for(int l=0; l<numD;l++) {
				for(int k=0; k < sz; k++) {
				  double dlength = dlengths[l][paths[i][pd][k]];
				  dlength *= (1.0 + (e * (origdemands[i])) /
							  vcapacity[paths[i][pd][k]]);
				  dlengths[l][paths[i][pd][k]] = dlength;
				}
			  }
			  // then update just the length set for this demand 
			  // to allow for congested flow
			  bool newupdate = false;
			  if(newupdate==true) {
				int l=i;
				for(int m=0;m<recdemands.size();m++) {
				  l=recdemands[m];
				  for(int k=0; k < sz; k++) {
					double dlength = dlengths[l][paths[i][pd][k]];
					//if( (vcapacity[paths[i][pd][k]] - weights [paths[i][pd][k]] )
					//	< 0) {
					if(minfree == 0) {
					  
					  double val = weights [paths[i][pd][k]] - vcapacity[paths[i][pd][k]];
					  // Rprintf("val=%lg\n",val);
					  dlength *= (1.0 + (e * (penalty)) /
								  vcapacity[paths[i][pd][k]]);
					  dlengths[l][paths[i][pd][k]] = dlength;
					}
				  }
				}
			  }
			}
		  } // end if(demand > 0.0)
			if(gammajustgamma) {
			  int pd=pathrecord3[i];
			  ve=paths[i][pd].end();
			  for(vi=paths[i][pd].begin(); vi != ve; vi++) {
				weights3[*vi] -= origdemands[i];
			  }
			  pd = findshortestpath(paths[i],vlengths2);
			  ve=paths[i][pd].end();
			  pathrecord3[i]=pd;
			  for(vi=paths[i][pd].begin(); vi != ve; vi++) {
				weights3[*vi] += origdemands[i];
			  }

			  double penalty=0;
			  double minfree=DBL_MAX;
			  double mincap=demand2;
			  ve=paths[i][pd].end();
			  for(vi=paths[i][pd].begin(); vi != ve; vi++) {
				if(mincap > vcapacity[*vi])
				  mincap = vcapacity[*vi];
				
				if(minfree > vcapacity[*vi] - weights3[*vi] - origdemands[i]) {
				  minfree = vcapacity[*vi] - weights3[*vi] - origdemands[i];
				}
			  }
			  demand2 -= mincap;
			  if(minfree < 0) {
				//penalty=-minfree;
				//penalty=mincap;
				penalty=origdemands[i];
				minfree=0;
			  }	


			  double gamma=DBL_MAX;
			  for(int jl=0; jl< M; j++) {
				double tmp = 1.0 - weights3[jl]/vcapacity[jl];
				if( gamma > tmp ) {
				  gamma = tmp;
				}
			  }
			  //			  Rprintf("gamma=%lg\n",gamma);
			  ve=paths[i][pd].end();
			  for(vi=paths[i][pd].begin(); vi != ve; vi++) {
				vlengths2[*vi]  *= (1.0 + (e * (weights3[j]))/vcapacity[*vi]);
			  }
			} //end if(true)
		  //end New
		  ///!!!! this needs to move earlier???
		  pathflows[i][p] += mincap;
		  D=0.0;
		  for(int k=0;k<M;k++) {
			D+=vcapacity[k]*vlengths[k];
		  }
		}
		if(countunderv.size() < countunder+1) {
		  countunderv.resize(countunder+1);
		}
		countunderv[countunder]++;

	  }// end for each demand


	  
	  // this is just used to calculate gamma

	  //if(sort_order)
	  //std::sort(orderv.begin(),orderv.end(),less_cost());
	  //	  for(int m=0;m<numD;m++) {
	  //Rprintf("%ld, %lg,",orderv[m].first,orderv[m].second);
	  //}
	  
	  //Rprintf("\n");

	  std::vector<double> gweights(M,0.0);
	  if(gammaindividual) {

		for(int j=0; j<numD;j++) {
		  i= demand_index[j];
		  ve=paths[i][pathrecord[i]].end();
		  for(vi=paths[i][pathrecord[i]].begin(); vi != ve; vi++) {
			/*
			  if(minfree > vcapacity[*vi] - demand) {
			  minfree = vcapacity[*vi] - demand;
			  smallestcap = vcapacity[*vi];
			  }*/
			gweights[*vi] += origdemands[i];
		  }
		}
	  }	  

	  std::vector<double> gweights2(M,0.0);	  
	  if(gammanormal) {

		for(int j=0; j<numD;j++) {
		  i= demand_index[j];
		  ve=paths[i][pathrecord2[i]].end();
		  for(vi=paths[i][pathrecord2[i]].begin(); vi != ve; vi++) {
			/*
			  if(minfree > vcapacity[*vi] - demand) {
			  minfree = vcapacity[*vi] - demand;
			  smallestcap = vcapacity[*vi];
			  }*/
			gweights2[*vi] += origdemands[i];
		  }
		}
	  }	  


	  std::vector<double> gweights3(M,0.0);	  
	  if(gammajustgamma) {
		for(int j=0; j<numD;j++) {
		  i= demand_index[j];
		  ve=paths[i][pathrecord3[i]].end();
		  for(vi=paths[i][pathrecord3[i]].begin(); vi != ve; vi++) {
			gweights3[*vi] += origdemands[i];
		  }
		}
	  }
	  
	  phases++;
	  totalphases++;
	  if(calcgamma && underdemand) {
		double gamma = DBL_MAX;
		double gamma2 = DBL_MAX;
		double gamma3 = DBL_MAX;
		if(gammaindividual) {
		  for(int j=0; j< M; j++) {
			double tmp = 1.0 - gweights[j]/vcapacity[j];
			if( gamma > tmp ) {
			  gamma = tmp;
			}
		  }
		  if(gamma > bestgamma) {
			bestgamma = gamma;
			bestpaths = pathrecord;
			//Rprintf("bestgamma so far %lg\n",bestgamma);
			//Rprintf("best paths: ");
			//for(int j=0; j<numD; j++)
			//Rprintf("%d,",bestpaths[j]);
			//Rprintf("\n");
		  }
		  if(gamma > bestgammaindividual)
			bestgammaindividual= gamma;
		  if(gamma >= intgamma * 0.9999999999) {
			countindividual++;
		  }
		} else
		  gamma = -DBL_MAX;
		if(gammanormal) {
		  for(int j=0; j< M; j++) {
			double tmp = 1.0 - gweights2[j]/vcapacity[j];
			if( gamma2 > tmp ) {
			  gamma2 = tmp;
			}
		  }
		  if(gamma2 > bestgamma) {
			bestgamma = gamma2;
			bestpaths = pathrecord2;
			//Rprintf("bestgamma so far %lg\n",bestgamma);
			//Rprintf("best paths: ");
			//for(int j=0; j<numD; j++)
			//Rprintf("%d,",bestpaths[j]);
			//Rprintf("\n");
		  }
		  if(gamma2 > bestgammanormal)
			bestgammanormal= gamma2;
		  if(gamma2 >= intgamma * 0.9999999999) {
			countnormal++;
		  }
		  
		} else
		  gamma2= -DBL_MAX;
		if(gammajustgamma) {
		  for(int j=0; j< M; j++) {
			double tmp = 1.0 - gweights3[j]/vcapacity[j];
			if( gamma3 > tmp ) {
			  gamma3 = tmp;
			}
		  }
		  //Rprintf("gamma3=%lg\n",gamma3);
		  if(gamma3 > bestgamma) {
			bestgamma = gamma3;
			bestpaths = pathrecord3;
			//Rprintf("bestgamma so far %lg\n",bestgamma);
			//Rprintf("best paths: ");
			//for(int j=0; j<numD; j++)
			//Rprintf("%d,",bestpaths[j]);
			//Rprintf("\n");
		  }
		  if(gamma3 > bestgammajustgamma)
			bestgammajustgamma = gamma3;

		  if(gamma3 >= intgamma * 0.9999999999) {
			countjustgamma++;
		  }
			
		} else
		  gamma3 = -DBL_MAX;
		if(gamma2 >= intgamma * 0.9999999999 || 
		   gamma >= intgamma * 0.9999999999 ||
		   gamma3 >= intgamma * 0.9999999999 ) {
		  countgamma++;
		  //Rprintf("gamma %lg, matches int gamma of %lg\n",gamma, intgamma);
		}
		
	  }

	}

	Rprintf("\n");
	for(i=0;i<countunderv.size();i++) {
	  Rprintf("%ld,",countunderv[i]);
	}
	Rprintf("; normal=%ld, individual=%ld, gamma=%ld\n",
			countnormal, countindividual,countjustgamma);
	Rprintf("; normal=%lg, individual=%lg, gamma=%lg\n",
			bestgammanormal, bestgammaindividual,bestgammajustgamma);
	
	UNPROTECT(1);
	
	Rcpp::List retlist;
	retlist.push_back(pathCount,"pathcount");
	retlist.push_back(totalphases,"totalphases");
	retlist.push_back(countgamma,"countgamma");
	retlist.push_back(bestgamma,"bestgamma");
	retlist.push_back(bestpaths,"bestpaths");
	retlist.push_back(pathflows,"pathflows");
	retlist.push_back(vlengths,"vlengths");
	return Rcpp::wrap(retlist);
  }



}

