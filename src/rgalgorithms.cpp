#include <boost/lexical_cast.hpp>
#include <R.h>
#include <Rcpp.h>
#include <Rinternals.h>
#include "rgalgorithms.h"
#undef __BEGIN_DECLS
#undef __END_DECLS
#ifdef __cplusplus
# define __BEGIN_DECLS extern "C" {
# define __END_DECLS }
#else
# define __BEGIN_DECLS /* empty */
# define __END_DECLS /* empty */
#endif



void make_heap(std::vector< std::pair<double, std::vector<int> > >::iterator,
	       std::vector< std::pair<double, std::vector<int> > >::iterator,
	       bool operator()(std::pair<double, std::vector<int> > a,
			       std::pair<double, std::vector<int> > b)
	       );


struct less_demand {
public:
  bool operator()(const std::pair<double, std::vector<int> > a,
		  const std::pair<double, std::vector<int> > b) {
    return (a.first > b.first);
  }
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


__BEGIN_DECLS

void rg_cversion() {
  Rprintf("Version 2.0");
}


class rg_demandi {
public:
  double flow;
  long source;
  long sink;
  double demand;
  std::map< std::vector< long > ,double > path_flow_map;
};


void igraph_rg_fleischer_max_concurrent_flow(std::vector<long> vedges,
					     std::vector<double> capacities,
					     std::vector<rg_demandi> &demands,
					     std::vector<double> &lengths,
					     double e,
					     long N,
					     long &totalphases);



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
  PROTECT(count = Rf_allocMatrix(INTSXP,N,N));
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
  bool updateflow=Rf_asLogical(Rupdateflow);
  std::string title("");
  if ( Rf_isLogical(Rprogress) ) {
    progress=Rf_asLogical(Rprogress);
  } else if ( Rf_isString(Rprogress) ) {
    title = std::string(CHAR(STRING_ELT(Rprogress,0)));
    progress=true;
  }
  Graph_rg gdual(num_verts_in, num_edges_in, R_edges_in, R_weights_in,capacities_in);
	
  double* dem_in = REAL(demands_in);
  int* dem_sources_in = INTEGER(demands_sources_in);
  int* dem_sinks_in = INTEGER(demands_sinks_in);
  int num_dem = Rf_asInteger(num_demands_in);
  std::vector<rg_demand> demands(num_dem);
  for(int i=0 ; i<num_dem; i++) {
    demands[i].demand = dem_in[i];
    demands[i].flow = 0;
    demands[i].source = vertex(dem_sources_in[i],gdual);
    demands[i].sink = vertex(dem_sinks_in[i],gdual);
  }
  int N = num_vertices(gdual);
  int m = num_edges(gdual);
  double e = Rf_asReal(Re);
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
  int doubreq =  2*int(ceil(1.0/e * log(m/(1-e))/log(1+e)));
	
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

  PROTECT(cmdsxp = Rf_allocVector(STRSXP,1));
	
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
      SET_STRING_ELT(cmdsxp,0,Rf_mkChar(cmd));
      cmdexpr = PROTECT(R_ParseVector(cmdsxp, -1, &status, R_NilValue));
      for(int i = 0; i < Rf_length(cmdexpr); i++)
	ansxp = Rf_eval(VECTOR_ELT(cmdexpr, i), env);
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
  PROTECT(duallenxp = Rf_allocVector(REALSXP,m));
  PROTECT(flowsxp = Rf_allocVector(REALSXP,num_dem));
  PROTECT(lenxp = Rf_allocVector(INTSXP,3));
  PROTECT(pmapkeyxp = Rf_allocVector(VECSXP,num_dem_p_flows * 2));
  PROTECT(pmapvalxp = Rf_allocVector(REALSXP,num_dem_p_flows));
  PROTECT(retlistxp = Rf_allocVector(VECSXP,5));

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
      PROTECT(demsxp = Rf_allocVector(STRSXP,1));
      SET_STRING_ELT(demsxp,0,Rf_mkChar(sdem.c_str()));
      // create the path string
      PROTECT(pathsxp = Rf_allocVector(STRSXP,1));
      SET_STRING_ELT(pathsxp,0,Rf_mkChar(spath.c_str()));
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
  random_shuffle(trypath.begin(),trypath.end());

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
  //return (std::pair<int,double>(minp,min));
	
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
  std::vector<double> vcapacity = Rcpp::as< std::vector<double> >(Rcapacity);
  std::vector<double> vdemands = Rcpp::as< std::vector<double> >(Rdemands);
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

  int doubreq =  2*int(ceil(1.0/e * log(M/(1-e))/log(1+e)));
  int updatepb = int(ceil(doubreq / 100.0));

  double D=0.0;
  for(int i=0;i<M;i++) {
    D+=vcapacity[i]*vlengths[i];
  }
  char cmd[256];
  SEXP cmdsxp, cmdexpr, ansxp = R_NilValue;
  ParseStatus status;

  PROTECT(cmdsxp = Rf_allocVector(STRSXP,1));
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
      SET_STRING_ELT(cmdsxp,0,Rf_mkChar(cmd));
      cmdexpr = PROTECT(R_ParseVector(cmdsxp, -1, &status, R_NilValue));
      for(int i = 0; i < Rf_length(cmdexpr); i++)
	ansxp = Rf_eval(VECTOR_ELT(cmdexpr, i), env);
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
  std::vector<double> vcapacity = Rcpp::as< std::vector<double> >(Rcapacity);
  std::vector<double> vdemands = Rcpp::as< std::vector<double> >(Rdemands);
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

  int doubreq =  2*int(ceil(1.0/e * log(M/(1-e))/log(1+e)));
  int updatepb = int(ceil(doubreq / 100.0));

  double D=0.0;
  for(int i=0;i<M;i++) {
    D+=vcapacity[i]*vlengths[i];
  }
  char cmd[256];
  SEXP cmdsxp, cmdexpr, ansxp = R_NilValue;
  ParseStatus status;

  PROTECT(cmdsxp = Rf_allocVector(STRSXP,1));
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
      SET_STRING_ELT(cmdsxp,0,Rf_mkChar(cmd));
      cmdexpr = PROTECT(R_ParseVector(cmdsxp, -1, &status, R_NilValue));
      for(int i = 0; i < Rf_length(cmdexpr); i++)
	ansxp = Rf_eval(VECTOR_ELT(cmdexpr, i), env);
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





// only used to inspect the cost of a path while investigating the aglorithm
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
					SEXP Rpermutation,
					SEXP Rdeltaf
					) {
  using namespace boost;
  bool progress=false;
  bool updateflow=Rcpp::as<bool>(Rupdateflow);
  std::string title("");
  if ( Rf_isLogical(Rprogress) ) {
    progress=Rcpp::as<bool>(Rprogress);
  } else if ( Rf_isString(Rprogress) ) {
    title = std::string(CHAR(STRING_ELT(Rprogress,0)));
    progress=true;
  }
  Graph_rg gdual(num_verts_in, num_edges_in, R_edges_in, R_weights_in,capacities_in);

  std::vector<int> permutation = Rcpp::as< std::vector<int> > (Rpermutation);
  double* dem_in = REAL(demands_in);
  int* dem_sources_in = INTEGER(demands_sources_in);
  int* dem_sinks_in = INTEGER(demands_sinks_in);
  int num_dem = Rcpp::as<long>(num_demands_in);
  double deltaf = Rcpp::as<double>(Rdeltaf);
  std::vector<rg_demand> demands(num_dem);
  for(int i=0 ; i<num_dem; i++) {
    demands[i].demand = dem_in[i];
    demands[i].flow = 0;
    demands[i].source = vertex(dem_sources_in[i],gdual);
    demands[i].sink = vertex(dem_sinks_in[i],gdual);
  }
  int N = num_vertices(gdual);
  int m = num_edges(gdual);
  double e = Rcpp::as<double>(Re);
  int num_dem_ed_flows =0;
  int num_dem_p_flows =0;
  double delta = pow(double(m) / (1.0 - e),-1.0/e) * deltaf;

  graph_traits < Graph_rg >::edge_iterator ei, eend;

#ifdef RECORDPATH
  std::vector< std::vector< std::vector<int> > > paths;

  std::vector< std::vector< int > > pathCount;
  int sz = demands.size();
  paths.resize(sz);
  pathCount.resize(sz);
#endif


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
  int doubreq =  2*int(ceil(1.0/e * log(m/(1-e))/log(1+e)));
	
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

  PROTECT(cmdsxp = Rf_allocVector(STRSXP,1));
	
  char cmd[256];
  //ansxp = Rf_eval(cmdexpr,env);

  double a = pow(2.0,-130);
  double b =1;
  double c = a +b;

  std::vector<int> demand_index(num_dem);
  if(permutation[0] >= 0) {
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

  if(permutation[0] == -2) {
    // sort demands on lowest cost path first
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
      SET_STRING_ELT(cmdsxp,0,Rf_mkChar(cmd));
      cmdexpr = PROTECT(R_ParseVector(cmdsxp, -1, &status, R_NilValue));
      for(int i = 0; i < Rf_length(cmdexpr); i++)
	ansxp = Rf_eval(VECTOR_ELT(cmdexpr, i), env);
      //Rprintf("D=%lg\n",D);
      UNPROTECT(1);
    }	  
    if(permutation[0] == -2) {
      // sort demands on lowest cost path first
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
    if(permutation[0] == -1) {
      random_shuffle(demand_index.begin(),demand_index.end());
    }

    for(int j=0; j<num_dem;j++) {
      int i;
      if(permutation[0] == -2)
	i = costpair[j].first;
      else
	i= demand_index[j];

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

#ifdef RECORDPATH
	std::vector<int> tmppath(0);
	tmppath.push_back(f);
#endif		  
	f = p;
	p = penult[p];
	double w =get(edge_weight,gdual,ed.first);
		  
	while(f != source) {

#ifdef RECORDPATH
	  tmppath.push_back(f);
#endif
	  ed = edge(p,f,gdual);
	  double cap =get(edge_capacity,gdual,ed.first);
	  double w =get(edge_weight,gdual,ed.first);
	  mincap = cap < mincap? cap : mincap;
	  f = p;
	  p = penult[p];
	}
#ifdef RECORDPATH
	tmppath.push_back(f);

	// keep this here. it is useful to inspect what is happening.
	// can be removed later
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
#endif

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

#ifdef RECORDPATH
  for(int i=0;i<paths.size();i++) {
    Rprintf("Demand %ld\n",i+1);
    for(int j=0;j<paths[i].size();j++) {
      for(int k=paths[i][j].size()-1;k>=0;k--) {
	Rprintf("%ld,",paths[i][j][k]+1);
      }
      Rprintf(": %ld\n",pathCount[i][j]);
    }
  }
#endif

  // need to return:
  // dual lengths
  // demand[].flow
  // demand[].[pathmap] key and values

  // set up return expressions
  SEXP duallenxp, flowsxp, lenxp, pmapkeyxp, pmapvalxp, retlistxp;
  PROTECT(duallenxp = Rf_allocVector(REALSXP,m));
  PROTECT(flowsxp = Rf_allocVector(REALSXP,num_dem));
  PROTECT(lenxp = Rf_allocVector(INTSXP,3));
  PROTECT(pmapkeyxp = Rf_allocVector(VECSXP,num_dem_p_flows * 2));
  PROTECT(pmapvalxp = Rf_allocVector(REALSXP,num_dem_p_flows));
  PROTECT(retlistxp = Rf_allocVector(VECSXP,5));

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
      PROTECT(demsxp = Rf_allocVector(STRSXP,1));
      SET_STRING_ELT(demsxp,0,Rf_mkChar(sdem.c_str()));
      // create the path string
      PROTECT(pathsxp = Rf_allocVector(STRSXP,1));
      SET_STRING_ELT(pathsxp,0,Rf_mkChar(spath.c_str()));
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

SEXP rg_karakostas_max_concurrent_flow_c(SEXP num_verts_in,
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
					 SEXP Rpermutation,
					 SEXP Rdeltaf
					 ) {
  using namespace boost;
  bool progress=false;
  bool updateflow=Rcpp::as<bool>(Rupdateflow);
  std::string title("");
  if ( Rf_isLogical(Rprogress) ) {
    progress=Rcpp::as<bool>(Rprogress);
  } else if ( Rf_isString(Rprogress) ) {
    title = std::string(CHAR(STRING_ELT(Rprogress,0)));
    progress=true;
  }
  Graph_rg gdual(num_verts_in, num_edges_in, R_edges_in, R_weights_in,capacities_in);

  std::vector<int> permutation = Rcpp::as< std::vector<int> > (Rpermutation);
  double* dem_in = REAL(demands_in);
  int* dem_sources_in = INTEGER(demands_sources_in);
  int* dem_sinks_in = INTEGER(demands_sinks_in);
  int num_dem = Rcpp::as<long>(num_demands_in);
  double deltaf = Rcpp::as<double>(Rdeltaf);
  std::vector<rg_demand> demands(num_dem);
  for(int i=0 ; i<num_dem; i++) {
    demands[i].demand = dem_in[i];
    demands[i].flow = 0;
    demands[i].source = vertex(dem_sources_in[i],gdual);
    demands[i].sink = vertex(dem_sinks_in[i],gdual);
  }
  int N = num_vertices(gdual);
  int m = num_edges(gdual);
  double e = Rcpp::as<double>(Re);
  int num_dem_ed_flows =0;
  int num_dem_p_flows =0;
  double delta = pow(double(m) / (1.0 - e),-1.0/e) * deltaf;

  graph_traits < Graph_rg >::edge_iterator ei, eend;

#ifdef RECORDPATH
  std::vector< std::vector< std::vector<int> > > paths;

  std::vector< std::vector< int > > pathCount;
  int sz = demands.size();
  paths.resize(sz);
  pathCount.resize(sz);
#endif


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
  int doubreq =  2*int(ceil(1.0/e * log(m/(1-e))/log(1+e)));
	
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

  PROTECT(cmdsxp = Rf_allocVector(STRSXP,1));
	
  char cmd[256];
  //ansxp = Rf_eval(cmdexpr,env);

  double a = pow(2.0,-130);
  double b =1;
  double c = a +b;

  std::vector<int> demand_index(num_dem);
  if(permutation[0] >= 0) {
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

  if(permutation[0] == -2) {
    // sort demands on lowest cost path first
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
      SET_STRING_ELT(cmdsxp,0,Rf_mkChar(cmd));
      cmdexpr = PROTECT(R_ParseVector(cmdsxp, -1, &status, R_NilValue));
      for(int i = 0; i < Rf_length(cmdexpr); i++)
	ansxp = Rf_eval(VECTOR_ELT(cmdexpr, i), env);
      //Rprintf("D=%lg\n",D);
      UNPROTECT(1);
    }	  
    if(permutation[0] == -2) {
      // sort demands on lowest cost path first
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
    if(permutation[0] == -1) {
      random_shuffle(demand_index.begin(),demand_index.end());
    }

    for(int j=0; j<num_dem;j++) {
      int i;
      if(permutation[0] == -2)
	i = costpair[j].first;
      else
	i= demand_index[j];

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

#ifdef RECORDPATH
	std::vector<int> tmppath(0);
	tmppath.push_back(f);
#endif		  
	f = p;
	p = penult[p];
	double w =get(edge_weight,gdual,ed.first);
		  
	while(f != source) {

#ifdef RECORDPATH
	  tmppath.push_back(f);
#endif
	  ed = edge(p,f,gdual);
	  double cap =get(edge_capacity,gdual,ed.first);
	  double w =get(edge_weight,gdual,ed.first);
	  mincap = cap < mincap? cap : mincap;
	  f = p;
	  p = penult[p];
	}
#ifdef RECORDPATH
	tmppath.push_back(f);

	// keep this here. it is useful to inspect what is happening.
	// can be removed later
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
#endif

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

#ifdef RECORDPATH
  for(int i=0;i<paths.size();i++) {
    Rprintf("Demand %ld\n",i+1);
    for(int j=0;j<paths[i].size();j++) {
      for(int k=paths[i][j].size()-1;k>=0;k--) {
	Rprintf("%ld,",paths[i][j][k]+1);
      }
      Rprintf(": %ld\n",pathCount[i][j]);
    }
  }
#endif

  // need to return:
  // dual lengths
  // demand[].flow
  // demand[].[pathmap] key and values

  // set up return expressions
  SEXP duallenxp, flowsxp, lenxp, pmapkeyxp, pmapvalxp, retlistxp;
  PROTECT(duallenxp = Rf_allocVector(REALSXP,m));
  PROTECT(flowsxp = Rf_allocVector(REALSXP,num_dem));
  PROTECT(lenxp = Rf_allocVector(INTSXP,3));
  PROTECT(pmapkeyxp = Rf_allocVector(VECSXP,num_dem_p_flows * 2));
  PROTECT(pmapvalxp = Rf_allocVector(REALSXP,num_dem_p_flows));
  PROTECT(retlistxp = Rf_allocVector(VECSXP,5));

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
      PROTECT(demsxp = Rf_allocVector(STRSXP,1));
      SET_STRING_ELT(demsxp,0,Rf_mkChar(sdem.c_str()));
      // create the path string
      PROTECT(pathsxp = Rf_allocVector(STRSXP,1));
      SET_STRING_ELT(pathsxp,0,Rf_mkChar(spath.c_str()));
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

  
  

SEXP rg_fleischer_max_concurrent_flow_c_boost(SEXP num_verts_in,
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
					      SEXP Rpermutation,
					      SEXP Rdeltaf
					      ) {
  using namespace boost;
  bool progress=false;
  bool updateflow=Rcpp::as<bool>(Rupdateflow);
  std::string title("");
  if ( Rf_isLogical(Rprogress) ) {
    progress=Rcpp::as<bool>(Rprogress);
  } else if ( Rf_isString(Rprogress) ) {
    title = std::string(CHAR(STRING_ELT(Rprogress,0)));
    progress=true;
  }
  Graph_rg gdual(num_verts_in, num_edges_in, R_edges_in, R_weights_in,capacities_in);
	
  std::vector<long> vedges = Rcpp::as< std::vector<long> > (R_edges_in);

  std::vector<double> capacities = Rcpp::as< std::vector<double> > (capacities_in);
  int N = Rcpp::as<int>(num_verts_in);
  int m = vedges.size()/2;

  std::vector<int> permutation = Rcpp::as< std::vector<int> > (Rpermutation);
  double* dem_in = REAL(demands_in);
  int* dem_sources_in = INTEGER(demands_sources_in);
  int* dem_sinks_in = INTEGER(demands_sinks_in);
  int num_dem = Rcpp::as<long>(num_demands_in);
  double deltaf = Rcpp::as<double>(Rdeltaf);
  std::vector<rg_demand> demands(num_dem);

  for(int i=0 ; i<num_dem; i++) {
    demands[i].demand = dem_in[i];
    demands[i].flow = 0;
    demands[i].source = vertex(dem_sources_in[i],gdual);
    demands[i].sink = vertex(dem_sinks_in[i],gdual);
  }


  double e = Rcpp::as<double>(Re);
  int num_dem_ed_flows =0;
  int num_dem_p_flows =0;
  double delta = pow(double(m) / (1.0 - e),-1.0/e) * deltaf;

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

  int doubreq =  2*int(ceil(1.0/e * log(m/(1-e))/log(1+e)));
	
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

  PROTECT(cmdsxp = Rf_allocVector(STRSXP,1));
	
  char cmd[256];
  //ansxp = Rf_eval(cmdexpr,env);

  double a = pow(2.0,-130);
  double b =1;
  double c = a +b;

  std::vector<int> demand_index(num_dem);
  if(permutation[0] >= 0) {
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

  if(permutation[0] == -2) {
    // sort demands on lowest cost path first
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
      SET_STRING_ELT(cmdsxp,0,Rf_mkChar(cmd));
      cmdexpr = PROTECT(R_ParseVector(cmdsxp, -1, &status, R_NilValue));
      for(int i = 0; i < Rf_length(cmdexpr); i++)
	ansxp = Rf_eval(VECTOR_ELT(cmdexpr, i), env);
      //Rprintf("D=%lg\n",D);
      UNPROTECT(1);
    }	  
    if(permutation[0] == -2) {
      // sort demands on lowest cost path first
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
    if(permutation[0] == -1) {
      random_shuffle(demand_index.begin(),demand_index.end());
    }

    for(int j=0; j<num_dem;j++) {
      int i;
      if(permutation[0] == -2)
	i = costpair[j].first;
      else
	i= demand_index[j];

      rg_demand demand=demands[i];
		
      rgVertex source = demand.source;
      rgVertex sink = demand.sink;
      rgVertex f,p;

      //iterations
      while( D < 1.0 && demand.demand > 0) {
	//while( Di < 1.0 && demandi.demand > 0) {
		  
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


  // need to return:
  // dual lengths
  // demand[].flow
  // demand[].[pathmap] key and values

  // set up return expressions
  std::vector<double> lengths(m);

  int* edges_in = INTEGER(R_edges_in);
  for (int i = 0;
       i < m ; 
       i++, edges_in += 2) {
    std::pair<rgEdge, bool> ed = edge(vertex(*edges_in,gdual),vertex(*(edges_in+1),gdual), gdual);
    lengths[i] = get(edge_weight,gdual,ed.first);
  }
  Rcpp::List retdemandsi;
  for(int i=0 ; i<num_dem; i++) {
    Rcpp::List demand;
    demand.push_back(demands[i].demand,"demand");
    demand.push_back(demands[i].flow,"flow");
    demand.push_back(boost::lexical_cast<std::string>(demands[i].source+1),
			     "source");
    demand.push_back(boost::lexical_cast<std::string>(demands[i].sink+1)
			     ,"sink");
    Rcpp::List paths;
    std::map<const std::vector<rgVertex>, double>::iterator pmi, pmend;
    for(pmi=demands[i].path_flow_map.begin();
	pmi != demands[i].path_flow_map.end(); pmi++) {
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
      paths.push_back(pmi->second,spath.c_str());
    }
    demand.push_back(paths,"paths");
    retdemandsi.push_back(demand,boost::lexical_cast<std::string> (i+1));

  }
  

  
  Rcpp::List retlist;
  retlist.push_back(retdemandsi,"demands");
  retlist.push_back(totalphases,"totalphases");
  retlist.push_back(lengths,"lengths");
  UNPROTECT(1);
  
  return wrap(retlist);
  
}



inline std::vector<double> 
calcgraphweights(std::vector< std::vector< std::vector<int> > > &paths,
		 std::vector<int> &pathrecord,
		 std::vector<double> &demands,
		 int M) {
  int numD=demands.size();
  std::vector<double> weights(M,0.0);
  std::vector<int>::iterator vi,ve;
	
  for(int i=0; i<numD;i++) {
    ve=paths[i][pathrecord[i]].end();
    for(vi=paths[i][pathrecord[i]].begin(); vi != ve; vi++) {
      weights[*vi] += demands[i];
    }
  }
  return(weights);
}

inline double calcgamma(std::vector<double> weights,
			std::vector<double> vcapacity) {
  int M = weights.size();
  double gamma = DBL_MAX;
  for(int j=0; j< M; j++) {
    double tmp = 1.0 - weights[j]/vcapacity[j];
    if( gamma > tmp ) {
      gamma = tmp;
    }
  }
  return(gamma);
}

double calcbeta_restricted(std::vector< std::vector< std::vector<int> > > &paths,
			   std::vector<double> &lengths,
			   std::vector<double> &demands,
			   std::vector<double> &capacity
			   ) {
  double alpha=0.0;
  long numD=paths.size();
  long M=lengths.size();
  for(int i=0; i< numD; i++) {
    std::pair<int,double> tmp=findshortestpathcost(paths[i],lengths);
    alpha+=demands[i] * tmp.second;
  }
  double D=0;
  for(int i=0;i<M;i++) {
    D+=capacity[i]*lengths[i];
  }

  double beta = D/alpha;
  return(beta);

}
double calcbeta_int(std::vector< std::vector< std::vector<int> > > &paths,
		    std::vector<double> &lengths,
		    std::vector<double> &demands,
		    std::vector<double> &capacity,
		    std::vector<int> pathsselected
		    ) {
  double alpha=0.0;
  long numD=paths.size();
  long M=lengths.size();
  for(int i=0; i< numD; i++) {
    double cost=0;
    for(int j=0; j< paths[i][pathsselected[i]].size(); j++) {
      cost += lengths[paths[i][pathsselected[i]][j]];
    }
    alpha+=demands[i] * cost;
  }
  double D=0;
  for(int i=0;i<M;i++) {
    D+=capacity[i]*lengths[i];
  }

  double beta = D/alpha;
  return(beta);

}
SEXP rg_max_concurrent_flow_int_c
(
 SEXP RRdemandpaths,
 SEXP Rdemands,
 SEXP Rcapacity,
 SEXP Re,
 SEXP ReInternal,
 SEXP Rprogress,
 SEXP pb,
 SEXP Rintgamma,
 SEXP env,
 SEXP Rpermutation,
 SEXP Rdeltaf
 ) {

  bool progress=Rcpp::as<bool>(Rprogress);
  double e = Rcpp::as<double>(Re);
  double eInternal = Rcpp::as<double>(ReInternal);
  double intgamma = Rcpp::as<double>(Rintgamma);
  std::vector<double> vcapacity = Rcpp::as< std::vector<double> > (Rcapacity);
  std::vector<double> vdemands = Rcpp::as< std::vector<double> > (Rdemands);
  std::vector<double> origdemands(vdemands);	
  std::vector<int> permutation = Rcpp::as< std::vector<int> > (Rpermutation);
  double deltaf = Rcpp::as<double>(Rdeltaf);

  int numD = vdemands.size();
  std::vector<double>::iterator vecdit;
  int M = vcapacity.size();
  //decode Rdemandpaths into paths
  std::vector<int>::iterator vi,vj,ve;
  std::vector< std::vector< std::vector<int> > > paths;
  std::vector< std::vector< int > > pathCount;
  std::vector< std::vector< std::vector<int> > >::iterator dit;
  std::vector< std::vector<int> >::iterator pit;
  std::vector<int>::iterator eit;
  std::vector<double> gammavals;
  std::vector<double> lambdavals;
  std::vector<double> betavals;
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

  double delta = pow(double(M) / (1.0 - e),-1.0/e) * deltaf;
  std::vector< std::vector<double> > dlengths(numD);
  for(int j=0; j<numD; j++) {
    dlengths[j]=std::vector<double>(M);
    for(int i=0;i<M;i++) {
      dlengths[j][i]=delta / vcapacity[j];
    }
  }
  std::vector<double> vlengths(M);
  for(int i=0;i<M;i++) {
    vlengths[i]=delta / vcapacity[i];
  }

  //int doubreq =  2*int(ceil(1.0/eInternal * log(M/(1-eInternal))/log(1+eInternal)));
  int doubreq =  2*int(ceil(log(1.0/delta)/log(1+eInternal)));
  int updatepb = int(ceil(doubreq / 100.0));

  double D=0.0;
  for(int i=0;i<M;i++) {
    D+=vcapacity[i]*vlengths[i];
  }
  char cmd[256];
  SEXP cmdsxp, cmdexpr, ansxp = R_NilValue;
  ParseStatus status;

  PROTECT(cmdsxp = Rf_allocVector(STRSXP,1));
  int i=0;
  int countgamma=0;
  double bestgamma= -DBL_MAX;

  std::vector<int> demand_index(numD);
  if(permutation.size()==numD) {
    demand_index = permutation;
  } else {
    for(int i=0; i<numD; i++) {
      demand_index[i]=i;
    }
  }

  std::vector<int> bestpaths;
  bool gammanormal=true;
  bool gammaindividual=false;//true;
  bool gammajustgamma=false;//true;
  // set to iteratively sort order from remaining demands in each step
  // must be set to false unless using permutation "less" (option -2)
  // suggest always leave false here and set to true below if needed
  bool sort_order=false;
  std::vector<int> countunderv(1);

  int countnormal=0;
  int countindividual=0;
  int countjustgamma=0;

  countunderv[0]=0;

  std::vector< std::pair<int, double> > costpair(numD);

  if(permutation[0] == -2) {
    // sort demands on lowest cost path first
    // set this true here if it is needed
    sort_order=false;
    std::vector< std::pair<int,double> >::iterator pit,pend;
    pend=costpair.end();
    int i=0;
    for(pit=costpair.begin(); pit != costpair.end();pit++, i++) {
      (*pit).first=i;
      std::pair<long, double> p = findshortestpathcost(paths[i],vlengths);
      (*pit).second = p.second;
    }
	  
    std::sort(costpair.begin(),costpair.end(),less_cost());
    for(int m=0;m<numD;m++) {
      demand_index[m]=costpair[m].first;
    }
  }
	
  if(permutation[0] == -1) {
    random_shuffle(demand_index.begin(),demand_index.end());
  }
	
  if(delta == 0) {
    Rprintf("Error delta=0\n");
    // fix this just to "bomb out"
    D=1.0;
  }

  std::vector<int> pathdiffcount(1,0);
  std::vector<int> phasepathdiffcount(1,0);
	
  bool first_phase=true;
  std::vector<int> prevphasepathrecord(numD);

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
      SET_STRING_ELT(cmdsxp,0,Rf_mkChar(cmd));
      cmdexpr = PROTECT(R_ParseVector(cmdsxp, -1, &status, R_NilValue));
      for(int i = 0; i < Rf_length(cmdexpr); i++)
	ansxp = Rf_eval(VECTOR_ELT(cmdexpr, i), env);
      //Rprintf("D=%lg\n",D);
      UNPROTECT(1);
    }
	  
    bool underdemand=true;

    // note fix for random rotation of demand order

	  
    std::vector<int> pathrecord(numD);
    std::vector<int> prevpathrecord(numD);
    std::vector<double> weights(M,0.0);
    std::vector<int> recdemands(0);


    if(permutation[0] == -2) {
      // need to do this as it may be resized later if sort_order is true
      costpair.resize(numD);
      // sort demands on lowest cost path first
      std::vector< std::pair<int,double> >::iterator pit,pend;
      pend=costpair.end();
      int i=0;
      for(pit=costpair.begin(); pit != costpair.end();i++,pit++) {
	//(*pit).first=i;
	i=(*pit).first;
	std::pair<long, double> p = findshortestpathcost(paths[i],vlengths);
	(*pit).second = p.second;
      }
		
      std::sort(costpair.begin(),costpair.end(),less_cost());
      // inserted for testing - need to remove
      random_shuffle(costpair.begin(),costpair.end());
      for(int m=0;m<numD;m++) {
	demand_index[m]=costpair[m].first;
      }
    }
	  
    if(permutation[0] == -1) {
      random_shuffle(demand_index.begin(),demand_index.end());
    }

    for(int j=0; j<numD;j++,i++) {
		
      if(sort_order) {
	std::vector< std::pair<int,double> >::iterator pit,pend;
	pend=costpair.end();
	for(pit=costpair.begin(); pit != pend;pit++) {
	  long it=(*pit).first;
	  std::pair<long, double> p = findshortestpathcost(paths[it],vlengths);
	  (*pit).second = p.second;
	}
		  
	std::sort(costpair.begin(),costpair.end(),less_cost());
	i = costpair[0].first;
	costpair.erase(costpair.begin());
      } else {
	i= demand_index[j];
      }
      recdemands.push_back(i);
      double demand = vdemands[i];
      if(D >= 1.0)
	underdemand=false;

      int countunder=0;
      bool first_time=true;
      while( D < 1.0 && demand > 0.0) {
	// could use this
	/*std::pair<long, double> ppair = 
	  findshortestpathcostopt(paths[i],
	  vlengths,
	  vcapacity,
	  weights,
	  vdemands[i]);*/
	std::pair<long, double> ppair = 
	  findshortestpathcost(paths[i],
			       vlengths);
	/*
	  Rprintf("%ld, %lg, %ld, %lg, %lg\n",
	  tmppair.first,
	  tmppair.second,
	  ppair.first,
	  ppair.second,tmppair.second - ppair.second);*/

	int p = ppair.first;
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
		  
	// this needs more thought - this only records
	// on the last time
	if(!first_time) {
	  prevpathrecord[i]=pathrecord[i];
	} else {
	  first_time=false;
	  prevpathrecord[i]=p;
	}
		  
	pathrecord[i]=p;


	if(demand <=0) {
	  // while this records only the first time
	  // which one is best?
	  //if(first_time) {
	  for(vi=paths[i][p].begin(); vi != ve; vi++) {
	    weights[*vi] += origdemands[i];
	  }
	} 
	int sz = paths[i][p].size();
	for(int k=0; k < sz; k++) {
	  double length = vlengths[paths[i][p][k]];
	  length *= (1.0 + (eInternal * mincap) /
		     vcapacity[paths[i][p][k]]);
	  vlengths[paths[i][p][k]] = length;
	}

	// New
	if(demand > 0.0) {
	  countunder++;
	} 
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

	  

    phases++;
    totalphases++;

    int numdiff=0;

    //compare paths
    for(int j=0; j<numD; j++) {
      if(prevpathrecord[j] != pathrecord[j]) {
	numdiff++;
      }
    }
    if(pathdiffcount.size()<numdiff+1)
      pathdiffcount.resize(numdiff+1);
    pathdiffcount[numdiff]++;

    numdiff=0;
    if(first_phase) {
      for(int j=0; j<numD; j++) {
	prevphasepathrecord[j]=pathrecord[j];
      }
      first_phase=false;
    } else {
      //Rprintf("pathrecord\n");
		
      for(int j=0; j<numD; j++) {
	//Rprintf("%ld,",pathrecord[j]);
	if(prevpathrecord[j] != prevphasepathrecord[j]) {
	  numdiff++;
	}
      }
      //Rprintf("\n");
      if(phasepathdiffcount.size()<numdiff+1)
	phasepathdiffcount.resize(numdiff+1);
      phasepathdiffcount[numdiff]++;
		
    }

    betavals.push_back(calcbeta_restricted(paths,vlengths,vdemands,vcapacity));
    lambdavals.push_back(calcbeta_int(paths,vlengths,vdemands,vcapacity,pathrecord));
    std::vector<double> gweights;
    bool using_normal_gamma=true;
    bool using_prev_path_gamma=false;
    bool using_new_path_gamma=false;
    double gamma;
    if(using_normal_gamma) {
      gweights = 
	calcgraphweights(paths,pathrecord,origdemands,M);
		
      gamma = calcgamma(gweights,vcapacity);
      //gammavals.push_back(gamma);
		
      if(gamma > bestgamma) {
	bestgamma = gamma;
	bestpaths = pathrecord;
      }
      if(gamma >= intgamma * 0.99999)
	countgamma++;
      gammavals.push_back(gamma);
    }

    if(using_prev_path_gamma) {
      gweights = 
	calcgraphweights(paths,prevpathrecord,origdemands,M);
		
      gamma = calcgamma(gweights,vcapacity);
		
      if(gamma > bestgamma) {
	bestgamma = gamma;
	bestpaths = pathrecord;
      }
      if(gamma >= intgamma * 0.99999)
	countgamma++;
      gammavals.push_back(gamma);
    }

    if(using_new_path_gamma) {
      std::vector<int> newpath(numD);
      for(int j=0;j<numD;j++) {
	if(drand48() > 0.85)
	  newpath[j]=pathrecord[j];
	else
	  //newpath[j]=prevphasepathrecord[j];
	  newpath[j]=bestpaths[j];
	prevphasepathrecord[j]=pathrecord[j];
      }
      gweights = 
	calcgraphweights(paths,newpath,origdemands,M);
      /*Rprintf("newpath\n");
		  
	for(int j=0; j<numD; j++) {
	Rprintf("%ld,",newpath[j]);
	}
	Rprintf("\n");*/
		
      gamma = calcgamma(gweights,vcapacity);
		
      if(gamma > bestgamma) {
	bestgamma = gamma;
	bestpaths = newpath;
      }
		
		
      if(gamma >= intgamma * 0.99999)
	countgamma++;
      gammavals.push_back(gamma);
    }
    //Rprintf("gamma %lg, matches int gamma of %lg\n",gamma, intgamma);
  }

  if(progress != false)
    Rprintf("\n");

  if(intgamma < DBL_MAX) {
    for(i=0;i<countunderv.size();i++) {
      Rprintf("%ld,",countunderv[i]);
    }
    Rprintf("; countgamma=%ld\n",
	    countgamma);
    Rprintf("******** bestgamma=%lg, intgamma=%lg, ratio=%lg\n",
	    bestgamma,intgamma,bestgamma/intgamma);
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
  retlist.push_back(pathdiffcount,"pathdiffcount");
  retlist.push_back(phasepathdiffcount,"phasepathdiffcount");
  retlist.push_back(gammavals,"gammavals");
  retlist.push_back(betavals,"betavals");
  retlist.push_back(lambdavals,"lambdavals");
  return Rcpp::wrap(retlist);
}


SEXP rg_test_every_path_inner
(
 SEXP RRdemandpaths,
 SEXP Rdemands,
 SEXP Rcapacity,
 SEXP Rprogress,
 SEXP pb,
 SEXP Rrecordlen,
 SEXP env
 ) {
	
  int recordlen = Rcpp::as<double>(Rrecordlen);
  bool progress=Rcpp::as<bool>(Rprogress);
  std::vector<double> vcapacity = Rcpp::as< std::vector<double> > (Rcapacity);
  std::vector<double> vdemands = Rcpp::as< std::vector<double> > (Rdemands);
	
  int L = vdemands.size();
  std::vector<double>::iterator vecdit;
  int M = vcapacity.size();
  //decode Rdemandpaths into paths
  std::vector<int>::iterator vi,ve;
  std::vector<double>::iterator vw,vc;
  std::vector< std::vector< std::vector<int> > > paths;
  std::vector< std::vector< std::vector<int> > >::iterator dit;
  std::vector< std::vector<int> >::iterator pit;
  std::vector<int>::iterator eit;

  std::vector< std::pair<double, std::vector<int> > > recorddemands(recordlen);
  std::vector< std::pair<double, std::vector<int> > >::iterator rit;
	
  std::vector<int> nullpath(0);
  for(rit = recorddemands.begin(); rit != recorddemands.end(); rit++) {
    *rit = std::pair<double, std::vector<int> >((double)-DBL_MAX, nullpath);
  }
	
  make_heap(recorddemands.begin(),recorddemands.end(),less_demand());

  Rcpp::List demands(RRdemandpaths);
  int sz = demands.size();
  long pathcomb=1;
  paths.resize(sz);
  for(int i=0;i<sz;i++) {
    Rcpp::List demandpaths((SEXPREC*)demands[i]);
    int sz = demandpaths.size();
    paths[i].resize(sz);
    pathcomb *= paths[i].size();
    for(int j=0;j<sz;j++) {
      paths[i][j]= Rcpp::as< std::vector<int> >(demandpaths[j]);
      //paths[i][j]= RcppVector<int>((SEXPREC*)demandpaths[j]).stlVector();
    }
  }
	
  std::vector< std::vector<double> > pathflows(paths.size());
  for(int i=0;i<L;i++) {
    pathflows[i] = std::vector<double>(paths[i].size(),0.0);
  }

  bool finished=false;
	
  //comput L

  // check this starts with all 0
  std::vector<int> pathptr(L,0);
  std::vector<int> bestpathptr;
  long count=0;
  double bestgamma = -DBL_MAX;

  // needed for progress bar
  char cmd[256];
  SEXP cmdsxp, cmdexpr, ansxp = R_NilValue;
  ParseStatus status;
  PROTECT(cmdsxp = Rf_allocVector(STRSXP,1));
  long updatepb = pathcomb / 100.0;

  while(!finished) {

    if(progress != false && count % updatepb == 0) {
      sprintf(cmd,"setTxtProgressBar(pb,%ld)",count);
      SET_STRING_ELT(cmdsxp,0,Rf_mkChar(cmd));
      cmdexpr = PROTECT(R_ParseVector(cmdsxp, -1, &status, R_NilValue));
      for(int i = 0; i < Rf_length(cmdexpr); i++)
	ansxp = Rf_eval(VECTOR_ELT(cmdexpr, i), env);
      UNPROTECT(1);
    }
	  


    std::vector<double> weights(M,0);
	
    for(int i=0; i<L; i++){
      //Rprintf("pathptr %d\n",pathptr[i]);
      std::vector<int> path=paths[i][pathptr[i]];
      int j;
      for(j=0, vi=path.begin();vi != path.end(); j++,vi++) {
	weights[*vi] += vdemands[i];
      }
    }
    //Rprintf("\n");

    double gamma = (double)DBL_MAX;

    for(vw=weights.begin(), vc=vcapacity.begin();
	vw != weights.end(); 
	vw++, vc++) {
      double tmp = (*vc - *vw)/ (*vc);
      if (tmp < gamma) 
	gamma = tmp;
    }
	  
    if(gamma > bestgamma) {
      bestgamma = gamma;
      //Rprintf("best gamma so far=%g\n",bestgamma);
    }
    if(gamma > recorddemands.front().first) {
      pop_heap(recorddemands.begin(),recorddemands.end(),less_demand());
      recorddemands.pop_back();
      recorddemands.push_back(std::pair<double, std::vector<int> >
			      (gamma,pathptr));
      push_heap(recorddemands.begin(),recorddemands.end(),less_demand());

		
    }
    bool increment=true;
    bool incrementnext=false;
    for(int i=0; i<L ; i++) {
      if(paths[i].size() > 1) {
	if( (pathptr[i] + 1) % paths[i].size() == 0)
	  incrementnext = true;
	else 
	  incrementnext = false;
	if(increment)
	  pathptr[i] = (pathptr[i] + 1) % paths[i].size();
	if(incrementnext && increment)
	  increment=true;
	else
	  increment=false;
      }

    }
    if(increment )
      finished=true;
    count++;
    // just for testing
    //if(count > 0)
    //finished=true;
		
		
  }

  sort_heap(recorddemands.begin(),recorddemands.end(),less_demand());
  std::vector< std::vector<int> >  pathptrs(recordlen);
  std::vector< double > gammas(recordlen);
  for(int i=0;i<recordlen;i++) {
    //Rprintf("%lg, ",recorddemands[i].first);
    pathptrs[i]=recorddemands[i].second;
    gammas[i]=recorddemands[i].first;
  }

  UNPROTECT(1);
  //Rprintf(" end \n");
  //Rprintf("best gamma %lg\n",bestgamma);
  Rcpp::List retlist;
  int var=0;
  retlist.push_back(pathptrs,"pathptrs");
  retlist.push_back(gammas,"gammas");
  return Rcpp::wrap(retlist);
	
}



SEXP rg_fleischer_max_concurrent_flow_c_igraph(SEXP num_verts_in,
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
					       SEXP Rpermutation,
					       SEXP Rdeltaf
					       ) {

  bool progress=false;
  bool updateflow=Rcpp::as<bool>(Rupdateflow);
  std::string title("");
  if ( Rf_isLogical(Rprogress) ) {
    progress=Rcpp::as<bool>(Rprogress);
  } else if ( Rf_isString(Rprogress) ) {
    title = std::string(CHAR(STRING_ELT(Rprogress,0)));
    progress=true;
  }
  std::vector<long> vedges = Rcpp::as< std::vector<long> > (R_edges_in);
  int m = vedges.size()/2;
  std::vector<rg_demandi>::iterator vi,ve;
    

    
    
  std::vector<double> capacities = Rcpp::as< std::vector<double> > (capacities_in);
  int N = Rcpp::as<int>(num_verts_in);
  double* dem_in = REAL(demands_in);
  int* dem_sources_in = INTEGER(demands_sources_in);
  int* dem_sinks_in = INTEGER(demands_sinks_in);
  int num_dem = Rcpp::as<int>(num_demands_in);
  double deltaf = Rcpp::as<double>(Rdeltaf);

  std::vector<rg_demandi> demandsi(num_dem);
  for(int i=0 ; i<num_dem; i++) {
    demandsi[i].demand = dem_in[i];
    demandsi[i].flow = 0;
    demandsi[i].source = dem_sources_in[i];
    demandsi[i].sink = dem_sinks_in[i];
  }

  double e = Rcpp::as<double>(Re);
  long totalphases;
  std::vector<double> lengths;
  
  igraph_rg_fleischer_max_concurrent_flow(vedges, capacities, demandsi,
					  lengths,e,N,totalphases);
  
  Rcpp::List retdemandsi;
  for(int i=0 ; i<num_dem; i++) {
    Rcpp::List demand;
    demand.push_back(demandsi[i].demand,"demand");
    demand.push_back(demandsi[i].flow,"flow");
    demand.push_back(boost::lexical_cast<std::string>(demandsi[i].source+1),
			     "source");
    demand.push_back(boost::lexical_cast<std::string>(demandsi[i].sink+1)
			     ,"sink");
    Rcpp::List paths;
    std::map< std::vector< long > ,double >::iterator mit;
    for(mit=demandsi[i].path_flow_map.begin();
	mit != demandsi[i].path_flow_map.end();
	mit++) {
      std::vector<long> path = mit->first;
      std::string spath;
      std::vector<long>::iterator k=path.begin();
      spath = boost::lexical_cast<std::string>( *k + 1);
      k++;
      for(;
	  k != path.end();
	  k++) {
	spath += "|";
	spath += boost::lexical_cast<std::string>( *k + 1);
      }
      paths.push_back(mit->second,spath.c_str());
    }
    demand.push_back(paths,"paths");
    retdemandsi.push_back(demand,boost::lexical_cast<std::string> (i+1));
  }
  
  Rcpp::List retlist;
  retlist.push_back(retdemandsi,"demands");
  retlist.push_back(totalphases,"totalphases");
  retlist.push_back(lengths,"lengths");
  
  return wrap(retlist);
}

__END_DECLS
