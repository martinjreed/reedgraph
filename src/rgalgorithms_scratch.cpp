#include "rgalgorithms.h"



  
extern "C" {



  SEXP rg_fleischer_max_concurrent_flow_with_int_c(SEXP num_verts_in,
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
	bool updateflow=asLogical(Rupdateflow);
	std::string title("");
	if ( isLogical(Rprogress) ) {
	  progress=asLogical(Rprogress);
	} else if ( isString(Rprogress) ) {
	  title = std::string(CHAR(STRING_ELT(Rprogress,0)));
	  progress=true;
	}
	Graph_rg gdual(num_verts_in, num_edges_in, R_edges_in, R_weights_in,capacities_in);

	std::vector<int> permutation = Rcpp::as< std::vector<int> > (Rpermutation);
	double* dem_in = REAL(demands_in);
	int* dem_sources_in = INTEGER(demands_sources_in);
	int* dem_sinks_in = INTEGER(demands_sinks_in);
	int num_dem = asInteger(num_demands_in);
	double deltaf = asReal(Rdeltaf);
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

	PROTECT(cmdsxp = allocVector(STRSXP,1));
	
	char cmd[256];
	//ansxp = eval(cmdexpr,env);

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
		SET_STRING_ELT(cmdsxp,0,mkChar(cmd));
		cmdexpr = PROTECT(R_ParseVector(cmdsxp, -1, &status, R_NilValue));
		for(int i = 0; i < length(cmdexpr); i++)
		  ansxp = eval(VECTOR_ELT(cmdexpr, i), env);
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

}
