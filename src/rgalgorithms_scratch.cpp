#include "rgalgorithms.h"



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
  
extern "C" {


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
	std::vector<double> vcapacity = RcppVector<double>(Rcapacity).stlVector();
	std::vector<double> vdemands = RcppVector<double>(Rdemands).stlVector();
	
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
	PROTECT(cmdsxp = allocVector(STRSXP,1));
	long updatepb = pathcomb / 100.0;

	while(!finished) {

	  if(progress != false && count % updatepb == 0) {
		sprintf(cmd,"setTxtProgressBar(pb,%ld)",count);
		SET_STRING_ELT(cmdsxp,0,mkChar(cmd));
		cmdexpr = PROTECT(R_ParseVector(cmdsxp, -1, &status, R_NilValue));
		for(int i = 0; i < length(cmdexpr); i++)
		  ansxp = eval(VECTOR_ELT(cmdexpr, i), env);
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

}
