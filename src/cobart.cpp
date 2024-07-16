#include <ctime>
#include "common.h"
#include "tree.h"
#include "treefuns.h"
#include "info.h"
#include "bartfuns.h"
#include "bd.h"
#include "heterbart.h"

double rtrunnormGen(double mu, double  sigma, double a, double b){
  
  double psi_lower = R::pnorm((a-mu)/sigma,0,1,1,0);
  double psi_higher = R::pnorm((b-mu)/sigma,0,1,1,0);
  double out = mu+sigma*R::qnorm(R::runif(0,1)*(psi_higher-psi_lower)+psi_lower,0,1,1,0);
  return(out);
}

double extractZByCategory(Rcpp::IntegerVector y, Rcpp::NumericVector z, int category, double gamma,bool max) {
  Rcpp::NumericVector result;
  
  // Loop through y and z to extract values for the specified category
  for (int i = 0; i < y.size(); ++i) {
    if (y[i] == category) {
      result.push_back(z[i]);
    }
  }
  result.push_back(gamma);
  double extreme = result[0];
  for (int i = 0; i < result.size(); ++i ){
    if(max && result[i] > extreme){ extreme = result[i]; }
    else if (!max && result[i] < extreme){ extreme = result[i]; }
  } 
  
  return extreme;
}

#ifndef NoRcpp

#define TRDRAW(a, b) trdraw(a, b)
#define TEDRAW(a, b) tedraw(a, b)
#define TRGAMMADRAW(a, b) trgammadraw(a, b)

RcppExport SEXP cobart(
    SEXP cintr,            //number of observations in training data
    SEXP cip,            //dimension of x
    SEXP cinp,           //number of observations in test data
    SEXP cix,            //x, train,  pxn (transposed so rows are contiguous in memory)
    SEXP ciy,            //y, train,  nx1
    SEXP cixp,           //x, test, pxnp (transposed so rows are contiguous in memory)
    SEXP cim,            //number of trees
    SEXP cinc,           //number of cut points
    SEXP cind,           //number of kept draws (except for thinnning ..)
    SEXP ciburn,         //number of burn-in draws skipped
    SEXP cipower,
    SEXP cibase,
    SEXP citau,
    SEXP cioldgamma,  // initial value of cutoff values, obtained by proportional odds logistic regression
    SEXP cioldz, // initial value of latent continuous values, obtained by proportional odds logistic regression
    SEXP ciXinfo,
    SEXP cioffset
    //SEXP cinprintevery
)
{
  
  //--------------------------------------------------
  //process args
  size_t n = Rcpp::as<int>(cintr);
  size_t p = Rcpp::as<int>(cip);
  size_t np = Rcpp::as<int>(cinp);
  Rcpp::NumericVector  xv(cix);
  double *ix = &xv[0];
  Rcpp::IntegerVector  yv(ciy); 
  //int *iy = &yv[0];
  Rcpp::NumericVector  xpv(cixp);
  double *ixp = &xpv[0];
  size_t m = Rcpp::as<int>(cim);
  //size_t nc = Rcpp::as<int>(_inc);
  Rcpp::IntegerVector _nc(cinc);
  int *numcut = &_nc[0];
  size_t nd = Rcpp::as<int>(cind);
  size_t burn = Rcpp::as<int>(ciburn);
  double mybeta = Rcpp::as<double>(cipower);
  double alpha = Rcpp::as<double>(cibase);
  double tau = Rcpp::as<double>(citau);
  double offset = Rcpp::as<double>(cioffset);
  Rcpp::NumericVector  coldgamma(cioldgamma);
  Rcpp::NumericVector  coldz(cioldz);
  size_t ngamma = coldgamma.size();
  
  Rcpp::NumericVector  oldz(n);
  Rcpp::NumericVector  oldgamma(ngamma);
  
  for(size_t i=0;i<n;i++) {oldz[i] = coldz[i];}
  for(size_t i=0;i<ngamma;i++) {oldgamma[i] = coldgamma[i];}
  
  double *iz = &oldz[0];
  Rcpp::NumericMatrix Xinfo(ciXinfo);
  Rcpp::NumericMatrix trdraw(nd,n);
  Rcpp::NumericMatrix tedraw(nd,np);

  Rcpp::NumericMatrix trgammadraw(nd,ngamma);
  
  //size_t printevery = Rcpp::as<int>(cinprintevery);
  //random number generation
  arn gen;
  
  heterbart bm(m);
  
  if(Xinfo.size()>0) {
    xinfo _xi;
    _xi.resize(p);
    for(size_t i=0;i<p;i++) {
      _xi[i].resize(numcut[i]);
      //Rcpp::IntegerVector cutpts(Xinfo[i]);
      for(size_t j=0;j<(size_t)numcut[i];j++) _xi[i][j]=Xinfo(i, j);
    }
    bm.setxinfo(_xi);
  }
#else
  
#define TRDRAW(a, b) trdraw[a][b]
#define TEDRAW(a, b) tedraw[a][b]
#define TRGAMMADRAW(a, b) trgammadraw[a][b]
  
  void cobart(
      size_t n,            //number of observations in training data
      size_t p,		//dimension of x
      size_t np,		//number of observations in test data
      double* ix,		//x, train,  pxn (transposed so rows are contiguous in memory)
      int* iy,		//y, train,  nx1
      double* ixp,		//x, test, pxnp (transposed so rows are contiguous in memory)
      size_t m,		//number of trees
      int *numcut,		//number of cut points
      size_t nd,		//number of kept draws (except for thinnning ..)
      size_t burn,		//number of burn-in draws skipped
      double mybeta,
      double alpha,
      unsigned int n1, // additional parameters needed to call from C++
      unsigned int n2,
      /*
       double* trmean,
       double* temean,
       */
      double tau,
      double offset,
      size_t ngamma = oldgamma.size();
  size_t printevery = Rcpp::as<int>(_inprintevery);
  
      Rcpp::NumericVector  coldgamma,
      Rcpp::NumericVector  coldz,
      Rcpp::NumericVector  oldz(n);
  Rcpp::NumericVector  oldgamma(ngamma);
  
  for(size_t i=0;i<n;i++) {oldz[i] = coldz[i];}
  for(size_t i=2;i<ngamma;i++) {oldgamma[i] = coldgamma[i];}
  
      double* iz,
      
  double* _trdraw,
  double* _tedraw,
  double* _trgammadraw
  )
  {
    
    //return data structures (using C++)
    std::vector<double*> trdraw(nd);
    std::vector<double*> tedraw(nd);
    std::vector<double*> trgammadraw(nd);
    
    for(size_t i=0; i<nd; ++i) trdraw[i]=&_trdraw[i*n];
    for(size_t i=0; i<nd; ++i) tedraw[i]=&_tedraw[i*np];
    for(size_t i=0; i<nd; ++i) trgammadraw[i]=&_trgammadraw[i*ngamma];
    
    //random number generation
    arn gen(n1, n2);
    
    heterbart bm(m);
#endif
    
    std::stringstream treess;  //string stream to write trees to
    treess.precision(10);
    treess << nd << " " << m << " " << p << endl;
    
    //printf("*****Into main of lbart\n");
    
    //--------------------------------------------------
    //print args
    //printf("*****Data:\n");
    //printf("data:n,p,np: %zu, %zu, %zu\n",n,p,np);
    //printf("y1,yn: %d, %d\n",iy[0],iy[n-1]);
    //printf("x1,x[n*p]: %lf, %lf\n",ix[0],ix[n*p-1]);
    //if(np) printf("xp1,xp[np*p]: %lf, %lf\n",ixp[0],ixp[np*p-1]);
    //printf("*****Number of Trees: %zu\n",m);
    //printf("*****Number of Cut Points: %d ... %d\n", numcut[0], numcut[p-1]);
    //printf("*****burn and ndpost: %zu, %zu\n",burn,nd);
    //   printf("*****Prior:\nbeta,alpha,tau,nu,lambda: %lf,%lf,%lf,%lf,%lf\n",
    //printf("*****Prior:\nbeta,alpha,tau: %lf,%lf,%lf\n",
    //       mybeta,alpha,tau);
    
    //--------------------------------------------------
    
    double *svec = new double[n];
    for(size_t i=0;i<n;i++) svec[i]=1;
    
    bm.setprior(alpha,mybeta,tau);
    bm.setdata(p,n,ix,iz,numcut);
    
    //--------------------------------------------------
    //mcmc
    //printf("\nMCMC\n");
    //size_t index;
    size_t trcnt=0; //count kept train draws
    size_t tecnt=0; //count kept train draws
    
    //time_t tp;
    //int time1 = time(&tp);
    xinfo& xi = bm.getxinfo();
    //for(size_t i=0;i<n;i++) trmean[i]=0.0;
    //for(size_t i=0;i<np;i++) temean[i]=0.0;
    
    double* fhattest=0; 
    if(np) { fhattest = new double[np]; }
    
    Rcpp::NumericVector cutoff1(n);
    Rcpp::NumericVector cutoff2(n);
    
    Rcpp::NumericVector oldg(ngamma);
    for(size_t k=0; k<ngamma; k++) {
      oldg[k] = oldgamma[k];
    }
    
    double acnum = 0;
    
    for(size_t i=0;i<(nd+burn);i++) {
      
      bm.draw(svec,gen);
      
      //if(i%printevery==0) printf("done %zu (out of %zu)\n",i,nd+burn);
      //if(i%printevery==0) printf("oldgamma %f \n", oldgamma[2]);
      
      double sigmaH = 0.5/(ngamma-1);
      
      // draw gamma
      for(size_t k=2; k < ngamma-1; k++) {
        oldg[k] = rtrunnormGen(oldgamma[k], sigmaH, oldg[k-1], oldgamma[k+1]);
      }
      //printf("oldg %lf\n",oldg[0]);
      //printf("oldg %lf\n",oldg[1]);
      //printf("oldg %lf\n",oldg[2]);
      //printf("oldg %lf\n",oldg[3]);
      // compute the acceptance ratio R
      double product1 = 1.0;
      double product2 = 1.0;
      
      for(size_t k=0; k < n; k++) {
        double temp1 = R::pnorm(oldg[yv[k]]-bm.f(k)-offset, 0, 1, true, false) - R::pnorm(oldg[yv[k]-1]-bm.f(k)-offset, 0, 1, true, false);  
        double temp2 = R::pnorm(oldgamma[yv[k]]-bm.f(k)-offset, 0, 1, true, false) - R::pnorm(oldgamma[yv[k]-1]-bm.f(k)-offset, 0, 1, true, false);  
        //printf("temp2 %lf\n",temp2);
        product1 *= temp1/temp2; 
      }
      //printf("product1 %lf\n",product1);
      
      for(size_t k=2; k < ngamma-1; k++) {
        double temp3 = R::pnorm((oldgamma[k+1]-oldgamma[k])/sigmaH, 0, 1, true, false) - R::pnorm((oldg[k-1]-oldgamma[k])/sigmaH, 0, 1, true, false);  
        double temp4 = R::pnorm((oldg[k+1]-oldg[k])/sigmaH, 0, 1, true, false) - R::pnorm((oldgamma[k-1]-oldg[k])/sigmaH, 0, 1, true, false);   
        product2 *= temp3/temp4; 
      }
      //printf("product2 %lf\n",product2);
      double acr = product2 * product1;
      //printf("acr %lf\n",acr);
      double uniformRandom = R::runif(0, 1);
      
      if(uniformRandom <= acr ) {
        acnum += 1;
        for(size_t k=2; k < ngamma-1; k++) {
          oldgamma[k] = oldg[k];
        }
        //printf("accept!!! \n");
        //printf("uniform and acr: %lf, %lf\n", uniformRandom, acr);
        //printf("accept number: %lf\n", acnum);
      }
      
      //double acratio = acnum/(i+1);
      //if(i%printevery==0 || i == (burn + nd - 1)) {
        //printf("accept ratio: %lf\n", acratio);
        //printf("oldgamma: %lf\n", oldgamma[0]);
        //printf("oldgamma: %lf\n", oldgamma[1]);
        //printf("oldgamma: %lf\n", oldgamma[2]);
        //printf("oldgamma: %lf\n", oldgamma[3]);
        //printf("oldz: %lf\n", oldgamma[0]);
        //printf("oldz: %lf\n", oldgamma[1]);
        //printf("oldz: %lf\n", oldgamma[2]);
        //printf("oldz: %lf\n", oldgamma[3]);
      //}
      
      // draw z
      for (size_t k=0; k<n; k++){
        cutoff1[k] = oldgamma[yv[k]-1];
        //printf("cutoff1 %f \n",cutoff1[k]);
        cutoff2[k] = oldgamma[yv[k]];
        //printf("cutoff2 %f \n",cutoff2[k]);
      }
      
      for(size_t k=0; k<n; k++) {
        oldz[k] = rtrunnormGen(bm.f(k)+offset, 1, cutoff1[k], cutoff2[k]) - offset;
        //trmean[k]+=bm.f[k];
      }
    
      if(i>=burn) {
        for(size_t k=0;k < n; k++) TRDRAW(trcnt,k)=bm.f(k);
        for(size_t k=0;k < ngamma;k++) TRGAMMADRAW(trcnt,k)=oldgamma[k];
        trcnt+=1;
        
        if(np) {
          bm.predict(p,np,ixp,fhattest);
          //index=tecnt*np;
          //for(size_t k=0;k<np;k++) tedraw[index+k]=fhattest[k];
          for(size_t k=0;k<np;k++) TEDRAW(tecnt,k)=fhattest[k];
          tecnt+=1;
        }
      }
      
      for(size_t j=0;j<m;j++) {
        treess << bm.gettree(j);
      }
      
    }
    
    //int time2 = time(&tp);
    //printf("time: %ds\n",time2-time1);
    //   for(size_t k=0;k<n;k++) trmean[k]/=nd;
    //   for(size_t k=0;k<np;k++) temean[k]/=temecnt;
    //printf("check counts\n");
    //printf("trcnt,tecnt: %zu,%zu\n",trcnt,tecnt);
    //printf("trcnt,tecnt,temecnt,treedrawscnt: %zu,%zu,%zu,%zu\n",trcnt,tecnt,temecnt,treedrawscnt);
    //--------------------------------------------------
    
    if(fhattest) delete[] fhattest;
    delete[] svec;
    
#ifndef NoRcpp
    //--------------------------------------------------
    //return
    Rcpp::List ret;
    //   ret["yhat.train.mean"]=trmean;
    ret["yhat.train"]=trdraw;
    //ret["yhat.test.mean"]=temean;
    ret["yhat.test"]=tedraw;
    ret["gammahat"]=trgammadraw;
    
    Rcpp::List xiret(xi.size());
    for(size_t i=0;i<xi.size();i++) {
      Rcpp::NumericVector vtemp(xi[i].size());
      std::copy(xi[i].begin(),xi[i].end(),vtemp.begin());
      xiret[i] = Rcpp::NumericVector(vtemp);
    }
    
    Rcpp::List treesL;
    treesL["nkeeptreedraws"] = Rcpp::wrap<int>(nd); //in trees
    treesL["ntree"] = Rcpp::wrap<int>(m); //in trees
    treesL["numx"] = Rcpp::wrap<int>(p); //in cutpoints
    treesL["cutpoints"] = xiret;
    treesL["trees"]=Rcpp::CharacterVector(treess.str());
       //if(treesaslists) treesL["lists"]=list_of_lists;
    ret["treedraws"] = treesL;
    
    return ret;
#else
    
#endif
    
  }
  