
#include <Rcpp.h>

using namespace Rcpp;

// [[Rcpp::export]]
double sample_m(double n) {
     IntegerVector pool = Rcpp::seq_len(n+1);
     std::random_shuffle(pool.begin(), pool.end());
     double samp = pool[0]-1.0;
     return samp;
}




// [[Rcpp::export]]
Rcpp::NumericVector reject(double n) {
    double temp = 0.0;
    double mm = 0.0;
    while(temp<1.0){
      mm = sample_m(n);
      double pm0 = mm/n;
      int lg = 1;
      double dmm = R::dbinom(mm,n,pm0,lg);
      Rcpp::NumericVector ru  = Rcpp::runif(1,0.0,1.0);
      double lru  = std::log(ru[0]);
      if(dmm <= lru){
        temp = temp + 1.0;
      }
    }
      Rcpp::NumericVector mmm(2);
      mmm[0] = mm;
      mmm[1] = n - mm;  
    return mmm;
  }


// [[Rcpp::export]]
Rcpp::NumericVector gibbs(NumericVector n, NumericVector ks){
  int ioff = 0;
  int counti=0;
  int kse = ks.size();
  for(int tt=0; tt < kse; tt++){
    counti = ks[tt]-1;
    for(int ii=0; ii < counti;ii++){
      int beginn = ioff+ii;
      int endd = ioff+counti;
      double N = n[beginn] + n[endd];
      if(N > 0){
        Rcpp::NumericVector ww = reject(N);
        n[beginn] = ww[0];
        n[endd] = ww[1];
      }
    }
    ioff = ioff + ks[tt];
  }
  return n;
}






// [[Rcpp::export]]
Rcpp::NumericVector burnin(int burn, NumericVector init, NumericVector ks, NumericVector Ns){
  Rcpp::NumericVector temp1 = init;
  Rcpp::NumericVector temp2;

  for(int zz=0; zz < burn; zz++){
      temp2 = gibbs(temp1,ks);
      temp1 = temp2;
  }
  return temp2;
}


// [[Rcpp::export]]
Rcpp::NumericMatrix gen_chain(NumericVector ks, NumericVector Ns, int batchsize, int thin, NumericVector startvec){
  Rcpp::NumericVector tmp1 = startvec;
  Rcpp::NumericVector tmp2 = startvec;  
  int KS = sum(ks);
  Rcpp::NumericMatrix   genat(batchsize, KS);
  for(int jj=0;jj < KS; jj++){
    genat(0,jj) = startvec[jj];
  }
  
  for(int www=1; www < batchsize; www++){
    for(int ddd=0; ddd < thin; ddd++){
      tmp2 = gibbs(tmp1,ks);
      tmp1 = gibbs(tmp2,ks);
    }
     for(int jj=0;jj < KS; jj++){
      genat(www,jj) = tmp2[jj];
    }
  }
  return genat;
}





// [[Rcpp::export]]
double SCR(NumericVector Q0, NumericVector chd, NumericVector NN){

  NumericVector Q = pnorm(Q0);
  
 double Vt = Q[0];
 double Gt = Q[1];
 double Vr = Q[2];
 double Gr = Q[3];
 double a  = Q[4];
 double b  = Q[5];
 double Vtx = Q[6];
 double Gtx = Q[7];
 double Vrx = Q[8];
 double Grx = Q[9];

  Rcpp::NumericVector ee(15);

	ee[0] = Vt + (1-Vt)*Gt*a + (1-Vt)*(1-Gt)*b*a;
	ee[1] = (1-Vt)*Gt*(1-a) + (1-Vt)*(1-Gt)*b*(1-a);
	ee[2] = (1-Vt)*(1-Gt)*(1-b);

	ee[3] = (1-Vr)*Gr*a + (1-Vr)*(1-Gr)*b*a;
	ee[4] = Vr + (1-Vr)*Gr*(1-a) + (1-Vr)*(1-Gr)*b*(1-a);
	ee[5] = (1-Vr)*(1-Gr)*(1-b);

	ee[6] = Vtx + (1-Vtx)*Gtx*a + (1-Vtx)*(1-Gtx)*b*a;
	ee[7] = (1-Vtx)*Gtx*(1-a) + (1-Vtx)*(1-Gtx)*b*(1-a);
	ee[8] = (1-Vtx)*(1-Gtx)*(1-b);

	ee[9] = (1-Vrx)*Grx*a + (1-Vrx)*(1-Grx)*b*a;
	ee[10] = Vrx + (1-Vrx)*Grx*(1-a) + (1-Vrx)*(1-Grx)*b*(1-a);
	ee[11] = (1-Vrx)*(1-Grx)*(1-b);

	ee[12] = b*a;
	ee[13] = b*(1-a);
	ee[14] = (1-b);
	
	ee = ee*NN;
	
  double LL = 0.0;
  
  for(int ii=0; ii<15; ii++){
    if(chd[ii] > 0.0){
      LL = LL + chd[ii]*(log(chd[ii])-log(ee[ii]));
    }
  }
  return LL;
}




// [[Rcpp::export]]
double LT(NumericVector Q0, NumericVector chd, NumericVector NN){

 NumericVector Q = pnorm(Q0);
  
 double Vt = Q[0];
 double Gt = Q[1];
 double Vr = Q[2];
 double Gr = Q[3];
 double a  = Q[4];
 double b  = Q[5];
 double Vtx = Q[6];
 double Gtx = Q[7];
 double Vrx = Q[8];
 double Grx = Q[9];

 Rcpp::NumericVector ee(15);

  ee[0] = Gt*Vt + Gt*(1-Vt)*a + (1-Gt)*b*a;
	ee[1] = Gt*(1-Vt)*(1-a) + (1-Gt)*b*(1-a);
	ee[2] = (1-Gt)*(1-b);

	ee[3] = Gr*(1-Vr)*a + (1-Gr)*b*a;
	ee[4] = Gr*Vr +Gr*(1-Vr)*(1-a) + (1-Gr)*b*(1-a);
	ee[5] = (1-Gr)*(1-b);

	ee[6] = Gtx*Vtx + Gtx*(1-Vtx)*a + (1-Gtx)*b*a;
	ee[7] = Gtx*(1-Vtx)*(1-a) + (1-Gtx)*b*(1-a);
	ee[8] = (1-Gtx)*(1-b);

	ee[9] = Grx*(1-Vrx)*a + (1-Grx)*b*a;
	ee[10] = Grx*Vrx +Grx*(1-Vrx)*(1-a) + (1-Grx)*b*(1-a);
	ee[11] = (1-Grx)*(1-b);

	ee[12] = b*a;
	ee[13] = b*(1-a);
	ee[14] = (1-b);
	
 	ee = ee*NN;
	
  double LL = 0.0;
  
  for(int ii=0; ii<15; ii++){
    if(chd[ii] > 0.0){
      LL = LL + chd[ii]*(log(chd[ii])-log(ee[ii]));
    }
  }
  return LL;
}







// [[Rcpp::export]]
double SCR_V(NumericVector Q0, NumericVector chd, NumericVector NN){

  NumericVector Q = pnorm(Q0);
  
 double Vt = Q[0];
 double Gt = Q[1];
 double Vr = Q[2];
 double Gr = Q[3];
 double a  = Q[4];
 double b  = Q[5];
 double Vtx = Q[0];
 double Gtx = Q[6];
 double Vrx = Q[2];
 double Grx = Q[7];

  Rcpp::NumericVector ee(15);

	ee[0] = Vt + (1-Vt)*Gt*a + (1-Vt)*(1-Gt)*b*a;
	ee[1] = (1-Vt)*Gt*(1-a) + (1-Vt)*(1-Gt)*b*(1-a);
	ee[2] = (1-Vt)*(1-Gt)*(1-b);

	ee[3] = (1-Vr)*Gr*a + (1-Vr)*(1-Gr)*b*a;
	ee[4] = Vr + (1-Vr)*Gr*(1-a) + (1-Vr)*(1-Gr)*b*(1-a);
	ee[5] = (1-Vr)*(1-Gr)*(1-b);

	ee[6] = Vtx + (1-Vtx)*Gtx*a + (1-Vtx)*(1-Gtx)*b*a;
	ee[7] = (1-Vtx)*Gtx*(1-a) + (1-Vtx)*(1-Gtx)*b*(1-a);
	ee[8] = (1-Vtx)*(1-Gtx)*(1-b);

	ee[9] = (1-Vrx)*Grx*a + (1-Vrx)*(1-Grx)*b*a;
	ee[10] = Vrx + (1-Vrx)*Grx*(1-a) + (1-Vrx)*(1-Grx)*b*(1-a);
	ee[11] = (1-Vrx)*(1-Grx)*(1-b);

	ee[12] = b*a;
	ee[13] = b*(1-a);
	ee[14] = (1-b);
	
	ee = ee*NN;
	
  double LL = 0.0;
  
  for(int ii=0; ii<15; ii++){
    if(chd[ii] > 0.0){
      LL = LL + chd[ii]*(log(chd[ii])-log(ee[ii]));
    }
  }
  return LL;
}




// [[Rcpp::export]]
double SCR_G(NumericVector Q0, NumericVector chd, NumericVector NN){

  NumericVector Q = pnorm(Q0);
  
 double Vt = Q[0];
 double Gt = Q[1];
 double Vr = Q[2];
 double Gr = Q[3];
 double a  = Q[4];
 double b  = Q[5];
 double Vtx = Q[6];
 double Gtx = Q[1];
 double Vrx = Q[7];
 double Grx = Q[2];

  Rcpp::NumericVector ee(15);

	ee[0] = Vt + (1-Vt)*Gt*a + (1-Vt)*(1-Gt)*b*a;
	ee[1] = (1-Vt)*Gt*(1-a) + (1-Vt)*(1-Gt)*b*(1-a);
	ee[2] = (1-Vt)*(1-Gt)*(1-b);

	ee[3] = (1-Vr)*Gr*a + (1-Vr)*(1-Gr)*b*a;
	ee[4] = Vr + (1-Vr)*Gr*(1-a) + (1-Vr)*(1-Gr)*b*(1-a);
	ee[5] = (1-Vr)*(1-Gr)*(1-b);

	ee[6] = Vtx + (1-Vtx)*Gtx*a + (1-Vtx)*(1-Gtx)*b*a;
	ee[7] = (1-Vtx)*Gtx*(1-a) + (1-Vtx)*(1-Gtx)*b*(1-a);
	ee[8] = (1-Vtx)*(1-Gtx)*(1-b);

	ee[9] = (1-Vrx)*Grx*a + (1-Vrx)*(1-Grx)*b*a;
	ee[10] = Vrx + (1-Vrx)*Grx*(1-a) + (1-Vrx)*(1-Grx)*b*(1-a);
	ee[11] = (1-Vrx)*(1-Grx)*(1-b);

	ee[12] = b*a;
	ee[13] = b*(1-a);
	ee[14] = (1-b);
	
	ee = ee*NN;
	
  double LL = 0.0;
  
  for(int ii=0; ii<15; ii++){
    if(chd[ii] > 0.0){
      LL = LL + chd[ii]*(log(chd[ii])-log(ee[ii]));
    }
  }
  return LL;
}






