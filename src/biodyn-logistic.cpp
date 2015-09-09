// Implements the logistic function
//PELLA , J. J., 1967 A study of methods to estimate the Schaefer model parameters with special reference to the yellowfin tuna fishery in the eastern tropical Pacific ocean. University of Washington, Seattle.
//PELLA , J. J., and P. K. T OMLINSON , 1969 A generalized stock production model. Bulletin of the Inter-American Tropical Tuna Commission 13: 419-496.
//PRAGER , M. H., 1994 A suite of extensions to a nonequilibrium surplus-production model. U. S. Fishery Bulletin 92: 374-389.

#include "../inst/include/biodyn-logistic.h"

logistic::logistic(FLQuant F, FLQuant C, FLQuant B, FLQuant params){
  hvt=F;
  ctc=C;
  stk=B;
  
  par=params;
  }
     
void logistic::nr(double tolVal, int niter){
    for (int i=1; i<=ctc.get_niter(); i++){
      for (int yr=1; yr<=ctc.get_nyear()-1; yr++){
        int  iters=0;
        double val=1;
        while (iters<niter){
         iters++;
         double func = ctc(1,yr,1,1,1,i)-yield(hvt(1,yr,1,1,1,i),stk(1,yr,1,1,1,i),par(1,1,1,1,1,i),par(2,1,1,1,1,i));
         double grad = gradY(hvt(1,yr,1,1,1,i),ctc(1,yr,1,1,1,i),stk(1,yr,1,1,1,i),par(1,1,1,1,1,i),par(2,1,1,1,1,i));
         val  = NewRhap(func,grad);
        
         hvt(1,yr,1,1,1,i)-=val;
         }
       
   stk(1,yr+1,1,1,1,i)=stock(hvt(1,yr,1,1,1,i),ctc(1,yr,1,1,1,i),
                             stk(1,yr,1,1,1,i),par(1,1,1,1,1,i),
                             par(2,1,1,1,1,i));
   }}}
     

double logistic::gradY(double F, double C, double  B, double r, double K){

  double expr1  = r/K;
  double expr2  = F/expr1;
  double expr3  = expr1 * B;
  double expr4  = r - F;
  double expr5  = exp(expr4);
  double expr7  = expr3 * (1 - expr5);
  double expr9  = 1 - expr7/expr4;
  double expr10 = log(expr9);

  return(-(1/expr1*expr10-expr2*((expr3*expr5/expr4 + expr7/(expr4*expr4))/expr9)));
  }
  
  
 logistic::~logistic(void){};    

 
//[[Rcpp::export]]  
FLQuant stockCPP(FLQuant F,FLQuant C,FLQuant B, FLQuant params, double tol=1e-10, int niter=200) {
	logistic nr(F,C,B,params);
	
  nr.nr(tol,niter);
	
	return(nr.stk);}

//[[Rcpp::export]]  
FLQuant fCPP(FLQuant F, FLQuant C, FLQuant B, FLQuant params, double tol=1e-10, int niter=200) {
	
	logistic nr(F,C,B,params);
	
	nr.nr(tol,niter);
	
	return(nr.hvt);}
