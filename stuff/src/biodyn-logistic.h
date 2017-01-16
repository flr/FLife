// Implements the logistic function
//PELLA , J. J., 1967 A study of methods to estimate the Schaefer model parameters with special reference to the yellowfin tuna fishery in the eastern tropical Pacific ocean. University of Washington, Seattle.
//PELLA , J. J., and P. K. T OMLINSON , 1969 A generalized stock production model. Bulletin of the Inter-American Tropical Tuna Commission 13: 419-496.
//PRAGER , M. H., 1994 A suite of extensions to a nonequilibrium surplus-production model. U. S. Fishery Bulletin 92: 374-389.


#ifndef _INC_logistic
#define _INC_logistic

#include <FLasher.h>

#define const_nan 0.0
#define pow(x) x*x

#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define flq(x,i,j,k,l,m,n) x(MIN(1,x.get_nquant()),MIN(1,x.get_nyear()),MIN(1,x.get_munit()),MIN(1,x.get_nseason()),MIN(1,x.get_narea()),MIN(1, x.get_niter()))
  
                                 
class logistic
{
public:        
   logistic(FLQuant,FLQuant,FLQuant,FLQuant);
   
   FLQuant ctc;
   FLQuant stk;
   FLQuant hvt;

   void nr(double tolVal=1e-5,int niter=100);
   
  ~logistic(void);      
 
protected: 
  FLQuant par;
   
  // reparameterisation of logistic
  inline double alpha(double r, double F) {return(r-F);}
  inline double beta( double r, double K) {return(r/K);}

  // estimation
  inline double NewRhap(double x, double func, double grad) {return(x-func/grad);}
  
  //yield functions
  double gradY(double F, double C, double  B, double r, double K);
  inline double yield(   double F,double B,double r,double K){ 
     return(F/(r/K)*log(1-(r/K)*B*(1-exp((r-F)))/(r-F)));}
  inline double gradMinY(double F, double C, double  B, double r, double K){
     return(2*(C-yield(F,B,r,K))*gradY(F,C,B,r,K));}

  // F functions
  double gradF(double F, double C, double  B, double r, double K);
  inline double f(double F, double C, double  B, double r, double K){
     return(C*(r/K)/(log(r/K*B*(exp(r-F)-1)/(r-F)+1)));}
  inline double gradMinF(double F, double C, double  B, double r, double K){
     return(2*(F-f(F,C,B,r,K))*gradF(F,C,B,r,K));}

  // estimate stock size
  inline double stock(double F, double C, double  B, double r, double K){
    return(alpha(r,F)*B*exp(alpha(r,F))/(alpha(r,F)+beta(r,K)*B*(exp(alpha(r,F))-1)));}
  };
   
FLQuant     fCPP(FLQuant, FLQuant, FLQuant, FLQuant);
FLQuant stockCPP(FLQuant, FLQuant, FLQuant, FLQuant);
#endif /* _INC_logistic */
