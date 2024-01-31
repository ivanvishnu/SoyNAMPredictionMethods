#include <Rcpp.h>
using namespace Rcpp;
//[[Rcpp::export]]
SEXP LD(NumericVector p1,NumericVector q1,NumericVector a, NumericMatrix pro,NumericMatrix D,NumericMatrix D1,NumericMatrix R2){
	
int m = D.ncol();
int n = a.size();
int s = pro.ncol();
double pab = 0;

for (int j = 0; j < m-1; j++) {
	 		
 for (int i = a(j); i < n; i++) { 
	
  for (int k = 0; k< s; k++) {     	
	if(pro(j,k)==1 && pro(i,k)==1){
		pab += 1;
	}
  }	
  pab=pab/s;
  D(i,j)=pab-(p1(j)*p1(i));

  if(D(i,j) >0){ 

      double PQ1 = p1(i)*q1(j);
      double PQ2 = q1(i)*p1(j); 
      double Dmax = std::min(PQ1,PQ2);
      D1(i,j) = D(i,j)/Dmax ;
  } 
  if(D(i,j) <0){
     double PQ1 = p1(i)*p1(j);
     double PQ2 = q1(i)*q1(j); 
     double Dmax = std::min(PQ1,PQ2); 
     D1(i,j) = D(i,j)/Dmax ;
  } 
  if(D(i,j)==0){ 
    D1(i,j) = D(i,j);
  }
 
  R2(i,j)=0;
  if(p1(j) !=0 && p1(i)!=0 && q1(j) !=0 && q1(i)!=0){
	R2(i,j)= (D(i,j)*D(i,j))/(p1(j)*p1(i)*q1(j)*q1(i));
  }
  
  pab = 0;} 
  
 }
  return List::create(Named("D")=D,Named("D1")=D1,Named("R2")=R2);
}




