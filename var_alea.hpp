#ifndef VAR_ALEA_HPP
#define VAR_ALEA_HPP
#include <cmath>
#include <ctime>
#include "mt19937.h"

void init_alea(unsigned seed = static_cast<unsigned>(std::time(0))) {
	init_genrand(seed);
};

double Pi = 3.141592653;

template <typename T>
struct var_alea {
  typedef T result_type;
  var_alea() : value(0) {};
  var_alea(T value) : value(value) {};
  virtual ~var_alea() {};
  virtual T operator()() = 0;
  virtual double pdf (double x) = 0;
  T current() const { return value; };
protected:
  T value;
};

struct uniform : public var_alea<double>
{
  uniform(double left = 0, double right = 1)
    : left(left), right(right), size(right-left), genrand(genrand_real3) {};
  double operator()() {
    return value = left + size * genrand();
  };
  double pdf(double x){
    if (x>left&&x<right){return 1/(right-left);}
    else {return 0;}
  };
private:
   double left, right, size;
   double (*genrand)(void);
};

struct expo : public var_alea<double>
{
  expo(double lambda) : inv_lambda(1./lambda),lambda(lambda), U(0,1) {};
  double operator()() {
    return value = - inv_lambda * log(U());
  };
  double pdf (double x){
   if (x<=0){return 0;}
   else {return lambda* exp ( -lambda *x);}
  };
private:
   double inv_lambda, lambda;
   uniform U;
};

struct gaussian : public var_alea<double>
{
  gaussian(double mean = 0, double std = 1)
    : mean(mean), std(std), flag(true), unif(-1,1) {};
  double operator()() {
    flag = !flag;
    if (!flag) {
      do {
	U = unif(); V = unif();
	R2 = U*U + V*V;
      } while (R2 > 1);
      rac = sqrt(-2 * log(R2) / R2);
      return value = mean + std * U * rac;
    } else
      return value = mean + std * V * rac;
  };
  double pdf (double x){return exp(-(x-mean)*(x-mean)/(2*std))/(sqrt(2*Pi)*std);};
  double SimpleMonteCarlo (double x, int n){
	double A=0;
	uniform U(0,1);
	for (int i =0; i<n; i++){ A+= exp(-pow((log(U())-mean+x),2)/(2*std))/(sqrt(2*Pi)*std*(U.current()));};
	return A/n;};
  double MonteCarloAnthetique (double x, int n){
        double A=0;
	uniform U(0,1);
	if (x<mean){
	for (int i =0; i<n/2; i++)
		{double y=U(); A+= exp(-pow((log(y)-mean+x),2)/(2*std))/(sqrt(2*Pi)*std*y)+exp(-pow((log(1-y)-mean+x),2)/(2*std))/(sqrt(2*Pi)*std*(1-y));};}
	else {for (int i =0; i<n/2; i++)
		{double y=U(); A+=1+(x-mean)*exp(-pow(y*(x-mean),2)/(2*std))/(sqrt(2*Pi)*std)+(x-mean)*exp(-pow((1-y)*(x-mean),2)/(2*std))/(sqrt(2*Pi)*std);};}
	return A/n;};
  double moyenne(){return mean;};
  double var(){return std;};
private:
  double mean, std, U, V, R2, rac;
  bool flag;
  uniform unif;
};


struct chi_deux : public var_alea<double>
{
  chi_deux(int n) : n(n), G(0,1) {};
  double operator()() {
    value = 0;
    for (int j = 0; j < n; j++) value += G()*G.current();
    return value;
  };
  double pdf (double x){ return  5.0*pow(x,n/2-1)*exp (-x/2)*pow(0.5,n/2)/ tgamma (n/2);};
private:
  int n;
  gaussian G;
};


double gamma (double x, double y){ return sqrt (x*x-y*y);};


struct t_student : public var_alea <double>
{
   t_student (int n) : n(n), G(0,1), X(n) {};
   double operator ()(){
   return G()/sqrt(X()/n);
   };
   double pdf (double x) { return tgamma((n+1)/2)*exp (-(n+1)*log(1+x*x/n)/2)/(sqrt(n*Pi)*tgamma(n/2));};
   double SimpleMonteCarlo (double x, int m)
	{ uniform U(0,1);
	  double A=0;
	  for (int i=0; i<m;i++) {double y=U();A+=tgamma((n+1)/2)*exp(-(n+1)*log(1+pow(log(y*exp(x)),2)/n)/2)/(sqrt(n*Pi)*y*tgamma(n/2));};
	  return A/m;};
   double MonteCarloAnthetique (double x, int m){
        double A=0;
	uniform U(0,1);
	if (x<this-> moyenne()){
	for (int i =0; i<m/2; i++)
		{double y=U(); A+=tgamma((n+1)/2)*exp (-(n+1)*log(1+pow(x+log(y),2)/n)/2)/(sqrt(n*Pi)*y*tgamma(n/2))+tgamma((n+1)/2)*exp(-(n+1)*log(1+pow(x+log(1-y),2)/n)/2)/(sqrt(n*Pi)*(1-y)*tgamma(n/2));};}
	else {for (int i =0; i<m/2; i++)
		{double y=U(); A+=1+tgamma((n+1)/2)*x*exp(-(n+1)*log(1+x*y*y*x/n)/2)/(sqrt(n*Pi)*tgamma(n/2))+tgamma((n+1)/2)*x*exp(-(n+1)*log(1+x*(1-y)*(1-y)*x/n)/2)/(sqrt(n*Pi)*tgamma(n/2));};}
	return A/m;};
  double moyenne () {return 0;};
  double var () {return n/(n-2) ;};

private :
    int n;
    gaussian G;
    chi_deux X;
};


struct inverse_gaussian : public var_alea<double>
{
  inverse_gaussian(double lambda, double mu)
    : lambda(lambda), mu(mu), Y(1), U(0,1) {};
  double operator()() {
    double Z = mu + 0.5*mu*mu/lambda*Y();
    double rac = sqrt(Z*Z - mu*mu);
    return value = (U() < mu/(mu+Z+rac)) ? Z+rac : Z-rac;
  };
  double pdf (double x){if (x>0) {return  lambda*exp(-pow(lambda-mu*x,2)/(2*mu*x))/(sqrt(2*Pi*mu*x)*x) ;} else {return 0;}};
private:
  double lambda, mu;
  chi_deux Y;
  uniform U;
};

struct normal_inverse_gaussian : public var_alea<double>
{
  normal_inverse_gaussian(double alpha, double beta, double mu, double delta)
    : alpha(alpha), beta(beta), mu(mu), delta(delta), G(0,1),
      Y(delta/gamma(alpha,beta), delta*delta)  {};
  double operator()() {
    double y_ = Y();
    return value = mu + beta*y_ * sqrt(y_) * G();
  };
  double pdf (double x) { return delta*alpha*exp(delta*gamma(alpha,beta)+beta*(x-mu))*5/(Pi*sqrt(delta*delta+(x-mu)*(x-mu)));}; // jamais utilisé
  double SimpleMonteCarlo (double x, int n)
	{ gaussian GA(0,1);
	  inverse_gaussian YIG( delta*gamma(alpha,beta), pow(gamma(alpha,beta),2));
	  uniform U(0,1);
	  double A=0;
	  for (int i=0; i<n;i++) {double y=U();A+=GA.SimpleMonteCarlo((x-(mu-beta*log(y)))/sqrt(-log(y)),2*n)*YIG.pdf(-log(y))/y;};
	  return A/n;};
/*  double MonteCarloAnthetique (double x, int n)
	{ gaussian GA(0,1);
	  inverse_gaussian YIG( delta*gamma(alpha,beta), pow(gamma(alpha,beta),2));
	  uniform U(0,1);
	  double A=0;
	  for (int i=0; i<n/2;i++) {double y=U();A+=GA.cdfMonteCarlo((x-(mu-beta*log(y)))/sqrt(-log(y)),2*n)*YIG.pdf(-log(y))/y+GA.cdfMonteCarlo((x-(mu-beta*log(1-y)))/sqrt(-log(1-y)),2*n)*YIG.pdf(-log(1-y))/(1-y);};
	  return (A/2)/(n/2);}; */ 
  double moyenne () {return mu + delta*beta/gamma(alpha, beta);};
  double var () {return delta*alpha*alpha/(pow(gamma(alpha,beta),3)) ;};
 
private:
  double alpha, beta, mu, delta;
  gaussian G;
  inverse_gaussian Y;
};

template < typename T>
void remplirTableau ( T X, double A[], double B[],int n , int m) {
//  int n = sizeof(A)/sizeof(double);
  double mu = X.moyenne();
  double variance = X.var();
  for (int i=0; i<=n; i++){ A[i] = mu - 2* variance + 4*variance*i/n ;};
  for (int i=0; i<=n; i++){ B[i] =  X.SimpleMonteCarlo (A[i] , m) ; }
};

template < typename T>
void remplirTableauAnthetique ( T X, double A[], double B[], int n, int m) {
//  int n = sizeof(A)/sizeof(double);
  double mu = X.moyenne();
  double variance = X.var();
  for (int i=0; i<=n; i++){ A[i] = mu - 2* variance + 4*variance*i/n ;};
  for (int i=0; i<=n; i++){ B[i] =  X.MonteCarloAnthetique (A[i] , m) ; }
};


double inverseTableau (double y, double A[], double B[], int m){
//  int  m = sizeof(A)/sizeof(double);
  double a=B[0];
  int i=0;
  if (y<=B[0]) {return (A[1]-A[0])*B[0]*log(y/B[0])/(B[1]-B[0])+A[0];} //la CDF est approximée avec y=B[0]exp((x-A[0])*(B[1]-B[0])/((A[1]-A[0])*B[0])), prolongement C1
  else if (y>B[m]) {return A[m]-(1-B[m])*(A[m]-A[m-1])*log((1-y)/(1-B[m]))/(B[m]-B[m-1]) ;} //la CDF est approximée avec la fonction y=1-(1-B[n])exp((-x+A[n])*(B[n]-B[n-1])/((1-B[n])*(A[n]-A[n-1])))
  else while (y>a){i++;a=B[i];}
  double u=(y-B[i])/(B[i+1]-B[i]); // on peut supprimer cette ligne pour remplacer u dans la suivante.
  return u*(A[i+1]-A[i])+A[i];};

 normal_inverse_gaussian fonction (double a, double alpha, double beta)
	{ normal_inverse_gaussian Z(alpha/a, beta/a, -beta*pow(gamma (alpha,beta),2)/(a*pow(alpha,2)),pow(gamma(alpha,beta),3)/(a*pow(alpha,2)));
	  return Z;};

/*double fonction_correlation (double a){
	return sqrt(1-a*a)/a;}; 

double calculCDONIG(double x, double C, double a, double alpha, double beta, double A[], double B[], double m){
	normal_inverse_gaussian T = (1/fonction_correlation (a), alpha, beta);
	normal_inverse_gaussian Tinfini = (1, alpha, beta);
	return (1-Tinfini((C-sqrt(1-a*a)*T.inverseTableau(x,A,B,m))/a  */
#endif	
