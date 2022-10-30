
#ifndef VARSOLVER_H
#define VARSOLVER_H


#include <math.h>
#include <Eigen/Dense>
#include "constants.h"


class VarSolver
{
public:
  void set_mass_z(double m1, double m2, double m3, double z1, double z2, double z3);
  VarSolver(int _nmax, int _n34_max, int _J ):nmax(_nmax),n34_max(_n34_max),J(_J){
    if ((n34_max > nmax) || (nmax > n0max_arr) || (n34_max > nmax_arr))
      throw std::runtime_error("Check bounds!");
    tab_vib = new int[(nmax+1)*(nmax+1)*(nmax+1)*(2*nmax+1)];
    totmem += (nmax+1)*(nmax+1)*(nmax+1)*(2*nmax+1)*sizeof(int);
    quannummat = new int[(nmax+1)*(nmax+1)*(nmax+1)*(2*nmax+1)] [4];
    totmem += (nmax+1)*(nmax+1)*(nmax+1)*(2*nmax+1)*sizeof(int)*4;

    pol0 =  new double[n0arr_dim];
    pol3 =  new double[nmax_arr+1];
    pol4 =  new double[nmax_arr+1];
    pol0_der =  new double[n0arr_dim];
    pol0_secder =  new double[n0arr_dim];
    pol3_der =  new double[nmax_arr+1];
    pol3_secder =  new double[nmax_arr+1];
    pol4_der =  new double[nmax_arr+1];
    pol4_secder =  new double[nmax_arr+1];

    Mm12 = Eigen::MatrixXd::Zero(9,9);
    X0_vec = Eigen::MatrixXd::Zero(9,1);

    totmem += ((n0arr_dim)+(nmax_arr+1) + (nmax_arr+1))*3*sizeof(double);


  };

  ~VarSolver()
  {
    delete [] tab_vib;
    delete [] quannummat;

    delete [] pol0;
    delete [] pol0_der;
    delete [] pol0_secder;
  //  delete [] pol2;
    delete [] pol3;
    delete [] pol3_der;
    delete [] pol3_secder;

    delete [] pol4;
    delete [] pol4_der;
    delete [] pol4_secder;
  }
  void (*setup_pot)(void) = nullptr;
  void (*extpot)(double &, double &, double &, double &) = nullptr;
  double pot(double * sval)//r1,r2,s3x,s3y
  {
    double res;
    double ang = asin(sqrt(sval[2]*sval[2]+sval[3]*sval[3]+sval[4]*sval[4]));
    double r1 = sval[0]/constants::atob;
    double r2 = sval[1]/constants::atob;
    (*extpot)(r1,r2,ang,res );
    return res;
  }
  //double pot_s3()
  void harmonic_solver();

  void print_normal_coords(std::ofstream & ,constants::matr_type);

  void compute_s (double * xval, double * sval);

  void fill_zet();

  void fill_tabvib();

  void print_mem();

  void check_tabvib();

  void fill_hamiltonian();

  double  mulin( double qq1, double qq2, double qq3, double qq4 );

  void diagonalize();

  void prepare_for_business();




private:

  double totmem = 0;
  double m[3] ;// {11.996709*em, 15.990525*em, 15.990525*em};
  double Z[3];
  int nmax;
  int n34_max;
  int J;
  int Dim = 0;
  const int n0max_arr = 30;
  const int nmax_arr = 30;
  //const int lmax_arr = n0max_arr;
  //const int kmax = n0max_arr-1;
  const int n0arr_dim = (n0max_arr*n0max_arr+3*n0max_arr)/2 + 1;

  static const int Nquad = 16;
  static double pts[Nquad];
  static double wei[Nquad];
  static double pts_lag[Nquad];
  static double wei_lag[Nquad];
  // double lam1;
  // double lam2;
  // double lam3;
  // double lam4;
  int nmodes = 4;
  double lambdas[4];

//  int tab_vib[nmax+1][nmax+1][nmax+1][2*nmax+1];
  int * tab_vib;
  int (* quannummat)[4];

  double * pol0;
  double * pol3;
  double * pol4;
  double * pol0_der;
  double * pol0_secder;
  double * pol3_der;
  double * pol3_secder;
  double * pol4_der;
  double * pol4_secder;

  Eigen::MatrixXd Mm12;
  Eigen::VectorXd X0_vec;
  Eigen::MatrixXd Hamiltonian_full;

  bool mm12_filled = false;
  bool masses_filled = false;
  bool zet_filled = false;
  bool harmonic_solved = false;



  double zet[4][4][3];
//  int Dim = 0;

  Eigen::MatrixXd evec;

  int levi(int i, int j, int k)
  {
    if((i==j)||(i==k)||(k==j)) return 0;
     else return ( (i-j)*(j-k)*(k-i)/2 );
  }

  int cor(int bet, int i){
    return 3*(i-1)+bet;
  }
  int del(int a, int b)
  {
    return (a == b) ? 1 : 0;
  }

  int tab_vib_idx(int i, int j, int k, int l)
  {
    return (nmax+1)*(nmax+1)*(nmax+1)*l + (nmax+1)*(nmax+1)*k + (nmax+1)* j +i;
  }

  int lag_ind(int n, int l)
  {
    int k = n-1;
    return  (k*k+3*k)/2 + 1 + (l+n)/2 ;
  }

  double pi2_dl0(double qq0, double qq3, double qq4, int n0, int l, int n3, int n4 );
  //llef = lrig+1
  void pix_dlp1(double & Re,double & Im, double qq0, double qq3, double qq4, int n0, int l, int n3, int n4 );

  //llef = lrig-1
  void pix_dlm1(double & Re,double & Im, double qq0, double qq3, double qq4, int n0, int l, int n3, int n4);


  void piy_dlm1(double & Re,double & Im, double qq0, double qq3, double qq4, int n0, int l, int n3, int n4 );


  void piy_dlp1(double & Re,double & Im, double qq0, double qq3, double qq4, int n0, int l, int n3, int n4 );
  //Real
  double Pix (int llef, int lrig);


  //Imag
  double Piy (int llef, int lrig);


  double Pixx (int llef, int lrig);


  double Piyy (int llef, int lrig);

  void  fill_lag(double x, double * wfarr, double *wfder,double * wfsecder);
  void fill_wf(double x, double * wfarr, double * wfder, double * wfsecder);






};




#endif
