#include <math.h>
#include <Eigen/Dense>
//#include <Eigen/SVD>
#include <Spectra/SymEigsSolver.h>

#include <iostream>
#include <fstream>
#include <iomanip>
#include <gsl/gsl_sf_gamma.h>
#include "VarSolver.h"



double VarSolver::pts[VarSolver::Nquad] = {-4.6887389393058184e+00, -3.8694479048601227e+00, -3.1769991619799560e+00, -2.5462021578474814e+00, -1.9517879909162540e+00, -1.3802585391988808e+00, -8.2295144914465589e-01, -2.7348104613815245e-01, 2.7348104613815245e-01, 8.2295144914465589e-01, 1.3802585391988808e+00, 1.9517879909162540e+00, 2.5462021578474814e+00, 3.1769991619799560e+00, 3.8694479048601227e+00, 4.6887389393058184e+00};
double VarSolver::wei[VarSolver::Nquad] = {2.6548074740111822e-10, 2.3209808448652107e-07, 2.7118600925378815e-05, 9.3228400862418053e-04, 1.2880311535509974e-02, 8.3810041398985829e-02, 2.8064745852853368e-01, 5.0792947901661374e-01, 5.0792947901661374e-01, 2.8064745852853368e-01, 8.3810041398985829e-02, 1.2880311535509974e-02, 9.3228400862418053e-04, 2.7118600925378815e-05, 2.3209808448652107e-07, 2.6548074740111822e-10};

double VarSolver::pts_lag[VarSolver::Nquad] = {8.7649410478927800e-02, 4.6269632891508100e-01, 1.1410577748312300e+00, 2.1292836450983800e+00, 3.4370866338932100e+00, 5.0780186145497700e+00, 7.0703385350482300e+00, 9.4383143363919400e+00, 1.2214223368866200e+01, 1.5441527368781600e+01, 1.9180156856753100e+01, 2.3515905693991900e+01, 2.8578729742882100e+01, 3.4583398702286600e+01, 4.1940452647688300e+01, 5.1701160339543300e+01};
double VarSolver::wei_lag[VarSolver::Nquad] ={2.0615171495782900e-01, 3.3105785495066000e-01, 2.6579577764289200e-01, 1.3629693429889600e-01, 4.7328928715557700e-02, 1.1299900055889500e-02, 1.8490709598650400e-03, 2.0427191334507600e-04, 1.4844586517905300e-05, 6.8283196322608000e-07, 1.8810249400935900e-08, 2.8623501570998700e-10, 2.1270794660234200e-12, 6.2979668760995700e-15, 5.0504735990499800e-18, 4.1614623616631900e-22};


void VarSolver::set_mass_z(double m1, double m2, double m3, double z1, double z2, double z3)
{
  masses_filled = true;
  Z[0] = z1;
  Z[1] = z2;
  Z[2] = z3;

  m[0] = m1;
  m[1] = m2;
  m[2] = m3;

  X0_vec(2) = z1;
  X0_vec(5) = z2;
  X0_vec(8) = z3;
}

// double VarSolver::s1(double * xval)
// {
//   double pot = compute_s1(xval[0], xval[1], xval[2],
//               xval[3], xval[4],xval[5],
//                 xval[6], xval[7], xval[8]);
//
//   return pot;
// }
//
// double VarSolver::s2(double * xval)
// {
//   double pot = compute_s2(xval[0], xval[1], xval[2],
//               xval[3], xval[4],xval[5],
//                 xval[6], xval[7], xval[8]);
//
//   return pot;
// }
//
// double VarSolver::s3x(double *xval)
// {
//   double x1 = xval[0];
//   double y1 = xval[1];
//   double z1 = xval[2];
//   double x2 = xval[3];
//   double y2 = xval[4];
//   double z2 = xval[5];
//   double x3 = xval[6];
//   double y3 = xval[7];
//   double z3 = xval[8];
// //        i       j     k
// //  r2 = x3-x1, y3-y1, z3-z1
// //  r1 = x2-x1, y2-y1, z2-z1;
// //  double
//   double r2xr1_x = (y3-y1)*(z2-z1)-(z3-z1)*(y2-y1);
//   double r1 = sqrt((x2-x1)*(x2-x1) + (y2-y1)*(y2-y1) + (z2-z1)*(z2-z1));
//   double r2 = sqrt((x3-x1)*(x3-x1) + (y3-y1)*(y3-y1) + (z3-z1)*(z3-z1));
//   return -r2xr1_x/r1/r2;
//
// }
//
// double VarSolver::s3y(double *xval)
// {
//   double x1 = xval[0];
//   double y1 = xval[1];
//   double z1 = xval[2];
//   double x2 = xval[3];
//   double y2 = xval[4];
//   double z2 = xval[5];
//   double x3 = xval[6];
//   double y3 = xval[7];
//   double z3 = xval[8];
// //        i       j     k
// //  r2 = x3-x1, y3-y1, z3-z1
// //  r1 = x2-x1, y2-y1, z2-z1;
// //  double
//   double r2xr1_y = -((x3-x1)*(z2-z1)-(z3-z1)*(x2-x1));
//   double r1 = sqrt((x2-x1)*(x2-x1) + (y2-y1)*(y2-y1) + (z2-z1)*(z2-z1));
//   double r2 = sqrt((x3-x1)*(x3-x1) + (y3-y1)*(y3-y1) + (z3-z1)*(z3-z1));
//   return r2xr1_y/r1/r2;
//
// }

void VarSolver::compute_s (double * xval, double * sval)
{
  double x1 = xval[0];
  double y1 = xval[1];
  double z1 = xval[2];
  double x2 = xval[3];
  double y2 = xval[4];
  double z2 = xval[5];
  double x3 = xval[6];
  double y3 = xval[7];
  double z3 = xval[8];
  //        i       j     k
  //  r2 = x3-x1, y3-y1, z3-z1
  //  r1 = x2-x1, y2-y1, z2-z1;
  //  double
  double r2xr1_x = (y3-y1)*(z2-z1)-(z3-z1)*(y2-y1);
  double r2xr1_y = -((x3-x1)*(z2-z1)-(z3-z1)*(x2-x1));
  //double r2xr1_z = 0.0;
  // if (full_s3)
  //   r2xr1_z = (x3-x1)*(y2-y1) - (y3-y1)*(x2-x1);
  double r1 = sqrt((x2-x1)*(x2-x1) + (y2-y1)*(y2-y1) + (z2-z1)*(z2-z1));
  double r2 = sqrt((x3-x1)*(x3-x1) + (y3-y1)*(y3-y1) + (z3-z1)*(z3-z1));
  // s2x = ;
  // s2y = ;
  // s1 = r1;
  // s2 = r2;

  sval[0] = r1;
  sval[1] = r2;
  sval[2] = -r2xr1_x/r1/r2;
  sval[3] = r2xr1_y/r1/r2;
  sval[4] = ((x3-x1)*(y2-y1) - (y3-y1)*(x2-x1))/r1/r2;


  // if (full_s3)
  // {

  // }
  // else
  //   sval[4] = 0.0;

//return ;
}

void VarSolver::harmonic_solver()
{
  if (!masses_filled)
    throw std::runtime_error("Fill masses first!");

  harmonic_solved = true;


  Eigen::MatrixXd G = Eigen::MatrixXd::Zero(9,9);
  Eigen::MatrixXd F = Eigen::MatrixXd::Zero(9,9);
//  Eigen::MatrixXd Mm12 = Eigen::MatrixXd::Zero(9,9);

  Eigen::MatrixXd Gtrue = Eigen::MatrixXd::Zero(4,4);
  Eigen::MatrixXd Ftrue = Eigen::MatrixXd::Zero(4,4);
  Eigen::MatrixXd B = Eigen::MatrixXd::Zero(4,9);

  double step = 1E-4;
 double xmin[9] = {0.0,0.0,Z[0],0.0,0.0,Z[1],0.0,0.0,Z[2]};
  double currx[9];

  double smin[5] = {Z[2],Z[2],0.0,0.0,0.0};
  double currs[5];


  double xtest [] = {0.4,0.2,-0.5,-0.32,0.543,Z[1],0.42,-0.123,Z[2]};

  // double s3_test = s3_wrapper(xtest);
  // double s3x_test = s3x_wrapper(xtest);
  // double s3y_test = s3y(xtest);
  // std::cout << "s3 " << (s3_test) << " " << (asin(sqrt(s3x_test*s3x_test+s3y_test*s3y_test)))<< "\n";


  for (int i = 0; i < 9; i++) {
    G(i,i) = 1.0/m[i/3];
    Mm12(i,i) = 1.0/sqrt(m[i/3]);
    mm12_filled = true;
  }
  for (int i = 0; i < 4; ++i)
    for (int j= 0; j < 4;++j)
    {
      if( ((i==2)&&(j==2)) || ((i==3)&&(j==3)))
        step = 1E-3;
        else
          step=1E-4;
      if (i != j)
      {
        memcpy(currs,smin, sizeof smin);
        currs[i] += step;
        currs[j] += step;
        Ftrue(i,j) += pot(currs);//+ +
        currs[i] -= 2*step;
        Ftrue(i,j) -= pot(currs);//- +
        currs[j] -= 2*step;
        Ftrue(i,j) += pot(currs);//- -
        currs[i] += 2*step;
        Ftrue(i,j) -= pot(currs);//+ -
        Ftrue(i,j) *= 1/4.0/step/step;
      }
      else
      {

        memcpy(currs,smin, sizeof smin);
        Ftrue(i,j) -= 2.0*pot(currs);
        // if ((i==2) && (j==2))
        //  std::cout << "mat "<< 2.0*pot(currs) << "\n";
        currs[i] += step;
        Ftrue(i,j) += pot(currs);
        // if ((i==2) && (j==2))
        //  std::cout << "mat "<< pot(currs) << "\n";
        currs[i] -= 2*step;
        Ftrue(i,j) += pot(currs);
        // if ((i==2) && (j==2))
        // { std::cout << "mat "<< pot(currs) << "\n";
        //  std::cout << Ftrue(i,j) << "\n";}
        Ftrue(i,j) *= 1/step/step;

      }

    //  F(i,j) *= 1/sqrt(m[i/3])/sqrt(m[j/3]);

    }

    step = 1E-2;
    //double sval[5];

    for(int i = 0; i < 9; ++i)
    {
      memcpy(currx,xmin, sizeof xmin);
      currx[i] += step;
      compute_s(currx,currs);
      for (int k = 0; k < 4; ++k)
        B(k,i) += 1/2.0*currs[k];
      currx[i] -= 2*step;
      compute_s(currx,currs);
      for (int k = 0; k < 4; ++k)
      {
        B(k,i) += -1/2.0*currs[k];
        B(k,i) *= 1/step;
      }

    }


    Gtrue = B*G*B.transpose();

    Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd> es(Ftrue,Gtrue.inverse());
  //  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(Gtrue);
       std:: cout << "\n"<<Ftrue <<"\n";
        std:: cout << "\n"<<Gtrue <<"\n";
        std::cout << Gtrue*Ftrue <<"\n";


    std::cout << "Harmonic frequencies:" << std::endl << std::setprecision(9) << (es.eigenvalues())*constants::htocm  << std::endl;
    std::cout << "Harmonic frequencies:" << "\n";
    for(int i =0; i < nmodes; ++i)
    {
      std::cout << sqrt(es.eigenvalues()(i))*constants::htocm <<"\n";
    }

//     std::cout <<"Matrix" <<"\n";
//     for (int i =0; i < 9; ++i)
//     {
//       std::cout << es.eigenvectors()(i,5) << std::endl;
//     }
     Eigen::MatrixXd evecs = es.eigenvectors();
    // lam1 = sqrt(es.eigenvalues()(0));
     //lam2 = sqrt(es.eigenvalues()(1));
     //lam3 = sqrt(es.eigenvalues()(2));
     //lam4 = sqrt(es.eigenvalues()(3));
//
//     //std::cout <<evec.transpose()*evec << std::endl;
//   //std::cout << evec.transpose()*F*evec*htocm <<std::endl;
//
// //  std::cout << evec.transpose()*G.inverse()*evec <<std::endl;
     evec = (evecs.inverse()*B*Mm12).transpose();
     for (int i = 0; i < nmodes; ++i)
      lambdas[i] = sqrt(es.eigenvalues()(i));

    // std::cout<< (evecs.inverse()*B*Mm12).transpose();
  //   std::cout <<evec.transpose()*evec << std::endl;
}

void VarSolver::fill_zet()
{
  // ev[0] = lam1;
  // ev[1] = lam2;
  // ev[2] = lam3;
  // ev[3] = lam4;
  zet_filled = true;

  for (int k = 1; k <= 4; k++)
    for (int l = 1; l <= 4; l++)
      for (int alpha = 1; alpha <= 3; ++alpha)
      {
        zet[k-1][l-1][alpha-1] = 0.0;
        for (int i = 1; i <= 3; i++)
          for(int bet = 1; bet <= 3; bet++)
            for (int gam = 1; gam <= 3; gam++)
            {
               zet[k-1][l-1][alpha-1] += sqrt(lambdas[l-1]/lambdas[k-1])*levi(alpha,bet,gam)*evec(cor(bet,i)-1,k-1)*evec(cor(gam,i)-1,l-1);
            }
      }

}


void VarSolver::fill_tabvib()
{

//  std::ofstream ofile("tabvib.txt");
  int ctr = 0, n4;
  for (int curN = 0; curN <= nmax; ++curN)
    for(int n0 = 0; n0 <= curN; ++n0)
      for (int n3 = 0; n3 <= curN-n0; ++n3)
        for(int l = -n0; l <= n0; l+=2)
        {
      //  printf("fsf");
          if (abs(l)>J)
            continue;


          n4 = curN - n0-n3;

          if((n3 > n34_max) || (n4>n34_max))
            continue;
          std::cout << n0 << " " << n3 << " "<< n4 << " " << l << " = " << ctr  << std::endl;
          //ofile << n0 << " " << n3 << " "<< n4 << " " << l <<"\n";//<< " = " << ctr << std::endl;

          tab_vib[tab_vib_idx(n0,n3,n4,l+n0)] = ctr;

          quannummat[ctr][0] = n0;
          quannummat[ctr][1] = n3;
          quannummat[ctr][2] = n4;
          quannummat[ctr][3] = l;

          ++ctr;
          ++Dim;
        //  printf("[%d][%d][%d] = %d\n",n1,n2,n3,ctr );

        }

      //  ofile.close();
      //--Dim;
}

void VarSolver::print_mem()
{
  std::cout << "Total memory usage: " << totmem/1024.0/1024.0 << " MB" << std::endl;
}


void VarSolver::check_tabvib()
{

//  std::ofstream ofile("tabvib.txt");
  int ctr = 0, n4;
  for (int curN = 0; curN <= nmax; ++curN)
    for(int n0 = 0; n0 <= curN; ++n0)
      for (int n3 = 0; n3 <= curN-n0; ++n3)
        for(int l = -n0; l <= n0; l+=2)
        {
      //  printf("fsf");
          if (abs(l)>J)
            continue;


          n4 = curN - n0-n3;

          if((n3 > n34_max) || (n4>n34_max))
            continue;

          int var = tab_vib[tab_vib_idx(n0,n3,n4,l+n0)];


        //  std::cout << n0 << " " << n3 << " "<< n4 << " " << l << " = " << var  << std::endl;
          //ofile << n0 << " " << n3 << " "<< n4 << " " << l <<"\n";//<< " = " << ctr << std::endl;
          //std::cout << quannummat[var][0] <<" " << quannummat[var][1] << " " << quannummat[var][2] << " "<<quannummat[var][3]<< std::endl;

          if ((var != ctr) || (quannummat[var][0] != n0) || (quannummat[var][1] != n3) || (quannummat[var][2] != n4) || (quannummat[var][3] != l))
            throw std::runtime_error("Tabvib check failed!");
        //  tab_vib[tab_vib_idx(n0,n3,n4,l+n0)] = ctr;

          // quannummat[ctr][0] = n0;
          // quannummat[ctr][1] = n3;
          // quannummat[ctr][2] = n4;
          // quannummat[ctr][3] = l;

          ++ctr;
        //  ++Dim;
        //  printf("[%d][%d][%d] = %d\n",n1,n2,n3,ctr );

        }

      std::cout << "Tabvib check passed" << std::endl;

      //  ofile.close();
      //--Dim;
}
// void fill_hamiltonian(Eigen::MatrixXd &evec,Eigen::MatrixXd & Hamiltonian)
// {
//     double ** pol0 =  new double* [Nquad];
//     for (int i =0; i < Nquad; ++i)
//       pol0[i] = new double[n0arr_dim];
//
//     double ** pol3 = new double *[Nquad];
//     for (int i =0; i < Nquad; ++i)
//       pol3[i] =  new double[nmax_arr+1];
//
//       double ** pol4 = new double* [Nquad];
//       for (int i =0; i < Nquad; ++i)
//         pol4[i] =  new double[nmax_arr+1];
//
//         Eigen::MatrixXd Ham = Eigen::MatrixXd::Zero(Dim,Dim);
//         Eigen::MatrixXd Vmat = Eigen::MatrixXd::Zero(Dim,Dim);
//         Eigen::MatrixXd Hvibrot = Eigen::MatrixXd::Zero(Dim,Dim);
//         Eigen::MatrixXd Mumat = Eigen::MatrixXd::Zero(Dim,Dim);
//
//
//     //double * pol3 =  new double[nmax_arr+1];
//     //double * pol4 =  new double[nmax_arr+1];
//     double * pol0_der =  new double[n0arr_dim];
//     double * pol0_secder =  new double[n0arr_dim];
//     double * pol3_der =  new double[nmax_arr+1];
//     double * pol3_secder =  new double[nmax_arr+1];
//     double * pol4_der =  new double[nmax_arr+1];
//     double * pol4_secder =  new double[nmax_arr+1];
//
//
//     std::cout << "n0dim: " << n0arr_dim << std::endl;
//
//     std::ofstream ofile("potval.txt");
//
//
//     for (int qqq = 0; qqq < Nquad; ++qqq)
//     {
//       double Q0 = sqrt(pts_lag[qqq]);
//       fill_lag(Q0, pol0[qqq],pol0_der, pol0_secder);
//       double Q3 = pts[qqq];
//
//       fill_wf(Q3, pol3[qqq], pol3_der, pol3_secder);    //double Q3 = pts[qqq];
//       fill_wf(Q3, pol4[qqq], pol4_der, pol4_secder);
//     }
//     double * potvec = new double[Nquad*Nquad*Nquad];
//
//     std::ifstream herm_file("herm_vals.txt");
//     std::ifstream lag_file("lag_vals.txt");
//     // for (int i = 0; i < Nquad; i++)
//     //   for(int n = 0; n <= 16; ++n)
//     //   {
//     //     herm_file >> pol3[i][n];
//     //     pol4[i][n] = pol3[i][n];
//     //   }
//
//       // for (int i = 0; i < Nquad; i++)
//       //   for(int n = 0; n <= 16; n+=2)
//       //     {
//       //     lag_file >> pol0[i][lag_ind(n,0)];
//       //   //  pol4[i][n] = pol3[i][n];
//       //   }
//
//
//
//     herm_file.close();
//     lag_file.close();
//
//     for (int q0 = 0; q0 < Nquad; ++q0)
//     {
//         double Q0 = sqrt(pts_lag[q0]);
//         for (int q3 = 0; q3 < Nquad; ++q3)
//         {
//             double Q3 = pts[q3];
//         for (int q4 = 0; q4 < Nquad; ++q4)
//         {
//             double Q4 = pts[q4];
//
//             double Q1 = Q0*sin(0);
//             double Q2 = Q0*cos(0);
//             Eigen::VectorXd Qvec(4);
//           //  Qvec << 0,0,0,0,0,Q1/sqrt(lam1),Q2/sqrt(lam2), Q3/sqrt(lam3), Q4/sqrt(lam4);
//               Qvec << Q1/sqrt(lam1),Q2/sqrt(lam2), Q3/sqrt(lam3), Q4/sqrt(lam4);
//             Eigen::VectorXd Xvec;
//             Xvec = evec*Qvec;
//             double pot = V(Xvec(0)/sqrt(m[0])+X[0], Xvec(1)/sqrt(m[0])+Y[0], Xvec(2)/sqrt(m[0])+Z[0],
//                           Xvec(3)/sqrt(m[1])+X[1], Xvec(4)/sqrt(m[1])+Y[1], Xvec(5)/sqrt(m[1])+Z[1],
//                           Xvec(6)/sqrt(m[2])+X[2], Xvec(7)/sqrt(m[2])+Y[2], Xvec(8)/sqrt(m[2])+Z[2]);
//
//             potvec[pot_idx(q0,q3,q4)] = pot;
//             ofile << std::setprecision(14) <<std::scientific <<pot << "\n";
//           }
//         }
//       }
//
//     ofile.close();
//
//
//     for (int currig = 0; currig <= nmax; ++currig)
//     {
//       std::cout << "currig " << currig << "\n";
//       for (int n0rig = 0; n0rig <= currig; n0rig++)
//         for(int lrig = -n0rig; lrig <= n0rig; lrig+=2)
//           for (int n3rig = 0; n3rig <= currig - n0rig; ++n3rig)
//           {
//
//             if (abs(lrig)>J)
//               continue;
//
//               // if((n0rig+lrig) %2 !=0)
//               //   continue;
//             int n4rig = currig - n0rig -n3rig;
//
//             if((n3rig > n34_max) || (n4rig>n34_max))
//               continue;
//               for (int curlef = 0; curlef <= nmax; curlef++)
//                 for(int n0lef = 0; n0lef <= curlef; n0lef++)
//                   for(int llef = -n0lef; llef <= n0lef; llef+=2)
//                   //for(int n2lef = 0; n2lef <= curlef-n1lef; ++n2lef)
//                     for(int n3lef = 0; n3lef <= curlef-n0lef; ++n3lef)
//
//                     {
//                       if (abs(llef)>J)
//                         continue;
//                       if (abs(llef-lrig) > 1)
//                         continue;
//
//                       // if((n0lef+llef) %2 !=0)
//                       //   continue;
//
//                        int n4lef = curlef - n0lef -n3lef;
//
//                        if((n3lef > n34_max) || (n4lef>n34_max))
//                          continue;
//
//                        int i = tab_vib[n0lef][n3lef][n4lef][llef+n0lef];
//                        int j = tab_vib[n0rig][n3rig][n4rig][lrig+n0rig];
//                        // if (j > i)
//                        //  continue;
//
//                         double kin0 = lam1/2.0*(1/2.0*(n0lef+llef)*del(n0lef,n0rig)*del(llef,lrig)+(n0lef-llef+2)/2.0*del(n0lef,n0rig)*del(llef,lrig) + sqrt((n0lef+llef)/2.0)*sqrt((n0lef-llef)/2.0)*del(llef,lrig)*del(n0lef-2,n0rig) + sqrt((n0lef+llef+2)/2.0)*sqrt((n0lef-llef+2)/2.0)*del(llef,lrig)*del(n0lef+2,n0rig))*del(n3lef,n3rig)*del(n4lef,n4rig);
//                         double kin3 =  -lam3/4*(sqrt(n3lef*(n3lef-1))*del(n3lef,n3rig+2) - (n3rig+1)*del(n3lef,n3rig)-n3rig*del(n3lef,n3rig)+sqrt((n3lef+1)*(n3lef+2))*del(n3lef+2,n3rig))*del(n0lef,n0rig)*del(n4lef,n4rig);
//                         double kin4 =  -lam4/4*(sqrt(n4lef*(n4lef-1))*del(n4lef,n4rig+2) - (n4rig+1)*del(n4lef,n4rig)-n4rig*del(n4lef,n4rig)+sqrt((n4lef+1)*(n4lef+2))*del(n4lef+2,n4rig))*del(n0lef,n0rig)*del(n3lef,n3rig);
//                        Ham(i,j) += del(llef,lrig)*(kin0 + kin3+kin4);
//
//                        double zetst = 1.0;//zet[3][0][0];
//                        // double res = zetst*zetst*( lam1/lam4*(n4lef+1/2.0)*( (n0lef+1) + (n0lef+1) ) +
//                        //                          1.0/2.0*((n0lef-1) - (n0lef+1) ) -
//                        //                        (+1/2)*(-(n0lef-1)+(n0lef+1) ) +
//                        //                      lam4/lam1*(n4lef+1/2.0)*((n0lef+1) + (n0lef+1) )   );
//                       // double res1 = lam1/lam4*( (n4lef+1/2.0)*( (n0lef+1) + (n0lef+1) )*del(n4rig,n4lef)*del(n0lef,n0rig)+
//                       //                         (n4lef+1/2.0)*( -sqrt((n0lef+llef+2)*(n0lef-llef+2)) )*del(n4rig,n4lef)*del(n0lef+2,n0rig) +
//                       //                       (n4lef+1/2.0)*( -sqrt((n0lef+llef+2)*(n0lef-llef+2)) )*del(n4rig,n4lef)*del(n0lef+2,n0rig))
//                       // //  Hvibrot(i,j) += res*del(llef,lrig)*del(n3lef,n3rig)*del(n4rig,n4lef)*del(n0lef,n0rig);
//                             double res1 =   lam1/lam4*  (  (n4lef+1/2.0)*del(n4rig,n4lef) + 1/2.0*sqrt((n4lef+1)*(n4lef+2))*del(n4lef+2,n4rig) + 1/2.0*sqrt(n4lef*(n4lef-1))*del(n4lef,n4rig+2)  )*
//                             (2*(n0lef+1)*del(n0lef,n0rig)+ sqrt((n0lef+llef+2)*(n0lef-llef+2))*del(n0lef+2,n0rig) + sqrt(n0lef*n0lef-llef*llef)*del(n0lef,n0rig+2) );
//
//
//                             //first (-) = i^2
//                             double res2 = -(  1/2.0*del(n4rig,n4lef) - 1/2.0*sqrt((n4lef+1)*(n4lef+2))*del(n4lef+2,n4rig) + 1/2.0*sqrt(n4lef*(n4lef-1))*del(n4lef,n4rig+2)  )*
//                             (((llef-1) -(llef+1)  )*del(n0lef,n0rig)+ sqrt((n0lef+llef+2)*(n0lef-llef+2))*del(n0lef+2,n0rig) -  sqrt(n0lef*n0lef-llef*llef)*del(n0lef,n0rig+2) );
//
//                             //first (-) = i^2
//                             double res3 = -(  -1/2.0*del(n4rig,n4lef) - 1/2.0*sqrt((n4lef+1)*(n4lef+2))*del(n4lef+2,n4rig) + 1/2.0*sqrt(n4lef*(n4lef-1))*del(n4lef,n4rig+2)  )*
//                             ((-(llef-1) +(llef+1)  )*del(n0lef,n0rig)+ sqrt((n0lef+llef+2)*(n0lef-llef+2))*del(n0lef+2,n0rig) -  sqrt(n0lef*n0lef-llef*llef)*del(n0lef,n0rig+2) );
//
//                             double res4 = lam4/lam1*  (  (n4lef+1/2.0)*del(n4rig,n4lef) - 1/2.0*sqrt((n4lef+1)*(n4lef+2))*del(n4lef+2,n4rig) - 1/2.0*sqrt(n4lef*(n4lef-1))*del(n4lef,n4rig+2)  )*
//                             (2*(n0lef+1)*del(n0lef,n0rig)- sqrt((n0lef+llef+2)*(n0lef-llef+2))*del(n0lef+2,n0rig) - sqrt(n0lef*n0lef-llef*llef)*del(n0lef,n0rig+2) );
//
//                             Hvibrot(i,j) += zetst*zetst*(res1-res2-res3+res4)*del(llef,lrig)*del(n3lef,n3rig);
//
//                             for (int q0 = 0; q0 < Nquad; ++q0)
//                             {
//                                 double Q0 = sqrt(pts_lag[q0]);
//                                 for (int q3 = 0; q3 < Nquad; ++q3)
//                                 {
//                                     double Q3 = pts[q3];
//                                 for (int q4 = 0; q4 < Nquad; ++q4)
//                                 {
//                                     double Q4 = pts[q4];
//
//                                     double Q1 = Q0*sin(0);
//                                     double Q2 = Q0*cos(0);
//
//                                     Eigen::VectorXd Qvec(4);
//                                   //  Qvec << 0,0,0,0,0,Q1/sqrt(lam1),Q2/sqrt(lam2), Q3/sqrt(lam3), Q4/sqrt(lam4);
//                                       Qvec << Q1/sqrt(lam1),Q2/sqrt(lam2), Q3/sqrt(lam3), Q4/sqrt(lam4);
//                                     // Eigen::VectorXd Qvec(9);
//                                     // Qvec << 0,0,0,0,0,Q1/sqrt(lam1),Q2/sqrt(lam2), Q3/sqrt(lam3), Q4/sqrt(lam4);
//                                     Eigen::VectorXd Xvec;
//                                     Xvec = evec*Qvec;
//                                     // double pot = V(Xvec(0)/sqrt(m[0])+X[0], Xvec(1)/sqrt(m[0])+Y[0], Xvec(2)/sqrt(m[0])+Z[0],
//                                     //               Xvec(3)/sqrt(m[1])+X[1], Xvec(4)/sqrt(m[1])+Y[1], Xvec(5)/sqrt(m[1])+Z[1],
//                                     //               Xvec(6)/sqrt(m[2])+X[2], Xvec(7)/sqrt(m[2])+Y[2], Xvec(8)/sqrt(m[2])+Z[2]);
//                                     double mu = mulin(Xvec,Q1/sqrt(lam1),Q2/sqrt(lam2), Q3/sqrt(lam3), Q4/sqrt(lam4),evec );
//                                     double pot = potvec[pot_idx(q0,q3,q4)];
//                                     Vmat(i,j) += del(llef,lrig)*pol0[q0][lag_ind(n0lef,llef)]*pol0[q0][lag_ind(n0rig,lrig)]*pol3[q3][n3lef]*pol3[q3][n3rig]*pol4[q4][n4lef]*pol4[q4][n4rig]*wei_lag[q0]*wei[q3]*wei[q4]*(pot)*M_PI;
//                                     Mumat(i,j) += del(llef,lrig)*pol0[q0][lag_ind(n0lef,llef)]*pol0[q0][lag_ind(n0rig,lrig)]*pol3[q3][n3lef]*pol3[q3][n3rig]*pol4[q4][n4lef]*pol4[q4][n4rig]*wei_lag[q0]*wei[q3]*wei[q4]*(mu)*M_PI;
//
//                                 }
//                               }
//                             }
//
//
//
//
//
//   }
//   }
//   }
//
//   std::ofstream potfile("vmat.txt");
//
//   for(int i =0; i < Dim; ++i)
//     for (int j = 0; j < Dim;++j)
//     {
//       if (fabs(Vmat(i,j)) > 1E-5)
//       {
//         potfile << i+1 << " " << j+1 << " "<<std::setprecision(12)<<Vmat(i,j)<<"\n";
//       }
//     }
//
//   potfile.close();
//       // std::cout << "result: " << res <<std::endl;
//
//
//   std::cout << "Hamiltonian matrix dimension: " << Dim << "\n";
//
//   Hamiltonian = (Ham+Vmat +1/4.0*Mumat*Hvibrot);
//
//   delete [] pol0_der;
//   delete [] pol0_secder;
// //  delete [] pol2;
//   //delete [] pol3;
//   delete [] pol3_der;
//   delete [] pol3_secder;
//
// //  delete [] pol4;
//   delete [] pol4_der;
//   delete [] pol4_secder;
//
//   for (int i = 0; i < Nquad; ++i)
//   {
//     delete [] pol0[i];
//     delete [] pol3[i];
//     delete [] pol4[i];
//   }
//   delete [] pol0;
//   delete [] pol3;
//   delete [] pol4;
//
//   delete [] potvec;
//
//   //std::cout << Mumat*Hvibrot-(Mumat*Hvibrot).transpose() << std::endl;
//
// }

void VarSolver::fill_hamiltonian()
{
  if (!masses_filled)
    throw std::runtime_error("Fill masses first!");
  if (!mm12_filled)
    throw std::runtime_error("Fill Mm12 matrix first!");

  if (!zet_filled)
    throw std::runtime_error("Fill zeta matrix first!");

  double lam1 = lambdas[0];
  double lam2 = lambdas[1];
  double lam3 = lambdas[2];
  double lam4 = lambdas[3];
  //double xarray = new double[9];
  Eigen::MatrixXd Ham = Eigen::MatrixXd::Zero(Dim,Dim);
  Eigen::MatrixXd Hrot_re = Eigen::MatrixXd::Zero(Dim,Dim);
  Eigen::MatrixXd Hrot_im = Eigen::MatrixXd::Zero(Dim,Dim);
  //Hamiltonian_full =
  //double res = 0;
  for (int q0 = 0; q0 < Nquad; ++q0)
  {
      double Q0 = sqrt(pts_lag[q0]);
      fill_lag(Q0, pol0,pol0_der, pol0_secder);
    //  fill_lag(Q0, pol0);
    //  fill_lag_der(Q0, pol0_der);
    //  fill_lag_secder(Q0, pol0_secder);


      std:: cout << "q1: " << q0 << "\n";
  //  for (int q2 = 0; q2 < Nquad; ++q2)
      for (int q3 = 0; q3 < Nquad; ++q3)
      {

          double Q3 = pts[q3];
        // fill_wf(Q3, pol3);
        // fill_wf_der(Q3, pol3_der);
        // fill_wf_secder(Q3, pol3_secder);
        fill_wf(Q3, pol3, pol3_der, pol3_secder);
      for (int q4 = 0; q4 < Nquad; ++q4)
      {



      //  double Q2 = pts[q2];

        double Q4 = pts[q4];

        //fill_wf(Q2, pol2);

        // fill_wf(Q4, pol4);
        // fill_wf_der(Q4, pol4_der);
        // fill_wf_secder(Q4, pol4_secder);
        fill_wf(Q4, pol4, pol4_der, pol4_secder);



        double Q1 = Q0*sin(0);
        double Q2 = Q0*cos(0);
        Eigen::VectorXd Qvec(4);
        Qvec << Q1/sqrt(lam1),Q2/sqrt(lam2), Q3/sqrt(lam3), Q4/sqrt(lam4);
        Eigen::VectorXd Xvec;
        Xvec = evec*Qvec;
        //double pot = 0;
        double sval[5];
        Xvec = Mm12*Xvec+X0_vec;
        compute_s (Xvec.data(), sval);
        double pot_val = pot(sval);
        // double pot = V(Xvec(0)/sqrt(m[0])+X[0], Xvec(1)/sqrt(m[0])+Y[0], Xvec(2)/sqrt(m[0])+Z[0],
        //               Xvec(3)/sqrt(m[1])+X[1], Xvec(4)/sqrt(m[1])+Y[1], Xvec(5)/sqrt(m[1])+Z[1],
        //               Xvec(6)/sqrt(m[2])+X[2], Xvec(7)/sqrt(m[2])+Y[2], Xvec(8)/sqrt(m[2])+Z[2]);



      //  double U = 0.0;

        double mu = mulin(Q1/sqrt(lam1),Q2/sqrt(lam2), Q3/sqrt(lam3), Q4/sqrt(lam4));



        for (int currig = 0; currig <= nmax; ++currig)
          for (int n0rig = 0; n0rig <= currig; ++n0rig)
            for(int lrig = -n0rig; lrig <= n0rig; lrig+=2)
              for (int n3rig = 0; n3rig <= currig - n0rig; ++n3rig)
              {
                if (abs(lrig)>J)
                  continue;
                int n4rig = currig - n0rig -n3rig;

                if((n3rig > n34_max) || (n4rig>n34_max))
                  continue;



                double pixxyy_dl0 = pi2_dl0(Q0,Q3,Q4, n0rig, lrig, n3rig, n4rig);
                double pix_dlm1_re, pix_dlm1_im,pix_dlp1_re, pix_dlp1_im,piy_dlm1_re, piy_dlm1_im,piy_dlp1_re, piy_dlp1_im;
                double pixx_dlm2_re, pixx_dlm2_im,pixx_dlp2_re, pixx_dlp2_im,piyy_dlm2_re, piyy_dlm2_im,piyy_dlp2_re, piyy_dlp2_im;
                pix_dlm1(pix_dlm1_re, pix_dlm1_im, Q0,Q3,Q4, n0rig, lrig, n3rig, n4rig );
                pix_dlp1(pix_dlp1_re, pix_dlp1_im, Q0,Q3,Q4, n0rig, lrig, n3rig, n4rig );
                piy_dlm1(piy_dlm1_re, piy_dlm1_im, Q0,Q3,Q4, n0rig, lrig, n3rig, n4rig );
                piy_dlp1(piy_dlp1_re, piy_dlp1_im, Q0,Q3,Q4, n0rig, lrig, n3rig, n4rig );



                for (int curlef = 0; curlef <= nmax; curlef++)
                  for(int n0lef = 0; n0lef <= curlef; ++n0lef)
                    for(int llef = -n0lef; llef <= n0lef; llef+=2)
                    //for(int n2lef = 0; n2lef <= curlef-n1lef; ++n2lef)
                      for(int n3lef = 0; n3lef <= curlef-n0lef; ++n3lef)

                      {
                        if (abs(llef)>J)
                          continue;
                        if (abs(llef-lrig) > 1)
                          continue;

                         int n4lef = curlef - n0lef -n3lef;

                         if((n3lef > n34_max) || (n4lef>n34_max))
                           continue;

                         int i = tab_vib[tab_vib_idx(n0lef,n3lef,n4lef,llef+n0lef)];
                         int j = tab_vib[tab_vib_idx(n0rig,n3rig,n4rig,lrig+n0rig)];
                         if (j > i)
                          continue;

                          double Pix_re = Pix(llef,lrig);
                          double Piy_im = Piy(llef,lrig);
                          double Pixx_re = Pixx(llef,lrig);
                          double Piyy_re = Piyy(llef,lrig);

                        //  double rot_re = Pixx_re + Piyy_re + pixxyy;
                          double rot_re =  pixxyy_dl0*del(llef,lrig);
                          double rot_im = 0.0;


//
// if ((llef+1) == lrig)
// {
//   rot_re +=  piy_dlm1_re;
//   rot_im +=  piy_dlm1_im;
// }
// if ((llef-1) == lrig)
// {
//   rot_re +=  piy_dlp1_re;
//   rot_im +=  piy_dlp1_im;
// }
// //
// if ((llef+1) == lrig)
// {
//   rot_re +=  0.0;
//   rot_im +=  Pix_im;
// }
// if ((llef-1) == lrig)
// {
//   rot_re +=  0.0;
//   rot_im +=  Pix_im;
// }

                          if ((llef+1) == lrig)
                          {
                            //rot_re += -2.0*pix_dlm1_re*Pix_re
                            rot_re -=  -2.0*piy_dlm1_im*Piy_im + 2.0*pix_dlm1_re*Pix_re;
                            rot_im -=  2.0*pix_dlm1_im*Pix_re + 2.0*piy_dlm1_re*Piy_im;
                          }
                          if ((llef-1) == lrig)
                          {
                            rot_re -=  -2.0*piy_dlp1_im*Piy_im + 2.0*pix_dlp1_re*Pix_re;
                            rot_im -=  2.0*pix_dlp1_im*Pix_re + 2.0*piy_dlp1_re*Piy_im;
                          }


                          //P_i = -i*\sqrt{lam_i}*\frac{\partial}{\partial Q_i}
                         if ((q0==0) && (q3==0) && (q4 ==0))
                         {
                          // double kin1 =  -lam1/4*(sqrt(n1rig*(n1rig-1))*del(n1lef,n1rig-2) - (n1rig+1)*del(n1lef,n1rig)-n1rig*del(n1lef,n1rig)+sqrt((n1rig+1)*(n1rig+2))*del(n1lef,n1rig+2))*del(n2lef,n2rig)*del(n3lef,n3rig)*del(n4lef,n4rig);
                           //double kin2 =  -lam2/4*(sqrt(n2rig*(n2rig-1))*del(n2lef,n2rig-2) - (n2rig+1)*del(n2lef,n2rig)-n2rig*del(n2lef,n2rig)+sqrt((n2rig+1)*(n2rig+2))*del(n2lef,n2rig+2))*del(n1lef,n1rig)*del(n3lef,n3rig)*del(n4lef,n4rig);

                           double kin0 = lam1/2.0*(1/2.0*(n0lef+llef)*del(n0lef,n0rig)*del(llef,lrig)+(n0lef-llef+2)/2.0*del(n0lef,n0rig)*del(llef,lrig) + sqrt((n0rig+lrig)/2.0)*sqrt((n0rig-lrig)/2.0)*del(llef,lrig)*del(n0lef+2,n0rig) + sqrt((n0rig+lrig+2)/2.0)*sqrt((n0rig-lrig+2)/2.0)*del(llef,lrig)*del(n0lef-2,n0rig))*del(n3lef,n3rig)*del(n4lef,n4rig);
                           double kin3 =  -lam3/4*(sqrt(n3rig*(n3rig-1))*del(n3lef,n3rig-2) - (n3rig+1)*del(n3lef,n3rig)-n3rig*del(n3lef,n3rig)+sqrt((n3rig+1)*(n3rig+2))*del(n3lef,n3rig+2))*del(n0lef,n0rig)*del(n4lef,n4rig);
                          double kin4 =  -lam4/4*(sqrt(n4rig*(n4rig-1))*del(n4lef,n4rig-2) - (n4rig+1)*del(n4lef,n4rig)-n4rig*del(n4lef,n4rig)+sqrt((n4rig+1)*(n4rig+2))*del(n4lef,n4rig+2))*del(n0lef,n0rig)*del(n3lef,n3rig);


                          Ham(i,j) += del(llef,lrig)*(kin0 + kin3+kin4);
                        //Ham(i,j) += del(llef,lrig)*(kin0 + kin3);

                        }

                      //  Ham(i,j) +=-del(llef,lrig)*1.0/2.0*lam4*pol0[lag_ind(n0lef,llef)]*pol0[lag_ind(n0rig,lrig)]*pol3[n3lef]*pol3[n3rig]*pol4[n4lef]*pol4_secder[n4rig]*wei_lag[q0]*wei[q3]*wei[q4]*M_PI;;
                         //std::cout << i << " " << j << "\n";
                        // std::cout << n1lef << " " << n2lef << " " << n3lef << "\n";


                       //Pix^2+Piy^2
                       Hrot_re(i,j) +=  pol0[lag_ind(n0lef,llef)]*pol0[lag_ind(n0rig,lrig)]*pol3[n3lef]*pol3[n3rig]*pol4[n4lef]*pol4[n4rig]*wei_lag[q0]*wei[q3]*wei[q4]*(1/2.0*(Pixx_re + Piyy_re)*mu)*M_PI;;
                       //pix^2+piy^2 -2*(Pix*pix+Piy*piy)
                       Hrot_re(i,j) += 1.0/2.0*mu*rot_re*wei_lag[q0]*wei[q3]*wei[q4]*pol0[lag_ind(n0lef,llef)]*pol3[n3lef]*pol4[n4lef]*1/2.0;

                       Hrot_im(i,j) += 1/2.0*mu*rot_im*wei_lag[q0]*wei[q3]*wei[q4]*pol0[lag_ind(n0lef,llef)]*pol3[n3lef]*pol4[n4lef]*1/2.0;
                       Ham(i,j) += del(llef,lrig)*pol0[lag_ind(n0lef,llef)]*pol0[lag_ind(n0rig,lrig)]*pol3[n3lef]*pol3[n3rig]*pol4[n4lef]*pol4[n4rig]*wei_lag[q0]*wei[q3]*wei[q4]*(pot_val)*M_PI;
                      // Piz_mat(i,j) += piz_val* wei[q1]*wei[q2]*wei[q3]*wei[q4]*pol1[n1lef+2]*pol2[n2lef+2]*pol3[n3lef+2]*pol4[n4lef+2];

                        // Hrot(i,j) +=  1/2.0*muten(1,1)*(-1/4.0*pixxx(n1rig,n2rig,n3rig,pol1,pol2,pol3))*wei[q1]*wei[q2]*wei[q3]*pol1[n1lef]*pol2[n2lef]*pol3[n3lef];
                      }
          }

        }
      }

      // res += pol3[10]*pol3[10]*wei[q0];
      // std::cout << std::scientific << pol0[lag_ind(22,0)] <<std::endl;

    }

    //delete [] xarray;
    Hamiltonian_full = Ham+Hrot_re;
}

double VarSolver::mulin( double qq1, double qq2, double qq3, double qq4 )
{
  //double t5,t7,t11,t21,t25;

  //double I0 = 0.7385233635433291033824945e5;
  double I0 = m[0]*Z[0]*Z[0] + m[1]*Z[1]*Z[1]+m[2]*Z[2]*Z[2];
  // t5 = qq2 * qq2;
  // t7 = qq4 * qq4;
  // t11 = qq1 * qq1;
  // t21 = qq3 * qq3;
  // t25 = 0.2310416607557373010179999e-17 * qq1 + 0.8178386727341119478034002e-17 * qq2 + 0.6255050799999999999999999e-14 * qq4 + 0.1058579445801688866157176e4 * qq3 + 0.5968818480950025175449127e-40 * t5 + 0.3491521880488099255985865e-34 * t7 + I0 + 0.3372415075403443046611281e-40 * qq1 * qq2 + 0.4763582389507596975712964e-41 * t11 + 0.4365126522568433857932602e-20 * qq1 * qq3 + 0.2579314028378363880763709e-37 * qq1 * qq4 + 0.1545162577976832206901634e-19 * qq2 * qq3 + 0.9130226793875112048020382e-37 * qq2 * qq4 + 0.1000000000000000000000001e1 * t21 + 0.1181782023977027542937910e-16 * qq3 * qq4;

  double qqarr [4] = {qq1,qq2,qq3,qq4};
  double res = 0.0;
  for (int i = 1; i <= 4; i++) {
    for (int I = 1; I <=3; ++I)
      res += sqrt(m[I-1])*Z[I-1]*evec(cor(3,I)-1,i-1)*qqarr[i-1];
        //res += m[I-1]*Z[I-1]*evec(cor(3,I)-1,4+i)*qqarr[i-1];
  }
  return 1.0/( (I0+res)*(I0+res)/I0 );

//  return 1.0/(res);

}

double VarSolver::pi2_dl0(double qq0, double qq3, double qq4, int n0, int l, int n3, int n4)
{

  double dfd0 = pol3[n3]*pol4[n4]*pol0_der[lag_ind(n0,l)];
  double dfd0d0 = pol3[n3]*pol4[n4]*pol0_secder[lag_ind(n0,l)];
  double hh = pol0[lag_ind(n0,l)]*pol3[n3]*pol4[n4];
  double dfd3d3 = pol0[lag_ind(n0,l)]*pol3_secder[n3]*pol4[n4];
  double dfd3d4 = pol0[lag_ind(n0,l)]*pol3_der[n3]*pol4_der[n4];
  double dfd4d4 = pol0[lag_ind(n0,l)]*pol3[n3]*pol4_secder[n4];
  double dfd0d3 = pol0_der[lag_ind(n0,l)]*pol3_der[n3]*pol4[n4];
  double dfd0d4 = pol0_der[lag_ind(n0,l)]*pol3[n3]*pol4_der[n4];
  double dfd3 = pol0[lag_ind(n0,l)]*pol3_der[n3]*pol4[n4];
  double dfd4 = pol0[lag_ind(n0,l)]*pol3[n3]*pol4_der[n4];

//  double zetprod, minsing = -1.0;//p^2 = -1* h^2*diff...
  double finval = 0.0;
  double minus_2pi = -2.0*M_PI;


        //xx:
        //030+030, 030+300, 300+300, 300+030
      //  if ((k==0)&& (p==3) && (m==0) && (n==3))
        {
          finval += dfd4d4*qq0*qq0/2.0*zet[0][3][0]*zet[0][3][0];
        }
        //if ((k==0)&& (p==3) && (m==3) && (n==0))
        {
          finval += (qq0/2.0*dfd0 + dfd0d4*qq4*qq0/2.0)*zet[0][3][0]*zet[3][0][0];
        }
      //  if ((k==3)&& (p==0) && (m==3) && (n==0))
        {
          finval += ( dfd0d0*qq4*qq4/2.0 - l*l/2.0/qq0/qq0*hh*qq4*qq4 + dfd0*qq4*qq4/2.0/qq0)*zet[3][0][0]*zet[3][0][0];
        }
      //  if ((k==3)&& (p==0) && (m==0) && (n==3))
        {
          finval += (qq4*dfd4 + dfd0d4*qq4*qq0/2.0)*zet[3][0][0]*zet[0][3][0];
        }

        return minus_2pi*finval*2.0;


  }

  //llef = lrig+1
  void VarSolver::pix_dlp1(double & Re,double & Im, double qq0, double qq3, double qq4, int n0, int l, int n3, int n4 )
  {

  //  double dfd3 = pol0[lag_ind(n0,l)]*pol3_der[n3]*pol4[n4];
    double dfd4 = pol0[lag_ind(n0,l)]*pol3[n3]*pol4_der[n4];
    double dfd0 = pol3[n3]*pol4[n4]*pol0_der[lag_ind(n0,l)];
    double hh = pol0[lag_ind(n0,l)]*pol3[n3]*pol4[n4];

    double t1;
    t1 = qq0 * qq0;
    Re = (0.5928783638484331759343003e1 * t1 * dfd4 - 0.1664692962823733752077279e1 * dfd0 * qq0 * qq4 + 0.1664692962823733752077279e1 * hh * l * qq4) / qq0;;

    Im = 0.0;
  }

  //llef = lrig-1
  void VarSolver::pix_dlm1(double & Re,double & Im, double qq0, double qq3, double qq4, int n0, int l, int n3, int n4)
  {

  //  double dfd3 = pol0[lag_ind(n0,l)]*pol3_der[n3]*pol4[n4];
    double dfd4 = pol0[lag_ind(n0,l)]*pol3[n3]*pol4_der[n4];
    double dfd0 = pol3[n3]*pol4[n4]*pol0_der[lag_ind(n0,l)];
    double hh = pol0[lag_ind(n0,l)]*pol3[n3]*pol4[n4];

    double t1;
    t1 = qq0 * qq0;
    Re = (-0.5928783638484331759343003e1 * dfd4 * t1 + 0.1664692962823733752077279e1 * dfd0 * qq0 * qq4 + 0.1664692962823733752077279e1 * hh * l * qq4) / qq0;

    Im = 0.0;
  }

  void VarSolver::piy_dlm1(double & Re,double & Im, double qq0, double qq3, double qq4, int n0, int l, int n3, int n4 )
  {
    //double dfd3 = pol0[lag_ind(n0,l)]*pol3_der[n3]*pol4[n4];
    double dfd4 = pol0[lag_ind(n0,l)]*pol3[n3]*pol4_der[n4];
    double dfd0 = pol3[n3]*pol4[n4]*pol0_der[lag_ind(n0,l)];
    double hh = pol0[lag_ind(n0,l)]*pol3[n3]*pol4[n4];

    double t1;
    t1 = qq0 * qq0;
    Im = (-0.5928783638484331759343006e1 * t1 * dfd4 + 0.1664692962823733752077281e1 * qq0 * qq4 * dfd0 + 0.1664692962823733752077281e1 * hh * l * qq4) / qq0;

    Re = 0.0;

  }

  void VarSolver::piy_dlp1(double & Re,double & Im, double qq0, double qq3, double qq4, int n0, int l, int n3, int n4 )
  {
    double dfd3 = pol0[lag_ind(n0,l)]*pol3_der[n3]*pol4[n4];
    double dfd4 = pol0[lag_ind(n0,l)]*pol3[n3]*pol4_der[n4];
    double dfd0 = pol3[n3]*pol4[n4]*pol0_der[lag_ind(n0,l)];
    double hh = pol0[lag_ind(n0,l)]*pol3[n3]*pol4[n4];

    double t1;
    t1 = qq0 * qq0;
    Im = (-0.5928783638484331759343006e1 * t1 * dfd4 + 0.1664692962823733752077281e1 * qq0 * qq4 * dfd0 - 0.1664692962823733752077281e1 * hh * l * qq4) / qq0;

    Re = 0.0;

  }

  //Real
  double VarSolver::Pix (int llef, int lrig)
  {
      if( (llef + 1) == lrig)
      {
        return 1/2.0*sqrt((J+lrig)*(J-lrig+1));
      }
      if ((llef-1) == lrig)
      {
        return 1/2.0*sqrt((J-lrig)*(J+lrig+1));
      }
      return 0.0;
  }

  //Imag
  double VarSolver::Piy (int llef, int lrig)
  {
      if( (llef + 1) == lrig)
      {
        return -1/2.0*sqrt((J+lrig)*(J-lrig+1));
      }
      if ((llef-1) == lrig)
      {
        return 1/2.0*sqrt((J-lrig)*(J+lrig+1));
      }
      return 0.0;
  }

  double VarSolver::Pixx (int llef, int lrig)
  {
    if (llef == lrig) {
      return 1/2.0*(J*(J+1) - llef*llef);
  //  return 1;
    }
    if( (llef + 2) == lrig)
    {
      return -1/4.0*sqrt((J-lrig+1)*(J-lrig+2)*(J+lrig)*(J+lrig-1));
    }
    if ((llef-2) == lrig)
    {
      return -1/4.0*sqrt((J+lrig+1)*(J+lrig+2)*(J-lrig)*(J-lrig-1));
    }
    return 0.0;
  }

  double VarSolver::Piyy (int llef, int lrig)
  {
    if (llef == lrig) {
      return 1/2.0*(J*(J+1) - llef*llef);
    //return 1;
    }
    if( (llef + 2) == lrig)
    {
      return 1/4.0*sqrt((J-lrig+1)*(J-lrig+2)*(J+lrig)*(J+lrig-1));
    }
    if ((llef-2) == lrig)
    {
      return 1/4.0*sqrt((J+lrig+1)*(J+lrig+2)*(J-lrig)*(J-lrig-1));
    }
    return 0.0;
  }

  void VarSolver::fill_wf(double x, double * wfarr, double * wfder, double * wfsecder)
  {

    // wfarr[0] = 0.0; //-2
    // wfarr[1] = 0.0; //-1;

    wfarr[0] = 1.0/constants::pi4;
    wfder[0] = -x/constants::pi4;
    wfsecder[0] = -0.751125544464943e0 + 0.751125544464943e0 * x * x;

    if (nmax == 0)
      return;
    wfarr[1] = sqrt(2.0)*x/constants::pi4;
    wfder[1] = sqrt(2.0)/constants::pi4*(1-x*x);
    wfsecder[1] = 0.10622519320272e1 * x*x*x - 0.318675579608161e1 * x;
    for (int i = 1; i < nmax; ++i)
    {
      wfarr[i+1] = (2.0*x/sqrt(2.0*(i+1))*wfarr[i] - i/sqrt(i*(i+1))*wfarr[i-1]);
      wfder[i+1] = sqrt(2.0*(i+1))*wfarr[i] - x*wfarr[i+1];
      wfsecder[i+1] =  sqrt(2.0*(i+1))*wfder[i]  - wfarr[i+1] - x*wfder[i+1];
    }
  }

  void VarSolver::fill_lag(double x, double * wfarr, double *wfder,double * wfsecder)
  {
    double xl = 1;
    double xlm1 = 0.0;
    double xlm2 = 0.0;
    for (int l = 0; l <= nmax; l++) {

      if (l == 1)
        xlm1 = 1.0;
      if (l==2)
        xlm2 = 1.0;


      double zerconst = 1.0/constants::pi2*sqrt(1.0/gsl_sf_fact(l));
      double secconst = 1.0/constants::pi2*sqrt(1.0/gsl_sf_fact(l+1));
      wfarr[lag_ind(l,l)] = xl*zerconst;
      wfarr[lag_ind(l,-l)] =  wfarr[lag_ind(l,l)];
      wfder[lag_ind(l,l)] = zerconst*(l*xlm1 - x*xl);
      wfder[lag_ind(l,-l)] =  wfder[lag_ind(l,l)];
    //  wfsecder[lag_ind(l,l)] = l*(l-1)*xlm2/pi2*sqrt(1.0/gsl_sf_fact(l)) - wfarr[lag_ind(l,l)] + x*x*wfarr[lag_ind(l,l)];
      wfsecder[lag_ind(l,l)] = zerconst*(l*(l-1)*xlm2 + x*x*xl - (2*l+1)*xl );
      wfsecder[lag_ind(l,-l)] =  wfsecder[lag_ind(l,l)];
      if (l == nmax)
        continue;
      wfarr[lag_ind(l+2,l)] = (1+l-x*x)*xl*secconst;
      wfarr[lag_ind(l+2,-l)] =  wfarr[lag_ind(l+2,l)];
      wfder[lag_ind(l+2,l)] = l*(1+l-x*x)*xlm1/constants::pi2*sqrt(1.0/gsl_sf_fact(l+1)) + (-2*x)*xl/constants::pi2*sqrt(1.0/gsl_sf_fact(l+1)) - x*wfarr[lag_ind(l+2,l)];
      wfder[lag_ind(l+2,-l)] =  wfder[lag_ind(l+2,l)];
      wfsecder[lag_ind(l+2,l)] = -secconst*(l*(1-l*l)*xlm2 - 3.0*(l+2)*x*x*xl+3*(l+1)*(l+1)*xl+x*x*x*x*xl);
      wfsecder[lag_ind(l+2,-l)] =  wfsecder[lag_ind(l+2,l)];

      for(int n = l+2; n < nmax; n+=2)
      {
        if ((n+2) > nmax)
          continue;
        int k = (n-l)/2;
        wfarr[lag_ind(n+2,l)] =  sqrt(k+1)/sqrt(l+k+1)*(2*k+1+l-x*x)/(k+1.0)*wfarr[lag_ind(n,l)]-(k+l)/(k+1.0)*wfarr[lag_ind(n-2,l)]*sqrt(k*(k+1))/sqrt((l+k)*(l+k+1));
        wfarr[lag_ind(n+2,-l)] =  wfarr[lag_ind(n+2,l)];

        wfder[lag_ind(n+2,l)] =  sqrt(k+1)/sqrt(l+k+1)*(-2*x)/(k+1.0)*wfarr[lag_ind(n,l)] + sqrt(k+1)/sqrt(l+k+1)*(2*k+1+l-x*x)/(k+1.0)*wfder[lag_ind(n,l)]-(k+l)/(k+1.0)*wfder[lag_ind(n-2,l)]*sqrt(k*(k+1))/sqrt((l+k)*(l+k+1));
        wfder[lag_ind(n+2,-l)] =  wfder[lag_ind(n+2,l)];

        wfsecder[lag_ind(n+2,l)] =  sqrt(k+1)/sqrt(l+k+1)*(-2)/(k+1.0)*wfarr[lag_ind(n,l)] + sqrt(k+1)/sqrt(l+k+1)*(-2*x)/(k+1.0)*wfder[lag_ind(n,l)] +
                                    sqrt(k+1)/sqrt(l+k+1)*(2*k+1+l-x*x)/(k+1.0)*wfsecder[lag_ind(n,l)] + sqrt(k+1)/sqrt(l+k+1)*(-2*x)/(k+1.0)*wfder[lag_ind(n,l)]-
                                    (k+l)/(k+1.0)*wfsecder[lag_ind(n-2,l)]*sqrt(k*(k+1))/sqrt((l+k)*(l+k+1));
        wfsecder[lag_ind(n+2,-l)] =  wfsecder[lag_ind(n+2,l)];
      }


      xl *= x;
      xlm1 *= x;
      xlm2 *= x;
    }
  }


  void VarSolver::diagonalize()
  {
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(Hamiltonian_full);
    std::cout << "The eigenvalues of A are:" << std::endl << std::setprecision(9) << es.eigenvalues()*constants::htocm  << std::endl;
  }

  void VarSolver::print_normal_coords(std::ofstream & file,constants::matr_type type)
  {
    if (!harmonic_solved)
      throw std::runtime_error("Solve the harmonic problem first!");
    Eigen::Matrix4d lammat = Eigen::Matrix4d::Zero();
    lammat(0,0) = sqrt(lambdas[0]);
    lammat(1,1) = sqrt(lambdas[1]);
    lammat(2,2) = sqrt(lambdas[2]);
    lammat(3,3) = sqrt(lambdas[3]);

    file << std::setprecision(10) << std::scientific;
    for (int i =0; i < 9; ++i)
    {
      for (int j = 0; j < 4; ++j)
      {
        if (type==constants::TRANSFORMED)
          file << (Mm12*evec*lammat.inverse())(i,j) << " ";
        else if (type == constants::PURE)
          file << evec(i,j) << " ";
        else if (type == constants::FOR_AAT)
          file << (2*Mm12*evec*lammat)(i,j) << " ";
        else
          throw std::runtime_error("Unknown type!");
      }
      file << "\n";
    }

  }

  void VarSolver::prepare_for_business()
  {
    if (setup_pot != nullptr)
      setup_pot();
    else
      std::cout << "REMINDER! Your potential is probably not set up!" << std::endl;

    if (extpot == nullptr)
      throw std::runtime_error("There is no potential! Check input");
    //Eigen::MatrixXd evec(9,4);
    harmonic_solver();
  //  std::cout << asin(0.0) << std::endl;
    fill_zet();
    fill_tabvib();
    print_mem();
    check_tabvib();
  }
