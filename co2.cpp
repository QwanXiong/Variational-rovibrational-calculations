
#include "VarSolver.h"
#include <iostream>
#include <fstream>

extern "C"
{

extern void __attribute__((stdcall)) co2potlongrange_(double &,double &,double &,double &);
extern void __attribute__((stdcall)) setupco2longrange_( void);

}


int main(int argc, char const *argv[]) {

  VarSolver solver(10,10,0);
  solver.set_mass_z(11.996709*constants::em,15.990525*constants::em,15.990525*constants::em,0.0, -2.192119941, 2.192119941);
  solver.setup_pot = &setupco2longrange_;
  solver.extpot = & co2potlongrange_;


  solver.prepare_for_business();
  solver.fill_hamiltonian();
  solver.diagonalize();

    std::ofstream normfile("normal_coords.txt");
    solver.print_normal_coords(normfile,constants::FOR_AAT);

    normfile.close();

  return 0;
}
