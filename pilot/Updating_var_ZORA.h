#pragma once

#include <iostream>
#include <vector>
#include <cstdlib>
#include <numeric>




#include "../api/Printer"
#include "../api/Timer"
#include "../api/Gaussians"
#include "../api/Plotter"
#include "../api/MWFunctions"
#include "../api/MWOperators"


void Update_SCF_Variables2(MultiResolutionAnalysis<3> &MRA, ABGVOperator<3> &D, std::vector<mrcpp::CompFunction<3>> &Psi_2c, CompFunction<3> &V_electron_nucleus ,std::vector<std::vector<mrcpp::CompFunction<3> *>> &Nabla_Psi_2c,CompFunction<3> &K_tree, std::vector<mrcpp::CompFunction<3> *> &Nabla_K_tree, CompFunction<3> &K_inverse, CompFunction<3> &V, CompFunction<3> &J, FunctionTree<3> &K_real, mrcpp::PoissonOperator &P){
    // Updating the Nabla_Psi_2c
    std::cout << "--------------UPDATING GRADIENT--------------" << '\n';
    Nabla_Psi_2c[0] = mrcpp::gradient(D,Psi_2c[0]);
    Nabla_Psi_2c[1] = mrcpp::gradient(D,Psi_2c[1]);
    // Print  values
    std::cout << "Nabla_Psi_2c[0] = " << Nabla_Psi_2c[0][0]->getSquareNorm() << "   " <<  Nabla_Psi_2c[0][1]->getSquareNorm() << "   " <<  Nabla_Psi_2c[0][2]->getSquareNorm() << '\n';
    std::cout << "Nabla_Psi_2c[1] = " << Nabla_Psi_2c[1][0]->getSquareNorm() << "   " <<  Nabla_Psi_2c[1][1]->getSquareNorm() << "   " <<  Nabla_Psi_2c[1][2]->getSquareNorm() << '\n';
    std::cout << '\n';

    // Same format, print the pointers
    std::cout << "Nabla_Psi_2c[0][0] = " << Nabla_Psi_2c[0][0] << "   " << Nabla_Psi_2c[0][1] << "   " << Nabla_Psi_2c[0][2] << '\n';
    std::cout << "Nabla_Psi_2c[1][0] = " << Nabla_Psi_2c[1][0] << "   " << Nabla_Psi_2c[1][1] << "   " << Nabla_Psi_2c[1][2] << '\n';

    std::cout << "------------------- DONE -------------------" << '\n' << '\n';
   
    // Create the Rho_t and Rho_b, then Rho = Rho_t + Rho_b
    std::cout << "--------------UPDATING RHO------------------" << '\n';
    CompFunction<3> Rho_t(MRA);
    CompFunction<3> Rho_b(MRA);
    CompFunction<3> Rho(MRA);


    std::cout << "Computing density for TOP component..." << '\n';
    make_density_local(Rho_t, Psi_2c[0] , MRA, building_precision);
    std::cout << "Norm of Rho_top = " << Rho_t.getSquareNorm() << '\n';

  
    std::cout << "Computing density for BOTTOM component..." << '\n';
    make_density_local(Rho_b, Psi_2c[1] , MRA, building_precision);
    std::cout << "Norm of Rho_bottom = " << Rho_b.getSquareNorm() << '\n' << '\n';
  

    mrcpp::add(Rho, 1.0, Rho_t, 1.0, Rho_b, building_precision, false);
    std::cout << "Rho is real = " << Rho.isreal() << '\n';
    std::cout << "Rho is complex = " << Rho.iscomplex() << '\n';
    std::cout << "Norm of Rho = " << Rho.getSquareNorm() << '\n';

    if ( num_cycle !=0  && (Rho_b.getSquareNorm() == 0.0 || Rho_t.getSquareNorm() == 0.0)) {
        std::cerr << "Error: Rho_t or Rho_b has zero norm!" << '\n';
        exit(1);
        return;
    }
    std::cout << "---------------- DONE -------------------" << '\n' << '\n';

    // Calculate the J operator
    std::cout << "--------------UPDATING J-----------------" << '\n';
    
    apply(building_precision, J, P, Rho);
    std::cout << "Applied the Poisson operator to the density" << '\n';
    J.rescale(ComplexDouble(2.0,0.0));  


    std::cout << "J is real = " << J.isreal() << '\n';
    std::cout << "J is complex = " << J.iscomplex() << '\n';
    std::cout << "---------------- DONE -------------------" << '\n' << '\n';


    // Compute the full potential
    std::cout << "--------------UPDATING V------------------" << '\n';
    mrcpp::add(V, 1.0, V_electron_nucleus,1.0, J, building_precision, false);
    std::cout << "V norm = " << V.getSquareNorm() << '\n';
    std::cout << "---------------- DONE -------------------" << '\n' << '\n';




    // Compute the new K tree
    std::cout << "-------------UPDATING K_inv---------------" << '\n';
    double constant = 2.0 * m * c * c;
    double prefactor = 1/constant;


    // Copy V in K_inverse
    
    K_inverse.add(ComplexDouble(-constant,0.0),V ); // K_inv = -2mc^2 +V
    std::cout << "K = -mc^2 *V, Norm =" << K_inverse.getSquareNorm() << '\n'; 
    K_inverse.rescale(-prefactor);                  // K_inv = 1 - V/2mc^2

    std::cout << "K_inverse norm = " << K_inverse.getSquareNorm() << '\n';
    std::cout << "K_inverse is real = " << K_inverse.isreal() << '\n';
    std::cout << "K_inverse is complex = " << K_inverse.iscomplex() << '\n';
    std::cout << "K_inv Function Tree:" << '\n';
    std::cout << *(K_inverse.CompD[0]) << '\n';
    std::cout << "---------------- DONE -------------------" << '\n' << '\n';





    // Compute the inverse of K_inverse (that is just the Kappa operator):
    std::cout << "--------------UPDATING K_tree------------------" << '\n';
    
    std::function<double(double)> Reciprocal = [] (const double &x) -> double {
        return  (1.0 / x) ;
    };
    
    if (K_inverse.CompD[0] == nullptr) {
        std::cerr << "Error: K_inverse.CompD[0] is null!" << '\n';
        return;
    }
    
    std::cout << "Now about to map the K_inverse function to K_real" << '\n';
    FunctionTree<3> K_real_inp(MRA);

    


    K_inverse.CompD[0]->deep_copy(&K_real_inp);
    
    map(building_precision, K_real, K_real_inp ,Reciprocal);

    std::cout << "Putting the K_real in the tree" << '\n';
    K_tree.setReal(&K_real,0);

    std::cout << "K_real norm = " << K_real.getSquareNorm() << '\n';
    std::cout << "K_tree norm = " << K_tree.getSquareNorm() << '\n';
    std::cout << "K_tree is real = " << K_tree.isreal() << '\n';
    std::cout << "K_tree is complex = " << K_tree.iscomplex() << '\n';

    std::cout << "---------------- DONE -------------------" << '\n' << '\n';




    // Compute the Nabla_K_tree
    std::cout << "--------------UPDATING Nabla_K_tree------------------" << '\n';
    Nabla_K_tree = mrcpp::gradient(D,K_tree);

    // Print the pointers of the Nabla_K_tree
    std::cout << "Nabla_K_tree[0] = " << Nabla_K_tree[0] << '\n';
    std::cout << "Nabla_K_tree[1] = " << Nabla_K_tree[1] << '\n';
    std::cout << "Nabla_K_tree[2] = " << Nabla_K_tree[2] << '\n' <<'\n';
    // Now we print the norms
    std::cout << "Nabla_K_tree[0] norm = " << Nabla_K_tree[0]->getSquareNorm() << '\n';
    std::cout << "Nabla_K_tree[1] norm = " << Nabla_K_tree[1]->getSquareNorm() << '\n';
    std::cout << "Nabla_K_tree[2] norm = " << Nabla_K_tree[2]->getSquareNorm() << '\n' <<'\n';

    std::cout << "Nabla_K_tree[0] is real = " << Nabla_K_tree[0]->isreal() << '\n';
    std::cout << "Nabla_K_tree[0] is complex = " << Nabla_K_tree[0]->iscomplex() << '\n';

    std::cout << "-------------- DONE UPDATING SCF VARIABLES -------------------" << '\n' << '\n';

}