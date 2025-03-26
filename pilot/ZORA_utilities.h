//void compute_V_phi_with_spin_orbit_coupling(){
//}
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






// FOR NOW I PLAY IT SAFE, A SIMPLE ADDITION MAY BE JUST TO ADD IN THE ARGUMENTS THE ABGV OPERATOR, NOT TO REDEFINE IT EVERY TIME



/*
*
*   CHANGE IN TERM 1 AND 2: THE FUNCTION DOT ALREADY CONJUGATES THE COMPFUNCTION BRA!!!!!!!!! UNCONJUGATE THEM AS IT IS DONE IN THE FUNCTION
*
*/


// IN principle, one could skip the step of recomputing all the gradient of the \Psi * K bu simpli adding \Nabla \Psi  K + \Psi  \Nabla K, but this would be longer, tho requiring less memory
ComplexDouble compute_Term1_T_ZORA(MultiResolutionAnalysis<3> &MRA, std::vector<std::vector<mrcpp::CompFunction<3>*>> &Nabla_Psi_2c, mrcpp::CompFunction<3> &K_tree,  std::vector<mrcpp::CompFunction<3>> &Psi_2c){
    // Nabla(\Psi K) = \Nabla \Psi  K + \Psi  \Nabla K

    bool debug = false;


    CompFunction<3> Psi_t_K(MRA);
    CompFunction<3> Psi_b_K(MRA);

    std::vector<mrcpp::CompFunction<3> *> Nabla_Psi_t = Nabla_Psi_2c[0];
    std::vector<mrcpp::CompFunction<3> *> Nabla_Psi_b = Nabla_Psi_2c[1];

    // Compute Psi_top * K
    mrcpp::multiply(Psi_t_K, Psi_2c[0], K_tree,building_precision, false, false, false);
    // Compute Psi_bottom * K
    mrcpp::multiply(Psi_b_K, Psi_2c[1], K_tree,building_precision, false, false, false);

    if (debug){
        std::cout << '\n';
        std::cout << "Psi_t_K = " << Psi_t_K.getSquareNorm() << '\n';
        std::cout << "Psi_b_K = " << Psi_b_K.getSquareNorm() << '\n' << '\n';
        std::cout << "Nabla_Psi_t[0] = " << Nabla_Psi_t[0]->getSquareNorm() << '\n';
        std::cout << "Nabla_Psi_t[1] = " << Nabla_Psi_t[1]->getSquareNorm() << '\n';
        std::cout << "Nabla_Psi_t[2] = " << Nabla_Psi_t[2]->getSquareNorm() << '\n';
        std::cout << "Nabla_Psi_b[0] = " << Nabla_Psi_b[0]->getSquareNorm() << '\n';
        std::cout << "Nabla_Psi_b[1] = " << Nabla_Psi_b[1]->getSquareNorm() << '\n';
        std::cout << "Nabla_Psi_b[2] = " << Nabla_Psi_b[2]->getSquareNorm() << '\n';
        std::cout << '\n';
        std::cout << "Nabla_Psi_2c[0][0] = " << Nabla_Psi_2c[0][0]->getSquareNorm() << '\n';
        std::cout << "Nabla_Psi_2c[0][1] = " << Nabla_Psi_2c[0][1]->getSquareNorm() << '\n';
        std::cout << "Nabla_Psi_2c[0][2] = " << Nabla_Psi_2c[0][2]->getSquareNorm() << '\n';
        std::cout << "Nabla_Psi_2c[1][0] = " << Nabla_Psi_2c[1][0]->getSquareNorm() << '\n';
        std::cout << "Nabla_Psi_2c[1][1] = " << Nabla_Psi_2c[1][1]->getSquareNorm() << '\n';
        std::cout << "Nabla_Psi_2c[1][2] = " << Nabla_Psi_2c[1][2]->getSquareNorm() << '\n';
        std::cout << '\n';
        std::cout << '\n';
    }

    // Operatr ABGV
    mrcpp::ABGVOperator<3> D(MRA, 0.0, 0.0);

    // Gradient of Psi_top * K
    auto Nabla_Psi_t_K = mrcpp::gradient(D,Psi_t_K);
    // Gradient of Psi_bottom * K
    auto Nabla_Psi_b_K = mrcpp::gradient(D,Psi_b_K);

    // REMEMBER THE MINUS SIGN!!! (it comes from the integration by parts)
    ComplexDouble Top_contribute = -1.0*(dot(*Nabla_Psi_t_K[0],*Nabla_Psi_t[0]) + dot(*Nabla_Psi_t_K[1],*Nabla_Psi_t[1]) + dot(*Nabla_Psi_t_K[2],*Nabla_Psi_t[2]));
    ComplexDouble Bottom_contribute = -1.0*(dot(*Nabla_Psi_b_K[0],*Nabla_Psi_b[0]) + dot(*Nabla_Psi_b_K[1],*Nabla_Psi_b[1]) + dot(*Nabla_Psi_b_K[2],*Nabla_Psi_b[2]));
    
    if (debug){
        std::cout << "Top_contribute = " << Top_contribute << '\n';
        std::cout << "Bottom_contribute = " << Bottom_contribute << '\n';
    }

    return Top_contribute + Bottom_contribute;
}



ComplexDouble compute_Term2_T_ZORA(MultiResolutionAnalysis<3> &MRA, std::vector<std::vector<mrcpp::CompFunction<3>*>> &Nabla_Psi_2c, mrcpp::CompFunction<3> &K_tree, std::vector<mrcpp::CompFunction<3> *> &Nabla_K_tree,  std::vector<mrcpp::CompFunction<3>> &Psi_2c){
    std::vector<mrcpp::CompFunction<3>> Psi_t_Nabla_K; 
    std::vector<mrcpp::CompFunction<3>> Psi_b_Nabla_K;

    CompFunction<3> Psi_Nabla_k_tmp(MRA);
    CompFunction<3> Psi_element(MRA);
    CompFunction<3> Nabla_K_element(MRA);
    
    // Operatr ABGV
    mrcpp::ABGVOperator<3> D(MRA, 0.0, 0.0);


    // Gradient of K
    auto Nabla_K = Nabla_K_tree;

    // Gradient of the 2 components of Psi
    // PAY ATTENTION!!! The gradient obtained is actually a pointer, so later on we'll have to dereference it
    auto Nabla_Psi_t = Nabla_Psi_2c[0];
    auto Nabla_Psi_b = Nabla_Psi_2c[1];

    // Compute Psi_top * Nabla_K
    for (int i=0;i<3;i++){
        mrcpp::multiply(Psi_Nabla_k_tmp, Psi_2c[0], *Nabla_K[i],building_precision, false, false, false);   // <---- Here we dereference the pointer
        Psi_t_Nabla_K.push_back(Psi_Nabla_k_tmp);
    }
    
    // Compute Psi_bottom * Nabla_K
    for (int i=0;i<3;i++){
        mrcpp::multiply(Psi_Nabla_k_tmp, Psi_2c[1], *Nabla_K[i],building_precision, false, false, true);   // <---- Here we dereference the pointer
        Psi_b_Nabla_K.push_back(Psi_Nabla_k_tmp);
    }

    std::vector<ComplexDouble> Kinetic_term_1;


    ComplexDouble Top_contribution  = dot(Psi_t_Nabla_K[0],*Nabla_Psi_t[0]) + dot(Psi_t_Nabla_K[1],*Nabla_Psi_t[1]) + dot(Psi_t_Nabla_K[2],*Nabla_Psi_t[2]);
    ComplexDouble Bottom_contribution  = dot(Psi_b_Nabla_K[0],*Nabla_Psi_b[0]) + dot(Psi_b_Nabla_K[1],*Nabla_Psi_b[1]) + dot(Psi_b_Nabla_K[2],*Nabla_Psi_b[2]);


    return Top_contribution + Bottom_contribution;
}









void compute_rotor( std::vector<mrcpp::CompFunction<3>*> &Rotor , mrcpp::ABGVOperator<3> &D,  std::vector<mrcpp::CompFunction<3>*> &K_Nabla_Psi, MultiResolutionAnalysis<3> &MRA){
    Eigen::Matrix<int, 3, 2> Curl_coef;
    Curl_coef << 1, 2,
                 2, 0,
                 0, 1;
                      
    CompFunction<3> Deriv_1(MRA);
    CompFunction<3> Deriv_2(MRA);


    // For synthax:
    //void mrcpp::apply<3>(mrcpp::CompFunction<3> &out, mrcpp::DerivativeOperator<3> &oper, mrcpp::CompFunction<3> &inp, int dir, ComplexDouble (*metric)[4] = (ComplexDouble (*)[4])nullptr)
    // mrcpp::apply(deriv_Psi1, D, Psi_in[s], j); // Derivative of Psi with respect to j
    for (int i=0; i<3; i++){
        mrcpp::apply<3>(Deriv_1, D, *(K_Nabla_Psi[Curl_coef(i,1)]), Curl_coef(i,0));
        mrcpp::apply<3>(Deriv_2, D, *(K_Nabla_Psi[Curl_coef(i,0)]), Curl_coef(i,1));
        mrcpp::add<3>(*Rotor[i], 1.0, Deriv_1, -1.0, Deriv_2, building_precision, false);
    }


}


void compute_sigma_cdot_spinor(MultiResolutionAnalysis<3> &MRA  , std::vector<mrcpp::CompFunction<3>*> &Curl_top, std::vector<mrcpp::CompFunction<3>*> &Curl_bottom, CompFunction<3> &SOC_Psi_t, CompFunction<3> &SOC_Psi_b){
    // Define the sigma matrices in a vector:
    std::vector<Eigen::Matrix2cd> sigma(3);
    sigma[0] << 0, 1,
                1, 0; // Pauli matrix sigma_x
    sigma[1] << 0, -std::complex<double>(0, 1),
                std::complex<double>(0, 1), 0; // Pauli matrix sigma_y
    sigma[2] << 1, 0,
                0, -1; // Pauli matrix sigma_z


    // I'll devide the xyz components of the top and bottoom components of the spinor in 3 spinors, one for each direction
    std::vector<mrcpp::CompFunction<3>> Spinor_component_x(2);
    std::vector<mrcpp::CompFunction<3>> Spinor_component_y(2);
    std::vector<mrcpp::CompFunction<3>> Spinor_component_z(2);

    Spinor_component_x[0] = *Curl_top[0];
    Spinor_component_x[1] = *Curl_bottom[0];

    Spinor_component_y[0] = *Curl_top[1];
    Spinor_component_y[1] = *Curl_bottom[1];

    Spinor_component_z[0] = *Curl_top[2];
    Spinor_component_z[1] = *Curl_bottom[2];


    // Now i compute the scalar product between the rotor and the sigma matrix vector
    // I'll do one by one
    // Compute SOC_i as the matrix-vector product of sigma_i and Spinor_component_i
    std::vector<mrcpp::CompFunction<3>*> SOC_x(2, new mrcpp::CompFunction<3>(MRA));
    std::vector<mrcpp::CompFunction<3>*> SOC_y(2, new mrcpp::CompFunction<3>(MRA));
    std::vector<mrcpp::CompFunction<3>*> SOC_z(2, new mrcpp::CompFunction<3>(MRA));

    // In general, for each matrix vector(spinor) multiplication we need 4 scalar multiplication and 2 additions

    // Sigma_x * Spinor_component_x
    mrcpp::add(*SOC_x[0], sigma[0](0,0), Spinor_component_x[0], sigma[0](0,1), Spinor_component_x[1], building_precision, false);
    mrcpp::add(*SOC_x[1], sigma[0](1,0), Spinor_component_x[0], sigma[0](1,1), Spinor_component_x[1], building_precision, false);

    // Sigma_y * Spinor_component_y
    mrcpp::add(*SOC_y[0], sigma[1](0,0), Spinor_component_y[0], sigma[1](0,1), Spinor_component_y[1], building_precision, false);
    mrcpp::add(*SOC_y[1], sigma[1](1,0), Spinor_component_y[0], sigma[1](1,1), Spinor_component_y[1], building_precision, false);

    // Sigma_z * Spinor_component_z
    mrcpp::add(*SOC_z[0], sigma[2](0,0), Spinor_component_z[0], sigma[2](0,1), Spinor_component_z[1], building_precision, false);
    mrcpp::add(*SOC_z[1], sigma[2](1,0), Spinor_component_z[0], sigma[2](1,1), Spinor_component_z[1], building_precision, false);


    // Now i sum the y and z components
    std::vector<mrcpp::CompFunction<3>> SOC_yz(2); // This is the sum of the y and z components, temporary, to be used in the final sum
    mrcpp::add(SOC_yz[0], 1.0, *SOC_y[0], 1.0, *SOC_z[0],building_precision, false);
    mrcpp::add(SOC_yz[1], 1.0, *SOC_y[1], 1.0, *SOC_z[1],building_precision, false);

    // Finally i sum the x component to the yz component to the total SOC
    mrcpp::add(SOC_Psi_t, 1.0, *SOC_x[0], 1.0, SOC_yz[0],building_precision, false);
    mrcpp::add(SOC_Psi_b, 1.0, *SOC_x[1], 1.0, SOC_yz[1],building_precision, false);

}





ComplexDouble compute_Term3_T_ZORA(MultiResolutionAnalysis<3> &MRA, std::vector<std::vector<mrcpp::CompFunction<3>*>> &Nabla_Psi_2c, mrcpp::CompFunction<3> &K_tree, std::vector<mrcpp::CompFunction<3> *> &Nabla_K_tree,  std::vector<mrcpp::CompFunction<3>> &Psi_2c){
    // Compute the product K * (\Nabla \Psi) for top and bottom components
    bool debug = false;


    std::vector<mrcpp::CompFunction<3>*> K_Nabla_Psi_top(3, new mrcpp::CompFunction<3>(MRA));
    std::vector<mrcpp::CompFunction<3>*> K_Nabla_Psi_bottom(3, new mrcpp::CompFunction<3>(MRA));

    std::vector<mrcpp::CompFunction<3> *> Nabla_Psi_top = Nabla_Psi_2c[0];
    std::vector<mrcpp::CompFunction<3> *> Nabla_Psi_bottom = Nabla_Psi_2c[1];
    if (debug){
        std::cout << '\n';
        std::cout << "Nabla_Psi_top[0] = " << Nabla_Psi_top[0]->getSquareNorm() << '\n';
        std::cout << "Nabla_Psi_top[1] = " << Nabla_Psi_top[1]->getSquareNorm() << '\n';
        std::cout << "Nabla_Psi_top[2] = " << Nabla_Psi_top[2]->getSquareNorm() << '\n';   
        std::cout << "Nabla_Psi_bottom[0] = " << Nabla_Psi_bottom[0]->getSquareNorm() << '\n';
        std::cout << "Nabla_Psi_bottom[1] = " << Nabla_Psi_bottom[1]->getSquareNorm() << '\n';
        std::cout << "Nabla_Psi_bottom[2] = " << Nabla_Psi_bottom[2]->getSquareNorm() << '\n';
        std::cout << '\n';
    }

    // For synthax:
    // void mrcpp::multiply<3>(mrcpp::CompFunction<3> &out, mrcpp::CompFunction<3> inp_a, mrcpp::CompFunction<3> inp_b, double prec, bool absPrec, bool useMaxNorms, bool conjugate)
    std::cout << '\n' << "Now computing K * Nabla_Psi" << '\n';
    //mrcpp::multiply(tmp, K_tree, *Nabla_Psi_top[1] ,building_precision, false, false, false);   

    for (int i = 0; i<2; i++){
        mrcpp::multiply( *K_Nabla_Psi_top[i] , K_tree, *Nabla_Psi_top[i] ,building_precision, false, false, false);   
        mrcpp::multiply( *K_Nabla_Psi_bottom[i], K_tree, *Nabla_Psi_bottom[i] ,building_precision, false, false, false);   
    }

    if (debug){
        std::cout << "K_Nabla_Psi_top[0] = " << K_Nabla_Psi_top[0]->getSquareNorm() << '\n';
        std::cout << "K_Nabla_Psi_top[1] = " << K_Nabla_Psi_top[1]->getSquareNorm() << '\n';
        std::cout << "K_Nabla_Psi_top[2] = " << K_Nabla_Psi_top[2]->getSquareNorm() << '\n';
        std::cout << "K_Nabla_Psi_bottom[0] = " << K_Nabla_Psi_bottom[0]->getSquareNorm() << '\n';
        std::cout << "K_Nabla_Psi_bottom[1] = " << K_Nabla_Psi_bottom[1]->getSquareNorm() << '\n';
        std::cout << "K_Nabla_Psi_bottom[2] = " << K_Nabla_Psi_bottom[2]->getSquareNorm() << '\n';
    }

    
    // Compute the rotor of K * (\Nabla \Psi) for top and bottom components
    std::vector<mrcpp::CompFunction<3>*> rotor_K_Nabla_Psi_top(3, new mrcpp::CompFunction<3>(MRA));
    std::vector<mrcpp::CompFunction<3>*> rotor_K_Nabla_Psi_bottom(3, new mrcpp::CompFunction<3>(MRA));

    mrcpp::ABGVOperator<3> D(MRA, 0.0, 0.0);
    std::cout << "Now computing rotor of the top component" << '\n';
    compute_rotor(rotor_K_Nabla_Psi_top, D, K_Nabla_Psi_top, MRA);
    std::cout << "Now computing rotor of the bottom component" << '\n';
    compute_rotor(rotor_K_Nabla_Psi_bottom, D, K_Nabla_Psi_bottom, MRA);
    // Compute the scalar product between the rotor and the sigma matrix vector
    CompFunction<3> SOC_Psi_t(MRA);
    CompFunction<3> SOC_Psi_b(MRA);
    std::cout << "Now computing sigma_cdot_spinor" << '\n';
    compute_sigma_cdot_spinor(MRA, rotor_K_Nabla_Psi_top, rotor_K_Nabla_Psi_bottom, SOC_Psi_t, SOC_Psi_b);

    // Now i do the innetr produc with the bra
    std::cout << "Now computing the inner product top" << '\n';
    ComplexDouble Top_contribution = dot(Psi_2c[0], SOC_Psi_t);
    std::cout << "Now computing the inner product bottom" << '\n';
    ComplexDouble Bottom_contribution = dot(Psi_2c[1], SOC_Psi_b);

    return ComplexDouble(0,1)*(Top_contribution + Bottom_contribution);
}







ComplexDouble compute_energy_ZORA(MultiResolutionAnalysis<3> &MRA, std::vector<std::vector<mrcpp::CompFunction<3>*>> &Nabla_Psi_2c, mrcpp::CompFunction<3> &K_tree, std::vector<mrcpp::CompFunction<3> *> &Nabla_K_tree,  std::vector<mrcpp::CompFunction<3>> &Psi_2c, CompFunction<3> &V){
    ComplexDouble Term1 = compute_Term1_T_ZORA(MRA, Nabla_Psi_2c, K_tree, Psi_2c);
    ComplexDouble Term2 = compute_Term2_T_ZORA(MRA, Nabla_Psi_2c, K_tree, Nabla_K_tree, Psi_2c);
    ComplexDouble Term3 = compute_Term3_T_ZORA(MRA, Nabla_Psi_2c, K_tree, Nabla_K_tree, Psi_2c);

    // Now we do <Psi | V | Psi>
    mrcpp::CompFunction<3> psi_top__V(MRA);
    mrcpp::CompFunction<3> psi_bottom__V(MRA);

    mrcpp::multiply(psi_top__V, V, Psi_2c[0],building_precision, false, false, false);
    mrcpp::multiply(psi_bottom__V, V, Psi_2c[1],building_precision, false, false, false);

    ComplexDouble potential_energy = dot(Psi_2c[0],psi_top__V) + dot(Psi_2c[1],psi_bottom__V);

    ComplexDouble kinetic_energy = -0.5*m* (Term1 + Term2 + Term3) ;

    return kinetic_energy + potential_energy;



}
    



void compute_term_A_Propagator(MultiResolutionAnalysis<3> &MRA, double E_n, std::vector<mrcpp::CompFunction<3>> &Psi_in, CompFunction<3> &V,CompFunction<3> &term_A_top,CompFunction<3> &term_A_bottom, std::vector<mrcpp::CompFunction<3>> &V_Psi){
    double const_factor = -E_n / (c*c);
    
    mrcpp::multiply(building_precision, term_A_top, 1.0, Psi_in[0], V);
    mrcpp::multiply(building_precision, term_A_bottom, 1.0, Psi_in[1], V);
    V_Psi.push_back(term_A_top);
    V_Psi.push_back(term_A_bottom);

    term_A_top.rescale(const_factor);
    term_A_bottom.rescale(const_factor);
}


void compute_term_B_Propagator(MultiResolutionAnalysis<3> &MRA, std::vector<mrcpp::CompFunction<3>> &Psi_in, std::vector<std::vector<mrcpp::CompFunction<3> *>> &Nabla_Psi_2c, std::vector<mrcpp::CompFunction<3> *> &Nabla_K_tree , CompFunction<3> &K_tree, CompFunction<3> &K_inverse, CompFunction<3> &term_B_top,CompFunction<3> &term_B_bottom){
    std::vector<mrcpp::CompFunction<3>*> Dot_product(2, new mrcpp::CompFunction<3>(MRA));
    for (int i = 0; i < 3; i++){
        mrcpp::multiply(*Dot_product[0], *Nabla_K_tree[i], *Nabla_Psi_2c[0][i],building_precision, false, false, false);
        mrcpp::multiply(*Dot_product[1], *Nabla_K_tree[i], *Nabla_Psi_2c[1][i],building_precision, false, false, false);
    }    


    // Times K INVERSE
    mrcpp::multiply(term_B_top, K_inverse, *Dot_product[0],building_precision, false, false, false);
    mrcpp::multiply(term_B_bottom, K_inverse, *Dot_product[1],building_precision, false, false, false);

}


void compute_term_C_Propagator(MultiResolutionAnalysis<3> &MRA, std::vector<std::vector<mrcpp::CompFunction<3> *>> &Nabla_Psi_2c , CompFunction<3> &K_tree, CompFunction<3> &K_inverse, CompFunction<3> &term_C_top,CompFunction<3> &term_C_bottom){
    // Where I'll put the rotor of K * (\Nabla \Psi)
    debug = true;
    std::vector<mrcpp::CompFunction<3> *>  Curl_t(3, new mrcpp::CompFunction(MRA));
    std::vector<mrcpp::CompFunction<3> *>  Curl_b(3, new mrcpp::CompFunction(MRA));

    // Let's compute K * (\Nabla \Psi)
    std::vector<mrcpp::CompFunction<3> *> K_Nabla_Psi_t(3, new mrcpp::CompFunction(MRA));
    std::vector<mrcpp::CompFunction<3> *> K_Nabla_Psi_b(3, new mrcpp::CompFunction(MRA));
    std::cout << " Correctly defined the variables" << '\n';
    for (int i = 0; i<3; i++){
        mrcpp::multiply(*K_Nabla_Psi_t[i], K_tree, *Nabla_Psi_2c[0][i],building_precision, false, false, false);
        mrcpp::multiply(*K_Nabla_Psi_b[i], K_tree, *Nabla_Psi_2c[1][i],building_precision, false, false, false);
    }
    std::cout << " Correctly computed the product K * Grad(Psi)" << '\n';
    
    ABGVOperator<3> D(MRA, 0.0, 0.0);
    compute_rotor(Curl_t, D, K_Nabla_Psi_t, MRA);
    compute_rotor(Curl_b, D, K_Nabla_Psi_b, MRA);
    std::cout << " Correctly computed the rotor" << '\n';

    // Now we compute the scalar product between the rotor and the sigma matrix vector
    CompFunction<3> SOC_Psi_t(MRA);
    CompFunction<3> SOC_Psi_b(MRA);
    compute_sigma_cdot_spinor(MRA, Curl_t, Curl_b, SOC_Psi_t, SOC_Psi_b);
    std::cout << " Correctly computed the sigma scalar product" << '\n';

    // Now the actual top and bottom components for C are the i * K_inverse * SOC_Psi

    // Times K INVERSE
    mrcpp::multiply(term_C_top, K_inverse, SOC_Psi_t,building_precision, false, false, false);
    mrcpp::multiply(term_C_bottom, K_inverse, SOC_Psi_b,building_precision, false, false, false);
    std::cout << " Correctly computed the term C" << '\n';
    // Times i 
    term_C_top.rescale(ComplexDouble(0,1));
    term_C_bottom.rescale(ComplexDouble(0,1));
    std::cout << " Correctly rescaled the term C" << '\n';
}


void compute_term_D_Propagator(std::vector<mrcpp::CompFunction<3>> &V_Psi , CompFunction<3> &K_inverse, CompFunction<3> &term_D_top,CompFunction<3> &term_D_bottom){
    mrcpp::multiply(building_precision, term_D_top, -2*m, K_inverse, V_Psi[0]);
    mrcpp::multiply(building_precision, term_D_bottom, -2*m, K_inverse, V_Psi[1]);
}






void apply_Helmholtz_ZORA(MultiResolutionAnalysis<3> &MRA, std::vector<mrcpp::CompFunction<3>> &Psi_in, std::vector<std::vector<mrcpp::CompFunction<3> *>> &Nabla_Psi_2c, std::vector<mrcpp::CompFunction<3> *> &Nabla_K_tree,  CompFunction<3> &V, CompFunction<3> &K_tree, CompFunction<3> &K_inverse, double &E_n, std::vector<mrcpp::CompFunction<3>> &Psi_out){
    // Build the Helmholtz operator
    double mu = std::sqrt(-2*E_n);
    mrcpp::HelmholtzOperator Helm(MRA, mu, building_precision);

    // we will apply the Helmholtz operator to the 2 components of the spinor separately
    // Psi_out = Helm(-term_1 + K_inverted * (Kinetic_term_2 + Kinetic_term_3 + 2*m*V)  )

    /*
     * Where the term_1 = (energy/c^2) * V * Psi_in
     * K_inverted = (1 - V / (2*m*c^2) = 1/K_r
     * Kinetic_term_2 =  \Nabla(K) \cdot \Nabla(\Psi_in)
     * Kinetic_term_3 = i * \sigma \cdot (\Nabla \cross K (\Nabla (\Psi_in)))
     *     
     * 
    */

    // Term a) -E^n/c^2 * V * Psi^n
    CompFunction<3> term_A_top(MRA);
    CompFunction<3> term_A_bottom(MRA);

    // Term b) K_invinverted * (\Nabla K \cdot \Nabla \Psi^n)
    CompFunction<3> term_B_top(MRA);
    CompFunction<3> term_B_bottom(MRA);

    // Term c) K_invinverted * (i * \sigma \cdot (\Nabla \cross K (\Nabla (\Psi^n))))
    CompFunction<3> term_C_top(MRA);
    CompFunction<3> term_C_bottom(MRA);

    // Term d) K_invinverted * 2*m*V
    CompFunction<3> term_D_top(MRA);
    CompFunction<3> term_D_bottom(MRA);

    // Compute the term A will also evaluate the potential energy times the wavefunction, we'll need it later in term d.
    std::vector<mrcpp::CompFunction<3>> V_Psi;
    std::cout << "Now computing term A" << '\n';
    compute_term_A_Propagator(MRA, E_n, Psi_in, V, term_A_top, term_A_bottom, V_Psi);

    // Compute the term B
    std::cout << "Now computing term B" << '\n';
    compute_term_B_Propagator(MRA, Psi_in, Nabla_Psi_2c, Nabla_K_tree, K_tree, K_inverse, term_B_top, term_B_bottom);
    
    // Compute the term C
    std::cout << "Now computing term C" << '\n';
    compute_term_C_Propagator(MRA, Nabla_Psi_2c, K_tree, K_inverse, term_C_top, term_C_bottom);

    // Compute the term D
    std::cout << "Now computing term D" << '\n';
    compute_term_D_Propagator(V_Psi, K_inverse, term_D_top, term_D_bottom);

    // Now we sum all the terms, we have a total of 4 terms to sum. We'll take them by pairs and sum them
    CompFunction<3> first_pair_top(MRA);
    CompFunction<3> first_pair_bottom(MRA);
    CompFunction<3> second_pair_top(MRA);
    CompFunction<3> second_pair_bottom(MRA);
    std::cout << "Now computing the additions" << '\n';
    mrcpp::add(first_pair_top, 1.0, term_A_top, 1.0, term_B_top, building_precision, false);
    mrcpp::add(first_pair_bottom, 1.0, term_A_bottom, 1.0, term_B_bottom, building_precision, false);
    mrcpp::add(second_pair_top, 1.0, term_C_top, 1.0, term_D_top, building_precision, false);
    mrcpp::add(second_pair_bottom, 1.0, term_C_bottom, 1.0, term_D_bottom, building_precision, false);

    // Now we sum the two pairs
    std::vector<mrcpp::CompFunction<3>> Psi_to_be_convoluted(2);

    std::cout << "Now computing the final addition" << '\n';
    mrcpp::add(Psi_to_be_convoluted[0], 1.0, first_pair_top, 1.0, second_pair_top, building_precision, false);
    mrcpp::add(Psi_to_be_convoluted[1], 1.0, first_pair_bottom, 1.0, second_pair_bottom, building_precision, false);

    // Now we apply the Helmholtz operator
    std::cout << "Now applying the Helmholtz operator" << '\n';
    mrcpp::apply(building_precision, Psi_out[0], Helm, Psi_to_be_convoluted[0]);
    mrcpp::apply(building_precision, Psi_out[1], Helm, Psi_to_be_convoluted[1]);

}



