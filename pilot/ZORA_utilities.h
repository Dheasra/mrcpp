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

    // < \Nabla (\Psi * K) | \Nabla \Psi > 
    ComplexDouble Top_contribute = dot(*Nabla_Psi_t_K[0],*Nabla_Psi_t[0]) + dot(*Nabla_Psi_t_K[1],*Nabla_Psi_t[1]) + dot(*Nabla_Psi_t_K[2],*Nabla_Psi_t[2]);
    ComplexDouble Bottom_contribute = dot(*Nabla_Psi_b_K[0],*Nabla_Psi_b[0]) + dot(*Nabla_Psi_b_K[1],*Nabla_Psi_b[1]) + dot(*Nabla_Psi_b_K[2],*Nabla_Psi_b[2]);
    
    if (debug){
        std::cout << "Top_contribute = " << Top_contribute << '\n';
        std::cout << "Bottom_contribute = " << Bottom_contribute << '\n';
    }
    // REMEMBER THE MINUS SIGN!!! (it comes from the integration by parts)
    return -(Top_contribute + Bottom_contribute);
}



ComplexDouble compute_Term2_T_ZORA(MultiResolutionAnalysis<3> &MRA, std::vector<std::vector<mrcpp::CompFunction<3>*>> &Nabla_Psi_2c, mrcpp::CompFunction<3> &K_tree, std::vector<mrcpp::CompFunction<3> *> &Nabla_K_tree,  std::vector<mrcpp::CompFunction<3>> &Psi_2c){
    /// NEW_APPROACH: I first do all the \nabla(K) * \nabla(\Psi) and then multiply by the Psi from the left as a bra
    debug=false;
    std::vector<mrcpp::CompFunction<3>> Nabla_K_Nabla_Psi_top(3, MRA);
    std::vector<mrcpp::CompFunction<3>> Nabla_K_Nabla_Psi_bottom(3, MRA);

    // \Naabla(K) * |\Nabla(\Psi)>
    for (int i=0;i<3;i++){
        mrcpp::multiply(Nabla_K_Nabla_Psi_top[i], *Nabla_K_tree[i], *Nabla_Psi_2c[0][i],building_precision, false, false, false);   // <---- Here we dereference the pointer
    }
    for (int i=0;i<3;i++){
        mrcpp::multiply(Nabla_K_Nabla_Psi_bottom[i], *Nabla_K_tree[i], *Nabla_Psi_2c[1][i],building_precision, false, false, false);   // <---- Here we dereference the pointer
    }

    // Debug
    if (debug){
        std::cout << '\n';
        std::cout << "Nabla_K_Nabla_Psi_top[0] = " << Nabla_K_Nabla_Psi_top[0].getSquareNorm() << '\n';
        std::cout << "Nabla_K_Nabla_Psi_top[1] = " << Nabla_K_Nabla_Psi_top[1].getSquareNorm() << '\n';
        std::cout << "Nabla_K_Nabla_Psi_top[2] = " << Nabla_K_Nabla_Psi_top[2].getSquareNorm() << '\n';
        std::cout << "Nabla_K_Nabla_Psi_bottom[0] = " << Nabla_K_Nabla_Psi_bottom[0].getSquareNorm() << '\n';
        std::cout << "Nabla_K_Nabla_Psi_bottom[1] = " << Nabla_K_Nabla_Psi_bottom[1].getSquareNorm() << '\n';
        std::cout << "Nabla_K_Nabla_Psi_bottom[2] = " << Nabla_K_Nabla_Psi_bottom[2].getSquareNorm() << '\n';
        std::cout << '\n';
    }


    ComplexDouble Top_contribution  = dot(Psi_2c[0],Nabla_K_Nabla_Psi_top[0]) + dot(Psi_2c[0],Nabla_K_Nabla_Psi_top[1]) + dot(Psi_2c[0],Nabla_K_Nabla_Psi_top[2]);
    ComplexDouble Bottom_contribution  = dot(Psi_2c[1],Nabla_K_Nabla_Psi_bottom[0]) + dot(Psi_2c[1],Nabla_K_Nabla_Psi_bottom[1]) + dot(Psi_2c[1],Nabla_K_Nabla_Psi_bottom[2]);

    return Top_contribution + Bottom_contribution;
}








/* 
void compute_rotor( std::vector<mrcpp::CompFunction<3>*> &Rotor , mrcpp::ABGVOperator<3> &D,  std::vector<mrcpp::CompFunction<3>*> &K_Nabla_Psi, MultiResolutionAnalysis<3> &MRA){
    Eigen::Matrix<int, 3, 2> Curl_coef;
    //          [A][B]
    Curl_coef << 1, 2,
                 2, 0,
                 0, 1;
                      
    CompFunction<3> Deriv_1(MRA);
    CompFunction<3> Deriv_2(MRA);

    // For synthax:
    //void mrcpp::apply<3>(mrcpp::CompFunction<3> &out, mrcpp::DerivativeOperator<3> &oper, mrcpp::CompFunction<3> &inp, int dir, ComplexDouble (*metric)[4] = (ComplexDouble (*)[4])nullptr)
    // mrcpp::apply(deriv_Psi1, D, Psi_in[s], j); // Derivative of Psi with respect to j
    for (int i=0; i<3; i++){
        //              Term 1     Derivative             K_Nabla_Psi[A]                Direction = [B]
        mrcpp::apply<3>(Deriv_1,      D       ,   *(K_Nabla_Psi[Curl_coef(i,1)])  ,     Curl_coef(i,0));    // \dv{F_A}{B}
        mrcpp::apply<3>(Deriv_2,      D       ,   *(K_Nabla_Psi[Curl_coef(i,0)])  ,     Curl_coef(i,1));    // \dv{F_B}{A}
        mrcpp::add<3>(*Rotor[i], 1.0, Deriv_1, -1.0, Deriv_2, building_precision, false); // \dv{F_A}{B} - \dv{F_B}{A}
    }
} */

void Cross_Product_Compfunct( std::vector<mrcpp::CompFunction<3>*> &Cross , std::vector<mrcpp::CompFunction<3>*> &A, std::vector<mrcpp::CompFunction<3>*> &B, MultiResolutionAnalysis<3> &MRA){
    Eigen::Matrix<int, 3, 2> Cross_coef;
    //          [A][B]
    Cross_coef << 1, 2,
                 2, 0,
                 0, 1;
                      
    CompFunction<3> Term_1(MRA);
    CompFunction<3> Term_2(MRA);

    if (debug){
        std::cout << "Cross_Product_Compfunct" << '\n';
        std::cout << "A[0] = " << A[0]->getSquareNorm() << '\n';
        std::cout << "A[1] = " << A[1]->getSquareNorm() << '\n';
        std::cout << "A[2] = " << A[2]->getSquareNorm() << '\n';
        std::cout << "B[0] = " << B[0]->getSquareNorm() << '\n';
        std::cout << "B[1] = " << B[1]->getSquareNorm() << '\n';
        std::cout << "B[2] = " << B[2]->getSquareNorm() << '\n';    
    }


    for (int i=0; i<3; i++){
        mrcpp::multiply(Term_1, *A[Cross_coef(i,0)], *B[Cross_coef(i,1)],building_precision, false, false, false); // A_i * B_j
        mrcpp::multiply(Term_2, *A[Cross_coef(i,1)], *B[Cross_coef(i,0)],building_precision, false, false, false); // A_j * B_i
        mrcpp::add<3>(*Cross[i], 1.0, Term_1 , -1.0, Term_2, building_precision, false); // \dv{F_A}{B} - \dv{F_B}{A}
    }

}




void compute_sigma_cdot_spinor(MultiResolutionAnalysis<3> &MRA  , std::vector<mrcpp::CompFunction<3>*> &Curl_top, std::vector<mrcpp::CompFunction<3>*> &Curl_bottom, CompFunction<3> &SOC_Psi_t, CompFunction<3> &SOC_Psi_b){
    // Define the sigma matrices in a vector:
    std::vector<Eigen::Matrix2cd> sigma(3);
    sigma[0] << 0, std::complex<double>(1,0),
                std::complex<double>(1,0), 0; // Pauli matrix sigma_x
    sigma[1] << 0, -std::complex<double>(0, 1),
                std::complex<double>(0, 1), 0; // Pauli matrix sigma_y
    sigma[2] << std::complex<double>(1,0), 0,
                0, std::complex<double>(-1,0); // Pauli matrix sigma_z


    // I'll devide the xyz components of the top and bottoom components of the spinor in 3 spinors, one for each direction
    std::vector<mrcpp::CompFunction<3>> Spinor_component_x(2, MRA); // <-- 2 component of X
    std::vector<mrcpp::CompFunction<3>> Spinor_component_y(2, MRA); // <-- 2 component of Y
    std::vector<mrcpp::CompFunction<3>> Spinor_component_z(2, MRA); // <-- 2 component of Z

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

    std::vector<mrcpp::CompFunction<3>> SOC_yz(2, MRA); // This is the sum of the y and z components, temporary, to be used in the final sum
    mrcpp::add(SOC_yz[0], ComplexDouble(0.,-1.) , *Curl_bottom[1]   ,   +1.0, *Curl_top[2],building_precision, false);
    mrcpp::add(SOC_yz[1], ComplexDouble(1.0,0) , *Curl_top[1]       ,   -1.0, *Curl_bottom[2],building_precision, false);

    // Finally i sum the x component to the yz component to the total SOC
    mrcpp::add(SOC_Psi_t, 1.0, *Curl_bottom[0], 1.0, SOC_yz[0],building_precision, false);
    mrcpp::add(SOC_Psi_b, 1.0, *Curl_top[0], 1.0, SOC_yz[1],building_precision, false);

    if (debug){
        std::cout << "> NOW COMPUTING THE SIGMA CDOT SPINOR" << '\n';
        // print the sigma matrices:
        //std::cout << "Sigma_x = [" << sigma[0](0,0) << ", " << sigma[0](0,1) << "]" << '\n';
        //std::cout << "           " << sigma[0](1,0) << ", " << sigma[0](1,1) << "]" << '\n';
        //std::cout << "Sigma_y = [" << sigma[1](0,0) << ", " << sigma[1](0,1) << "]" << '\n';
        //std::cout << "           " << sigma[1](1,0) << ", " << sigma[1](1,1) << "]" << '\n';
        //std::cout << "Sigma_z = [" << sigma[2](0,0) << ", " << sigma[2](0,1) << "]" << '\n';
        //std::cout << "           " << sigma[2](1,0) << ", " << sigma[2](1,1) << "]" << '\n'; 
        
        std::cout << '\n';
        std::cout << "Curl_top[0] = " << Curl_top[0]->getSquareNorm() << '\n';
        std::cout << "Curl_top[1] = " << Curl_top[1]->getSquareNorm() << '\n';
        std::cout << "Curl_top[2] = " << Curl_top[2]->getSquareNorm() << '\n';
        std::cout << "Curl_bottom[0] = " << Curl_bottom[0]->getSquareNorm() << '\n';
        std::cout << "Curl_bottom[1] = " << Curl_bottom[1]->getSquareNorm() << '\n';
        std::cout << "Curl_bottom[2] = " << Curl_bottom[2]->getSquareNorm() << '\n';

        std::cout << '\n';
        //std::cout << "Spinor_component_x[0] = " << Spinor_component_x[0].getSquareNorm() << '\n';
        //std::cout << "Spinor_component_y[0] = " << Spinor_component_y[0].getSquareNorm() << '\n';
        //std::cout << "Spinor_component_z[0] = " << Spinor_component_z[0].getSquareNorm() << '\n';
        //std::cout << "Spinor_component_x[1] = " << Spinor_component_x[1].getSquareNorm() << '\n';
        //std::cout << "Spinor_component_y[1] = " << Spinor_component_y[1].getSquareNorm() << '\n';
        //std::cout << "Spinor_component_z[1] = " << Spinor_component_z[1].getSquareNorm() << '\n';
        std::cout << "SOC_yz[0] = " << SOC_yz[0].getSquareNorm() << '\n';
        std::cout << "SOC_yz[1] = " << SOC_yz[1].getSquareNorm() << '\n';
        //std::cout << "SOC_x[0] = " << SOC_x[0]->getSquareNorm() << '\n';
        //std::cout << "SOC_y[0] = " << SOC_y[0]->getSquareNorm() << '\n';
        //std::cout << "SOC_z[0] = " << SOC_z[0]->getSquareNorm() << '\n';
        //std::cout << "SOC_x[1] = " << SOC_x[1]->getSquareNorm() << '\n';
        //std::cout << "SOC_y[1] = " << SOC_y[1]->getSquareNorm() << '\n';
        //std::cout << "SOC_z[1] = " << SOC_z[1]->getSquareNorm() << '\n';
        std::cout << '\n';
        std::cout << "SOC_Psi_t = " << SOC_Psi_t.getSquareNorm() << '\n';
        std::cout << "SOC_Psi_b = " << SOC_Psi_b.getSquareNorm() << '\n';
    }

}





ComplexDouble compute_Term3_T_ZORA(MultiResolutionAnalysis<3> &MRA, std::vector<std::vector<mrcpp::CompFunction<3>*>> &Nabla_Psi_2c, mrcpp::CompFunction<3> &K_tree, std::vector<mrcpp::CompFunction<3> *> &Nabla_K_tree,  std::vector<mrcpp::CompFunction<3>> &Psi_2c){
    // Compute the product K * (\Nabla \Psi) for top and bottom components
    bool debug = true;

    std::vector<mrcpp::CompFunction<3>*> K_Nabla_Psi_top(3, new mrcpp::CompFunction<3>(MRA));
    std::vector<mrcpp::CompFunction<3>*> K_Nabla_Psi_bottom(3, new mrcpp::CompFunction<3>(MRA));

    std::vector<mrcpp::CompFunction<3> *> Nabla_Psi_top = Nabla_Psi_2c[0];
    std::vector<mrcpp::CompFunction<3> *> Nabla_Psi_bottom = Nabla_Psi_2c[1];


    // For synthax:
    // void mrcpp::multiply<3>(mrcpp::CompFunction<3> &out, mrcpp::CompFunction<3> inp_a, mrcpp::CompFunction<3> inp_b, double prec, bool absPrec, bool useMaxNorms, bool conjugate)
    //std::cout << '\n' << "Now computing K * Nabla_Psi" << '\n';
    //mrcpp::multiply(tmp, K_tree, *Nabla_Psi_top[1] ,building_precision, false, false, false);   
/*
    for (int i = 0; i<2; i++){
        mrcpp::multiply( *K_Nabla_Psi_top[i] , K_tree, *Nabla_Psi_top[i] ,building_precision, false, false, false);   
        mrcpp::multiply( *K_Nabla_Psi_bottom[i], K_tree, *Nabla_Psi_bottom[i] ,building_precision, false, false, false);   
    }
*/
    /* if (debug){
        std::cout << "K_Nabla_Psi_top[0] = " << K_Nabla_Psi_top[0]->getSquareNorm() << '\n';
        std::cout << "K_Nabla_Psi_top[1] = " << K_Nabla_Psi_top[1]->getSquareNorm() << '\n';
        std::cout << "K_Nabla_Psi_top[2] = " << K_Nabla_Psi_top[2]->getSquareNorm() << '\n';
        std::cout << "K_Nabla_Psi_bottom[0] = " << K_Nabla_Psi_bottom[0]->getSquareNorm() << '\n';
        std::cout << "K_Nabla_Psi_bottom[1] = " << K_Nabla_Psi_bottom[1]->getSquareNorm() << '\n';
        std::cout << "K_Nabla_Psi_bottom[2] = " << K_Nabla_Psi_bottom[2]->getSquareNorm() << '\n';
    }
 */
    
    // Compute the rotor of K * (\Nabla \Psi) for top and bottom components
    //std::vector<mrcpp::CompFunction<3>*> rotor_K_Nabla_Psi_top(3, new mrcpp::CompFunction<3>(MRA));
    //std::vector<mrcpp::CompFunction<3>*> rotor_K_Nabla_Psi_bottom(3, new mrcpp::CompFunction<3>(MRA));

    mrcpp::ABGVOperator<3> D(MRA, 0.0, 0.0);

    // Compute the scalar product between the rotor and the sigma matrix vector
    CompFunction<3> SOC_Psi_t(MRA);
    CompFunction<3> SOC_Psi_b(MRA);

    std::vector<mrcpp::CompFunction<3>*> cross_top(3, new mrcpp::CompFunction<3>(MRA));
    std::vector<mrcpp::CompFunction<3>*> cross_bottom(3, new mrcpp::CompFunction<3>(MRA));
    Cross_Product_Compfunct(cross_top, Nabla_K_tree, Nabla_Psi_top, MRA);
    Cross_Product_Compfunct(cross_bottom, Nabla_K_tree, Nabla_Psi_bottom, MRA);

    compute_sigma_cdot_spinor(MRA,cross_top , cross_bottom , SOC_Psi_t, SOC_Psi_b);
    

    // Now i do the innetr produc with the bra
    //std::cout << "Now computing the inner product top" << '\n';
    ComplexDouble Top_contribution = dot(Psi_2c[0], SOC_Psi_t);
    //std::cout << "Now computing the inner product bottom" << '\n';
    ComplexDouble Bottom_contribution = dot(Psi_2c[1], SOC_Psi_b);

    if (debug){

        std::cout << '\n';
        std::cout << "TERM 3 COMPUTATION" << '\n';
        std::cout << "Nabla_Psi_top[0] = " << Nabla_Psi_top[0]->getSquareNorm() << '\n';
        std::cout << "Nabla_Psi_top[1] = " << Nabla_Psi_top[1]->getSquareNorm() << '\n';
        std::cout << "Nabla_Psi_top[2] = " << Nabla_Psi_top[2]->getSquareNorm() << '\n';   
        std::cout << '\n';
        std::cout << "Nabla_Psi_bottom[0] = " << Nabla_Psi_bottom[0]->getSquareNorm() << '\n';
        std::cout << "Nabla_Psi_bottom[1] = " << Nabla_Psi_bottom[1]->getSquareNorm() << '\n';
        std::cout << "Nabla_Psi_bottom[2] = " << Nabla_Psi_bottom[2]->getSquareNorm() << '\n';
        std::cout << '\n';
        std::cout << "K_tree = " << K_tree.getSquareNorm() << '\n';
        std::cout << '\n';
        std::cout << "Nabla_K_tree[0] = " << Nabla_K_tree[0]->getSquareNorm() << '\n';
        std::cout << "Nabla_K_tree[1] = " << Nabla_K_tree[1]->getSquareNorm() << '\n';
        std::cout << "Nabla_K_tree[2] = " << Nabla_K_tree[2]->getSquareNorm() << '\n';
        std::cout << '\n';
        std::cout << "Cross top[0] = " << cross_top[0]->getSquareNorm() << '\t' << "Real: " << cross_top[0]->isreal() << "; Complex: " << cross_top[0]->iscomplex() << '\n';
        std::cout << "Cross top[1] = " << cross_top[1]->getSquareNorm() << '\t' << "Real: " << cross_top[1]->isreal() << "; Complex: " << cross_top[1]->iscomplex() << '\n';
        std::cout << "Cross top[2] = " << cross_top[2]->getSquareNorm() << '\t' << "Real: " << cross_top[2]->isreal() << "; Complex: " << cross_top[2]->iscomplex() << '\n';
        std::cout << "Cross bottom[0] = " << cross_bottom[0]->getSquareNorm() << '\t' << "Real: " << cross_bottom[0]->isreal() << "; Complex: " << cross_bottom[0]->iscomplex() << '\n';
        std::cout << "Cross bottom[1] = " << cross_bottom[1]->getSquareNorm() << '\t' << "Real: " << cross_bottom[1]->isreal() << "; Complex: " << cross_bottom[1]->iscomplex() << '\n';
        std::cout << "Cross bottom[2] = " << cross_bottom[2]->getSquareNorm() << '\t' << "Real: " << cross_bottom[2]->isreal() << "; Complex: " << cross_bottom[2]->iscomplex() << '\n';
        std::cout << '\n';
        std::cout << "SOC_Psi_t = " << SOC_Psi_t.getSquareNorm() << '\n';
        std::cout << "SOC_Psi_b = " << SOC_Psi_b.getSquareNorm() << '\n';
        std::cout << '\n';
        std::cout << "Top_contribution = " << Top_contribution << '\n';
        std::cout << "Bottom_contribution = " << Bottom_contribution << '\n';


        std::cout << '\n';

    }

    return ComplexDouble(0,1)*(Top_contribution + Bottom_contribution);
}







double compute_energy_ZORA(MultiResolutionAnalysis<3> &MRA, std::vector<std::vector<mrcpp::CompFunction<3>*>> &Nabla_Psi_2c, mrcpp::CompFunction<3> &K_tree, std::vector<mrcpp::CompFunction<3> *> &Nabla_K_tree,  std::vector<mrcpp::CompFunction<3>> &Psi_2c, CompFunction<3> &V){
    ComplexDouble Term1 = compute_Term1_T_ZORA(MRA, Nabla_Psi_2c, K_tree, Psi_2c);
    ComplexDouble Term2 = compute_Term2_T_ZORA(MRA, Nabla_Psi_2c, K_tree, Nabla_K_tree, Psi_2c);
    ComplexDouble Term3 = compute_Term3_T_ZORA(MRA, Nabla_Psi_2c, K_tree, Nabla_K_tree, Psi_2c);

    // Now we do <Psi | V | Psi>
    mrcpp::CompFunction<3> psi_top__V(MRA);
    mrcpp::CompFunction<3> psi_bottom__V(MRA);

    mrcpp::multiply(psi_top__V, V, Psi_2c[0],building_precision, false, false, false);
    mrcpp::multiply(psi_bottom__V, V, Psi_2c[1],building_precision, false, false, false);

    ComplexDouble potential_energy = dot(Psi_2c[0],psi_top__V) + dot(Psi_2c[1],psi_bottom__V);

    double kin_prefactor = -1.0/(2.0*m);
    ComplexDouble kinetic_energy = kin_prefactor * (Term1 + Term2 + Term3) ;

    ComplexDouble total_energy = kinetic_energy + potential_energy;
    std::cout << "--------------------------------------------------" << '\n';
    std::cout << "Kinetic energy TERM 1 ="  << kin_prefactor * Term1 << '\n';
    std::cout << "Kinetic energy TERM 2 ="  << kin_prefactor * Term2 << '\n';
    std::cout << "Kinetic energy TERM 3 ="  << kin_prefactor * Term3 << '\n';
    std::cout << "Potential energy ="  << potential_energy << '\n'<< '\n';

    std::cout << "Total energy   = "  << total_energy << '\n';
    std::cout << "Energy modulus = "  << std::abs(total_energy) << '\n';
    std::cout << "--------------------------------------------------" << '\n';


    

    // IF the enery is above zero it stops the program
    if (total_energy.real() > 0){
        std::cout << "Energy is above zero, stopping the program" << '\n';
        exit(1);
    }

    return total_energy.real();

}
    



void compute_term_A_Propagator(MultiResolutionAnalysis<3> &MRA, double E_n, std::vector<mrcpp::CompFunction<3>> &Psi_in, CompFunction<3> &V,CompFunction<3> &term_A_top,CompFunction<3> &term_A_bottom){
    double const_factor = -E_n / (c*c);
    
    std::cout << "--------------------------------------------------" << '\n';
    std::cout << " Input in A:" << '\n';
    std::cout << "  Psi_in[0] = " << Psi_in[0].getSquareNorm() << '\n';
    std::cout << "  Psi_in[1] = " << Psi_in[1].getSquareNorm() << '\n';
    std::cout << "  V = " << V.getSquareNorm() << '\n';


    mrcpp::multiply(building_precision, term_A_top, 1.0, V, Psi_in[0]);
    mrcpp::multiply(building_precision, term_A_bottom, 1.0, V, Psi_in[1]);
    term_A_top.rescale(const_factor);
    term_A_bottom.rescale(const_factor);

}


void compute_term_B_Propagator(MultiResolutionAnalysis<3> &MRA, std::vector<mrcpp::CompFunction<3>> &Psi_in, std::vector<std::vector<mrcpp::CompFunction<3> *>> &Nabla_Psi_2c, std::vector<mrcpp::CompFunction<3> *> &Nabla_K_tree , CompFunction<3> &K_tree, CompFunction<3> &K_inverse, CompFunction<3> &term_B_top,CompFunction<3> &term_B_bottom){
    std::vector<mrcpp::CompFunction<3>> Nabla_K_Nabla_Psi_top(3, MRA);
    std::vector<mrcpp::CompFunction<3>> Nabla_K_Nabla_Psi_bottom(3, MRA);

    // EXPLICIT RUTE
    mrcpp::CompFunction<3> DOT_TOP_x(MRA); // ----> X component for the dot product
    mrcpp::CompFunction<3> DOT_TOP_y(MRA); // ----> Y component for the dot product
    mrcpp::CompFunction<3> DOT_TOP_z(MRA); // ----> Z component for the dot product
    mrcpp::CompFunction<3> DOT_BOTTOM_x(MRA); // ----> X component for the dot product
    mrcpp::CompFunction<3> DOT_BOTTOM_y(MRA); // ----> Y component for the dot product
    mrcpp::CompFunction<3> DOT_BOTTOM_z(MRA); // ----> Z component for the dot product


    mrcpp::multiply(DOT_TOP_x, *Nabla_K_tree[0], *Nabla_Psi_2c[0][0],building_precision, false, false, false);
    mrcpp::multiply(DOT_TOP_y, *Nabla_K_tree[1], *Nabla_Psi_2c[0][1],building_precision, false, false, false);
    mrcpp::multiply(DOT_TOP_z, *Nabla_K_tree[2], *Nabla_Psi_2c[0][2],building_precision, false, false, false);
    mrcpp::multiply(DOT_BOTTOM_x, *Nabla_K_tree[0], *Nabla_Psi_2c[1][0],building_precision, false, false, false);
    mrcpp::multiply(DOT_BOTTOM_y, *Nabla_K_tree[1], *Nabla_Psi_2c[1][1],building_precision, false, false, false);
    mrcpp::multiply(DOT_BOTTOM_z, *Nabla_K_tree[2], *Nabla_Psi_2c[1][2],building_precision, false, false, false);

    // Now i sum the y and z components
    mrcpp::CompFunction<3> DOT_TOP_yz(MRA);
    mrcpp::CompFunction<3> DOT_BOTTOM_yz(MRA);
    mrcpp::CompFunction<3> DOT_TOP(MRA);
    mrcpp::CompFunction<3> DOT_BOTTOM(MRA);

    mrcpp::add(DOT_TOP_yz, 1.0, DOT_TOP_y, 1.0, DOT_TOP_z,building_precision, false);
    mrcpp::add(DOT_BOTTOM_yz, 1.0, DOT_BOTTOM_y, 1.0, DOT_BOTTOM_z,building_precision, false);

    // Finally i sum the x component to the yz component to the total SOC
    mrcpp::add(DOT_TOP, 1.0, DOT_TOP_x, 1.0, DOT_TOP_yz,building_precision, false);
    mrcpp::add(DOT_BOTTOM, 1.0, DOT_BOTTOM_x, 1.0, DOT_BOTTOM_yz,building_precision, false);

    // Multiply by K_inverse
    mrcpp::multiply(term_B_top, K_inverse, DOT_TOP,building_precision, false, false, false);
    mrcpp::multiply(term_B_bottom, K_inverse, DOT_BOTTOM,building_precision, false, false, false);

/* 
With this i was getting : term_B_top = 7.442213823661e-09	term_B_bottom = 0.000000000000e+00

    for (int i=0;i<3;i++){
        mrcpp::multiply(Nabla_K_Nabla_Psi_top[i], *Nabla_K_tree[i], *Nabla_Psi_2c[0][i],building_precision, false, false, false);   // <---- Here we dereference the pointer
    }
    for (int i=0;i<3;i++){
        mrcpp::multiply(Nabla_K_Nabla_Psi_bottom[i], *Nabla_K_tree[i], *Nabla_Psi_2c[1][i],building_precision, false, false, false);   // <---- Here we dereference the pointer
    }

    // Times K INVERSE
    CompFunction<3> Top_x(MRA);
    CompFunction<3> Top_y(MRA);
    CompFunction<3> Top_z(MRA);
    CompFunction<3> Bottom_x(MRA);
    CompFunction<3> Bottom_y(MRA);
    CompFunction<3> Bottom_z(MRA);

    CompFunction<3> Top_yz(MRA);
    CompFunction<3> Bottom_yz(MRA);

    mrcpp::multiply(Top_x, K_inverse, Nabla_K_Nabla_Psi_top[0],building_precision, false, false, false);
    mrcpp::multiply(Top_y, K_inverse, Nabla_K_Nabla_Psi_top[1],building_precision, false, false, false);
    mrcpp::multiply(Top_z, K_inverse, Nabla_K_Nabla_Psi_top[2],building_precision, false, false, false);

    mrcpp::multiply(Bottom_x, K_inverse, Nabla_K_Nabla_Psi_bottom[0],building_precision, false, false, false);
    mrcpp::multiply(Bottom_y, K_inverse, Nabla_K_Nabla_Psi_bottom[1],building_precision, false, false, false);
    mrcpp::multiply(Bottom_z, K_inverse, Nabla_K_Nabla_Psi_bottom[2],building_precision, false, false, false);

    // Now i sum the y and z components
    mrcpp::add(Top_yz, 1.0, Top_y, 1.0, Top_z,building_precision, false);
    mrcpp::add(Bottom_yz, 1.0, Bottom_y, 1.0, Bottom_z,building_precision, false);

    // Finally i sum the x component to the yz component to the total SOC
    mrcpp::add(term_B_top, 1.0, Top_x, 1.0, Top_yz,building_precision, false);
    mrcpp::add(term_B_bottom, 1.0, Bottom_x, 1.0, Bottom_yz,building_precision, false); */

}


void compute_term_C_Propagator(MultiResolutionAnalysis<3> &MRA, std::vector<std::vector<mrcpp::CompFunction<3> *>> &Nabla_Psi_2c , std::vector<mrcpp::CompFunction<3> *> &Nabla_K_tree, CompFunction<3> &K_inverse, CompFunction<3> &term_C_top,CompFunction<3> &term_C_bottom){
    // Where I'll put the rotor of K * (\Nabla \Psi)
    debug = true;

    // Now we compute the scalar product between the rotor and the sigma matrix vector
    CompFunction<3> SOC_Psi_t(MRA);
    CompFunction<3> SOC_Psi_b(MRA);
    
    std::vector<mrcpp::CompFunction<3>*> cross_top(3, new mrcpp::CompFunction<3>(MRA));
    std::vector<mrcpp::CompFunction<3>*> cross_bottom(3, new mrcpp::CompFunction<3>(MRA));
    
    std::vector<mrcpp::CompFunction<3> *> Nabla_Psi_top(3, new mrcpp::CompFunction<3>(MRA));
    std::vector<mrcpp::CompFunction<3> *> Nabla_Psi_bottom(3, new mrcpp::CompFunction<3>(MRA));

    for (int i=0;i<3;i++){
        Nabla_Psi_top[i] = Nabla_Psi_2c[0][i];
        Nabla_Psi_bottom[i] = Nabla_Psi_2c[1][i];
    }
 
    Cross_Product_Compfunct(cross_top, Nabla_K_tree, Nabla_Psi_top, MRA);
    Cross_Product_Compfunct(cross_bottom, Nabla_K_tree, Nabla_Psi_bottom, MRA);



    compute_sigma_cdot_spinor(MRA, cross_top, cross_bottom, SOC_Psi_t, SOC_Psi_b);
    //std::cout << " Correctly computed the sigma scalar product" << '\n';

    // Now the actual top and bottom components for C are the i * K_inverse * SOC_Psi

    // Times K INVERSE
    mrcpp::multiply(term_C_top, K_inverse, SOC_Psi_t,building_precision, false, false, false);
    mrcpp::multiply(term_C_bottom, K_inverse, SOC_Psi_b,building_precision, false, false, false);
    //std::cout << " Correctly computed the term C" << '\n';
    // Times i 
    std::cout << "NORM PRE RESCALE:" << '\n';
    std::cout << "term_C_top = " << term_C_top.getSquareNorm() << '\n';
    std::cout << "term_C_bottom = " << term_C_bottom.getSquareNorm() << '\n';
    term_C_top.rescale(ComplexDouble(0,1));
    term_C_bottom.rescale(ComplexDouble(0,1));
    //std::cout << " Correctly rescaled the term C" << '\n';
}


void compute_term_D_Propagator(MultiResolutionAnalysis<3> &MRA,std::vector<mrcpp::CompFunction<3>> &Psi_in ,CompFunction<3> &V , CompFunction<3> &K_inverse, CompFunction<3> &term_D_top,CompFunction<3> &term_D_bottom){
    // Term_D = (-2m/K) * V * Psi
    double factor = -2*m;
    CompFunction<3> V_Psi_t(MRA);
    CompFunction<3> V_Psi_b(MRA);

    mrcpp::multiply(building_precision, V_Psi_t, factor, V, Psi_in[0]);
    mrcpp::multiply(building_precision, V_Psi_b, factor, V, Psi_in[1]);
    
    std::cout << "NORM INPUT:" << '\n';
    std::cout << "V_Psi_top = " << V_Psi_t.getSquareNorm() << '\n';
    std::cout << "V_Psi_bottom = " << V_Psi_b.getSquareNorm() << '\n';
    
    mrcpp::multiply(building_precision, term_D_top ,   1, K_inverse,  V_Psi_t);
    mrcpp::multiply(building_precision, term_D_bottom ,    1, K_inverse, V_Psi_b);
    //mrcpp::multiply(building_precision, K_inverse, 1, term_D_top , V_Psi[0]);
    //mrcpp::multiply(building_precision, K_inverse, 1, term_D_bottom , V_Psi[1]);
    term_D_top.rescale(factor);
    term_D_bottom.rescale(factor);
}






void apply_Helmholtz_ZORA(MultiResolutionAnalysis<3> &MRA, std::vector<mrcpp::CompFunction<3>> &Psi_in, std::vector<std::vector<mrcpp::CompFunction<3> *>> &Nabla_Psi_2c, std::vector<mrcpp::CompFunction<3> *> &Nabla_K_tree,  CompFunction<3> &V, CompFunction<3> &K_tree, CompFunction<3> &K_inverse, double &E_n, std::vector<mrcpp::CompFunction<3>> &Psi_out){
    // Build the Helmholtz operator
    double mu = std::sqrt(-2*E_n);
    std::cout << "--------------------------------------------------" << '\n';
    std::cout << "Energy = " << E_n << '\n';
    std::cout << "Mu = " << mu << '\n';
    std::cout << "--------------------------------------------------" << '\n';

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
    
    //std::cout << "Now computing term A" << '\n';
    compute_term_A_Propagator(MRA, E_n, Psi_in, V, term_A_top, term_A_bottom);

    // Compute the term B
    //std::cout << "Now computing term B" << '\n';
    compute_term_B_Propagator(MRA, Psi_in, Nabla_Psi_2c, Nabla_K_tree, K_tree, K_inverse, term_B_top, term_B_bottom);
    
    // Compute the term C
    //std::cout << "Now computing term C" << '\n';
    compute_term_C_Propagator(MRA, Nabla_Psi_2c, Nabla_K_tree, K_inverse, term_C_top, term_C_bottom);

    // Compute the term D
    //std::cout << "Now computing term D" << '\n';
    compute_term_D_Propagator(MRA, Psi_in, V, K_inverse, term_D_top, term_D_bottom);

    // Now we sum all the terms, we have a total of 4 terms to sum. We'll take them by pairs and sum them
    CompFunction<3> first_pair_top(MRA);
    CompFunction<3> first_pair_bottom(MRA);
    CompFunction<3> second_pair_top(MRA);
    CompFunction<3> second_pair_bottom(MRA);

    // DEBUG
    std::cout << '\n' << "Norm of the components" << '\n';
    std::cout << "term_A_top = " << term_A_top.getSquareNorm() << '\t' << "term_A_bottom = " << term_A_bottom.getSquareNorm() << '\n';
    std::cout << "term_B_top = " << term_B_top.getSquareNorm() << '\t' << "term_B_bottom = " << term_B_bottom.getSquareNorm() << '\n';
    std::cout << "term_C_top = " << term_C_top.getSquareNorm() << '\t' << "term_C_bottom = " << term_C_bottom.getSquareNorm() << '\n';
    std::cout << "term_D_top = " << term_D_top.getSquareNorm() << '\t' << "term_D_bottom = " << term_D_bottom.getSquareNorm() << '\n';




    //std::cout << "Now computing the additions" << '\n';
    mrcpp::add(first_pair_top, 1.0, term_A_top, 1.0, term_B_top, building_precision, false);
    mrcpp::add(first_pair_bottom, 1.0, term_A_bottom, 1.0, term_B_bottom, building_precision, false);
    mrcpp::add(second_pair_top, 1.0, term_C_top, 1.0, term_D_top, building_precision, false);
    mrcpp::add(second_pair_bottom, 1.0, term_C_bottom, 1.0, term_D_bottom, building_precision, false);

    // SECOND DEBUG
    std::cout << '\n' << "Norm of the components after the addition" << '\n';
    std::cout << "first_pair_top = " << first_pair_top.getSquareNorm() << '\t' << "first_pair_bottom = " << first_pair_bottom.getSquareNorm() << '\n';
    std::cout << "second_pair_top = " << second_pair_top.getSquareNorm() << '\t' << "second_pair_bottom = " << second_pair_bottom.getSquareNorm() << '\n';


    // Now we sum the two pairs
    std::vector<mrcpp::CompFunction<3>> Psi_to_be_convoluted(2,MRA);

    //std::cout << "Now computing the final addition" << '\n';
    mrcpp::CompFunction<3> tmp1(MRA);
    mrcpp::CompFunction<3> tmp2(MRA);

//    mrcpp::add(Psi_to_be_convoluted[0], 1.0, first_pair_top, 1.0, second_pair_top, building_precision, false);
//    mrcpp::add(Psi_to_be_convoluted[1], 1.0, first_pair_bottom, 1.0, second_pair_bottom, building_precision, false);

    mrcpp::add(tmp1, 1.0, first_pair_top, 1.0, second_pair_top, building_precision, false);
    mrcpp::add(tmp2, 1.0, first_pair_bottom, 1.0, second_pair_bottom, building_precision, false);



    Psi_to_be_convoluted[0] = tmp1;
    Psi_to_be_convoluted[1] = tmp2;

    // LAST DEBUG
    std::cout << '\n' << "Norm of the components after the final addition" << '\n';
    std::cout << "Psi_to_be_convoluted[0] = " << Psi_to_be_convoluted[0].getSquareNorm() << '\n' << "Psi_to_be_convoluted[1] = " << Psi_to_be_convoluted[1].getSquareNorm() << '\n';
    //std::cout << "Real = " << Psi_to_be_convoluted[0].isreal() << '\t'  << "Complex = " << Psi_to_be_convoluted[0].iscomplex() << '\n';


    // Now we apply the Helmholtz operator
    std::cout << "Now applying the Helmholtz operator" << '\n';




    mrcpp::apply(building_precision, Psi_out[0], Helm, Psi_to_be_convoluted[0]);
    mrcpp::apply(building_precision, Psi_out[1], Helm, Psi_to_be_convoluted[1]);

}



void Renormalize_Spinor(std::vector<mrcpp::CompFunction<3>> &Psi_2c){
    double norm = std::sqrt(Psi_2c[0].getSquareNorm() + Psi_2c[1].getSquareNorm());
    Psi_2c[0].rescale(1.0/norm);
    Psi_2c[1].rescale(1.0/norm);

    norm = std::sqrt(Psi_2c[0].getSquareNorm() + Psi_2c[1].getSquareNorm());
    if (norm < 0.999 || norm > 1.001){
        std::cerr << "ERROR : Norm is not 1.0, it is " << norm << '\n';
    }
}




void Update_Nabla_Psi(MultiResolutionAnalysis<3> &MRA, std::vector<mrcpp::CompFunction<3>> &Psi_2c, std::vector<std::vector<mrcpp::CompFunction<3> *>> &Nabla_Psi_2c){
    // Operatr ABGV
    mrcpp::ABGVOperator<3> D(MRA, 0.0, 0.0);

    Nabla_Psi_2c[0] = mrcpp::gradient(D,Psi_2c[0]);
    Nabla_Psi_2c[1] = mrcpp::gradient(D,Psi_2c[1]);

}
