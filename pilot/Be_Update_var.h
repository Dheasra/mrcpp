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






void make_density_spinor(Spinor_CompFunction &Psi_2c, MultiResolutionAnalysis<3> &MRA, FunctionTree<3,double> &Rho){
    // This function will create the density for the spinor
    // The density is given by:
    // Rho = |Psi_2c[0]|^2 + |Psi_2c[1]|^2
    // We need to compute the norm of each component and then sum them up

    CompFunction<3> Rho_top(MRA);
    CompFunction<3> Rho_bottom(MRA);

    make_density_local(Rho_top, Psi_2c[0], MRA, building_precision);
    make_density_local(Rho_bottom, Psi_2c[1], MRA, building_precision);

    add(building_precision, Rho, 1.0, *(Rho_top.CompD[0]), 1.0, *(Rho_bottom.CompD[0]));
}

void make_density_spinor_separate_components(Spinor_CompFunction &Psi_2c, MultiResolutionAnalysis<3> &MRA, FunctionTree<3,double> &Rho_t, FunctionTree<3,double> &Rho_b){
    // This function will create the density for the spinor
    // The density is given by:
    // Rho = |Psi_2c[0]|^2 + |Psi_2c[1]|^2
    // We need to compute the norm of each component and then sum them up

    CompFunction<3> Rho_top(MRA);
    CompFunction<3> Rho_bottom(MRA);

    make_density_local(Rho_top, Psi_2c[0], MRA, building_precision);
    make_density_local(Rho_bottom, Psi_2c[1], MRA, building_precision);

    Rho_top.CompD[0]->deep_copy(&Rho_t);
    Rho_bottom.CompD[0]->deep_copy(&Rho_b);

}



void Update_Nabla_Psi_2c(Spinor_gradients_pointer_list &Nabla_Psi_list, Orbital_Pointer_List &Psi_list, MultiResolutionAnalysis<3> &MRA){
    // This function will update the Nabla_Psi_2c vector with the new Psi_2c vector
    // The Nabla_Psi_2c vector is a vector of vectors of CompFunction, where each vector corresponds to a component of the spinor
    // The first vector corresponds to the top component and the second vector corresponds to the bottom component

    //Define the derivsative operator
    ABGVOperator D(MRA, 1.0, 1.0);

    for (int i=0; i<4; i++){
        Nabla_Psi_list[i][0] = gradient(D, (*Psi_list[i])[0]);
        Nabla_Psi_list[i][1] = gradient(D, (*Psi_list[i])[1]);
        // Debug print
        //std::cout << "Psi[" << i << "] gradient TOP (x,y,z) " << Nabla_Psi_list[i][0][0]->getSquareNorm() << '\t' << Nabla_Psi_list[i][0][1]->getSquareNorm() << '\t' << Nabla_Psi_list[i][0][2]->getSquareNorm() << '\n';
        //std::cout << "Psi[" << i << "] gradient BOT (x,y,z) " << Nabla_Psi_list[i][1][0]->getSquareNorm() << '\t' << Nabla_Psi_list[i][1][1]->getSquareNorm() << '\t' << Nabla_Psi_list[i][1][2]->getSquareNorm() << '\n';
    }
    
}



void K_Psi(std::vector<mrcpp::CompFunction<3>> &K_ij_T, std::vector<mrcpp::CompFunction<3>> &K_ij_B, Orbital_Pointer_List &Psi_k, MultiResolutionAnalysis<3> &MRA, PoissonOperator &P){
    // THIS FUNCTION WILL EVALUATE ALL THE K |PSI> FOR THE 4 SPINORS

    // Variable declaration
   

    // I'll need these because make_density_local only works with CompFunction, making them some temporary variables
    std::vector<CompFunction<3>> Overlap_t_C(4);
    std::vector<CompFunction<3>> Overlap_b_C(4);

    std::vector<CompFunction<3>> Overlap_tot(4);
    std::vector<CompFunction<3>> Convoluted(4);

    std::vector<CompFunction<3>> K_i__j_TOP(4); // the 4 top components of K|Psi_i>
    std::vector<CompFunction<3>> K_i__j_BOT(4); // the 4 top components of K|Psi_i>


    // To sum them over
    std::vector<mrcpp::FunctionTree<3, double> *> K_i__j_TOP_RealFuncTree(4);
    std::vector<mrcpp::FunctionTree<3, double> *> K_i__j_BOT_RealFuncTree(4);

    // Same for the comlex function
    std::vector<mrcpp::FunctionTree<3, ComplexDouble> *> K_i__j_TOP_ComplexFuncTree(4);
    std::vector<mrcpp::FunctionTree<3, ComplexDouble> *> K_i__j_BOT_ComplexFuncTree(4);

    bool fetch_complex = false;
    
    
    for (int i=0; i<4; i++){
        for (int j=0; j<4; j++){

            // Take Psi^*_j Psi_i, for top and bottom components
            // OL_ij(r) = \Psi^*_j(r) * \Psi_i(r)
            multiply(Overlap_t_C[j], (*Psi_k[j])[0], (*Psi_k[i])[0], building_precision, false, false, true);
            multiply(Overlap_b_C[j], (*Psi_k[j])[1], (*Psi_k[i])[1], building_precision, false, false, true);
            std::cout << "---------------------------------------------------------------------------" << '\n';
            std::cout << "[i] = " << i << '\n';
            std::cout << "Overlap_t_C[" << j << "] = " << Overlap_t_C[j].getSquareNorm() << '\n';
            std::cout << "Overlap_b_C[" << j << "] = " << Overlap_b_C[j].getSquareNorm() << '\n';
            
            // Add them together 
            // OL_tot_ij(r) = [\PsiT^*_j(r) * \PsiT_i(r)] + [\PsiB^*_j(r) * \PsiB_i(r)]
            add(Overlap_tot[j], 1.0, Overlap_t_C[j], 1.0, Overlap_b_C[j], building_precision);
            std::cout << "Overlap_tot[" << j << "] = " << Overlap_tot[j].getSquareNorm() << '\n';
            
            // Convolute the density
            // Conv_ij(r) = \int OL_tot_ij(r') / \abs{r-r'} dr'
            apply(building_precision, Convoluted[j], P, Overlap_tot[j]);
            std::cout << "Convoluted[" << j << "] = " << Convoluted[j].getSquareNorm() << '\n';

            // Now i multiply each integral by the corresponding spinor
            multiply(K_i__j_TOP[j], Convoluted[j], (*Psi_k[j])[0], building_precision, false, false, false);
            multiply(K_i__j_BOT[j], Convoluted[j], (*Psi_k[j])[1], building_precision, false, false, false);

            std::cout << "K_["<<i << j << "] = " << K_i__j_TOP[j].getSquareNorm() << '\t' << K_i__j_BOT[j].getSquareNorm() << '\n';

            
            // Now i put the corresponding tree pointer in the vector to sum them up later
            if (K_i__j_TOP[j].iscomplex()){
                K_i__j_TOP_ComplexFuncTree[j] = K_i__j_TOP[j].CompC[0];
                fetch_complex = true;
            }
            else if (K_i__j_TOP[j].isreal()){
                K_i__j_TOP_RealFuncTree[j] = K_i__j_TOP[j].CompD[0];
                fetch_complex = false;
            }
            else{
                //std::cout << "K_i__j_TOP[j] is neither complex nor real" << '\n';
            }

            if (K_i__j_BOT[j].iscomplex()){
                K_i__j_BOT_ComplexFuncTree[j] = K_i__j_BOT[j].CompC[0];
                fetch_complex = true;
            }
            else if (K_i__j_TOP[j].isreal()){
                K_i__j_BOT_RealFuncTree[j] = K_i__j_BOT[j].CompD[0];
                fetch_complex = false;
            }
            else{
                //std::cout << "K_i__j_TOP[j] is neither complex nor real" << '\n';
            }
        }   

        //std::cout << '\n';
        //std::cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << '\n'; 
        //std::cout << "K_" << i << "0 top / bottom = " << K_i__j_TOP[0].getSquareNorm() << '\t' << K_i__j_BOT[0].getSquareNorm() << '\n';
        //std::cout << "K_" << i << "1 top / bottom = " << K_i__j_TOP[1].getSquareNorm() << '\t' << K_i__j_BOT[1].getSquareNorm() << '\n';
        //std::cout << "K_" << i << "2 top / bottom = " << K_i__j_TOP[2].getSquareNorm() << '\t' << K_i__j_BOT[2].getSquareNorm() << '\n';
        //std::cout << "K_" << i << "3 top / bottom = " << K_i__j_TOP[3].getSquareNorm() << '\t' << K_i__j_BOT[3].getSquareNorm() << '\n' << '\n';
        

        
        if (fetch_complex){
            // sum all the entries of the vector of FunctionTree getting the total density for top and bottom
            mrcpp::add(building_precision,  *K_ij_T[i].CompC[0], K_i__j_TOP_ComplexFuncTree);
            mrcpp::add(building_precision,  *K_ij_B[i].CompC[0], K_i__j_BOT_ComplexFuncTree);
        }
        else{
            // sum all the entries of the vector of FunctionTree getting the total density for top and bottom
            mrcpp::add(building_precision,  *K_ij_T[i].CompD[0], K_i__j_TOP_RealFuncTree);
            std::cout << "K_Psi_k[" << i << "][0] = " << K_ij_T[i].getSquareNorm() << '\n'; 
            mrcpp::add(building_precision,  *K_ij_B[i].CompD[0], K_i__j_BOT_RealFuncTree);

            
            std::cout << "K_Psi_k[" << i << "][1] = " << K_ij_B[i].getSquareNorm() << '\n' << '\n';
            std::cout << '\n';
            
            

        }
        
    }

    //Print the norms for K top and bottom
    std::cout << "K_ij[0] TOP AND BOT = " << K_ij_T[0].getSquareNorm() << '\t' << K_ij_B[0].getSquareNorm() << '\n';
    std::cout << "K_ij[1] TOP AND BOT = " << K_ij_T[1].getSquareNorm() << '\t' << K_ij_B[1].getSquareNorm() << '\n';
    std::cout << "K_ij[2] TOP AND BOT = " << K_ij_T[2].getSquareNorm() << '\t' << K_ij_B[2].getSquareNorm() << '\n';
    std::cout << "K_ij[3] TOP AND BOT = " << K_ij_T[3].getSquareNorm() << '\t' << K_ij_B[3].getSquareNorm() << '\n';




}










void J_Psi(std::vector<std::vector<mrcpp::CompFunction<3>*>> &J_Psi_k, Orbital_Pointer_List &Psi_k, CompFunction<3> Mean_field_potential, MultiResolutionAnalysis<3> &MRA, PoissonOperator &P){
    // THIS FUNCTION WILL EVALUATE ALL THE J |PSI> FOR THE 4 SPINORS
    std::cout << "J_Psi function called" << '\n';
    // Variable declaration
    std::vector<FunctionTree<3, double>*> Rho_t(4); // the 4 top components of the density
    std::vector<FunctionTree<3, double>*> Rho_b(4); // the 4 bottom components of the density

    FunctionTree<3, double> Rho_top_TOT(MRA);       // This will hold the sum of the top components
    FunctionTree<3, double> Rho_bottom_TOT(MRA);    // This will hold the sum of the bottom components
    FunctionTree<3, double> Rho_TOTAL(MRA);         // Will hold the sum of the previous two

    // I'll need these because make_density_local only works with CompFunction, making them some temporary variables
    std::vector<CompFunction<3>> Rho_t_C(4);
    std::vector<CompFunction<3>> Rho_b_C(4);
    
    for (int k=0; k<4; k++){
        // Store in Rho_t_C and Rho_b_C the density of the top and bottom components
        make_density_local(Rho_t_C[k], (*Psi_k[k])[0], MRA, building_precision);
        make_density_local(Rho_b_C[k], (*Psi_k[k])[1], MRA, building_precision);

        // Transfer all the data to the Rho_t and Rho_b vectors of FunctionTree, so to be able to sum them up
        Rho_t[k] = Rho_t_C[k].CompD[0];
        Rho_b[k] = Rho_b_C[k].CompD[0];
                        
        //std::cout << "Rho_t_C " << k << " = " << Rho_t_C[k].CompD[0]->getSquareNorm() << '\n';
        //std::cout << "Rho_b_C " << k << " = " << Rho_b_C[k].CompD[0]->getSquareNorm() << '\n';    
    }

    // sum all the entries of the vector of FunctionTree getting the total density for top and bottom
    mrcpp::add(building_precision, Rho_top_TOT, Rho_t);
    mrcpp::add(building_precision, Rho_bottom_TOT, Rho_b);
    // And sum the two components to get the total density
    mrcpp::add(building_precision, Rho_TOTAL, 1.0, Rho_top_TOT, 1.0, Rho_bottom_TOT);
    std::cout << "Rho_top_TOT = " << Rho_top_TOT.getSquareNorm() << '\n';
    
    // Convolute the density   
    mrcpp::apply(building_precision, *Mean_field_potential.CompD[0], P, Rho_TOTAL);
    std::cout << "Mean_field_potential = " << Mean_field_potential.getSquareNorm() << '\n';
    

    // Set every elemeent as the mean field potential times the respective spinors 
    for (int k=0; k<4; k++){
        // J_psi_k = V_mean_field * Psi_k 
        mrcpp::multiply(*(J_Psi_k[k][0]), Mean_field_potential, (*Psi_k[k])[0], building_precision);
        mrcpp::multiply(*(J_Psi_k[k][1]), Mean_field_potential, (*Psi_k[k])[1], building_precision);
        
    }
  
    
}










void Update_SCF_Variables(MultiResolutionAnalysis<3> &MRA, ABGVOperator<3> &D, std::vector<mrcpp::CompFunction<3>> &Psi_2c, CompFunction<3> &V_electron_nucleus ,std::vector<std::vector<mrcpp::CompFunction<3> *>> &Nabla_Psi_2c,CompFunction<3> &K_tree, std::vector<mrcpp::CompFunction<3> *> &Nabla_K_tree, CompFunction<3> &K_inverse, CompFunction<3> &V, CompFunction<3> &J, mrcpp::PoissonOperator &P){
    // Updating the Nabla_Psi_2c
    //std::cout << "--------------UPDATING GRADIENT--------------" << '\n';
    Nabla_Psi_2c[0] = mrcpp::gradient(D,Psi_2c[0]);
    Nabla_Psi_2c[1] = mrcpp::gradient(D,Psi_2c[1]);
    // Print  values
    if (debug){
    std::cout << "Nabla_Psi_2c[0] = " << Nabla_Psi_2c[0][0]->getSquareNorm() << "   " <<  Nabla_Psi_2c[0][1]->getSquareNorm() << "   " <<  Nabla_Psi_2c[0][2]->getSquareNorm() << '\n';
    std::cout << "Nabla_Psi_2c[1] = " << Nabla_Psi_2c[1][0]->getSquareNorm() << "   " <<  Nabla_Psi_2c[1][1]->getSquareNorm() << "   " <<  Nabla_Psi_2c[1][2]->getSquareNorm() << '\n';
    std::cout << '\n';
    
    // Same format, print the pointers
    std::cout << "Nabla_Psi_2c[0][0] = " << Nabla_Psi_2c[0][0] << "   " << Nabla_Psi_2c[0][1] << "   " << Nabla_Psi_2c[0][2] << '\n';
    std::cout << "Nabla_Psi_2c[1][0] = " << Nabla_Psi_2c[1][0] << "   " << Nabla_Psi_2c[1][1] << "   " << Nabla_Psi_2c[1][2] << '\n';
    }
    //std::cout << "------------------- DONE -------------------" << '\n' << '\n';
    
    // Create the Rho_t and Rho_b, then Rho = Rho_t + Rho_b
    //std::cout << "--------------UPDATING RHO------------------" << '\n';
    CompFunction<3> Rho_t(MRA);
    CompFunction<3> Rho_b(MRA);
    CompFunction<3> Rho(MRA);


    //std::cout << "Computing density for TOP component..." << '\n';
    make_density_local(Rho_t, Psi_2c[0] , MRA, building_precision);
    //std::cout << "Norm of Rho_top = " << Rho_t.getSquareNorm() << '\n';

  
    //std::cout << "Computing density for BOTTOM component..." << '\n';
    make_density_local(Rho_b, Psi_2c[1] , MRA, building_precision);
    //std::cout << "Norm of Rho_bottom = " << Rho_b.getSquareNorm() << '\n' << '\n';
  

    mrcpp::add(Rho, 1.0, Rho_t, 1.0, Rho_b, building_precision, false);
    
    if ( num_cycle !=0  && (Rho_b.getSquareNorm() == 0.0 || Rho_t.getSquareNorm() == 0.0)) {
        std::cerr << "Error: Rho_t or Rho_b has zero norm!" << '\n';
        exit(1);
        return;
    }

    //std::cout << "---------------- DONE -------------------" << '\n' << '\n';

    // Calculate the J operator
    //std::cout << "--------------UPDATING J-----------------" << '\n';
    //std::cout << "Rho status" << *(Rho.CompD[0]) << '\n';
    mrcpp::apply(building_precision, J, P, Rho);
    //std::cout << "Applied the Poisson operator to the density" << '\n';
    double two_pi = (2.0 * M_PI);
    J.rescale(ComplexDouble(2.0 ,0.0));  
    //std::cout << "---------------- DONE -------------------" << '\n' << '\n';


    // Compute the full potential
    //std::cout << "--------------UPDATING V------------------" << '\n';
    mrcpp::add(V, 1.0, V_electron_nucleus,0.5, J, building_precision, false); // V = V_electron_nucleus + 0.5 J = V_ce(r)+J(r)
    //std::cout << "V norm = " << V.getSquareNorm() << '\n';
    //std::cout << "---------------- DONE -------------------" << '\n' << '\n';





    // Compute the inverse of K_inverse (that is just the Kappa operator):
    //std::cout << "--------------UPDATING K_tree------------------" << '\n';
    double constant = 2.0 * m * c * c;
    double prefactor = 1/constant;


    // Define and project a the function whcih is one everywhere
    std::function<double(const Coord<3> &x)> one = [] (const mrcpp::Coord<3> &r) -> double {
        return 1.0;
    };
    mrcpp::CompFunction<3> one_tree(MRA);
    mrcpp::project(one_tree, one, building_precision);

    // Function to calculate Kappa
    std::function<double(const Coord<3> &x)> K_function_compute = [&V_electron_nucleus, &J, prefactor, &V] (const mrcpp::Coord<3> &r) -> double {
        double VAL = 1.0 - (prefactor * V.CompD[0]->evalf(r));
        return (1.0 / VAL) -1.0;   
    };
    mrcpp::project(K_tree, K_function_compute, building_precision);
    //std::cout << "K_tree norm NO MAP= " << K_tree.getSquareNorm() << '\n';
 
    K_tree.add(ComplexDouble(1.0,0.0), one_tree); // K_tree = K_tree + 1

    
    if (debug){
    std::cout << "K_tree norm = " << K_tree.getSquareNorm() << '\n';
    std::cout << "K_tree is real = " << K_tree.isreal() << '\n';
    std::cout << "K_tree is complex = " << K_tree.iscomplex() << '\n';
    }
    //std::cout << "---------------- DONE -------------------" << '\n' << '\n';





    // Compute the new K tree
    //std::cout << "-------------UPDATING K_inv---------------" << '\n';


    // K_inverse
    std::function<double(const Coord<3> &x)> K_inverse_function = [&K_tree, &V, prefactor, &V_electron_nucleus, &J] (const mrcpp::Coord<3> &r) -> double {
        return - (prefactor * (V.CompD[0]->evalf(r)));
        //return 1.0 / K_tree.CompD[0]->evalf(r);   
    };

    mrcpp::project(K_inverse, K_inverse_function, building_precision);
    //add: K_inverse + 1
    K_inverse.add(ComplexDouble(1.0,0.0), one_tree);
    

    
    if (debug){
    std::cout << "K_inverse norm = " << K_inverse.getSquareNorm() << '\n';
    std::cout << "K_inverse is real = " << K_inverse.isreal() << '\n';
    std::cout << "K_inverse is complex = " << K_inverse.iscomplex() << '\n';
    std::cout << "K_inv Function Tree:" << '\n';
    std::cout << *(K_inverse.CompD[0]) << '\n';
    }
    



    //std::cout << "---------------- DONE -------------------" << '\n' << '\n';



    



    // Compute the Nabla_K_tree
    //std::cout << "--------------UPDATING Nabla_K_tree------------------" << '\n';
    Nabla_K_tree = mrcpp::gradient(D,K_tree);

    if (debug){
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
    }
    //std::cout << "-------------- DONE UPDATING SCF VARIABLES -------------------" << '\n' << '\n';

  
}






ComplexDouble compute_Term1_T_ZORA_MEL(MultiResolutionAnalysis<3> &MRA, std::vector<std::vector<mrcpp::CompFunction<3>*>> &Nabla_Psi_2c, mrcpp::CompFunction<3> &K_tree,std::vector<mrcpp::CompFunction<3>> &Psi_2c_Bra){
    // Nabla(\Psi K) = \Nabla \Psi  K + \Psi  \Nabla K
    CompFunction<3> Psi_t_K(MRA);
    CompFunction<3> Psi_b_K(MRA);

    std::vector<mrcpp::CompFunction<3> *> Nabla_Psi_t = Nabla_Psi_2c[0];
    std::vector<mrcpp::CompFunction<3> *> Nabla_Psi_b = Nabla_Psi_2c[1];


    // Compute Psi_top * K
    mrcpp::multiply<3>(Psi_t_K, Psi_2c_Bra[0], K_tree,building_precision);
    // Compute Psi_bottom * K
    mrcpp::multiply<3>(Psi_b_K, Psi_2c_Bra[1], K_tree,building_precision);


    // Operatr ABGV
    mrcpp::ABGVOperator<3> D(MRA, 0.0, 0.0);

    // Gradient of | Psi_top * K >
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


ComplexDouble compute_Term2_T_ZORA_MEL(MultiResolutionAnalysis<3> &MRA, std::vector<std::vector<mrcpp::CompFunction<3>*>> &Nabla_Psi_2c, mrcpp::CompFunction<3> &K_tree, std::vector<mrcpp::CompFunction<3> *> &Nabla_K_tree,  std::vector<mrcpp::CompFunction<3>> &Psi_bra){
    /// NEW_APPROACH: I first do all the \nabla(K) * \nabla(\Psi) and then multiply by the Psi from the left as a bra
    std::vector<mrcpp::CompFunction<3>> Nabla_K_Nabla_Psi_top(3, MRA);
    std::vector<mrcpp::CompFunction<3>> Nabla_K_Nabla_Psi_bottom(3, MRA);

    // \Naabla(K) * |\Nabla(\Psi)>
    /*
    for (int i=0;i<3;i++){
        mrcpp::multiply(Nabla_K_Nabla_Psi_top[i], *Nabla_K_tree[i], *Nabla_Psi_2c[0][i],building_precision, false, false, false);   // <---- Here we dereference the pointer
    }
    for (int i=0;i<3;i++){
        mrcpp::multiply(Nabla_K_Nabla_Psi_bottom[i], *Nabla_K_tree[i], *Nabla_Psi_2c[1][i],building_precision, false, false, false);   // <---- Here we dereference the pointer
    }

    */

    CompFunction<3> Nabla_K_Nabla_Psi_XTOP(MRA);
    CompFunction<3> Nabla_K_Nabla_Psi_YTOP(MRA);
    CompFunction<3> Nabla_K_Nabla_Psi_ZTOP(MRA);
    CompFunction<3> Nabla_K_Nabla_Psi_XBOTTOM(MRA);
    CompFunction<3> Nabla_K_Nabla_Psi_YBOTTOM(MRA);
    CompFunction<3> Nabla_K_Nabla_Psi_ZBOTTOM(MRA);

    // Let's do explicitely what the loops do:
    mrcpp::multiply(Nabla_K_Nabla_Psi_XTOP, *Nabla_K_tree[0], *Nabla_Psi_2c[0][0],building_precision, false, false, false);   
    mrcpp::multiply(Nabla_K_Nabla_Psi_YTOP, *Nabla_K_tree[1], *Nabla_Psi_2c[0][1],building_precision, false, false, false);
    mrcpp::multiply(Nabla_K_Nabla_Psi_ZTOP, *Nabla_K_tree[2], *Nabla_Psi_2c[0][2],building_precision, false, false, false);
    mrcpp::multiply(Nabla_K_Nabla_Psi_XBOTTOM, *Nabla_K_tree[0], *Nabla_Psi_2c[1][0],building_precision, false, false, false);
    mrcpp::multiply(Nabla_K_Nabla_Psi_YBOTTOM, *Nabla_K_tree[1], *Nabla_Psi_2c[1][1],building_precision, false, false, false);
    mrcpp::multiply(Nabla_K_Nabla_Psi_ZBOTTOM, *Nabla_K_tree[2], *Nabla_Psi_2c[1][2],building_precision, false, false, false);


    //ComplexDouble Top_contribution  = dot(Psi_2c[0],Nabla_K_Nabla_Psi_top[0]) + dot(Psi_2c[0],Nabla_K_Nabla_Psi_top[1]) + dot(Psi_2c[0],Nabla_K_Nabla_Psi_top[2]);
    //ComplexDouble Bottom_contribution  = dot(Psi_2c[1],Nabla_K_Nabla_Psi_bottom[0]) + dot(Psi_2c[1],Nabla_K_Nabla_Psi_bottom[1]) + dot(Psi_2c[1],Nabla_K_Nabla_Psi_bottom[2]);

    ComplexDouble Top_contribution  = dot(Psi_bra[0],Nabla_K_Nabla_Psi_XTOP) + dot(Psi_bra[0],Nabla_K_Nabla_Psi_YTOP) + dot(Psi_bra[0],Nabla_K_Nabla_Psi_ZTOP);
    ComplexDouble Bottom_contribution  = dot(Psi_bra[1],Nabla_K_Nabla_Psi_XBOTTOM) + dot(Psi_bra[1],Nabla_K_Nabla_Psi_YBOTTOM) + dot(Psi_bra[1],Nabla_K_Nabla_Psi_ZBOTTOM);


    return Top_contribution + Bottom_contribution;
}



ComplexDouble compute_Term3_T_ZORA_MEL(MultiResolutionAnalysis<3> &MRA, std::vector<std::vector<mrcpp::CompFunction<3>*>> &Nabla_Psi_2c, mrcpp::CompFunction<3> &K_tree, std::vector<mrcpp::CompFunction<3> *> &Nabla_K_tree,  std::vector<mrcpp::CompFunction<3>> &Psi_2c_bra){
    // Compute the product K * (\Nabla \Psi) for top and bottom components
    std::vector<mrcpp::CompFunction<3>*> K_Nabla_Psi_top(3, new mrcpp::CompFunction<3>(MRA));
    std::vector<mrcpp::CompFunction<3>*> K_Nabla_Psi_bottom(3, new mrcpp::CompFunction<3>(MRA));

    std::vector<mrcpp::CompFunction<3> *> Nabla_Psi_top = Nabla_Psi_2c[0];
    std::vector<mrcpp::CompFunction<3> *> Nabla_Psi_bottom = Nabla_Psi_2c[1];

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
    
    ComplexDouble Top_contribution = dot(Psi_2c_bra[0], SOC_Psi_t);
    ComplexDouble Bottom_contribution = dot(Psi_2c_bra[1], SOC_Psi_b);

    return ComplexDouble(0,1)*(Top_contribution + Bottom_contribution);
}


ComplexDouble F_matrix_element(int row, int col,MultiResolutionAnalysis<3> &MRA, Spinor_gradients_pointer_list &Nabla_Psi_List, mrcpp::CompFunction<3> &Kappa_tree, std::vector<mrcpp::CompFunction<3> *> &Nabla_Kappa_tree,  Orbital_Pointer_List Psi_List, CompFunction<3> &V,  std::vector<mrcpp::CompFunction<3>> J_List_T, std::vector<mrcpp::CompFunction<3>> J_List_B, std::vector<mrcpp::CompFunction<3>> K_List_T, std::vector<mrcpp::CompFunction<3>> K_List_B, bool Energy_Orbital){
    
    if (verbose){std::cout << '\t' << "Computing Kinetic Term 1..." << '\n';}
    ComplexDouble Term1 = compute_Term1_T_ZORA_MEL(MRA, Nabla_Psi_List[row], Kappa_tree,  (*Psi_List[col]));
    if (verbose){std::cout << '\t' << "Computing Kinetic Term 2..." << '\n';}
    ComplexDouble Term2 = compute_Term2_T_ZORA_MEL(MRA, Nabla_Psi_List[row], Kappa_tree, Nabla_Kappa_tree, (*Psi_List[col]));
    if (verbose){std::cout << '\t' << "Computing Kinetic Term 3..." << '\n';}
    ComplexDouble Term3 = compute_Term3_T_ZORA_MEL(MRA, Nabla_Psi_List[row], Kappa_tree, Nabla_Kappa_tree, (*Psi_List[col]));
    if (verbose){std::cout << '\t' << "Kinetic Terms computed!" << '\n';}
    // Now we do <Psi | V | Psi>
    mrcpp::CompFunction<3> psi_top__V(MRA);
    mrcpp::CompFunction<3> psi_bottom__V(MRA);

    // V*|Psi>
    //std::cout << "V is real = " << V.isreal() << '\n';
    //std::cout << "V is complex = " << V.iscomplex() << '\n';
    //std::cout << "Psi_List[row][0] is real = " << (*Psi_List[row])[0].isreal() << '\n';
    //std::cout << "Psi_List[row][0] is complex = " << (*Psi_List[row])[0].iscomplex() << '\n';
    mrcpp::multiply(psi_top__V, V, (*Psi_List[row])[0],building_precision, false, false, false);
    mrcpp::multiply(psi_bottom__V, V, (*Psi_List[row])[1],building_precision, false, false, false);

    if (verbose){std::cout << '\t' << "Computing Potential Term..." << '\n' <<'\n';}

    ComplexDouble monoel_pot_energy = dot((*Psi_List[col])[0],psi_top__V) + dot((*Psi_List[col])[1],psi_bottom__V);

    // BIELECTRONIC PART
    ComplexDouble Coulomb_repulsion = dot((*Psi_List[col])[0],J_List_T[row]) + dot((*Psi_List[col])[1],J_List_B[row]);
    ComplexDouble Exchange_attraction = dot((*Psi_List[col])[0],K_List_T[row]) + dot((*Psi_List[col])[1],K_List_B[row]);
    ComplexDouble biel_energy = Coulomb_repulsion - Exchange_attraction;


    double kin_prefactor = -1.0/(2.0*m);
    ComplexDouble kinetic_energy = kin_prefactor * (Term1 + Term2 + Term3) ;

    double real_Z = double(Z);
    ComplexDouble Fock_ij =  kinetic_energy + monoel_pot_energy + biel_energy;
    
    if (Energy_Orbital){
        std::cout << "Element [" << row << "][" << col << "] of the Fock matrix:" << '\n'; 
        std::cout << "Kinetic energy = " << kinetic_energy.real() << '\n';
        std::cout << '\t' << "term 1 = " << kin_prefactor * Term1.real() << '\n';
        std::cout << '\t' << "term 2 = " << kin_prefactor * Term2.real() << '\n';
        std::cout << '\t' << "term 3 = " << kin_prefactor * Term3.real() << '\n';
        std::cout << "Mono-electronic potential energy = " << monoel_pot_energy.real() << '\n';
        std::cout << "Coulomb repulsion = " << 0.5 * Coulomb_repulsion.real() << '\n';
        std::cout << "Exchange attraction = " << 0.5 * Exchange_attraction.real() << '\n';
        std::cout << "Bielectronic energy = " << 0.5 * biel_energy.real() << '\n';
        std::cout << "tot = " << kinetic_energy+monoel_pot_energy + 0.5*biel_energy << '\n' << '\n';

        
        return Fock_ij - 0.5*biel_energy;
        
    }else{
        return Fock_ij;
    }

}


void Update_Fock(std::vector<double> &En, std::vector<std::vector<std::complex<double>>> &F_matrix, MultiResolutionAnalysis<3> &MRA, Spinor_gradients_pointer_list &Nabla_Psi_List, mrcpp::CompFunction<3> &Kappa_tree, std::vector<mrcpp::CompFunction<3> *> &Nabla_Kappa_tree, Orbital_Pointer_List &Psi_List, CompFunction<3> &V, std::vector<mrcpp::CompFunction<3>> J_List_T, std::vector<mrcpp::CompFunction<3>> J_List_B, std::vector<mrcpp::CompFunction<3>> K_List_T, std::vector<mrcpp::CompFunction<3>> K_List_B){
    // This function will update the Fock matrix with the new values of the orbitals
    if (verbose){std::cout << "Updating Fock matrix..." << '\n';}
    ComplexDouble tmp_ComplexEnergy;
    for (int i=0; i<4; i++){
        for (int j=i; j<4; j++){
            F_matrix[i][j] = F_matrix_element(i,j,MRA,Nabla_Psi_List,Kappa_tree,Nabla_Kappa_tree,Psi_List,V, J_List_T, J_List_B, K_List_T, K_List_B, false); // Fock matrix element
            F_matrix[j][i] = std::conj(F_matrix[i][j]); // Fock matrix is hermitian   
            if (i == j){
                std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++" << '\n';
                std::cout << "Computing energy for orbital " << i << '\n';
                std::cout << "F_matrix[" << i << "][" << j << "] = " << '\n';
                tmp_ComplexEnergy = F_matrix_element(i,j,MRA,Nabla_Psi_List,Kappa_tree,Nabla_Kappa_tree,Psi_List,V,J_List_T, J_List_B, K_List_T, K_List_B, true); // Energy of the orbital
                En[i] = tmp_ComplexEnergy.real();
            }
        }
    }
    
    if (verbose){std::cout << "Fock matrix updated!" << '\n';}

}



void Update_Kappa(CompFunction<3> &Kappa_tree, CompFunction<3> &Kappa_inverted_tree, std::vector<mrcpp::CompFunction<3> *> &Nabla_Kappa_tree, MultiResolutionAnalysis<3> &MRA, CompFunction<3> &V, CompFunction<3> &Mean_field_potential){
    // This function will update the Kappa operator with the new values of the orbitals
    if (verbose){std::cout << "Updating Kappa operator..." << '\n';}
    CompFunction<3> Tot_potential(MRA);
    mrcpp::add(Tot_potential, 1.0, V, 1.0, Mean_field_potential, building_precision, false); // Kappa = V + 0.5 * Mean_field_potential
    

    //std::cout << "--------------UPDATING K_tree------------------" << '\n';
    double constant = 2.0 * m * c * c;
    double prefactor = 1/constant;


    // Define and project a the function whcih is one everywhere
    std::function<double(const Coord<3> &x)> one = [] (const mrcpp::Coord<3> &r) -> double {
        return 1.0;
    };
    mrcpp::CompFunction<3> one_tree(MRA);
    mrcpp::project(one_tree, one, building_precision);

    // Function to calculate Kappa
    std::function<double(const Coord<3> &x)> K_function_compute = [&Tot_potential, prefactor] (const mrcpp::Coord<3> &r) -> double {
        double VAL = 1.0 - (prefactor * Tot_potential.CompD[0]->evalf(r));
        return (1.0 / VAL) -1.0;   
    };
    mrcpp::project(Kappa_tree, K_function_compute, building_precision);
    //std::cout << "K_tree norm NO MAP= " << K_tree.getSquareNorm() << '\n';
 
    Kappa_tree.add(ComplexDouble(1.0,0.0), one_tree); // K_tree = K_tree + 1

    
    if (debug){
    std::cout << "K_tree norm = " << Kappa_tree.getSquareNorm() << '\n';
    std::cout << "K_tree is real = " << Kappa_tree.isreal() << '\n';
    std::cout << "K_tree is complex = " << Kappa_tree.iscomplex() << '\n';
    }
    //std::cout << "---------------- DONE -------------------" << '\n' << '\n';

    if (verbose){std::cout << "Kappa operator updated!" << '\n';}



    //std::cout << "-------------UPDATING K_inv---------------" << '\n';
    // K_inverse
    std::function<double(const Coord<3> &x)> K_inverse_function = [&Tot_potential, prefactor] (const mrcpp::Coord<3> &r) -> double {
        return - (prefactor * (Tot_potential.CompD[0]->evalf(r)));
        //return 1.0 / K_tree.CompD[0]->evalf(r);   
    };

    mrcpp::project(Kappa_inverted_tree, K_inverse_function, building_precision);
    //add: K_inverse + 1
    Kappa_inverted_tree.add(ComplexDouble(1.0,0.0), one_tree);
    

    
    if (debug){
    std::cout << "K_inverse norm = " << Kappa_inverted_tree.getSquareNorm() << '\n';
    std::cout << "K_inverse is real = " << Kappa_inverted_tree.isreal() << '\n';
    std::cout << "K_inverse is complex = " << Kappa_inverted_tree.iscomplex() << '\n';
    std::cout << "K_inv Function Tree:" << '\n';
    std::cout << *(Kappa_inverted_tree.CompD[0]) << '\n';
    }

    //std::cout << "---------------- DONE -------------------" << '\n' << '\n';



    //std::cout << "--------------UPDATING Nabla_K_tree------------------" << '\n';
    ABGVOperator<3> D(MRA, 0.0, 0.0);
    Nabla_Kappa_tree = mrcpp::gradient(D,Kappa_tree);

    if (debug){
    // Print the pointers of the Nabla_K_tree
    std::cout << "Nabla_K_tree[0] = " << Nabla_Kappa_tree[0] << '\n';
    std::cout << "Nabla_K_tree[1] = " << Nabla_Kappa_tree[1] << '\n';
    std::cout << "Nabla_K_tree[2] = " << Nabla_Kappa_tree[2] << '\n' <<'\n';
    // Now we print the norms
    std::cout << "Nabla_K_tree[0] norm = " << Nabla_Kappa_tree[0]->getSquareNorm() << '\n';
    std::cout << "Nabla_K_tree[1] norm = " << Nabla_Kappa_tree[1]->getSquareNorm() << '\n';
    std::cout << "Nabla_K_tree[2] norm = " << Nabla_Kappa_tree[2]->getSquareNorm() << '\n' <<'\n';

    std::cout << "Nabla_K_tree[0] is real = " << Nabla_Kappa_tree[0]->isreal() << '\n';
    std::cout << "Nabla_K_tree[0] is complex = " << Nabla_Kappa_tree[0]->iscomplex() << '\n';
    }
    //std::cout << "-------------- DONE UPDATING SCF VARIABLES -------------------" << '\n' << '\n';


}



