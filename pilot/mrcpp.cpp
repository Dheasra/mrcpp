#include <cstdlib>
#include <numeric>
#include <iostream>
#include <vector>


#include "../api/Printer"
#include "../api/Timer"
#include "../api/Gaussians"
#include "../api/Plotter"
#include "../api/MWFunctions"
#include "../api/MWOperators"



using namespace mrcpp;

// DEBUG
bool debug = false;
bool verbose = false;

// GLOBAL VARIABLES
int Z =2;
//double c = 137.035999146; // --> This is the speed of light in atomic units
double c = 137.035999177; 
double m = 1.0;
int n_electrons;


int order;
int MaxLevel;
double building_precision;
double epsilon;
int Relativity;
int num_cycle =  0;     //    -> SCF cycle counter 



#include "my_utilities.h"
//#include "ZORA_utilities.h"

#include "He_atom_ZORA_utils.h"
#include "He_atom_utils.h"
#include "Updating_var_ZORA.h"
#include "Smeared_potential.h"
//#include "ZORA_utilities.h"

/*
 * ==================================================================================================================================
 * This program will introduce relativistic effects in the H atom, in detail, with the ZORA approximation.
 * 
 * In this framework we can say that the Hamiltonian is given by:
 *      H_{ZORA} = T_{ZORA} + V
 * V Is the same as the H atom in general. While T_{ZORA} is different. From classical to relativistic:
 *      T = - p^2/2m -> T_{ZORA} = p \cdot K(r) p
 * Where this K term is given by:
 *     K(r) = [1-V/(2mc^2)]^{-1}
 * Remember that each p is an operator, so the product is not trivial.
 * ==================================================================================================================================
*/


// ==================================================================================================================================
// ================================================ START OF THE MAIN FUNCTION ======================================================
// ==================================================================================================================================

int main(int argc, char **argv) {
    Timer timer;

    // Initialize printing
    int printlevel = 3;
    Printer::init(printlevel);
    print::environment(0);
    print::header(0, "MRCPP pilot code");
    
// ==================================================================================================================================
// ================= Start of the parameters reading ================================================================================
// ==================================================================================================================================
    auto parameters = readParameters("Input_Parameters.txt");

    // Assegna i valori ai parametri corrispondenti
    order = static_cast<int>(parameters["order"]);
    MaxLevel = static_cast<int>(parameters["MaxLevel"]);
    building_precision = parameters["building_precision"];
    epsilon = parameters["epsilon"];
    Z = static_cast<int>(parameters["Z"]);
    // Set the charge as neutra by:
    Z = 2;
    n_electrons = 2; // --> As we are dealing with the H type atom
    

    std::cout << "Parameters read from file: " << '\n' << '\n';
    std::cout << " Order =" << '\t'<< '\t' << order << '\n';
    std::cout << " MaxLevel =" << '\t'<< '\t' << MaxLevel << '\n';
    std::cout << " Building precision =" << '\t' << building_precision << '\n';
    std::cout << " Epsilon =" << '\t'<< '\t' << epsilon << '\n';
    std::cout << " Z =" << '\t'<< '\t' << '\t' << Z << '\n';
    std::cout << " Number of electrons =" << '\t' << n_electrons << '\n';
    
    
    
    
    
    
    
    
    
    
    std::cout << '\n' << '\n' << "************************************************************" << '\n';
    




    // ==================================================================================================================================
    // ================= End of the parameters reading ==================================================================================
    // ==================================================================================================================================

    // 1) Build the Box and MRA
    std::array<int,2> box = {-20,20};
    mrcpp::BoundingBox<3> bb(box);
    mrcpp::MultiResolutionAnalysis mra(bb, order, MaxLevel);
    // Also, we create the trees for each function we'll be dealing with
    mrcpp::CompFunction<3> VPhi_tree_top(mra);  // -> This is the tree that will hold the function
    mrcpp::CompFunction<3> VPhi_tree_bottom(mra);  // -> This is the tree that will hold the function
    mrcpp::CompFunction<3> Trial_function_tree_top(mra); // -> This is the tree that will hold the output of the convolution
    mrcpp::CompFunction<3> Trial_function_tree_bottom(mra); // -> This is the tree that will hold the output of the convolution
    mrcpp::CompFunction<3> core_el_tree(mra); // -> This is the tree that will hold the [nucleus-electron] potential
    mrcpp::CompFunction<3> J_function_tree(mra);  // -> This is the tree that will hold <b a | b a> = <\psi | J - K | \psi>



    mrcpp::CompFunction<3> Potential_tree(mra); // -> This is the tree that will hold the total potential
    // Then the 2 auxiliary trees
    mrcpp::CompFunction<3> Normalized_difference_tree(mra);
    mrcpp::CompFunction<3> Phi_n_copy_tree(mra);




    
    // ==================================================================================================================================
    // ================= My function will be now a 1s H atom orbital in 3D, given by the following lambda function ======================
    // ==================================================================================================================================

    // 2) We define declare operator G^\mu as the Helmltz Convolutional Operator:      [G^\mu]
    double mu = 1.0;



    // 3) Define a trial function as a Gaussian 
    mrcpp::Coord<3> pos = {0,0,0};
    std::array<int,3> pow = {0,0,0};
    
    // Define the trial function as a Slater-type orbital
    std::function<double(const Coord<3> &x)> slater = [] (const mrcpp::Coord<3> &r) -> double {
        auto R = std::sqrt(r[0]*r[0] + r[1]*r[1] + r[2]*r[2]);
        double lambda = std::sqrt(1. - (1/(c*c)) );
        double alpha = 1.69;
        return exp(-alpha*R);
    };
    mrcpp::CompFunction He_NR_Solution(mra);
    mrcpp::project(He_NR_Solution, slater, building_precision);

    // Define another function that is zero everywhere
    std::function<double(const Coord<3> &x)> zero = [] (const mrcpp::Coord<3> &r) -> double {
        return 0.0;
    };
    

    // 4) Let's DEFINE now the POTENTIAL for CORE ELECTRON function V_ee(r) = -1/r
    std::function<double(const Coord<3> &x)> V_ce = [] (const mrcpp::Coord<3> &r) -> double {
        auto R = std::sqrt(r[0]*r[0] + r[1]*r[1] + r[2]*r[2]);
        return -Z/R;
    };

    // Here is the smeared potential
    // Define the smeared potential
    std::function<double(const Coord<3> &x)> smeared_potential = [=] (const mrcpp::Coord<3> &r) -> double {
        auto R = std::sqrt(r[0]*r[0] + r[1]*r[1] + r[2]*r[2]);
        return -SmearedPotential::coulomb_HFYGB(R, Z, building_precision);
    };



    // We now project the potential on the tree
    //mrcpp::project(core_el_tree, V_ce, building_precision);
    mrcpp::project(core_el_tree, smeared_potential, building_precision);


    // 5) PROJECT the function ON the TREE
    mrcpp::project(Trial_function_tree_top, slater, building_precision); // I project the trial-function on the tree
    mrcpp::project(Trial_function_tree_bottom, zero, building_precision); // I project the trial-function on the tree

    double Psi_norm = std::sqrt(Trial_function_tree_top.getSquareNorm() + Trial_function_tree_bottom.getSquareNorm());
    
    // Renormalize the spinor components
    Trial_function_tree_top.rescale(1.0/Psi_norm);
    Trial_function_tree_bottom.rescale(1.0/Psi_norm);

    // Create a vector of CompFunction<3> to hold the components
    std::vector<mrcpp::CompFunction<3>> Psi_2c(2,mra);
    Psi_2c[0] = Trial_function_tree_top;
    Psi_2c[1] = Trial_function_tree_bottom;

    
    if (debug){ 

    // Befofe starting the SCF cycle, let's print some debugging information:
    std::cout << "Here is some debugging information: " << '\n';
    std::cout << "************************************" << '\n';
    std::cout << '\n'<< " Psi_trial (top) square norm = " << Psi_2c[0].getSquareNorm() << '\n';
    std::cout << '\n'<< " Psi_trial (bottom) square norm = " << Psi_2c[1].getSquareNorm() << '\n' << '\n';

    }
    


    std::cout << '\n'<< " Psi_trial (top) square norm = " << Psi_2c[0].getSquareNorm() << '\n';
    std::cout << '\n'<< " Psi_trial (bottom) square norm = " << Psi_2c[1].getSquareNorm() << '\n' << '\n';


   
      /*
     *-----------------------------------
     *-----------------------------------
     *         BEGIN SCF CYCLE
     *-----------------------------------
     *-----------------------------------
    */

    // Parameters for the SCF
    double norm_diff = 1;   //    -> Norm of the difference in the 2 consecutive iterations (initialized as 1 to begin the while loop)
    double E_tot;               //    -> Energy of the system
    double E_Psi;               //    -> Energy of the orbital
    
    // We now define all the trees that will be used to compute the enrgy and the SCF cycle
    mrcpp::CompFunction<3> K_tree(mra);  // -> This is the tree that will hold the K function
    mrcpp::CompFunction<3> K_inverted_tree(mra);  // -> This is the tree that will hold the K^-1 function
    std::vector<mrcpp::CompFunction<3> *> Nabla_K_tree(3, new mrcpp::CompFunction<3>(mra)); // -> This is the tree that will hold the gradient of K
    
    // Gradient of K
    mrcpp::ABGVOperator<3> D(mra, 0.0, 0.0); // deine the ABGV operator
    // Using this to calculate the gradient of K, in the one el case is a constant.
    //Nabla_K_tree = mrcpp::gradient(D, K_tree);
    std::vector<std::vector<mrcpp::CompFunction<3>*>> Nabla_Psi_2c(2, std::vector<mrcpp::CompFunction<3>*>(3, new mrcpp::CompFunction<3>(mra))); // -> This is the tree that will hold the gradient of Psi_2c

    
    // ==================================================================================================================================
    
    
    // Create a vector to hold both Nabla_Psi_t and Nabla_Psi_b
    
    
    std::cout << "************************************************************" << '\n';
    
    
    std::cout << '\n' << '\n';
    // A few utilities variables for the SCF cycle
    std::vector<mrcpp::CompFunction<3>> Psi_2c_next(2, mra);    // -> This will hold the next iteration of the spinor
    mrcpp::CompFunction<3> Psi_2c_diff_top(mra);                // -> This will hold the difference between the 2 spinors for the TOP component
    mrcpp::CompFunction<3> Psi_2c_diff_bottom(mra);             // -> This will hold the difference between the 2 spinors for the BOTTOM component
    
    
    int max_cycle = 9;
    std::cout << "Initializing variables for the SCF cycle...";
    mrcpp::PoissonOperator P(mra, building_precision);
    

    
    Update_SCF_Variables(mra, D, Psi_2c, core_el_tree, Nabla_Psi_2c, K_tree, Nabla_K_tree, K_inverted_tree, Potential_tree, J_function_tree, P);
    std::cout << " done!" << '\n' << '\n';
 // ==================================================================================================================================

 

    std::vector<double> E_values;
    std::vector<double> norm_diff_values;
    std::vector<int> cycle_values;

    


    while (norm_diff > epsilon) {
        num_cycle++;    
        std::cout << "************************************************************" << '\n';
        std::cout << "> CYCLE " << num_cycle << " STARDED:" << '\n' << '\n';

        // Compute the energy
        if (verbose) {
            std::cout << "$ STARTING TO COMPUTE THE ENERGY..." << '\n'<< '\n';
        }
        compute_energy_ZORA(mra, Nabla_Psi_2c, K_tree, Nabla_K_tree, Psi_2c, core_el_tree, J_function_tree, E_tot, E_Psi);
        

        // Apply the Helmholtz operator

        if (verbose) {
            std::cout << "$ STARTING TO COMPUTE NEXT ITERATION'S WAVEFUNCTION..." << '\n'<< '\n';
        }

        
        apply_Helmholtz_ZORA(mra, Psi_2c, Nabla_Psi_2c, Nabla_K_tree, Potential_tree, K_tree, K_inverted_tree, E_Psi , Psi_2c_next);
        

        
        if (debug){
        std::cout << "BEFORE RENORMALIZATION" << '\n';
        std::cout << '\t' << "Psi_2c_next[0] square norm = " << Psi_2c_next[0].getSquareNorm() << '\n';
        std::cout << '\t' << "Psi_2c_next[1] square norm = " << Psi_2c_next[1].getSquareNorm() << '\n';
        }
        if (verbose) {
            std::cout << "$ RENORMALIZING THE WAVEFUNCTION..." << '\n' << '\n';
        }
        // Renormalize the spinor Psi_2c_next
        Renormalize_Spinor(Psi_2c_next);
        
        if (debug){
        std::cout << "AFTER RENORMALIZATION" << '\n';
        std::cout << '\t' << "Psi_2c_next[0] square norm = " << Psi_2c_next[0].getSquareNorm() << '\n';
        std::cout << '\t' << "Psi_2c_next[1] square norm = " << Psi_2c_next[1].getSquareNorm() << '\n';
        
        // Debug, check Psi_2c norm:
        std::cout << "Norm of the spinor" << '\n';
        std::cout << '\t' << "Psi_2c[0] square norm = " << Psi_2c[0].getSquareNorm() << '\n';
        std::cout << '\t' << "Psi_2c[1] square norm = " << Psi_2c[1].getSquareNorm() << '\n';
        }
        

      
        if (verbose) {
            std::cout << "$ COMPUTING THE DIFFERENCE BETWEEN THE 2 ITERATIONS..." << '\n' << '\n';
        }
        // Compute the difference between the 2 spinors
        mrcpp::add(Psi_2c_diff_top, 1.0, Psi_2c[0] , -1.0, Psi_2c_next[0] ,building_precision, false);
        mrcpp::add(Psi_2c_diff_bottom, 1.0, Psi_2c[1] , -1.0, Psi_2c_next[1] ,building_precision, false);

        if (debug){
        std::cout << "Norm of the difference" << '\n';
        std::cout << '\t' << "Psi_2c_diff[0] square norm = " << Psi_2c_diff_top.getSquareNorm() << '\n';
        std::cout << '\t' << "Psi_2c_diff[1] square norm = " << Psi_2c_diff_bottom.getSquareNorm() << '\n';
        }

        

        norm_diff = std::sqrt(Psi_2c_diff_top.getSquareNorm() + Psi_2c_diff_bottom.getSquareNorm());
        std::cout << "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@" << '\n';
        std::cout << "@ Norm of the difference = " << norm_diff << " @" << '\n';
        std::cout << "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@" << '\n';
        std::cout << '\n';


        if (verbose) {
            std::cout << "$ UPDATING THE VARIABLES FOR THE NEXT ITERATION..." << '\n' << '\n';
        }
        
        
       
        std::cout << "SCF CYCLE COMPLETED!" << '\n' << '\n';


        
        Update_SCF_Variables(mra, D, Psi_2c_next, core_el_tree, Nabla_Psi_2c, K_tree, Nabla_K_tree, K_inverted_tree, Potential_tree, J_function_tree, P);
        std::cout<< "SCF variables updated!" << '\n' << '\n';
        // Update the Psi_2c
        Psi_2c.swap(Psi_2c_next);
        //Psi_2c_next.clear();
        //Psi_2c_diff.clear();

        // Store all the values for the SCF cycle
        E_values.push_back(E_tot);
        norm_diff_values.push_back(norm_diff);
        cycle_values.push_back(num_cycle);

        std::cout << "Last cleanup, MOVING ON TO THE NEXT CYCLE!" << '\n' << '\n';



        if (num_cycle >= max_cycle){
            std::cout << "The SCF cycle did not converge" << '\n';
            print::footer(0, timer, 2);
            exit(0);
        }
        


    }



    std::cout << "************************************************************" << '\n';


    std::cout << '\n' << '\n';
    std::cout << "-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-" << '\n';
    std::cout << "Your calculation is finished, HURRAY. Here are some cute cats:"<<'\n' << '\n';
    std::cout << "=^.^=" << '\t' <<"=^.^=" << '\t' <<"=^.^=" << '\t' <<"=^.^=" << '\t' <<"=^.^=" << '\t' <<"=^.^=" << '\t' <<"=^.^=" << '\t' << '\n' << '\n';

    std::cout << "Here is a summary of the SCF cycle:" << '\n' << '\n';
    std::cout << "Cycle" << '\t' << "Energy" << '\t' << '\t' << '\t' << "Norm. Diff." << '\n';
    for (int i = 0; i < E_values.size(); i++){
        std::cout << cycle_values[i] << '\t' << E_values[i] << '\t' << norm_diff_values[i] << '\n';
    }
    std::cout << '\n' << '\n' << "Non relativistic energuy for He atom" << '\n';
    std::cout << "*************************************************************" << '\n';
    compute_He_energy(mra, *He_NR_Solution.CompD[0], *core_el_tree.CompD[0], mu, E_tot);
    std::cout << "*************************************************************" << '\n';

    print::footer(0, timer, 2);
    return 0;
}