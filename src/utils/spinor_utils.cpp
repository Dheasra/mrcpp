#include "spinor_utils.h"
#include "utils/CompFunction.h"
#include "utils/mpi_utils.h"
#include "utils/parallel.h"
// #include "FunctionTreeVector.h"

#include <complex>
#include <iostream>

// remove below
#include "Printer.h"
// till here

using namespace std::complex_literals;

namespace mrcpp {

    CompFunction<3> apply_alpha(const CompFunction<3> &inp, int alpha, bool conjugate) {
        // Implementation of applying Alpha matrices to spinor functions
        // This function will modify 'out' based on the Alpha matrix specified by 'Alpha'
        // and the input function 'inp'.
        // The 'prec' parameter is used for precision control.
        // The 'conjugate' parameter indicates whether to apply conjugation.
        ComplexDouble comp_i(0.0, 1.0); // Define the imaginary unit
        CompFunction<3> out(inp);
        switch (alpha) {
        case 0:
            //Identity, base case, nothing to apply
            break;
        case 1:
            // Apply Alpha-X matrix
            for (int i = 0; i < inp.Ncomp(); i = i + 2) {
                if (inp.isreal() == 1) {
                    out.defreal(); // Safety catch for later operations on out. If the input is real, the output should also be defined as real
                    // Warning: shallow copy
                    out.CompD[i] = inp.CompD[i+1];
                    out.CompD[i+1] = inp.CompD[i];
                    // out.setReal(inp.CompD[i+1], i);
                    // out.setReal(inp.CompD[i], i+1);

                } else {
                    out.defcomplex(); // Safety catch for later operations on out. If the input is complex, the output should also be defined as complex
                    // Warning: shallow copy
                    out.CompC[i] = inp.CompC[i+1];
                    out.CompC[i+1] = inp.CompC[i];
                    // out.setCplx(inp.CompC[i+1], i);
                    // out.setCplx(inp.CompC[i], i+1);
                }
                // coefficient multiplication needn't be done separately for real and complex cases, since the coefficient is purely imaginary, so it will just be multiplied to the complex part of the function, even if the function is defined as real. However, we need to make sure that the output function is defined as complex in this case, otherwise we might run into issues later on when trying to multiply it by a complex coefficient. 
                out.func_ptr->data.c1[i] = inp.func_ptr->data.c1[i+1];
                out.func_ptr->data.c1[i+1] = inp.func_ptr->data.c1[i];
            }
            break;
        case 2:
            // Apply Alpha-Y matrix
            // std::cout << "apply Alpha Y tut0 " << out.Ncomp() << " " << inp.Ncomp() << '\n';
            for (int i = 0; i < inp.Ncomp(); i = i + 2) {
                // std::cout << "Applying Alpha-Y matrix " << i << " " << inp.isreal() << std::endl;
                if (inp.isreal() == 1) {
                    out.defreal(); // Safety catch for later operations on out. If the input is real, the output should also be defined as real
                    // Warning: shallow copy
                    out.setReal(inp.CompD[i+1], i);
                    out.setReal(inp.CompD[i], i+1);
                    // out.CompD[i] = inp.CompD[i+1];
                    // out.CompD[i+1] = inp.CompD[i];

                } else {
                    out.defcomplex(); // Safety catch for later operations on out. If the input is complex, the output should also be defined as complex
                    // Warning: shallow copy
                    out.setCplx(inp.CompC[i+1], i);
                    out.setCplx(inp.CompC[i], i+1);
                    // out.CompC[i] = inp.CompC[i+1];
                    // out.CompC[i+1] = inp.CompC[i];
                }
                // coefficient multiplication needn't be done separately for real and complex cases, since the coefficient is purely imaginary, so it will just be multiplied to the complex part of the function, even if the function is defined as real. However, we need to make sure that the output function is defined as complex in this case, otherwise we might run into issues later on when trying to multiply it by a complex coefficient.
                out.func_ptr->data.c1[i] = inp.func_ptr->data.c1[i+1] * (-1.0)*comp_i;
                out.func_ptr->data.c1[i+1] = inp.func_ptr->data.c1[i] * comp_i;
            }
            break;
        case 3:
            // Apply Alpha-Z matrix
            for (int i = 0; i < inp.Ncomp(); i = i + 2) {
                if (inp.isreal() == 1) {
                    out.defreal(); // Safety catch for later operations on out. If the input is real, the output should also be defined as real
                    // Warning: shallow copy
                    out.setReal(inp.CompD[i], i); 
                    out.setReal(inp.CompD[i+1], i+1); 
                    // out.CompD[i] = inp.CompD[i+1];
                    // out.CompD[i+1] = inp.CompD[i];

                } else {
                    out.defcomplex(); // Safety catch for later operations on out. If the input is complex, the output should also be defined as complex
                    // Warning: shallow copy
                    out.setCplx(inp.CompC[i], i); 
                    out.setCplx(inp.CompC[i+1], i+1); 
                    // out.CompC[i] = inp.CompC[i+1]; 
                    // out.CompC[i+1] = inp.CompC[i];
                }
                out.func_ptr->data.c1[i] = inp.func_ptr->data.c1[i];
                // out.CompC[i+1] = inp.CompC[i+1]; 
                // out.CompC[i+1]->rescale(-1.0);
                out.func_ptr->data.c1[i+1] = inp.func_ptr->data.c1[i+1] * (-1.0);
                // out.func_ptr->data.c1[i+1] *= -1.0;
            }
            break;
        default:
            std::cerr << "Invalid Alpha matrix index, values must be 0,1,2,3. Current value: " << alpha << std::endl;
        }
        return out;
    }
    
    // /*
    //  * @brief: shuffles the indices of a spinor, simulating the application of a Dirac matrix to it
    //  * pauli represents the index of the Dirac matrices.
    //  * For scalar operators, it is unused.
    //  * For 2 component (Weyl/Pauli) spinors, pauli = 0,1,2,3 corresponds to indentiy, sigma_x, y and z respectively.
    //  * 
    // */
    //TODO: changer ça pour que ça retourne un CompFunction et retier out? Ou faire une provision si out et inp sont identiques?
    void apply_Pauli(CompFunction<3> &out, CompFunction<3> &inp, int pauli, double prec, bool conjugate) { //NOTE: assumes 2-component spinors for now
        // Implementation of applying Pauli matrices to spinor functions
        // This function will modify 'out' based on the Pauli matrix specified by 'pauli'
        // and the input function 'inp'.
        // The 'prec' parameter is used for precision control.
        // The 'conjugate' parameter indicates whether to apply conjugation.
        ComplexDouble comp_i(0.0, 1.0); // Define the imaginary unit

        switch (pauli) {
        case 0:
            //Identity, base case, nothing to apply, just copy the input to the output if they are not the same function
            if (&out != &inp) {
                // out.deep_copy(inp);
                out = inp;
            }
            break;
        case 1:
            // MSG_INFO("case X "<< inp.Ncomp() << " " << out.Ncomp());
            // Apply Pauli-X matrix
            for (int i = 0; i < inp.Ncomp(); i = i + 2) {
                if (inp.isreal() == 1) {
                    out.defreal(); // Safety catch for later operations on out. If the input is real, the output should also be defined as real
                    out.func_ptr->data.iscomplex = 0;
                    // Warning: shallow copy
                    // out.CompD[i] = inp.CompD[i+1];
                    // out.CompD[i+1] = inp.CompD[i];
                    // out.setReal(new FunctionTree<3, double> (*inp.CompD[i+1]), i);
                    // out.setReal(new FunctionTree<3, double> (*inp.CompD[i]), i+1);
                    // MSG_INFO(i << " comp addresses input" << &inp.CompD[i] << " " << &inp.CompD[i+1] << " comp addresses output" << &out.CompD[i] << " " << &out.CompD[i+1]);
                    inp.CompD[i]->deep_copy(out.CompD[i+1]);
                    inp.CompD[i+1]->deep_copy(out.CompD[i]);
                    // out.setReal(inp.CompD[i+1], i);
                    // out.setReal(inp.CompD[i], i+1);

                } else {
                    out.defcomplex(); // Safety catch for later operations on out. If the input is complex, the output should also be defined as complex
                    out.func_ptr->data.isreal = 0;
                    // Warning: shallow copy
                    // out.CompC[i] = inp.CompC[i+1];
                    // out.CompC[i+1] = inp.CompC[i];
                    inp.CompC[i]->deep_copy(out.CompC[i+1]);
                    inp.CompC[i+1]->deep_copy(out.CompC[i]);
                    // out.setCplx(inp.CompC[i+1], i);
                    // out.setCplx(inp.CompC[i], i+1);
                }
                // coefficient multiplication needn't be done separately for real and complex cases, since the coefficient is purely imaginary, so it will just be multiplied to the complex part of the function, even if the function is defined as real. However, we need to make sure that the output function is defined as complex in this case, otherwise we might run into issues later on when trying to multiply it by a complex coefficient. 
                out.func_ptr->data.c1[i] = inp.func_ptr->data.c1[i+1];
                out.func_ptr->data.c1[i+1] = inp.func_ptr->data.c1[i];
            }
            break;
        case 2: //WARNING: We may need to enforce out to be complex in this case, rather than just multiplying the whole phase by i
            // MSG_INFO("Y ");
            // Apply Pauli-Y matrix
            // std::cout << "apply Pauli Y tut0 " << out.Ncomp() << " " << inp.Ncomp() << '\n';
            for (int i = 0; i < inp.Ncomp(); i = i + 2) {
                // std::cout << "Applying Pauli-Y matrix " << i << " " << inp.isreal() << std::endl;
                if (inp.isreal() == 1) {
                    out.defreal(); // Safety catch for later operations on out. If the input is real, the output should also be defined as real
                    out.func_ptr->data.iscomplex = 0;
                    // Warning: shallow copy
                    // out.setReal(inp.CompD[i+1], i);
                    // out.setReal(inp.CompD[i], i+1);
                    inp.CompD[i]->deep_copy(out.CompD[i+1]);
                    inp.CompD[i+1]->deep_copy(out.CompD[i]);
                    // out.CompD[i] = inp.CompD[i+1];
                    // out.CompD[i+1] = inp.CompD[i];

                } else {
                    out.defcomplex(); // Safety catch for later operations on out. If the input is complex, the output should also be defined as complex
                    out.func_ptr->data.isreal = 0;
                    // Warning: shallow copy
                    // out.setCplx(inp.CompC[i+1], i);
                    // out.setCplx(inp.CompC[i], i+1);
                    inp.CompC[i]->deep_copy(out.CompC[i+1]);
                    inp.CompC[i+1]->deep_copy(out.CompC[i]);
                    // out.CompC[i] = inp.CompC[i+1];
                    // out.CompC[i+1] = inp.CompC[i];
                }
                // coefficient multiplication needn't be done separately for real and complex cases, since the coefficient is purely imaginary, so it will just be multiplied to the complex part of the function, even if the function is defined as real. However, we need to make sure that the output function is defined as complex in this case, otherwise we might run into issues later on when trying to multiply it by a complex coefficient.
                out.func_ptr->data.c1[i] = inp.func_ptr->data.c1[i+1] * (-1.0)*comp_i;
                out.func_ptr->data.c1[i+1] = inp.func_ptr->data.c1[i] * comp_i;
            }
            break;
        case 3:
            // MSG_INFO("Z ");
            // Apply Pauli-Z matrix
            for (int i = 0; i < inp.Ncomp(); i = i + 2) {
                if (inp.isreal() == 1) {
                    out.defreal(); // Safety catch for later operations on out. If the input is real, the output should also be defined as real
                    out.func_ptr->data.iscomplex = 0;
                    // Warning: shallow copy
                    // out.setReal(inp.CompD[i], i); 
                    // out.setReal(inp.CompD[i+1], i+1); 
                    inp.CompD[i]->deep_copy(out.CompD[i]);
                    inp.CompD[i+1]->deep_copy(out.CompD[i+1]);
                    // out.CompD[i] = inp.CompD[i+1];
                    // out.CompD[i+1] = inp.CompD[i];

                } else {
                    out.defcomplex(); // Safety catch for later operations on out. If the input is complex, the output should also be defined as complex
                    out.func_ptr->data.isreal = 0;
                    // Warning: shallow copy
                    // out.setCplx(inp.CompC[i], i); 
                    // out.setCplx(inp.CompC[i+1], i+1); 
                    inp.CompC[i]->deep_copy(out.CompC[i]);
                    inp.CompC[i+1]->deep_copy(out.CompC[i+1]);
                    // out.CompC[i] = inp.CompC[i+1]; 
                    // out.CompC[i+1] = inp.CompC[i];
                }
                out.func_ptr->data.c1[i] = inp.func_ptr->data.c1[i];
                // out.CompC[i+1] = inp.CompC[i+1]; 
                // out.CompC[i+1]->rescale(-1.0);
                out.func_ptr->data.c1[i+1] = inp.func_ptr->data.c1[i+1] * (-1.0);
                // out.func_ptr->data.c1[i+1] *= -1.0;
            }
            break;
        default:
            std::cerr << "Invalid Pauli matrix index, values must be 0,1,2,3. Current value: " << pauli << std::endl;
        }
    }
    
    void normalize_spinor(CompFunction<3> &inp, double prec) {
        // Implementation of normalization for spinor functions
        // This function normalizes the input spinor function 'inp' in place.
 
        double norm = inp.norm();
        if (norm < prec) {
            std::cerr << "normalize_spinor: Norm is too small, cannot normalize." << std::endl;
            return;
        }
        // inp.rescale(1.0 / norm);
        if (inp.isreal() == 1) {
            for (int i = 0; i < inp.Ncomp(); i++) {
                inp.CompD[i]->rescale(1.0 / norm); //Rescaling each component in place
                inp.func_ptr->data.c1[i] = 1.0; // Resetting the overall multiplicative factor to 1 after normalization, since the components have already been rescaled.
                // std::cout << "normalize_spinor: Component " << i << " normalized." << std::endl;
            }
        } else {
            for (int i = 0; i < inp.Ncomp(); i++) {
                inp.CompC[i]->rescale(1.0 / norm); //Rescaling each component in place
                inp.func_ptr->data.c1[i] = ComplexDouble(1.0,0.0); // Resetting the overall multiplicative factor to 1 after normalization, since the components have already been rescaled.
                // std::cout << "normalize_spinor: Component " << i << " normalized." << std::endl;
            }
        }
    }

    // @brief Disjoining (filetering out) scalar paired orbitals in the vector of spinors Phi
    // The elements of Phi that have n1[spin]==1 are moved to the output vector, while the others remain in Phi, with ownership transferred as needed.
    // @param Phi: vector of spinors to be disjoined
    // @param spin: index of the spin component to filter by (0 for alpha, 1 for beta) (For scalar and 2C calculations). For 4C spinors, we also have 2 is alpha small component and 3 is beta small component. 
    // CompFunctionVector disjoin(CompFunctionVector &Phi, int spin) {
    //     CompFunctionVector out;
    //     CompFunctionVector tmp;
    //     for (auto &i : Phi) {
    //         if (i.func_ptr->data.n1[spin] == 1) { //checking if the element's spin is the desired one and transferring it to out (with ownership)
    //             if (i.getRank() % mrcpp::mpi::wrk_size != out.size() % mrcpp::mpi::wrk_size) { 
    //                 // need to send orbital from owner to new owner
    //                 if (mrcpp::mpi::my_func(i)) { mrcpp::mpi::send_function(i, out.size() % mrcpp::mpi::wrk_size, i.getRank(), mrcpp::mpi::comm_wrk); }
    //                 if (mrcpp::mpi::my_func(out.size())) { mrcpp::mpi::recv_function(i, i.getRank() % mrcpp::mpi::wrk_size, i.getRank(), mrcpp::mpi::comm_wrk); }
    //             }
    //             i.setRank(out.size());
    //             out.push_back(i);
    //         } else { //otherwise transferring it to tmp, also with ownership.
    //             if (i.getRank() % mrcpp::mpi::wrk_size != tmp.size() % mrcpp::mpi::wrk_size) {
    //                 // need to send orbital from owner to new owner
    //                 if (mrcpp::mpi::my_func(i)) { mrcpp::mpi::send_function(i, tmp.size() % mrcpp::mpi::wrk_size, i.getRank(), mrcpp::mpi::comm_wrk); }
    //                 if (mrcpp::mpi::my_func(tmp.size())) { mrcpp::mpi::recv_function(i, i.getRank() % mrcpp::mpi::wrk_size, i.getRank(), mrcpp::mpi::comm_wrk); }
    //             }
    //             i.setRank(tmp.size());
    //             tmp.push_back(i);
    //         }
    //     }
    //     Phi.clear();
    //     Phi = tmp;
    //     return out;
    // }
    
    // CompFunctionVector adjoin(CompFunctionVector &Phi_a, CompFunctionVector &Phi_b) {
    //     CompFunctionVector out;
    //     for (auto &phi : Phi_a) {
    //         if (phi.getRank() % mrcpp::mpi::wrk_size != out.size() % mrcpp::mpi::wrk_size) {
    //             // need to send orbital from owner to new owner
    //             if (mrcpp::mpi::my_func(phi)) { mrcpp::mpi::send_function(phi, out.size() % mrcpp::mpi::wrk_size, phi.getRank(), mrcpp::mpi::comm_wrk); }
    //             if (mrcpp::mpi::my_func(out.size())) { mrcpp::mpi::recv_function(phi, phi.getRank() % mrcpp::mpi::wrk_size, phi.getRank(), mrcpp::mpi::comm_wrk); }
    //         }
    //         phi.setRank(out.size());
    //         out.push_back(phi);
    //     }
    //     for (auto &phi : Phi_b) {
    //         if (phi.getRank() % mrcpp::mpi::wrk_size != out.size() % mrcpp::mpi::wrk_size) {
    //             // need to send orbital from owner to new owner
    //             if (mrcpp::mpi::my_func(phi)) { mrcpp::mpi::send_function(phi, out.size() % mrcpp::mpi::wrk_size, phi.getRank(), mrcpp::mpi::comm_wrk); }
    //             if (mrcpp::mpi::my_func(out.size())) { mrcpp::mpi::recv_function(phi, phi.getRank() % mrcpp::mpi::wrk_size, phi.getRank(), mrcpp::mpi::comm_wrk); }
    //         }
    //         phi.setRank(out.size());
    //         out.push_back(phi);
    //     }
    //     Phi_a.clear();
    //     Phi_b.clear();
    //     return out;
    // }
}
