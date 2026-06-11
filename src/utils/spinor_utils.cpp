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

    /** @brief in-place application of Pauli/Gamma matrices to a CompFunction 
     *  @param inp: input (and output) CompFunction, passed by reference, will be modified
     *  @param index: integer index of the gamma matrix to be applied. 0 is the identity. 5 is not yet implemented
     *  
     * WARNING: 4C/Gamma matrices are not implemented. Only 2C/Pauli matrices are, currently.
     */
    template<int D> void apply_gamma(CompFunction<D> &inp, int index) {
        ComplexDouble comp_i(0.0, 1.0); // Define the imaginary unit for convenience later
        switch (index) {
        case 0:
            //Identity, base case, nothing to apply
            break;
        case 1:
            // Apply sigma_X matrix
            // Basically amounts to swapping first and second components
            for (int i = 0; i < inp.Ncomp(); i = i + 2) {
                //swapping trees
                std::swap(inp.CompD[i], inp.CompD[i+1]);
                std::swap(inp.CompC[i], inp.CompC[i+1]);
                //swapping tree metadata
                std::swap(inp.func_ptr->data.Nchunks[i], inp.func_ptr->data.Nchunks[i+1]);
                //swapping prefactors
                std::swap(inp.func_ptr->data.c1[i], inp.func_ptr->data.c1[i+1]);
            }
            break;
        case 2:
            // Apply sigma_Y matrix
            // Amounts to swapping first and second components, and multiplying former by i and the latter by -i
            for (int i = 0; i < inp.Ncomp(); i = i + 2) {
                //swapping trees
                std::swap(inp.CompD[i], inp.CompD[i+1]);
                std::swap(inp.CompC[i], inp.CompC[i+1]);
                //swapping tree metadata
                std::swap(inp.func_ptr->data.Nchunks[i], inp.func_ptr->data.Nchunks[i+1]);
                //swapping prefactors
                std::swap(inp.func_ptr->data.c1[i], inp.func_ptr->data.c1[i+1]);
                //multiplying the 1st and 2nd prefactors by the complex unit ±i
                inp.func_ptr->data.c1[i] *= (-1.0*comp_i);
                inp.func_ptr->data.c1[i+1] *= comp_i;
            }
            break;
        case 3:
            // Apply Alpha-Z matrix
            // nothing to do except multiply the second element's prefactor by -1
            for (int i = 0; i < inp.Ncomp(); i = i + 2) {
                inp.func_ptr->data.c1[i+1] *= -1.0;
            }
            break;
        default:
            //identity case again. Nothing to apply.
            break;
        }
    }
    
    /*
     * @brief: shuffles the indices of a spinor, simulating the application of a Dirac matrix to it
     * pauli represents the index of the Dirac matrices.
     * For scalar operators, it is unused.
     * For 2 component (Weyl/Pauli) spinors, pauli = 0,1,2,3 corresponds to indentiy, sigma_x, y and z respectively.
     * 
    */
    //TODO: Add a provision in case inp and out are identical 
    void apply_Pauli(CompFunction<3> &out, CompFunction<3> &inp, int pauli, double prec, bool conjugate) { //NOTE: assumes 2-component spinors for now
        // Implementation of applying Pauli matrices to spinor functions
        // This function will modify 'out' based on the Pauli matrix specified by 'pauli'
        // and the input function 'inp'.
        // The 'prec' parameter is used for precision control.
        // The 'conjugate' parameter indicates whether to apply conjugation.
        ComplexDouble comp_i = {0.0, 1.0}; // Define the imaginary unit

        switch (pauli) {
        case 0:
            //Identity, base case, nothing to apply, just copy the input to the output if they are not the same function
            if (&out != &inp) {
                // out.deep_copy(inp);
                out = inp;
            }
            break;
        case 1:
            // Apply Pauli-X matrix
            for (int i = 0; i < inp.Ncomp(); i = i + 2) {
                if (inp.isreal() == 1) {
                    out.defreal(); // Safety catch for later operations on out. If the input is real, the output should also be defined as real
                    out.func_ptr->data.iscomplex = 0;
                    inp.CompD[i]->deep_copy(out.CompD[i+1]);
                    inp.CompD[i+1]->deep_copy(out.CompD[i]);
                } else {
                    out.defcomplex(); // Safety catch for later operations on out. If the input is complex, the output should also be defined as complex
                    out.func_ptr->data.isreal = 0;
                    inp.CompC[i]->deep_copy(out.CompC[i+1]);
                    inp.CompC[i+1]->deep_copy(out.CompC[i]);;
                }
                // coefficient multiplication needn't be done separately for real and complex cases, since the coefficient is purely imaginary, so it will just be multiplied to the complex part of the function, even if the function is defined as real. However, we need to make sure that the output function is defined as complex in this case, otherwise we might run into issues later on when trying to multiply it by a complex coefficient. 
                out.func_ptr->data.c1[i] = inp.func_ptr->data.c1[i+1];
                out.func_ptr->data.c1[i+1] = inp.func_ptr->data.c1[i];
            }
            break;
        case 2: //WARNING: We may need to enforce out to be complex in this case, rather than just multiplying the whole phase by i
            // Apply Pauli-Y matrix
            for (int i = 0; i < inp.Ncomp(); i = i + 2) {
                if (inp.isreal() == 1) {
                    out.defreal(); // Safety catch for later operations on out. If the input is real, the output should also be defined as real
                    out.func_ptr->data.iscomplex = 0;
                    inp.CompD[i]->deep_copy(out.CompD[i+1]);
                    inp.CompD[i+1]->deep_copy(out.CompD[i]);


                } else {
                    out.defcomplex(); // Safety catch for later operations on out. If the input is complex, the output should also be defined as complex
                    out.func_ptr->data.isreal = 0;
                    inp.CompC[i]->deep_copy(out.CompC[i+1]);
                    inp.CompC[i+1]->deep_copy(out.CompC[i]);
                }
                // coefficient multiplication needn't be done separately for real and complex cases, since the coefficient is purely imaginary, so it will just be multiplied to the complex part of the function, even if the function is defined as real. However, we need to make sure that the output function is defined as complex in this case, otherwise we might run into issues later on when trying to multiply it by a complex coefficient.
                out.func_ptr->data.c1[i] = inp.func_ptr->data.c1[i+1] * (-1.0)*comp_i;
                out.func_ptr->data.c1[i+1] = inp.func_ptr->data.c1[i] * comp_i;
            }
            break;
        case 3:
            // Apply Pauli-Z matrix
            for (int i = 0; i < inp.Ncomp(); i = i + 2) {
                if (inp.isreal() == 1) {
                    out.defreal(); // Safety catch for later operations on out. If the input is real, the output should also be defined as real
                    out.func_ptr->data.iscomplex = 0;
                    inp.CompD[i]->deep_copy(out.CompD[i]);
                    inp.CompD[i+1]->deep_copy(out.CompD[i+1]);
                } else {
                    out.defcomplex(); // Safety catch for later operations on out. If the input is complex, the output should also be defined as complex
                    out.func_ptr->data.isreal = 0;
                    inp.CompC[i]->deep_copy(out.CompC[i]);
                    inp.CompC[i+1]->deep_copy(out.CompC[i+1]);
                }
                out.func_ptr->data.c1[i] = inp.func_ptr->data.c1[i];
                out.func_ptr->data.c1[i+1] = inp.func_ptr->data.c1[i+1] * (-1.0);
            }
            break;
        default:
            // std::cerr << "Invalid Pauli matrix index, values must be 0,1,2,3. Current value: " << pauli << std::endl;
            MSG_ABORT("Invalid Pauli matrix index, values must be 0,1,2,3. Current value: " << pauli);
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

    template void apply_gamma(CompFunction<3> &inp, int index);
}
