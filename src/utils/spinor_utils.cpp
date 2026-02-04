#include "spinor_utils.h"
#include "utils/CompFunction.h"
#include "utils/mpi_utils.h"
#include "utils/parallel.h"
// #include "FunctionTreeVector.h"

#include <complex>
#include <iostream>
using namespace std::complex_literals;

namespace mrcpp {
    
    // /*
    //  * @brief: shuffles the indices of a spinor, simulating the application of a Dirac matrix to it
    //  * pauli represents the index of the Dirac matrices.
    //  * For scalar operators, it is unused.
    //  * For 2 component (Weyl/Pauli) spinors, pauli = 0,1,2,3 corresponds to indentiy, sigma_x, y and z respectively.
    //  * 
    // */
    void apply_Pauli(CompFunction<3> &out, const CompFunction<3> &inp, int pauli, double prec, bool conjugate) { //NOTE: assumes 2-component spinors for now
        // Implementation of applying Pauli matrices to spinor functions
        // This function will modify 'out' based on the Pauli matrix specified by 'pauli'
        // and the input function 'inp'.
        // The 'prec' parameter is used for precision control.
        // The 'conjugate' parameter indicates whether to apply conjugation.
        ComplexDouble comp_i(0.0, 1.0); // Define the imaginary unit
        switch (pauli) {
        case 0:
            // Apply Pauli-X matrix
            for (int i = 0; i < inp.Ncomp(); i = i + 2) {
                out.setReal(inp.CompD[i+1], i); //WARNING: No copy might create some issues down the line
                out.setCplx(inp.CompC[i+1], i); //WARNING: No copy might create some issues down the line
                out.func_ptr->data.c1[i] *= inp.func_ptr->data.c1[i+1];
                // out.CompC[i] = inp.CompC[i+1]; //WARNING: No copy might create some issues down the line
                // out.CompC[i+1] = inp.CompC[i];
                out.setReal(inp.CompD[i], i+1); //WARNING: No copy might create some issues down the line
                out.setCplx(inp.CompC[i], i+1); //WARNING: No copy might create some issues down the line
                out.func_ptr->data.c1[i+1] *= inp.func_ptr->data.c1[i];
            }
            break;
        case 1:
            // Apply Pauli-Y matrix
            // std::cout << "apply Pauli Y tut0 " << out.Ncomp() << " " << inp.Ncomp() << '\n';
            for (int i = 0; i < inp.Ncomp(); i = i + 2) {
                // ComplexDouble dotut0 = dot(out, out);
                std::cout << "Applying Pauli-Y matrix " << i << std::endl;
                // out.CompC[i] = inp.CompC[i+1]; //WARNING: No copy might create some issues down the line
                out.setReal(inp.CompD[i+1], i); //WARNING: No copy might create some issues down the line
                out.setCplx(inp.CompC[i+1], i); //WARNING: No copy might create some issues down the line
                // std::cout << "apply Pauli Y tut1" << '\n'; 
                // out.CompC[i]->rescale(-1.0i); //Là j'utilise le rescale des functiontree, pas celui de CompFunction TODO: créer un rescale qui ne change que un seul terme
                out.func_ptr->data.c1[i] *= inp.func_ptr->data.c1[i+1] * (-1.0)*comp_i;
                // out.func_ptr->data.c1[i] *= comp_i;
                // std::cout << "apply Pauli Y tut2" << '\n';
                // out.CompC[i+1] = inp.CompC[i];
                out.setReal(inp.CompD[i], i+1); //WARNING: No copy might create some issues down the line
                out.setCplx(inp.CompC[i], i+1); //WARNING: No copy might create some issues down the line
                // std::cout << "apply Pauli Y tut3" << '\n';
                // out.CompC[i+1]->rescale(1.0i);
                out.func_ptr->data.c1[i+1] *= inp.func_ptr->data.c1[i] * comp_i;
                // out.func_ptr->data.c1[i+1] *= -1.0*comp_i;
                // std::cout << "apply Pauli Y tut4 pouet" << '\n';

                // rescale(1.0, out.CompC[i], -1.0i, inp.CompC[i+1]); //WARNING: No copy might create some issues down the line
                // out.CompC[i+1]->multiply(1.0i, inp.CompC[i]);
                // ComplexDouble dotuta = dot(inp, inp);
                // std::cout << "apply Pauli Y tut5" << dotuta << '\n';
                // ComplexDouble dotutb = dot(out, out);
                // std::cout << "apply Pauli Y tut6" << dotutb << '\n';
                // ComplexDouble dotutc = dot(out, inp);
                // std::cout << "apply Pauli Y tut7" << dotutc << '\n';
            }
            break;
        case 2:
            // Apply Pauli-Z matrix
            for (int i = 0; i < inp.Ncomp(); i = i + 2) {
                // out.CompC[i]->multiply(-1.0i, inp.CompC[i+1]); //WARNING: No copy might create some issues down the line
                out.setReal(inp.CompD[i], i); //WARNING: No copy might create some issues down the line
                out.setCplx(inp.CompC[i], i); //WARNING: No copy might create some issues down the line
                // out.CompC[i] = inp.CompC[i]; //WARNING: No copy might create some issues down the line
                out.func_ptr->data.c1[i] *= inp.func_ptr->data.c1[i];
                // out.CompC[i+1] = inp.CompC[i+1]; //WARNING: No copy might create some issues down the line
                // out.CompC[i+1]->rescale(-1.0);
                out.func_ptr->data.c1[i+1] *= inp.func_ptr->data.c1[i+1] * (-1.0);
                // out.func_ptr->data.c1[i+1] *= -1.0;
            }
            break;
        default:
            std::cerr << "Invalid Pauli matrix index, values must be 0,1,2,3. Current value: " << pauli << std::endl;
        }
    }
    
    // void dot_spinor(CompFunction<3> &out, const CompFunction<3> &inp_a, const CompFunction<3> &inp_b, double prec, bool conjugate) {
    //     // Implementation of the dot product for spinor functions
    //     // This function computes the dot product of two spinor functions and stores the result in 'out'.
    //     if (inp_a.Ncomp() != inp_b.Ncomp()) {
    //         std::cerr << "dot_spinor: Input functions must have the same number of components." << std::endl;
    //         return;
    //     }
    //     if (out.Ncomp() != 1) {
    //         std::cerr << "dot_spinor: Output function must have exactly one component." << std::endl;
    //         return;
    //     }
    //     // Compute the dot product
    //     // ComplexDouble result = 0.0;
    //     for (int i = 0; i < inp_a.Ncomp(); ++i) {
    //         auto tmp_add = 
    //         // ComplexDouble val_a = inp_a.CompC[i]->integrate();
    //         // ComplexDouble val_b = conjugate ? std::conj(inp_b.CompC[i]->integrate()) : inp_b.CompC[i]->integrate();
    //         // result += val_a * val_b;
    //     }
    //     out.CompC[0]->setValue(result);
    // }

    void normalize_spinor(CompFunction<3> &inp, double prec) {
        // Implementation of normalization for spinor functions
        // This function normalizes the input spinor function 'inp' in place.
        // if (inp.Ncomp() != 2 && inp.Ncomp() != 4) {
        //     std::cerr << "normalize_spinor: Input function must have 2 or 4 components." << std::endl;
        //     return;
        // }
        double norm = inp.norm();
        if (norm < prec) {
            std::cerr << "normalize_spinor: Norm is too small, cannot normalize." << std::endl;
            return;
        }
        // inp.rescale(1.0 / norm);
        for (int i = 0; i < inp.Ncomp(); i++) {
            inp.func_ptr->data.c1[i] *= 1.0 / norm; // Normalize each component
            // inp.CompC[i]->rescale(1.0 / norm); //WARNING: No copy might create some issues down the line
            // inp.CompC[i]->func_ptr->data.c1[i] *= 1.0 / norm; // Normalize each component
            // std::cout << "normalize_spinor: Component " << i << " normalized." << std::endl;
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
