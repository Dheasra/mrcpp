#include "spinor_utils.h"

#include "CompFunction.h"
#include "FunctionTreeVector.h"

#include <complex>
#include <iostream>

namespace mrcpp {
    void apply_Pauli(CompFunction<3> &out, const CompFunction<3> &inp, int pauli, double prec, bool conjugate) {
        // Implementation of applying Pauli matrices to spinor functions
        // This function will modify 'out' based on the Pauli matrix specified by 'pauli'
        // and the input function 'inp'.
        // The 'prec' parameter is used for precision control.
        // The 'conjugate' parameter indicates whether to apply conjugation.
        if (inp.Ncomp() != 2 && inp.Ncomp() != 4) {
            std::cerr << "apply_Pauli: Input function must have 2 or 4 components." << std::endl;
            return;
        }
        if (out.Ncomp() != 2 && out.Ncomp() != 4) {
            std::cerr << "apply_Pauli: Output function must have 2 or 4 components." << std::endl;
            return;
        }
        switch (pauli)
        {
        case 0:
            // Apply Pauli-X matrix
            for (int i = 0; i < inp.Ncomp(); i = i + 2) {
                out.CompC[i] = inp.CompC[i+1]; //WARNING: No copy might create some issues down the line
                out.CompC[i+1] = inp.CompC[i];
            }
            break;
        case 1:
            // Apply Pauli-Y matrix
            for (int i = 0; i < inp.Ncomp(); i = i + 2) {
                out.CompC[i] = inp.CompC[i+1]; //WARNING: No copy might create some issues down the line
                out.CompC[i]->rescale(-1.0i);
                out.CompC[i+1] = inp.CompC[i];
                out.CompC[i+1]->rescale(1.0i);

                // rescale(1.0, out.CompC[i], -1.0i, inp.CompC[i+1]); //WARNING: No copy might create some issues down the line
                // out.CompC[i+1]->multiply(1.0i, inp.CompC[i]);
            }
            break;
        case 2:
            // Apply Pauli-Z matrix
            for (int i = 0; i < inp.Ncomp(); i = i + 2) {
                // out.CompC[i]->multiply(-1.0i, inp.CompC[i+1]); //WARNING: No copy might create some issues down the line
                out.CompC[i+1] = inp.CompC[i+1]; //WARNING: No copy might create some issues down the line
                out.CompC[i+1]->rescale(-1.0);
            }
            break;
        default:
            std::cerr << "Invalid Pauli matrix index: " << pauli << std::endl;
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
        inp.rescale(1.0 / norm);
    }
}