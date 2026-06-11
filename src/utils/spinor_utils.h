#pragma once

// #include "mpi_utils.h"
#include "trees/FunctionTreeVector.h"
#include "utils/CompFunction.h"
// using namespace std::complex_literals;

namespace mrcpp {
    //NOTE: These functions assume that the input functions are spinors, i.e. they have either 2 or 4 components and are complex valued.

    template<int D> void apply_gamma(CompFunction<D> &inp, int index);

    // orginial implementation of apply_alpha, may be more efficient because everything is done by reference, but it doesn't work in mrchem
    void apply_Pauli(CompFunction<3> &out, CompFunction<3> &inp, int pauli, double prec = -1.0, bool conjugate = false);

    // void dot_spinor(CompFunction<3> &out, const CompFunction<3> &inp_a, const CompFunction<3> &inp_b, double prec = -1.0, bool conjugate = false);

    /* @brief Normalization of spinor functions. 
     * The function is normalized in place, i.e. the input function is modified. 
     * The norm is computed as the sum of the norms of the components, i.e. ||Psi||^2 = ||Psi_1||^2 + ||Psi_2||^2 for a 2-component spinor. 
     * For a 4-component spinor, the norm is computed as ||Psi||^2 = ||Psi_1||^2 + ||Psi_2||^2 + ||Psi_3||^2 + ||Psi_4||^2. 
     * The function is normalized by rescaling each component by the same factor, i.e. Psi_i -> Psi_i / ||Psi||. 
     * @param inp The input spinor function to be normalized. It is modified in place.
     * @param prec The precision for the normalization. If the norm is smaller than prec, the norm is considered to be zero.
     * NOTE: This function rescales the components of the input function rather than its overall scaling coefficient.
     *       This makes it slower but allows to reset the scaling coefficient to a reasonable value.
     */
    void normalize_spinor(CompFunction<3> &inp, double prec = -1.0);

    // std::array<std::array<ComplexDouble>> compute_overlap(const CompFunction<3> &bra, const CompFunction<3> &ket, double prec = -1.0, bool conjugate = false);

    // void compute_overlap_kramer_partner(ComplexDouble &overlap, const CompFunction<3> &bra, const CompFunction<3> &ket, double prec = -1.0, bool conjugate = false);
    // CompFunctionVector disjoin(CompFunctionVector &Phi, int spin);
    // CompFunctionVector adjoin(CompFunctionVector &Phi_a, CompFunctionVector &Phi_b);
}