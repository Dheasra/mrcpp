#pragma once

#include "mpi_utils.h"
#include "trees/FunctionTreeVector.h"
#include "utils/CompFunction.h"

namespace mrcpp {
    //WARNING: These functions assume that the input functions are spinors, i.e. they have either 2 or 4 components and are complex valued.

    void apply_Pauli(CompFunction<3> &out, const CompFunction<3> &inp, int pauli, double prec = -1.0, bool conjugate = false);

    // void dot_spinor(CompFunction<3> &out, const CompFunction<3> &inp_a, const CompFunction<3> &inp_b, double prec = -1.0, bool conjugate = false);

    void normalize_spinor(CompFunction<3> &inp, double prec = -1.0);

    // std::array<std::array<ComplexDouble>> compute_overlap(const CompFunction<3> &bra, const CompFunction<3> &ket, double prec = -1.0, bool conjugate = false);

    // void compute_overlap_kramer_partner(ComplexDouble &overlap, const CompFunction<3> &bra, const CompFunction<3> &ket, double prec = -1.0, bool conjugate = false);
}