/*
 * MRCPP, a numerical library based on multiresolution analysis and
 * the multiwavelet basis which provide low-scaling algorithms as well as
 * rigorous error control in numerical computations.
 * Copyright (C) 2021 Stig Rune Jensen, Jonas Juselius, Luca Frediani and contributors.
 *
 * This file is part of MRCPP.
 *
 * MRCPP is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * MRCPP is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with MRCPP.  If not, see <https://www.gnu.org/licenses/>.
 *
 * For information on the complete list of contributors to MRCPP, see:
 * <https://mrcpp.readthedocs.io/>
 */

/*
 *
 *
 *  \date Jul 18, 2009
 *  \author Jonas Juselius <jonas.juselius@uit.no> \n
 *          CTCC, University of Tromsø
 *
 * \breif
 */

#include "SchrodingerEvolution_CrossCorrelation.h"

#include <fstream>

#include "MRCPP/config.h"
#include "MRCPP/constants.h"

#include "utils/Printer.h"
#include "utils/details.h"

using namespace Eigen;

namespace mrcpp {

SchrodingerEvolution_CrossCorrelation::SchrodingerEvolution_CrossCorrelation(int amount, int k, int t)
    : type(t), order(k), amount(amount)
{
    if (this->order < 1 or this->order > MaxOrder) MSG_ABORT("Invalid cross correlation order: " << this->order);
    switch (this->type) {
        case (Interpol):
            MSG_ERROR("Not implemented yet filter type: " << this->type);
        case (Legendre):
            break;
        default:
            MSG_ERROR("Unknown filter type: " << this->type);
    }

    setCCCPath(details::find_filters());

    readCCCBin();
}

/*
SchrodingerEvolution_CrossCorrelation::SchrodingerEvolution_CrossCorrelation(int t, const MatrixXd &L, const MatrixXd &R)
    : type(t), order(L.cols() / 2 - 1)
{
    if (this->order < 1 or this->order > MaxOrder) MSG_ABORT("Invalid cross correlation order, " << this->order);
    if (R.cols() != L.cols()) MSG_ABORT("Right and Left cross correlation have different order!");
    switch (this->type) {
        case (Interpol):
        case (Legendre):
            break;
        default:
            MSG_ERROR("Unknown filter type: " << this->type);
    }

    this->Left = L;
    this->Right = R;
}
*/

void SchrodingerEvolution_CrossCorrelation::setCCCPath(const std::string &lib) {
    switch (this->type) {
        case (Interpol):
            MSG_ERROR("Not implemented yet filter type: " << this->type);
            break;
        case (Legendre):
            this->path = lib + "/Schrodinger_evolution_cross_correlation_coefficients_Legendre_scaling_type.bin";
            break;
        default:
            MSG_ERROR("Invalid CrossCorrelation type");
    }
}

void SchrodingerEvolution_CrossCorrelation::readCCCBin()
{
    std::ifstream input_file(this->path.c_str(), std::ios::binary);

    if (not input_file) MSG_ABORT("Could not open cross correlation: " << this->path);

    // Read the text length
    int text_length;
    input_file.read(reinterpret_cast<char*>(&text_length), sizeof(text_length));

    // Read the Unicode characters
    std::vector<char32_t> unicode_chars(text_length);
    input_file.read(reinterpret_cast<char*>(unicode_chars.data()), sizeof(char32_t) * text_length);

    // Read the amount of matrices
    int K;
    input_file.read(reinterpret_cast<char*>(&K), sizeof(K));

    // Read the size/order of each matrix
    int order;
    input_file.read(reinterpret_cast<char*>(&order), sizeof(order));

    // Read the matrices
    std::vector<Eigen::MatrixXd> C_even(K, Eigen::MatrixXd(order, order));
    auto data_amount = order * order * sizeof(double);
    for (auto& matrix : C_even) input_file.read(reinterpret_cast<char*>(matrix.data()), data_amount);

    // Print the text length
    std::cout << text_length << std::endl;
    
    // Print the text
    for (char32_t c : unicode_chars) {
        std::wcout << static_cast<wchar_t>(c);
    }

    // Print the matrices
    std::cout << std::endl;
    std::cout << "----------------------------------" << std::endl;
    for (auto& matrix : C_even)
    {
        std::cout << matrix  << std::endl;
        std::cout << "----------------------------------" << std::endl;
    }

    //TODO: create matrix containing the appropriate amount of coefficients

    for (int k = 0; k < this->amount; k++)
    {
        Matrix.push_back( C_even[k].block(0, 0, this->order + 1, this->order + 1).eval() );
    }

/*
    int K = this->order + 1;
    this->Left = MatrixXd::Zero(K * K, 2 * K);
    this->Right = MatrixXd::Zero(K * K, 2 * K);
    double dL[2 * K];
    double dR[2 * K];
    for (int i = 0; i < K * K; i++) {
        L_fis.read((char *)dL, sizeof(double) * 2 * K);
        R_fis.read((char *)dR, sizeof(double) * 2 * K);
        for (int j = 0; j < 2 * K; j++) {
            if (std::abs(dL[j]) < MachinePrec) dL[j] = 0.0;
            if (std::abs(dR[j]) < MachinePrec) dR[j] = 0.0;
            this->Left(i, j) = dL[j];
            this->Right(i, j) = dR[j];
        }
    }
*/
    C_even.clear();
    input_file.close();
}

} // namespace mrcpp
