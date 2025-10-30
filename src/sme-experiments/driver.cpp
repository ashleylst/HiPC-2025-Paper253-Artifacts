#include "mm3x3t3xN_sme.h"
#include "mm3x3t3xN.hpp"
#include "sme_helpers.hpp"

#include <algorithm>
#include <array>
#include <cassert>
#include <complex>
#include <cstdint>
#include <functional>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <ranges>
#include <sstream>
#include <string>
#include <vector>

#include <cblas.h>

void complex_iota(std::vector<std::complex<double>>& matrix)
{
    using std::complex_literals::operator""i;
    std::ranges::iota(matrix, 0);
    std::ranges::copy(
            std::views::transform(matrix,
            [](std::complex<double> a)
                {return (a*2.0)*0.3 + (a*2.0+1.0)*0.3i;}),
            matrix.begin());
}

void print_cmat(
        std::vector<std::complex<double>>& matrix,
        std::size_t ldA,
        std::size_t padding = 0,
        bool col_major = false)
{
    std::size_t M = matrix.size()/ldA;
    std::size_t N = ldA;

    std::string bleftmargin(padding, ' ');

    for (std::size_t i = 0; i < M; i++)
    {
        std::cout << bleftmargin;
        for(std::size_t j = 0; j < N; j++)
        {
            std::stringstream numberstr;
            if(col_major)
            {
                numberstr << matrix[(j*M+i)].real() << "+" << matrix[(j*M+i)].imag() << "i";
            }
            else
            {
                numberstr << matrix[(i*N+j)].real() << "+" << matrix[(i*N+j)].imag() << "i";
            }
            std::cout << "  " << std::setw(10) << numberstr.str();
        }
        std::cout << "\n";
    }

}


int main(int argc, char *argv[])
{
    if ( (argc != 3) && (argc != 4) )
    {
        std::cout << argv[0] << " N <cpp|blas|sme_nega1|sme_nega8|sme_negb> [roi_iterations]\n";
        return -1;
    }
    std::uint64_t N = std::stoul(argv[1]);
    std::string method{argv[2]};
    std::uint64_t roi_iterations = 1;
    if (argc == 4)
    {
        roi_iterations = std::stoul(argv[3]);
    }

    std::size_t vlen = get_sme_vlen();

    std::vector<std::complex<double>> A(3*3);
    std::vector<std::complex<double>> B(3*N);
    std::vector<std::complex<double>> C(3*N);

    complex_iota(A);
    complex_iota(B);

    std::string bleftmargin(35,' ');
    std::cout << bleftmargin << "B:\n";

    print_cmat(B, N, 35);

    std::cout << "A:\n";

    print_cmat(A, 3, 0, true);

    
    auto mm3x3t3xn_blas = [](
            std::complex<double>* a,
            std::complex<double>* b,
            std::complex<double>* c,
            std::size_t N)
    {
        std::complex<double> one{1.0,0.0};
        cblas_zgemm(CBLAS_ORDER::CblasRowMajor,
                    CBLAS_TRANSPOSE::CblasTrans, CBLAS_TRANSPOSE::CblasNoTrans,
                    3, N, 3, &one,
                    a, 3,
                    b, N, &one,
                    c, N);
    };

    std::ranges::fill(C, 0.0);

    std::function<void()> implementation;
    if ("cpp" == method)
    {
        implementation=std::bind(mm3x3t3xN_cpp, A.data(), B.data(), C.data(), N);
    }
    else if ("blas" == method)
    {
        implementation=std::bind(mm3x3t3xn_blas, A.data(), B.data(), C.data(), N);
    }
    else if ("sme_negb" == method)
    {
        implementation=std::bind(mm3x3t3xn_negb_1xtile, A.data(), B.data(), C.data(), N);
    }
    else if ("sme_nega1" == method)
    {
        implementation=std::bind(mm3x3t3xn_nega_1xtile, A.data(), B.data(), C.data(), N);
    }
    else if ("sme_nega8" == method)
    {
        implementation=std::bind(mm3x3t3xn_nega_8xtile, A.data(), B.data(), C.data(), N);
    }
    else
    {
        std::string msg("Unknown method: ");
        msg += method;
        throw std::runtime_error(msg);
    }

    for (std::size_t i = 0; i < roi_iterations; i++)
    {
        implementation();
    }

    std::cout << bleftmargin << "C:\n";

    print_cmat(C, N, 35);

    return 0;
}
