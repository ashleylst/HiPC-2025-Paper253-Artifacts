#include <algorithm>
#include <array>
#include <cassert>
#include <complex>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <sstream>
#include <vector>

#include <cblas.h>

void mvmh_block_test(std::complex<double> *eta, std::complex<double> *D, std::complex<double> *phi,
                const int blen){
  //memset(eta, 0, 3 * blen * sizeof(std::complex<T>));
  for(int col = 0; col < 3; col++){
    for(int row = 0; row < 3; row++){
      for(int j = 0; j < blen; j++){
        eta[blen*row + j] += conj(D[col*3 + row]) * phi[blen*col + j];
      }
    }
  }
}

void mvm_block_test(std::complex<double> *eta, std::complex<double> *D, std::complex<double> *phi,
                     const int blen){
  //memset(eta, 0, 3 * blen * sizeof(std::complex<T>));
    for(int col = 0; col < 3; col++){
      for(int row = 0; row < 3; row++){
        for(int j = 0; j < blen; j++){
          eta[blen*col + j] += D[col*3 + row] * phi[blen*row + j];
        }
      }
    }
}



int main(int argc, char *argv[])
{
    if (argc != 2)
    {
        std::cout << argv[0] << " N\n";
        return -1;
    }
    std::uint64_t N = std::stoul(argv[1]);

//    std::size_t vlen = get_sme_vlen();

    std::vector<double> A(3*3*2);
    std::vector<double> B(3*N*2);
    std::vector<double> C(3*N*2);

    std::iota(A.begin(), A.end(), 0.0);

    std::iota(B.begin(), B.end(), 0.0);
    std::transform(B.begin(), B.end(), B.begin(), [](double a){return a*0.3;});

    std::string bleftmargin = "                                   ";
    std::cout << bleftmargin << "B:\n";

    for (std::size_t i = 0; i < 3; i++)
    {
        std::cout << bleftmargin;
        for(std::size_t j = 0; j < N; j++)
        {
            std::stringstream numberstr;
            numberstr << B[(i*N+j)*2] << "+" << B[(i*N+j)*2+1] << "i";
            std::cout << "  " << std::setw(10) << numberstr.str();
        }
        std::cout << "\n";
    }

    std::cout << "A:\n";

    for (std::size_t i = 0; i < 3; i++)
    {
        for(std::size_t j = 0; j < 3; j++)
        {
            std::stringstream numberstr;
            numberstr << A[(j*3+i)*2] << "+" << A[(j*3+i)*2+1] << "i";
            std::cout << "  " << std::setw(10) << numberstr.str();
        }
        std::cout << "\n";
    }

    std::complex<double> one{1.0,0.0};
    std::fill(C.begin(), C.end(), 0.0);
    cblas_zgemm(CBLAS_ORDER::CblasRowMajor,
                CBLAS_TRANSPOSE::CblasNoTrans, CBLAS_TRANSPOSE::CblasNoTrans,
                3, N, 3, &one,
                A.data(), 3,
                B.data(), N, &one,
                C.data(), N);

    std::cout << bleftmargin << "C (reference):\n";
    for (std::size_t i = 0; i < 3; i++)
    {
        std::cout << bleftmargin;
        for(std::size_t j = 0; j < N; j++)
        {
            std::stringstream numberstr;
            numberstr << C[(i*N+j)*2] << "+" << C[(i*N+j)*2+1] << "i";
            std::cout << "  " << std::setw(10) << numberstr.str();
        }
        std::cout << "\n";
    }

    std::fill(C.begin(), C.end(), 0.0);

    mvm_block_test((std::complex<double>*)C.data(), (std::complex<double>*)A.data(), (std::complex<double>*)B.data(), N);

    std::cout << bleftmargin << "C (ours):\n";

    for (std::size_t i = 0; i < 3; i++)
    {
        std::cout << bleftmargin;
        for(std::size_t j = 0; j < N; j++)
        {
            std::stringstream numberstr;
            numberstr << C[(i*N+j)*2] << "+" << C[(i*N+j)*2+1] << "i";
            std::cout << "  " << std::setw(10) << numberstr.str();
        }
        std::cout << "\n";
    }

    return 0;
}

