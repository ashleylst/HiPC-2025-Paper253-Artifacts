#include <cstdlib>
#include <complex>


extern "C" void mm3x3t3xn_negb_1xtile(
        const std::complex<double>* a,
        const std::complex<double>* b,
        std::complex<double> *c,
        std::size_t N);
extern "C" void mm3x3t3xn_nega_1xtile(
        const std::complex<double>* a,
        const std::complex<double>* b,
        std::complex<double> *c,
        std::size_t N);
extern "C" void mm3x3t3xn_nega_8xtile(
        const std::complex<double>* a,
        const std::complex<double>* b,
        std::complex<double> *c,
        std::size_t N);
