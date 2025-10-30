#include <complex>

void mm3x3t3xN_cpp(
        const std::complex<double> *A,
        const std::complex<double> *B,
        std::complex<double> *C,
        const std::size_t N)
{
    //memset(eta, 0, 3 * N * sizeof(std::complex<T>));
    for(int col = 0; col < 3; col++)
    {
        for(int row = 0; row < 3; row++)
        {
            for(int j = 0; j < N; j++)
            {
                C[N*col + j] += A[row*3 + col] * B[N*row + j];
            }
        }
    }
}
