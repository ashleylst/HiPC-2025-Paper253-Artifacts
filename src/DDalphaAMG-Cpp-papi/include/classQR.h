// #ifndef QR_HEADER
// #define QR_HEADER
// extern "C" {
// extern void zgeqrf_(int *, int *, complex_double *, int *, complex_double *, complex_double *,
//                     int *, int *); // computes HH reflectors and R for a QR decomposition
// extern void zungqr_(int *, int *, int *, complex_double *, int *, complex_double *,
//                     complex_double *, int *, int *); // recover Q from HH reflectors
// }
// class Orthogonalization {
//
// public:
//   Orthogonalization();
//   Orthogonalization(int, int, int, int, Level *, Thread *);
//
//   Orthogonalization(const Orthogonalization &) = delete;            // copy constructor
//   Orthogonalization &operator=(const Orthogonalization &) = delete; // assignment operator
//
//   Orthogonalization(const Orthogonalization &&) = delete;            // move constructor
//   Orthogonalization &operator=(const Orthogonalization &&) = delete; // move assignment operator
//
//   ~Orthogonalization(); // default destructor
//
//   //   void allocateMemory();
//   void orthogonalize(complex_double **);
//
//   Level *level;
//   Thread *threading;
//
//   complex_double *coefficients;
//   complex_double *tmp;
//   complex_double *work;
//
//   double norm;
//
//   int type;
//   int vectorSize;
//   int numberOfVectors;
//   int beginOrthogonalization;
//   int startThreading;
//   int endThreading;
//   int lwork;
//   int info;
//
//
// private:
//   void allocateMemory();
//   void gramSchmidt(complex_double **);
//   void modifiedGramSchmidt(complex_double **);
//   void householderQR(complex_double **);
//   void denseQR(complex_double **);
//   void gramSchmidtOnAggregates(complex_double **V);
// };
// #endif
