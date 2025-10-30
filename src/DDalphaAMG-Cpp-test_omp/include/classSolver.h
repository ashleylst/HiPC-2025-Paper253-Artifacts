#ifndef SOLVER_HEADER
#define SOLVER_HEADER

// TODO: compute_core_start_end can probably be moved to constructors

/**
 *  Solver class.
 */
class Solver {

public:
  Solver(); // default constructor
  Solver(operator_struct<double> *, int, int, int, double, int, int, int,
         void (*precond)(vector_double, vector_double, vector_double, int, Level *),
         void (*evalOp)(vector_double, vector_double, operator_struct<double> *, Level *), Level *);

  Solver(const Solver &) = delete;            // copy constructor
  Solver &operator=(const Solver &) = delete; // assignment operator

  Solver(const Solver &&) = delete;            // move constructor
  Solver &operator=(const Solver &&) = delete; // move assignment operator

  ~Solver(); // default destructor

  virtual void allocateMemory();

  void (*preconditioner)(vector_double, vector_double, vector_double, const int, Level *);
  void (*applyOperator)(vector_double, vector_double, operator_struct<double> *, Level *); // TODO: This should rather be an operator* overload (maybe?)
  //vector_double solve();

  Level *level;

  vector_double x, b, r, w; // solution, residual and one buffer vector
  vector_double *V;    /**< non-preconditioned subspace basis vector*/
  vector_double *Z;    /**< preconditioned subspace basis vector*/

  //complex_double **H;  /**< projected subspace matrix (typically something like V* A V) */

  config_double *D, *clover;
  operator_struct<double> *matrix;

  double tolerance;
  double norm;

  int problemSize;
  int restartNumber;
  int restartLength;
  int initialGuessIsZero = 1;
  int timing;
  int verbose;
  int vectorLength;
  int vectorStart = 0;
  int vectorEnd = problemSize;
};


/**
 *  GMRES class, inherit from Solver class.
 */
class Gmres : public Solver {

public:
  Gmres(); // default constructor
  Gmres(operator_struct<double> *, int, int, int, double, const int, const int, const int, const int,
        int, void (*precond)(vector_double, vector_double, vector_double, const int, Level *),
        void (*evalOp)(vector_double, vector_double, operator_struct<double> *, Level *), Level *);

  Gmres(const Gmres &) = delete;            // copy constructor
  Gmres &operator=(const Gmres &) = delete; // copy assignment operator

  Gmres(const Gmres &&) = delete;            // move constructor
  Gmres &operator=(const Gmres &&) = delete; // move assignment operator

  ~Gmres(); // default destructor

  void allocateMemory() override;
  void initialize(operator_struct<double> *op, int resLen, int resNum, int probSize, double tol,
                  int t, int precKind, int tim, int print, int vecLen,
                  void (*precond)(vector_double, vector_double, vector_double, int, Level *),
                  void (*evalOp)(vector_double, vector_double, operator_struct<double> *, Level *),
                  Level *lev);

  int arnoldiStep(complex_double *Hp);
  void qrUpdate(complex_double *Hp);
  void computeSolution(int, complex_double *Hp);
  vector_double solve();

  int j = -1;
  complex_double *sinus;
  complex_double *cosinus;
  complex_double *gamma;
  complex_double *y;

private:
  complex_double beta;
  int type;
  int preconditionerKind;
  int finished = 0, iter = 0, res;
  double normR0 = 1, gammaJPlusOne = 1, elapsedTime = 0;
};


/**
 *  Schwarz class, inherit from Solver class.
 */
class Schwarz : public Solver {

public:
  Schwarz() = default;

  Schwarz(const Schwarz &);            // copy constructor
  Schwarz &operator=(const Schwarz &); // copy assignment operator

  Schwarz(const Schwarz &&) = delete;            // move constructor
  Schwarz &operator=(const Schwarz &&) = delete; // move assignment operator

  ~Schwarz();

  void allocateMemory() override;
  void initialize(void (*evalOp)(vector_double, vector_double, operator_struct<double> *, Level *),
                  Level *); /**< initialize the Schwarz class */

  inline int connectLink(int, int, int, int, int, int, int *, int *);
  void defineLayout();
  void boundaryUpdate();
  void setupOddEvenOperator();
  void setupFineOperator(operator_struct<double> *);
  void setupIntermediateOperator();
  void setupCoarseOperator();

  void solve();

  complex_double **odd_even_buffer;
  complex_double **local_minres_buffer;

  int type;
  int number_of_colors;
  int number_of_cycles;
  int number_of_blocks;
  int has_residual;
  int *index_table[4];
  int *odd_even_index_table[4];
  int direction_length[4];
  int direction_length_even[4];
  int direction_length_odd[4];
  int block_odd_even_offset;
  int number_of_even_block_sites;
  int number_of_odd_block_sites;
  int number_of_aggregates;
  int block_vector_size;
  int number_of_block_sites;
  int block_boundary_length[9];
  int **block_list;
  int *block_list_length;

  block_struct *block;
  operator_struct<double> op; // TODO: not a pointer?

private:
  // members for general schwarz
  int k;

  // function pointers
  void (*blockOperator)(vector_double eta, vector_double phi, int start, Schwarz *s, Level *level);
  void (Schwarz::*boundaryOperator)(vector_double eta, vector_double phi, int dir, int k);
  void (Schwarz::*blockSolve)(vector_double phi, vector_double eta, vector_double latest_iter,
                              int start);

  // bounday operators
  void blockBoundaryOperator(vector_double eta, vector_double phi, int dir, int k);
  void coarseBlockBoundaryOperator(vector_double eta, vector_double phi, int dir, int k);

  // block solvers
  void blockSolveOddeven(vector_double, vector_double, vector_double, int);
  void localMinres(vector_double, vector_double, vector_double, int);

  // members for additive schwarz
  void additiveSolve();

  int nb_thread_start, nb_thread_end;
  vector_double Dphi, latestIter, latestIter2;

  // members for two color schwarz
  void twoColorSolve();

  int initialResidual;
  int residualCommunication;
  void (*communicate[2])(vector_double, int, int, comm_struct<double> *, Level *);
  int blockThreadStart[8];
  int blockThreadEnd[8];
  int start;
  int end;

  // members for sixteen color schwarz
  void sixteenColorSolve();
};
#endif
