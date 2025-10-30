// enumerations
enum { _EVEN, _ODD };
enum { _NO_DEFAULT_SET, _DEFAULT_SET };
enum { _NO_REORDERING, _REORDER };
enum { _ADD, _COPY };
enum { _ORDINARY, _SCHWARZ, _ODDEVEN };
enum { _RES, _NO_RES };
enum { _STANDARD, _LIME }; // formats
enum { _READ, _WRITE };
enum { _NO_SHIFT };
enum { _BTWN_ORTH = 20 };
enum { _GLOBAL_FGMRES, _K_CYCLE, _COARSE_GMRES, _SMOOTHER }; // GMRES types
enum { _ADDITIVE = 1, _TWO_COLOR = 2, _SIXTEEN_COLOR = 3 }; // Schwarz types, correspond to g.method
enum { _FINE, _INTERMEDIATE, _COARSEST };                   // level/grid types
enum { _ZERO, _ONES, _KTH_UNIT, _RANDOM, _INDEX };          // vector initialization types
enum { _LEX_ORDER, _BLOCK_ORDER };                          // vector orderings
enum {
  _DENSE_QR,
  _GRAM_SCHMID,
  _MODIFIED_GRAM_SCHMIDT,
  _HOUSEHOLDER_QR,
  _GRAM_SCHMID_ON_AGGREGATES
}; // orthogonalization types TODO: GramSchmidt enum currently doubled with profiling enum -> no T
   // in Schmidt for now
enum { _T, _Z, _Y, _X };       // time-space dimensions
enum { _NEGATIVE, _POSITIVE }; // relative directions
enum { _COARSE_GLOBAL };
enum { _FULL_SYSTEM, _EVEN_SITES, _ODD_SITES };
enum { _LEFT, _RIGHT, _NOTHING };
enum { _PERIODIC, _ANTIPERIODIC, _TWISTED, _DIRICHLET };
enum {
  _GIP,
  _PIP,
  _LA2,
  _LA6,
  _LA8,
  _LA,
  _CPY,
  _SET,
  _PR,
  _SC,
  _NC,
  _SM,
  _OP_COMM,
  _OP_IDLE,
  _ALLR,
  _GD_COMM,
  _GD_IDLE,
  _GRAM_SCHMIDT,
  _GRAM_SCHMIDT_ON_AGGREGATES,
  _SM1,
  _SM2,
  _SM3,
  _SM4,
  _SMALL1,
  _SMALL2,
  _NUM_PROF
}; // _NUM_PROF has always to be the last constant
enum { _VTS = 20 };
enum {
  _TRCKD_VAL,
  _STP_TIME,
  _SLV_ITER,
  _SLV_TIME,
  _CRS_ITER,
  _CRS_TIME,
  _SLV_ERR,
  _CGNR_ERR,
  _NUM_OPTB
};
