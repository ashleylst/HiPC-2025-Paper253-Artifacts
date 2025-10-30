
#include "main.h"

inline void mapGlobalIndexToLocal(int *process, int *localIndex, int globalIndex,
                                  int localVectorSize) {
  // NOTE: Currently only works for lexicographical ordering
  *process = globalIndex / localVectorSize;
  *localIndex = globalIndex % localVectorSize;
}

inline void mapLocalIndexToGlobal(int *globalIndex, int process, int localIndex,
                                  int localVectorSize) {
  // NOTE: Currently only works for lexicographical ordering
  *globalIndex = process * localVectorSize + localIndex;
}

void Vector::print() const { print(0, globalLength - 1); }

void Vector::print(const int &start, const int &end) const {

  int firstProcess, lastProcess;
  int startOnFirstProcess, endOnLastProcess;
  int globalIndex;
  int printStart = 0, printEnd = 0; // Setting to zero ensures that idle processors dont print stuff

  mapGlobalIndexToLocal(&firstProcess, &startOnFirstProcess, start, localLength);
  mapGlobalIndexToLocal(&lastProcess, &endOnLastProcess, end, localLength);

  if (my_rank == firstProcess) {
    printStart = startOnFirstProcess;
    printEnd = localLength - 1;
  }
  if ((my_rank > firstProcess) && (my_rank < lastProcess)) {
    printStart = 0;
    printEnd = localLength - 1;
  }
  if (my_rank == lastProcess) {
    printStart = 0;
    printEnd = endOnLastProcess;
  }
  for (int i = printStart; i <= printEnd; i++) {
    mapLocalIndexToGlobal(&globalIndex, my_rank, i, localLength);
    printf("global index: %d, on rank: %d =  %lf%+lfi\n", globalIndex, my_rank, real(data[i]),
           imag(data[i]));
  }
}

Vector::Vector(const int &size) : Vector(size, _ZERO) {}

Vector::Vector(const int &size, const int &type) : Vector(size, type, 0) {}

Vector::Vector(const int &size, const int &type, const int &k) {

  int cores;
  MPI_Comm_size(MPI_COMM_WORLD, &cores);

  globalLength = size;
  localLength = size / cores;
  data.resize(localLength);

  Vector::set(type, k);
}

Vector Vector::operator+(const Vector &vec2) const {

  ASSERT(globalLength == vec2.globalLength);

  // quick return
  if (vec2.isZero)
    return *this;
  if (isZero)
    return vec2;

  Vector tmp(globalLength);
  tmp.data = data + vec2.data;
  tmp.isZero = 0;

  return tmp;
}

Vector Vector::operator-(const Vector &vec2) const {

  ASSERT(globalLength == vec2.globalLength);

  // quick return
  if (vec2.isZero)
    return *this;

  Vector tmp(globalLength);
  tmp.data = data - vec2.data;
  tmp.isZero = 0;

  return tmp;
}

Vector Vector::operator*(const std::complex<double> &scalar) const {

  // quick return
  if (isZero || scalar == 0.0)
    return Vector(globalLength);

  Vector tmp(globalLength);
  tmp.data = scalar * data;
  tmp.isZero = 0;

  return tmp;
}

// Operators to enable left multiplication with a scalar
Vector operator*(const int &scalar, const Vector &vec) {
  return vec * (std::complex<double>)(scalar);
}
Vector operator*(const double &scalar, const Vector &vec) {
  return vec * (std::complex<double>)(scalar);
}
Vector operator*(const std::complex<double> &scalar, const Vector &vec) { return vec * scalar; }

std::complex<double>
Vector::operator*(const Vector &vec2) const { // NOTE: This implements the inner product

  ASSERT(globalLength == vec2.globalLength);

  // quick returns
  if (isZero || vec2.isZero)
    return 0;

  std::complex<double> localInnerProduct = 0;
  std::complex<double> globalInnerProduct = 0;

  Vector localIPValues(globalLength);
  Vector thisVectorConj = *this;

  for (int i = 0; i < localLength; i++)
    thisVectorConj.data[i] = conj(data[i]);

  localIPValues.data = thisVectorConj.data * vec2.data;
  localInnerProduct = localIPValues.data.sum();
  MPI_Allreduce(&localInnerProduct, &globalInnerProduct, 1, MPI_COMPLEX_double, MPI_SUM,
                MPI_COMM_WORLD);

  return globalInnerProduct;
}

void Vector::set(const int &type) {

  if (type == _KTH_UNIT)
    error0("Incomplete call to Vector::set()! type == _KTH_UNIT requires second argument!\n");

  Vector::set(type, 0);
}

void Vector::set(const int &type, const int &k) {

  switch (type) {
  case _ZERO:
    isZero = 1;
    break; // NOTE: std::valarray comes already zero-initialized
  case _ONES:
    data = std::complex<double>(1.0, 0.0);
    break;
  case _KTH_UNIT: {
    int targetRank;
    int localIndex;
    mapGlobalIndexToLocal(&targetRank, &localIndex, k - 1, localLength);
    if (my_rank == targetRank)
      data[localIndex] = std::complex<double>(1.0, 0.0);
    break;
  }
  case _RANDOM: {
    std::random_device seed;
    std::mt19937 randomInt(seed());
    std::uniform_real_distribution<> unit(-1.0, 1.0);
    for (int i = 0; i < localLength; i++)
      data[i] = std::complex<double>(unit(randomInt), unit(randomInt));
    break;
  }
  case _INDEX: {
    int targetRank;
    int localIndex;
    for (int i = 0; i < globalLength; i++) {
      mapGlobalIndexToLocal(&targetRank, &localIndex, i, localLength);
      if (my_rank == targetRank)
        data[localIndex] = std::complex<double>(i, i);
    }
    break;
  }
  default:
    error0("Unknown vector type!\n");
    break;
  }
}

double Vector::getNorm() const {

  if (isZero)
    return 0;

  double localNorm = 0;
  double globalNorm = 0;

  for (int i = 0; i < localLength; i++)
    localNorm += norm(data[i]);

  MPI_Allreduce(&localNorm, &globalNorm, 1, MPI_double, MPI_SUM, MPI_COMM_WORLD);

  return sqrt(globalNorm);
}

void Vector::normalize() {

  if (isZero)
    error0("Trying to normalize a zero vector!\n");

  data /= getNorm();
}

Vector Vector::changeOrdering(const int newOrdering, int *translationTable) {

  if (ordering == newOrdering || isZero) {
    warning0("Vector already in requested ordering or is zero.\n");
    return *this;
  }

  Vector tmp(globalLength);

  if (newOrdering == _BLOCK_ORDER) {
    for (int i = 0; i < localLength / 12; i++)
      for (int k = 0; k < 12; k++)
        tmp.data[12 * translationTable[i] + k] = data[12 * i + k];
    tmp.ordering = _BLOCK_ORDER;
  }

  if (newOrdering == _LEX_ORDER) {
    for (int i = 0; i < localLength / 12; i++)
      for (int k = 0; k < 12; k++)
        tmp.data[12 * i + k] = data[12 * translationTable[i] + k];
    tmp.ordering = _LEX_ORDER;
  }

  return tmp;
}
