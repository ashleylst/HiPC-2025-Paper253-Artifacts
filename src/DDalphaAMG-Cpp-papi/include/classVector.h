#ifndef VECTOR_HEADER
#define VECTOR_HEADER

class Vector {

public:
  Vector() = delete; // No point in constructing a vector without size
  Vector(const int &size);
  Vector(const int &size,
         const int &type); // type = ( _ZERO | _ONES | _KTH_UNIT | _RANDOM | _INDEX )
  Vector(const int &size, const int &type, const int &k); // k required if tpye == _KTH_UNIT

  ~Vector() = default;

  // NOTE: Currently WE are not managing resources (std::valarray is!), so no copy und move
  // semantics required and destructor = default

  Vector operator+(const Vector &vec2) const;
  Vector operator-(const Vector &vec2) const;
  Vector operator*(const int &scalar) const { return (*this) * (std::complex<double>)(scalar); }
  Vector operator*(const double &scalar) const { return (*this) * (std::complex<double>)(scalar); }
  Vector operator*(const std::complex<double> &scalar) const;
  std::complex<double> operator*(const Vector &vec2) const; // inner product

  void normalize();
  Vector changeOrdering(const int newOrdering, int *translationTable);
  double getNorm() const;
  void set(const int &type);
  void set(const int &type, const int &position);
  void print() const;
  void print(const int &start, const int &end) const;

  int globalLength;
  int localLength;
  int isZero = 0;
  int ordering = _LEX_ORDER;

  std::valarray<std::complex<double>> data;
};

// Operators to enable left multiplication with a scalar
Vector operator*(const int &scalar, const Vector &vec);
Vector operator*(const double &scalar, const Vector &vec);
Vector operator*(const std::complex<double> &scalar, const Vector &vec);

#endif
