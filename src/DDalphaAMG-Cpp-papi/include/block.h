//
// Created by shiting on 2024-11-07.
//

#ifndef DDALPHAAMG_BLOCK_H
#define DDALPHAAMG_BLOCK_H

template<typename T>
class BlockVector{
public:
  explicit BlockVector(std::complex<T> *input){
    this->v = input;
  };

  std::complex<T> *v;
};

/**
 * 12 of the first v | 12 of the second v | ...
 * @tparam T precision
 */
template<typename T>
class BlockVecPerSite : public BlockVector<T>{
public:
  BlockVecPerSite(std::complex<T> *in, int nv): BlockVector<T>(in){
    this->nv = nv;
  }
  void transform(int len, int blen);
  void detransform(int len, int blen);

  int nv;
};

/**
 * 1 of the first v | 1 of the second | ...
 * @tparam T precision
 */
template<typename T>
class BlockVecPerEntry : public BlockVector<T>{
public:
  BlockVecPerEntry(std::complex<T> *in): BlockVector<T>(in){};
  void transform(int len, int blen);
  void detransform(int len, int blen);
};

template <typename T> void BlockVecPerEntry<T>::detransform(int len, int blen) {
  auto *tmp = new std::complex<T>[blen * len];

  for (int i = 0; i < len; i++) {
    for (int j = 0; j < blen; j++) {
        tmp[j*len + i] = this->v[i*blen + j];
    }
  }

  for (int i = 0; i < len * blen; i++) {
    this->v[i] = tmp[i];
  }

  delete[] tmp;
}

template <typename T> void BlockVecPerEntry<T>::transform(int len, const int blen) {
  auto *tmp = new std::complex<T>[blen * len];

  for (int i = 0; i < len; i++) {
    for (int j = 0; j < blen; j++) {
        tmp[i*blen + j] = this->v[j*len + i];
    }
  }

  for (int i = 0; i < len * blen; i++) {
      this->v[i] = tmp[i];
  }

  //todo: change here for optimal
  delete[] tmp;
}

template <typename T> void BlockVecPerSite<T>::detransform(int len, const int blen) {
  auto *tmp = new std::complex<T>[blen*len];

  for (int i = 0; i < len/nv; i++) {
    for(int j = 0; j < blen; j++){
      for (int k = 0; k < nv; k++) {
        tmp[j*len + i*nv + k] = this->v[i*blen*nv + j*nv + k];
      }
    }
  }

  for(int i = 0; i < blen * len; i++){
    this->v[i] = tmp[i];
  }

  delete[] tmp;
}

template <typename T>
void BlockVecPerSite<T>::transform(int len, const int blen) {
  auto *tmp = new std::complex<T>[blen * len];

  for(int i = 0; i < len/nv; i++){
    for(int j = 0; j < blen; j++){
      for (int k = 0; k < nv; ++k) {
        tmp[i*blen*nv + j*nv + k] = this->v[j*len + i*nv + k];
      }
    }
  }

  for(int i = 0; i < blen * len; i++){
    this->v[i] = tmp[i];
  }

  delete[] tmp;
}
#endif // DDALPHAAMG_BLOCK_H
