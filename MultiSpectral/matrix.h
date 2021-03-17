#ifndef _MATRIX_H_
#define _MATRIX_H_

#include "Coregistration.h"

#define endll endl << endl // double end line definition

typedef double FLOAT;      // double precision
//typedef float  FLOAT;    // single precision

class LocalMatrix {

public:

  // constructor / deconstructor
	LocalMatrix();                                                  // init empty 0x0 matrix
	LocalMatrix(const int32_t m,const int32_t n);                   // init empty mxn matrix
	LocalMatrix(const int32_t m,const int32_t n,const FLOAT* val_); // init mxn matrix with values from array 'val'
	LocalMatrix(const LocalMatrix &M);                                   // creates deepcopy of M
    ~LocalMatrix();

  // assignment operator, copies contents of M
	LocalMatrix& operator= (const LocalMatrix &M);

  // copies submatrix of M into array 'val', default values copy whole row/column/matrix
  void getData(FLOAT* val_,int32_t i1=0,int32_t j1=0,int32_t i2=-1,int32_t j2=-1);

  // set or get submatrices of current matrix
  LocalMatrix getMat(int32_t i1,int32_t j1,int32_t i2=-1,int32_t j2=-1);
  void   setMat(const LocalMatrix &M,const int32_t i,const int32_t j);

  // set sub-matrix to scalar (default 0), -1 as end replaces whole row/column/matrix
  void setVal(FLOAT s,int32_t i1=0,int32_t j1=0,int32_t i2=-1,int32_t j2=-1);

  // set (part of) diagonal to scalar, -1 as end replaces whole diagonal
  void setDiag(FLOAT s,int32_t i1=0,int32_t i2=-1);

  // clear matrix
  void zero();
  
  // extract columns with given index
  LocalMatrix extractCols (std::vector<int> idx);

  // create identity matrix
  static LocalMatrix eye (const int32_t m);
  void          eye ();

  // create diagonal matrix with nx1 or 1xn matrix M as elements
  static LocalMatrix diag(const LocalMatrix &M);
  
  // returns the m-by-n matrix whose elements are taken column-wise from M
  static LocalMatrix reshape(const LocalMatrix &M,int32_t m,int32_t n);

  // create 3x3 rotation matrices (convention: http://en.wikipedia.org/wiki/Rotation_matrix)
  static LocalMatrix rotMatX(const FLOAT &angle);
  static LocalMatrix rotMatY(const FLOAT &angle);
  static LocalMatrix rotMatZ(const FLOAT &angle);

  // simple arithmetic operations
  LocalMatrix  operator+ (const LocalMatrix &M); // add matrix
  LocalMatrix  operator- (const LocalMatrix &M); // subtract matrix
  LocalMatrix  operator* (const LocalMatrix &M); // multiply with matrix
  LocalMatrix  operator* (const FLOAT &s);  // multiply with scalar
  LocalMatrix  operator/ (const LocalMatrix &M); // divide elementwise by matrix (or vector)
  LocalMatrix  operator/ (const FLOAT &s);  // divide by scalar
  LocalMatrix  operator- ();                // negative matrix
  LocalMatrix  operator~ ();                // transpose
  FLOAT   l2norm ();                   // euclidean norm (vectors) / frobenius norm (matrices)
  FLOAT   mean ();                     // mean of all elements in matrix

  // complex arithmetic operations
  static LocalMatrix cross (const LocalMatrix &a, const LocalMatrix &b);    // cross product of two vectors
  static LocalMatrix inv (const LocalMatrix &M);                       // invert matrix M
  bool   inv ();                                             // invert this matrix
  FLOAT  det ();                                             // returns determinant of matrix
  bool   solve (const LocalMatrix &M,FLOAT eps=1e-20);            // solve linear system M*x=B, replaces *this and M
  bool   lu(int32_t *idx, FLOAT &d, FLOAT eps=1e-20);        // replace *this by lower upper decomposition
  void   svd(LocalMatrix &U, LocalMatrix &W, LocalMatrix &V);                 // singular value decomposition *this = U*diag(W)*V^T

  // print matrix to stream
  friend std::ostream& operator<< (std::ostream& out,const LocalMatrix& M);

  // direct data access
  FLOAT   **val;
  int32_t   m,n;

private:

  void allocateMemory (const int32_t m_,const int32_t n_);
  void releaseMemory ();
  inline FLOAT pythag(FLOAT a,FLOAT b);

};


#endif