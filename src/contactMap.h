/*
 * contactMap.h
 *
 *  Created on: 4 Feb 2018
 *      Author: yanlin
 */

#ifndef CONTACTMAP_H_
#define CONTACTMAP_H_
#include "inlineFunc.h"
#include "readOptions.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
using namespace std;

#define MIN(a, b) ((a) > (b)) ? (b) : (a)
#define MAX(a, b) ((a) > (b)) ? (a) : (b)

template <class T>
class MyMatrix
{
public:
  int firstRow, firstCol, lastRow, lastCol; // inclusive
  T **data;

  MyMatrix(){};
  MyMatrix(int fr, int fc, int lr, int lc){};
  inline int startRow()
  {
    return firstRow;
  }
  inline int endRow()
  {
    return lastRow + 1;
  }
  inline int startCol(const int &row)
  {
    return -999;
  }
  inline int endCol(const int &row)
  {
    return -999;
  }
  inline T get(const int &i, const int &j){};
  inline void set(const int &i, const int &j, const T &x){};
  inline void add(const int &i, const int &j, const T &x){};
  void clear(){};
};

template <class T>
class UpperDiag : public MyMatrix<T>
{
public:
  UpperDiag()
  {
  }

  inline int startCol(const int &row)
  {
    return MAX(row, this->firstCol);
  }
  inline int endCol(const int &row)
  {
    return MIN(this->lastCol + 1, row + bandSize);
  }

  UpperDiag(int fr, int fc, int lr, int lc)
  {

    this->firstRow = fr;
    this->firstCol = fc;
    this->lastRow = lr;
    this->lastCol = lc;
    long int nalloc = 0;
    this->data = (T **)malloc(
        (this->lastRow - this->firstRow + 1) * sizeof(T *));
    for (int i = this->firstRow; i <= this->lastRow; i++)
    {
      int s = (endCol(i) - startCol(i));
      this->data[i - this->firstRow] = (T *)malloc(s * sizeof(T));
      for(int j=0;j<s;j++)
      this->data[i - this->firstRow][j]=0;
      // memset(this->data[i - this->firstRow], 0, s * sizeof(T));
      nalloc += s * sizeof(T);
    }
  }

  inline int increase(int i, int j)
  {
    int tmp;
    if (j < i)
    {
      tmp = i;
      i = j;
      j = tmp;
    }
    this->data[i - this->firstRow][j - startCol(i)]++;
    return 0;
  }

  inline T get(const int &i, const int &j)
  {
    //#ifdef SAFE
    //    if (i > j || i < this->startRow() || i >= this->endRow()
    //        || j < this->startCol(i) || j >= this->endCol(i)) {
    //      fprintf(stderr,
    //          "ERROR: Getting element outside matrix i=%d j=%d %d %d %d %d\n",
    //          i, j, this->startRow(), this->endRow(), this->startCol(i),
    //          this->endCol(i));
    //      while (1)
    //        ;
    //    }
    //#endif

    if (j < i)
      return this->data[j - this->firstRow][i - startCol(j)];
    return this->data[i - this->firstRow][j - startCol(i)];
  }

  void set(const int &i, const int &j, const T &x)
  {

// #ifdef SAFE
//     if (i > j || i < this->startRow() || i >= this->endRow() || j < this->startCol(i) || j >= this->endCol(i))
//     {
//       fprintf(stderr,
//               "ERROR: Setting element outside matrix i=%d j=%d %d %d %d %d\n",
//               i, j, this->startRow(), this->endRow(), this->startCol(i),
//               this->endCol(i));
//       while (1)
//         ;
//     }
// #endif

    this->data[i - this->firstRow][j - startCol(i)] = x;
  }

  inline void add(const int &i, const int &j, const T &x)
  {

// #ifdef SAFE
//     if (i > j || i < this->startRow() || i >= this->endRow() || j < this->startCol(i) || j >= this->endCol(i))
//     {
//       fprintf(stderr,
//               "ERROR: Setting element outside matrix i=%d j=%d %d %d %d %d\n",
//               i, j, this->startRow(), this->endRow(), this->startCol(i),
//               this->endCol(i));
//       while (1)
//         ;
//     }
// #endif

    this->data[i - this->firstRow][j - startCol(i)] += x;
  }

  inline void clear()
  {
    for (int i = this->firstRow; i <= this->lastRow; i++)
    {
      free(this->data[i - this->firstRow]);
    }
    free(this->data);
  }
};

extern UpperDiag<float> O;
extern UpperDiag<float> testO;
extern UpperDiag<float> OoverBias;
extern UpperDiag<float> T;
extern UpperDiag<char> horizontalBoundary;
extern UpperDiag<char> verticalBoundary;

#endif /* CONTACTMAP_H_ */
