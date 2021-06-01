/***************************************************************************
 * libRSF - A Robust Sensor Fusion Library
 *
 * Copyright (C) 2019 Chair of Automation Technology / TU Chemnitz
 * For more information see https://www.tu-chemnitz.de/etit/proaut/self-tuning
 *
 * libRSF is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * libRSF is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with libRSF.  If not, see <http://www.gnu.org/licenses/>.
 *
 * Author: Tim Pfeifer (tim.pfeifer@etit.tu-chemnitz.de)
 ***************************************************************************/

/**
 * @file VectorMath.h
 * @author Tim Pfeifer
 * @date 07.03.2019
 * @brief Derived vector type and some simple helper function.
 * @copyright GNU Public License.
 *
 */

#ifndef VECTORMATH_H
#define VECTORMATH_H

#include <algorithm>

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <ceres/ceres.h>

namespace libRSF
{
  /** define vector and matrix types for libRSF */
  typedef Eigen::Matrix<double, Eigen::Dynamic, 1> Vector;
  typedef Eigen::Matrix<double, 1, 1> Vector1;
  typedef Eigen::Matrix<double, 2, 1> Vector2;
  typedef Eigen::Matrix<double, 3, 1> Vector3;
  typedef Eigen::Matrix<double, 4, 1> Vector4;
  typedef Eigen::Matrix<double, 6, 1> Vector6;
  typedef Eigen::Matrix<double, 9, 1> Vector9;
  typedef Eigen::Matrix<double, 12, 1> Vector12;
  typedef Eigen::Matrix<double, 15, 1> Vector15;

  typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> Matrix;
  typedef Eigen::Matrix<double, 1, 1, Eigen::RowMajor> Matrix11;
  typedef Eigen::Matrix<double, 2, 2, Eigen::RowMajor> Matrix22;
  typedef Eigen::Matrix<double, 3, 3, Eigen::RowMajor> Matrix33;
  typedef Eigen::Matrix<double, 4, 4, Eigen::RowMajor> Matrix44;

  typedef Eigen::Matrix<double, 3, Eigen::Dynamic, Eigen::RowMajor> Matrix3X;

  /** math with vector/matrix class */
  template <int Dimension, typename T>
  Eigen::Matrix<T, Dimension, Dimension> SquareRoot(Eigen::Matrix<T, Dimension, Dimension> A)
  {
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix<T, Dimension, Dimension>> ES(A);
    return ES.operatorSqrt();
  }

  template <int Dimension, typename T>
  Eigen::Matrix<T, Dimension, Dimension> InverseSquareRoot(Eigen::Matrix<T, Dimension, Dimension> A)
  {
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix<T, Dimension, Dimension>> ES(A);
    return ES.operatorInverseSqrt();
  }

  /** vector math with raw pointers (compatible to ceres jet type) */
  template <int Dimension, typename T1, typename T2>
  void VectorDifference(const T1* const Vector1, const T2* const Vector2,  T1* Difference)
  {
    Eigen::Map<const Eigen::Matrix<T1, Dimension, 1> > V1(Vector1);
    Eigen::Map<const Eigen::Matrix<T2, Dimension, 1> > V2(Vector2);
    Eigen::Map<Eigen::Matrix<T1, Dimension, 1> > Diff(Difference);

    Diff = V1 - V2;
  }

  template <int Dimension, typename T>
  T VectorLength(const T* const Vector)
  {
    Eigen::Map<const Eigen::Matrix<T, Dimension, 1> > V(Vector);
    T SquaredNorm = V.squaredNorm();

    /** for stability of the derivation */
    if(SquaredNorm < T(1e-10))
      SquaredNorm += 1e-10;

    return ceres::sqrt(SquaredNorm);
  }

  template <int Dimension, typename T1, typename T2>
  T1 VectorDistance(const T1* const Vector1, const T2* const Vector2)
  {
    T1 Difference[Dimension];
    VectorDifference<Dimension, T1, T2>(Vector1, Vector2, Difference);
    return VectorLength<Dimension, T1>(Difference);
  }

  void RemoveColumn (Eigen::MatrixXd& Matrix, unsigned int ColToRemove);

  double Median(std::vector<double> &V);

  double Median(Vector V);

  double MAD(Vector V);
}

#endif // VECTORMATH_H
