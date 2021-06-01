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

#include "VectorMath.h"

namespace libRSF
{
  void RemoveColumn(Eigen::MatrixXd& Matrix, unsigned int ColToRemove)
  {
      unsigned int numRows = Matrix.rows();
      unsigned int numCols = Matrix.cols()-1;

      if( ColToRemove < numCols )
          Matrix.block(0,ColToRemove,numRows,numCols-ColToRemove) = Matrix.rightCols(numCols-ColToRemove);

      Matrix.conservativeResize(numRows,numCols);
  }

  double Median(std::vector<double> &V)
  {
      size_t n = V.size() / 2;
      std::nth_element(V.begin(), V.begin()+n, V.end());
      return V[n];
  }

  double Median(Vector V)
  {
    std::vector<double> Vec(V.data(), V.data() + V.rows() * V.cols());
    return Median(Vec);
  }

  double MAD(Vector V)
  {
    return Median((V.array()- Median(V)).abs().matrix());
  }
}
