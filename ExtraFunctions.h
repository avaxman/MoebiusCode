//
//  ExtraFunctions.h
//  testigl
//
//  Created by Amir Vaxman on 18/12/14.
//  Copyright (c) 2014 Amir Vaxman. All rights reserved.
//

#ifndef testigl_ExtraFunctions_h
#define testigl_ExtraFunctions_h

#include <iostream>
#include <Eigen/Dense>
#include <Eigen/Sparse>

using namespace Eigen;
using namespace std;

typedef std::complex<double> Complex;

Vector2cd GetSingleMobius(const VectorXcd& OrigVc, const VectorXcd& X);

VectorXcd SolveComplexSytem(const SparseMatrix<Complex>& A, const VectorXcd& b);

void GetComplexMobiusCoeffs(Complex& a, Complex& b, Complex& c, Complex& d, const Vector3cd& z, const Vector3cd& w);

RowVector4d Rot2Quat(Matrix3d& R);


MatrixXd GetCenters(const MatrixXd& V, const MatrixXi& D, const MatrixXi& F);





#endif
