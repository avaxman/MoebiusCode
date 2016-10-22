//
//  Deform3D.h
//  testigl
//
//  Created by Amir Vaxman on 19/08/14.
//  Copyright (c) 2014 Amir Vaxman. All rights reserved.
//

#ifndef testigl_Deform3D_h
#define testigl_Deform3D_h


#include <Eigen/Dense>
#include <Eigen/Sparse>
#include "DeformTraits.h"
#include "PrescribeEdgeJumps.h"
#include "QuadConstSolver.h"
#include <hedra/EigenSolverWrapper.h>

using namespace Eigen;
using namespace std;

void ComputeCR(const MatrixXd& Vq, const MatrixXi& D, const MatrixXi& F, const MatrixXi& QuadVertexIndices, MatrixXd& ECR, MatrixXd& FCR);


class MoebiusDeformation3D{
public:
    
    MatrixXi tF;     //triangulated mesh (only for comparisons and visualization)
    VectorXi FromFace;  //the original quad
    MatrixXi EdgeF;  //the edge mesh for MC/IAP error etc. visualization
    MatrixXd EdgeOrigV;
    MatrixXd EdgeDeformV;
    MatrixXd EdgeInterpV;
    
    MatrixXi QuadVertexIndices;  //a list of all four consecutive vertices, corresponding to edge-mesh, that on which the cross ratio is defined
    
    MatrixXi D;  //face degrees
    MatrixXi F;  //general mesh
    MatrixXi E2V,E2F,F2E, E2Fi;
    MatrixXi ExtE2V;  //extended with diagonals
    MatrixXd F2ESigns;
    VectorXi InnerEdges;  //indices of edges which are non-boundary
    
    MatrixXi FaceCornerPairs;   //faces (f,g, i,k) and vertices around every edge for compatibility and smoothness, like in the paper
    MatrixXi CornerPairs;  //neighboring corner pairs across edges for AMAP/MC/IAP comparisons
    
    //Raw 3D positions
    MatrixXd OrigV;    //original vertex positions
    MatrixXd DeformV;  //current deformed mesh
    MatrixXd InterpV;  //current interpolated mesh
    
    //Quaternionic Positions
    MatrixXd OrigVq;
    MatrixXd DeformVq;
    MatrixXd InterpVq;
    
    //reciprocal corner variables for deformations
    int NumCorners;  //# of all corners in the mesh by order
    MatrixXd DeformX;
    VectorXi CornerOffset;  //where does every face begin in the corner list
    
    //of positional handles
    VectorXi ConstIndices;
    MatrixXd QuatConstPoses;
    
    //edge-based cross ratios
    MatrixXd OrigECR;
    MatrixXd DeformECR;
    MatrixXd InterpECR;
    
    //face-based cross ratios (if not triangular, otherwise 0)
    MatrixXd OrigFCR;
    MatrixXd DeformFCR;
    MatrixXd InterpFCR;
    
    MatrixXd EdgeCompCoeffs;
    
    VectorXd DeformConvErrors; //last process convergence errors
    VectorXd InterpConvErrors;
    
    //optimization operators
    DeformTraitsCornerVars3D DeformTraits;
    QuadConstSolver<DeformTraitsCornerVars3D> DeformSolver;
    
    PrescribeEdgeJumps3D InterpTraits;
    QuadConstSolver<PrescribeEdgeJumps3D> InterpSolver;
    
    hedra::optimization::EigenSolverWrapper<Eigen::SparseQR<Eigen::SparseMatrix<double, Eigen::ColMajor>, Eigen::COLAMDOrdering<int> > > d0Solver;
    
    Vector4i DeformGlobalVertices;    //four points defining the global transformation
    MatrixXd DeformGlobalMoebius;  //the global mobius (c,d)
    
    SparseMatrix<double> d0;
    SparseMatrix<double> d0NoFirst;
    MatrixXd RhsAdd;  //-d0*First OrigVq;
    
    void SetupMesh(const MatrixXd& InV, const MatrixXi& InD, const MatrixXi& InF);
    void InitDeformation(const VectorXi& InConstIndices, bool isExactMC, bool isExactIAP, double RigidRatio);
    void UpdateDeformation(const MatrixXd& ConstPoses, int MaxIterations);
    
    void SetupInterpolation(bool isExactMC, bool isExactIAP, bool isFromDeformed);
    void Interpolate(double t, int NumIterations);
    

};


#endif
