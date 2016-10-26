//
//  Deform2D.h
//  testigl
//
//  Created by Amir Vaxman on 14/08/14.
//  Copyright (c) 2014 Amir Vaxman. All rights reserved.
//

#ifndef MoebiusCode_Deform2D_h
#define MoebiusCode_Deform2D_h

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <hedra/Moebius2DEdgeDeviationTraits.h>
#include <hedra/EigenSolverWrapper.h>
#include <hedra/LMSolver.h>
#include "PrescribeEdgeJumps.h"
#include "QuadConstSolver.h"

using namespace Eigen;
using namespace std;

typedef std::complex<double> Complex;


void Complex2Coords(const VectorXcd& CV, MatrixXd& V);
void Coords2Complex(const MatrixXd& V, VectorXcd& CV);
void ComputeCR(const VectorXcd& Vc, const MatrixXi& D, const MatrixXi& F, const MatrixXi& QuadVertexIndices, VectorXcd& ECR, VectorXcd& FCR);


class MoebiusDeformation2D{
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
    
    //Raw 3D positions
    MatrixXd OrigV;    //original vertex positions
    MatrixXd DeformV;  //current deformed mesh
    MatrixXd InterpV;  //current interpolated mesh
    
    //Complex Positions
    VectorXcd OrigVc;
    VectorXcd DeformVc;
    VectorXcd InterpVc;
    
    
    //Edge deviation method variables
    VectorXcd DeformY;  //candidate reciprocals of the Mobius transformation.
    VectorXcd DeformE;  //errors of vertex reciprocals vs. actual transformed vertices
   
       
    //for inversion control and other usages
    SparseMatrix<Complex> d0;

    //of positional handles
    VectorXi ConstIndices;
    VectorXcd ComplexConstPoses;

    
    //edge-based cross ratios
    VectorXcd OrigECR;
    VectorXcd DeformECR;
    VectorXcd InterpECR;
    
    //face-based cross ratios (if not triangular, otherwise 0
    VectorXcd OrigFCR;
    VectorXcd DeformFCR;
    VectorXcd InterpFCR;
    
    VectorXd DeformConvErrors; //last process convergence errors
    VectorXd InterpConvErrors;
    
    //optimization operators
    //optimization operators
    hedra::optimization::EigenSolverWrapper<Eigen::SimplicialLLT<Eigen::SparseMatrix<double> > > DeformLinearSolver;
    hedra::optimization::Moebius2DEdgeDeviationTraits DeformTraits;
    hedra::optimization::LMSolver<hedra::optimization::EigenSolverWrapper<Eigen::SimplicialLLT<Eigen::SparseMatrix<double> > >,hedra::optimization::Moebius2DEdgeDeviationTraits> DeformSolver;


    PrescribeEdgeJumps2D InterpTraits;
    QuadConstSolver<PrescribeEdgeJumps2D> InterpSolver;
    
    Vector3i DeformGlobalVertices;    //three points going to others
    Vector2cd DeformGlobalMoebius;  //the global mobius (c,d)
    

    void SetupMesh(const MatrixXd& InV, const MatrixXi& InD, const MatrixXi& InF);
    void InitDeformation(const VectorXi& InConstIndices, bool isExactMC, bool isExactIAP, double RigidRatio);
    void UpdateDeformation(const MatrixXd& ConstPoses, int MaxIterations,bool isExactMC,bool isExactIAP);
    
    void SetupInterpolation(bool isExactMC, bool isExactIAP);
    void Interpolate(double t, int NumIterations);
    
};



#endif
