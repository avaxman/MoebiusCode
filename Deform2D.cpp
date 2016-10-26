//
//  Deform2D.cpp
//  testigl
//
//  Created by Amir Vaxman on 14/08/14.
//  Copyright (c) 2014 Amir Vaxman. All rights reserved.
//

#include "Deform2D.h"
#include <hedra/polygonal_edge_topology.h>
#include <hedra/Moebius2DEdgeDeviationTraits.h>
#include <hedra/EigenSolverWrapper.h>
#include <hedra/LMSolver.h>
#include <hedra/triangulate_mesh.h>
#include <hedra/check_traits.h>
#include <igl/colon.h>
#include <igl/setdiff.h>
#include <igl/slice.h>
#include <igl/slice_into.h>
#include <igl/speye.h>
#include "PrescribeEdgeJumps.h"
#include "ExtraFunctions.h"
#include <set>
#include <iostream>

#define PRINT_OUT 0


void Complex2Coords(const VectorXcd& CV, MatrixXd& V)
{
    V.resize(CV.rows(), 3);
    for (int i=0;i<CV.rows();i++)
        V.row(i)<<CV(i).real(), CV(i).imag(),0.0;
}


void Coords2Complex(const MatrixXd& V, VectorXcd& CV)
{
    CV.resize(V.rows());
    for (int i=0;i<V.rows();i++)
        CV(i)=std::complex<double>(V(i,0), V(i,1));
}



void ConstructEFi(const MatrixXi& D, const MatrixXi& F2E, const MatrixXi& E2F, MatrixXi& E2Fi, MatrixXcd& F2ESigns, VectorXi& InnerEdges)
{
    //cout<<"F: "<<F<<endl;
    //cout<<"F2E: "<<F2E<<endl;
    //cout<<"E2F: "<<E2F<<endl;
    std::vector<int> InnerEdgesVec;
    E2Fi=MatrixXi::Constant(E2F.rows(), 2,-1);
    F2ESigns=MatrixXcd::Zero(F2E.rows(),F2E.cols());
    for (int i=0;i<E2F.rows();i++)
        for (int k=0;k<2;k++){
            if (E2F(i,k)==-1)
                continue;
            
            for (int j=0;j<D(E2F(i,k));j++)
                if (F2E(E2F(i,k),j)==i)
                    E2Fi(i,k)=j;
        }
    
    //cout<<"E2Fi: "<<E2Fi<<endl;
    
    //doing edge signs
    for (int i=0;i<E2F.rows();i++){
        if (E2Fi(i,0)!=-1) F2ESigns(E2F(i,0),E2Fi(i,0))=1.0;
        if (E2Fi(i,1)!=-1) F2ESigns(E2F(i,1),E2Fi(i,1))=-1.0;
        if ((E2F(i,0)!=-1)&&(E2F(i,1)!=-1))
            InnerEdgesVec.push_back(i);
    }
    
    //cout<<"F2ESigns :"<<F2ESigns<<endl;
    
    InnerEdges.resize(InnerEdgesVec.size());
    for (int i=0;i<InnerEdgesVec.size();i++)
        InnerEdges(i)=InnerEdgesVec[i];
    
    //cout<<"InnerEdges: "<<InnerEdges<<endl;
}


//this computes the cross ratio on every edge, and on every face of the mesh
void ComputeCR(const VectorXcd& Vc, const MatrixXi& D, const MatrixXi& F, const MatrixXi& QuadVertexIndices, VectorXcd& ECR, VectorXcd& FCR)
{
    
    ECR.resize(QuadVertexIndices.rows());
    
    for (int i=0;i<QuadVertexIndices.rows();i++){
        Complex Zi=Vc(QuadVertexIndices(i,0));
        Complex Zj=Vc(QuadVertexIndices(i,1));
        Complex Zk=Vc(QuadVertexIndices(i,2));
        Complex Zl=Vc(QuadVertexIndices(i,3));
        
        ECR(i)=((Zj-Zi)*(Zl-Zk))/((Zi-Zl)*(Zk-Zj));
    }
    
    int NumFCR=F.col(0).sum()-3*D.rows();
    FCR.resize(NumFCR); FCR.setZero();
    int CurrFCR=0;
    for (int i=0;i<D.rows();i++){
        for (int j=0;j<D(i)-3;j++){
            Complex Zi=Vc(F(i,j));
            Complex Zj=Vc(F(i,j+1));
            Complex Zk=Vc(F(i,j+2));
            Complex Zl=Vc(F(i,j+3));
            
            FCR(CurrFCR++)=((Zj-Zi)*(Zl-Zk))/((Zi-Zl)*(Zk-Zj));
        }
    }
    
}





void MoebiusDeformation2D::SetupMesh(const MatrixXd& InV,  const MatrixXi& InD, const MatrixXi& InF)
{
    F=InF;
    D=InD;
    hedra::triangulate_mesh(D, F, tF, FromFace);
    std::cout<<"Finished Triangulation"<<std::endl;
    OrigV=InV;
    DeformV=InV;
    InterpV=InV;
    
    DeformConvErrors.resize(1);
    DeformConvErrors(0)=0.0;
    
    InterpConvErrors.resize(1);
    InterpConvErrors(0)=0.0;
    
    //complex representations
    Coords2Complex(OrigV,OrigVc);
    DeformVc=OrigVc;
    InterpVc=OrigVc;
    
    
    hedra::polygonal_edge_topology(D, F, E2V, F2E, E2F, E2Fi, F2ESigns, InnerEdges);
    //igl::general_edge_topology(D, F, E2V,F2E,E2F);
    //ConstructEFi(D, F2E, E2F, E2Fi, F2ESigns, InnerEdges);
    
    //creating full edge list
    vector<pair<int,int>>  Diagonals;
    for (int i=0;i<F.rows();i++)
        for (int j=0;j<D(i);j++)
            for (int k=j+2;k<D(i);k++)
                Diagonals.push_back(pair<int,int>(F(i,j),F(i,k)));
    
    ExtE2V.resize(E2V.rows()+Diagonals.size(),2);
    ExtE2V.block(0,0,E2V.rows(),2)=E2V;
    for (int i=0;i<Diagonals.size();i++)
        ExtE2V.row(E2V.rows()+i)<<Diagonals[i].first, Diagonals[i].second;
    
    
    
    DeformY=VectorXcd::Ones(DeformVc.rows());
    DeformE=VectorXcd::Ones(ExtE2V.rows());
    
    //setting up the full differential operators
    std::vector<Triplet<Complex> > d0tris;
    d0.resize(ExtE2V.rows(),OrigV.rows());
    for (int i=0;i<ExtE2V.rows();i++)
    {
        d0tris.push_back(Triplet<Complex>(i,ExtE2V(i,0),-1.0));
        d0tris.push_back(Triplet<Complex>(i,ExtE2V(i,1),1.0));
    }
    d0.setFromTriplets(d0tris.begin(), d0tris.end());
    
    //building edge mesh
    
    //center vertices
    EdgeOrigV.resize(OrigV.rows()+F.rows(),3);
    
    MatrixXd Centers(F.rows(),3); Centers.setZero();
    for (int i=0;i<F.rows();i++)
        for (int j=0;j<D(i);j++)
            Centers.row(i)+=OrigV.row(F(i,j))/(double)D(i);
    
    
    EdgeOrigV<<OrigV, Centers;
    
    EdgeDeformV=EdgeOrigV;
    EdgeInterpV=EdgeOrigV;
    
    EdgeF.resize(2*InnerEdges.rows(),3);
    QuadVertexIndices.resize(InnerEdges.size(),4);
    for (int i=0;i<InnerEdges.rows();i++){
        int f=E2F(InnerEdges(i),0);
        int g=E2F(InnerEdges(i),1);
        int vi=E2V(InnerEdges(i),0);
        int vk=E2V(InnerEdges(i),1);
        int vj=F(g,1+(E2Fi(InnerEdges(i),1)+2)%(F(g,0)));
        int vl=F(f,1+(E2Fi(InnerEdges(i),0)+2)%(F(f,0)));
        int vf=OrigV.rows()+E2F(InnerEdges(i),0);
        int vg=OrigV.rows()+E2F(InnerEdges(i),1);
        
        QuadVertexIndices.row(i)<<vi,vj,vk,vl;
        EdgeF.row(i)<<vi,vk,vf;
        EdgeF.row(i+InnerEdges.rows())<<vk,vi,vg;
    }
    
    //cout<<"QuadVertexIndices: "<<QuadVertexIndices<<endl;
    
    
    FaceCornerPairs.resize(InnerEdges.size(),4);
    for (int i=0;i<InnerEdges.size();i++)
        FaceCornerPairs.row(i)<<E2F.row(InnerEdges(i)), E2V(InnerEdges(i),1), E2V(InnerEdges(i),0);
    
    //cout<<"FaceCornerPairs: "<<FaceCornerPairs<<endl;
    
    
    ComputeCR(OrigVc, D, F, QuadVertexIndices, OrigECR, OrigFCR);
    DeformECR=OrigECR;
    InterpECR=OrigECR;
    DeformFCR=OrigFCR;
    InterpFCR=OrigFCR;
    
    std::cout<<"Finished regular init"<<std::endl;
    
    DeformTraits.init(OrigVc, D, F, E2V, false, false);
    
    std::cout<<"Finished DeformTraits.init"<<std::endl;
    DeformSolver.init(&DeformLinearSolver, &DeformTraits, 150);
    
    std::cout<<"Finished DeformSolver.init"<<std::endl;
    
}




void MoebiusDeformation2D::InitDeformation(const VectorXi& InConstIndices, bool isExactMC, bool isExactIAP, double RigidRatio)
{
 
    ConstIndices=InConstIndices;
    
    DeformTraits.init(OrigVc, D, F, E2V, isExactMC, isExactIAP, ConstIndices);
    DeformTraits.rigidRatio=RigidRatio;
    DeformSolver.init(&DeformLinearSolver, &DeformTraits, 150);
    
    //checking traits
    if (isExactMC || isExactIAP){
        DeformTraits.initSolution=VectorXcd::Random(OrigVc.rows()+DeformY.rows()+DeformE.rows());
        DeformTraits.complexConstPoses=VectorXcd::Random(InConstIndices.size());
        DeformTraits.smoothFactor=100.0;
        DeformTraits.posFactor=10.0;
        hedra::optimization::check_traits<hedra::optimization::Moebius2DEdgeDeviationTraits>(DeformTraits);
     }else{
         DeformTraits.initSolution=VectorXcd::Random(OrigVc.rows()+DeformY.rows());
         DeformTraits.complexConstPoses=VectorXcd::Random(InConstIndices.size());
         DeformTraits.smoothFactor=100.0;
         DeformTraits.posFactor=10.0;
         hedra::optimization::check_traits<hedra::optimization::Moebius2DEdgeDeviationTraits>(DeformTraits);
     }
    
    DeformSolver.init(&DeformLinearSolver, &DeformTraits, 150);
    
}




void MoebiusDeformation2D::UpdateDeformation(const MatrixXd& ConstPoses, int MaxIterations,bool isExactMC,bool isExactIAP)
{
    
    Coords2Complex(ConstPoses, ComplexConstPoses);
    DeformTraits.complexConstPoses=ComplexConstPoses;
    DeformTraits.smoothFactor=100.0;
    DeformTraits.posFactor=10.0;
    DeformSolver.solve(true);
    DeformVc=DeformTraits.finalPositions;
    DeformY=DeformTraits.finalY;
    DeformE=DeformTraits.finalE;

    Complex2Coords(DeformVc, DeformV);
    ComputeCR(DeformVc, D, F,QuadVertexIndices, DeformECR, DeformFCR);
    
    EdgeDeformV<<DeformV, GetCenters(DeformV, D, F);
    
}


//must be called before interpolation when the deformed mesh has changed.
void MoebiusDeformation2D::SetupInterpolation(bool isExactMC, bool isExactIAP)
{
    
    InterpTraits.Initialize(OrigVc,  D, F, FaceCornerPairs,isExactMC, isExactIAP);
    
    //checking traits
    //InterpTraits.InitSolution=VectorXcd::Random(2*F.rows());
    //InterpTraits.SmoothFactor=10.0;
    //InterpTraits.PresJumps=VectorXcd::Random(FaceCornerPairs.rows());
    //CheckTraits<PrescribeEdgeJumps2D>(InterpTraits, 4*F.rows());
    
    
    InterpSolver.Initialize(&InterpTraits);
    
    //computing the global Moebius transformations between three points
    VectorXd Distances(OrigVc.rows());
    for (int i=0;i<OrigVc.rows();i++)
        Distances(i)=abs(OrigVc(i)-DeformVc(i));
    
    VectorXi SortedIndices;
    VectorXd SortedDistances;
    
    igl::sort(Distances, 1,true, SortedDistances, SortedIndices);
    
    DeformGlobalVertices=SortedIndices.head(3);
    
    Vector3cd z,w;
    for (int i=0;i<3;i++){
        z(i)=OrigVc(DeformGlobalVertices(i));
        w(i)=DeformVc(DeformGlobalVertices(i));
    }
    
    Complex a,b,c,d;
    GetComplexMobiusCoeffs(a, b, c, d, z,w);
    
    //signing properly
    if (d.real()<0.0){
        c=-c; d=-d;
    }
    
    DeformGlobalMoebius<<c,d;
    
    //testing deformglobalmoebius
    Vector3cd TestGlobalPoints;
    Vector2cd TestGlobalVectors;  //z1-z0, z2-z0
    TestGlobalPoints(0)=DeformVc(DeformGlobalVertices(0));
    
    Vector3cd CornerValues=(DeformGlobalMoebius(0)*z).array()+DeformGlobalMoebius(1);
    
    TestGlobalVectors<<(z(1)-z(0))/(CornerValues(0)*CornerValues(1)), (z(2)-z(0))/(CornerValues(0)*CornerValues(2));
    TestGlobalPoints(1)=TestGlobalPoints(0)+TestGlobalVectors(0);
    TestGlobalPoints(2)=TestGlobalPoints(0)+TestGlobalVectors(1);
    //cout<<"TestGlobalPoints: "<<TestGlobalPoints<<endl;
    //cout<<"w: "<<w<<endl;
    
}



void MoebiusDeformation2D::Interpolate(double t, int NumIterations)
{
    
    cout<<"Solving for integrable Edge Jump Coefficiengs"<<endl;
    cout<<"+++++++++++++++++++++++++++++++++++++++++++++"<<endl;
    
    VectorXcd PresJumps=exp(log(DeformECR.cwiseQuotient(OrigECR).array())*t/2.0);   //this assumes that the difference is small enough for logarithem to always be in the principal branch;
    InterpTraits.PresJumps=PresJumps;
    
    
    //creating matrix with current PresJump Values
    vector<Triplet<Complex> > JumpMatTris(8*FaceCornerPairs.size());
    for (int i=0;i<FaceCornerPairs.rows();i++){
        JumpMatTris[8*i]=Triplet<Complex>(2*i,2*FaceCornerPairs(i,0), OrigVc(FaceCornerPairs(i,2))*PresJumps(i));
        JumpMatTris[8*i+1]=Triplet<Complex>(2*i,2*FaceCornerPairs(i,0)+1, PresJumps(i));
        JumpMatTris[8*i+2]=Triplet<Complex>(2*i,2*FaceCornerPairs(i,1), -OrigVc(FaceCornerPairs(i,2)));
        JumpMatTris[8*i+3]=Triplet<Complex>(2*i,2*FaceCornerPairs(i,1)+1, -1.0);
        
        JumpMatTris[8*i+4]=Triplet<Complex>(2*i+1,2*FaceCornerPairs(i,0), OrigVc(FaceCornerPairs(i,3)));
        JumpMatTris[8*i+5]=Triplet<Complex>(2*i+1,2*FaceCornerPairs(i,0)+1, 1.0);
        JumpMatTris[8*i+6]=Triplet<Complex>(2*i+1,2*FaceCornerPairs(i,1), -OrigVc(FaceCornerPairs(i,3))*PresJumps(i));
        JumpMatTris[8*i+7]=Triplet<Complex>(2*i+1,2*FaceCornerPairs(i,1)+1, -PresJumps(i));
    }
    
    //can be in advance.
    SparseMatrix<Complex> JumpMat(2*FaceCornerPairs.rows(),2*F.rows());
    JumpMat.reserve(JumpMatTris.size());
    JumpMat.setFromTriplets(JumpMatTris.begin(), JumpMatTris.end());
    
    VectorXcd torhs(2*F.rows()); torhs.setZero();
    torhs(1)=1.0;
    VectorXcd Rhs=-JumpMat*torhs;
    
    VectorXi VarIndices(2*F.rows()-2);
    for (int i=0;i<VarIndices.size();i++)
        VarIndices(i)=i+2;
    
    SparseMatrix<Complex> JumpMatVar(2*FaceCornerPairs.rows(),2*(F.rows()-1));
    igl::slice(JumpMat, VarIndices, 2, JumpMatVar);
    
    VectorXcd InitialMobCoeffs=SolveComplexSytem(JumpMatVar.adjoint()*JumpMatVar, JumpMatVar.adjoint()*Rhs);
    
    cout<<"Initial MobCoeff from Jump Error: "<<(JumpMatVar*InitialMobCoeffs-Rhs).lpNorm<Infinity>()<<endl;
    
    VectorXcd InitSolution(2*F.rows());
    InitSolution<<0.0, 1.0, InitialMobCoeffs;
    
    VectorXd InitSolutionReal(2*InitSolution.size());
    InitSolutionReal<<InitSolution.real(), InitSolution.imag();
    
    InterpTraits.InitSolution=InitSolution;
    InterpTraits.SmoothFactor=10.0;
    InterpTraits.PresJumps=PresJumps;
    VectorXd IntegratedSolutionReal=InterpSolver.Solve(InitSolutionReal, 1.0, NumIterations);
    
    VectorXcd IntegratedSolution(IntegratedSolutionReal.size()/2);
    IntegratedSolution.real()=IntegratedSolutionReal.head(IntegratedSolutionReal.size()/2);
    IntegratedSolution.imag()=IntegratedSolutionReal.tail(IntegratedSolutionReal.size()/2);
    
    InterpConvErrors=InterpSolver.ConvErrors;
    
    
    //reconstructing the current interpolated mesh
    VectorXcd FaceEdgeVectors(E2F.rows());
    for (int i=0;i<E2F.rows();i++){
        FaceEdgeVectors(i)=0.0;
        if (E2F(i,0)!=-1){
            Complex c1=IntegratedSolution(2*E2F(i,0));
            Complex d1=IntegratedSolution(2*E2F(i,0)+1);
            FaceEdgeVectors(i)=(OrigVc(E2V(i,1))-OrigVc(E2V(i,0)))/((c1*OrigVc(E2V(i,0))+d1)*(c1*OrigVc(E2V(i,1))+d1));
        }
        if (E2F(i,1)!=-1){
            Complex c2=IntegratedSolution(2*E2F(i,1));
            Complex d2=IntegratedSolution(2*E2F(i,1)+1);
            FaceEdgeVectors(i)+=(OrigVc(E2V(i,1))-OrigVc(E2V(i,0)))/((c2*OrigVc(E2V(i,0))+d2)*(c2*OrigVc(E2V(i,1))+d2));
        }
        
        if ((E2F(i,0)!=-1)&&(E2F(i,1)!=-1))
            FaceEdgeVectors(i)/=2.0;
    }
    
    
    VarIndices.resize(OrigVc.rows()-1);
    for (int i=1;i<OrigVc.rows();i++)
        VarIndices(i-1)=i;
    
    torhs.resize(OrigVc.rows()); torhs.setZero();
    torhs(0)=OrigVc(0);
    Rhs=-d0*torhs;
    
    
    SparseMatrix<Complex> d0LocalVar(E2F.rows(), OrigV.rows()-1);
    igl::slice(d0, VarIndices, 2, d0LocalVar);
    
    
    VectorXcd CurrInterpVc(OrigVc.rows());
    CurrInterpVc<<OrigVc(0), SolveComplexSytem(d0LocalVar.adjoint()*d0LocalVar, d0LocalVar.adjoint()*(Rhs+FaceEdgeVectors));
    
    cout<<"Reconstruction Error: "<<(d0*CurrInterpVc-FaceEdgeVectors).lpNorm<Infinity>()<<endl;
    
    Vector3cd ReconPoints;
    ReconPoints<<CurrInterpVc(DeformGlobalVertices(0)),CurrInterpVc(DeformGlobalVertices(1)),CurrInterpVc(DeformGlobalVertices(2));
    
    Vector3cd OrigPoints;
    OrigPoints<<OrigVc(DeformGlobalVertices(0)),OrigVc(DeformGlobalVertices(1)),OrigVc(DeformGlobalVertices(2));
    
    
    
    Vector2cd InterpGlobalMoebius; InterpGlobalMoebius<<t*DeformGlobalMoebius(0), pow(DeformGlobalMoebius(1),t);
    //cout<<"InterpGlobalMoebius: "<<InterpGlobalMoebius<<endl;
    Vector3cd CornerValues=(InterpGlobalMoebius(0)*OrigPoints).array()+InterpGlobalMoebius(1);
    
    //interpolating between the t=0 and t=1 chosen global points. First point is linearly interpolated, and the rest are vectors
    Vector3cd TargetGlobalPoints;
    Vector2cd TargetGlobalVectors;  //z1-z0, z2-z0
    TargetGlobalPoints(0)=DeformVc(DeformGlobalVertices(0))*t+OrigVc(DeformGlobalVertices(0))*(1-t);
    
    TargetGlobalVectors<<(OrigPoints(1)-OrigPoints(0))/(CornerValues(0)*CornerValues(1)), (OrigPoints(2)-OrigPoints(0))/(CornerValues(0)*CornerValues(2));
    TargetGlobalPoints(1)=TargetGlobalPoints(0)+TargetGlobalVectors(0);
    TargetGlobalPoints(2)=TargetGlobalPoints(0)+TargetGlobalVectors(1);
    
    
    Complex a,b,c,d;
    GetComplexMobiusCoeffs(a, b, c, d, ReconPoints,TargetGlobalPoints);
    
    for (int i=0;i<CurrInterpVc.rows();i++)
        InterpVc(i)=(a*CurrInterpVc(i)+b)/(c*CurrInterpVc(i)+d);
    
    Complex2Coords(InterpVc, InterpV);
    ComputeCR(InterpVc, D, F, QuadVertexIndices, InterpECR, InterpFCR);
 
    EdgeInterpV<<InterpV, GetCenters(InterpV, D, F);
    
    InterpConvErrors=InterpSolver.ConvErrors;
}

