//
//  Deform3D.cpp
//  testigl
//
//  Created by Amir Vaxman on 19/08/14.
//  Copyright (c) 2014 Amir Vaxman. All rights reserved.
//

#include "Deform3D.h"
#include <igl/edge_topology.h>
#include <hedra/triangulate_mesh.h>
#include <hedra/polygonal_edge_topology.h>
#include <hedra/EigenSolverWrapper.h>
#include <hedra/Moebius3DCornerVarsTraits.h>
#include <hedra/quaternionic_operations.h>
#include <hedra/check_traits.h>
#include <hedra/LMSolver.h>
#include <hedra/triangulate_mesh.h>
#include <igl/colon.h>
#include <igl/setdiff.h>
#include <igl/slice.h>
#include <igl/slice_into.h>
#include <igl/speye.h>
#include <igl/colon.h>
#include "ExtraFunctions.h"
#include "PrescribeEdgeJumps.h"
#include "UnitSphereMoebiusTraits.h"
#include <set>


double ComputePlanarity(const RowVector3d& z1, const RowVector3d& z2,const RowVector3d& z3,const RowVector3d& z4)
{
    Matrix3d M;
    M<<z2-z1, z3-z1, z4-z1;
    double Planarity=M.determinant();
    Planarity/=((z1-z2).norm()*(z1-z3).norm()*(z1-z4).norm());
    return Planarity;
}



void Quat2Vector(const MatrixXd& QVec, VectorXd& RawVec)
{
    RawVec.resize(QVec.size());
    for (int i=0;i<QVec.rows();i++)
        RawVec.segment(4*i,4)=QConj(QVec.row(i)).transpose();
    
}

//this computes the cross ratio on every edge, and on every face of the mesh
void ComputeCR(const MatrixXd& Vq, const MatrixXi& D, const MatrixXi& F, const MatrixXi& QuadVertexIndices, MatrixXd& ECR, MatrixXd& FCR)
{
    
    ECR.resize(QuadVertexIndices.rows(),4);
    
    for (int i=0;i<QuadVertexIndices.rows();i++){
        RowVector4d qi=Vq.row(QuadVertexIndices(i,0));
        RowVector4d qj=Vq.row(QuadVertexIndices(i,1));
        RowVector4d qk=Vq.row(QuadVertexIndices(i,2));
        RowVector4d ql=Vq.row(QuadVertexIndices(i,3));
        
        ECR.row(i)=QMult(QMult(qj-qi, QInv(qk-qj)),QMult(ql-qk, QInv(qi-ql)));
        
    }
    
    int NumFCR=F.col(0).sum()-3*D.rows();
    FCR.resize(NumFCR,4); FCR.setZero();
    int CurrFCR=0;
    for (int i=0;i<D.rows();i++){
        for (int j=0;j<D(i)-3;j++){
            RowVector4d qi=Vq.row(F(i,j));
            RowVector4d qj=Vq.row(F(i,j+1));
            RowVector4d qk=Vq.row(F(i,j+2));
            RowVector4d ql=Vq.row(F(i,j+3));
            
            FCR.row(CurrFCR++)=QMult(QMult(qj-qi, QInv(qk-qj)),QMult(ql-qk, QInv(qi-ql)));
        }
    }
    
 //   cout<<"FCR: "<<FCR<<endl;
    
}




void MoebiusDeformation3D::SetupMesh(const MatrixXd& InV, const MatrixXi& InD, const MatrixXi& InF)
{
    F=InF;
    D=InD;
    hedra::triangulate_mesh(D, F, tF, FromFace);

    OrigV=InV;
    DeformV=InV;
    InterpV=InV;
    
    DeformConvErrors.resize(1);
    DeformConvErrors(0)=0.0;
    
    InterpConvErrors.resize(1);
    InterpConvErrors(0)=0.0;
    
    //complex representations
    Coords2Quat(OrigV,OrigVq);
    DeformVq=OrigVq;
    InterpVq=OrigVq;
    
    
    hedra::polygonal_edge_topology(D, F, E2V, F2E, E2F, E2Fi, F2ESigns, InnerEdges);
    //hedra::polygonal_edge_topology(D, F,  E2V,F2E,E2F);
    //ConstructEFi(D, F2E, E2F, E2Fi, F2ESigns, InnerEdges);
    
    //creating full edge list
    vector<pair<int,int>>  Diagonals;
    for (int i=0;i<F.rows();i++)
        for (int j=0;j<D(i);j++)
            for (int k=j+1;k<D(i);k++)
                Diagonals.push_back(pair<int,int>(F(i,j),F(i,k)));
    
    ExtE2V.resize(E2V.rows()+Diagonals.size(),2);
    ExtE2V.block(0,0,E2V.rows(),2)=E2V;
    for (int i=0;i<Diagonals.size();i++)
        ExtE2V.row(E2V.rows()+i)<<Diagonals[i].first, Diagonals[i].second;
    
    
    NumCorners=D.sum();
    DeformX.resize(NumCorners,4);
    DeformX.setZero();
    DeformX.col(0).setConstant(1.0);
    CornerOffset.resize(F.rows());
    CornerOffset(0)=0;
    
    for (int i=1;i<D.rows();i++)
        CornerOffset(i)=CornerOffset(i-1)+D(i-1);
    
    
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
        int vj=F(g,(E2Fi(InnerEdges(i),1)+2)%D(g));
        int vl=F(f,(E2Fi(InnerEdges(i),0)+2)%D(f));
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
    
    //constructing corner pairs
    MatrixXi AdjCorners(OrigV.rows(),12);
    VectorXi Valences(OrigV.rows()); Valences.setZero();
    for (int i=0;i<D.rows();i++){
        for (int j=0;j<D(i);j++){
            AdjCorners(F(i,j),Valences(F(i,j)))=CornerOffset(i)+j;
            Valences(F(i,j))++;
        }
    }
    
    //cout<<"F: "<<F<<endl;
    //cout<<"AdjCorners: "<<AdjCorners<<endl;
    
    vector<pair<int,int> > CornerPairList;
    for (int i=0;i<OrigV.rows();i++)
        for (int j=0;j<Valences(i)-1;j++)
            CornerPairList.push_back(pair<int,int>(AdjCorners(i,j),AdjCorners(i,(j+1)%Valences(i))));
    
    CornerPairs.resize(CornerPairList.size(),2);
    for (int i=0;i<CornerPairList.size();i++)
        CornerPairs.row(i)<<CornerPairList[i].first, CornerPairList[i].second;
    
    //cout<<"CornerPairs: "<<CornerPairs<<endl;
    
    
    ComputeCR(OrigVq, D, F, QuadVertexIndices, OrigECR, OrigFCR);
    DeformECR=OrigECR;
    InterpECR=OrigECR;
    DeformFCR=OrigFCR;
    InterpFCR=OrigFCR;
    
    //setting up the full differential operators
    std::vector<Triplet<double> > d0tris;
    d0.resize(E2V.rows(),OrigV.rows());
    VectorXi I,J;
    VectorXd S;
    for (int i=0;i<E2V.rows();i++)
    {
        d0tris.push_back(Triplet<double>(i,E2V(i,0),-1.0));
        d0tris.push_back(Triplet<double>(i,E2V(i,1),1.0));
    }
    d0.setFromTriplets(d0tris.begin(), d0tris.end());

    //setuping differential oepratot
    VectorXi LocalVarIndices(OrigVq.rows()-1);
    for (int i=0;i<OrigVq.rows()-1;i++)
        LocalVarIndices(i)=i+1;
    
    
    d0NoFirst.resize(E2V.rows(), OrigV.rows()-1);
    igl::slice(d0, LocalVarIndices, 2, d0NoFirst);
    SparseMatrix<double> d0td0=d0NoFirst.transpose()*d0NoFirst;
    int NumMembers=0;
    for (int k=0; k<d0td0.outerSize(); ++k)
        for (SparseMatrix<double>::InnerIterator it(d0td0,k); it; ++it)
            if (it.row()<=it.col())
                NumMembers++;
         
    I.resize(NumMembers);
    J.resize(NumMembers);
    S.resize(NumMembers);
    int Counter=0;
    for (int k=0; k<d0td0.outerSize(); ++k){
        for (SparseMatrix<double>::InnerIterator it(d0td0,k); it; ++it){
            if (it.row()<=it.col()){
                I(Counter)=it.row();
                J(Counter)=it.col();
                S(Counter++)=it.value();
            }
        }
    }
    
    //cout<<"I:"<<I<<endl;
    //cout<<"J:"<<J<<endl;
    
    d0Solver.analyze(I,J, false);
    if(!d0Solver.factorize(S, false))
        // decomposition failed
        cout<<"Solver Failed to factorize! "<<endl;
    
    DeformTraits.init(OrigV, D, F, false);
    DeformSolver.init(&DeformLinearSolver, &DeformTraits, 150);
}

//only effectively works for quad meshes

typedef pair<int,int> IntPair;

struct compare {
    bool operator()(const pair<int,int>& a, const pair<int,int>& b) {
        if (a.first < b.first) return true;
        else if ( (a.first == b.first) && (a.second < b.second) ) return true;
        else return false;
    }
};



void MoebiusDeformation3D::InitDeformation(const VectorXi& InConstIndices, bool isExactMC, bool isExactIAP, double RigidRatio)
{

    ConstIndices=InConstIndices;
    
    DeformTraits.init(OrigV, D, F, isExactMC, ConstIndices);
    DeformTraits.rigidRatio=RigidRatio;
    DeformSolver.init(&DeformLinearSolver, &DeformTraits, 150);
    
    //checking traits
    /*DeformTraits.initSolution=VectorXd::Random(4*NumCorners+3*OrigVq.rows());
    DeformTraits.constPoses=MatrixXd::Random(ConstIndices.size(),3);
    DeformTraits.smoothFactor=100.0;
    DeformTraits.posFactor=10.0;
    hedra::optimization::check_traits<hedra::optimization::MoebiusCornerVarsTraits>(DeformTraits);*/
 
}



void MoebiusDeformation3D::UpdateDeformation(const MatrixXd& ConstPoses, int MaxIterations)
{
    
    DeformTraits.constPoses=ConstPoses;
    DeformTraits.smoothFactor=100.0;
    DeformTraits.posFactor=10.0;
    DeformSolver.solve(true);
    DeformV=DeformTraits.finalPositions;
    DeformX=DeformTraits.finalX;
  
    Coords2Quat(DeformV, DeformVq);
    ComputeCR(DeformVq, D, F, QuadVertexIndices, DeformECR, DeformFCR);
    MatrixXd Centers(F.rows(),3); Centers.setZero();
    for (int i=0;i<F.rows();i++)
        for (int j=0;j<D(i);j++)
            Centers.row(i)+=DeformV.row(F(i,j))/(double)D(i);
    
    EdgeDeformV<<DeformV, Centers;
}



//assuming there are only four points!!
RowVector4d FindSphere(const MatrixXd& z){
    
    Matrix3d A;
    RowVector3d zMeans=z.colwise().mean();
    
    MatrixXd zCentered=z.rowwise()-zMeans;
    
    A<< (zCentered.col(0).cwiseProduct(z.col(0))).mean(),2*(zCentered.col(1).cwiseProduct(z.col(0))).mean(),2*(zCentered.col(2).cwiseProduct(z.col(0))).mean(),
        0.0, (zCentered.col(1).cwiseProduct(z.col(1))).mean(),2*(zCentered.col(2).cwiseProduct(z.col(1))).mean(),
        0.0, 0.0, (zCentered.col(2).cwiseProduct(z.col(2))).mean();
    
    Matrix3d AtA=A+A.transpose();
    
    Vector3d b;
    VectorXd zsqrnorm=((z.rowwise()).squaredNorm()).transpose();
    for (int i=0;i<3;i++)
        b(i)=(zsqrnorm.cwiseProduct(zCentered.col(i))).mean();
    
    RowVector4d Circle;
    
    Circle.tail(3)=(AtA.colPivHouseholderQr().solve(b)).transpose();
    Circle(0)=(((z.rowwise()-Circle.tail(3)).rowwise()).norm()).mean();
    
    //cout<<"Testing Find Sphere Function: "<<endl;
    //cout<<"Distances to center: "<<(z.rowwise()-Circle.tail(3)).rowwise().norm()<<endl;
    //cout<<"radius: "<<Circle(0)<<endl;
    
    
    return Circle;
    
}

void ComposeMobius(const RowVector4d& a1, const RowVector4d& b1, const RowVector4d& c1, const RowVector4d& d1, const RowVector4d& a2, const RowVector4d& b2, const RowVector4d& c2, const RowVector4d& d2, RowVector4d& a, RowVector4d& b, RowVector4d& c, RowVector4d& d){
    
    RowVector4d ta, tb, tc, td;  //so to be able to use same variables
    
    ta=QMult(a1, a2)+QMult(b1, c2);
    tb=QMult(a1, b2)+QMult(b1, d2);
    tc=QMult(c1, a2)+QMult(d1, c2);
    td=QMult(c1, b2)+QMult(d1, d2);
    
    a=ta;
    b=tb;
    c=tc;
    d=td;
}


void GetMobiusCoeffs(RowVector4d& a, RowVector4d& b, RowVector4d& c, RowVector4d& d, const MatrixXd& z, const MatrixXd& w)
{
    
    //finding spheres
    RowVector4d Sz=FindSphere(z);
    RowVector4d Sw=FindSphere(w);
    
    
    MatrixXd zunit=(z.rowwise()-Sz.tail(3))/Sz(0);
    MatrixXd wunit=(w.rowwise()-Sw.tail(3))/Sw(0);
    
    RowVector4d UnitQuat; UnitQuat<<1.0,0.0,0.0,0.0;
    RowVector4d ZeroQuat; ZeroQuat.setZero();
    
    MatrixXd qzu(z.rows(),4); qzu.setZero(); qzu.block(0,1,z.rows(),3)=zunit;
    MatrixXd qwu(w.rows(),4); qwu.setZero(); qwu.block(0,1,w.rows(),3)=wunit;
    
    MatrixXd qz(z.rows(),4); qz.setZero(); qz.block(0,1,z.rows(),3)=z;
    MatrixXd qw(w.rows(),4); qw.setZero(); qw.block(0,1,w.rows(),3)=w;

    VectorXd InitSolution(16);
    UnitSphereMoebiusTraitsTraits usmt;
    /*InitSolution=VectorXd::Random(16);
    usmt.InitSolution=InitSolution;
    usmt.Initialize(qzu,qwu,InitSolution);
    CheckTraits<UnitSphereMoebiusTraitsTraits>(usmt, 16);*/

    InitSolution<<UnitQuat.transpose(), ZeroQuat.transpose(), ZeroQuat.transpose(), UnitQuat.transpose();
    usmt.InitSolution=InitSolution;
    usmt.Initialize(qzu,qwu,InitSolution);
    
    QuadConstSolver<UnitSphereMoebiusTraitsTraits> qcs;
    qcs.Initialize(&usmt);
    VectorXd Solution=qcs.Solve(InitSolution, 1.0, 1000);
    RowVector4d ua,ub,uc,ud;
    ua<<Solution.segment(0,4).transpose();
    ub<<Solution.segment(4,4).transpose();
    uc<<Solution.segment(8,4).transpose();
    ud<<Solution.segment(12,4).transpose();
    
    cout<<"uaud-ubuc="<<QMult(QConj(ua), ud)-QMult(QConj(ub),uc)<<endl;
    
    //testing center
    if (b.norm()>d.norm()){
        RowVector4d ta=uc,tb=ud,tc=-ua,td=-ub;
        ua=ta; ub=tb; uc=tc; ud=td;
    }
    
    
    double rz=Sz(0);
    RowVector4d cz; cz<<0.0, Sz.tail(3);
    
    double rw=Sw(0);
    RowVector4d cw; cw<<0.0, Sw.tail(3);
    

    //composing the transformations
    ComposeMobius( ua,ub,uc,ud,UnitQuat/rz, -cz/rz, ZeroQuat, UnitQuat, a, b, c, d);
    ComposeMobius(rw*UnitQuat, cw, ZeroQuat, UnitQuat, a, b, c, d, a, b, c, d);
 
    
    //double det=sqrt(abs((QMult(QConj(a), d)-QMult(QConj(b),c))(0)));
    double det=sqrt(abs((QMult(a, QConj(d))+QMult(b,QConj(c)))(0)));
    a/=det;
    b/=det;
    c/=det;
    d/=det;
    cout<<"ad-bc="<<QMult(QConj(a), d)-QMult(QConj(b),c)<<endl;
    cout<<"ad+bc="<<QMult(a, QConj(d))+QMult(b,QConj(c))<<endl;
    
    cout<<"Testing Final Mobius Coeffs: "<<endl;
    MatrixXd Transz=qz;
    for (int i=0;i<z.rows();i++){
        Transz.row(i)=QMult(QMult(a,qz.row(i))+b,QInv(QMult(c,qz.row(i))+d));
        cout<<"(aq+b)*inv(cq+d)-w: "<<Transz.row(i)-qw.row(i)<<endl;
    }
    
    
    /*for (int i=1;i<z.rows();i++){
        RowVector4d Xi=QInv(QMult(c,qz.row(0))+d);
        RowVector4d Xj=QInv(QMult(c,qz.row(i))+d);
        RowVector4d TransVecq=QMult(QMult(QConj(Xj), qz.row(i)-qz.row(0)), Xi);
        cout<<"Xj*qij*Xi:"<<TransVecq<<endl;
        cout<<"Actual vectors: "<<qw.row(i)-qw.row(0)<<endl;
    }
    cout<<"stop here"<<endl;*/
}

inline int sign(double x){
    if (x==0.0) return 0;
    return (x>0.0 ? 1 : -1);
}

void MoebiusDeformation3D::SetupInterpolation(bool isExactMC, bool isExactIAP, bool isFromDeformed)
{
    RowVector4d UnitQuat; UnitQuat<<1.0,0.0,0.0,0.0;
    RowVector4d ZeroQuat; ZeroQuat.setZero();
    
    MatrixXd Mattorhs(OrigVq.rows(),4); Mattorhs.setZero();
    Mattorhs.row(0)=OrigVq.row(0);
    RhsAdd=-d0*Mattorhs;
    
    InterpTraits.Initialize(OrigVq,  D, F, FaceCornerPairs,isExactMC);
    
    
    //checking traits
    /*InterpTraits.InitSolution=VectorXd::Random(4*2*F.rows());
    InterpTraits.PresFactor=10.0;
    InterpTraits.PresJumps=VectorXd::Random(4*FaceCornerPairs.rows());
    CheckTraits<PrescribeEdgeJumps3D>(InterpTraits, 4*2*F.rows());*/
    
    InterpSolver.Initialize(&InterpTraits);
    

    cout<<"Computing the global Mobius transformations between points"<<endl;
    cout<<"=========================================================="<<endl;
    
    //RowVector4d OrigMean=OrigVq.colwise().mean();
    //RowVector4d DeformMean=DeformVq.colwise().mean();
    //choosing the least moving four points
    VectorXd Distances(OrigVq.rows());
    for (int i=0;i<OrigVq.rows();i++)
        Distances(i)=(DeformVq.row(i)-OrigVq.row(i)-(DeformVq.row(0)-OrigVq.row(0))).norm();
    
    VectorXi SortedIndices;
    VectorXd SortedDistances;
    
    igl::sort(Distances, 1,true, SortedDistances, SortedIndices);
    
    DeformGlobalVertices.head(3)=SortedIndices.head(3);
    //finding the nearest nonplanar point
    for (int i=3;i<SortedIndices.size();i++){
        RowVector4d Fourth=OrigVq.row(SortedIndices(i));
        double SignOrig=ComputePlanarity(OrigV.row(DeformGlobalVertices(0)),OrigV.row(DeformGlobalVertices(1)),OrigV.row(DeformGlobalVertices(2)),OrigV.row(SortedIndices(i)));
        double SignDeform=ComputePlanarity(DeformV.row(DeformGlobalVertices(0)).tail(3),DeformV.row(DeformGlobalVertices(1)).tail(3),DeformV.row(DeformGlobalVertices(2)).tail(3),DeformV.row(SortedIndices(i)).tail(3));
        if ((abs(SignOrig)<10e-3)||(sign(SignOrig)!=sign(SignDeform)))
            continue;
        
        DeformGlobalVertices(3)=SortedIndices(i);
        break;
    }
    
    
    MatrixXd z(4,3),w(4,3),cw(4,3);
    for (int i=0;i<4;i++){
        z.row(i)=OrigVq.row(DeformGlobalVertices(i)).tail(3);
        w.row(i)=DeformVq.row(DeformGlobalVertices(i)).tail(3)-OrigVq.row(0).tail(3);
    }
    
    RowVector4d a,b,c,d;

    GetMobiusCoeffs(a,b,c,d, z, w);
    
    //signing properly
    if (d(0)<0.0){
        a=-a; b=-b;c=-c; d=-d;
    }
    
    DeformGlobalMoebius.resize(2,4); DeformGlobalMoebius<<c,d;  //TODO: the inside out testing
    //cout<<"det of t=0->t=1 transformation: "<<QMult(QConj(a),d)-QMult(QConj(b),c)<<endl;
    
    //testing entire mesh
    for (int i=1;i<OrigVq.rows();i++){
        RowVector4d Xi=QInv(QMult(c,OrigVq.row(0))+d);
        RowVector4d Xj=QInv(QMult(c,OrigVq.row(i))+d);
        RowVector4d TransVecq=QMult(QMult(QConj(Xj), OrigVq.row(i)-OrigVq.row(0)), Xi);
        cout<<"Xj*qij*Xi:"<<TransVecq<<endl;
        cout<<"Actual vectors: "<<DeformVq.row(i)-DeformVq.row(0)<<endl;
    }
    
    
    cout<<"Retrieving the Edge Compatibility Coefficients"<<endl;
    cout<<"=============================================="<<endl;
    
    if (isFromDeformed){
        
        //Updating edge coefficient
        EdgeCompCoeffs.resize(FaceCornerPairs.rows(),4);
        MatrixXd Coherence(InnerEdges.size(),4);
        MatrixXd EdgeCoherence(InnerEdges.size(),4);
        for (int i=0;i<InnerEdges.rows();i++){
            int f1=E2F(InnerEdges(i),0);
            int f2=E2F(InnerEdges(i),1);
            int f1i=(E2Fi(InnerEdges(i),0)+1)%F(f1,0);
            int f1j=E2Fi(InnerEdges(i),0);
            int f2j=(E2Fi(InnerEdges(i),1)+1)%F(f2,0);
            int f2i=E2Fi(InnerEdges(i),1);
            RowVector4d zij=OrigVq.row(E2V(InnerEdges(i),1))-OrigVq.row(E2V(InnerEdges(i),0));
            
            RowVector4d G1i=QInv(DeformX.row(CornerOffset(f1)+f1i));
            RowVector4d G2i=QInv(DeformX.row(CornerOffset(f2)+f2i));
            RowVector4d G1j=QInv(DeformX.row(CornerOffset(f1)+f1j));
            RowVector4d G2j=QInv(DeformX.row(CornerOffset(f2)+f2j));
            
            EdgeCoherence.row(i)=QMult(QConj(QInv(G1i)),QMult(zij, QInv(G1j)))-QMult(QConj(QInv(G2i)),QMult(zij, QInv(G2j)));
            
            EdgeCompCoeffs.row(i)=QMult(QInv(G2i),G1i);
            Coherence.row(i)=QMult(QInv(zij),QMult(QMult(QInv(G1j),G2j),zij))-QConj(EdgeCompCoeffs.row(i));
            
        }
        
        cout<<"Comp. Coefficient Coherence:"<<Coherence.lpNorm<Infinity>()<<endl;
        cout<<"Edge Coherence: "<<EdgeCoherence.lpNorm<Infinity>()<<endl;
        return;
    }
    
}


void MoebiusDeformation3D::Interpolate(double t, int NumIterations)
{
    
    RowVector4d UnitQuat; UnitQuat<<1.0,0.0,0.0,0.0;
    RowVector4d ZeroQuat; ZeroQuat.setZero();
    
    //creating gt->c,d matrix for initial solution
    
    cout<<"Linear Initial Interpolation solution"<<endl;
    cout<<"+++++++++++++++++++++++++++++++++++++"<<endl;
    
    MatrixXd gt(InnerEdges.size(),4);
    for (int i=0;i<InnerEdges.size();i++)
        gt.row(i)=QExp(QLog(EdgeCompCoeffs.row(i))*t);
    
    //cout<<"gt:"<<gt<<endl;
    
    //vector<Triplet<RowVector4d> > cdtoGTris;
    VectorXi cdtoGRowIndices(8*FaceCornerPairs.rows()), cdtoGColIndices(8*FaceCornerPairs.rows());
    MatrixXd cdtoGValues(8*FaceCornerPairs.rows(),4);
    //creating matrix with current PresJump Values
    
    //vector<Triplet<RowVector4d> > CompMatTris;
    VectorXi CompMatRowIndices(4*FaceCornerPairs.rows()), CompMatColIndices(4*FaceCornerPairs.rows());
    MatrixXd CompMatValues(4*FaceCornerPairs.rows(),4);
    for (int i=0;i<FaceCornerPairs.rows();i++){
        
        RowVector4d zi=OrigVq.row(FaceCornerPairs(i,2));
        RowVector4d zj=OrigVq.row(FaceCornerPairs(i,3));
        
        //G variables order: G1i, G1j, G2i, G2j
        
        //for conj(G)
        //G2i*g=G1i   (actually conj(g)*conj(G2i)-conj(G1i) )
        //G1j*e*conj(g)*inv(e)=G2j  (actually conj(inv(e))*g*conj(e)*conj(G2j)-conj(G1j) )
        CompMatRowIndices(4*i)=2*i; CompMatColIndices(4*i)=4*i+2; CompMatValues.row(4*i)= QConj(gt.row(i));
        CompMatRowIndices(4*i+1)=2*i; CompMatColIndices(4*i+1)=4*i; CompMatValues.row(4*i+1)= -UnitQuat;
        CompMatRowIndices(4*i+2)=2*i+1; CompMatColIndices(4*i+2)=4*i+1; CompMatValues.row(4*i+2)= QMult(QMult(QInv(QConj(zj-zi)),gt.row(i)), QConj(zj-zi));
        CompMatRowIndices(4*i+3)=2*i+1; CompMatColIndices(4*i+3)=4*i+3; CompMatValues.row(4*i+3)= -UnitQuat;
        
        
        //derivative conj(G1i), conj(G1j) to conj(c1),conj(d1)
        cdtoGRowIndices(8*i)=4*i; cdtoGColIndices(8*i)=2*FaceCornerPairs(i,0); cdtoGValues.row(8*i)=QConj(zi);
        cdtoGRowIndices(8*i+1)=4*i; cdtoGColIndices(8*i+1)=2*FaceCornerPairs(i,0)+1; cdtoGValues.row(8*i+1)=UnitQuat;
        cdtoGRowIndices(8*i+2)=4*i+1; cdtoGColIndices(8*i+2)=2*FaceCornerPairs(i,0); cdtoGValues.row(8*i+2)=QConj(zj);
        cdtoGRowIndices(8*i+3)=4*i+1; cdtoGColIndices(8*i+3)=2*FaceCornerPairs(i,0)+1; cdtoGValues.row(8*i+3)=UnitQuat;
        
        //derivative conj(Gi2), conj(G2j) to conj(c2),conj(d2)
        cdtoGRowIndices(8*i+4)=4*i+2; cdtoGColIndices(8*i+4)=2*FaceCornerPairs(i,1); cdtoGValues.row(8*i+4)=QConj(zi);
        cdtoGRowIndices(8*i+5)=4*i+2; cdtoGColIndices(8*i+5)=2*FaceCornerPairs(i,1)+1; cdtoGValues.row(8*i+5)=UnitQuat;
        cdtoGRowIndices(8*i+6)=4*i+3; cdtoGColIndices(8*i+6)=2*FaceCornerPairs(i,1); cdtoGValues.row(8*i+6)=QConj(zj);
        cdtoGRowIndices(8*i+7)=4*i+3; cdtoGColIndices(8*i+7)=2*FaceCornerPairs(i,1)+1; cdtoGValues.row(8*i+7)=UnitQuat;
    }
    
    SparseMatrix<double> cdtoG;
    quat2SparseMatrix(cdtoGRowIndices, cdtoGColIndices, cdtoGValues,cdtoG,  4*FaceCornerPairs.rows(), 2*F.rows());

    SparseMatrix<double> CompMat;
    quat2SparseMatrix(CompMatRowIndices, CompMatColIndices, CompMatValues, CompMat, 2*FaceCornerPairs.rows(), 4*FaceCornerPairs.rows());
    
    CompMat=CompMat*cdtoG*quatConjMat((int)(2*F.rows()));
    
    VectorXd torhs(8*F.rows()); torhs.setZero();
    torhs.head(8)<<ZeroQuat.transpose(), UnitQuat.transpose();
    VectorXd Rhs=-CompMat*torhs;
    
    VectorXi LocalVarIndices(8*F.rows()-8);
    for (int i=0;i<8*(F.rows()-1);i++)
        LocalVarIndices(i)=i+8;
    
    SparseMatrix<double> CompMatVar;
    igl::slice(CompMat, LocalVarIndices, 2, CompMatVar);
    
    
    VectorXd InitialMobCoeffs=hedra::optimization::EigenSingleSolveWrapper<Eigen::SparseQR<Eigen::SparseMatrix<double, Eigen::ColMajor>, Eigen::COLAMDOrdering<int> > >(CompMatVar.adjoint()*CompMatVar,CompMatVar.adjoint()*Rhs);
    
    cout<<"Initial solution g^t->(c,d) Error: "<<(CompMatVar*InitialMobCoeffs-Rhs).lpNorm<Infinity>()<<endl;
    
    //cout<<"InitialMobCoeffs: "<<InitialMobCoeffs<<endl;
    
    //testing imagiary of coeffs
    VectorXd cdImagValue(F.rows());
    for (int i=0;i<F.rows()-1;i++){
        RowVector4d c=InitialMobCoeffs.segment(8*i,4);
        RowVector4d d=InitialMobCoeffs.segment(8*i+4,4);
        cdImagValue(i)=QMult(c,QConj(d))(0);
    }
    
    cout<<"Initial c,d imaginarity: "<<cdImagValue.lpNorm<Infinity>()<<endl;


    VectorXd gtvec(4*gt.rows());
    for (int i=0;i<gt.rows();i++)
        gtvec.segment(4*i,4)=gt.row(i).transpose();

    
    //VectorXd InitJumpSolution(4*(InnerEdges.size()+2*F.rows()));
    //InitJumpSolution<<gtvec, ZeroQuat.transpose(), UnitQuat.transpose(), InitialMobCoeffs;
    VectorXd InitJumpSolution(4*2*F.rows());
    InitJumpSolution<<ZeroQuat.transpose(), UnitQuat.transpose(), InitialMobCoeffs;
    
    cout<<"Optimizing for Exact Interpolation solution"<<endl;
    cout<<"+++++++++++++++++++++++++++++++++++++++++++"<<endl;
    
    
    InterpTraits.InitSolution=InitJumpSolution;
    InterpTraits.PresFactor=10.0;
    InterpTraits.PresJumps.resize(4*EdgeCompCoeffs.rows());
    for (int i=0;i<EdgeCompCoeffs.rows();i++)
        InterpTraits.PresJumps.segment(4*i,4)=gt.row(i).transpose();
    VectorXd IntegratedSolution=InterpSolver.Solve(InitJumpSolution, 1.0, NumIterations);
    
    //reconstructing the current interpolated mesh
    MatrixXd FaceEdgeVectors(E2F.rows(),4);
    FaceEdgeVectors.setZero();
    for (int i=0;i<E2F.rows();i++){
        if (E2F(i,0)!=-1){
            RowVector4d c1=IntegratedSolution.segment(4*(2*E2F(i,0)),4).transpose();
            RowVector4d d1=IntegratedSolution.segment(4*(2*E2F(i,0)+1),4).transpose();
            //cout<<"Imaginarity of c,d: "<<QMult(c1,QConj(d1))<<endl;
            RowVector4d X1=QInv(QMult(c1,OrigVq.row(E2V(i,0)))+d1);
            RowVector4d X2=QInv(QMult(c1,OrigVq.row(E2V(i,1)))+d1);
            FaceEdgeVectors.row(i)=QMult(QConj(X2),QMult(OrigVq.row(E2V(i,1))-OrigVq.row(E2V(i,0)),X1));
        }
        if (E2F(i,1)!=-1){
            RowVector4d c2=IntegratedSolution.segment(4*(2*E2F(i,1)),4).transpose();
            RowVector4d d2=IntegratedSolution.segment(4*(2*E2F(i,1)+1),4).transpose();
            //cout<<"Imaginarity of c,d: "<<QMult(c2,QConj(d2))<<endl;
            RowVector4d X1=QInv(QMult(c2,OrigVq.row(E2V(i,0)))+d2);
            RowVector4d X2=QInv(QMult(c2,OrigVq.row(E2V(i,1)))+d2);
            FaceEdgeVectors.row(i)+=QMult(QConj(X2),QMult(OrigVq.row(E2V(i,1))-OrigVq.row(E2V(i,0)),X1));

        }
        
        if ((E2F(i,0)!=-1)&&(E2F(i,1)!=-1))
            FaceEdgeVectors.row(i)/=2.0;
    }
    
  
    
    //cout<<"FaceEdgeVectors: "<<FaceEdgeVectors<<endl;
    
    cout<<"Reconstructing an Integrable Mesh"<<endl;
    cout<<"+++++++++++++++++++++++++++++++++"<<endl;
    
    MatrixXd d0Rhs=d0NoFirst.adjoint()*(RhsAdd+FaceEdgeVectors);
    MatrixXd CurrInterpVqVar;
    d0Solver.solve(d0Rhs,CurrInterpVqVar);
    MatrixXd CurrInterpVq(OrigVq.rows(),4);
    CurrInterpVq.row(0)=OrigVq.row(0);
    CurrInterpVq.block(1,0,OrigVq.rows()-1,4)=CurrInterpVqVar;
    
  
    cout<<"Reconstruction Error: "<<(d0*CurrInterpVq-FaceEdgeVectors).lpNorm<Infinity>()<<endl;
    
    //cout<<"Reconstructed FaceEdgeVectors: "<<FaceEdgeVectors<<endl;
    
    cout<<"Final Global Mobius Transformation"<<endl;
    cout<<"++++++++++++++++++++++++++++++++++"<<endl;
    
    cout<<"Global Mobius Transformation: "<<DeformGlobalMoebius<<endl;
    

    cout<<"Transforming reconstructed -> Original"<<endl;
    
    MatrixXd z(4,3),cw(4,3);
    cout<<"DeformGlobalVertices:"<<DeformGlobalVertices<<endl;
    for (int i=0;i<4;i++){
        z.row(i)=OrigVq.row(DeformGlobalVertices(i)).tail(3);
        cw.row(i)=CurrInterpVq.row(DeformGlobalVertices(i)).tail(3);
    }
    
    //transforming reconstruction -> t=0
    RowVector4d a2zero, b2zero, c2zero, d2zero;
    GetMobiusCoeffs(a2zero, b2zero, c2zero, d2zero, cw, z);
    //cout<<"det of 2zero transformation: "<<QMult(QConj(a2zero),d2zero)-QMult(QConj(b2zero),c2zero)<<endl;
    
    /*for (int i=0;i<d0.rows();i++){
        RowVector4d Xi2zero=QInv(QMult(c2zero, CurrInterpVq.row(E2V(i,0)))+d2zero);
        RowVector4d Xj2zero=QInv(QMult(c2zero, CurrInterpVq.row(E2V(i,1)))+d2zero);
        
        FaceEdgeVectors.row(i)=QMult(QMult(QConj(Xj2zero), FaceEdgeVectors.row(i)), Xi2zero);
        
    }*/
    for (int i=0;i<OrigVq.rows();i++)
        CurrInterpVq.row(i)=QMult(QMult(a2zero, CurrInterpVq.row(i))+b2zero,QInv(QMult(c2zero, CurrInterpVq.row(i))+d2zero));
    
    //cout<<"FaceEdgeVectors: "<<FaceEdgeVectors<<endl;
    cout<<"Imaginarity of 2zero c,d: "<<QMult(c2zero,QConj(d2zero))<<endl;
    cout<<"c2zero, d2zero: "<<c2zero<<endl<<d2zero<<endl;
    
    /*d0Rhs=d0NoFirst.adjoint()*(RhsAdd+FaceEdgeVectors);
    d0Solver.solve(d0Rhs,CurrInterpVqVar);
    CurrInterpVq.row(0)=OrigVq.row(0);
    CurrInterpVq.block(1,0,OrigVq.rows()-1,4)=CurrInterpVqVar;
    
    cout<<"Reconstruction Error: "<<(d0*CurrInterpVq-FaceEdgeVectors).lpNorm<Infinity>()<<endl;*/

    //interpolating t=0->1 c,d into t=0,t
    cout<<"Transforming original -> interpolated"<<endl;
    RowVector4d ct, dt;
    
    dt=QExp(QLog(DeformGlobalMoebius.row(1))*t);
    ct=t*QMult(DeformGlobalMoebius.row(0),QConj(QExp(QLog(DeformGlobalMoebius.row(1))*(1-t))));
    
    
    FaceEdgeVectors=d0*CurrInterpVq;
    //interpolating edge vectors and then re-integrating NEED TO MAKE THIS BETTER
    /*cout<<"CurrInterpVq:"<<CurrInterpVq<<endl;
    cout<<"FaceEdgeVectors: "<<FaceEdgeVectors<<endl;
    //testing entire mesh
    cout<<"LOCAL TESTING::::::"<<endl;
    for (int i=1;i<OrigVq.rows();i++){
        RowVector4d Xi=QInv(QMult(ct,CurrInterpVq.row(0))+dt);
        RowVector4d Xj=QInv(QMult(ct,CurrInterpVq.row(i))+dt);
        RowVector4d TransVecq=QMult(QMult(QConj(Xj), CurrInterpVq.row(i)-CurrInterpVq.row(0)), Xi);
        cout<<"Xj*qij*Xi:"<<TransVecq<<endl;
        cout<<"Actual vectors: "<<DeformVq.row(i)-DeformVq.row(0)<<endl;
    }*/
    
    
    for (int i=0;i<d0.rows();i++){
        
        RowVector4d Xit=QInv(QMult(ct, CurrInterpVq.row(E2V(i,0)))+dt);
        RowVector4d Xjt=QInv(QMult(ct, CurrInterpVq.row(E2V(i,1)))+dt);
        
        FaceEdgeVectors.row(i)=QMult(QMult(QConj(Xjt), FaceEdgeVectors.row(i)), Xit);

    }
    
    
    //cout<<"FaceEdgeVectors: "<<FaceEdgeVectors<<endl;
    cout<<"Imaginarity of interpolated ct,dt: "<<QMult(ct,QConj(dt))<<endl;
    cout<<"ct, dt: "<<ct<<endl<<dt<<endl;
    cout<<"DeformGlobalMoebius: "<<DeformGlobalMoebius<<endl;
    
    d0Rhs=d0NoFirst.adjoint()*(RhsAdd+FaceEdgeVectors);
    d0Solver.solve(d0Rhs,CurrInterpVqVar);
    CurrInterpVq.row(0)=OrigVq.row(0);
    CurrInterpVq.block(1,0,OrigVq.rows()-1,4)=CurrInterpVqVar;
    
    for (int i=0;i<OrigVq.rows();i++)
        CurrInterpVq.row(i)+=t*(DeformVq.row(0)-OrigVq.row(0));
    
    cout<<"Reconstruction Error: "<<(d0*CurrInterpVq-FaceEdgeVectors).lpNorm<Infinity>()<<endl;

    InterpVq=CurrInterpVq;
    
    Quat2Coords(InterpVq, InterpV);
    ComputeCR(InterpVq, D, F, QuadVertexIndices, InterpECR, InterpFCR);

    MatrixXd Centers(F.rows(),3); Centers.setZero();
    for (int i=0;i<F.rows();i++)
        for (int j=0;j<D(i);j++)
            Centers.row(i)+=InterpV.row(F(i,j))/(double)D(i);
    
    EdgeInterpV<<InterpV, Centers;
    
    InterpConvErrors=InterpSolver.ConvErrors;
}



