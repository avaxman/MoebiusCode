//
//  ExtraFunctions.cpp
//  testigl
//
//  Created by Amir Vaxman on 18/12/14.
//  Copyright (c) 2014 Amir Vaxman. All rights reserved.
//

#include "ExtraFunctions.h"
#include <hedra/EigenSolverWrapper.h>


using namespace Eigen;
using namespace std;

typedef std::complex<double> Complex;

void GetComplexMobiusCoeffs(Complex& a, Complex& b, Complex& c, Complex& d, const Vector3cd& z, const Vector3cd& w)
{
    Vector3cd zw=z.cwiseProduct(w);
    Vector3cd ones=Vector3cd::Ones();
    Matrix3cd amat,bmat,cmat,dmat;
    amat.col(0)=zw; amat.col(1)=w; amat.col(2)=ones;
    bmat.col(0)=zw; bmat.col(1)=z; bmat.col(2)=w;
    cmat.col(0)=z; cmat.col(1)=w; cmat.col(2)=ones;
    dmat.col(0)=zw; dmat.col(1)=z; dmat.col(2)=ones;
    
    a=amat.determinant();
    b=bmat.determinant();
    c=cmat.determinant();
    d=dmat.determinant();
    
    Complex TotalDet=a*d-b*c;
    Complex Sign=(real(a)>0 ? 1.0 : -1.0);
    a/=sqrt(TotalDet)*Sign;
    b/=sqrt(TotalDet)*Sign;
    c/=sqrt(TotalDet)*Sign;
    d/=sqrt(TotalDet)*Sign;
    
}


Vector2cd GetSingleMobius(const VectorXcd& OrigVc, const VectorXcd& X)
{
    
    //estimating target global mobius transformation
    MatrixXcd MobMat(X.rows(),2);
    MobMat.col(0)=OrigVc;
    MobMat.col(1).setConstant(1.0);
    
    VectorXcd GFunc(X.size());
    for (int i=0;i<X.rows();i++)
        GFunc(i)=Complex(1.0)/X(i);
    
    Matrix2cd A=(MobMat.adjoint()*MobMat);
    Vector2cd b=(MobMat.adjoint()*GFunc);
    Vector2cd Result=A.colPivHouseholderQr().solve(b);
    
    cout<<"Total Target Mobius Error: "<<(MobMat*Result-GFunc).lpNorm<Infinity>()<<endl;
    return(Result);
    
}


VectorXcd SolveComplexSytem(const SparseMatrix<Complex>& A, const VectorXcd& b)
{
    //creating real matrices
    SparseMatrix<double> rA(A.rows()*2, A.cols()*2);
    
    vector<Triplet<double> > RealTris(A.nonZeros()*4);
    
    int Counter=0;
    for (int k=0; k<A.outerSize(); k++){
        for (SparseMatrix<Complex>::InnerIterator it(A,k); it; ++it){
            RealTris[Counter++]=Triplet<double>(it.row(),			it.col(),			it.value().real());
            RealTris[Counter++]=Triplet<double>(it.row(),			it.col()+A.cols(),	-it.value().imag());
            RealTris[Counter++]=Triplet<double>(it.row()+A.rows(),	it.col(),			it.value().imag());
            RealTris[Counter++]=Triplet<double>(it.row()+A.rows(),	it.col()+A.cols(),	it.value().real());
        }
    }
    rA.setFromTriplets(RealTris.begin(), RealTris.end());
    VectorXd rb(b.rows()*2);
    
    rb<<b.real(), b.imag();
    
    MatrixXd RawSolution=hedra::optimization::EigenSingleSolveWrapper<Eigen::SimplicialLDLT<Eigen::SparseMatrix<double> > >(rA, rb);
    
    int ComplexSize=A.cols();
    
    VectorXcd Solution(ComplexSize);
    for (int i=0;i<ComplexSize;i++)
        Solution[i]=Complex(RawSolution(i), RawSolution(ComplexSize+i));
    
    return Solution;
    
}


RowVector4d Rot2Quat(Matrix3d& R)
{
    
    double Rxx = R(0,0); double Rxy = R(0,1); double Rxz = R(0,2);
    double Ryx = R(1,0); double Ryy = R(1,1); double Ryz = R(1,2);
    double Rzx = R(2,0); double Rzy = R(2,1); double Rzz = R(2,2);
    
    double w = sqrt( R.trace() + 1 ) / 2;
    double x = sqrt( 1 + Rxx - Ryy - Rzz ) / 2;
    double y = sqrt( 1 + Ryy - Rxx - Rzz ) / 2;
    double z = sqrt( 1 + Rzz - Ryy - Rxx ) / 2;
    
    VectorXd TestVec(4); TestVec<<w,x,y,z;
    int MinIndex;
    TestVec.maxCoeff(&MinIndex);
    
    if( MinIndex == 0 ){
        x = ( Rzy - Ryz ) / (4*w);
        y = ( Rxz - Rzx ) / (4*w);
        z = ( Ryx - Rxy ) / (4*w);
    }
    
    if( MinIndex == 1 ){
        w = ( Rzy - Ryz ) / (4*x);
        y = ( Rxy + Ryx ) / (4*x);
        z = ( Rzx + Rxz ) / (4*x);
    }
    
    if( MinIndex == 2 ){
        w = ( Rxz - Rzx ) / (4*y);
        x = ( Rxy + Ryx ) / (4*y);
        z = ( Ryz + Rzy ) / (4*y);
    }
    
    if( MinIndex == 3 ){
        w = ( Ryx - Rxy ) / (4*z);
        x = ( Rzx + Rxz ) / (4*z);
        y = ( Ryz + Rzy ) / (4*z);
    }
    RowVector4d Result; Result<<w,x,y,z;
    
    return Result;
}

void TriangulateGeneralMesh(const MatrixXi& D, const MatrixXi& F, MatrixXi& tF, VectorXi& FromFace)
{
    vector<Vector3i> NewTriangles;
    vector<int> RawFromFace;
    
    
    //cout<<"F :"<<F<<endl;
    
    for (int i=0;i<D.rows();i++){
        //triangulating the face greedily
        for (int CurrIndex=1;CurrIndex<D(i)-1;CurrIndex++){
            Vector3i NewFace;
            NewFace<<F(i,0),F(i,CurrIndex),F(i,CurrIndex+1);
            RawFromFace.push_back(i);
            NewTriangles.push_back(NewFace);
        }
    }
    
    tF.resize(NewTriangles.size(),3);
    FromFace.resize(RawFromFace.size());
    for (int i=0;i<NewTriangles.size();i++){
        tF.row(i)=NewTriangles[i];
        FromFace(i)=RawFromFace[i];
    }
    
    //cout<<"tF :"<<tF<<endl;
    //cout<<"FromFace :"<<FromFace<<endl;
    
}

MatrixXd GetCenters(const MatrixXd& V, const MatrixXi& D, const MatrixXi& F)
{
    MatrixXd Centers(D.rows(),3); Centers.setZero();
    for (int i=0;i<D.rows();i++)
        for (int j=0;j<D(i);j++)
            Centers.row(i)+=V.row(F(i,j))/(double)D(i,0);
    
    return Centers;
    
}






