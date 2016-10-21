//
//  QuaternionOps.h
//  testigl
//
//  Created by Amir Vaxman on 22/08/14.
//  Copyright (c) 2014 Amir Vaxman. All rights reserved.
//

#ifndef __testigl__QuaternionOps__
#define __testigl__QuaternionOps__

#include <iostream>
#include <Eigen/Dense>
#include <Eigen/Sparse>

using namespace Eigen;
using namespace std;

inline RowVector4d QConj1(const RowVector4d& q)
{
    RowVector4d newq;
    newq<<q(0), -q.tail(3);
    return newq;
}

inline MatrixXd QConjn(const MatrixXd& q)
{
    MatrixXd newq(q.rows(),4);
    newq<<q.col(0), -q.block(0,1,q.rows(), q.cols()-1);
    return newq;
}

inline MatrixXd QMultn(const MatrixXd& q1, const MatrixXd& q2)
{
    MatrixXd newq(q1.rows(),4);
    VectorXd r1=q1.col(0);
    VectorXd r2=q2.col(0);
    MatrixXd v1=q1.block(0,1,q1.rows(), 3);
    MatrixXd v2=q2.block(0,1,q2.rows(), 3);
    MatrixXd r1mat; r1mat.resize(r1.rows(),3); r1mat<<r1, r1, r1;
    MatrixXd r2mat; r2mat.resize(r2.rows(),3); r2mat<<r2, r2, r2;
    newq.col(0)=r1.cwiseProduct(r2)-(v1.cwiseProduct(v2)).rowwise().sum();
    MatrixXd v1cv2; v1cv2.resize(v1.rows(),3);
    for (int i=0;i<v1.rows();i++){
        Vector3d vv1=v1.row(i);
        Vector3d vv2=v2.row(i);
        v1cv2.row(i)=vv1.cross(vv2);
        
    }
    newq.block(0,1,newq.rows(),newq.cols()-1)=v1.cwiseProduct(r2mat)+v2.cwiseProduct(r1mat)+v1cv2;
    return newq;
}

inline RowVector4d QMult1(const RowVector4d& q1, const RowVector4d& q2)
{
    RowVector4d newq;
    double r1=q1(0);
    double r2=q2(0);
    RowVector3d v1=q1.tail(3);
    RowVector3d v2=q2.tail(3);
    newq<<r1*r2-v1.dot(v2), r1*v2+r2*v1+v1.cross(v2);
    return newq;
}

inline RowVector4d QInv1(const RowVector4d& q)
{
    return(QConj1(q)/q.squaredNorm());
}

inline MatrixXd QInvn(const MatrixXd& q)
{
    MatrixXd newq;
    VectorXd qabs=q.rowwise().squaredNorm();
    MatrixXd qabsmat(qabs.rows(),4); qabsmat<<qabs,qabs,qabs,qabs;
    newq=QConjn(q).cwiseQuotient(qabsmat);
    return newq;
}

inline MatrixXd QLog(const MatrixXd& q)
{
    VectorXd nq=q.rowwise().norm();
    VectorXd nv=q.block(0,1,q.rows(),q.cols()-1).rowwise().norm();
    VectorXd acosqnq=acos((q.col(0).cwiseQuotient(nq)).array()).matrix().cwiseQuotient(nv);
    MatrixXd acosmat(acosqnq.rows(),3); acosmat<<acosqnq, acosqnq, acosqnq;
    MatrixXd logq(q.rows(),q.cols());
    logq<<log(nq.array()), q.block(0,1,q.rows(),q.cols()-1).cwiseProduct(acosmat);
    for (int i=0;i<logq.rows();i++)
        if (nv(i)<10e-6)
            logq.row(i)<<log(nq(i)),0.0,0.0,0.0;
    
    //cout<<"logq: "<<logq<<endl;
    return logq;
}

inline MatrixXd QExp(const MatrixXd& q)
{
    VectorXd nv=q.block(0,1,q.rows(),q.cols()-1).rowwise().norm();
    VectorXd exp1=exp(q.col(0).array());
    MatrixXd exp1mat(exp1.rows(),4); exp1mat<<exp1,exp1,exp1,exp1;
    MatrixXd expq(q.rows(),q.cols());
    VectorXd sinnv=(sin(nv.array()).matrix()).cwiseQuotient(nv);
    MatrixXd sinnvmat(sinnv.rows(),3); sinnvmat<<sinnv,sinnv,sinnv;
    expq<<cos(nv.array()), q.block(0,1,q.rows(),q.cols()-1).cwiseProduct(sinnvmat);
    expq=expq.cwiseProduct(exp1mat);
    for (int i=0;i<expq.rows();i++)
        if (nv(i)<10e-6)
            expq.row(i)<<exp1(i),0.0,0.0,0.0;
    
    return expq;
}



inline void Quat2SparseMatrix(const VectorXi& RowIndices, const VectorXi& ColIndices, const MatrixXd& Values, SparseMatrix<double>& Mat, const int m, const int n)
{
    /*[  rq, -vqx, -vqy, -vqz]
     [ vqx,   rq, -vqz,  vqy]
     [ vqy,  vqz,   rq, -vqx]
     [ vqz, -vqy,  vqx,   rq]*/
    Mat.resize(4*m,4*n);
    vector<Triplet<double> > RealTris(16*RowIndices.size());
    for (int i=0;i<RowIndices.size();i++){
        
        //r=rq*rp-<vq,vp>
        RealTris[16*i]=Triplet<double>(4*RowIndices(i),      4*ColIndices(i),     Values(i,0));
        RealTris[16*i+1]=Triplet<double>(4*RowIndices(i)+1,    4*ColIndices(i)+1,   Values(i,0));
        RealTris[16*i+2]=Triplet<double>(4*RowIndices(i)+2,    4*ColIndices(i)+2,   Values(i,0));
        RealTris[16*i+3]=Triplet<double>(4*RowIndices(i)+3,    4*ColIndices(i)+3,   Values(i,0));
        
        //v=rq*vp+vp*vq+(vq x vp)
        RealTris[16*i+4]=Triplet<double>(4*RowIndices(i),      4*ColIndices(i)+1,   -Values(i,1));
        RealTris[16*i+5]=Triplet<double>(4*RowIndices(i)+1,    4*ColIndices(i),     Values(i,1));
        RealTris[16*i+6]=Triplet<double>(4*RowIndices(i)+2,    4*ColIndices(i)+3,   -Values(i,1));
        RealTris[16*i+7]=Triplet<double>(4*RowIndices(i)+3,    4*ColIndices(i)+2,   Values(i,1));
        
        RealTris[16*i+8]=Triplet<double>(4*RowIndices(i),      4*ColIndices(i)+2,   -Values(i,2));
        RealTris[16*i+9]=Triplet<double>(4*RowIndices(i)+1,    4*ColIndices(i)+3,   Values(i,2));
        RealTris[16*i+10]=Triplet<double>(4*RowIndices(i)+2,    4*ColIndices(i),     Values(i,2));
        RealTris[16*i+11]=Triplet<double>(4*RowIndices(i)+3,    4*ColIndices(i)+1,   -Values(i,2));
        
        RealTris[16*i+12]=Triplet<double>(4*RowIndices(i),      4*ColIndices(i)+3,   -Values(i,3));
        RealTris[16*i+13]=Triplet<double>(4*RowIndices(i)+1,    4*ColIndices(i)+2,   -Values(i,3));
        RealTris[16*i+14]=Triplet<double>(4*RowIndices(i)+2,    4*ColIndices(i)+1,   Values(i,3));
        RealTris[16*i+15]=Triplet<double>(4*RowIndices(i)+3,    4*ColIndices(i),     Values(i,3));
    }
    
    Mat.setFromTriplets(RealTris.begin(), RealTris.end());
}

inline SparseMatrix<double> QuatConjMat(int m){
    SparseMatrix<double> Mat(4*m,4*m);
    vector<Triplet<double> > Tris(4*m);
    for (int i=0;i<m;i++){
        Tris[4*i]=Triplet<double>(4*i,4*i,1.0);
        Tris[4*i+1]=Triplet<double>(4*i+1,4*i+1,-1.0);
        Tris[4*i+2]=Triplet<double>(4*i+2,4*i+2,-1.0);
        Tris[4*i+3]=Triplet<double>(4*i+3,4*i+3,-1.0);
    }
    
    Mat.setFromTriplets(Tris.begin(), Tris.end());
    return Mat;
    
}




//deriving an expression a*X*b (or a*conj(X)*b) by X.

inline void GetQuatDerivativeIndices(VectorXi& Rows, VectorXi& Cols, const int CurrTriPos, Vector4i TriSkips, const int Row, const int Col)
{
    for (int i=0;i<4;i++)
        for (int j=0;j<4;j++){
            Rows(CurrTriPos+TriSkips(i)+j)=Row+i;
            Cols(CurrTriPos+TriSkips(i)+j)=Col+j;
        }
}


inline void GetQuatDerivativeValues(VectorXd& Values, const int CurrTriPos, Vector4i TriSkips, const RowVector4d& LeftCoeff, const RowVector4d& RightCoeff, const bool isConj, const bool add)
{
    //[  ra, -vax, -vay, -vaz]
    //[ vax,   ra, -vaz,  vay]
    //[ vay,  vaz,   ra, -vax]
    //[ vaz, -vay,  vax,   ra]
    
    Matrix4d a; a<<LeftCoeff(0), -LeftCoeff(1), -LeftCoeff(2), -LeftCoeff(3),
    LeftCoeff(1),  LeftCoeff(0), -LeftCoeff(3), LeftCoeff(2),
    LeftCoeff(2),  LeftCoeff(3), LeftCoeff(0), -LeftCoeff(1),
    LeftCoeff(3), -LeftCoeff(2), LeftCoeff(1), LeftCoeff(0);
    
    //cout<<"a: "<<a<<endl;
    
    //[  rb, -vbx, -vby, -vbz]
    //[ vbx,   rb,  vbz, -vby]
    //[ vby, -vbz,   rb,  vbx]
    //[ vbz,  vby, -vbx,   rb]
    
    Matrix4d b; b<<RightCoeff(0), -RightCoeff(1), -RightCoeff(2), -RightCoeff(3),
    RightCoeff(1),  RightCoeff(0),  RightCoeff(3), -RightCoeff(2),
    RightCoeff(2), -RightCoeff(3),  RightCoeff(0),  RightCoeff(1),
    RightCoeff(3),  RightCoeff(2), -RightCoeff(1),  RightCoeff(0);
    
    //cout<<"b: "<<b<<endl;
    
    if (!isConj){
        
        //not conjugate
        //[ ra*rb - vax*vbx - vay*vby - vaz*vbz, vay*vbz - rb*vax - ra*vbx - vaz*vby, vaz*vbx - rb*vay - vax*vbz - ra*vby, vax*vby - rb*vaz - ra*vbz - vay*vbx]
        //[ ra*vbx + rb*vax + vay*vbz - vaz*vby, ra*rb - vax*vbx + vay*vby + vaz*vbz, ra*vbz - rb*vaz - vax*vby - vay*vbx, rb*vay - ra*vby - vax*vbz - vaz*vbx]
        //[ ra*vby + rb*vay - vax*vbz + vaz*vbx, rb*vaz - ra*vbz - vax*vby - vay*vbx, ra*rb + vax*vbx - vay*vby + vaz*vbz, ra*vbx - rb*vax - vay*vbz - vaz*vby]
        //[ ra*vbz + rb*vaz + vax*vby - vay*vbx, ra*vby - rb*vay - vax*vbz - vaz*vbx, rb*vax - ra*vbx - vay*vbz - vaz*vby, ra*rb + vax*vbx + vay*vby - vaz*vbz]
        
        
        Matrix4d ab=a*b;
        
        //cout<<"ab: "<<ab<<endl;
        
        for (int i=0;i<4;i++)
            for (int j=0;j<4;j++)
                if (add)
                    Values(CurrTriPos+TriSkips(i)+j)+=ab(i,j);
                else
                    Values(CurrTriPos+TriSkips(i)+j)=ab(i,j);
        
        
    } else {
        
        Matrix4d ConjMat; ConjMat<<1.0,0.0,0.0,0.0,
        0.0,-1.0,0.0,0.0,
        0.0,0.0,-1.0,0.0,
        0.0,0.0,0.0,-1.0;
        Matrix4d abc=a*b*ConjMat;
        
        //cout<<"abc:"<<abc<<endl;
        
        for (int i=0;i<4;i++)
            for (int j=0;j<4;j++)
                if (add)
                    Values(CurrTriPos+TriSkips(i)+j)+=abc(i,j);
                else
                    Values(CurrTriPos+TriSkips(i)+j)=abc(i,j);
        
    }
}

inline void Quat2Coords(const MatrixXd& QV, MatrixXd& V)
{
    V.resize(QV.rows(), 3);
    V=QV.block(0,1,QV.rows(),QV.cols()-1);
}


inline void Coords2Quat(const MatrixXd& V, MatrixXd& QV)
{
    QV.resize(V.rows(),4); QV.setZero();
    QV.block(0,1,QV.rows(),QV.cols()-1)=V;
}





#endif /* defined(__testigl__QuaternionOps__) */
