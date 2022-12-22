//
//  PrescribeEdgeJumps.h
//  testigl
//
//  Created by Amir Vaxman on 22/12/14.
//  Copyright (c) 2014 Amir Vaxman. All rights reserved.
//

#ifndef MoebiusCode_PrescribeEdgeJumps_h
#define MoebiusCode_PrescribeEdgeJumps_h

#include <Eigen/Core>
#include <hedra/quaternionic_operations.h>
#include <hedra/quaternionic_derivatives.h>



class PrescribeEdgeJumps3D{
public:
    
    Eigen::MatrixXi D, F;
    
    Eigen::MatrixXd OrigVq;
    Eigen::MatrixXi FaceCornerIndices;  //rows of f1, f2, zi, zj (oriented)
    Eigen::VectorXd PresJumps;
    Eigen::VectorXd InitSolution;
    Eigen::VectorXd Firstcd;
    
    double PresFactor;
    double CloseFactor;
    double ConstTolerance;
    double SolutionSize;

    Eigen::RowVector4d UnitQuat;
    
    bool isExactMC;
    
    //intermediate variables
    Eigen::VectorXd PresVec;
    Eigen::VectorXd CloseVec;
    Eigen::VectorXd CompVec;
    Eigen::VectorXd ImagVec;
    Eigen::VectorXd FirstcdVec;
    Eigen::VectorXd MCVec;
    
    Eigen::VectorXd TotalVec;
    Eigen::VectorXd ConstVec;
    
    int PresTriOffset, PresRowOffset;
    int CloseTriOffset, CloseRowOffset;
    int CompTriOffset, CompRowOffset;
    int ImagTriOffset, ImagRowOffset;
    int FirstcdTriOffset, FirstcdRowOffset;
    int MCTriOffset, MCRowOffset;
    
    Eigen::VectorXi GradRows, GradCols;
    Eigen::VectorXd GradValues;
    
    void Initialize(const Eigen::MatrixXd& inOrigVq,
                    const Eigen::MatrixXi& inD,
                    const Eigen::MatrixXi& inF,
                    Eigen::MatrixXi& inFaceCornerIndices,
                    const bool& inisExactMC){
        
        using namespace Eigen;
        using namespace std;
        F=inF; D=inD;
        FaceCornerIndices=inFaceCornerIndices;
        OrigVq=inOrigVq;
        isExactMC=inisExactMC;
        
        SolutionSize=2*4*F.rows();
        
        UnitQuat<<1.0,0.0,0.0,0.0;
        
        PresVec.resize(4*2*FaceCornerIndices.rows());
        CloseVec.resize(SolutionSize);
        CompVec.resize(4*FaceCornerIndices.rows());
        ImagVec.resize(F.rows());
        FirstcdVec.resize(8);
        Firstcd.resize(8);
        Firstcd<<0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0;
        

        if (isExactMC)
            MCVec.resize(2*FaceCornerIndices.rows());
        else
            MCVec.resize(0);
        
        
        
        TotalVec.resize(PresVec.size()+CloseVec.size()+CompVec.size()+ImagVec.size()+FirstcdVec.size()+MCVec.size());
        ConstVec.resize(CompVec.size()+ImagVec.size()+FirstcdVec.size()+MCVec.size());
        
        ConstTolerance=10e-7;
        
        //***********************Creating Gradient Pattern**************************************/
        
        if (!isExactMC){
            GradRows.resize(128*FaceCornerIndices.rows()+SolutionSize+64*FaceCornerIndices.rows()+8*F.rows()+8);
            GradCols.resize(GradRows.size());
            GradValues.resize(GradRows.size());
        }  //not currently doing exact MC error reproduction
        
        /**************************Prescription Energy*******************************/
        PresTriOffset=0;
        PresRowOffset=0;
        Vector4i c1TriPoses; c1TriPoses<<0,16,32,48;
        Vector4i d1TriPoses; d1TriPoses<<4,20,36,52;
        Vector4i c2TriPoses; c2TriPoses<<8,24,40,56;
        Vector4i d2TriPoses; d2TriPoses<<12,28,44,60;
        int PresTriCounter=PresTriOffset;
        for (int i=0;i<FaceCornerIndices.rows();i++){
            
            int Colc1=4*(2*FaceCornerIndices(i,0));
            int Cold1=4*(2*FaceCornerIndices(i,0)+1);
            int Colc2=4*(2*FaceCornerIndices(i,1));
            int Cold2=4*(2*FaceCornerIndices(i,1)+1);
            int CurrRowOffset=PresRowOffset+4*2*i;
            
          
            hedra::quatDerivativeIndices(GradRows, GradCols, PresTriCounter, c1TriPoses, CurrRowOffset, Colc1);
            hedra::quatDerivativeIndices(GradRows, GradCols, PresTriCounter, d1TriPoses, CurrRowOffset, Cold1);
            hedra::quatDerivativeIndices(GradRows, GradCols, PresTriCounter, c2TriPoses, CurrRowOffset, Colc2);
            hedra::quatDerivativeIndices(GradRows, GradCols, PresTriCounter, d2TriPoses, CurrRowOffset, Cold2);
            
            CurrRowOffset=PresRowOffset+4*(2*i+1);
            
            hedra::quatDerivativeIndices(GradRows, GradCols, PresTriCounter+64, c1TriPoses, CurrRowOffset, Colc1);
            hedra::quatDerivativeIndices(GradRows, GradCols, PresTriCounter+64, d1TriPoses, CurrRowOffset, Cold1);
            hedra::quatDerivativeIndices(GradRows, GradCols, PresTriCounter+64, c2TriPoses, CurrRowOffset, Colc2);
            hedra::quatDerivativeIndices(GradRows, GradCols, PresTriCounter+64, d2TriPoses, CurrRowOffset, Cold2);
            
            PresTriCounter+=128;
        }
        
        /**********************Closeness Energy***************************************/
        CloseTriOffset=PresTriOffset+128*FaceCornerIndices.rows();
        CloseRowOffset=PresRowOffset+4*2*FaceCornerIndices.rows();
        
        for (int i=0;i<SolutionSize;i++){
            GradRows(CloseTriOffset+i)=CloseRowOffset+i;
            GradCols(CloseTriOffset+i)=i;
            GradValues(CloseTriOffset+i)=CloseFactor;
        }
        
        /************************Compatibility Constraint*****************************/
        CompTriOffset=CloseTriOffset+SolutionSize;
        CompRowOffset=CloseRowOffset+SolutionSize;
        int CompTriCounter=CompTriOffset;
        for (int i=0;i<FaceCornerIndices.rows();i++){
            
            int Colc1=4*(2*FaceCornerIndices(i,0));
            int Cold1=4*(2*FaceCornerIndices(i,0)+1);
            int Colc2=4*(2*FaceCornerIndices(i,1));
            int Cold2=4*(2*FaceCornerIndices(i,1)+1);
            int CurrRowOffset=CompRowOffset+4*i;
            
            hedra::quatDerivativeIndices(GradRows, GradCols, CompTriCounter, c1TriPoses, CurrRowOffset, Colc1);
            hedra::quatDerivativeIndices(GradRows, GradCols, CompTriCounter, d1TriPoses, CurrRowOffset, Cold1);
            hedra::quatDerivativeIndices(GradRows, GradCols, CompTriCounter, c2TriPoses, CurrRowOffset, Colc2);
            hedra::quatDerivativeIndices(GradRows, GradCols, CompTriCounter, d2TriPoses, CurrRowOffset, Cold2);
            
            CompTriCounter+=64;
        }
        
        /************************Imaginarity Constraint******************************/
        ImagTriOffset=CompTriOffset+64*FaceCornerIndices.rows();
        ImagRowOffset=CompRowOffset+4*FaceCornerIndices.rows();
        for (int i=0;i<F.rows();i++){
            for (int j=0;j<4;j++){
                GradRows(ImagTriOffset+4*2*i+j)=ImagRowOffset+i;
                GradCols(ImagTriOffset+4*2*i+j)=4*(2*i)+j;
                GradRows(ImagTriOffset+4*(2*i+1)+j)=ImagRowOffset+i;
                GradCols(ImagTriOffset+4*(2*i+1)+j)=4*(2*i+1)+j;
            }
        }
        
        /***************************First cd constraints*****************************/
        FirstcdTriOffset=ImagTriOffset+8*F.rows();
        FirstcdRowOffset=ImagRowOffset+F.rows();
        int FirstcdTriCounter=FirstcdTriOffset;
        for (int j=0;j<8;j++){
            GradRows(FirstcdTriOffset+j)=FirstcdRowOffset+j;
            GradCols(FirstcdTriOffset+j)=j;
            GradValues(FirstcdTriOffset+j)=1.0;
        }
        
        if (isExactMC){
            //not currently doing
        }
        
        /*MatrixXi GradIndexMat(GradRows.size(),3);
        for (int i=0;i<GradRows.size();i++)
            GradIndexMat(i)=i;
        GradIndexMat.col(1)=GradRows;
        GradIndexMat.col(2)=GradCols;
        cout<<"GradIndices: "<<GradIndexMat<<endl;
        
        cout<<"PresRowOffset: "<<PresRowOffset<<endl;
        cout<<"CloseRowOffset: "<<CloseRowOffset<<endl;
        cout<<"CompRowOffset: "<<CompRowOffset<<endl;
        cout<<"ImagRowOffset: "<<ImagRowOffset<<endl;
        cout<<"FirstcdRowOffset: "<<FirstcdRowOffset<<endl;*/
    }
    
    
    
    
    void UpdateEnergy(const Eigen::VectorXd& CurrSolution)
    {
        using namespace Eigen;
        using namespace std;
        
        for (int i=0;i<FaceCornerIndices.rows();i++){
            RowVector4d c1=CurrSolution.segment(4*2*FaceCornerIndices(i,0),4).transpose();
            RowVector4d d1=CurrSolution.segment(4*(2*FaceCornerIndices(i,0)+1),4).transpose();
            RowVector4d c2=CurrSolution.segment(4*2*FaceCornerIndices(i,1),4).transpose();
            RowVector4d d2=CurrSolution.segment(4*(2*FaceCornerIndices(i,1)+1),4).transpose();
            RowVector4d zi=OrigVq.row(FaceCornerIndices(i,2));
            RowVector4d zj=OrigVq.row(FaceCornerIndices(i,3));
            RowVector4d zij=zj-zi;
            
            RowVector4d G1i=QMult(c1,zi)+d1;
            RowVector4d G2i=QMult(c2,zi)+d2;
            RowVector4d G1j=QMult(c1,zj)+d1;
            RowVector4d G2j=QMult(c2,zj)+d2;
            
            RowVector4d g=PresJumps.segment(4*i,4);
            PresVec.segment(4*(2*i),4)=PresFactor*((QMult(G2i,g)-G1i).transpose());
            PresVec.segment(4*(2*i+1),4)=PresFactor*((QMult(QMult(G1j, zij), QMult(QConj(g), QInv(zij)))-G2j).transpose());
            
        }
        
        
        CloseVec<<CloseFactor*(CurrSolution-InitSolution);
      
        UpdateConstraints(CurrSolution);
        
        TotalVec<<PresVec, CloseVec, ConstVec;
        
    }
    
    
    void UpdateGradient(const Eigen::VectorXd& CurrSolution){
        
        using namespace Eigen;
        using namespace std;
        PresTriOffset=0;
        PresRowOffset=0;
        Vector4i c1TriPoses; c1TriPoses<<0,16,32,48;
        Vector4i d1TriPoses; d1TriPoses<<4,20,36,52;
        Vector4i c2TriPoses; c2TriPoses<<8,24,40,56;
        Vector4i d2TriPoses; d2TriPoses<<12,28,44,60;
        
        
        /**************************Prescription Gradient********************************************/
        int PresTriCounter=PresTriOffset;
        for (int i=0;i<FaceCornerIndices.rows();i++){
            
            RowVector4d zi=OrigVq.row(FaceCornerIndices(i,2));
            RowVector4d zj=OrigVq.row(FaceCornerIndices(i,3));
            RowVector4d zij=zj-zi;
            RowVector4d izij=QInv(zij);
            RowVector4d g=PresJumps.segment(4*i,4);
            
            
            //derivative of G2i*g-G1i=(c2*zi+d2)*g-c1*zi-d1
            
            hedra::quatDerivativeValues(GradValues, PresTriCounter, c1TriPoses, PresFactor*UnitQuat, -zi, false, false);
            hedra::quatDerivativeValues(GradValues, PresTriCounter, d1TriPoses, PresFactor*UnitQuat, -UnitQuat, false, false);
            hedra::quatDerivativeValues(GradValues, PresTriCounter, c2TriPoses, PresFactor*UnitQuat, QMult(zi,g), false, false);
            hedra::quatDerivativeValues(GradValues, PresTriCounter, d2TriPoses, PresFactor*UnitQuat, g, false, false);
            
           //derivatives of G1j*e*conj(g)*inv(e)-G2j=(c1*zj+d1)*zij*conj(g)*inv(zij)-c2*zj-d2
    
            hedra::quatDerivativeValues(GradValues, PresTriCounter+64, c1TriPoses, PresFactor*UnitQuat, QMult(QMult(zj, zij), QMult(QConj(g), izij)), false, false);
            hedra::quatDerivativeValues(GradValues, PresTriCounter+64, d1TriPoses, PresFactor*UnitQuat, QMult(zij, QMult(QConj(g), izij)), false, false);
            hedra::quatDerivativeValues(GradValues, PresTriCounter+64, c2TriPoses, PresFactor*UnitQuat, -zj, false, false);
            hedra::quatDerivativeValues(GradValues, PresTriCounter+64, d2TriPoses, PresFactor*UnitQuat, -UnitQuat, false, false);
            
            PresTriCounter+=128;
        }
        
        //closeness is constant
        
        /***************************Compatibility Gradient*****************************************/
        int CompTriCounter=CompTriOffset;
        for (int i=0;i<FaceCornerIndices.rows();i++){
            
            RowVector4d c1=CurrSolution.segment(4*2*FaceCornerIndices(i,0),4).transpose();
            RowVector4d d1=CurrSolution.segment(4*(2*FaceCornerIndices(i,0)+1),4).transpose();
            RowVector4d c2=CurrSolution.segment(4*2*FaceCornerIndices(i,1),4).transpose();
            RowVector4d d2=CurrSolution.segment(4*(2*FaceCornerIndices(i,1)+1),4).transpose();
            RowVector4d zi=OrigVq.row(FaceCornerIndices(i,2));
            RowVector4d zj=OrigVq.row(FaceCornerIndices(i,3));
            RowVector4d zij=zj-zi;
            RowVector4d izij=QInv(zij);
       
            
            RowVector4d G1i=QMult(c1,zi)+d1;
            RowVector4d G2i=QMult(c2,zi)+d2;
            RowVector4d G1j=QMult(c1,zj)+d1;
            RowVector4d G2j=QMult(c2,zj)+d2;
            
            
            
            //derivatives of Gi1*inv(zij)*conj(Gj1)-Gi2*inv(zij)*conj(Gj2)=
            //(c1*zi+d1)*inv(zij)*(conj(zj)*conj(c1)+conj(d1))-(c2*zi+d2)*inv(zij)*(conj(zj)*conj(c2)+conj(d2))
            
            
            //c1 derivatives
            hedra::quatDerivativeValues(GradValues, CompTriCounter, c1TriPoses, UnitQuat, QMult(zi,QMult(izij,QConj(G1j))), false, false);
            hedra::quatDerivativeValues(GradValues, CompTriCounter, c1TriPoses, QMult(G1i,QMult(izij,QConj(zj))), UnitQuat, true, true);
            
            //d1 derivatives
            hedra::quatDerivativeValues(GradValues, CompTriCounter, d1TriPoses, UnitQuat, QMult(izij,QConj(G1j)), false, false);
            hedra::quatDerivativeValues(GradValues, CompTriCounter, d1TriPoses, QMult(G1i,izij), UnitQuat, true, true);
            
            //c2 derivatives
            hedra::quatDerivativeValues(GradValues, CompTriCounter, c2TriPoses, -UnitQuat, QMult(zi,QMult(izij,QConj(G2j))), false, false);
            hedra::quatDerivativeValues(GradValues, CompTriCounter, c2TriPoses, QMult(G2i,QMult(izij,QConj(zj))), -UnitQuat, true, true);
            
            //d2 derivatives
            hedra::quatDerivativeValues(GradValues, CompTriCounter, d2TriPoses, -UnitQuat, QMult(izij,QConj(G2j)), false, false);
            hedra::quatDerivativeValues(GradValues, CompTriCounter, d2TriPoses, QMult(G2i,izij), -UnitQuat, true, true);
            
            CompTriCounter+=64;
        }
        
        /*****************************Imaginarity Gradient***********************/
        
        for (int i=0;i<F.rows();i++){
            
            //derivative of re(c*conj(d))=0=rc*rd+<Vc,Vd>
            RowVector4d c=CurrSolution.segment(4*2*i,4);
            RowVector4d d=CurrSolution.segment(4*(2*i+1),4);
            
            for (int j=0;j<4;j++){
                GradValues(ImagTriOffset+4*2*i+j)=d(j);
                GradValues(ImagTriOffset+4*(2*i+1)+j)=c(j);
            }
        }
        
        
        //firstcd is constant
        
    }
    
    void Reformulate(int CurrIter, int MaxIterations, const Eigen::VectorXd& CurrSolution, double PrevError)
    {
        using namespace Eigen;
        using namespace std;
        double rate=ConstVec.lpNorm<Infinity>()/PrevError;
        double ReduceRate=min(rate/2.0,1.0);
        
        InitSolution=CurrSolution;
        PresFactor*=0.75-0.25*(1.0-ReduceRate);
        
    }
    
    void UpdateConstraints(const Eigen::VectorXd& CurrSolution)
    {
        using namespace Eigen;
        using namespace std;
        for (int i=0;i<FaceCornerIndices.rows();i++){
            RowVector4d c1=CurrSolution.segment(4*2*FaceCornerIndices(i,0),4).transpose();
            RowVector4d d1=CurrSolution.segment(4*(2*FaceCornerIndices(i,0)+1),4).transpose();
            RowVector4d c2=CurrSolution.segment(4*2*FaceCornerIndices(i,1),4).transpose();
            RowVector4d d2=CurrSolution.segment(4*(2*FaceCornerIndices(i,1)+1),4).transpose();
            RowVector4d zi=OrigVq.row(FaceCornerIndices(i,2));
            RowVector4d zj=OrigVq.row(FaceCornerIndices(i,3));
            RowVector4d zij=zj-zi;
            
            RowVector4d G1i=QMult(c1,zi)+d1;
            RowVector4d G2i=QMult(c2,zi)+d2;
            RowVector4d G1j=QMult(c1,zj)+d1;
            RowVector4d G2j=QMult(c2,zj)+d2;
            
            CompVec.segment(4*i,4)=QMult(QMult(G1i, QInv(zij)), QConj(G1j))-QMult(QMult(G2i, QInv(zij)), QConj(G2j));
        }
        
        
        for (int i=0;i<F.rows();i++){
            RowVector4d c=CurrSolution.segment(4*2*i,4).transpose();
            RowVector4d d=CurrSolution.segment(4*(2*i+1),4).transpose();
            
            ImagVec(i)=QMult(c,QConj(d))(0);
        }
        
        FirstcdVec<<CurrSolution.segment(0,8)-Firstcd;
        
        /*if (isExactMC){
         MCVec(2*i)=G2i.squaredNorm()*g.squaredNorm()-G1i.squaredNorm();
         MCVec(2*i+1)=G1j.squaredNorm()*g.squaredNorm()-G2j.squaredNorm();
         
         }*/
        
        if (!isExactMC)
            ConstVec<<CompVec, ImagVec, FirstcdVec;
        else
            ConstVec<<CompVec, ImagVec, FirstcdVec, MCVec;
    
    }
    
    bool isTerminate()
    {
        return(ConstVec.lpNorm<Eigen::Infinity>() < ConstTolerance);
    }
};




#endif
