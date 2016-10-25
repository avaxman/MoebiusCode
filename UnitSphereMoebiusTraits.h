//
//  GetMobius3DCoeffs.h
//  testigl
//
//  Created by Amir Vaxman on 26/12/14.
//  Copyright (c) 2014 Amir Vaxman. All rights reserved.
//

#ifndef testigl_GetMobius3DCoeffs_h
#define testigl_GetMobius3DCoeffs_h



#include <hedra/QuaternionOps.h>


//both meshes have to be on the same unit sphere, 3 points interpolating, and the rest have the same radius
class UnitSphereMoebiusTraitsTraits{
public:
    
    Matrix4d SourceVq;
    Matrix4d TargetVq;
    
    VectorXd InitSolution;
    double SolutionSize;
    
    double CloseFactor;
    double ConstTolerance;
    double ValueTolerance;
 
    RowVector4d UnitQuat;
    
    VectorXd CloseVec;
    VectorXd ImagVec;
    VectorXd PointVec;
    VectorXd UnitVec;
    
    VectorXd TotalVec;
    VectorXd ConstVec;
    
    VectorXi GradRows, GradCols;
    VectorXd GradValues;
    
    double CloseTriOffset, CloseRowOffset;
    double ImagTriOffset, ImagRowOffset;
    double PosTriOffset, PosRowOffset;
    double UnitTriOffset, UnitRowOffset;
    
    void Initialize(const Matrix4d& inSourceVq, const Matrix4d& inTargetVq,  const VectorXd& inInitSolution){
        
        
        InitSolution=inInitSolution;
        SourceVq=inSourceVq;
        TargetVq=inTargetVq;
        
        UnitQuat<<1.0,0.0,0.0,0.0;
        SolutionSize=16;
    
        CloseFactor=10e-6;
        ConstTolerance=10e-11;
        
        ImagVec.resize(6);
        PointVec.resize(4*3);
        UnitVec.resize(1);
        CloseVec.resize(SolutionSize);
        
        TotalVec.resize(CloseVec.size()+ImagVec.size()+PointVec.size()+UnitVec.size());
        ConstVec.resize(ImagVec.size()+PointVec.size()+UnitVec.size());
        
        //createing gradient pattern
        GradRows.resize(SolutionSize+(64+8+8)+64*3+SolutionSize);
        GradCols.resize(GradRows.size());
        GradValues.resize(GradRows.size());
        
        
        /*******************Closeness energy*************/
        CloseTriOffset=0;
        CloseRowOffset=0;
        for (int i=0;i<SolutionSize;i++){
            GradRows(CloseTriOffset+i)=CloseRowOffset+i;
            GradCols(CloseTriOffset+i)=i;
            GradValues(CloseTriOffset+i)=CloseFactor;
        }
        
        /*******************Imaginarity Constraints************/
        ImagTriOffset=CloseTriOffset+SolutionSize;
        ImagRowOffset=CloseRowOffset+SolutionSize;
        
        //a*conj(c) in imH
        for (int j=0;j<4;j++){
            GradRows(ImagTriOffset+j)=ImagRowOffset;
            GradCols(ImagTriOffset+j)=j;
            GradRows(ImagTriOffset+j+4)=ImagRowOffset;
            GradCols(ImagTriOffset+j+4)=8+j;
        }
        
        //b*conj(d) in imH
        for (int j=0;j<4;j++){
            GradRows(ImagTriOffset+8+j)=ImagRowOffset+1;
            GradCols(ImagTriOffset+8+j)=j+4;
            GradRows(ImagTriOffset+8+j+4)=ImagRowOffset+1;
            GradCols(ImagTriOffset+8+j+4)=12+j;
        }
        
        //conj(a)*d-conj(b)*c = (1,0,0,0) derivatives
        Vector4i aTriPoses; aTriPoses<<0,16,32,48;
        Vector4i bTriPoses; bTriPoses<<4,20,36,52;
        Vector4i cTriPoses; cTriPoses<<8,24,40,56;
        Vector4i dTriPoses; dTriPoses<<12,28,44,60;
        hedra::quatDerivativeIndices(GradRows, GradCols, ImagTriOffset+16, aTriPoses, ImagRowOffset+2, 0);
        hedra::quatDerivativeIndices(GradRows, GradCols, ImagTriOffset+16, bTriPoses, ImagRowOffset+2, 4);
        hedra::quatDerivativeIndices(GradRows, GradCols, ImagTriOffset+16, cTriPoses, ImagRowOffset+2, 8);
        hedra::quatDerivativeIndices(GradRows, GradCols, ImagTriOffset+16, dTriPoses, ImagRowOffset+2, 12);
        
        /****************************Positional Constraints******************************/
        PosTriOffset=ImagTriOffset+64+8+8;
        PosRowOffset=ImagRowOffset+6;
        for (int i=0;i<3;i++){
            //a*q+b-w*c*q-w*d
            hedra::quatDerivativeIndices(GradRows, GradCols, PosTriOffset+64*i, aTriPoses, PosRowOffset+4*i, 0);
            hedra::quatDerivativeIndices(GradRows, GradCols, PosTriOffset+64*i, bTriPoses, PosRowOffset+4*i, 4);
            hedra::quatDerivativeIndices(GradRows, GradCols, PosTriOffset+64*i, cTriPoses, PosRowOffset+4*i, 8);
            hedra::quatDerivativeIndices(GradRows, GradCols, PosTriOffset+64*i, dTriPoses, PosRowOffset+4*i, 12);
            
        }
        
        /********************Unit Constraints***********************************************/
        UnitTriOffset=PosTriOffset+64*3;
        UnitRowOffset=PosRowOffset+4*3;
        
        for (int i=0;i<SolutionSize;i++){
            GradRows(UnitTriOffset+i)=UnitRowOffset;
            GradCols(UnitTriOffset+i)=i;
        }
    }
    
    
    void UpdateEnergy(const VectorXd& CurrSolution)
    {
        CloseVec<<CloseFactor*(CurrSolution-InitSolution);
        UpdateConstraints(CurrSolution);
        TotalVec<<CloseVec, ConstVec;
    }
    
    
    void UpdateGradient(const VectorXd& CurrSolution)
    {
        RowVector4d a=CurrSolution.segment(0,4);
        RowVector4d b=CurrSolution.segment(4,4);
        RowVector4d c=CurrSolution.segment(8,4);
        RowVector4d d=CurrSolution.segment(12,4);
        
        //closeness is constant
        
        /****************Imaginarity Gradient***********************/
        //a*conj(c)
        for (int j=0;j<4;j++){
            GradValues(ImagTriOffset+j)=c(j);
            GradValues(ImagTriOffset+j+4)=a(j);
        }
        
        //b*conj(d)
        for (int j=0;j<4;j++){
            GradValues(ImagTriOffset+8+j)=d(j);
            GradValues(ImagTriOffset+12+j)=b(j);
        }
        
        //conj(a)*d-conj(b)*c derivatives
        Vector4i aTriPoses; aTriPoses<<0,16,32,48;
        Vector4i bTriPoses; bTriPoses<<4,20,36,52;
        Vector4i cTriPoses; cTriPoses<<8,24,40,56;
        Vector4i dTriPoses; dTriPoses<<12,28,44,60;
        hedra::quatDerivativeValues(GradValues, ImagTriOffset+16, aTriPoses,  UnitQuat, d, true, false);
        hedra::quatDerivativeValues(GradValues, ImagTriOffset+16, bTriPoses, -UnitQuat, c, true, false);
        hedra::quatDerivativeValues(GradValues, ImagTriOffset+16, cTriPoses, -QConj(b), UnitQuat, false, false);
        hedra::quatDerivativeValues(GradValues, ImagTriOffset+16, dTriPoses, QConj(a), UnitQuat, false, false);

        
        /****************************Point Constraints******************************/
        for (int i=0;i<3;i++){
            //a*q+b-w*c*q-w*d
            RowVector4d q=SourceVq.row(i);
            RowVector4d w=TargetVq.row(i);
            hedra::quatDerivativeValues(GradValues, PosTriOffset+64*i, aTriPoses, UnitQuat, q, false, false);
            hedra::quatDerivativeValues(GradValues, PosTriOffset+64*i, bTriPoses, UnitQuat, UnitQuat, false, false);
            hedra::quatDerivativeValues(GradValues, PosTriOffset+64*i, cTriPoses, -w,q,  false, false);
            hedra::quatDerivativeValues(GradValues, PosTriOffset+64*i, dTriPoses, -w, UnitQuat, false, false);
        }
        
        /*****************************Unit Constraints******************************/
        

        
        //(diff(|aq+b|^2)/da=2*a*sum(q.^2)-2*QMult(b,q)
        RowVector4d diffa=2*a*SourceVq.row(3).squaredNorm()-2*QMult(b,SourceVq.row(3));
        //(diff(|aq+b|^2)/db=2*b+2*QMult(a,q)
        RowVector4d diffb=2*b+2*QMult(a,SourceVq.row(3));
        //(diff(|cq+d|^2)/dc=2*c*sum(q.^2)-2*QMult(d,q)
        RowVector4d diffc=2*c*SourceVq.row(3).squaredNorm()-2*QMult(d,SourceVq.row(3));
        //(diff(|cq+d|^2)/dd=2*d+2*QMult(c,q)
        RowVector4d diffd=2*d+2*QMult(c,SourceVq.row(3));
        
        for (int j=0;j<4;j++){
            GradValues(UnitTriOffset+j)=diffa(j);
            GradValues(UnitTriOffset+4+j)=diffb(j);
            GradValues(UnitTriOffset+8+j)=-diffc(j);
            GradValues(UnitTriOffset+12+j)=-diffd(j);
        }
    }
    
    void Reformulate(int NumIteration, int MaxIteration, VectorXd& CurrSolution, double PrevError)
    {
        InitSolution=CurrSolution;
    }
    
    
    void UpdateConstraints(const VectorXd& CurrSolution)
    {
        RowVector4d a=CurrSolution.segment(0,4);
        RowVector4d b=CurrSolution.segment(4,4);
        RowVector4d c=CurrSolution.segment(8,4);
        RowVector4d d=CurrSolution.segment(12,4);
        
        
        //a*conj(c) in imH
        RowVector4d ac=QMult(a,QConj(c));
        ImagVec(0)=ac(0);
        
        //b*conj(d) in imH
        RowVector4d bd=QMult(b,QConj(d));
        ImagVec(1)=bd(0);
        
        //conj(a)*d-conj(b)*c = UnitQuat
        RowVector4d adcb=QMult(QConj(a), d)-QMult(QConj(b), c)-UnitQuat;
        ImagVec.segment(2,4)=adcb.transpose();
        
        
        for (int i=0;i<3;i++){
            RowVector4d Res=QMult(a,SourceVq.row(i))+b-QMult(TargetVq.row(i),QMult(c,SourceVq.row(i))+d);
            PointVec.segment(4*i,4)=Res.transpose();
        }
        
        //unit-preserving
        UnitVec(0)=(QMult(a,SourceVq.row(3))+b).squaredNorm()-(QMult(c,SourceVq.row(3))+d).squaredNorm();
        
        
        ConstVec<<ImagVec, PointVec, UnitVec;
    }
    
    bool isTerminate(){
        return (ConstVec.lpNorm<Infinity>()<ConstTolerance);
        
    }
};


#endif
