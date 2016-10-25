//
//  PrescribeEdgeJumps.h
//  testigl
//
//  Created by Amir Vaxman on 22/12/14.
//  Copyright (c) 2014 Amir Vaxman. All rights reserved.
//

#ifndef testigl_PrescribeEdgeJumps_h
#define testigl_PrescribeEdgeJumps_h

#include <Eigen/core>
#include <hedra/QuaternionOps.h>
#include <hedra/quaternionic_derivatives.h>


class PrescribeEdgeJumps2D{
public:
    
    typedef std::complex<double> Complex;
    Eigen::MatrixXi D, F;
    Eigen::VectorXi InnerEdges;
    
    Eigen::VectorXcd OrigVc;
    Eigen::MatrixXi FaceCornerIndices;  //rows of f1, f2, zi, zj (oriented)
    Eigen::VectorXcd PresJumps;
    Eigen::Vector2cd Firstcd;
    Eigen::VectorXcd InitSolution;
    Eigen::VectorXcd CurrSolution;
    bool isExactMC, isExactIAP;
    
    double SmoothFactor;
    double CloseFactor;
    double ConstTolerance;
    int SolutionSize;
    
    
    Eigen::VectorXcd PresVec;
    Eigen::VectorXcd CompVec;
    Eigen::VectorXcd CloseVec;
    Eigen::Vector2cd FirstcdVec;
    Eigen::VectorXd TotalVec;
    Eigen::VectorXcd ConstVec;
    Eigen::VectorXd MCVec;
    Eigen::VectorXd IAPVec;
    
    Eigen::VectorXi GradRows;
    Eigen::VectorXi GradCols;
    Eigen::VectorXd GradValues;
    Eigen::VectorXi ComplexGradRows, ComplexGradCols;
    Eigen::VectorXcd ComplexGradValues;
    
    int PresTriOffset, PresRowOffset;
    int CloseTriOffset, CloseRowOffset;
    int CompTriOffset, CompRowOffset;
    int FirstcdTriOffset, FirstcdRowOffset;
    int MCTriOffset,MCRowOffset;
    int ComplexTriOffset, ComplexRowOffset;
    
    
    void Initialize(const Eigen::VectorXcd& inOrigVc, const Eigen::MatrixXi& inD, const Eigen::MatrixXi& inF, Eigen::MatrixXi& inFaceCornerIndices, const bool& inisExactMC, const bool& inisExactIAP){
        
        using namespace Eigen;
        using namespace std;
        
        F=inF; D=inD;
        FaceCornerIndices=inFaceCornerIndices;
        OrigVc=inOrigVc;
        isExactMC=inisExactMC;
        isExactIAP=inisExactIAP;
        
        CloseFactor=10e-6;
        ConstTolerance=10e-7;
        
        SolutionSize=2*F.rows();
        
        Firstcd<<Complex(0.0), Complex(1.0);
        
        CloseVec.resize(SolutionSize);
        PresVec.resize(2*FaceCornerIndices.rows());
        CompVec.resize(FaceCornerIndices.rows());
        ConstVec.resize(CompVec.size()+FirstcdVec.size());
        CurrSolution.resize(SolutionSize);
        
        if (isExactMC){
            MCVec.resize(2*FaceCornerIndices.rows());
            IAPVec.resize(0);
        } else if (isExactIAP){
            IAPVec.resize(2*FaceCornerIndices.rows());
            MCVec.resize(0);
        } else {
            IAPVec.resize(0);
            MCVec.resize(0);
        }
        
        if (!isExactMC && !isExactIAP)
            TotalVec.resize(2*(CloseVec.size()+PresVec.size()+CompVec.size()+FirstcdVec.size()));
        else
            TotalVec.resize(2*(CloseVec.size()+PresVec.size()+CompVec.size()+FirstcdVec.size())+MCVec.size());
        
        
        /************************Creating Gradient Pattern*******************/
        
        ComplexGradRows.resize(8*FaceCornerIndices.rows()+SolutionSize+4*FaceCornerIndices.rows()+2);
        ComplexGradCols.resize(ComplexGradRows.size());
        ComplexGradValues.resize(ComplexGradRows.size());
        
        if (!isExactMC && !isExactIAP){
            GradRows.resize(4*ComplexGradRows.size());
            GradCols.resize(4*ComplexGradRows.size());
            GradValues.resize(4*ComplexGradRows.size());
        } else {
            
            if (isExactMC){
                GradRows.resize(4*ComplexGradRows.size()+16*FaceCornerIndices.rows());
                GradCols.resize(4*ComplexGradRows.size()+16*FaceCornerIndices.rows());
                GradValues.resize(4*ComplexGradRows.size()+16*FaceCornerIndices.rows());
            } else{
                //currently out
                //GradRows.resize(4*ComplexGradRows.size()+E2V.rows());
                //GradCols.resize(4*CompelxGradRows.size()+E2V.rows());
                //GradValues.resize(4*ComplexGradRows.size()+E2V.rows());
            }
            
        }
        
        
        /******************Prescription Energy**********************/
        PresTriOffset=0;
        PresRowOffset=0;
        for (int i=0;i<FaceCornerIndices.rows();i++){
            ComplexGradRows(PresTriOffset+8*i)=PresRowOffset+2*i;
            ComplexGradCols(PresTriOffset+8*i)=2*FaceCornerIndices(i,0);
            
            ComplexGradRows(PresTriOffset+8*i+1)=PresRowOffset+2*i;
            ComplexGradCols(PresTriOffset+8*i+1)=2*FaceCornerIndices(i,0)+1;
            
            ComplexGradRows(PresTriOffset+8*i+2)=PresRowOffset+2*i;
            ComplexGradCols(PresTriOffset+8*i+2)=2*FaceCornerIndices(i,1);
            
            ComplexGradRows(PresTriOffset+8*i+3)=PresRowOffset+2*i;
            ComplexGradCols(PresTriOffset+8*i+3)=2*FaceCornerIndices(i,1)+1;
            
            ComplexGradRows(PresTriOffset+8*i+4)=PresRowOffset+2*i+1;
            ComplexGradCols(PresTriOffset+8*i+4)=2*FaceCornerIndices(i,0);
            
            ComplexGradRows(PresTriOffset+8*i+5)=PresRowOffset+2*i+1;
            ComplexGradCols(PresTriOffset+8*i+5)=2*FaceCornerIndices(i,0)+1;
            
            ComplexGradRows(PresTriOffset+8*i+6)=PresRowOffset+2*i+1;
            ComplexGradCols(PresTriOffset+8*i+6)=2*FaceCornerIndices(i,1);
            
            ComplexGradRows(PresTriOffset+8*i+7)=PresRowOffset+2*i+1;
            ComplexGradCols(PresTriOffset+8*i+7)=2*FaceCornerIndices(i,1)+1;

            
        }
        
        /*********************Closeness Energy**************************/
        CloseTriOffset=PresTriOffset+8*FaceCornerIndices.rows();
        CloseRowOffset=PresRowOffset+2*FaceCornerIndices.rows();
        for (int i=0;i<SolutionSize;i++){
            ComplexGradRows(CloseTriOffset+i)=CloseRowOffset+i;
            ComplexGradCols(CloseTriOffset+i)=i;
            ComplexGradValues(CloseTriOffset+i)=CloseFactor;
        }
        
        /**********************Comp Energy*****************************/
        CompTriOffset=CloseTriOffset+SolutionSize;
        CompRowOffset=CloseRowOffset+SolutionSize;
        for (int i=0;i<FaceCornerIndices.rows();i++){
            ComplexGradRows(CompTriOffset+4*i)=CompRowOffset+i;
            ComplexGradCols(CompTriOffset+4*i)=2*FaceCornerIndices(i,0);
            
            ComplexGradRows(CompTriOffset+4*i+1)=CompRowOffset+i;
            ComplexGradCols(CompTriOffset+4*i+1)=2*FaceCornerIndices(i,0)+1;
            
            ComplexGradRows(CompTriOffset+4*i+2)=CompRowOffset+i;
            ComplexGradCols(CompTriOffset+4*i+2)=2*FaceCornerIndices(i,1);
            
            ComplexGradRows(CompTriOffset+4*i+3)=CompRowOffset+i;
            ComplexGradCols(CompTriOffset+4*i+3)=2*FaceCornerIndices(i,1)+1;
            
        }
        
        /************************Firstcd Energy***********************/
        FirstcdTriOffset=CompTriOffset+4*FaceCornerIndices.rows();
        FirstcdRowOffset=CompRowOffset+FaceCornerIndices.rows();
        
        ComplexGradRows(FirstcdTriOffset)=FirstcdRowOffset;
        ComplexGradCols(FirstcdTriOffset)=0;
        ComplexGradValues(FirstcdTriOffset)=1.0;
        
        ComplexGradRows(FirstcdTriOffset+1)=FirstcdRowOffset+1;
        ComplexGradCols(FirstcdTriOffset+1)=1;
        ComplexGradValues(FirstcdTriOffset+1)=1.0;
        
        
        ComplexRowOffset=FirstcdRowOffset+2;
        ComplexTriOffset=FirstcdTriOffset+2;
        
        
        
        
        cout<<"ComplexGradRows: "<<ComplexGradRows<<endl;
        cout<<"ComplexGradCols: "<<ComplexGradCols<<endl;
        
        //creating the real-valued pattern and adding MC\IAP constraints
        //[Real -imag; imag real]
        for (int i=0;i<ComplexGradRows.size();i++){
            
            //real upper left
            GradRows(2*i)=ComplexGradRows(i);
            GradCols(2*i)=ComplexGradCols(i);
            GradValues(2*i)=ComplexGradValues(i).real();
            
            //-imag upper right
            GradRows(2*i+1)=ComplexGradRows(i);
            GradCols(2*i+1)=SolutionSize+ComplexGradCols(i);
            GradValues(2*i+1)=-ComplexGradValues(i).imag();
            
            //imag lower left
            GradRows(2*i+2*ComplexGradRows.size())=ComplexRowOffset+ComplexGradRows(i);
            GradCols(2*i+2*ComplexGradRows.size())=ComplexGradCols(i);
            GradValues(2*i+2*ComplexGradRows.size())=ComplexGradValues(i).imag();
            
            //real lower right
            GradRows(2*i+2*ComplexGradRows.size()+1)=ComplexRowOffset+ComplexGradRows(i);
            GradCols(2*i+2*ComplexGradRows.size()+1)=SolutionSize+ComplexGradCols(i);
            GradValues(2*i+2*ComplexGradRows.size()+1)=ComplexGradValues(i).real();
        }
        
        
        //currently not doing IAP
        if (isExactMC){
            MCRowOffset=ComplexRowOffset*2;
            MCTriOffset=4*ComplexGradRows.size();
            
            for (int i=0;i<FaceCornerIndices.rows();i++){
                for (int j=0;j<2;j++) {
                    
                    GradRows(CompTriOffset+16*i+8*j)=MCRowOffset+2*i+j;
                    GradCols(CompTriOffset+16*i+8*j)=2*FaceCornerIndices(i,0);
                    
                    GradRows(CompTriOffset+16*i+1+8*j)=MCRowOffset+2*i+j;
                    GradCols(CompTriOffset+16*i+1+8*j)=2*F.rows()+2*FaceCornerIndices(i,0);
                    
                    GradRows(CompTriOffset+16*i+2+8*j)=MCRowOffset+2*i+j;
                    GradCols(CompTriOffset+16*i+2+8*j)=2*FaceCornerIndices(i,0)+1;
                    
                    GradRows(CompTriOffset+16*i+3+8*j)=MCRowOffset+2*i+j;
                    GradCols(CompTriOffset+16*i+3+8*j)=2*F.rows()+2*FaceCornerIndices(i,0)+1;
                    
                    GradRows(CompTriOffset+16*i+4+8*j)=MCRowOffset+2*i+j;
                    GradCols(CompTriOffset+16*i+4+8*j)=2*FaceCornerIndices(i,1);
                    
                    GradRows(CompTriOffset+16*i+5+8*j)=MCRowOffset+2*i+j;
                    GradCols(CompTriOffset+16*i+5+8*j)=2*F.rows()+2*FaceCornerIndices(i,1);
                    
                    GradRows(CompTriOffset+16*i+6+8*j)=MCRowOffset+2*i+j;
                    GradCols(CompTriOffset+16*i+6+8*j)=2*FaceCornerIndices(i,1)+1;
                    
                    GradRows(CompTriOffset+16*i+6+8*j)=MCRowOffset+2*i+j;
                    GradCols(CompTriOffset+16*i+6+8*j)=2*F.rows()+2*FaceCornerIndices(i,1)+1;
                }
            }
            
        }
        
    }
    
    
    
    void UpdateEnergy(const Eigen::VectorXd& CurrSolutionReal){
        
        using namespace Eigen;
        using namespace std;
        
        CurrSolution.array().real()<<CurrSolutionReal.head(CurrSolutionReal.size()/2);
        CurrSolution.array().imag()<<CurrSolutionReal.tail(CurrSolutionReal.size()/2);
        
        for (int i=0;i<FaceCornerIndices.rows();i++){
            Complex c1=CurrSolution(2*FaceCornerIndices(i,0));
            Complex d1=CurrSolution(2*FaceCornerIndices(i,0)+1);
            Complex c2=CurrSolution(2*FaceCornerIndices(i,1));
            Complex d2=CurrSolution(2*FaceCornerIndices(i,1)+1);
            Complex zi=OrigVc(FaceCornerIndices(i,2));
            Complex zj=OrigVc(FaceCornerIndices(i,3));
            Complex g=PresJumps(i);
            PresVec(2*i)=SmoothFactor*((c1*zi+d1)*g-(c2*zi+d2));
            PresVec(2*i+1)=SmoothFactor*((c2*zj+d2)*g-(c1*zj+d1));
            
            CompVec(i)=(c1*zi+d1)*(c1*zj+d1)-(c2*zi+d2)*(c2*zj+d2);
            
        }
        
        CloseVec<<CloseFactor*(CurrSolution-InitSolution);
        
        UpdateConstraints(CurrSolutionReal);
        
        if (!isExactMC && !isExactIAP)
            TotalVec<<PresVec.real(), CloseVec.real(), CompVec.real(), FirstcdVec.real(),
                      PresVec.imag(), CloseVec.imag(), CompVec.imag(), FirstcdVec.imag();
        else
            TotalVec<<PresVec.real(), CloseVec.real(), CompVec.real(), FirstcdVec.real(),
                      PresVec.imag(), CloseVec.imag(), CompVec.imag(), FirstcdVec.imag(),  MCVec, IAPVec;
        
    }
    
    void UpdateGradient(const Eigen::VectorXd& CurrSolutionReal){
        
        using namespace Eigen;
        using namespace std;
        
        CurrSolution.array().real()<<CurrSolutionReal.head(CurrSolutionReal.size()/2);
        CurrSolution.array().imag()<<CurrSolutionReal.tail(CurrSolutionReal.size()/2);
        
        
        //Energy and compatibility
        for (int i=0;i<FaceCornerIndices.rows();i++){
            
            Complex c1=CurrSolution(2*FaceCornerIndices(i,0));
            Complex d1=CurrSolution(2*FaceCornerIndices(i,0)+1);
            Complex c2=CurrSolution(2*FaceCornerIndices(i,1));
            Complex d2=CurrSolution(2*FaceCornerIndices(i,1)+1);
            Complex zi=OrigVc(FaceCornerIndices(i,2));
            Complex zj=OrigVc(FaceCornerIndices(i,3));
            Complex g=PresJumps(i);

            ComplexGradValues(PresTriOffset+8*i)=SmoothFactor*zi*g;
            ComplexGradValues(PresTriOffset+8*i+1)=SmoothFactor*g;
            ComplexGradValues(PresTriOffset+8*i+2)=-SmoothFactor*zi;
            ComplexGradValues(PresTriOffset+8*i+3)=-SmoothFactor*1.0;
            
            ComplexGradValues(PresTriOffset+8*i+4)=-SmoothFactor*zj;
            ComplexGradValues(PresTriOffset+8*i+5)=-SmoothFactor*1.0;
            ComplexGradValues(PresTriOffset+8*i+6)=SmoothFactor*g*zj;
            ComplexGradValues(PresTriOffset+8*i+7)=SmoothFactor*g;
            
            
            ComplexGradValues(CompTriOffset+4*i)=zi*(c1*zj+d1)+zj*(c1*zi+d1);
            ComplexGradValues(CompTriOffset+4*i+1)=(c1*zi+d1)+(c1*zj+d1);
            ComplexGradValues(CompTriOffset+4*i+2)=-(zj*(c2*zi+d2)+zi*(c2*zj+d2));
            ComplexGradValues(CompTriOffset+4*i+3)=-((c2*zi+d2)+(c2*zj+d2));
            
        }
        
        
        //Firstcd and closeness are constant
        
        //updating real values from complex ones
        for (int i=0;i<ComplexGradRows.size();i++){
            GradValues(2*i)=   ComplexGradValues(i).real();
            GradValues(2*i+1)=-ComplexGradValues(i).imag();
            GradValues(2*i+2*ComplexGradRows.size())= ComplexGradValues(i).imag();
            GradValues(2*i+1+2*ComplexGradRows.size())= ComplexGradValues(i).real();
        }
        
        if (isExactMC){
            //[ 2*cx*zx^2 + 2*dx*zx + 2*cx*zy^2 + 2*dy*zy,
            //2*cy*zx^2 + 2*dy*zx + 2*cy*zy^2 - 2*dx*zy,
            //2*dx + 2*cx*zx - 2*cy*zy,
            //2*dy + 2*cx*zy + 2*cy*zx]
            for (int i=0;i<FaceCornerIndices.rows();i++){
                Complex c1=CurrSolution(2*FaceCornerIndices(i,0));
                Complex d1=CurrSolution(2*FaceCornerIndices(i,0)+1);
                Complex c2=CurrSolution(2*FaceCornerIndices(i,1));
                Complex d2=CurrSolution(2*FaceCornerIndices(i,1)+1);
                Complex zi=OrigVc(FaceCornerIndices(i,2));
                Complex zj=OrigVc(FaceCornerIndices(i,3));
                double gabs2=PresJumps(i).real()*PresJumps(i).real()+PresJumps(i).imag()*PresJumps(i).imag();
                
                //derivatives of |(c1*zi+d1)|^2*|g|^2-|(c2*zi+d2)|^2
                //c1.real() derivative
                GradValues(MCTriOffset+16*i)=(2*c1.real()*zi.real()*zi.real()  + 2*d1.real()*zi.real() + 2*c1.real()*zi.imag()*zi.imag() + 2*d1.imag()*zi.imag())*gabs2;
                //c1.imag() derivative
                GradValues(MCTriOffset+16*i+1)=(2*c1.imag()*zi.real()*zi.real() + 2*d1.imag()*zi.real() + 2*c1.imag()*zi.imag()*zi.imag() - 2*d1.real()*zi.imag())*gabs2,
                //d1.real() derivative
                GradValues(MCTriOffset+16*i+2)=(2*d1.real() + 2*c1.real()*zi.real() - 2*c1.imag()*zi.imag())*gabs2;
                //d1.real() derivative
                GradValues(MCTriOffset+16*i+3)=(2*d1.imag() + 2*c1.real()*zi.imag() + 2*c1.imag()*zi.real())*gabs2;
                
                //c2.real() derivative
                GradValues(MCTriOffset+16*i+4)=-(2*c2.real()*zi.real()*zi.real() + 2*d2.real()*zi.real() + 2*c2.real()*zi.imag()*zi.imag() + 2*d2.imag()*zi.imag());
                //c2.imag() derivative
                GradValues(MCTriOffset+16*i+5)=-(2*c2.imag()*zi.real()*zi.real()+ 2*d2.imag()*zi.real() + 2*c2.imag()*zi.imag()*zi.imag() - 2*d2.real()*zi.imag());
                //d2.real() derivative
                GradValues(MCTriOffset+16*i+6)=-(2*d2.real() + 2*c2.real()*zi.real() - 2*c2.imag()*zi.imag());
                //d2.real() derivative
                GradValues(MCTriOffset+16*i+7)=-(2*d2.imag() + 2*c2.real()*zi.imag() + 2*c2.imag()*zi.real());
                
                //derivatives of |(c1*zj+d1)|^2-|(c2*zj+d2)|^2 * |g|^2
                //c1.real() derivative
                GradValues(MCTriOffset+16*i+8)=2*c1.real()*zj.real()*zj.real() + 2*d1.real()*zj.real() + 2*c1.real()*zj.imag()*zj.imag() + 2*d1.imag()*zj.imag();
                //c1.imag() derivative
                GradValues(MCTriOffset+16*i+9)=2*c1.imag()*zj.real()*zj.real() + 2*d1.imag()*zj.real() + 2*c1.imag()*zj.imag()*zj.imag() - 2*d1.real()*zj.imag(),
                //d1.real() derivative
                GradValues(MCTriOffset+16*i+10)=2*d1.real() + 2*c1.real()*zj.real() - 2*c1.imag()*zj.imag();
                //d1.real() derivative
                GradValues(MCTriOffset+16*i+11)=2*d1.imag() + 2*c1.real()*zj.imag() + 2*c1.imag()*zj.real();
                
                //c2.real() derivative
                GradValues(MCTriOffset+16*i+12)=-(2*c2.real()*zj.real()*zj.real() + 2*d2.real()*zj.real() + 2*c2.real()*zj.imag()*zj.imag() + 2*d2.imag()*zj.imag())*gabs2;
                //c2.imag() derivative
                GradValues(MCTriOffset+16*i+13)=-(2*c2.imag()*zj.real()*zj.real() + 2*d2.imag()*zj.real() + 2*c2.imag()*zj.imag()*zj.imag() - 2*d2.real()*zj.imag())*gabs2;
                //d2.real() derivative
                GradValues(MCTriOffset+16*i+14)=-(2*d2.real() + 2*c2.real()*zj.real() - 2*c2.imag()*zj.imag())*gabs2;
                //d2.real() derivative
                GradValues(MCTriOffset+16*i+15)=-(2*d2.imag() + 2*c2.real()*zj.imag() + 2*c2.imag()*zj.real())*gabs2;
            }
        }
    }
    
    
    //only to be called after UpdateConstraints is updated!!
    void Reformulate(int CurrIter, int MaxIterations, const Eigen::VectorXd& CurrSolutionReal, double PrevError)
    {
        using namespace Eigen;
        using namespace std;
        VectorXcd CurrSolution(CurrSolutionReal.size()/2);
        CurrSolution.array().real()=CurrSolutionReal.head(CurrSolutionReal.size()/2);
        CurrSolution.array().imag()=CurrSolutionReal.tail(CurrSolutionReal.size()/2);
        InitSolution=CurrSolution;
        
        double rate=ConstVec.lpNorm<Infinity>()/PrevError;
        double ReduceRate=min(rate/2.0,1.0);
        
        SmoothFactor*=0.9-0.7*(1.0-ReduceRate);
        
    }
    
    void UpdateConstraints(const Eigen::VectorXd& CurrSolutionReal)
    {
        using namespace Eigen;
        using namespace std;
        VectorXcd CurrSolution(CurrSolutionReal.size()/2);
        CurrSolution.array().real()=CurrSolutionReal.head(CurrSolutionReal.size()/2);
        CurrSolution.array().imag()=CurrSolutionReal.tail(CurrSolutionReal.size()/2);
        
        for (int i=0;i<FaceCornerIndices.rows();i++){
            Complex c1=CurrSolution(2*FaceCornerIndices(i,0));
            Complex d1=CurrSolution(2*FaceCornerIndices(i,0)+1);
            Complex c2=CurrSolution(2*FaceCornerIndices(i,1));
            Complex d2=CurrSolution(2*FaceCornerIndices(i,1)+1);
            Complex zi=OrigVc(FaceCornerIndices(i,2));
            Complex zj=OrigVc(FaceCornerIndices(i,3));
            
            CompVec(i)=(c1*zi+d1)*(c1*zj+d1)-(c2*zj+d2)*(c2*zi+d2);
        }
        
        FirstcdVec=CurrSolution.head(2)-Firstcd;
        
        ConstVec<<CompVec, FirstcdVec;
        
    }
    
    
    //only to be called after UpdateConstraints is updated!!
    bool isTerminate(){
        return (ConstVec.lpNorm<Eigen::Infinity>()<ConstTolerance);
    }
};




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
