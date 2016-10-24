//
//  DeformTraits.h
//  testigl
//
//  Created by Amir Vaxman on 30/09/14.
//  Copyright (c) 2014 Amir Vaxman. All rights reserved.
//

#ifndef testigl_DeformTraits_h
#define testigl_DeformTraits_h

#include "QuaternionOps.h"

typedef std::complex<double> Complex;

class DeformTraitsEdgeDeviation2D{
public:
    VectorXcd OrigVc;
    MatrixXi E2V, D, F;
    VectorXi ConstIndices;//, VarIndices;
    VectorXcd ComplexConstPoses;
    VectorXcd InitSolution;
    
    double SmoothFactor;
    double CloseFactor;
    double ConstTolerance;
    double RigidRatio;
    
    bool isExactMC;
    bool isExactIAP;
    
    int SolutionSize;
    
    VectorXcd CurrSolution;  //the complex one
    VectorXcd CurrMobRecips;
    VectorXcd CurrLocations;
    VectorXcd CurrEdgeDeviations;
    VectorXcd CurrEdges;
    VectorXcd AMAPVec;
    VectorXcd MobVec;
    VectorXcd PosVec;
    VectorXcd RigidVec;
    VectorXcd CloseVec;
    VectorXcd DeviationVec;
    VectorXd MCVec;
    VectorXd IAPVec;
    VectorXcd FCR;
    
    SparseMatrix<Complex> d0;
    
    VectorXi ComplexGradRows;
    VectorXi ComplexGradCols;
    VectorXcd ComplexGradValues;
    
    
    VectorXi GradRows;
    VectorXi GradCols;
    VectorXd GradValues;
    
    //into the complex values
    int AMAPTriOffset, AMAPRowOffset;
    int RigidTriOffset, RigidRowOffset;
    int PosTriOffset, PosRowOffset;
    int MobTriOffset, MobRowOffset;
    int CloseTriOffset, CloseRowOffset;
    int DeviationTriOffset, DeviationRowOffset;
    int ComplexRowOffset, ComplexTriOffset;
    
    //into the real values
    int MCTriOffset, MCRowOffset;
    int IAPTriOffset, IAPRowOffset;
    
    VectorXd TotalVec;
    VectorXcd ConstVec;
    

    void Initialize(const VectorXcd& inOrigVc, const MatrixXi& inD, const MatrixXi& inF, const MatrixXi& inE2V, const VectorXi& inConstIndices, const VectorXcd& inFCR, bool inisExactMC, bool inisExactIAP, double inRigidRatio){
        
        F=inF; D=inD; E2V=inE2V;
        ConstIndices=inConstIndices; //VarIndices=inVarIndices;
        OrigVc=inOrigVc;
        RigidRatio=inRigidRatio;
        FCR=inFCR;
        
        isExactMC=inisExactMC;
        isExactIAP=inisExactIAP;
        
        
        SolutionSize=OrigVc.rows()+OrigVc.rows();
        
        if (isExactMC || isExactIAP)
            SolutionSize+=E2V.rows();
        
        CurrSolution.resize(SolutionSize);
        
        CloseFactor=10e-3;
        
        if (isExactMC || isExactIAP)
            ConstTolerance=10e-9;
        else
            ConstTolerance=10e-9;
        
        //creating difference operator
        d0.resize(E2V.rows(), OrigVc.rows());
        d0.reserve(2*E2V.rows());
        vector<Triplet<Complex> > d0Tris(2*E2V.rows());
        for (int i=0;i<E2V.rows();i++){
            d0Tris[2*i]=Triplet<Complex>(i,E2V(i,0), -1.0);
            d0Tris[2*i+1]=Triplet<Complex>(i,E2V(i,1), 1.0);
        }
        
        d0.setFromTriplets(d0Tris.begin(), d0Tris.end());
        
        //Allocating intermediate and output vectors
        AMAPVec.resize(E2V.rows());
        RigidVec.resize(E2V.rows());
        PosVec.resize(ConstIndices.size());
        CloseVec.resize(SolutionSize);
        MobVec.resize(D.sum()-3*D.rows());
        if (isExactMC || isExactIAP)
            DeviationVec.resize(E2V.rows());
        if (isExactMC){
            MCVec.resize(E2V.rows());
            IAPVec.resize(0);
        } else if (isExactIAP){
            IAPVec.resize(E2V.rows());
            MCVec.resize(0);
        } else {
            IAPVec.resize(0);
            MCVec.resize(0);
        }
        CurrMobRecips.resize(OrigVc.rows());
        CurrLocations.resize(OrigVc.rows());
        CurrEdgeDeviations.resize(E2V.rows());
        CurrEdges.resize(E2V.rows());
        
        if (!isExactMC && !isExactIAP)
            TotalVec.resize(2*(AMAPVec.size()+RigidVec.size()+CloseVec.size()+PosVec.size()+MobVec.size()));
        else
            TotalVec.resize(2*(AMAPVec.size()+RigidVec.size()+CloseVec.size()+PosVec.size()+MobVec.size()+DeviationVec.size())+MCVec.size()+IAPVec.size());
        
        
        if (!isExactMC && !isExactIAP)
            ConstVec.resize(PosVec.size()+MobVec.size());
        else
            ConstVec.resize(PosVec.size()+MobVec.size()+DeviationVec.size()+MCVec.size()+IAPVec.size());
                                   
        
        /**************************************************Creating Gradient Pattern*******************************************/
        
        if (!isExactMC && !isExactIAP){
            ComplexGradRows.resize(4*E2V.rows()+2*E2V.rows()+SolutionSize+ConstIndices.size()+4*FCR.size());
            ComplexGradCols.resize(ComplexGradRows.size());
            ComplexGradValues.resize(ComplexGradRows.size());
            
            GradRows.resize(4*ComplexGradRows.size());
            GradCols.resize(4*ComplexGradRows.size());
            GradValues.resize(4*ComplexGradRows.size());
        } else {
            ComplexGradRows.resize(4*E2V.rows()+2*E2V.rows()+SolutionSize+ConstIndices.size()+4*FCR.size()+5*E2V.rows());
            ComplexGradCols.resize(ComplexGradRows.size());
            ComplexGradValues.resize(ComplexGradRows.size());
            
            if (isExactMC){
                GradRows.resize(4*ComplexGradRows.size()+2*E2V.rows());
                GradCols.resize(4*ComplexGradRows.size()+2*E2V.rows());
                GradValues.resize(4*ComplexGradRows.size()+2*E2V.rows());
            } else{
                GradRows.resize(4*ComplexGradRows.size()+E2V.rows());
                GradCols.resize(4*ComplexGradRows.size()+E2V.rows());
                GradValues.resize(4*ComplexGradRows.size()+E2V.rows());
            }

        }
        

        
        /*********************AMAPEnergy**********************/
        //Xj*Xi*zij part
        AMAPTriOffset=0;
        AMAPRowOffset=0;
        for (int i=0;i<E2V.rows();i++){
            
            //the actual values are updated in the gradient
            ComplexGradRows(AMAPTriOffset+4*i)=AMAPRowOffset+i;
            ComplexGradCols(AMAPTriOffset+4*i)=E2V(i,0);
            
            ComplexGradRows(AMAPTriOffset+4*i+1)=AMAPRowOffset+i;
            ComplexGradCols(AMAPTriOffset+4*i+1)=E2V(i,1);
            
            ComplexGradRows(AMAPTriOffset+4*i+2)=AMAPRowOffset+i;
            ComplexGradCols(AMAPTriOffset+4*i+2)=OrigVc.rows()+E2V(i,0);
            
            ComplexGradRows(AMAPTriOffset+4*i+3)=AMAPRowOffset+i;
            ComplexGradCols(AMAPTriOffset+4*i+3)=OrigVc.rows()+E2V(i,1);
            
            
        }
        
        /*******************Rigidity Energy*******************/
        RigidTriOffset=AMAPTriOffset+4*E2V.rows();
        RigidRowOffset=AMAPRowOffset+E2V.rows();
        
        for (int i=0;i<E2V.rows();i++){
            ComplexGradRows(RigidTriOffset+2*i)=RigidRowOffset+i;
            ComplexGradCols(RigidTriOffset+2*i)=E2V(i,0);
            
            ComplexGradRows(RigidTriOffset+2*i+1)=RigidRowOffset+i;
            ComplexGradCols(RigidTriOffset+2*i+1)=E2V(i,1);

        }
        
        /************************Closeness Energy******************/
        CloseTriOffset=RigidTriOffset+2*E2V.rows();
        CloseRowOffset=RigidRowOffset+E2V.rows();
        for (int i=0;i<SolutionSize;i++){
            ComplexGradRows(CloseTriOffset+i)=CloseRowOffset+i;
            ComplexGradCols(CloseTriOffset+i)=i;
            ComplexGradValues(CloseTriOffset+i)=CloseFactor;
        }
        
        /******************Positional Constraints************/
        
        PosTriOffset=CloseTriOffset+SolutionSize;
        PosRowOffset=CloseRowOffset+SolutionSize;
        for (int i=0;i<ConstIndices.size();i++){
            ComplexGradRows(PosTriOffset+i)=PosRowOffset+i;
            ComplexGradCols(PosTriOffset+i)=OrigVc.rows()+ConstIndices(i);
            ComplexGradValues(PosTriOffset+i)=1.0;
        }
        
        /*****************Mobius-Equivalent Constraints********/
        
        //the actual values are updated in the gradient
        MobTriOffset=PosTriOffset+ConstIndices.size();
        MobRowOffset=PosRowOffset+ConstIndices.size();
        int MobTriCounter=0;
        for (int i=0;i<D.rows();i++){
            for (int j=0;j<D(i)-3;j++){
                for (int k=0;k<4;k++){
                    ComplexGradRows(MobTriOffset+4*MobTriCounter+k)=MobRowOffset+MobTriCounter;
                    ComplexGradCols(MobTriOffset+4*MobTriCounter+k)=OrigVc.rows()+F(i,j+k);
                }
                MobTriCounter++;
            }
        }
        
    
        ComplexRowOffset=MobRowOffset+FCR.size();
        ComplexTriOffset=MobTriOffset+4*FCR.size();

        
        if (isExactMC || isExactIAP){  //deviation constraints are included
            //actual values are computed in the gradient update
            DeviationTriOffset=ComplexTriOffset;
            DeviationRowOffset=ComplexRowOffset;
            for (int i=0;i<E2V.rows();i++){
                ComplexGradRows(DeviationTriOffset+5*i)=DeviationRowOffset+i;
                ComplexGradCols(DeviationTriOffset+5*i)=E2V(i,0);
                
                ComplexGradRows(DeviationTriOffset+5*i+1)=DeviationRowOffset+i;
                ComplexGradCols(DeviationTriOffset+5*i+1)=E2V(i,1);
                
                ComplexGradRows(DeviationTriOffset+5*i+2)=DeviationRowOffset+i;
                ComplexGradCols(DeviationTriOffset+5*i+2)=OrigVc.rows()+E2V(i,0);
                
                ComplexGradRows(DeviationTriOffset+5*i+3)=DeviationRowOffset+i;
                ComplexGradCols(DeviationTriOffset+5*i+3)=OrigVc.rows()+E2V(i,1);
                
                ComplexGradRows(DeviationTriOffset+5*i+4)=DeviationRowOffset+i;
                ComplexGradCols(DeviationTriOffset+5*i+4)=2*OrigVc.rows()+i;
            }
            
            ComplexRowOffset+=E2V.rows();
            ComplexTriOffset+=5*E2V.rows();
        }
        
        //cout<<"ComplexGradRows: "<<ComplexGradRows<<endl;
        
        
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
        
        
        
        /*************************Metric-Conformal and Intersection-Angle Preserving Constraints********************************/
        if (isExactMC){
            MCRowOffset=ComplexRowOffset*2;
            MCTriOffset=4*ComplexGradRows.size();
            
            
            //actual values are updated in the gradient function
            for (int i=0;i<E2V.rows();i++){
                GradRows(MCTriOffset+2*i)=MCRowOffset+i;
                GradCols(MCTriOffset+2*i)=2*OrigVc.rows()+i;  //real part of e_ij
                
                GradRows(MCTriOffset+2*i+1)=MCRowOffset+i;
                GradCols(MCTriOffset+2*i+1)=SolutionSize+2*OrigVc.rows()+i;  //imaginary part of e_ij
            }
        } else if (isExactIAP){
            IAPRowOffset=2*ComplexRowOffset;
            IAPTriOffset=4*ComplexGradRows.size();
            for (int i=0;i<E2V.rows();i++){
                GradRows(IAPTriOffset+i)=IAPRowOffset+i;
                GradCols(IAPTriOffset+i)=SolutionSize+2*OrigVc.rows()+i;  //imaginary part of e_ij
                GradValues(IAPTriOffset+i)=1.0;
            }
            
        }
        
        //cout<<"GradRows: "<<GradRows<<endl;
    }
    
    void UpdateEnergy(const VectorXd& CurrSolutionReal){
        
        CurrSolution.array().real()<<CurrSolutionReal.head(CurrSolutionReal.size()/2);
        CurrSolution.array().imag()<<CurrSolutionReal.tail(CurrSolutionReal.size()/2);
        
        CurrMobRecips<<CurrSolution.head(OrigVc.rows());
        CurrLocations<<CurrSolution.segment(OrigVc.rows(),OrigVc.rows());
        
        for (int i=0;i<E2V.rows();i++)
            AMAPVec(i)=SmoothFactor*(CurrMobRecips(E2V(i,0))*(OrigVc(E2V(i,1))-OrigVc(E2V(i,0)))*CurrMobRecips(E2V(i,1))-(CurrLocations(E2V(i,1))-CurrLocations(E2V(i,0))));
        
        RigidVec<<SmoothFactor*RigidRatio*d0*CurrMobRecips;
        CloseVec<<CloseFactor*(CurrSolution-InitSolution);
        
        //to update the solution
        UpdateConstraints(CurrSolutionReal);
        
        if (!isExactMC && !isExactIAP)
            TotalVec<<AMAPVec.real(), RigidVec.real(), CloseVec.real(), PosVec.real(), MobVec.real(),
                      AMAPVec.imag(), RigidVec.imag(), CloseVec.imag(), PosVec.imag(), MobVec.imag();
        
        else
            TotalVec<<AMAPVec.real(), RigidVec.real(), CloseVec.real(), PosVec.real(), MobVec.real(), DeviationVec.real(),
                      AMAPVec.imag(), RigidVec.imag(), CloseVec.imag(), PosVec.imag(), MobVec.imag(), DeviationVec.imag(),MCVec, IAPVec;
        
    }
    
    
    void UpdateGradient(const VectorXd& CurrSolutionReal){
        
        
        CurrSolution.array().real()<<CurrSolutionReal.head(CurrSolutionReal.size()/2);
        CurrSolution.array().imag()<<CurrSolutionReal.tail(CurrSolutionReal.size()/2);
        
        CurrMobRecips<<CurrSolution.head(OrigVc.rows());
        CurrLocations<<CurrSolution.segment(OrigVc.rows(),OrigVc.rows());
        
        
        /*********************AMAPEnergy**********************/
        for (int i=0;i<E2V.rows();i++){
            //Xi Derivative
            ComplexGradValues(4*i)=SmoothFactor*CurrMobRecips(E2V(i,1))*(OrigVc(E2V(i,1))-OrigVc(E2V(i,0)));
            //Xj Derivative
            ComplexGradValues(4*i+1)=SmoothFactor*CurrMobRecips(E2V(i,0))*(OrigVc(E2V(i,1))-OrigVc(E2V(i,0)));
            //zi derivative
            ComplexGradValues(4*i+2)=SmoothFactor;
            //zj derivative
            ComplexGradValues(4*i+3)=-SmoothFactor;
        }
        
        /*******************Rigidity Energy******************/
        
        for (int i=0;i<E2V.rows();i++){
            //zi derivative
            ComplexGradValues(RigidTriOffset+2*i)=-SmoothFactor*RigidRatio;
            //zj derivative
            ComplexGradValues(RigidTriOffset+2*i+1)=SmoothFactor*RigidRatio;
        }
        
        //Closeness gradient is constant
        
        //Positional gradient is constant
        
        /*****************Mobius-Equivalent Constraints********/
        
        //the actual values are updated in the gradient
        int MobCounter=0;
        for (int i=0;i<D.rows();i++){
            for (int j=0;j<D(i)-3;j++){
                Complex Zi=CurrLocations(F(i,j));
                Complex Zj=CurrLocations(F(i,j+1));
                Complex Zk=CurrLocations(F(i,j+2));
                Complex Zl=CurrLocations(F(i,j+3));
                //MobVec(Counter++)=(Zj-Zi)*(Zl-Zk)-FCR(MobCounter)*(Zi-Zl)*(Zk-Zj);
                
                //zi derivative
                ComplexGradValues(MobTriOffset+4*MobCounter)=-(Zl-Zk)-FCR(MobCounter)*(Zk-Zj);
                //zj derivative
                ComplexGradValues(MobTriOffset+4*MobCounter+1)=(Zl-Zk)+FCR(MobCounter)*(Zi-Zl);
                //zk derivative
                ComplexGradValues(MobTriOffset+4*MobCounter+2)=-(Zj-Zi)-FCR(MobCounter)*(Zi-Zl);
                //zl derivative
                ComplexGradValues(MobTriOffset+4*MobCounter+3)=(Zj-Zi)+FCR(MobCounter)*(Zk-Zj);
                
                MobCounter++;
                
            }
        }
        
        
        if (isExactMC || isExactIAP){
            /*****************Deviation constraints***************/
            CurrEdgeDeviations<<CurrSolution.segment(OrigVc.rows()+OrigVc.rows(),E2V.rows());
            for (int i=0;i<E2V.rows();i++){
                //Xi derivative
                ComplexGradValues(DeviationTriOffset+5*i)=CurrMobRecips(E2V(i,1))*(OrigVc(E2V(i,1))-OrigVc(E2V(i,0)));
                //Xj derivative
                ComplexGradValues(DeviationTriOffset+5*i+1)=CurrMobRecips(E2V(i,0))*(OrigVc(E2V(i,1))-OrigVc(E2V(i,0)));
    
                //zi derivative
                ComplexGradValues(DeviationTriOffset+5*i+2)=CurrEdgeDeviations(i);
                
                //zj derivative
                ComplexGradValues(DeviationTriOffset+5*i+3)=-CurrEdgeDeviations(i);
                
                //e_ij derivatve
                ComplexGradValues(DeviationTriOffset+5*i+4)=-(CurrLocations(E2V(i,1))-CurrLocations(E2V(i,0)));
            }
        }
        
        //updating real values from complex ones
        for (int i=0;i<ComplexGradRows.size();i++){
            GradValues(2*i)=   ComplexGradValues(i).real();
            GradValues(2*i+1)=-ComplexGradValues(i).imag();
            GradValues(2*i+2*ComplexGradRows.size())= ComplexGradValues(i).imag();
            GradValues(2*i+1+2*ComplexGradRows.size())= ComplexGradValues(i).real();
        }
        
        
        
        /*************************Metric-Conformal and Intersection-Angle Preserving Constraints********************************/
        if (isExactMC){
            for (int i=0;i<E2V.rows();i++){
                GradValues(MCTriOffset+2*i)=2.0*CurrEdgeDeviations(i).real();
                GradValues(MCTriOffset+2*i+1)=2.0*CurrEdgeDeviations(i).imag();
                
            }
        }
        
        //IAP constraints are constant
        
    }
    
    
    
    //should be called when constraints are updated!!!
    void Reformulate(int NumIteration, int MaxIteration, const VectorXd& CurrSolutionReal, double PrevError)
    {
        VectorXcd CurrSolution(CurrSolutionReal.size()/2);
        CurrSolution.array().real()=CurrSolutionReal.head(CurrSolutionReal.size()/2);
        CurrSolution.array().imag()=CurrSolutionReal.tail(CurrSolutionReal.size()/2);
        InitSolution=CurrSolution;

        double rate=ConstVec.lpNorm<Infinity>()/PrevError;
        double ReduceRate=min(rate/2.0,1.0);
        
        SmoothFactor*=0.9-0.7*(1.0-ReduceRate);

    }
    
    void UpdateConstraints(const VectorXd& CurrSolutionReal)
    {
        VectorXcd CurrSolution(CurrSolutionReal.size()/2);
        CurrSolution.array().real()=CurrSolutionReal.head(CurrSolutionReal.size()/2);
        CurrSolution.array().imag()=CurrSolutionReal.tail(CurrSolutionReal.size()/2);
        
        CurrMobRecips<<CurrSolution.head(OrigVc.rows());
        CurrLocations<<CurrSolution.segment(OrigVc.rows(),OrigVc.rows());
        
        for (int i=0;i<ConstIndices.size();i++)
            PosVec(i)=CurrLocations(ConstIndices(i))-ComplexConstPoses(i);
        
        int MobCounter=0;
        for (int i=0;i<D.rows();i++){
            for (int j=0;j<D(i)-3;j++){
                Complex Zi=CurrLocations(F(i,j));
                Complex Zj=CurrLocations(F(i,j+1));
                Complex Zk=CurrLocations(F(i,j+2));
                Complex Zl=CurrLocations(F(i,j+3));
                MobVec(MobCounter)=(Zj-Zi)*(Zl-Zk)-FCR(MobCounter)*(Zi-Zl)*(Zk-Zj);
                MobCounter++;
            }
        }
        
        if (!isExactMC && !isExactIAP){
            ConstVec<<PosVec, MobVec;
        } else{
            //Deviation constraints
            CurrEdgeDeviations<<CurrSolution.segment(OrigVc.rows()+OrigVc.rows(),E2V.rows());
            CurrEdges<<d0*CurrLocations;
            for (int i=0;i<E2V.rows();i++)
                DeviationVec(i)=CurrMobRecips(E2V(i,0))*(OrigVc(E2V(i,1))-OrigVc(E2V(i,0)))*CurrMobRecips(E2V(i,1))-CurrEdgeDeviations(i)*CurrEdges(i);
            
            if (isExactMC)
                MCVec<<CurrEdgeDeviations.array().abs2()-1;
            if (isExactIAP)
                IAPVec<<CurrEdgeDeviations.array().imag();
            
            ConstVec<<PosVec, MobVec, DeviationVec, MCVec.cast<Complex>(), IAPVec.cast<Complex>();
            //cout<<"MCVec: "<<MCVec<<endl;
        }

    }
    
    
    //should be called when constraints are updated!!
    bool isTerminate()
    {
        return (ConstVec.lpNorm<Infinity>()<ConstTolerance);
    }

};



class DeformTraitsCornerVars3D
{
public:
    
    MatrixXi D, F;
    VectorXi ConstIndices;
    MatrixXi CornerPairs;
    VectorXi CornerOffset;
    SparseMatrix<double> d0;
    VectorXd InitSolution;
    MatrixXd OrigVq;
    MatrixXd ConstPoses;
    
    MatrixXi EdgeCornerPairs;  //for rigidity energy
    MatrixXi EdgeCornerVertices;  //the vertices themselves, for compatibility
    
    int NumCorners;
    int NumEdgePairs;
    int SolutionSize;
    
    bool isExactMC;
    bool isExactIAP;
    
    double SmoothFactor;
    double CloseFactor;
    double ConstTolerance;
    double PosFactor;
    double RigidRatio;
    
    RowVector4d UnitQuat;
    
    
    //intermediate variables
    VectorXd CurrX;
    VectorXd CurrLocations;
    VectorXd AMAPVec;
    VectorXd RigidVec;
    VectorXd CompVec;
    VectorXd CloseVec;
    VectorXd PosVec;
    VectorXd MCVec;
    VectorXd IAPVec;
    
    VectorXd TotalVec;
    VectorXd ConstVec;
    
    //raw gradient triplets for pardiso
    VectorXi GradRows, GradCols;
    VectorXd GradValues;
    
    
    int AMAPTriOffset, AMAPRowOffset;
    int RigidTriOffset, RigidRowOffset;
    int CloseTriOffset, CloseRowOffset;
    int CompTriOffset, CompRowOffset;
    int PosTriOffset, PosRowOffset;
    int MCTriOffset, MCRowOffset;
    int IAPTriOffset, IAPRowOffset;

    
    void Initialize(const MatrixXd& inOrigVq,const MatrixXi& inD , const MatrixXi& inF ,const VectorXi& inConstIndices, const MatrixXi& inCornerPairs, const int inNumCorners, const VectorXi& inCornerOffset, const bool inExactMC, const bool inExactIAP, double inRigidRatio){
        
        
        OrigVq=inOrigVq; F=inF; D=inD;
        ConstIndices=inConstIndices;
        CornerPairs=inCornerPairs;
        NumCorners=inNumCorners;
        CornerOffset=inCornerOffset;
        
        UnitQuat<<1.0,0.0,0.0,0.0;

        //MobConstMat=inMobConstMat;
        isExactMC=inExactMC;
        isExactIAP=inExactIAP;
        RigidRatio=inRigidRatio;
        
        SolutionSize=4*NumCorners+3*OrigVq.rows();
        CurrX.resize(4*NumCorners);
        CurrLocations.resize(3*OrigVq.rows());
        
        VectorXi FaceNums=D;
        FaceNums=(FaceNums.cwiseAbs2()-FaceNums);
        NumEdgePairs=FaceNums.sum()/2;
        
        EdgeCornerPairs.resize(NumEdgePairs,2);
        EdgeCornerVertices.resize(NumEdgePairs,2);
        int CurrPair=0;
        for (int i=0;i<D.rows();i++)
            for (int j=0;j<D(i);j++)
                for (int k=j+1;k<D(i);k++){
                    EdgeCornerPairs.row(CurrPair)<<CornerOffset(i)+j, CornerOffset(i)+k;
                    EdgeCornerVertices.row(CurrPair++)<<F(i,j), F(i,k);
                }
        
        
        AMAPVec.resize(4*CornerPairs.rows());
        RigidVec.resize(4*EdgeCornerPairs.rows());
        CompVec.resize(4*EdgeCornerPairs.rows());
        CloseVec.resize(SolutionSize);
        PosVec.resize(3*ConstIndices.size());
        if (isExactMC){
            MCVec.resize(CornerPairs.rows());
            IAPVec.resize(0);
        } else if (isExactIAP){
            MCVec.resize(0);
            IAPVec.resize(CornerPairs.rows());
        } else{
            MCVec.resize(0);
            IAPVec.resize(0);
        }
        
        TotalVec.resize(AMAPVec.size()+RigidVec.size()+CloseVec.size()+CompVec.size()+PosVec.size()+MCVec.size()+IAPVec.size());
        ConstVec.resize(CompVec.size()+PosVec.size()+MCVec.size()+IAPVec.size());
        
        CloseFactor=10e-6;
        ConstTolerance=10e-7;
        
        //Constructing Gradient Pattern
        
        if (!isExactMC){  //not currently doing IAP
            GradRows.resize(2*4*CornerPairs.rows()+2*4*EdgeCornerPairs.rows()+SolutionSize+38*EdgeCornerPairs.rows()+3*ConstIndices.size());
            GradCols.resize(GradRows.size());
            GradValues.resize(GradRows.size());
        } else{
            GradRows.resize(2*4*CornerPairs.rows()+2*4*EdgeCornerPairs.rows()+SolutionSize+38*EdgeCornerPairs.rows()+3*ConstIndices.size()+2*4*CornerPairs.rows());
            GradCols.resize(GradRows.size());
            GradValues.resize(GradRows.size());
        }
        
        /*******************************AMAP Energy********************************************/
        AMAPTriOffset=0;
        AMAPRowOffset=0;
        for (int i=0;i<CornerPairs.rows();i++){
            for (int j=0;j<4;j++){
                GradRows(AMAPTriOffset+2*(4*i+j))=AMAPRowOffset+4*i+j;
                GradCols(AMAPTriOffset+2*(4*i+j))=4*CornerPairs(i,0)+j;
                GradRows(AMAPTriOffset+2*(4*i+j)+1)=AMAPRowOffset+4*i+j;
                GradCols(AMAPTriOffset+2*(4*i+j)+1)=4*CornerPairs(i,1)+j;
            }
        }
        
        /*******************************Rigidity Energy*****************************************/
        RigidTriOffset=AMAPTriOffset+2*4*CornerPairs.rows();
        RigidRowOffset=AMAPRowOffset+4*CornerPairs.rows();
        for (int i=0;i<EdgeCornerPairs.rows();i++){
            for (int j=0;j<4;j++){
                GradRows(RigidTriOffset+2*(4*i+j))=RigidRowOffset+4*i+j;
                GradCols(RigidTriOffset+2*(4*i+j))=4*EdgeCornerPairs(i,0)+j;
                GradRows(RigidTriOffset+2*(4*i+j)+1)=RigidRowOffset+4*i+j;
                GradCols(RigidTriOffset+2*(4*i+j)+1)=4*EdgeCornerPairs(i,1)+j;
            }
        }

        
        /****************************Closeness Energy*******************/
        CloseTriOffset=RigidTriOffset+2*4*EdgeCornerPairs.rows();
        CloseRowOffset=RigidRowOffset+4*EdgeCornerPairs.rows();
        for (int i=0;i<SolutionSize;i++){
            GradRows(CloseTriOffset+i)=CloseRowOffset+i;
            GradCols(CloseTriOffset+i)=i;
            GradValues(CloseTriOffset+i)=CloseFactor;
        }
        
        /****************************Compatibility Constraints*****************/
        CompTriOffset=CloseTriOffset+SolutionSize;
        CompRowOffset=CloseRowOffset+SolutionSize;
        int CompTriCounter=CompTriOffset;
        Vector4i XiTriPoses; XiTriPoses<<0,8,18,28;
        Vector4i XjTriPoses; XjTriPoses<<4,12,22,32;
        for (int i=0;i<EdgeCornerPairs.rows();i++){
            int ColCorneri=4*EdgeCornerPairs(i,0);
            int ColCornerj=4*EdgeCornerPairs(i,1);
            int CurrRowOffset=4*i;
            //derivative of Xi
            GetQuatDerivativeIndices(GradRows, GradCols, CompTriCounter, XiTriPoses, CompRowOffset+CurrRowOffset, ColCorneri);
            
            //Derivative of Xj
            GetQuatDerivativeIndices(GradRows, GradCols, CompTriCounter, XjTriPoses,CompRowOffset+CurrRowOffset, ColCornerj);
            
            for (int k=0;k<3;k++){
                //wj derivative
                GradRows(CompTriCounter+16+10*k)=CompRowOffset+CurrRowOffset+1+k;
                GradCols(CompTriCounter+16+10*k)=4*NumCorners+3*EdgeCornerVertices(i,1)+k;
                GradValues(CompTriCounter+16+10*k)=-1.0;
                
                //wi derivative
                GradRows(CompTriCounter+17+10*k)=CompRowOffset+CurrRowOffset+1+k;
                GradCols(CompTriCounter+17+10*k)=4*NumCorners+3*EdgeCornerVertices(i,0)+k;
                GradValues(CompTriCounter+17+10*k)=1.0;
            }
            
            CompTriCounter+=38;
        }
        
        /****************************Positional Constraints*******************/
        PosTriOffset=CompTriOffset+38*EdgeCornerPairs.rows();
        PosRowOffset=CompRowOffset+4*EdgeCornerPairs.rows();
        
        for (int i=0;i<ConstIndices.size();i++){
            for (int k=0;k<3;k++){
                GradRows(PosTriOffset+3*i+k)=PosRowOffset+3*i+k;
                GradCols(PosTriOffset+3*i+k)=4*NumCorners+3*ConstIndices(i)+k;
            }
        }
        
        
        /****************************Metric-Conformal Constraints*************/
        if (isExactMC){
            MCTriOffset=PosTriOffset+3*ConstIndices.size();
            MCRowOffset=PosRowOffset+3*ConstIndices.size();
            for (int i=0;i<CornerPairs.rows();i++){
                for (int j=0;j<4;j++){
                    GradRows(MCTriOffset+2*(4*i+j))=MCRowOffset+i;
                    GradCols(MCTriOffset+2*(4*i+j))=4*CornerPairs(i,0)+j;
                    GradRows(MCTriOffset+2*(4*i+j)+1)=MCRowOffset+i;
                    GradCols(MCTriOffset+2*(4*i+j)+1)=4*CornerPairs(i,1)+j;
                }
            }
        }
        
        if (isExactIAP){
            //not currently doing it
        }
        
        /*MatrixXi MatGrads(GradRows.size(),3);
        for (int i=0;i<GradRows.size();i++)
            MatGrads(i,0)=i;
        MatGrads.col(1)=GradRows;
        MatGrads.col(2)=GradCols;
        cout<<"AMAP Tri/Row Offsets: "<<AMAPTriOffset<<","<<AMAPRowOffset<<endl;
        cout<<"Rigid Tri/Row Offsets: "<<RigidTriOffset<<","<<RigidRowOffset<<endl;
        cout<<"Close Tri/row Offsets: "<<CloseTriOffset<<","<<CloseRowOffset<<endl;
        cout<<"Comp Tri/row Offsets: "<<CompTriOffset<<","<<CompRowOffset<<endl;
        cout<<"Pos Tri/row Offsets: "<<PosTriOffset<<","<<PosRowOffset<<endl;
        cout<<"MC Tri/row Offsets: "<<MCTriOffset<<","<<MCRowOffset<<endl;
        cout<<"GradRows and Cols"<<MatGrads<<endl;
        int kaka=8;*/
        
    }
    
    void UpdateEnergy(const VectorXd& CurrSolution){
        
        CurrX<<CurrSolution.head(4*NumCorners);
        CurrLocations<<CurrSolution.tail(3*OrigVq.rows());
        
        for (int i=0;i<CornerPairs.rows();i++)
            AMAPVec.segment(4*i,4)=(CurrX.segment(4*CornerPairs(i,1),4)-CurrX.segment(4*CornerPairs(i,0),4));
        
        AMAPVec.array()*=SmoothFactor;
        
        for (int i=0;i<EdgeCornerPairs.rows();i++)
            RigidVec.segment(4*i,4)=(CurrX.segment(4*EdgeCornerPairs(i,1),4)-CurrX.segment(4*EdgeCornerPairs(i,0),4));

        
        RigidVec.array()*=SmoothFactor*RigidRatio;
        
        CloseVec<<CloseFactor*(CurrSolution-InitSolution);
        
        UpdateConstraints(CurrSolution);
        
        TotalVec<<AMAPVec, RigidVec, CloseVec, ConstVec;

    }
    
    void UpdateGradient(const VectorXd& CurrSolution){
    
        CurrX<<CurrSolution.head(4*NumCorners);
        CurrLocations<<CurrSolution.tail(3*OrigVq.rows());
        
        /*******************************AMAP Energy********************************************/
        for (int i=0;i<CornerPairs.rows();i++){
            for (int j=0;j<4;j++){
                GradValues(AMAPTriOffset+2*(4*i+j))=-SmoothFactor;
                GradValues(AMAPTriOffset+2*(4*i+j)+1)=SmoothFactor;
            }
        }
        
        /*******************************Rigidity Energy*****************************************/
        for (int i=0;i<EdgeCornerPairs.rows();i++){
            for (int j=0;j<4;j++){
                GradValues(RigidTriOffset+2*(4*i+j))=-SmoothFactor*RigidRatio;
                GradValues(RigidTriOffset+2*(4*i+j)+1)=SmoothFactor*RigidRatio;
            }
        }
        
        //closeness energy is constant
        
        /****************************Compatibility Constraints*****************/
        int CompTriCounter=CompTriOffset;
        Vector4i XiTriPoses; XiTriPoses<<0,8,18,28;
        Vector4i XjTriPoses; XjTriPoses<<4,12,22,32;
        for (int i=0;i<EdgeCornerPairs.rows();i++){
            RowVector4d Xi=CurrX.segment(4*EdgeCornerPairs(i,0),4).transpose();
            RowVector4d Xj=CurrX.segment(4*EdgeCornerPairs(i,1),4).transpose();
            RowVector4d RightPart=QMult1(OrigVq.row(EdgeCornerVertices(i,1))-OrigVq.row(EdgeCornerVertices(i,0)),Xj);
            RowVector4d LeftPart=QMult1(QConj1(Xi),OrigVq.row(EdgeCornerVertices(i,1))-OrigVq.row(EdgeCornerVertices(i,0)));
            
  
            //derivative of Xi
            GetQuatDerivativeValues(GradValues, CompTriCounter, XiTriPoses, UnitQuat, RightPart, true, false);
            
            //Derivative of Xj
            GetQuatDerivativeValues(GradValues, CompTriCounter, XjTriPoses, LeftPart, UnitQuat, false, false);
            
            //the other compatibility values are constant
            CompTriCounter+=38;
            
            
        }
    
        
        
        /****************************Positional Constraints*******************/
        for (int i=0;i<ConstIndices.size();i++)
            for (int k=0;k<3;k++)
                GradValues(PosTriOffset+3*i+k)=PosFactor;
        
        /****************************Metric-Conformal Constraints*************/
        if (isExactMC){
            for (int i=0;i<CornerPairs.rows();i++){
                RowVector4d Xi=CurrX.segment(4*CornerPairs(i,0),4).transpose();
                RowVector4d Xj=CurrX.segment(4*CornerPairs(i,1),4).transpose();
                for (int k=0;k<4;k++){
                    GradValues(MCTriOffset+2*(4*i+k))=-2*Xi(k);
                    GradValues(MCTriOffset+2*(4*i+k)+1)=2*Xj(k);
                }
            }
        }

    }
    
    
    
    void Reformulate(int NumIteration, int MaxIteration, const VectorXd& CurrSolution, double PrevError)
    {
        //when error is halved, the smoothness is reduced by slowest, and when error change is zero, smoothness is halved.
        double rate=ConstVec.lpNorm<Infinity>()/PrevError;
        double ReduceRate=min(rate/2.0,1.0);
        
        if (ConstIndices.size()==OrigVq.rows())
            SmoothFactor*=0.9-0.2*(1.0-ReduceRate);
        else
            SmoothFactor*=0.9-0.7*(1.0-ReduceRate);
        
        InitSolution=CurrSolution;
        
        
    }
    
    void UpdateConstraints(const VectorXd& CurrSolution)
    {
        
        CurrX<<CurrSolution.head(4*NumCorners);
        CurrLocations<<CurrSolution.tail(3*OrigVq.rows());
        
        for (int i=0;i<EdgeCornerPairs.rows();i++){
            RowVector4d Xi=CurrX.segment(4*EdgeCornerPairs(i,0),4).transpose();
            RowVector4d Xj=CurrX.segment(4*EdgeCornerPairs(i,1),4).transpose();
            RowVector4d qi=OrigVq.row(EdgeCornerVertices(i,0));
            RowVector4d qj=OrigVq.row(EdgeCornerVertices(i,1));
            RowVector3d wi=CurrLocations.segment(3*EdgeCornerVertices(i,0),3).transpose();
            RowVector3d wj=CurrLocations.segment(3*EdgeCornerVertices(i,1),3).transpose();
            RowVector4d CurrEdgeVector=QMult1(QMult1(QConj1(Xi),qj-qi),Xj);
            CurrEdgeVector.tail(3)-=(wj-wi);
            CompVec.segment(4*i,4)<<CurrEdgeVector.transpose();
        
        }
        
        
        for (int i=0;i<ConstIndices.size();i++)
            PosVec.segment(3*i,3)<<PosFactor*(CurrLocations.segment(3*ConstIndices(i),3)-ConstPoses.row(i).transpose());
        
        
        if (!isExactMC){
            ConstVec<<CompVec, PosVec;
            
        } else {
            for (int i=0;i<CornerPairs.rows();i++)
                MCVec(i)=CurrX.segment(4*CornerPairs(i,1),4).squaredNorm()-CurrX.segment(4*CornerPairs(i,0),4).squaredNorm();
            ConstVec<<CompVec, PosVec, MCVec;
            //cout<<"Inside Const Error: "<<ConstVec.lpNorm<Infinity>()<<endl;
        }
        
    }
    
    //should be called when constraints are updated!!
    bool isTerminate()
    {
        return (ConstVec.lpNorm<Infinity>()<ConstTolerance);
    }

};



#endif
