//
//  InterpolateTraits.h
//  testigl
//
//  Created by Amir Vaxman on 30/09/14.
//  Copyright (c) 2014 Amir Vaxman. All rights reserved.
//

#ifndef MoebiusCode_InterpolateTraits_h
#define MoebiusCode_InterpolateTraits_h

#include "AuxSparse.h"
#include "MobiusFromShapes.h"
#include "InitialSolution3DTraits.h"
#include "QuadConstSolver.h"
#include "QuaternionOps.h"
#include "ExtraFunctions.h"

//currently not treating the const indices at all

class InterpolateTraits2D{
public:
    double OffWeight;
    double PosWeight;
    MatrixXi E2V, E2F, E2Fi, F;
    VectorXi InnerEdges, ConstIndices;
    
    SparseMatrix<Complex> d0;
    SparseMatrix<Complex> IDMat;
    
    VectorXcd SourceVc;
    
    SparseMatrix<Complex> gGlobal;
    
    
    //VectorXd OrigLengths, DeformLengths, PrescribeLengths;
    VectorXcd InitSolution;
    VectorXcd InterpEdgeDeviations;
    VectorXcd ConstPoses;
    Vector2cd Interpcd;
    
    double SmoothFactor;
    double CloseFactor;
    double ConstTolerance;
    double VertexTolerance;
    
    
    void Initialize(const VectorXcd& inSourceVc, const VectorXcd& inTargetVc, const VectorXcd& inInterpEdgeDeviations, const Vector2cd& inInterpcd, const MatrixXi& inF, const MatrixXi& inE2V, const VectorXi& inInnerEdges, const MatrixXi& inE2F, const MatrixXi& inE2Fi, const VectorXi& inConstIndices, const VectorXcd& inConstPoses, const SparseMatrix<Complex> ind0, double t, const VectorXcd& inInitSolution){
        
        F=inF; E2V=inE2V; E2F=inE2F; InnerEdges=inInnerEdges; E2F=inE2F; E2Fi=inE2Fi; ConstIndices=inConstIndices;
        InitSolution=inInitSolution;
        d0=ind0;
        SourceVc=inSourceVc;
        ConstPoses=inConstPoses;
        Interpcd=inInterpcd;
        InterpEdgeDeviations=inInterpEdgeDeviations;
        
        //InitialSolution2DTraits ist;
        //ist.Initialize(ConstIndices, SourceVc, TargetConstPoses, TargetConstX);
        //cout<<"Finding Base Mobius Transformation:"<<endl;
        //QuadConstSolver<InitialSolution2DTraits, Complex> qcs(ist);
        //TargetConstX=qcs.Solve(TargetConstX, 1.0, 25);
        
        //for (int i=0;i<ConstIndices.size();i++)
        //    InterpConstX(i)=exp(log(TargetConstX(i))*t);


        //cout<<"TargetConstX: "<<TargetConstX<<endl;
        
        //retrieving interpolated const positions by Poisson, keeping the first constant
        /*vector<pair<int,int> > ConstPairs;
        for (int i=0;i<ConstIndices.size();i++)
            for (int j=i+1;j<ConstIndices.size();j++)
                ConstPairs.push_back(pair<int,int>(i,j));*/
        
        
        //computing current deform error
        
       
        /*VectorXcd InterpConstDeviations(ConstPairs.size());
        MatrixXcd dConst(ConstPairs.size(), ConstIndices.size());
        dConst.setZero();
        VectorXcd MobEdgeVectors(ConstPairs.size());
        for (int i=0;i<ConstPairs.size();i++){
 
            MobEdgeVectors(i)=InterpConstX(ConstPairs[i].second)*(SourceConstPoses(ConstPairs[i].second)-SourceConstPoses(ConstPairs[i].first))*InterpConstX(ConstPairs[i].first);
            InterpConstDeviations(i)=TargetConstX(ConstPairs[i].second)*(SourceConstPoses(ConstPairs[i].second)-SourceConstPoses(ConstPairs[i].first))*TargetConstX(ConstPairs[i].first);
            InterpConstDeviations(i)/=(TargetConstPoses(ConstPairs[i].second)-TargetConstPoses(ConstPairs[i].first));
            //cout<<"InterpConstDeviations: "<<InterpConstDeviations(i)<<endl;
            InterpConstDeviations(i)=exp(log(InterpConstDeviations(i))*t);
            //cout<<"InterpConstDeviations: "<<InterpConstDeviations(i)<<endl;
            dConst(i,ConstPairs[i].first)=-InterpConstDeviations(i);
            dConst(i,ConstPairs[i].second)=InterpConstDeviations(i);
        }
        
        //cout<<"MobEdgeVectors: "<<MobEdgeVectors<<endl;
        //cout<<"Source Edge Vectors: "<<dConst*SourceConstPoses<<endl;
        //cout<<"InterpConstDeviations: "<<InterpConstDeviations<<endl;
        //cout<<"TargetEdgeDeviations: "<<TargetEdgeDeviations<<endl;
        //cout<<"Error with Target positions:"<<dConst*TargetConstPoses-MobEdgeVectors<<endl;
        VectorXcd torhs(ConstIndices.size()); torhs.setZero();
        torhs(0)=SourceConstPoses(0);
        VectorXcd rhs=MobEdgeVectors-dConst*torhs;
        MatrixXcd dConstSmall=dConst.block(0,1, ConstPairs.size(), ConstIndices.size()-1);
        VectorXcd ConstSolution = (dConstSmall.adjoint()*dConstSmall).llt().solve(dConstSmall.adjoint()*rhs);
        //InterpConstPoses(0)=torhs(0);
        //InterpConstPoses.tail(ConstIndices.size()-1)=ConstSolution;
        
        //cout<<"SourceConstPoses: "<<SourceConstPoses<<endl;
        //cout<<"TargetConstPoses: "<<TargetConstPoses<<endl;
        //cout<<"InterpConstPoses: "<<InterpConstPoses<<endl;*/
        
        //Interpolating edge deviations
        
        
        /*InterpEdgeDeviations.resize(TargetEdgeDeviations.size());
        for (int i=0;i<TargetEdgeDeviations.size();i++)
            InterpEdgeDeviations(i)=exp(log(TargetEdgeDeviations(i))*t);*/
    
        igl::speye(SourceVc.rows()+SourceVc.rows(), IDMat);
        
        SmoothFactor=100;
        CloseFactor=10e-5;
        ConstTolerance=10e-9;
        VertexTolerance=10e-5;
        
        
        
    }
    
    
    
    VectorXcd EvalEnergy(const VectorXcd& CurrSolution){
        
        //VectorXcd CurrSolution(CurrSolutionReal.size()/2);
        //CurrSolution.array().real()=CurrSolutionReal.head(CurrSolutionReal.size()/2);
        //CurrSolution.array().imag()=CurrSolutionReal.tail(CurrSolutionReal.size()/2);
        
        VectorXcd CurrMobRecips=CurrSolution.head(SourceVc.rows());
        VectorXcd CurrLocations=CurrSolution.segment(SourceVc.rows(),SourceVc.rows());
        
        VectorXcd CurrEdges=(d0*CurrLocations).cwiseProduct(InterpEdgeDeviations);
        VectorXcd CurrMobiusEdges(E2V.rows());
        for (int i=0;i<E2V.rows();i++)
            CurrMobiusEdges(i)=CurrMobRecips(E2V(i,0))*(SourceVc(E2V(i,1))-SourceVc(E2V(i,0)))*CurrMobRecips(E2V(i,1));
        
        VectorXcd EnergyVec=SmoothFactor*(CurrMobiusEdges-CurrEdges);
        VectorXcd CloseVec=CloseFactor*(CurrSolution-InitSolution);
        //VectorXcd MobVec=MobConstFactor*MobConstMat*CurrSolution;
        
        VectorXcd PosVec(ConstIndices.size());
        for (int i=0;i<ConstIndices.size();i++)
            PosVec(i)=CurrLocations(ConstIndices(i))-ConstPoses(i);
      

        VectorXcd ConstVec(CloseVec.size()+EnergyVec.size()+PosVec.size());
        ConstVec<<CloseVec, EnergyVec, PosVec;
        return ConstVec;

        
    }
    
    SparseMatrix<Complex> EvalGradient(const VectorXcd& CurrSolution){
        
        //VectorXcd CurrSolution(CurrSolutionReal.size()/2);
        //CurrSolution.array().real()=CurrSolutionReal.head(CurrSolutionReal.size()/2);
        //CurrSolution.array().imag()=CurrSolutionReal.tail(CurrSolutionReal.size()/2);
        
        VectorXcd CurrMobRecips=CurrSolution.head(SourceVc.rows());
        VectorXcd CurrLocations=CurrSolution.segment(SourceVc.rows(),SourceVc.rows());
        
        SparseMatrix<Complex> gCurrEdges=-SparseDiagonalMatrix<Complex>(InterpEdgeDeviations)*d0;
        SparseMatrix<Complex> gCurrMobiusEdges(E2V.rows(), SourceVc.rows());
        vector<Triplet<Complex> > gCurrMobiusEdgesTris;
        for (int i=0;i<E2V.rows();i++){
            gCurrMobiusEdgesTris.push_back(Triplet<Complex>(i,E2V(i,1),CurrMobRecips(E2V(i,0))*(SourceVc(E2V(i,1))-SourceVc(E2V(i,0)))));
            gCurrMobiusEdgesTris.push_back(Triplet<Complex>(i,E2V(i,0),CurrMobRecips(E2V(i,1))*(SourceVc(E2V(i,1))-SourceVc(E2V(i,0)))));
        }
        
        gCurrMobiusEdges.setFromTriplets(gCurrMobiusEdgesTris.begin(), gCurrMobiusEdgesTris.end());
        
        SparseMatrix<Complex> gAMAPEnergy=SparseBlock<Complex>(gCurrMobiusEdges, gCurrEdges,true);
        
        SparseMatrix<Complex> gPos(ConstIndices.size(), CurrSolution.size());
        vector<Triplet<Complex> > gPosTris;
        for (int i=0;i<ConstIndices.size();i++)
            gPosTris.push_back(Triplet<Complex>(i,SourceVc.rows()+ConstIndices(i),1.0));
        
        gPos.setFromTriplets(gPosTris.begin(), gPosTris.end());

        SparseMatrix<Complex> ComplexGradient=SparseBlock<Complex>(SparseBlock<Complex>(CloseFactor*IDMat,SmoothFactor*gAMAPEnergy, false),gPos,false);//, gPos, false);

        return ComplexGradient;
    }
    
    void Reformulate(int CurrIter, int MaxIterations, const VectorXcd& CurrSolution)
    {
        InitSolution=CurrSolution;
        SmoothFactor*=0.95;

    }
    
    double GetConstError(const VectorXcd& CurrSolution)
    {
        
        VectorXcd CurrMobRecips=CurrSolution.head(SourceVc.rows());
        VectorXcd CurrLocations=CurrSolution.segment(SourceVc.rows(),SourceVc.rows());
        
        VectorXcd CurrEdges=(d0*CurrLocations).cwiseProduct(InterpEdgeDeviations);
        VectorXcd CurrMobiusEdges(E2V.rows());
        for (int i=0;i<E2V.rows();i++)
            CurrMobiusEdges(i)=CurrMobRecips(E2V(i,0))*(SourceVc(E2V(i,1))-SourceVc(E2V(i,0)))*CurrMobRecips(E2V(i,1));
        
        VectorXcd EnergyVec=(CurrMobiusEdges-CurrEdges);
        
        VectorXcd PosVec(ConstIndices.size());
        for (int i=0;i<ConstIndices.size();i++)
            PosVec(i)=CurrLocations(ConstIndices(i))-ConstPoses(i);
        
        
        VectorXcd ConstVec(EnergyVec.size()+PosVec.size());
        ConstVec<<EnergyVec, PosVec;

        
        
        //return PosVec.lpNorm<Infinity>();
        return EnergyVec.lpNorm<Infinity>();
    }
    
    bool isTerminate(const VectorXcd& CurrSolution){
        //VectorXcd CurrSolution(CurrSolutionReal.size()/2);
        //CurrSolution.array().real()=CurrSolutionReal.head(CurrSolutionReal.size()/2);
        //CurrSolution.array().imag()=CurrSolutionReal.tail(CurrSolutionReal.size()/2);
        
        VectorXcd CurrLocations=CurrSolution.segment(SourceVc.rows(),SourceVc.rows());
        VectorXcd InitLocations=InitSolution.segment(SourceVc.rows(),SourceVc.rows());
        bool NoVerticesMoved=((CurrLocations-InitLocations).lpNorm<Infinity>()<VertexTolerance);
        cout<<"Vertex Movement: "<<(CurrLocations-InitLocations).lpNorm<Infinity>()<<endl;
        
        bool WithinConstTolerance=(GetConstError(CurrSolution) < ConstTolerance);
        return (WithinConstTolerance || NoVerticesMoved);

    }
};

class InterpolateTraits3D{
public:
    double OffWeight;
    double PosWeight;
    MatrixXi E2V, E2F, E2Fi, F;
    VectorXi InnerEdges, ConstIndices;
    
    
    VectorXd InitSolution;
    VectorXd InterpConstPoses;
    VectorXd InterpEdgeDeviations;
    MatrixXd SourceVq;
    
    MatrixXd ConstPoses;
    
    SparseMatrix<double> d0;
    SparseMatrix<double> IDMat;
    
    double SmoothFactor;
    double CloseFactor;
    double ConstTolerance;
    
    
    void Initialize(const MatrixXd& inSourceVq, const MatrixXd& TargetVq, const MatrixXd& TargetX, const MatrixXd& TargetEdgeDeviations, const MatrixXi& inF, const MatrixXi& inE2V, const VectorXi& inInnerEdges, const MatrixXi& inE2F, const MatrixXi& inE2Fi, const VectorXi& inConstIndices, const MatrixXd& inConstPoses, const SparseMatrix<double>& ind0, const double t, const VectorXd& inInitSolution){
        
        
        F=inF; E2V=inE2V; E2F=inE2F; InnerEdges=inInnerEdges; E2F=inE2F; E2Fi=inE2Fi; ConstIndices=inConstIndices;
        InitSolution=inInitSolution;
        SourceVq=inSourceVq;
        
        ConstPoses=inConstPoses;
        
        d0=ind0;
        
        MatrixXd SourceConstPoses(ConstIndices.size(),4);
        MatrixXd TargetConstPoses(ConstIndices.size(),4);
        InterpConstPoses.resize(3*ConstIndices.size());
        
        MatrixXd InterpConstX(ConstIndices.size(),4);
        VectorXd TargetConstX(4*ConstIndices.size());
        
        for (int i=0;i<ConstIndices.size();i++){
            SourceConstPoses.row(i)=SourceVq.row(ConstIndices(i));
            TargetConstPoses.row(i)=TargetVq.row(ConstIndices(i));
            TargetConstX.segment(4*i,4)=TargetX.row(ConstIndices(i)).transpose();
        }
        
        InitialSolution3DTraits ist;
        ist.Initialize(ConstIndices, SourceVq, TargetConstPoses, TargetConstX);
        QuadConstSolver<InitialSolution3DTraits, double> qcs(ist);
        TargetConstX=qcs.Solve(TargetConstX, 1.0, 25);
        
        for (int i=0;i<ConstIndices.size();i++)
            InterpConstX.row(i)=QExp(QLog(TargetConstX.segment(4*i,4).transpose())*t);
                                             
        
        //cout<<"TargetConstX: "<<TargetConstX<<endl;
        
        //retrieving interpolated const positions by Poisson, keeping the first constant
        vector<pair<int,int> > ConstPairs;
        for (int i=0;i<ConstIndices.size();i++)
            for (int j=i+1;j<ConstIndices.size();j++)
                ConstPairs.push_back(pair<int,int>(i,j));
        
        
        //computing current deform error
        MatrixXd InterpConstDeviations(ConstPairs.size(),4);
        MatrixXd dConst(ConstPairs.size(), ConstIndices.size());
        
        dConst.setZero();
                                             
        MatrixXd MobEdgeVectors(ConstPairs.size(),4);
        
        for (int i=0;i<ConstPairs.size();i++){
            
            MobEdgeVectors.row(i)=QMult1(QConj1(InterpConstX.row(ConstPairs[i].second)),QMult1(SourceConstPoses.row(ConstPairs[i].second)-SourceConstPoses.row(ConstPairs[i].first),InterpConstX.row(ConstPairs[i].first)));
            InterpConstDeviations.row(i)=QMult1(QConj1(TargetConstX.segment(4*ConstPairs[i].second,4).transpose()),QMult1(SourceConstPoses.row(ConstPairs[i].second)-SourceConstPoses.row(ConstPairs[i].first), TargetConstX.segment(4*ConstPairs[i].first,4).transpose()));
            InterpConstDeviations.row(i)=QMult1(InterpConstDeviations.row(i),QInv1(TargetConstPoses.row(ConstPairs[i].second)-TargetConstPoses.row(ConstPairs[i].first)));
            //cout<<"InterpConstDeviations: "<<InterpConstDeviations(i)<<endl;
            InterpConstDeviations.row(i)=QExp(QLog(InterpConstDeviations.row(i))*t);
            //cout<<"InterpConstDeviations: "<<InterpConstDeviations(i)<<endl;
            dConst(i,ConstPairs[i].first)=-1.0;
            dConst(i,ConstPairs[i].second)=1.0;
        }
        
        //cout<<"MobEdgeVectors: "<<MobEdgeVectors<<endl;
        //cout<<"Source Edge Vectors: "<<dConst*SourceConstPoses<<endl;
        //cout<<"InterpConstDeviations: "<<InterpConstDeviations<<endl;
        //cout<<"Error with Target positions:"<<dConst*TargetConstPoses-MobEdgeVectors<<endl;
        MatrixXd torhs(ConstIndices.size(),3); torhs.setZero();
        torhs.row(0)=SourceConstPoses.row(0).tail(3);
        MatrixXd rhs=MobEdgeVectors.block(0,1,MobEdgeVectors.rows(), 3)-dConst*torhs;
                                               
        MatrixXd dConstSmall=dConst.block(0,1, ConstPairs.size(), ConstIndices.size()-1);
                                               
        MatrixXd ConstSolution = (dConstSmall.adjoint()*dConstSmall).llt().solve(dConstSmall.adjoint()*rhs);
        
        InterpConstPoses.setZero();
        InterpConstPoses.segment(0,3)=torhs.row(0).transpose();
        for (int i=1;i<ConstIndices.size();i++)
            InterpConstPoses.segment(3*i,3)=ConstSolution.row(i-1);
        
        //cout<<"SourceConstPoses: "<<SourceConstPoses<<endl;
        //cout<<"TargetConstPoses: "<<TargetConstPoses<<endl;
        //cout<<"InterpConstPoses: "<<InterpConstPoses<<endl;
        
        //Interpolating edge deviations
        InterpEdgeDeviations.resize(TargetEdgeDeviations.size());
        for (int i=0;i<E2V.rows();i++)
            InterpEdgeDeviations.segment(4*i,4)=QExp(QLog(TargetEdgeDeviations.row(i))*t).transpose();
        
        //cout<<"TargetEdgeDeviations: "<<TargetEdgeDeviations<<endl;
        //cout<<"InterpEdgeDeviations: "<<InterpEdgeDeviations<<endl;
        
        igl::speye(4*SourceVq.rows()+3*SourceVq.rows(), IDMat);
        
        SmoothFactor=10;
        CloseFactor=10e-5;
        ConstTolerance=10e-7;
        
        //checking finite element gradient of a random choice
        /*VectorXd CurrSolution=InitSolution+100*VectorXd::Random(InitSolution.rows(), InitSolution.cols());
         VectorXd CurrConst=EvalEnergy(CurrSolution);
         SparseMatrix<double> FEGradient(CurrConst.size(), CurrSolution.size());
         vector<Triplet<double> > FEGradientTris;
         for (int i=0;i<CurrSolution.size();i++){
         VectorXd vh(CurrSolution.size()); vh.setZero(); vh(i)=10e-5;
         VectorXd CurrGradient=(EvalEnergy(CurrSolution+vh)-EvalEnergy(CurrSolution-vh))/(2*10e-5);
         //cout<<CurrGradient<<endl;
         for (int j=0;j<CurrGradient.size();j++)
         if (abs(CurrGradient(j))>10e-7)
         FEGradientTris.push_back(Triplet<double>(j,i,CurrGradient(j)));
         }
         
         FEGradient.setFromTriplets(FEGradientTris.begin(), FEGradientTris.end());
         SparseMatrix<double> OrigGradient=EvalGradient(CurrSolution);
         SparseMatrix<double> DiffMat=FEGradient-OrigGradient;
         double maxcoeff=-32767.0;
         int Maxi,Maxj;
         for (int k=0; k<DiffMat.outerSize(); ++k)
         for (SparseMatrix<double>::InnerIterator it(DiffMat,k); it; ++it){
         if (maxcoeff<abs(it.value())){
         maxcoeff=abs(it.value());
         Maxi=it.row();
         Maxj=it.col();
         
         }
         if (abs(it.value())>10e-6){
         cout<<"Gradient Discrepancy at: ("<<it.row()<<","<<it.col()<<") of "<<it.value()<<endl;
         cout<<"FE Gradient: "<<FEGradient.coeffRef(it.row(), it.col())<<endl;
         cout<<"Calculated gradient: "<<OrigGradient.coeffRef(it.row(),it.col())<<endl;
         }
         }
         cout<<"Maximum gradient difference: "<<maxcoeff<<endl;
         cout<<"At Location: ("<<Maxi<<","<<Maxj<<")"<<endl;
         //cout<<"E2V.rows()="<<E2V.rows()<<endl;
         //cout<<"OrigVq.rows()="<<OrigVq.rows()<<endl;
         //cout<<"VarIndices.rows()="<<VarIndices.rows()<<endl;*/


    }
    
    VectorXd EvalEnergy(const VectorXd& CurrSolution){
        
        VectorXd CurrX=CurrSolution.head(4*SourceVq.rows());
        VectorXd CurrLocations=CurrSolution.segment(4*SourceVq.rows(),3*SourceVq.rows());
        VectorXd ExtEdgeVectors(4*E2V.rows()); ExtEdgeVectors.setZero();
        VectorXd CurrWij=d0*CurrLocations;
        for (int i=0;i<E2V.rows();i++){
            RowVector4d ExtCurrWij; ExtCurrWij<<0,CurrWij.segment(3*i,3).transpose();
            ExtEdgeVectors.segment(4*i,4)=QMult1(InterpEdgeDeviations.segment(4*i,4).transpose(),ExtCurrWij).transpose();
        }
        
        VectorXd EnergyVec(4*E2V.rows());
        
        VectorXd MobEdgeVectors(4*E2V.rows());;
        
        for (int i=0;i<E2V.rows();i++){
            RowVector4d Xi=CurrX.segment(4*E2V(i,0),4).transpose();
            RowVector4d Xj=CurrX.segment(4*E2V(i,1),4).transpose();
            MobEdgeVectors.segment(4*i,4)=QMult1(QMult1(QConj1(Xj),SourceVq.row(E2V(i,1))-SourceVq.row(E2V(i,0))),Xi).transpose();
        }
        
        EnergyVec=MobEdgeVectors-ExtEdgeVectors;
        VectorXd CloseVec=(CurrSolution-InitSolution);
        
        VectorXd PosVec(3*ConstIndices.size());
        for (int i=0;i<ConstIndices.size();i++)
            PosVec.segment(3*i,3)=CurrLocations.segment(3*ConstIndices(i),3)-InterpConstPoses.segment(3*i,3);
        
        VectorXd BasicConstVec(EnergyVec.size()+CloseVec.size()+PosVec.size());
        BasicConstVec<<CloseFactor*CloseVec,SmoothFactor*EnergyVec, PosVec;
        
        //cout<<CloseVec.size()<<","<<EnergyVec.size()<<","<<PosVec.size()<<endl;
        
        return BasicConstVec;
        
    }
    
    SparseMatrix<double> EvalGradient(const VectorXd& CurrSolution){
        
        VectorXd CurrX=CurrSolution.head(4*SourceVq.rows());
        VectorXd CurrLocations=CurrSolution.segment(4*SourceVq.rows(),3*SourceVq.rows());
        VectorXd ExtEdgeVectors(4*E2V.rows()); ExtEdgeVectors.setZero();
        VectorXd CurrWij=d0*CurrLocations;
        for (int i=0;i<E2V.rows();i++){
            RowVector4d ExtCurrWij; ExtCurrWij<<0,CurrWij.segment(3*i,3).transpose();
            ExtEdgeVectors.segment(4*i,4)=QMult1(InterpEdgeDeviations.segment(4*i,4).transpose(),ExtCurrWij).transpose();
        }
        
        vector<Triplet<double> > gMobEdgesTris, gCurrEdgesTris;
        for (int i=0;i<E2V.rows();i++){
            int MatPosi=4*E2V(i,0);
            int MatPosj=4*E2V(i,1);
            int CurrRow=4*i;
            
            RowVector4d Xi=CurrX.segment(4*E2V(i,0),4).transpose();
            RowVector4d Xj=CurrX.segment(4*E2V(i,1),4).transpose();
            RowVector4d Cij=InterpEdgeDeviations.segment(4*i,4).transpose();
            
            RowVector4d ResttoRight=QMult1(SourceVq.row(E2V(i,1))-SourceVq.row(E2V(i,0)),Xi);
            RowVector4d ResttoLeft=QMult1(QConj1(Xj),SourceVq.row(E2V(i,1))-SourceVq.row(E2V(i,0)));
            
            //Derivative by Gj, using rest to right as "e":
            //[diff(qce',rq) diff(qce',vqx), diff(qce', vqy), diff(qce', vqz)]=
            
            
            //[  re,  vex,  vey,  vez]
            //[ vex,  -re, -vez,  vey]
            //[ vey,  vez,  -re, -vex]
            //[ vez, -vey,  vex,  -re]
            gMobEdgesTris.push_back(Triplet<double>(4*i,MatPosj,  ResttoRight(0)));
            gMobEdgesTris.push_back(Triplet<double>(4*i,MatPosj+1,ResttoRight(1)));
            gMobEdgesTris.push_back(Triplet<double>(4*i,MatPosj+2,ResttoRight(2)));
            gMobEdgesTris.push_back(Triplet<double>(4*i,MatPosj+3,ResttoRight(3)));
            
            gMobEdgesTris.push_back(Triplet<double>(4*i+1,MatPosj,   ResttoRight(1)));
            gMobEdgesTris.push_back(Triplet<double>(4*i+1,MatPosj+1,-ResttoRight(0)));
            gMobEdgesTris.push_back(Triplet<double>(4*i+1,MatPosj+2,-ResttoRight(3)));
            gMobEdgesTris.push_back(Triplet<double>(4*i+1,MatPosj+3, ResttoRight(2)));
            
            gMobEdgesTris.push_back(Triplet<double>(4*i+2,MatPosj,   ResttoRight(2)));
            gMobEdgesTris.push_back(Triplet<double>(4*i+2,MatPosj+1, ResttoRight(3)));
            gMobEdgesTris.push_back(Triplet<double>(4*i+2,MatPosj+2,-ResttoRight(0)));
            gMobEdgesTris.push_back(Triplet<double>(4*i+2,MatPosj+3,-ResttoRight(1)));
            
            gMobEdgesTris.push_back(Triplet<double>(4*i+3,MatPosj,   ResttoRight(3)));
            gMobEdgesTris.push_back(Triplet<double>(4*i+3,MatPosj+1,-ResttoRight(2)));
            gMobEdgesTris.push_back(Triplet<double>(4*i+3,MatPosj+2, ResttoRight(1)));
            gMobEdgesTris.push_back(Triplet<double>(4*i+3,MatPosj+3,-ResttoRight(0)));
            
            //Derivative by Gi, using rest to left as "e":
            //[  re, -vex, -vey, -vez]
            //[ vex,   re, -vez,  vey]
            //[ vey,  vez,   re, -vex]
            //[ vez, -vey,  vex,   re]
            
            gMobEdgesTris.push_back(Triplet<double>(4*i,MatPosi,  ResttoLeft(0)));
            gMobEdgesTris.push_back(Triplet<double>(4*i,MatPosi+1,-ResttoLeft(1)));
            gMobEdgesTris.push_back(Triplet<double>(4*i,MatPosi+2,-ResttoLeft(2)));
            gMobEdgesTris.push_back(Triplet<double>(4*i,MatPosi+3,-ResttoLeft(3)));
            
            gMobEdgesTris.push_back(Triplet<double>(4*i+1,MatPosi,   ResttoLeft(1)));
            gMobEdgesTris.push_back(Triplet<double>(4*i+1,MatPosi+1, ResttoLeft(0)));
            gMobEdgesTris.push_back(Triplet<double>(4*i+1,MatPosi+2,-ResttoLeft(3)));
            gMobEdgesTris.push_back(Triplet<double>(4*i+1,MatPosi+3, ResttoLeft(2)));
            
            gMobEdgesTris.push_back(Triplet<double>(4*i+2,MatPosi,   ResttoLeft(2)));
            gMobEdgesTris.push_back(Triplet<double>(4*i+2,MatPosi+1, ResttoLeft(3)));
            gMobEdgesTris.push_back(Triplet<double>(4*i+2,MatPosi+2, ResttoLeft(0)));
            gMobEdgesTris.push_back(Triplet<double>(4*i+2,MatPosi+3,-ResttoLeft(1)));
            
            gMobEdgesTris.push_back(Triplet<double>(4*i+3,MatPosi,   ResttoLeft(3)));
            gMobEdgesTris.push_back(Triplet<double>(4*i+3,MatPosi+1,-ResttoLeft(2)));
            gMobEdgesTris.push_back(Triplet<double>(4*i+3,MatPosi+2, ResttoLeft(1)));
            gMobEdgesTris.push_back(Triplet<double>(4*i+3,MatPosi+3, ResttoLeft(0)));
            
            //Derivative with relation to positions
            //Derivative by Wij (with proper sign), using rest to left (Cij) as "e":
            //[  re, -vex, -vey, -vez]
            //[ vex,   re, -vez,  vey]
            //[ vey,  vez,   re, -vex]
            //[ vez, -vey,  vex,   re]
            //gCurrEdgesTris.push_back(Triplet<double>(4*i,CurrX.size()+4*E2V(i,1),  Cij(0)));
            gCurrEdgesTris.push_back(Triplet<double>(4*i,CurrX.size()+3*E2V(i,1),  -Cij(1)));
            gCurrEdgesTris.push_back(Triplet<double>(4*i,CurrX.size()+3*E2V(i,1)+1,-Cij(2)));
            gCurrEdgesTris.push_back(Triplet<double>(4*i,CurrX.size()+3*E2V(i,1)+2,-Cij(3)));
            
            //gCurrEdgesTris.push_back(Triplet<double>(4*i+1,CurrX.size()+4*E2V(i,1),   Cij(1)));
            gCurrEdgesTris.push_back(Triplet<double>(4*i+1,CurrX.size()+3*E2V(i,1), Cij(0)));
            gCurrEdgesTris.push_back(Triplet<double>(4*i+1,CurrX.size()+3*E2V(i,1)+1,-Cij(3)));
            gCurrEdgesTris.push_back(Triplet<double>(4*i+1,CurrX.size()+3*E2V(i,1)+2, Cij(2)));
            
            //gCurrEdgesTris.push_back(Triplet<double>(4*i+2,CurrX.size()+4*E2V(i,1),   Cij(2)));
            gCurrEdgesTris.push_back(Triplet<double>(4*i+2,CurrX.size()+3*E2V(i,1), Cij(3)));
            gCurrEdgesTris.push_back(Triplet<double>(4*i+2,CurrX.size()+3*E2V(i,1)+1, Cij(0)));
            gCurrEdgesTris.push_back(Triplet<double>(4*i+2,CurrX.size()+3*E2V(i,1)+2,-Cij(1)));
            
            //gCurrEdgesTris.push_back(Triplet<double>(4*i+3,CurrX.size()+4*E2V(i,1),   Cij(3)));
            gCurrEdgesTris.push_back(Triplet<double>(4*i+3,CurrX.size()+3*E2V(i,1),-Cij(2)));
            gCurrEdgesTris.push_back(Triplet<double>(4*i+3,CurrX.size()+3*E2V(i,1)+1, Cij(1)));
            gCurrEdgesTris.push_back(Triplet<double>(4*i+3,CurrX.size()+3*E2V(i,1)+2,Cij(0)));
            
            
            //gCurrEdgesTris.push_back(Triplet<double>(4*i,CurrX.size()+4*E2V(i,0),  -Cij(0)));
            gCurrEdgesTris.push_back(Triplet<double>(4*i,CurrX.size()+3*E2V(i,0),Cij(1)));
            gCurrEdgesTris.push_back(Triplet<double>(4*i,CurrX.size()+3*E2V(i,0)+1,Cij(2)));
            gCurrEdgesTris.push_back(Triplet<double>(4*i,CurrX.size()+3*E2V(i,0)+2,Cij(3)));
            
            //gCurrEdgesTris.push_back(Triplet<double>(4*i+1,CurrX.size()+4*E2V(i,0),   -Cij(1)));
            gCurrEdgesTris.push_back(Triplet<double>(4*i+1,CurrX.size()+3*E2V(i,0), -Cij(0)));
            gCurrEdgesTris.push_back(Triplet<double>(4*i+1,CurrX.size()+3*E2V(i,0)+1,Cij(3)));
            gCurrEdgesTris.push_back(Triplet<double>(4*i+1,CurrX.size()+3*E2V(i,0)+2, -Cij(2)));
            
            //gCurrEdgesTris.push_back(Triplet<double>(4*i+2,CurrX.size()+4*E2V(i,0),   -Cij(2)));
            gCurrEdgesTris.push_back(Triplet<double>(4*i+2,CurrX.size()+3*E2V(i,0), -Cij(3)));
            gCurrEdgesTris.push_back(Triplet<double>(4*i+2,CurrX.size()+3*E2V(i,0)+1, -Cij(0)));
            gCurrEdgesTris.push_back(Triplet<double>(4*i+2,CurrX.size()+3*E2V(i,0)+2,Cij(1)));
            
            //gCurrEdgesTris.push_back(Triplet<double>(4*i+3,CurrX.size()+4*E2V(i,0),   -Cij(3)));
            gCurrEdgesTris.push_back(Triplet<double>(4*i+3,CurrX.size()+3*E2V(i,0),Cij(2)));
            gCurrEdgesTris.push_back(Triplet<double>(4*i+3,CurrX.size()+3*E2V(i,0)+1, -Cij(1)));
            gCurrEdgesTris.push_back(Triplet<double>(4*i+3,CurrX.size()+3*E2V(i,0)+2,-Cij(0)));
        }
        
        SparseMatrix<double> gMobEdges(4*E2V.rows(),CurrSolution.size());
        gMobEdges.setFromTriplets(gMobEdgesTris.begin(), gMobEdgesTris.end());
        
        SparseMatrix<double> gCurrEdges(4*E2V.rows(), CurrSolution.size());
        gCurrEdges.setFromTriplets(gCurrEdgesTris.begin(),gCurrEdgesTris.end());
        SparseMatrix<double> gAMAPEnergy=gMobEdges-gCurrEdges;
        
        SparseMatrix<double> gPos(3*ConstIndices.size(), CurrSolution.size());
        vector<Triplet<double> > gPosTris;
        for (int i=0;i<ConstIndices.size();i++)
            for (int j=0;j<3;j++)
                gPosTris.push_back(Triplet<double>(3*i+j, 4*SourceVq.rows()+3*ConstIndices(i)+j,1.0));
        
        gPos.setFromTriplets(gPosTris.begin(), gPosTris.end());
        
        SparseMatrix<double> BasicGradient=SparseBlock<double>(SparseBlock<double>(CloseFactor*IDMat, SmoothFactor*gAMAPEnergy,false),gPos,false);
        
        
        return BasicGradient;
        
    }
    
    void Reformulate(int CurrIter, int MaxIterations, const VectorXd& CurrSolution)
    {
       InitSolution=CurrSolution;
       SmoothFactor*=0.5;
    }

    double GetConstError(const VectorXd& CurrSolution)
    {
        VectorXd CurrLocations=CurrSolution.segment(4*SourceVq.rows(),3*SourceVq.rows());
        VectorXd PosVec(3*ConstIndices.size());
        for (int i=0;i<ConstIndices.size();i++)
            PosVec.segment(3*i,3)=CurrLocations.segment(3*ConstIndices(i),3)-InterpConstPoses.segment(3*i,3);
        
        return PosVec.lpNorm<Infinity>();

    }

    bool isTerminate(const VectorXd& CurrSolution){
       //VectorXcd CurrSolution(CurrSolutionReal.size()/2);
       //CurrSolution.array().real()=CurrSolutionReal.head(CurrSolutionReal.size()/2);
       //CurrSolution.array().imag()=CurrSolutionReal.tail(CurrSolutionReal.size()/2);
       
       VectorXd CurrLocations=CurrSolution.segment(4*SourceVq.rows(),3*SourceVq.rows());
       VectorXd InitLocations=InitSolution.segment(4*SourceVq.rows(),3*SourceVq.rows());
       bool NoVerticesMoved=((CurrLocations-InitLocations).lpNorm<Infinity>()<ConstTolerance);
       
       bool WithinConstTolerance=(GetConstError(CurrSolution) < ConstTolerance);
       return (WithinConstTolerance && NoVerticesMoved);
       
    }

};


/*class MobiusBoundedInterpolatationTraits3D{
    
    MatrixXi E2V, E2F, E2Fi, F, F2E;
    MatrixXd F2ESigns;
    VectorXi InnerEdges, ConstIndices, VarIndices;
    SparseMatrix<double> d0Var;
    VectorXd d0Rhs;
    vector<vector<int> > AdjFaces;
    vector<vector<int> > AdjFacesi;
    
    MatrixXd OrigVq, DeformVq;
    
    VectorXd InitSolution;
    
    VectorXd AbsOrigCR, AbsDeformCR, AbsInterpCR;
    VectorXd RealOrigCR, RealDeformCR, RealInterpCR;
    
    VectorXi VarLocIndices;

    SparseMatrix<double> IDMat;
    SparseMatrix<double> NormSumMatrix;
    SparseMatrix<double> CornerDiffLCRSmall, CornerDiffSmall, CornerDiffMat;
    
    double CloseFactor;
    double SmoothFactor;
    
    
    vector<pair<int, int> > CornerPairs;
    
    void Initialize(const MatrixXd& inOrigVq, const MatrixXd& inDeformVq, const MatrixXi& inF, MatrixXi& inF2E, MatrixXd& inF2ESigns, const MatrixXi& inE2V, const VectorXi& inInnerEdges, const MatrixXi& inE2F, const MatrixXi& inE2Fi, vector<vector<int> >& inAdjFaces, vector<vector<int> >& inAdjFacesi, const VectorXi& inConstIndices,  VectorXi& inVarIndices, SparseMatrix<double>& ind0var, VectorXd& ind0Rhs, VectorXd inAbsOrigCR, VectorXd& inAbsDeformCR, VectorXd& inRealOrigCR, VectorXd& inRealDeformCR, vector<pair<int,int> >& inCornerPairs,const double t, const VectorXd& inInitSolution){
        F=inF; E2V=inE2V; F2E=inF2E; F2ESigns=inF2ESigns; E2F=inE2F; InnerEdges=inInnerEdges; E2F=inE2F; E2Fi=inE2Fi; ConstIndices=inConstIndices;
        
        
        OrigVq=inOrigVq;
        DeformVq=inDeformVq;
        InitSolution=inInitSolution;
        AbsInterpCR.resize(InnerEdges.size());
        RealInterpCR.resize(InnerEdges.size());
        AbsOrigCR=inAbsOrigCR;
        AbsDeformCR=inAbsDeformCR;
        RealOrigCR=inRealOrigCR;
        RealDeformCR=inRealDeformCR;
        
        VarIndices=inVarIndices;
        d0Var=ind0var; d0Rhs=ind0Rhs;
        CornerPairs=inCornerPairs;
        AdjFaces=inAdjFaces;
        AdjFacesi=inAdjFacesi;
        
        for (int i=0;i<InnerEdges.rows();i++){
            AbsInterpCR(i)=exp(log(AbsOrigCR(i))*(1-t)+log(AbsDeformCR(i))*t);
            RealInterpCR(i)=exp(log(RealOrigCR(i))*(1-t)+log(RealDeformCR(i))*t);
        }
        
        //assuming CornerPairs is two per edge on the order of inneredges!!!
        
        vector<Triplet<double> > CornerDiffLCRTris, CornerDiffTris;
        //creating smoothness matrix for G values
        for (int i=0;i<InnerEdges.size();i++){
            CornerDiffLCRTris.push_back(Triplet<double>(2*i, CornerPairs[i].first,1.0));
            CornerDiffLCRTris.push_back(Triplet<double>(2*i+1, CornerPairs[i].second,-AbsInterpCR(i)));
            
            CornerDiffTris.push_back(Triplet<double>(2*i, CornerPairs[i].first,1.0));
            CornerDiffTris.push_back(Triplet<double>(2*i+1, CornerPairs[i].second,-1.0));
        }
        
        CornerDiffLCRSmall.resize(2*InnerEdges.size(), 3*F.rows());
        CornerDiffLCRSmall.setFromTriplets(CornerDiffLCRTris.begin(), CornerDiffLCRTris.end());
        
        CornerDiffSmall.resize(2*InnerEdges.size(), 3*F.rows());
        CornerDiffSmall.setFromTriplets(CornerDiffTris.begin(), CornerDiffTris.end());
        
        SparseMatrix<double> eye4;
        igl::speye(4, eye4);
        SparseMatrix<double> CornerDiffOnlyG=Eigen::kroneckerProduct(CornerDiffSmall,eye4);
        SparseMatrix<double> ZeroMat(CornerDiffOnlyG.rows(), 3*VarIndices.size());
        CornerDiffMat=SparseBlock(CornerDiffOnlyG, ZeroMat, true);

        
        VarLocIndices.resize(12*F.rows()+3*VarIndices.size());
        for (int i=0;i<12*F.rows();i++)
            VarLocIndices(i)=i;
        for (int i=0;i<VarIndices.size();i++)
            for (int j=0;j<3;j++)
                VarLocIndices(12*F.rows()+3*i+j)=12*F.rows()+3*VarIndices(i)+j;
        
        igl::speye(12*F.rows()+3*VarIndices.size(), IDMat);
        
        //constructing sum matrix
        vector<Triplet<double> > NormSumTris;
        for (int i=0;i<3*F.rows();i++)
            for(int j=0;j<4;j++)
                NormSumTris.push_back(Triplet<double>(i,4*i+j,1.0));
        
        NormSumMatrix.resize(3*F.rows(),12*F.rows()+3*VarIndices.size());
        NormSumMatrix.setFromTriplets(NormSumTris.begin(), NormSumTris.end());
        
        CloseFactor=10e-6;
        SmoothFactor=100;
        
    }
    
    VectorXd EvalEnergy(const VectorXd& CurrSolution){
    
        VectorXd CurrGInv=CurrSolution.head(12*F.rows());
        VectorXd CurrLocations=CurrSolution.tail(3*VarIndices.size());
        
        VectorXd CompVec(12*F.rows());
        VectorXd EdgeVectors=d0Var*CurrLocations-d0Rhs;
        for (int i=0;i<F.rows();i++){
            for (int j=0;j<3;j++){
                RowVector4d GiInv=CurrGInv.segment(4*(3*i+j),4).transpose();
                RowVector4d GjInv=CurrGInv.segment(4*(3*i+(j+1)%3),4).transpose();
                //cout<<"GiInv: "<<GiInv<<endl;
                //cout<<"GjInv: "<<GjInv<<endl;
                RowVector4d CurrEdgeVector=QMult(QMult(QConj(GjInv),OrigVq.row(F(i,(j+1)%3))-OrigVq.row(F(i,j))),GiInv);
                //cout<<"Transformed Edge Vector for ("<<i<<","<<j<<") is: ("<<CurrEdgeVector<<")"<<endl;
                //cout<<"Edge Vector is: ("<<EdgeVectors.segment(3*F2E(i,j),3)*F2ESigns(i,j)<<")"<<endl;
                CurrEdgeVector.tail(3)-=EdgeVectors.segment(3*F2E(i,j),3)*F2ESigns(i,j);
                CompVec.segment(4*(3*i+j),4)=CurrEdgeVector.transpose();
            }
        }
        
        VectorXd EnergyVec=SmoothFactor*CornerDiffMat*CurrSolution;
        
        VectorXd CloseVec=CloseFactor*(CurrSolution-InitSolution);
        VectorXd LCRVec=CornerDiffLCRSmall*NormSumMatrix*(CurrSolution.cwiseAbs2());
        
        /*VectorXd IAPVec(3*CornerPairs.size());
        for (int i=0;i<CornerPairs.size();i++){
            RowVector4d Gi=CurrGInv.segment(4*CornerPairs[i].first,4);
            RowVector4d Gj=CurrGInv.segment(4*CornerPairs[i].second,4);
            RowVector4d GiGj=QMult(QConj(Gj),Gi);
            IAPVec.segment(3*i,3)=GiGj.tail(3);
        }
        
        IAPVec=IAPVec*IAPConstFactor;*/ /*
        
        VectorXd ConstVec(CloseVec.size()+EnergyVec.size()+CompVec.size()+LCRVec.size());
        
        ConstVec<<CloseVec,EnergyVec, CompVec, LCRVec;//, RCRVec;
        //ConstVec=CompVec;
        
        
        return ConstVec;
    }

    SparseMatrix<double> EvalGradient(const VectorXd& CurrSolution){
        
        VectorXd CurrGInv=CurrSolution.head(12*F.rows());
        VectorXd CurrLocations=CurrSolution.tail(3*VarIndices.size());
        
        vector<Triplet<double> > gCompTris;
        for (int i=0;i<F.rows();i++){
            for (int j=0;j<3;j++){
                int MatPosi=4*(3*i+j);
                int MatPosj=4*(3*i+(j+1)%3);
                int CurrRow=4*(3*i+j);
                
                RowVector4d GiInv=CurrGInv.segment(4*(3*i+j),4).transpose();
                RowVector4d GjInv=CurrGInv.segment(4*(3*i+(j+1)%3),4).transpose();
                
                RowVector4d ResttoRight=QMult(OrigVq.row(F(i,(j+1)%3))-OrigVq.row(F(i,j)),GiInv);
                RowVector4d ResttoLeft=QMult(QConj(GjInv),OrigVq.row(F(i,(j+1)%3))-OrigVq.row(F(i,j)));
                
                
                //Derivative by Gj, using rest to right as "e":
                //[diff(qce',rq) diff(qce',vqx), diff(qce', vqy), diff(qce', vqz)]=
                
                
                //[  re,  vex,  vey,  vez]
                //[ vex,  -re, -vez,  vey]
                //[ vey,  vez,  -re, -vex]
                //[ vez, -vey,  vex,  -re]
                gCompTris.push_back(Triplet<double>(CurrRow,MatPosj,  ResttoRight(0)));
                gCompTris.push_back(Triplet<double>(CurrRow,MatPosj+1,ResttoRight(1)));
                gCompTris.push_back(Triplet<double>(CurrRow,MatPosj+2,ResttoRight(2)));
                gCompTris.push_back(Triplet<double>(CurrRow,MatPosj+3,ResttoRight(3)));
                
                gCompTris.push_back(Triplet<double>(CurrRow+1,MatPosj,   ResttoRight(1)));
                gCompTris.push_back(Triplet<double>(CurrRow+1,MatPosj+1,-ResttoRight(0)));
                gCompTris.push_back(Triplet<double>(CurrRow+1,MatPosj+2,-ResttoRight(3)));
                gCompTris.push_back(Triplet<double>(CurrRow+1,MatPosj+3, ResttoRight(2)));
                
                gCompTris.push_back(Triplet<double>(CurrRow+2,MatPosj,   ResttoRight(2)));
                gCompTris.push_back(Triplet<double>(CurrRow+2,MatPosj+1, ResttoRight(3)));
                gCompTris.push_back(Triplet<double>(CurrRow+2,MatPosj+2,-ResttoRight(0)));
                gCompTris.push_back(Triplet<double>(CurrRow+2,MatPosj+3,-ResttoRight(1)));
                
                gCompTris.push_back(Triplet<double>(CurrRow+3,MatPosj,   ResttoRight(3)));
                gCompTris.push_back(Triplet<double>(CurrRow+3,MatPosj+1,-ResttoRight(2)));
                gCompTris.push_back(Triplet<double>(CurrRow+3,MatPosj+2, ResttoRight(1)));
                gCompTris.push_back(Triplet<double>(CurrRow+3,MatPosj+3,-ResttoRight(0)));
                
                //Derivative by Gi, using rest to left as "e":
                //[  re, -vex, -vey, -vez]
                //[ vex,   re, -vez,  vey]
                //[ vey,  vez,   re, -vex]
                //[ vez, -vey,  vex,   re]
                
                gCompTris.push_back(Triplet<double>(CurrRow,MatPosi,  ResttoLeft(0)));
                gCompTris.push_back(Triplet<double>(CurrRow,MatPosi+1,-ResttoLeft(1)));
                gCompTris.push_back(Triplet<double>(CurrRow,MatPosi+2,-ResttoLeft(2)));
                gCompTris.push_back(Triplet<double>(CurrRow,MatPosi+3,-ResttoLeft(3)));
                
                gCompTris.push_back(Triplet<double>(CurrRow+1,MatPosi,   ResttoLeft(1)));
                gCompTris.push_back(Triplet<double>(CurrRow+1,MatPosi+1, ResttoLeft(0)));
                gCompTris.push_back(Triplet<double>(CurrRow+1,MatPosi+2,-ResttoLeft(3)));
                gCompTris.push_back(Triplet<double>(CurrRow+1,MatPosi+3, ResttoLeft(2)));
                
                gCompTris.push_back(Triplet<double>(CurrRow+2,MatPosi,   ResttoLeft(2)));
                gCompTris.push_back(Triplet<double>(CurrRow+2,MatPosi+1, ResttoLeft(3)));
                gCompTris.push_back(Triplet<double>(CurrRow+2,MatPosi+2, ResttoLeft(0)));
                gCompTris.push_back(Triplet<double>(CurrRow+2,MatPosi+3,-ResttoLeft(1)));
                
                gCompTris.push_back(Triplet<double>(CurrRow+3,MatPosi,   ResttoLeft(3)));
                gCompTris.push_back(Triplet<double>(CurrRow+3,MatPosi+1,-ResttoLeft(2)));
                gCompTris.push_back(Triplet<double>(CurrRow+3,MatPosi+2, ResttoLeft(1)));
                gCompTris.push_back(Triplet<double>(CurrRow+3,MatPosi+3, ResttoLeft(0)));
                
                //Derivative with relation to positions
                gCompTris.push_back(Triplet<double>(CurrRow+1,CurrGInv.size()+3*F(i,(j+1)%3),  -1.0));
                gCompTris.push_back(Triplet<double>(CurrRow+2,CurrGInv.size()+3*F(i,(j+1)%3)+1,-1.0));
                gCompTris.push_back(Triplet<double>(CurrRow+3,CurrGInv.size()+3*F(i,(j+1)%3)+2,-1.0));
                
                gCompTris.push_back(Triplet<double>(CurrRow+1,CurrGInv.size()+3*F(i,j),  1.0));
                gCompTris.push_back(Triplet<double>(CurrRow+2,CurrGInv.size()+3*F(i,j)+1,1.0));
                gCompTris.push_back(Triplet<double>(CurrRow+3,CurrGInv.size()+3*F(i,j)+2,1.0));
                
            }
        }
        
        SparseMatrix<double> gCompFull(12*F.rows(), 3*OrigVq.rows()+CurrGInv.size());
        SparseMatrix<double> gComp(12*F.rows(), 3*VarIndices.size()+CurrGInv.size());
        gCompFull.setFromTriplets(gCompTris.begin(),gCompTris.end());
        
        igl::slice(gCompFull,VarLocIndices, 2,gComp);
        
        SparseMatrix<double> gLCR=CornerDiffLCRSmall*NormSumMatrix*SparseDiagonalMatrix<double>(2*CurrSolution);
        
        SparseMatrix<double> Gradient=SparseBlock<double>(SparseBlock<double>(CloseFactor*IDMat, SmoothFactor*CornerDiffMat, false), SparseBlock<double>(gComp,gLCR, false),false);
        
        return Gradient;
        
    }

    
    void Reformulate(int CurrIter, int MaxIterations, VectorXd CurrSolution){
    }
    
    double GetConstError(VectorXd CurrError)
    {
        return CurrError.lpNorm<Infinity>();
    }

};*/


#endif
