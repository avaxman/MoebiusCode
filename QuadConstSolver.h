//
//  QuadConstSolver.h
//  testigl
//
//  Created by Amir Vaxman on 28/09/14.
//  Copyright (c) 2014 Amir Vaxman. All rights reserved.
//

#ifndef testigl_QuadConstSolver_h
#define testigl_QuadConstSolver_h

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <hedra/EigenSolverWrapper.h>
#include <iostream>
#include <igl/sortrows.h>

template<class QuadConstSolverTraits>
class QuadConstSolver{
private:
    QuadConstSolverTraits* SolverTraits;
    hedra::optimization::EigenSolverWrapper<Eigen::SimplicialLDLT<Eigen::SparseMatrix<double> > > LinearSolver;
public:
    
    Eigen::VectorXd ConvErrors;
    
    QuadConstSolver(){};
    
    //returns M^T*M by (I,J,S) representation
    //prerequisite: iI are sorted by rows (not necessary columns)
    void MultiplyAdjointMatrix(const Eigen::VectorXi& iI,
                               const Eigen::VectorXi& iJ,
                               const Eigen::VectorXd& iS,
                               Eigen::VectorXi& oI,
                               Eigen::VectorXi& oJ,
                               Eigen::VectorXd& oS)
    {
        int CurrTri=0;
        
        std::vector<int> oIlist;
        std::vector<int> oJlist;
        std::vector<double> oSlist;
        do{
            int CurrRow=iI(CurrTri);
            int NumCurrTris=0;
            while ((CurrTri+NumCurrTris<iI.size())&&(iI(CurrTri+NumCurrTris)==CurrRow))
                NumCurrTris++;
            
            for (int i=CurrTri;i<CurrTri+NumCurrTris;i++){
                for (int j=CurrTri;j<CurrTri+NumCurrTris;j++){
                    if (iJ(j)>=iJ(i)){
                        oIlist.push_back(iJ(i));
                        oJlist.push_back(iJ(j));
                        oSlist.push_back(iS(j)*iS(i));
                    }
                }
            }
            CurrTri+=NumCurrTris;
        }while (CurrTri!=iI.size());
        
        oI.resize(oIlist.size());
        oJ.resize(oJlist.size());
        oS.resize(oSlist.size());
        
        for (int i=0;i<oIlist.size();i++){
            oI(i)=oIlist[i];
            oJ(i)=oJlist[i];
            oS(i)=oSlist[i];
        }
    }
    
    //returns M^t*ivec by (I,J,S) representation
    void MultiplyAdjointVector(const Eigen::VectorXi& iI,
                               const Eigen::VectorXi& iJ,
                               const Eigen::VectorXd& iS,
                               const Eigen::VectorXd& iVec,
                               Eigen::VectorXd& oVec)
    {
        int MaxCol=iJ.maxCoeff()+1;
        
        oVec.resize(MaxCol); oVec.setZero();
        for (int i=0;i<iI.size();i++)
            oVec(iJ(i))+=iS(i)*iVec(iI(i));
        
    }
    
    
    void Initialize(QuadConstSolverTraits* st){
        SolverTraits=st;
        
        
        //analysing pattern
        Eigen::VectorXi I,J;
        Eigen::VectorXd S;
        
        MultiplyAdjointMatrix(st->GradRows, st->GradCols,Eigen::VectorXd::Ones(st->GradRows.size()),I,J,S);
    
        LinearSolver.analyze(I,J);
    }
    
    Eigen::VectorXd Solve(const Eigen::VectorXd& InitSolution,
                   double Inith,
                   int MaxIterations){
        
        using namespace Eigen;
        using namespace std;
        
        vector<double> ConvErrorsList;
 
        VectorXd CurrSolution=InitSolution;
        int CurrIter=0;
        
        SolverTraits->UpdateEnergy(InitSolution);
        SolverTraits->UpdateGradient(InitSolution);
        
        cout<<"Initial Total Error: "<<(SolverTraits->TotalVec).template lpNorm<Infinity>()<<endl;
        cout<<"Initial Const Error: "<<(SolverTraits->ConstVec).template lpNorm<Infinity>()<<endl;
        
        double ConstError=(SolverTraits->ConstVec).template lpNorm<Infinity>();
        ConvErrorsList.push_back(ConstError);
    
        bool isTerminate;
        double PrevConstError=32767.0;
        VectorXi MatRows, MatCols;
        VectorXd MatValues;
        VectorXd Rhs;
        VectorXd Direction;
        do{
            MultiplyAdjointMatrix(SolverTraits->GradRows, SolverTraits->GradCols,SolverTraits->GradValues, MatRows, MatCols, MatValues);
            MultiplyAdjointVector(SolverTraits->GradRows, SolverTraits->GradCols, SolverTraits->GradValues, -SolverTraits->TotalVec, Rhs);
            
            //solving to get the direction
            if(!LinearSolver.factorize(MatValues)) {
                // decomposition failed
                cout<<"Solver Failed to factorize! "<<endl;
                return CurrSolution;
            }
    
            LinearSolver.solve(Rhs,Direction);
            
            //doing a line search
            double h=Inith;
            double PrevTotalError=(SolverTraits->TotalVec).template lpNorm<Infinity>();
            VectorXd CurrSolutionTest;
            double TestTotalError;
            do{
                CurrSolutionTest=CurrSolution+h*Direction;
                SolverTraits->UpdateEnergy(CurrSolutionTest);
                TestTotalError=(SolverTraits->TotalVec).template lpNorm<Infinity>();
                h*=0.5;
            }while ((TestTotalError>PrevTotalError)&&(h>10e-8));

            CurrSolution=CurrSolutionTest;
            PrevConstError=ConstError;
            
            //assuming UpdateEnergy is already taken within the line search iteration
            ConstError=(SolverTraits->ConstVec).template lpNorm<Infinity>();
            cout<<"Total Error:"<<TestTotalError<<endl;
            isTerminate=SolverTraits->isTerminate();
            if (isTerminate)
                cout<<"Terminating because Constraints are met!"<<endl;
            cout<<"Const Error "<<ConstError<<endl;
            if ((h<10e-8)&&(isTerminate))
                break;
            CurrIter++;
            SolverTraits->Reformulate(CurrIter, MaxIterations, CurrSolution, ConstError);
            ConvErrorsList.push_back(ConstError);
            SolverTraits->UpdateGradient(CurrSolution);
        }while ((CurrIter<MaxIterations)&&(!isTerminate));
        
        SolverTraits->UpdateEnergy(CurrSolution);
        double FinalTotalError=(SolverTraits->TotalVec).template lpNorm<Infinity>();
        double FinalConstError=(SolverTraits->ConstVec).template lpNorm<Infinity>();
        cout<<"Final Const Error:"<<FinalConstError<<endl;
        cout<<"Final Total Error:"<<FinalTotalError<<endl;
        cout<<"Number of Iterations and max iterations:"<<CurrIter<<","<<MaxIterations<<endl;
        ConvErrors.resize(ConvErrorsList.size());
        for (int i=0;i<ConvErrorsList.size();i++)
            ConvErrors(i)=ConvErrorsList[i];
        
        return CurrSolution;
        
    }
};


#endif
