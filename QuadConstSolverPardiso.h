//
//  QuadConstSolver.h
//  testigl
//
//  Created by Amir Vaxman on 28/09/14.
//  Copyright (c) 2014 Amir Vaxman. All rights reserved.
//

#ifndef testigl_QuadConstSolver_pardiso_h
#define testigl_QuadConstSolver_pardiso_h

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include "PardisoSolver.h"
#include "AuxSparse.h"
#include <iostream>
#include <igl/sortrows.h>

using namespace igl;

template<class QuadConstSolverTraits>
class QuadConstSolverPardiso{
private:
    QuadConstSolverTraits* SolverTraits;
    PardisoSolver ps;
public:
    
    VectorXd ConvErrors;
    
    QuadConstSolverPardiso(){};
    
    //returns M^T*M by (I,J,S) representation
    //prerequisite: iI are sorted by rows (not necessary columns)
    void MultiplyAdjointMatrix(const VectorXi& iI, const VectorXi& iJ, const VectorXd& iS, VectorXi& oI, VectorXi& oJ, VectorXd& oS)
    {
        int CurrTri=0;
        
        vector<int> oIlist;
        vector<int> oJlist;
        vector<double> oSlist;
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
    void MultiplyAdjointVector(const VectorXi& iI, const VectorXi& iJ, const VectorXd& iS, const VectorXd& iVec, VectorXd& oVec)
    {
        int MaxCol=iJ.maxCoeff()+1;
        
        oVec.resize(MaxCol); oVec.setZero();
        for (int i=0;i<iI.size();i++)
            oVec(iJ(i))+=iS(i)*iVec(iI(i));
        
    }
    
    
    void Initialize(QuadConstSolverTraits* st){
        SolverTraits=st;
        
        
        //analysing pattern
        VectorXi I,J;
        VectorXd S;
        
        //st->GradRows.resize(4);
        //st->GradRows<<0,0,1,1;
        //st->GradCols.resize(4);
        //st->GradCols<<0,1,0,1;
        MultiplyAdjointMatrix(st->GradRows, st->GradCols,VectorXd::Ones(st->GradRows.size()),I,J,S);
        //testing multiplyadjoint
        /*SparseMatrix<double> M(st->GradRows.maxCoeff()+1,st->GradCols.maxCoeff()+1);
        
        //cout<<"Our I"<<I<<endl;
        //cout<<"Our J"<<J<<endl;
        //cout<<"Our S"<<S<<endl;
        
        vector<Triplet<double> > Tris;
        for (int i=0;i<st->GradRows.size();i++)
            Tris.push_back(Triplet<double>(st->GradRows(i), st->GradCols(i), 1.0));
        
        
        M.setFromTriplets(Tris.begin(), Tris.end());
        Tris.clear();
        
        SparseMatrix<double> MtM1=M.transpose()*M;
        SparseMatrix<double> MtM2(MtM1.rows(), MtM1.cols());
        for (int i=0;i<I.size();i++)
            Tris.push_back(Triplet<double>(I(i), J(i), S(i)));
        MtM2.setFromTriplets(Tris.begin(), Tris.end());
        SparseMatrix<double> DiffMat=MtM1-MtM2;
        for (int k=0; k<DiffMat.outerSize(); ++k)
            for (SparseMatrix<double>::InnerIterator it(DiffMat,k); it; ++it)
                if ((abs(it.value())>10e-6)&&(it.row()<=it.col())){
                    cout<<"Discrepancy at ("<<it.row()<<","<<it.col()<<"): "<<it.value()<<endl;
                    cout<<"MtM Our Values:"<<MtM2.coeffRef(it.row(),it.col())<<endl;
                    cout<<"MtM Evaluated:"<<MtM1.coeffRef(it.row(),it.col())<<endl;
                }*/
        
        ps.set_pattern(I,J,S);
        ps.analyze_pattern();
    }
    
    VectorXd Solve(const VectorXd& InitSolution, double Inith, int MaxIterations){
        
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
            
            //testing
            /*SparseMatrix<double> M(SolverTraits->GradRows.maxCoeff()+1,SolverTraits->GradCols.maxCoeff()+1);
            
            //cout<<"Our I"<<I<<endl;
            //cout<<"Our J"<<J<<endl;
            //cout<<"Our S"<<S<<endl;
            
            vector<Triplet<double> > Tris;
            for (int i=0;i<SolverTraits->GradRows.size();i++)
                Tris.push_back(Triplet<double>(SolverTraits->GradRows(i), SolverTraits->GradCols(i), SolverTraits->GradValues(i)));
            
            
            M.setFromTriplets(Tris.begin(), Tris.end());
            Tris.clear();
            
            SparseMatrix<double> MtM1=M.transpose()*M;
            SparseMatrix<double> MtM2(MtM1.rows(), MtM1.cols());
            for (int i=0;i<MatRows.size();i++)
                Tris.push_back(Triplet<double>(MatRows(i), MatCols(i), MatValues(i)));
            MtM2.setFromTriplets(Tris.begin(), Tris.end());
            SparseMatrix<double> DiffMat=MtM1-MtM2;
            for (int k=0; k<DiffMat.outerSize(); ++k)
                for (SparseMatrix<double>::InnerIterator it(DiffMat,k); it; ++it)
                    if ((abs(it.value())>10e-6)&&(it.row()<=it.col())){
                        cout<<"Discrepancy at ("<<it.row()<<","<<it.col()<<"): "<<it.value()<<endl;
                        cout<<"MtM Our Values:"<<MtM2.coeffRef(it.row(),it.col())<<endl;
                        cout<<"MtM Evaluated:"<<MtM1.coeffRef(it.row(),it.col())<<endl;
                    }
            
            //testing vector multiplication
            VectorXd CompVector=M.transpose()*(-SolverTraits->TotalVec);
            cout<<"Maximum vector multiplication difference: "<<(CompVector-Rhs).lpNorm<Infinity>()<<endl;*/

            
            //solving to get the direction
            
            ps.update_a(MatValues);
            if(!ps.factorize()) {
                // decomposition failed
                int kaka=8;
                cout<<"Solver Failed to factorize! "<<endl;
            }
    
            ps.solve(Rhs,Direction);
            
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
