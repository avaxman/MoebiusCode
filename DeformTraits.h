// This file is part of libhedra, a library for polyhedral mesh processing
//
// Copyright (C) 2016 Amir Vaxman <avaxman@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef HEDRA_MOEBIUS_EDGE_DEVIATION_2D_TRAITS_H
#define HEDRA_MOEBIUS_EDGE_DEVIATION_2D_TRAITS_H
#include <igl/igl_inline.h>
#include <Eigen/Core>
#include <string>
#include <vector>
#include <cstdio>
#include <set>


namespace hedra { namespace optimization {
    

    //This traits class implements the "unonstraintTraits" concept, for the energy and constraints of the complex system of deformation in [Vaxman et. al 2015] (see Tutorial).
    //TODO: fully integrate this as a constraint class
    
class Moebius2DEdgeDeviationTraits{
public:
    
    
    //concept requirements
    Eigen::VectorXi JRows, JCols;  //rows and column indices for the jacobian matrix
    Eigen::VectorXd JVals;         //values for the jacobian matrix.
    Eigen::VectorXd EVec;          //energy vector
    int xSize;                  //size of the solution

    typedef std::complex<double> Complex;
    
    Eigen::VectorXcd origVc;
    Eigen::MatrixXi EV, D, F;
    Eigen::VectorXi constIndices;//, VarIndices;
    Eigen::VectorXcd complexConstPoses;
    Eigen::VectorXcd initSolution;
    Eigen::VectorXcd finalPositions;  //final vertex positions
    Eigen::VectorXcd finalY;  //final vertex reciprocals
    Eigen::VectorXcd finalE;  //final edge deivations
    
    double smoothFactor;
    double closeFactor;
    double constTolerance;
    double rigidRatio;
    
    bool isExactMC;
    bool isExactIAP;
    
    
    //Complex variables
    Eigen::VectorXcd currSolution;
    Eigen::VectorXcd currMobRecips;
    Eigen::VectorXcd currLocations;
    Eigen::VectorXcd currEdgeDeviations;
    Eigen::VectorXcd currEdges;
    Eigen::VectorXcd AMAPVec;
    Eigen::VectorXcd mobVec;
    Eigen::VectorXcd posVec;
    Eigen::VectorXcd rigidVec;
    Eigen::VectorXcd closeVec;
    Eigen::VectorXcd deviationVec;
    
    //real variables
    Eigen::VectorXd MCVec;
    Eigen::VectorXd IAPVec;
    Eigen::VectorXcd FCR;
    
    Eigen::SparseMatrix<Complex> d0;
    
    Eigen::VectorXi complexJRows;
    Eigen::VectorXi complexJCols;
    Eigen::VectorXcd complexJVals;
    
    //into the complex values
    int AMAPTriOffset, AMAPRowOffset;
    int rigidTriOffset, RigidRowOffset;
    int posTriOffset, PosRowOffset;
    int mobTriOffset, mobRowOffset;
    int closeTriOffset, closeRowOffset;
    int deviationTriOffset, deviationRowOffset;
    int complexRowOffset, complexTriOffset;
    
    //into the real values
    int MCTriOffset, MCRowOffset;
    int IAPTriOffset, IAPRowOffset;
    
    Eigen::VectorXcd ConstVec;
    

    void Initialize(const Eigen::VectorXcd& inOrigVc,
                    const Eigen::MatrixXi& inD,
                    const Eigen::MatrixXi& inF,
                    const Eigen::MatrixXi& inE2V,
                    const Eigen::VectorXi& inConstIndices,
                    bool inisExactMC,
                    bool inisExactIAP,
                    double inRigidRatio){
        
        using namespace Eigen;
        using namespace std;
        
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
    
    void UpdateEnergy(const Eigen::VectorXd& CurrSolutionReal){
        
        using namespace Eigen;
        using namespace std;
        
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
    
    
    void UpdateGradient(const Eigen::VectorXd& CurrSolutionReal){
        
        using namespace Eigen;
        using namespace std;
        
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
    void Reformulate(int NumIteration, int MaxIteration, const Eigen::VectorXd& CurrSolutionReal, double PrevError)
    {
        Eigen::VectorXcd CurrSolution(CurrSolutionReal.size()/2);
        CurrSolution.array().real()=CurrSolutionReal.head(CurrSolutionReal.size()/2);
        CurrSolution.array().imag()=CurrSolutionReal.tail(CurrSolutionReal.size()/2);
        InitSolution=CurrSolution;

        double rate=ConstVec.lpNorm<Eigen::Infinity>()/PrevError;
        double ReduceRate=std::min(rate/2.0,1.0);
        
        SmoothFactor*=0.9-0.7*(1.0-ReduceRate);

    }
    
    void UpdateConstraints(const Eigen::VectorXd& CurrSolutionReal)
    {
        Eigen::VectorXcd CurrSolution(CurrSolutionReal.size()/2);
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
        return (ConstVec.lpNorm<Eigen::Infinity>()<ConstTolerance);
    }

};



#endif
