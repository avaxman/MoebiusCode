//#include <igl/readOFF.h>
#include <igl/viewer/Viewer.h>
#include <igl/unproject.h>
#include <igl/project.h>
#include <igl/jet.h>
#include "FileDialog.h"
#include "Deform3D.h"
#include "Deform2D.h"
#include "math.h"
#include "QuaternionOps.h"
#include "igl/serialize.h"
#include <igl/sortrows.h>
#include <igl/readOBJ.h>
#include <igl/writeOBJ.h>
#include <hedra/polygonal_read_OFF.h>
#include <hedra/triangulate_mesh.h>
#include <hedra/polygonal_edge_topology.h>
#include <igl/per_vertex_normals.h>
#include <algorithm>

using namespace Eigen;
using namespace std;

typedef enum {ORIGINAL, DEFORMATION, INTERPOLATION} EditingModes;
EditingModes EditingMode=DEFORMATION;

typedef enum {STANDARD, MC_ERROR, IA_ERROR, FACE_CONCYCLITY, FACE_MC, QC_ERROR, AREA_ERROR, SMALLK_ERROR} ViewingModes;
ViewingModes ViewingMode=STANDARD;

#define NUM_INTERPOLATION_STEPS 20
#define NUM_RIGID_STEPS 10

bool isExactIAP=false, isExactMC=false;
bool RecomputeInterpolation=true;  //after an external load or a deformation

bool Loadedtexture=false;
Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic> R, G, B;


igl::viewer::Viewer MobiusViewer;

// Mesh
Eigen::MatrixXd V;
Eigen::MatrixXi F;
Eigen::VectorXi D;
Eigen::MatrixXi tF;

double CurrWinZ;

//colormap for scalar value values
MatrixXd ColorMap;

// Per-face colors
MatrixXd FaceColors;

// pre-face scalar color [0,1]
VectorXd FaceValues;

MatrixXi FTC;  //for obj reading
MatrixXi FN;
MatrixXd TC;


bool Deforming=false;

std::vector<int> ConstIndices;
std::vector<Vector3d> ConstPoses;
int CurrActiveHandle;


Vector3d Mins, Maxs;
MatrixXd UV;
double uFreq, vFreq;

bool IS_COMPLEX;  //if the mesh if complex or quaternionic

double RigidRatio=0.5;       ///inversion control ratio
double MCSensitivity=0.01;  //ratio of abs(cr)
double IAPSensitivity=1.0;  //degrees difference
double QCSensitivity=0.5;   //abs(max_sing/min_sing)
double CircSensitivity=5.0;  //circle IA difference
double AreaSensitivity=0.1;
double SmallkSensitivity=0.25;  //Smallk error
double InterpScalar=0.0;


MoebiusDeformation2D md2;
MoebiusDeformation3D md3;


MatrixXd Scalar2RGB(const VectorXd ScalarValues)
{
    MatrixXd Result(ScalarValues.size(),3);
    for (int i=0;i<ScalarValues.size();i++){
        double CurrValue=ScalarValues(i);
        if (CurrValue<0.0) CurrValue=0.0;
        if (CurrValue>1.0) CurrValue=1.0;
        
        int Entry=floor(CurrValue*(ColorMap.rows()-1));
        if (Entry==ColorMap.rows()-1)
            Result.row(i)<<ColorMap.row(ColorMap.rows()-1);
        else{
            double Residue=CurrValue*(double)(ColorMap.rows()-1)-(double)(Entry);
            Result.row(i)=ColorMap.row(Entry)*(1-Residue)+ColorMap.row(Entry+1)*Residue;
            
        }
        
    }
    return Result;
}


void UpdateCurrentView()
{
    
    MobiusViewer.data.clear();
    MatrixXd LocalV, LocalEdgeV;
    switch(EditingMode){
        case ORIGINAL:{
            LocalV=(IS_COMPLEX ? md2.OrigV : md3.OrigV);
            LocalEdgeV=(IS_COMPLEX ? md2.EdgeOrigV : md3.EdgeOrigV);
            break;
        }
            
        case DEFORMATION:{
            LocalV=(IS_COMPLEX ? md2.DeformV : md3.DeformV);
            LocalEdgeV=(IS_COMPLEX ? md2.EdgeDeformV : md3.EdgeDeformV);
            break;
        }
        case INTERPOLATION:{
            LocalV=(IS_COMPLEX ? md2.InterpV : md3.InterpV);
            LocalEdgeV=(IS_COMPLEX ? md2.EdgeInterpV : md3.EdgeInterpV);
            break;
            
        }
    }
    
    cout<<"Viewing Mode in UpdateCurrentView(): <<"<<ViewingMode<<endl;
    
    if (ViewingMode==STANDARD){
        V=LocalV; tF=(IS_COMPLEX ? md2.tF : md3.tF);
        MobiusViewer.data.set_mesh(V,tF);
        if (Loadedtexture)
            MobiusViewer.data.set_texture(R,G,B);
        if (UV.rows()!=0){
            cout<<"setting UV"<<endl;
             MobiusViewer.data.set_uv(UV);
        }
        
    } else if ((ViewingMode==MC_ERROR)||(ViewingMode==IA_ERROR)){
        
        tF=(IS_COMPLEX ? md2.EdgeF : md3.EdgeF);
        V=LocalEdgeV;
        
        VectorXd AbsSecond,AbsFirst;
        VectorXd RealSecond, RealFirst;
        if (IS_COMPLEX){
            VectorXcd SecondCR;
            switch(EditingMode){
                case ORIGINAL: SecondCR=md2.OrigECR; break;
                case DEFORMATION:SecondCR=md2.DeformECR; break;
                case INTERPOLATION:SecondCR=md2.InterpECR; break;
                    
            }
            VectorXcd FirstCR=md2.OrigECR;
            
            AbsSecond=abs(SecondCR.array());
            AbsFirst=abs(FirstCR.array());
            
            RealSecond=SecondCR.real();
            RealFirst=FirstCR.real();
            
        } else {
            MatrixXd SecondCR;
            switch(EditingMode){
                case ORIGINAL: SecondCR=md3.OrigECR; break;
                case DEFORMATION:SecondCR=md3.DeformECR; break;
                case INTERPOLATION:SecondCR=md3.InterpECR; break;
                    
            }
            
            MatrixXd FirstCR=md3.OrigECR;
            
            AbsSecond=SecondCR.rowwise().norm();
            AbsFirst=FirstCR.rowwise().norm();
            
            RealSecond=SecondCR.col(0);
            RealFirst=FirstCR.col(0);
        }
        
        VectorXd AbsQuotient=(AbsSecond.cwiseQuotient(AbsFirst));
        VectorXd AngleFirst=acos(-RealFirst.cwiseQuotient(AbsFirst).array());
        VectorXd AngleSecond=acos(-RealSecond.cwiseQuotient(AbsSecond).array());
        
        //setting normal place for |cr2|/|cr1| from [1-MCSensitivity, 1+MCSensitivity]  -> [0, 1]
        AbsQuotient=((AbsQuotient.array()-1.0).cwiseAbs()/MCSensitivity).matrix();
        
        //setting normal place for anglediff from [-pi/180/IAPSensitivity, pi/180/IAPSensitivity] -> [0,1]
        VectorXd AngleDiff=(((AngleSecond-AngleFirst).cwiseAbs().array())/(M_PI/(180.0/IAPSensitivity))).matrix();
        
        
        FaceValues.resize(AbsQuotient.rows()*2);
        if (ViewingMode==MC_ERROR)
            FaceValues<<AbsQuotient,AbsQuotient;
        if (ViewingMode==IA_ERROR)
            FaceValues<<AngleDiff,AngleDiff;
        
        FaceColors=Scalar2RGB(FaceValues);
        
        MobiusViewer.data.set_mesh(V.cast<double>(),tF);
        MobiusViewer.data.set_colors(FaceColors.cast<double>());
        
    } else {  //face-based errors
        tF=(IS_COMPLEX ? md2.tF : md3.tF);
        V=LocalV;
        
        //Polygonal face-based concyclity and MC errors
        if ((ViewingMode==FACE_CONCYCLITY)||(ViewingMode==FACE_MC)){
            F=(IS_COMPLEX ? md2.F : md3.F);
            D=(IS_COMPLEX ? md2.D : md3.D);
            VectorXd AbsCR(F.rows()), RealCR(F.rows());
            AbsCR.setOnes(); RealCR.setConstant(-1);  //to get zero angle error for triangles
            if (IS_COMPLEX){
                VectorXcd FCR;
                switch(EditingMode){
                    case ORIGINAL: FCR=md2.OrigFCR; break;
                    case DEFORMATION:FCR=md2.DeformFCR; break;
                    case INTERPOLATION:FCR=md2.InterpFCR; break;
                        
                }
                
                int CurrFCR=0;
                for (int i=0;i<F.rows();i++){
                    if (D(i)==3)
                        continue;
                    RealCR(i)=0;
                    AbsCR(i)=0;
                    for (int j=0;j<D(i)-3;j++){
                        AbsCR(i)+=abs(FCR(CurrFCR))/(double)(D(i)-3);
                        RealCR(i)+=FCR(CurrFCR++).real()/(double)(D(i)-3);
                    }
                }
                
            } else {
                MatrixXd FCR;
                switch(EditingMode){
                    case ORIGINAL: FCR=md3.OrigFCR; break;
                    case DEFORMATION:FCR=md3.DeformFCR; break;
                    case INTERPOLATION:FCR=md3.InterpFCR; break;
                        
                }
                int CurrFCR=0;
                for (int i=0;i<F.rows();i++){
                    if (D(i)==3)
                        continue;
                    RealCR(i)=0;
                    AbsCR(i)=0;
                    for (int j=0;j<D(i)-3;j++){
                        AbsCR(i)+=(FCR.row(CurrFCR).norm())/(double)(D(i)-3);
                        RealCR(i)+=FCR(CurrFCR++,0)/(double)(D(i)-3);
                    }
                }
                
                
            }
            
            VectorXd Angle=acos(-RealCR.cwiseQuotient(AbsCR).array());
            FaceValues.resize(tF.rows());
            FaceValues.setZero();
            FaceColors.resize(FaceValues.rows(),3);
            if (ViewingMode==FACE_CONCYCLITY)
                for (int i=0;i<(IS_COMPLEX ? md2.tF.rows() : md3.tF.rows());i++)
                    FaceValues(i)=abs(Angle((IS_COMPLEX ? md2.FromFace(i) : md3.FromFace(i))))/(M_PI/(180.0/CircSensitivity));
            else  //FACE_MC
                for (int i=0;i<(IS_COMPLEX ? md2.tF.rows() : md3.tF.rows());i++)
                    FaceValues(i)=AbsCR((IS_COMPLEX ? md2.FromFace(i) : md3.FromFace(i)));
            
            FaceColors=Scalar2RGB(FaceValues);
            
        }
        
        if ((ViewingMode==QC_ERROR)||(ViewingMode==SMALLK_ERROR)){
            Vector2d TrueSings;
            FaceValues.resize(tF.rows());
            FaceValues.setZero();
            FaceColors.resize(FaceValues.rows(),3);
            if (!IS_COMPLEX){  //computing normals and affine maps
                for (int i=0;i<tF.rows();i++){
                    Vector3d p01=md3.OrigV.row(tF(i,1))-md3.OrigV.row(tF(i,0));
                    Vector3d p12=md3.OrigV.row(tF(i,2))-md3.OrigV.row(tF(i,1));
                    Vector3d p20=md3.OrigV.row(tF(i,0))-md3.OrigV.row(tF(i,2));
                    
                    Vector3d q01=V.row(tF(i,1))-V.row(tF(i,0));
                    Vector3d q12=V.row(tF(i,2))-V.row(tF(i,1));
                    Vector3d q20=V.row(tF(i,0))-V.row(tF(i,2));
                    
                    Vector3d np=(p12.cross(p01)).normalized();
                    Vector3d nq=(q12.cross(q01)).normalized();
                    
                    Matrix3d pMat; pMat<<p01, p12, np;
                    Matrix3d qMat; qMat<<q01, q12, nq;
                    
                    Matrix3d A=pMat*qMat.inverse();
                    
                    JacobiSVD<Matrix3d> svd(A);
                    Vector3d Sings=(svd.singularValues());
                    
                    TrueSings.setOnes();
                    int Count=0;
                    for (int i=0;i<3;i++)
                        if (abs(Sings(i)-1.0)>10e-5)
                            TrueSings(Count++)=abs(Sings(i));
                    
                    //cout<<"Sings: "<<Sings<<endl;
                    //cout<<"True Sings: "<<TrueSings<<endl;
                    double smax=max(TrueSings(0),TrueSings(1));
                    double smin=min(TrueSings(0), TrueSings(1));
                    double Bigk=smax/smin;
                    if (ViewingMode==QC_ERROR)
                        FaceValues(i)=abs(Bigk-1.0)/QCSensitivity;
                    else
                        FaceValues(i)=((Bigk-1.0)/(Bigk+1.0))/SmallkSensitivity;
                    
                }
            } else {
                for (int i=0;i<tF.rows();i++){
                    Vector2d p01=md2.OrigV.row(tF(i,1)).head(2)-md2.OrigV.row(tF(i,0)).head(2);
                    Vector2d p12=md2.OrigV.row(tF(i,2)).head(2)-md2.OrigV.row(tF(i,0)).head(2);
                    
                    Vector2d q01=V.row(tF(i,1)).head(2)-V.row(tF(i,0)).head(2);
                    Vector2d q12=V.row(tF(i,2)).head(2)-V.row(tF(i,0)).head(2);
                    
                    
                    Matrix2d pMat; pMat<<p01, p12;
                    Matrix2d qMat; qMat<<q01, q12;
                    
                    Matrix2d A=pMat*qMat.inverse();
                    
                    JacobiSVD<Matrix2d> svd(A);
                    Vector2d TrueSings=(svd.singularValues());
                    
                    if ((TrueSings(0)<=0.0)||(TrueSings(1)<=0.0))
                        cout<<"Non local Injectivity with values: "<<TrueSings<<endl;
                    
                    double smax=max(abs(TrueSings(0)),abs(TrueSings(1)));
                    double smin=min(abs(TrueSings(0)), abs(TrueSings(1)));
                    double Bigk=smax/smin;
                    if (ViewingMode==QC_ERROR)
                        FaceValues(i)=abs(Bigk-1.0)/QCSensitivity;
                    else
                        FaceValues(i)=((Bigk-1.0)/(Bigk+1.0))/SmallkSensitivity;
                    
                }
                
                
            }
        }
        
        if (ViewingMode==AREA_ERROR){
            FaceValues.resize(tF.rows());
            FaceValues.setZero();
            FaceColors.resize(FaceValues.rows(),3);
            for (int i=0;i<tF.rows();i++){
                if (!IS_COMPLEX){
                    Vector3d p01=md3.OrigV.row(tF(i,1))-md3.OrigV.row(tF(i,0));
                    Vector3d p12=md3.OrigV.row(tF(i,2))-md3.OrigV.row(tF(i,1));
                    
                    Vector3d q01=V.row(tF(i,1))-V.row(tF(i,0));
                    Vector3d q12=V.row(tF(i,2))-V.row(tF(i,1));
                    
                    double Areap=(p01.cross(p12)).norm()/2.0;
                    double Areaq=(q01.cross(q12)).norm()/2.0;
                    
                    FaceValues(i)=abs(Areaq/Areap-1.0)/AreaSensitivity;
                } else {
                    Vector2d p01=md2.OrigV.row(tF(i,1)).head(2)-md2.OrigV.row(tF(i,0)).head(2);
                    Vector2d p12=md2.OrigV.row(tF(i,2)).head(2)-md2.OrigV.row(tF(i,1)).head(2);
                    
                    Vector2d q01=V.row(tF(i,1)).head(2)-V.row(tF(i,0)).head(2);
                    Vector2d q12=V.row(tF(i,2)).head(2)-V.row(tF(i,1)).head(2);
                    
                    double Areap=abs(p01(0)*p12(1)-p01(1)*p12(0))/2.0;
                    double Areaq=abs(q01(0)*q12(1)-q01(1)*q12(0))/2.0;
                    
                    FaceValues(i)=abs(Areaq/Areap-1.0)/AreaSensitivity;
                }
                
            }
        }
        
        
        if (ViewingMode==SMALLK_ERROR)
            cout<<"Small k errors between ["<<FaceValues.minCoeff()*SmallkSensitivity<<","<<FaceValues.maxCoeff()*SmallkSensitivity<<"]"<<endl;
        
        FaceColors=Scalar2RGB(FaceValues);
        
        MobiusViewer.data.set_mesh(V.cast<double>(),tF);
        MobiusViewer.data.set_colors(FaceColors.cast<double>());
    }
    
    
    
    if ((EditingMode==DEFORMATION)&&(ConstIndices.size()!=0)){
        VectorXi EConstIndices(ConstIndices.size());
        for (int i=0;i<ConstIndices.size();i++)
            EConstIndices(i)=ConstIndices[i];
        
        MatrixXd EConstPoses(ConstIndices.size(),3);
        for (int i=0;i<ConstIndices.size();i++)
            EConstPoses.row(i)=ConstPoses[i];
        
        
        MatrixXd HandleColors(ConstIndices.size(),3);
        HandleColors.col(0)=VectorXd::Constant(ConstIndices.size(),1.0);
        HandleColors.col(1)=VectorXd::Constant(ConstIndices.size(),0.5);
        HandleColors.col(2)=VectorXd::Constant(ConstIndices.size(),0.0);
        HandleColors.row(CurrActiveHandle)<<0.0,1.0,0.5;
        
        MobiusViewer.data.set_points(EConstPoses.cast<double>(), HandleColors.cast<double>());
    }
    
    
    //drawing only original edges
    vector<pair<int,int> > OrigEdges;
    for (int i=0;i<(IS_COMPLEX ? md2.F.rows() : md3.F.rows());i++){
        int PolygonSize=(IS_COMPLEX ? md2.D(i) : md3.D(i));
        for (int j=0;j<(IS_COMPLEX ? md2.D(i) : md3.D(i));j++){
            int CurrVertex=(IS_COMPLEX ? md2.F(i,j) : md3.F(i,j));
            int NextVertex=(IS_COMPLEX ? md2.F(i,(j+1)%PolygonSize) : md3.F(i,(j+1)%PolygonSize));
            OrigEdges.push_back(pair<int,int> (CurrVertex, NextVertex));
        }
    }
    
    MatrixXd OrigEdges1(OrigEdges.size(),3);
    MatrixXd OrigEdges2(OrigEdges.size(),3);
    for (int i=0;i<OrigEdges.size();i++){
        OrigEdges1.row(i)=V.row(OrigEdges[i].first);
        OrigEdges2.row(i)=V.row(OrigEdges[i].second);
    }
    
    MatrixXd OrigEdgeColors(OrigEdges1.rows(),3);
    OrigEdgeColors.col(0)=VectorXd::Constant(OrigEdges1.rows(),0.25);
    OrigEdgeColors.col(1)=VectorXd::Constant(OrigEdges1.rows(),0.25);
    OrigEdgeColors.col(2)=VectorXd::Constant(OrigEdges1.rows(),0.0);
    
    //MobiusViewer.data.add_edges(OrigEdges1.cast<double>(),OrigEdges2.cast<double>(),OrigEdgeColors.cast<double>());
}

void SaveGeneralOff(const std::string str, const MatrixXd& V, const MatrixXi& D, const MatrixXi& F)
{
    /*ofstream FileName;
    FileName.open(str);
    
    FileName<<"OFF"<<endl<<V.rows()<<" "<<F.rows()<<" 0"<<endl;
    FileName<<V<<endl;
    MatrixXi FD(D.rows(), D.cols()+F.cols());
    FD<<D, F;
    FileName<<FD<<endl;
    FileName.close();*/
    
    //igl::writePolygonalOFF(str, V, D, F);
    
}


//Currently reading only triangle meshes!!!
void ReadGeneralObj(const std::string str, MatrixXd& V, MatrixXi& F)
{
    vector<vector<double> > N, vV, vUV, vTC;
    vector<vector<int> > vFTC, vFN, vF;
    MatrixXi TempF;
    
    MatrixXd RawUV;
    
    /*const std::string obj_file_name,
    std::vector<std::vector<Scalar > > & V,
    std::vector<std::vector<Scalar > > & TC,
    std::vector<std::vector<Scalar > > & N,
    std::vector<std::vector<Index > > & F,
    std::vector<std::vector<Index > > & FTC,
    std::vector<std::vector<Index > > & FN);*/
    
    igl::readOBJ(str, vV, vTC, N, vF, vFTC, vFN);
    
    igl::list_to_matrix(vV,V);
    igl::list_to_matrix(vF,TempF);
    igl::list_to_matrix(vTC, TC);
    igl::list_to_matrix(vFTC, FTC);
    igl::list_to_matrix(vFTC, FN);
    
    //reallocating UV to normal form
    UV.resize(3*TempF.rows(),2);
    for (int i=0;i<FTC.rows();i++)
        for (int j=0;j<FTC.cols();j++)
            UV.row(3*i+j)=TC.row(FTC(i,j));

    //cout<<"Current UV:"<<UV.row(0)<<endl;
    
    IS_COMPLEX=(((V.col(2).cwiseAbs()).array()).minCoeff()<10e-6);
    V.rowwise()-=V.colwise().mean();
    
    F.resize(TempF.rows(),4);
    F.block(0,1,TempF.rows(), 3)=TempF;
    F.col(0).setConstant(3);
    
    int MinIndex=F.minCoeff();
    F-=MatrixXi::Constant(F.rows(), F.cols(), MinIndex);
    F.col(0).array()+=MinIndex;
}

//Currently reading only triangle meshes!!!
void SaveOBJ(const std::string str, MatrixXd& V, MatrixXi& F)
{
    MatrixXd CN;
    
    igl::writeOBJ(str, V, F, CN, FN, UV, FTC);
}



void ReadGeneralOff(const std::string str, MatrixXd& V, VectorXi& D, MatrixXi& F)
{
    
    hedra::polygonal_read_OFF(str, V, D, F);

    //determining if complex or not
    IS_COMPLEX=(V.lpNorm<Infinity>()<10e-5);
    
    
    UV=V.block(0,0,V.rows(),2);
    
    //centering object
    V.rowwise()-=V.colwise().mean();

}



bool key_down(igl::viewer::Viewer& viewer, unsigned char key, int modifiers)
{
    switch(key)
    {
        case 'X':{
            char buffer[1024];
            get_save_file_path(buffer);
            if (EditingMode==INTERPOLATION)
                SaveGeneralOff(buffer, (IS_COMPLEX ? md2.InterpV : md3.InterpV), (IS_COMPLEX ? md2.D : md3.D), (IS_COMPLEX ? md2.F : md3.F));
            else
                SaveGeneralOff(buffer, (IS_COMPLEX ? md2.DeformV : md3.DeformV), (IS_COMPLEX ? md2.D : md3.D), (IS_COMPLEX ? md2.F : md3.F));
            
            break;
        }
            
        case 'P':{
            MatrixXd ConstPosesMat(ConstPoses.size(),3);
            for (int i=0;i<ConstPoses.size();i++)
                ConstPosesMat.row(i)=ConstPoses[i];
            
            VectorXi EConstIndices(ConstIndices.size());
            for (int i=0;i<ConstIndices.size();i++)
                EConstIndices(i)=ConstIndices[i];
            
            
            (IS_COMPLEX ? md2.InitDeformation(EConstIndices, isExactMC, isExactIAP, RigidRatio): md3.InitDeformation(EConstIndices, isExactMC, isExactIAP, RigidRatio));
            
            (IS_COMPLEX ? md2.UpdateDeformation(ConstPosesMat,1000, isExactMC, isExactIAP) : md3.UpdateDeformation(ConstPosesMat,1000));
            CurrActiveHandle=(int)(ConstIndices.size()-1);
            UpdateCurrentView();
            break;
        }

        case 'V':{
            UV.array()*=1.05;
            cout<<"Pressed V:"<<endl;
            UpdateCurrentView();
            break;
        }
            
        case 'B':{
            UV.array()/=1.05;
            UpdateCurrentView();
            break;
        }
            
        case 'N':
            CurrActiveHandle=(CurrActiveHandle+1)%ConstIndices.size();
            UpdateCurrentView();
            break;
            
        case 'M':
            CurrActiveHandle=(CurrActiveHandle+ConstIndices.size()-1)%ConstIndices.size();
            UpdateCurrentView();
            break;
            
        case 'D':  //delete current handle
            std::vector<int> TempConstIndices;
            std::vector<Vector3d> TempConstPoses;
            for (int i=0;i<CurrActiveHandle;i++){
                TempConstIndices.push_back(ConstIndices[i]);
                TempConstPoses.push_back(ConstPoses[i]);
            }
            for (int i=CurrActiveHandle+1;i<ConstIndices.size();i++){
                TempConstIndices.push_back(ConstIndices[i]);
                TempConstPoses.push_back(ConstPoses[i]);
            }
            ConstIndices=TempConstIndices;
            ConstPoses=TempConstPoses;
            if (ConstIndices.empty())
                CurrActiveHandle=-1;
            else
                CurrActiveHandle=0;
            UpdateCurrentView();
            break;
    }
    return true;
}

bool mouse_up(igl::viewer::Viewer& viewer, int button, int modifier)
{
    if (((igl::viewer::Viewer::MouseButton)button==igl::viewer::Viewer::MouseButton::Left))
        return false;
    
    Deforming=false;
    return true;
}

bool mouse_move(igl::viewer::Viewer& viewer, int mouse_x, int mouse_y)
{
    
    double x = viewer.current_mouse_x;
    double y = viewer.core.viewport(3) - viewer.current_mouse_y;
    Vector3f NewPos=igl::unproject<float>(Vector3f(x,y,CurrWinZ),
                                   viewer.core.view * viewer.core.model,
                                   viewer.core.proj,
                                   viewer.core.viewport);
    
    /*IGL_INLINE Eigen::Matrix<Scalar,3,1> igl::unproject(
                                                        const    Eigen::Matrix<Scalar,3,1>&  win,
                                                        const    Eigen::Matrix<Scalar,4,4>& model,
                                                        const    Eigen::Matrix<Scalar,4,4>& proj,
                                                        const    Eigen::Matrix<Scalar,4,1>&  viewport)*/
    
    if (!Deforming)
        return false;
    
    
    //just updating the deformation
    if (EditingMode==DEFORMATION){
        if (IS_COMPLEX)
            ConstPoses[CurrActiveHandle]<<(double)NewPos(0),(double)NewPos(1),0.0;
        else
            ConstPoses[CurrActiveHandle]=NewPos.cast<double>();
        
        MatrixXd ConstPosesMat(ConstPoses.size(),3);
        for (int i=0;i<ConstPoses.size();i++)
            ConstPosesMat.row(i)=ConstPoses[i];
        
        (IS_COMPLEX ? md2.UpdateDeformation(ConstPosesMat,200, isExactMC, isExactIAP) : md3.UpdateDeformation(ConstPosesMat,200));
        RecomputeInterpolation=true;
    }
    
    UpdateCurrentView();
    
    
    return true;
}


bool mouse_down(igl::viewer::Viewer& viewer, int button, int modifier)
{
    if ((igl::viewer::Viewer::MouseButton)button==igl::viewer::Viewer::MouseButton::Left)
        return false;
    
    //finding out which vertex
    int ClosestVertex=0;
    
    // Cast a ray in the view direction starting from the mouse position
    double x = viewer.current_mouse_x;
    double y = (double)viewer.core.viewport(3) - (double)viewer.current_mouse_y;
    double MinDistance=3276700.0;
    V=(IS_COMPLEX ? md2.DeformV : md3.DeformV);  //not deforming interpolations!
    for (int i=0;i<V.rows();i++){
        Vector3f WinCoords=igl::project<float>(V.row(i).cast<float>(),
                                        viewer.core.view * viewer.core.model,
                                        viewer.core.proj,
                                        viewer.core.viewport);
        double Distance=(WinCoords(0)-x)*(WinCoords(0)-x)+(WinCoords(1)-y)*(WinCoords(1)-y)+WinCoords(2)*WinCoords(2);
        
        if (Distance<MinDistance){
            ClosestVertex=i;
            MinDistance=Distance;
            CurrWinZ=WinCoords(2);
        }
    }
    
    if (EditingMode==DEFORMATION){
        
        //Creating a new handle, or bringing the focus to an old one.
        if ((igl::viewer::Viewer::MouseButton)button==igl::viewer::Viewer::MouseButton::Middle){
            
            //checking if handle already exists
            bool Found=false;
            for (int i=0;i<ConstIndices.size();i++){
                if (ConstIndices[i]==ClosestVertex){
                    Found=true;
                    CurrActiveHandle=i;
                    break;
                }
            }
            
            if (!Found){
                ConstIndices.push_back(ClosestVertex);
                ConstPoses.push_back(V.row(ClosestVertex));
                CurrActiveHandle=(int)(ConstIndices.size()-1);
            }
            
            cout << "Picked (vertex): "<<CurrActiveHandle<< ")" << endl;
        }
        
        //Starting to deform a given vertex
        if ((igl::viewer::Viewer::MouseButton)button==igl::viewer::Viewer::MouseButton::Right)
            Deforming=true;
    }
    
    
    if (EditingMode==DEFORMATION){
        VectorXi EConstIndices(ConstIndices.size());
        for (int i=0;i<ConstIndices.size();i++)
            EConstIndices(i)=ConstIndices[i];
        (IS_COMPLEX ? md2.InitDeformation(EConstIndices, isExactMC, isExactIAP, RigidRatio) : md3.InitDeformation(EConstIndices, isExactMC, isExactIAP, RigidRatio));
    }
    
    
    
    UpdateCurrentView();
    return true;
}

void SetEditingModeCallback(EditingModes value)
{
    EditingMode = value;
    cout<<"Editing Mode: "<<EditingMode<<endl;
    UpdateCurrentView();
    
}

EditingModes GetEditingModeCallback()
{
    return EditingMode;
    
}

void SetViewingModeCallback(ViewingModes value)
{
    ViewingMode = value;
    cout<<"Viewing Mode: "<<ViewingMode<<endl;
    UpdateCurrentView();
    
}

ViewingModes GetViewingModeCallback()
{
    return ViewingMode;
    
}

/*void nothing(){
 
 }*/

void SetInterpolationFactorCallback(double value)
{
    InterpScalar=value;
    cout<<"Current Interpolation Factor is "<<InterpScalar<<endl;
    if (RecomputeInterpolation){
        (IS_COMPLEX ? md2.SetupInterpolation(isExactMC, isExactIAP) : md3.SetupInterpolation(isExactMC, isExactIAP, true));
        RecomputeInterpolation=false;
    }
    (IS_COMPLEX ? md2.Interpolate(InterpScalar, 500) : md3.Interpolate(InterpScalar, 500));
    UpdateCurrentView();
    
}

double GetInterpolationFactorCallback()
{
    return(InterpScalar);
    
}


void SetRigidRatioCallback(double value)
{
    RigidRatio=value;
    UpdateCurrentView();
    
}



double GetRigidRatioCallback()
{
    return RigidRatio;
    
}


void SetMCSenseCallback(double value)
{
    MCSensitivity=value;
    UpdateCurrentView();
}



double GetMCSenseCallback()
{
    return MCSensitivity;
    
}


void SetIAPSenseCallback(double value)
{
    IAPSensitivity=value;
    UpdateCurrentView();
}

double GetIAPSenseCallback()
{
    return IAPSensitivity;
    
}



void SetQCSenseCallback(double value)
{
    QCSensitivity=value;
    UpdateCurrentView();
}



double GetQCSenseCallback()
{
    return QCSensitivity;
}


void SetCircSenseCallback(double value)
{
    CircSensitivity=value;
    UpdateCurrentView();
}


double GetCircSenseCallback()
{
    return CircSensitivity;
    
}



double GetAreaSenseCallback()
{
    return AreaSensitivity;
    
}

void SetAreaSenseCallback(double value)
{
    AreaSensitivity=value;
    UpdateCurrentView();
}

double GetSmallkCallback()
{
    return SmallkSensitivity;
    
}

void SetSmallkCallback(double value)
{
    SmallkSensitivity=value;
    UpdateCurrentView();
}




}

bool init(igl::viewer::Viewer& viewer)
{
    using namespace nanogui;
    
    viewer.ngui->addGroup("Algorithm Options");
    
    viewer.ngui->addVariable<EditingModes>("Editing Mode", SetEditingModeCallback, GetEditingModeCallback)->setItems({"Original","Deformation","Interpolation"});
    
    viewer.ngui->addVariable<ViewingModes>("Viewing Mode", SetViewingModeCallback, GetViewingModeCallback)->setItems({"Mesh View","Edge MC Error","Edge IA Error", "Concyclity", "Face MC Error", "Face QC Error", "Face Area Error", "Small k Error"});
    
    viewer.ngui->addVariable<double>("Interpolation Step",SetInterpolationFactorCallback, GetInterpolationFactorCallback);
    viewer.ngui->addVariable<double>("Rigidity Ratio", SetRigidRatioCallback, GetRigidRatioCallback);
    
    viewer.ngui->addVariable<double>("MC Sensitivity",SetMCSenseCallback, GetMCSenseCallback);
    viewer.ngui->addVariable<double>("IAP Sensitivity", SetIAPSenseCallback, GetIAPSenseCallback);
    viewer.ngui->addVariable<double>("QC Sensitivity", SetQCSenseCallback, GetQCSenseCallback);
    viewer.ngui->addVariable<double>("Circularity Sensitivity", SetCircSenseCallback, GetCircSenseCallback);
    viewer.ngui->addVariable<double>("Area Sensitivity", SetAreaSenseCallback, GetAreaSenseCallback);
    viewer.ngui->addVariable<double>("Small k Sensitivity", SetSmallkCallback, GetSmallkCallback);
    
    viewer.ngui->addVariable<bool>("Exact MC", isExactMC);
    viewer.ngui->addVariable<bool>("Exact IAP", isExactIAP);
    
    viewer.screen->performLayout();

    return false;
}

void SetupColormap(){
    ColorMap.resize(33, 3);
    ColorMap<<59, 76, 192,
    68, 90, 204,
    77, 104, 215,
    87 ,117 ,225,
    98 ,130 ,234,
    108 ,142 ,241,
    119 ,154 ,247,
    130 ,165 ,251,
    141 ,176 ,254,
    152 ,185 ,255,
    163 ,194 ,255,
    174 ,201 ,253,
    184 ,208 ,249,
    194 ,213 ,244,
    204 ,217 ,238,
    213 ,219 ,230,
    221, 221 ,221,
    229 ,216 ,209,
    236 ,211 ,197,
    241 ,204 ,185,
    245 ,196 ,173,
    247 ,187 ,160,
    247 ,177 ,148,
    247 ,166 ,135,
    244 ,154 ,123,
    241 ,141 ,111,
    236 ,127 ,99,
    229 ,112 ,88,
    222 ,96 ,77,
    213 ,80 ,66,
    203 ,62 ,56,
    192 ,40 ,47,
    180 ,4, 38;
    
    ColorMap/=255.0;
}





int main(int argc, char *argv[])
{
    char buffer[1024];
    get_open_file_path(buffer);
    std::string buffstr(buffer);
    if(strstr(buffer, ".obj") == NULL) { //an off file
        ReadGeneralOff(buffer, V, D, F);
        //F.array()-=1;
        MobiusViewer.data.set_uv(UV);
    }else
        ReadGeneralObj(buffer, V, F);
    
    (IS_COMPLEX ? md2.SetupMesh(V, D, F): md3.SetupMesh(V, D, F));
    
    if (IS_COMPLEX)
        cout<<"***************Complex Mesh**************"<<endl;
    else
        cout<<"***************Quaternionic Mesh*********"<<endl;
    
    
    SetupColormap();
    
    tF=(IS_COMPLEX ? md2.tF : md3.tF);
    
    Mins=V.colwise().minCoeff();
    Maxs=V.colwise().maxCoeff();
    
    //glPointSize(2.0);
    
    Vector3d Spans=Mins-Maxs;
    MobiusViewer.core.background_color<<0.5,0.5,0.5,1.0;
    
    UpdateCurrentView();
    
    
    
    uFreq=vFreq=1.0;
    
    // Show mesh
    MobiusViewer.callback_mouse_down = &mouse_down;
    MobiusViewer.callback_mouse_up = &mouse_up;
    MobiusViewer.callback_mouse_move = &mouse_move;
    MobiusViewer.callback_key_down = &key_down;
    MobiusViewer.callback_init = &init;
    MobiusViewer.launch();
    
}
