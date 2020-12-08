#ifndef TTSP_3D_H
#define TTSP_3D_H

#include <vector>
#include <utility>
#include <string>
#include <filesystem>
//#include "T_mesh.h"
#include "BasicDataStructure.h"

using namespace std;

class HexQuality
{
public:
	HexQuality();


	void BuildInitialEdges();
	void InitialConnect();
	void InitializeMesh(string fn);
	void SetSharpFeature_1(double tol);
	void SetSharpFeature_Manual(string fn);

	void OutputCM(string fn);

	void run_MeshQualityImprove(int mode, int flag_sharp, double tol_sharp, int opt_par1, double opt_par2, string fn);

	//quality
	int IsPillowNeeded();
	void Pillow(int nlayer = 2);
	void Smoothing(int nStep, double stepSize);//move interior points towards mass center

	int SmoothingPoint(int pid, double stepSize);
	int LaplaceSmoothingPoint(int pid, double stepSize);
	int SmoothingPointBoundary(int pid, double stepSize);
	int SmoothingPointBoundarySharp(int pid, double stepSize);
	void GetHexVolAndCenter(int eid, double& vol, double center[3]);
	void GetQuadInfo(int fcid, int ploc, double& area, double center[3], double nm[3]);
	void GetHexMinJacob(int eid, double& minJacob, double GaussPos[3]);
	void GetGaussPoint(int ng, vector<double>& Gpt, vector<double>& wght);
	void JacobEval(int eid, double u, double v, double w, double& detJ);
	void JacobEval_Scale(int eid, double u, double v, double w, double& detJ);
	void GlobalMinJacob(double& minJacob_glb, int& min_pos, vector<int>& nBadEle);
	void Optimizing(int nStep, double stepSize);
	int OptimizingElement(int eid, double stepSize);
	int OptimizeElement(int eid, double stepSize);
	double AverageEdgeLength(int eid);
	void JacobEval_Grad(int eid, int ploc, double Jacob_Grad[3]);

	void Smoothing_Adapt(int nSize, double stepSize);
	int SmoothingPoint_Adapt(int lev, int pid, double stepSize);
	int SmoothingPointBoundary_Adapt(int lev, int pid, double stepSize);
	void GetHexMinJacob_Adapt(int lev, int eid, double& minJacob_ele, double GaussPos[3]);
	void GetHexVolAndCenter_Adapt(int lev, int eid, double& vol, double center[3]);
	void GetQuadInfo_Adapt(int lev, int fcid, int ploc, double& area, double center[3], double nm[3]);
	//void JacobEval_Adapt(int eid, double u, double v, double w, double& detJ);
	void JacobEval_Scale_Adapt(int lev, int eid, double u, double v, double w, double detJ[2]);
	void GlobalMinJacob_Adapt(double& minJacob_glb, array<int, 2>& min_pos, vector<array<int, 2>>& BadEle);
	void BadElementFlag_Adapt(const vector<array<int,2>>& rfid);


	void Optimizing_1(int nStep, int pnitrMax, double stepSize);
	double MaxEdgeLength();
	double OptimizePoint(int pid, int nitrMax, double stepSize);
	void TranslateScale(int pid, double trs[3], vector<double>& scl);
	void TranslateScale_Reverse(int pid, double trs[3], const vector<double>& scl);
	void GetAdvanceDirection(int pid, double grad[3], double dir[3], double& f);
	double GetStepSize(int pid, double grad[3], double dir[3], double keta);
	void objGrad(int eid, int iloc, double delta, double& obj, double grad[3], double grad2[3][3]);
	void GetJacobMat(int eid, int iloc, double Jmat[3][3], double DJmat[3][3][3], double& detJ, double grad[3]);

	void Optimizing_glb(int nStep, double stepSize);

	//quality Laplace
	void LaplaceSmoothing(int nstep = 100);
	void LaplaceSmooth_Interior(int pid);
	void LaplaceSmooth_Boundary_NonSharp(int pid);
	void LaplaceSmooth_Boundary_Sharp(int pid);

	bool ComputeElementNormal_Quad(int cnct[4], array<double, 3>& normal_quad);
	bool CrossProduct(double vector_1[3], double vector_2[3], array<double, 3>& result);
	bool DotProduct(vector<array<double, 3>> check_normal_surface_from_edge, double & angle);

	//vector<int> paid;//active control points
	vector<int> haid;//active hex elements
	vector<int> faid;//active faces
	vector<int> eaid;//active edges

private:
	vector<Vertex3D> cp;//control points
	vector<Vertex3D> cp_smooth; // for smoothing, store control points after smoothing
	vector<array<double, 3>> cpa;
	vector<Element3D> tmesh;//elements in T-mesh
	vector<Edge3D> tmedge;
	vector<Face3D> tmface;
	unsigned int npt_old;
	unsigned int nel_old;
	unsigned int nfc_old;
	unsigned int ned_old;
	vector<vector<Vertex3D>> hcp;
	vector<vector<Element3D>> hmesh;
	vector<vector<Edge3D>> hedge;
	vector<vector<Face3D>> hface;
	vector<vector<double>> kvec;

	vector<BPatch3D> bsp;

	vector<Element3D> tmtmp;//elements in T-mesh
	vector<Edge3D> tmedtmp;
	vector<Face3D> tmfctmp;

	vector<MFNode> mp;//meshfree nodes
	vector<array<double, 3>> bzcp;//Bezier control points
	//vector<array<double, 3>> bzcp_c01;//Bezier control points
	vector<int> bzcp_c01;

	vector<array<double, 3>> ptri;
	vector<array<int, 3>> etri;

	//used for solution
	double dmrg[3][2];
	double nmpl[3];
	double acoef;
	double dmlen[3];
};

#endif