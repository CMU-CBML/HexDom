#ifndef BASIC_DATA_STRUCTURE_H
#define BASIC_DATA_STRUCTURE_H

#include <vector>
#include <array>
#include <string>
#include <fstream>
#include <iostream>

#include <Eigen/Dense>

using namespace std;
using namespace Eigen;

//const int phys_dim(3);

/////////////////////////////////

class ELCS2D
{
public:
	int rot;//0, 1, 2, 3
	double u[2];
};

class ELCS3D
{
public:
	int rot[3];//1 for +u and -1 for -u, 2 for +v and -2 for -v, 3 for +w and -3 for -w
	double u[3];
};

class Edge
{
public:
	Edge();
	int act;
	int pt[2];//two ends
	int lev;
	vector<int> face;//faces sharing this edge
	double len;//parametric length
	int pn[2][2];//previous edge (or face) and next edge (or face), 1st flag (0 for edge, 1 for face, 2 for XP, 3 for end), 2nd ID
	int prt;
	int chd[2];
	int midpt;//middle point that splits the edge into two
	bool operator==(const Edge& ed);
};

class Edge3D
{
public:
	Edge3D();
	int act;
	int id_act;
	int pt[2];//two ends
	int lev;
	int type;//0 for regular, 1 for boundary, 2 for extraordinary
	int sharp;
	vector<int> face;//faces sharing this edge
	vector<int> hex;
	double len;//parametric length
	int pn[2][2];//previous edge (or face) and next edge (or face), 1st flag (0 for edge, 1 for face, 2 for hex, 3 for XP, 4 for end), 2nd ID
	int prt;
	int chd[2];
	int midpt;//middle point that splits the edge into two
	bool operator==(const Edge3D& ed);

	//used for coupling meshfree
	int mpflag;//0 not assigned, 1 assigned
	vector<int> mpid;//ID of meshfree nodes

	vector<int> bzid;
	int c0flag;
	int c0flag_b;

	int ped_id;
};

class Face3D
{
public:
	Face3D();
	int act;
	int id_act;
	int lev;
	int type;//0 for interior, 1 for boundary
	int ctpt;
	int cnct[4];//four corners, could be in reverse order
	int edge[4];
	int pn[4][2];//relation of 4 edges with neighbor, 1st flag (0 for edge, 1 for face, 2 for hex, 3 for XP, 4 for end), 2nd ID
	vector<int> Tedge;
	vector<int> hex;//solids sharing this face
	//int pn[2][2];//previous edge (or face) and next edge (or face), 1st flag (0 for edge, 1 for face, 2 for XP, 3 for end), 2nd ID
	int prt;
	vector<int> chd;
	bool operator==(const Face3D& fc);

	//used for coupling meshfree
	int mpflag;//0 not assigned, 1 assigned
	vector<int> mpid;//ID of meshfree nodes

	vector<int> bzid;
	int c0flag;

	int pfc_id;
};

class Vertex2D
{
public:
	double coor[3];
	double coortmp[3];
	int uved[2];//global u-v direction indicated by two edges
	//PointIPP ipp;
	int update;//0 not updated, 1 updated, 2 later coordinates, 3 later update tbf and tc
	double kitvU[4];//knot interval in local definition
	double kitvV[4];
	double kitvUtmp[4];
	double kitvVtmp[4];
	int type;//0 for regular (including corner and boundary), 1 for T-junction, 2 for extraordinary
	int trun;
	int rfc;//reference face
	int aff;
	vector<int> tbf;//truncated basis functions
	vector<double> tc;//truncation coefficients
	vector<int> edge;//edges that connect to this vertex
	vector<int> face;//faces that share this vertex

	Vertex2D();
	//bool operator==(const Vertex2D& v);
};

class Vertex3D
{
public:
	double coor[3];
	double coortmp[3];
	double wght;
	int uved[3];//global u-v-w direction indicated by three edges
	//PointIPP ipp;
	double kitvU[4];//knot interval in local definition
	double kitvV[4];
	double kitvW[4];
	int rhx;//reference hex
	Matrix3d rot_ref;//a rotation wrt the reference hex
	int type;//0 for regular, 1 for boundary, 2 for T-junction, 3 for extraordinary, 12 for boundary T-junc, 13 for boundary XP
	int trun;
	int aff;
	int update;
	int lev;
	int act;
	int sharp;//0 for non-sharp, 1 for sharp edge, 2 for sharp corner
	int bcxp;//boundary extraordinary point, 0 false, 1 for symmetric, 2 for unsymmetric
	vector<array<int,2>> supph;
	vector<int> supp;
	//vector<int> r1;//one ring neighborhood at the same level
	vector<int> tbf;//truncated basis functions fot type 2
	vector<double> tc;//truncation coefficients for type 2

	int truntmp;
	vector<int> tbftmp;
	vector<double> tctmp;

	vector<int> edge;//edges that connect to this vertex
	vector<int> face;//faces that share this vertex
	vector<int> hex;
	Vertex3D();
	bool operator==(const Vertex3D& v);

	vector<int> chd;
	vector<double> coef;

	//used for coupling meshfree
	int mpflag;//0 not assigned, 1 assigned
	int mpid;//ID of meshfree node

	int bzid;
	int c0flag;
	int c0flag_b;

	int smth;
	int pvt_id;
};

class Element2D
{
public:
	int cnct[4];
	//vector<vector<int>> edge;
	int edge[4];
	//vector<vector<int>> edge_act;
	int act;
	//int aff;
	int type;//0 for square, 1 for rectangular, 2 for non-corner boundary, 3 for corner boundary, 4 for extraordinary, 5 for invalid
	int lev;
	int prt;
	vector<int> chd;
	vector<array<double,2>> chd_o;
	//vector<array<int,2>> Tjunc;//[edge pos,id]
	vector<int> node;
	vector<ELCS2D> lcs;//local coordinate system for each node
	vector<int> IEN;
	vector<array<double,5>> patch_ku;
	vector<array<double,5>> patch_kv;
	vector<vector<double>> bemat;//for irregular element

	vector<int> IENtmp;
	vector<array<double,5>> patch_kutmp;
	vector<array<double,5>> patch_kvtmp;

	//vector<int> IENb;//global index of Bezier control points

	Element2D();
	bool operator==(const Element2D& e);
	//Element& operator=(const Element& e);
};

class Element3D//one element at most contains one extraordinary edge
{
public:
	int cnct[8];
	int edge[12];
	int face[6];
	int act;
	int id_act;
	int type;//0 for regular, 1 for boundary, 2 for extraordinary, 3 for T-edge
	int lev;
	int prt;
	int ref_flag;//0 for no refine, 1 for subdivide into 8, 21 for subdv into 4 in u direction, 22 (4,v), 23 (4,w), 31 for subdv into 2 in u, 32 (2,v), 33 (2,w)
	int trun;
	int ghost;
	int bc_lev;
	int bc;
	int smth;
	double Jacobian[8];
	double ConditionNumb[8];

	vector<int> chd;
	vector<array<double,3>> chd_o;
	//vector<array<int,2>> TjuncF;//[face pos, id]
	//vector<array<int,2>> TjuncE;//[edge pos, id]
	vector<int> IEN;//order in the desired way, in each hierarchy!
	vector<int> node;//all points, including T-junctions
	vector<Matrix4d> lcs;//local coordinate system for each node
	vector<array<double,5>> patch_ku;
	vector<array<double,5>> patch_kv;
	vector<array<double,5>> patch_kw;
	vector<vector<double>> bemat;//Bezier matrix, "C2" splines to C0 Bezier

	double kvlen[3];
	int fcnb[6][2];//[][0] for neighbor hex id, [][1] for face loc id in that hex
	int ednb[12][2];//for regular element only 
	int vtnb[8][2];//for regular element only
	double pkv[3][8];//patch knot vectors in three directions

	double dm[3][2];
	vector<vector<double>> tmat;
	vector<array<int, 2>> IEN_act;

	vector<array<double, 3>> bzpt;//tmp use only

	vector<int> IENtmp;
	vector<array<double, 5>> patch_kutmp;
	vector<array<double, 5>> patch_kvtmp;
	vector<array<double, 5>> patch_kwtmp;
	array<Matrix3d,6> nbrot;//rotation matrix for each face neighbor

	vector<int> IENb;//global index of Bezier control points
	vector<int> IENb_loc;//paired with IENb in transition element when only C0 Bezier added
	int bzflag;//0 for non-Bezier patch, 1 for Bezier
	vector<int> IENc01;//C0 and C1 splines
	//vector<int> IENc01_b;//boundary C0 and C1 splines
	//vector<int> IENc0_loc;
	//vector<int> IENc0;//C0 splines
	//vector<int> IENc1;//C1 splines
	vector<int> IENc2;//C2 splines with truncation
	vector<vector<double>> cmat;//Bezier extraction matrix, "C0", "C1" and "C2" splines to C0 Bezier

	Element3D();
	void Clear();
	void Initialize();

	int pbd_id;
	int jacobFlag;
};

class BezierElement2D
{
public:
	BezierElement2D(int p=3);
	int order;
	int degree;
	int nbf;
	int prt;
	vector<array<double,3>> pts;
	//double pts[16][3];
	//double pts4[25][3];
	vector<int> IEN;
	vector<vector<double>> cmat;
	//vector<array<double,16>> cmat;
	//vector<array<double,25>> cmat4;
	void BezierPolyn(double u, vector<double>& Nu, vector<double>& dNdu) const;
	void Basis(double u, double v, vector<double>& Nt, vector<array<double,2>>& dNdt) const;
	//void Basis(double u, double v, double Nx[16], double dNdx[16][2]) const;
	//void Basis4(double u, double v, double Nx[25], double dNdx[25][2]) const;
	void Para2Phys(double u, double v, double pt[3]);
	//void Para2Phys4(double u, double v, double pt[3]);
	void SurfPointNormal(double u, double v, array<double,3>& pt, array<double,3>& nm) const;
	//void SurfPointNormal4(double u, double v, array<double,3>& pt, array<double,3>& nm) const;
};

class BezierElement3D
{
public:
	BezierElement3D(int p=3);
	int degree;
	int order;
	int nbf;
	int type;//0 for interior and 1 for boundary, for visualization purpose
	//int prt;
	int prt[2];
	int trun;
	vector<int> bfc;
	vector<array<double,3>> pts;
	vector<int> IEN;
	vector<vector<double>> cmat;
	void BezierPolyn(double u, vector<double>& Nu, vector<double>& dNdu) const;
	void Basis(double u, double v, double w, vector<double>& Nt, vector<array<double,3>>& dNdt) const;
	void Para2Phys(double u, double v, double w, double pt[3]) const;
	void GeomMap(const array<double,3>& u, array<double,3>& pt) const;
	//void GeomMap(const array<double, 3>& u, array<double, 3>& pt);
	void SurfPointNormal(double u, double v, double w, array<double,3>& pt, array<double,3>& nm) const;

	vector<array<double, 3>> cnpt;

	//for meshfree
	int rkpm;
	vector<int> mid;
	vector<array<double, 3>> mp;
	vector<double> ra;
	int bc[6];//0 for interior non-coupling face, 1 for boundary face, 2 for coupling interface
	double bcval[6];

	//for Bezier coupling
	int bzcouple;//1 for coupling element, 0 for no-coupling
	int bzflag;//0 for spline element, 1 for Bezier element
	int bcflag;
	vector<int> IENb;
};


class BPatch3D
{
public:
	BPatch3D();
	int deg[3];
	int npt[3];
	vector<vector<double>> kv;
	vector<double> wght;
	vector<array<double, 3>> cp;
	vector<int> pid;
	vector<array<int, 3>> ele;
	vector<vector<int>> IEN;//local in the patch
	void Initialize(const vector<vector<double>>& kv_in);
	void BuildElement();
	void GlobalRefine(int dir);
	void BezierExtract(int eid, vector<vector<double>>& cmat, vector<int>& IENg);
};



class MFNode
{
public:
	MFNode();
	array<double, 3> coor;
	double a;
};

//class PointIPP
//{
//public:
//	int index[2];//index space
//	double pm[2];//parameter space
//	//double coor[3];//physical coordinates
//	bool operator==(const PointIPP& pt);
//};

//bool ComparePointIPPu(PointIPP ipp1, PointIPP ipp2)
//{
//	return (ipp1.pm[0] < ipp2.pm[0]);
//}
//
//bool ComparePointIPPv(PointIPP ipp1, PointIPP ipp2)
//{
//	return (ipp1.pm[1] < ipp2.pm[1]);
//}

class RegularPatchBasis
{
public:
	RegularPatchBasis(){for(int i=0;i<4;i++) {val[i]=0.; Dval[i]=0.; D2val[i]=0.;}};
	double val[4];
	double Dval[4];
	double D2val[4];
	void Evaluate(double x,double a=1.,double b=0.);
	void Clear();
};


//mesh format conversion

void Raw2Vtk_hex(string fn);

void Rawn2Vtk_hex(string fn);

void Raw2Vtk_hex_delete(string fn, string fnd);

void ReadVtk_hex(string fn, vector<array<double, 3>>& pts, vector<array<int, 8>>& cnct);

void WriteVtk_hex(string fn, const vector<array<double, 3>>& pts, const vector<array<int, 8>>& cnct);

void ReadRaw_hex(string fn, vector<array<double, 3>>& pts, vector<array<int, 8>>& cnct);

void WriteRaw_hex(string fn, const vector<array<double, 3>>& pts, const vector<array<int, 8>>& cnct);

void Vtk2Raw_hex(string fn);

void ReadRaw_tri(string fn, vector<array<double, 3>>& pts, vector<array<int, 3>>& cnct);

void WriteRaw_tri(string fn, const vector<array<double, 3>>& pts, const vector<array<int, 3>>& cnct);

void ReadVtk_tri(string fn, vector<array<double, 3>>& pts, vector<array<int, 3>>& cnct);

void WriteVtk_tri(string fn, const vector<array<double, 3>>& pts, const vector<array<int, 3>>& cnct);

void RepairConnect_tri(int npt, vector<array<int, 3>>& cnct, int ref);

#endif