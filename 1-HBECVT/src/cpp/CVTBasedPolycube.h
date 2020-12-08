#pragma once
#include <vector>
#include <fstream>
#include "rawmesh.h"
#include "GeometricMethod.h"
#include <algorithm>
#include <sstream>
#include <iostream>
#include "Eigen/Dense"
#include "Eigen/LU"
#include "Eigen/SparseCore"
#include "Eigen/SparseCholesky"
#include "Eigen/IterativeLinearSolvers"

using namespace Eigen;
using namespace std;

const double EPSILON = 1e-6;



class CVTBasedPolycube :
	public RawMesh
{

public:
	CVTBasedPolycube(void);
	CVTBasedPolycube(string inputFileName);
	~CVTBasedPolycube(void);

	bool ReadKTri(const char * filename);

	static const int NUM_CLUSTER = 6;
	//static const int NUM_DIRECT_NEI = 20;
	
	static const int NUM_NEI_CLUSTER = 10; // why the value is larger than NUM_CLUSTER?
	
	static const int NUM_RING = 2;

	static const int INITIALIZATION_METHOD = 0;

	double WEIGHT_LENGTH_EWCVT;

	static const int OCTREE_MAX_LEVEL = 5;//5 before

	static const int OCTREE_MIN_LEVEL = 3;//3 before

	//static const double OCTREE_TOL_ERR;

	double OCTREE_TOL_ERR;

	static const double DOMAIN_SIZE;

	static const int MANUAL_MODIFICATION = 0;

	static const int READ_IN_MAPPING = 0;//1 before

	static const int READ_IN_PARAHEX_TORUS = 1;//1 before

	static const int READ_IN_UNIFORMHEX = 0;

	static const int HAS_SPECIAL_ELEMENT = 1;//0 before

	static const int READ_IN_CORNERS = 0;//0 before

	static const int MARECHAL_METHOD = 1;//1 before

	static const int CURVE_SMOOTH = 0;

	static const int CURVE_FITTING = 0;

	struct CVTElement 
	{

		int indexCluster;
		int index;

		int indexPatch;//12/18/2014

		vector<int> directNei;
		int numDirectNei;

		vector<int> neiElements;
		int numNeiElements;

		int numNeiCluster;
		int numNeiElementEachCluster[NUM_NEI_CLUSTER];
		int indexNeiClusters[NUM_NEI_CLUSTER];

		double normal[3];

	};

	struct Centroid
	{
		int index;

		int numElements;

		double normal[3];

	};

	struct Polycube
	{

		int index;
		int indexCluster;

		int numCorner;
		vector<int> cornerPoint;

		vector<vector<int> > boundaryEdge;
		
		////////////////////////////////////////////Parametric mapping step
		int numIntVert;
		vector<int> InterVert;//interior vertices
		///////////////////////////////////////////

		int numElements;
		vector<int> element;

	};

	//For propagation in hex meshing step
	struct PropagationElement
	{
		int px;
		int nx;
		int py;
		int ny;
		int pz;
		int nz;
	};

	////////////////////////////////////////////////////////////////////////////////////
	//For adaptive hex meshing

	struct	Edge
	{
		int		NodeIndex;
		int		Twin;
	};

	struct Node
	{
		double	Coords[3];
		bool		BoundaryFlag;
		//char		Type; // 0: Regular Node; 1:Face TJunction; 2:Edge TJunction; 3:Partial Extraordinary Node; 4: Extraordinary Node; 5: Newly Inserted ExtraNodes
		//double	KnotVector[3][5];
		//double	KnotInterval[3][4];
		//unsigned int ParaPos[3];
		double ParaPos[3];
		vector<int> NeighborElement;
		//vector<Edge> NeighborEdge;
	};

	struct Junction
	{
		char	Position;
		int		NodeIndex;
	};

	struct ExtraNode
	{
		int		NodeIndex;
		unsigned char Pos[2]; // Pos[0] indicate the node which the extranode associated with; Pos[1] indicate the position spot around that node;
	};

	struct	NodeInfor
	{
		int		NodeIndex;
		double	UVWknots[3];
		char	RelationFlag[3];
	};

	struct Quad
	{
		int NodeIndex[4];
	};

	struct Element
	{
		int		NodeIndex[8];
		//unsigned int	EdgeInterval[3];
		//char		OctreeLevel;
		int OctreeID;
		bool		BoundaryFlag;
		vector<Junction> Junctions;
		//vector<ExtraNode> ExtraNodes;
		//vector<NodeInfor>	NodesList;  // RelationFlag shows the knot vector relation with the knot interval;
	};
	////////////////////////////////////////////////////////////////////////////////////
	//End for adaptive hex meshing

	///////////////////////////////////////////////////////////////////////////////////
	//For curve fitting
	struct XYZ
	{
		double x;
		double y;
		double z;
	};
	///////////////////////////////////////////////////////////////////////////////////

	string inputName;

	vector<CVTElement> elementArray;
	vector<Centroid> generators;

	int numberCornerPoints;
	vector<int> cornerPoints;

	RawMesh* polycubePara;

	RawMesh parametricHex;
	RawMesh realDomainHex;

	RawMesh adaptiveOctreeHex;

	RawMesh adaptiveFinalHex;

	vector<PropagationElement> propElement;

	GeometricMethod GM;

	vector<Polycube> polycubePatch;
	int NUM_POLYCUBE_PATCH; // total number of separated polycube patches

	/////////////////////////////////////////////
	//////FOR ADAPITVE HEX MESHING
	int nNode;
	int nElement;
	vector<Node> Nodes;
	vector<Element> Elements;
	////////////////////////////////////////////


	bool Initialization(const char *inputName);
	
	bool InitializeElement(void);
	bool InitializeSegments(void);
	bool SearchNei(void);

	//classic CVT method

	bool ClassicalCVT(void);

	//Edge-weighted CVT method

	bool EdgeWeightedCVT(void);

	//post processing
	bool PostProcessing(const char *inputManualName);

	//Polycube Construction
	bool InitializePolycube(void);

	//Hex meshing for the parametric domain
	bool HexMeshParametricDomain(const char *outputName);

	//Hex meshing for the real domain
	bool HexMeshRealDomain(const char *outputName);

	//Adaptive Hex meshing for both parametric domain and real physical domain
	bool AdaptiveHexMeshing(const char *outputName);

	//Output results

	ofstream outputCentroids;
	ofstream outputCornerPoints;
	ofstream outputSpecialElements;

	bool OutputResults(void);
	bool OutputPatchesVTK(const char *outputName);
	bool WriteKFileBeforePostProcessing(const char * outputName);
	bool OutputPatchesVTKPara(const char *outputName);

	//me added
	bool OutputPatchesVTK_PatchID(const char *outputName);
	bool ModifySpecialElement(const char *inputManualName);

private:

	///////////////////////////////////////////////////////////////////////////
	//Parameters for octree construction

	double cellSize;
	double origCood[3];
	int gridSize;
	int leafNum;
	int levelRes[10];
	int numVoxels;
	int numGrids;
	int octreeDepth;
	int octreeCellNum;
	int voxelSize;
	vector<char> refineFlagArray;
	vector<int> cutArray;

	//restore the information of hex meshes
	//vector<vector<double> > hexVerts;
	vector<Node> hexVerts; //Restore both Para and Phys information
	vector<int> hexIdx;
	vector<vector<int> > hexIdxM; // Marechal's method

	int numVerts;
	int numElems;
	vector<int> vtxIdxArray;
	vector<int> dualVertexOnce;

	//////////////////////////////////////////////////////////////////////////

	bool GetNeighboringElementByElement(int elementID); // more than three

	bool GetDirectNeighboringElementByElement(int elementID); // three neighboring elements for each triangle

	bool GetNeighboringElementInRings(int ringNumber = 1);

	bool InitializeGeneratorsByInput(void);
	bool InitializeGeneratorsByRandom(void);
	//bool ReadInputMesh(const char *inputName);

	double GetNormalDist(const Centroid &currentGenerator, const CVTElement & currentElement);
	double NormalDist(double veca[], double vecb[]);
	double GetEWDist(const Centroid &currentGenerator, const CVTElement & currentElement);

	int GetShortestNormalDist(CVTElement &currentElement);
	int GetShortestEWDist(CVTElement &currentElement);

	bool DataTransfer(CVTElement &currentElement, int newIndex);
	bool UpdateGenerator(Centroid &currentGenerator, CVTElement &currentElement, int inOut);

	bool NormalizeGeneratorNormal(Centroid &currentGenerator);
	bool RecalculateCentroids(void);

	bool IsCounted(CVTElement &currentElement, int indexClusterNeiElement, int &position);

	bool IsBoundaryElement(CVTElement &currentElement);

	//For post processing
	bool EnforceLabelConnectivity(void);
	bool CheckLabelConnectivity(void);

	bool EnforceBoundaryConnectivity(void);

	

	bool SmoothBoundaryCurve(void);
	bool CurveFittingBoundaryCurve(void);
	//For Curve Fitting
	void SplineKnots(vector<int> &u, int n, int t);
	void SplineCurve(const vector<XYZ> &inp, int n, const vector<int> &knots, int t, vector<XYZ> &outp, const vector<double> &ratio);
	void SplinePoint(const vector<int> &u, int n, int t, double v, const vector<XYZ> &control, XYZ &output);

	void SplineCurve(const vector<XYZ> &inp, int n, const vector<int> &knots, int t, vector<XYZ> &outp, int res);
	//void SplinePoint(const vector<int> &u, int n, int t, double v, const vector<XYZ> &control, vector<XYZ> &output);

	double SplineBlend(int k, int t, const vector<int> &u, double v);
	//End for curve fitting

	bool ModifyWrongBoundaryElements(void);

	bool EdgeFlipTwoElements(int vertexID, int elementID);

	bool CheckBoundaryConnectivity(void);

	bool CheckPolycubeValidity(void);

	//For polycube construction
	bool SearchCornerandBoundary(void);
	
	bool ParametricMapping(void);

	bool ParametricMappingCorner(void);
	bool ParametricMappingCornerByInput(void);
	bool ParametricMappingCornerByInput_ID(void);//Xiaodong

	bool CornerMappingMaxMin(Polycube &currentPolycubePatch);
	bool CornerMappingMaxMin_ManualAdjust(Polycube &currentPolycubePatch);//Xiaodong added
	bool RegularizeCornerPoints(void);

	bool ParametricMappingEdge(void);

	bool ParametricMappingInterior(void);
	bool InteriorMappingOnePatch(Polycube &currentPolycubePatch);

	bool MappingOnePatchWeighted(Polycube &currentPolycubePatch);
	bool MappingOnePatchPostProcessing(Polycube &currentPolycubePatch);
	bool FlipCheckEachPatch(Polycube &currentPolycubePatch);

	//double SignedTriAreaElement(int elementID, int indexCluster);
	bool SignedTriAreaElement(int elementID, int indexCluster);
	double TriAreaElement(int elementID);

	bool IsCornerPoint(int vertexID);
	bool IsBoundaryPoint(int vertexID);
	bool IsInteriorPoint(int vertexID);

	bool IsBoundaryEdge(int vertexIDone, int vertexIDtwo);

	int FindOneCornerPoint(Polycube &currentPolycubePatch);

	int SearchCCWVertex(int vertexID, int elementID);

	//Hex meshing for the parametric domain

	bool InitializeOctree(void);
	bool ConstructeOctree(void);
	//bool HexMeshOctree(const char* file_name);
	bool HexMeshOctree(void);

	int GetDepth(int res);
	int GetOctreeNum(int depth);
	int GetLevel(int octree_id);
	bool CheckAdaptation(int octree_id);
	void RefineBrothers(int octree_id, int *octree_idx);
	int XYZToOctreeidx(int x, int y, int z, int level, int direction = 0);
	void OctreeidxToXYZ(int octree_id, int &x, int &y, int &z, int level);
	int Child(int octree_id, int level, int i);

	bool DeleteOutsideMesh(const char *outputName);
	bool IsOutOfSurface(double dP_[3], RawMesh &surfaceMesh, int iTime=0);
	int	IntersectAxis(double *dP1_, double *dP2_, double *dP3_, double *dObject_, double *dIntersectP_, char cAxis_='x');

	bool CheckUniformHexValidity(void);

	//Hex meshing for the real domain

	bool HexMeshProjectionBoundary(void);
	bool Projection(double paraPosition[3], double *realPosition);

	bool HexMeshProjectionInterior(void);
	bool FindPropagationBoundElements(int vertexID);
	bool Propagation(int vertexID);

	bool Pillowing(void); // Only for paraHex

	//Adaptive Hex meshing process

	bool AdaptiveOctreeConstruction(void);

	double TriArea(double v0[3], double v1[3], double v2[3]);
	double QuadArea(double v0[3], double v1[3], double v2[3], double v3[3], double *MassCenter);
	double GetHexVolume(double p[8][3], double *MassCenter);

	bool AdaptiveHexMeshExtraction(const char *outputName);
	bool AdaptiveHexMeshExtractionMarechal(const char *outputName);

	bool CalculateCenterPoint(const Element currentElement, Node &tempNode);

	bool IsLeafCell(int x, int y, int z, int level, int direction = 0);
	int IsRegularNode(int x, int y, int z, int level);
	int IsIrregularNode(int x, int y, int z, int level, vector<int> &eight_cells);
	int IsSharedByEightCells(int x, int y, int z, int level, vector<int> &vec_eight_cell);

	int IsTransitionNode(int x, int y, int z, int level, int direction = 0);

	void TwentyOneVertices(int octree_id, double *vertex_cood, int i, int direction = 0);

	void AddOneHexRegular(int x, int y, int z, int level, int num_hexa);
	void AddOneHexIrregular(int x, int y, int z, int level, int num_hexa, vector<int> &eight_cells);

	void AddThirteenHexTypeOneFirst(int x, int y, int z, int level, int num_hexa, int direction = 0);
	void AddThirteenHexTypeOneSecond(int x, int y, int z, int level, int num_hexa, int direction = 0);
	void AddThirteenHexTypeOneThird(int x, int y, int z, int level, int num_hexa, int direction = 0);
	void AddThirteenHexTypeOneFourth(int x, int y, int z, int level, int num_hexa, int direction = 0);
	void AddThirteenHexTypeOneFifth(int x, int y, int z, int level, int num_hexa, int direction = 0);
	void AddThirteenHexTypeOneSixth(int x, int y, int z, int level, int num_hexa, int direction = 0);
	void AddThirteenHexTypeOneSeventh(int x, int y, int z, int level, int num_hexa, int direction = 0);
	void AddThirteenHexTypeOneEighth(int x, int y, int z, int level, int num_hexa, int direction = 0);

	void AddFiveHexTypeOneFirst(int x, int y, int z, int level, int num_hexa, int direction = 0);
	void AddFourHexTypeOneSecond(int x, int y, int z, int level, int num_hexa, int direction = 0);
	void AddFourHexTypeOneThird(int x, int y, int z, int level, int num_hexa, int direction = 0);
	void AddFiveHexTypeOneFourth(int x, int y, int z, int level, int num_hexa, int direction = 0);
	void AddFiveHexTypeOneFifth(int x, int y, int z, int level, int num_hexa, int direction = 0);
	void AddFourHexTypeOneSixth(int x, int y, int z, int level, int num_hexa, int direction = 0);
	void AddFourHexTypeOneSeventh(int x, int y, int z, int level, int num_hexa, int direction = 0);
	void AddFourHexTypeOneEighth(int x, int y, int z, int level, int num_hexa, int direction = 0);
	void AddFourHexTypeOneNinth(int x, int y, int z, int level, int num_hexa, int direction = 0);
	void AddFourHexTypeOneTenth(int x, int y, int z, int level, int num_hexa, int direction = 0);
	void AddFiveHexTypeOneEleventh(int x, int y, int z, int level, int num_hexa, int direction = 0);
	void AddFiveHexTypeOneTwelfth(int x, int y, int z, int level, int num_hexa, int direction = 0);
	void AddFiveHexTypeOneThirteenth(int x, int y, int z, int level, int num_hexa, int direction = 0);
	void AddFourHexTypeOneFourteenth(int x, int y, int z, int level, int num_hexa, int direction = 0);

	//Marechal's method
	
	int IsTransitionNodeMarechal(int x, int y, int z, int level, int direction = 0);

	void TwentyOneVerticesMarechal(int octree_id, double *vertex_cood, int i, int direction = 0, int para_phy = 0);

	void AddOneHexRegularMarechal(int x, int y, int z, int level, int num_hexa);
	void AddOneHexIrregularMarechal(int x, int y, int z, int level, int num_hexa, vector<int> &eight_cells);

	void AddThirteenHexTypeOneFirstMarechal(int x, int y, int z, int level, int num_hexa, int direction = 0);
	void AddThirteenHexTypeOneSecondMarechal(int x, int y, int z, int level, int num_hexa, int direction = 0);
	void AddThirteenHexTypeOneThirdMarechal(int x, int y, int z, int level, int num_hexa, int direction = 0);
	void AddThirteenHexTypeOneFourthMarechal(int x, int y, int z, int level, int num_hexa, int direction = 0);
	void AddThirteenHexTypeOneFifthMarechal(int x, int y, int z, int level, int num_hexa, int direction = 0);
	void AddThirteenHexTypeOneSixthMarechal(int x, int y, int z, int level, int num_hexa, int direction = 0);
	void AddThirteenHexTypeOneSeventhMarechal(int x, int y, int z, int level, int num_hexa, int direction = 0);
	void AddThirteenHexTypeOneEighthMarechal(int x, int y, int z, int level, int num_hexa, int direction = 0);

	void AddFiveHexTypeOneFirstMarechal(int x, int y, int z, int level, int num_hexa, int direction = 0);
	void AddFourHexTypeOneSecondMarechal(int x, int y, int z, int level, int num_hexa, int direction = 0);
	void AddFourHexTypeOneThirdMarechal(int x, int y, int z, int level, int num_hexa, int direction = 0);
	void AddFiveHexTypeOneFourthMarechal(int x, int y, int z, int level, int num_hexa, int direction = 0);
	void AddFiveHexTypeOneFifthMarechal(int x, int y, int z, int level, int num_hexa, int direction = 0);
	void AddFourHexTypeOneSixthMarechal(int x, int y, int z, int level, int num_hexa, int direction = 0);
	void AddFourHexTypeOneSeventhMarechal(int x, int y, int z, int level, int num_hexa, int direction = 0);
	void AddFourHexTypeOneEighthMarechal(int x, int y, int z, int level, int num_hexa, int direction = 0);
	void AddFourHexTypeOneNinthMarechal(int x, int y, int z, int level, int num_hexa, int direction = 0);
	void AddFourHexTypeOneTenthMarechal(int x, int y, int z, int level, int num_hexa, int direction = 0);
	void AddFiveHexTypeOneEleventhMarechal(int x, int y, int z, int level, int num_hexa, int direction = 0);
	void AddFiveHexTypeOneTwelfthMarechal(int x, int y, int z, int level, int num_hexa, int direction = 0);
	void AddFiveHexTypeOneThirteenthMarechal(int x, int y, int z, int level, int num_hexa, int direction = 0);
	void AddFourHexTypeOneFourteenthMarechal(int x, int y, int z, int level, int num_hexa, int direction = 0);


};

bool AddProgress(double progress);

