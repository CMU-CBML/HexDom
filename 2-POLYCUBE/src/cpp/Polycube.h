#pragma once
#include <vector>
#include <fstream>
#include "rawmesh.h"

#include "Eigen/Dense"
#include "Eigen/LU"
#include "Eigen/SparseCore"
#include "Eigen/SparseCholesky"
#include "Eigen/IterativeLinearSolvers"

using namespace Eigen;
using namespace std;

const double EPSILON = 1e-6; // 1e-8 is too small sometimes

class PolycubeConstruction :
	public RawMesh
{
public:
	PolycubeConstruction(void);
	PolycubeConstruction(string inputFileName);

	~PolycubeConstruction(void);

	bool ReadKTri(const char * filename);
	
	static const int NUM_CLUSTER = 6;
	static const int NUM_NEI_CLUSTER = 12;
	static const int NUM_RING = 0;
	static const int INITIALIZATION_METHOD = 0;
	static const int INITIALIZATION_RANDOM = 0;

	static const int SET_INITIAL_SKELETONPOINTS = 0;

	static const double WEIGHT_LENGTH_EWCVT;
	
	static const int ASSIGN_PART_ELEMENT = 1;
	static const int ASSIGN_PART_ELEMENT_K = 1;

	static const int ASSIGN_CLUSTER_ELEMENT= 0; //Not used for k files

	static const int OUTPUT_CORNER_POINTS = 1;
	static const int READ_IN_PARAHEX_TORUS = 0;  //
	static const int READ_IN_UNIFORMHEX = 0;
	static const int READ_IN_MAPPING = 0;

	static const int OCTREE_MAX_LEVEL = 0;
	static const int OCTREE_MIN_LEVEL = 0;

	static const int CURVE_SMOOTH = 0;
	static const int CURVE_FITTING = 0;

	static const int READ_WRITE_K_BEFORE_POSTPROCESSING = 1; // READ == 1, WRITE == 0
	static const int READ_WRITE_K_INDEXPATCH = 1; // READ == 1, WRITE == 0
	///////////////////////////////////////////////////////////////////////////////////
	//For curve fitting
	struct XYZ
	{
		double x;
		double y;
		double z;
	};
	///////////////////////////////////////////////////////////////////////////////////

	struct CVTElement 
	{
		int indexCluster;
		int index;

		int indexPatch;//12/18/2014
		int indexSkeletonPoint;
		int indexPart;

		vector<int> directNei;
		int numDirectNei;

		vector<int> neiElements;
		int numNeiElements;

		int numNeiCluster;
		int numNeiElementEachCluster[NUM_NEI_CLUSTER];
		int indexNeiClusters[NUM_NEI_CLUSTER];

		double normal[3];
		double normalNew[3];
	};

	struct SkeletonPoint
	{
		int index;

		double localCoordX[3];
		double localCoordY[3];
		double localCoordZ[3];
		double transMatrix[3][3];
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
		vector<int> InterVert;
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

	string inputName;

	RawMesh curveSkeleton;
	vector<CVTElement> elementArray;
	vector<SkeletonPoint> skeletonPoint;
	vector<bool> skeletonPointUsed;
	vector<Centroid> generators;

	int numberCornerPoints;
	vector<int> cornerPoints;
	vector<Polycube> polycubePatch;
	int NUM_POLYCUBE_PATCH; // total number of separated polycube patches
	RawMesh* polycubePara;

	RawMesh parametricHex;
	RawMesh realDomainHex;
	vector<PropagationElement> propElement;

	RawMesh adaptiveOctreeHex;
	RawMesh adaptiveFinalHex;

	ofstream outputCornerPoints;

	/////////////////////////////////////////////
	//////FOR ADAPITVE HEX MESHING
	int nNode;
	int nElement;
	vector<Node> Nodes;
	vector<Element> Elements;
	////////////////////////////////////////////

	bool Initialization(void);
	bool InitializeSkeletonPoints(void);
	//bool InitializeElement(void);
	bool InitializeElement(const char * filename);
	bool InitializeSegments(void);
	bool SearchNei(void);
	//Edge-weighted CVT method
	bool EdgeWeightedCVT(void);

	//post processing
	bool PostProcessing(void);
	//Polycube Construction
	bool InitializePolycube(void);
	//Hex meshing for the parametric domain
	bool HexMeshParametricDomain(const char *outputName);
	//Hex meshing for the real domain
	bool HexMeshRealDomain(const char *outputName);

	//Hex meshing through octree subdivision
	bool HexMeshOctreeSubdivision(const char *outputName);

	bool OutputPatchesVTK(const char *outputName);
	bool OutputPatchesVTKBifurcation(const char *outputName);
	bool OutputPatchesVTKPara(const char *outputName);

	////////////////////////////////////////////////////////////////
	//Read and Write K files
	bool ReadKFileInitial(const char *inputName);//indexPart
	bool ReadKFileBeforePostProcessing(const char *inputName);//indexCluster
	bool WriteKFileBeforePostProcessing(const char *outputName);
	bool ReadKFileIndexPatch(const char *outputName); //indexPatch
	bool WriteKFileIndexPatch(const char *outputName); //indexPatch

	//bool ReadKFileMapping(const char *inputName);//indexPatch
	bool WriteKFileMapping(const char *outputName);

	bool ReadKFileHexMeshPara(const char *inputName);
	bool OutputExtraInformation(const char * outputmeshName);
	bool CreateInitialPolycube(const char *outputmeshName);
	//indexHexPart
	bool WriteKFileHexMeshPara(const char *outputName);

	bool ReadKFileHexMeshPhy(const char *inputName); //indexHexPart
	bool WriteKFileHexMeshPhy(const char *outputName);
	////////////////////////////////////////////////////////////////

	//For octree construction
	int levelRes[10];

private:

	bool RecursiveLocalCoordinateSystem(int start_vertex);
	bool CalculateLocalCoordinateSystem(double vz[3], double vxpre[3], double vx[3], double vy[3]);
	bool CalculateRotationMatrix(SkeletonPoint &sp);
	bool OutputLocalCoordinateSystem(void);

	

	bool AssignPartToElement();
	bool AssignClusterToElement();

	bool GetDirectNeighboringElementByElement(int elementID); // three neighboring elements for each triangle
	bool GetNeighboringElementInRings(int ringNumber = 1);
	double GetPhysicalDist(int index, const CVTElement &currentElement);

	bool InitializeGeneratorsByInput(void);
	double GetNormalDist(const Centroid &currentGenerator, const CVTElement & currentElement);
	double NormalDist(double veca[], double vecb[]);
	bool NormalizeGeneratorNormal(Centroid &currentGenerator);

	bool IsCounted(CVTElement &currentElement, int indexClusterNeiElement, int &position);

	bool IsBoundaryElement(CVTElement &currentElement);

	//For Edge-weighted CVT
	double GetEWDist(const Centroid &currentGenerator, const CVTElement & currentElement);
	int GetShortestEWDist(CVTElement &currentElement);
	bool DataTransfer(CVTElement &currentElement, int newIndex);

	//For post processing
	bool EnforceLabelConnectivity(void);
	bool CheckLabelConnectivity(void);
	bool EnforceBoundaryConnectivity(void);
	bool CheckBoundaryConnectivity(void);

	//For polycube construction
	bool SearchCornerandBoundary(void);
	bool IsCornerPoint(int vertexID);
	bool IsBoundaryPoint(int vertexID);
	bool IsBoundaryEdge(int vertexIDone, int vertexIDtwo);
	int FindOneCornerPoint(Polycube &currentPolycubePatch);
	int SearchCCWVertex(int vertexID, int elementID);

	bool ModifyWrongBoundaryElements(void);
	bool EdgeFlipTwoElements(int vertexID, int elementID);
	bool CornerMappingMaxMin(Polycube &currentPolycubePatch);
	bool ParametricMapping(void);
	bool ParametricMappingCorner(void);
	bool ParametricMappingCornerByInput(void);
	bool ParametricMappingEdge(void);
	bool ParametricMappingEdgeGeneral(void);
	bool ParametricMappingInterior(void);
	bool InteriorMappingOnePatch(Polycube &currentPolycubePatch);
	bool MappingOnePatchWeighted(Polycube &currentPolycubePatch);
	double TriAreaElement(int elementID);
	bool FlipCheckEachPatch(Polycube &currentPolycubePatch);
	bool SignedTriAreaElement(int elementID, int indexCluster);
	bool MappingOnePatchPostProcessing(Polycube &currentPolycubePatch);
	bool IsInteriorPoint(int vertexID);

	bool SmoothBoundaryCurve(void);
	bool CurveFittingBoundaryCurve(void);
	//For Curve Fitting
	void SplineKnots(vector<int> &u, int n, int t);
	void SplineCurve(const vector<XYZ> &inp, int n, const vector<int> &knots, int t, vector<XYZ> &outp, const vector<double> &ratio);
	void SplinePoint(const vector<int> &u, int n, int t, double v, const vector<XYZ> &control, XYZ &output);

	void SplineCurve(const vector<XYZ> &inp, int n, const vector<int> &knots, int t, vector<XYZ> &outp, int res);
	double SplineBlend(int k, int t, const vector<int> &u, double v);
	//End for curve fitting

	//Hex meshing for the real domain
	bool HexMeshProjectionBoundary(void);
	bool Projection(double paraPosition[3], double *realPosition);
	bool HexMeshProjectionInterior(void);
	bool FindPropagationBoundElements(int vertexID);
	bool Propagation(int vertexID);

	//Hex meshing through octree subdivision
	double TriArea(double v0[3], double v1[3], double v2[3]);
	double QuadArea(double v0[3], double v1[3], double v2[3], double v3[3], double *MassCenter);
	double GetHexVolume(double p[8][3], double *MassCenter);

	int GetLevel(int octree_id);
	int XYZToOctreeidx(int x, int y, int z, int level, int direction = 0);
	void OctreeidxToXYZ(int octree_id, int &x, int &y, int &z, int level);



	//by yu
	//bool CreateInitialPolycube(void);
};

