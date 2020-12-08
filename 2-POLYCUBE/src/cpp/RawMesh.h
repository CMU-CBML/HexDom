//**************************************************************
//**  Basic mesh class						                  **
//**  Developed By Xinghua Liang (liangxh@cmu.edu)            **
//**  Version: 0.1_2010-11-11					              **
//**************************************************************

#ifndef _RAWMESH_
#define _RAWMESH_

#define _CRT_SECURE_NO_DEPRECATE

#include <cstdio>
#include <cstring>
#include <vector>
#include <string>
#include <cmath>
#include <omp.h>
#include <iostream>
#include <fstream>
#include <sstream>


#include "Eigen/Dense"


using namespace std;
using namespace Eigen;

class RawMesh
{
public:
	//
	//Define a new structure for element property
	//
	typedef int ElementType;

	struct ElementProperty
	{
		//The type of element
		ElementType		elementType;
		//Number of vertries in one element;
		int				vertexNumber;
		//Number of edges in one element;
		int				edgeNumber;
		//Number of faces in one element;
		int				faceNumber;
		//Number of vertrices in each face
		int				faceVertexNumber;
		//Element type of the face
		ElementType		faceType;
	};

	struct MeshInfo
	{
		double minCorner[3];
		double maxCorner[3];
		double center[3];
		double biggestSize;
		double smallestSize;

		bool isUseVertexNormal;
		bool isUseVertexColor;
		bool isUseElementNormal;
		bool isUseElementColor;
	};

	//Constants in this class
	const ElementType	NA;
	const ElementType	POINT;
	const ElementType	LINE;
	const ElementType	TRIANGLE;
	const ElementType	TETRAHEDRON;
	const ElementType	QUADRILATERAL;
	const ElementType	HEXAHEDRON;
	const ElementType	HEXAGON;
	const ElementType	POLYGON;

	double const	MAX;
	double const	MIN;

	//*****************************
	//** Variabls for this class **
	//*****************************
	//==Essential components for vertex==
	//"vertex" is used to store the coordinates of all vertries. vertex[i][0] is the x-coordinate of i-th vertex
	double	(*vertex)[3];
	//"vertexNumber" is the total number of vertries in the mesh
	int		vertexNumber;
	//"vertexSign" is used to give more information for each vertex
	int		*vertexSign;
	//store the color for each vertex
	double	(*vertexColor)[3];
	//store the normal for each vertex
	double	(*vertexNormal)[3];
	//store the curvature for each vertex
	double	*vertexCurvature;


	//==Essential components for element==
	//"element" is used to store the connectivity of all elements. element[i][0] is the number of the 1st vertex of i-th vertex
	int				**element;
	//"elementNumber" is the total number of elements in the mesh
	int				elementNumber;
	//"elementSign" is used to give more information for each element
	int				*elementSign;
	//store the color for each element
	double			(*elementColor)[3];
	//store the normal for each element
	double			(*elementNormal)[3];
	//The property of element
	ElementProperty elementProperty;

	//==Relation between vertex, edge, face, element
	//edges connect to a vertex
	int		**edgeValence;
	int		*edgeValenceNumber;
	//elements using this vertex
	int		**elementValence;
	int		*elementValenceNumber;

	//==Mesh information
	int		meshID;
	string	meshName;
	string	version;
	MeshInfo	meshInfo;

	//==Some default value for color and normal
	double	defaultColor[3];
	double	defaultNormal[3];

	//==parallel threshold number
	int		thresholdForParallel;
	int		thresholdForParallel_Light;
	int		thresholdForParallel_Heavy;


	//******************************
	//** Functions for this class **
	//******************************
	RawMesh();
	RawMesh(const RawMesh& inMesh);
	virtual ~RawMesh();

	virtual RawMesh& operator = (RawMesh const &inMesh);

	//Mesh operations
	virtual void	AddVertex(int addedNumber, double (*addedVertex)[3], int *addedVertexSign=NULL);
	virtual void	AddVertex(int addedNumber, double (*addedVertex)[3], int *addedVertexSign, double (*addedVertexColor)[3]);
	virtual void	AddVertex(vector <double> &addedVertex, vector <int> &addedVertexSign);
	virtual void	AddVertex(vector <double> &addedVertex);
	virtual void	AddElement(int addedNumber, int **addedElement, int *addedElementSign=NULL);
	virtual void	AddElement(vector <int> &addedElement, vector <int> &addedElementSign);
	virtual void	AddElement(vector <int> &addedElement);

	virtual void	ApplyColorMapToVertex(int mapType=0);

	virtual void	ChangeElementDirection();
	virtual void	CopyMesh(RawMesh* &mesh);
	virtual bool	ComputeNormal(int normalType = 1);	//normalType: 1 - vertex normal; 2 - element normal.
	virtual bool	ComputeVertexNormal();
	virtual bool	ComputeElementNormal();
	virtual bool	ComputeVertexCurvature(int curvatureType = 1, bool isUseExtendQuadric=true);	//1: max curv; 2: min curv; 3: mean curv; 4: Guassian curv.
	virtual void	ClearOldMesh();
	virtual bool	CreateNewMesh(ElementType const myType, int const myVertexNumber, int const myElementNumber);

	virtual void	DeleteUnusedVertex();
	virtual void	DeleteElement(bool *badElement);
	virtual void	DeleteElement(vector <int> &deletedElement);
	virtual void	DeleteDuplicatedPoint(int startPoint=0);
	virtual void	DeleteDuplicatedElement(int startPoint=0);
	virtual void	DeleteDomainContainTheElement(int iElem);
	virtual void	DeleteDomainContainTheElement(int iElem, int loopNum);

	virtual void	 EmptyCurrentMesh();
	virtual RawMesh* ExtractSurface();
	virtual RawMesh* ExtractFrame();
	virtual RawMesh* ExtractCrossSection(const char axis, double value[3]);

	virtual	bool	GetElementValence()	{return InitiateElementValence();}	//Use InitiateElementValence() if possible
	virtual void	GetElementEdge(int elementID, int edgeID, int &point_1, int &point_2);
	virtual void	GetElementFace(int elementID, int faceID, int *vertexIndex);
	virtual bool	GetEdgeValence()	{return InitiateEdgeValence();}		//Use InitiateEdgeValence() if possible
	virtual void	GetMeshInfo();
	virtual void	GetNewIndexOrder(int elementID, int vertexID, int *newInd);
	virtual void	GetNewIndexOrder(int elementID, int faceIndex, int vertexID, int *newInd);
	virtual int		GetNeighboringElementByEdge(int elementID, int point_1, int point_2){return -1;};
	virtual int		GetNeighboringElementByFace(int elementID, int face[]);
	virtual ElementProperty GetElementProperty(const ElementType myElementType);

	virtual bool	InitiateElementValence();
	virtual bool	InitiateEdgeValence();
	virtual bool	IsBoundaryFace(int *vertexIndex);
	virtual bool	IsBoundaryFace(int elementID, int faceIndex);
	virtual	bool	IsEmpty();
	virtual bool	IsEdge(const int vertex_1, const int vertex_2);
	virtual bool	IsInElement(const int vertexID, const int elementID);
	virtual int		IsInSameElement(const int vertex_1, const int vertex_2);

	virtual	void	Move(double x, double y, double z);
	virtual void	Move(double pos[3]);
	virtual void	Merge(RawMesh &addedMesh);
	
	virtual void	PushVertexSign();
	virtual void	PushElementSign();
	virtual void	PushElementType();
	virtual void	PushAndSetElementType(ElementType iType);
	virtual void	PopVertexSign();
	virtual void	PopElementSign();
	virtual void	PopElementType();

	virtual void	ReplaceVertexInElement(int elementID, int originalVertexID, int replacedVertexID);
	virtual void	Rotate(double xRad, double yRad, double zRad, bool rotateAtCenter=false);
	virtual void	Rotate(double rad[3], bool rotateAtCenter=false);

	virtual void	SetElementType(const ElementType myElementType);
	virtual void	SetBoundaryVertexSign(int sign=1);
	virtual void	SetInnerVertexSign(int sign=0);

	virtual	void	Transform(int elementID, double zoom[3], double rotate[3], double move[3]);

	virtual	void	UpdateMesh();

	virtual void	Zoom(double xRatio, double yRatio, double zRatio, bool zoomAtCenter=false);
	virtual void	Zoom(double ratio[3], bool zoomAtCenter=false);	

	//File operations
	virtual bool	Read(const char * filename);
	virtual	bool	Write(const char * filename);

	//____________________________________________________________________________________________________
	// modified by Kangkang Hu

	void Smooth(int iLoop);
	void SmoothSurface(int iLoop);

	//void MoveToSurface(RawMesh &surface);

	//____________________________________________________________________________________________________


	//Functions for memory operations
	bool	InitiateMatrix(double ** &matrix, int const dim1, int const dim2, double initialValue=0.0);
	bool	InitiateMatrix(int ** &matrix, int const dim1, int const dim2, int initialValue=0);
	bool	InitiateMatrix(bool ** &matrix, int const dim1, int const dim2, bool initialValue=false);
	bool	InitiateMatrix(double * &matrix, int const dim, double initialValue=0.0);
	bool	InitiateMatrix(int * &matrix, int const dim, int initialValue=0);
	bool	InitiateMatrix(bool * &matrix, int const dim, bool initialValue=false);
	bool	InitiateMatrix(double (* &matrix)[3], int const dim, double initialValue=0.0);

	void	FreeMatrix(double ** &matrix);
	void	FreeMatrix(int ** &matrix);
	void	FreeMatrix(bool ** &matrix);
	void	FreeMatrix(double * &matrix);
	void	FreeMatrix(int * &matrix);
	void	FreeMatrix(bool * &matrix);
	void	FreeMatrix(double (* &matrix)[3]);
	//End

	//Testing functions
	void	SortElementByZ(void);
	void	SortElementByZ(int start, int end);
	//End

protected:
	vector <int> copyOfVertexSign;
	vector <int> copyOfElementSign;
	vector <ElementType> copyOfElementType;

	//Functions for file operation
	bool	Read_RawMesh(const char * filename);
	bool	Read_Mesh(const char * filename);
	bool	Read_RawFile(const char * filename);
	bool	Read_Ply(const char * filename);
	bool	Read_Inp(const char * filename);
	bool	Read_Plt(const char * filename);
	bool	Read_Vtk(const char * filename);
	
	//Math function
	void	CrossProduct(double vector_1[3], double vector_2[3], double result[3]);
	void	CrossProduct(double point_1[3], double point_2[3], double point_3[3], double result[3]);
	void	Normalize(double vector[3]);

	//For Debug and Error
	void	OutputError(const char *className, const char *message);
	void	OutputError(const char *className, const char *message, int index);

private:
	bool	ReadRaw_tet(const char *filename, const char *type);
	bool	ReadRaw_tri(const char *filename, const char *type);
	bool	ReadRaw_quad(const char *filename, const char *type);
	bool	ReadRaw_hex(const char *filename, const char *type);
	bool	ReadRaw_line(const char *filename, const char *type);
	bool	ReadRaw_ReadVertex(FILE *&input, const char * type, bool withSign = true);
	bool	Read_RawMesh_OldVersion(const char * filename);

	bool	Write_RawMesh(const char *filename);
	bool	Write_RawMesh_OldVersion(const char *filename);
	bool	Write_Mesh(const char *filename);
	bool	Write_RawFile(const char *filename);

	bool	Write_Ply(const char *filename);
	bool	Write_Plt(const char *filename);
	bool	Write_Inp(const char *filename);
	bool	Write_Vtk(const char *filename);

	bool	WriteRaw_tet(const char *filename, const char *type);
	bool	WriteRaw_tri(const char *filename, const char *type);
	bool	WriteRaw_quad(const char *filename, const char *type);
	bool	WriteRaw_hex(const char *filename, const char *type);
	bool	WriteRaw_line(const char *filename, const char *type);
	bool	WriteRaw_WriteVertex(FILE *&output, const char * type, bool withSign = true);

	//Function for Normal operation
	bool	ComputeVertexNormal_Point();
	bool	ComputeVertexNormal_Line();
	bool	ComputeVertexNormal_Tri();
	bool	ComputeVertexNormal_Quad();
	bool	ComputeVertexNormal_Tet();
	bool	ComputeVertexNormal_Hex();

	bool	ComputeElementNormal_Point();
	bool	ComputeElementNormal_Line();
	bool	ComputeElementNormal_Tri();
	bool	ComputeElementNormal_Quad();

	bool	InitiateElementValenceMember();
	bool	InitiateEdgeValenceMember(const int &sumNum);

	void	FreeElementValence();
	void	FreeEdgeValence();


	

};

#endif _RAWMESH_
