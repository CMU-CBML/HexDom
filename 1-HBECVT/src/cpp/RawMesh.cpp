#include "RawMesh.h"

#define _CRT_SECURE_NO_DEPRECATE


RawMesh::RawMesh(): meshName("Default.rawm"), MAX(1.0E+40), MIN(-1.0E+40),
					NA(0), POINT(1), LINE(2), TRIANGLE(3), TETRAHEDRON(4),
					QUADRILATERAL(5), HEXAHEDRON(6), HEXAGON(7), POLYGON(8)
{
	//Initial values for vertex components
	vertexNumber = 0;
	vertex = NULL;
	vertexSign = NULL;
	vertexNormal = NULL;
	vertexColor = NULL;
	vertexCurvature = NULL;
	//boundaryFace = NULL;

	//Initial values for element components
	elementNumber = 0;
	element = NULL;
	elementSign = NULL;
	elementColor = NULL;
	elementNormal = NULL;
	SetElementType(NA);
	
	//Initial values for other components
	edgeValence = NULL;
	edgeValenceNumber = NULL;
	elementValence = NULL;
	elementValenceNumber = NULL;
	
	//Initial values for default color and normal
	defaultColor[0] = 1.0f;
	defaultColor[1] = 0.7f;
	defaultColor[2] = 0.0f;
	defaultNormal[0] = 0.0f;
	defaultNormal[1] = 0.0f;
	defaultNormal[2] = 0.0f;

	meshInfo.biggestSize = 1.0f;
	meshInfo.center[0] = 0.0f;
	meshInfo.center[1] = 0.0f;
	meshInfo.center[2] = 0.0f;
	meshInfo.maxCorner[0] = 0.0f;
	meshInfo.maxCorner[1] = 0.0f;
	meshInfo.maxCorner[2] = 0.0f;
	meshInfo.minCorner[0] = 0.0f;	
	meshInfo.minCorner[1] = 0.0f;
	meshInfo.minCorner[2] = 0.0f;
	meshInfo.isUseVertexNormal = false;
	meshInfo.isUseVertexColor = false;
	meshInfo.isUseElementNormal = false;
	meshInfo.isUseElementColor = false;

	version = "2.0.20120515";

	thresholdForParallel = 10000;
	thresholdForParallel_Light = 1000000;
	thresholdForParallel_Heavy = 100;

	copyOfElementType.reserve(3);
}


RawMesh::~RawMesh()
{
	ClearOldMesh();

	meshName.clear();
	version.clear();
}

RawMesh::RawMesh(const RawMesh& inMesh): meshName("Default.rawm"), MAX(1.0E+40), MIN(-1.0E+40),
										 NA(0), POINT(1), LINE(2), TRIANGLE(3), TETRAHEDRON(4),
										 QUADRILATERAL(5), HEXAHEDRON(6), HEXAGON(7), POLYGON(8)
{
	if (this == &inMesh)
		return;

	int i, j;

	this->ClearOldMesh();
	this->SetElementType(inMesh.elementProperty.elementType);

	this->vertexNumber = inMesh.vertexNumber;
	this->elementNumber = inMesh.elementNumber;
	this->meshName = inMesh.meshName;

	if (inMesh.vertexNumber > 0)
	{
		this->InitiateMatrix(this->vertex, inMesh.vertexNumber);
		this->InitiateMatrix(this->vertexSign, inMesh.vertexNumber);

		for (i=0; i<inMesh.vertexNumber; ++i)
		{
			this->vertex[i][0] = inMesh.vertex[i][0];
			this->vertex[i][1] = inMesh.vertex[i][1];
			this->vertex[i][2] = inMesh.vertex[i][2];
			this->vertexSign[i] = inMesh.vertexSign[i];
		}

		if (inMesh.vertexColor != NULL)	//Not finished
		{
			this->InitiateMatrix(this->vertexColor, inMesh.vertexNumber);

			for (i=0; i<inMesh.vertexNumber; ++i)
			{
				this->vertexColor[i][0] = inMesh.vertexColor[i][0];
				this->vertexColor[i][1] = inMesh.vertexColor[i][1];
				this->vertexColor[i][2] = inMesh.vertexColor[i][2];
			}
		}
	}

	if (inMesh.elementNumber > 0)
	{
		int elementVertexNumber = inMesh.elementProperty.vertexNumber;

		this->InitiateMatrix(this->element, inMesh.elementNumber, elementVertexNumber);
		this->InitiateMatrix(this->elementSign, inMesh.elementNumber);

		for (i=0; i<inMesh.elementNumber; ++i)
		{
			for (j=0; j<elementVertexNumber; ++j)
				this->element[i][j] = inMesh.element[i][j];
			this->elementSign[i] = inMesh.elementSign[i];
		}

		if (inMesh.elementColor != NULL)
		{
			this->InitiateMatrix(this->elementColor, inMesh.elementNumber);

			for (i=0; i<inMesh.elementNumber; ++i)
			{
				this->elementColor[i][0] = inMesh.elementColor[i][0];
				this->elementColor[i][1] = inMesh.elementColor[i][1];
				this->elementColor[i][2] = inMesh.elementColor[i][2];
			}
		}
	}

	return;
}

RawMesh& RawMesh::operator = (RawMesh const &inMesh)
{
	if (this == &inMesh)
		return *this;

	int i, j;

	this->ClearOldMesh();
	this->SetElementType(inMesh.elementProperty.elementType);

	this->vertexNumber = inMesh.vertexNumber;
	this->elementNumber = inMesh.elementNumber;
	this->meshName = inMesh.meshName;

	if (inMesh.vertexNumber > 0)
	{
		this->InitiateMatrix(this->vertex, inMesh.vertexNumber);
		this->InitiateMatrix(this->vertexSign, inMesh.vertexNumber);

		for (i=0; i<inMesh.vertexNumber; ++i)
		{
			this->vertex[i][0] = inMesh.vertex[i][0];
			this->vertex[i][1] = inMesh.vertex[i][1];
			this->vertex[i][2] = inMesh.vertex[i][2];
			this->vertexSign[i] = inMesh.vertexSign[i];
		}

		if (inMesh.vertexColor != NULL)
		{
			this->InitiateMatrix(this->vertexColor, inMesh.vertexNumber);

			for (i=0; i<inMesh.vertexNumber; ++i)
			{
				this->vertexColor[i][0] = inMesh.vertexColor[i][0];
				this->vertexColor[i][1] = inMesh.vertexColor[i][1];
				this->vertexColor[i][2] = inMesh.vertexColor[i][2];
			}
		}
	}

	if (inMesh.elementNumber > 0)
	{
		int elementVertexNumber = inMesh.elementProperty.vertexNumber;

		this->InitiateMatrix(this->element, inMesh.elementNumber, elementVertexNumber);
		this->InitiateMatrix(this->elementSign, inMesh.elementNumber);

		for (i=0; i<inMesh.elementNumber; ++i)
		{
			for (j=0; j<elementVertexNumber; ++j)
				this->element[i][j] = inMesh.element[i][j];
			this->elementSign[i] = inMesh.elementSign[i];
		}

		if (inMesh.elementColor != NULL)
		{
			this->InitiateMatrix(this->elementColor, inMesh.elementNumber);

			for (i=0; i<inMesh.elementNumber; ++i)
			{
				this->elementColor[i][0] = inMesh.elementColor[i][0];
				this->elementColor[i][1] = inMesh.elementColor[i][1];
				this->elementColor[i][2] = inMesh.elementColor[i][2];
			}
		}
	}

	return *this;
}

void RawMesh::ChangeElementDirection()
{
	int iTmp;
	if (elementProperty.elementType == TRIANGLE)
	{
		for (int i=0; i<elementNumber; ++i)
		{
			iTmp = element[i][1];
			element[i][1] = element[i][2];
			element[i][2] = iTmp;
		}
	}
	else if (elementProperty.elementType == QUADRILATERAL)
	{
		for (int i=0; i<elementNumber; ++i)
		{
			iTmp = element[i][1];
			element[i][1] = element[i][3];
			element[i][3] = iTmp;
		}
	}
	else
	{
		OutputError("RawMesh::ChangeElementDirection", "No implementation for this element type");
	}

	return;
}

void RawMesh::CopyMesh(RawMesh* &mesh)
{
	if (mesh != NULL)
	{
		mesh->EmptyCurrentMesh();
		//delete mesh;
	}

	mesh = new RawMesh;
	(*mesh) = (*this);
}


inline bool	RawMesh::IsEmpty()
{
	if (vertexNumber == 0 || vertex == NULL || elementNumber == 0 || element == NULL)
		return true;

	return false;
}

inline bool RawMesh::IsInElement(const int vertexID, const int elementID)
{
	for (int i=0; i<elementProperty.vertexNumber; ++i)
		if (element[elementID][i] == vertexID)
			return true;

	return false;
}

int RawMesh::IsInSameElement(const int vertex_1, const int vertex_2)
{
	if (elementValence == NULL)
		GetElementValence();

	int i, j;
	for (i=0; i<elementValenceNumber[vertex_1]; ++i)
		for (j=0; j<elementValenceNumber[vertex_2]; ++j)
			if (elementValence[vertex_1][i] == elementValence[vertex_2][j])
				return elementValence[vertex_1][i];

	return -1;
}

bool RawMesh::IsBoundaryFace(int *vertexIndex)
{
	if (elementValence == NULL)
		GetElementValence();

	int i, iValenceNum, iIndex, j, k, iCount, jCount, iElem, faceVertexNumber, elementVertexNumber;
	iIndex = vertexIndex[0];
	iValenceNum = elementValenceNumber[iIndex];
	elementVertexNumber = elementProperty.vertexNumber;
	if (elementProperty.faceType == TRIANGLE)
		faceVertexNumber = 3;
	else if (elementProperty.faceType == QUADRILATERAL)
		faceVertexNumber = 4;
	else if (elementProperty.faceType == HEXAGON)
		faceVertexNumber = 6;
	else
		faceVertexNumber = GetElementProperty(elementProperty.faceType).vertexNumber;

	for (i=0; i<faceVertexNumber; ++i)
	{
		if (vertexIndex[i] < 0)
			return false;
	}

	iCount = 0;
	for(i=0; i<iValenceNum; ++i)
	{
		iElem = elementValence[iIndex][i];
		
		jCount = 0;
		for (j=0; j<faceVertexNumber; ++j)
			for (k=0; k<elementVertexNumber; ++k)
				if (element[iElem][k] == vertexIndex[j])
					{jCount++; break;}

		if (jCount == faceVertexNumber)
			++iCount;
	}

	if (iCount == 1)
		return true;
	else
		return false;
}

inline bool RawMesh::IsBoundaryFace(int elementID, int faceIndex)
{
	int face[20];
	GetElementFace(elementID, faceIndex, face);

	return IsBoundaryFace(face);
}

bool RawMesh::IsEdge(const int vertex_1, const int vertex_2)
{
	int i, imax;

	if(NULL == edgeValence || NULL == edgeValenceNumber)
		GetEdgeValence();

	if(vertex_1 == vertex_2)
	{
		OutputError("RawMesh::IsEdge", "The two points are the same");
		return false;		
	}

	imax = edgeValenceNumber[vertex_1];
	for(i = 0; i < imax; ++i)
		if(edgeValence[vertex_1][i] == vertex_2)
			return true;

	return false;	
}

void RawMesh::UpdateMesh()
{
	//Clear connectivity information
	FreeMatrix(edgeValence);
	FreeMatrix(edgeValenceNumber);
	FreeMatrix(elementValence);
	FreeMatrix(elementValenceNumber);

	FreeMatrix(elementNormal);
	FreeMatrix(vertexNormal);
	FreeMatrix(vertexCurvature);
}

void RawMesh::EmptyCurrentMesh()
{
	//Clear vertex information
	vertexNumber = 0;
	FreeMatrix(vertex);
	FreeMatrix(vertexSign);
	FreeMatrix(vertexColor);
	FreeMatrix(vertexNormal);
	FreeMatrix(vertexCurvature);

	//Clear element information
	elementNumber = 0;
	FreeMatrix(element);
	FreeMatrix(elementSign);
	FreeMatrix(elementColor);
	FreeMatrix(elementNormal);

	//Clear connectivity information
	FreeMatrix(edgeValence);
	FreeMatrix(edgeValenceNumber);
	FreeMatrix(elementValence);
	FreeMatrix(elementValenceNumber);

	SetElementType(NA);

	meshInfo.biggestSize = 0.0f;
	meshInfo.center[0] = 0.0f;
	meshInfo.center[1] = 0.0f;
	meshInfo.center[2] = 0.0f;
	meshInfo.maxCorner[0] = 0.0f;
	meshInfo.maxCorner[1] = 0.0f;
	meshInfo.maxCorner[2] = 0.0f;
	meshInfo.minCorner[0] = 0.0f;	
	meshInfo.minCorner[1] = 0.0f;
	meshInfo.minCorner[2] = 0.0f;
}

void RawMesh::AddVertex(int addedNumber, double (*addedVertex)[3], int *addedVertexSign)
{
	double (*newVertex)[3]=NULL;
	int *newVertexSign=NULL;
	int i, j, newVertexNum;

	if (addedNumber == 0)
		return;

	//add vertexes
	newVertexNum = vertexNumber + addedNumber;
	if (!InitiateMatrix(newVertex, newVertexNum) || !InitiateMatrix(newVertexSign, newVertexNum, 0))
	{
		OutputError("RawMesh::AddVertex", "Fail to initiate newVertex");
		return;
	}

	for (i=0; i<vertexNumber; ++i)
	{
		for (j=0; j<3; ++j)
			newVertex[i][j] = vertex[i][j];
		newVertexSign[i] = vertexSign[i];
	}
	for (i=0; i<addedNumber; ++i)
	{
		for (j=0; j<3; j++)
			newVertex[vertexNumber + i][j] = addedVertex[i][j];
	}
	if (addedVertexSign != NULL)
	{
		for (i=0; i<addedNumber; ++i)
			newVertexSign[vertexNumber + i] = addedVertexSign[i];
	}
	
	FreeMatrix(vertex);
	FreeMatrix(vertexSign);
	vertex = newVertex;
	vertexSign = newVertexSign;
	vertexNumber = newVertexNum;
	
	UpdateMesh();
	return;
}

void RawMesh::AddVertex(int addedNumber, double (*addedVertex)[3], int *addedVertexSign, double (*addedVertexColor)[3])
{
	double (*newVertex)[3]=NULL, (*newVertexColor)[3]=NULL;
	int *newVertexSign=NULL;
	int i, j, newVertexNum;

	if (addedNumber == 0)
		return;

	//add vertexes
	newVertexNum = vertexNumber + addedNumber;
	if (!InitiateMatrix(newVertex, newVertexNum) || !InitiateMatrix(newVertexSign, newVertexNum, 0) || !InitiateMatrix(newVertexColor, newVertexNum, 0.0))
	{
		OutputError("RawMesh::AddVertex", "Fail to initiate newVertex");
		return;
	}

	for (i=0; i<vertexNumber; ++i)
	{
		for (j=0; j<3; ++j)
		{
			newVertex[i][j] = vertex[i][j];
			newVertexColor[i][j] = vertexColor[i][j];
		}
		newVertexSign[i] = vertexSign[i];
	}
	for (i=0; i<addedNumber; ++i)
	{
		for (j=0; j<3; j++)
			newVertex[vertexNumber + i][j] = addedVertex[i][j];
	}
	if (addedVertexSign != NULL)
	{
		for (i=0; i<addedNumber; ++i)
			newVertexSign[vertexNumber + i] = addedVertexSign[i];
	}
	if (addedVertexSign != NULL)
	{
		for (i=0; i<addedNumber; ++i)
		{
			for (j=0; j<3; j++)
				newVertexColor[vertexNumber + i][j] = addedVertexColor[i][j];
		}
	}
	
	FreeMatrix(vertex);
	FreeMatrix(vertexSign);
	FreeMatrix(vertexColor);
	vertex = newVertex;
	vertexSign = newVertexSign;
	vertexColor = newVertexColor;
	vertexNumber = newVertexNum;
	
	UpdateMesh();
	return;
}

void RawMesh::AddVertex(vector <double> &addedVertex, vector <int> &addedVertexSign)
{
	int i, j, addedNumber, newVertexNum;

	addedNumber = addedVertex.size()/3;
	if (addedNumber == 0)
		return;

	//add vertexes
	double (*newVertex)[3]=NULL;
	int *newVertexSign=NULL;
	newVertexNum = vertexNumber + addedNumber;
	if (!InitiateMatrix(newVertex, newVertexNum) || !InitiateMatrix(newVertexSign, newVertexNum, 0))
	{
		OutputError("RawMesh::AddVertex", "Fail to initiate newVertex");
		return;
	}

	for (i=0; i<vertexNumber; ++i)
	{
		for (j=0; j<3; ++j)
			newVertex[i][j] = vertex[i][j];
		newVertexSign[i] = vertexSign[i];
	}
	for (i=0; i<addedNumber; ++i)
	{
		for (j=0; j<3; ++j)
			newVertex[vertexNumber + i][j] = addedVertex[3*i + j];
		newVertexSign[vertexNumber + i] = addedVertexSign[i];
	}
	
	FreeMatrix(vertex);
	FreeMatrix(vertexSign);
	vertex = newVertex;
	vertexSign = newVertexSign;
	vertexNumber = newVertexNum;
	UpdateMesh();	

	addedVertex.clear();
	addedVertexSign.clear();
	
	return;
}

void RawMesh::AddVertex(vector <double> &addedVertex)
{
	int i, j, addedNumber, newVertexNum;

	addedNumber = addedVertex.size()/3;
	if (addedNumber == 0)
		return;

	//add vertexes
	double (*newVertex)[3]=NULL;
	int *newVertexSign=NULL;
	newVertexNum = vertexNumber + addedNumber;
	if (!InitiateMatrix(newVertex, newVertexNum) || !InitiateMatrix(newVertexSign, newVertexNum, 0))
	{
		OutputError("RawMesh::AddVertex", "Fail to initiate newVertex");
		return;
	}

	for (i=0; i<vertexNumber; ++i)
	{
		for (j=0; j<3; ++j)
			newVertex[i][j] = vertex[i][j];
		newVertexSign[i] = vertexSign[i];
	}
	for (i=0; i<addedNumber; ++i)
	{
		for (j=0; j<3; ++j)
			newVertex[vertexNumber + i][j] = addedVertex[3*i + j];
	}
	
	FreeMatrix(vertex);
	FreeMatrix(vertexSign);
	vertex = newVertex;
	vertexSign = newVertexSign;
	vertexNumber = newVertexNum;
	UpdateMesh();	

	addedVertex.clear();	
	return;
}

void RawMesh::AddElement(int addedNumber, int **addedElement, int *addedElementSign)
{
	int **newElement=NULL, *newElementSign=NULL;
	int i, j, newElementNumber, elementVertexNumber;

	if (addedNumber == 0)
		return;

	newElementNumber = elementNumber + addedNumber;
	if (!InitiateMatrix(newElement, newElementNumber, elementProperty.vertexNumber) ||
		!InitiateMatrix(newElementSign, newElementNumber, 0))
	{
		OutputError("RawMesh::AddElement", "Fail to initiate newElement");
		return;
	}

	elementVertexNumber = elementProperty.vertexNumber;
	for (i=0; i<elementNumber; ++i)
	{
		for (j=0; j<elementVertexNumber; ++j)
			newElement[i][j] = element[i][j];
		newElementSign[i] = elementSign[i];
	}
	for (i=0; i<addedNumber; ++i)
	{
		for (j=0; j<elementVertexNumber; ++j)
			newElement[elementNumber + i][j] = addedElement[i][j];
	}
	if (addedElementSign != NULL)
	{
		for (i=0; i<addedNumber; ++i)
			newElementSign[elementNumber + i] = addedElementSign[i];
	}
	
	FreeMatrix(element);
	FreeMatrix(elementSign);
	element = newElement;
	elementSign = newElementSign;
	elementNumber = newElementNumber;
	
	return;
}

void RawMesh::AddElement(vector <int> &addedElement, vector <int> &addedElementSign)
{
	int i, j, addedNumber, newElementNumber, elementVertexNumber;

	if (addedElement.size() == 0)
		return;

	addedNumber = addedElement.size()/elementProperty.vertexNumber;

	int **newElement=NULL, *newElementSign=NULL;
	newElementNumber = elementNumber + addedNumber;
	if (!InitiateMatrix(newElement, newElementNumber, elementProperty.vertexNumber) ||
		!InitiateMatrix(newElementSign, newElementNumber, 0))
	{
		OutputError("RawMesh::AddElement", "Fail to initiate newElement");
		return;
	}

	elementVertexNumber = elementProperty.vertexNumber;
	for (i=0; i<elementNumber; ++i)
	{
		for (j=0; j<elementVertexNumber; ++j)
			newElement[i][j] = element[i][j];
		newElementSign[i] = elementSign[i];
	}
	for (i=0; i<addedNumber; ++i)
	{
		for (j=0; j<elementVertexNumber; ++j)
			newElement[elementNumber + i][j] = addedElement[i*elementVertexNumber + j];
		newElementSign[elementNumber + i] = addedElementSign[i];
	}

	FreeMatrix(element);
	FreeMatrix(elementSign);
	element = newElement;
	elementSign = newElementSign;
	elementNumber = newElementNumber;

	addedElement.clear();
	return;
}

void RawMesh::AddElement(vector <int> &addedElement)
{
	int i, j, addedNumber, newElementNumber, elementVertexNumber;

	if (addedElement.size() == 0)
		return;

	addedNumber = addedElement.size()/elementProperty.vertexNumber;

	int **newElement=NULL, *newElementSign=NULL;
	newElementNumber = elementNumber + addedNumber;
	if (!InitiateMatrix(newElement, newElementNumber, elementProperty.vertexNumber) ||
		!InitiateMatrix(newElementSign, newElementNumber, 0))
	{
		OutputError("RawMesh::AddElement", "Fail to initiate newElement");
		return;
	}

	elementVertexNumber = elementProperty.vertexNumber;
	for (i=0; i<elementNumber; ++i)
	{
		for (j=0; j<elementVertexNumber; ++j)
			newElement[i][j] = element[i][j];
		newElementSign[i] = elementSign[i];
	}
	for (i=0; i<addedNumber; ++i)
	{
		for (j=0; j<elementVertexNumber; ++j)
			newElement[elementNumber + i][j] = addedElement[i*elementVertexNumber + j];
	}

	FreeMatrix(element);
	FreeMatrix(elementSign);
	element = newElement;
	elementSign = newElementSign;
	elementNumber = newElementNumber;

	addedElement.clear();
	return;
}

void RawMesh::ApplyColorMapToVertex(int mapType)
{
	//double minColor[3]={1.0, 0.7, 0.0}, maxColor[3]={0.6, 0.8, 1.0};
	//double minColor[3]={1.0, 1.0, 1.0}, maxColor[3]={0.0, 0.0, 1.0};
	//double minColor[3]={1.0, 0.0, 0.0}, midColor[3] = {0.0, 1.0, 0.0}, maxColor[3]={0.0, 0.0, 1.0};
	double colorMap[3][3][3]={{{0.0, 0.0, 1.0}, {0.0, 1.0, 0.0}, {1.0, 0.0, 0.0}},
							 {{1.0, 0.7, 0.0}, {0.8, 0.75, 0.5}, {0.6, 0.8, 1.0}},
							 {{1.0, 1.0, 1.0}, {0.5, 0.5, 1.0}, {0.0, 0.0, 1.0}}};

	int i;
	double minColor[3], midColor[3], maxColor[3];

	if (mapType >= 3)
		mapType = 0;
	for (i=0; i<3; ++i)
	{
		minColor[i] = colorMap[mapType][0][i];
		midColor[i] = colorMap[mapType][1][i];
		maxColor[i] = colorMap[mapType][2][i];
	}

	int minSign=0, maxSign=0, midSign;

	for (i=0; i<vertexNumber; ++i)
	{
		if (vertexSign[i] < minSign)
			minSign = vertexSign[i];
		else if (vertexSign[i] > maxSign)
			maxSign = vertexSign[i];
	}
	midSign = (maxSign + minSign)/2;

	InitiateMatrix(vertexColor, vertexNumber);
	for (i=0; i<vertexNumber; ++i)
	{
		if (vertexSign[i] <= midSign)
		{
			vertexColor[i][0] = minColor[0] + (midColor[0] - minColor[0])*(vertexSign[i]-minSign)/(midSign - minSign);
			vertexColor[i][1] = minColor[1] + (midColor[1] - minColor[1])*(vertexSign[i]-minSign)/(midSign - minSign);
			vertexColor[i][2] = minColor[2] + (midColor[2] - minColor[2])*(vertexSign[i]-minSign)/(midSign - minSign);
		}
		else
		{
			vertexColor[i][0] = midColor[0] + (maxColor[0] - midColor[0])*(vertexSign[i]-midSign)/(maxSign - midSign);
			vertexColor[i][1] = midColor[1] + (maxColor[1] - midColor[1])*(vertexSign[i]-midSign)/(maxSign - midSign);
			vertexColor[i][2] = midColor[2] + (maxColor[2] - midColor[2])*(vertexSign[i]-midSign)/(maxSign - midSign);
		}
	}

	printf("ApplyColorMapToVertex::Colormap type: %d; Max: %d; Min: %d\n", mapType, maxSign, minSign);

	return;
}

void RawMesh::ClearOldMesh()
{
	meshName.clear();
	meshName += "Default.rawm";

	meshInfo.isUseVertexNormal = false;
	meshInfo.isUseVertexColor = false;
	meshInfo.isUseElementNormal = false;
	meshInfo.isUseElementColor = false;

	EmptyCurrentMesh();
}

bool RawMesh::CreateNewMesh(ElementType const myType, int const myVertexNumber, int const myElementNumber)
{
	SetElementType(myType);

	vertexNumber = myVertexNumber;
	if (! InitiateMatrix(vertex, vertexNumber))
		return false;
	if (! InitiateMatrix(vertexSign, vertexNumber, 0))
		return false;

	elementNumber = myElementNumber;
	if (! InitiateMatrix(element, elementNumber, elementProperty.vertexNumber))
		return false;
	if (! InitiateMatrix(elementSign, elementNumber, 0))
		return false;

	return true;
}

void RawMesh::Transform(int elementID, double zoom[3], double rotate[3], double move[3])
{
	int i;
	double dTmp;
	//zoom
	for (i=0; i<elementProperty.vertexNumber; ++i)
	{
		vertex[element[elementID][i]][0] *= zoom[0];
		vertex[element[elementID][i]][1] *= zoom[1];
		vertex[element[elementID][i]][2] *= zoom[2];
	}

	//Fotate
	for (i=0; i<elementProperty.vertexNumber; ++i)
	{
		//x-axis
		dTmp = vertex[element[elementID][i]][1]*cos(rotate[0]) - vertex[element[elementID][i]][2]*sin(rotate[0]);
		vertex[element[elementID][i]][2] = vertex[element[elementID][i]][1]*sin(rotate[0]) + vertex[element[elementID][i]][2]*cos(rotate[0]);
		vertex[element[elementID][i]][1] = dTmp;

		//y-axis
		dTmp = vertex[element[elementID][i]][0]*cos(rotate[1]) + vertex[element[elementID][i]][2]*sin(rotate[1]);
		vertex[element[elementID][i]][2] = -vertex[element[elementID][i]][0]*sin(rotate[1]) + vertex[element[elementID][i]][2]*cos(rotate[1]);
		vertex[element[elementID][i]][0] = dTmp;

		//z-axis
		dTmp = vertex[element[elementID][i]][0]*cos(rotate[2]) - vertex[element[elementID][i]][1]*sin(rotate[2]);
		vertex[element[elementID][i]][1] = vertex[element[elementID][i]][0]*sin(rotate[2]) + vertex[element[elementID][i]][1]*cos(rotate[2]);
		vertex[element[elementID][i]][0] = dTmp;
	}

	//translate
	for (i=0; i<elementProperty.vertexNumber; ++i)
	{
		vertex[element[elementID][i]][0] += move[0];
		vertex[element[elementID][i]][1] += move[1];
		vertex[element[elementID][i]][2] += move[2];
	}

	return;
}

void RawMesh::Rotate(double rad[3], bool rotateAtCenter)
{
	GetMeshInfo();

	if (rotateAtCenter)
		Move(-meshInfo.center[0], -meshInfo.center[1], -meshInfo.center[2]);
	//zoom
	double dTmp;
	for (int i=0; i<vertexNumber; ++i)
	{
		//x-axis
		dTmp = vertex[i][1]*cos(rad[0]) - vertex[i][2]*sin(rad[0]);
		vertex[i][2] = vertex[i][1]*sin(rad[0]) + vertex[i][2]*cos(rad[0]);
		vertex[i][1] = dTmp;

		//y-axis
		dTmp = vertex[i][0]*cos(rad[1]) + vertex[i][2]*sin(rad[1]);
		vertex[i][2] = -vertex[i][0]*sin(rad[1]) + vertex[i][2]*cos(rad[1]);
		vertex[i][0] = dTmp;

		//z-axis
		dTmp = vertex[i][0]*cos(rad[2]) - vertex[i][1]*sin(rad[2]);
		vertex[i][1] = vertex[i][0]*sin(rad[2]) + vertex[i][1]*cos(rad[2]);
		vertex[i][0] = dTmp;
	}

	if (rotateAtCenter)
		Move(meshInfo.center);
	return;
}

inline void RawMesh::Rotate(double xRad, double yRad, double zRad, bool rotateAtCenter)
{
	double rad[3];

	rad[0]= xRad; rad[1] = yRad; rad[2] = zRad;
	Rotate(rad, rotateAtCenter);

	return;
}

void RawMesh::Zoom(double ratio[3], bool zoomAtCenter)
{
	GetMeshInfo();

	if (zoomAtCenter)
		Move(-meshInfo.center[0], -meshInfo.center[1], -meshInfo.center[2]);
	//zoom
	for (int i=0; i<vertexNumber; ++i)
	{
		vertex[i][0] *= ratio[0];
		vertex[i][1] *= ratio[1];
		vertex[i][2] *= ratio[2];
	}

	if (zoomAtCenter)
		Move(meshInfo.center);
	return;
}

void RawMesh::Zoom(double xRatio, double yRatio, double zRatio, bool zoomAtCenter)
{
	GetMeshInfo();

	if (zoomAtCenter)
		Move(-meshInfo.center[0], -meshInfo.center[1], -meshInfo.center[2]);
	//zoom
	for (int i=0; i<vertexNumber; ++i)
	{
		vertex[i][0] *= xRatio;
		vertex[i][1] *= yRatio;
		vertex[i][2] *= zRatio;
	}

	if (zoomAtCenter)
		Move(meshInfo.center);
	return;
}

void RawMesh::SetElementType(const ElementType myElementType)
{
	if (myElementType == NA)
	{
		elementProperty.elementType = NA;
		elementProperty.vertexNumber = 0;
		elementProperty.edgeNumber = 0;
		elementProperty.faceNumber = 0;
		elementProperty.faceVertexNumber = 0;
		elementProperty.faceType = NA;
	}
	else if (myElementType == POINT)
	{
		elementProperty.elementType = POINT;
		elementProperty.vertexNumber = 1;
		elementProperty.edgeNumber = 0;
		elementProperty.faceNumber = 0;
		elementProperty.faceVertexNumber = 0;
		elementProperty.faceType = NA;
	}
	else if (myElementType == LINE)
	{
		elementProperty.elementType = LINE;
		elementProperty.vertexNumber = 2;
		elementProperty.edgeNumber = 1;
		elementProperty.faceNumber = 0;
		elementProperty.faceVertexNumber = 0;
		elementProperty.faceType = NA;
	}
	else if (myElementType == TRIANGLE)
	{
		elementProperty.elementType = TRIANGLE;
		elementProperty.vertexNumber = 3;
		elementProperty.edgeNumber = 3;
		elementProperty.faceNumber = 1;
		elementProperty.faceVertexNumber = 3;
		elementProperty.faceType = TRIANGLE;
	}
	else if (myElementType == TETRAHEDRON)
	{
		elementProperty.elementType = TETRAHEDRON;
		elementProperty.vertexNumber = 4;
		elementProperty.edgeNumber = 6;
		elementProperty.faceNumber = 4;
		elementProperty.faceVertexNumber = 3;
		elementProperty.faceType = TRIANGLE;
	}
	else if (myElementType == QUADRILATERAL)
	{
		elementProperty.elementType = QUADRILATERAL;
		elementProperty.vertexNumber = 4;
		elementProperty.edgeNumber = 4;
		elementProperty.faceNumber = 1;
		elementProperty.faceVertexNumber = 4;
		elementProperty.faceType = QUADRILATERAL;
	}
	else if (myElementType == HEXAHEDRON)
	{
		elementProperty.elementType = HEXAHEDRON;
		elementProperty.vertexNumber = 8;
		elementProperty.edgeNumber = 12;
		elementProperty.faceNumber = 6;
		elementProperty.faceVertexNumber = 4;
		elementProperty.faceType = QUADRILATERAL;
	}
	else if (myElementType == HEXAGON)
	{
		elementProperty.elementType = HEXAGON;
		elementProperty.vertexNumber = 6;
		elementProperty.edgeNumber = 6;
		elementProperty.faceNumber = 1;
		elementProperty.faceVertexNumber = 6;
		elementProperty.faceType = HEXAGON;
	}
	else if (myElementType == POLYGON)	//Not Finish!!!
	{
		elementProperty.elementType = POLYGON;
		elementProperty.vertexNumber = 0;
		elementProperty.edgeNumber = 0;
		elementProperty.faceNumber = 0;
		elementProperty.faceVertexNumber = 0;
		elementProperty.faceType = POLYGON;
	}
	else
	{
		elementProperty.elementType = NA;
		elementProperty.vertexNumber = 0;
		elementProperty.edgeNumber = 0;
		elementProperty.faceNumber = 0;
		elementProperty.faceType = NA;
	}

	return;
}

RawMesh::ElementProperty RawMesh::GetElementProperty(const ElementType myElementType)
{
	ElementProperty *myEementProperty = new ElementProperty;
	if (myElementType == NA)
	{
		myEementProperty->elementType = NA;
		myEementProperty->vertexNumber = 0;
		myEementProperty->edgeNumber = 0;
		myEementProperty->faceNumber = 0;
		myEementProperty->faceType = NA;
	}
	if (myElementType == POINT)
	{
		myEementProperty->elementType = POINT;
		myEementProperty->vertexNumber = 1;
		myEementProperty->edgeNumber = 0;
		myEementProperty->faceNumber = 0;
		myEementProperty->faceType = NA;
	}
	if (myElementType == LINE)
	{
		myEementProperty->elementType = LINE;
		myEementProperty->vertexNumber = 2;
		myEementProperty->edgeNumber = 1;
		myEementProperty->faceNumber = 0;
		myEementProperty->faceType = NA;
	}
	else if (myElementType == TRIANGLE)
	{
		myEementProperty->elementType = TRIANGLE;
		myEementProperty->vertexNumber = 3;
		myEementProperty->edgeNumber = 3;
		myEementProperty->faceNumber = 1;
		myEementProperty->faceType = TRIANGLE;
	}
	else if (myElementType == TETRAHEDRON)
	{
		myEementProperty->elementType = TETRAHEDRON;
		myEementProperty->vertexNumber = 4;
		myEementProperty->edgeNumber = 6;
		myEementProperty->faceNumber = 4;
		myEementProperty->faceType = TRIANGLE;
	}
	else if (myElementType == QUADRILATERAL)
	{
		myEementProperty->elementType = QUADRILATERAL;
		myEementProperty->vertexNumber = 4;
		myEementProperty->edgeNumber = 4;
		myEementProperty->faceNumber = 1;
		myEementProperty->faceType = QUADRILATERAL;
	}
	else if (myElementType == HEXAHEDRON)
	{
		myEementProperty->elementType = HEXAHEDRON;
		myEementProperty->vertexNumber = 8;
		myEementProperty->edgeNumber = 12;
		myEementProperty->faceNumber = 6;
		myEementProperty->faceType = QUADRILATERAL;
	}
	else if (myElementType == HEXAGON)
	{
		myEementProperty->elementType = HEXAGON;
		myEementProperty->vertexNumber = 6;
		myEementProperty->edgeNumber = 6;
		myEementProperty->faceNumber = 1;
		myEementProperty->faceType = HEXAGON;
	}
	else if (myElementType == POLYGON)	//Not Finish!!!
	{
		myEementProperty->elementType = POLYGON;
		myEementProperty->vertexNumber = 0;
		myEementProperty->edgeNumber = 0;
		myEementProperty->faceNumber = 0;
		myEementProperty->faceType = POLYGON;
	}
	else
	{
		myEementProperty->elementType = NA;
		myEementProperty->vertexNumber = 0;
		myEementProperty->edgeNumber = 0;
		myEementProperty->faceNumber = 0;
		myEementProperty->faceType = NA;
	}

	return *myEementProperty;
}

inline bool RawMesh::InitiateMatrix(double ** &matrix, int const dim1, int const dim2, double initialValue)
{
	if (matrix != NULL)
		FreeMatrix(matrix);

	double *tmp = NULL;
	int i, iSum;
	iSum = dim1*dim2;
	tmp = new double [iSum];
	matrix = new double *[dim1];
	if (tmp == NULL || matrix == NULL)
		return false;

	for (i=0; i<iSum; ++i)
		tmp[i] = initialValue;

	for (i=0; i<dim1; ++i)
		matrix[i] = tmp + dim2*i;

	return true;
}

inline bool RawMesh::InitiateMatrix(int ** &matrix, int const dim1, int const dim2, int initialValue)
{
	if (matrix != NULL)
		FreeMatrix(matrix);

	int *tmp = NULL;
	int i, iSum;
	iSum = dim1*dim2;
	tmp = new int [iSum];
	matrix = new int *[dim1];
	if (tmp == NULL || matrix == NULL)
		return false;

	for (i=0; i<iSum; ++i)
		tmp[i] = initialValue;

	for (i=0; i<dim1; ++i)
		matrix[i] = tmp + dim2*i;

	return true;
}


bool RawMesh::InitiateElementValenceMember()
{
	int i, i_tmp;

	if(vertexNumber == 0)
		return false;

	FreeElementValence();

	if((elementValence = new int *[vertexNumber]) == NULL)
		return false;
	int *tmp;
	i_tmp = elementNumber * elementProperty.vertexNumber;
	if((tmp = new int [i_tmp]) == NULL)
		return false;
	for(i = 0; i < i_tmp; ++i)
		tmp[i] = -1;
	elementValence[0] = tmp;

	if((elementValenceNumber = new int [vertexNumber]) == NULL)
		return false;
	for(i = 0; i < vertexNumber; ++i)
		elementValenceNumber[i] = 0;

	return true;
}

bool RawMesh::InitiateEdgeValenceMember(const int &sumNum)
{
	int i;

	if(vertexNumber == 0)
		return false;

	FreeEdgeValence();

	if((edgeValence = new int *[vertexNumber]) == NULL)
		return false;
	int *tmp;
	if((tmp = new int [sumNum]) == NULL)
		return false;
	edgeValence[0] = tmp;
	for(i = 0; i < sumNum; ++i)
			tmp[i] = -1;

	if((edgeValenceNumber = new int [vertexNumber]) == NULL)
		return false;
	for(i = 0; i < vertexNumber; ++i)
		edgeValenceNumber[i] = 0;

	return true;
}


inline bool RawMesh::InitiateMatrix(bool ** &matrix, int const dim1, int const dim2, bool initialValue)
{
	if (matrix != NULL)
		FreeMatrix(matrix);

	bool *tmp = NULL;
	int i, iSum;
	iSum = dim1*dim2;
	tmp = new bool [iSum];
	matrix = new bool *[dim1];
	if (tmp == NULL || matrix == NULL)
		return false;

	for (i=0; i<iSum; ++i)
		tmp[i] = initialValue;

	for (i=0; i<dim1; ++i)
		matrix[i] = tmp + dim2*i;

	return true;
}

inline bool RawMesh::InitiateMatrix(double * &matrix, int const dim, double initialValue)
{
	if (matrix != NULL)
		FreeMatrix(matrix);

	matrix = new double [dim];
	if (matrix == NULL)
		return false;

	for (int i=0; i<dim; ++i)
		matrix[i] = initialValue;

	return true;
}

inline bool RawMesh::InitiateMatrix(bool * &matrix, int const dim, bool initialValue)
{
	if (matrix != NULL)
		FreeMatrix(matrix);

	matrix = new bool [dim];
	if (matrix == NULL)
		return false;

	for (int i=0; i<dim; ++i)
		matrix[i] = initialValue;

	return true;
}

inline bool RawMesh::InitiateMatrix(int * &matrix, int const dim, int initialValue)
{
	if (matrix != NULL)
		FreeMatrix(matrix);

	matrix = new int [dim];
	if (matrix == NULL)
		return false;

	for (int i=0; i<dim; ++i)
		matrix[i] = initialValue;

	return true;
}

inline bool RawMesh::InitiateMatrix(double (* &matrix)[3], int const dim, double initialValue)
{
	if (matrix != NULL)
		FreeMatrix(matrix);

	matrix = new double [dim][3];
	if (matrix == NULL)
		return false;

	for (int i=0; i<dim; ++i)
	{
		matrix[i][0] = initialValue;
		matrix[i][1] = initialValue;
		matrix[i][2] = initialValue;
	}

	return true;
}

inline void RawMesh::FreeMatrix(double ** &matrix)
{
	if (matrix != NULL)
	{
		delete [] matrix[0];
		delete [] matrix;
		matrix = NULL;
	}

	return;	
}

inline void RawMesh::FreeMatrix(int ** &matrix)
{
	if (matrix != NULL)
	{
		delete [] matrix[0];
		delete [] matrix;
		matrix = NULL;
	}

	return;	
}

inline void RawMesh::FreeMatrix(bool ** &matrix)
{
	if (matrix != NULL)
	{
		delete [] matrix[0];
		delete [] matrix;
		matrix = NULL;
	}

	return;	
}

inline void RawMesh::FreeMatrix(bool * &matrix)
{
	if (matrix != NULL)
	{
		delete [] matrix;
		matrix = NULL;
	}

	return;	
}

inline void RawMesh::FreeMatrix(double * &matrix)
{
	if (matrix != NULL)
	{
		delete [] matrix;
		matrix = NULL;
	}

	return;	
}

inline void RawMesh::FreeMatrix(int * &matrix)
{
	if (matrix != NULL)
	{
		delete [] matrix;
		matrix = NULL;
	}

	return;	
}


inline void RawMesh::FreeMatrix(double (* &matrix)[3])
{
	if (matrix != NULL)
	{
		delete [] matrix;
		matrix = NULL;
	}

	return;	
}

void RawMesh::FreeElementValence()
{
	if(elementValence != NULL)
	{
		delete [] elementValence[0];
		delete [] elementValence;
		elementValence = NULL;
	}

	if(elementValenceNumber != NULL)
	{
		delete [] elementValenceNumber;
		elementValenceNumber = NULL;
	}

	return;
}

void RawMesh::FreeEdgeValence()
{
	if(edgeValence != NULL)
	{
		delete [] edgeValence[0];
		delete [] edgeValence;
		edgeValence = NULL;
	}

	if(edgeValenceNumber != NULL)
	{
		delete [] edgeValenceNumber;
		edgeValenceNumber = NULL;
	}

	return;
}


bool RawMesh::Read(const char * filename)
{
	string::size_type fileLength;
	bool readSuccess;
	string filenameInUpperCase;
	char * cTmp;

	ClearOldMesh();
	meshName.clear();
	meshName += filename;
	fileLength = meshName.size();

	cTmp = _strupr(_strdup(filename));
	filenameInUpperCase += cTmp;

	if		(filenameInUpperCase.find(".RAWM") == (fileLength - 5))
		readSuccess = Read_RawMesh(filename);
	else if	(filenameInUpperCase.find(".MESH") == (fileLength - 5))
		readSuccess = Read_Mesh(filename);
	else if (filenameInUpperCase.find(".RAW") < fileLength)
		readSuccess = Read_RawFile(filename);
	else if (filenameInUpperCase.find(".PLY") == (fileLength - 4))
		readSuccess = Read_Ply(filename);
	else if (filenameInUpperCase.find(".PLT") == (fileLength - 4))
		readSuccess = Read_Plt(filename);
	else if (filenameInUpperCase.find(".INP") == (fileLength - 4))
		readSuccess = Read_Inp(filename);
	else if (filenameInUpperCase.find(".VTK") == (fileLength - 4))
		readSuccess = Read_Vtk(filename);
	else
		readSuccess = false;

	//Clear
	filenameInUpperCase.clear();
	delete [] cTmp;

	return readSuccess;
}

bool RawMesh::Read_RawMesh(const char * filename)
{
	FILE *input;
	int i, j, elementType;
	int elemInd, vSign, eSign;
	double x, y, z;

	input = NULL;
	input = fopen(filename, "r");
	if(input == NULL)
	{
		//printf("Open file Error!\n");
		return false;
	}

	//Check if Rawm is the old version.
	char title[20], tmpVersion[20];
	fscanf(input, "%s %s\n", title, tmpVersion);
	if (strcmp(_strupr(title), "VERSION") != 0)
	{
		fclose(input);
		version = "1.0";
		return Read_RawMesh_OldVersion(filename);		
	}
	version = tmpVersion;
	//End check

	fscanf(input,"%d %d %d\n", &vertexNumber, &elementNumber, &elementType);
	fscanf(input, "%d %d %d %d\n", &meshInfo.isUseVertexNormal, &meshInfo.isUseVertexColor, &meshInfo.isUseElementNormal, &meshInfo.isUseElementColor);

	if (! CreateNewMesh((ElementType) elementType, vertexNumber, elementNumber))
		return false;	

	//Read vertex information
	for (i=0; i<vertexNumber; ++i)
	{
		fscanf(input, "%lf %lf %lf %d", &x, &y, &z, &vSign);
		vertex[i][0] = x;
		vertex[i][1] = y;
		vertex[i][2] = z;
		vertexSign[i] = vSign;
	}

	if (meshInfo.isUseVertexNormal)
	{
		InitiateMatrix(vertexNormal, vertexNumber);

		double nx, ny, nz;
		for (i=0; i<vertexNumber; ++i)
		{
			fscanf(input, "%lf %lf %lf", &nx, &ny, &nz);
			vertexNormal[i][0] = nx;
			vertexNormal[i][1] = ny;
			vertexNormal[i][2] = nz;
		}
	}

	if (meshInfo.isUseVertexColor)
	{
		InitiateMatrix(vertexColor, vertexNumber);

		double r, g, b;
		for (i=0; i<vertexNumber; ++i)
		{
			fscanf(input, "%lf %lf %lf", &r, &g, &b);
			vertexColor[i][0] = r;
			vertexColor[i][1] = g;
			vertexColor[i][2] = b;
		}
	}
	//End of reading vertex information

	//Read element information
	for(i=0; i<elementNumber; ++i)
	{
		for (j=0; j<elementProperty.vertexNumber; ++j)
		{
			fscanf(input, "%d ", &elemInd);
			element[i][j] = elemInd;
		}

		fscanf(input, "%d\n", &eSign);
		elementSign[i] = eSign;
	}

	if (meshInfo.isUseElementNormal)
	{
		InitiateMatrix(elementNormal, elementNumber);

		double nx, ny, nz;
		for (i=0; i<elementNumber; ++i)
		{
			fscanf(input, "%lf %lf %lf", &nx, &ny, &nz);
			elementNormal[i][0] = nx;
			elementNormal[i][1] = ny;
			elementNormal[i][2] = nz;
		}
	}

	if (meshInfo.isUseElementColor)
	{
		InitiateMatrix(elementColor, elementNumber);

		double r, g, b;
		for (i=0; i<elementNumber; ++i)
		{
			fscanf(input, "%lf %lf %lf", &r, &g, &b);
			elementColor[i][0] = r;
			elementColor[i][1] = g;
			elementColor[i][2] = b;
		}
	}
	//End of reading element information

	fclose(input);

	return true;
}

bool RawMesh::Read_RawMesh_OldVersion(const char * filename)
{
	FILE *input;
	int i, j, elementType;
	int elemInd, vSign, eSign;
	double x, y, z;

	input = NULL;
	input = fopen(filename, "r");
	if(input == NULL)
	{
		//printf("Open file Error!\n");
		return false;
	}

	fscanf(input,"%d %d %d\n", &vertexNumber, &elementNumber, &elementType);

	if (! CreateNewMesh((ElementType) elementType, vertexNumber, elementNumber))
		return false;	

	for (i=0; i<vertexNumber; ++i)
	{
		fscanf(input, "%lf %lf %lf %d", &x, &y, &z, &vSign);
		vertex[i][0] = x;
		vertex[i][1] = y;
		vertex[i][2] = z;
		vertexSign[i] = vSign;
	}

	for (i=0; i<elementNumber; ++i)
	{
		for (j=0; j<elementProperty.vertexNumber; ++j)
		{
			fscanf(input, "%d ", &elemInd);
			element[i][j] = elemInd;
		}

		fscanf(input, "%d\n", &eSign);
		elementSign[i] = eSign;
	}

	fclose(input);

	return true;
}

bool RawMesh::Read_RawFile(const char * filename)
{
	int i, fdot, flen;
	string::size_type fileLength;
	bool readSuccess;
	char * fileType;
	string filenameInUpperCase;
	char * cTmp;

	fileLength = meshName.size();

	//Find the suffix of the file
	flen = (int) fileLength;
	for(i = flen-1; i >= 0; i--)
		if(filename[i] == '.') break;
	if(i < 0)
		return false;

	fdot = i+1;
	fileType = new char[flen - fdot + 1];
	for(i = fdot; i <= flen; ++i)
		fileType[i-fdot] = filename[i];

	cTmp = _strupr(_strdup(filename));
	filenameInUpperCase += cTmp;

	//Read RawFile
	if		(filenameInUpperCase.find("_TRI") < fileLength)
		readSuccess = ReadRaw_tri(filename, fileType);
	else if	(filenameInUpperCase.find("_TET") < fileLength)
		readSuccess = ReadRaw_tet(filename, fileType);
	else if (filenameInUpperCase.find("_QUAD") < fileLength)
		readSuccess = ReadRaw_quad(filename, fileType);
	else if (filenameInUpperCase.find("_HEX") < fileLength)
		readSuccess = ReadRaw_hex(filename, fileType);
	else if (filenameInUpperCase.find("_LINE") < fileLength)
		readSuccess = ReadRaw_line(filename, fileType);
	else
		readSuccess = false;

	if (readSuccess)
	{
		for (i=0; i<elementNumber; ++i)
			elementSign[i] = 0;
	}

	//Clear
	filenameInUpperCase.clear();
	delete [] cTmp;
	delete [] fileType;

	return readSuccess;
}

bool RawMesh::Read_Ply(const char * filename)
{
	FILE *input;
	int i, j, elementType;
	int elemInd, vSign, eSign;
	double x, y, z, tmp;

	input = NULL;
	input = fopen(filename, "r");
	if(input == NULL)
	{
		//printf("Open file Error!\n");
		return false;
	}

	char str[6][100], vertCount=0, elemCount=0;
	fscanf(input, "%s\n", str[0]);
	fscanf(input, "%s %s %s\n", str[0], str[1], str[2]);
	fscanf(input, "%s ", str[0]);
	fgets(str[1], 100, input);
	if (strcmp(str[0], "element") == 0)
		sscanf(str[1], "%s %d", str[2], &vertexNumber);
	else
		fscanf(input, "%s %s %d\n", str[0], str[1], &vertexNumber);
	fscanf(input, "%s %s %s\n", str[0], str[1], str[2]);
	while (strcmp(str[0], "property") == 0)
	{
		++vertCount;
		fscanf(input, "%s %s %s\n", str[0], str[1], str[2]);
	}
	elementNumber = atoi(str[2]);
	fscanf(input, "%s", str[0]);
	fgets(str[1], 100, input);
	while (strcmp(str[0], "property") == 0)
	{
		++elemCount;
		fscanf(input, "%s", str[0]);
		fgets(str[1], 100, input);
	}

	if (! CreateNewMesh(TRIANGLE, vertexNumber, elementNumber))
		return false;	

	if (vertCount == 3)
	{
		for (i=0; i<vertexNumber; ++i)
		{
		
			fscanf(input, "%lf %lf %lf\n", &x, &y, &z);
			vertex[i][0] = x;
			vertex[i][1] = y;
			vertex[i][2] = z;
			vertexSign[i] = 1;
		}
	}
	else if (vertCount == 4)
	{
		for (i=0; i<vertexNumber; ++i)
		{
		
			fscanf(input, "%lf %lf %lf %lf\n", &x, &y, &z, &tmp);
			vertex[i][0] = x;
			vertex[i][1] = y;
			vertex[i][2] = z;
			vertexSign[i] = (int) (tmp);
		}
	}
	else if (vertCount == 5)
	{
		for (i=0; i<vertexNumber; ++i)
		{
		
			fscanf(input, "%lf %lf %lf %lf %lf\n", &x, &y, &z, &tmp, &tmp);
			vertex[i][0] = x;
			vertex[i][1] = y;
			vertex[i][2] = z;
			vertexSign[i] = 1;
		}
	}

	if (elemCount == 1)
	{
		for (i=0; i<elementNumber; ++i)
		{
			fscanf(input, "%d ", &elemInd);	//Pass
			for (j=0; j<elementProperty.vertexNumber; ++j)
			{
				fscanf(input, "%d ", &elemInd);
				element[i][j] = elemInd;
			}

			elementSign[i] = 0;
		}
	}
	else if (elemCount == 2)
	{
		for (i=0; i<elementNumber; ++i)
		{
			fscanf(input, "%d ", &elemInd);	//Pass
			for (j=0; j<elementProperty.vertexNumber; ++j)
			{
				fscanf(input, "%d ", &elemInd);
				element[i][j] = elemInd;
			}

			fscanf(input, "%d\n", &eSign);
			elementSign[i] = eSign;
		}
	}

	fclose(input);

	return true;
}

bool RawMesh::Read_Vtk(const char * filename)
{
	FILE *input;
	int i, j, elementType;
	int v[8];

	ifstream infile;
	string oneLine,oneWrod, formerLine;

	infile.open(filename);
	if (!infile)
	{
		cerr << "error: unable to open the file "<<ends;
		cerr <<filename <<endl;
		return false;
	}
	getline(infile,oneLine);
	getline(infile,oneLine);
	getline(infile,oneLine);
	getline(infile,oneLine);
	getline(infile,oneLine);

	istringstream strStream(oneLine);
	strStream >> oneWrod >> vertexNumber >> oneWrod;

	for (i = 0; i < vertexNumber; ++i)
	{
		getline(infile,oneLine);
		istringstream st(oneLine);
		st >> vertex[i][0] >> vertex[i][1] >> vertex[i][2];
	}

	getline(infile,oneLine);
	strStream.clear();
	strStream.str(oneLine);
	strStream >> oneWrod >> elementNumber;
	getline(infile,oneLine);
	strStream.clear();
	strStream.str(oneLine);
	strStream >> elementType;

	switch(elementType)
	{
	case 2:
		{
			if (! CreateNewMesh(LINE, vertexNumber, elementNumber))
				return false;
			break;
		}
	case 3:
		{
			if (! CreateNewMesh(TRIANGLE, vertexNumber, elementNumber))
				return false;
			break;
		}
	case 4:
		{
			if (! CreateNewMesh(QUADRILATERAL, vertexNumber, elementNumber))
				return false;
			break;
		}
	case 8:
		{
			if (! CreateNewMesh(HEXAHEDRON, vertexNumber, elementNumber))
				return false;
			break;
		}
	}

	switch(elementType)
	{
	case 2:
		{
			for (i = 0; i < elementNumber; ++i)
			{
				getline(infile,oneLine);
				istringstream st(oneLine);
				for (j = 0; j < 2; ++j)
				{
					st>>element[i][j];
				}
			}
			break;
		}
	case 3:
		{
			for (i = 0; i < elementNumber; ++i)
			{
				getline(infile,oneLine);
				istringstream st(oneLine);
				for (j = 0; j < 3; ++j)
				{
					st>>element[i][j];
				}
			}
			break;
		}
	case 4:
		{
			for (i = 0; i < elementNumber; ++i)
			{
				getline(infile,oneLine);
				istringstream st(oneLine);
				for (j = 0; j < 4; ++j)
				{
					st>>element[i][j];
				}
			}
			break;
		}
	case 8:
		{
			for (i = 0; i < elementNumber; ++i)
			{
				getline(infile,oneLine);
				istringstream st(oneLine);
				for (j = 0; j < 8; ++j)
				{
					st>>element[i][j];
				}
			}
			break;
		}
	}

	infile.close(); 


	return true;
}

bool RawMesh::Read_Inp(const char * filename)
{
	FILE *input;
	int i, j, elementType;
	int elemInd, vSign, eSign;
	double x, y, z, tmp;
	 
	ifstream infile;
	string oneLine,oneWrod, formerLine;
	istringstream strStream(oneLine);
	infile.open(filename);
	if (!infile)
	{
		cerr << "error: unable to open the file "<<ends;
		cerr <<filename <<endl;
		return false;
	}

	while(oneWrod.compare("*Node"))
	{
		getline(infile,oneLine);
		istringstream st(oneLine);
		st >> oneWrod;
	}
	
	while(oneWrod.compare("*Element,"))
	{
		formerLine = oneLine;
		getline(infile,oneLine);
		istringstream st(oneLine);
		st >> oneWrod;
	}
	strStream.clear();
	strStream.str(oneLine);	
	strStream >> oneWrod ;
	strStream >> oneWrod ;
	if (0 == oneWrod.compare(5,4,"CPS3"))
	{
		elementType = 3;
	}
	else if (0 == oneWrod.compare(5,3,"S4R"))
	{
		elementType = 5;
	}
	else if (0 == oneWrod.compare(5,2,"S3"))
	{
		elementType = 3;
	}
	else if (0 == oneWrod.compare(5,4,"C3D8") || 0 == oneWrod.compare(5,5,"C3D20") || 0 == oneWrod.compare(5,5,"C3D27")
		|| 0 == oneWrod.compare(6,4,"C3D8")|| 0 == oneWrod.compare(7,4,"C3D8"))
	{
		elementType = 6;
	}
	else if (0 == oneWrod.compare(5,4,"F2D2"))
	{
		elementType = 2;
	}
	else if (0 == oneWrod.compare(5,5,"C3D10") || 0 == oneWrod.compare(5,4,"C3D4")|| 0 == oneWrod.compare(6,4,"C3D4"))
	{
		elementType = 4;
	}
	strStream.clear();
	strStream.str(formerLine);
	strStream >> oneWrod;
	oneWrod.erase(oneWrod.size()-1,1);
	strStream.clear();
	strStream.str(oneWrod);
	strStream >> vertexNumber;
	

	while(oneWrod.compare("*End"))
	{
		formerLine = oneLine;
		getline(infile,oneLine);
		istringstream st(oneLine);
		st >> oneWrod;
	}
	strStream.clear();
	strStream.str(formerLine);
	strStream >> oneWrod;
	oneWrod.erase(oneWrod.size()-1,1);
	strStream.clear();
	strStream.str(oneWrod);
	strStream >> elementNumber;

	infile.close();

	switch(elementType)
	{
	case 2:
		{
			if (! CreateNewMesh(LINE, vertexNumber, elementNumber))
				return false;
			break;
		}
	case 3:
		{
			if (! CreateNewMesh(TRIANGLE, vertexNumber, elementNumber))
				return false;
			break;
		}
	case 4:
		{
			if (! CreateNewMesh(TETRAHEDRON, vertexNumber, elementNumber))
				return false;
			break;
		}
	case 5:
		{
			if (! CreateNewMesh(QUADRILATERAL, vertexNumber, elementNumber))
				return false;
			break;
		}
	case 8:
		{
			if (! CreateNewMesh(HEXAHEDRON, vertexNumber, elementNumber))
				return false;
			break;
		}
	}

	infile.open(filename);
	
	while(oneWrod.compare("*Node"))
	{
		getline(infile,oneLine);
		istringstream st(oneLine);
		st >> oneWrod;
	}

	for (i = 0; i < vertexNumber; ++i)
	{
		getline(infile,oneLine);
		istringstream st(oneLine);
		st >> oneWrod;
		for (j = 0; j < 3; ++j)
		{
			st>>oneWrod;
			if (j < 2)
				oneWrod.erase(oneWrod.size()-1,1);
			istringstream str2Num(oneWrod);
			str2Num >> vertex[i][j];

		}
		
	}
	getline(infile,oneLine);
	 
	switch(elementType)
	{
	case 2:
		{
			for (i = 0; i < elementNumber; ++i)
			{
				getline(infile,oneLine);
				istringstream st(oneLine);
				st >> oneWrod;
				for (j = 0; j < 2; ++j)
				{
					st>>oneWrod;
					if (j < 1)
						oneWrod.erase(oneWrod.size()-1,1);
					istringstream str2Num(oneWrod);
					str2Num >> element[i][j];
					element[i][j]--;
				}
			}
			break;
		}
	case 3:
		{
			for (i = 0; i < elementNumber; ++i)
			{
				getline(infile,oneLine);
				istringstream st(oneLine);
				st >> oneWrod;
				for (j = 0; j < 3; ++j)
				{
					st>>oneWrod;
					if (j < 2)
						oneWrod.erase(oneWrod.size()-1,1);
					istringstream str2Num(oneWrod);
					str2Num >> element[i][j];
					element[i][j]--;
				}
			}
			break;
		}
	case 4:
		{
			for (i = 0; i < elementNumber; ++i)
			{
				getline(infile,oneLine);
				istringstream st(oneLine);
				st >> oneWrod;
				for (j = 0; j < 4; ++j)
				{
					st>>oneWrod;
					if (j < 3)
						oneWrod.erase(oneWrod.size()-1,1);
					istringstream str2Num(oneWrod);
					str2Num >> element[i][j];
					element[i][j]--;
				}
			}
			break;
		}
	case 8:
		{
			for (i = 0; i < elementNumber; ++i)
			{
				getline(infile,oneLine);
				istringstream st(oneLine);
				st >> oneWrod;
				for (j = 0; j < 8; ++j)
				{
					st>>oneWrod;
					if (j < 7)
						oneWrod.erase(oneWrod.size()-1,1);
					istringstream str2Num(oneWrod);
					str2Num >> element[i][j];
					element[i][j]--;
				}
			}
			break;
		}
	}

	infile.close(); 


	return true;
}

bool RawMesh::Read_Plt(const char * filename)
{
	FILE *input;
	int i, j, elementType;
	int v[8];
	 
	ifstream infile;
	string oneLine,oneWrod, formerLine;
	
	infile.open(filename);
	if (!infile)
	{
		cerr << "error: unable to open the file "<<ends;
		cerr <<filename <<endl;
		return false;
	}
	getline(infile,oneLine);
	getline(infile,oneLine);
	istringstream strStream(oneLine);
	strStream >> oneWrod >> oneWrod >> oneWrod >> vertexNumber >> oneWrod >> oneWrod >> oneWrod >> elementNumber;
	getline(infile,oneLine);
	strStream.clear();
	strStream.str(oneLine);
	strStream >> oneWrod >> oneWrod;
	oneWrod.erase(0,9);
	if (0 == oneWrod.compare("FETRIANGLE"))
	{
		elementType = 3;
	}
	else if (0 == oneWrod.compare("FEQUADRILATERAL"))
	{
		elementType = 5;
	}
	else if (0 == oneWrod.compare("FETETRAHEDRAL"))
	{
		elementType = 8;
	}
	else if (0 == oneWrod.compare("FEBRICK"))
	{
		elementType = 8;
	}
	else if (0 == oneWrod.compare("FELINESEG"))
	{
		elementType = 2;
	}
	

	switch(elementType)
	{
	case 2:
		{
			if (! CreateNewMesh(LINE, vertexNumber, elementNumber))
				return false;
			break;
		}
	case 3:
		{
			if (! CreateNewMesh(TRIANGLE, vertexNumber, elementNumber))
				return false;
			break;
		}
	case 4:
		{
			if (! CreateNewMesh(QUADRILATERAL, vertexNumber, elementNumber))
				return false;
			break;
		}
	case 8:
		{
			if (! CreateNewMesh(HEXAHEDRON, vertexNumber, elementNumber))
				return false;
			break;
		}
	}


	for (i = 0; i < vertexNumber; ++i)
	{
		getline(infile,oneLine);
		istringstream st(oneLine);
		st >> vertex[i][0] >> vertex[i][1] >> vertex[i][2];
		
	}
	getline(infile,oneLine);

	switch(elementType)
	{
	case 2:
		{
			for (i = 0; i < elementNumber; ++i)
			{
				getline(infile,oneLine);
				istringstream st(oneLine);
				for (j = 0; j < 2; ++j)
				{
					st>>v[j];
					element[i][j] = v[j]-1;
				}
			}
			break;
		}
	case 3:
		{
			for (i = 0; i < elementNumber; ++i)
			{
				getline(infile,oneLine);
				istringstream st(oneLine);
				for (j = 0; j < 3; ++j)
				{
					st>>v[j];
					element[i][j] = v[j]-1;
				}
			}
			break;
		}
	case 4:
		{
			for (i = 0; i < elementNumber; ++i)
			{
				getline(infile,oneLine);
				istringstream st(oneLine);
				for (j = 0; j < 4; ++j)
				{
					st>>v[j];
					element[i][j] = v[j]-1;
				}
			}
			break;
		}
	case 8:
		{
			for (i = 0; i < elementNumber; ++i)
			{
				getline(infile,oneLine);
				istringstream st(oneLine);
				for (j = 0; j < 8; ++j)
				{
					st>>v[j];
					element[i][j] = v[j]-1;
				}
			}
			break;
		}
	}

	infile.close(); 


	return true;
}

bool  RawMesh::ReadRaw_ReadVertex(FILE *&input, const char * type, bool withSign)
{
	int i, sign;
	double x, y, z, nx, ny, nz, r, g, b;

	if (withSign)
	{
		if(_stricmp(type, "raw") == 0)
			for(i = 0; i < vertexNumber; ++i)
			{
				fscanf(input, "%lf %lf %lf %d\n", &x, &y, &z, &sign);
				vertex[i][0] = x;
				vertex[i][1] = y;
				vertex[i][2] = z;
				vertexSign[i] = sign;
			}
		else if(_stricmp(type, "rawn") == 0)
		{
			if(! InitiateMatrix(vertexNormal, vertexNumber))
				return false;

			for(i = 0; i < vertexNumber; ++i)
			{
				fscanf(input, "%lf %lf %lf %lf %lf %lf %d\n", &x, &y, &z, &nx, &ny, &nz, &sign);
				vertex[i][0] = x;
				vertex[i][1] = y;
				vertex[i][2] = z;
				vertexNormal[i][0] = nx;
				vertexNormal[i][1] = ny;
				vertexNormal[i][2] = nz;
				vertexSign[i] = sign;
			}
		}
		else if(_stricmp(type, "rawc") == 0)
		{
			if(! InitiateMatrix(vertexColor, vertexNumber))
				return false;

			for(i = 0; i < vertexNumber; ++i)
			{
				fscanf(input, "%lf %lf %lf %lf %lf %lf %d\n", &x, &y, &z, &r, &g, &b, &sign);
				vertex[i][0] = x;
				vertex[i][1] = y;
				vertex[i][2] = z;
				vertexColor[i][0] = r;
				vertexColor[i][1] = g;
				vertexColor[i][2] = b;
				vertexSign[i] = sign;
			}
		}
		else if(_stricmp(type, "rawnc") == 0)
		{
			if(! InitiateMatrix(vertexNormal, vertexNumber))
				return false;
			if(! InitiateMatrix(vertexColor, vertexNumber))
				return false;

			for(i = 0; i < vertexNumber; ++i)
			{
				fscanf(input, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %d\n", &x, &y, &z, &nx, &ny, &nz, &r, &g, &b, &sign);
				vertex[i][0] = x;
				vertex[i][1] = y;
				vertex[i][2] = z;
				vertexNormal[i][0] = nx;
				vertexNormal[i][1] = ny;
				vertexNormal[i][2] = nz;
				vertexColor[i][0] = r;
				vertexColor[i][1] = g;
				vertexColor[i][2] = b;
				vertexSign[i] = sign;
			}
		}
		else
			return false;
	}
	else
	{
		if(_stricmp(type, "raw") == 0)
			for(i = 0; i < vertexNumber; ++i)
			{
				fscanf(input, "%lf %lf %lf\n", &x, &y, &z);
				vertex[i][0] = x;
				vertex[i][1] = y;
				vertex[i][2] = z;
			}
		else if(_stricmp(type, "rawn") == 0)
		{
			if(! InitiateMatrix(vertexNormal, vertexNumber))
				return false;

			for(i = 0; i < vertexNumber; ++i)
			{
				fscanf(input, "%lf %lf %lf %lf %lf %lf\n", &x, &y, &z, &nx, &ny, &nz);
				vertex[i][0] = x;
				vertex[i][1] = y;
				vertex[i][2] = z;
				vertexNormal[i][0] = nx;
				vertexNormal[i][1] = ny;
				vertexNormal[i][2] = nz;
			}
		}
		else if(_stricmp(type, "rawc") == 0)
		{
			if(! InitiateMatrix(vertexColor, vertexNumber))
				return false;

			for(i = 0; i < vertexNumber; ++i)
			{
				fscanf(input, "%lf %lf %lf %lf %lf %lf\n", &x, &y, &z, &r, &g, &b);
				vertex[i][0] = x;
				vertex[i][1] = y;
				vertex[i][2] = z;
				vertexColor[i][0] = r;
				vertexColor[i][1] = g;
				vertexColor[i][2] = b;
			}
		}
		else if(_stricmp(type, "rawnc") == 0)
		{
			if(! InitiateMatrix(vertexNormal, vertexNumber))
				return false;
			if(! InitiateMatrix(vertexColor, vertexNumber))
				return false;

			for(i = 0; i < vertexNumber; ++i)
			{
				fscanf(input, "%lf %lf %lf %lf %lf %lf %lf %lf %lf\n", &x, &y, &z, &nx, &ny, &nz, &r, &g, &b);
				vertex[i][0] = x;
				vertex[i][1] = y;
				vertex[i][2] = z;
				vertexNormal[i][0] = nx;
				vertexNormal[i][1] = ny;
				vertexNormal[i][2] = nz;
				vertexColor[i][0] = r;
				vertexColor[i][1] = g;
				vertexColor[i][2] = b;
			}
		}
		else
			return false;
	}

	return true;
}


bool RawMesh::ReadRaw_tet(const char *filename, const char *type)
{
	FILE *input;
	int i;
	int v0, v1, v2, v3;

	input = NULL;
	input = fopen(filename, "r");
	if (input == NULL)
	{
		//printf("Open file Error!\n");
		return false;
	}

	fscanf(input,"%d %d\n", &vertexNumber, &elementNumber);

	if (! CreateNewMesh(TETRAHEDRON, vertexNumber, elementNumber))
		return false;	

	if (! ReadRaw_ReadVertex(input, type))
		return false;

	for(i = 0; i < elementNumber; ++i)
	{
		fscanf(input, "%d %d %d %d\n", &v0, &v1, &v2, &v3);
		element[i][0] = v0;
		element[i][1] = v1;
		element[i][2] = v2;
		element[i][3] = v3;
	}

	fclose(input);

	return true;
}

bool RawMesh::ReadRaw_tri(const char *filename, const char *type)
{
	FILE *input;
	int i;
	int v0, v1, v2;

	input = NULL;
	input = fopen(filename, "r");
	if(input == NULL)
	{
		//printf("Open file Error!\n");
		return false;
	}

	fscanf(input,"%d %d\n", &vertexNumber, &elementNumber);

	if (! CreateNewMesh(TRIANGLE, vertexNumber, elementNumber))
		return false;	

	if (! ReadRaw_ReadVertex(input, type, false))
		return false;

	for(i = 0; i < elementNumber; ++i)
	{
		fscanf(input, "%d %d %d\n", &v0, &v1, &v2);
		element[i][0] = v0;
		element[i][1] = v1;
		element[i][2] = v2;
	}

	fclose(input);

	return true;
}


bool RawMesh::ReadRaw_quad(const char *filename, const char *type)
{
	FILE *input;
	int i;
	int v0, v1, v2, v3;


	input = NULL;
	input = fopen(filename, "r");
	if(input == NULL)
	{
		//printf("Open file Error!\n");
		return false;
	}

	fscanf(input,"%d %d\n", &vertexNumber, &elementNumber);

	if (! CreateNewMesh(QUADRILATERAL, vertexNumber, elementNumber))
		return false;	

	if (! ReadRaw_ReadVertex(input, type, false))
		return false;

	for(i = 0; i < elementNumber; ++i)
	{
		fscanf(input, "%d %d %d %d\n", &v0, &v1, &v2, &v3);
		element[i][0] = v0;
		element[i][1] = v1;
		element[i][2] = v2;
		element[i][3] = v3;
	}

	fclose(input);

	return true;
}

bool RawMesh::ReadRaw_line(const char *filename, const char *type)
{
	FILE *input;
	int i;
	int v0, v1;


	input = NULL;
	input = fopen(filename, "r");
	if(input == NULL)
	{
		//printf("Open file Error!\n");
		return false;
	}

	fscanf(input,"%d %d\n", &vertexNumber, &elementNumber);

	if (! CreateNewMesh(LINE, vertexNumber, elementNumber))
		return false;	

	if (! ReadRaw_ReadVertex(input, type))
		return false;

	for(i = 0; i < elementNumber; ++i)
	{
		fscanf(input, "%d %d\n", &v0, &v1);
		element[i][0] = v0;
		element[i][1] = v1;
	}

	fclose(input);

	return true;
}

bool RawMesh::ReadRaw_hex(const char *filename, const char *type)
{
	FILE *input;
	int i;
	int v0, v1, v2, v3, v4, v5, v6, v7;


	input = NULL;
	input = fopen(filename, "r");
	if(input == NULL)
	{
		//printf("Open file Error!\n");
		return false;
	}

	fscanf(input,"%d %d\n", &vertexNumber, &elementNumber);

	if (! CreateNewMesh(HEXAHEDRON, vertexNumber, elementNumber))
		return false;	

	if (! ReadRaw_ReadVertex(input, type))
		return false;

	for(i = 0; i < elementNumber; ++i)
	{
		fscanf(input, "%d %d %d %d %d %d %d %d\n", &v0, &v1, &v2, &v3, &v4, &v5, &v6, &v7);
		element[i][0] = v0;
		element[i][1] = v1;
		element[i][2] = v2;
		element[i][3] = v3;
		element[i][4] = v4;
		element[i][5] = v5;
		element[i][6] = v6;
		element[i][7] = v7;
	}

	fclose(input);

	return true;
}

bool RawMesh::Read_Mesh(const char *filename)
{
	FILE *input;
	int i, j;
	int v0,nHex, nQuad, nTet, nTri, nLine, nPyr, sign;
	double x, y, z;


	input = NULL;
	input = fopen(filename, "r");
	if(input == NULL)
	{
		//printf("Open file Error!\n");
		return false;
	}

	fscanf(input,"%d %d %d %d %d %d %d\n", &vertexNumber, &nHex, &nPyr, &nTet, &nQuad, &nTri, &nLine);

	if (nHex != 0)
	{
		elementNumber = nHex;
		SetElementType(HEXAHEDRON);
	}
	else if (nPyr != 0)
	{
		elementNumber = nPyr;
		return false;
	}
	else if (nTet != 0)
	{
		elementNumber = nTet;
		SetElementType(TETRAHEDRON);
	}
	else if (nQuad != 0)
	{
		elementNumber = nQuad;
		SetElementType(QUADRILATERAL);
	}
	else if (nTri != 0)
	{
		elementNumber = nTri;
		SetElementType(TRIANGLE);
	}
	else
	{
		elementNumber = nLine;
		SetElementType(LINE);
	}


	if(! CreateNewMesh(elementProperty.elementType, vertexNumber, elementNumber))
		return false;	

	for(i = 0; i < vertexNumber; ++i)
	{
		fscanf(input, "%lf %lf %lf %d\n", &x, &y, &z, &sign);
		vertex[i][0] = x;
		vertex[i][1] = y;
		vertex[i][2] = z;
		vertexSign[i] = sign;
	}


	for(i = 0; i < elementNumber; ++i)
	{
		for (j=0; j<elementProperty.vertexNumber; j++)
		{
			fscanf(input, "%d", &v0);
			element[i][j] = v0;
		}
		fscanf(input, "\n");

		elementSign[i] = 0;
	}

	fclose(input);

	return true;
}


bool RawMesh::Write(const char *filename)
{
	string::size_type fileLength;
	bool writeSuccess;
	string sTmpName;
	

	if(vertexNumber == 0 || vertex == NULL || elementNumber == 0 || element == NULL)
	{
		//printf("Error! EMPTY Vertex or Element.\n");
		return false;
	}

	char * cTmp;

	cTmp = _strupr(_strdup(filename));
	sTmpName += cTmp;
	fileLength = sTmpName.size();

	//Write files
	if		(sTmpName.find(".RAWM") == (fileLength - 5))
		writeSuccess = Write_RawMesh(filename);
	else if (sTmpName.find(".MESH") == (fileLength - 5))
		writeSuccess = Write_Mesh(filename);
	else if (sTmpName.find(".RAW") < fileLength)
		writeSuccess = Write_RawFile(filename);
	else if (sTmpName.find(".PLY") < fileLength)
		writeSuccess = Write_Ply(filename);
	else if (sTmpName.find(".PLT") < fileLength)
		writeSuccess = Write_Plt(filename);
	else if (sTmpName.find(".INP") < fileLength)
		writeSuccess = Write_Inp(filename);
	else if (sTmpName.find(".VTK") < fileLength)
		writeSuccess = Write_Vtk(filename);

	else
		writeSuccess = false;

	//Clear
	sTmpName.clear();
	delete [] cTmp;

	return writeSuccess;
}

bool RawMesh::Write_RawMesh(const char *filename)
{
	version = "2.0.20120515";
	if (version.compare("1.0") == 0)
	{
		return Write_RawMesh_OldVersion(filename);
	}

	FILE *output;
	int i, j;
	
	output = fopen(filename,"w");

	fprintf(output, "Version %s\n", version.c_str());
	fprintf(output, "%d %d %d\n", vertexNumber, elementNumber, elementProperty.elementType);
	fprintf(output, "%d %d %d %d\n", meshInfo.isUseVertexNormal, meshInfo.isUseVertexColor, meshInfo.isUseElementNormal, meshInfo.isUseElementColor);

	for (i=0; i<vertexNumber; ++i)
		fprintf(output, "%lf %lf %lf %d\n", vertex[i][0], vertex[i][1], vertex[i][2], vertexSign[i]);

	if (meshInfo.isUseVertexNormal)
	{
		if (vertexNormal == NULL)
			ComputeVertexNormal();

		for (i=0; i<vertexNumber; ++i)
			fprintf(output, "%lf %lf %lf\n", vertexNormal[i][0], vertexNormal[i][1], vertexNormal[i][2]);
	}

	if (meshInfo.isUseVertexColor)
	{
		if (vertexColor != NULL)
		{
			for (i=0; i<vertexNumber; ++i)
				fprintf(output, "%lf %lf %lf\n", vertexColor[i][0], vertexColor[i][1], vertexColor[i][2]);
		}
		else
		{
			for (i=0; i<vertexNumber; ++i)
				fprintf(output, "%lf %lf %lf\n", defaultColor[0], defaultColor[1], defaultColor[2]);
		}
	}

	for (i=0; i<elementNumber; ++i)
	{
		for (j=0; j<elementProperty.vertexNumber; ++j)
			fprintf(output, "%d ", element[i][j]);
		fprintf(output, "%d\n", elementSign[i]);
	}

	if (meshInfo.isUseElementNormal)
	{
		if (elementNormal == NULL)
			ComputeElementNormal();

		for (i=0; i<elementNumber; ++i)
			fprintf(output, "%lf %lf %lf\n", elementNormal[i][0], elementNormal[i][1], elementNormal[i][2]);
	}

	if (meshInfo.isUseElementColor)
	{
		if (elementColor != NULL)
		{
			for (i=0; i<elementNumber; ++i)
				fprintf(output, "%lf %lf %lf\n", elementColor[i][0], elementColor[i][1], elementColor[i][2]);
		}
		else
		{
			for (i=0; i<elementNumber; ++i)
				fprintf(output, "%lf %lf %lf\n", defaultColor[0], defaultColor[1], defaultColor[2]);
		}
	}

	fclose(output);
	return true;
}

bool RawMesh::Write_RawMesh_OldVersion(const char *filename)
{
	FILE *output;
	int i, j;
	
	output = fopen(filename,"w");

	fprintf(output, "%d %d %d\n", vertexNumber, elementNumber, elementProperty.elementType);

	for (i=0; i<vertexNumber; ++i)
		fprintf(output, "%lf %lf %lf %d\n", vertex[i][0], vertex[i][1], vertex[i][2], vertexSign[i]);

	for (i=0; i<elementNumber; ++i)
	{
		for (j=0; j<elementProperty.vertexNumber; ++j)
			fprintf(output, "%d ", element[i][j]);
		fprintf(output, "%d\n", elementSign[i]);
	}

	fclose(output);
	return true;
}

bool RawMesh::Write_RawFile(const char *filename)
{
	int i, fileLength, fdot;
	bool writeSuccess;
	char *fileType;

	//Find the suffix of the file
	fileLength = strlen(filename);
	for(i = fileLength-1; i >= 0; i--)
		if(filename[i] == '.') break;
	if(i < 0)
	{
		//printf("Filename Error!\n");
		return false;
	}
	fdot = i+1;
	fileType = new char[fileLength-fdot+1];
	for(i = fdot; i <= fileLength; ++i)
		fileType[i-fdot] = filename[i];

	//Write files
	if		(elementProperty.elementType == TRIANGLE)
		writeSuccess = WriteRaw_tri(filename, fileType);
	else if	(elementProperty.elementType == TETRAHEDRON)
		writeSuccess = WriteRaw_tet(filename, fileType);
	else if	(elementProperty.elementType == QUADRILATERAL)
		writeSuccess = WriteRaw_quad(filename, fileType);
	else if	(elementProperty.elementType == HEXAHEDRON)
		writeSuccess = WriteRaw_hex(filename, fileType);
	else if	(elementProperty.elementType == LINE)
		writeSuccess = WriteRaw_line(filename, fileType);
	else
		writeSuccess = false;

	//Clear
	delete [] fileType;

	return writeSuccess;
}


bool RawMesh::WriteRaw_WriteVertex(FILE *&output, const char * type, bool withSign)
{
	int i;

	if (withSign)
	{
		if(_stricmp(type, "raw") == 0)
			for(i = 0; i < vertexNumber; ++i)
				fprintf(output, "%lf %lf %lf %d\n", vertex[i][0], vertex[i][1], vertex[i][2], vertexSign[i]);
		else if(_stricmp(type, "rawn") == 0)
		{
			if(NULL == vertexNormal)
				ComputeVertexNormal();
			if(vertexNormal != NULL)
				for(i = 0; i < vertexNumber; ++i)
					fprintf(output, "%lf %lf %lf %lf %lf %lf %d\n", vertex[i][0], vertex[i][1], vertex[i][2],
						vertexNormal[i][0], vertexNormal[i][1], vertexNormal[i][2], vertexSign[i]);
			else
				for(i = 0; i < vertexNumber; ++i)
					fprintf(output, "%lf %lf %lf %lf %lf %lf %d\n", vertex[i][0], vertex[i][1], vertex[i][2],
						defaultNormal[0], defaultNormal[1], defaultNormal[2], vertexSign[i]);
		}
		else if(_stricmp(type, "rawc") == 0)
		{
			if(vertexColor != NULL)
				for(i = 0; i < vertexNumber; ++i)
					fprintf(output, "%lf %lf %lf %lf %lf %lf %d\n", vertex[i][0], vertex[i][1], vertex[i][2],
						vertexColor[i][0], vertexColor[i][1], vertexColor[i][2], vertexSign[i]);
			else
				for(i = 0; i < vertexNumber; ++i)
					fprintf(output, "%lf %lf %lf %lf %lf %lf %d\n", vertex[i][0], vertex[i][1], vertex[i][2],
						defaultColor[0], defaultColor[1], defaultColor[2], vertexSign[i]);
		}
		else if(_stricmp(type, "rawnc") == 0)
		{
			if(NULL == vertexNormal)
				ComputeVertexNormal();
			if(vertexNormal != NULL && vertexColor != NULL)
				for(i = 0; i < vertexNumber; ++i)
					fprintf(output, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %d\n", vertex[i][0], vertex[i][1], vertex[i][2],
						vertexNormal[i][0], vertexNormal[i][1], vertexNormal[i][2], vertexColor[i][0], vertexColor[i][1], vertexColor[i][2], vertexSign[i]);
			else if(vertexNormal != NULL && vertexColor == NULL)
				for(i = 0; i < vertexNumber; ++i)
					fprintf(output, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %d\n", vertex[i][0], vertex[i][1], vertex[i][2],
						vertexNormal[i][0], vertexNormal[i][1], vertexNormal[i][2], defaultColor[0], defaultColor[1], defaultColor[2], vertexSign[i]);
			else if(vertexNormal == NULL && vertexColor != NULL)
				for(i = 0; i < vertexNumber; ++i)
					fprintf(output, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %d\n", vertex[i][0], vertex[i][1], vertex[i][2],
						defaultNormal[0], defaultNormal[1], defaultNormal[2], vertexColor[i][0], vertexColor[i][1], vertexColor[i][2], vertexSign[i]);
			else
				for(i = 0; i < vertexNumber; ++i)
					fprintf(output, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %d\n", vertex[i][0], vertex[i][1], vertex[i][2],
						defaultNormal[0], defaultNormal[1], defaultNormal[2], defaultColor[0], defaultColor[1], defaultColor[2], vertexSign[i]);
		}
		else
			return false;
	}
	else
	{
		if(_stricmp(type, "raw") == 0)
			for(i = 0; i < vertexNumber; ++i)
				fprintf(output, "%lf %lf %lf\n", vertex[i][0], vertex[i][1], vertex[i][2]);
		else if(_stricmp(type, "rawn") == 0)
		{
			if(NULL == vertexNormal)
				ComputeVertexNormal();
			if(vertexNormal != NULL)
				for(i = 0; i < vertexNumber; ++i)
					fprintf(output, "%lf %lf %lf %lf %lf %lf\n", vertex[i][0], vertex[i][1], vertex[i][2],
						vertexNormal[i][0], vertexNormal[i][1], vertexNormal[i][2]);
			else
				for(i = 0; i < vertexNumber; ++i)
					fprintf(output, "%lf %lf %lf %lf %lf %lf\n", vertex[i][0], vertex[i][1], vertex[i][2],
						defaultNormal[0], defaultNormal[1], defaultNormal[2]);
		}
		else if(_stricmp(type, "rawc") == 0)
		{
			if(vertexColor != NULL)
				for(i = 0; i < vertexNumber; ++i)
					fprintf(output, "%lf %lf %lf %lf %lf %lf\n", vertex[i][0], vertex[i][1], vertex[i][2],
						vertexColor[i][0], vertexColor[i][1], vertexColor[i][2]);
			else
				for(i = 0; i < vertexNumber; ++i)
					fprintf(output, "%lf %lf %lf %lf %lf %lf\n", vertex[i][0], vertex[i][1], vertex[i][2],
						defaultColor[0], defaultColor[1], defaultColor[2]);
		}
		else if(_stricmp(type, "rawnc") == 0)
		{
			if(NULL == vertexNormal)
				ComputeVertexNormal();
			if(vertexNormal != NULL && vertexColor != NULL)
				for(i = 0; i < vertexNumber; ++i)
					fprintf(output, "%lf %lf %lf %lf %lf %lf %lf %lf %lf\n", vertex[i][0], vertex[i][1], vertex[i][2],
						vertexNormal[i][0], vertexNormal[i][1], vertexNormal[i][2], vertexColor[i][0], vertexColor[i][1], vertexColor[i][2]);
			else if(vertexNormal != NULL && vertexColor == NULL)
				for(i = 0; i < vertexNumber; ++i)
					fprintf(output, "%lf %lf %lf %lf %lf %lf %lf %lf %lf\n", vertex[i][0], vertex[i][1], vertex[i][2],
						vertexNormal[i][0], vertexNormal[i][1], vertexNormal[i][2], defaultColor[0], defaultColor[1], defaultColor[2]);
			else if(vertexNormal == NULL && vertexColor != NULL)
				for(i = 0; i < vertexNumber; ++i)
					fprintf(output, "%lf %lf %lf %lf %lf %lf %lf %lf %lf\n", vertex[i][0], vertex[i][1], vertex[i][2],
						defaultNormal[0], defaultNormal[1], defaultNormal[2], vertexColor[i][0], vertexColor[i][1], vertexColor[i][2]);
			else
				for(i = 0; i < vertexNumber; ++i)
					fprintf(output, "%lf %lf %lf %lf %lf %lf %lf %lf %lf\n", vertex[i][0], vertex[i][1], vertex[i][2],
						defaultNormal[0], defaultNormal[1], defaultNormal[2], defaultColor[0], defaultColor[1], defaultColor[2]);
		}
		else
			return false;
	}

	return true;
}

bool RawMesh::WriteRaw_tri(const char *filename, const char *type)
{
	FILE *output;
	int i;

	output = 0;
	output = fopen(filename, "w");
	if(output == 0)
	{
		//printf("Open file Error!\n");
		return false;
	}

	fprintf(output, "%d %d\n", vertexNumber, elementNumber);

	WriteRaw_WriteVertex(output, type, false);

	for(i = 0; i < elementNumber; ++i)
		fprintf(output, "%d %d %d\n", element[i][0], element[i][1], element[i][2]);

	fclose(output);

	return true;
}

bool RawMesh::WriteRaw_tet(const char *filename, const char *type)
{
	FILE *output;
	int i;

	output = NULL;
	output = fopen(filename, "w");
	if(output == NULL)
	{
		//printf("Open file Error!\n");
		return false;
	}

	fprintf(output, "%d %d\n", vertexNumber, elementNumber);

	WriteRaw_WriteVertex(output, type);

	for(i = 0; i < elementNumber; ++i)
		fprintf(output, "%d %d %d %d\n", element[i][0], element[i][1], element[i][2], element[i][3]);

	fclose(output);

	return true;
}

bool RawMesh::WriteRaw_quad(const char *filename, const char *type)
{
	FILE *output;
	int i;

	output = NULL;
	output = fopen(filename, "w");
	if(output == NULL)
	{
		//printf("Open file Error!\n");
		return false;
	}

	fprintf(output, "%d %d\n", vertexNumber, elementNumber);

	WriteRaw_WriteVertex(output, type, false);

	for(i = 0; i < elementNumber; ++i)
		fprintf(output, "%d %d %d %d\n", element[i][0], element[i][1], element[i][2], element[i][3]);

	fclose(output);

	return true;
}

bool RawMesh::WriteRaw_hex(const char *filename, const char *type)
{
	FILE *output;
	int i;

	output = NULL;
	output = fopen(filename, "w");
	if(output == NULL)
	{
		//printf("Open file Error!\n");
		return false;
	}

	fprintf(output, "%d %d\n", vertexNumber, elementNumber);

	WriteRaw_WriteVertex(output, type);

	for(i = 0; i < elementNumber; ++i)
		fprintf(output, "%d %d %d %d %d %d %d %d\n", element[i][0], element[i][1], element[i][2], element[i][3], 
			element[i][4], element[i][5], element[i][6], element[i][7]);

	fclose(output);

	return true;
}

bool RawMesh::WriteRaw_line(const char *filename, const char *type)
{
	FILE *output;
	int i;

	output = NULL;
	output = fopen(filename, "w");
	if(output == NULL)
	{
		//printf("Open file Error!\n");
		return false;
	}

	fprintf(output, "%d %d\n", vertexNumber, elementNumber);

	WriteRaw_WriteVertex(output, type);

	for(i = 0; i < elementNumber; ++i)
		fprintf(output, "%d %d\n", element[i][0], element[i][1]);

	fclose(output);

	return true;
}

bool RawMesh::Write_Mesh(const char *filename)
{
	FILE *output;
	int i, j, nHex=0, nQuad=0, nTet=0, nTri=0, nLine=0, nPyr=0;

	output = NULL;
	output = fopen(filename, "w");
	if(output == NULL)
	{
		//printf("Open file Error!\n");
		return false;
	}

	if (elementProperty.elementType == HEXAHEDRON)
		nHex = elementNumber;
	else if (elementProperty.elementType == TETRAHEDRON)
		nTet = elementNumber;
	else if (elementProperty.elementType == QUADRILATERAL)
		nQuad = elementNumber;
	else if (elementProperty.elementType == TRIANGLE)
		nTri = elementNumber;
	else if (elementProperty.elementType == LINE)
		nLine = elementNumber;

	fprintf(output, "%d %d %d %d %d %d %d\n", vertexNumber, nHex, nPyr, nTet, nQuad, nTri, nLine);

	for (i=0; i<vertexNumber; ++i)
	{
		fprintf(output, "%.9lf %.9lf %.9lf %d\n", vertex[i][0], vertex[i][1], vertex[i][2], vertexSign[i]);
	}
	for (i=0; i<elementNumber; ++i)
	{
		for(j=0; j<elementProperty.vertexNumber; j++)
		{
			fprintf(output, "%d ", element[i][j]);
		}
		fprintf(output, "\n");
	}

	fclose(output);
	return true;
}

bool RawMesh::Write_Inp(const char *filename)
{
	string inputFileName = filename;

	fstream output(inputFileName, fstream::out);
	if (!output)
	{
		cout << "open input file error!!!" <<endl;
		return false;
	}
	int i, j;

	output << "*Heading" << "\n";
	output << "** Job name: TSpline Model name: Input"<<"\n";
	output << "**"<<"\n";
	output << "** PARTS"<<"\n";
	output << "**"<<"\n";
	output << "*Part, name=Input"<<"\n";
	output << "*Node"<<"\n";
	
	output.precision(9);
	for (i=0; i<vertexNumber; ++i)
	{
		output << i+1<<", ";

		output <<vertex[i][0]<<", "<<vertex[i][1]<<", "<<vertex[i][2]<<"\n";

	}
	output.unsetf(ostream::floatfield);
	
	if (elementProperty.elementType == HEXAHEDRON)
		output << "*Element, type=C3D8"<<"\n";
	else if (elementProperty.elementType == QUADRILATERAL)
		output << "*Element, type=S4R"<<"\n";
	else if (elementProperty.elementType == TETRAHEDRON)
		output << "*Element, type=C3D4"<<"\n";
	else if (elementProperty.elementType == TRIANGLE)
		output << "*Element, type=CPS3"<<"\n";
	else if (elementProperty.elementType == LINE)
		output << "*Element, type=F2D2"<<"\n";

	
	for (i=0; i<elementNumber; ++i)
	{
		output << i+1<<", ";
		for(j=0; j<elementProperty.vertexNumber-1; j++)
		{
			output<< element[i][j]+1<<", ";
		}
		output<< element[i][j]+1<<"\n";
	}
	output << "*End Part"<<"\n";
	output.close();


	return true;
}



bool RawMesh::Write_Plt(const char *filename)
{
	string inputFileName = filename;

	fstream output(inputFileName, fstream::out);
	if (!output)
	{
		cout << "open input file error!!!" <<endl;
		return false;
	}
	int i, j;

	output << "VARIABLES = \"X\" \"Y\" \"Z\"" << "\n";
	output << "ZONE N = "<<vertexNumber << " , E = "<<elementNumber<<"\n";
	output << "DATAPACKING=POINT, ZONETYPE=";

	switch(elementProperty.vertexNumber)
	{
	case 2:
		output <<"FELINESEG"<<"\n";
		break;
	case 3:
		output <<"FETRIANGLE"<<"\n";
		break;
	case 4:
		output <<"FETETRAHEDRON,"<<"\n";
		break;
	case 5:
		output <<"FEQUADRILATERAL"<<"\n";
		break;
	case 8:
		output <<"FEBRICK"<<"\n";
		break;
	}
	 

	output.precision(9);
	for (i=0; i<vertexNumber; ++i)
		output <<vertex[i][0]<<"\t"<<vertex[i][1]<<"\t"<<vertex[i][2]<<"\n";
	output.unsetf(ostream::floatfield);

	
	for (i=0; i<elementNumber; ++i)
	{
		for(j=0; j<elementProperty.vertexNumber; j++)
		{
			output<< element[i][j]+1<<" ";
		}
		output<<"\n";
	}
	output.close();


	return true;
}

bool RawMesh::Write_Ply(const char *filename)

{
	string inputFileName = filename;

	fstream output(inputFileName, fstream::out);
	if (!output)
	{
		cout << "open input file error!!!" <<endl;
		return false;
	}
	if (elementProperty.vertexNumber>4)
		cout << "error !!!\n solid mesh cannot be write to ply format...\n";
	int i, j;

	output << "ply" << "\n";
	output << "format ascii 1.0" << "\n";
	output << "comment File exported by Rhinoceros Version 4.0" << "\n";
	output << "element vertex " <<vertexNumber << "\n";
	output << "property float x" << "\n";
	output << "property float y" << "\n";
	output << "property float z" << "\n";
	output << "element face " <<elementNumber<< "\n";
	output << "property list uchar uint vertex_index" << "\n";
	output << "end_header" << "\n";


	output.precision(9);
	for (i=0; i<vertexNumber; ++i)
		output <<vertex[i][0]<<"\t"<<vertex[i][1]<<"\t"<<vertex[i][2]<<"\n";
	output.unsetf(ostream::floatfield);
	

	for (i=0; i<elementNumber; ++i)
	{
		output << elementProperty.vertexNumber<<" ";
		for(j=0; j<elementProperty.vertexNumber; j++)
		{
			output<< element[i][j]<<" ";
		}
		output<<"\n";
	}
	output.close();


	return true;
}

bool RawMesh::Write_Vtk(const char *filename)

{
	string inputFileName = filename;

	fstream output(inputFileName, fstream::out);
	if (!output)
	{
		cout << "open input file error!!!" <<endl;
		return false;
	}
	
	int i, j;

	output << "# vtk DataFile Version 3.1 " << "\n";
	output << "for LSEConsole" << "\n";
	output << "ASCII" << "\n";
	output << "DATASET UNSTRUCTURED_GRID" << "\n";
	output << "POINTS " << vertexNumber<<" FLOAT"<<"\n";
	

	output.precision(9);
	for (i=0; i<vertexNumber; ++i)
		output <<vertex[i][0]<<"\t"<<vertex[i][1]<<"\t"<<vertex[i][2]<<"\n";
	output.unsetf(ostream::floatfield);
	output << "CELLS " << elementNumber <<" "<<(elementProperty.vertexNumber+1) * elementNumber<<"\n";

	for (i=0; i<elementNumber; ++i)
	{
		output << elementProperty.vertexNumber<<" ";
		for(j=0; j<elementProperty.vertexNumber; j++)
		{
			output<< element[i][j]<<" ";
		}
		output<<"\n";
	}
	output <<"CELL_TYPES "<<elementNumber<<endl;
	for (i=0; i<elementNumber; ++i)
	{
		if (elementProperty.elementType == NA)
		{
			continue;
		}
		else if (elementProperty.elementType == POINT)
		{
			output << 1 <<endl;
		}
		else if (elementProperty.elementType == LINE)
		{
			output << 3 <<endl;
		}
		else if (elementProperty.elementType == TRIANGLE)
		{
			output << 5 <<endl;
		}
		else if (elementProperty.elementType == TETRAHEDRON)
		{
			output << 10 <<endl;
		}
		else if (elementProperty.elementType == QUADRILATERAL)
		{
			output << 9 <<endl;
		}
		else if (elementProperty.elementType == HEXAHEDRON)
		{
			output << 12 <<endl;
		}
		else if (elementProperty.elementType == POLYGON)
		{
			output << 7 <<endl;
		}
		else if (elementProperty.elementType == HEXAGON)
		{
			continue;
		}
		
	}
	output.close();


	return true;
}

bool RawMesh::InitiateElementValence()
{
	//calculate the vertex valence
	int i, j, vertexIndex, sum, elemVertNum;

	if (IsEmpty())
	{
		OutputError("RawMesh::GetElementValence","Empty Vertex or Element");
		return false;
	}

	//allocate space to edgeValence and edgeValenceNum
	if(! InitiateElementValenceMember())
	{
		OutputError("RawMesh::GetElementValence","Failed to initiate elementValence");
		return false;
	}

	//Begin
	elemVertNum = elementProperty.vertexNumber;
	for(i = 0; i < elementNumber; ++i)
		for(j = 0; j < elemVertNum; ++j)
		{
			vertexIndex = element[i][j];
			if (vertexIndex != -1)
				elementValenceNumber[vertexIndex]++;
		}
	

	//set where the elementValence[i] point to
	sum = 0;
	for(i = 1; i < vertexNumber; ++i)
	{
		elementValence[i] = elementValence[i-1] + elementValenceNumber[i-1];
		sum += elementValenceNumber[i-1];
		elementValenceNumber[i-1] = 0;
	}
	sum += elementValenceNumber[vertexNumber-1];
	elementValenceNumber[vertexNumber-1] = 0;

	//Store
	for(i = 0; i < elementNumber; ++i)
	{
		for(j = 0; j < elemVertNum; ++j)
		{
			vertexIndex = element[i][j];

			if (i == 2850)
				i = i;

			if (vertexIndex != -1)
				elementValence[vertexIndex][elementValenceNumber[vertexIndex]++] = i;
		}
	}
	//End

	if(sum != elementNumber*elemVertNum)
		return false;

	return true;
}


bool RawMesh::InitiateEdgeValence()
{
	if(IsEmpty())
	{
		OutputError("RawMesh::GetEdgeValence","Empty Vertex or Element");
		return false;
	}

	int iNum, i, j, k, iCount, point_1, point_2, elementEdgeNum;

	elementEdgeNum = elementProperty.edgeNumber;

	iNum = elementNumber * elementEdgeNum*2;
	if (!InitiateEdgeValenceMember(iNum))
	{
		OutputError("RawMesh::GetEdgeValence","Failed to initiate edgeValence and edgeValenceNum");
		return false;
	}

	//Store the edge valence information
	for (i=0; i<elementNumber; ++i)
	{
		for (j=0; j<elementEdgeNum; ++j)
		{
			GetElementEdge(i, j, point_1, point_2);
			if (point_1 < 0 || point_2 < 0)
				continue;

			edgeValenceNumber[point_1]++;
			edgeValenceNumber[point_2]++;
		}
	}

	for (i=1; i<vertexNumber; i++)
		edgeValence[i] = edgeValence[i-1] + edgeValenceNumber[i-1];

	for (i=0; i<elementNumber; ++i)
	{
		for (j=0; j<elementEdgeNum; ++j)
		{
			GetElementEdge(i, j, point_1, point_2);
			if (point_1 < 0 || point_2 < 0)
				continue;

			for (k=0; k<edgeValenceNumber[point_1]; ++k)
				if (edgeValence[point_1][k] == -1)
				{
					edgeValence[point_1][k] = point_2;
					break;
				}
			for (k=0; k<edgeValenceNumber[point_2]; ++k)
				if (edgeValence[point_2][k] == -1)
				{
					edgeValence[point_2][k] = point_1;
					break;
				}
		}
	}

	//Delete some unused space, restore the edge valence information, and then reduce used memory
	int *newEdgeValenceNum=NULL;
	if (!InitiateMatrix(newEdgeValenceNum, vertexNumber))
	{
		OutputError("RawMesh::GetEdgeValence","Failed to initiate newEdgeValenceNum");
		return false;
	}

	iNum = 0;
	for (i=0; i<vertexNumber; ++i)
	{
		newEdgeValenceNum[i] = edgeValenceNumber[i];

		for (j=0; j<edgeValenceNumber[i]; ++j)
		{
			if (edgeValence[i][j] == -1)
				continue;

			for (k=j+1; k<edgeValenceNumber[i]; ++k)
				if (edgeValence[i][k] == edgeValence[i][j])
				{
					edgeValence[i][k] = -1;
					newEdgeValenceNum[i]--;
				}
		}

		iNum += newEdgeValenceNum[i];
	}

	int **newEdgeValence;
	newEdgeValence = new int *[vertexNumber];
	newEdgeValence[0] = new int [iNum];
	if (newEdgeValence == NULL || newEdgeValence[0] == NULL)
	{
		OutputError("RawMesh::GetEdgeValence","Failed to initiate newEdgeValence");
		return false;
	}
	for (i=1; i<vertexNumber; i++)
		newEdgeValence[i] = newEdgeValence[i-1] + newEdgeValenceNum[i-1];

	for (i=0; i<vertexNumber; i++)
	{
		iCount = 0;
		for (j=0; j<edgeValenceNumber[i]; j++)
			if (edgeValence[i][j] != -1)
				newEdgeValence[i][iCount++] = edgeValence[i][j];
	}

	delete [] edgeValence[0];
	delete [] edgeValence;
	delete [] edgeValenceNumber;
	edgeValence = newEdgeValence;
	edgeValenceNumber = newEdgeValenceNum;	
	
	return true;
}

void RawMesh::GetElementEdge(int elementID, int edgeID, int &point_1, int &point_2)
{
	if (elementProperty.elementType == HEXAHEDRON)
	{
		if (edgeID == 0)
			{point_1 = element[elementID][0]; point_2 = element[elementID][1];}
		else if (edgeID == 1)
			{point_1 = element[elementID][1]; point_2 = element[elementID][2];}
		else if (edgeID == 2)
			{point_1 = element[elementID][2]; point_2 = element[elementID][3];}
		else if (edgeID == 3)
			{point_1 = element[elementID][3]; point_2 = element[elementID][0];}
		else if (edgeID == 4)
			{point_1 = element[elementID][0]; point_2 = element[elementID][4];}
		else if (edgeID == 5)
			{point_1 = element[elementID][1]; point_2 = element[elementID][5];}
		else if (edgeID == 6)
			{point_1 = element[elementID][2]; point_2 = element[elementID][6];}
		else if (edgeID == 7)
			{point_1 = element[elementID][3]; point_2 = element[elementID][7];}
		else if (edgeID == 8)
			{point_1 = element[elementID][4]; point_2 = element[elementID][5];}
		else if (edgeID == 9)
			{point_1 = element[elementID][5]; point_2 = element[elementID][6];}
		else if (edgeID == 10)
			{point_1 = element[elementID][6]; point_2 = element[elementID][7];}
		else if (edgeID == 11)
			{point_1 = element[elementID][7]; point_2 = element[elementID][4];}
		else
			OutputError("RawMesh::GetElementEdge", "Wrong edgeID", edgeID);
	}
	else if (elementProperty.elementType == QUADRILATERAL)
	{
		if (edgeID == 0)
			{point_1 = element[elementID][0]; point_2 = element[elementID][1];}
		else if (edgeID == 1)
			{point_1 = element[elementID][1]; point_2 = element[elementID][2];}
		else if (edgeID == 2)
			{point_1 = element[elementID][2]; point_2 = element[elementID][3];}
		else if (edgeID == 3)
			{point_1 = element[elementID][3]; point_2 = element[elementID][0];}
		else
			OutputError("RawMesh::GetElementEdge", "Wrong edgeID", edgeID);
	}
	else if (elementProperty.elementType == TETRAHEDRON)
	{
		if (edgeID == 0)
			{point_1 = element[elementID][0]; point_2 = element[elementID][1];}
		else if (edgeID == 1)
			{point_1 = element[elementID][0]; point_2 = element[elementID][2];}
		else if (edgeID == 2)
			{point_1 = element[elementID][0]; point_2 = element[elementID][3];}
		else if (edgeID == 3)
			{point_1 = element[elementID][1]; point_2 = element[elementID][2];}
		else if (edgeID == 4)
			{point_1 = element[elementID][1]; point_2 = element[elementID][3];}
		else if (edgeID == 5)
			{point_1 = element[elementID][2]; point_2 = element[elementID][3];}
		else
			OutputError("RawMesh::GetElementEdge", "Wrong edgeID", edgeID);
	}
	else if (elementProperty.elementType == TRIANGLE)
	{
		if (edgeID == 0)
			{point_1 = element[elementID][0]; point_2 = element[elementID][1];}
		else if (edgeID == 1)
			{point_1 = element[elementID][1]; point_2 = element[elementID][2];}
		else if (edgeID == 2)
			{point_1 = element[elementID][2]; point_2 = element[elementID][0];}
		else
			OutputError("RawMesh::GetElementEdge", "Wrong edgeID", edgeID);
	}
	else if (elementProperty.elementType == HEXAGON)
	{
		if (edgeID == 0)
			{point_1 = element[elementID][0]; point_2 = element[elementID][1];}
		else if (edgeID == 1)
			{point_1 = element[elementID][1]; point_2 = element[elementID][2];}
		else if (edgeID == 2)
			{point_1 = element[elementID][2]; point_2 = element[elementID][3];}
		else if (edgeID == 3)
			{point_1 = element[elementID][3]; point_2 = element[elementID][4];}
		else if (edgeID == 4)
			{point_1 = element[elementID][4]; point_2 = element[elementID][5];}
		else if (edgeID == 5)
			{point_1 = element[elementID][5]; point_2 = element[elementID][0];}
		else
			OutputError("RawMesh::GetElementEdge", "Wrong edgeID", edgeID);

		if (point_1 >= 0 && point_2 == -1)
			point_2 = element[elementID][0];
		else if (point_1 == -1 && point_2 >= 0)
			point_2 = -1;
	}
	else if (elementProperty.elementType == LINE)
	{
		if (edgeID == 0)
			{point_1 = element[elementID][0]; point_2 = element[elementID][1];}
		else
			OutputError("RawMesh::GetElementEdge", "Wrong edgeID", edgeID);
	}
	else
	{
		point_1 = -1; point_2 = -1;
		OutputError("RawMesh::GetElementEdge", "No correct element type");
	}

	return;
}

void RawMesh::GetElementFace(int elementID, int faceID, int *vertexIndex)
{
	assert(elementID < elementNumber);

	if (elementProperty.elementType == HEXAHEDRON)
	{
		if (faceID == 0)
			{vertexIndex[0] = element[elementID][0]; vertexIndex[1] = element[elementID][3]; vertexIndex[2] = element[elementID][2]; vertexIndex[3] = element[elementID][1];}
		else if (faceID == 1)
			{vertexIndex[0] = element[elementID][0]; vertexIndex[1] = element[elementID][1]; vertexIndex[2] = element[elementID][5]; vertexIndex[3] = element[elementID][4];}
		else if (faceID == 2)
			{vertexIndex[0] = element[elementID][1]; vertexIndex[1] = element[elementID][2]; vertexIndex[2] = element[elementID][6]; vertexIndex[3] = element[elementID][5];}
		else if (faceID == 3)
			{vertexIndex[0] = element[elementID][2]; vertexIndex[1] = element[elementID][3]; vertexIndex[2] = element[elementID][7]; vertexIndex[3] = element[elementID][6];}
		else if (faceID == 4)
			{vertexIndex[0] = element[elementID][3]; vertexIndex[1] = element[elementID][0]; vertexIndex[2] = element[elementID][4]; vertexIndex[3] = element[elementID][7];}
		else if (faceID == 5)
			{vertexIndex[0] = element[elementID][4]; vertexIndex[1] = element[elementID][5]; vertexIndex[2] = element[elementID][6]; vertexIndex[3] = element[elementID][7];}
		else
			OutputError("RawMesh::GetElementFace", "Wrong face", faceID);
	}
	else if (elementProperty.elementType == TETRAHEDRON)
	{
		if (faceID == 0)
			{vertexIndex[0] = element[elementID][0]; vertexIndex[1] = element[elementID][2]; vertexIndex[2] = element[elementID][1];}
		else if (faceID == 1)
			{vertexIndex[0] = element[elementID][0]; vertexIndex[1] = element[elementID][1]; vertexIndex[2] = element[elementID][3];}
		else if (faceID == 2)
			{vertexIndex[0] = element[elementID][0]; vertexIndex[1] = element[elementID][3]; vertexIndex[2] = element[elementID][2];}
		else if (faceID == 3)
			{vertexIndex[0] = element[elementID][1]; vertexIndex[1] = element[elementID][2]; vertexIndex[2] = element[elementID][3];}
		else
			OutputError("RawMesh::GetElementFace", "Wrong face", faceID);
	}
	else if (elementProperty.elementType == QUADRILATERAL)
	{
		if (faceID == 0)
			{vertexIndex[0] = element[elementID][0]; vertexIndex[1] = element[elementID][1]; vertexIndex[2] = element[elementID][2]; vertexIndex[3] = element[elementID][3];}
		else
			OutputError("RawMesh::GetElementFace", "Wrong face", faceID);
	}
	else if (elementProperty.elementType == TRIANGLE)
	{
		if (faceID == 0)
			{vertexIndex[0] = element[elementID][0]; vertexIndex[1] = element[elementID][1]; vertexIndex[2] = element[elementID][2];}
		else
			OutputError("RawMesh::GetElementFace", "Wrong face", faceID);
	}
	else
		OutputError("RawMesh::GetElementFace", " No such face in the element");

	return;
}

void RawMesh::GetNewIndexOrder(int elementID, int vertexID, int *newInd)
{
	if (elementProperty.elementType == QUADRILATERAL)
	{
		if (element[elementID][0] == vertexID)
		{
			newInd[0] = element[elementID][0];
			newInd[1] = element[elementID][1];
			newInd[2] = element[elementID][2];
			newInd[3] = element[elementID][3];
		}
		else if (element[elementID][1] == vertexID)
		{
			newInd[0] = element[elementID][1];
			newInd[1] = element[elementID][2];
			newInd[2] = element[elementID][3];
			newInd[3] = element[elementID][0];
		}
		else if (element[elementID][2] == vertexID)
		{
			newInd[0] = element[elementID][2];
			newInd[1] = element[elementID][3];
			newInd[2] = element[elementID][0];
			newInd[3] = element[elementID][1];
		}
		else if (element[elementID][3] == vertexID)
		{
			newInd[0] = element[elementID][3];
			newInd[1] = element[elementID][0];
			newInd[2] = element[elementID][1];
			newInd[3] = element[elementID][2];
		}
		else
		{
			OutputError("RawMesh::GetNewIndexOrder", "No such vertex ID %d in the element. ID", elementID);
		}
	}
	else if (elementProperty.elementType == HEXAHEDRON)
	{
		if (element[elementID][0] == vertexID)
		{
			newInd[0] = element[elementID][0];
			newInd[1] = element[elementID][1];
			newInd[2] = element[elementID][2];
			newInd[3] = element[elementID][3];
			newInd[4] = element[elementID][4];
			newInd[5] = element[elementID][5];
			newInd[6] = element[elementID][6];
			newInd[7] = element[elementID][7];
		}
		else if (element[elementID][1] == vertexID)
		{
			newInd[0] = element[elementID][1];
			newInd[1] = element[elementID][2];
			newInd[2] = element[elementID][3];
			newInd[3] = element[elementID][0];
			newInd[4] = element[elementID][5];
			newInd[5] = element[elementID][6];
			newInd[6] = element[elementID][7];
			newInd[7] = element[elementID][4];
		}
		else if (element[elementID][2] == vertexID)
		{
			newInd[0] = element[elementID][2];
			newInd[1] = element[elementID][3];
			newInd[2] = element[elementID][0];
			newInd[3] = element[elementID][1];
			newInd[4] = element[elementID][6];
			newInd[5] = element[elementID][7];
			newInd[6] = element[elementID][4];
			newInd[7] = element[elementID][5];
		}
		else if (element[elementID][3] == vertexID)
		{
			newInd[0] = element[elementID][3];
			newInd[1] = element[elementID][0];
			newInd[2] = element[elementID][1];
			newInd[3] = element[elementID][2];
			newInd[4] = element[elementID][7];
			newInd[5] = element[elementID][4];
			newInd[6] = element[elementID][5];
			newInd[7] = element[elementID][6];
		}
		else if (element[elementID][4] == vertexID)
		{
			newInd[0] = element[elementID][4];
			newInd[1] = element[elementID][5];
			newInd[2] = element[elementID][1];
			newInd[3] = element[elementID][0];
			newInd[4] = element[elementID][7];
			newInd[5] = element[elementID][6];
			newInd[6] = element[elementID][2];
			newInd[7] = element[elementID][3];
		}
		else if (element[elementID][5] == vertexID)
		{
			newInd[0] = element[elementID][5];
			newInd[1] = element[elementID][6];
			newInd[2] = element[elementID][2];
			newInd[3] = element[elementID][1];
			newInd[4] = element[elementID][4];
			newInd[5] = element[elementID][7];
			newInd[6] = element[elementID][3];
			newInd[7] = element[elementID][0];

		}
		else if (element[elementID][6] == vertexID)
		{
			newInd[0] = element[elementID][6];
			newInd[1] = element[elementID][7];
			newInd[2] = element[elementID][3];
			newInd[3] = element[elementID][2];
			newInd[4] = element[elementID][5];
			newInd[5] = element[elementID][4];
			newInd[6] = element[elementID][0];
			newInd[7] = element[elementID][1];
		}
		else if (element[elementID][7] == vertexID)
		{
			newInd[0] = element[elementID][7];
			newInd[1] = element[elementID][4];
			newInd[2] = element[elementID][0];
			newInd[3] = element[elementID][3];
			newInd[4] = element[elementID][6];
			newInd[5] = element[elementID][5];
			newInd[6] = element[elementID][1];
			newInd[7] = element[elementID][2];
		}
		else
		{
			OutputError("RawMesh::GetNewIndexOrder", "No such vertex ID %d in the element. ID", elementID);
		}
	}

	return;
}

void RawMesh::GetNewIndexOrder(int elementID, int faceIndex, int vertexID, int *newInd)
{
	if (elementProperty.elementType == QUADRILATERAL)
	{
		if (element[elementID][0] == vertexID)
		{
			newInd[0] = element[elementID][0];
			newInd[1] = element[elementID][1];
			newInd[2] = element[elementID][2];
			newInd[3] = element[elementID][3];
		}
		else if (element[elementID][1] == vertexID)
		{
			newInd[0] = element[elementID][1];
			newInd[1] = element[elementID][2];
			newInd[2] = element[elementID][3];
			newInd[3] = element[elementID][0];
		}
		else if (element[elementID][2] == vertexID)
		{
			newInd[0] = element[elementID][2];
			newInd[1] = element[elementID][3];
			newInd[2] = element[elementID][0];
			newInd[3] = element[elementID][1];
		}
		else if (element[elementID][3] == vertexID)
		{
			newInd[0] = element[elementID][3];
			newInd[1] = element[elementID][0];
			newInd[2] = element[elementID][1];
			newInd[3] = element[elementID][2];
		}
		else
		{
			OutputError("RawMesh::GetNewIndexOrder", "No such vertex ID %d in the element. ID", elementID);
		}
	}
	else if (elementProperty.elementType == HEXAHEDRON)
	{
		//face 0
		if (faceIndex == 0)
		{
			if (element[elementID][0] == vertexID)
			{
				newInd[0] = element[elementID][0];
				newInd[1] = element[elementID][1];
				newInd[2] = element[elementID][2];
				newInd[3] = element[elementID][3];
				newInd[4] = element[elementID][4];
				newInd[5] = element[elementID][5];
				newInd[6] = element[elementID][6];
				newInd[7] = element[elementID][7];
			}
			else if (element[elementID][1] == vertexID)
			{
				newInd[0] = element[elementID][1];
				newInd[1] = element[elementID][2];
				newInd[2] = element[elementID][3];
				newInd[3] = element[elementID][0];
				newInd[4] = element[elementID][5];
				newInd[5] = element[elementID][6];
				newInd[6] = element[elementID][7];
				newInd[7] = element[elementID][4];
			}
			else if (element[elementID][2] == vertexID)
			{
				newInd[0] = element[elementID][2];
				newInd[1] = element[elementID][3];
				newInd[2] = element[elementID][0];
				newInd[3] = element[elementID][1];
				newInd[4] = element[elementID][6];
				newInd[5] = element[elementID][7];
				newInd[6] = element[elementID][4];
				newInd[7] = element[elementID][5];
			}
			else if (element[elementID][3] == vertexID)
			{
				newInd[0] = element[elementID][3];
				newInd[1] = element[elementID][0];
				newInd[2] = element[elementID][1];
				newInd[3] = element[elementID][2];
				newInd[4] = element[elementID][7];
				newInd[5] = element[elementID][4];
				newInd[6] = element[elementID][5];
				newInd[7] = element[elementID][6];
			}
			else
			{
				OutputError("RawMesh::GetNewIndexOrder", "No such vertex ID %d in the element. ID", elementID);
			}
		}
		//face 1
		else if (faceIndex == 1)
		{
			if (element[elementID][0] == vertexID)
			{
				newInd[0] = element[elementID][0];
				newInd[1] = element[elementID][4];
				newInd[2] = element[elementID][5];
				newInd[3] = element[elementID][1];
				newInd[4] = element[elementID][3];
				newInd[5] = element[elementID][7];
				newInd[6] = element[elementID][6];
				newInd[7] = element[elementID][2];
			}
			else if (element[elementID][4] == vertexID)
			{
				newInd[0] = element[elementID][4];
				newInd[1] = element[elementID][5];
				newInd[2] = element[elementID][1];
				newInd[3] = element[elementID][0];
				newInd[4] = element[elementID][7];
				newInd[5] = element[elementID][6];
				newInd[6] = element[elementID][2];
				newInd[7] = element[elementID][3];
			}
			else if (element[elementID][5] == vertexID)
			{
				newInd[0] = element[elementID][5];
				newInd[1] = element[elementID][1];
				newInd[2] = element[elementID][0];
				newInd[3] = element[elementID][4];
				newInd[4] = element[elementID][6];
				newInd[5] = element[elementID][2];
				newInd[6] = element[elementID][3];
				newInd[7] = element[elementID][7];
			}
			else if (element[elementID][1] == vertexID)
			{
				newInd[0] = element[elementID][1];
				newInd[1] = element[elementID][0];
				newInd[2] = element[elementID][4];
				newInd[3] = element[elementID][5];
				newInd[4] = element[elementID][2];
				newInd[5] = element[elementID][3];
				newInd[6] = element[elementID][7];
				newInd[7] = element[elementID][6];
			}
			else
			{
				OutputError("RawMesh::GetNewIndexOrder", "No such vertex ID %d in the element. ID", elementID);
			}
		}
		//face 2
		else if (faceIndex == 2)
		{
			if (element[elementID][1] == vertexID)
			{
				newInd[0] = element[elementID][1];
				newInd[1] = element[elementID][5];
				newInd[2] = element[elementID][6];
				newInd[3] = element[elementID][2];
				newInd[4] = element[elementID][0];
				newInd[5] = element[elementID][4];
				newInd[6] = element[elementID][7];
				newInd[7] = element[elementID][3];
			}
			else if (element[elementID][5] == vertexID)
			{
				newInd[0] = element[elementID][5];
				newInd[1] = element[elementID][6];
				newInd[2] = element[elementID][2];
				newInd[3] = element[elementID][1];
				newInd[4] = element[elementID][4];
				newInd[5] = element[elementID][7];
				newInd[6] = element[elementID][3];
				newInd[7] = element[elementID][0];
			}
			else if (element[elementID][6] == vertexID)
			{
				newInd[0] = element[elementID][6];
				newInd[1] = element[elementID][2];
				newInd[2] = element[elementID][1];
				newInd[3] = element[elementID][5];
				newInd[4] = element[elementID][7];
				newInd[5] = element[elementID][3];
				newInd[6] = element[elementID][0];
				newInd[7] = element[elementID][4];
			}
			else if (element[elementID][2] == vertexID)
			{
				newInd[0] = element[elementID][2];
				newInd[1] = element[elementID][1];
				newInd[2] = element[elementID][5];
				newInd[3] = element[elementID][6];
				newInd[4] = element[elementID][3];
				newInd[5] = element[elementID][0];
				newInd[6] = element[elementID][4];
				newInd[7] = element[elementID][7];
			}
			else
			{
				OutputError("RawMesh::GetNewIndexOrder", "No such vertex ID %d in the element. ID", elementID);
			}
		}
		//face 3
		else if (faceIndex == 3)
		{
			if (element[elementID][2] == vertexID)
			{
				newInd[0] = element[elementID][2];
				newInd[1] = element[elementID][6];
				newInd[2] = element[elementID][7];
				newInd[3] = element[elementID][3];
				newInd[4] = element[elementID][1];
				newInd[5] = element[elementID][5];
				newInd[6] = element[elementID][4];
				newInd[7] = element[elementID][0];
			}
			else if (element[elementID][6] == vertexID)
			{
				newInd[0] = element[elementID][6];
				newInd[1] = element[elementID][7];
				newInd[2] = element[elementID][3];
				newInd[3] = element[elementID][2];
				newInd[4] = element[elementID][5];
				newInd[5] = element[elementID][4];
				newInd[6] = element[elementID][0];
				newInd[7] = element[elementID][1];
			}
			else if (element[elementID][7] == vertexID)
			{
				newInd[0] = element[elementID][7];
				newInd[1] = element[elementID][3];
				newInd[2] = element[elementID][2];
				newInd[3] = element[elementID][6];
				newInd[4] = element[elementID][4];
				newInd[5] = element[elementID][0];
				newInd[6] = element[elementID][1];
				newInd[7] = element[elementID][5];
			}
			else if (element[elementID][3] == vertexID)
			{
				newInd[0] = element[elementID][3];
				newInd[1] = element[elementID][2];
				newInd[2] = element[elementID][6];
				newInd[3] = element[elementID][7];
				newInd[4] = element[elementID][0];
				newInd[5] = element[elementID][1];
				newInd[6] = element[elementID][5];
				newInd[7] = element[elementID][4];
			}
			else
			{
				OutputError("RawMesh::GetNewIndexOrder", "No such vertex ID %d in the element. ID", elementID);
			}
		}
		//face 4
		else if (faceIndex == 4)
		{
			if (element[elementID][0] == vertexID)
			{
				newInd[0] = element[elementID][0];
				newInd[1] = element[elementID][3];
				newInd[2] = element[elementID][7];
				newInd[3] = element[elementID][4];
				newInd[4] = element[elementID][1];
				newInd[5] = element[elementID][2];
				newInd[6] = element[elementID][6];
				newInd[7] = element[elementID][5];
			}
			else if (element[elementID][3] == vertexID)
			{
				newInd[0] = element[elementID][3];
				newInd[1] = element[elementID][7];
				newInd[2] = element[elementID][4];
				newInd[3] = element[elementID][0];
				newInd[4] = element[elementID][2];
				newInd[5] = element[elementID][6];
				newInd[6] = element[elementID][5];
				newInd[7] = element[elementID][1];
			}
			else if (element[elementID][7] == vertexID)
			{
				newInd[0] = element[elementID][7];
				newInd[1] = element[elementID][4];
				newInd[2] = element[elementID][0];
				newInd[3] = element[elementID][3];
				newInd[4] = element[elementID][6];
				newInd[5] = element[elementID][5];
				newInd[6] = element[elementID][1];
				newInd[7] = element[elementID][2];
			}
			else if (element[elementID][4] == vertexID)
			{
				newInd[0] = element[elementID][4];
				newInd[1] = element[elementID][0];
				newInd[2] = element[elementID][3];
				newInd[3] = element[elementID][7];
				newInd[4] = element[elementID][5];
				newInd[5] = element[elementID][1];
				newInd[6] = element[elementID][2];
				newInd[7] = element[elementID][6];
			}
			else
			{
				OutputError("RawMesh::GetNewIndexOrder", "No such vertex ID %d in the element. ID", elementID);
			}
		}
		//face 5
		else if (faceIndex == 5)
		{
			if (element[elementID][4] == vertexID)
			{
				newInd[0] = element[elementID][4];
				newInd[1] = element[elementID][7];
				newInd[2] = element[elementID][6];
				newInd[3] = element[elementID][5];
				newInd[4] = element[elementID][0];
				newInd[5] = element[elementID][3];
				newInd[6] = element[elementID][2];
				newInd[7] = element[elementID][1];
			}
			else if (element[elementID][5] == vertexID)
			{
				newInd[0] = element[elementID][5];
				newInd[1] = element[elementID][4];
				newInd[2] = element[elementID][7];
				newInd[3] = element[elementID][6];
				newInd[4] = element[elementID][1];
				newInd[5] = element[elementID][0];
				newInd[6] = element[elementID][3];
				newInd[7] = element[elementID][2];
			}
			else if (element[elementID][6] == vertexID)
			{
				newInd[0] = element[elementID][6];
				newInd[1] = element[elementID][5];
				newInd[2] = element[elementID][4];
				newInd[3] = element[elementID][7];
				newInd[4] = element[elementID][2];
				newInd[5] = element[elementID][1];
				newInd[6] = element[elementID][0];
				newInd[7] = element[elementID][3];
			}
			else if (element[elementID][7] == vertexID)
			{
				newInd[0] = element[elementID][7];
				newInd[1] = element[elementID][6];
				newInd[2] = element[elementID][5];
				newInd[3] = element[elementID][4];
				newInd[4] = element[elementID][3];
				newInd[5] = element[elementID][2];
				newInd[6] = element[elementID][1];
				newInd[7] = element[elementID][0];
			}
			else
			{
				OutputError("RawMesh::GetNewIndexOrder", "No such vertex ID %d in the element. ID", elementID);
			}
		}
		else
		{
			OutputError("RawMesh::GetNewIndexOrder", "No such vertex ID %d in the element. ID", elementID);
		}
	}

	return;
}

int  RawMesh::GetNeighboringElementByFace(int elementID, int face[])
{
	if (elementValence == NULL)
		GetElementValence();

	int i, j, k, iVert, iVert_2, iElem=-1;
	bool isTheSameElement=true;

	iVert = face[0];
	for (i=0; i<elementValenceNumber[iVert]; ++i)
	{
		iElem = elementValence[iVert][i];
		if (iElem == elementID)
			continue;

		isTheSameElement = true;
		for (j=1; j<elementProperty.faceVertexNumber; ++j)
		{
			iVert_2 = face[j];
			for (k=0; k<elementValenceNumber[iVert_2]; ++k)
			{
				if (elementValence[iVert_2][k] == iElem)
					break;
			}

			if (k >= elementValenceNumber[iVert_2])
			{
				isTheSameElement = false;
				break;
			}
		}

		if (isTheSameElement)
			break;
	}
	if (iElem == elementID)
		isTheSameElement = false;

	if (isTheSameElement)
		return iElem;
	return -1;
}

void RawMesh::GetMeshInfo()
{
	int i, j;

	for (i=0; i<3; ++i)
	{
		meshInfo.minCorner[i] = MAX;
		meshInfo.maxCorner[i] = MIN;
	}

	for (i=0; i<3; ++i)
	{
		for (j=0; j<vertexNumber; j++)
		{
			meshInfo.minCorner[i] = meshInfo.minCorner[i] > vertex[j][i] ? vertex[j][i] : meshInfo.minCorner[i];
			meshInfo.maxCorner[i] = meshInfo.maxCorner[i] < vertex[j][i] ? vertex[j][i] : meshInfo.maxCorner[i];
		}
	}

	for (i=0; i<3; ++i)
	{
		meshInfo.center[i] = (meshInfo.maxCorner[i] + meshInfo.minCorner[i])/2;
	}

    meshInfo.biggestSize = (meshInfo.maxCorner[0] - meshInfo.minCorner[0] > meshInfo.maxCorner[1] - meshInfo.minCorner[1] ?
							meshInfo.maxCorner[0] - meshInfo.minCorner[0] : meshInfo.maxCorner[1] - meshInfo.minCorner[1]);
    meshInfo.biggestSize = (meshInfo.maxCorner[2] - meshInfo.minCorner[2] > meshInfo.biggestSize ?
							meshInfo.maxCorner[2] - meshInfo.minCorner[2] : meshInfo.biggestSize);

	meshInfo.smallestSize = (meshInfo.maxCorner[0] - meshInfo.minCorner[0] < meshInfo.maxCorner[1] - meshInfo.minCorner[1] ?
							meshInfo.maxCorner[0] - meshInfo.minCorner[0] : meshInfo.maxCorner[1] - meshInfo.minCorner[1]);
    meshInfo.smallestSize = (meshInfo.maxCorner[2] - meshInfo.minCorner[2] < meshInfo.smallestSize ?
							meshInfo.maxCorner[2] - meshInfo.minCorner[2] : meshInfo.smallestSize);

	return;
}

//normalType: 1 - vertex normal; 2 - element normal.
bool RawMesh::ComputeNormal(int normalType)
{
	if (normalType == 1)
		return ComputeVertexNormal();
	else if (normalType == 2)
		return ComputeElementNormal();

	return false;
}


bool RawMesh::ComputeVertexNormal()
{
	if (vertexNormal != NULL)
		FreeMatrix(vertexNormal);

	if (elementProperty.elementType == TRIANGLE)
		return ComputeVertexNormal_Tri();
	else if (elementProperty.elementType == QUADRILATERAL)
		return ComputeVertexNormal_Quad();
	else if (elementProperty.elementType == POINT)
		return ComputeVertexNormal_Point();
	else if (elementProperty.elementType == LINE)
		return ComputeVertexNormal_Line();
	else if (elementProperty.elementType == TETRAHEDRON)
		return ComputeVertexNormal_Tet();
	else if (elementProperty.elementType == HEXAHEDRON)
		return ComputeVertexNormal_Hex();

	return false;
}

bool RawMesh::ComputeVertexNormal_Point()
{
	int i;

	if(IsEmpty())
	{
		OutputError("RawMesh::ComputeVertexNormal_Point", "Mesh is empty");
		return false;
	}

	if(!InitiateMatrix(vertexNormal, vertexNumber))
	{
		OutputError("RawMesh::ComputeVertexNormal_Point", "Failed to initiate the normal");
		return false;
	}

	if (!GetElementValence())
	{
		OutputError("RawMesh::ComputeVertexNormal_Point", "Failed to get Element Valence");
		return false;
	}

	// calculate the normal for each vertex
	for(i = 0; i < vertexNumber; i++)
	{
		vertexNormal[i][0] = 0.0;
		vertexNormal[i][1] = 0.0;
		vertexNormal[i][2] = 0.0;
	}

	return true;
}

bool RawMesh::ComputeVertexNormal_Line()
{
	int i, j, k;
	int v0, v1, v2;
	double mv0[3], mv1[3];
	double norm[3];

	if(IsEmpty())
	{
		OutputError("RawMesh::ComputeVertexNormal_Line", "Mesh is empty");
		return false;
	}

	if(!InitiateMatrix(vertexNormal, vertexNumber))
	{
		OutputError("RawMesh::ComputeVertexNormal_Line", "Failed to initiate the normal");
		return false;
	}

	if (!GetEdgeValence())
	{
		OutputError("RawMesh::ComputeVertexNormal_Line", "Failed to get Element Valence");
		return false;
	}

	// calculate the normal for each vertex
	mv0[0] = 0.0; mv0[1] = 0.0; mv0[2] = 1.0;
	for(i = 0; i < vertexNumber; i++)
	{
		vertexNormal[i][0] = 0.0;
		vertexNormal[i][1] = 0.0;
		vertexNormal[i][2] = 0.0;

		for(j = 0; j < edgeValenceNumber[i]; j++)
		{
			k = edgeValence[i][j];
			mv1[0] = vertex[k][0] - vertex[i][0];
			mv1[1] = vertex[k][1] - vertex[i][1];
			mv1[2] = vertex[k][2] - vertex[i][2];

			CrossProduct(mv0, mv1, norm);
			Normalize(norm);
			for(k = 0; k < 3; k++)
				vertexNormal[i][k] += norm[k];
		}

		Normalize(vertexNormal[i]);
	}

	return true;
}

bool RawMesh::ComputeVertexNormal_Tri()
{
	int i, j, k, itri;
	int v0, v1, v2;
	double mv0[3], mv1[3], mv2[3];
	double norm[3];

	if(IsEmpty())
	{
		OutputError("RawMesh::ComputeVertexNormal_Tri", "Mesh is empty");
		return false;
	}

	if(!InitiateMatrix(vertexNormal, vertexNumber))
	{
		OutputError("RawMesh::ComputeVertexNormal_Tri", "Failed to initiate the normal");
		return false;
	}

	if (!GetElementValence())
	{
		OutputError("RawMesh::ComputeVertexNormal_Tri", "Failed to get Element Valence");
		return false;
	}

	// calculate the normal for each vertex
	for(i = 0; i < vertexNumber; i++)
	{
		vertexNormal[i][0] = 0.0;
		vertexNormal[i][1] = 0.0;
		vertexNormal[i][2] = 0.0;

		for(j = 0; j < elementValenceNumber[i]; j++)
		{
			itri = elementValence[i][j];
			v0 = element[itri][0];	v1 = element[itri][1];
			v2 = element[itri][2];	
			for(k = 0; k < 3; k++)
			{
				mv0[k] = vertex[v0][k];	mv1[k] = vertex[v1][k];
				mv2[k] = vertex[v2][k];	
			}

			if(i == v0)
				CrossProduct(mv0, mv1, mv2, norm);
			if(i == v1)
				CrossProduct(mv1, mv2, mv0, norm);
			if(i == v2)
				CrossProduct(mv2, mv0, mv1, norm);

			Normalize(norm);
			for(k = 0; k < 3; k++)
				vertexNormal[i][k] += norm[k];
		}

		Normalize(vertexNormal[i]);
	}

	return true;
}

bool RawMesh::ComputeVertexNormal_Quad()
{
	int i, j, k, iquad;
	int v0, v1, v2, v3;
	double mv0[3], mv1[3], mv2[3], mv3[3];
	double norm[3];

	if(IsEmpty())
	{
		OutputError("RawMesh::ComputeVertexNormal_Quad", "Mesh is empty");
		return false;
	}

	if(!InitiateMatrix(vertexNormal, vertexNumber))
	{
		OutputError("RawMesh::ComputeVertexNormal_Quad", "Failed to initiate the normal");
		return false;
	}

	if (!GetElementValence())
	{
		OutputError("RawMesh::ComputeVertexNormal_Quad", "Failed to get Element Valence");
		return false;
	}

	// calculate the normal for each vertex
	for(i = 0; i < vertexNumber; i++)
	{
		vertexNormal[i][0] = 0.0;
		vertexNormal[i][1] = 0.0;
		vertexNormal[i][2] = 0.0;

		for(j = 0; j < elementValenceNumber[i]; j++)
		{
			iquad = elementValence[i][j];
			v0 = element[iquad][0];	v1 = element[iquad][1];
			v2 = element[iquad][2];	v3 = element[iquad][3];
			for(k = 0; k < 3; k++)
			{
				mv0[k] = vertex[v0][k];	mv1[k] = vertex[v1][k];
				mv2[k] = vertex[v2][k];	mv3[k] = vertex[v3][k];
			}

			if(i == v0)
				CrossProduct(mv0, mv1, mv3, norm);
			else if(i == v1)
				CrossProduct(mv1, mv2, mv0, norm);
			else if(i == v2)
				CrossProduct(mv2, mv3, mv1, norm);
			else if(i == v3)
				CrossProduct(mv3, mv0, mv2, norm);

			Normalize(norm);
			for(k = 0; k < 3; k++)
				vertexNormal[i][k] += norm[k];
		}

		Normalize(vertexNormal[i]);
	}

	return true;
}

bool RawMesh::ComputeVertexNormal_Tet()
{
	return false;
}

bool RawMesh::ComputeVertexNormal_Hex()
{
	return false;
}

bool RawMesh::ComputeElementNormal()
{
	if (elementNormal != NULL)
		FreeMatrix(elementNormal);

	if (elementProperty.elementType == TRIANGLE)
		return ComputeElementNormal_Tri();
	else if (elementProperty.elementType == QUADRILATERAL)
		return ComputeElementNormal_Quad();
	else if (elementProperty.elementType == POINT)
		return ComputeElementNormal_Point();
	else if (elementProperty.elementType == LINE)
		return ComputeElementNormal_Line();

	return false;
}

bool RawMesh::ComputeElementNormal_Point()
{
	if (!InitiateMatrix(elementNormal, elementNumber))
		return false;

	double vector_1[3], vector_2[3];

	for (int i=0; i<elementNumber; ++i)
	{
		elementNormal[i][0] = 0.0;
		elementNormal[i][1] = 0.0;
		elementNormal[i][2] = 0.0;
	}

	return true;
}

bool RawMesh::ComputeElementNormal_Line()
{
	if (!InitiateMatrix(elementNormal, elementNumber))
		return false;

	double vector_1[3], vector_2[3];

	for (int i=0; i<elementNumber; ++i)
	{
		elementNormal[i][0] = 0.0;
		elementNormal[i][1] = 0.0;
		elementNormal[i][2] = 0.0;
	}

	return true;
}

bool RawMesh::ComputeElementNormal_Tri()
{
	if (!InitiateMatrix(elementNormal, elementNumber))
		return false;

	double vector_1[3], vector_2[3];

	for (int i=0; i<elementNumber; ++i)
	{
		vector_1[0] = vertex[element[i][1]][0] - vertex[element[i][0]][0];
		vector_1[1] = vertex[element[i][1]][1] - vertex[element[i][0]][1];
		vector_1[2] = vertex[element[i][1]][2] - vertex[element[i][0]][2];

		vector_2[0] = vertex[element[i][2]][0] - vertex[element[i][0]][0];
		vector_2[1] = vertex[element[i][2]][1] - vertex[element[i][0]][1];
		vector_2[2] = vertex[element[i][2]][2] - vertex[element[i][0]][2];

		CrossProduct(vector_1, vector_2, elementNormal[i]);
	}

	return true;
}

bool RawMesh::ComputeElementNormal_Quad()
{
	if (!InitiateMatrix(elementNormal, elementNumber))
		return false;

	double vector_1[3], vector_2[3];

	for (int i=0; i<elementNumber; ++i)
	{
		vector_1[0] = vertex[element[i][1]][0] - vertex[element[i][0]][0];
		vector_1[1] = vertex[element[i][1]][1] - vertex[element[i][0]][1];
		vector_1[2] = vertex[element[i][1]][2] - vertex[element[i][0]][2];

		vector_2[0] = vertex[element[i][3]][0] - vertex[element[i][0]][0];
		vector_2[1] = vertex[element[i][3]][1] - vertex[element[i][0]][1];
		vector_2[2] = vertex[element[i][3]][2] - vertex[element[i][0]][2];

		CrossProduct(vector_1, vector_2, elementNormal[i]);
	}

	return true;
}

//CurvatureType: 1 - max curv; 2 - min curv; 3 - mean curv; 4 - Guassian curv.
bool RawMesh::ComputeVertexCurvature(int curvatureType, bool isUseExtendQuadric)
{
	if (!InitiateMatrix(vertexCurvature, vertexNumber))
		return false;
	
	if (!ComputeVertexNormal() || !InitiateEdgeValence())
		return false;

	int i, j, jSize;
	Vector3d n1, n2, n3, vX, vXX;
	VectorXd vRight, vCo(5);
	Matrix3d mI, mRotate;
	MatrixXd mLeft, mCo, mLeft_T;
	mI.setIdentity();

	if (!isUseExtendQuadric)
	{
		for (i=0; i<vertexNumber; ++i)
		{
			//Simple quadric fix
			//Obtain the rotation matrix from global coordinate to local coordinate
			n3<<vertexNormal[i][0], vertexNormal[i][1], vertexNormal[i][2];
			n3.normalize();
			vX<<1.0, 0.0, 0.0;
			if (n3.dot(vX) > 0.99)		//Check if n3 equal to vX
				vX<<0.0, 1.0, 0.0;
			n1 = (mI - n3*n3.transpose())*vX;
			n1.normalize();
			n2 = n3.cross(n1);
			n2.normalize();
			mRotate<<n1, n2, n3;

			//Calculate the coefficient matrix
			jSize = edgeValenceNumber[i];
			mLeft.resize(jSize, 3);
			vRight.resize(jSize);
			for (j=0; j<jSize; ++j)
			{
				vX<< vertex[edgeValence[i][j]][0] - vertex[i][0],
						vertex[edgeValence[i][j]][1] - vertex[i][1],
						vertex[edgeValence[i][j]][2] - vertex[i][2];
				vXX = mRotate*vX;

				mLeft(j, 0) = vXX(0)*vXX(0);
				mLeft(j, 1) = vXX(0)*vXX(1);
				mLeft(j, 2) = vXX(1)*vXX(1);
				vRight(j) = vXX(2);
			}
			vX = (mLeft.transpose()*mLeft).inverse()*mLeft.transpose()*vRight;

			//Obtain curvature
			double k1, k2;
			k1 = (vX(0) + vX(2) + sqrt((vX(0) - vX(2))*(vX(0) - vX(2)) + vX(1)*vX(1)))/2;
			k2 = (vX(0) + vX(2) - sqrt((vX(0) - vX(2))*(vX(0) - vX(2)) + vX(1)*vX(1)))/2;
			if		(curvatureType == 1)
				vertexCurvature[i] = k1 > k2 ? k1 : k2;	//fabs(k1) > fabs(k2) ? fabs(k1) : fabs(k2);
			else if (curvatureType == 2)
				vertexCurvature[i] = k2 < k1 ? k2 : k1;	//fabs(k1) < fabs(k2) ? fabs(k1) : fabs(k2);
			else if (curvatureType == 3)
				vertexCurvature[i] = vX(0) + vX(2);
			else if (curvatureType == 4)
				vertexCurvature[i] = 4*vX(0)*vX(2) - vX(1)*vX(1);
			else
				vertexCurvature[i] = 0.0;
		}
	}
	else
	{
		vector <int> vEdgeValence;
		for (i=0; i<vertexNumber; ++i)
		{
			//Extended quadric fix
			n3<<vertexNormal[i][0], vertexNormal[i][1], vertexNormal[i][2];
			n3.normalize();
			int iLoop=0;
			while(++iLoop <= 100)
			{
				//Obtain the rotation matrix from global coordinate to local coordinate
				vX<<1.0, 0.0, 0.0;
				if (n3.dot(vX) > 0.99985)		//Check if n3 equal to vX
					vX<<0.0, 1.0, 0.0;
				n1 = (mI - n3*n3.transpose())*vX;
				n1.normalize();
				n2 = n3.cross(n1);
				n2.normalize();
				mRotate<<n1, n2, n3;

				//Calculate the coefficient matrix
				jSize = edgeValenceNumber[i];
				vEdgeValence.clear();
				for (j=0; j<jSize; ++j)
					vEdgeValence.push_back(edgeValence[i][j]);
				while (jSize < 5)
				{
					//printf("%d less than 5.\n", i);

					//mLeft.resize(jSize, 3);
					//vRight.resize(jSize);
					//for (j=0; j<jSize; ++j)
					//{
					//	vX<< vertex[edgeValence[i][j]][0] - vertex[i][0],
					//			vertex[edgeValence[i][j]][1] - vertex[i][1],
					//			vertex[edgeValence[i][j]][2] - vertex[i][2];
					//	vXX = mRotate*vX;

					//	mLeft(j, 0) = vXX(0)*vXX(0);
					//	mLeft(j, 1) = vXX(0)*vXX(1);
					//	mLeft(j, 2) = vXX(1)*vXX(1);
					//	vRight(j) = vXX(2);
					//}
					//vX = (mLeft.transpose()*mLeft).inverse()*mLeft.transpose()*vRight;
					//vCo<< vX(0), vX(1), vX(2), 1.0, 1.0;
					//break;

					int k, ii, jNode, kNode, kSize;
					for (j=0; j<jSize; ++j)
					{
						jNode = vEdgeValence[j];
						for (k=0; k<edgeValenceNumber[jNode]; ++k)
						{
							kNode = edgeValence[jNode][k];
							kSize = vEdgeValence.size();
							for (ii=0; ii<kSize; ++ii)
							{
								if (kNode == vEdgeValence[ii])
									break;
							}
							if (ii >= jSize)
								vEdgeValence.push_back(kNode);
						}
					}
					jSize = vEdgeValence.size();
				}
				mLeft.resize(jSize, 5);
				mLeft_T.resize(5, jSize);
				mCo.resize(5, 5);
				vRight.resize(jSize);
				for (j=0; j<jSize; ++j)
				{
					vX<< vertex[vEdgeValence[j]][0] - vertex[i][0],
						 vertex[vEdgeValence[j]][1] - vertex[i][1],
						 vertex[vEdgeValence[j]][2] - vertex[i][2];
					vXX = mRotate*vX;

					mLeft(j, 0) = vXX(0)*vXX(0);
					mLeft(j, 1) = vXX(0)*vXX(1);
					mLeft(j, 2) = vXX(1)*vXX(1);
					mLeft(j, 3) = vXX(0);
					mLeft(j, 4) = vXX(1);
					vRight(j) = vXX(2);
				}
				mLeft_T = mLeft.transpose();
				mCo = mLeft_T*mLeft;
				vCo = (mCo.inverse()*mLeft_T)*vRight;
				n2<< -vCo(3), -vCo(4), 1;
				n2.normalize();
				n2 = mRotate.inverse()*n2;
				n2.normalize();

				if (n3.dot(n2) > 0.99985)	//angle < 1 degree
					break;
				else
					n3 = n2;
			}
			vX<< vCo(0), vCo(1), vCo(2);
			//printf("%d ", i);
			//*/

			//Obtain curvature
			double k1, k2;
			k1 = (vX(0) + vX(2) + sqrt((vX(0) - vX(2))*(vX(0) - vX(2)) + vX(1)*vX(1)))/2;
			k2 = (vX(0) + vX(2) - sqrt((vX(0) - vX(2))*(vX(0) - vX(2)) + vX(1)*vX(1)))/2;
			if		(curvatureType == 1)
				vertexCurvature[i] = k1 > k2 ? k1 : k2;	//fabs(k1) > fabs(k2) ? fabs(k1) : fabs(k2);
			else if (curvatureType == 2)
				vertexCurvature[i] = k2 < k1 ? k2 : k1;	//fabs(k1) < fabs(k2) ? fabs(k1) : fabs(k2);
			else if (curvatureType == 3)
				vertexCurvature[i] = vX(0) + vX(2);
			else if (curvatureType == 4)
				vertexCurvature[i] = 4*vX(0)*vX(2) - vX(1)*vX(1);
			else
				vertexCurvature[i] = 0.0;
		}
	}

	return true;
}

void RawMesh::CrossProduct(double vector_1[3], double vector_2[3], double result[3])
{
	result[0] = vector_1[1]*vector_2[2] - vector_1[2]*vector_2[1];
	result[1] = vector_1[2]*vector_2[0] - vector_1[0]*vector_2[2];
	result[2] = vector_1[0]*vector_2[1] - vector_1[1]*vector_2[0];

	return;
}

void RawMesh::CrossProduct(double point_1[3], double point_2[3], double point_3[3], double result[3])
{
	double vector_1[3], vector_2[3];

	vector_1[0] = point_2[0] - point_1[0];
	vector_1[1] = point_2[1] - point_1[1];
	vector_1[2] = point_2[2] - point_1[2];

	vector_2[0] = point_3[0] - point_1[0];
	vector_2[1] = point_3[1] - point_1[1];
	vector_2[2] = point_3[2] - point_1[2];

	result[0] = vector_1[1]*vector_2[2] - vector_1[2]*vector_2[1];
	result[1] = vector_1[2]*vector_2[0] - vector_1[0]*vector_2[2];
	result[2] = vector_1[0]*vector_2[1] - vector_1[1]*vector_2[0];

	return;
}

void RawMesh::DeleteUnusedVertex()
{
	int *newVertexOrder=NULL;
	bool *usedVertex=NULL;
	int i, j, iCount, elementVertexNumber;

	if (!InitiateMatrix(newVertexOrder, vertexNumber, -1) || !InitiateMatrix(usedVertex, vertexNumber, false))
	{
		OutputError("RawMesh::DeleteUnusedVertex", "Fail to initiate newVertexOrder or usedVertex");
		return;
	}

	elementVertexNumber = elementProperty.vertexNumber;
	for (i=0; i<elementNumber; ++i)
		for (j=0; j<elementVertexNumber; ++j)
		{
			if (element[i][j] == -1)
				continue;
			usedVertex[element[i][j]] = true;
		}

	int newVertexNum = 0;
	for (i=0; i<vertexNumber; ++i)
		if (usedVertex[i])
			++newVertexNum;

	if (newVertexNum == vertexNumber)
	{
		FreeMatrix(newVertexOrder);
		FreeMatrix(usedVertex);
		return;
	}

	if (newVertexNum == 0)
	{
		FreeMatrix(vertex);
		FreeMatrix(vertexSign);
		FreeMatrix(vertexColor);
		vertexNumber = 0;
	}
	else
	{
		double (*newVertex)[3]=NULL;
		int *newVertexSign=NULL;
		if (!InitiateMatrix(newVertex, newVertexNum) || !InitiateMatrix(newVertexSign, newVertexNum))
		{
			OutputError("RawMesh::DeleteUnusedVertex", "Fail to initiate newVertexOrder or usedVertex");
			return;
		}

		iCount = 0;
		for (i=0; i<vertexNumber; i++)
		{
			if(usedVertex[i])
			{
				for (j=0; j<3; j++)
					newVertex[iCount][j] = vertex[i][j];
				newVertexSign[iCount] = vertexSign[i];
				newVertexOrder[i] = iCount++;
			}
		}
		for (i=0; i<elementNumber; i++)
			for (j=0; j<elementVertexNumber; j++)
			{
				if (element[i][j] == -1)
					continue;
				element[i][j] = newVertexOrder[element[i][j]];
			}

		//Adding VertexColor
		if (vertexColor != NULL)
		{
			double (*newVertexColor)[3]=NULL;

			if (!InitiateMatrix(newVertexColor, newVertexNum))
			{
				OutputError("RawMesh::DeleteUnusedVertex", "Fail to initiate newVertexOrder or usedVertex");
				return;
			}
			
			iCount = 0;
			for (i=0; i<vertexNumber; i++)
			{
				if(usedVertex[i])
				{
					for (j=0; j<3; j++)
						newVertexColor[iCount][j] = vertexColor[i][j];
					iCount++;
				}
			}

			FreeMatrix(vertexColor);
			vertexColor = newVertexColor;
		}
		//End of For VertexColor

		FreeMatrix(vertex);
		FreeMatrix(vertexSign);
		vertex = newVertex;
		vertexNumber = newVertexNum;
		vertexSign = newVertexSign;
	}

	FreeMatrix(newVertexOrder);
	FreeMatrix(usedVertex);
	UpdateMesh();

	return;
}

void RawMesh::DeleteElement(bool *badElement)
{
	int i, j, iCount;
	int **newElement=NULL, *newElementSign=NULL;

	if (badElement == NULL)
		return;

	iCount = 0;
	for (i=0; i<elementNumber; ++i)
		if (!badElement[i])
			++iCount;

	if (iCount == elementNumber)
		return;
	if (iCount == 0)
	{
		FreeMatrix(element);
		FreeMatrix(elementSign);
		FreeMatrix(elementColor);
		elementNumber = 0;

		UpdateMesh();
		return;
	}

	if (!InitiateMatrix(newElement, iCount, elementProperty.vertexNumber) ||
		!InitiateMatrix(newElementSign, iCount, 0))
	{
		OutputError("RawMesh::DeleteElement", "Fail to initiate newElement");
		return;
	}

	iCount = 0;
	for (i=0; i<elementNumber; ++i)
		if (!badElement[i])
		{
			for (j=0; j<elementProperty.vertexNumber; ++j)
				newElement[iCount][j] = element[i][j];
			newElementSign[iCount] = elementSign[i];
			++iCount;
		}

	//Keep new element color information
	if (elementColor != NULL)
	{
		double (*newColor)[3];
		newColor = NULL;
		InitiateMatrix(newColor, iCount);

		iCount = 0;
		for (i=0; i<elementNumber; ++i)
			if (!badElement[i])
			{
				newColor[iCount][0] = elementColor[i][0];
				newColor[iCount][1] = elementColor[i][1];
				newColor[iCount][2] = elementColor[i][2];
				++iCount;
			}

		FreeMatrix(elementColor);
		elementColor = newColor;
	}
	//End

	FreeMatrix(element);
	FreeMatrix(elementSign);
	element = newElement;
	elementSign = newElementSign;
	elementNumber = iCount;

	UpdateMesh();
	return;
}

void RawMesh::DeleteElement(vector <int> &deletedElement)
{
	int i, j, iCount, iSize;

	iSize = deletedElement.size();
	if (iSize == 0)
		return;

	bool *badElement=NULL;
	if (!InitiateMatrix(badElement, elementNumber, false))
	{
		OutputError("RawMesh::DeleteElement", "Fail to initiate newElement");
		return;
	}

	for (i=0; i<iSize; ++i)
		badElement[deletedElement[i]] = true;

	iSize = 0;
	for (i=0; i<elementNumber; ++i)
		if (badElement[i])
			++iSize;

	if (iSize == elementNumber)
	{
		FreeMatrix(element);
		FreeMatrix(elementSign);
		FreeMatrix(elementColor);
		elementNumber = 0;
		
		UpdateMesh();
		return;
	}

	int **newElement=NULL, *newElementSign=NULL;
	if (!InitiateMatrix(newElement, elementNumber - iSize, elementProperty.vertexNumber) ||
		!InitiateMatrix(newElementSign, elementNumber - iSize, 0))
	{
		OutputError("RawMesh::DeleteElement", "Fail to initiate newElement");
		return;
	}

	iCount = 0;
	for (i=0; i<elementNumber; ++i)
		if (!badElement[i])
		{
			for (j=0; j<elementProperty.vertexNumber; ++j)
				newElement[iCount][j] = element[i][j];
			newElementSign[iCount] = elementSign[i];
			++iCount;
		}

	//Keep new element color information
	if (elementColor != NULL)
	{
		double (*newColor)[3];
		newColor = NULL;
		InitiateMatrix(newColor, elementNumber - iSize);

		iCount = 0;
		for (i=0; i<elementNumber; ++i)
			if (!badElement[i])
			{
				newColor[iCount][0] = elementColor[i][0];
				newColor[iCount][1] = elementColor[i][1];
				newColor[iCount][2] = elementColor[i][2];
				++iCount;
			}

		FreeMatrix(elementColor);
		elementColor = newColor;
	}
	//End

	FreeMatrix(element);
	FreeMatrix(elementSign);
	element = newElement;
	elementSign = newElementSign;
	elementNumber = iCount;

	UpdateMesh();
	return;
}

void RawMesh::DeleteDomainContainTheElement(int iElem)
{
	GetElementValence();
	PushElementSign();

	bool* badElem=NULL;
	InitiateMatrix(badElem, elementNumber, false);

	int i, j, k, iVert;
	bool isQuit=false;

	for (i=0; i<elementNumber; ++i)
		elementSign[i] = 0;
	elementSign[iElem] = 1;
	badElem[iElem] = true;

	int iCount = 0;
	while(!isQuit)	// && ++iCount < 100)
	{
		isQuit = true;
		for (i=0; i<elementNumber; i++)
		{
			if (!(elementSign[i]==1))
				continue;

			isQuit = false;

			for (k=0; k<elementProperty.vertexNumber; ++k)
			{
				iVert = element[i][k];
				for (j=0; j<elementValenceNumber[iVert]; ++j)
				{
					if (elementSign[elementValence[iVert][j]] == 0)
					{
						elementSign[elementValence[iVert][j]] = 1;
						badElem[elementValence[iVert][j]] = true;
					}
				}
			}
			elementSign[i] = 2;
		}
	}

	PopElementSign();
	DeleteElement(badElem);
	FreeMatrix(badElem);	

	return;
}

void RawMesh::DeleteDomainContainTheElement(int iElem, int loopNum)
{
	if (loopNum < 0)
		return;

	GetElementValence();
	PushElementSign();

	bool* badElem=NULL;
	InitiateMatrix(badElem, elementNumber, false);

	int i, j, k, iVert, iSize, ii;
	bool isQuit=false;
	vector <int> loopElem;

	for (i=0; i<elementNumber; ++i)
		elementSign[i] = 0;
	elementSign[iElem] = 1;
	badElem[iElem] = true;
	loopElem.push_back(iElem);

	int iCount = 0;
	while(!isQuit && ++iCount <= loopNum)
	{
		isQuit = true;
		iSize = loopElem.size();
		for (i=0; i<iSize; i++)
		{
			ii = loopElem[i];
			if (!(elementSign[ii]==1))
				continue;

			isQuit = false;

			for (k=0; k<elementProperty.vertexNumber; ++k)
			{
				iVert = element[ii][k];
				for (j=0; j<elementValenceNumber[iVert]; ++j)
				{
					if (elementSign[elementValence[iVert][j]] == 0)
					{
						elementSign[elementValence[iVert][j]] = 1;
						badElem[elementValence[iVert][j]] = true;
						loopElem.push_back(elementValence[iVert][j]);
					}
				}
			}
			elementSign[ii] = 2;
		}
	}

	PopElementSign();
	DeleteElement(badElem);
	FreeMatrix(badElem);

	return;
}

void RawMesh::DeleteDuplicatedPoint(int startPoint)
{
	int *newVertexOrder=NULL, *vertexOrder=NULL;
	bool *usedVertex=NULL;
	int i, j, iCount, p1, p2;
	double ERR, length;

	if (!InitiateMatrix(newVertexOrder, vertexNumber, -1) ||
		!InitiateMatrix(vertexOrder, vertexNumber) ||
		!InitiateMatrix(usedVertex, vertexNumber, true))
	{
		OutputError("RawMesh::DeleteDuplicatedPoint", "Fail to initiate newVertexOrder");
		return;
	}

	for (i=0; i<vertexNumber; ++i)
		vertexOrder[i] = i;

	//Find the ERR
	ERR = MAX;
	for (i=0; i<elementNumber; ++i)
	{
		for (j=0; j<elementProperty.edgeNumber; ++j)
		{
			GetElementEdge(i, j, p1, p2);
			if (p1 == -1 || p2 == -1)
				continue;

			length = (vertex[p1][0] - vertex[p2][0])*(vertex[p1][0] - vertex[p2][0]);
			length += (vertex[p1][1] - vertex[p2][1])*(vertex[p1][1] - vertex[p2][1]);
			length += (vertex[p1][2] - vertex[p2][2])*(vertex[p1][2] - vertex[p2][2]);
			length = sqrt(length);

			if (length <ERR)
				ERR = length;
		}
	}
	ERR = ERR/100;
	if (ERR < 1.0e-6)
		ERR = 1.0e-6;
	//End

	#pragma omp parallel for if (vertexNumber > thresholdForParallel) private(j)
	for (i=0; i<vertexNumber; ++i)
	{
		if (!usedVertex[i])
			continue;

		j = i + 1;
		j = (j < startPoint ? startPoint : j);
		for (j; j<vertexNumber; ++j)
		{
			if (fabs(vertex[j][0]-vertex[i][0]) > ERR)
				continue;
			if (fabs(vertex[j][1]-vertex[i][1]) > ERR)
				continue;
			if (fabs(vertex[j][2]-vertex[i][2]) > ERR)
				continue;

			usedVertex[j] = false;
			vertexOrder[j] = i;
		}
	}

	int newVertexNum = 0;
	for (i=0; i<vertexNumber; ++i)
		if (usedVertex[i])
			++newVertexNum;

	if (newVertexNum == vertexNumber)
	{
		FreeMatrix(vertexOrder);
		FreeMatrix(newVertexOrder);
		FreeMatrix(usedVertex);
		return;
	}

	if (newVertexNum == 0)
	{
		FreeMatrix(vertex);
		FreeMatrix(vertexSign);
		vertexNumber = 0;
	}
	else
	{
		double (*newVertex)[3]=NULL;
		int *newVertexSign=NULL;

		if (!InitiateMatrix(newVertex, newVertexNum) || !InitiateMatrix(newVertexSign, newVertexNum))
		{
			OutputError("RawMesh::DeleteDuplicatedPoint", "Fail to initiate newVertex");
			return;
		}

		iCount = 0;
		for (i=0; i<vertexNumber; ++i)
		{
			if(usedVertex[i])
			{
				for (j=0; j<3; ++j)
					newVertex[iCount][j] = vertex[i][j];
				newVertexSign[iCount] = vertexSign[i];
				newVertexOrder[i] = iCount++;
			}
		}
		for (i=0; i<elementNumber; ++i)
			for (j=0; j<elementProperty.vertexNumber; ++j)
			{
				if (element[i][j] == -1)
					continue;
				element[i][j] = newVertexOrder[vertexOrder[element[i][j]]];
			}

		//Keep color
		if (vertexColor != NULL)
		{
			double (*newVertexColor)[3] = NULL;
			InitiateMatrix(newVertexColor, newVertexNum);

			iCount = 0;
			for (i=0; i<vertexNumber; ++i)
			{
				if(usedVertex[i])
				{
					for (j=0; j<3; ++j)
						newVertexColor[iCount][j] = vertexColor[i][j];
					iCount++;
				}
			}

			FreeMatrix(vertexColor);
			vertexColor = newVertexColor;
		}
		//End of keeping color

		FreeMatrix(vertex);
		FreeMatrix(vertexSign);
		vertex = newVertex;
		vertexNumber = newVertexNum;
		vertexSign = newVertexSign;
	}

	FreeMatrix(vertexOrder);
	FreeMatrix(newVertexOrder);
	FreeMatrix(usedVertex);
	UpdateMesh();
	return;
}

void RawMesh::DeleteDuplicatedElement(int startPoint)
{
	bool *badElem=NULL;
	InitiateMatrix(badElem, elementNumber, false);

	int i, j, k;
	
	#pragma omp parallel for if (elementNumber > thresholdForParallel) private(j, k)
	for (i=0; i<elementNumber; ++i)
	{
		if (badElem[i])
			continue;

		j = i + 1;
		j = (j < startPoint ? startPoint : j);
		for (j; j<elementNumber; ++j)
		{
			for (k=0; k<elementProperty.vertexNumber; ++k)
			{
				if (!IsInElement(element[i][k], j))
					break;
			}

			if (k == elementProperty.vertexNumber)
				badElem[j] = true;
		}
	}

	DeleteElement(badElem);
	DeleteUnusedVertex();
	
	FreeMatrix(badElem);
	return;
}

void RawMesh::Move(double x, double y, double z)
{
	int i;

	if(vertexNumber == 0)
		return;

	for (i=0; i<vertexNumber; ++i)
	{
		vertex[i][0] += x;
		vertex[i][1] += y;
		vertex[i][2] += z;
	}

	return;
}

void RawMesh::Move(double pos[3])
{
	Move(pos[0], pos[1], pos[2]);
	return;
}

void RawMesh::Merge(RawMesh &addedMesh)
{
	AddVertex(addedMesh.vertexNumber, addedMesh.vertex, addedMesh.vertexSign);
	AddElement(addedMesh.elementNumber, addedMesh.element, addedMesh.elementSign);

	int i, j, elementVertexNumber, oldVertexNumber, oldElementNumber;
	elementVertexNumber = addedMesh.elementProperty.vertexNumber;
	oldVertexNumber = vertexNumber - addedMesh.vertexNumber;
	oldElementNumber = elementNumber - addedMesh.elementNumber;

	for (i=oldElementNumber; i<elementNumber; ++i)
		for (j=0; j<elementVertexNumber; ++j)
			element[i][j] += oldVertexNumber;
	
	return;
}

void RawMesh::Normalize(double vector[3])
{
	double sum;
	
	sum = sqrt(vector[0]*vector[0] + vector[1]*vector[1] + vector[2]*vector[2]);
	if (sum == 0.0)
		return;
	vector[0] /= sum; vector[1] /= sum; vector[2] /= sum;

	return;
}

void RawMesh::OutputError(const char *className, const char *message)
{
	fprintf(stderr, "Error. %s: %s.\n", className, message);
}

void RawMesh::OutputError(const char *className, const char *message, int index)
{
	fprintf(stderr, "Error. %s: %s: .\n", className, message, index);
}

void  RawMesh::ReplaceVertexInElement(int elementID, int originalVertexID, int replacedVertexID)
{
	for (int i=0; i<elementProperty.vertexNumber; ++i)
	{
		if (element[elementID][i] == originalVertexID)
		{
			element[elementID][i] = replacedVertexID;
			return;
		}
	}

	return;
}

void RawMesh::PushVertexSign()
{
	copyOfVertexSign.clear();
	for (int i=0; i<vertexNumber; ++i)
		copyOfVertexSign.push_back(vertexSign[i]);

	return;
}

void RawMesh::PushElementSign()
{
	copyOfElementSign.clear();
	for (int i=0; i<elementNumber; ++i)
		copyOfElementSign.push_back(elementSign[i]);

	return;
}

void RawMesh::PopVertexSign()
{
	int i, iSize;
	
	iSize = copyOfVertexSign.size();
	if (iSize != vertexNumber)
	{
		OutputError("RawMesh::PopVertexSign", "Size is not correct");
		
		iSize = (iSize>vertexNumber ? vertexNumber : iSize);
	}

	for (i=0; i<iSize; ++i)
		vertexSign[i] = copyOfVertexSign[i];

	return;
}

void RawMesh::PopElementSign()
{
	int i, iSize;
	
	iSize = copyOfElementSign.size();
	if (iSize != elementNumber)
	{
		OutputError("RawMesh::PopElementSign", "Size is not correct.\n");
		
		iSize = (iSize>elementNumber ? elementNumber : iSize);
	}

	for (i=0; i<iSize; ++i)
		elementSign[i] = copyOfElementSign[i];

	return;
}

void RawMesh::PushElementType()
{
	copyOfElementType.push_back(elementProperty.elementType);
	
	return;
}

void RawMesh::PushAndSetElementType(ElementType iType)
{
	copyOfElementType.push_back(elementProperty.elementType);
	SetElementType(iType);
	
	return;
}

void RawMesh::PopElementType()
{
	ElementType i=copyOfElementType.back();
	copyOfElementType.pop_back();
	SetElementType(i);

	return;
}

RawMesh* RawMesh::ExtractSurface()
{
	int i, j, k, faceVertexNumber;
	int *vertexIndex=NULL;

	faceVertexNumber = GetElementProperty(elementProperty.faceType).vertexNumber;
	InitiateMatrix(vertexIndex,faceVertexNumber);

	vector <int> newFace, newFaceSign, newColor;
	for (i=0; i<elementNumber; ++i)
	{
		//if (elementSign[i] < 0)
		//	continue;

		for (j=0; j<elementProperty.faceNumber; ++j)
		{
			GetElementFace(i, j, vertexIndex);
			if (IsBoundaryFace(vertexIndex))
			{
				for (k=0; k<faceVertexNumber; ++k)
					newFace.push_back(vertexIndex[k]);
				newFaceSign.push_back(elementSign[i]);

				if (elementColor != NULL)
					newColor.push_back(i);
			}
		}
	}

	FreeMatrix(vertexIndex);

	RawMesh *surfaceMesh = new RawMesh;
	surfaceMesh->SetElementType(elementProperty.faceType);
	if (vertexColor == NULL)
		surfaceMesh->AddVertex(vertexNumber, vertex, vertexSign);
	else
		surfaceMesh->AddVertex(vertexNumber, vertex, vertexSign, vertexColor);
	surfaceMesh->AddElement(newFace, newFaceSign);
	
	int iSize = newColor.size();
	if (iSize != 0)
	{
		InitiateMatrix(surfaceMesh->elementColor, surfaceMesh->elementNumber);
		for (i=0; i<iSize; ++i)
		{
			surfaceMesh->elementColor[i][0] = elementColor[newColor[i]][0];
			surfaceMesh->elementColor[i][1] = elementColor[newColor[i]][1];
			surfaceMesh->elementColor[i][2] = elementColor[newColor[i]][2];
		}
	}
	//surfaceMesh->DeleteUnusedVertex();

	return surfaceMesh;
}

RawMesh* RawMesh::ExtractFrame()
{
	int i, j, k, faceVertexNumber;
	int *vertexIndex=NULL;

	faceVertexNumber = GetElementProperty(elementProperty.faceType).vertexNumber;
	InitiateMatrix(vertexIndex,faceVertexNumber);

	vector <int> newFace, newFaceSign, newColor;
	for (i=0; i<elementNumber; ++i)
	{
		//if (elementSign[i] < 0)
		//	continue;

		for (j=0; j<elementProperty.faceNumber; ++j)
		{
			GetElementFace(i, j, vertexIndex);
			for (k=0; k<faceVertexNumber; ++k)
				newFace.push_back(vertexIndex[k]);
			newFaceSign.push_back(elementSign[i]);

			if (elementColor != NULL)
				newColor.push_back(i);
		}
	}

	FreeMatrix(vertexIndex);

	RawMesh *surfaceMesh = new RawMesh;
	surfaceMesh->SetElementType(elementProperty.faceType);
	if (vertexColor == NULL)
		surfaceMesh->AddVertex(vertexNumber, vertex, vertexSign);
	else
		surfaceMesh->AddVertex(vertexNumber, vertex, vertexSign, vertexColor);
	surfaceMesh->AddElement(newFace, newFaceSign);

	int iSize = newColor.size();
	if (iSize != 0)
	{
		InitiateMatrix(surfaceMesh->elementColor, surfaceMesh->elementNumber);
		for (i=0; i<iSize; ++i)
		{
			surfaceMesh->elementColor[i][0] = elementColor[newColor[i]][0];
			surfaceMesh->elementColor[i][1] = elementColor[newColor[i]][1];
			surfaceMesh->elementColor[i][2] = elementColor[newColor[i]][2];
		}
	}
	//surfaceMesh->DeleteUnusedVertex();

	return surfaceMesh;
}

RawMesh* RawMesh::ExtractCrossSection(const char axis, double value[3])
{
	RawMesh* csMesh = new RawMesh;
	(*csMesh) = *this;

	double center[3];
	int i, j, iCount;
	vector <int> dElem;

	for (i=0; i<elementNumber; ++i)
	{
		iCount = 0;
		center[0] = 0.0; center[1] = 0.0; center[2] = 0.0;
		for (j=0; j<elementProperty.vertexNumber; ++j)
		{
			if (element[i][j] == -1)
				continue;

			++iCount;
			center[0] += vertex[element[i][j]][0];
			center[1] += vertex[element[i][j]][1];
			center[2] += vertex[element[i][j]][2];
		}
		center[0] /= iCount;
		center[1] /= iCount;
		center[2] /= iCount;

		if (axis == 'x')
		{
			if (center[0] < value[0])
				dElem.push_back(i);
		}
		else if (axis == 'X')
		{
			if (center[0] >= value[0])
				dElem.push_back(i);
		}
		else if (axis == 'y')
		{
			if (center[1] < value[1])
				dElem.push_back(i);
		}
		else if (axis == 'Y')
		{
			if (center[1] >= value[1])
				dElem.push_back(i);
		}
		else if (axis == 'z')
		{
			if (center[2] < value[2])
				dElem.push_back(i);
		}
		else if (axis == 'Z')
		{
			if (center[2] >= value[2] - 2.0)
				dElem.push_back(i);
		}
	}

	csMesh->DeleteElement(dElem);
	//csMesh->DeleteUnusedVertex();
	dElem.clear();

	return csMesh;
}

void RawMesh::SetBoundaryVertexSign(int sign)
{
	int i, j, k;
	int *vertexIndex=NULL;

	InitiateMatrix(vertexIndex,elementProperty.faceVertexNumber);

	for (i=0; i<vertexNumber; ++i)
		vertexSign[i] = 0;

	for (i=0; i<elementNumber; ++i)
	{
		//if (elementSign[i] < 0)
		//	continue;

		for (j=0; j<elementProperty.faceNumber; ++j)
		{
			GetElementFace(i, j, vertexIndex);
			if (IsBoundaryFace(vertexIndex))
			{
				for (k=0; k<elementProperty.faceVertexNumber; ++k)
					vertexSign[vertexIndex[k]] = sign;
			}
		}
	}

	FreeMatrix(vertexIndex);
	return;
}

void RawMesh::SetInnerVertexSign(int sign)
{
	int i, j, k;
	int *vertexIndex=NULL;
	bool *boundaryVertex=NULL;

	InitiateMatrix(vertexIndex,elementProperty.faceVertexNumber);
	InitiateMatrix(boundaryVertex, vertexNumber, false);

	for (i=0; i<elementNumber; ++i)
	{
		//if (elementSign[i] < 0)
		//	continue;

		for (j=0; j<elementProperty.faceNumber; ++j)
		{
			GetElementFace(i, j, vertexIndex);
			if (IsBoundaryFace(vertexIndex))
			{
				for (k=0; k<elementProperty.faceVertexNumber; ++k)
					boundaryVertex[vertexIndex[k]] = true;
			}
		}
	}

	for (i=0; i<vertexNumber; ++i)
		if (!boundaryVertex[i])
			vertexSign[i] = 0;

	FreeMatrix(vertexIndex);
	FreeMatrix(boundaryVertex);
	return;
}

void RawMesh::SortElementByZ(void)
{
	SortElementByZ(0, elementNumber-1);

	return;
}

void RawMesh::SortElementByZ(int start, int end)
{
	if (start >= end)
		return;

	int left=start-1, right=end, j, tmp;
	while(true)
	{
		while(left < end && vertex[element[++left][0]][2] < vertex[element[end][0]][2]);
		while(right > 0 && vertex[element[--right][0]][2] > vertex[element[end][0]][2]);
		if (left >= right)
			break;

		for (j=0; j<elementProperty.vertexNumber; ++j)
		{
			tmp = element[left][j];
			element[left][j] = element[right][j];
			element[right][j] = tmp;
		}
	}
	for (j=0; j<elementProperty.vertexNumber; ++j)
	{
			tmp = element[left][j];
			element[left][j] = element[end][j];
			element[end][j] = tmp;
	}

	SortElementByZ(start, left -1);
	SortElementByZ(left+1, end);

	return;
}



/////////////////////////////////////////////////////////////////////////////////////////////////
//Modified by Kangkang Hu, 2014

void RawMesh::Smooth(int iLoop)
{
	int i, j, nVert, iCount, ii;
	double center[3];

	PushVertexSign();
	InitiateEdgeValence();

	for (ii=0; ii<iLoop; ++ii)
		for (i=0; i<vertexNumber; ++i)
		{
			if (vertexSign[i] > 0)
				continue;

			center[0] = 0.0; center[1] = 0.0; center[2] = 0.0;
			iCount = 0;
			for (j=0; j<edgeValenceNumber[i]; ++j)
			{
				nVert = edgeValence[i][j];
				//if (vertexSign[nVert] == vertexSign[i])
				//	continue;

				center[0] += vertex[nVert][0];
				center[1] += vertex[nVert][1];
				center[2] += vertex[nVert][2];
				++iCount;
			}

			if (iCount > 0)
			{
				vertex[i][0] += (center[0]/iCount - vertex[i][0])*0.02;
				vertex[i][1] += (center[1]/iCount - vertex[i][1])*0.02;
				vertex[i][2] += (center[2]/iCount - vertex[i][2])*0.02;
			}
		}

		PopVertexSign();
		return;
}

void RawMesh::SmoothSurface(int iLoop)
{
	int i, j, nVert, iCount, ii;
	double center[3];

	PushVertexSign();
	InitiateEdgeValence();

	for (ii=0; ii<iLoop; ++ii)
		for (i=0; i<vertexNumber; ++i)
		{
			if (vertexSign[i] == 0)
				continue;

			center[0] = 0.0; center[1] = 0.0; center[2] = 0.0;
			iCount = 0;
			for (j=0; j<edgeValenceNumber[i]; ++j)
			{
				nVert = edgeValence[i][j];
				if (vertexSign[nVert] == 0)
					continue;

				center[0] += vertex[nVert][0];
				center[1] += vertex[nVert][1];
				center[2] += vertex[nVert][2];
				++iCount;
			}

			if (iCount > 0)
			{
				vertex[i][0] += (center[0]/iCount - vertex[i][0])*0.02;
				vertex[i][1] += (center[1]/iCount - vertex[i][1])*0.02;
				vertex[i][2] += (center[2]/iCount - vertex[i][2])*0.02;
			}
		}

		PopVertexSign();
		return;
}


