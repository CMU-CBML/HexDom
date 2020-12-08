#include "Polycube.h"
#include "StaticVars.h"
#include <sstream>

const double PolycubeConstruction::WEIGHT_LENGTH_EWCVT = 1.0f; 

PolycubeConstruction::PolycubeConstruction(void)
{

}


PolycubeConstruction::PolycubeConstruction(string inputFileName)
{
	ReadKTri(inputFileName.c_str());

	InitiateEdgeValence();
	//////////////////////////////////////////////////////////////////////

	InitiateElementValence();


	//Write("check_tri.vtk");
	//Write("check_tri.raw");
	elementArray.resize(elementNumber);
	InitializeElement(inputFileName.c_str());
	//getchar();


}


PolycubeConstruction::~PolycubeConstruction(void)
{

}
bool PolycubeConstruction::ReadKTri(const char * filename)
{
	FILE *input;
	int i, j, elementType;
	int v[8];

	ifstream infile;
	string oneLine, oneWrod, formerLine;

	infile.open(filename);
	if (infile.is_open())
	{
		while (1)
		{
			getline(infile, oneLine);

			if (oneLine.find("ELEMENT_SHELL") != std::string::npos)
			{
				break;
			}
		}

		elementNumber = 0;
		while (1)
		{
			getline(infile, oneLine);

			if (oneLine.find("NODE") != std::string::npos)
			{
				break;
			}
			elementNumber++;
		}
		vertexNumber = 0;
		while (1)
		{
			getline(infile, oneLine);

			if (oneLine.find("*") != std::string::npos)
			{
				break;
			}
			vertexNumber++;
		}
		cout << "Triangle inforamtion: vertex number " << vertexNumber << " element number " << elementNumber << endl;
		if (!CreateNewMesh(TRIANGLE, vertexNumber, elementNumber)) {
			return false;
		}


		infile.clear();
		infile.seekg(0);

		while (1)
		{
			getline(infile, oneLine);

			if (oneLine.find("ELEMENT_SHELL") != std::string::npos)
			{
				break;
			}
		}

		for (i = 0; i < elementNumber; ++i)
		{

			getline(infile, oneLine);
			std::istringstream iss(oneLine); // string stream
			std::string token;
			size_t pos = -1;
			int location_assign_value = 0;
			while (iss >> token) {
				//cout << token << endl;
				while ((pos = token.rfind(',')) != std::string::npos) {
					token.erase(pos, 1);
				}
				//cout << token << endl;
				//getchar();
				/*std::cout << token << '\n';
				getchar();*/
				if (location_assign_value >= 2 && location_assign_value <= 4)
				{
					std::istringstream(token) >> element[i][location_assign_value - 2];
					
				}
				location_assign_value++;
			}
		}
		for (i = 0; i < elementNumber; ++i)
		{
			element[i][0] = element[i][0] - 1;
			element[i][1] = element[i][1] - 1;
			element[i][2] = element[i][2] - 1;
			/*cout << element[i][0] << " "<<element[i][1] << " " << element[i][2] << endl;
			getchar();*/
		}

		while (1)
		{
			getline(infile, oneLine);

			if (oneLine.find("NODE") != std::string::npos)
			{
				break;
			}

		}

		for (i = 0; i < vertexNumber; ++i)
		{

			getline(infile, oneLine);
			std::istringstream iss(oneLine); // string stream
			std::string token;
			size_t pos = -1;
			int location_assign_value = 0;
			while (iss >> token) {
				while ((pos = token.rfind(',')) != std::string::npos) {
					token.erase(pos, 1);
				}
				/*std::cout << token << '\n';
				getchar();*/
				if (location_assign_value >= 1 && location_assign_value <= 3)
				{
					std::istringstream(token) >> vertex[i][location_assign_value - 1];
					/*cout << vertex[i][location_assign_value - 1] << endl;
					getchar();*/
				}
				location_assign_value++;
			}
		}


		//getline(infile, oneLine);
		//istringstream strStream(oneLine);
		//strStream >> oneWrod >> vertexNumber >> oneWrod;
		////cout << vertexNumber << endl;

		////vector<array<double, 3>> pts;
		////pts.resize(vertexNumber);
		//vector<vector<double> > pts(vertexNumber, vector<double>(3));

		//for (i = 0; i < vertexNumber; ++i)
		//{
		//	infile >> pts[i][0] >> pts[i][1] >> pts[i][2];
		//	//cout << pts[i][0] << endl;
		//}

		//while (1)
		//{
		//	getline(infile, oneLine);
		//	if (oneLine.find("CELLS") != std::string::npos) {
		//		break;
		//	}
		//}
		//strStream.clear();
		//strStream.str(oneLine);
		//strStream >> oneWrod >> elementNumber >> oneWrod;

		//if (!CreateNewMesh(HEXAHEDRON, vertexNumber, elementNumber)) {
		//	return false;
		//}

		//for (i = 0; i < vertexNumber; ++i)
		//{
		//	vertex[i][0] = pts[i][0];
		//	vertex[i][1] = pts[i][1];
		//	vertex[i][2] = pts[i][2];
		//	//cout << vertex[i][0] << endl;
		//}
		////cout << elementNumber << endl;
		//for (i = 0; i < elementNumber; ++i)
		//{
		//	infile >> oneWrod >> element[i][0] >> element[i][1] >> element[i][2] >> element[i][3] >>
		//		element[i][4] >> element[i][5] >> element[i][6] >> element[i][7];
		//	//cout << element[i][0] << endl;
		//}
		
		infile.close();
	}
	else
	{
		cerr << "Cannot open " << filename << "!\n";
	}
	return true;
}

bool PolycubeConstruction::Initialization()
{
	
	//InitializeSkeletonPoints();
	//SearchNei();

	return true;

}

bool PolycubeConstruction::InitializeSkeletonPoints()
{

	curveSkeleton.InitiateEdgeValence();
	curveSkeleton.InitiateElementValence();

	int startVertex = 0;
	int i, j, index, firstVertex, prevVertex, secondVertex, countTotal;
	double tempDouble[3], tempX[3], tempY[3], tempZ[3], tempXPre[3];

	//First initialization method, may not be correct
	//for (i = 0; i < curveSkeleton.vertexNumber; i++)
	//{
	//	if (curveSkeleton.elementValenceNumber[i] == 1) //elementValence and edgeValence are all useful!
	//	{
	//		
	//		index = curveSkeleton.elementValence[i][0];

	//		if (curveSkeleton.element[index][0] == i)
	//		{
	//			startVertex = i;
	//		}

	//	}
	//}
	
	//Second initialization method, should be correct
	for (i = 0; i < curveSkeleton.vertexNumber; i++)
	{
		if (curveSkeleton.edgeValenceNumber[i] == 1) //elementValence and edgeValence are all useful!
		{

			index = curveSkeleton.edgeValence[i][0];

			if (skeletonPointUsed[index] == false)
			{
				startVertex = i;
			}

		}
	}

	if (SET_INITIAL_SKELETONPOINTS == 1)
	{
		//set the initial curve-skeleton point manually
		//startVertex = 329; //for cactus model
		//startVertex = 663; //for ant model
		startVertex = 439; //for humanB model
	}

	cout<<"The index of the startvertex of the curve-skeleton is "<<startVertex<<endl;

	//countTotal = 0;
	//firstVertex = startVertex;	

	//while (countTotal != curveSkeleton.vertexNumber)
	//{

	//	prevVertex = -1;
	//	secondVertex = -1;

	//	//for (i = 0; i < curveSkeleton.elementValenceNumber[firstVertex]; i++)
	//	//{

	//	//	index = curveSkeleton.elementValence[firstVertex][i];			
	//	//	if (curveSkeleton.element[index][0] == firstVertex)
	//	//	{
	//	//		secondVertex = curveSkeleton.element[index][1];
	//	//	}
	//	//	else if (curveSkeleton.element[index][1] == firstVertex)
	//	//	{
	//	//		prevVertex = curveSkeleton.element[index][0];
	//	//	}

	//	//}

	//	if (curveSkeleton.edgeValenceNumber[firstVertex] == 1)
	//	{
	//		index = curveSkeleton.edgeValence[firstVertex][0];
	//		if (skeletonPointUsed[index] == false)
	//		{
	//			secondVertex = index;
	//		}
	//		else
	//		{
	//			prevVertex = index;
	//		}
	//	}
	//	else
	//	{
	//		for (i = 0; i < curveSkeleton.edgeValenceNumber[firstVertex]; i++)
	//		{
	//			index = curveSkeleton.edgeValence[firstVertex][i];
	//			if (skeletonPointUsed[index] == false)
	//			{
	//				secondVertex = index;
	//			}
	//			else
	//			{
	//				prevVertex = index;
	//			}
	//		}
	//	}

	//	skeletonPoint[firstVertex].index = firstVertex;

	//	for (i = 0; i < 3; i++)
	//	{
	//		tempDouble[i] = curveSkeleton.vertex[secondVertex][i] - curveSkeleton.vertex[firstVertex][i];
	//	}
	//	Normalize(tempDouble);
	//	for (i = 0; i < 3; i++)
	//	{
	//		skeletonPoint[firstVertex].localCoordZ[i] = tempDouble[i];
	//	}

	//	if (prevVertex == -1)
	//	{
	//		for (i = 0; i < 3; i++)
	//		{
	//			tempZ[i] = skeletonPoint[firstVertex].localCoordZ[i];
	//		}
	//		tempXPre[0] = 1.0f; tempXPre[1] = 0.0f; tempXPre[2] = 0.0f;
	//		CalculateLocalCoordinateSystem(tempZ, tempXPre, tempX, tempY);
	//		for (i = 0; i < 3; i++)
	//		{
	//			skeletonPoint[firstVertex].localCoordX[i] = tempX[i];
	//			skeletonPoint[firstVertex].localCoordY[i] = tempY[i];
	//		}
	//	}
	//	else if (secondVertex == -1)
	//	{
	//		for (i = 0; i < 3; i++)
	//		{
	//			skeletonPoint[firstVertex].localCoordX[i] = skeletonPoint[prevVertex].localCoordX[i];
	//			skeletonPoint[firstVertex].localCoordY[i] = skeletonPoint[prevVertex].localCoordY[i];
	//			skeletonPoint[firstVertex].localCoordZ[i] = skeletonPoint[prevVertex].localCoordZ[i];
	//		}
	//	}
	//	else
	//	{
	//		for (i = 0; i < 3; i++)
	//		{
	//			tempZ[i] = skeletonPoint[firstVertex].localCoordZ[i];
	//			tempXPre[i] = skeletonPoint[prevVertex].localCoordX[i];
	//		}
	//		CalculateLocalCoordinateSystem(tempZ, tempXPre, tempX, tempY);
	//		for (i = 0; i < 3; i++)
	//		{
	//			skeletonPoint[firstVertex].localCoordX[i] = tempX[i];
	//			skeletonPoint[firstVertex].localCoordY[i] = tempY[i];
	//		}
	//	}

	//	CalculateRotationMatrix(skeletonPoint[firstVertex]);

	//	skeletonPointUsed[firstVertex] = true;
	//	firstVertex = secondVertex;

	//	countTotal++;

	//}

	RecursiveLocalCoordinateSystem(startVertex);

	OutputLocalCoordinateSystem();

	return true;

}

bool PolycubeConstruction::RecursiveLocalCoordinateSystem(int start_vertex)
{

	//int startVertex = 0;
	int i, j, index, firstVertex, prevVertex, secondVertex, countTotal;
	double tempDouble[3], tempX[3], tempY[3], tempZ[3], tempXPre[3];
	bool stopCriteria;

	stopCriteria = false;
	firstVertex = start_vertex;

	while (stopCriteria == false)
	{

		prevVertex = -1;
		secondVertex = -1;

		if (curveSkeleton.edgeValenceNumber[firstVertex] == 1)
		{
			index = curveSkeleton.edgeValence[firstVertex][0];
			if (skeletonPointUsed[index] == false)
			{
				secondVertex = index;
			}
			else
			{
				prevVertex = index;
			}
		}
		else
		{
			for (i = 0; i < curveSkeleton.edgeValenceNumber[firstVertex]; i++)
			{
				index = curveSkeleton.edgeValence[firstVertex][i];
				if (skeletonPointUsed[index] == false)
				{
					secondVertex = index;
				}
				else
				{
					prevVertex = index;
				}
			}
		}

		skeletonPoint[firstVertex].index = firstVertex;

		for (i = 0; i < 3; i++)
		{
			tempDouble[i] = curveSkeleton.vertex[secondVertex][i] - curveSkeleton.vertex[firstVertex][i];
		}
		Normalize(tempDouble);
		for (i = 0; i < 3; i++)
		{
			skeletonPoint[firstVertex].localCoordZ[i] = tempDouble[i];
		}

		if (prevVertex == -1)
		{
			for (i = 0; i < 3; i++)
			{
				tempZ[i] = skeletonPoint[firstVertex].localCoordZ[i];
			}
			tempXPre[0] = 1.0f; tempXPre[1] = 0.0f; tempXPre[2] = 0.0f;
			CalculateLocalCoordinateSystem(tempZ, tempXPre, tempX, tempY);
			for (i = 0; i < 3; i++)
			{
				skeletonPoint[firstVertex].localCoordX[i] = tempX[i];
				skeletonPoint[firstVertex].localCoordY[i] = tempY[i];
			}
		}
		else if (secondVertex == -1)
		{
			for (i = 0; i < 3; i++)
			{
				skeletonPoint[firstVertex].localCoordX[i] = skeletonPoint[prevVertex].localCoordX[i];
				skeletonPoint[firstVertex].localCoordY[i] = skeletonPoint[prevVertex].localCoordY[i];
				skeletonPoint[firstVertex].localCoordZ[i] = skeletonPoint[prevVertex].localCoordZ[i];
			}
		}
		else
		{
			for (i = 0; i < 3; i++)
			{
				tempZ[i] = skeletonPoint[firstVertex].localCoordZ[i];
				tempXPre[i] = skeletonPoint[prevVertex].localCoordX[i];
			}
			CalculateLocalCoordinateSystem(tempZ, tempXPre, tempX, tempY);
			for (i = 0; i < 3; i++)
			{
				skeletonPoint[firstVertex].localCoordX[i] = tempX[i];
				skeletonPoint[firstVertex].localCoordY[i] = tempY[i];
			}
		}

		CalculateRotationMatrix(skeletonPoint[firstVertex]);

		skeletonPointUsed[firstVertex] = true;

		if (curveSkeleton.edgeValenceNumber[firstVertex] == 1)
		{

			index = curveSkeleton.edgeValence[firstVertex][0];
			if (skeletonPointUsed[index] == true)
			{
				stopCriteria = true;
			}

		}

		firstVertex = secondVertex;
	}

	//check all curve-skeleton points
	bool stopCriteriaAll = true;
	start_vertex = -1;
	for (i = 0; i < curveSkeleton.vertexNumber; i++)
	{
		if (skeletonPointUsed[i] == false)
		{

			stopCriteriaAll = false;
			break;

		}
	}

	if (stopCriteriaAll == false)
	{
		for (i = 0; i < curveSkeleton.vertexNumber; i++)
		{
			if (skeletonPointUsed[i] == false)
			{
				for (j = 0; j < curveSkeleton.edgeValenceNumber[i]; j++)
				{
					index = curveSkeleton.edgeValence[i][j];
					if (skeletonPointUsed[index] == true)
					{
						start_vertex = i;
						break;
					}
				}

				if (start_vertex != -1)
				{
					break;
				}
			}
		}

		RecursiveLocalCoordinateSystem(start_vertex);
		if (start_vertex == -1)
		{
			cout<<"Error in RecursiveLocalCoordinateSystem(start_vertex)!"<<endl;
		}
	}
	else
	{
		return true;
	}

}

//Given vz and vxpre, n=vz*vxpre, tv=n*vz, vx=xpre-tv, vy=xz x vx
bool PolycubeConstruction::CalculateLocalCoordinateSystem(double vz[3], double vxpre[3], double vx[3], double vy[3])
{

	int i;
	double tempLength;
	double tempArray[3];

	tempLength = vxpre[0]*vz[0] + vxpre[1]*vz[1] + vxpre[2]*vz[2];
	for (i = 0; i < 3; i++)
	{
		tempArray[i] = tempLength * vz[i];
	}

	for (i = 0; i < 3; i++)
	{
		vx[i] = vxpre[i] - tempArray[i];
	}
	Normalize(vx);

	CrossProduct(vz, vx, vy);
	Normalize(vy);

	return true;

}

bool PolycubeConstruction::CalculateRotationMatrix(SkeletonPoint &sp)
{

	sp.transMatrix[0][0] = sp.localCoordX[0]*global_CoordSys[0][0] + sp.localCoordX[1]*global_CoordSys[0][1] + sp.localCoordX[2]*global_CoordSys[0][2];
	sp.transMatrix[0][1] = sp.localCoordX[0]*global_CoordSys[1][0] + sp.localCoordX[1]*global_CoordSys[1][1] + sp.localCoordX[2]*global_CoordSys[1][2];
	sp.transMatrix[0][2] = sp.localCoordX[0]*global_CoordSys[2][0] + sp.localCoordX[1]*global_CoordSys[2][1] + sp.localCoordX[2]*global_CoordSys[2][2];

	sp.transMatrix[1][0] = sp.localCoordY[0]*global_CoordSys[0][0] + sp.localCoordY[1]*global_CoordSys[0][1] + sp.localCoordY[2]*global_CoordSys[0][2];
	sp.transMatrix[1][1] = sp.localCoordY[0]*global_CoordSys[1][0] + sp.localCoordY[1]*global_CoordSys[1][1] + sp.localCoordY[2]*global_CoordSys[1][2];
	sp.transMatrix[1][2] = sp.localCoordY[0]*global_CoordSys[2][0] + sp.localCoordY[1]*global_CoordSys[2][1] + sp.localCoordY[2]*global_CoordSys[2][2];

	sp.transMatrix[2][0] = sp.localCoordZ[0]*global_CoordSys[0][0] + sp.localCoordZ[1]*global_CoordSys[0][1] + sp.localCoordZ[2]*global_CoordSys[0][2];
	sp.transMatrix[2][1] = sp.localCoordZ[0]*global_CoordSys[1][0] + sp.localCoordZ[1]*global_CoordSys[1][1] + sp.localCoordZ[2]*global_CoordSys[1][2];
	sp.transMatrix[2][2] = sp.localCoordZ[0]*global_CoordSys[2][0] + sp.localCoordZ[1]*global_CoordSys[2][1] + sp.localCoordZ[2]*global_CoordSys[2][2];

	return true;
}

bool PolycubeConstruction::OutputLocalCoordinateSystem()
{

	string inputFileName = "local_coordSys.txt";

	fstream output(inputFileName, fstream::out);
	if (!output)
	{
		cout << "open input file error!!!" <<endl;
		return false;
	}

	int i;

	output << curveSkeleton.vertexNumber << "\n";

	output.precision(9);
	for (i = 0; i < curveSkeleton.vertexNumber; i++)
	{
		output << curveSkeleton.vertex[i][0]<<"\t"<<curveSkeleton.vertex[i][1]<<"\t"<<curveSkeleton.vertex[i][2]<<"\t"<<
			skeletonPoint[i].localCoordX[0]<<"\t"<<skeletonPoint[i].localCoordX[1]<<"\t"<<skeletonPoint[i].localCoordX[2]<<"\t"<<
			skeletonPoint[i].localCoordY[0]<<"\t"<<skeletonPoint[i].localCoordY[1]<<"\t"<<skeletonPoint[i].localCoordY[2]<<"\t"<<
			skeletonPoint[i].localCoordZ[0]<<"\t"<<skeletonPoint[i].localCoordZ[1]<<"\t"<<skeletonPoint[i].localCoordZ[2]<<"\n";
	}

	output.close();

	return true;

}


bool PolycubeConstruction::InitializeElement(const char * filename)
{

	double tempVec[3];
	int i, j;

	double tempDist, dist; // Voronoi Diagram in physical space
	int indexNearestSkeleton;

	ComputeNormal(2); //compute normal of each element
	if (elementValence == NULL)
	{
		InitiateElementValence();
	}
	

	string tempString;
	tempString = filename;
	ReadKFileInitial(tempString.c_str());
	ReadKFileIndexPatch(tempString.c_str());

	for (i = 0; i < elementNumber; i++)
	{

		// Assign the index to each element
		elementArray[i].index = i;


		for (j = 0; j < 3; j++)
		{

			tempVec[j] = elementNormal[i][j];

		}

		Normalize(tempVec);

		for (j = 0; j < 3; j++)
		{

			elementArray[i].normal[j] = tempVec[j];

		}

		GetDirectNeighboringElementByElement(i);

		//Find all neighboring elements

		elementArray[i].numNeiElements = 1;

		elementArray[i].neiElements.push_back(i);

		//End

		////Assign the triangle to the nearest skeleton point and update the new normal vector
		//tempDist = 10000000.0f;
		//indexNearestSkeleton = 0;

		//for (j = 0; j < curveSkeleton.vertexNumber; j++)
		//{
		//	dist = GetPhysicalDist(j, elementArray[i]);

		//	if (dist < tempDist)
		//	{
		//		indexNearestSkeleton = j;
		//		tempDist = dist;
		//	}
		//}

		//elementArray[i].indexSkeletonPoint = indexNearestSkeleton;

		////Calculate new normal for each triangle through a rotation matrix
		//int tempIndex = elementArray[i].indexSkeletonPoint;
		//for (j = 0; j < 3; j++)
		//{
		//	elementArray[i].normalNew[j] = skeletonPoint[tempIndex].transMatrix[j][0] * elementArray[i].normal[0] + skeletonPoint[tempIndex].transMatrix[j][1] * elementArray[i].normal[1]
		//	+ skeletonPoint[tempIndex].transMatrix[j][2] * elementArray[i].normal[2]; // be careful about i and j here
		//}

	}

	GetNeighboringElementInRings();
	//getchar();
	return true;

}

bool PolycubeConstruction::AssignPartToElement()
{

	string tempName;
	string oneLine;
	tempName = inputName + "_elementSign.txt";

	ifstream myFile(tempName);

	int tempPart;
	int i = 0;

	if (myFile.is_open())
	{

		while (getline(myFile, oneLine))
		{

			istringstream st(oneLine);

			st >> tempPart; //>>, not <<

			elementArray[i].indexPart = tempPart;

			i++;
		}

		myFile.close();

	}
	else
	{
		cout<<"Unable to open file!";
	}


	return true;

}

bool PolycubeConstruction::AssignClusterToElement()
{

	string tempName;
	string oneLine;
	tempName = inputName + "_elementSignCluster.txt";

	ifstream myFile(tempName);

	int tempPart;
	int i = 0;

	if (myFile.is_open())
	{

		while (getline(myFile, oneLine))
		{

			istringstream st(oneLine);

			st >> tempPart; //>>, not <<

			elementArray[i].indexCluster = tempPart;

			i++;
		}

		myFile.close();

	}
	else
	{
		cout<<"Unable to open file!";
	}


	return true;

}

double PolycubeConstruction::GetPhysicalDist(int index, const CVTElement &currentElement)
{

	double dist = 0.0f;

	double vec1[3] = {0.};
	double vec2[3] = {0.};

	int i;
	int tempIndex = currentElement.index;

	for (i = 0; i < 3; i++)
	{
		vec1[i] = curveSkeleton.vertex[index][i];

		for (int j = 0; j < 3; j++)
		{
			vec2[i] += vertex[element[tempIndex][j]][i];
		}
		vec2[i] /= 3.0f;
	}

	for (i = 0; i < 3; i++)
	{
		dist = dist + (vec1[i]-vec2[i]) * (vec1[i]-vec2[i]);
	}

	dist = sqrt(dist);

	return dist;

}

bool PolycubeConstruction::GetDirectNeighboringElementByElement(int elementID)
{

	int i, iVert, j, jElem, k, kSize;

	vector<int> tempElements;

	for (i = 0; i < 3; i++)
	{

		iVert = element[elementID][i];

		for (j = 0; j < elementValenceNumber[iVert]; j++)
		{

			jElem = elementValence[iVert][j];

			kSize = tempElements.size();

			for (k = 0; k < kSize; k++)
			{

				if (jElem == elementID || jElem == tempElements[k])
				{

					break;

				}

			}
			if (k == kSize && jElem != elementID) //Very important when kSize = 0
			{

				tempElements.push_back(jElem);

			}

		}

	}

	//Find three direct neighboring elements
	kSize = 0;

	for (i = 0; i < tempElements.size(); i++)
	{

		for (j = 0; j < 3; j++)
		{

			iVert = element[tempElements[i]][j];

			for (k = 0; k < 3; k++)
			{

				if (iVert == element[elementID][k])
				{
					kSize++;
				}

			}

		}

		if (kSize == 2)
		{

			elementArray[elementID].directNei.push_back(tempElements[i]);

		}

		kSize = 0;

	}

	elementArray[elementID].numDirectNei = elementArray[elementID].directNei.size();

	return true;

}


//Recursive function to find all neighboring elements in N rings
bool PolycubeConstruction::GetNeighboringElementInRings(int ringNumber)
{

	int i, j, iSize, k, kk, kSize;

	int tempVertex, tempElement;
	//vector<int> tempNodes;

	for (i = 0; i < elementNumber; i++)
	{

		vector<int> tempNodes;

		iSize = elementArray[i].neiElements.size();

		for (j = 0; j < iSize; j++)
		{

			for (k = 0; k < 3; k++)
			{

				tempVertex = element[elementArray[i].neiElements[j]][k];

				kSize = tempNodes.size();

				for (kk = 0; kk < kSize; kk++)
				{

					if (tempVertex == tempNodes[kk])
					{
						break;
					}

				}
				if (kk == kSize)
				{
					tempNodes.push_back(tempVertex);
				}

			}

		}

		//Next step
		iSize = tempNodes.size();

		for (j = 0; j < iSize; j++)
		{

			tempVertex = tempNodes[j];

			for (k = 0; k < elementValenceNumber[tempVertex]; k++)
			{

				tempElement = elementValence[tempVertex][k];

				kSize = elementArray[i].neiElements.size();

				for (kk = 0; kk < kSize; kk++)
				{

					if (tempElement == elementArray[i].neiElements[kk])
					{
						break;
					}

				}
				if (kk == kSize)
				{
					elementArray[i].neiElements.push_back(tempElement);
				}

			}

		}

		elementArray[i].numNeiElements = elementArray[i].neiElements.size(); //numNeiElements
	}

	if (ringNumber < NUM_RING)
	{
		ringNumber++;
		GetNeighboringElementInRings(ringNumber);

	}
	else
	{

		return true;
	}
	//return true;

}


bool PolycubeConstruction::InitializeSegments()
{

	int i, j, k;

	InitializeGeneratorsByInput();

	// Initialize the segments
	double tempDist, dist;
	int indexNearestGenerator;

	for (i = 0; i < elementNumber; i++)
	{

		tempDist = 10000000.0f;
		indexNearestGenerator = 0;

		for (j = 0; j < NUM_CLUSTER; j++)
		{

			dist = GetNormalDist(generators[j], elementArray[i]);

			if (dist < tempDist)
			{

				indexNearestGenerator = j;
				tempDist = dist;

			}

		}

		// Initialize more info for each element
		elementArray[i].numNeiCluster = 1;

		for (k = 0; k < NUM_NEI_CLUSTER; k++)
		{

			elementArray[i].numNeiElementEachCluster[k] = 0;
			elementArray[i].indexNeiClusters[k] = -1;

		}

		elementArray[i].indexCluster = indexNearestGenerator;

		elementArray[i].indexNeiClusters[0] = elementArray[i].indexCluster;
		elementArray[i].numNeiElementEachCluster[0] = 1;

		generators[indexNearestGenerator].numElements++;

	}

	//Recalculate the generators

	int indexCluster;

	for (i = 0; i < NUM_CLUSTER; i++)
	{

		for (j = 0; j < 3; j++)
		{

			generators[i].normal[j] = 0.f;

		}

	}

	for (i = 0; i < elementNumber; i++)
	{

		indexCluster = elementArray[i].indexCluster;

		for (j = 0; j < 3; j++)
		{

			//generators[indexCluster].normal[j] += elementArray[i].normal[j];
			generators[indexCluster].normal[j] += elementArray[i].normalNew[j];

		}

	}

	for (i = 0; i < NUM_CLUSTER; i++)
	{

		for (j = 0; j < 3; j++)
		{

			generators[i].normal[j] /= generators[i].numElements;

		}

		NormalizeGeneratorNormal(generators[i]);

	}

	//end

	return true;

}


bool PolycubeConstruction::InitializeGeneratorsByInput()
{

	int i, j;

	for (i = 0; i < NUM_CLUSTER; i++)
	{

		generators[i].index = i;

		generators[i].numElements = 0;

		for (j = 0; j < 3; j++)
		{

			generators[i].normal[j] = prin_normal[i][j];

			//outputCentroids<<generators[i].normal[j]<<",";

		}

		//outputCentroids<<endl;

	}

	//outputCentroids<<endl;

	return true;

}

double PolycubeConstruction::GetNormalDist(const Centroid &currentGenerator, const CVTElement & currentElement)
{

	double dist = 0.0f;

	double vec1[3] = {0.};
	double vec2[3] = {0.};

	int i;

	for (i = 0; i < 3; i++)
	{

		vec1[i] = currentGenerator.normal[i];
		//vec2[i] = currentElement.normal[i];
		vec2[i] = currentElement.normalNew[i];

	}

	dist = NormalDist(vec1, vec2);

	dist = dist / 2.0; // normalize the distance to [0, 1]

	return dist;

}

double PolycubeConstruction::NormalDist(double veca[], double vecb[])
{

	double dist = 0.;

	int i = 0;

	for (i = 0; i < 3; i++)
	{

		dist = dist + (veca[i] - vecb[i]) * (veca[i] - vecb[i]);

	}

	dist = sqrt(dist);

	return dist;

}

bool PolycubeConstruction::NormalizeGeneratorNormal(Centroid &currentGenerator)
{

	double sum = 0.;

	int i;

	for (i = 0; i < 3; i++)
	{
		sum += currentGenerator.normal[i]*currentGenerator.normal[i];
	}

	sum = sqrt(sum);

	for (i = 0; i < 3; i++)
	{
		currentGenerator.normal[i] /= sum;
	}

	return true;

}

bool PolycubeConstruction::SearchNei()
{

	int neiIndex;
	int neiClusterIndex;
	int i, j;
	int position;

	for (i = 0; i < elementNumber; i++)
	{

		for (j = 0; j < elementArray[i].numNeiElements; j++)
		{

			neiIndex = elementArray[i].neiElements[j];

			if (neiIndex != i)
			{

				neiClusterIndex = elementArray[neiIndex].indexCluster;

				IsCounted(elementArray[i], neiClusterIndex, position);

			}

		}

	}

	//Verification

	int error = 0;
	int numNei = 0;

	for (i = 0; i < elementNumber; i++)
	{

		numNei = 0;

		for (j = 0; j < elementArray[i].numNeiCluster; j++)
		{

			numNei = numNei + elementArray[i].numNeiElementEachCluster[j];

		}

		if (numNei != elementArray[i].numNeiElements)
		{

			error++;

		}
	}

	printf("%d errors occured in searching the neighbor pixels\n", error);

	return true;

}

bool PolycubeConstruction::IsCounted(CVTElement &currentElement, int indexClusterNeiElement, int &position)
{

	int i;

	bool counted = false; // Nice method!!! Learn how to use it!
	int numNeiCluster;

	for (i = 0; i < currentElement.numNeiCluster; i++)
	{

		if (currentElement.indexNeiClusters[i] == indexClusterNeiElement)
		{

			counted = true;
			currentElement.numNeiElementEachCluster[i]++;
			position = i;

			break;

		}

	}

	if (counted == false)
	{

		currentElement.numNeiCluster++;
		numNeiCluster = currentElement.numNeiCluster;

		if (numNeiCluster > NUM_NEI_CLUSTER)
		{

			printf("Warning! the predefined number of neighbor clusters is not OK!\n");
			printf("element index: %d\n", currentElement.index);
			exit(0);

		}
		else
		{

			currentElement.indexNeiClusters[numNeiCluster - 1] = indexClusterNeiElement;

			currentElement.numNeiElementEachCluster[numNeiCluster - 1]++;

			position = numNeiCluster - 1;

		}


	}

	return counted;

}


bool PolycubeConstruction::OutputPatchesVTK(const char *outputName)
{

	string inputFileName = outputName;

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

	//start new stuff
	output <<"CELL_DATA "<<elementNumber<<endl;
	output <<"SCALARS patches float"<<endl;
	output <<"LOOKUP_TABLE default"<<endl;

	for (i = 0; i < elementNumber; i++)
	{
		
		float tempFloat = (float) elementArray[i].indexPatch;
		output << tempFloat <<endl;

	}

	//end new stuff

	output.close();


	return true;

}

bool PolycubeConstruction::OutputPatchesVTKPara(const char *outputName)
{

	string inputFileName = outputName;

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
	output << "POINTS " << polycubePara->vertexNumber<<" FLOAT"<<"\n";


	output.precision(9);
	for (i=0; i<polycubePara->vertexNumber; ++i)
		output <<polycubePara->vertex[i][0]<<"\t"<<polycubePara->vertex[i][1]<<"\t"<<polycubePara->vertex[i][2]<<"\n";
	output.unsetf(ostream::floatfield);
	output << "CELLS " << polycubePara->elementNumber <<" "<<(polycubePara->elementProperty.vertexNumber+1) * polycubePara->elementNumber<<"\n";

	for (i=0; i<polycubePara->elementNumber; ++i)
	{
		output << polycubePara->elementProperty.vertexNumber<<" ";
		for(j=0; j<polycubePara->elementProperty.vertexNumber; j++)
		{
			output<< polycubePara->element[i][j]<<" ";
		}
		output<<"\n";
	}
	output <<"CELL_TYPES "<<polycubePara->elementNumber<<endl;
	for (i=0; i<polycubePara->elementNumber; ++i)
	{
		if (polycubePara->elementProperty.elementType == NA)
		{
			continue;
		}
		else if (polycubePara->elementProperty.elementType == POINT)
		{
			output << 1 <<endl;
		}
		else if (polycubePara->elementProperty.elementType == LINE)
		{
			output << 3 <<endl;
		}
		else if (polycubePara->elementProperty.elementType == TRIANGLE)
		{
			output << 5 <<endl;
		}
		else if (polycubePara->elementProperty.elementType == TETRAHEDRON)
		{
			output << 10 <<endl;
		}
		else if (polycubePara->elementProperty.elementType == QUADRILATERAL)
		{
			output << 9 <<endl;
		}
		else if (polycubePara->elementProperty.elementType == HEXAHEDRON)
		{
			output << 12 <<endl;
		}
		else if (polycubePara->elementProperty.elementType == POLYGON)
		{
			output << 7 <<endl;
		}
		else if (polycubePara->elementProperty.elementType == HEXAGON)
		{
			continue;
		}

	}

	//start new stuff
	output <<"CELL_DATA "<<polycubePara->elementNumber<<endl;
	output <<"SCALARS patches float"<<endl;
	output <<"LOOKUP_TABLE default"<<endl;

	for (i = 0; i < polycubePara->elementNumber; i++)
	{

		float tempFloat = (float) elementArray[i].indexPatch;
		output << tempFloat <<endl;

	}

	//end new stuff

	output.close();


	return true;

}

bool PolycubeConstruction::OutputPatchesVTKBifurcation(const char *outputName)
{

	string inputFileName = outputName;

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

	//start new stuff
	output <<"CELL_DATA "<<elementNumber<<endl;
	output <<"SCALARS patches float"<<endl;
	output <<"LOOKUP_TABLE default"<<endl;

	for (i = 0; i < elementNumber; i++)
	{

		float tempFloat = (float) elementArray[i].indexPatch;
		output << tempFloat <<endl;

	}

	//end new stuff

	output.close();


	return true;

}

bool PolycubeConstruction::EdgeWeightedCVT()
{

	cout<<"**********************************************************"<<endl;
	cout<<"The Edge-weighted CVT algorithm is running......"<<endl;

	int i, j, k;

	int numTransfer = elementNumber;
	int threshold = 3;

	int newNearestGenerator;
	int oldNearestGenerator;

	int step = 0;

	while(numTransfer >= threshold)
	{

		numTransfer = 0;

		for (i = 0; i < elementNumber; i++)
		{

			if (IsBoundaryElement(elementArray[i]))
				//if(true)
			{

				newNearestGenerator = GetShortestEWDist(elementArray[i]);

				oldNearestGenerator = elementArray[i].indexCluster;

				if (newNearestGenerator != oldNearestGenerator)
				{

					DataTransfer(elementArray[i], newNearestGenerator);

					//UpdateGenerator(generators[newNearestGenerator], elementArray[i], 1);

					//UpdateGenerator(generators[oldNearestGenerator], elementArray[i], 0);

					numTransfer++;

				}

			}

		}

		for (j = 0; j < NUM_CLUSTER; j++)
		{

			for (k = 0; k < 3; k++)
			{

				//outputCentroids<<generators[j].normal[k]<<",";

			}

			//outputCentroids<<endl;

		}

		//outputCentroids<<endl;

		step++;
		printf("step %d    %d\n", step, numTransfer);
		if(step > 100)
		{
			break;
		}

	}

	printf("Edge-weighted CVT is done!\n");

	return true;

}

bool PolycubeConstruction::DataTransfer(CVTElement &currentElement, int newIndex)
{

	int oldIndex = currentElement.indexCluster;

	int i;
	int newIndexPosition;
	int oldIndexPosition;

	bool oldIndexCounted = false;

	IsCounted(currentElement, newIndex, newIndexPosition);

	for (i = 0; i < currentElement.numNeiCluster; i++)
	{

		if (currentElement.indexNeiClusters[i] == oldIndex)
		{

			currentElement.numNeiElementEachCluster[i]--;

			oldIndexPosition = i;
			oldIndexCounted = true;

			break;

		}

	}

	if(oldIndexCounted == false)
	{
		printf("warning! can not find the cluster for the current pixel in the indexNeiClusters array!\n");
		exit(0);
	}

	oldIndexCounted = false;

	int j;
	int indexNeiData;

	for (i = 0; i < currentElement.numNeiElements; i++)
	{

		indexNeiData = currentElement.neiElements[i];

		if (indexNeiData != currentElement.index)
		{

			IsCounted(elementArray[indexNeiData], newIndex, newIndexPosition);

			for (j = 0; j < elementArray[indexNeiData].numNeiCluster; j++)
			{

				if (elementArray[indexNeiData].indexNeiClusters[j] == oldIndex)
				{

					elementArray[indexNeiData].numNeiElementEachCluster[j]--;

					oldIndexPosition = j;
					oldIndexCounted = true;


					//test

					if (elementArray[indexNeiData].numNeiElementEachCluster[j] == -1)
					{

						j = j;

					}

					//end test

					break;

				}

			}

			if(oldIndexCounted == false)
			{
				printf("warning! can not find the cluster for the current pixel in the indexNeiClusters array!\n");
				exit(0);
			}

		}

	}

	currentElement.indexCluster = newIndex;

	return true;
}

int PolycubeConstruction::GetShortestEWDist(CVTElement &currentElement)
{

	double shortestEWDist = 0.;
	double tempshortestEWDist = 0.;

	int indexCluster_shortestEWDist = currentElement.indexCluster;
	int currentCluster = currentElement.indexCluster;
	int indexNeiCluster;

	int i;

	shortestEWDist = GetEWDist(generators[currentCluster], currentElement);

	for (i = 0; i < currentElement.numNeiCluster; i++)
	{

		if (currentElement.numNeiElementEachCluster[i] != 0)
		{

			if (currentCluster != currentElement.indexNeiClusters[i])
			{
				indexNeiCluster = currentElement.indexNeiClusters[i];
				tempshortestEWDist = GetEWDist(generators[indexNeiCluster], currentElement);

				if (tempshortestEWDist < shortestEWDist)
				{
					shortestEWDist = tempshortestEWDist;
					indexCluster_shortestEWDist = indexNeiCluster;
				}

			}

		}

	}

	return indexCluster_shortestEWDist;

}

double PolycubeConstruction::GetEWDist(const Centroid &currentGenerator, const CVTElement & currentElement)
{

	double EWDist = 0.;

	int j;
	int indexNeiCluster = currentGenerator.index;
	int currentCluster = currentElement.indexCluster;

	EWDist = GetNormalDist(currentGenerator, currentElement);

	EWDist = EWDist * EWDist;

	for (j = 0; j < currentElement.numNeiCluster; j++)
	{

		if (currentElement.indexNeiClusters[j] == indexNeiCluster)
		{

			break;

		}

	}

	double normalizedEWPart = 0.;

	if (currentCluster == indexNeiCluster)
	{

		normalizedEWPart = 2 * WEIGHT_LENGTH_EWCVT * (currentElement.numNeiElements - currentElement.numNeiElementEachCluster[j]);

		normalizedEWPart /= currentElement.numNeiElements;

		EWDist += normalizedEWPart;

	}
	else
	{

		normalizedEWPart = 2 * WEIGHT_LENGTH_EWCVT * (currentElement.numNeiElements - currentElement.numNeiElementEachCluster[j] - 1);

		normalizedEWPart /= currentElement.numNeiElements;

		EWDist += normalizedEWPart;

	}

	EWDist = sqrt(EWDist);

	return EWDist;

}


bool PolycubeConstruction::PostProcessing()
{
	
	////////////////////////////////////////////////////////////////////////
	if (ASSIGN_CLUSTER_ELEMENT == 1)
	{
		AssignClusterToElement();
	}
	////////////////////////////////////////////////////////////////////////

	if (READ_WRITE_K_BEFORE_POSTPROCESSING == 1)
	{
		string tempString = inputName + "_initial_read.k";
		ReadKFileBeforePostProcessing(tempString.c_str());
	}
	else if (READ_WRITE_K_BEFORE_POSTPROCESSING == 0)
	{
		string tempString = inputName + "_initial_write.k";
		WriteKFileBeforePostProcessing(tempString.c_str());
	}

	EnforceLabelConnectivity();
	EnforceBoundaryConnectivity();

	return true;

}

bool PolycubeConstruction::IsBoundaryElement(CVTElement &currentElement)
{

	int i;
	int neiIndex;

	bool isOnBoundary = false;

	for (i = 0; i < currentElement.numDirectNei; i++)
	{

		neiIndex = currentElement.directNei[i];

		if (currentElement.indexCluster != elementArray[neiIndex].indexCluster)
		{

			isOnBoundary = true;
			break;

		}

	}

	return isOnBoundary;

}

bool PolycubeConstruction::EnforceLabelConnectivity()
{

	cout<<"Label connectivity needs to be enforced!"<<endl;

	int i, j, k, index;
	int label = 0;

	//const int MAX_NUM_CLUSTER = (int) sqrt(elementNumber); //Error may occur if number of clusters is vergy large

	vector<vector<int> > tempCluster(elementNumber);

	for (i = 0; i < elementNumber; i++)
	{

		tempCluster[i].push_back(-1);
		tempCluster[i].push_back(elementArray[i].indexCluster);

	}

	for (i = 0; i < elementNumber; i++)
	{

		if (0 > tempCluster[i][0])
		{

			tempCluster[i][0] = label;

			int count = 1;

			vector<int> tempIndex;

			tempIndex.push_back(i);

			for (int c = 0; c < count; c++)
			{

				for (int n = 0; n < elementArray[tempIndex[c]].numDirectNei; n++)
				{

					index = elementArray[tempIndex[c]].directNei[n];

					if (0 > tempCluster[index][0] && elementArray[i].indexCluster == elementArray[index].indexCluster)
					{

						tempCluster[index][0] = label;

						tempIndex.push_back(index);

						count++;

					}

				}

			}

			label++;

		}

	}

	const int MAX_NUM_CLUSTER = label; // total number of seperated clusters

	vector<CVTElement> elementArrayCopy;

	elementArrayCopy = elementArray;

	for (i = 0; i < elementNumber; i++)
	{

		elementArrayCopy[i].indexCluster = tempCluster[i][0];

	}

	vector<vector<int> > elementInCluster;
	vector<vector<int> > neighborCluster;
	int indexCluster, indexClusterNei;

	elementInCluster.resize(MAX_NUM_CLUSTER);
	neighborCluster.resize(MAX_NUM_CLUSTER, vector<int>(MAX_NUM_CLUSTER, 0));

	for (i = 0; i < elementNumber; i++)
	{

		//Search for all elements in each cluster
		elementInCluster[tempCluster[i][0]].push_back(i);

		if (IsBoundaryElement(elementArrayCopy[i]))
		{

			for (j = 0; j < elementArrayCopy[i].numDirectNei; j++)
			{

				index = elementArrayCopy[i].directNei[j];

				if (elementArrayCopy[i].indexCluster != elementArrayCopy[index].indexCluster)
				{

					indexCluster = elementArrayCopy[i].indexCluster;
					indexClusterNei = elementArrayCopy[index].indexCluster;

					neighborCluster[indexCluster][indexClusterNei]++;

				}

			}

		}

	}

	int kSize = 0;
	int maxClusterIndex;
	int indexNei;

	for (i = 0; i < MAX_NUM_CLUSTER; i++)
	{

		kSize = 0;

		for (j = 0; j < MAX_NUM_CLUSTER; j++)
		{

			if (neighborCluster[i][j] > 0)
			{

				kSize++;

			}

		}

		if (kSize < 4)
		{

			maxClusterIndex = 0;

			for (j = 0; j < MAX_NUM_CLUSTER; j++)
			{

				if (neighborCluster[i][j] > neighborCluster[i][maxClusterIndex])
				{

					maxClusterIndex = j;

				}

			}

			indexNei = elementInCluster[maxClusterIndex][0];

			for (k = 0; k < elementInCluster[i].size(); k++)
			{

				index = elementInCluster[i][k];

				elementArray[index].indexCluster = tempCluster[indexNei][1];

			}


		}


	}

	if (CheckLabelConnectivity() == false)
	{

		EnforceLabelConnectivity();

	}
	else
	{

		cout<<"Label connectivity is now correct!"<<endl;

		return true;

	}

	//return true;

}

bool PolycubeConstruction::CheckLabelConnectivity()
{

	int i, j, k, index;
	int label = 0;

	vector<vector<int> > tempCluster(elementNumber);

	for (i = 0; i < elementNumber; i++)
	{

		tempCluster[i].push_back(-1);
		tempCluster[i].push_back(elementArray[i].indexCluster);

	}

	for (i = 0; i < elementNumber; i++)
	{

		if (0 > tempCluster[i][0])
		{

			tempCluster[i][0] = label;

			int count = 1;

			vector<int> tempIndex;

			tempIndex.push_back(i);

			for (int c = 0; c < count; c++)
			{

				for (int n = 0; n < elementArray[tempIndex[c]].numDirectNei; n++)
				{

					index = elementArray[tempIndex[c]].directNei[n];

					if (0 > tempCluster[index][0] && elementArray[i].indexCluster == elementArray[index].indexCluster)
					{

						tempCluster[index][0] = label;

						tempIndex.push_back(index);

						count++;

					}

				}

			}

			label++;

		}

	}

	const int MAX_NUM_CLUSTER = label; // total number of seperated clusters

	vector<CVTElement> elementArrayCopy;

	elementArrayCopy = elementArray;

	for (i = 0; i < elementNumber; i++)
	{

		elementArrayCopy[i].indexCluster = tempCluster[i][0];

	}

	vector<vector<int> > elementInCluster;
	vector<vector<int> > neighborCluster;
	int indexCluster, indexClusterNei;

	elementInCluster.resize(MAX_NUM_CLUSTER);
	neighborCluster.resize(MAX_NUM_CLUSTER, vector<int>(MAX_NUM_CLUSTER, 0));

	for (i = 0; i < elementNumber; i++)
	{

		//Search for all elements in each cluster
		elementInCluster[tempCluster[i][0]].push_back(i);

		if (IsBoundaryElement(elementArrayCopy[i]))
		{

			for (j = 0; j < elementArrayCopy[i].numDirectNei; j++)
			{

				index = elementArrayCopy[i].directNei[j];

				if (elementArrayCopy[i].indexCluster != elementArrayCopy[index].indexCluster)
				{

					indexCluster = elementArrayCopy[i].indexCluster;
					indexClusterNei = elementArrayCopy[index].indexCluster;

					neighborCluster[indexCluster][indexClusterNei]++;

				}

			}

		}

	}

	int kSize = 0;
	int maxClusterIndex;
	int indexNei;

	for (i = 0; i < MAX_NUM_CLUSTER; i++)
	{

		kSize = 0;

		for (j = 0; j < MAX_NUM_CLUSTER; j++)
		{

			if (neighborCluster[i][j] > 0)
			{

				kSize++;

			}

		}

		if (kSize < 4)
		{

			return false;

		}


	}


	return true;

}

bool PolycubeConstruction::EnforceBoundaryConnectivity()
{

	cout<<"Boundary connectivity needs to be enforced!"<<endl;

	int i, j, k, index;
	int label = 0;

	vector<vector<int> > tempCluster(elementNumber);

	for (i = 0; i < elementNumber; i++)
	{

		tempCluster[i].push_back(-1);
		tempCluster[i].push_back(elementArray[i].indexCluster);

	}

	for (i = 0; i < elementNumber; i++)
	{

		if (0 > tempCluster[i][0])
		{

			tempCluster[i][0] = label;

			int count = 1;

			vector<int> tempIndex;

			tempIndex.push_back(i);

			for (int c = 0; c < count; c++)
			{

				for (int n = 0; n < elementArray[tempIndex[c]].numDirectNei; n++)
				{

					index = elementArray[tempIndex[c]].directNei[n];

					if (0 > tempCluster[index][0] && elementArray[i].indexCluster == elementArray[index].indexCluster)
					{

						tempCluster[index][0] = label;

						tempIndex.push_back(index);

						count++;

					}

				}

			}

			label++;

		}

	}

	int tempNeiIndex[2] = {-1};

	for (i = 0; i < elementNumber; i++)
	{

		int kSize = 0;

		for (j = 0; j < elementArray[i].numDirectNei; j++)
		{

			if (tempCluster[i][0] != tempCluster[elementArray[i].directNei[j]][0])
			{

				tempNeiIndex[kSize] = elementArray[i].directNei[j];

				kSize++;
			}

		}

		if (kSize > 1 && tempCluster[tempNeiIndex[0]][0] == tempCluster[tempNeiIndex[1]][0])
		{

			elementArray[i].indexCluster = elementArray[tempNeiIndex[0]].indexCluster;

			tempCluster[i][0] = tempCluster[tempNeiIndex[0]][0];
			tempCluster[i][1] = tempCluster[tempNeiIndex[0]][1];

		}

	}


	if (CheckBoundaryConnectivity() == false)
	{

		EnforceBoundaryConnectivity();

		//return false;

	}
	else
	{

		cout<<"Label connectivity is now correct!"<<endl;

		return true;

	}


}

bool PolycubeConstruction::CheckBoundaryConnectivity()
{

	int i, j, k, index;
	int label = 0;

	vector<vector<int> > tempCluster(elementNumber);

	for (i = 0; i < elementNumber; i++)
	{

		tempCluster[i].push_back(-1);
		tempCluster[i].push_back(elementArray[i].indexCluster);

	}

	for (i = 0; i < elementNumber; i++)
	{

		if (0 > tempCluster[i][0])
		{

			tempCluster[i][0] = label;

			int count = 1;

			vector<int> tempIndex;

			tempIndex.push_back(i);

			for (int c = 0; c < count; c++)
			{

				for (int n = 0; n < elementArray[tempIndex[c]].numDirectNei; n++)
				{

					index = elementArray[tempIndex[c]].directNei[n];

					if (0 > tempCluster[index][0] && elementArray[i].indexCluster == elementArray[index].indexCluster)
					{

						tempCluster[index][0] = label;

						tempIndex.push_back(index);

						count++;

					}

				}

			}

			label++;

		}

	}

	int tempNeiIndex[2] = {-1};

	for (i = 0; i < elementNumber; i++)
	{

		int kSize = 0;

		for (j = 0; j < elementArray[i].numDirectNei; j++)
		{

			if (tempCluster[i][0] != tempCluster[elementArray[i].directNei[j]][0])
			{

				tempNeiIndex[kSize] = elementArray[i].directNei[j];

				kSize++;

			}

		}

		if (kSize > 1 && tempCluster[tempNeiIndex[0]][0] == tempCluster[tempNeiIndex[1]][0])
		{

			return false;

		}

	}

	return true;

}

bool PolycubeConstruction::InitializePolycube()
{
	
	
	string tempString;
	
	int labelMAX = -1;
	int label;
	int i, index;
	int j, k;
		

	cout<<"**********************************************************"<<endl;
	cout<<"Initializing the polycube structure!"<<endl;

	label = 0;

	vector<vector<int> > tempCluster(elementNumber);

	for (i = 0; i < elementNumber; i++)
	{

		tempCluster[i].push_back(-1);
		tempCluster[i].push_back(elementArray[i].indexCluster);

	}

	for (i = 0; i < elementNumber; i++)
	{

		if (0 > tempCluster[i][0])
		{

			tempCluster[i][0] = label;

			int count = 1;

			vector<int> tempIndex;

			tempIndex.push_back(i);
			for (int c = 0; c < count; c++)
			{

				for (int n = 0; n < elementArray[tempIndex[c]].numDirectNei; n++)
				{

					index = elementArray[tempIndex[c]].directNei[n];

					if (0 > tempCluster[index][0] && elementArray[i].indexCluster == elementArray[index].indexCluster
						&& elementArray[i].indexPart == elementArray[index].indexPart) //For bifurcation
					{

						tempCluster[index][0] = label;

						tempIndex.push_back(index);

						count++;

					}

				}

			}

			label++;

		}

	}

	NUM_POLYCUBE_PATCH = label; // total number of separated polycube patches

	polycubePatch.resize(NUM_POLYCUBE_PATCH);
	for (i = 0; i < NUM_POLYCUBE_PATCH; i++)
	{

		polycubePatch[i].index = i;
		polycubePatch[i].numCorner = 0;
		polycubePatch[i].numElements = 0;

	}

	for (i = 0; i < elementNumber; i++)
	{

		index = tempCluster[i][0];
		polycubePatch[index].element.push_back(i);
		polycubePatch[index].numElements++;
		///////////////////////////////////////////////////
		////assign each element to one patch, 12/18/2014
		elementArray[i].indexPatch = index;
		///////////////////////////////////////////////////

	}

	for (i = 0; i < NUM_POLYCUBE_PATCH; i++)
	{

		index = polycubePatch[i].element[0];
		polycubePatch[i].indexCluster = tempCluster[index][1];

	}

	//OutputPatchesVTKBifurcation("test.vtk");

		
	
		

	//getchar();


	//cout << "=====1";
	SearchCornerandBoundary();
	//cout << "=====1";
	//ModifyWrongBoundaryElements();
	//cout << "=====1";
	/*if (CURVE_SMOOTH == 1)
	{
		
		SmoothBoundaryCurve();
		if (CURVE_FITTING == 1)
		{
			CurveFittingBoundaryCurve();
		}

		OutputPatchesVTK("smooth_boundary_result.vtk");

		string tempString;
		tempString = inputName + "_indexPatch_writeFitting.k";
		WriteKFileIndexPatch(tempString.c_str());
	}*/




	//ParametricMappingCorner();
	//CreateInitialPolycube();
	//ParametricMapping();
	
	return true;

}

bool PolycubeConstruction::CornerMappingMaxMin(Polycube &currentPolycubePatch)
{

	int indexCluster;
	int index;
	int i;
	double minMax = 0.f;

	indexCluster = currentPolycubePatch.indexCluster;

	switch (indexCluster)
	{
	case 0:

		minMax = -1000000.f;
		for (i = 0; i < currentPolycubePatch.numCorner; i++)
		{

			index = currentPolycubePatch.cornerPoint[i];
			if (minMax < polycubePara->vertex[index][0])
			{

				minMax = polycubePara->vertex[index][0];

			}

		}

		for (i = 0; i < currentPolycubePatch.numCorner; i++)
		{

			index = currentPolycubePatch.cornerPoint[i];
			polycubePara->vertex[index][0] = minMax;

		}

		break;

	case 1:

		minMax = 1000000.f;
		for (i = 0; i < currentPolycubePatch.numCorner; i++)
		{

			index = currentPolycubePatch.cornerPoint[i];
			if (minMax > polycubePara->vertex[index][0])
			{

				minMax = polycubePara->vertex[index][0];

			}

		}

		for (i = 0; i < currentPolycubePatch.numCorner; i++)
		{

			index = currentPolycubePatch.cornerPoint[i];
			polycubePara->vertex[index][0] = minMax;

		}

		break;

	case 2:

		minMax = -1000000.f;
		for (i = 0; i < currentPolycubePatch.numCorner; i++)
		{

			index = currentPolycubePatch.cornerPoint[i];
			if (minMax < polycubePara->vertex[index][1])
			{

				minMax = polycubePara->vertex[index][1];

			}

		}

		for (i = 0; i < currentPolycubePatch.numCorner; i++)
		{

			index = currentPolycubePatch.cornerPoint[i];
			polycubePara->vertex[index][1] = minMax;

		}

		break;

	case 3:

		minMax = 1000000.f;
		for (i = 0; i < currentPolycubePatch.numCorner; i++)
		{

			index = currentPolycubePatch.cornerPoint[i];
			if (minMax > polycubePara->vertex[index][1])
			{

				minMax = polycubePara->vertex[index][1];

			}

		}

		for (i = 0; i < currentPolycubePatch.numCorner; i++)
		{

			index = currentPolycubePatch.cornerPoint[i];
			polycubePara->vertex[index][1] = minMax;

		}

		break;

	case 4:

		minMax = -1000000.f;
		for (i = 0; i < currentPolycubePatch.numCorner; i++)
		{

			index = currentPolycubePatch.cornerPoint[i];
			if (minMax < polycubePara->vertex[index][2])
			{

				minMax = polycubePara->vertex[index][2];

			}

		}

		for (i = 0; i < currentPolycubePatch.numCorner; i++)
		{

			index = currentPolycubePatch.cornerPoint[i];
			polycubePara->vertex[index][2] = minMax;

		}

		break;

	case 5:

		minMax = 1000000.f;
		for (i = 0; i < currentPolycubePatch.numCorner; i++)
		{

			index = currentPolycubePatch.cornerPoint[i];
			if (minMax > polycubePara->vertex[index][2])
			{

				minMax = polycubePara->vertex[index][2];

			}

		}

		for (i = 0; i < currentPolycubePatch.numCorner; i++)
		{

			index = currentPolycubePatch.cornerPoint[i];
			polycubePara->vertex[index][2] = minMax;

		}

		break;

	default:

		cout << "Error in CornerMappingMaxMin()!!!" << endl;
		break;
	}

	return true;
}

bool PolycubeConstruction::ParametricMappingCorner()
{
	polycubePara = NULL;

	CopyMesh(polycubePara);
	int i, j;

	for (j = 0; j < 1; j++)
	{
		for (i = 0; i < NUM_POLYCUBE_PATCH; i++)
		{

			CornerMappingMaxMin(polycubePatch[i]);

		}
	}
	

	string tempString;
	tempString = inputName + "_mapping_write_test.k";
	WriteKFileMapping(tempString.c_str());

	//int index;

	//double maxXYZ[3] = { -1000000.0 }, minXYZ[3] = { 1000000.0 };

	//double maxDimension = -1000000.0f;

	//for (i = 0; i < vertexNumber; i++)
	//{

	//	if (IsCornerPoint(i))
	//	{

	//		cornerPoints.push_back(i);

	//	}

	//}

	//numberCornerPoints = cornerPoints.size();

	//for (i = 0; i < numberCornerPoints; i++)
	//{

	//	index = cornerPoints[i];

	//	for (j = 0; j < 3; j++)
	//	{

	//		if (maxXYZ[j] < polycubePara->vertex[index][j])
	//		{

	//			maxXYZ[j] = polycubePara->vertex[index][j];

	//		}

	//		if (minXYZ[j] > polycubePara->vertex[index][j])
	//		{

	//			minXYZ[j] = polycubePara->vertex[index][j];

	//		}

	//	}

	//}

	//for (i = 0; i < 3; i++)
	//{

	//	if (maxDimension < (maxXYZ[i] - minXYZ[i]))
	//	{

	//		maxDimension = maxXYZ[i] - minXYZ[i];

	//	}

	//}

	//for (i = 0; i < numberCornerPoints; i++)
	//{

	//	index = cornerPoints[i];

	//	for (j = 0; j < 3; j++)
	//	{

	//		//int tempValue = (int) 0.f + (polycubePara->vertex[index][j]-minXYZ[j]) / (maxXYZ[j]-minXYZ[j]) * 64.0f + 0.5f;
	//		int tempValue = (int) 0.f + (polycubePara->vertex[index][j] - minXYZ[j]) / maxDimension * DOMAIN_SIZE + 0.5f;

	//		polycubePara->vertex[index][j] = tempValue;

	//		outputCornerPoints << polycubePara->vertex[index][j] << ", ";
	//	}

	//	outputCornerPoints << endl;

	//}

	return true;

}

bool PolycubeConstruction::SmoothBoundaryCurve()
{

	int i, j, k, p;
	int index;

	for (i = 0; i < NUM_POLYCUBE_PATCH; i++)
	{

		for (j = 0; j < polycubePatch[i].numCorner; j++)
		{

			for (k = 0; k < polycubePatch[i].boundaryEdge[j].size(); k++)
			{

				int count = 0;
				vector<int> neiElements;
				int indexCluster_one, indexCluster_two;

				indexCluster_one = polycubePatch[i].indexCluster;

				index = polycubePatch[i].boundaryEdge[j][k];

				if (IsCornerPoint(index) == true)
				{
					continue;
				}

				for (p = 0; p < elementValenceNumber[index]; p++)
				{
					int tempIndex = elementValence[index][p];

					if (elementArray[tempIndex].indexCluster == indexCluster_one)
					{
						count++;
						neiElements.push_back(tempIndex);
					}
					else
					{
						indexCluster_two = elementArray[tempIndex].indexCluster;
					}
				}

				if (count == 2)
				{

					int fourNodes[4];

					fourNodes[0] = polycubePatch[i].boundaryEdge[j][k];
					fourNodes[1] = polycubePatch[i].boundaryEdge[j][k+1];
					fourNodes[2] = polycubePatch[i].boundaryEdge[j][k-1];

					for (p = 0; p < 3; p++)
					{

						if (element[neiElements[0]][p] != fourNodes[0] && element[neiElements[0]][p] != fourNodes[1] && element[neiElements[0]][p] != fourNodes[2])
						{
							fourNodes[3] = element[neiElements[0]][p];
							break;
						}

					}

					double tempCenter[3];

					for (p = 0; p < 3; p++)
					{
						tempCenter[p] = 0.f;
					}

					for (p = 0; p < 3; p++)
					{
						tempCenter[p] += vertex[fourNodes[1]][p] + vertex[fourNodes[2]][p];
						tempCenter[p] += vertex[fourNodes[0]][p] + vertex[fourNodes[3]][p];
					}

					for (p = 0; p < 3; p++)
					{
						vertex[fourNodes[0]][p] = 0.25 * tempCenter[p]; // calculate the center point, assign it to the node!
					}

				}

			}

		}

	}

	return true;

}

bool PolycubeConstruction::CurveFittingBoundaryCurve()
{

	int i, j, k, l;
	int index, indexPrev;
	int edgeSize;

	double totalLength, tempLength;

	vector<int> edgeValueVec;
	int edgeValue;
	bool used_edge;

	for (i = 0; i < NUM_POLYCUBE_PATCH; i++)
	{

		for (j = 0; j < polycubePatch[i].boundaryEdge.size(); j++)
		{

			edgeSize = polycubePatch[i].boundaryEdge[j].size();

			//////////////////////////////////////////////////////////////
			//Check if used already
			edgeValue = 0;
			used_edge = false;
			for (k = 0; k < edgeSize; k++)
			{
				edgeValue += polycubePatch[i].boundaryEdge[j][k];
			}

			for (k = 0; k < edgeValueVec.size(); k++)
			{

				if (edgeValue == edgeValueVec[k])
				{
					used_edge = true;
					break;
				}

			}

			if (used_edge == true)
			{
				continue;
			}
			else
			{
				edgeValueVec.push_back(edgeValue);
			}
			//////////////////////////////////////////////////////////////

			int N = edgeSize-1;
			vector<double> lengthRatio(edgeSize, 0.0);
			//XYZ inputNodes[N+1];

			vector<XYZ> inputNodes;
			inputNodes.resize(N+1);
			vector<XYZ> outputNodes(edgeSize);

			const int RESOLUTION = 1000;

			vector<XYZ> outputRes(RESOLUTION);

			for (k = 0; k < edgeSize; k++)
			{
				index = polycubePatch[i].boundaryEdge[j][k];
				inputNodes[k].x = vertex[index][0];
				inputNodes[k].y = vertex[index][1];
				inputNodes[k].z = vertex[index][2];
			}

			int T = 4;

			vector<int> knots;
			knots.resize(N+T+1);

			SplineKnots(knots, N, T);

			totalLength = 0.f;
			for (k = 1; k < edgeSize; k++)
			{

				index = polycubePatch[i].boundaryEdge[j][k];
				indexPrev = polycubePatch[i].boundaryEdge[j][k-1];

				tempLength = 0.f;

				for (l = 0; l < 3; l++)
				{
					tempLength += pow(vertex[index][l]-vertex[indexPrev][l], 2);
				}

				tempLength = sqrt(tempLength);
				totalLength += tempLength;

			}

			tempLength = 0.f;
			for (k = 1; k < edgeSize; k++)
			{
				double tempLengthOne = 0.f;
				index = polycubePatch[i].boundaryEdge[j][k];
				indexPrev = polycubePatch[i].boundaryEdge[j][k-1];

				for (l = 0; l < 3; l++)
				{
					tempLengthOne += pow(vertex[index][l]-vertex[indexPrev][l], 2);
				}

				tempLengthOne = sqrt(tempLengthOne);

				tempLength += tempLengthOne;

				lengthRatio[k] = tempLength / totalLength;
			}

			//SplineCurve(inputNodes, N, knots, T, outputNodes, lengthRatio);
			//for (k = 0; k < edgeSize; k++)
			//{
			//	index = polycubePatch[i].boundaryEdge[j][k];
			//	vertex[index][0] = outputNodes[k].x;
			//	vertex[index][1] = outputNodes[k].y;
			//	vertex[index][2] = outputNodes[k].z;
			//}
			
			SplineCurve(inputNodes, N, knots, T, outputRes, RESOLUTION);

			/*for (k = 0; k < edgeSize - 1; k++)
			{
				index = polycubePatch[i].boundaryEdge[j][k];

				int tempInt = (int) (lengthRatio[k] * RESOLUTION);

				vertex[index][0] = outputRes[tempInt].x;
				vertex[index][1] = outputRes[tempInt].y;
				vertex[index][2] = outputRes[tempInt].z;
			}

			index = polycubePatch[i].boundaryEdge[j][edgeSize-1];

			vertex[index][0] = outputRes[RESOLUTION-1].x;
			vertex[index][1] = outputRes[RESOLUTION-1].y;
			vertex[index][2] = outputRes[RESOLUTION-1].z;*/

			double tempMin = 10000000.0f;
			int tempCount;
			for (k = 0; k < edgeSize; k++)
			{
				tempMin = 10000000.0f;
				index = polycubePatch[i].boundaryEdge[j][k];
				for (l = 0; l < RESOLUTION; l++)
				{

					double tempDouble = 0;
					tempDouble += pow(vertex[index][0]-outputRes[l].x, 2);
					tempDouble += pow(vertex[index][1]-outputRes[l].y, 2);
					tempDouble += pow(vertex[index][2]-outputRes[l].z, 2);

					tempDouble = sqrt(tempDouble);

					if (tempDouble < tempMin)
					{
						tempMin = tempDouble;
						tempCount = l;
					}

				}

				vertex[index][0] = outputRes[tempCount].x;
				vertex[index][1] = outputRes[tempCount].y;
				vertex[index][2] = outputRes[tempCount].z;

			}

		}

	}

	return true;

}

void PolycubeConstruction::SplineKnots(vector<int> &u, int n, int t)
{

	int j;

	for (j = 0; j <= n+t; j++)
	{

		if (j < t)
		{
			u[j] = 0;
		}
		else if (j <= n)
		{
			u[j] = j - t + 1;
		}
		else if (j > n)
		{
			u[j] = n - t + 2;
		}

	}

}

void PolycubeConstruction::SplineCurve(const vector<XYZ> &inp, int n, const vector<int> &knots, int t, vector<XYZ> &outp, const vector<double> &ratio)
{

	int i;
	int numNodes;
	XYZ tempOutp;

	numNodes = ratio.size();

	double tempDouble;

	//for (i = 0; i < numNodes; i++)
	//{
	//	tempDouble = (n - t + 2) * ratio[i];
	//	SplinePoint(knots, n, t, tempDouble, inp, tempOutp);
	//	outp[i] = tempOutp;
	//}

	for (i = 0; i < numNodes-1; i++)
	{
		tempDouble = (n - t + 2) * ratio[i];
		SplinePoint(knots, n, t, tempDouble, inp, tempOutp);
		outp[i] = tempOutp;
	}

	outp[numNodes-1] = inp[n];
}

void PolycubeConstruction::SplineCurve(const vector<XYZ> &inp, int n, const vector<int> &knots, int t, vector<XYZ> &outp, int res)
{

	int i; 

	double interval, increment;

	interval = 0.f;
	increment = (n - t + 2) / (double)(res-1);

	for (i = 0; i < res-1; i++)
	{

		SplinePoint(knots, n, t, interval, inp, outp[i]);

		interval += increment;
	}

	outp[res-1] = inp[n];

}

void PolycubeConstruction::SplinePoint(const vector<int> &u, int n, int t, double v, const vector<XYZ> &control, XYZ &output)
{

	int k;
	double b;

	output.x = 0.f;
	output.y = 0.f;
	output.z = 0.f;

	for (k = 0; k <= n; k++)
	{
		b = SplineBlend(k, t, u, v);

		output.x += control[k].x * b;
		output.y += control[k].y * b;
		output.z += control[k].z * b;
	}

}

//void CVTBasedPolycube::SplinePoint(const vector<int> &u, int n, int t, double v, const vector<XYZ> &control, vector<XYZ> &output)
//{
//
//	int k;
//	double b;
//
//	output.x = 0.f;
//	output.y = 0.f;
//	output.z = 0.f;
//
//	for (k = 0; k <= n; k++)
//	{
//		b = SplineBlend(k, t, u, v);
//
//		output.x += control[k].x * b;
//		output.y += control[k].y * b;
//		output.z += control[k].z * b;
//	}
//
//}

double PolycubeConstruction::SplineBlend(int k, int t, const vector<int> &u, double v)
{

	double value;

	if (t == 1) 
	{
		if ((u[k] <= v) && (v < u[k+1]))
			value = 1;
		else
			value = 0;
	} 
	else 
	{
		if ((u[k+t-1] == u[k]) && (u[k+t] == u[k+1]))
			value = 0;
		else if (u[k+t-1] == u[k]) 
			value = (u[k+t] - v) / (u[k+t] - u[k+1]) * SplineBlend(k+1,t-1,u,v);
		else if (u[k+t] == u[k+1])
			value = (v - u[k]) / (u[k+t-1] - u[k]) * SplineBlend(k,t-1,u,v);
		else
			value = (v - u[k]) / (u[k+t-1] - u[k]) * SplineBlend(k,t-1,u,v) + 
			(u[k+t] - v) / (u[k+t] - u[k+1]) * SplineBlend(k+1,t-1,u,v);
	}

	return(value);

}

bool PolycubeConstruction::SearchCornerandBoundary()
{

	InitiateElementValence();

	int i, j;
	int index;
	int vertexCCW;

	//vector<int> cornerPoints; //It is defined as a variable of the class

	for (i = 0; i < vertexNumber; i++)
	{

		if (IsCornerPoint(i))
		{

			cornerPoints.push_back(i);

		}

	}
	/*for (int loopi = 0; loopi < cornerPoints.size(); loopi++)
	{
		cout << cornerPoints[loopi] << ", ";
	}
*/

	string file_temp = "test_tri.vtk";
	//Write(file_temp.c_str());
	//cout << NUM_POLYCUBE_PATCH << endl;
	numberCornerPoints = cornerPoints.size(); //04142015

	int firstCorner;
	int tempIndex;

	for (i = 0; i < NUM_POLYCUBE_PATCH; i++)
	{

		vector<int> tempCornerPoints;

		for (j = 0; j < polycubePatch[i].numElements; j++)
		{

			index = polycubePatch[i].element[j];

			for (int k = 0; k < 3; k++)
			{

				tempIndex = element[index][k];

				if (IsCornerPoint(tempIndex))
				{

					bool kCount = true;

					for (int kk = 0; kk < tempCornerPoints.size(); kk++)
					{

						if (tempIndex == tempCornerPoints[kk])
						{

							kCount = false;

						}

					}

					if (kCount == true)
					{

						tempCornerPoints.push_back(tempIndex);

					}

				}

			}

		}

		polycubePatch[i].numCorner = tempCornerPoints.size();

		//polycubePatch[i].numCorner = 1;

		firstCorner = FindOneCornerPoint(polycubePatch[i]);

		polycubePatch[i].cornerPoint.push_back(firstCorner);

		int count = 1;

		vector<int> tempNodes;

		tempNodes.push_back(firstCorner);

		for (int c = 0; c < count; c++)
		{

			for (int n = 0; n < elementValenceNumber[tempNodes[c]]; n++)
			{

				index = elementValence[tempNodes[c]][n];

				vertexCCW = SearchCCWVertex(tempNodes[c], index);

				//if (elementArray[index].indexCluster == polycubePatch[i].indexCluster && IsBoundaryPoint(vertexCCW) && IsBoundaryEdge(tempNodes[c], vertexCCW))
				if (elementArray[index].indexPatch == polycubePatch[i].index && IsBoundaryPoint(vertexCCW) && IsBoundaryEdge(tempNodes[c], vertexCCW))
				{

					tempNodes.push_back(vertexCCW);

					count++;

				}
				//else if (elementArray[index].indexCluster == polycubePatch[i].indexCluster && IsCornerPoint(vertexCCW) && IsBoundaryEdge(tempNodes[c], vertexCCW))
				else if (elementArray[index].indexPatch == polycubePatch[i].index && IsCornerPoint(vertexCCW) && IsBoundaryEdge(tempNodes[c], vertexCCW))
				{

					if (vertexCCW == firstCorner)
					{

						tempNodes.push_back(vertexCCW);
						polycubePatch[i].boundaryEdge.push_back(tempNodes);

						tempNodes.clear();

						if (polycubePatch[i].cornerPoint.size() == polycubePatch[i].numCorner)
						{

							count = -1;

							break;
						}
						else
						{					

							for (int ii = 0; ii < tempCornerPoints.size(); ii++)
							{

								bool tempCount = false;

								tempIndex = tempCornerPoints[ii];

								for (int jj = 0; jj < polycubePatch[i].cornerPoint.size(); jj++)
								{

									if (tempIndex == polycubePatch[i].cornerPoint[jj])
									{

										tempCount = true;

									}

								}

								if (tempCount == false)
								{

									firstCorner = tempIndex;

									break;

								}

							}

							polycubePatch[i].cornerPoint.push_back(firstCorner);

							tempNodes.push_back(firstCorner);

							count = 1;

							c = -1;

							break;

						}

					}
					else
					{

						tempNodes.push_back(vertexCCW);

						polycubePatch[i].cornerPoint.push_back(vertexCCW);

						//polycubePatch[i].numCorner++;

						polycubePatch[i].boundaryEdge.push_back(tempNodes);

						tempNodes.clear();

						tempNodes.push_back(vertexCCW);

						count = 1;

						c = -1; // Very important

						break;

					}

				}

			}

		}

	}

	return true;

}

bool PolycubeConstruction::IsCornerPoint(int vertexID)
{

	//int i, j, index, indexCluster;
	//int kSize = 0;

	////bool counted = true;

	//vector<int> neiCluster;

	//for (i = 0; i < elementValenceNumber[vertexID]; i++)
	//{

	//	bool counted = true;

	//	index = elementValence[vertexID][i];

	//	indexCluster = elementArray[index].indexCluster;

	//	for (j = 0; j < neiCluster.size(); j++)
	//	{

	//		if (neiCluster[j] == indexCluster)
	//		{

	//			counted = false;

	//		}

	//	}

	//	if (counted == true)
	//	{

	//		neiCluster.push_back(indexCluster);

	//		kSize++;

	//	}

	//}

	//if (kSize == 3)
	//{

	//	return true;

	//}
	//else
	//{

	//	return false;

	//}

	//For Bifurcation Cases, but may be also more general
	int i, j, index, indexPatch;
	int kSize = 0;

	vector<int> neiPatch;

	for (i = 0; i < elementValenceNumber[vertexID]; i++)
	{

		bool counted = true;

		index = elementValence[vertexID][i];

		indexPatch = elementArray[index].indexPatch;

		for (j = 0; j < neiPatch.size(); j++)
		{

			if (neiPatch[j] == indexPatch)
			{

				counted = false;

			}

		}

		if (counted == true)
		{

			neiPatch.push_back(indexPatch);

			kSize++;

		}

	}

	if (kSize >= 3)
	{

		return true;

	}
	else
	{

		return false;

	}

}

int PolycubeConstruction::FindOneCornerPoint(Polycube &currentPolycubePatch)
{

	int i, j;
	int index;
	int vertexID;

	for (i = 0; i < currentPolycubePatch.numElements; i++)
	{

		index = currentPolycubePatch.element[i];

		for (j = 0; j < 3; j++)
		{	

			vertexID = element[index][j];

			if (IsCornerPoint(vertexID))
			{

				return vertexID;

			}

		}

	}

}

bool PolycubeConstruction::IsBoundaryPoint(int vertexID)
{

	//int i, j, index, indexCluster;
	//int kSize = 0;

	////bool counted = true;

	//vector<int> neiCluster;

	//for (i = 0; i < elementValenceNumber[vertexID]; i++)
	//{

	//	bool counted = true;

	//	index = elementValence[vertexID][i];

	//	indexCluster = elementArray[index].indexCluster;

	//	for (j = 0; j < neiCluster.size(); j++)
	//	{

	//		if (neiCluster[j] == indexCluster)
	//		{

	//			counted = false;

	//		}

	//	}

	//	if (counted == true)
	//	{

	//		neiCluster.push_back(indexCluster);

	//		kSize++;

	//	}

	//}

	//if (kSize == 2)
	//{

	//	return true;

	//}
	//else
	//{

	//	return false;

	//}

	//For Bifurcation Cases, but may be also more general 
	int i, j, index, indexPatch;
	int kSize = 0;

	vector<int> neiPatch;

	for (i = 0; i < elementValenceNumber[vertexID]; i++)
	{

		bool counted = true;

		index = elementValence[vertexID][i];

		indexPatch = elementArray[index].indexPatch;

		for (j = 0; j < neiPatch.size(); j++)
		{

			if (neiPatch[j] == indexPatch)
			{

				counted = false;

			}

		}

		if (counted == true)
		{

			neiPatch.push_back(indexPatch);

			kSize++;

		}

	}

	if (kSize == 2)
	{

		return true;

	}
	else
	{

		return false;

	}

}

bool PolycubeConstruction::IsBoundaryEdge(int vertexIDone, int vertexIDtwo)
{

	int twoElement[2] = {0};

	int i, j, index;

	int count = 0;

	for (i = 0; i < elementValenceNumber[vertexIDone]; i++)
	{

		index = elementValence[vertexIDone][i];

		for (j = 0; j < 3; j++)
		{

			if (element[index][j] == vertexIDtwo)
			{

				twoElement[count] = index;

				count++;

			}

		}

	}

	if (count != 2)
	{

		cout<<"Error in IsboundaryEdge!!!"<<endl;

		exit;

	}

	//if (elementArray[twoElement[0]].indexCluster != elementArray[twoElement[1]].indexCluster)
	//{

	//	return true;

	//}
	//else if (elementArray[twoElement[0]].indexPart != elementArray[twoElement[1]].indexPart)//For bifurcation cases
	//{
	//	return true;
	//}
	//else
	//{

	//	return false;

	//}


	//For Bifurcation Cases, but may be also more general 
	if (elementArray[twoElement[0]].indexPatch != elementArray[twoElement[1]].indexPatch)
	{
		return true;
	}
	else
	{
		return false;
	}

}


int PolycubeConstruction::SearchCCWVertex(int vertexID, int elementID)
{

	int i;
	int vertexCCW;
	int position = 0;

	for (i = 0; i < 3; i++)
	{

		if (element[elementID][i] == vertexID)
		{

			position = i;

			break;
		}

	}

	if (position == 2)
	{
		position = 0;

		vertexCCW = element[elementID][position];

		return vertexCCW;
	}
	else
	{

		position++;

		vertexCCW = element[elementID][position];

		return vertexCCW;

	}

}

bool PolycubeConstruction::ModifyWrongBoundaryElements()
{

	int i, j, k;
	int index;

	for (i = 0; i < NUM_POLYCUBE_PATCH; i++)
	{

		for (j = 0; j < polycubePatch[i].numCorner; j++)
		{

			int count = 0;
			vector<int> neiElements;

			index = polycubePatch[i].cornerPoint[j];

			for (k = 0; k < elementValenceNumber[index]; k++)
			{

				int tempIndex = elementValence[index][k];

				/*if (elementArray[tempIndex].indexCluster == polycubePatch[i].indexCluster)
				{

					count++;
					neiElements.push_back(tempIndex);

				}*/

				if (elementArray[tempIndex].indexPatch == polycubePatch[i].index)
				{

					count++;
					neiElements.push_back(tempIndex);

				}
			}

			if (count == 1)
			{

				if (neiElements.size() != 1)
				{
					cout <<"ERROR!!! neiElement.size != 1"<<endl;
				}

				EdgeFlipTwoElements(index, neiElements[0]);

			}

		}

	}

	//Reinitialization
	////////////////////////////////////////////////////////////////////////////
	
	///////////////////////////////////////////////////////////////////////////

	return true;

}

//This function will modify the original vertex and element of the mesh, be careful
bool PolycubeConstruction::EdgeFlipTwoElements(int vertexID, int elementID)
{

	int i, j, indexVertex;
	int indexElement;
	int fourVertices[4];
	int twoElements[2];

	fourVertices[0] = vertexID;

	twoElements[0] = elementID;

	for (i = 0; i < 3; i++)
	{

		if (element[elementID][i] == vertexID)
		{

			fourVertices[1] = element[elementID][(i+1)%3];

			fourVertices[2] = element[elementID][(i+2)%3];

			break;

		}

	}

	indexVertex = fourVertices[1];

	for (i = 0; i < elementValenceNumber[indexVertex]; i++)
	{

		bool count = true;

		for (j = 0; j < 3; j++)
		{

			if (element[elementValence[indexVertex][i]][j] == indexVertex)
			{

				if (fourVertices[2] == element[elementValence[indexVertex][i]][(j+3-1)%3])
				{

					fourVertices[3] = element[elementValence[indexVertex][i]][(j+1)%3];

					twoElements[1] = elementValence[indexVertex][i];

					count = false;

					break;

				}

			}

		}

		if (count == false)
		{

			break;

		}

	}

	indexElement = twoElements[0];

	element[indexElement][0] = fourVertices[0];
	element[indexElement][1] = fourVertices[1];
	element[indexElement][2] = fourVertices[3];

	indexElement = twoElements[1];

	element[indexElement][0] = fourVertices[3];
	element[indexElement][1] = fourVertices[2];
	element[indexElement][2] = fourVertices[0];

	return true;

}

bool PolycubeConstruction::ParametricMapping()
{
	
	polycubePara = NULL;

	CopyMesh(polycubePara);

	if (OUTPUT_CORNER_POINTS == 1)
	{
		///////////////////////////////////////////////////////////////OUT PUT CORNER POINTS
		outputCornerPoints.open(inputName + "_Output_CornerPoints.txt");

		for (int i = 0; i < numberCornerPoints; i++)
		{
			int index = cornerPoints[i];
			outputCornerPoints<<index<<" ";

			for (int j = 0; j < 3; j++)
			{
				outputCornerPoints << polycubePara->vertex[index][j]<<" ";
			}

			outputCornerPoints<<endl;
		}

		outputCornerPoints.close();
		//by yu
		outputCornerPoints.open(inputName + "_Output_CornerPoints_forK.txt");

		for (int i = 0; i < numberCornerPoints; i++)
		{
			int index = cornerPoints[i];
			outputCornerPoints<<index+1<<",";
		}

		outputCornerPoints.close();
		///////////////////////////////////////////////////////////////
	}

	if (READ_IN_MAPPING == 1)
	{
		polycubePara->Read("mapping_tri.raw");
	}
	else
	{

		ParametricMappingCornerByInput();
		ParametricMappingEdgeGeneral();

		ParametricMappingInterior();

		polycubePara->Write("mapping_tri.raw");
		polycubePara->Write("mapping_tri.inp");

		string tempString;
		tempString = inputName + "_mapping_write.k";
		WriteKFileMapping(tempString.c_str());

	}

	return true;

}

bool PolycubeConstruction::ParametricMappingCornerByInput()
{

	string tempName;
	string oneLine;
	tempName = inputName + "_Output_CornerPoints.txt";

	ifstream myFile(tempName);

	//int tempX, tempY, tempZ;
	double tempX, tempY, tempZ;

	int i = 0, j, index;

	//for (i = 0; i < vertexNumber; i++)
	//{

	//	if (IsCornerPoint(i))
	//	{

	//		cornerPoints.push_back(i);

	//	}

	//}

	//numberCornerPoints = cornerPoints.size();

	i = 0;

	if (myFile.is_open())
	{

		while (getline(myFile, oneLine))
		{

			istringstream st(oneLine);

			st >> j >> tempX >> tempY >> tempZ; //>>, not <<

			index = cornerPoints[i];

			index = j; // Added on 08/11/2015

			if (j != index)
			{
				cout <<"Error in ParametricMappingCornerByInput()!!"<<endl;
			}

			polycubePara->vertex[index][0] = tempX;
			polycubePara->vertex[index][1] = tempY;
			polycubePara->vertex[index][2] = tempZ;

			//printf("%lf,%lf,%lf\n",polycubePara->vertex[index][0],polycubePara->vertex[index][1],polycubePara->vertex[index][2]);

			i++;
		}

		myFile.close();

	}
	else
	{
		cout<<"Unable to open file!";
	}

	return true;

}

bool PolycubeConstruction::ParametricMappingEdge()
{

	int i, j, k, l;
	int index, indexPrev, indexStart, indexEnd;
	int edgeSize;
	int direction;

	int tempIntStart, tempIntEnd;

	double totalLength, tempLength;

	for (i = 0; i < NUM_POLYCUBE_PATCH; i++)
	{

		for (j = 0; j < polycubePatch[i].boundaryEdge.size(); j++)
		{

			edgeSize = polycubePatch[i].boundaryEdge[j].size();

			indexStart = polycubePatch[i].boundaryEdge[j][0];
			indexEnd = polycubePatch[i].boundaryEdge[j][edgeSize-1];

			/////////////////////////////////////
			//For DEBUG
			if (j == 3 && i == 2)
			{

				j = j;

			}
			////////////////////////////////////

			for (k = 0; k < 3; k++)
			{

				tempIntStart = (int) polycubePara->vertex[indexStart][k] + 0.5f;
				tempIntEnd = (int) polycubePara->vertex[indexEnd][k] + 0.5f;

				if (tempIntStart != tempIntEnd)
				{

					direction = k;

					break;

				}

			}

			totalLength = 0.f;

			for (k = 1; k < edgeSize-1; k++) // k < edgeSize is very important!
			{

				index = polycubePatch[i].boundaryEdge[j][k];
				indexPrev = polycubePatch[i].boundaryEdge[j][k-1];

				tempLength = 0.f;

				for (l = 0; l < 3; l++)
				{

					tempLength += pow(vertex[index][l]-vertex[indexPrev][l], 2);

				}

				tempLength = sqrt(tempLength);

				totalLength += tempLength;

				for (l = 0; l < 3; l++)
				{

					if (l != direction)
					{

						polycubePara->vertex[index][l] = polycubePara->vertex[indexStart][l];

					}
					else
					{

						polycubePara->vertex[index][l] = totalLength;

					}

				}

			}

			////////////////////////////////////////////////////////////
			//Final total Length
			index = polycubePatch[i].boundaryEdge[j][edgeSize-1];
			indexPrev = polycubePatch[i].boundaryEdge[j][edgeSize-2];
			tempLength = 0.f;

			for (l = 0; l < 3; l++)
			{

				tempLength += pow(vertex[index][l]-vertex[indexPrev][l], 2);

			}

			tempLength = sqrt(tempLength);

			totalLength += tempLength;
			///////////////////////////////////////////////////////////

			for (k = 1; k < edgeSize-1; k++)
			{

				index = polycubePatch[i].boundaryEdge[j][k];

				for (l = 0; l < 3; l++)
				{

					if (l == direction)
					{

						///////////////////////////////////////////////////////////////////////////////////////
						//////In fact, the following if else conditions are the same! 12/18/2014
						if (polycubePara->vertex[indexStart][l] < polycubePara->vertex[indexEnd][l])
						{

							polycubePara->vertex[index][l] = polycubePara->vertex[indexStart][l] + polycubePara->vertex[index][l] / totalLength * (polycubePara->vertex[indexEnd][l]-polycubePara->vertex[indexStart][l]);

						}
						else
						{

							polycubePara->vertex[index][l] = polycubePara->vertex[indexStart][l] - polycubePara->vertex[index][l] / totalLength * (polycubePara->vertex[indexStart][l]-polycubePara->vertex[indexEnd][l]);

						}

					}

				}

			}

		}

	}


	return true;

}

bool PolycubeConstruction::ParametricMappingEdgeGeneral()
{

	int i, j, k, l;
	int index, indexPrev, indexStart, indexEnd;
	int edgeSize;
	double direction[3];

	int tempIntStart, tempIntEnd;

	double totalLength, tempLength, totalLengthPara;

	for (i = 0; i < NUM_POLYCUBE_PATCH; i++)
	{

		for (j = 0; j < polycubePatch[i].boundaryEdge.size(); j++)
		{

			edgeSize = polycubePatch[i].boundaryEdge[j].size();

			indexStart = polycubePatch[i].boundaryEdge[j][0];
			indexEnd = polycubePatch[i].boundaryEdge[j][edgeSize-1];

			totalLengthPara = 0.0f; // length in the parametric domain
			for (k = 0; k < 3; k++)
			{

				//tempIntStart = (int) polycubePara->vertex[indexStart][k] + 0.5f;
				//tempIntEnd = (int) polycubePara->vertex[indexEnd][k] + 0.5f;

				//if (tempIntStart != tempIntEnd)
				//{
				//	direction = k;
				//	break;
				//}

				direction[k] = polycubePara->vertex[indexEnd][k] - polycubePara->vertex[indexStart][k];
				totalLengthPara += pow(polycubePara->vertex[indexEnd][k] - polycubePara->vertex[indexStart][k], 2);
			}

			totalLengthPara = sqrt(totalLengthPara);
			Normalize(direction);

			totalLength = 0.f;

			for (k = 1; k < edgeSize-1; k++) // k < edgeSize is very important!
			{

				index = polycubePatch[i].boundaryEdge[j][k];
				indexPrev = polycubePatch[i].boundaryEdge[j][k-1];

				tempLength = 0.f;

				for (l = 0; l < 3; l++)
				{

					tempLength += pow(vertex[index][l]-vertex[indexPrev][l], 2);

				}

				tempLength = sqrt(tempLength);

				totalLength += tempLength;

				for (l = 0; l < 3; l++)
				{

					/*if (l != direction)
					{

						polycubePara->vertex[index][l] = polycubePara->vertex[indexStart][l];

					}
					else
					{

						polycubePara->vertex[index][l] = totalLength;

					}*/

					polycubePara->vertex[index][l] = totalLength;

				}

			}

			////////////////////////////////////////////////////////////
			//Final total Length
			index = polycubePatch[i].boundaryEdge[j][edgeSize-1];
			indexPrev = polycubePatch[i].boundaryEdge[j][edgeSize-2];
			tempLength = 0.f;

			for (l = 0; l < 3; l++)
			{

				tempLength += pow(vertex[index][l]-vertex[indexPrev][l], 2);

			}

			tempLength = sqrt(tempLength);

			totalLength += tempLength;
			///////////////////////////////////////////////////////////

			for (k = 1; k < edgeSize-1; k++)
			{

				index = polycubePatch[i].boundaryEdge[j][k];

				for (l = 0; l < 3; l++)
				{

					//if (l == direction)
					//{

					//	///////////////////////////////////////////////////////////////////////////////////////
					//	//////In fact, the following if else conditions are the same! 12/18/2014
					//	if (polycubePara->vertex[indexStart][l] < polycubePara->vertex[indexEnd][l])
					//	{

					//		polycubePara->vertex[index][l] = polycubePara->vertex[indexStart][l] + polycubePara->vertex[index][l] / totalLength * (polycubePara->vertex[indexEnd][l]-polycubePara->vertex[indexStart][l]);

					//	}
					//	else
					//	{

					//		polycubePara->vertex[index][l] = polycubePara->vertex[indexStart][l] - polycubePara->vertex[index][l] / totalLength * (polycubePara->vertex[indexStart][l]-polycubePara->vertex[indexEnd][l]);

					//	}

					//}

					polycubePara->vertex[index][l] = polycubePara->vertex[indexStart][l] + polycubePara->vertex[index][l] / totalLength * totalLengthPara * direction[l];

				}

			}

		}

	}


	return true;

}

bool PolycubeConstruction::ParametricMappingInterior()
{

	int i, j, k, l;
	int index, indexCluster;

	for (i = 0; i < NUM_POLYCUBE_PATCH; i++)
	{

		polycubePatch[i].numIntVert = 0;

		vector<int> InteriorList;

		for (j = 0; j < polycubePatch[i].numElements; j++)
		{

			index = polycubePatch[i].element[j];

			for (k = 0; k < 3; k++)
			{

				if (IsInteriorPoint(element[index][k]))
				{

					bool count = true;

					for (l = 0; l < InteriorList.size(); l++)
					{

						if (InteriorList[l] == element[index][k])
						{

							count = false;
							break;

						}

					}

					if (count == true)
					{

						InteriorList.push_back(element[index][k]);
						polycubePatch[i].numIntVert++;

					}

				}

			}

		}

		polycubePatch[i].InterVert = InteriorList;

	}

	for (i = 0; i < NUM_POLYCUBE_PATCH; i++)
	{
		//i=47;
		InteriorMappingOnePatch(polycubePatch[i]);
		cout<<"Mapping patch ID: "<<i<<endl;

	}

	return true;

}

bool PolycubeConstruction::IsInteriorPoint(int vertexID)
{

	//int i, j, index, indexCluster;
	//int kSize = 0;

	////bool counted = true;

	//vector<int> neiCluster;

	//for (i = 0; i < elementValenceNumber[vertexID]; i++)
	//{

	//	bool counted = true;

	//	index = elementValence[vertexID][i];

	//	indexCluster = elementArray[index].indexCluster;

	//	for (j = 0; j < neiCluster.size(); j++)
	//	{

	//		if (neiCluster[j] == indexCluster)
	//		{

	//			counted = false;

	//		}

	//	}

	//	if (counted == true)
	//	{

	//		neiCluster.push_back(indexCluster);

	//		kSize++;

	//	}

	//}

	//if (kSize == 1)
	//{

	//	return true;

	//}
	//else
	//{

	//	return false;

	//}

	//For Bifurcation Cases, but may be also more general
	int i, j, index, indexPatch;
	int kSize = 0;

	vector<int> neiPatch;

	for (i = 0; i < elementValenceNumber[vertexID]; i++)
	{

		bool counted = true;

		index = elementValence[vertexID][i];

		indexPatch = elementArray[index].indexPatch;

		for (j = 0; j < neiPatch.size(); j++)
		{

			if (neiPatch[j] == indexPatch)
			{

				counted = false;

			}

		}

		if (counted == true)
		{

			neiPatch.push_back(indexPatch);

			kSize++;

		}

	}

	if (kSize == 1)
	{

		return true;

	}
	else
	{

		return false;

	}

}

bool PolycubeConstruction::InteriorMappingOnePatch(Polycube &currentPolycubePatch)
{

	int i, j, k, p, q;
	int index, indexCluster, iConst, nSize;
	int nCount, nIndex, nNeighbor, iIndex[2];
	double weight, tempLength, vecta[2][3];

	

	vector<int> InteriorList;
	InteriorList.clear();

	nCount = currentPolycubePatch.numIntVert;

	//MatrixXd mW = MatrixXd::Zero(nCount, nCount);
	//MatrixXd vB = MatrixXd::Zero(nCount, 2);
	//VectorXd x_solve0 = VectorXd::Zero(nCount);
	//VectorXd x_solve1 = VectorXd::Zero(nCount);

	Eigen::SparseMatrix<double, ColMajor> mW(nCount, nCount);
	Eigen::VectorXd vB0(nCount), vB1(nCount), x_solve0(nCount), x_solve1(nCount), vBy(nCount), x_solvey(nCount);
	mW.setZero();
	vB0.setZero();
	vB1.setZero();
	x_solve0.setZero();
	x_solve1.setZero();

	indexCluster = currentPolycubePatch.indexCluster;  
	double theta=1.0;
	//iConst = patch_const[indexCluster];
	for (i = 0; i < 3; i++)
	{
		bool planarDir = true;
		for (j = 0; j < currentPolycubePatch.numCorner; j++)
		{
			index = currentPolycubePatch.cornerPoint[j];
			//printf("xyz:=%d, index:= %d  %lf %lf\n", i,index, polycubePara->vertex[index][i],polycubePara->vertex[currentPolycubePatch.cornerPoint[0]][i]);

			if (fabs(polycubePara->vertex[index][i]-polycubePara->vertex[currentPolycubePatch.cornerPoint[0]][i]) > EPSILON)
			{
				planarDir = false;
			}
		}

		if (planarDir == true)
		{
			iConst = i; 
			theta=1.0;
			printf("it is a plane\n");
			break;
		}


	    //by yu
		if (planarDir == false)    
		{
			iConst = 1;
			theta=1.0/sqrt(2.0);
		}
		//by yu

	}


	


	sort(currentPolycubePatch.InterVert.begin(), currentPolycubePatch.InterVert.end());

	InteriorList = currentPolycubePatch.InterVert;

	nSize = currentPolycubePatch.numIntVert;

	for (j = 0; j < nSize; j++)
	{

		nNeighbor = elementValenceNumber[InteriorList[j]];

		for (k = 0; k < nNeighbor; k++)
		{

			for (p = 0; p < 3; p++)
			{
				//printf("%d %d %d\n", InteriorList[j],elementValence[InteriorList[j]][k],element[elementValence[InteriorList[j]][k]][p]);
				////printf("%d %d %d", InteriorList[j+1],elementValence[InteriorList[j]][k+3],element[elementValence[InteriorList[j]][k]][p+1]);
				////printf("%d %d %d\n", InteriorList[j+1],elementValence[InteriorList[j]][k+4],element[elementValence[InteriorList[j]][k]][p+2]);
				if (element[elementValence[InteriorList[j]][k]][p] == InteriorList[j])
				{

					nIndex = element[elementValence[InteriorList[j]][k]][(p+1)%3];

					iIndex[0] = element[elementValence[InteriorList[j]][k]][(p+2)%3];

					break;

				}

			}

			for (i = 0; i < nNeighbor; i++)
			{
				
				bool count = true;

				for (p = 0; p < 3; p++)
				{
					//printf("%d %d %d\n", InteriorList[j],elementValence[InteriorList[j]][i],element[elementValence[InteriorList[j]][i]][p]);
					if (element[elementValence[InteriorList[j]][i]][p] == InteriorList[j])
					{

						if (nIndex == element[elementValence[InteriorList[j]][i]][(p+3-1)%3])
						{

							iIndex[1] = element[elementValence[InteriorList[j]][i]][(p+1)%3];
							count = false;

							break;

						}

					}
					

				}
				/*printf("%d %d %d\n", InteriorList[j],elementValence[InteriorList[j]][i],element[elementValence[InteriorList[j]][i]][p]);*/
				if (count == false)
				{

					break;

				}

			}

			weight = 0;

			for (p = 0; p < 2; p++)
			{

				for (q = 0; q < 3; q++)
				{

					vecta[0][q] = vertex[InteriorList[j]][q] - vertex[iIndex[p]][q];

					vecta[1][q] = vertex[nIndex][q] - vertex[iIndex[p]][q];

				}

				tempLength = sqrt(pow(vecta[0][1]*vecta[1][2]-vecta[0][2]*vecta[1][1],2) + 
					pow(vecta[0][2]*vecta[1][0]-vecta[0][0]*vecta[1][2],2) + pow(vecta[0][0]*vecta[1][1]-vecta[0][1]*vecta[1][0],2));

				tempLength = (vecta[0][0]*vecta[1][0] + vecta[0][1]*vecta[1][1] + vecta[0][2]*vecta[1][2])/tempLength;

				weight += tempLength;

			}

			//mW(j, j) -= weight;
			mW.coeffRef(j, j) -= weight;

			if (IsInteriorPoint(nIndex))
			{

				//double distance; // int or double???
				int distance;

				for (i = 0; i < InteriorList.size(); i++)
				{

					if (InteriorList[i] == nIndex)
					{

						distance = i;
						break;

					}

				}

				//mW(j, distance) += weight;
				mW.coeffRef(j, distance) += weight;

			}
			else
			{

				/*vB(j, 0) -= weight * polycubePara->vertex[nIndex][(iConst+1)%3];
				vB(j, 1) -= weight * polycubePara->vertex[nIndex][(iConst+2)%3];*/

				vB0[j] -= weight * polycubePara->vertex[nIndex][(iConst+1)%3];
				vB1[j] -= weight * polycubePara->vertex[nIndex][(iConst+2)%3];
				if(fabs(1.0/sqrt(2.0)-theta)<EPSILON)
				{
					vBy[j] -= weight * polycubePara->vertex[nIndex][(iConst)%3];
				}
				

			}

		}

	}

	//Eigen::FullPivLU<MatrixXd> lu(mW);
	//x_solve0 = lu.solve(vB.topLeftCorner(nCount, 1));
	//x_solve1 = lu.solve(vB.topRightCorner(nCount, 1));

	//Eigen::SimplicialLLT<SparseMatrix<double> > lu;
	//lu.compute(mW);

	Eigen::ConjugateGradient<SparseMatrix<double> > lu;
	lu.compute(mW);

	x_solve0 = lu.solve(vB0);
	x_solve1 = lu.solve(vB1);
	if(fabs(1.0/sqrt(2.0)-theta)<EPSILON)
	{
		printf("calculate");
		x_solvey = lu.solve(vBy);
	}
	

	if (lu.info() != Eigen::Success) 
	{
		cout <<"solve surface error in InteriorMappingOnePatch()!\n";
		return false;
	}

	////////////////////////////////////////////////////////////////////////////////

	int oneCorner = currentPolycubePatch.cornerPoint[0];
	int secCorner = currentPolycubePatch.cornerPoint[2];

	for(j = 0; j < nSize; j++)
	{
	

		
		polycubePara->vertex[InteriorList[j]][(iConst+1)%3] = x_solve0(j);
		polycubePara->vertex[InteriorList[j]][(iConst+2)%3] = x_solve1(j);

		

		if(fabs(1.0/sqrt(2.0)-theta)>EPSILON)
         {
          polycubePara->vertex[InteriorList[j]][iConst] = polycubePara->vertex[oneCorner][iConst];
         } else
         {
			//polycubePara->vertex[InteriorList[j]][(iConst)%3] = x_solve0(j);
           //polycubePara->vertex[InteriorList[j]][(iConst+1)%3] = polycubePara->vertex[oneCorner][(iConst+1)%3];
			// polycubePara->vertex[InteriorList[j]][(iConst)%3] = polycubePara->vertex[oneCorner][(iConst)%3];
			 //printf("%d, %lf,%lf,%lf\n",InteriorList[j],x_solve0(j),x_solve1(j),x_solvey(j));
			 //printf("%d, %lf,%lf,%lf\n",InteriorList[j],vertex[InteriorList[j]][(iConst)%3],vertex[oneCorner][(iConst)%3]);
			 polycubePara->vertex[InteriorList[j]][(iConst)%3] = polycubePara->vertex[oneCorner][(iConst)%3]+(vertex[InteriorList[j]][(iConst)%3]-vertex[oneCorner][(iConst)%3])/fabs(vertex[oneCorner][(iConst)%3]-vertex[secCorner][(iConst)%3])*fabs(polycubePara->vertex[oneCorner][(iConst)%3]-polycubePara->vertex[secCorner][(iConst)%3]);
          //polycubePara->vertex[InteriorList[j]][iConst] =polycubePara->vertex[oneCorner][iConst]-fabs(polycubePara->vertex[InteriorList[j]][(iConst+2)%3]-polycubePara->vertex[oneCorner][(iConst+2)%3]);
         }


	}

	/////////////////////////////////////////////////////////////

	//MappingOnePatchWeighted(currentPolycubePatch);
	//MappingOnePatchPostProcessing(currentPolycubePatch);

	int stepNumber = 0;

	while(FlipCheckEachPatch(currentPolycubePatch) == false)
	{

		cout<<"Polycube patch "<<currentPolycubePatch.index<<" has flipped elements!"<<"Cluster No. "<<currentPolycubePatch.indexCluster<<endl;
		
		if (stepNumber < 0)
		{
			MappingOnePatchWeighted(currentPolycubePatch);
		}
		if(fabs(1.0/sqrt(2.0)-theta)>EPSILON)
		{
		MappingOnePatchPostProcessing(currentPolycubePatch);
		}
		stepNumber++;

		cout<<"Polycube patch "<<currentPolycubePatch.index<<" still has flipped elements after "<<stepNumber<<" iterations!"<<endl;

		if (stepNumber > 0)
		{
			
			cout<<"Polycube patch "<<currentPolycubePatch.index<<" still has flipped elements after 1000 iterations!"<<endl;
			break;

		}

	}	

	//MappingOnePatchWeighted(currentPolycubePatch);
	///////////////////////////////////////////////////////////

	return true;

}

bool PolycubeConstruction::FlipCheckEachPatch(Polycube &currentPolycubePatch)
{

	int elementID, indexCluster;

	bool trueFalse = true;

	int i;

	indexCluster = currentPolycubePatch.indexCluster;

	for (i = 0; i < currentPolycubePatch.numElements; i++)
	{

		elementID = currentPolycubePatch.element[i];

		trueFalse = SignedTriAreaElement(elementID, indexCluster);

		if (trueFalse == false)
		{

			return false;

		}

	}

	return true;

}

bool PolycubeConstruction::SignedTriAreaElement(int elementID, int indexCluster)
{

	double point_1[3], point_2[3], point_3[3];
	double result[3];
	int i;

	for (i = 0; i < 3; i++)
	{

		point_1[i] = polycubePara->vertex[element[elementID][0]][i];

		point_2[i] = polycubePara->vertex[element[elementID][1]][i];

		point_3[i] = polycubePara->vertex[element[elementID][2]][i];

	}

	CrossProduct(point_1, point_2, point_3, result);

	switch (indexCluster)
	{
	case 0:

		//return 0.5 * result[0];
		if (result[0] < 0.f)
		{
			return false;
		}

		break;

	case 1:

		//return 0.5 * result[0];
		if (result[0] > 0.f)
		{
			return false;
		}

		break;

	case 2:

		//return 0.5 * result[1];
		if (result[1] < 0.f)
		{
			return false;
		}

		break;

	case 3:

		//return 0.5 * result[1];
		if (result[1] > 0.f)
		{
			return false;
		}

		break;

	case 4:

		//return 0.5 * result[2];
		if (result[2] < 0.f)
		{
			return false;
		}

		break;

	case 5:

		//return 0.5 * result[2];
		if (result[2] > 0.f)
		{
			return false;
		}

		break;

	default:

		cout<<"Error in SignedTriAreaElement()!!!"<<endl;
		return false;
		break;
		exit;

	}

	return true;

}

bool PolycubeConstruction::MappingOnePatchWeighted(Polycube &currentPolycubePatch)
{

	int i, j, k, p, q;
	int index, indexCluster, iConst, nSize;
	int nCount, nIndex, nNeighbor, iIndex[2];
	double weight, tempLength, vecta[2][3];
	int twoElements[2];

	vector<int> InteriorList;
	InteriorList.clear();

	nCount = currentPolycubePatch.numIntVert;

	//MatrixXd mW = MatrixXd::Zero(nCount, nCount);
	//MatrixXd vB = MatrixXd::Zero(nCount, 2);
	//VectorXd x_solve0 = VectorXd::Zero(nCount);
	//VectorXd x_solve1 = VectorXd::Zero(nCount);

	Eigen::SparseMatrix<double, ColMajor> mW(nCount, nCount);
	Eigen::VectorXd vB0(nCount), vB1(nCount), x_solve0(nCount), x_solve1(nCount);
	mW.setZero();
	vB0.setZero();
	vB1.setZero();
	x_solve0.setZero();
	x_solve1.setZero();

	indexCluster = currentPolycubePatch.indexCluster;

	//iConst = patch_const[indexCluster];
	for (i = 0; i < 3; i++)
	{
		bool planarDir = true;
		for (j = 0; j < currentPolycubePatch.numCorner; j++)
		{
			index = currentPolycubePatch.cornerPoint[j];

			if (fabs(polycubePara->vertex[index][i]-polycubePara->vertex[currentPolycubePatch.cornerPoint[0]][i]) > EPSILON)
			{
				planarDir = false;
			}
		}

		if (planarDir == true)
		{
			iConst = i;
			break;
		}

	}

	sort(currentPolycubePatch.InterVert.begin(), currentPolycubePatch.InterVert.end());

	InteriorList = currentPolycubePatch.InterVert;

	nSize = currentPolycubePatch.numIntVert;

	for (j = 0; j < nSize; j++)
	{

		nNeighbor = elementValenceNumber[InteriorList[j]];

		for (k = 0; k < nNeighbor; k++)
		{

			for (p = 0; p < 3; p++)
			{

				if (element[elementValence[InteriorList[j]][k]][p] == InteriorList[j])
				{

					nIndex = element[elementValence[InteriorList[j]][k]][(p+1)%3];

					iIndex[0] = element[elementValence[InteriorList[j]][k]][(p+2)%3];

					twoElements[0] = elementValence[InteriorList[j]][k]; //for weighted mapping

					break;

				}

			}

			for (i = 0; i < nNeighbor; i++)
			{

				bool count = true;

				for (p = 0; p < 3; p++)
				{

					if (element[elementValence[InteriorList[j]][i]][p] == InteriorList[j])
					{

						if (nIndex == element[elementValence[InteriorList[j]][i]][(p+3-1)%3])
						{

							iIndex[1] = element[elementValence[InteriorList[j]][i]][(p+1)%3];

							twoElements[1] = elementValence[InteriorList[j]][i]; //for weighted mapping

							count = false;

							break;

						}

					}

				}

				if (count == false)
				{

					break;

				}

			}

			weight = 0;

			for (p = 0; p < 2; p++)
			{

				for (q = 0; q < 3; q++)
				{

					vecta[0][q] = polycubePara->vertex[InteriorList[j]][q] - polycubePara->vertex[iIndex[p]][q];

					vecta[1][q] = polycubePara->vertex[nIndex][q] - polycubePara->vertex[iIndex[p]][q];

				}

				tempLength = sqrt(pow(vecta[0][1]*vecta[1][2]-vecta[0][2]*vecta[1][1],2) + 
					pow(vecta[0][2]*vecta[1][0]-vecta[0][0]*vecta[1][2],2) + pow(vecta[0][0]*vecta[1][1]-vecta[0][1]*vecta[1][0],2));

				tempLength = (vecta[0][0]*vecta[1][0] + vecta[0][1]*vecta[1][1] + vecta[0][2]*vecta[1][2])/tempLength;

				//double triArea = abs(TriAreaElement(twoElements[p]));
				double triArea = fabs(TriAreaElement(twoElements[p]));
				
				weight += triArea * tempLength;

			}

			/*mW(j, j) -= weight;*/
			mW.coeffRef(j, j) -= weight;

			if (IsInteriorPoint(nIndex))
			{

				//double distance; // int or double???
				int distance;

				for (i = 0; i < InteriorList.size(); i++)
				{

					if (InteriorList[i] == nIndex)
					{

						distance = i;
						break;

					}

				}

				/*mW(j, distance) += weight;*/
				mW.coeffRef(j, distance) += weight;

			}
			else
			{

				/*vB(j, 0) -= weight * polycubePara->vertex[nIndex][(iConst+1)%3];
				vB(j, 1) -= weight * polycubePara->vertex[nIndex][(iConst+2)%3];*/
				vB0[j] -= weight * polycubePara->vertex[nIndex][(iConst+1)%3];
				vB1[j] -= weight * polycubePara->vertex[nIndex][(iConst+2)%3];

			}

		}

	}

	//Eigen::FullPivLU<MatrixXd> lu(mW);
	//x_solve0 = lu.solve(vB.topLeftCorner(nCount, 1));
	//x_solve1 = lu.solve(vB.topRightCorner(nCount, 1));

	//Eigen::SimplicialLLT<SparseMatrix<double> > lu;
	//lu.compute(mW);

	Eigen::ConjugateGradient<SparseMatrix<double> > lu;
	lu.compute(mW);

	x_solve0 = lu.solve(vB0);
	x_solve1 = lu.solve(vB1);

	if (lu.info() != Eigen::Success) 
	{
		cout <<"solve surface error in MappingOnePatchWeighted()!\n";
		return false;
	}

	int oneCorner = currentPolycubePatch.cornerPoint[0];

	for(j = 0; j < nSize; j++)
	{

		polycubePara->vertex[InteriorList[j]][iConst] = polycubePara->vertex[oneCorner][iConst];
		polycubePara->vertex[InteriorList[j]][(iConst+1)%3] = x_solve0(j);
		polycubePara->vertex[InteriorList[j]][(iConst+2)%3] = x_solve1(j);

	}

	return true;

}

double PolycubeConstruction::TriAreaElement(int elementID)
{

	double v0[3], v1[3], v2[3];

	int i;

	for (i = 0; i < 3; i++)
	{

		v0[i] = polycubePara->vertex[element[elementID][0]][i];

		v1[i] = polycubePara->vertex[element[elementID][1]][i];

		v2[i] = polycubePara->vertex[element[elementID][2]][i];

	}

	Vector3d n;

	Vector3d a(v1[0]-v0[0], v1[1]-v0[1], v1[2]-v0[2]);

	Vector3d b(v2[0]-v0[0], v2[1]-v0[1], v2[2]-v0[2]);

	n = a.cross(b);

	return 0.5 * n.norm();

}

bool PolycubeConstruction::MappingOnePatchPostProcessing(Polycube &currentPolycubePatch)
{

	int i, j, k, p, q;
	int index, indexCluster, iConst, nSize;
	int nCount, nIndex, nNeighbor, iIndex[2];
	double weight, tempLength, vecta[2][3];

	vector<int> InteriorList;
	InteriorList.clear();

	nCount = currentPolycubePatch.numIntVert;

	//MatrixXd mW = MatrixXd::Zero(nCount, nCount);
	//MatrixXd vB = MatrixXd::Zero(nCount, 2);
	//VectorXd x_solve0 = VectorXd::Zero(nCount);
	//VectorXd x_solve1 = VectorXd::Zero(nCount);

	Eigen::SparseMatrix<double, ColMajor> mW(nCount, nCount);
	Eigen::VectorXd vB0(nCount), vB1(nCount), x_solve0(nCount), x_solve1(nCount);
	mW.setZero();
	vB0.setZero();
	vB1.setZero();
	x_solve0.setZero();
	x_solve1.setZero();

	indexCluster = currentPolycubePatch.indexCluster;

	//iConst = patch_const[indexCluster];
	for (i = 0; i < 3; i++)
	{
		bool planarDir = true;
		for (j = 0; j < currentPolycubePatch.numCorner; j++)
		{
			index = currentPolycubePatch.cornerPoint[j];

			if (fabs(polycubePara->vertex[index][i]-polycubePara->vertex[currentPolycubePatch.cornerPoint[0]][i]) > EPSILON)
			{
				planarDir = false;
			}
		}

		if (planarDir == true)
		{
			iConst = i;
			break;
		}

	}

	sort(currentPolycubePatch.InterVert.begin(), currentPolycubePatch.InterVert.end());

	InteriorList = currentPolycubePatch.InterVert;

	nSize = currentPolycubePatch.numIntVert;

	for (j = 0; j < nSize; j++)
	{

		nNeighbor = elementValenceNumber[InteriorList[j]];

		for (k = 0; k < nNeighbor; k++)
		{

			for (p = 0; p < 3; p++)
			{

				if (element[elementValence[InteriorList[j]][k]][p] == InteriorList[j])
				{

					nIndex = element[elementValence[InteriorList[j]][k]][(p+1)%3];

					iIndex[0] = element[elementValence[InteriorList[j]][k]][(p+2)%3];

					break;

				}

			}

			for (i = 0; i < nNeighbor; i++)
			{

				bool count = true;

				for (p = 0; p < 3; p++)
				{

					if (element[elementValence[InteriorList[j]][i]][p] == InteriorList[j])
					{

						if (nIndex == element[elementValence[InteriorList[j]][i]][(p+3-1)%3])
						{

							iIndex[1] = element[elementValence[InteriorList[j]][i]][(p+1)%3];
							count = false;

							break;

						}

					}

				}

				if (count == false)
				{

					break;

				}

			}

			weight = 0;

			for (p = 0; p < 2; p++)
			{

				for (q = 0; q < 3; q++)
				{

					vecta[0][q] = polycubePara->vertex[InteriorList[j]][q] - polycubePara->vertex[iIndex[p]][q];

					vecta[1][q] = polycubePara->vertex[nIndex][q] - polycubePara->vertex[iIndex[p]][q];

				}

				tempLength = sqrt(pow(vecta[0][1]*vecta[1][2]-vecta[0][2]*vecta[1][1],2) + 
					pow(vecta[0][2]*vecta[1][0]-vecta[0][0]*vecta[1][2],2) + pow(vecta[0][0]*vecta[1][1]-vecta[0][1]*vecta[1][0],2));

				tempLength = (vecta[0][0]*vecta[1][0] + vecta[0][1]*vecta[1][1] + vecta[0][2]*vecta[1][2])/tempLength;

				weight += tempLength;

			}

			/*mW(j, j) -= weight;*/
			mW.coeffRef(j, j) -= weight;

			if (IsInteriorPoint(nIndex))
			{

				//double distance; // int or double???
				int distance;

				for (i = 0; i < InteriorList.size(); i++)
				{

					if (InteriorList[i] == nIndex)
					{

						distance = i;
						break;

					}

				}

				/*mW(j, distance) += weight;*/
				mW.coeffRef(j, distance) += weight;

			}
			else
			{

				/*vB(j, 0) -= weight * polycubePara->vertex[nIndex][(iConst+1)%3];
				vB(j, 1) -= weight * polycubePara->vertex[nIndex][(iConst+2)%3];*/
				vB0[j] -= weight * polycubePara->vertex[nIndex][(iConst+1)%3];
				vB1[j] -= weight * polycubePara->vertex[nIndex][(iConst+2)%3];

			}

		}

	}

	//Eigen::FullPivLU<MatrixXd> lu(mW);
	//x_solve0 = lu.solve(vB.topLeftCorner(nCount, 1));
	//x_solve1 = lu.solve(vB.topRightCorner(nCount, 1));

	//Eigen::SimplicialLLT<SparseMatrix<double> > lu;
	//lu.compute(mW);

	Eigen::ConjugateGradient<SparseMatrix<double> > lu;
	lu.compute(mW);

	x_solve0 = lu.solve(vB0);
	x_solve1 = lu.solve(vB1);

	if (lu.info() != Eigen::Success) 
	{
		cout <<"solve surface error in MappingOnePatchPostProcessing()!\n";
		return false;
	}

	int oneCorner = currentPolycubePatch.cornerPoint[0];

	for(j = 0; j < nSize; j++)
	{

		polycubePara->vertex[InteriorList[j]][iConst] = polycubePara->vertex[oneCorner][iConst];
		polycubePara->vertex[InteriorList[j]][(iConst+1)%3] = x_solve0(j);
		polycubePara->vertex[InteriorList[j]][(iConst+2)%3] = x_solve1(j);

	}

	return true;

}

////////////////////////////////////////////////////////////////////////////////////////////////////////
//For hex meshing in parametric domain

bool PolycubeConstruction::HexMeshParametricDomain(const char *outputName)
{

	/*InitializeOctree();
	ConstructeOctree();
	HexMeshOctree();
	DeleteOutsideMesh(outputName);*/

	//if (READ_IN_PARAHEX_TORUS == 1)
	//{
	//	string tempName;
	//	tempName = inputName + "_paraHex_hex.raw";
	//	RawMesh tempParametricHex;
	//	int index;
	//	tempParametricHex.Read(tempName.c_str());

	//	//int count = tempParametricHex.elementNumber;

	//	for (int i = 0; i < tempParametricHex.elementNumber; i++)
	//	{

	//		tempParametricHex.elementSign[i] = -1;

	//		double centeOne[3] = {0.f, 0.f, 0.f};
	//		for (int j = 0; j < 8; j++)
	//		{
	//			index = tempParametricHex.element[i][j];

	//			centeOne[0] += 0.125 * tempParametricHex.vertex[index][0];
	//			centeOne[1] += 0.125 * tempParametricHex.vertex[index][1];
	//			centeOne[2] += 0.125 * tempParametricHex.vertex[index][2];
	//		}


	//		//Do remember this FUCKING MISTAKE!!!
	//		//for (int j = 0; j < parametricHex.elementNumber; j++)
	//		//{

	//		//	if (fabs(tempParametricHex.vertex[i][0]-parametricHex.vertex[j][0]) < EPSILON &&
	//		//		fabs(tempParametricHex.vertex[i][1]-parametricHex.vertex[j][1]) < EPSILON && 
	//		//		fabs(tempParametricHex.vertex[i][2]-parametricHex.vertex[j][2]) < EPSILON && 
	//		//		fabs(tempParametricHex.vertex[i][3]-parametricHex.vertex[j][3]) < EPSILON &&
	//		//		fabs(tempParametricHex.vertex[i][4]-parametricHex.vertex[j][4]) < EPSILON &&
	//		//		fabs(tempParametricHex.vertex[i][5]-parametricHex.vertex[j][5]) < EPSILON &&
	//		//		fabs(tempParametricHex.vertex[i][6]-parametricHex.vertex[j][6]) < EPSILON &&
	//		//		fabs(tempParametricHex.vertex[i][7]-parametricHex.vertex[j][7]) < EPSILON)
	//		//	{

	//		//		tempParametricHex.elementSign[i] = parametricHex.elementSign[j];
	//		//		break;

	//		//	}

	//		//}

	//		for (int k = 0; k < parametricHex.elementNumber; k++)
	//		{

	//			double tempCenter[3] = {0.f, 0.f, 0.f};

	//			for (int j = 0; j < 8; j++)
	//			{
	//				index = parametricHex.element[k][j];

	//				tempCenter[0] += 0.125 * parametricHex.vertex[index][0];
	//				tempCenter[1] += 0.125 * parametricHex.vertex[index][1];
	//				tempCenter[2] += 0.125 * parametricHex.vertex[index][2];
	//			}

	//			if (fabs(tempCenter[0]-centeOne[0]) < EPSILON && fabs(tempCenter[1]-centeOne[1]) < EPSILON && fabs(tempCenter[2]-centeOne[2]) < EPSILON)
	//			{

	//				tempParametricHex.elementSign[i] = parametricHex.elementSign[k];

	//				//count--;

	//				break;

	//			}

	//		}

	//	}

	//	//cout << count<<endl;

	//	parametricHex = tempParametricHex;

	//	parametricHex.Write(outputName);

	//}

	if (READ_IN_PARAHEX_TORUS == 1)
	{
		
		string tempName;
		tempName = inputName + "_paraHex_hex.raw";
		RawMesh tempParametricHex;
		int index;
		tempParametricHex.Read(tempName.c_str());

		parametricHex = tempParametricHex;
		parametricHex.Write(outputName);

		/*string tempString;
		tempString = inputName + "_paraHex_read.k";
		ReadKFileHexMeshPara(tempString.c_str());*/
	}

	return true;

}

// For hex mesh generation in real domain
bool PolycubeConstruction::HexMeshRealDomain(const char *outputName)
{

	parametricHex.SetBoundaryVertexSign(1);

	//////////////////////////////////////////////////////////////////////
	////Useful for propagation step
	parametricHex.InitiateEdgeValence();
	//////////////////////////////////////////////////////////////////////

	parametricHex.InitiateElementValence();

	//CheckUniformHexValidity(); // Not valid for unusual polycube structures

	realDomainHex = parametricHex;

	if (READ_IN_UNIFORMHEX == 0)
	{

		realDomainHex.SetBoundaryVertexSign(1);
		realDomainHex.InitiateEdgeValence();
		realDomainHex.InitiateElementValence();

		propElement.resize(parametricHex.vertexNumber);

		HexMeshProjectionBoundary();

		HexMeshProjectionInterior();

		//Pillowing();

		//realDomainHex.Write("PillowNoSmooth_hex.raw");
		string tempName;	
		realDomainHex.Write(outputName);
		tempName = inputName + "_realHex_hex.inp";

		realDomainHex.Write(tempName.c_str());

		realDomainHex.Smooth(10000);



		tempName = inputName + "_NoPillowSmooth_hex.raw";

		realDomainHex.Write(tempName.c_str());

		tempName = inputName + "_NoPillowSmooth_hex.vtk";

		realDomainHex.Write(tempName.c_str());
		tempName = inputName + "_NoPillowSmooth_hex.inp";

		realDomainHex.Write(tempName.c_str());

	}
	else
	{

		realDomainHex.Read("uniform_hex.raw");
		realDomainHex.SetBoundaryVertexSign(1);
		realDomainHex.InitiateEdgeValence();
		realDomainHex.InitiateElementValence();

		for (int i = 0; i < realDomainHex.elementNumber; i++)
		{
			realDomainHex.elementSign[i] = parametricHex.elementSign[i];
		}
		realDomainHex.Smooth(10000);
		string tempName;	
		tempName = inputName + "_NoPillowSmooth_hex_NEW.raw";

		realDomainHex.Write(tempName.c_str());

		tempName = inputName + "_NoPillowSmooth_hex_NEW.vtk";

		realDomainHex.Write(tempName.c_str());
		tempName = inputName + "_NoPillowSmooth_hex_NEW.inp";

		realDomainHex.Write(tempName.c_str());

	}


	return true;

}

bool PolycubeConstruction::HexMeshProjectionBoundary()
{

	int i, j;

	double paraPos[3], realPos[3];

	for (i = 0; i < parametricHex.vertexNumber; i++)
	{

		if (parametricHex.vertexSign[i] != 0)
		{

			for (j = 0; j < 3; j++)
			{

				paraPos[j] = parametricHex.vertex[i][j];

			}

			Projection(paraPos, realPos);

			for (j = 0; j < 3; j++)
			{

				realDomainHex.vertex[i][j] = realPos[j];

			}

		}
		else
		{

			//for (j = 0; j < 3; j++)
			//{

			//	realDomainHex.vertex[i][j] = 0.f;

			//}

			FindPropagationBoundElements(i);

			//outputCentroids<<i<<": "<<propElement[i].px<<" "<<propElement[i].nx<<" "<<propElement[i].py<<" "<<propElement[i].ny<<" "<<propElement[i].pz<<" "<<propElement[i].nz<<endl;

		}

	}

	return true;

}

bool PolycubeConstruction::FindPropagationBoundElements(int vertexID)
{

	int i, j;
	int nVert;

	for (i = 0; i < parametricHex.edgeValenceNumber[vertexID]; i++)
	{

		nVert = parametricHex.edgeValence[vertexID][i];

		if ((parametricHex.vertex[nVert][0]-parametricHex.vertex[vertexID][0]) > EPSILON)
		{
			propElement[vertexID].px = nVert;
		}
		else if ((parametricHex.vertex[nVert][0]-parametricHex.vertex[vertexID][0]) < -EPSILON)
		{
			propElement[vertexID].nx = nVert;
		}
		else if ((parametricHex.vertex[nVert][1]-parametricHex.vertex[vertexID][1]) > EPSILON)
		{
			propElement[vertexID].py = nVert;
		}
		else if ((parametricHex.vertex[nVert][1]-parametricHex.vertex[vertexID][1]) < -EPSILON)
		{
			propElement[vertexID].ny = nVert;
		}
		else if ((parametricHex.vertex[nVert][2]-parametricHex.vertex[vertexID][2]) > EPSILON)
		{
			propElement[vertexID].pz = nVert;
		}
		else
		{
			propElement[vertexID].nz = nVert;
		}

	}

	int index;

	while (!parametricHex.vertexSign[propElement[vertexID].px])
	{

		index = propElement[vertexID].px;

		for (j = 0; j < parametricHex.edgeValenceNumber[index]; j++)
		{
			nVert = parametricHex.edgeValence[index][j];

			if ((parametricHex.vertex[nVert][0]-parametricHex.vertex[index][0]) > EPSILON)
			{

				propElement[vertexID].px = nVert;

				break;

			}
		}

	}

	while (!parametricHex.vertexSign[propElement[vertexID].nx])
	{

		index = propElement[vertexID].nx;

		for (j = 0; j < parametricHex.edgeValenceNumber[index]; j++)
		{
			nVert = parametricHex.edgeValence[index][j];

			if ((parametricHex.vertex[nVert][0]-parametricHex.vertex[index][0]) < -EPSILON)
			{

				propElement[vertexID].nx = nVert;

				break;

			}
		}

	}

	while (!parametricHex.vertexSign[propElement[vertexID].py])
	{

		index = propElement[vertexID].py;

		for (j = 0; j < parametricHex.edgeValenceNumber[index]; j++)
		{
			nVert = parametricHex.edgeValence[index][j];

			if ((parametricHex.vertex[nVert][1]-parametricHex.vertex[index][1]) > EPSILON)
			{

				propElement[vertexID].py = nVert;

				break;

			}
		}

	}

	while (!parametricHex.vertexSign[propElement[vertexID].ny])
	{

		index = propElement[vertexID].ny;

		for (j = 0; j < parametricHex.edgeValenceNumber[index]; j++)
		{
			nVert = parametricHex.edgeValence[index][j];

			if ((parametricHex.vertex[nVert][1]-parametricHex.vertex[index][1]) < -EPSILON)
			{

				propElement[vertexID].ny = nVert;

				break;

			}
		}

	}

	while (!parametricHex.vertexSign[propElement[vertexID].pz])
	{

		index = propElement[vertexID].pz;

		for (j = 0; j < parametricHex.edgeValenceNumber[index]; j++)
		{
			nVert = parametricHex.edgeValence[index][j];

			if ((parametricHex.vertex[nVert][2]-parametricHex.vertex[index][2]) > EPSILON)
			{

				propElement[vertexID].pz = nVert;

				break;

			}
		}

	}

	while (!parametricHex.vertexSign[propElement[vertexID].nz])
	{

		index = propElement[vertexID].nz;

		for (j = 0; j < parametricHex.edgeValenceNumber[index]; j++)
		{
			nVert = parametricHex.edgeValence[index][j];

			if ((parametricHex.vertex[nVert][2]-parametricHex.vertex[index][2]) < -EPSILON)
			{

				propElement[vertexID].nz = nVert;

				break;

			}
		}

	}

	return true;

}

bool PolycubeConstruction::Projection(double paraPosition[3], double *realPosition)
{

	int i, j, k, p, q, index;
	
	int patchIndex;
	bool onArcsFlag, inBoxFlag;

	double dTemp, minCoords[3], maxCoords[3], dP[3][3], dPara[3][3], dSign[3];
	double dVec[3];

	for (i = 0; i < 3; i++)
	{
		realPosition[i] = 0.f;
	}

	onArcsFlag = false;
	
	for (i = 0; i < NUM_POLYCUBE_PATCH; i++)
	{

		bool counted = false;

		for (j = 0; j < polycubePatch[i].numCorner; j++)
		{

			p = 0;

			index = polycubePatch[i].boundaryEdge[j][1];

			int indexStart = polycubePatch[i].boundaryEdge[j][0];
			int indexEnd = polycubePatch[i].boundaryEdge[j][polycubePatch[i].boundaryEdge[j].size()-1];

			for (k = 0; k < 3; k++)
			{

				if (fabs(paraPosition[k] - polycubePara->vertex[index][k]) < EPSILON)
				{
					p++;
				}
				else
				{
					q = k;
				}
			}

			//Very important for the if condition
			if (p == 2&& 
				(
				(paraPosition[q] >= polycubePara->vertex[indexStart][q] && paraPosition[q] <= polycubePara->vertex[indexEnd][q])
				||
				(paraPosition[q] <= polycubePara->vertex[indexStart][q] && paraPosition[q] >= polycubePara->vertex[indexEnd][q])
				)
				)
			{

				patchIndex = j;

				onArcsFlag = true;

				for (k = 0; k < polycubePatch[i].boundaryEdge[j].size()-1; k++)
				{

					int tempIndexOne = polycubePatch[i].boundaryEdge[j][k];
					int tempIndexTwo = polycubePatch[i].boundaryEdge[j][k+1];

					dTemp = (paraPosition[q] - polycubePara->vertex[tempIndexOne][q]) / (polycubePara->vertex[tempIndexTwo][q] - polycubePara->vertex[tempIndexOne][q]);

					if (dTemp > -EPSILON && dTemp < 1 + EPSILON)
					{

						for (int l = 0; l < 3; l++)
						{

							realPosition[l] = (1-dTemp)*vertex[tempIndexOne][l] + dTemp*vertex[tempIndexTwo][l];

						}

						counted = true;

						break;

					}

				}

				break;

			}

		}

		///////May be Just break; is OK enough? Not enough
		if (counted == true)
		{
			break;
		}

	}

	int indexCluster;
	int iConst;
	int indexPatch;
	int theta;
	if (!onArcsFlag)
	{

		for (i = 0; i < elementNumber; i++)
		{

			/*indexCluster = elementArray[i].indexCluster;
			iConst = patch_const[indexCluster];*/
			
			indexPatch = elementArray[i].indexPatch;

			for (j = 0; j < 3; j++)
			{
				bool planarDir = true;
				for (k = 0; k < polycubePatch[indexPatch].numCorner; k++)
				{
					index = polycubePatch[indexPatch].cornerPoint[k];

					if (fabs(polycubePara->vertex[index][j]-polycubePara->vertex[polycubePatch[indexPatch].cornerPoint[0]][j]) > EPSILON)
					{
						planarDir = false;
					}
				}

				if (planarDir == true)
				{
					iConst = j;
					break;
				}

				//by yu
		        if (planarDir == false)    
		        {
					iConst = 1;
					theta=1.0/sqrt(2.0);
				}
				//by yu

			}
		


			////////////////////////////////////////////////////////////////////////////////////////
			////Reduce computational cost
			//if (fabs(paraPosition[iConst] - polycubePara->vertex[element[i][0]][iConst]) > EPSILON)
			//{
			//	continue;
			//}
			///////////////////////////////////////////////////////////////////////////////////////

			inBoxFlag = true;

			for (j = 0; j < 3; j++)
			{

				for (k = 0; k < 3; k++)
				{

					dP[j][k] = vertex[element[i][j]][k];
					dPara[j][k] = polycubePara->vertex[element[i][j]][k];

				}

			}

			for (j = 0; j < 3; j++)
			{

				if (j == 0)
				{

					for (k = 0; k < 3; k++)
					{
						minCoords[k] = dPara[j][k];
						maxCoords[k] = dPara[j][k];
					}

				}
				else
				{

					for (k = 0; k < 3; k++)
					{
						minCoords[k] = min(minCoords[k], dPara[j][k]);
						maxCoords[k] = max(maxCoords[k], dPara[j][k]);
					}

				}
			}

			for (j = 0; j < 3; j++)
			{

				if (!(paraPosition[j] > minCoords[j] - EPSILON && paraPosition[j] < maxCoords[j] + EPSILON))
				{

					inBoxFlag = false;
					break;

				}

			}

			if (inBoxFlag)
			{

				//dSign is the barycentric coordinates;
				dVec[0] = (dPara[0][0]-dPara[2][0])*(dPara[1][1]-dPara[2][1])-(dPara[1][0]-dPara[2][0])*(dPara[0][1]-dPara[2][1]);
				dVec[1] = (dPara[0][2]-dPara[2][2])*(dPara[1][1]-dPara[2][1])-(dPara[1][2]-dPara[2][2])*(dPara[0][1]-dPara[2][1]);
				dVec[2] = (dPara[0][0]-dPara[2][0])*(dPara[1][2]-dPara[2][2])-(dPara[1][0]-dPara[2][0])*(dPara[0][2]-dPara[2][2]);

				if(fabs(dVec[0]) > EPSILON)
				{

					dSign[0] = ((dPara[1][1]-dPara[2][1])*(paraPosition[0]-dPara[2][0])+(dPara[2][0]-dPara[1][0])*(paraPosition[1]-dPara[2][1]))/dVec[0];
					dSign[1] = ((dPara[2][1]-dPara[0][1])*(paraPosition[0]-dPara[2][0])+(dPara[0][0]-dPara[2][0])*(paraPosition[1]-dPara[2][1]))/dVec[0];
					dSign[2] = 1 - dSign[0] - dSign[1];

				}
				else if(fabs(dVec[1]) > EPSILON)
				{

					dSign[0] = ((dPara[1][1]-dPara[2][1])*(paraPosition[2]-dPara[2][2])+(dPara[2][2]-dPara[1][2])*(paraPosition[1]-dPara[2][1]))/dVec[1];
					dSign[1] = ((dPara[2][1]-dPara[0][1])*(paraPosition[2]-dPara[2][2])+(dPara[0][2]-dPara[2][2])*(paraPosition[1]-dPara[2][1]))/dVec[1];
					dSign[2] = 1 - dSign[0] - dSign[1];

				}
				else if(fabs(dVec[2]) > EPSILON)
				{

					dSign[0] = ((dPara[1][2]-dPara[2][2])*(paraPosition[0]-dPara[2][0])+(dPara[2][0]-dPara[1][0])*(paraPosition[2]-dPara[2][2]))/dVec[2];
					dSign[1] = ((dPara[2][2]-dPara[0][2])*(paraPosition[0]-dPara[2][0])+(dPara[0][0]-dPara[2][0])*(paraPosition[2]-dPara[2][2]))/dVec[2];
					dSign[2] = 1 - dSign[0] - dSign[1];

				}
				else
				{

					dSign[0] = 0.33;
					dSign[1] = 0.33;
					dSign[2] = 1 - dSign[0] - dSign[1];

				}

				if(dSign[0] > -1e-3 && dSign[1] > -1e-3 && dSign[2] > -1e-3)
				{

					for(j = 0; j < 3; j++)
					{
						realPosition[j] = dSign[0]*dP[0][j] + dSign[1]*dP[1][j] + dSign[2]*dP[2][j];
					}

					return true;

				}

			}

		}

	}

	//return true;
	return false;

}

bool PolycubeConstruction::Propagation(int vertexID)
{

	double dTemp;

	int tempIndexOne, tempIndexTwo;

	int i;
	//if (vertexID == 41833)
	//{
	//	vertexID = vertexID;
	//}

	for (i = 0; i < 3; i++)
	{
		realDomainHex.vertex[vertexID][i] = 0.f;
	}


	tempIndexOne = propElement[vertexID].nx;
	tempIndexTwo = propElement[vertexID].px;

	dTemp = (parametricHex.vertex[vertexID][0] - parametricHex.vertex[tempIndexOne][0]) / (parametricHex.vertex[tempIndexTwo][0] - parametricHex.vertex[tempIndexOne][0]);
	//dTemp = 0.5;

	if (dTemp > -EPSILON && dTemp < 1 + EPSILON)
	{

		realDomainHex.vertex[vertexID][0] = (1-dTemp) * realDomainHex.vertex[tempIndexOne][0] + dTemp * realDomainHex.vertex[tempIndexTwo][0];

		//for (i = 0; i < 3; i++)
		//{

		//	realDomainHex.vertex[vertexID][i] += 1.0/3.0 * ((1-dTemp) * realDomainHex.vertex[tempIndexOne][i] + dTemp * realDomainHex.vertex[tempIndexTwo][i]);

		//}

	}

	tempIndexOne = propElement[vertexID].ny;
	tempIndexTwo = propElement[vertexID].py;

	dTemp = (parametricHex.vertex[vertexID][1] - parametricHex.vertex[tempIndexOne][1]) / (parametricHex.vertex[tempIndexTwo][1] - parametricHex.vertex[tempIndexOne][1]);

	if (dTemp > -EPSILON && dTemp < 1 + EPSILON)
	{

		realDomainHex.vertex[vertexID][1] = (1-dTemp) * realDomainHex.vertex[tempIndexOne][1] + dTemp * realDomainHex.vertex[tempIndexTwo][1];

		//for (i = 0; i < 3; i++)
		//{

		//	realDomainHex.vertex[vertexID][i] += 1.0/3.0 * ((1-dTemp) * realDomainHex.vertex[tempIndexOne][i] + dTemp * realDomainHex.vertex[tempIndexTwo][i]);

		//}

	}

	tempIndexOne = propElement[vertexID].nz;
	tempIndexTwo = propElement[vertexID].pz;

	dTemp = (parametricHex.vertex[vertexID][2] - parametricHex.vertex[tempIndexOne][2]) / (parametricHex.vertex[tempIndexTwo][2] - parametricHex.vertex[tempIndexOne][2]);

	if (dTemp > -EPSILON && dTemp < 1 + EPSILON)
	{

		realDomainHex.vertex[vertexID][2] = (1-dTemp) * realDomainHex.vertex[tempIndexOne][2] + dTemp * realDomainHex.vertex[tempIndexTwo][2];

		//for (i = 0; i < 3; i++)
		//{

		//	realDomainHex.vertex[vertexID][i] += 1.0/3.0 * ((1-dTemp) * realDomainHex.vertex[tempIndexOne][i] + dTemp * realDomainHex.vertex[tempIndexTwo][i]);

		//}

	}

	return true;

}

bool PolycubeConstruction::HexMeshProjectionInterior()
{

	int i;

	for (i = 0; i < parametricHex.vertexNumber; i++)
	{

		if (parametricHex.vertexSign[i] != 0)
		{
			continue;
		}

		Propagation(i);

	}

	return true;

}

bool PolycubeConstruction::HexMeshOctreeSubdivision(const char *outputName)
{
	
	for (int i = 0; i <= OCTREE_MAX_LEVEL; i++)
	{
		levelRes[i] = (1 << i);
	}
	
	int i, j, k, p, q, nNeighbor;
	Node TempNode, NewNode[19];
	Element TempElement, NewElement[8];
	bool SubFlag, eFlag, BalanceFlag;
	int iIndex, iTemp, octree_id, octree_level;
	vector<int> ElementList;
	Junction NewJunction;
	double TParaPos[3];
	double PointPos[3], TempCoords[3], dtemp, t_area, vPos8[8][3];
	int NumNode, NumElement, NewNodeIndex[19], NewElementIndex[8];

	for (i = 0; i < realDomainHex.vertexNumber; i++)
	{
		//printf("%d Sign:%d\n",i,realDomainHex.vertexSign[i]);
		if (realDomainHex.vertexSign[i] == 1)
		{
			TempNode.BoundaryFlag = true;
		}
		else
		{
			TempNode.BoundaryFlag = false;
		}

		for (j = 0; j < 3; j++)
		{
			TempNode.Coords[j] = realDomainHex.vertex[i][j];
			TempNode.ParaPos[j] = parametricHex.vertex[i][j];
		}

		Nodes.push_back(TempNode);

		//Nodes[i].NeighborElement.push_back(0);

		for (j = 0; j < realDomainHex.elementValenceNumber[i]; j++)
		{

			int index = realDomainHex.elementValence[i][j];

			Nodes[i].NeighborElement.push_back(index);

		}

	}

	for (i = 0; i < realDomainHex.elementNumber; i++)
	{

		//TempElement.OctreeID = realDomainHex.elementSign[i];

		TempElement.OctreeID = 0;

		TempElement.BoundaryFlag = false;

		for (j = 0; j < 8; j++)
		{

			TempElement.NodeIndex[j] = realDomainHex.element[i][j];

			int index = TempElement.NodeIndex[j];
			if (realDomainHex.vertexSign[index] == 1)
			{
				TempElement.BoundaryFlag = true;
			}

		}

		Elements.push_back(TempElement);

	}

	int MAX_LEVEL_TOR = 1;

	while (MAX_LEVEL_TOR <= OCTREE_MAX_LEVEL)
	{

		for (i = 0; i < Elements.size(); i++)
		{


			octree_id = Elements[i].OctreeID;
			octree_level = GetLevel(octree_id);
			/*if (octree_level >= (OCTREE_MAX_LEVEL))
			{
				continue;
			}*/
			if (octree_level >= (MAX_LEVEL_TOR))
			{
				continue;
			}

			SubFlag = true;

			if (SubFlag)
			{

				int octree_idx[8];

				NumNode = Nodes.size();
				NumElement = Elements.size();
				TempElement = Elements[i];
				Elements[i].Junctions.clear();
				NewElementIndex[0] = i;

				for (j = 0; j < 7; j++)
				{

					NewElementIndex[j+1] = NumElement + j;

				}

				//Generate 19 New Nodes
				for (j = 0; j < 19; j++)
				{

					NewNodeIndex[j] = NumNode + j;
					NewNode[j].NeighborElement.clear();

					if (j == 13)
					{
						NewNode[j].BoundaryFlag = false;
					}
					else
					{
						NewNode[j].BoundaryFlag = true;
						for (k = 0; k < NewNodeRelation[j][0]; k++)
						{
							if(!Nodes[TempElement.NodeIndex[NewNodeRelation[j][k+1]]].BoundaryFlag)
							{
								NewNode[j].BoundaryFlag = false;
								break;
							}
						}
					}

					//printf("%d %d\n", j, NewNode[j].BoundaryFlag);

					for (k = 0; k < 3; k++)
					{
						//iIndex = 0;
						NewNode[j].Coords[k] = 0.0f;
						NewNode[j].ParaPos[k] = 0.0f;

						for(p = 0; p < NewNodeRelation[j][0]; p++)
						{
							NewNode[j].Coords[k] += Nodes[TempElement.NodeIndex[NewNodeRelation[j][p+1]]].Coords[k];
							NewNode[j].ParaPos[k] += Nodes[TempElement.NodeIndex[NewNodeRelation[j][p+1]]].ParaPos[k];
						}
						NewNode[j].Coords[k] /= NewNodeRelation[j][0];
						NewNode[j].ParaPos[k] /= NewNodeRelation[j][0];

					}

					if (NewNode[j].BoundaryFlag)
					{
						Projection(NewNode[j].ParaPos, NewNode[j].Coords);
					}
				}

				//smoothing the local mesh

				if(TempElement.BoundaryFlag)
				{
					for(j = 0; j < 6; j++)
					{
						if(!NewNode[FACENODES[j][8]-8].BoundaryFlag && (Nodes[TempElement.NodeIndex[FACENODES[j][0]]].BoundaryFlag ||
							Nodes[TempElement.NodeIndex[FACENODES[j][2]]].BoundaryFlag || Nodes[TempElement.NodeIndex[FACENODES[j][4]]].BoundaryFlag || Nodes[TempElement.NodeIndex[FACENODES[j][6]]].BoundaryFlag))
						{
							for(k = 0; k < 3; k++)
							{
								PointPos[k] = 0.0f;
							}
							t_area = 0.0f;
							for(k = 0; k < 4; k++)
							{
								dtemp = QuadArea(Nodes[TempElement.NodeIndex[FACENODES[j][k*2]]].Coords, NewNode[FACENODES[j][k*2+1]-8].Coords, NewNode[FACENODES[j][8]-8].Coords, NewNode[FACENODES[j][(k+3)%4*2+1]-8].Coords, TempCoords);
								t_area += dtemp;
								for(p = 0; p < 3; p++)
								{
									PointPos[p] += TempCoords[p]*dtemp;
								}
							}
							if(t_area < EPSILON)
								continue;
							for(k = 0; k < 3; k++)
							{
								PointPos[k] /= t_area;
								NewNode[FACENODES[j][8]-8].Coords[k] = PointPos[k];
							}
						}
					}
					for(j = 0; j < 3; j++)
					{
						PointPos[j] = 0.0f;
					}
					t_area = 0.0f;
					for(j = 0; j < 8; j++)
					{
						for(k = 0; k < 8; k++)
						{
							if(NewElementNIndex[j][k*2])
							{
								for(p = 0; p < 3; p++)
									vPos8[k][p] = NewNode[NewElementNIndex[j][k*2+1]].Coords[p];
							}
							else
							{
								for(p = 0; p < 3; p++)
									vPos8[k][p] = Nodes[TempElement.NodeIndex[NewElementNIndex[j][k*2+1]]].Coords[p];
							}
						}
						dtemp = GetHexVolume(vPos8, TempCoords);
						t_area += dtemp;
						for(p = 0; p < 3; p++)
						{
							PointPos[p] += TempCoords[p]*dtemp;
						}
					}
					if(t_area > EPSILON)
					{
						for(k = 0; k < 3; k++)
						{
							PointPos[k] /= t_area;
							NewNode[13].Coords[k] = PointPos[k];
						}
					}
				}


				//End smoothing the local mesh

				//Check whether the new nodes are already exist or not
				for (j = 0; j < 19; j++)
				{

					if (NewNodeRelation[j][0] == 8)
					{
						continue;
					}

					for (k = 0; k < TempElement.Junctions.size(); k++)
					{

						eFlag = true;
						for(p = 0; p < 3; p++)
						{
							//if(!(NewNode[j].ParaPos[p] == Nodes[TempElement.Junctions[k].NodeIndex].ParaPos[p]))
							if(!(fabs(NewNode[j].ParaPos[p] - Nodes[TempElement.Junctions[k].NodeIndex].ParaPos[p]) < EPSILON))
							{
								eFlag = false;
								break;
							}
						}

						if(eFlag)
						{
							NewNodeIndex[j] = TempElement.Junctions[k].NodeIndex;
							for(p = 0; p < NewNodeRelation[j][0]; p++)
							{
								if(p == 0)
								{
									for(q = 0; q < Nodes[NewNodeIndex[j]].NeighborElement.size(); q++)
									{
										if(Nodes[NewNodeIndex[j]].NeighborElement[q] == i)
										{
											Nodes[NewNodeIndex[j]].NeighborElement[q] = NewElementIndex[NewNodeRelation[j][p+1]];
										}
									}
								}
								else
								{
									Nodes[NewNodeIndex[j]].NeighborElement.push_back(NewElementIndex[NewNodeRelation[j][p+1]]);
								}
							}
							TempElement.Junctions.erase(TempElement.Junctions.begin() + k);
							//if(Nodes[NewNodeIndex[j]].NeighborElement.size() == 8)
								//Nodes[NewNodeIndex[j]].Type = 0;
							for(p = j+1; p < 19; p++)
							{
								NewNodeIndex[p] -= 1;
							}
							break;
						}

					}

				}

				for(j = 0; j < 19; j++)
				{
					if(NewNodeIndex[j] < NumNode)
						continue;
					if(j == 13)
					{
						for(k = 0; k < 8; k++)
							NewNode[j].NeighborElement.push_back(NewElementIndex[k]);
					}
					else
					{
						for(k = 0; k < NewNodeRelation[j][0]; k++)
						{
							NewNode[j].NeighborElement.push_back(NewElementIndex[NewNodeRelation[j][k+1]]);
						}
						ElementList.clear();
						for(k = 0; k < NewNodeRelation[j][0]; k++)
						{
							iIndex = TempElement.NodeIndex[NewNodeRelation[j][k+1]];
							nNeighbor = Nodes[iIndex].NeighborElement.size();
							for(p = 0; p < nNeighbor; p++)
							{
								eFlag = false;
								iTemp = Nodes[iIndex].NeighborElement[p];
								for(q = 0; q < ElementList.size(); q++)
								{
									if(ElementList[q] == iTemp)
									{
										eFlag = true;
										break;
									}
								}
								if(iTemp == i || eFlag)
									continue;

								for (int l = 0; l < 19; l++)
								{

									double ninteenNode[3] = {0.f, 0.f, 0.f};

									for (int m = 0;  m < NewNodeRelation[l][0]; m++)
									{
										ninteenNode[0] += Nodes[Elements[iTemp].NodeIndex[NewNodeRelation[l][m+1]]].ParaPos[0];
										ninteenNode[1] += Nodes[Elements[iTemp].NodeIndex[NewNodeRelation[l][m+1]]].ParaPos[1];
										ninteenNode[2] += Nodes[Elements[iTemp].NodeIndex[NewNodeRelation[l][m+1]]].ParaPos[2];
									}

									ninteenNode[0] /= NewNodeRelation[l][0];
									ninteenNode[1] /= NewNodeRelation[l][0];
									ninteenNode[2] /= NewNodeRelation[l][0];

									if (fabs(NewNode[j].ParaPos[0] - ninteenNode[0]) < EPSILON && fabs(NewNode[j].ParaPos[1] - ninteenNode[1]) < EPSILON && fabs(NewNode[j].ParaPos[2] - ninteenNode[2]) < EPSILON)
									{
										ElementList.push_back(iTemp);
										NewNode[j].NeighborElement.push_back(iTemp);
										NewJunction.NodeIndex = NewNodeIndex[j];
										Elements[iTemp].Junctions.push_back(NewJunction);

										break;

									}

								}

								///*if(NewNode[j].ParaPos[0] >= Nodes[Elements[iTemp].NodeIndex[0]].ParaPos[0]
								//&& NewNode[j].ParaPos[1] >= Nodes[Elements[iTemp].NodeIndex[0]].ParaPos[1]
								//&& NewNode[j].ParaPos[2] >= Nodes[Elements[iTemp].NodeIndex[0]].ParaPos[2]
								//&& NewNode[j].ParaPos[0] <= Nodes[Elements[iTemp].NodeIndex[6]].ParaPos[0]
								//&& NewNode[j].ParaPos[1] <= Nodes[Elements[iTemp].NodeIndex[6]].ParaPos[1]
								//&& NewNode[j].ParaPos[2] <= Nodes[Elements[iTemp].NodeIndex[6]].ParaPos[2])*/
								//if((NewNode[j].ParaPos[0] > Nodes[Elements[iTemp].NodeIndex[0]].ParaPos[0] || fabs(NewNode[j].ParaPos[0] - Nodes[Elements[iTemp].NodeIndex[0]].ParaPos[0]) < EPSILON)
								//&& (NewNode[j].ParaPos[1] > Nodes[Elements[iTemp].NodeIndex[0]].ParaPos[1] || fabs(NewNode[j].ParaPos[1] - Nodes[Elements[iTemp].NodeIndex[0]].ParaPos[1]) < EPSILON)
								//&& (NewNode[j].ParaPos[2] > Nodes[Elements[iTemp].NodeIndex[0]].ParaPos[2] || fabs(NewNode[j].ParaPos[2] - Nodes[Elements[iTemp].NodeIndex[0]].ParaPos[2]) < EPSILON)
								//&& (NewNode[j].ParaPos[0] < Nodes[Elements[iTemp].NodeIndex[6]].ParaPos[0] || fabs(NewNode[j].ParaPos[0] - Nodes[Elements[iTemp].NodeIndex[6]].ParaPos[0]) < EPSILON)
								//&& (NewNode[j].ParaPos[1] < Nodes[Elements[iTemp].NodeIndex[6]].ParaPos[1] || fabs(NewNode[j].ParaPos[1] - Nodes[Elements[iTemp].NodeIndex[6]].ParaPos[1]) < EPSILON)
								//&& (NewNode[j].ParaPos[2] < Nodes[Elements[iTemp].NodeIndex[6]].ParaPos[2] || fabs(NewNode[j].ParaPos[2] - Nodes[Elements[iTemp].NodeIndex[6]].ParaPos[2]) < EPSILON))
								//{
								//	ElementList.push_back(iTemp);
								//	NewNode[j].NeighborElement.push_back(iTemp);
								//	NewJunction.NodeIndex = NewNodeIndex[j];
								//	Elements[iTemp].Junctions.push_back(NewJunction);
								//}

							}
							//if(eFlag)
							//	break;
						}
					}
					Nodes.push_back(NewNode[j]);

				}

				// update old nodes' neighbors;
				for(j = 0; j < 8; j++)
				{
					nNeighbor = Nodes[TempElement.NodeIndex[j]].NeighborElement.size();
					for(k = 0; k < nNeighbor; k++)
					{
						if(Nodes[TempElement.NodeIndex[j]].NeighborElement[k] == i)
						{
							Nodes[TempElement.NodeIndex[j]].NeighborElement[k] = NewElementIndex[j];
							break;
						}
					}
				}

				octree_id = TempElement.OctreeID;
				octree_level = GetLevel(octree_id);
				int xx, yy, zz;
				OctreeidxToXYZ(octree_id, xx, yy, zz, octree_level);
			
				octree_id = XYZToOctreeidx(xx*2, yy*2, zz*2, octree_level+1);
				NewElement[0].OctreeID = octree_id;

				octree_id = XYZToOctreeidx(xx*2 + 1, yy*2, zz*2, octree_level+1);
				NewElement[1].OctreeID = octree_id;

				octree_id = XYZToOctreeidx(xx*2 + 1, yy*2 + 1, zz*2, octree_level+1);
				NewElement[2].OctreeID = octree_id;

				octree_id = XYZToOctreeidx(xx*2, yy*2 + 1, zz*2, octree_level+1);
				NewElement[3].OctreeID = octree_id;

				octree_id = XYZToOctreeidx(xx*2, yy*2, zz*2 + 1, octree_level+1);
				NewElement[4].OctreeID = octree_id;

				octree_id = XYZToOctreeidx(xx*2 + 1, yy*2, zz*2 + 1, octree_level+1);
				NewElement[5].OctreeID = octree_id;

				octree_id = XYZToOctreeidx(xx*2 + 1, yy*2 + 1, zz*2 + 1, octree_level+1);
				NewElement[6].OctreeID = octree_id;

				octree_id = XYZToOctreeidx(xx*2, yy*2 + 1, zz*2 + 1, octree_level+1);
				NewElement[7].OctreeID = octree_id;

				// generate new elements;
				for(j = 0; j < 8; j++)
				{
					NewElement[j].Junctions.clear();
				
					//NewElement[j].OctreeLevel = TempElement.OctreeLevel + 1;
					//octree_id = TempElement.OctreeID;
					//octree_level = GetLevel(octree_id);
					//int xx, yy, zz;
					//OctreeidxToXYZ(octree_id, xx, yy, zz, octree_level);
					//octree_id = XYZToOctreeidx(xx*2, yy*2, zz*2, octree_level+1);

					//NewElement[j].OctreeID = octree_id;
				
					NewElement[j].BoundaryFlag = Nodes[TempElement.NodeIndex[j]].BoundaryFlag;
					for(k = 0; k < 8; k++)
					{
						if(NewElementNIndex[j][k*2])
						{
							NewElement[j].NodeIndex[k] = NewNodeIndex[NewElementNIndex[j][k*2+1]];
						}
						else
						{
							NewElement[j].NodeIndex[k] = TempElement.NodeIndex[NewElementNIndex[j][k*2+1]];
						}
					}
				}

				for(j = 0; j < TempElement.Junctions.size(); j++)
				{
					iIndex = TempElement.Junctions[j].NodeIndex;
					eFlag = true;
					for(k = 0; k < 8; k++)
					{
						for (int l = 0; l < 19; l++)
						{

							double ninteenNode[3] = {0.f, 0.f, 0.f};

							for (int m = 0;  m < NewNodeRelation[l][0]; m++)
							{
								ninteenNode[0] += Nodes[NewElement[k].NodeIndex[NewNodeRelation[l][m+1]]].ParaPos[0];
								ninteenNode[1] += Nodes[NewElement[k].NodeIndex[NewNodeRelation[l][m+1]]].ParaPos[1];
								ninteenNode[2] += Nodes[NewElement[k].NodeIndex[NewNodeRelation[l][m+1]]].ParaPos[2];
							}

							ninteenNode[0] /= NewNodeRelation[l][0];
							ninteenNode[1] /= NewNodeRelation[l][0];
							ninteenNode[2] /= NewNodeRelation[l][0];

							if (fabs(Nodes[iIndex].ParaPos[0] - ninteenNode[0]) < EPSILON && fabs(Nodes[iIndex].ParaPos[1] - ninteenNode[1]) < EPSILON && fabs(Nodes[iIndex].ParaPos[2] - ninteenNode[2]) < EPSILON)
							{
								NewElement[k].Junctions.push_back(TempElement.Junctions[j]);
								if(eFlag)
								{
									eFlag = false;
									for(p = 0; p < Nodes[iIndex].NeighborElement.size(); p++)
									{
										if(Nodes[iIndex].NeighborElement[p] == i)
										{
											Nodes[iIndex].NeighborElement[p] = NewElementIndex[k];
											break;
										}
									}
								}
								else
								{
									Nodes[iIndex].NeighborElement.push_back(NewElementIndex[k]);
								}

								break;

							}

						}

						///*if(Nodes[iIndex].ParaPos[0] >= Nodes[NewElement[k].NodeIndex[0]].ParaPos[0]
						//&& Nodes[iIndex].ParaPos[1] >= Nodes[NewElement[k].NodeIndex[0]].ParaPos[1]
						//&& Nodes[iIndex].ParaPos[2] >= Nodes[NewElement[k].NodeIndex[0]].ParaPos[2]
						//&& Nodes[iIndex].ParaPos[0] <= Nodes[NewElement[k].NodeIndex[6]].ParaPos[0]
						//&& Nodes[iIndex].ParaPos[1] <= Nodes[NewElement[k].NodeIndex[6]].ParaPos[1]
						//&& Nodes[iIndex].ParaPos[2] <= Nodes[NewElement[k].NodeIndex[6]].ParaPos[2])*/
						//if((Nodes[iIndex].ParaPos[0] > Nodes[NewElement[k].NodeIndex[0]].ParaPos[0] || fabs(Nodes[iIndex].ParaPos[0] - Nodes[NewElement[k].NodeIndex[0]].ParaPos[0]) < EPSILON)
						//&& (Nodes[iIndex].ParaPos[1] > Nodes[NewElement[k].NodeIndex[0]].ParaPos[1] || fabs(Nodes[iIndex].ParaPos[1] - Nodes[NewElement[k].NodeIndex[0]].ParaPos[1]) < EPSILON)
						//&& (Nodes[iIndex].ParaPos[2] > Nodes[NewElement[k].NodeIndex[0]].ParaPos[2] || fabs(Nodes[iIndex].ParaPos[2] - Nodes[NewElement[k].NodeIndex[0]].ParaPos[2]) < EPSILON)
						//&& (Nodes[iIndex].ParaPos[0] < Nodes[NewElement[k].NodeIndex[6]].ParaPos[0] || fabs(Nodes[iIndex].ParaPos[0] - Nodes[NewElement[k].NodeIndex[6]].ParaPos[0]) < EPSILON)
						//&& (Nodes[iIndex].ParaPos[1] < Nodes[NewElement[k].NodeIndex[6]].ParaPos[1] || fabs(Nodes[iIndex].ParaPos[1] - Nodes[NewElement[k].NodeIndex[6]].ParaPos[1]) < EPSILON)
						//&& (Nodes[iIndex].ParaPos[2] < Nodes[NewElement[k].NodeIndex[6]].ParaPos[2] || fabs(Nodes[iIndex].ParaPos[2] - Nodes[NewElement[k].NodeIndex[6]].ParaPos[2]) < EPSILON))
						//{
						//	NewElement[k].Junctions.push_back(TempElement.Junctions[j]);
						//	if(eFlag)
						//	{
						//		eFlag = false;
						//		for(p = 0; p < Nodes[iIndex].NeighborElement.size(); p++)
						//		{
						//			if(Nodes[iIndex].NeighborElement[p] == i)
						//			{
						//				Nodes[iIndex].NeighborElement[p] = NewElementIndex[k];
						//				break;
						//			}
						//		}
						//	}
						//	else
						//	{
						//		Nodes[iIndex].NeighborElement.push_back(NewElementIndex[k]);
						//	}
						//	//break;
						//}
					}
				}
				Elements[i] = NewElement[0];
				for(j = 0; j < 7; j++)
					Elements.push_back(NewElement[j+1]);
				//i--;
				i = -1;

			}

		}

		MAX_LEVEL_TOR++;

	}

	///////////////////////////////////////////////////////////////////////////////////////
	nNode = Nodes.size();
	nElement = Elements.size();

	adaptiveOctreeHex.CreateNewMesh(adaptiveOctreeHex.HEXAHEDRON, nNode, nElement);

	for (i = 0; i < nNode; i++)
	{

		for (j = 0; j < 3; j++)
		{

			adaptiveOctreeHex.vertex[i][j] = Nodes[i].Coords[j];

		}
		
	}

	for (i = 0; i < nElement; i++)
	{

		for (j = 0; j < 8; j++)
		{

			adaptiveOctreeHex.element[i][j] = Elements[i].NodeIndex[j];

		}

	}
	//adaptiveOctreeHex.Smooth(10);
	string tempName;

	tempName = inputName + "_AdaptiveOctreePhys_hex.inp";

	adaptiveOctreeHex.Write(tempName.c_str());
	tempName = inputName + "_AdaptiveOctreePhys_hex.vtk";
	adaptiveOctreeHex.Write(tempName.c_str());
	tempName = inputName + "_AdaptiveOctreePhys_hex.raw";
	adaptiveOctreeHex.Write(tempName.c_str());
	for (i = 0; i < nNode; i++)
	{

		for (j = 0; j < 3; j++)
		{

			adaptiveOctreeHex.vertex[i][j] = Nodes[i].ParaPos[j];

		}

	}

	tempName = inputName + "_AdaptiveOctreePara_hex.inp";

	adaptiveOctreeHex.Write(tempName.c_str());
	tempName = inputName + "_AdaptiveOctreePara_hex.vtk";
	adaptiveOctreeHex.Write(tempName.c_str());
	tempName = inputName + "_AdaptiveOctreePara_hex.raw";
	adaptiveOctreeHex.Write(tempName.c_str());


	return true;

}


double		PolycubeConstruction::TriArea(double v0[3], double v1[3], double v2[3])
{
	Vector3d n;
	Vector3d a(v1[0]-v0[0], v1[1]-v0[1], v1[2]-v0[2]);
	Vector3d b(v2[0]-v0[0], v2[1]-v0[1], v2[2]-v0[2]);
	n = a.cross(b);
	return n.norm();
}

double		PolycubeConstruction::QuadArea(double v0[3], double v1[3], double v2[3], double v3[3], double *MassCenter)
{
	int i, j;
	//vector<Tri> TempTriList;
	double t_area = 0.0f, area;

	Vector3d TempNormal, va, vb;

	for(i = 0; i < 3; i++)
	{
		MassCenter[i] = 0.0f;
	}

	for(i = 0; i < 2; i++)
	{
		if(i == 0)
		{
			area = TriArea(v0, v1, v2);
		}
		else
		{
			area = TriArea(v0, v2, v3);
		}
		t_area += area;
	}
	for(j = 0; j < 3; j++)
	{
		MassCenter[j] += (v0[j] + v1[j] + v2[j] + v3[j])/4;
	}
	return t_area;
}
double		PolycubeConstruction::GetHexVolume(double p[8][3], double *MassCenter)
{	
	int i, j;
	double u, v, w, su[3], sv[3], sw[3], cuv[3], volume;

	for(i = 0; i < 3; i++)
	{
		MassCenter[i] = 0.0f;
	}
	for(i = 0; i < 8; i++)
	{
		for(j = 0; j < 3; j++)
		{
			MassCenter[j] += p[i][j];
		}
	}
	for(i = 0; i < 3; i++)
	{
		MassCenter[i] /= 8;
	}

	volume = 0.0f;
	for(i = 0; i < 8; i++) 
	{
		u = 0.5f - (double)sqrt(3.0f)/6.0f;
		v = u;
		w = u;
		if(i == 1 || i == 3 || i == 5 || i == 7) u = 0.5f + (double)sqrt(3.0f)/6.0f;
		if(i == 2 || i == 3 || i == 6 || i == 7) v = 0.5f + (double)sqrt(3.0f)/6.0f;
		if(i == 4 || i == 5 || i == 6 || i == 7) w = 0.5f + (double)sqrt(3.0f)/6.0f;
		for(j = 0; j < 3; j++) 
		{
			su[j] = (1-v)*(1-w)*(p[1][j] - p[0][j]) + v*(1-w)*(p[2][j] - p[3][j]) + (1-v)*w*(p[5][j] - p[4][j]) + v*w*(p[6][j] - p[7][j]);
			sv[j] = (1-u)*(1-w)*(p[3][j] - p[0][j]) + u*(1-w)*(p[2][j] - p[1][j]) + (1-u)*w*(p[7][j] - p[4][j]) + u*w*(p[6][j] - p[5][j]);
			sw[j] = (1-u)*(1-v)*(p[4][j] - p[0][j]) + u*(1-v)*(p[5][j] - p[1][j]) + (1-u)*v*(p[7][j] - p[3][j]) + u*v*(p[6][j] - p[2][j]);
		}
		cuv[0] = su[1]*sv[2] - su[2]*sv[1];
		cuv[1] = su[2]*sv[0] - su[0]*sv[2];
		cuv[2] = su[0]*sv[1] - su[1]*sv[0];

		volume += (cuv[0]*sw[0]+cuv[1]*sw[1]+cuv[2]*sw[2])/8.0f;
	}
	return volume;
}

int PolycubeConstruction::GetLevel(int octree_id)
{
	int num = 0;
	int i = 0;
	while (1)
	{
		num += (1 << (i * 3));
		if (num > octree_id)
			break;
		i++;
	}
	return i;
}

void PolycubeConstruction::OctreeidxToXYZ(int octree_id, int &x, int &y, int &z, int level)
{
	int idx;
	int lres;//level resolution

	idx = octree_id - level_id[level];
	lres = levelRes[level];

	x = idx % lres;
	y = (idx / lres) % lres;
	z = idx / (lres * lres);

}

//Given X, Y and Z, gain the octree cell idex
//Here, X, Y and Z are all local, related to level
int PolycubeConstruction::XYZToOctreeidx(int x, int y, int z, int level, int direction)
{
	int lres;
	lres = levelRes[level];
	int xx, yy, zz;
	xx = x; yy = y; zz = z;

	if (x < 0 || y < 0 || z < 0 || x >= lres || y >= lres || z >= lres)
	{
		return -1;
	} 
	else
	{
		switch (direction)
		{
		case 0:
			x = xx; y = yy; z = zz;
			break;
		case 1:
			x = zz; y = xx; z = yy;
			break;
		case 2:
			x = yy; y = zz; z = xx;
			break;
		}

		return level_id[level] + x + y * lres + z * lres *lres;
	}

}

bool PolycubeConstruction::ReadKFileInitial(const char *inputName)
{
	int i, j, tempElementSign;
	string oneLine;

	string inputFileName = inputName;
	fstream input(inputFileName, fstream::in);
	if (!input)
	{
		cout << "ReadKFileInitial open input file error!!!" << endl;
		return false;
	}

	getline(input, oneLine);
	getline(input, oneLine);
	getline(input, oneLine);
	getline(input, oneLine);

	for (i = 0; i < elementNumber; i++)
	{
		getline(input, oneLine);
		std::istringstream iss(oneLine); // string stream
		std::string token;
		size_t pos = -1;
		int location_assign_value = 0;
		while (iss >> token) {
			//cout << token << endl;
			while ((pos = token.rfind(',')) != std::string::npos) {
				token.erase(pos, 1);
			}
			//cout << token << endl;
			//getchar();
			/*std::cout << token << '\n';
			getchar();*/
			if (location_assign_value ==1)
			{
				std::istringstream(token) >> tempElementSign;
			}
			location_assign_value++;
		}
		
		elementArray[i].indexPart = tempElementSign - 1;
	}

	input.close();

	return true;

}




bool PolycubeConstruction::ReadKFileBeforePostProcessing(const char *inputName)
{
	int i, j, tempElementSign;
	string oneLine;

	string inputFileName = inputName;
	fstream input(inputFileName, fstream::in);
	if (!input)
	{
		cout << "ReadKFileBeforePostProcessing open input file error!!!" << endl;
		return false;
	}

	getline(input, oneLine);
	getline(input, oneLine);
	getline(input, oneLine);
	getline(input, oneLine);

	for (i = 0; i < elementNumber; i++)
	{
		input >> j >> tempElementSign >> j >> j >> j >> j >> j >> j >> j >> j;
		elementArray[i].indexCluster = tempElementSign - 1;
	}

	input.close();

	return true;

}

bool PolycubeConstruction::WriteKFileBeforePostProcessing(const char *outputName)
{

	int i, j;

	string outputFileName = outputName;

	fstream output(outputFileName, fstream::out);
	if (!output)
	{
		cout << "WriteKFileBeforePostProcessing open input file error!!!" << endl;
		return false;
	}

	output << "$# LS-DYNA Keyword file created by LS-PrePost 4.2 (Beta)" << "\n";
	output << "$# Created on Jun-10-2015" << "\n";

	output << "*KEYWORD" << "\n";
	output << "*ELEMENT_SHELL" << "\n";
	for (i = 0; i < elementNumber; i++)
	{
		output <<i+1<<","<<elementArray[i].indexCluster+1<<","<<element[i][0]+1<<","<<element[i][1]+1<<","<<element[i][2]+1<<","<<element[i][2]+1<<"\n";
	}

	output << "*NODE" << "\n";
	for (i = 0; i < vertexNumber; i++)
	{
		output <<i+1<<","<<vertex[i][0]<<","<<vertex[i][1]<<","<<vertex[i][2]<<"\n";
	}

	output << "*END" << "\n";
	
	output.close();
	return true;

}

bool PolycubeConstruction::ReadKFileIndexPatch(const char *inputName)
{
	int i, j, tempElementSign;
	string oneLine;

	string inputFileName = inputName;
	fstream input(inputFileName, fstream::in);
	if (!input)
	{
		cout << "ReadKFileIndexPatch open input file error!!!" << endl;
		return false;
	}

	getline(input, oneLine);
	getline(input, oneLine);
	getline(input, oneLine);
	getline(input, oneLine);

	for (i = 0; i < elementNumber; i++)
	{
		getline(input, oneLine);
		std::istringstream iss(oneLine); // string stream
		std::string token;
		size_t pos = -1;
		int location_assign_value = 0;
		while (iss >> token) {
			//cout << token << endl;
			while ((pos = token.rfind(',')) != std::string::npos) {
				token.erase(pos, 1);
			}
			//cout << token << endl;
			//getchar();
			/*std::cout << token << '\n';
			getchar();*/
			if (location_assign_value == 1)
			{
				std::istringstream(token) >> tempElementSign;
			}
			location_assign_value++;
		}

		elementArray[i].indexPatch = tempElementSign - 1;
	}

	input.close();

	return true;

}

bool PolycubeConstruction::WriteKFileIndexPatch(const char *outputName)
{

	int i, j;

	string outputFileName = outputName;

	fstream output(outputFileName, fstream::out);
	if (!output)
	{
		cout << "WriteKFileIndexPatch open input file error!!!" << endl;
		return false;
	}

	output << "$# LS-DYNA Keyword file " << "\n";
	output << "$# Created " << "\n";

	output << "*KEYWORD" << "\n";
	output << "*ELEMENT_SHELL" << "\n";
	for (i = 0; i < elementNumber; i++)
	{
		output <<i+1<<","<<elementArray[i].indexPatch+1<<","<<element[i][0]+1<<","<<element[i][1]+1<<","<<element[i][2]+1<<","<<element[i][2]+1<<"\n";
		//cout << element[i][0] << "," << element[i][1] << "," << element[i][2] << "\n";
	}

	output << "*NODE" << "\n";
	for (i = 0; i < vertexNumber; i++)
	{
		output <<i+1<<","<<vertex[i][0]<<","<<vertex[i][1]<<","<<vertex[i][2]<<"\n";
		//cout<< vertex[i][0] << "," <<vertex[i][1] << "," << vertex[i][2] << "\n";
	}
	
	output << "*END" << "\n";

	output.close();
	return true;

}

bool PolycubeConstruction::WriteKFileMapping(const char *outputName)
{

	int i, j;

	string outputFileName = outputName;

	fstream output(outputFileName, fstream::out);
	if (!output)
	{
		cout << "WriteKFileMapping open input file error!!!" << endl;
		return false;
	}

	output << "$# LS-DYNA Keyword file created by LS-PrePost 4.2 (Beta)" << "\n";
	output << "$# Created on Jun-10-2015" << "\n";

	output << "*KEYWORD" << "\n";
	output << "*ELEMENT_SHELL" << "\n";
	for (i = 0; i < elementNumber; i++)
	{
		output <<i+1<<","<<elementArray[i].indexPatch+1<<","<<element[i][0]+1<<","<<element[i][1]+1<<","<<element[i][2]+1<<","<<element[i][2]+1<<"\n";
	}

	output << "*NODE" << "\n";
	for (i = 0; i < vertexNumber; i++)
	{
		output <<i+1<<","<<polycubePara->vertex[i][0]<<","<<polycubePara->vertex[i][1]<<","<<polycubePara->vertex[i][2]<<"\n";
	}

	output << "*END" << "\n";

	output.close();

	return true;
}

bool PolycubeConstruction::ReadKFileHexMeshPara(const char *inputName)
{

	
	return true;

}

//by yu
bool PolycubeConstruction::OutputExtraInformation(const char *outputmeshName)
{

	int i, j;

	//save polycube to new vector.

	

	string tempString = outputmeshName;

	if (OUTPUT_CORNER_POINTS == 1)
	{
		///////////////////////////////////////////////////////////////OUT PUT CORNER POINTS
		outputCornerPoints.open(tempString + "_corner.txt");

		for (int i = 0; i < numberCornerPoints; i++)
		{
			int index = cornerPoints[i];
			outputCornerPoints << index << " ";

			for (int j = 0; j < 3; j++)
			{
				outputCornerPoints << vertex[index][j] << " ";
			}

			outputCornerPoints << endl;
		}

		outputCornerPoints.close();
		//by yu

		///////////////////////////////////////////////////////////////
	}

	//cout << "test" << endl;

	//export edge

	string outputFileName;

	outputFileName = tempString + "_edge.txt";

	vector<vector <int>> polycube_edge_vector;

	fstream output2(outputFileName, fstream::out);
	if (!output2)
	{
		cout << "open input file error!!!" << endl;
		return false;
	}
	//no need+1
	//edit for hex dominant
	for (i = 0; i < NUM_POLYCUBE_PATCH; i++)
	{
		if (polycubePatch[i].numCorner!=4)
		{
			cout << i+1 << ", in LS prepost with number of corner: "<< polycubePatch[i].numCorner<<endl;  //ls prepost
		}
		
		vector <int> polycube_edge_vector_temp;
		for (j = 0; j < polycubePatch[i].numCorner - 1; j++)
		{
			polycube_edge_vector_temp.clear();
			if (polycubePatch[i].cornerPoint[j]<polycubePatch[i].cornerPoint[j + 1])
			{
				output2 << polycubePatch[i].cornerPoint[j] << "," << polycubePatch[i].cornerPoint[j + 1] << "\n";
				polycube_edge_vector_temp.push_back(polycubePatch[i].cornerPoint[j]);
				polycube_edge_vector_temp.push_back(polycubePatch[i].cornerPoint[j + 1]);
			}
			else {
				output2 << polycubePatch[i].cornerPoint[j + 1] << "," << polycubePatch[i].cornerPoint[j] << "\n";
				polycube_edge_vector_temp.push_back(polycubePatch[i].cornerPoint[j + 1]);
				polycube_edge_vector_temp.push_back(polycubePatch[i].cornerPoint[j]);
			}
			polycube_edge_vector.push_back(polycube_edge_vector_temp);

		}
		polycube_edge_vector_temp.clear();

		if (polycubePatch[i].cornerPoint[0]<polycubePatch[i].cornerPoint[polycubePatch[i].numCorner - 1])
		{
			output2 << polycubePatch[i].cornerPoint[0] << "," << polycubePatch[i].cornerPoint[polycubePatch[i].numCorner - 1] << "\n";
			polycube_edge_vector_temp.push_back(polycubePatch[i].cornerPoint[0]);
			polycube_edge_vector_temp.push_back(polycubePatch[i].cornerPoint[polycubePatch[i].numCorner - 1]);
		}
		else {
			output2 << polycubePatch[i].cornerPoint[polycubePatch[i].numCorner - 1] << "," << polycubePatch[i].cornerPoint[0] << "\n";
			polycube_edge_vector_temp.push_back(polycubePatch[i].cornerPoint[polycubePatch[i].numCorner - 1]);
			polycube_edge_vector_temp.push_back(polycubePatch[i].cornerPoint[0]);
		}
		polycube_edge_vector.push_back(polycube_edge_vector_temp);
	}

	output2.close();



	//output edge vtk
	RawMesh polycube_edge;

	int polycube_edge_nn = vertexNumber;
	int polycube_edge_ne = polycube_edge_vector.size();

	polycube_edge.CreateNewMesh(polycube_edge.LINE, polycube_edge_nn, polycube_edge_ne);

	for (i = 0; i < polycube_edge_nn; i++)
	{

		for (j = 0; j < 3; j++)
		{

			polycube_edge.vertex[i][j] = vertex[i][j];

		}

	}

	for (i = 0; i < polycube_edge_ne; i++)
	{


		for (j = 0; j < 2; j++)
		{
			polycube_edge.element[i][j] = polycube_edge_vector[i][j];
		}

	}

	polycube_edge.DeleteUnusedVertex();

	if (true) //flag_yu
	{

		outputFileName = tempString + "_polycube_edge.vtk";
		polycube_edge.Write(outputFileName.c_str());
	}

	////export face



	//outputFileName = tempString + "_face.txt";

	//fstream output3(outputFileName, fstream::out);
	//if (!output3)
	//{
	//	cout << "open input file error!!!" << endl;
	//	return false;
	//}
	////no need+1
	//int elementKnumber = 0;
	//for (i = 0; i < NUM_POLYCUBE_PATCH; i++)
	//{

	//	for (int k = 0; k < polycubePatch[i].numCorner; k++)
	//	{
	//		output3 <<  polycubePatch[i].cornerPoint[k] << ",";

	//	}
	//	output3 << "\n";
	//	elementKnumber++;
	//}


	//output3.close();





}

bool PolycubeConstruction::CreateInitialPolycube(const char *outputmeshName)
{

	int i, j;

	//save polycube to new vector.

	RawMesh polycube_initial;

	int polycube_initial_nn = vertexNumber;
	int polycube_initial_ne = NUM_POLYCUBE_PATCH;
	
	polycube_initial.CreateNewMesh(polycube_initial.QUADRILATERAL, polycube_initial_nn, polycube_initial_ne);

	for (i = 0; i < polycube_initial_nn; i++)
	{

		for (j = 0; j < 3; j++)
		{

			polycube_initial.vertex[i][j] = vertex[i][j];

		}

	}
	//getchar();

	vector<vector<int> > polycube_surface_triangle_patch;
	
	for (i = 0; i < polycube_initial_ne; i++)
	{

		if (polycubePatch[i].numCorner!=4)
		{
			if (polycubePatch[i].numCorner == 3)
			{
				vector<int>  polycube_surface_triangle_patch_temp;
				for (j = 0; j < 3; j++)
				{
					polycube_surface_triangle_patch_temp.push_back(polycubePatch[i].cornerPoint[j]);
					polycube_initial.element[i][j] = polycubePatch[i].cornerPoint[j];
				}
				polycube_initial.element[i][3] = polycubePatch[i].cornerPoint[2];
				polycube_surface_triangle_patch.push_back(polycube_surface_triangle_patch_temp);
				continue;
			}
			else
			{
				cout << "Not equal to 4 and 3" << endl;
				continue;
			}
			
		}

		for (j = 0; j < 4; j++)
		{
			polycube_initial.element[i][j] =polycubePatch[i].cornerPoint[j];
		}		

	}

	polycube_initial.DeleteUnusedVertex();
	
	string tempString= outputmeshName;
	string outputFileName;
	if (false) //flag_yu
	{
		outputFileName = tempString + "_without_hybrid_.vtk";
		polycube_initial.Write(outputFileName.c_str());
	}
	

	

	outputFileName = tempString + ".k";

	fstream output(outputFileName, fstream::out);
	if (!output)
	{
		cout << "open input file error!!!" << endl;
		return false;
	}

	output << "$# LS-DYNA Keyword file created by LS-PrePost 4.2 (Beta)" << "\n";
	output << "$# Created on Jun-10-2015" << "\n";

	output << "*KEYWORD" << "\n";
	output << "*ELEMENT_SHELL" << "\n";
	int elementKnumber = 0;
	for (i = 0; i < polycube_initial.elementNumber; i++)
	{	
		output << elementKnumber + 1 << "," << elementKnumber + 1 ;
		for (j = 0; j < 4; j++)
		{
			output << "," << polycube_initial.element[i][j] + 1;


		}
		elementKnumber++;
		output << "\n";
	}

	/*for (i = 0; i < polycube_surface_triangle_patch.size(); i++)
	{
		output << elementKnumber + 1 << "," << elementKnumber + 1;
		for (j = 0; j < 3; j++)
		{
			output << "," << polycube_surface_triangle_patch[i][j] + 1;


		}
		output << "," << polycube_surface_triangle_patch[i][2] + 1;
		elementKnumber++;
		output << "\n";
	}*/

	output << "*NODE" << "\n";

	for (i = 0; i < polycube_initial.vertexNumber; i++)
		{
			output<<i+1;

			for (j = 0; j < 3; j++)
			{
				output<<"," << polycube_initial.vertex[i][j]<<" ";
			}

			output<<endl;
		}

	output << "*END" << "\n";
	
	output.close();


	


	
	

}
