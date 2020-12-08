#include "physical_polycube.h"
#include "StaticVars.h"
#include "CellQueue.h"
#include <sstream>
#include <array>

const double PhysicalPolycube::WEIGHT_LENGTH_EWCVT = 1.0f; 
const double PhysicalPolycube::DOMAIN_SIZE = 8.0f; //64.0f

PhysicalPolycube::PhysicalPolycube(void)
{

}

PhysicalPolycube::~PhysicalPolycube(void)
{

}

bool PhysicalPolycube::Initialization(string inputFileName, string polycubefilename)
{
	int i, j;
	ReadKTri(inputFileName.c_str()); //Read(tempName.c_str()); //in order to compatiable
	InitiateEdgeValence();
	InitiateElementValence();
	if (Debug_yu == true)
	{
		Write("test_1_tri.vtk");
	}
	surface_mesh_tri_.CreateNewMesh(surface_mesh_tri_.TRIANGLE, vertexNumber, elementNumber);

	for (i = 0; i < vertexNumber; i++)
	{

		for (j = 0; j < 3; j++)
		{

			surface_mesh_tri_.vertex[i][j] = vertex[i][j];

		}

	}

	for (i = 0; i < elementNumber; i++)
	{

		for (j = 0; j < 3; j++)
		{

			surface_mesh_tri_.element[i][j] = element[i][j];
			//cout<< element[i][j]<<" ";

		}
		//cout << endl;

	}

	surface_mesh_tri_.InitiateEdgeValence();
	//////////////////////////////////////////////////////////////////////

	surface_mesh_tri_.InitiateElementValence();
	if (Debug_yu == true)
	{
		cout << vertexNumber<< elementNumber << endl;
		cout << surface_mesh_tri_.vertexNumber << surface_mesh_tri_.elementNumber << endl;
		Write("check_surface_tri_original.vtk");
		surface_mesh_tri_.Write("check_surface_tri_in_surface_mesh_tri_.vtk");
	}

	

	inputName = inputFileName;

	string tempName;
	tempName = polycubefilename;


	ReadKHex(tempName.c_str());
	/*polycube_structure_hex_.InitiateEdgeValence();
	polycube_structure_hex_.InitiateElementValence();*/

	//cout << "run 0 here" << endl;
	




	if (Debug_yu==true)
	{
		tempName = "polycube_hex.vtk";
		polycube_structure_hex_.Write(tempName.c_str());
	}
	
	
	//getchar();
	//cout << polycube_structure_hex_.elementProperty.edgeNumber << endl;  //12
	//cout << "run 0_0 here" << endl;
	//getchar();
	elementArray.resize(elementNumber);
	InitializeElement();

	//SearchNei();

	return true;

}

bool PhysicalPolycube::ReadKHex(const char * filename)
{
	FILE *input;
	int i, j, elementType;
	int v[8];
	int elementNumber_hex;
	int vertexNumber_hex;
	ifstream infile;
	string oneLine, oneWrod, formerLine;

	infile.open(filename);
	if (infile.is_open())
	{
		while (1)
		{
			getline(infile, oneLine);

			if (oneLine.find("ELEMENT_SOLID") != std::string::npos)
			{
				break;
			}
		}

		elementNumber_hex = 0;
		while (1)
		{
			getline(infile, oneLine);

			if (oneLine.find("NODE") != std::string::npos)
			{
				break;
			}
			elementNumber_hex++;
		}
		vertexNumber_hex = 0;
		while (1)
		{
			getline(infile, oneLine);

			if (oneLine.find("*") != std::string::npos)
			{
				break;
			}
			vertexNumber_hex++;
		}
		cout << "Triangle inforamtion: vertex number " << vertexNumber_hex << " element number " << elementNumber_hex << endl;
		if (!polycube_structure_hex_.CreateNewMesh(polycube_structure_hex_.HEXAHEDRON, vertexNumber_hex, elementNumber_hex)) {
			return false;
		}


		infile.clear();
		infile.seekg(0);

		while (1)
		{
			getline(infile, oneLine);

			if (oneLine.find("ELEMENT_SOLID") != std::string::npos)
			{
				break;
			}
		}

		for (i = 0; i < elementNumber_hex; ++i)
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
				if (location_assign_value >= 2 && location_assign_value <= 9)
				{
					std::istringstream(token) >> polycube_structure_hex_.element[i][location_assign_value - 2];

				}
				location_assign_value++;
			}
		}
		for (i = 0; i < elementNumber_hex; ++i)
		{
			for (int loopj = 0; loopj < 8; loopj++)
			{
				polycube_structure_hex_.element[i][loopj] = polycube_structure_hex_.element[i][loopj] - 1;
			}
			
			/*cout << polycube_structure_hex_.element[i][0] << " "<<polycube_structure_hex_.element[i][1] << " " << polycube_structure_hex_.element[i][2] << endl;
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

		for (i = 0; i < vertexNumber_hex; ++i)
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
					std::istringstream(token) >> polycube_structure_hex_.vertex[i][location_assign_value - 1];
					/*cout << vertex[i][location_assign_value - 1] << endl;
					getchar();*/
				}
				location_assign_value++;
			}
		}


		//getline(infile, oneLine);
		//istringstream strStream(oneLine);
		//strStream >> oneWrod >> vertexNumber_hex >> oneWrod;
		////cout << vertexNumber_hex << endl;

		////vector<array<double, 3>> pts;
		////pts.resize(vertexNumber_hex);
		//vector<vector<double> > pts(vertexNumber_hex, vector<double>(3));

		//for (i = 0; i < vertexNumber_hex; ++i)
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
		//strStream >> oneWrod >> elementNumber_hex >> oneWrod;

		//if (!CreateNewMesh(HEXAHEDRON, vertexNumber_hex, polycube_structure_hex_.elementNumber_hex)) {
		//	return false;
		//}

		//for (i = 0; i < vertexNumber_hex; ++i)
		//{
		//	vertex[i][0] = pts[i][0];
		//	vertex[i][1] = pts[i][1];
		//	vertex[i][2] = pts[i][2];
		//	//cout << vertex[i][0] << endl;
		//}
		////cout << elementNumber_hex << endl;
		//for (i = 0; i < elementNumber_hex; ++i)
		//{
		//	infile >> oneWrod >> polycube_structure_hex_.element[i][0] >> polycube_structure_hex_.element[i][1] >> polycube_structure_hex_.element[i][2] >> polycube_structure_hex_.element[i][3] >>
		//		polycube_structure_hex_.element[i][4] >> polycube_structure_hex_.element[i][5] >> polycube_structure_hex_.element[i][6] >> polycube_structure_hex_.element[i][7];
		//	//cout << polycube_structure_hex_.element[i][0] << endl;
		//}

		infile.close();
	}
	else
	{
		cerr << "Cannot open " << filename << "!\n";
	}
	return true;
}

bool PhysicalPolycube::ReadKTri(const char * filename)
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

bool PhysicalPolycube::InitializeSkeletonPoints()
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

bool PhysicalPolycube::RecursiveLocalCoordinateSystem(int start_vertex)
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
bool PhysicalPolycube::CalculateLocalCoordinateSystem(double vz[3], double vxpre[3], double vx[3], double vy[3])
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

bool PhysicalPolycube::CalculateRotationMatrix(SkeletonPoint &sp)
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

bool PhysicalPolycube::OutputLocalCoordinateSystem()
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


bool PhysicalPolycube::InitializeElement()
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

		//Assign the triangle to the nearest skeleton point and update the new normal vector
		tempDist = 10000000.0f;
		indexNearestSkeleton = 0;

		for (j = 0; j < curveSkeleton.vertexNumber; j++)
		{
			dist = GetPhysicalDist(j, elementArray[i]);

			if (dist < tempDist)
			{
				indexNearestSkeleton = j;
				tempDist = dist;
			}
		}

		elementArray[i].indexSkeletonPoint = indexNearestSkeleton;

		//Calculate new normal for each triangle through a rotation matrix
		int tempIndex = elementArray[i].indexSkeletonPoint;
		//for (j = 0; j < 3; j++)
		//{
		//	elementArray[i].normalNew[j] = skeletonPoint[tempIndex].transMatrix[j][0] * elementArray[i].normal[0] + skeletonPoint[tempIndex].transMatrix[j][1] * elementArray[i].normal[1]
		//	+ skeletonPoint[tempIndex].transMatrix[j][2] * elementArray[i].normal[2]; // be careful about i and j here
		//}

	}

	GetNeighboringElementInRings();

	return true;

}

bool PhysicalPolycube::AssignPartToElement()
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

bool PhysicalPolycube::AssignClusterToElement()
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

			elementArray[i].indexPatch = tempPart;

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

double PhysicalPolycube::GetPhysicalDist(int index, const CVTElement &currentElement)
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

bool PhysicalPolycube::GetDirectNeighboringElementByElement(int elementID)
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
bool PhysicalPolycube::GetNeighboringElementInRings(int ringNumber)
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


bool PhysicalPolycube::InitializeSegments()
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

		elementArray[i].indexPatch = indexNearestGenerator;

		elementArray[i].indexNeiClusters[0] = elementArray[i].indexPatch;
		elementArray[i].numNeiElementEachCluster[0] = 1;

		generators[indexNearestGenerator].numElements++;

	}

	//Recalculate the generators

	int indexPatch;

	for (i = 0; i < NUM_CLUSTER; i++)
	{

		for (j = 0; j < 3; j++)
		{

			generators[i].normal[j] = 0.f;

		}

	}

	for (i = 0; i < elementNumber; i++)
	{

		indexPatch = elementArray[i].indexPatch;

		for (j = 0; j < 3; j++)
		{

			//generators[indexCluster].normal[j] += elementArray[i].normal[j];
			generators[indexPatch].normal[j] += elementArray[i].normalNew[j];

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


bool PhysicalPolycube::InitializeGeneratorsByInput()
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

double PhysicalPolycube::GetNormalDist(const Centroid &currentGenerator, const CVTElement & currentElement)
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

double PhysicalPolycube::NormalDist(double veca[], double vecb[])
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

bool PhysicalPolycube::NormalizeGeneratorNormal(Centroid &currentGenerator)
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

bool PhysicalPolycube::SearchNei()
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

bool PhysicalPolycube::IsCounted(CVTElement &currentElement, int indexClusterNeiElement, int &position)
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


bool PhysicalPolycube::OutputPatchesVTK(const char *outputName)
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

bool PhysicalPolycube::OutputPatchesVTKPara(const char *outputName)
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




bool PhysicalPolycube::OutputPatchesVTKBifurcation(const char *outputName)
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

bool PhysicalPolycube::EdgeWeightedCVT()
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

bool PhysicalPolycube::DataTransfer(CVTElement &currentElement, int newIndex)
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

int PhysicalPolycube::GetShortestEWDist(CVTElement &currentElement)
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

double PhysicalPolycube::GetEWDist(const Centroid &currentGenerator, const CVTElement & currentElement)
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


bool PhysicalPolycube::PreProcessing()
{
	string tempName;
	int i, j, k, faceVertexNumber;
	int *vertexIndex = NULL;
	
	//0. extract surface
	RawMesh * exract_surface_hex;
	exract_surface_hex=  polycube_structure_hex_.ExtractSurface();
	if (Debug_yu==true)
	{
		tempName = "extract_surface_quad.raw";
		exract_surface_hex->Write(tempName.c_str());
	}
	/*tempName = "extract_surface_quad.raw";
	exract_surface_hex->Write(tempName.c_str());	*/

	faceVertexNumber = polycube_structure_hex_.GetElementProperty(polycube_structure_hex_.elementProperty.faceType).vertexNumber;
	//cout << faceVertexNumber << endl;
	polycube_structure_hex_.InitiateMatrix(vertexIndex, faceVertexNumber);

	//1. initial polycube element structure
	polycube_element_.resize(polycube_structure_hex_.elementNumber);
	

	//2. get node information
	for (i = 0; i<polycube_structure_hex_.elementNumber; ++i)
	{
		
		polycube_element_[i].octree_ID = 0;


		for (j = 0; j < polycube_structure_hex_.elementProperty.vertexNumber; ++j)
		{
			polycube_element_[i].node_index[j] = polycube_structure_hex_.element[i][j];
		}
		

		//3. get face information

		
		for (j = 0; j < polycube_structure_hex_.elementProperty.faceNumber; ++j)
			{
				polycube_structure_hex_.GetElementFace(i, j, vertexIndex);
				for (k = 0; k < faceVertexNumber; ++k)
				{
					polycube_element_[i].face_id[j][k] = vertexIndex[k];
				}
				//4. get boundary information
				//5. get patch information
				if (polycube_structure_hex_.IsBoundaryFace(vertexIndex))
				{
					polycube_element_[i].boundary_sign[j] = 1;
					polycube_element_[i].patch_sign[j] = -1;
					//ToDO match the patch number Done
					//6. Obtain the surface relationship between \textbf{mesh A} and \textbf{mesh B}
					for (int inner_j = 0; inner_j < NUM_POLYCUBE_PATCH; inner_j++)
					{
						int sum_inner = 0;
						for (int inner_i = 0; inner_i < faceVertexNumber; inner_i++)
						{
							
							int vertex_number_on_face = polycubePatch[inner_j].numCorner;
							
							
							for (int inner_k = 0; inner_k < vertex_number_on_face; inner_k++)
							{
								double sum_lsq = 0;
								for (int loop_coordinate = 0; loop_coordinate < 3; loop_coordinate++)
								{
									sum_lsq= sum_lsq+pow(polycube_structure_hex_.vertex[vertexIndex[inner_i]][loop_coordinate] - vertex[polycubePatch[inner_j].cornerPoint[inner_k]][loop_coordinate], 2);
									
								}
								
								sum_lsq = sqrt(sum_lsq);
								//cout << sum_lsq << endl;
								if (sum_lsq<EPSILON)
								{
									sum_inner++;
								}
							}

						}
						if (sum_inner == faceVertexNumber)
						{
							polycube_element_[i].patch_sign[j] = inner_j;
							break;

						}
					}
					
				}
				else
				{
					polycube_element_[i].boundary_sign[j] = 0;
					polycube_element_[i].patch_sign[j] = -1;
				}
			}
		
	}

	FreeMatrix(vertexIndex);																

	//3. get global face ID 
	int loop_yu = 0;
	for (i = 0; i < polycube_structure_hex_.elementNumber; ++i)
	{
		for (j = 0; j < polycube_structure_hex_.elementProperty.faceNumber; ++j)
		{
			polycube_element_[i].face_global[j] = loop_yu;
			loop_yu++;
		}
	}
	//3.2 reduce ID

	for (i = 0; i < polycube_structure_hex_.elementNumber*6-1; ++i)  //compare face
	{
		for (j = i+1; j < polycube_structure_hex_.elementNumber*6; ++j)
		{

			int polycube_structure_hex_A = i / 6;
			int polycube_structure_hex_A_face = i % 6;
			int polycube_structure_hex_B = j / 6;
			int polycube_structure_hex_B_face = j % 6;
			std::vector<int> temp_A(polycube_element_[polycube_structure_hex_A].face_id[polycube_structure_hex_A_face], polycube_element_[polycube_structure_hex_A].face_id[polycube_structure_hex_A_face]+4);
			std::vector<int> temp_B(polycube_element_[polycube_structure_hex_B].face_id[polycube_structure_hex_B_face], polycube_element_[polycube_structure_hex_B].face_id[polycube_structure_hex_B_face]+4);
			
			/*if (polycube_structure_hex_A == 0 && polycube_structure_hex_B == 2)
			{
				for (int temp_i = 0; temp_i < temp_A.size(); ++temp_i)
					std::cout << temp_A[temp_i] << ' ';
				cout << endl;
				for (int temp_i = 0; temp_i < temp_B.size(); ++temp_i)
					std::cout << temp_B[temp_i] << ' ';
				cout << endl;
				getchar();
			}*/
			//getchar();
			std::sort(temp_A.begin(), temp_A.end());
			std::sort(temp_B.begin(), temp_B.end());
			if (temp_A == temp_B)
			{
				//cout << "yes" << endl;
				if (polycube_element_[polycube_structure_hex_A].face_global[polycube_structure_hex_A_face]>= polycube_element_[polycube_structure_hex_B].face_global[polycube_structure_hex_B_face])
				{
					polycube_element_[polycube_structure_hex_A].face_global[polycube_structure_hex_A_face] = polycube_element_[polycube_structure_hex_B].face_global[polycube_structure_hex_B_face];
				}
				else
				{
					polycube_element_[polycube_structure_hex_B].face_global[polycube_structure_hex_B_face]= polycube_element_[polycube_structure_hex_A].face_global[polycube_structure_hex_A_face] ;
				}
				std::vector<int> temp_compare_index{ i, j };
				//cout << temp_compare_index.size()<<"===="<<endl;
				match_surface_index.push_back(temp_compare_index);
			}
			
		}
	}
	if (Debug_yu==true)
	{
		cout << match_surface_index.size() << match_surface_index[0].size() << endl;
		for (int loop_i = 0; loop_i < match_surface_index.size(); loop_i++)
		{
			for (int loop_j = 0; loop_j < match_surface_index[0].size(); loop_j++)
			{
				cout << match_surface_index[loop_i][loop_j] << " ";
			}
			cout << endl;
		}
	}
	

	//getchar();

	return true;

}

bool PhysicalPolycube::IsBoundaryElement(CVTElement &currentElement)
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

bool PhysicalPolycube::EnforceLabelConnectivity()
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

bool PhysicalPolycube::CheckLabelConnectivity()
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

bool PhysicalPolycube::EnforceBoundaryConnectivity()
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

bool PhysicalPolycube::CheckBoundaryConnectivity()
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

bool PhysicalPolycube::InitializePolycube(string inputFileName)
{
	
	
	string tempString;
	tempString = inputFileName;

	ReadKFileIndexPatch(tempString.c_str());
	int labelMAX = -1;
	int label;
	int i, index;
	int j, k;


	/*cout << "**********************************************************" << endl;
	cout << "Initializing the polycube structure!" << endl;*/

	label = 0;

	vector<vector<int> > tempCluster(elementNumber);

	for (i = 0; i < elementNumber; i++)
	{

		tempCluster[i].push_back(-1);
		tempCluster[i].push_back(elementArray[i].indexPatch);

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

					if (0 > tempCluster[index][0] && elementArray[i].indexPatch == elementArray[index].indexPatch
) //For bifurcation
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
		polycubePatch[i].indexPatch = tempCluster[index][1];

	}

	//OutputPatchesVTKBifurcation("test.vtk");

	if (Debug_yu==true)
	{
		tempString = inputName + "_indexPatch_write.k";
		WriteKFileIndexPatch(tempString.c_str());
		tempString = inputName + "_segmentation_colors.vtk";
		OutputPatchesVTK(tempString.c_str());
	}
	
	InitiateElementValence();

	polycube_structure_hex_.SetBoundaryVertexSign(1);


	/*for (int loopj = 0; loopj < vertexNumber; loopj++)
	{
		if (IsCornerPoint(loopj))
		{
			cout << loopj << endl;
		}
	}*/
	/*cout << "polycube_boundary" << endl;
	for (int loopi = 0; loopi < polycube_structure_hex_.vertexNumber; loopi++)
	{
		if (polycube_structure_hex_.vertexSign[loopi])
		{
			cout << loopi << ", ";
		}
	}

	cout << "patch corner" << endl;
	for (int loopi = 0; loopi < vertexNumber; loopi++)
	{
		if (IsCornerPoint(loopi))
		{
			cout << loopi << ", ";
		}
	}*/

	//cout << "start" << endl;
	for (int loopi = 0; loopi < polycube_structure_hex_.vertexNumber; loopi++)
	{
		if (polycube_structure_hex_.vertexSign[loopi])
		{
			double sum_lsq = INT_MAX;
			int loopj_index = -1;
			//cout << sum_lsq << endl;
			for (int loopj = 0; loopj < vertexNumber; loopj++)
			{
				if (IsCornerPoint(loopj))
				{
					double new_sum_lsq = 0;

					for (int loop_coordinate = 0; loop_coordinate < 3; loop_coordinate++)
					{
						new_sum_lsq = new_sum_lsq + pow(polycube_structure_hex_.vertex[loopi][loop_coordinate] - vertex[loopj][loop_coordinate], 2);

					}

					new_sum_lsq = sqrt(new_sum_lsq);
					//cout << sum_lsq << endl;
					if (new_sum_lsq<=sum_lsq)
					{
						loopj_index = loopj;
						sum_lsq = new_sum_lsq;
					}
					//cornerPoints.push_back(i);

				}

			}
			if (loopj_index<0)
			{
				cout << "check maybe problem" << endl;
			}
			for (int loop_coordinate = 0; loop_coordinate < 3; loop_coordinate++)
			{
				polycube_structure_hex_.vertex[loopi][loop_coordinate] = vertex[loopj_index][loop_coordinate];

			}

			//cout << loopi << ", " << loopj_index << endl;
		
		
		}
	}
	//cout << "end" << endl;


	if (Debug_yu == true)
	{
		string tempName =  "polycube_hex_edit.vtk";
		polycube_structure_hex_.Write(tempName.c_str());
	}

	

	SearchCornerandBoundary();

	
	ModifyWrongBoundaryElements();

	
	
	if (Debug_yu == true)
	{
		cout << "run here physical_polycube cpp" << endl;
	}


	//ParametricMappingCorner();
	//CreateInitialPolycube();
	//ParametricMapping();
	
	return true;

}

bool PhysicalPolycube::CornerMappingMaxMin(Polycube &currentPolycubePatch)
{

	int indexPatch;
	int index;
	int i;
	double minMax = 0.f;

	indexPatch = currentPolycubePatch.indexPatch;

	switch (indexPatch)
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

bool PhysicalPolycube::ParametricMappingCorner()
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

bool PhysicalPolycube::SmoothBoundaryCurve()
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
				int indexPatch_one, indexPatch_two;

				indexPatch_one = polycubePatch[i].indexPatch;

				index = polycubePatch[i].boundaryEdge[j][k];

				if (IsCornerPoint(index) == true)
				{
					continue;
				}

				for (p = 0; p < elementValenceNumber[index]; p++)
				{
					int tempIndex = elementValence[index][p];

					if (elementArray[tempIndex].indexPatch == indexPatch_one)
					{
						count++;
						neiElements.push_back(tempIndex);
					}
					else
					{
						indexPatch_two = elementArray[tempIndex].indexPatch;
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

bool PhysicalPolycube::CurveFittingBoundaryCurve()
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

void PhysicalPolycube::SplineKnots(vector<int> &u, int n, int t)
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

void PhysicalPolycube::SplineCurve(const vector<XYZ> &inp, int n, const vector<int> &knots, int t, vector<XYZ> &outp, const vector<double> &ratio)
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

void PhysicalPolycube::SplineCurve(const vector<XYZ> &inp, int n, const vector<int> &knots, int t, vector<XYZ> &outp, int res)
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

void PhysicalPolycube::SplinePoint(const vector<int> &u, int n, int t, double v, const vector<XYZ> &control, XYZ &output)
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

//void PhysicalPolycube::SplinePoint(const vector<int> &u, int n, int t, double v, const vector<XYZ> &control, vector<XYZ> &output)
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

double PhysicalPolycube::SplineBlend(int k, int t, const vector<int> &u, double v)
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

bool PhysicalPolycube::SearchCornerandBoundary()
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

				//if (elementArray[index].indexPatch == polycubePatch[i].indexPatch && IsBoundaryPoint(vertexCCW) && IsBoundaryEdge(tempNodes[c], vertexCCW))
				if (elementArray[index].indexPatch == polycubePatch[i].index && IsBoundaryPoint(vertexCCW) && IsBoundaryEdge(tempNodes[c], vertexCCW))
				{

					tempNodes.push_back(vertexCCW);

					count++;

				}
				//else if (elementArray[index].indexPatch == polycubePatch[i].indexPatch && IsCornerPoint(vertexCCW) && IsBoundaryEdge(tempNodes[c], vertexCCW))
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

bool PhysicalPolycube::IsCornerPoint(int vertexID)
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

int PhysicalPolycube::FindOneCornerPoint(Polycube &currentPolycubePatch)
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

bool PhysicalPolycube::IsBoundaryPoint(int vertexID)
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

bool PhysicalPolycube::IsBoundaryEdge(int vertexIDone, int vertexIDtwo)
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


int PhysicalPolycube::SearchCCWVertex(int vertexID, int elementID)
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

bool PhysicalPolycube::ModifyWrongBoundaryElements()
{

	int i, j, k;
	int index;
	
	if (Debug_yu == true)
	{
		Write("test_2_0_tri.vtk");
	}

	for (i = 0; i < NUM_POLYCUBE_PATCH; i++)
	{
		//cout << endl;
		for (j = 0; j < polycubePatch[i].numCorner; j++)
		{

			int count = 0;
			vector<int> neiElements;

			index = polycubePatch[i].cornerPoint[j];
			//cout << index<< ", ";
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
			//following is for lS prepost index from 1
			/*cout << "node index"<<index+1 << endl;
			cout << "element ID: ";*/
			/*for (int loop_yu = 0; loop_yu < neiElements.size(); loop_yu++)
			{
				cout << neiElements[loop_yu]+1 << ", ";
			}*/
			//cout <<  endl;

			//getchar();

			if (count == 1)
			{
				//cout << "happens here" << endl;
				if (neiElements.size() != 1)
				{
					cout <<"ERROR!!! neiElement.size != 1"<<endl;
				}

				//EdgeFlipTwoElements(index, neiElements[0]);
				/*cout << "problem area" << endl;
				cout<<index<< ", ";
				for (int loop_yu = 0; loop_yu < neiElements.size(); loop_yu++)
				{
					cout << neiElements[loop_yu] << ", " ;
				}
				cout<< endl;
				cout << "problem area" << endl;*/
			}

		}

	}
	
	//Reinitialization
	////////////////////////////////////////////////////////////////////////////
	if (Debug_yu == true)
	{
		Write("test_2_tri.vtk");
	}

	InitiateElementValence();

	

	if (Debug_yu == true)
	{
		cout << "run here ModifyWrongBoundaryElements" << endl;
		cout << NUM_POLYCUBE_PATCH << endl;
	}

	InitializeElement();
	///////////////////////////////////////////////////////////////////////////
	
	return true;

}

//This function will modify the original vertex and element of the mesh, be careful
bool PhysicalPolycube::EdgeFlipTwoElements(int vertexID, int elementID)
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

bool PhysicalPolycube::ParametricMapping(int patch_ID, const vector<vector<int> > &corner_coordinate)
{
	
	

	for (int loopi = 0; loopi < corner_coordinate.size(); loopi++)
	{
		for (int loopj = 1; loopj < corner_coordinate[loopi].size(); loopj++)
		{
			//cout << corner_coordinate[loopi][loopj] << ", ";
			polycubePara->vertex[corner_coordinate[loopi][0]][loopj - 1] = corner_coordinate[loopi][loopj];
		}
		//cout << endl;
	}
	
	//To Do given polycubepara
	/*polycubePara->vertex[index][0] = tempX;
	polycubePara->vertex[index][1] = tempY;
	polycubePara->vertex[index][2] = tempZ;*/

		//ParametricMappingCornerByInput();
		ParametricMappingEdgeGeneral(patch_ID);

		ParametricMappingInterior(patch_ID);

		/*polycubePara->Write("mapping_tri.raw");
		polycubePara->Write("mapping_tri.inp");

		string tempString;
		tempString = inputName + "_mapping_write.k";
		WriteKFileMapping(tempString.c_str());*/



	return true;

}

bool PhysicalPolycube::ParametricMappingCornerByInput()
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

bool PhysicalPolycube::ParametricMappingEdge()
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


bool PhysicalPolycube::ParametricMappingEdgeGeneral(int patch_ID)
{

	int i, j, k, l;
	int index, indexPrev, indexStart, indexEnd;
	int edgeSize;
	double direction[3];

	int tempIntStart, tempIntEnd;

	double totalLength, tempLength, totalLengthPara;
	i = patch_ID; 
	

		for (j = 0; j < polycubePatch[i].boundaryEdge.size(); j++)
		{

			edgeSize = polycubePatch[i].boundaryEdge[j].size();

			indexStart = polycubePatch[i].boundaryEdge[j][0];
			indexEnd = polycubePatch[i].boundaryEdge[j][edgeSize - 1];

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

			for (k = 1; k < edgeSize - 1; k++) // k < edgeSize is very important!
			{

				index = polycubePatch[i].boundaryEdge[j][k];
				indexPrev = polycubePatch[i].boundaryEdge[j][k - 1];

				tempLength = 0.f;

				for (l = 0; l < 3; l++)
				{

					tempLength += pow(vertex[index][l] - vertex[indexPrev][l], 2);

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
			index = polycubePatch[i].boundaryEdge[j][edgeSize - 1];
			indexPrev = polycubePatch[i].boundaryEdge[j][edgeSize - 2];
			tempLength = 0.f;

			for (l = 0; l < 3; l++)
			{

				tempLength += pow(vertex[index][l] - vertex[indexPrev][l], 2);

			}

			tempLength = sqrt(tempLength);

			totalLength += tempLength;
			///////////////////////////////////////////////////////////

			for (k = 1; k < edgeSize - 1; k++)
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

	


	return true;

}




bool PhysicalPolycube::ParametricMappingInterior(int patch_ID)
{

	int i, j, k, l;
	int index, indexCluster;
	i = patch_ID;
	
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

	
		InteriorMappingOnePatch(polycubePatch[i]);
		

	return true;

}

bool PhysicalPolycube::IsInteriorPoint(int vertexID)
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

bool PhysicalPolycube::InteriorMappingOnePatch(Polycube &currentPolycubePatch)
{

	int i, j, k, p, q;
	int index, indexPatch, iConst, nSize;
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

	indexPatch = currentPolycubePatch.indexPatch;  
	double theta=1.0;
	//iConst = patch_const[indexPatch];
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
			//printf("it is a plane\n");
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
		cout << endl << theta << endl;
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
          
			 polycubePara->vertex[InteriorList[j]][(iConst) % 3] = x_solvey(j);
			 //polycubePara->vertex[InteriorList[j]][iConst] =polycubePara->vertex[oneCorner][iConst]-fabs(polycubePara->vertex[InteriorList[j]][(iConst+2)%3]-polycubePara->vertex[oneCorner][(iConst+2)%3]);
         }


	}

	/////////////////////////////////////////////////////////////

	//MappingOnePatchWeighted(currentPolycubePatch);
	//MappingOnePatchPostProcessing(currentPolycubePatch);

	int stepNumber = 0;

	while(FlipCheckEachPatch(currentPolycubePatch) == false)
	{

		//cout<<"Polycube patch "<<currentPolycubePatch.index<<" has flipped elements!"<<"Cluster No. "<<currentPolycubePatch.indexPatch<<endl;
		
		if (stepNumber < 0)
		{
			MappingOnePatchWeighted(currentPolycubePatch);
		}
		if(fabs(1.0/sqrt(2.0)-theta)>EPSILON)
		{
		MappingOnePatchPostProcessing(currentPolycubePatch);
		}
		stepNumber++;

		//cout<<"Polycube patch "<<currentPolycubePatch.index<<" still has flipped elements after "<<stepNumber<<" iterations!"<<endl;

		if (stepNumber > 0)
		{
			
			//cout<<"Polycube patch "<<currentPolycubePatch.index<<" still has flipped elements after 1000 iterations!"<<endl;
			break;

		}

	}	

	//MappingOnePatchWeighted(currentPolycubePatch);
	///////////////////////////////////////////////////////////

	return true;

}

bool PhysicalPolycube::FlipCheckEachPatch(Polycube &currentPolycubePatch)
{

	int elementID, indexPatch;

	bool trueFalse = true;

	int i;

	indexPatch = currentPolycubePatch.indexPatch;

	for (i = 0; i < currentPolycubePatch.numElements; i++)
	{

		elementID = currentPolycubePatch.element[i];

		trueFalse = SignedTriAreaElement(elementID, indexPatch);

		if (trueFalse == false)
		{

			return false;

		}

	}

	return true;

}

bool PhysicalPolycube::SignedTriAreaElement(int elementID, int indexCluster)
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

		//cout<<"Error in SignedTriAreaElement()!!!"<<endl;
		return false;
		break;
		exit;

	}

	return true;

}

bool PhysicalPolycube::MappingOnePatchWeighted(Polycube &currentPolycubePatch)
{

	int i, j, k, p, q;
	int index, indexPatch, iConst, nSize;
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

	indexPatch = currentPolycubePatch.indexPatch;

	//iConst = patch_const[indexPatch];
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

double PhysicalPolycube::TriAreaElement(int elementID)
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

bool PhysicalPolycube::MappingOnePatchPostProcessing(Polycube &currentPolycubePatch)
{

	int i, j, k, p, q;
	int index, indexPatch, iConst, nSize;
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

	indexPatch = currentPolycubePatch.indexPatch;

	//iConst = patch_const[indexPatch];
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

bool PhysicalPolycube::InitializeOctree()
{

	voxelSize = 1 << OCTREE_MAX_LEVEL;
	gridSize = voxelSize + 1;
	numVoxels = voxelSize * voxelSize * voxelSize;
	numGrids = gridSize * gridSize * gridSize;

	

	for (int i = 0; i < 3; i++)
	{
		origCood[i] = 0.f;
	}

	//cellSize = 64.0 / voxelSize;
	cellSize = DOMAIN_SIZE / voxelSize;

	leafNum = 0;
	octreeDepth = GetDepth(voxelSize);
	octreeCellNum = GetOctreeNum(octreeDepth);

	refineFlagArray.resize(octreeCellNum);

	cutArray.resize(numVoxels);

	for (int i = 0; i <= octreeDepth; i++)
	{
		levelRes[i] = (1 << i);
	}


	//cout <<" cellSize:  "<< cellSize << endl;
	//cout <<" gridSize:  "<< gridSize << endl;;
	//cout <<" leafNum :  "<< leafNum << endl;;
	//cout <<" numVoxel:  "<< numVoxels << endl;;
	//cout <<" numGrids:  "<< numGrids << endl;;
	//cout <<" octreeDe:  "<< octreeDepth << endl;;
	//cout <<" octreeCe:  "<< octreeCellNum << endl;;
	//cout <<" voxelSiz:  "<< voxelSize << endl;;


	return true;

}

bool PhysicalPolycube::ConstructeOctree()
{

	CellQueue prev_queue, cur_queue;

	int octree_id, level;

	int start_level, end_level;

	vector<int> cell_sign;

	start_level = OCTREE_MIN_LEVEL;
	//end_level = OCTREE_MAX_LEVEL;
	end_level = OCTREE_MIN_LEVEL + 1;

	octree_id = 0;

	leafNum = 0;

	prev_queue.Add(octree_id);

	while (prev_queue.Empty() == 0)//while prev_queue is not empty
	{

		while (prev_queue.Get(octree_id) >= 0)
		{

			cell_sign.push_back(octree_id);
			level = GetLevel(octree_id);

			if (
				(((CheckAdaptation(octree_id) == false) && level < end_level) || level <= start_level)
				&& refineFlagArray[octree_id] != 1
				)

			{

				int octree_idx[8];
				if (octree_id == 0)
				{

					refineFlagArray[octree_id] = 1;
					cur_queue.Add(octree_id);

				}
				else
				{

					RefineBrothers(octree_id, octree_idx);

					for (int ii = 0; ii < 8; ii++)
					{
						cur_queue.Add(octree_idx[ii]);
					}

				}

			}

		}

		while (cur_queue.Get(octree_id) >= 0)
		{

			level = GetLevel(octree_id);
			for (int i = 0; i < 8; i++)
			{
				prev_queue.Add(Child(octree_id, level, i));
			}

		}

	}

	for (int i = 0; i < cell_sign.size(); i++)
	{

		if (refineFlagArray[cell_sign[i]] != 1)
		{

			cutArray[leafNum] = cell_sign[i];
			leafNum++;

		}

	}

	return true;

}

void PhysicalPolycube::RefineBrothers(int octree_id, int *octree_idx)
{

	int level, parent_level;
	int x, y, z;
	int parent_x, parent_y, parent_z;
	int parent_idx;

	level = GetLevel(octree_id);
	OctreeidxToXYZ(octree_id, x, y, z, level);
	parent_x = x / 2;
	parent_y = y / 2;
	parent_z = z / 2;
	parent_level = level - 1;
	//the usage of parent_idx?
	parent_idx = XYZToOctreeidx(parent_x, parent_y, parent_z, parent_level);

	octree_idx[0] = XYZToOctreeidx(parent_x * 2, parent_y * 2, parent_z * 2, parent_level + 1);
	octree_idx[1] = XYZToOctreeidx(parent_x * 2 + 1, parent_y * 2, parent_z * 2, parent_level + 1);
	octree_idx[2] = XYZToOctreeidx(parent_x * 2, parent_y * 2 + 1, parent_z * 2, parent_level + 1);
	octree_idx[3] = XYZToOctreeidx(parent_x * 2 + 1, parent_y * 2 + 1, parent_z * 2, parent_level + 1);
	octree_idx[4] = XYZToOctreeidx(parent_x * 2, parent_y * 2, parent_z * 2 + 1, parent_level + 1);
	octree_idx[5] = XYZToOctreeidx(parent_x * 2 + 1, parent_y * 2, parent_z * 2 + 1, parent_level + 1);
	octree_idx[6] = XYZToOctreeidx(parent_x * 2, parent_y * 2 + 1, parent_z * 2 + 1, parent_level + 1);
	octree_idx[7] = XYZToOctreeidx(parent_x * 2 + 1, parent_y * 2 + 1, parent_z * 2 + 1, parent_level + 1);

	for (int i = 0; i < 8; i++)
	{
		//octreeArray[octree_idx[i]].refine_flag = 1;
		refineFlagArray[octree_idx[i]] = 1;
	}

}

int PhysicalPolycube::IntersectAxis(double *dP1_, double *dP2_, double *dP3_, double *dObject_, double *dIntersectP_, char cAxis_/* ='x' */)
{

	double dCoff[3], dSign[3];
	double ERR = 1.0E-6;


	if (dP1_[0] < dObject_[0] && dP2_[0] < dObject_[0] && dP3_[0] < dObject_[0])
		return 0;
	if (!(cAxis_ == 'x') && !(cAxis_ == 'X'))
	{
		if (dP1_[0] > dObject_[0] && dP2_[0] > dObject_[0] && dP3_[0] > dObject_[0])
			return 0;
	}

	if (dP1_[1] < dObject_[1] && dP2_[1] < dObject_[1] && dP3_[1] < dObject_[1])
		return 0;
	if (!(cAxis_ == 'y') && !(cAxis_ == 'Y'))
	{
		if (dP1_[1] > dObject_[1] && dP2_[1] > dObject_[1] && dP3_[1] > dObject_[1])
			return 0;
	}

	if (dP1_[2] < dObject_[2] && dP2_[2] < dObject_[2] && dP3_[2] < dObject_[2])
		return 0;
	if (!(cAxis_ == 'z') && !(cAxis_ == 'Z'))
	{
		if (dP1_[2] > dObject_[2] && dP2_[2] > dObject_[2] && dP3_[2] > dObject_[2])
			return 0;
	}


	dCoff[0] = (dP2_[1] - dP1_[1])*(dP3_[2] - dP1_[2]) - (dP3_[1] - dP1_[1])*(dP2_[2] - dP1_[2]);
	dCoff[1] = (dP2_[2] - dP1_[2])*(dP3_[0] - dP1_[0]) - (dP3_[2] - dP1_[2])*(dP2_[0] - dP1_[0]);
	dCoff[2] = (dP2_[0] - dP1_[0])*(dP3_[1] - dP1_[1]) - (dP3_[0] - dP1_[0])*(dP2_[1] - dP1_[1]);

	if (cAxis_ == 'X' || cAxis_ == 'x')
	{
		if (fabs(dCoff[0]) < ERR)
			if (dP1_[0] == dObject_[0])
				return 5;
			else
				return 0;
		dIntersectP_[0] = -(dCoff[1] * (dObject_[1] - dP1_[1]) + dCoff[2] * (dObject_[2] - dP1_[2])) / dCoff[0] + dP1_[0];
		dIntersectP_[1] = dObject_[1];
		dIntersectP_[2] = dObject_[2];
		if (dIntersectP_[0] <= dObject_[0])
			return 0;
	}
	else if (cAxis_ == 'Y' || cAxis_ == 'y')
	{
		if (fabs(dCoff[1]) < ERR)
			if (dP1_[1] == dObject_[1])
				return 5;
			else
				return 0;
		dIntersectP_[0] = dObject_[0];
		dIntersectP_[1] = -(dCoff[0] * (dObject_[0] - dP1_[0]) + dCoff[2] * (dObject_[2] - dP1_[2])) / dCoff[1] + dP1_[1];
		dIntersectP_[2] = dObject_[2];
		if (dIntersectP_[1] <= dObject_[1])
			return 0;
	}
	else if (cAxis_ == 'Z' || cAxis_ == 'z')
	{
		if (fabs(dCoff[2]) < ERR)
			if (dP1_[2] == dObject_[2])
				return 5;
			else
				return 0;
		dIntersectP_[0] = dObject_[0];
		dIntersectP_[1] = dObject_[1];
		dIntersectP_[2] = -(dCoff[1] * (dObject_[1] - dP1_[1]) + dCoff[0] * (dObject_[0] - dP1_[0])) / dCoff[2] + dP1_[2];
		if (dIntersectP_[2] <= dObject_[2])
			return 0;
	}
	else
	{
		printf("Please indicates the axis direction.\n");
		return 0;
	}

	//if ((dIntersectP_[0] - dObject_[0]) < ERR && (dIntersectP_[1] - dObject_[1]) < ERR && (dIntersectP_[2] - dObject_[2]) < ERR)
	//	return 0;

	dSign[0] = (dP2_[1] - dIntersectP_[1])*(dP3_[2] - dIntersectP_[2]) - (dP3_[1] - dIntersectP_[1])*(dP2_[2] - dIntersectP_[2])
		+ (dP2_[2] - dIntersectP_[2])*(dP3_[0] - dIntersectP_[0]) - (dP3_[2] - dIntersectP_[2])*(dP2_[0] - dIntersectP_[0])
		+ (dP2_[0] - dIntersectP_[0])*(dP3_[1] - dIntersectP_[1]) - (dP3_[0] - dIntersectP_[0])*(dP2_[1] - dIntersectP_[1]);
	dSign[1] = (dP3_[1] - dIntersectP_[1])*(dP1_[2] - dIntersectP_[2]) - (dP1_[1] - dIntersectP_[1])*(dP3_[2] - dIntersectP_[2])
		+ (dP3_[2] - dIntersectP_[2])*(dP1_[0] - dIntersectP_[0]) - (dP1_[2] - dIntersectP_[2])*(dP3_[0] - dIntersectP_[0])
		+ (dP3_[0] - dIntersectP_[0])*(dP1_[1] - dIntersectP_[1]) - (dP1_[0] - dIntersectP_[0])*(dP3_[1] - dIntersectP_[1]);
	dSign[2] = (dP1_[1] - dIntersectP_[1])*(dP2_[2] - dIntersectP_[2]) - (dP2_[1] - dIntersectP_[1])*(dP1_[2] - dIntersectP_[2])
		+ (dP1_[2] - dIntersectP_[2])*(dP2_[0] - dIntersectP_[0]) - (dP2_[2] - dIntersectP_[2])*(dP1_[0] - dIntersectP_[0])
		+ (dP1_[0] - dIntersectP_[0])*(dP2_[1] - dIntersectP_[1]) - (dP2_[0] - dIntersectP_[0])*(dP1_[1] - dIntersectP_[1]);

	if ((dSign[0]<0 && dSign[1]<0 && dSign[2]<0) || (dSign[0]>0 && dSign[1]>0 && dSign[2]>0))
		return 5;
	else if (fabs(dSign[0]) < ERR || fabs(dSign[1]) < ERR || fabs(dSign[2]) < ERR)
	{
		if (fabs(dSign[0]) < ERR && fabs(dSign[1]) < ERR)
			return 3;
		else if (fabs(dSign[0]) < ERR && fabs(dSign[2]) < ERR)
			return 2;
		else if (fabs(dSign[1]) < ERR&& fabs(dSign[2]) < ERR)
			return 1;
		else if ((fabs(dSign[0]) < ERR && (dSign[1] * dSign[2])>0) || (fabs(dSign[1]) < ERR && (dSign[2] * dSign[0])>0) || (fabs(dSign[2]) < ERR && (dSign[1] * dSign[0])>0))
			return 4;
		else
			return 0;
	}
	else
		return 0;

}
int PhysicalPolycube::Child(int octree_id, int level, int i)
{
	int x, y, z;
	int ret_idx = 0;
	OctreeidxToXYZ(octree_id, x, y, z, level);

	switch (i)
	{
	case 0:
		ret_idx = XYZToOctreeidx(x * 2, y * 2, z * 2, level + 1);
		break;
	case 1:
		ret_idx = XYZToOctreeidx(x * 2 + 1, y * 2, z * 2, level + 1);
		break;
	case 2:
		ret_idx = XYZToOctreeidx(x * 2, y * 2 + 1, z * 2, level + 1);
		break;
	case 3:
		ret_idx = XYZToOctreeidx(x * 2 + 1, y * 2 + 1, z * 2, level + 1);
		break;
	case 4:
		ret_idx = XYZToOctreeidx(x * 2, y * 2, z * 2 + 1, level + 1);
		break;
	case 5:
		ret_idx = XYZToOctreeidx(x * 2 + 1, y * 2, z * 2 + 1, level + 1);
		break;
	case 6:
		ret_idx = XYZToOctreeidx(x * 2, y * 2 + 1, z * 2 + 1, level + 1);
		break;
	case 7:
		ret_idx = XYZToOctreeidx(x * 2 + 1, y * 2 + 1, z * 2 + 1, level + 1);
		break;
	}

	return ret_idx;

}


bool PhysicalPolycube::CheckAdaptation(int octree_id)
{


	return false;

}
bool PhysicalPolycube::HexMeshOctree()
{

	parametricHex.CreateNewMesh(parametricHex.HEXAHEDRON, leafNum * 8, leafNum);
	
	//cout << "leafNum"<<leafNum << endl;
	//getchar();

	int octree_id;
	int level;
	int cell_size;
	int x, y, z;

	for (int i = 0; i < leafNum * 8; i += 8)
	{
		//i=0,8,16,24,...
		octree_id = cutArray[(i / 8)];//be careful, understand why i/8
		//cout << octree_id << endl;
		level = GetLevel(octree_id);
		cell_size = voxelSize / (1 << level);
		OctreeidxToXYZ(octree_id, x, y, z, level);

		//parametricHex.vertex[i][0] = origCood[0] + x * cell_size * cellSize;
		//parametricHex.vertex[i][1] = origCood[1] + y * cell_size * cellSize;
		//parametricHex.vertex[i][2] = origCood[2] + z * cell_size * cellSize;
		//parametricHex.vertex[i + 1][0] = origCood[0] + (x + 1) * cell_size * cellSize;
		//parametricHex.vertex[i + 1][1] = origCood[1] + y * cell_size * cellSize;
		//parametricHex.vertex[i + 1][2] = origCood[2] + z * cell_size * cellSize;
		//parametricHex.vertex[i + 2][0] = origCood[0] + (x + 1) * cell_size * cellSize;
		//parametricHex.vertex[i + 2][1] = origCood[1] + y * cell_size * cellSize;
		//parametricHex.vertex[i + 2][2] = origCood[2] + (z + 1) * cell_size * cellSize;
		//parametricHex.vertex[i + 3][0] = origCood[0] + x * cell_size * cellSize;
		//parametricHex.vertex[i + 3][1] = origCood[1] + y * cell_size * cellSize;
		//parametricHex.vertex[i + 3][2] = origCood[2] + (z + 1) * cell_size * cellSize;
		//parametricHex.vertex[i + 4][0] = origCood[0] + x * cell_size * cellSize;
		//parametricHex.vertex[i + 4][1] = origCood[1] + (y + 1) * cell_size * cellSize;
		//parametricHex.vertex[i + 4][2] = origCood[2] + z * cell_size * cellSize;
		//parametricHex.vertex[i + 5][0] = origCood[0] + (x + 1) * cell_size * cellSize;
		//parametricHex.vertex[i + 5][1] = origCood[1] + (y + 1) * cell_size * cellSize;
		//parametricHex.vertex[i + 5][2] = origCood[2] + z * cell_size * cellSize;
		//parametricHex.vertex[i + 6][0] = origCood[0] + (x + 1) * cell_size * cellSize;
		//parametricHex.vertex[i + 6][1] = origCood[1] + (y + 1) * cell_size * cellSize;
		//parametricHex.vertex[i + 6][2] = origCood[2] + (z + 1) * cell_size * cellSize;
		//parametricHex.vertex[i + 7][0] = origCood[0] + x * cell_size * cellSize;
		//parametricHex.vertex[i + 7][1] = origCood[1] + (y + 1) * cell_size * cellSize;
		//parametricHex.vertex[i + 7][2] = origCood[2] + (z + 1) * cell_size * cellSize;

		parametricHex.vertex[i][0] = origCood[0] + x * cell_size * cellSize;
		parametricHex.vertex[i][1] = origCood[1] + y * cell_size * cellSize;
		parametricHex.vertex[i][2] = origCood[2] + z * cell_size * cellSize;

		parametricHex.vertex[i + 1][0] = origCood[0] + (x + 1) * cell_size * cellSize;
		parametricHex.vertex[i + 1][1] = origCood[1] + y * cell_size * cellSize;
		parametricHex.vertex[i + 1][2] = origCood[2] + z * cell_size * cellSize;

		parametricHex.vertex[i + 2][0] = origCood[0] + (x + 1) * cell_size * cellSize;
		parametricHex.vertex[i + 2][1] = origCood[1] + (y + 1) * cell_size * cellSize;
		parametricHex.vertex[i + 2][2] = origCood[2] + z * cell_size * cellSize;

		parametricHex.vertex[i + 3][0] = origCood[0] + x * cell_size * cellSize;
		parametricHex.vertex[i + 3][1] = origCood[1] + (y + 1) * cell_size * cellSize;
		parametricHex.vertex[i + 3][2] = origCood[2] + z * cell_size * cellSize;

		parametricHex.vertex[i + 4][0] = origCood[0] + x * cell_size * cellSize;
		parametricHex.vertex[i + 4][1] = origCood[1] + y * cell_size * cellSize;
		parametricHex.vertex[i + 4][2] = origCood[2] + (z + 1) * cell_size * cellSize;

		parametricHex.vertex[i + 5][0] = origCood[0] + (x + 1) * cell_size * cellSize;
		parametricHex.vertex[i + 5][1] = origCood[1] + y * cell_size * cellSize;
		parametricHex.vertex[i + 5][2] = origCood[2] + (z + 1) * cell_size * cellSize;

		parametricHex.vertex[i + 6][0] = origCood[0] + (x + 1) * cell_size * cellSize;
		parametricHex.vertex[i + 6][1] = origCood[1] + (y + 1) * cell_size * cellSize;
		parametricHex.vertex[i + 6][2] = origCood[2] + (z + 1) * cell_size * cellSize;

		parametricHex.vertex[i + 7][0] = origCood[0] + x * cell_size * cellSize;
		parametricHex.vertex[i + 7][1] = origCood[1] + (y + 1) * cell_size * cellSize;
		parametricHex.vertex[i + 7][2] = origCood[2] + (z + 1) * cell_size * cellSize;

	}

	for (int i = 0; i < leafNum; i++)
	{
		//parametricHex.element[i][0] = i * 8;
		//parametricHex.element[i][3] = i * 8 + 1;
		//parametricHex.element[i][2] = i * 8 + 2;
		//parametricHex.element[i][1] = i * 8 + 3;
		//parametricHex.element[i][4] = i * 8 + 4;
		//parametricHex.element[i][7] = i * 8 + 5;
		//parametricHex.element[i][6] = i * 8 + 6;
		//parametricHex.element[i][5] = i * 8 + 7;

		parametricHex.element[i][0] = i * 8;
		parametricHex.element[i][1] = i * 8 + 1;
		parametricHex.element[i][2] = i * 8 + 2;
		parametricHex.element[i][3] = i * 8 + 3;
		parametricHex.element[i][4] = i * 8 + 4;
		parametricHex.element[i][5] = i * 8 + 5;
		parametricHex.element[i][6] = i * 8 + 6;
		parametricHex.element[i][7] = i * 8 + 7;

		parametricHex.elementSign[i] = cutArray[i];

	}


	
	//getchar();
	return true;

}



bool PhysicalPolycube::MappingRealPara(string hexfilename, int hex_number)
{
	parametricHex.SetBoundaryVertexSign(1);

	parametricHex.InitiateEdgeValence();
	//////////////////////////////////////////////////////////////////////

	parametricHex.InitiateElementValence();


	

	realDomainHex = parametricHex;
	
	realDomainHex.SetBoundaryVertexSign(1);
	realDomainHex.InitiateEdgeValence();
	realDomainHex.InitiateElementValence();

	propElement.resize(parametricHex.vertexNumber);


	int para_coord_sys[8][3] =
	{
		{ 0	,0,	0  },
		{8	,0,	0  },
		{8	,8,	0  },
		{0	,8,	0  },
		{0	,0,	8  },
		{8	,0,	8  },
		{8	,8,	8  },
		{0	,8,	8  },
	};

	int para_index[6][4] =
	{
		{ 0	,3,2,1 },
		{ 0	,1,5,4 },
		{ 1,2,6,5 },
		{ 2	,3,7,6 },
		{ 3,0,4,7 },
		{ 4,5,6,7 },
	};

	int  EDGENODES[12][2] = { { 0, 1 },{ 1, 2 },{ 2, 3 },{ 0, 3 },{ 4, 5 },{ 5, 6 },{ 6, 7 },{ 4, 7 },{ 0, 4 },{ 1, 5 },{ 2, 6 },{ 3, 7 } };
	int i;

	vector<vector<double> > real_hex_total_coordinate;
	vector<vector<int> > real_hex_total_connection;
	vector<vector<int> > smooth_index;

	real_hex_coordinate.clear();
	vector<vector<vector<int >>> real_hex_connection;

	for (i = 0; i < polycube_structure_hex_.elementNumber; ++i)
	{

		cout << "loop: " << i << endl;
		polycubePara = NULL;

		surface_mesh_tri_.CopyMesh(polycubePara);

		vector<int> patch_list_mapping;

		//initial realDomainHex = parametricHex;
		for (int loopi = 0; loopi < realDomainHex.vertexNumber; loopi++)
		{
			for (int loopj = 0; loopj < 3; loopj++)
			{
				realDomainHex.vertex[loopi][loopj] = parametricHex.vertex[loopi][loopj];
			}
		}
		realDomainHex = parametricHex;

		//first assign corner
		unit_cube_index = { -1,-1, -1,-1, -1,-1, -1,-1 };
		//unit_cube_index(8, -1);
		for (int loopi = 0; loopi < 6; loopi++)
		{
			//find element matching information.
			vector<vector<int> > corner_coordinate(4);

			for (int loopj = 0; loopj < 4; loopj++)
			{
				int para_coord_sys_index = para_index[loopi][loopj];


				int temp_index_polycube = polycube_element_[i].face_id[loopi][loopj];//polycube coordinate hex
				for (int loop_unitcube = 0; loop_unitcube < parametricHex.vertexNumber; loop_unitcube++)
				{
					double sum_unit_search = 0;
					for (int loopk = 0; loopk < 3; loopk++)
					{
						sum_unit_search = sum_unit_search + pow(parametricHex.vertex[loop_unitcube][loopk] - para_coord_sys[para_index[loopi][loopj]][loopk], 2);
					}
					sum_unit_search = sqrt(sum_unit_search);
					if (sum_unit_search < EPSILON)
					{
						unit_cube_index[para_coord_sys_index] = loop_unitcube;
						break;

					}
				}

				for (int loopk = 0; loopk < 3; loopk++)
				{
					realDomainHex.vertex[unit_cube_index[para_coord_sys_index]][loopk] = polycube_structure_hex_.vertex[temp_index_polycube][loopk];
				}

			}
		}

		
		

		//unit_cube_index cornerpoint;
		/*for (int i_internal_temp = 0; i_internal_temp < 6; ++i_internal_temp)
		{
			for (int j_internal_temp = 0; j_internal_temp < 4; ++j_internal_temp)
			{
				std::cout << unit_cube_index[para_index[i_internal_temp][j_internal_temp]] << ", ";
			}
			std::cout << endl << endl;
		}*/
			

		//getchar();
		//create edge;
		int unit_edge[12][2];
		for (int loopi = 0; loopi < 12; loopi++)
		{
			for (int loopj = 0; loopj < 2; loopj++)
			{
				unit_edge[loopi][loopj] = unit_cube_index[EDGENODES[loopi][loopj]];
				//cout << unit_edge[loopi][loopj] << ", " ;
			}
			//cout <<  endl;
		}

		vector<vector <int>> interior_point_on_edge;
		//linear interpolation on edge.

		
		



		for (int loopj = 0; loopj < 12; loopj++)
		{
			vector <int> temp_interior_point_on_edge;
			for (int loopi = 0; loopi < parametricHex.vertexNumber; loopi++)
			{

				if ((parametricHex.vertexSign[loopi] == 1) && (!std::count(unit_cube_index.begin(), unit_cube_index.end(), loopi)))//cout means whether appear
				{



					Vector3f V0(parametricHex.vertex[loopi][0] - parametricHex.vertex[unit_edge[loopj][0]][0],
						parametricHex.vertex[loopi][1] - parametricHex.vertex[unit_edge[loopj][0]][1],
						parametricHex.vertex[loopi][2] - parametricHex.vertex[unit_edge[loopj][0]][2]);

					Vector3f V1(parametricHex.vertex[unit_edge[loopj][1]][0] - parametricHex.vertex[unit_edge[loopj][0]][0],
						parametricHex.vertex[unit_edge[loopj][1]][1] - parametricHex.vertex[unit_edge[loopj][0]][1],
						parametricHex.vertex[unit_edge[loopj][1]][2] - parametricHex.vertex[unit_edge[loopj][0]][2]);



					/*double x_segment = (parametricHex.vertex[loopi][0] - parametricHex.vertex[unit_edge[loopj][0]][0])
						/(parametricHex.vertex[unit_edge[loopj][1]][0] - parametricHex.vertex[unit_edge[loopj][0]][0]);
					double y_segment = (parametricHex.vertex[loopi][1] - parametricHex.vertex[unit_edge[loopj][0]][1])
						/(parametricHex.vertex[unit_edge[loopj][1]][1] - parametricHex.vertex[unit_edge[loopj][0]][1]);
					double z_segment = (parametricHex.vertex[loopi][2] - parametricHex.vertex[unit_edge[loopj][0]][2])
						/(parametricHex.vertex[unit_edge[loopj][1]][2] - parametricHex.vertex[unit_edge[loopj][0]][2]);
					cout << x_segment;*/
					if (abs(V0.cross(V1).norm()) < EPSILON)
					{
						double ratio_para_phys;
						double para_length;
						double x_segment = (parametricHex.vertex[unit_edge[loopj][1]][0] - parametricHex.vertex[unit_edge[loopj][0]][0]);
						double y_segment = (parametricHex.vertex[unit_edge[loopj][1]][1] - parametricHex.vertex[unit_edge[loopj][0]][1]);
						double z_segment = (parametricHex.vertex[unit_edge[loopj][1]][2] - parametricHex.vertex[unit_edge[loopj][0]][2]);
						if (abs(x_segment) > EPSILON)
						{
							para_length = x_segment;
							ratio_para_phys = (parametricHex.vertex[loopi][0] - parametricHex.vertex[unit_edge[loopj][0]][0])
								/ para_length;
						}
						else if (abs(y_segment) > EPSILON)
						{
							para_length = y_segment;
							ratio_para_phys = (parametricHex.vertex[loopi][1] - parametricHex.vertex[unit_edge[loopj][0]][1])
								/ para_length;
						}
						else if (abs(z_segment) > EPSILON)
						{
							para_length = z_segment;
							ratio_para_phys = (parametricHex.vertex[loopi][2] - parametricHex.vertex[unit_edge[loopj][0]][2])
								/ para_length;
						}
						else
						{
							cout << "two points too close" << endl;
						}


						double x_segment_physical = (realDomainHex.vertex[unit_edge[loopj][1]][0] - realDomainHex.vertex[unit_edge[loopj][0]][0]);
						double y_segment_physical = (realDomainHex.vertex[unit_edge[loopj][1]][1] - realDomainHex.vertex[unit_edge[loopj][0]][1]);
						double z_segment_physical = (realDomainHex.vertex[unit_edge[loopj][1]][2] - realDomainHex.vertex[unit_edge[loopj][0]][2]);


						realDomainHex.vertex[loopi][0] = realDomainHex.vertex[unit_edge[loopj][0]][0] + ratio_para_phys * x_segment_physical;
						realDomainHex.vertex[loopi][1] = realDomainHex.vertex[unit_edge[loopj][0]][1] + ratio_para_phys * y_segment_physical;
						realDomainHex.vertex[loopi][2] = realDomainHex.vertex[unit_edge[loopj][0]][2] + ratio_para_phys * z_segment_physical;

						temp_interior_point_on_edge.push_back(loopi);
						//cout << loopi << ", ";
					}
				}


			}
			//cout << endl;
			interior_point_on_edge.push_back(temp_interior_point_on_edge);
		}
		
		/*if (Debug_yu)
		{
			string new_string = std::string(4 - std::to_string(i).length(), '0') + std::to_string(i);

			string tempName = inputName + "_initial_" + std::to_string(OCTREE_MAX_LEVEL) + "_patch_hex." + new_string + ".vtk";

			realDomainHex.Write(tempName.c_str());
		}*/

		//for (int loopi_edge = 0; loopi_edge < interior_point_on_edge.size(); loopi_edge++)
		//{
		//	for (int loopj_edge = 0; loopj_edge < interior_point_on_edge[loopi_edge].size(); loopj_edge++)
		//	{
		//		cout << interior_point_on_edge[loopi_edge][loopj_edge] << ", ";
		//	}
		//	cout << endl << endl;
		//}

		////getchar();


		/*for (int loopi = 0; loopi < unit_cube_corner_index.size()-1; loopi++)
		{
			for (int loopj = loopi+1; loopj < unit_cube_corner_index.size(); loopj++)
			{
				int index_unit_1 = unit_cube_corner_index[loopi];
				int index_unit_2 = unit_cube_corner_index[loopj];

			}
		}*/

		//Next step interpolate on the edge. 

		/*if (i==4)
		{


		string test_name = inputName + "physical_test.vtk";

		realDomainHex.Write(test_name.c_str());
		}*/

		//
		//Project for boundary surface
		//find correspoinding face index;


		int boundary_stop = hex_number;  
		
		for (int loopi = 0; loopi < 6; loopi++)
		{
			if (-1 != polycube_element_[i].patch_sign[loopi])
			{
				//find element matching information.
				vector<vector<int> > corner_coordinate(4);
				for (int loopj = 0; loopj < 4; loopj++)
				{
					int temp_index_polycube = polycube_element_[i].face_id[loopi][loopj];

					for (int loopk = 0; loopk < 4; loopk++)
					{
						int temp_index_physical = polycubePatch[polycube_element_[i].patch_sign[loopi]].cornerPoint[loopk];
						double sum_lsq = 0;
						for (int loop_coordinate = 0; loop_coordinate < 3; loop_coordinate++)
						{
							sum_lsq = sum_lsq + pow(polycube_structure_hex_.vertex[temp_index_polycube][loop_coordinate]
								- vertex[temp_index_physical][loop_coordinate], 2);
						}
						sum_lsq = sqrt(sum_lsq);
						if (sum_lsq < EPSILON)
						{
							corner_coordinate[loopj].push_back(temp_index_physical);

							break;

						}
					}
					for (int loopk = 0; loopk < 3; loopk++)
					{
						corner_coordinate[loopj].push_back(para_coord_sys[para_index[loopi][loopj]][loopk]);
					}
				}

				//cout << sum_lsq << endl;
				//cout << corner_coordinate << endl;
				ParametricMapping(polycube_element_[i].patch_sign[loopi], corner_coordinate);


				patch_list_mapping.push_back(polycube_element_[i].patch_sign[loopi]);


				//HexMeshProjectionInterior();

				parametricHex.SetBoundaryVertexSign(0);

				for (int loopj = 0; loopj < parametricHex.vertexNumber; ++loopj)
				{
					if (equation_plane(corner_coordinate, parametricHex.vertex[loopj]) == 1)
					{
						parametricHex.vertexSign[loopj] = 2;
					}

				}



				HexMeshProjectionBoundary(polycube_element_[i].patch_sign[loopi]);

				//ToDo create hex


			}
			else
			{
				/*struct PropagationElement
				{
					int px;
					int nx;
					int py;
					int ny;
					int pz;
					int nz;
				};*/

			}
			//polycubePara->Write("mapping_tri.vtk");

		}

		
		/*if (i==4)
		{
			string test_name = inputName + "_2_test.vtk";

			realDomainHex.Write(test_name.c_str());
			getchar();
		}*/


		string tempName;
		//linear interpolation interiro surface
		//first separate inner surface with boundary
		//following must repeat other wise will wrong.
		parametricHex.SetBoundaryVertexSign(1);
		realDomainHex.SetBoundaryVertexSign(1);


		


		for (int loopi = 0; loopi < 6; loopi++)
		{
			if (-1 != polycube_element_[i].patch_sign[loopi])
			{
				//find element matching information.
				vector<vector<int> > corner_coordinate(4);
				for (int loopj = 0; loopj < 4; loopj++)
				{
					int temp_index_polycube = polycube_element_[i].face_id[loopi][loopj];

					for (int loopk = 0; loopk < 4; loopk++)
					{
						int temp_index_physical = polycubePatch[polycube_element_[i].patch_sign[loopi]].cornerPoint[loopk];
						double sum_lsq = 0;
						for (int loop_coordinate = 0; loop_coordinate < 3; loop_coordinate++)
						{
							sum_lsq = sum_lsq + pow(polycube_structure_hex_.vertex[temp_index_polycube][loop_coordinate]
								- vertex[temp_index_physical][loop_coordinate], 2);
						}
						sum_lsq = sqrt(sum_lsq);
						if (sum_lsq < EPSILON)
						{
							corner_coordinate[loopj].push_back(temp_index_physical);

							break;

						}
					}
					for (int loopk = 0; loopk < 3; loopk++)
					{
						corner_coordinate[loopj].push_back(para_coord_sys[para_index[loopi][loopj]][loopk]);
					}
				}


				for (int loopj = 0; loopj < parametricHex.vertexNumber; ++loopj)
				{
					if (equation_plane(corner_coordinate, parametricHex.vertex[loopj]) == 1)
					{
						parametricHex.vertexSign[loopj] = 2;
						realDomainHex.vertexSign[loopj] = 2;
						//cout << loopj << ", ";
					}

				}

				//cout <<endl<< endl;

			}
		}
		//now corner becomes 1, boundary surface not corner is 2.

		

		//getchar();
		/*if (i == 4)
		{
			for (int loopi = 0; loopi < parametricHex.vertexNumber; loopi++)
			{

				if (parametricHex.vertexSign[loopi] == 0)
				{

					cout << loopi << ", ";
				}


			}
			getchar();

			for (int loopi = 0; loopi < parametricHex.vertexNumber; loopi++)
			{

				if (parametricHex.vertexSign[loopi] == 1)
				{

					cout << loopi << ", ";
				}


			}
			getchar();

			for (int loopi = 0; loopi < parametricHex.vertexNumber; loopi++)
			{

				if (parametricHex.vertexSign[loopi] == 2)
				{

					cout << loopi << ", ";
				}


			}
			getchar();
		}

*/

		for (int loopi = 0; loopi < interior_point_on_edge.size(); ++loopi)
		{
			for (int loopj = 0; loopj < interior_point_on_edge[loopi].size(); ++loopj)
			{
				parametricHex.vertexSign[interior_point_on_edge[loopi][loopj]] = 2;
				realDomainHex.vertexSign[interior_point_on_edge[loopi][loopj]] = 2;
			}
		}



		/*tempName = inputName + "_" + std::to_string(OCTREE_MAX_LEVEL) + "_patch_hex.vtk";

		realDomainHex.Write(tempName.c_str());*/



		//getchar();

		//remove interior corner point.
		for (int loopi = 0; loopi < unit_cube_index.size(); loopi++)
		{
			parametricHex.vertexSign[unit_cube_index[loopi]] = 2;
			realDomainHex.vertexSign[unit_cube_index[loopi]] = 2;
		}
		

		for (int loopi = 0; loopi < parametricHex.vertexNumber; loopi++)
		{

			if (parametricHex.vertexSign[loopi] == 1)
			{
				realDomainHex.vertexSign[loopi] == 1;
				FindPropagationBoundElementsSurface(loopi);
				
				/*if (i==31)
				{
					cout << loopi << ", ";
				}*/
				
			}


		}
		
		/*if (i==31)
		{
			cout << endl << "testtest" << endl;
			getchar();
		}*/
		//getchar();
		/*cout << "test_now" << endl;
		getchar();*/
		//parametricHex.SetBoundaryVertexSign(1);
		for (int loopi = 0; loopi < parametricHex.vertexNumber; loopi++)
		{

			if (parametricHex.vertexSign[loopi] == 1)
			{

				/*for (int loopj = 0; loopj < 3; loopj++)
				{
				realDomainHex.vertex[loopi][loopj] = 0.f;
				}*/
				realDomainHex.vertexSign[loopi] == 1;
				
				//cout << loopi << ", ";
			}


		}


		
		//3


		realDomainHex.SmoothInnerSurface(10000);
		for (int loopi = 0; loopi < parametricHex.vertexNumber; loopi++)
		{

			if (parametricHex.vertexSign[loopi] == 1)
			{
				PropagationSurface(loopi);
				//cout << loopi << ", ";
			}


		}
		//realDomainHex.SmoothInnerSurface(10000);
		vector<int> smooth_index_temp;

		for (int loopi = 0; loopi < parametricHex.vertexNumber; loopi++)
		{
			smooth_index_temp.push_back(parametricHex.vertexSign[loopi]);

		}

		if (smooth_index_temp.size()!= parametricHex.vertexNumber)
		{
			cout << "define boundary type error";
		}

		smooth_index.push_back(smooth_index_temp);

		//realDomainHex.SmoothInnerSurface(10000);

		//if (i == 31)
		//{


		//	for (int loopi = 0; loopi < unit_cube_index.size(); loopi++)
		//	{
		//		//find element matching information.
		//		//vector<vector<int> > corner_coordinate(4);

		//		cout << unit_cube_index[loopi] << ", ";


		//		for (int loopk = 0; loopk < 3; loopk++)
		//		{
		//			cout << realDomainHex.vertex[unit_cube_index[loopi]][loopk] << ", ";
		//		}


		//		cout << endl;
		//	}
		//	getchar();
		//}
		//



		parametricHex.SetBoundaryVertexSign(1);

		//linear interiro point in the volume


		realDomainHex.SetBoundaryVertexSign(1);

		std::string new_string;

		/*realDomainHex.Smooth(10000);

		new_string = std::string(4 - std::to_string(i).length(), '0') + std::to_string(i);

		tempName = inputName + "_before_0_octree_" + std::to_string(OCTREE_MAX_LEVEL) + "_patch_hex." + new_string + ".vtk";

		realDomainHex.Write(tempName.c_str());*/

		for (int loopi = 0; loopi < parametricHex.vertexNumber; loopi++)
		{

			if (parametricHex.vertexSign[loopi] == 0)
			{
				FindPropagationBoundElements(loopi);

			}


		}

		HexMeshProjectionInterior();

		/*if (i == 14)
		{
			int array_temp[3] = { 20,23,40 };
			for (int loopi = 0; loopi < 3; loopi++) {
				cout << propElement[array_temp[loopi]].nx << ", ";
				cout << propElement[array_temp[loopi]].px << ", ";
				cout << propElement[array_temp[loopi]].ny << ", ";
				cout << propElement[array_temp[loopi]].py << ", ";
				cout << propElement[array_temp[loopi]].nz << ", ";
				cout << propElement[array_temp[loopi]].pz << ", ";
				cout << endl;
			}
			getchar();
		}*/
		//cout << "here 2" << endl;
		///Modified pxnx if the point is outside.
		
		//remove  value close to zero

		for (int loopi = 0; loopi < realDomainHex.vertexNumber; loopi++)
		{
			for (int loopj = 0; loopj < 3; loopj++)
			{
				if (abs(realDomainHex.vertex[loopi][loopj])<EPSILON)
				{
					realDomainHex.vertex[loopi][loopj] = 0.0;
				}
			}

			
		}

		
		if (Debug_yu)
		{
			new_string = std::string(4 - std::to_string(i).length(), '0') + std::to_string(i);

			tempName = inputName + "_before_octree_" + std::to_string(OCTREE_MAX_LEVEL) + "_patch_hex." + new_string + ".vtk";

			realDomainHex.Write(tempName.c_str());
		}

		

		
		vector<vector<double >> real_hex_coordinate_temp;
		for (int loopi = 0; loopi < realDomainHex.vertexNumber; loopi++)
		{
			vector<double> inner_temp;
			
			inner_temp.assign(realDomainHex.vertex[loopi], realDomainHex.vertex[loopi] + 3);
			real_hex_coordinate_temp.push_back(inner_temp);
		}
		real_hex_coordinate.push_back(real_hex_coordinate_temp);
		vector<vector<int >> real_hex_connection_temp;
		for (int loopi = 0; loopi < realDomainHex.elementNumber; loopi++)
		{
			vector<int> inner_temp_connection;
			inner_temp_connection.assign(realDomainHex.element[loopi], realDomainHex.element[loopi] + 8);


			real_hex_connection_temp.push_back(inner_temp_connection);
		}
		real_hex_connection.push_back(real_hex_connection_temp);

	}


	
	//////////////test



	if (Debug_yu == true)
	{
		cout << "run here MappingReakPara before merging surface" << endl;
	}

	if (smooth_index.size()!= polycube_structure_hex_.elementNumber)
	{
		cout << "error for polycube number." << endl;
	}
	

	//final step
	vector<vector<int> > match_surface_index_incl_bnd;
	match_surface_index_incl_bnd.resize(match_surface_index.size());

	

	for (int loop_i = 0; loop_i < match_surface_index.size(); loop_i++)
	{
		
			int temp_A = match_surface_index[loop_i][0] / 6;   // polycube element ID
			
			int temp_A_face = match_surface_index[loop_i][0] % 6;

			int temp_B = match_surface_index[loop_i][1] / 6;
			


			int temp_B_face = match_surface_index[loop_i][1] % 6;

			//int temp_A_face_odd_even = temp_A_face % 2;
			int temp_A_face_next;
			if (temp_A_face==0)
			{
				temp_A_face_next = 5;
			}
			else if (temp_A_face == 1)
			{
				temp_A_face_next = 3;
			}
			else if (temp_A_face == 2)
			{
				temp_A_face_next = 4;
			}
			else if (temp_A_face == 3)
			{
				temp_A_face_next = 1;
			}
			else if (temp_A_face == 4)
			{
				temp_A_face_next = 2;
			}
			else if (temp_A_face == 5)
			{
				temp_A_face_next = 0;
			}
			int boundary_surface_number_A = 0;
			for (int loop_internal = 0; loop_internal < polycube_element_[temp_A].boundary_sign.size(); loop_internal++)
			{
				if (loop_internal== temp_A_face)
				{
					continue;
				}

				if (loop_internal == temp_A_face_next)
				{
					continue;
				}


				if (polycube_element_[temp_A].boundary_sign[loop_internal]>0)
				{
					boundary_surface_number_A++;
				}
			}

			

			
			int temp_B_face_next;
			if (temp_B_face == 0)
			{
				temp_B_face_next = 5;
			}
			else if (temp_B_face == 1)
			{
				temp_B_face_next = 3;
			}
			else if (temp_B_face == 2)
			{
				temp_B_face_next = 4;
			}
			else if (temp_B_face == 3)
			{
				temp_B_face_next = 1;
			}
			else if (temp_B_face == 4)
			{
				temp_B_face_next = 2;
			}
			else if (temp_B_face == 5)
			{
				temp_B_face_next = 0;
			}
			int boundary_surface_number_B = 0;
			for (int loop_internal = 0; loop_internal < polycube_element_[temp_B].boundary_sign.size(); loop_internal++)
			{
				if (loop_internal == temp_B_face)
				{
					continue;
				}

				if (loop_internal == temp_B_face_next)
				{
					continue;
				}


				if (polycube_element_[temp_B].boundary_sign[loop_internal]>0)
				{
					boundary_surface_number_B++;
				}
			}
			match_surface_index_incl_bnd[loop_i].push_back(match_surface_index[loop_i][0]);
			match_surface_index_incl_bnd[loop_i].push_back(match_surface_index[loop_i][1]);
			match_surface_index_incl_bnd[loop_i].push_back(boundary_surface_number_A);
			match_surface_index_incl_bnd[loop_i].push_back(boundary_surface_number_B);
			match_surface_index_incl_bnd[loop_i].push_back(-1); //used to mark whether done or not.
		//cout << endl;
	}

	//unit_cube_index cornerpoint;
	
	unit_cube_index_vector.clear();

	/*for (int loopi = 0; loopi < match_surface_index_incl_bnd.size(); loopi++)
	{
		for (int loopj = 0; loopj < match_surface_index_incl_bnd[loopi].size(); loopj++)
		{
			cout << match_surface_index_incl_bnd[loopi][loopj] << ", ";
		}
		cout << endl;	
	}*/
	//getchar();
	
	
	for (int i_internal_temp = 0; i_internal_temp < 6; ++i_internal_temp)
	{
		vector <int > unit_cube_index_vector_temp;
		for (int j_internal_temp = 0; j_internal_temp < 4; ++j_internal_temp)
		{
			unit_cube_index_vector_temp.push_back(unit_cube_index[para_index[i_internal_temp][j_internal_temp]] );
		}
		unit_cube_index_vector.push_back(unit_cube_index_vector_temp);
	}

	int start_boundary_number = 4;
	bool verify_to_stop = false;

	int hex_element_ID_stop_merge = polycube_structure_hex_.elementNumber - hex_number;  //for hex dominant heli 5, another 3
	int hex_element_ID_stop_output = polycube_structure_hex_.elementNumber - hex_number;  //for hex dominant

	while (true)
	{

		for (int loopi = 0; loopi < match_surface_index_incl_bnd.size(); loopi++)
		{
			
			int temp_A = match_surface_index_incl_bnd[loopi][0] / 6;   // polycube element ID

			int temp_B = match_surface_index_incl_bnd[loopi][1] / 6;
			if (temp_A>polycube_structure_hex_.elementNumber - hex_element_ID_stop_merge)
			{
				match_surface_index_incl_bnd[loopi][4] = 1;
				continue;
			}
			if (temp_B>polycube_structure_hex_.elementNumber - hex_element_ID_stop_merge)
			{
				match_surface_index_incl_bnd[loopi][4] = 1;
				continue;
			}	

			if (match_surface_index_incl_bnd[loopi][2] >= match_surface_index_incl_bnd[loopi][3])
			{
				if (match_surface_index_incl_bnd[loopi][2] == start_boundary_number && match_surface_index_incl_bnd[loopi][4] < 0)
				{
					MergeInternalSurface(match_surface_index_incl_bnd[loopi][0], match_surface_index_incl_bnd[loopi][1]);
					match_surface_index_incl_bnd[loopi][4] = 1;
					//cout << loopi << endl;
				}
				//MergeInternalSurface(1, 23);  //copy 1 to 23
			}
			else
			{
				if (match_surface_index_incl_bnd[loopi][3] == start_boundary_number && match_surface_index_incl_bnd[loopi][4] < 0)
				{
					MergeInternalSurface(match_surface_index_incl_bnd[loopi][1], match_surface_index_incl_bnd[loopi][0]);
					match_surface_index_incl_bnd[loopi][4] = 1;
					//cout << loopi << endl;
				}
			}
			
		}
		start_boundary_number--;
		
		for (int loopi = 0; loopi < match_surface_index_incl_bnd.size(); loopi++)
		{
			if (match_surface_index_incl_bnd[loopi][4]==1)
			{
				verify_to_stop = true;
			}
			else
			{
				verify_to_stop = false;
				break;
			}
		}

		if (verify_to_stop)
		{
			break;
		}
		if (start_boundary_number<0)
		{
			cout << "must have problem in merging" << endl;
			break;
		}
	}
	
	//After merging
	//MergeInternalSurface(1, 23);  //copy 1 to 23

	int coordinate_index_offset;
	real_hex_total_coordinate.clear();
	real_hex_total_connection.clear();
	string tempName;

	if (Debug_yu == true)
	{
		cout << "run here MappingReakPara before edit volume linear interpolation" << endl;
	}

	for (int loopk = 0; loopk < real_hex_coordinate.size()- hex_element_ID_stop_output; loopk++)
	{
		for (int loopi = 0; loopi < real_hex_coordinate[loopk].size(); loopi++)
		{
			for (int loopj = 0; loopj < 3; loopj++)
			{
				realDomainHex.vertex[loopi][loopj] = real_hex_coordinate[loopk][loopi][loopj];
			}
		}

		//deal with interior surface.




		/*parametricHex.SetBoundaryVertexSign(1);
		realDomainHex.SetBoundaryVertexSign(1);

		for (int loopi = 0; loopi < smooth_index[loopk].size(); loopi++)
		{
			realDomainHex.vertexSign[loopi] = smooth_index[loopk][loopi];
			parametricHex.vertexSign[loopi] = smooth_index[loopk][loopi];
		}*/
		/*if (loopk == 24)
		{

			cout << loopk << " here "<<endl;
			for (int loopi = 0; loopi < parametricHex.vertexNumber; loopi++)
			{
				if (realDomainHex.vertexSign[loopi] == 1)
				{
					cout << loopi << ", ";
				}
			}

			cout << loopk << " here " << endl;
			for (int loopi = 0; loopi < parametricHex.vertexNumber; loopi++)
			{
				if (realDomainHex.vertexSign[loopi] == 2)
				{
					cout << loopi << ", ";
				}
			}
			cout << endl;
			cout << realDomainHex.vertex[49][0] << endl;
			realDomainHex.SmoothInnerSurface(4000);
			cout << realDomainHex.vertex[49][0] << endl;

			getchar();
		}
*/

		parametricHex.SetBoundaryVertexSign(1);
		realDomainHex.SetBoundaryVertexSign(1);

		for (int loopi = 0; loopi < realDomainHex.vertexNumber; loopi++)
		{

			if (parametricHex.vertexSign[loopi] == 0)
			{
				FindPropagationBoundElements(loopi);

			}


		}
		HexMeshProjectionInterior();

		std::array<int, 3> coord_index{ 0,1,2 };
		int loop_count = 0;
		while (true)
		{
			int show_problem_interior = -1;
			//cout << "Here===========" << show_problem_interior<< endl;;
			
			for (int loopi = 0; loopi < parametricHex.vertexNumber; loopi++)
			{

				if (parametricHex.vertexSign[loopi] != 0)
				{
					continue;
				}

				
				//outputCentroids<<i<<": "<<propElement[i].px<<" "<<propElement[i].nx<<" "<<propElement[i].py<<" "<<propElement[i].ny<<" "<<propElement[i].pz<<" "<<propElement[i].nz<<endl;


				show_problem_interior = inner_angle_between_lines(loopi);//-1 no problem, 0 is nx px, 1 is ny,py
																		 //cout << loopi << "  " << show_problem_interior << endl;
																		 //getchar();
				if (show_problem_interior >= 0)
				{
					//cout << "need to change Propagation direction!!";
					break;
				}


				//Propagation(i);

			}


			if (show_problem_interior >= 0)
			{
				loop_count++;
				if (loop_count >= 3)
				{
					break;
				}
				int temp_x = coord_index[0];


				int temp_y = coord_index[1];

				int temp_z = coord_index[2];

				coord_index[0] = temp_z;
				coord_index[1] = temp_x;
				coord_index[2] = temp_y;
				for (int loopi = 0; loopi < parametricHex.vertexNumber; loopi++)
				{

					if (parametricHex.vertexSign[loopi] != 0)
					{
						continue;
					}



					/*cout << "coord_index:  " << coord_index[0]
					<< coord_index[1]
					<< coord_index[2] << endl;*/

					Propagation_wthout_initial(loopi, coord_index);

					/*cout << "test___"<<loopi<<", " << propElement[loopi].nx
					<< ", " << propElement[loopi].px
					<< ", " << propElement[loopi].ny
					<< ", " << propElement[loopi].py
					<< ", " << propElement[loopi].nz
					<< ", " << propElement[loopi].pz<<endl;*/


					//getchar();
				}
				
			}
			else
			{
				break;
			}

			
		}
		std::string new_string;
		if (Debug_yu == true)
		{
			cout << loopk << endl;

			new_string = std::string(4 - std::to_string(loopk).length(), '0') + std::to_string(loopk);

			tempName = inputName + "_octree_" + std::to_string(OCTREE_MAX_LEVEL) + "_patch_hex." + new_string + ".vtk";

			realDomainHex.Write(tempName.c_str());
		}

		if (loop_count>=3)
		{
			cout << "block " << loopk << " or its neighbour has irregular shape which will reduce mesh quality"<< endl;
			/*realDomainHex.SetBoundaryVertexSign(1);
			realDomainHex.Smooth(100000);*/
			//for (int loopi_local = 0; loopi_local < realDomainHex.vertexNumber; loopi_local++)
			//{
			//	if (realDomainHex.vertexSign[loopi_local])
			//	{
			//		cout << loopi_local << ", " ;
			//	}
			//}

			//cout << endl;
			////realDomainHex.Smooth(100000);
			//getchar();
		}

		for (int loopi = 0; loopi < realDomainHex.vertexNumber; loopi++)
		{
			for (int loopj = 0; loopj < 3; loopj++)
			{
				real_hex_coordinate[loopk][loopi][loopj]= realDomainHex.vertex[loopi][loopj];
			}
		}
		//following is combination
		coordinate_index_offset = real_hex_total_coordinate.size();

		for (int loopi = 0; loopi < realDomainHex.vertexNumber; loopi++)
		{
			vector<double> inner_temp;
			inner_temp.assign(realDomainHex.vertex[loopi], realDomainHex.vertex[loopi] + 3);
			real_hex_total_coordinate.push_back(inner_temp);
		}
		
		for (int loopi = 0; loopi < realDomainHex.elementNumber; loopi++)
		{
			vector<int> inner_temp_connection;
			inner_temp_connection.assign(realDomainHex.element[loopi], realDomainHex.element[loopi] + 8);

			for (int loopj = 0; loopj < inner_temp_connection.size(); loopj++)
			{
				inner_temp_connection[loopj] = inner_temp_connection[loopj] + coordinate_index_offset;
			}


			real_hex_total_connection.push_back(inner_temp_connection);
		}
	}

	/*string test_name = inputName + "parametric_test.vtk";

	parametricHex.Write(test_name.c_str());*/

	nNode = real_hex_total_coordinate.size();
	nElement = real_hex_total_connection.size();
	RawMesh real_domain_hex_total;
	real_domain_hex_total.CreateNewMesh(real_domain_hex_total.HEXAHEDRON, nNode, nElement);

	for (int loopi = 0; loopi < nNode; loopi++)
	{

		for (int loopj = 0; loopj < 3; loopj++)
		{

			real_domain_hex_total.vertex[loopi][loopj] = real_hex_total_coordinate[loopi][loopj];

		}

	}

	for (int loopi = 0; loopi < nElement; loopi++)
	{

		for (int loopj = 0; loopj < 8; loopj++)
		{

			real_domain_hex_total.element[loopi][loopj] = real_hex_total_connection[loopi][loopj];

		}

	}

	real_domain_hex_total.InitiateElementValence();
	real_domain_hex_total.InitiateEdgeValence();
	real_domain_hex_total.DeleteUnusedVertex();
	real_domain_hex_total.DeleteDuplicatedPoint();


	//real_domain_hex_total.Smooth(10);

	tempName = hexfilename;

	real_domain_hex_total.Write(tempName.c_str());
	
	return true;
}

int PhysicalPolycube::MergeInternalSurface(int FromID, int ToID)
{
	//parametric_surface_total[loopi][loopj][loopk]   //connection unit cube with octree loopi faceID, loopj element
	//vector<vector <int >> unit_cube_index_vector;  //corner point.
	//vector<vector<vector<double >>> real_hex_coordinate;// each polycube coordinate first polycube element then matrix.

	
	int temp_A = FromID / 6;   // polycube element ID
	int temp_A_face = FromID % 6;
	int temp_B = ToID / 6;
	int temp_B_face = ToID % 6;
	

	
	double normal_face[6][3] =
	{
		{ 0	,0,-1 },
		{ 0	,-1,0 },
		{ 1,0,0 },
		{ 0,1,0 },
		{ -1,0,0 },
		{ 0,0,1 },
	};

	vector <int> vertex_index_from;
	
	for (int loopj = 0; loopj < parametric_surface_total[temp_A_face].size(); loopj++)
	{
		for (int loopk = 0; loopk < 4; loopk++)
		{
			vertex_index_from.push_back(parametric_surface_total[temp_A_face][loopj][loopk]);
		}
	}
	sort(vertex_index_from.begin(), vertex_index_from.end());
	vertex_index_from.erase(unique(vertex_index_from.begin(), vertex_index_from.end()), vertex_index_from.end());

	vector <int> vertex_index_to;

	for (int loopj = 0; loopj < parametric_surface_total[temp_B_face].size(); loopj++)
	{
		for (int loopk = 0; loopk < 4; loopk++)
		{
			vertex_index_to.push_back(parametric_surface_total[temp_B_face][loopj][loopk]);
		}
	}
	sort(vertex_index_to.begin(), vertex_index_to.end());
	vertex_index_to.erase(unique(vertex_index_to.begin(), vertex_index_to.end()), vertex_index_to.end());

	/*for (int loopi = 0; loopi < vertex_index_from.size(); loopi++)
	{
		cout << vertex_index_from[loopi] << ", ";
	}
	cout << endl;
	for (int loopi = 0; loopi < vertex_index_to.size(); loopi++)
	{
		cout << vertex_index_to[loopi] << ", ";
	}
	cout << endl;*/

	vector<vector <double >> vertex_change_from;
	for (int loopi = 0; loopi < parametricHex.vertexNumber; loopi++)
	{
		
		vector <double > vertex_change_from_temp;
		for (int loopj = 0; loopj < 3; loopj++)
		{

			vertex_change_from_temp.push_back(parametricHex.vertex[loopi][loopj]);

		}
		vertex_change_from.push_back(vertex_change_from_temp);

	}

	int index;
	double centeOne[3] = { 4,4,4 };

	Vector3d Norm0, Norm1;
	Matrix3d Rt;
	Vector3d V0, V1;

	Norm1 << normal_face[temp_A_face][0], normal_face[temp_A_face][1], normal_face[temp_A_face][2];
	Norm0 << normal_face[temp_B_face][0], normal_face[temp_B_face][1], normal_face[temp_B_face][2];
	Rt = rotation_matrix(Norm0, Norm1);
	for (int loopi = 0; loopi < vertex_index_from.size(); loopi++)
	{			
		index = vertex_index_from[loopi];
		
		for (int loopj = 0; loopj < 3; loopj++)
		{
			V0 << vertex_change_from[index][0] - centeOne[0],
				  vertex_change_from[index][1] - centeOne[1],
				  vertex_change_from[index][2] - centeOne[2];
			//cout << V0 << endl;

			V1 = Rt*V0;
			//cout << V1 << endl;
			vertex_change_from[index][0] = V1[0] + centeOne[0];
			vertex_change_from[index][1] = V1[1] + centeOne[1];
			vertex_change_from[index][2] = V1[2] + centeOne[2];

		}
	}

	//calculate surface center
	centeOne[0]	   =0.0;
	centeOne[1]	   =0.0;
	centeOne[2]	   =0.0;

	
	for (int loopj = 0; loopj < 4; loopj++)
	{
		index = unit_cube_index_vector[temp_A_face][loopj];

		
		
		centeOne[0]=centeOne[0]+vertex_change_from[index][0] ;
		centeOne[1]=centeOne[1]+vertex_change_from[index][1] ;
		centeOne[2]=centeOne[2]+vertex_change_from[index][2] ;
		
	}
	centeOne[0]=centeOne[0]/4.0;
	centeOne[1]=centeOne[1]/4.0;
	centeOne[2]=centeOne[2]/4.0;
	
	
	
	Norm1 << normal_face[temp_B_face][0], normal_face[temp_B_face][1], normal_face[temp_B_face][2];
	Norm0 = -1 * Norm1;
	Rt = rotation_matrix(Norm0, Norm1);
	for (int loopi = 0; loopi < vertex_index_from.size(); loopi++)
	{
		index = vertex_index_from[loopi];

		for (int loopj = 0; loopj < 3; loopj++)
		{
			V0 << vertex_change_from[index][0] - centeOne[0],
				vertex_change_from[index][1] - centeOne[1],
				vertex_change_from[index][2] - centeOne[2];
			//cout << V0 << endl;

			V1 = Rt * V0;
			//cout << V1 << endl;
			vertex_change_from[index][0] = V1[0] + centeOne[0];
			vertex_change_from[index][1] = V1[1] + centeOne[1];
			vertex_change_from[index][2] = V1[2] + centeOne[2];

		}
	}
	vector <vector <int >> interface_match_vector;
	int corner_index_yu;

	//following  loop
	while (true)
	{


		interface_match_vector.clear();
		for (int loopi = 0; loopi < vertex_index_from.size(); loopi++)
		{
			vector <int > interface_match_vector_temp;
			interface_match_vector_temp.resize(2);
			int index_from = vertex_index_from[loopi];
			if (index_from == unit_cube_index_vector[temp_A_face][0])
			{
				corner_index_yu = loopi;
			}
			for (int loopj = 0; loopj < vertex_index_to.size(); loopj++)
			{
				int index_to = vertex_index_to[loopj];
				double sum_lsq = 0;
				for (int loop_coordinate = 0; loop_coordinate < 3; loop_coordinate++)
				{
					sum_lsq = sum_lsq + pow(vertex_change_from[index_from][loop_coordinate]
						- parametricHex.vertex[index_to][loop_coordinate], 2);
				}
				sum_lsq = sqrt(sum_lsq);
				if (sum_lsq < EPSILON)
				{
					interface_match_vector_temp[0] = index_from;
					interface_match_vector_temp[1] = index_to;
					break;

				}
			}
			interface_match_vector.push_back(interface_match_vector_temp);
		}

		double sum_lsq_exterior = 0;
		for (int loop_coordinate = 0; loop_coordinate < 3; loop_coordinate++)
		{
			sum_lsq_exterior = sum_lsq_exterior + pow(real_hex_coordinate[temp_A][interface_match_vector[corner_index_yu][0]][loop_coordinate]
				- real_hex_coordinate[temp_B][interface_match_vector[corner_index_yu][1]][loop_coordinate], 2);
		}
		sum_lsq_exterior = sqrt(sum_lsq_exterior);
		if (sum_lsq_exterior < EPSILON)
		{
			
			break;

		}
		else
		{
			Norm1 << vertex_change_from[unit_cube_index_vector[temp_A_face][0]][0] - centeOne[0],
				     vertex_change_from[unit_cube_index_vector[temp_A_face][0]][1] - centeOne[1],
				     vertex_change_from[unit_cube_index_vector[temp_A_face][0]][2] - centeOne[2];

			Norm0 << vertex_change_from[unit_cube_index_vector[temp_A_face][1]][0] - centeOne[0],
				    vertex_change_from[unit_cube_index_vector[temp_A_face][1]][1] - centeOne[1],
				    vertex_change_from[unit_cube_index_vector[temp_A_face][1]][2] - centeOne[2];


			Rt = rotation_matrix(Norm0, Norm1);
			for (int loopi = 0; loopi < vertex_index_from.size(); loopi++)
			{
				index = vertex_index_from[loopi];

				for (int loopj = 0; loopj < 3; loopj++)
				{
					V0 << vertex_change_from[index][0] - centeOne[0],
						  vertex_change_from[index][1] - centeOne[1],
						  vertex_change_from[index][2] - centeOne[2];
					//cout << V0 << endl;

					V1 = Rt * V0;
					//cout << V1 << endl;
					vertex_change_from[index][0] = V1[0] + centeOne[0];
					vertex_change_from[index][1] = V1[1] + centeOne[1];
					vertex_change_from[index][2] = V1[2] + centeOne[2];

				}
			}
		}
	}


	//cout << "the first corner is at rows: " << corner_index_yu << endl;
	//cout << temp_B << " " << temp_A << endl;
	for (int loopi = 0; loopi < interface_match_vector.size(); loopi++)
	{
		for (int loopj = 0; loopj < 3; loopj++)
		{
			
			real_hex_coordinate[temp_B][interface_match_vector[loopi][1]][loopj] = real_hex_coordinate[temp_A][interface_match_vector[loopi][0]][loopj];

		}
		
	}


	//cout << "test" << endl;
	
	/*RawMesh Hex_check;
	Hex_check.CreateNewMesh(Hex_check.HEXAHEDRON, parametricHex.vertexNumber, parametricHex.elementNumber);

	for (int i = 0; i < parametricHex.vertexNumber; i++)
	{

		for (int j = 0; j < 3; j++)
		{

			Hex_check.vertex[i][j] = vertex_change_from[i][j];

		}

	}

	for (int i = 0; i < parametricHex.elementNumber; i++)
	{

		for (int j = 0; j < 3; j++)
		{

			Hex_check.element[i][j] = parametricHex.element[i][j];

		}

	}


	Hex_check.Write("test_hex.vtk");*/

	//temp_change_parametricHex






	//getchar();
	return 1;
}

Matrix3d PhysicalPolycube::rotation_matrix(Vector3d vector0, Vector3d vector1)
{
	return Matrix3d(Quaterniond::FromTwoVectors(vector0, vector1));
}

int PhysicalPolycube::inner_angle_between_lines(int vertexID)
{
	int tempIndexOne, tempIndexTwo;
	double angle;
	double pi_value = 3.14159265358979323846;
	int return_number = -1;
	//			cout << vertexID << ": ";
	Eigen::Vector3d a_line;
	Eigen::Vector3d b_line;

	tempIndexOne = propElement[vertexID].nx;
	tempIndexTwo = propElement[vertexID].px;
	
	a_line[0] = realDomainHex.vertex[tempIndexOne][0] - realDomainHex.vertex[vertexID][0];
	a_line[1] = realDomainHex.vertex[tempIndexOne][1] - realDomainHex.vertex[vertexID][1];
	a_line[2] = realDomainHex.vertex[tempIndexOne][2] - realDomainHex.vertex[vertexID][2];
	b_line[0] = realDomainHex.vertex[tempIndexTwo][0] - realDomainHex.vertex[vertexID][0];
	b_line[1] = realDomainHex.vertex[tempIndexTwo][1] - realDomainHex.vertex[vertexID][1];
	b_line[2] = realDomainHex.vertex[tempIndexTwo][2] - realDomainHex.vertex[vertexID][2];



	//cout << "debug:realDomainHex.vertex   " << realDomainHex.vertex[vertexID][0]<<" "
	// <<realDomainHex.vertex[vertexID][1]<<" "
	// <<realDomainHex.vertex[vertexID][2]<<endl;
	/*cout<<  a_line[0] << endl;
	cout<<a_line[1] <<endl;
	cout<<a_line[2] <<endl;
	cout<<b_line[0] <<endl;
	cout<<b_line[1] <<endl;
	cout<<b_line[2] <<endl;*/
	//cout<< "debug:realDomainHex.nx px   " << realDomainHex.vertex[tempIndexOne][0]	   
	//<<	realDomainHex.vertex[tempIndexOne][1]  <<" "
	//<<	realDomainHex.vertex[tempIndexOne][2]  <<" "
	//<<	realDomainHex.vertex[tempIndexTwo][0]  <<" "
	//<<	realDomainHex.vertex[tempIndexTwo][1]  <<" "
	//<<	realDomainHex.vertex[tempIndexTwo][2]  <<" "
	// << "   debug" << endl;


	angle = std::atan2(a_line.cross(b_line).norm(), a_line.dot(b_line));
	//cout << angle << " " << endl;

	if (angle<pi_value/2)
	{
		return_number = 0;
	}

	tempIndexOne = propElement[vertexID].ny;
	tempIndexTwo = propElement[vertexID].py;

	a_line[0] = realDomainHex.vertex[tempIndexOne][0] - realDomainHex.vertex[vertexID][0];
	a_line[1] = realDomainHex.vertex[tempIndexOne][1] - realDomainHex.vertex[vertexID][1];
	a_line[2] = realDomainHex.vertex[tempIndexOne][2] - realDomainHex.vertex[vertexID][2];
	b_line[0] = realDomainHex.vertex[tempIndexTwo][0] - realDomainHex.vertex[vertexID][0];
	b_line[1] = realDomainHex.vertex[tempIndexTwo][1] - realDomainHex.vertex[vertexID][1];
	b_line[2] = realDomainHex.vertex[tempIndexTwo][2] - realDomainHex.vertex[vertexID][2];

	//cout << "debug:realDomainHex.ny py   " << realDomainHex.vertex[tempIndexOne][0]
	//	<< realDomainHex.vertex[tempIndexOne][1] << " "
	//	<< realDomainHex.vertex[tempIndexOne][2] << " "
	//	<< realDomainHex.vertex[tempIndexTwo][0] << " "
	//	<< realDomainHex.vertex[tempIndexTwo][1] << " "
	//	<< realDomainHex.vertex[tempIndexTwo][2] << " "
	//	<< "   debug" << endl;
	//
	angle = std::atan2(a_line.cross(b_line).norm(), a_line.dot(b_line));
	//cout << angle << " " << endl;
	if (angle<pi_value / 2)
	{
		return_number = 1;
	}

	tempIndexOne = propElement[vertexID].nz;
	tempIndexTwo = propElement[vertexID].pz;

	a_line[0] = realDomainHex.vertex[tempIndexOne][0] - realDomainHex.vertex[vertexID][0];
	a_line[1] = realDomainHex.vertex[tempIndexOne][1] - realDomainHex.vertex[vertexID][1];
	a_line[2] = realDomainHex.vertex[tempIndexOne][2] - realDomainHex.vertex[vertexID][2];
	b_line[0] = realDomainHex.vertex[tempIndexTwo][0] - realDomainHex.vertex[vertexID][0];
	b_line[1] = realDomainHex.vertex[tempIndexTwo][1] - realDomainHex.vertex[vertexID][1];
	b_line[2] = realDomainHex.vertex[tempIndexTwo][2] - realDomainHex.vertex[vertexID][2];

	//cout << "debug:realDomainHex.nz pz   " << realDomainHex.vertex[tempIndexOne][0]
	//	<< realDomainHex.vertex[tempIndexOne][1] << " "
	//	<< realDomainHex.vertex[tempIndexOne][2] << " "
	//	<< realDomainHex.vertex[tempIndexTwo][0] << " "
	//	<< realDomainHex.vertex[tempIndexTwo][1] << " "
	//	<< realDomainHex.vertex[tempIndexTwo][2] << " "
	//	<< "   debug" << endl;

	angle = std::atan2(a_line.cross(b_line).norm(), a_line.dot(b_line));
	//cout << angle << " "<<endl;

	if (angle<pi_value / 2)
	{
		return_number = 2;
	}

	//cout << "debug:realDomainHex.vertex   " << realDomainHex.vertex[vertexID][0] << " "
	//	<< realDomainHex.vertex[vertexID][1] << " "
	//	<< realDomainHex.vertex[vertexID][2] << endl;
	//cout << endl;

	return return_number;
}

int PhysicalPolycube::equation_plane(vector<vector<int> > &corner_coordinate, double point[3])
{
	double x1 = corner_coordinate[0][1];
	double y1 = corner_coordinate[0][2];
	double z1 = corner_coordinate[0][3];
	double x2 = corner_coordinate[1][1];
	double y2 = corner_coordinate[1][2];
	double z2 = corner_coordinate[1][3];
	double x3 = corner_coordinate[2][1];
	double y3 = corner_coordinate[2][2];
	double z3 = corner_coordinate[2][3];
	double x = point[0];
	double y = point[1];
	double z = point[2];


	double a1 = x2 - x1;
	double b1 = y2 - y1;
	double c1 = z2 - z1;
	double a2 = x3 - x1;
	double b2 = y3 - y1;
	double c2 = z3 - z1;
	double a = b1 * c2 - b2 * c1;
	double b = a2 * c1 - a1 * c2;
	double c = a1 * b2 - b1 * a2;
	double d = (-a * x1 - b * y1 - c * z1);

	// equation of plane is: a*x + b*y + c*z = 0 #  

	// checking if the 4th point satisfies  
	// the above equation  
	if (abs(a * x + b * y + c * z + d - 0)<EPSILON)
		return 1;
	else
		return 0;

}

bool PhysicalPolycube::CleanMesh(const char *outputName)
{

	cout << "***** Start removing outside elements! *****" << endl;

	int i, j, k;

	bool *outVertex = new bool[parametricHex.vertexNumber];

	polycubePara->InitiateElementValence();

#pragma omp parallel for
	for (i = 0; i < parametricHex.elementNumber; i++)
	{

		double elementCenter[3] = { 0.f };

		for (j = 0; j < 8; j++)
		{

			for (k = 0; k < 3; k++)
			{

				elementCenter[k] += parametricHex.vertex[parametricHex.element[i][j]][k] / 8.0f;

			}

		}

		if (IsOutOfSurface(elementCenter, polycubePara[0])) //polycubePara[0], not polycubePara.
		{

			for (j = 0; j < 8; j++)
			{

				int index = parametricHex.element[i][j];

				outVertex[index] = true;

			}

		}
		else
		{

			for (j = 0; j < 8; j++)
			{

				int index = parametricHex.element[i][j];

				outVertex[index] = false;

			}

		}

	}

	bool *badElement = new bool[parametricHex.elementNumber];

#pragma omp parallel for private(j, k)
	for (i = 0; i < parametricHex.elementNumber; i++)
	{

		badElement[i] = false;

		for (j = 0; j < parametricHex.elementProperty.vertexNumber; ++j)
		{
			if (parametricHex.element[i][j] == -1)
				continue;

			if (outVertex[parametricHex.element[i][j]])
			{

				badElement[i] = true;

				for (k = 0; k < parametricHex.elementProperty.vertexNumber; ++k)
				{
					if (parametricHex.element[i][k] == -1)
						continue;
					parametricHex.vertexSign[parametricHex.element[i][k]] = 1;
				}

				break;

			}
		}

	}

	parametricHex.DeleteElement(badElement);

	parametricHex.DeleteUnusedVertex();
	parametricHex.DeleteDuplicatedPoint();

	//parametricHex.Write(outputName);

	return true;

}


bool PhysicalPolycube::IsOutOfSurface(double dP_[3], RawMesh &surfaceMesh, int iTime)
{

	if (iTime > 10)
		return 0;

	int i, iIns;
	double dInsP[3];
	bool atPoint;
	int dCount;

	//double ERR = 0.0013456;

	dCount = 0;
	atPoint = false;

	for (i = 0; i<surfaceMesh.elementNumber; i++)
	{
		iIns = IntersectAxis(surfaceMesh.vertex[surfaceMesh.element[i][0]], surfaceMesh.vertex[surfaceMesh.element[i][1]], surfaceMesh.vertex[surfaceMesh.element[i][2]], dP_, dInsP, 'x');

		if (iIns == 0)
			continue;
		else if (iIns == 5)
			dCount += 2;
		else if (iIns == 4)
			dCount += 1;
		else
			atPoint = true;
		/*
		{
		if (intersectAxis(rf.vertex[rf.element[i][0]], rf.vertex[rf.element[i][1]], rf.vertex[rf.element[i][2]], dP_, dInsP, 'x')>0)
		dCount += 1.0/rf.elementValenceNum[rf.element[i][-iIns]];
		else if (intersectAxis(rf.vertex[rf.element[i][0]], rf.vertex[rf.element[i][1]], rf.vertex[rf.element[i][2]], dP_, dInsP, 'x')>0)
		dCount += 1.0/rf.elementValenceNum[rf.element[i][-iIns]];
		}
		//*/
	}

	if (atPoint)
	{
		double dTP[3];
		bool ii, ij;
		iTime++;
		//printf("Times: %d\n", iTime);
		dTP[0] = dP_[0] + (rand() % 10)*GM.ERR; dTP[1] = dP_[1] + (rand() % 10)*GM.ERR; dTP[2] = dP_[2] + (rand() % 10)*GM.ERR;
		ii = IsOutOfSurface(dTP, surfaceMesh, iTime);
		dTP[0] = dP_[0] - (rand() % 10)*GM.ERR; dTP[1] = dP_[1] - (rand() % 10)*GM.ERR; dTP[2] = dP_[2] - (rand() % 10)*GM.ERR;
		ij = IsOutOfSurface(dTP, surfaceMesh, iTime);
		if (ii != ij)
			return false;
		else
			return ii;
	}

	if (((int)dCount) % 4 == 0)
		return  true;
	else
		return false;

}

bool PhysicalPolycube::HexMeshParametricDomain()
{
	
	InitializeOctree();//For each cube, define

	ConstructeOctree(); //connectivity

	HexMeshOctree(); // coordinate

	


	/*cout << "test======" << endl;
	cout << parametricHex.vertex[1][0]<< endl;
	cout << realDomainHex.vertex[1][0]<< endl;

	realDomainHex.vertex[1][0] = 1111;

	cout << parametricHex.vertex[1][0] << endl;
	cout << realDomainHex.vertex[1][0] << endl;*/
	//cout << "test end======" << endl;
	//Above is to create hex mesh for unit cube. 





	//To do , 

	//To do
	//CleanMesh(outputName);

	string file_name = "octree_another" + std::to_string(OCTREE_MAX_LEVEL) + "_test_hex.vtk";
	parametricHex.DeleteDuplicatedPoint();
	parametricHex.InitiateElementValence();
	parametricHex.InitiateEdgeValence();

	//parametricHex.Write(file_name.c_str());

	/*RawMesh *parametric_surface_hex;
	parametric_surface_hex = parametricHex.ExtractSurface();

	string tempName = "parametric_extract_surface_quad.raw";
	parametric_surface_hex->Write(tempName.c_str());*/

	int i, j, k, faceVertexNumber;
	int *vertexIndex = NULL;
	InitiateMatrix(vertexIndex, 4);
	
	
	
	parametric_surface_total.resize(6);
	for ( i = 0; i<parametricHex.elementNumber; ++i)
	{
		//if (elementSign[i] < 0)
		//	continue;

		for (j = 0; j<6; ++j)
		{
			parametricHex.GetElementFace(i, j, vertexIndex);
			if (parametricHex.IsBoundaryFace(vertexIndex))
			{
				vector <int> parametric_surface_temp;
				for (k = 0; k < 4; ++k)
				{
					parametric_surface_temp.push_back(vertexIndex[k]);
				}
				parametric_surface_total[j].push_back(parametric_surface_temp);
			}
		}
	}
	/*for (int loopi = 0; loopi < 6; loopi++)
	{
		for (int loopj = 0; loopj < parametric_surface_total[loopi].size(); loopj++)
		{
			for (int loopk = 0; loopk < 4; loopk++)
			{
				cout << parametric_surface_total[loopi][loopj][loopk] << ", ";
			}
		}
		cout << endl << endl;
	}*/
	//getchar();
	return true;

}




// For hex mesh generation in real domain
//bool PhysicalPolycube::HexMeshRealDomain(const char *outputName)
//{
//
//	parametricHex.SetBoundaryVertexSign(1);
//
//	//////////////////////////////////////////////////////////////////////
//	////Useful for propagation step
//	parametricHex.InitiateEdgeValence();
//	//////////////////////////////////////////////////////////////////////
//
//	parametricHex.InitiateElementValence();
//
//	//CheckUniformHexValidity(); // Not valid for unusual polycube structures
//
//	realDomainHex = parametricHex;
//
//
//
//
//	if (READ_IN_UNIFORMHEX == 0)
//	{
//
//		realDomainHex.SetBoundaryVertexSign(1);
//		realDomainHex.InitiateEdgeValence();
//		realDomainHex.InitiateElementValence();
//
//		propElement.resize(parametricHex.vertexNumber);
//
//		HexMeshProjectionBoundary();
//
//		HexMeshProjectionInterior();
//
//		//Pillowing();
//
//		//realDomainHex.Write("PillowNoSmooth_hex.raw");
//		string tempName;	
//		realDomainHex.Write(outputName);
//		tempName = inputName + "_realHex_hex.inp";
//
//		realDomainHex.Write(tempName.c_str());
//
//		realDomainHex.Smooth(10000);
//
//
//
//		tempName = inputName + "_NoPillowSmooth_hex.raw";
//
//		realDomainHex.Write(tempName.c_str());
//
//		tempName = inputName + "_NoPillowSmooth_hex.vtk";
//
//		realDomainHex.Write(tempName.c_str());
//		tempName = inputName + "_NoPillowSmooth_hex.inp";
//
//		realDomainHex.Write(tempName.c_str());
//
//	}
//	else
//	{
//
//		realDomainHex.Read("uniform_hex.raw");
//		realDomainHex.SetBoundaryVertexSign(1);
//		realDomainHex.InitiateEdgeValence();
//		realDomainHex.InitiateElementValence();
//
//		for (int i = 0; i < realDomainHex.elementNumber; i++)
//		{
//			realDomainHex.elementSign[i] = parametricHex.elementSign[i];
//		}
//		realDomainHex.Smooth(10000);
//		string tempName;	
//		tempName = inputName + "_NoPillowSmooth_hex_NEW.raw";
//
//		realDomainHex.Write(tempName.c_str());
//
//		tempName = inputName + "_NoPillowSmooth_hex_NEW.vtk";
//
//		realDomainHex.Write(tempName.c_str());
//		tempName = inputName + "_NoPillowSmooth_hex_NEW.inp";
//
//		realDomainHex.Write(tempName.c_str());
//
//	}
//
//
//	return true;
//
//}

bool PhysicalPolycube::HexMeshProjectionBoundary(int patch_ID)
{

	int i, j;

	double paraPos[3], realPos[3];

	for (i = 0; i < parametricHex.vertexNumber; i++)
	{

		if (parametricHex.vertexSign[i] == 2)
		{

			for (j = 0; j < 3; j++)
			{

				paraPos[j] = parametricHex.vertex[i][j];

			}
			//cout << parametricHex.vertexSign[i] << endl;
			Projection(paraPos, realPos, patch_ID);

			for (j = 0; j < 3; j++)
			{

				realDomainHex.vertex[i][j] = realPos[j];

			}

		}
		

	}

	return true;

}
bool PhysicalPolycube::FindPropagationBoundElementsSurface(int vertexID)
{

	int i, j;
	int nVert;
	int five_direction_on_surface_flag[6]= { 0, 0, 0, 0, 0, 0 };


	for (i = 0; i < parametricHex.edgeValenceNumber[vertexID]; i++)
	{
		nVert = parametricHex.edgeValence[vertexID][i];

		if ((parametricHex.vertex[nVert][0] - parametricHex.vertex[vertexID][0]) > EPSILON)
		{
			propElement[vertexID].px = nVert; 
			five_direction_on_surface_flag[0] = 1;
		}
		else if ((parametricHex.vertex[nVert][0] - parametricHex.vertex[vertexID][0]) < -EPSILON)
		{
			propElement[vertexID].nx = nVert;
			five_direction_on_surface_flag[1] = 1;

		}
		else if ((parametricHex.vertex[nVert][1] - parametricHex.vertex[vertexID][1]) > EPSILON)
		{
			propElement[vertexID].py = nVert;
			five_direction_on_surface_flag[2] = 1;

		}
		else if ((parametricHex.vertex[nVert][1] - parametricHex.vertex[vertexID][1]) < -EPSILON)
		{
			propElement[vertexID].ny = nVert;
			five_direction_on_surface_flag[3] = 1;

		}
		else if ((parametricHex.vertex[nVert][2] - parametricHex.vertex[vertexID][2]) > EPSILON)
		{
			propElement[vertexID].pz = nVert;
			five_direction_on_surface_flag[4] = 1;

		}
		else
		{
			propElement[vertexID].nz = nVert;
			five_direction_on_surface_flag[5] = 1;
		}

	}

	if (0 == five_direction_on_surface_flag[0])
	{
		propElement[vertexID].px = vertexID;
	}
	else if (0 == five_direction_on_surface_flag[1])
	{
		propElement[vertexID].nx = vertexID;
	}
	else if (five_direction_on_surface_flag[2] == 0)
	{
		propElement[vertexID].py = vertexID;
	}
	else if (five_direction_on_surface_flag[3] == 0)
	{
		propElement[vertexID].ny = vertexID;

	}
	else if (five_direction_on_surface_flag[4] == 0)
	{
		propElement[vertexID].pz = vertexID;


	}
	else if (five_direction_on_surface_flag[5] == 0)
	{
		propElement[vertexID].nz = vertexID;
	}

	/*cout << "testtest" << endl;
	for (int loopi = 0; loopi < 6; loopi++)
	{
		cout << five_direction_on_surface_flag[loopi] << endl;

	}*/

	int index;
	int jump_out_endless_loop = 0;
	//cout << numGrids << endl;
	while ((parametricHex.vertexSign[propElement[vertexID].px]!=2)&& (five_direction_on_surface_flag[0] != 0))
	{
		jump_out_endless_loop++;

		//To do future, if there is interior line, need to rewrite.
		if (jump_out_endless_loop> numGrids)
		{
			propElement[vertexID].px = vertexID;
			break;
		}

		index = propElement[vertexID].px;

		for (j = 0; j < parametricHex.edgeValenceNumber[index]; j++)
		{
			nVert = parametricHex.edgeValence[index][j];

			if ((parametricHex.vertex[nVert][0] - parametricHex.vertex[index][0]) > EPSILON)
			{

				propElement[vertexID].px = nVert;

				break;

			}
		}

	}
	jump_out_endless_loop = 0;

	while ((parametricHex.vertexSign[propElement[vertexID].nx]!= 2) && (five_direction_on_surface_flag[1] != 0))
	{
		jump_out_endless_loop++;

		//To do future, if there is interior line, need to rewrite.
		if (jump_out_endless_loop> numGrids)
		{
			propElement[vertexID].nx = vertexID;
			break;
		}

		index = propElement[vertexID].nx;

		for (j = 0; j < parametricHex.edgeValenceNumber[index]; j++)
		{
			nVert = parametricHex.edgeValence[index][j];

			if ((parametricHex.vertex[nVert][0] - parametricHex.vertex[index][0]) < -EPSILON)
			{

				propElement[vertexID].nx = nVert;

				break;

			}
		}

	}
	jump_out_endless_loop = 0;

	while ((parametricHex.vertexSign[propElement[vertexID].py]!= 2) && (five_direction_on_surface_flag[2] != 0))
	{
		jump_out_endless_loop++;

		//To do future, if there is interior line, need to rewrite.
		if (jump_out_endless_loop> numGrids)
		{
			propElement[vertexID].py = vertexID;
			break;
		}

		index = propElement[vertexID].py;

		for (j = 0; j < parametricHex.edgeValenceNumber[index]; j++)
		{
			nVert = parametricHex.edgeValence[index][j];

			if ((parametricHex.vertex[nVert][1] - parametricHex.vertex[index][1]) > EPSILON)
			{

				propElement[vertexID].py = nVert;

				break;

			}
		}

	}
	jump_out_endless_loop = 0;

	while ((parametricHex.vertexSign[propElement[vertexID].ny] != 2) && (five_direction_on_surface_flag[3] != 0))
	{

		jump_out_endless_loop++;

		//To do future, if there is interior line, need to rewrite.
		if (jump_out_endless_loop> numGrids)
		{
			propElement[vertexID].ny = vertexID;
			break;
		}


		index = propElement[vertexID].ny;

		for (j = 0; j < parametricHex.edgeValenceNumber[index]; j++)
		{
			nVert = parametricHex.edgeValence[index][j];

			if ((parametricHex.vertex[nVert][1] - parametricHex.vertex[index][1]) < -EPSILON)
			{

				propElement[vertexID].ny = nVert;

				break;

			}
		}

	}
	jump_out_endless_loop = 0;

	while ((parametricHex.vertexSign[propElement[vertexID].pz] != 2) && (five_direction_on_surface_flag[4] != 0))
	{
		jump_out_endless_loop++;

		//To do future, if there is interior line, need to rewrite.
		if (jump_out_endless_loop> numGrids)
		{
			propElement[vertexID].pz = vertexID;
			break;
		}

		index = propElement[vertexID].pz;

		for (j = 0; j < parametricHex.edgeValenceNumber[index]; j++)
		{
			nVert = parametricHex.edgeValence[index][j];

			if ((parametricHex.vertex[nVert][2] - parametricHex.vertex[index][2]) > EPSILON)
			{

				propElement[vertexID].pz = nVert;

				break;

			}
		}

	}
	jump_out_endless_loop = 0;

	while ((parametricHex.vertexSign[propElement[vertexID].nz] != 2) && (five_direction_on_surface_flag[5] != 0))
	{
		jump_out_endless_loop++;

		//To do future, if there is interior line, need to rewrite.
		if (jump_out_endless_loop> numGrids)
		{
			propElement[vertexID].nz = vertexID;
			break;
		}

		index = propElement[vertexID].nz;

		for (j = 0; j < parametricHex.edgeValenceNumber[index]; j++)
		{
			nVert = parametricHex.edgeValence[index][j];

			if ((parametricHex.vertex[nVert][2] - parametricHex.vertex[index][2]) < -EPSILON)
			{

				propElement[vertexID].nz = nVert;

				break;

			}
		}

	}
	jump_out_endless_loop = 0;

	return true;

}


bool PhysicalPolycube::FindPropagationBoundElements(int vertexID)
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

bool PhysicalPolycube::Projection(double paraPosition[3], double *realPosition, int patch_ID)
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
	i = patch_ID;
	
	

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
		

	

	//int indexPatch;
	int iConst;
	int indexPatch;
	int theta;
	if (!onArcsFlag)
	{

		for (i = 0; i < elementNumber; i++)
		{

			
			indexPatch = elementArray[i].indexPatch;

			if (indexPatch==patch_ID)
			{

			

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
				dVec[0] = (dPara[0][0] - dPara[2][0])*(dPara[1][1] - dPara[2][1]) - (dPara[1][0] - dPara[2][0])*(dPara[0][1] - dPara[2][1]);
				dVec[1] = (dPara[0][2] - dPara[2][2])*(dPara[1][1] - dPara[2][1]) - (dPara[1][2] - dPara[2][2])*(dPara[0][1] - dPara[2][1]);
				dVec[2] = (dPara[0][0] - dPara[2][0])*(dPara[1][2] - dPara[2][2]) - (dPara[1][0] - dPara[2][0])*(dPara[0][2] - dPara[2][2]);

				if (fabs(dVec[0]) > EPSILON)
				{

					dSign[0] = ((dPara[1][1] - dPara[2][1])*(paraPosition[0] - dPara[2][0]) + (dPara[2][0] - dPara[1][0])*(paraPosition[1] - dPara[2][1])) / dVec[0];
					dSign[1] = ((dPara[2][1] - dPara[0][1])*(paraPosition[0] - dPara[2][0]) + (dPara[0][0] - dPara[2][0])*(paraPosition[1] - dPara[2][1])) / dVec[0];
					dSign[2] = 1 - dSign[0] - dSign[1];

				}
				else if (fabs(dVec[1]) > EPSILON)
				{

					dSign[0] = ((dPara[1][1] - dPara[2][1])*(paraPosition[2] - dPara[2][2]) + (dPara[2][2] - dPara[1][2])*(paraPosition[1] - dPara[2][1])) / dVec[1];
					dSign[1] = ((dPara[2][1] - dPara[0][1])*(paraPosition[2] - dPara[2][2]) + (dPara[0][2] - dPara[2][2])*(paraPosition[1] - dPara[2][1])) / dVec[1];
					dSign[2] = 1 - dSign[0] - dSign[1];

				}
				else if (fabs(dVec[2]) > EPSILON)
				{

					dSign[0] = ((dPara[1][2] - dPara[2][2])*(paraPosition[0] - dPara[2][0]) + (dPara[2][0] - dPara[1][0])*(paraPosition[2] - dPara[2][2])) / dVec[2];
					dSign[1] = ((dPara[2][2] - dPara[0][2])*(paraPosition[0] - dPara[2][0]) + (dPara[0][0] - dPara[2][0])*(paraPosition[2] - dPara[2][2])) / dVec[2];
					dSign[2] = 1 - dSign[0] - dSign[1];

				}
				else
				{

					dSign[0] = 0.33;
					dSign[1] = 0.33;
					dSign[2] = 1 - dSign[0] - dSign[1];

				}

				if (dSign[0] > -1e-3 && dSign[1] > -1e-3 && dSign[2] > -1e-3)
				{

					for (j = 0; j < 3; j++)
					{
						realPosition[j] = dSign[0] * dP[0][j] + dSign[1] * dP[1][j] + dSign[2] * dP[2][j];
					}

					return true;

				}
			}

			}

		}

	}

	//return true;
	return false;

}



bool PhysicalPolycube::PropagationSurface(int vertexID)
{

	double dTemp;

	int tempIndexOne, tempIndexTwo;

	int i;
	//if (vertexID == 41833)
	//{
	//	vertexID = vertexID;
	//}

	/*for (i = 0; i < 3; i++)
	{
		realDomainHex.vertex[vertexID][i] = 0.f;
	}*/

	






	tempIndexOne = propElement[vertexID].nx;
	tempIndexTwo = propElement[vertexID].px;

	dTemp = (parametricHex.vertex[vertexID][0] - parametricHex.vertex[tempIndexOne][0]) / (parametricHex.vertex[tempIndexTwo][0] - parametricHex.vertex[tempIndexOne][0]);
	//dTemp = 0.5;

	if (dTemp > -EPSILON && dTemp < 1 + EPSILON)
	{

		realDomainHex.vertex[vertexID][0] = (1 - dTemp) * realDomainHex.vertex[tempIndexOne][0] + dTemp * realDomainHex.vertex[tempIndexTwo][0];

	}

	tempIndexOne = propElement[vertexID].ny;
	tempIndexTwo = propElement[vertexID].py;

	dTemp = (parametricHex.vertex[vertexID][1] - parametricHex.vertex[tempIndexOne][1]) / (parametricHex.vertex[tempIndexTwo][1] - parametricHex.vertex[tempIndexOne][1]);

	if (dTemp > -EPSILON && dTemp < 1 + EPSILON)
	{

		realDomainHex.vertex[vertexID][1] = (1 - dTemp) * realDomainHex.vertex[tempIndexOne][1] + dTemp * realDomainHex.vertex[tempIndexTwo][1];

	}

	tempIndexOne = propElement[vertexID].nz;
	tempIndexTwo = propElement[vertexID].pz;

	dTemp = (parametricHex.vertex[vertexID][2] - parametricHex.vertex[tempIndexOne][2]) / (parametricHex.vertex[tempIndexTwo][2] - parametricHex.vertex[tempIndexOne][2]);

	if (dTemp > -EPSILON && dTemp < 1 + EPSILON)
	{

		realDomainHex.vertex[vertexID][2] = (1 - dTemp) * realDomainHex.vertex[tempIndexOne][2] + dTemp * realDomainHex.vertex[tempIndexTwo][2];

	}


	//edit for surface interpolation.

	//double dTemp_1, dTemp_2;
	//if ((propElement[vertexID].nx == vertexID) || (propElement[vertexID].px == vertexID))
	//{
	//	tempIndexOne = propElement[vertexID].ny;
	//	tempIndexTwo = propElement[vertexID].py;

	//	dTemp_1 = realDomainHex.vertex[tempIndexTwo][0] - realDomainHex.vertex[tempIndexOne][0];

	//	tempIndexOne = propElement[vertexID].nz;
	//	tempIndexTwo = propElement[vertexID].pz;

	//	dTemp_2 = realDomainHex.vertex[tempIndexTwo][0] - realDomainHex.vertex[tempIndexOne][0];

	//	if (abs(dTemp_1)>abs(dTemp_2))
	//	{
	//		tempIndexOne = propElement[vertexID].ny;
	//		tempIndexTwo = propElement[vertexID].py;
	//		dTemp = (parametricHex.vertex[vertexID][1] - parametricHex.vertex[tempIndexOne][1]) / (parametricHex.vertex[tempIndexTwo][1] - parametricHex.vertex[tempIndexOne][1]);
	//	}
	//	else
	//	{
	//		tempIndexOne = propElement[vertexID].nz;
	//		tempIndexTwo = propElement[vertexID].pz;
	//		dTemp = (parametricHex.vertex[vertexID][2] - parametricHex.vertex[tempIndexOne][2]) / (parametricHex.vertex[tempIndexTwo][2] - parametricHex.vertex[tempIndexOne][2]);

	//	}



	//	if (dTemp > -EPSILON && dTemp < 1 + EPSILON)
	//	{

	//		realDomainHex.vertex[vertexID][0] = (1 - dTemp) * realDomainHex.vertex[tempIndexOne][0] + dTemp * realDomainHex.vertex[tempIndexTwo][0];

	//	}

	//}
	//else if ((propElement[vertexID].ny == vertexID) || (propElement[vertexID].py == vertexID))
	//{
	//	tempIndexOne = propElement[vertexID].nx;
	//	tempIndexTwo = propElement[vertexID].px;

	//	dTemp_1 = realDomainHex.vertex[tempIndexTwo][1] - realDomainHex.vertex[tempIndexOne][1];

	//	tempIndexOne = propElement[vertexID].nz;
	//	tempIndexTwo = propElement[vertexID].pz;

	//	dTemp_2 = realDomainHex.vertex[tempIndexTwo][1] - realDomainHex.vertex[tempIndexOne][1];

	//	if (abs(dTemp_1)>abs(dTemp_2))
	//	{
	//		tempIndexOne = propElement[vertexID].nx;
	//		tempIndexTwo = propElement[vertexID].px;

	//		dTemp = (parametricHex.vertex[vertexID][0] - parametricHex.vertex[tempIndexOne][0]) / (parametricHex.vertex[tempIndexTwo][0] - parametricHex.vertex[tempIndexOne][0]);

	//	}
	//	else
	//	{
	//		tempIndexOne = propElement[vertexID].nz;
	//		tempIndexTwo = propElement[vertexID].pz;
	//		dTemp = (parametricHex.vertex[vertexID][2] - parametricHex.vertex[tempIndexOne][2]) / (parametricHex.vertex[tempIndexTwo][2] - parametricHex.vertex[tempIndexOne][2]);

	//	}



	//	if (dTemp > -EPSILON && dTemp < 1 + EPSILON)
	//	{

	//		realDomainHex.vertex[vertexID][1] = (1 - dTemp) * realDomainHex.vertex[tempIndexOne][1] + dTemp * realDomainHex.vertex[tempIndexTwo][1];

	//	}

	//}
	//else if ((propElement[vertexID].nz == vertexID) || (propElement[vertexID].pz == vertexID))
	//{
	//	tempIndexOne = propElement[vertexID].nx;
	//	tempIndexTwo = propElement[vertexID].px;

	//	dTemp_1 = realDomainHex.vertex[tempIndexTwo][2] - realDomainHex.vertex[tempIndexOne][2];

	//	tempIndexOne = propElement[vertexID].ny;
	//	tempIndexTwo = propElement[vertexID].py;

	//	dTemp_2 = realDomainHex.vertex[tempIndexTwo][2] - realDomainHex.vertex[tempIndexOne][2];

	//	if (abs(dTemp_1)>abs(dTemp_2))
	//	{
	//		tempIndexOne = propElement[vertexID].nx;
	//		tempIndexTwo = propElement[vertexID].px;

	//		dTemp = (parametricHex.vertex[vertexID][0] - parametricHex.vertex[tempIndexOne][0]) / (parametricHex.vertex[tempIndexTwo][0] - parametricHex.vertex[tempIndexOne][0]);

	//	}
	//	else
	//	{
	//		tempIndexOne = propElement[vertexID].ny;
	//		tempIndexTwo = propElement[vertexID].py;

	//		dTemp = (parametricHex.vertex[vertexID][1] - parametricHex.vertex[tempIndexOne][1]) / (parametricHex.vertex[tempIndexTwo][1] - parametricHex.vertex[tempIndexOne][1]);

	//	}



	//	if (dTemp > -EPSILON && dTemp < 1 + EPSILON)
	//	{

	//		realDomainHex.vertex[vertexID][2] = (1 - dTemp) * realDomainHex.vertex[tempIndexOne][2] + dTemp * realDomainHex.vertex[tempIndexTwo][2];

	//	}

	//}




	return true;

}


bool PhysicalPolycube::Propagation(int vertexID)
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

	/*if (vertexID==20)
	{
		cout << "vertexID==20: " << dTemp << endl;
	}*/

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

	/*if (vertexID == 20)
	{
		cout << "vertexID==20: " << dTemp << endl;
	}
*/

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

	/*if (vertexID == 20)
	{
		cout << "vertexID==20: " << dTemp << endl;
	}

*/
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


bool PhysicalPolycube::Propagation_wthout_initial(int vertexID, array<int, 3>& coord_index)
{

	double dTemp;

	int tempIndexOne, tempIndexTwo;

	int i;
	//if (vertexID == 41833)
	//{
	//	vertexID = vertexID;
	//}



	tempIndexOne = propElement[vertexID].nx;
	tempIndexTwo = propElement[vertexID].px;

	dTemp = (parametricHex.vertex[vertexID][0] - parametricHex.vertex[tempIndexOne][0]) / (parametricHex.vertex[tempIndexTwo][0] - parametricHex.vertex[tempIndexOne][0]);
	//dTemp = 0.5;

	/*if (vertexID==20)
	{
	cout << "vertexID==20: " << dTemp << endl;
	}*/

	if (dTemp > -EPSILON && dTemp < 1 + EPSILON)
	{

		realDomainHex.vertex[vertexID][coord_index[0]] = (1 - dTemp) * realDomainHex.vertex[tempIndexOne][coord_index[0]] + dTemp * realDomainHex.vertex[tempIndexTwo][coord_index[0]];

		//for (i = 0; i < 3; i++)
		//{

		//	realDomainHex.vertex[vertexID][i] += 1.0/3.0 * ((1-dTemp) * realDomainHex.vertex[tempIndexOne][i] + dTemp * realDomainHex.vertex[tempIndexTwo][i]);

		//}

	}

	tempIndexOne = propElement[vertexID].ny;
	tempIndexTwo = propElement[vertexID].py;

	dTemp = (parametricHex.vertex[vertexID][1] - parametricHex.vertex[tempIndexOne][1]) / (parametricHex.vertex[tempIndexTwo][1] - parametricHex.vertex[tempIndexOne][1]);

	/*if (vertexID == 20)
	{
	cout << "vertexID==20: " << dTemp << endl;
	}
	*/

	if (dTemp > -EPSILON && dTemp < 1 + EPSILON)
	{

		realDomainHex.vertex[vertexID][coord_index[1]] = (1 - dTemp) * realDomainHex.vertex[tempIndexOne][coord_index[1]] + dTemp * realDomainHex.vertex[tempIndexTwo][coord_index[1]];

		//for (i = 0; i < 3; i++)
		//{

		//	realDomainHex.vertex[vertexID][i] += 1.0/3.0 * ((1-dTemp) * realDomainHex.vertex[tempIndexOne][i] + dTemp * realDomainHex.vertex[tempIndexTwo][i]);

		//}

	}

	tempIndexOne = propElement[vertexID].nz;
	tempIndexTwo = propElement[vertexID].pz;

	dTemp = (parametricHex.vertex[vertexID][2] - parametricHex.vertex[tempIndexOne][2]) / (parametricHex.vertex[tempIndexTwo][2] - parametricHex.vertex[tempIndexOne][2]);

	/*if (vertexID == 20)
	{
	cout << "vertexID==20: " << dTemp << endl;
	}

	*/
	if (dTemp > -EPSILON && dTemp < 1 + EPSILON)
	{

		realDomainHex.vertex[vertexID][coord_index[2]] = (1 - dTemp) * realDomainHex.vertex[tempIndexOne][coord_index[2]] + dTemp * realDomainHex.vertex[tempIndexTwo][coord_index[2]];

		//for (i = 0; i < 3; i++)
		//{

		//	realDomainHex.vertex[vertexID][i] += 1.0/3.0 * ((1-dTemp) * realDomainHex.vertex[tempIndexOne][i] + dTemp * realDomainHex.vertex[tempIndexTwo][i]);

		//}

	}

	return true;

}

bool PhysicalPolycube::HexMeshProjectionInterior()
{

	int i;

	

	for (i = 0; i < parametricHex.vertexNumber; i++)
	{

		if (parametricHex.vertexSign[i] != 0)
		{
			continue;
		}
		

		//outputCentroids<<i<<": "<<propElement[i].px<<" "<<propElement[i].nx<<" "<<propElement[i].py<<" "<<propElement[i].ny<<" "<<propElement[i].pz<<" "<<propElement[i].nz<<endl;




		Propagation(i);

	}

	return true;

}

bool PhysicalPolycube::HexMeshOctreeSubdivision(const char *outputName)
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
					int patch_ID = 0;   //Just make it run future rewrite e
					if (NewNode[j].BoundaryFlag)
					{
						Projection(NewNode[j].ParaPos, NewNode[j].Coords, patch_ID);
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


double		PhysicalPolycube::TriArea(double v0[3], double v1[3], double v2[3])
{
	Vector3d n;
	Vector3d a(v1[0]-v0[0], v1[1]-v0[1], v1[2]-v0[2]);
	Vector3d b(v2[0]-v0[0], v2[1]-v0[1], v2[2]-v0[2]);
	n = a.cross(b);
	return n.norm();
}

double		PhysicalPolycube::QuadArea(double v0[3], double v1[3], double v2[3], double v3[3], double *MassCenter)
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
double		PhysicalPolycube::GetHexVolume(double p[8][3], double *MassCenter)
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


int PhysicalPolycube::GetDepth(int res)
{
	int i = 0;
	while (1)
	{
		if (res <= (1 << i))//take care of this!
		{
			break;
		}
		i++;
	}
	if (res != (1 << i))
	{
		cout << "Unsupported resolution:" << res;
	}
	return i;
}

int PhysicalPolycube::GetOctreeNum(int depth)
{
	int num = 0;
	for (int i = 0; i <= depth; i++)
	{
		num += (1 << (i * 3));
	}
	return num;
}

int PhysicalPolycube::GetLevel(int octree_id)
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

void PhysicalPolycube::OctreeidxToXYZ(int octree_id, int &x, int &y, int &z, int level)
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
int PhysicalPolycube::XYZToOctreeidx(int x, int y, int z, int level, int direction)
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

bool PhysicalPolycube::ReadKFileInitial(const char *inputName)
{
	int i, j, tempElementSign;
	string oneLine;

	string inputFileName = inputName;
	fstream input(inputFileName, fstream::in);
	if (!input)
	{
		cout << "open input file error!!!" << endl;
		return false;
	}

	getline(input, oneLine);
	getline(input, oneLine);
	getline(input, oneLine);
	getline(input, oneLine);

	for (i = 0; i < elementNumber; i++)
	{
		input >> j >> tempElementSign >> j >> j >> j >> j >> j >> j >> j >> j;
		elementArray[i].indexPart = tempElementSign - 1;
	}

	input.close();

	return true;

}




bool PhysicalPolycube::ReadKFileBeforePostProcessing(const char *inputName)
{
	int i, j, tempElementSign;
	string oneLine;

	string inputFileName = inputName;
	fstream input(inputFileName, fstream::in);
	if (!input)
	{
		cout << "open input file error!!!" << endl;
		return false;
	}

	getline(input, oneLine);
	getline(input, oneLine);
	getline(input, oneLine);
	getline(input, oneLine);

	for (i = 0; i < elementNumber; i++)
	{
		input >> j >> tempElementSign >> j >> j >> j >> j >> j >> j >> j >> j;
		elementArray[i].indexPatch = tempElementSign - 1;
	}

	input.close();

	return true;

}

bool PhysicalPolycube::WriteKFileBeforePostProcessing(const char *outputName)
{

	int i, j;

	string outputFileName = outputName;

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
	for (i = 0; i < elementNumber; i++)
	{
		output <<i+1<<","<<elementArray[i].indexPatch+1<<","<<element[i][0]+1<<","<<element[i][1]+1<<","<<element[i][2]+1<<","<<element[i][2]+1<<"\n";
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

bool PhysicalPolycube::ReadKFileIndexPatch(const char *inputName)
{
	int i, j, tempElementSign;
	string oneLine;

	string inputFileName = inputName;
	fstream input(inputFileName, fstream::in);
	if (!input)
	{
		cout << "open input file error!!!" << endl;
		return false;
	}

	getline(input, oneLine);
	getline(input, oneLine);
	getline(input, oneLine);
	getline(input, oneLine);

	for (i = 0; i < elementNumber; i++)
	{
		input >> j >> tempElementSign >> j >> j >> j >> j >> j >> j >> j >> j;
		elementArray[i].indexPatch = tempElementSign - 1;
	}

	input.close();

	return true;

}

bool PhysicalPolycube::WriteKFileIndexPatch(const char *outputName)
{

	int i, j;

	string outputFileName = outputName;

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
	for (i = 0; i < elementNumber; i++)
	{
		output <<i+1<<","<<elementArray[i].indexPatch+1<<","<<element[i][0]+1<<","<<element[i][1]+1<<","<<element[i][2]+1<<","<<element[i][2]+1<<"\n";
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

bool PhysicalPolycube::WriteKFileMapping(const char *outputName)
{

	int i, j;

	string outputFileName = outputName;

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

bool PhysicalPolycube::ReadKFileHexMeshPara(const char *inputName)
{

	
	return true;

}

//by yu

bool PhysicalPolycube::CreateInitialPolycube()
{
	int i, j;

	string outputFileName = "InitialPolycube.k";
	

	outputFileName =  inputName+"_"+outputFileName;
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
	for (i = 0; i < NUM_POLYCUBE_PATCH; i++)
	{	
		for (int k = 0; k < polycubePatch[i].numCorner-2; k++)
		{
			output << elementKnumber + 1 << "," << polycubePatch[i].index + 1<<"," << polycubePatch[i].cornerPoint[0]+1;
			for (j = 1; j < 3; j++)
			{
				output << "," << polycubePatch[i].cornerPoint[k+j] + 1;
			}
			output << "," << polycubePatch[i].cornerPoint[k + 2] + 1;
			elementKnumber++;
			output << "\n";
		}
	}

	output << "*NODE" << "\n";

	for (int i = 0; i < numberCornerPoints; i++)
		{
			int index = cornerPoints[i];
			output<<index+1;

			for (int j = 0; j < 3; j++)
			{
				output<<"," << vertex[index][j]<<" ";
			}

			output<<endl;
		}

	output << "*END" << "\n";
	
	output.close();


	if (OUTPUT_CORNER_POINTS == 1)
	{
		///////////////////////////////////////////////////////////////OUT PUT CORNER POINTS
		outputCornerPoints.open(inputName + "_corner.txt");

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

	cout << "test" << endl;

	//export edge

	outputFileName = "edge.txt";


	outputFileName = inputName + "_" + outputFileName;
	fstream output2(outputFileName, fstream::out);
	if (!output2)
	{
		cout << "open input file error!!!" << endl;
		return false;
	}
	//no need+1
	
	for (i = 0; i < NUM_POLYCUBE_PATCH; i++)
	{
		
			for (j = 0; j < polycubePatch[i].numCorner-1; j++)
			{
				if (polycubePatch[i].cornerPoint[j]<polycubePatch[i].cornerPoint[j + 1])
				{
					output2 << polycubePatch[i].cornerPoint[j]  << "," << polycubePatch[i].cornerPoint[j + 1]  << "\n";
				}
				else {
					output2 << polycubePatch[i].cornerPoint[j+1]  << "," << polycubePatch[i].cornerPoint[j]  << "\n";
				}
				
			}

			if (polycubePatch[i].cornerPoint[0]<polycubePatch[i].cornerPoint[polycubePatch[i].numCorner - 1])
			{
				output2 << polycubePatch[i].cornerPoint[0] << "," << polycubePatch[i].cornerPoint[polycubePatch[i].numCorner - 1]  << "\n";
			}
			else {
				output2 << polycubePatch[i].cornerPoint[polycubePatch[i].numCorner - 1]  << "," << polycubePatch[i].cornerPoint[0] << "\n";
			}
		
	}

	output2.close();


	//export edge

	outputFileName = "face.txt";


	outputFileName = inputName + "_" + outputFileName;
	fstream output3(outputFileName, fstream::out);
	if (!output3)
	{
		cout << "open input file error!!!" << endl;
		return false;
	}
	//no need+1
	elementKnumber = 0;
	for (i = 0; i < NUM_POLYCUBE_PATCH; i++)
	{
		output3 << elementKnumber + 1 ;
		for (int k = 0; k < polycubePatch[i].numCorner; k++)
		{
			output3 << "," << polycubePatch[i].cornerPoint[k] + 1;
			
		}
		output3 << "\n";
	    elementKnumber++;
	}
	

	output3.close();


	return true;
	

}
