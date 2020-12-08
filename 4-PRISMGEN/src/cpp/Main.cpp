#include <stdio.h>
#include <time.h>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include "PrismGen.h"
#include "cxxopts.hpp"
#include "atlstr.h" 





using namespace std;

int main(int argc, char* argv[])
{

	string inputmeshName;
	//inputmeshName = "heli_triangle_1.k";
	
	//inputmeshName = "heli_triangle_2.k";
	//inputmeshName = "heli_triangle_3.k";

	/*inputmeshName = "heli_bottom_1.k";
	inputmeshName = "heli_bottom_2.k";
	inputmeshName = "heli_bottom_3.k";
	inputmeshName = "heli_bottom_4.k";
	inputmeshName = "heli_triangle_4.k";
	inputmeshName = "heli_triangle_5.k";*/
	//inputmeshName = "heli_tet_initial_read.k";
	//inputmeshName = "heli_triangle_1.k";
	//inputmeshName = "pfm_tri_1.k";
	//inputmeshName = "pfm_tri_2.k";
	//inputmeshName = "pfm_tri_3.k";
	//inputmeshName = "heli_triangle_1.k";
	//string inputmeshName = "Honda";

	string fn_in;
	string fn_out;
	string corner_points;
	int octree_level_default = 2;
	
	bool flip_triangle_for_inclined_plane;
	try
	{
		cxxopts::Options options(argv[0], "CMU Create Segmentation");
		options
			.positional_help("[optional args]")
			.show_positional_help();

		options.add_options("General")
			("h,help", "Print help")
			("i,input", "The input LS-DYNA Keyword file with segmentation information", cxxopts::value<std::string>(fn_in))
			("o,output", "The output LS-DYNA Keyword file with boundary surface information of polycube", cxxopts::value<std::string>(fn_out))
			("f,flip_triangle_for_inclined_plane", "0-No change, 1-change", cxxopts::value<bool>(flip_triangle_for_inclined_plane))
			("s,octree_subdivision", "Octree level, default is 2", cxxopts::value<int>(octree_level_default))
			//("c,corner_points", "0-No output file of corner point, edge and face of polycube structure, 1-Output", cxxopts::value<string>(corner_points))
#ifdef CXXOPTS_USE_UNICODE
			("unicode", u8"A help option with non-ascii: ��. Here the size of the"
				" string should be correct")
#endif
			;

		auto result = options.parse(argc, argv);
		if (result.count("help"))
		{
			cout << options.help({ "General" }) << endl;
			exit(0);
		}
		
		if (result.count("input"))
		{
		}
		else
		{
			cout << "Need to assign the input file name" << endl;
			exit(0);
		}
		
		if (!result.count("output"))
		{
			cout << "Need to assign the output file name" << endl;
			exit(0);
		}

		/*if (!result.count("corner_points"))
		{
			cout << "Need to assign the corner_points file name" << endl;
			exit(0);
		}
*/
	}
	catch (const cxxopts::OptionException& e)
	{
		std::cout << "error parsing options: " << e.what() << std::endl;
		exit(1);
	}


	inputmeshName = fn_in;

	string inputFullName;

	string outputmeshName;

	//inputFullName = inputmeshName + ".k";

	PrismGenClass PGC(inputmeshName);

	PGC.Initialization();
	//outputmeshName = inputmeshName + "_output_initial.vtk";
	//PGC.OutputPatchesVTK(outputmeshName.c_str());

	//PGC.EdgeWeightedCVT();
	PGC.Flip_yu = flip_triangle_for_inclined_plane;
	PGC.PostProcessing();


	/*outputmeshName = inputmeshName + "_patch_output_CleanUp.vtk";
	PGC.OutputPatchesVTK(outputmeshName.c_str());*/
	


	outputmeshName = inputmeshName + "_tri.raw";
	PGC.mapping_file = outputmeshName;
	PGC.InitializePolycube();
	
	PGC.polycubePara->Write(outputmeshName.c_str());


	outputmeshName = inputmeshName + "_output_ModifiedBifur.vtk";
	//PGC.OutputPatchesVTKBifurcation(outputmeshName.c_str()); // Be careful, the results are different with OutputPatchesVTK()!
	//getchar();

	outputmeshName = inputmeshName + "_output_paraMapping.vtk";
	cout << "here" << endl;
	PGC.OutputPatchesVTKPara(outputmeshName.c_str());
	int octree_level_yu = octree_level_default;
	outputmeshName = inputmeshName + "_paraHexNew_hex.raw";
	
	PGC.HexMeshParametricDomain(outputmeshName.c_str(), octree_level_yu);

	outputmeshName = fn_out;
	PGC.HexMeshRealDomain(outputmeshName.c_str());

	outputmeshName = inputmeshName + "_finalHex_hex.raw";
	//PGC.HexMeshOctreeSubdivision(outputmeshName.c_str());

	return 0;

}