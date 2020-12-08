#include <stdio.h>
#include <time.h>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>

#include "physical_polycube.h"
#include "Main.h"
#include "cxxopts.hpp"
#include "atlstr.h" 
using namespace std;

int main(int argc, char* argv[])
{
	bool debug_yu = true; //used for paper picture .

	string fn_in;
	string fn_out;

	fn_in = "heli_dominant_V2";
	fn_out = "heli_hex.vtk";


	string fn_manual_file;
	bool flag_manual = false;
	string fn_out_k;
	string fn_polycube="heli_dominant_only_prism_te_polycube_structure.k";

	//fn_polycube = "PART_FOR_HEX_MESHING_Polycube_V5.k";
	string fn_out_vtk;
	int octree_level_default = 2;
	double value_lambda = 0;
	bool flag_value_lambda = false;
	int hex_number = 0;
	try
	{
		cxxopts::Options options(argv[0], "CMU Create Segmentation");
		options
			.positional_help("[optional args]")
			.show_positional_help();

		options.add_options("General")
			("h,help", "Print help")
			("i,input", "The input LS-DYNA Keyword file with segmentation information", cxxopts::value<std::string>(fn_in))
			("p,polycube_file", "The output LS-DYNA Keyword file with boundary surface information of polycube", cxxopts::value<std::string>(fn_polycube))
			("o,output", "The output LS-DYNA Keyword file with unstructured hex mesh", cxxopts::value<std::string>(fn_out))
			("s,octree_subdivision", "Octree level, default is 2", cxxopts::value<int>(octree_level_default))
			("n,hex_number", "hex_number", cxxopts::value<int>(hex_number))

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
			char *buffer = strdup(fn_in.c_str());
			char *fn_wo_extension_char = buffer;
			PathRemoveExtensionA(fn_wo_extension_char);
			stringstream ss_in;
			ss_in << fn_wo_extension_char;
			ss_in >> fn_in;
		}
		else
		{
			cout << "Need to assign the input LS-DYNA Keyword file name" << endl;
			exit(0);
		}
		if (!result.count("polycube_file"))
		{
			cout << "Need to assign the polycube LS-DYNA Keyword file name" << endl;
			exit(0);
		}
		if (!result.count("hex_number"))
		{
			cout << "Need to assign the hex number" << endl;
			exit(0);
		}
		if (!result.count("output"))
		{
			cout << "Need to assign the output LS-DYNA Keyword file name" << endl;
			exit(0);
		}

		

	}
	catch (const cxxopts::OptionException& e)
	{
		std::cout << "error parsing options: " << e.what() << std::endl;
		exit(1);
	}





	string inputmeshName = fn_in;

	cout << inputmeshName << endl;

	//getchar();

	string outputmeshName;

	PhysicalPolycube a_physical_polycube_class;	

	string inputFullName = inputmeshName + ".k";//input file, triangle mesh

	a_physical_polycube_class.Initialization(inputFullName, fn_polycube);
	
	a_physical_polycube_class.OCTREE_MAX_LEVEL = octree_level_default + 1;

	a_physical_polycube_class.OCTREE_MIN_LEVEL = octree_level_default - 1;
	
	if (2==(a_physical_polycube_class.OCTREE_MAX_LEVEL- a_physical_polycube_class.OCTREE_MIN_LEVEL))
	{
		cout << "octree level is: " << a_physical_polycube_class.OCTREE_MAX_LEVEL - 1 << endl;
	}

	a_physical_polycube_class.InitializePolycube(inputFullName);
	
	if (debug_yu)
	{
		cout << "run 1 here" << endl;
	}
	//cout << "run 1 here" << endl;
	a_physical_polycube_class.PreProcessing();
	

	if (debug_yu)
	{
		cout << "run 2 here" << endl;
	}
	//

	a_physical_polycube_class.HexMeshParametricDomain();

	if (debug_yu)
	{
		cout << "run 3 here" << endl;
	}

	//cout << "run 3 here" << endl;

	//output unit cube
	/*string file_name = "octree_outside_" + std::to_string(a_physical_polycube_class.OCTREE_MAX_LEVEL-1) + "_test_hex.vtk";
	a_physical_polycube_class.parametricHex.DeleteDuplicatedPoint();
	a_physical_polycube_class.parametricHex.Write(file_name.c_str());*/

	a_physical_polycube_class.MappingRealPara(fn_out, hex_number);


	//cout << "run 4 here" << endl;
	

	//outputmeshName = inputmeshName + "_output_initial.vtk";
	//a_physical_polycube_class.OutputPatchesVTK(outputmeshName.c_str());

	//a_physical_polycube_class.EdgeWeightedCVT();
	
	//a_physical_polycube_class.PostProcessing();
	

	
	//outputmeshName = inputmeshName + "_output_ModifiedBifur.vtk";
	//a_physical_polycube_class.OutputPatchesVTKBifurcation(outputmeshName.c_str()); // Be careful, the results are different with OutputPatchesVTK()!
	//outputmeshName = inputmeshName + "_output_paraMapping.vtk";
	//a_physical_polycube_class.OutputPatchesVTKPara(outputmeshName.c_str());

	//outputmeshName = inputmeshName + "_paraHexNew_hex.raw";
	//a_physical_polycube_class.HexMeshParametricDomain(outputmeshName.c_str());

	//outputmeshName = inputmeshName + "_realHex_hex.vtk";
	//a_physical_polycube_class.HexMeshRealDomain(outputmeshName.c_str());

	//outputmeshName = inputmeshName + "_finalHex_hex.raw";
	//a_physical_polycube_class.HexMeshOctreeSubdivision(outputmeshName.c_str());

	return 0;

}