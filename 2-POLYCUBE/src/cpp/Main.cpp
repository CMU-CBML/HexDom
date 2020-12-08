#include <stdio.h>
#include <time.h>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include "cxxopts.hpp"
#include "Polycube.h"
#include "atlstr.h" 

using namespace std;
//int argc, char* argv[]
int main(int argc, char* argv[])
{

	string fn_in;
	string fn_out;
	bool flag_value_output_corner = false;

	//replace following try when debug.
	/*fn_in = "heli_dominant_only_prism";
	fn_in = "con_surface_tri_patch";
	fn_in = "heli_dominant_V2";
	fn_out = "heli_dominant_only_prism_polycube_structure";
	fn_out = "con_surface_tri_patch_polycube_structure";
	fn_out = "con_surface_tri_patch_polycube_structure";
	fn_out = "heli_dominant_V2_polycube_structure";
	flag_value_output_corner = 1;*/


	try
	{
		cxxopts::Options options(argv[0], "CMU Polycube Construction");
		options
			.positional_help("[optional args]")
			.show_positional_help();

		options.add_options("General")
			("h,help", "Print help")
			("i,input", "The input LS-DYNA Keyword file with segmentation information", cxxopts::value<std::string>(fn_in))
			("o,output", "The output LS-DYNA Keyword file with boundary surface information of polycube", cxxopts::value<std::string>(fn_out))
			("c,corner_edge_face_output_flag", "0-No output file of corner point, edge and face of polycube structure, 1-Output", cxxopts::value<bool>(flag_value_output_corner))
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
		if (!result.count("output"))
		{
			cout << "Need to assign the output LS-DYNA Keyword file name" << endl;
			exit(0);
		}
		else
		{
			char *buffer = strdup(fn_out.c_str());
			char *fn_wo_extension_char = buffer;
			PathRemoveExtensionA(fn_wo_extension_char);
			stringstream ss_in;
			ss_in << fn_wo_extension_char;
			ss_in >> fn_out;
		}


	}
	catch (const cxxopts::OptionException& e)
	{
		std::cout << "error parsing options: " << e.what() << std::endl;
		exit(1);
	}


	string inputmeshName = fn_in;
	string inputFullName;
	inputFullName = inputmeshName + ".k";

	PolycubeConstruction Pcube(inputFullName.c_str());


	bool Debug_yu = true;
	
	if (Debug_yu)
	{
		cout << inputFullName << endl;
	}
	//Pcube.Initialization();

	Pcube.InitializePolycube();

	string outputmeshName = fn_out;
	string outputFullName = outputmeshName + ".k";

	string tempString;
	if (Debug_yu)
	{
		//tempString = outputmeshName + "_write.k";
		Pcube.WriteKFileIndexPatch(tempString.c_str());
		tempString = outputmeshName + "_segmentation_colors.vtk";
		//Pcube.OutputPatchesVTK(tempString.c_str());
	}
	
	cout << "here"<< Pcube.vertexNumber;
	Pcube.CreateInitialPolycube(outputmeshName.c_str());
	cout << "here";
	if (flag_value_output_corner)
	{
		Pcube.OutputExtraInformation(outputmeshName.c_str());
	}
	

	return 0;

}