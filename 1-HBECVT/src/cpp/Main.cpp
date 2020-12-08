#include <stdio.h>
#include <time.h>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include "cxxopts.hpp"
#include "CVTBasedPolycube.h"
#include "atlstr.h" 
using namespace std;


int main(int argc, char* argv[])
{

	bool debug_yu = true; //used for paper picture .

	string fn_in;
	string fn_out;
	string fn_manual_file;
	bool flag_manual = false;
	string fn_out_k;
	string fn_out_vtk;
	double value_lambda = 0;
	bool flag_value_lambda = false;

	try
	{
		cxxopts::Options options(argv[0], "CMU Create Segmentation");
		options
			.positional_help("[optional args]")
			.show_positional_help();

		options.add_options("General")
			("h,help", "Print help")
			("i,input", "The input LS-DYNA Keyword file name", cxxopts::value<std::string>(fn_in))
			//("m,manual", "Need the manual information, default is 0", cxxopts::value<int>(manual_seg))
			("m,manual-file", "The file including manual information", cxxopts::value<std::string>(fn_manual_file))
			("o,output", "The output LS-DYNA Keyword file name", cxxopts::value<std::string>(fn_out))
			/*("o,output","Output path", cxxopts::value<std::string>(fn_out))*/
			("l,lambda-weight", "Weight", cxxopts::value<double>(value_lambda))
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

		if (result.count("manual-file"))
		{
			flag_manual = true;
			

		}
		if (result.count("lambda-weight"))
		{
			flag_value_lambda=true;
		}

	}
	catch (const cxxopts::OptionException& e)
	{
		std::cout << "error parsing options: " << e.what() << std::endl;
		exit(1);
	}


	string inputmeshName = fn_in;




	string inputFullName;

	string outputmeshName;


	CVTBasedPolycube polycube;

	inputFullName = inputmeshName + ".k";//input file, triangle mesh

	polycube.Initialization(inputFullName.c_str());//generate initial segmentation results

	polycube.ClassicalCVT();
	if (debug_yu)
	{
		outputmeshName = inputmeshName + "_output_CVT.vtk";
		polycube.OutputPatchesVTK(outputmeshName.c_str());
	}
	
	//cout << value_lambda<<"xxxxxx";
	polycube.WEIGHT_LENGTH_EWCVT = value_lambda;
	
	

	polycube.EdgeWeightedCVT();
	if (debug_yu)
	{
		outputmeshName = inputmeshName + "_output_EWCVT.vtk";
		polycube.OutputPatchesVTK(outputmeshName.c_str());
	}

	if (flag_manual)
	{
		cout << fn_manual_file << endl;
		polycube.ModifySpecialElement(fn_manual_file.c_str());
	}
	

	if (debug_yu)
	{
		outputmeshName = inputmeshName + "_output_EWCVT_edit.vtk";
		polycube.OutputPatchesVTK(outputmeshName.c_str());
	}

	outputmeshName = fn_out;//output initial segmentation results
	polycube.WriteKFileBeforePostProcessing(outputmeshName.c_str());
	/*}*/
		

	return 0;

}