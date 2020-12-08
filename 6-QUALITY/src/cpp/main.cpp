#include <iostream>
#include "hex_quality.h"
#include "cxxopts.hpp"
using namespace std;

void Commandline(int argc, char** argv);

int main(int argc, char** argv)
{
	Commandline(argc, argv);
}

void Commandline(int argc, char** argv)
{
	HexQuality app;

	try
	{
		cxxopts::Options options(argv[0], "CMU Solid Software");
		options
			.positional_help("[optional args]")
			.show_positional_help();

		int flag_quality_option = 0;
		int flag_sharp = 0;
		double tol_sharp;
		int quality_par1;
		double quality_par2;

		string fn_in;
		string fn_out;

		options.add_options("General Settings")
			("h,help", "Print help")
			("s,sharp", "0-No sharp feature, 1-Automatic sharp feature, 2-Manual sharp feature", cxxopts::value<int>(flag_sharp))
			("t,stol", "Tolerance for automatically detecting sharp feature", cxxopts::value<double>(tol_sharp))
			("I,input", "Input file", cxxopts::value<std::string>(fn_in))
#ifdef CXXOPTS_USE_UNICODE
			("unicode", u8"A help option with non-ascii: ��. Here the size of the"
				" string should be correct")
#endif
			;
		options.add_options("Mesh Quality Improvement")
			("m,method",
				"Improvement methods: 0-Laplacian Smoothing interior points (Give iternation number -n)\; 1-Pillowing (Give pillow layer number -n)\; 2-Smoothing (Give iteration number -n and smooth step -p)\; 3-Optimization (Give iteration number -n and optimization step -p)", cxxopts::value<int>(flag_quality_option))
				("n,number", "Pillowing layer number, Smoothing and Optimization number of steps", cxxopts::value<int>(quality_par1))
			("p,parameter", "Smoothing / Optimization step size", cxxopts::value<double>(quality_par2))
#ifdef CXXOPTS_USE_UNICODE
			("unicode", u8"A help option with non-ascii: ��. Here the size of the"
				" string should be correct")
#endif
			;
			

		auto result = options.parse(argc, argv);
		if (result.count("help"))
		{
			cout << options.help({ "General Settings", "Mesh Quality Improvement" }) << endl;
			exit(0);
		}


		app.run_MeshQualityImprove(flag_quality_option, flag_sharp, tol_sharp, quality_par1, quality_par2, fn_in);//mesh quality improvement

			
	}
	catch (const cxxopts::OptionException& e)
	{
		std::cout << "error parsing options: " << e.what() << std::endl;
		exit(1);
	}
}

