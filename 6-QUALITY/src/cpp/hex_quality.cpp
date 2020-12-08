#include "hex_quality.h"
#include <fstream>
#include <iostream>
#include <algorithm>
#include <cmath>
#include <map>
#include <sstream>
#include <iomanip>

typedef unsigned int uint;

HexQuality::HexQuality()
{
	//cp.clear();
	//tmesh.clear();
}

void HexQuality::BuildInitialEdges()
{
	int edloc[12][2] = { {0,1},{1,2},{2,3},{3,0},{0,4},{1,5},{2,6},{3,7},{4,5},{5,6},{6,7},{7,4} };
	int fcloc[6][4] = { {0,3,2,1},{0,1,5,4},{1,2,6,5},{2,3,7,6},{0,4,7,3},{4,5,6,7} };
	int fced[6][4] = { {3,2,1,0},{0,5,8,4},{1,6,9,5},{2,7,10,6},{4,11,7,3},{8,9,10,11} };

	uint i, j, k;
	tmedge.clear();
	tmface.clear();
	//point-hex relation
	for (i = 0; i < tmesh.size(); i++)
	{
		for (j = 0; j < 8; j++)
		{
			cp[tmesh[i].cnct[j]].hex.push_back(i);
		}
	}
	//construct edges
	for (i = 0; i < tmesh.size(); i++)
	{
		vector<int> nb;
		for (j = 0; j < 8; j++)
		{
			for (k = 0; k < cp[tmesh[i].cnct[j]].hex.size(); k++)
			{
				int eid(cp[tmesh[i].cnct[j]].hex[k]);
				if (eid < i)
				{
					vector<int>::iterator it = find(nb.begin(), nb.end(), eid);
					if (it == nb.end())
					{
						nb.push_back(eid);
					}
				}
			}
		}
		for (j = 0; j < 12; j++)//edge
		{
			Edge3D edtmp;
			edtmp.pt[0] = tmesh[i].cnct[edloc[j][0]];
			edtmp.pt[1] = tmesh[i].cnct[edloc[j][1]];
			int flag(-1);
			for (k = 0; k < nb.size(); k++)
			{
				for (int k0 = 0; k0 < 12; k0++)
				{
					if (edtmp == tmedge[tmesh[nb[k]].edge[k0]])
					{
						flag = tmesh[nb[k]].edge[k0]; break;
					}
				}
				if (flag != -1) break;
			}
			if (flag != -1)
			{
				tmesh[i].edge[j] = flag;
			}
			else
			{
				tmedge.push_back(edtmp);
				tmesh[i].edge[j] = tmedge.size() - 1;
			}
		}
		for (j = 0; j < 6; j++)//face
		{
			Face3D fctmp;
			for (k = 0; k < 4; k++)
			{
				fctmp.cnct[k] = tmesh[i].cnct[fcloc[j][k]];
				fctmp.edge[k] = tmesh[i].edge[fced[j][k]];
			}
			int flag(-1);
			for (k = 0; k < nb.size(); k++)
			{
				for (int k0 = 0; k0 < 6; k0++)
				{
					if (fctmp == tmface[tmesh[nb[k]].face[k0]])
					{
						flag = tmesh[nb[k]].face[k0]; break;
					}
				}
				if (flag != -1) break;
			}
			if (flag != -1)
			{
				tmesh[i].face[j] = flag;
			}
			else
			{
				tmface.push_back(fctmp);
				tmesh[i].face[j] = tmface.size() - 1;
			}
		}
	}

	for (i = 0; i < tmesh.size(); i++)
	{
		for (j = 0; j < 12; j++)
		{
			tmedge[tmesh[i].edge[j]].hex.push_back(i);
		}
		for (j = 0; j < 6; j++)
		{
			tmface[tmesh[i].face[j]].hex.push_back(i);
		}
	}
	for (i = 0; i < tmface.size(); i++)
	{
		for (j = 0; j < 4; j++)
		{
			cp[tmface[i].cnct[j]].face.push_back(i);
			tmedge[tmface[i].edge[j]].face.push_back(i);
		}
	}
	for (i = 0; i < tmedge.size(); i++)
	{
		for (j = 0; j < 2; j++)
		{
			cp[tmedge[i].pt[j]].edge.push_back(i);
		}
	}
}

void HexQuality::InitialConnect()
{
	uint i, j;
	tmedge.clear();
	tmface.clear();
	BuildInitialEdges();

	////construct edges
	//for(i=0; i<tmesh.size(); i++)
	//{
	//	for(j=0; j<4; j++)
	//	{
	//		Edge3D edtmp;
	//		edtmp.pt[0]=tmesh[i].cnct[j];
	//		edtmp.pt[1]=tmesh[i].cnct[(j+1)%4];
	//		vector<Edge3D>::iterator it=find(tmedge.begin(),tmedge.end(),edtmp);
	//		int edid(it-tmedge.begin());
	//		if(it==tmedge.end())
	//		{
	//			tmedge.push_back(edtmp);
	//		}
	//		tmesh[i].edge[j]=edid;
	//		//tmedge[edid].hex.push_back(i);
	//	}
	//	for(j=0; j<4; j++)
	//	{
	//		Edge3D edtmp;
	//		edtmp.pt[0]=tmesh[i].cnct[j];
	//		edtmp.pt[1]=tmesh[i].cnct[j+4];
	//		vector<Edge3D>::iterator it=find(tmedge.begin(),tmedge.end(),edtmp);
	//		int edid(it-tmedge.begin());
	//		if(it==tmedge.end())
	//		{
	//			tmedge.push_back(edtmp);
	//		}
	//		tmesh[i].edge[j+4]=edid;
	//		//tmedge[edid].hex.push_back(i);
	//	}
	//	for(j=0; j<4; j++)
	//	{
	//		Edge3D edtmp;
	//		edtmp.pt[0]=tmesh[i].cnct[j+4];
	//		edtmp.pt[1]=tmesh[i].cnct[(j+1)%4+4];
	//		vector<Edge3D>::iterator it=find(tmedge.begin(),tmedge.end(),edtmp);
	//		int edid(it-tmedge.begin());
	//		if(it==tmedge.end())
	//		{
	//			tmedge.push_back(edtmp);
	//		}
	//		tmesh[i].edge[j+8]=edid;
	//		//tmedge[edid].hex.push_back(i);
	//	}
	//}
	////construct faces
	//for(i=0; i<tmesh.size(); i++)
	//{
	//	//one bottom face
	//	Face3D fc1;
	//	for(j=0; j<4; j++)
	//	{
	//		fc1.cnct[j]=tmesh[i].cnct[j];
	//		fc1.edge[j]=tmesh[i].edge[j];
	//	}
	//	vector<Face3D>::iterator it1=find(tmface.begin(),tmface.end(),fc1);
	//	int fc1id(it1-tmface.begin());
	//	if(it1==tmface.end())
	//	{
	//		tmface.push_back(fc1);
	//	}
	//	tmesh[i].face[0]=fc1id;
	//	//tmface[fc1id].hex.push_back(i);
	//	//4 side faces
	//	for(j=0; j<4; j++)
	//	{
	//		Face3D fc;
	//		for(int k=0; k<4; k++)
	//		{
	//			fc.cnct[k]=tmesh[i].cnct[fc_cnct[j][k]];
	//			fc.edge[k]=tmesh[i].edge[ed_cnct[j][k]];
	//		}
	//		vector<Face3D>::iterator it=find(tmface.begin(),tmface.end(),fc);
	//		int fcid(it-tmface.begin());
	//		if(it==tmface.end())
	//		{
	//			tmface.push_back(fc);
	//		}
	//		tmesh[i].face[j+1]=fcid;
	//		//tmface[fcid].hex.push_back(i);
	//	}
	//	//one top face
	//	Face3D fc2;
	//	for(j=0; j<4; j++)
	//	{
	//		fc2.cnct[j]=tmesh[i].cnct[j+4];
	//		fc2.edge[j]=tmesh[i].edge[j+8];
	//	}
	//	vector<Face3D>::iterator it2=find(tmface.begin(),tmface.end(),fc2);
	//	int fc2id(it2-tmface.begin());
	//	if(it2==tmface.end())
	//	{
	//		tmface.push_back(fc2);
	//	}
	//	tmesh[i].face[5]=fc2id;
	//	//tmface[fc2id].hex.push_back(i);
	//}
	////vertex-to-hex, edge-to-hex, face-to-hex
	//for(i=0; i<tmesh.size(); i++)
	//{
	//	for(j=0; j<8; j++)
	//	{
	//		cp[tmesh[i].cnct[j]].hex.push_back(i);
	//	}
	//	for(j=0; j<12; j++)
	//	{
	//		tmedge[tmesh[i].edge[j]].hex.push_back(i);
	//	}
	//	for(j=0; j<6; j++)
	//	{
	//		tmface[tmesh[i].face[j]].hex.push_back(i);
	//	}
	//}
	////vertex-to-face, edge-to-face
	//for(i=0; i<tmface.size(); i++)
	//{
	//	for(j=0; j<4; j++)
	//	{
	//		cp[tmface[i].cnct[j]].face.push_back(i);
	//		tmedge[tmface[i].edge[j]].face.push_back(i);
	//	}
	//}
	////vertex-to-edge
	//for(i=0; i<tmedge.size(); i++)
	//{
	//	for(j=0; j<2; j++)
	//	{
	//		cp[tmedge[i].pt[j]].edge.push_back(i);
	//	}
	//}

	//find BC face, edge, vertex

	//int ed0[6][4] = { { 4, 5, 6, 7 }, { 1, 3, 9, 11 }, { 0, 2, 8, 10 }, { 1, 3, 9, 11 }, { 0, 2, 8, 10 }, { 4, 5, 6, 7 } };//order could be wrong, but doesn't matter
	for (i = 0; i < tmface.size(); i++)
	{
		if (tmface[i].hex.size() == 1)
		{
			tmface[i].type = 1;
			tmesh[tmface[i].hex[0]].type = 1;
			for (j = 0; j < 4; j++)
			{
				cp[tmface[i].cnct[j]].type = 1;
				tmedge[tmface[i].edge[j]].type = 1;
			}
			//set zero length edges
			//int hexid(tmface[i].hex[0]);
			//int* it=find(tmesh[hexid].face,tmesh[hexid].face+6,i);
			//int fc_loc(it-tmesh[hexid].face);
			//for(j=0; j<4; j++)
			//{
			//	tmedge[tmesh[hexid].edge[ed0[fc_loc][j]]].len=0.;
			//}
		}
	}
	//additional boundary elements
	for (i = 0; i < tmesh.size(); i++)
	{
		if (tmesh[i].type != 1)
		{
			for (j = 0; j < 8; j++)
			{
				if (cp[tmesh[i].cnct[j]].type == 1)
				{
					tmesh[i].type = 1;
					break;
				}
			}
		}
	}
	//find extraordinary edges and vertices
	for (i = 0; i < tmedge.size(); i++)
	{
		if (tmedge[i].type != 1 && tmedge[i].hex.size() != 4)
		{
			tmedge[i].type = 2;
			if (cp[tmedge[i].pt[0]].type != 1)
				cp[tmedge[i].pt[0]].type = 3;
			if (cp[tmedge[i].pt[1]].type != 1)
				cp[tmedge[i].pt[1]].type = 3;
		}
	}
	//boundary
	//for(i=0; i<cp.size(); i++)
	//{
	//	if(cp[i].type==1)//not consider surface extraordinary points yet
	//	{
	//		//int val(0);
	//		//for(j=0; j<cp[i].face.size(); j++)
	//		//{
	//		//	if(tmface[cp[i].face[j]].type==1) val++;
	//		//}
	//		//if(val==3 || val>4) cp[i].type=13;
	//	}
	//}
	//find irregular elements
	for (i = 0; i < tmesh.size(); i++)
	{
		if (tmesh[i].type != 1)
		{
			for (j = 0; j < 12; j++)
			{
				if (tmedge[tmesh[i].edge[j]].type == 2)
				{
					tmesh[i].type = 2;
					break;
				}
			}
			//additional
			for (j = 0; j < 8; j++)
			{
				if (cp[tmesh[i].cnct[j]].type == 3)
				{
					tmesh[i].type = 2;
					break;
				}
			}
		}
	}

	//boundry extraordinary points
	for (i = 0; i < cp.size(); i++)
	{
		if (cp[i].type == 1)
		{
			int count(0);
			for (j = 0; j < cp[i].edge.size(); j++)
			{
				if (tmedge[cp[i].edge[j]].type == 2) count++;
			}
			if (count == 1) cp[i].bcxp = 1;
			else if (count > 1) cp[i].bcxp = 2;
		}
	}

	//check 3D EP
	int n3d(0);
	for (i = 0; i < cp.size(); i++)
	{
		if (cp[i].type == 3)
		{
			int ned(0);
			for (j = 0; j < cp[i].edge.size(); j++)
			{
				if (tmedge[cp[i].edge[j]].type == 2)
				{
					ned++;
				}
			}
			if (ned > 2)
			{
				n3d++;
			}
		}
	}
	cout << "# 3D EP: " << n3d << "\n";
	//getchar();
}

void HexQuality::InitializeMesh(string fn)
{
	//read hex vtk
	string fname(fn + ".vtk"), stmp;
	int npts, neles, itmp;
	ifstream fin;
	fin.open(fname);
	if (fin.is_open())
	{
		for (int i = 0; i < 4; i++) getline(fin, stmp);//skip lines
		fin >> stmp >> npts >> stmp;
		cp.resize(npts);
		for (int i = 0; i < npts; i++)
		{
			cp[i].act = 1;
			fin >> cp[i].coor[0] >> cp[i].coor[1] >> cp[i].coor[2];
			cp[i].coortmp[0] = cp[i].coor[0];
			cp[i].coortmp[1] = cp[i].coor[1];
			cp[i].coortmp[2] = cp[i].coor[2];
		}
		getline(fin, stmp);
		fin >> stmp >> neles >> itmp;
		tmesh.resize(neles);
		for (int i = 0; i < neles; i++)
		{
			tmesh[i].act = 1;
			fin >> itmp >> tmesh[i].cnct[0] >> tmesh[i].cnct[1] >> tmesh[i].cnct[2] >> tmesh[i].cnct[3] >>
				tmesh[i].cnct[4] >> tmesh[i].cnct[5] >> tmesh[i].cnct[6] >> tmesh[i].cnct[7];
		}
		fin.close();
	}
	else
	{
		cerr << "Cannot open " << fname << "!\n";
	}
	InitialConnect();

	hmesh.push_back(tmesh);
	hcp.push_back(cp);
	hface.push_back(tmface);
	hedge.push_back(tmedge);
}

void HexQuality::SetSharpFeature_1(double tol)
{
	int fc_dir[6][3] = { { 0, 3, 2 }, { 0, 1, 5 }, { 1, 2, 6 }, { 3, 7, 6 }, { 0, 4, 7 }, { 4, 5, 6 } };
	uint i, j, k;
	for (i = 0; i < tmedge.size(); i++)
	{
		if (tmedge[i].type == 1)
		{
			vector<array<double, 3>> vec;
			for (j = 0; j < tmedge[i].hex.size(); j++)
			{
				if (tmesh[tmedge[i].hex[j]].type == 1)
				{
					int hxid(tmedge[i].hex[j]);
					for (k = 0; k < 6; k++)
					{
						if (tmface[tmesh[hxid].face[k]].type == 1)
						{
							vector<int>::iterator it = find(tmedge[i].face.begin(), tmedge[i].face.end(), tmesh[hxid].face[k]);
							if (it != tmedge[i].face.end())
							{
								array<double, 3> tmp1, tmp2;
								for (int dof = 0; dof < 3; dof++)
								{
									tmp1[dof] = cp[tmesh[hxid].cnct[fc_dir[k][1]]].coor[dof] - cp[tmesh[hxid].cnct[fc_dir[k][0]]].coor[dof];
									tmp2[dof] = cp[tmesh[hxid].cnct[fc_dir[k][2]]].coor[dof] - cp[tmesh[hxid].cnct[fc_dir[k][1]]].coor[dof];
								}
								array<double, 3> tmp3 = { tmp1[1] * tmp2[2] - tmp1[2] * tmp2[1], -tmp1[0] * tmp2[2] + tmp1[2] * tmp2[0], tmp1[0] * tmp2[1] - tmp1[1] * tmp2[0] };
								double dis = sqrt(tmp3[0] * tmp3[0] + tmp3[1] * tmp3[1] + tmp3[2] * tmp3[2]);
								tmp3[0] /= dis; tmp3[1] /= dis; tmp3[2] /= dis;
								vec.push_back(tmp3);
							}
						}
					}
				}
			}
			if (vec.size() == 2)
			{
				double ang(vec[0][0] * vec[1][0] + vec[0][1] * vec[1][1] + vec[0][2] * vec[1][2]);
				if (ang < tol)
				{
					tmedge[i].sharp = 1;
					cp[tmedge[i].pt[0]].sharp = 1;
					cp[tmedge[i].pt[1]].sharp = 1;
				}
			}
			else
			{
				cerr << "Something wrong in determining sharp edge!\n";
				getchar();
			}
		}
	}

	//int cp_sharp[] = { 1179 };//manual control, heli and heli coarse
	////int cp_sharp[] = { 1543 };//manual control, heli dense loc1
	//int ncpshp(1);
	//for (int i = 0; i < ncpshp; i++)
	//{
	//	cp[cp_sharp[i]].sharp = 1;
	//}
	//for (i = 0; i < tmedge.size(); i++)
	//{
	//	if (tmedge[i].type == 1 && tmedge[i].sharp == 0)
	//	{
	//		if (cp[tmedge[i].pt[0]].sharp == 1 && cp[tmedge[i].pt[1]].sharp == 1 &&
	//			(tmedge[i].pt[0] == cp_sharp[0] || tmedge[i].pt[1] == cp_sharp[0]))//heli coarse
	//		{
	//			tmedge[i].sharp = 1;
	//		}
	//	}
	//}

	////int cp_desharp[] = { 1216 };//muanlly remove certain sharp, heli coarse
	//int cp_desharp[] = { 1580 };//muanlly remove certain sharp, heli coarse loc
	////int cp_desharp[] = { 12003,12022,12016 };//muanlly remove certain sharp, heli loc
	//int ndeshp(1);
	//for (int i = 0; i < ndeshp; i++)
	//{
	//	cp[cp_desharp[i]].sharp = 0;
	//	for (uint j = 0; j < cp[cp_desharp[i]].edge.size(); j++)
	//	{
	//		tmedge[cp[cp_desharp[i]].edge[j]].sharp = 0;
	//	}
	//}

	////int cp_desharp[] = { 7418,7419 };//muanlly remove certain sharp, navair coarse, no pillow
	//int cp_desharp[] = { 13109,13110 };//muanlly remove certain sharp, navair coarse
	////int cp_desharp[] = { 30994,36209 };//muanlly remove certain sharp, navair
	//int ndeshp(2);
	//for (int i = 0; i < ndeshp; i++)
	//{
	//	cp[cp_desharp[i]].sharp = 0;
	//	for (uint j = 0; j < cp[cp_desharp[i]].edge.size(); j++)
	//	{
	//		tmedge[cp[cp_desharp[i]].edge[j]].sharp = 0;
	//	}
	//}

	for (i = 0; i < cp.size(); i++)
	{
		if (cp[i].sharp == 1)
		{
			int nshp(0);
			for (j = 0; j < cp[i].edge.size(); j++)
			{
				if (tmedge[cp[i].edge[j]].sharp == 1)
				{
					nshp++;
				}
			}
			if (nshp >= 3) cp[i].sharp = 2;//sharp corner
			if (nshp == 1) cp[i].sharp = 2;//
		}
	}
	//for (i = 0; i < tmedge.size(); i++)
	//{
	//	if (tmedge[i].sharp == 1)
	//	{
	//		if (cp[tmedge[i].pt[0]].sharp == 0 || cp[tmedge[i].pt[1]].sharp == 0)
	//		{
	//			tmedge[i].sharp = 0;
	//		}
	//	}
	//}
}

void HexQuality::SetSharpFeature_Manual(string fn)
{
	ifstream fin;
	fin.open(fn);
	vector<int> pshp;
	//if (fin.is_open())
	//{
	//	string oneline;
	//	int itmp;
	//	while (getline(fin, oneline))
	//	{
	//		istringstream is(oneline);
	//		is >> itmp;
	//		pshp.push_back(itmp);
	//	}
	//	//int nshp;
	//	//string stmp;
	//	//fin >> stmp >> nshp;
	//	//pshp.resize(nshp);
	//	//for (int i = 0; i < nshp; i++)
	//	//{
	//	//	fin >> pshp[i];
	//	//}
	//	fin.close();
	//}
	//else
	//{
	//	cerr << "Can't open " << fn << "!\n";
	//	return;
	//}

	vector<int> itmp_vec;
	if (fin.is_open())
	{
		string oneline;
		int itmp;
		while (getline(fin, oneline))
		{
			istringstream is(oneline);
			while (is >> itmp)
				itmp_vec.push_back(itmp);
			pshp.push_back(itmp_vec.back());
		}
		//int nshp;
		//string stmp;
		//fin >> stmp >> nshp;
		//pshp.resize(nshp);
		//for (int i = 0; i < nshp; i++)
		//{
		//	fin >> pshp[i];
		//}
		fin.close();
	}
	else
	{
		cerr << "Can't open " << fn << "!\n";
		return;
	}

	for (uint i = 0; i < pshp.size(); i++)
	{
		if (cp[pshp[i]].type == 1)
		{
			cp[pshp[i]].sharp = 1;
		}
	}

	for (uint i = 0; i < tmedge.size(); i++)
	{
		if (tmedge[i].type == 1 && cp[tmedge[i].pt[0]].sharp == 1 && cp[tmedge[i].pt[1]].sharp == 1)
		{
			tmedge[i].sharp = 1;
		}
	}

	for (uint i = 0; i < cp.size(); i++)
	{
		if (cp[i].sharp == 1)
		{
			int nshp(0);
			for (uint j = 0; j < cp[i].edge.size(); j++)
			{
				if (tmedge[cp[i].edge[j]].sharp == 1)
				{
					nshp++;
				}
			}
			if (nshp >= 3 || nshp == 1) cp[i].sharp = 2;//sharp corner
			else if (nshp == 0) cp[i].sharp = 0;
		}
	}
}

void HexQuality::OutputCM(string fn)
{
	//string fname(fn + "_CM.vtk");
	string fname(fn + ".vtk");
	ofstream fout;
	fout.open(fname.c_str());
	if (fout.is_open())
	{
		fout << "# vtk DataFile Version 2.0\nSquare plate test\nASCII\nDATASET UNSTRUCTURED_GRID\n";
		fout << "POINTS " << cp.size() << " float\n";
		for (uint i = 0; i < cp.size(); i++)
		{
			fout << cp[i].coor[0] << " " << cp[i].coor[1] << " " << cp[i].coor[2] << "\n";
		}
		/*fout << "\nCELLS " << tmesh.size() << " " << 9 * tmesh.size() << '\n';
		for (uint i = 0; i<tmesh.size(); i++)
		{
			fout << "8 ";
			for (int j = 0; j<8; j++)
			{
				fout << tmesh[i].cnct[j] << ' ';
			}
			fout << '\n';
		}
		fout << "\nCELL_TYPES " << tmesh.size() << '\n';
		for (uint i = 0; i<tmesh.size(); i++)
		{
			fout << "12\n";
		}*/
		int neles(0);
		for (uint i = 0; i < tmesh.size(); i++)
		{
			//if (tmesh[i].type == 1) neles++;
			if (tmesh[i].act == 1) neles++;
		}
		fout << "\nCELLS " << neles << " " << 9 * neles << '\n';
		for (uint i = 0; i < tmesh.size(); i++)
		{
			//if (tmesh[i].type == 1)
			if (tmesh[i].act == 1)
			{
				fout << "8 ";
				for (int j = 0; j < 8; j++)
				{
					fout << tmesh[i].cnct[j] << ' ';
				}
				fout << '\n';
			}
		}
		fout << "\nCELL_TYPES " << neles << '\n';
		for (uint i = 0; i < neles; i++)
		{
			fout << "12\n";
		}

		fout.close();
	}
	else
	{
		cout << "Cannot open " << fname << "!\n";
	}

	cout << "Done output control mesh!\n"; //getchar();
}


void HexQuality::run_MeshQualityImprove(int mode, int flag_sharp, double tol_sharp, int opt_par1, double opt_par2, string fn)
{
	string fld;
	string fn_out;

	std::filesystem::path path_in = fn;
	fld = path_in.parent_path().string();

	size_t lastindex = fn.find_last_of(".");
	fn_out = fn.substr(0, lastindex);
	cout << "Input Mesh: " << fn << endl;



	InitializeMesh(fn_out);

	//tt3.LaplaceSmoothing(50);
	//tt3.OutputCM(fn_out);
	//cout << "done lap!\n";
	//getchar();

	switch (flag_sharp)
	{
	default:
		cout << "Sharp feature OFF" << endl;
		break;
	case 1:
		cout << "Automatically detect sharp feature ON: TOL = " << tol_sharp << endl;
		SetSharpFeature_1(tol_sharp);
		break;
	case 2:
		cout << "Manually apply sharp feature ON" << endl;
		SetSharpFeature_Manual(fld + "\\sharp.txt");
		break;
	}

	////tt3.SetSharpFeature_1();//automatic way may not work
	//string fn1(fld + "navair_sharp");
	//tt3.OutputEdge(fn1);
	//tt3.OutputCM(fn1);
	//cout << "done setting sharp feature\n";
	//getchar();

	//tt3.CheckJacobBEXT();

	//mesh quality improvement includes three parts (Pillow, Smoothing and Optimizizing) and they should be run one by one

	switch (mode)
	{
	case 0: //Laplace Smoothing
		LaplaceSmoothing(opt_par1);
		OutputCM(fn_out + "_lap");
		cout << "done lap!\n";
		break;
	case 1: //Pillowing
		///tt3.Pillow(1);//run Pillow first if needed, which outputs an mesh with an outer layer (input argument "1" means one layer) added to the input mesh.
					  ///the name of the outpue mesh can be changed within the function definition
		Pillow(opt_par1);
		OutputCM(fn_out + "_pillow");
		break;
	case 2: //Smoothing
		///tt3.Smoothing(50, 0.4);//then run Smoothing on the pillowed mesh by specifying the
		///					   //number of steps (e.g. 350) and the step size (e.g. 0.5), which
		///					   //iteratively moves each point in the mesh to its mass center.
		///					   //The input "fn_in" should be changed to the pillowed mesh or leave it be if smoothing is applied without pillowing
		///					   //It outputs a smoothed mesh, see the commented code "OutputCM" in the definition
		Smoothing(opt_par1, opt_par2);
		//Smoothing_Angran(opt_par1, opt_par2);
		OutputCM(fn_out + "_smooth");
		break;
	case 3: //Optimizing
		///tt3.Optimizing(50, 0.1);//last run Optimizing if there is still negative Jacobian elements after smoothing
		///                        //specifying the number of steps (e.g. 100) and the step size (e.g. 0.01),
		///						//which in each iteration quality improves the element of the worst Jacobian.
		///						//The input "fn_in" should be changed to the smoothed mesh.
		///						//It outputs an optimized mesh, check the commented code "OutputCM" in the function definition.
		Optimizing(opt_par1, opt_par2);
		OutputCM(fn_out + "_opt");
		break;
	}

}


int HexQuality::IsPillowNeeded()
{
	int pillow_flag(0);
	uint i, j;
	for (i = 0; i<tmface.size(); i++)
	{
		if (tmface[i].hex.size() == 1)
		{
			tmface[i].type = 1;
			tmesh[tmface[i].hex[0]].type = 1;
			for (j = 0; j<4; j++)
			{
				cp[tmface[i].cnct[j]].type = 1;
				tmedge[tmface[i].edge[j]].type = 1;
			}
		}
	}
	for (i = 0; i < tmesh.size(); i++)//boundary elements
	{
		if (tmesh[i].type == 0)
		{
			for (j = 0; j < 8; j++)
			{
				if (cp[tmesh[i].cnct[j]].type == 1)
				{
					pillow_flag = 2;
					tmesh[i].type = 1;
					break;
				}
			}
		}
		if (pillow_flag == 2) break;
	}
	//find extraordinary edges and vertices
	for (i = 0; i<tmedge.size(); i++)
	{
		if (tmedge[i].type == 1 && (tmedge[i].hex.size() == 1 || tmedge[i].hex.size() > 2))
		{
			pillow_flag = 2; break;
		}
	}
	if (pillow_flag == 0)
	{
		for (i = 0; i < cp.size(); i++)
		{
			if (cp[i].type == 0)
			{
				int nb[2] = { 0,0 };
				for (j = 0; j < cp[i].hex.size(); j++)
				{
					if (tmesh[cp[i].hex[j]].type == 0)
					{
						nb[0]++;
					}
					else
					{
						nb[1]++;
					}
				}
				if (nb[0] == nb[1])
				{
					pillow_flag = 1; break;
				}
			}
		}
	}
	return pillow_flag;
}

void HexQuality::Pillow(int nlayer)
{
	
	cout << "Pillowing" << endl;
	//input is a hex mesh
	int fcloc[6][4] = { {0,3,2,1},{0,1,5,4},{1,2,6,5},{2,3,7,6},{0,4,7,3},{4,5,6,7} };
	//initial layer
	vector<int> pts0;
	vector<int> pid(cp.size(), -1);
	vector<array<int, 4>> layer0;
	array<int, 4> fctmp;
	for (uint i = 0; i < cp.size(); i++)
	{
		if (cp[i].type == 1)//boundary point
		{
			pts0.push_back(i);
			pid[i] = pts0.size() - 1;
		}
	}
	for (uint i = 0; i < tmesh.size(); i++)
	{
		if (tmesh[i].type == 1)//boundary element
		{
			for (int j = 0; j < 6; j++)
			{
				if (tmface[tmesh[i].face[j]].type == 1)//boundary face
				{
					for (int k = 0; k < 4; k++)
					{
						fctmp[k] = pid[tmesh[i].cnct[fcloc[j][k]]];
					}
					layer0.push_back(fctmp);
				}
			}
		}
	}

	vector<vector<int>> p2f(pts0.size());
	vector<double> dist;
	vector<array<double, 3>> dir;
	for (uint i = 0; i < layer0.size(); i++)
	{
		for (int j = 0; j < 4; j++)
		{
			p2f[layer0[i][j]].push_back(i);
		}
	}
	for (uint i = 0; i < pts0.size(); i++)
	{
		double dtmp, dmin(1.e5);
		array<double, 3> vtmp;
		for (uint j = 0; j < cp[pts0[i]].edge.size(); j++)
		{
			if (tmedge[cp[pts0[i]].edge[j]].type == 1)
			{
				int edpt[2] = { tmedge[cp[pts0[i]].edge[j]].pt[0],tmedge[cp[pts0[i]].edge[j]].pt[1] };
				vtmp[0] = cp[edpt[1]].coor[0] - cp[edpt[0]].coor[0];
				vtmp[1] = cp[edpt[1]].coor[1] - cp[edpt[0]].coor[1];
				vtmp[2] = cp[edpt[1]].coor[2] - cp[edpt[0]].coor[2];
				dtmp = sqrt(vtmp[0] * vtmp[0] + vtmp[1] * vtmp[1] + vtmp[2] * vtmp[2]);
				if (dtmp < dmin) dmin = dtmp;
				//if (dtmp < 1.e-10)
				//{
				//	cout << "zero edge length!\n"; getchar();
				//}
			}
		}
		dist.push_back(dmin);
		array<double, 3> nm = { 0.,0.,0. };
		for (uint j = 0; j < p2f[i].size(); j++)
		{
			int fcid(p2f[i][j]);
			array<int,4>::iterator it = find(layer0[fcid].begin(),layer0[fcid].end(),i);
			int loc(it - layer0[fcid].begin());
			if (loc == 4)
			{
				cerr << "can't find!\n"; getchar();
			}
			int edpt[3] = { pts0[layer0[fcid][loc]],pts0[layer0[fcid][(loc + 1) % 4]],pts0[layer0[fcid][(loc + 3) % 4]] };
			double v1[3] = { cp[edpt[1]].coor[0] - cp[edpt[0]].coor[0],cp[edpt[1]].coor[1] - cp[edpt[0]].coor[1],cp[edpt[1]].coor[2] - cp[edpt[0]].coor[2] };
			double v2[3] = { cp[edpt[2]].coor[0] - cp[edpt[0]].coor[0],cp[edpt[2]].coor[1] - cp[edpt[0]].coor[1],cp[edpt[2]].coor[2] - cp[edpt[0]].coor[2] };
			double v3[3] = { v1[1] * v2[2] - v1[2] * v2[1],v1[2] * v2[0] - v1[0] * v2[2], v1[0] * v2[1] - v1[1] * v2[0] };
			dtmp = sqrt(v3[0] * v3[0] + v3[1] * v3[1] + v3[2] * v3[2]);
			v3[0] /= dtmp; v3[1] /= dtmp; v3[2] /= dtmp;
			nm[0] += v3[0];
			nm[1] += v3[1];
			nm[2] += v3[2];
		}
		nm[0] /= double(p2f[i].size());
		nm[1] /= double(p2f[i].size());
		nm[2] /= double(p2f[i].size());
		dir.push_back(nm);
	}
	
	for (int il = 0; il < nlayer; il++)//now only consider nlayer==1
	{
		vector<int> pts1(pts0.size());
		for (uint i = 0; i < pts0.size(); i++)
		{
			Vertex3D ptmp;
			ptmp.act = 1;
			ptmp.coor[0] = cp[pts0[i]].coor[0];
			ptmp.coor[1] = cp[pts0[i]].coor[1];
			ptmp.coor[2] = cp[pts0[i]].coor[2];
			cp.push_back(ptmp);
			pts1[i] = cp.size() - 1;
		}
		for (uint i = 0; i < layer0.size(); i++)
		{
			Element3D etmp;
			etmp.act = 1;
			for (int j = 0; j < 4; j++)
			{
				etmp.cnct[j] = pts0[layer0[i][j]];
				etmp.cnct[j + 4] = pts1[layer0[i][j]];
			}
			tmesh.push_back(etmp);
		}
		//for (uint i = 0; i < pts0.size(); i++)
		//{
		//	pts0[i] = pts1[i];
		//}
	}

	double eta(0.2);
	double smin(1.e5);
	for (uint i = 0; i < pts0.size(); i++)
	{
		double step(eta*dist[i]);
		cp[pts0[i]].coor[0] -= step*dir[i][0];
		cp[pts0[i]].coor[1] -= step*dir[i][1];
		cp[pts0[i]].coor[2] -= step*dir[i][2];
		if (dist[i] < smin) smin = dist[i];
	}
	cout << "min edge len: " << smin << "\n";

	cout << "Done pillowing!\n";
	//OutputCM("../io/hex_input/pillow/cube_coarse");
	//OutputCM("../io/hex_input/pillow/fertility1");
	//OutputCM("../io/hex_input/pillow/rockerarm1");
	//OutputCM("../io/hex_input/pillow/honda1m");
	//OutputCM("../io/hex_input/pillow/honda2");
	//OutputCM("../io/hex_input/pillow/honda2_1");
	//OutputCM("../io/hex_input/pillow/rod");
	//OutputCM("../io/hex_input/pillow/navair_coarse");
}

void HexQuality::Smoothing(int nSize, double stepSize)
{
	cout << "\nSmoothing...\n";

	//int nLapStep(50);
	//for (int it = 0; it < nLapStep; it++)
	//{
	//	for (uint i = 0; i < cp.size(); i++)
	//	{
	//		if (cp[i].type != 1)//interior points only
	//		{
	//			int tmp = LaplaceSmoothingPoint(i, stepSize);
	//		}
	//	}
	//	double minJacob_glb;
	//	int min_pos;
	//	vector<int> BadEle;
	//	GlobalMinJacob(minJacob_glb, min_pos, BadEle);
	//}
	////OutputCM("../io/NAVAIR_GEM/output/smooth1/navair3_0");
	////getchar();

	for (int it = 0; it < nSize; it++)
	{
		cout << "it: " << it << "\n";
		int flag(0);

		double minJacob_glb0;
		int min_pos0(0);
		vector<int> BadEle0;
		GlobalMinJacob(minJacob_glb0, min_pos0, BadEle0);

		//for (uint i0 = 0; i0 < BadEle0.size(); i0++)
		//{
		//	//cout << BadEle0[i0] << " ";
		//	tmesh[BadEle0[i0]].trun = 1;
		//	//int eid(BadEle[i0]);
		//	//for (int i = 0; i < 8; i++)
		//	//{
		//	//	if (cp[tmesh[eid].cnct[i]].type != 1)
		//	//	{
		//	//		int tmp = SmoothingPoint(tmesh[eid].cnct[i], stepSize);
		//	//		if (tmp == 1) flag = 1;
		//	//	}
		//	//}
		//}
		//break;

		for (uint i = 0; i < cp.size(); i++)
		{
			if (cp[i].type != 1)//interior points
			{
				int tmp = SmoothingPoint(i, stepSize);
				if (tmp == 1) flag = 1;
			}
			else if (cp[i].sharp == 0)
			{
				int tmp = SmoothingPointBoundary(i, stepSize);
				if (tmp == 1) flag = 1;
			}
			else if (cp[i].sharp == 1)
			{
				int tmp = SmoothingPointBoundarySharp(i, stepSize);
				if (tmp == 1) flag = 1;
			}
		}
		double minJacob_glb;
		int min_pos;
		vector<int> BadEle;
		GlobalMinJacob(minJacob_glb, min_pos, BadEle);

		for (uint i = 0; i < tmesh.size(); i++) tmesh[i].trun = 0;
		for (uint i = 0; i < BadEle.size(); i++) tmesh[BadEle[i]].trun = 1;

		//if (flag == 0 || minJacob_glb0 > minJacob_glb)
		//{
		//	cout << "After " << it << "steps, no interior point moved!\n";
		//	break;
		//}
		
	}
	cout << "Done smoothing!\n";
}

int HexQuality::LaplaceSmoothingPoint(int pid, double stepSize)
{
	double ptmp[3] = { 0.,0.,0., };
	for (uint i = 0; i < cp[pid].edge.size(); i++)
	{
		int edid(cp[pid].edge[i]);
		int pnb(tmedge[edid].pt[0]);
		if (pnb == pid) pnb = tmedge[edid].pt[1];
		ptmp[0] += cp[pnb].coor[0];
		ptmp[1] += cp[pnb].coor[1];
		ptmp[2] += cp[pnb].coor[2];
	}
	cp[pid].coor[0] = ptmp[0] / cp[pid].edge.size();
	cp[pid].coor[1] = ptmp[1] / cp[pid].edge.size();
	cp[pid].coor[2] = ptmp[2] / cp[pid].edge.size();
	
	return 1;
}

int HexQuality::SmoothingPoint(int pid, double stepSize)
{
	const double eps(1.e-6);
	double vol_all(0.), center[3] = { 0.,0.,0. }, avg[3] = { 0.,0.,0. }, dir[3];
	double pold[3] = { cp[pid].coor[0],cp[pid].coor[1],cp[pid].coor[2] };

	for (uint i = 0; i < cp[pid].hex.size(); i++)
	{
		double e_vol, e_center[3];
		GetHexVolAndCenter(cp[pid].hex[i], e_vol, e_center);
		vol_all += e_vol;
		center[0] += e_center[0] * e_vol;
		center[1] += e_center[1] * e_vol;
		center[2] += e_center[2] * e_vol;
		avg[0] += e_center[0];
		avg[1] += e_center[1];
		avg[2] += e_center[2];
	}
	avg[0] /= double(cp[pid].hex.size());
	avg[1] /= double(cp[pid].hex.size());
	avg[2] /= double(cp[pid].hex.size());

	if (vol_all < eps)
	{
		dir[0] = avg[0] - cp[pid].coor[0];
		dir[1] = avg[1] - cp[pid].coor[1];
		dir[2] = avg[2] - cp[pid].coor[2];
	}
	else
	{
		center[0] /= vol_all;
		center[1] /= vol_all;
		center[2] /= vol_all;
		dir[0] = center[0] - cp[pid].coor[0];
		dir[1] = center[1] - cp[pid].coor[1];
		dir[2] = center[2] - cp[pid].coor[2];
	}

	//before smoothing
	double minJacob_0(1.e5), minJacob_1(1.e5), tmpJacob, GaussPos[3];
	int nBad0(0), nBad1(0);
	for (uint i = 0; i < cp[pid].hex.size(); i++)
	{
		GetHexMinJacob(cp[pid].hex[i], tmpJacob, GaussPos);
		if (tmpJacob < minJacob_0) minJacob_0 = tmpJacob;
		if (tmpJacob < 0.) nBad0++;
	}
	//double edlen = AverageEdgeLength(pid);
	//double stepSize_scale(stepSize*edlen);
	double stepSize_scale(stepSize);
	//after smoothing
	cp[pid].coor[0] = cp[pid].coor[0] + stepSize_scale*dir[0];
	cp[pid].coor[1] = cp[pid].coor[1] + stepSize_scale*dir[1];
	cp[pid].coor[2] = cp[pid].coor[2] + stepSize_scale*dir[2];
	for (uint i = 0; i < cp[pid].hex.size(); i++)
	{
		GetHexMinJacob(cp[pid].hex[i], tmpJacob, GaussPos);
		if (tmpJacob < minJacob_1) minJacob_1 = tmpJacob;
		if (tmpJacob < 0.) nBad1++;
	}
	if ((minJacob_0 < 0. && minJacob_1 < minJacob_0) || (nBad1 > nBad0))
	{
		cp[pid].coor[0] = pold[0];
		cp[pid].coor[1] = pold[1];
		cp[pid].coor[2] = pold[2];
		return 0;
	}
	return 1;
}

int HexQuality::SmoothingPointBoundary(int pid, double stepSize)
{
	const double eps(1.e-6);
	double area_all(0.), center[3] = { 0.,0.,0. }, nm[3] = { 0., 0., 0. };
	int nbf(0);
	for (uint i = 0; i < cp[pid].face.size(); i++)
	{
		if (tmface[cp[pid].face[i]].type == 1)
		{
			double e_area, e_center[3], e_nm[3];
			int* it = find(tmface[cp[pid].face[i]].cnct, tmface[cp[pid].face[i]].cnct + 4, pid);
			int loc(it - tmface[cp[pid].face[i]].cnct);
			GetQuadInfo(cp[pid].face[i], loc, e_area, e_center, e_nm);
			area_all += e_area;
			center[0] += e_center[0] * e_area;
			center[1] += e_center[1] * e_area;
			center[2] += e_center[2] * e_area;
			nm[0] += e_nm[0] * e_area;
			nm[1] += e_nm[1] * e_area;
			nm[2] += e_nm[2] * e_area;
			//nm[0] += e_nm[0];
			//nm[1] += e_nm[1];
			//nm[2] += e_nm[2];
			nbf++;
		}
	}
	if (area_all < eps)
	{
		cerr << "Element with almost zero area!\n";
		return 0;
	}
	center[0] /= area_all;	center[1] /= area_all;	center[2] /= area_all;
	nm[0] /= area_all;	nm[1] /= area_all;	nm[2] /= area_all;
	//nm[0] /= double(nbf);	nm[1] /= double(nbf);	nm[2] /= double(nbf);
	double dst = sqrt(nm[0] * nm[0] + nm[1] * nm[1] + nm[2] * nm[2]);
	nm[0] /= dst; nm[1] /= dst; nm[2] /= dst;
	double dir[3] = { center[0] - cp[pid].coor[0],center[1] - cp[pid].coor[1], center[2] - cp[pid].coor[2] };
	dst = dir[0] * nm[0] + dir[1] * nm[1] + dir[2] * nm[2];
	dir[0] = dir[0] - dst*nm[0];
	dir[1] = dir[1] - dst*nm[1];
	dir[2] = dir[2] - dst*nm[2];

	//dst = dir[0] * nm[0] + dir[1] * nm[1] + dir[2] * nm[2];
	//cout << "inner product: " << dst << "\n";
	//getchar();

	double pold[3] = { cp[pid].coor[0],cp[pid].coor[1],cp[pid].coor[2] };

	//before smoothing
	double minJacob_0(1.e5), minJacob_1(1.e5), tmpJacob, GaussPos[3];
	int nBad0(0), nBad1(0);
	for (uint i = 0; i < cp[pid].hex.size(); i++)
	{
		GetHexMinJacob(cp[pid].hex[i], tmpJacob, GaussPos);
		if (tmpJacob < minJacob_0) minJacob_0 = tmpJacob;
		if (tmpJacob < 0.) nBad0++;
	}
	//double edlen = AverageEdgeLength(pid);
	//double stepSize_scale(stepSize*edlen);
	double stepSize_scale(stepSize);
	//after smoothing
	cp[pid].coor[0] = cp[pid].coor[0] + stepSize_scale*dir[0];
	cp[pid].coor[1] = cp[pid].coor[1] + stepSize_scale*dir[1];
	cp[pid].coor[2] = cp[pid].coor[2] + stepSize_scale*dir[2];
	for (uint i = 0; i < cp[pid].hex.size(); i++)
	{
		GetHexMinJacob(cp[pid].hex[i], tmpJacob, GaussPos);
		if (tmpJacob < minJacob_1) minJacob_1 = tmpJacob;
		if (tmpJacob < 0.) nBad1++;
	}
	if ((minJacob_0 < 0. && minJacob_1 < minJacob_0) || nBad1 > nBad0)
	{
		cp[pid].coor[0] = pold[0];
		cp[pid].coor[1] = pold[1];
		cp[pid].coor[2] = pold[2];
		return 0;
	}
	return 1;
}

int HexQuality::SmoothingPointBoundarySharp(int pid, double stepSize)
{
	const double eps(1.e-6);
	double len_all(0.), center[3] = { 0.,0.,0. }, ldir[3] = { 0., 0., 0. };
	int ned(0);
	for (uint i = 0; i < cp[pid].edge.size(); i++)
	{
		if (tmedge[cp[pid].edge[i]].sharp == 1)
		{
			int edid(cp[pid].edge[i]);
			int edpt[2] = { tmedge[edid].pt[0],tmedge[edid].pt[1] };
			if (edpt[1] == pid)
			{
				edpt[0] = tmedge[edid].pt[1];
				edpt[1] = tmedge[edid].pt[0];
			}
			if (ned == 1)
			{
				int itmp(edpt[0]);
				edpt[0] = edpt[1];
				edpt[1] = itmp;
			}
			double vec[3] = { cp[edpt[1]].coor[0] - cp[edpt[0]].coor[0],cp[edpt[1]].coor[1] - cp[edpt[0]].coor[1],
				cp[edpt[1]].coor[2] - cp[edpt[0]].coor[2] };
			double edlen = sqrt(vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2]);
			double midpt[3] = { (cp[edpt[1]].coor[0] + cp[edpt[0]].coor[0]) / 2.,(cp[edpt[1]].coor[1] + cp[edpt[0]].coor[1]) / 2.,
				(cp[edpt[1]].coor[2] + cp[edpt[0]].coor[2]) / 2. };
			len_all += edlen;
			center[0] += midpt[0] * edlen;
			center[1] += midpt[1] * edlen;
			center[2] += midpt[2] * edlen;
			ldir[0] += vec[0] * edlen;
			ldir[1] += vec[1] * edlen;
			ldir[2] += vec[2] * edlen;
			ned++;
		}
	}
	if (len_all < eps)
	{
		cerr << "Edge with almost zero length!\n";
		return 0;
	}
	center[0] /= len_all;	center[1] /= len_all;	center[2] /= len_all;
	ldir[0] /= len_all;	ldir[1] /= len_all;	ldir[2] /= len_all;
	//nm[0] /= double(nbf);	nm[1] /= double(nbf);	nm[2] /= double(nbf);
	double dst = sqrt(ldir[0] * ldir[0] + ldir[1] * ldir[1] + ldir[2] * ldir[2]);
	ldir[0] /= dst; ldir[1] /= dst; ldir[2] /= dst;
	double dir[3] = { center[0] - cp[pid].coor[0],center[1] - cp[pid].coor[1], center[2] - cp[pid].coor[2] };
	dst = dir[0] * ldir[0] + dir[1] * ldir[1] + dir[2] * ldir[2];
	dir[0] = dst*ldir[0];
	dir[1] = dst*ldir[1];
	dir[2] = dst*ldir[2];

	//dst = dir[0] * nm[0] + dir[1] * nm[1] + dir[2] * nm[2];
	//cout << "inner product: " << dst << "\n";
	//getchar();

	double pold[3] = { cp[pid].coor[0],cp[pid].coor[1],cp[pid].coor[2] };

	//before smoothing
	double minJacob_0(1.e5), minJacob_1(1.e5), tmpJacob, GaussPos[3];
	for (uint i = 0; i < cp[pid].hex.size(); i++)
	{
		GetHexMinJacob(cp[pid].hex[i], tmpJacob, GaussPos);
		if (tmpJacob < minJacob_0) minJacob_0 = tmpJacob;
	}
	//double edlen = AverageEdgeLength(pid);
	//double stepSize_scale(stepSize*edlen);
	double stepSize_scale(stepSize);
	//after smoothing
	cp[pid].coor[0] = cp[pid].coor[0] + stepSize_scale*dir[0];
	cp[pid].coor[1] = cp[pid].coor[1] + stepSize_scale*dir[1];
	cp[pid].coor[2] = cp[pid].coor[2] + stepSize_scale*dir[2];
	for (uint i = 0; i < cp[pid].hex.size(); i++)
	{
		GetHexMinJacob(cp[pid].hex[i], tmpJacob, GaussPos);
		if (tmpJacob < minJacob_1) minJacob_1 = tmpJacob;
	}
	if (minJacob_1 < minJacob_0)
	{
		cp[pid].coor[0] = pold[0];
		cp[pid].coor[1] = pold[1];
		cp[pid].coor[2] = pold[2];
		return 0;
	}
	return 1;
}

void HexQuality::GetHexMinJacob(int eid, double& minJacob_ele, double GaussPos[3])
{
	vector<double> Gpt, wght;
	GetGaussPoint(2, Gpt, wght);
	Gpt[0] = 0.; Gpt[1] = 1.;
	double detJ;
	minJacob_ele = 1.e5;
	GaussPos[0] = 0.; GaussPos[1] = 0.; GaussPos[2] = 0.;
	for (uint i = 0; i < Gpt.size(); i++)
	{
		for (uint j = 0; j < Gpt.size(); j++)
		{
			for (uint k = 0; k < Gpt.size(); k++)
			{
				//JacobEval(eid, Gpt[i], Gpt[j], Gpt[k], detJ);
				JacobEval_Scale(eid, Gpt[i], Gpt[j], Gpt[k], detJ);
				if (detJ < minJacob_ele)
				{
					minJacob_ele = detJ;
					GaussPos[0] = Gpt[i];
					GaussPos[1] = Gpt[j];
					GaussPos[2] = Gpt[k];
				}
			}
		}
	}
}

void HexQuality::GetHexVolAndCenter(int eid, double& vol, double center[3])
{
	vol = 0.;
	center[0] = 0.; center[1] = 0.; center[2] = 0.;
	for (int i = 0; i < 8; i++)
	{
		center[0] += cp[tmesh[eid].cnct[i]].coor[0];
		center[1] += cp[tmesh[eid].cnct[i]].coor[1];
		center[2] += cp[tmesh[eid].cnct[i]].coor[2];
	}
	center[0] /= 8.; center[1] /= 8.; center[2] /= 8.;

	vector<double> Gpt, wght;
	GetGaussPoint(2, Gpt, wght);
	double detJ;
	for (uint i = 0; i < Gpt.size(); i++)
	{
		for (uint j = 0; j < Gpt.size(); j++)
		{
			for (uint k = 0; k < Gpt.size(); k++)
			{
				JacobEval(eid, Gpt[i], Gpt[j], Gpt[k], detJ);
				vol += wght[i] * wght[j] * wght[k] * detJ;
			}
		}
	}
}

void HexQuality::GetQuadInfo(int fcid, int ploc, double& area, double center[3], double nm[3])
{
	//face is counter-clock-wise oriented
	area = 0.;
	center[0] = 0.; center[1] = 0.; center[2] = 0.;
	nm[0] = 0.; nm[1] = 0.; nm[2] = 0.;
	for (int i = 0; i < 4; i++)
	{
		center[0] += cp[tmface[fcid].cnct[i]].coor[0];
		center[1] += cp[tmface[fcid].cnct[i]].coor[1];
		center[2] += cp[tmface[fcid].cnct[i]].coor[2];
	}
	center[0] /= 4.; center[1] /= 4.; center[2] /= 4.;

	int tp[2][3] = { { tmface[fcid].cnct[ploc],tmface[fcid].cnct[(ploc + 1) % 4],tmface[fcid].cnct[(ploc + 3) % 4] },
	{ tmface[fcid].cnct[(ploc + 2) % 4],tmface[fcid].cnct[(ploc + 3) % 4],tmface[fcid].cnct[(ploc + 1) % 4] } };
	double v0[2][3] = { {cp[tp[0][1]].coor[0] - cp[tp[0][0]].coor[0],cp[tp[0][1]].coor[1] - cp[tp[0][0]].coor[1],cp[tp[0][1]].coor[2] - cp[tp[0][0]].coor[2] },
	{ cp[tp[0][2]].coor[0] - cp[tp[0][0]].coor[0],cp[tp[0][2]].coor[1] - cp[tp[0][0]].coor[1],cp[tp[0][2]].coor[2] - cp[tp[0][0]].coor[2] } };
	double v1[2][3] = { { cp[tp[1][1]].coor[0] - cp[tp[1][0]].coor[0],cp[tp[1][1]].coor[1] - cp[tp[1][0]].coor[1],cp[tp[1][1]].coor[2] - cp[tp[1][0]].coor[2] },
	{ cp[tp[1][2]].coor[0] - cp[tp[1][0]].coor[0],cp[tp[1][2]].coor[1] - cp[tp[1][0]].coor[1],cp[tp[1][2]].coor[2] - cp[tp[1][0]].coor[2] } };
	double vc0[3] = { v0[0][1] * v0[1][2] - v0[0][2] * v0[1][1],v0[0][2] * v0[1][0] - v0[0][0] * v0[1][2],v0[0][0] * v0[1][1] - v0[0][1] * v0[1][0] };
	double vc1[3] = { v1[0][1] * v1[1][2] - v1[0][2] * v1[1][1],v1[0][2] * v1[1][0] - v1[0][0] * v1[1][2],v1[0][0] * v1[1][1] - v1[0][1] * v1[1][0] };
	double dst0 = sqrt(vc0[0] * vc0[0] + vc0[1] * vc0[1] + vc0[2] * vc0[2]);
	double dst1 = sqrt(vc1[0] * vc1[0] + vc1[1] * vc1[1] + vc1[2] * vc1[2]);

	area = (dst0 + dst1) / 2.;
	nm[0] = vc0[0] / dst0; nm[1] = vc0[1] / dst0; nm[2] = vc0[2] / dst0;
}

void HexQuality::GetGaussPoint(int ng, vector<double>& Gpt, vector<double>& wght)
{
	//already map [-1,1] to [0,1]
	Gpt.clear();
	wght.clear();
	switch (ng)
	{
	case 2:
	{
		Gpt.resize(ng);
		wght.resize(ng);
		Gpt[0] = 0.2113248654051871;			Gpt[1] = 0.7886751345948129;
		wght[0] = 1.;			wght[1] = 1.;
		break;
	}
	case 3:
	{
		Gpt.resize(ng);
		wght.resize(ng);
		Gpt[0] = 0.1127016653792583;			Gpt[1] = 0.5;			Gpt[2] = 0.8872983346207417;
		wght[0] = 0.5555555555555556;			wght[1] = 0.8888888888888889;			wght[2] = 0.5555555555555556;
		break;
	}
	case 4:
	{
		Gpt.resize(ng);
		wght.resize(ng);
		Gpt[0] = 0.06943184420297371;			Gpt[1] = 0.33000947820757187;			Gpt[2] = 0.6699905217924281;			Gpt[3] = 0.9305681557970262;
		wght[0] = 0.3478548451374539;			wght[1] = 0.6521451548625461;			wght[2] = 0.6521451548625461;			wght[3] = 0.3478548451374539;
		break;
	}
	case 5:
	{
		Gpt.resize(ng);
		wght.resize(ng);
		Gpt[0] = 0.046910077030668;			Gpt[1] = 0.2307653449471585;			Gpt[2] = 0.5;			Gpt[3] = 0.7692346550528415;  Gpt[4] = 0.953089922969332;
		wght[0] = 0.2369268850561891;			wght[1] = 0.4786286704993665;			wght[2] = 0.5688888888888889;			wght[3] = 0.4786286704993665; wght[4] = 0.2369268850561891;
		break;
	}
	default:
	{
		Gpt.resize(2);
		wght.resize(2);
		Gpt[0] = 0.2113248654051871;			Gpt[1] = 0.7886751345948129;
		wght[0] = 1.;			wght[1] = 1.;
		break;
	}
	}
}

void HexQuality::JacobEval(int eid, double u, double v, double w, double& detJ)
{
	double Nu[2] = { 1. - u,u };
	double Nv[2] = { 1. - v,v };
	double Nw[2] = { 1. - w,w };
	double dNdu[2] = { -1.,1. };
	double dNdv[2] = { -1.,1. };
	double dNdw[2] = { -1.,1. };
	double dNdt[8][3];
	int ploc[8] = { 0,1,3,2,4,5,7,6 };
	int loc(0);
	for (int i = 0; i < 2; i++)
	{
		for (int j = 0; j < 2; j++)
		{
			for (int k = 0; k < 2; k++)
			{
				dNdt[ploc[loc]][0] = dNdu[k] * Nv[j] * Nw[i];
				dNdt[ploc[loc]][1] = Nu[k] * dNdv[j] * Nw[i];
				dNdt[ploc[loc]][2] = Nu[k] * Nv[j] * dNdw[i];
				loc++;
			}
		}
	}
	double dxdt[3][3] = { {0.,0.,0.},{ 0.,0.,0. },{ 0.,0.,0. } };
	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			for (int k = 0; k < 8; k++)
			{
				dxdt[i][j] += cp[tmesh[eid].cnct[k]].coor[i] * dNdt[k][j];
			}
		}
	}
	detJ = dxdt[0][0] * (dxdt[1][1] * dxdt[2][2] - dxdt[1][2] * dxdt[2][1]) -
		dxdt[0][1] * (dxdt[1][0] * dxdt[2][2] - dxdt[1][2] * dxdt[2][0]) +
		dxdt[0][2] * (dxdt[1][0] * dxdt[2][1] - dxdt[1][1] * dxdt[2][0]);
	detJ *= 0.125;
}

void HexQuality::JacobEval_Scale(int eid, double u, double v, double w, double& detJ)
{
	double Nu[2] = { 1. - u,u };
	double Nv[2] = { 1. - v,v };
	double Nw[2] = { 1. - w,w };
	double dNdu[2] = { -1.,1. };
	double dNdv[2] = { -1.,1. };
	double dNdw[2] = { -1.,1. };
	double dNdt[8][3];
	int ploc[8] = { 0,1,3,2,4,5,7,6 };
	int loc(0);
	for (int i = 0; i < 2; i++)
	{
		for (int j = 0; j < 2; j++)
		{
			for (int k = 0; k < 2; k++)
			{
				dNdt[ploc[loc]][0] = dNdu[k] * Nv[j] * Nw[i];
				dNdt[ploc[loc]][1] = Nu[k] * dNdv[j] * Nw[i];
				dNdt[ploc[loc]][2] = Nu[k] * Nv[j] * dNdw[i];
				loc++;
			}
		}
	}
	double dxdt[3][3] = { { 0.,0.,0. },{ 0.,0.,0. },{ 0.,0.,0. } };
	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			for (int k = 0; k < 8; k++)
			{
				dxdt[i][j] += cp[tmesh[eid].cnct[k]].coor[i] * dNdt[k][j];
			}
		}
	}
	detJ = dxdt[0][0] * (dxdt[1][1] * dxdt[2][2] - dxdt[1][2] * dxdt[2][1]) -
		dxdt[0][1] * (dxdt[1][0] * dxdt[2][2] - dxdt[1][2] * dxdt[2][0]) +
		dxdt[0][2] * (dxdt[1][0] * dxdt[2][1] - dxdt[1][1] * dxdt[2][0]);
	double nm[3] = { sqrt(dxdt[0][0] * dxdt[0][0] + dxdt[1][0] * dxdt[1][0] + dxdt[2][0] * dxdt[2][0]),
		sqrt(dxdt[0][1] * dxdt[0][1] + dxdt[1][1] * dxdt[1][1] + dxdt[2][1] * dxdt[2][1]), 
		sqrt(dxdt[0][2] * dxdt[0][2] + dxdt[1][2] * dxdt[1][2] + dxdt[2][2] * dxdt[2][2]) };
	double tmp(nm[0] * nm[1] * nm[2]);
	if (tmp > 1.e-12)
	{
		detJ /= tmp;
	}
	//double eps(1.e-12);
	//if (detJ < eps)
	//{
	//	double sumtmp(0.);
	//	for (int i = 0; i < 3; i++)
	//	{
	//		for (int j = 0; j < 3; j++)
	//		{
	//			sumtmp += dxdt[i][j] * dxdt[i][j];
	//		}
	//	}
	//	detJ = sqrt(sumtmp);
	//}
}

void HexQuality::GlobalMinJacob(double& minJacob_glb, int& min_pos, vector<int>& BadEle)
{
	minJacob_glb = 1.e5;
	min_pos = -1;
	BadEle.clear();
	const double eps(1.e-10);
	//const double eps(1.e-2);
	for (uint eid = 0; eid < tmesh.size(); eid++)
	{
		double minJacob_ele, GaussPos[3];
		GetHexMinJacob(eid, minJacob_ele, GaussPos);
		if (minJacob_ele < eps)
		{
			BadEle.push_back(eid);
		}
		if (minJacob_ele < minJacob_glb)
		{
			minJacob_glb = minJacob_ele;
			min_pos = eid;
		}
	}
	cout << "minJacob eid nBad: " << minJacob_glb << " " << min_pos << " " << BadEle.size() << "\n";
}

void HexQuality::Optimizing(int nSize, double stepSize)
{
	//OutputCM("../io/hex_input/smooth/fertility_b");
	cout << "\nOptimizing...\n";
	for (int it = 0; it < nSize; it++)
	{
		cout << it << "\n";
		int flag(0);
		double minJacob_glb;
		int min_pos;
		vector<int> BadEle;
		GlobalMinJacob(minJacob_glb, min_pos, BadEle);

		if (BadEle.size() > 0)
		{
			for (uint i = 0; i < BadEle.size(); i++)
			{
				//int tmp = OptimizingElement(BadEle[i], stepSize);
				int tmp = OptimizeElement(BadEle[i], stepSize);
				if (tmp == 1) flag = 1;
			}
		}
		else
		{
			//int tmp = OptimizingElement(min_pos, stepSize);
			int tmp = OptimizeElement(min_pos, stepSize);
			if (tmp == 1) flag = 1;
		}

		double minJacob_glb1;
		int min_pos1;
		vector<int> BadEle1;
		GlobalMinJacob(minJacob_glb1, min_pos1, BadEle1);

		//if (flag == 0 /*|| fabs(minJacob_glb - minJacob_glb1) < 1.e-6*/)
		//{
		//	cout << "After " << it << "steps, no interior point moved!\n";
		//	break;
		//}

		for (uint i = 0; i < tmesh.size(); i++) tmesh[i].trun = 0;
		for (uint i = 0; i < BadEle.size(); i++) tmesh[BadEle[i]].trun = 1;

	}
	cout << "Done optimizing!\n";
	//OutputCM("../io/hex_input/optimize/fertility1");
	//OutputCM("../io/hex_input/optimize/rockerarm1");
	//OutputCM("../io/hex_input/optimize/honda2_2");
	//OutputCM("../io/NAVAIR_GEM/output/optimize/navair3_0");
}

int HexQuality::OptimizeElement(int eid, double stepSize)
{
	double uc[8][3] = { { 0.,0.,0. },{ 1.,0.,0. },{ 1.,1.,0. },{ 0.,1.,0. },
	{ 0.,0.,1. },{ 1.,0.,1. },{ 1.,1.,1. },{ 0.,1.,1. } };
	int cned[8][4] = { { 0,1,3,4 },{ 1,2,0,5 },{ 2,3,1,6 },{ 3,0,2,7 },{ 4,7,5,0 },{ 5,4,6,1 },{ 6,5,7,2 },{ 7,6,4,3 } };

	//find which corner has min Jacob
	double jacob0(1.e5), detJtmp;
	int ploc(-1);
	for (int i = 0; i < 8; i++)
	{
		JacobEval_Scale(eid, uc[i][0], uc[i][1], uc[i][2], detJtmp);
		if (detJtmp < jacob0 && cp[tmesh[eid].cnct[i]].type != 1)
		{
			jacob0 = detJtmp;
			ploc = i;
		}
	}

	int flag(0);
	for (int ic = 0; ic < 4; ic++)
	{
		int i(cned[ploc][ic]);
		if (cp[tmesh[eid].cnct[i]].type != 1)
		{
			//before
			//JacobEval_Scale(eid, uc[i][0], uc[i][1], uc[i][2], jacob0);
			jacob0 = 1.e5;
			int nBad0(0), nBad1(0);
			for (uint j = 0; j < cp[tmesh[eid].cnct[i]].hex.size(); j++)
			{
				int hxid(cp[tmesh[eid].cnct[i]].hex[j]);
				double utmp[3];
				GetHexMinJacob(hxid, detJtmp, utmp);
				if (detJtmp < jacob0)
				{
					jacob0 = detJtmp;
				}
				if (detJtmp < 0.) nBad0++;
			}
			//cout << "jacob before: " << jacob0 << "\n";
			double pold[3] = { cp[tmesh[eid].cnct[i]].coor[0],cp[tmesh[eid].cnct[i]].coor[1],cp[tmesh[eid].cnct[i]].coor[2] };
			double detJ_grad[3];
			JacobEval_Grad(eid, i, detJ_grad);
			//cout << detJ_grad[0] << " " << detJ_grad[1] << " " << detJ_grad[2] << "\n";
			//getchar();
			double edlen = AverageEdgeLength(tmesh[eid].cnct[i]);
			double stepScale(edlen*stepSize);
			cp[tmesh[eid].cnct[i]].coor[0] += stepScale*detJ_grad[0];
			cp[tmesh[eid].cnct[i]].coor[1] += stepScale*detJ_grad[1];
			cp[tmesh[eid].cnct[i]].coor[2] += stepScale*detJ_grad[2];
			//after
			double jacob1(1.e5);
			for (uint j = 0; j < cp[tmesh[eid].cnct[i]].hex.size(); j++)
			{
				int hxid(cp[tmesh[eid].cnct[i]].hex[j]);
				double utmp[3];
				GetHexMinJacob(hxid, detJtmp, utmp);
				if (detJtmp < jacob1)
				{
					jacob1 = detJtmp;
				}
				if (detJtmp < 0.)
				{
					nBad1++;
				}
			}
			//JacobEval_Scale(eid, uc[i][0], uc[i][1], uc[i][2], jacob1);
			//cout << "jacob after: " << jacob1 << "\n";
			//getchar();
			//if (tmesh[eid].cnct[i] == 1765)
			//{
			//	cout << "jacob before: " << jacob0 << "\n";
			//	cout << "jacob after: " << jacob1 << "\n";
			//	getchar();
			//}
			if ((jacob0 < 0. && jacob1 < jacob0) || nBad1 > nBad0)
			{
				cp[tmesh[eid].cnct[i]].coor[0] = pold[0];
				cp[tmesh[eid].cnct[i]].coor[1] = pold[1];
				cp[tmesh[eid].cnct[i]].coor[2] = pold[2];
			}
			else
			{
				flag = 1;
			}
		}
	}
	
	return flag;
}

int HexQuality::OptimizingElement(int eid, double stepSize)
{
	//double stepscale[8];
	//for (int i = 0; i < 8; i++)
	//{
	//	double edlen = 0.;
	//	for (uint j = 0; j < cp[tmesh[eid].cnct[i]].edge.size(); j++)
	//	{
	//		int edid(cp[tmesh[eid].cnct[i]].edge[j]);
	//		int edpid[2] = { tmedge[edid].pt[0],tmedge[edid].pt[1] };
	//		edlen += sqrt((cp[edpid[1]].coor[0] - cp[edpid[0]].coor[0])*(cp[edpid[1]].coor[0] - cp[edpid[0]].coor[0]) +
	//			(cp[edpid[1]].coor[1] - cp[edpid[0]].coor[1])*(cp[edpid[1]].coor[1] - cp[edpid[0]].coor[1]) + 
	//			(cp[edpid[1]].coor[2] - cp[edpid[0]].coor[2])*(cp[edpid[1]].coor[2] - cp[edpid[0]].coor[2]));
	//	}
	//	edlen /= double(cp[tmesh[eid].cnct[i]].edge.size());
	//	stepscale[i] = stepSize*edlen;
	//}

	double uc[8][3] = { { 0.,0.,0. },{ 1.,0.,0. },{ 1.,1.,0. },{ 0.,1.,0. },
	{ 0.,0.,1. },{ 1.,0.,1. },{ 1.,1.,1. },{ 0.,1.,1. } };
	int flag(0);
	for (int i = 0; i < 8; i++)
	{
		if (cp[tmesh[eid].cnct[i]].type != 1)
		{
			//before
			double jacob0(1.e5), jacob1(1.e5), detJtmp;
			//JacobEval_Scale(eid, uc[i][0], uc[i][1], uc[i][2], jacob0);
			for (uint j = 0; j < cp[tmesh[eid].cnct[i]].hex.size(); j++)
			{
				int hxid(cp[tmesh[eid].cnct[i]].hex[j]);
				double utmp[3];
				GetHexMinJacob(hxid, detJtmp, utmp);
				if (detJtmp < jacob0)
				{
					jacob0 = detJtmp;
				}
			}
			//cout << "jacob before: " << jacob0 << "\n";
			double pold[3] = { cp[tmesh[eid].cnct[i]].coor[0],cp[tmesh[eid].cnct[i]].coor[1],cp[tmesh[eid].cnct[i]].coor[2] };
			//after
			double detJ_grad[3];
			JacobEval_Grad(eid, i, detJ_grad);
			//cout << detJ_grad[0] << " " << detJ_grad[1] << " " << detJ_grad[2] << "\n";
			//getchar();
			double edlen = AverageEdgeLength(tmesh[eid].cnct[i]);
			double stepScale(edlen*stepSize);
			cp[tmesh[eid].cnct[i]].coor[0] += stepScale*detJ_grad[0];
			cp[tmesh[eid].cnct[i]].coor[1] += stepScale*detJ_grad[1];
			cp[tmesh[eid].cnct[i]].coor[2] += stepScale*detJ_grad[2];
			for (uint j = 0; j < cp[tmesh[eid].cnct[i]].hex.size(); j++)
			{
				int hxid(cp[tmesh[eid].cnct[i]].hex[j]);
				double utmp[3];
				GetHexMinJacob(hxid, detJtmp, utmp);
				if (detJtmp < jacob1)
				{
					jacob1 = detJtmp;
				}
			}
			//JacobEval_Scale(eid, uc[i][0], uc[i][1], uc[i][2], jacob1);
			//cout << "jacob after: " << jacob1 << "\n";
			//getchar();
			if (jacob1 < jacob0)
			{
				cp[tmesh[eid].cnct[i]].coor[0] = pold[0];
				cp[tmesh[eid].cnct[i]].coor[1] = pold[1];
				cp[tmesh[eid].cnct[i]].coor[2] = pold[2];
			}
			else
			{
				flag = 1;
			}
		}
	}
	return flag;
}

double HexQuality::AverageEdgeLength(int pid)
{
	//int edpt[12][2] = { {0,1},{1,2},{2,3},{3,0},{0,4},{1,5},{2,6},{3,7},{4,5},{5,6},{6,7},{7,4} };
	double len(0.), dtmp;
	for (int i = 0; i < cp[pid].edge.size(); i++)
	{
		int edid(cp[pid].edge[i]);
		int edpt[2] = { tmedge[edid].pt[0],tmedge[edid].pt[1] };
		double vtmp[3] = { cp[edpt[1]].coor[0] - cp[edpt[0]].coor[0],cp[edpt[1]].coor[1] - cp[edpt[0]].coor[1],
			cp[edpt[1]].coor[2] - cp[edpt[0]].coor[2] };
		dtmp = sqrt(vtmp[0] * vtmp[0] + vtmp[1] * vtmp[1] + vtmp[2] * vtmp[2]);
		len += dtmp;
	}
	len /= double(cp[pid].edge.size());
	return len;
}

void HexQuality::JacobEval_Grad(int eid, int ploc, double grad[3])
{
	int cned[8][3] = { {1,3,4},{2,0,5},{3,1,6},{0,2,7},{7,5,0},{4,6,1},{5,7,2},{6,4,3} };
	//double a[3][3] = { {cp[tmesh[eid].cnct[cned[ploc][0]]].coor[0]-cp[tmesh[eid].cnct[ploc]].coor[0],
	//	cp[tmesh[eid].cnct[cned[ploc][1]]].coor[0] - cp[tmesh[eid].cnct[ploc]].coor[0], 
	//	cp[tmesh[eid].cnct[cned[ploc][2]]].coor[0] - cp[tmesh[eid].cnct[ploc]].coor[0]},
	//	{ cp[tmesh[eid].cnct[cned[ploc][0]]].coor[1] - cp[tmesh[eid].cnct[ploc]].coor[1],
	//	cp[tmesh[eid].cnct[cned[ploc][1]]].coor[1] - cp[tmesh[eid].cnct[ploc]].coor[1],
	//	cp[tmesh[eid].cnct[cned[ploc][2]]].coor[1] - cp[tmesh[eid].cnct[ploc]].coor[1] },
	//	{ cp[tmesh[eid].cnct[cned[ploc][0]]].coor[2] - cp[tmesh[eid].cnct[ploc]].coor[2],
	//	cp[tmesh[eid].cnct[cned[ploc][1]]].coor[2] - cp[tmesh[eid].cnct[ploc]].coor[2],
	//	cp[tmesh[eid].cnct[cned[ploc][2]]].coor[2] - cp[tmesh[eid].cnct[ploc]].coor[2] } };
	////double a[3][3] = { {-3,2,-5},{-1,0,-2},{3,-4,1} };
	//double adj[3][3] = { {a[1][1] * a[2][2] - a[1][2] * a[2][1],a[0][2] * a[2][1] - a[0][1] * a[2][2], a[0][1] * a[1][2] - a[0][2] * a[1][1] },
	//{ a[1][2] * a[2][0] - a[1][0] * a[2][2],a[0][0] * a[2][2] - a[0][2] * a[2][0], a[0][2] * a[1][0] - a[0][0] * a[1][2] }, 
	//{ a[1][0] * a[2][1] - a[1][1] * a[2][0],a[0][1] * a[2][0] - a[0][0] * a[2][1], a[0][0] * a[1][1] - a[0][1] * a[1][0] } };
	//Jacob_Grad[0] = -adj[0][0] - adj[1][0] - adj[2][0];
	//Jacob_Grad[1] = -adj[0][1] - adj[1][1] - adj[2][1];
	//Jacob_Grad[2] = -adj[0][2] - adj[1][2] - adj[2][2];
	//double nm = sqrt(Jacob_Grad[0] * Jacob_Grad[0] + Jacob_Grad[1] * Jacob_Grad[1] + Jacob_Grad[2] * Jacob_Grad[2]);
	//if (nm > 1.e-8)
	//{
	//	Jacob_Grad[0] /= nm;
	//	Jacob_Grad[1] /= nm;
	//	Jacob_Grad[2] /= nm;
	//}
	////cout << adj[0][0] << "\t" << adj[0][1] << "\t" << adj[0][2] << "\n";
	////cout << adj[1][0] << "\t" << adj[1][1] << "\t" << adj[1][2] << "\n";
	////cout << adj[2][0] << "\t" << adj[2][1] << "\t" << adj[2][2] << "\n";
	////getchar();

	double XSIETA[8][3] = { {0.,0.,0.},{ 1.,0.,0. },{ 1.,1.,0. },{ 0.,1.,0. },
	{ 0.,0.,1. },{ 1.,0.,1. },{ 1.,1.,1. },{ 0.,1.,1. }, };
	double jacob;
	double g00, g01, g02, g10, g11, g12, g20, g21, g22;
	double phi[8], phi_xsi[8], phi_eta[8], phi_zeta[8];

	double xsi, eta, zeta;
	//bool flag = true;
	//for (i = 0; i < 8; i++)
	//{
		xsi = XSIETA[ploc][0];
		eta = XSIETA[ploc][1];
		zeta = XSIETA[ploc][2];

		phi[0] = (1 - xsi)*(1 - eta)*(1 - zeta); 	phi[1] = xsi*(1 - eta)*(1 - zeta);
		phi[2] = xsi*eta*(1 - zeta);			phi[3] = (1 - xsi)*eta*(1 - zeta);
		phi[4] = (1 - xsi)*(1 - eta)*zeta;		phi[5] = xsi*(1 - eta)*zeta;
		phi[6] = xsi*eta*zeta;				phi[7] = (1 - xsi)*eta*zeta;

		phi_xsi[0] = -(1 - eta)*(1 - zeta);		phi_xsi[1] = (1 - eta)*(1 - zeta);
		phi_xsi[2] = eta*(1 - zeta);			phi_xsi[3] = -eta*(1 - zeta);
		phi_xsi[4] = -(1 - eta)*zeta;			phi_xsi[5] = (1 - eta)*zeta;
		phi_xsi[6] = eta*zeta;				phi_xsi[7] = -eta*zeta;

		phi_eta[0] = -(1 - xsi)*(1 - zeta);		phi_eta[1] = -xsi*(1 - zeta);
		phi_eta[2] = xsi*(1 - zeta);			phi_eta[3] = (1 - xsi)*(1 - zeta);
		phi_eta[4] = -(1 - xsi)*zeta;			phi_eta[5] = -xsi*zeta;
		phi_eta[6] = xsi*zeta;				phi_eta[7] = (1 - xsi)*zeta;

		phi_zeta[0] = -(1 - xsi)*(1 - eta);		phi_zeta[1] = -xsi*(1 - eta);
		phi_zeta[2] = -xsi*eta;				phi_zeta[3] = -(1 - xsi)*eta;
		phi_zeta[4] = (1 - xsi)*(1 - eta);		phi_zeta[5] = xsi*(1 - eta);
		phi_zeta[6] = xsi*eta;				phi_zeta[7] = (1 - xsi)*eta;

		g00 = 0.0;	g01 = 0.0;	g02 = 0.0;	g10 = 0.0;	g11 = 0.0;	g12 = 0.0;
		g20 = 0.0; g21 = 0.0; g22 = 0.0;
		for (int j = 0; j < 8; j++)
		{
			g00 += cp[tmesh[eid].cnct[j]].coor[0] * phi_xsi[j];
			g01 += cp[tmesh[eid].cnct[j]].coor[0] * phi_eta[j];
			g02 += cp[tmesh[eid].cnct[j]].coor[0] * phi_zeta[j];
			g10 += cp[tmesh[eid].cnct[j]].coor[1] * phi_xsi[j];
			g11 += cp[tmesh[eid].cnct[j]].coor[1] * phi_eta[j];
			g12 += cp[tmesh[eid].cnct[j]].coor[1] * phi_zeta[j];
			g20 += cp[tmesh[eid].cnct[j]].coor[2] * phi_xsi[j];
			g21 += cp[tmesh[eid].cnct[j]].coor[2] * phi_eta[j];
			g22 += cp[tmesh[eid].cnct[j]].coor[2] * phi_zeta[j];
			//cout << phi_zeta[j] << "\n";
			//cout << cp[tmesh[eid].cnct[j]].coor[0] << " " << cp[tmesh[eid].cnct[j]].coor[1] << " " <<
			//	cp[tmesh[eid].cnct[j]].coor[2] << "\n";
		}
		//cout << "\n";
		//cout << g00 << "\t" << g01 << "\t" << g02 << "\n" << g10 << "\t" << g11 << "\t" << g12 << "\n" << g20 << "\t" << g21 << "\t" << g22 << "\n";
		//jacob = g00*g11*g22 + g10*g21*g02 + g20*g01*g12 - g02*g20*g11 - g01*g10*g22 - g00*g12*g21;
		//cout << "jacob: " << jacob << "\n";
		//getchar();

		// normalize
		double dtemp = 0., vtmp[3] = { 0.,0.,0. };
		//for (int j = 0; j < 3; j++)
		//{
		//	vtmp[0] = cp[tmesh[eid].cnct[cned[ploc][j]]].coor[0] - cp[tmesh[eid].cnct[ploc]].coor[0];
		//	vtmp[1] = cp[tmesh[eid].cnct[cned[ploc][j]]].coor[1] - cp[tmesh[eid].cnct[ploc]].coor[1];
		//	vtmp[2] = cp[tmesh[eid].cnct[cned[ploc][j]]].coor[2] - cp[tmesh[eid].cnct[ploc]].coor[2];
		//	dtemp += sqrt(vtmp[0] * vtmp[0] + vtmp[1] * vtmp[1] + vtmp[2] * vtmp[2]);
		//}
		//dtemp /= 3.;
		//if (dtemp > 1.e-6)
		//	jacob /= pow(dtemp, 3);

		double tempgrad[3];
		//if (flag || jacob < minJacob)
		{
			tempgrad[0] = phi_xsi[ploc] * g11*g22 + g10*g21*phi_zeta[ploc] + g20*phi_eta[ploc] * g12
				- g20*phi_zeta[ploc] * g11 - phi_eta[ploc] * g10*g22 - phi_xsi[ploc] * g12*g21;

			tempgrad[1] = g00*phi_eta[ploc] * g22 + phi_xsi[ploc] * g21*g02 + g20*g01*phi_zeta[ploc]
				- g02*g20*phi_eta[ploc] - g01*phi_xsi[ploc] * g22 - g00*phi_zeta[ploc] * g21;

			tempgrad[2] = g00*g11*phi_zeta[ploc] + g10*phi_eta[ploc] * g02 + phi_xsi[ploc] * g01*g12
				- g02*phi_xsi[ploc] * g11 - g01*g10*phi_zeta[ploc] - g00*g12*phi_eta[ploc];

			dtemp = tempgrad[0] * tempgrad[0] + tempgrad[1] * tempgrad[1] + tempgrad[2] * tempgrad[2];
			//if (fabs(dtemp) > 1.e-6)
			{
				//minJacob = jacob[i];
				grad[0] = tempgrad[0];
				grad[1] = tempgrad[1];
				grad[2] = tempgrad[2];
				//flag = false;
			}
		}
	//}

	dtemp = sqrt(grad[0] * grad[0] + grad[1] * grad[1] + grad[2] * grad[2]);
	if (dtemp > 1.e-12)
	{
		for (int i = 0; i < 3; i++)
			grad[i] /= dtemp;
	}
}

//void HexQuality::GlobalMinJacobBEXT(double& minJacob_glb, int& min_pos, vector<int>& BadEle)
//{
//	minJacob_glb = 1.e5;
//	min_pos = -1;
//	BadEle.clear();
//	const double eps(1.e-10);
//	for (uint eid = 0; eid < tmesh.size(); eid++)
//	{
//		double minJacob_ele, GaussPos[3];
//		GetHexMinJacobBEXT(eid, minJacob_ele, GaussPos);
//		if (minJacob_ele < eps)
//		{
//			BadEle.push_back(eid);
//		}
//		if (minJacob_ele < minJacob_glb)
//		{
//			minJacob_glb = minJacob_ele;
//			min_pos = eid;
//		}
//	}
//	cout << "minJacob eid nBad: " << minJacob_glb << " " << min_pos << " " << BadEle.size() << "\n";
//}
//
//void HexQuality::GetHexMinJacobBEXT(int eid, double& minJacob_ele, double GaussPos[3])
//{
//	vector<double> Gpt, wght;
//	GetGaussPoint(2, Gpt, wght);
//	Gpt[0] = 0.; Gpt[1] = 1.;
//	double detJ;
//	minJacob_ele = 1.e5;
//	GaussPos[0] = 0.; GaussPos[1] = 0.; GaussPos[2] = 0.;
//	for (uint i = 0; i < Gpt.size(); i++)
//	{
//		for (uint j = 0; j < Gpt.size(); j++)
//		{
//			for (uint k = 0; k < Gpt.size(); k++)
//			{
//				//JacobEval(eid, Gpt[i], Gpt[j], Gpt[k], detJ);
//				JacobEval_ScaleBEXT(eid, Gpt[i], Gpt[j], Gpt[k], detJ);
//				if (detJ < minJacob_ele)
//				{
//					minJacob_ele = detJ;
//					GaussPos[0] = Gpt[i];
//					GaussPos[1] = Gpt[j];
//					GaussPos[2] = Gpt[k];
//				}
//			}
//		}
//	}
//}

//void HexQuality::JacobEval_ScaleBEXT(int eid, double u, double v, double w, double& detJ)
//{
//	double Nu[4] = { (1. - u)*(1. - u)*(1. - u), 3.*(1. - u)*(1. - u)*u, 3.*(1. - u)*u*u, u*u*u };
//	double Nv[4] = { (1. - v)*(1. - v)*(1. - v), 3.*(1. - v)*(1. - v)*v, 3.*(1. - v)*v*v, v*v*v };
//	double Nw[4] = { (1. - w)*(1. - w)*(1. - w), 3.*(1. - w)*(1. - w)*w, 3.*(1. - w)*w*w, w*w*w };
//
//	double dNdu[4] = { -3.*(1. - u)*(1. - u), 3. - 12.*u + 9.*u*u, 3.*(2. - 3.*u)*u, 3.*u*u };
//	double dNdv[4] = { -3.*(1. - v)*(1. - v), 3. - 12.*v + 9.*v*v, 3.*(2. - 3.*v)*v, 3.*v*v };
//	double dNdw[4] = { -3.*(1. - w)*(1. - w), 3. - 12.*w + 9.*w*w, 3.*(2. - 3.*w)*w, 3.*w*w };
//
//	double dNdt[64][3];
//	double Nx_bz[64];
//	double dNdx_bz[64][3];
//
//	vector<double> Nx;
//	vector<array<double, 3>> dNdx;
//
//	Nx.resize(bzmesh[eid].IEN.size());
//	dNdx.resize(bzmesh[eid].IEN.size());
//
//	int i, j, k, a, b, c, loc;
//	loc = 0;
//	for (i = 0; i<4; i++) {
//		for (j = 0; j<4; j++) {
//			for (k = 0; k < 4; k++) {
//				Nx_bz[loc] = Nu[k] * Nv[j] * Nw[i];
//				dNdt[loc][0] = dNdu[k] * Nv[j] * Nw[i];
//				dNdt[loc][1] = Nu[k] * dNdv[j] * Nw[i];
//				dNdt[loc][2] = Nu[k] * Nv[j] * dNdw[i];
//				loc++;
//			}
//		}
//	}
//
//	for (i = 0; i < bzmesh[eid].IEN.size; i++)
//	{
//		for (j = 0; j < 64; j++)
//		{
//			Nx[i] += bzmesh[eid].cmat[i][j] * Nx_bz[j];
//				for (k = 0; k < 3; k++)
//					dNdx[i] += bzmesh[eid].cmat[i][j] * dNdx_bz[j][k];
//		}
//	}
//
//
//
//	double dxdt[3][3] = { { 0.,0.,0. },{ 0.,0.,0. },{ 0.,0.,0. } };
//	for (i = 0; i < bzmesh[eid].IEN.size; i++)
//		for (a = 0; a<3; a++)
//			for (b = 0; b<3; b++)
//				dxdt[a][b] += cp[bzmesh[eid].IEN[i]].coor[a] * dNdt[i][b];
//
//
//	for (int i = 0; i < 3; i++)
//	{
//		for (int j = 0; j < 3; j++)
//		{
//			for (int k = 0; k < 64; k++)
//			{
//				dxdt[i][j] += cp[tmesh[eid].cnct[k]].coor[i] * dNdt[k][j];
//			}
//		}
//	}
//	detJ = dxdt[0][0] * (dxdt[1][1] * dxdt[2][2] - dxdt[1][2] * dxdt[2][1]) -
//		dxdt[0][1] * (dxdt[1][0] * dxdt[2][2] - dxdt[1][2] * dxdt[2][0]) +
//		dxdt[0][2] * (dxdt[1][0] * dxdt[2][1] - dxdt[1][1] * dxdt[2][0]);
//	double nm[3] = { sqrt(dxdt[0][0] * dxdt[0][0] + dxdt[1][0] * dxdt[1][0] + dxdt[2][0] * dxdt[2][0]),
//		sqrt(dxdt[0][1] * dxdt[0][1] + dxdt[1][1] * dxdt[1][1] + dxdt[2][1] * dxdt[2][1]),
//		sqrt(dxdt[0][2] * dxdt[0][2] + dxdt[1][2] * dxdt[1][2] + dxdt[2][2] * dxdt[2][2]) };
//	double tmp(nm[0] * nm[1] * nm[2]);
//	if (tmp > 1.e-12)
//	{
//		detJ /= tmp;
//	}
//}

//////////////////////////////////////////////////////////////////////////////////////

void HexQuality::Smoothing_Adapt(int nSize, double stepSize)
{
	cout << "\nSmoothing...\n";

	vector<array<int, 2>> psm;
	for (uint i = 0; i < hcp.size(); i++)
	{
		for (uint j = 0; j < hcp[i].size(); j++)
		{
			hcp[i][j].smth = 0;
		}
	}
	//for (uint i = 0; i < hmesh.size(); i++)
	//{
	//	for (uint j = 0; j < hmesh[i].size(); j++)
	//	{
	//		//if (hmesh[i][j].act == 1)
	//		if (hmesh[i][j].act == 1 && i == hmesh.size() - 1)
	//		{
	//			for (int k = 0; k < 8; k++)
	//			{
	//				hcp[i][hmesh[i][j].cnct[k]].smth = 1;
	//			}
	//		}
	//	}
	//}
	//for (uint i = 0; i < hcp.size(); i++)
	//{
	//	for (uint j = 0; j < hcp[i].size(); j++)
	//	{
	//		if (hcp[i][j].smth == 1)
	//		{
	//			array<int, 2> itmp = { i,j };
	//			psm.push_back(itmp);
	//		}
	//	}
	//}

	int il(hcp.size() - 1);
	vector<int> flag(hcp[il].size(), 1);
	for (uint j = 0; j < hface[il].size(); j++)
	{
		if (hface[il][j].hex.size() == 1)
		{
			for (int k = 0; k < 4; k++)
			{
				flag[hface[il][j].cnct[k]] = 0;
			}
		}
	}

	vector<array<int, 2>> psmb;
	for (uint j = 0; j < hcp[il].size(); j++)
	{
		if (flag[j] == 1)
		{
			array<int, 2> itmp = { il,j };
			psm.push_back(itmp);
		}
		else
		{
			array<int, 2> itmp = { il,j };
			psmb.push_back(itmp);
		}
	}
	

	cout << "\n# points to be smoothed: "<<psm.size()<<"\n";

	for (int it = 0; it < nSize; it++)
	{
		cout << "it: " << it << "\n";
		int flag(0);

		//double minJacob_glb0;
		//array<int,2> min_pos0;
		//vector<array<int, 2>> BadEle0;
		//GlobalMinJacob_Adapt(minJacob_glb0, min_pos0, BadEle0);

		//for (uint i0 = 0; i0 < BadEle0.size(); i0++)
		//{
		//	//cout << BadEle0[i0] << " ";
		//	tmesh[BadEle0[i0]].trun = 1;
		//	//int eid(BadEle[i0]);
		//	//for (int i = 0; i < 8; i++)
		//	//{
		//	//	if (cp[tmesh[eid].cnct[i]].type != 1)
		//	//	{
		//	//		int tmp = SmoothingPoint(tmesh[eid].cnct[i], stepSize);
		//	//		if (tmp == 1) flag = 1;
		//	//	}
		//	//}
		//}
		//break;

		for (uint i = 0; i < psm.size(); i++)
		{
			if (hcp[psm[i][0]][psm[i][1]].type != 1)//interior points only
			{
				int tmp = SmoothingPoint_Adapt(psm[i][0], psm[i][1], stepSize);
				if (tmp == 1) flag = 1;
			}
			//else 
			//{
			//	int tmp = SmoothingPointBoundary(i, stepSize);
			//	if (tmp == 1) flag = 1;
			//}
		}
		for (uint i = 0; i < psmb.size(); i++)
		{
			{
				int tmp = SmoothingPointBoundary_Adapt(psmb[i][0], psmb[i][1], stepSize);
			}
			//else 
			//{
			//	int tmp = SmoothingPointBoundary(i, stepSize);
			//	if (tmp == 1) flag = 1;
			//}
		}

		double minJacob_glb;
		array<int,2> min_pos;
		vector<array<int,2>> BadEle;
		GlobalMinJacob_Adapt(minJacob_glb, min_pos, BadEle);

		//if (flag == 0 || minJacob_glb0 > minJacob_glb)
		//{
		//	cout << "After " << it << "steps, no interior point moved!\n";
		//	break;
		//}

	}
	cout << "Done smoothing!\n";

	//OutputCM("../io/hex_input/smooth/fertility2");
	//OutputCM("../io/hex_input/smooth/cube_coarse");
	//OutputCM("../io/hex_input/smooth/rockerarm2");
	//OutputCM("../io/hex_input/smooth/navair_coarse");
	//OutputCM("../io/hex_input/smooth/honda1");
	//OutputCM("../io/hex_input/smooth/honda2");
	//OutputCM("../io/hex_input/smooth/honda1_dense");
	//OutputCM("../io/hex_input/smooth/honda1m");
	//OutputCM("../io/hex_input/smooth/honda2_2");
	//OutputCM("../io/hex_input/smooth/heli_dense_loc1");
}

int HexQuality::SmoothingPoint_Adapt(int lev, int pid, double stepSize)
{
	const double eps(1.e-8);
	double vol_all(0.), center[3] = { 0.,0.,0., };
	for (uint i = 0; i < hcp[lev][pid].hex.size(); i++)
	{
		double e_vol, e_center[3];
		GetHexVolAndCenter_Adapt(lev,hcp[lev][pid].hex[i], e_vol, e_center);
		vol_all += e_vol;
		center[0] += e_center[0] * e_vol;
		center[1] += e_center[1] * e_vol;
		center[2] += e_center[2] * e_vol;
	}
	if (vol_all < eps)
	{
		//cerr << "Element with almost zero volume!\n";
		return 0;
	}
	center[0] /= vol_all;
	center[1] /= vol_all;
	center[2] /= vol_all;
	double dir[3] = { center[0] - hcp[lev][pid].coor[0],center[1] - hcp[lev][pid].coor[1], center[2] - hcp[lev][pid].coor[2] };
	double pold[3] = { hcp[lev][pid].coor[0],hcp[lev][pid].coor[1],hcp[lev][pid].coor[2] };

	//cout << dir[0] << " " << dir[1] << " " << dir[2] << "\n";
	//getchar();

	//before smoothing
	double minJacob_0(1.e5), minJacob_1(1.e5), tmpJacob, GaussPos[3];
	for (uint i = 0; i < hcp[lev][pid].hex.size(); i++)
	{
		GetHexMinJacob_Adapt(lev, hcp[lev][pid].hex[i], tmpJacob, GaussPos);
		if (tmpJacob < minJacob_0) minJacob_0 = tmpJacob;
	}
	//double edlen = AverageEdgeLength(pid);
	//double stepSize_scale(stepSize*edlen);
	double stepSize_scale(stepSize);
	//after smoothing
	hcp[lev][pid].coor[0] = hcp[lev][pid].coor[0] + stepSize_scale*dir[0];
	hcp[lev][pid].coor[1] = hcp[lev][pid].coor[1] + stepSize_scale*dir[1];
	hcp[lev][pid].coor[2] = hcp[lev][pid].coor[2] + stepSize_scale*dir[2];
	for (uint i = 0; i < hcp[lev][pid].hex.size(); i++)
	{
		GetHexMinJacob_Adapt(lev, hcp[lev][pid].hex[i], tmpJacob, GaussPos);
		if (tmpJacob < minJacob_1) minJacob_1 = tmpJacob;
	}
	if (minJacob_1 < minJacob_0)
	{
		hcp[lev][pid].coor[0] = pold[0];
		hcp[lev][pid].coor[1] = pold[1];
		hcp[lev][pid].coor[2] = pold[2];
		return 0;
	}
	return 1;
}

int HexQuality::SmoothingPointBoundary_Adapt(int lev, int pid, double stepSize)
{
	const double eps(1.e-8);
	double area_all(0.), center[3] = { 0.,0.,0. }, nm[3] = { 0., 0., 0. };
	int nbf(0);
	for (uint i = 0; i < hcp[lev][pid].face.size(); i++)
	{
		if (hface[lev][hcp[lev][pid].face[i]].hex.size() == 1)
		{
			double e_area, e_center[3], e_nm[3];
			int* it = find(hface[lev][hcp[lev][pid].face[i]].cnct, hface[lev][hcp[lev][pid].face[i]].cnct + 4, pid);
			int loc(it - hface[lev][hcp[lev][pid].face[i]].cnct);
			GetQuadInfo_Adapt(lev, hcp[lev][pid].face[i], loc, e_area, e_center, e_nm);
			area_all += e_area;
			center[0] += e_center[0] * e_area;
			center[1] += e_center[1] * e_area;
			center[2] += e_center[2] * e_area;
			nm[0] += e_nm[0] * e_area;
			nm[1] += e_nm[1] * e_area;
			nm[2] += e_nm[2] * e_area;
			//nm[0] += e_nm[0];
			//nm[1] += e_nm[1];
			//nm[2] += e_nm[2];
			nbf++;
		}
	}
	if (area_all < eps)
	{
		cerr << "Element with almost zero area!\n";
		return 0;
	}
	center[0] /= area_all;	center[1] /= area_all;	center[2] /= area_all;
	nm[0] /= area_all;	nm[1] /= area_all;	nm[2] /= area_all;
	//nm[0] /= double(nbf);	nm[1] /= double(nbf);	nm[2] /= double(nbf);
	double dst = sqrt(nm[0] * nm[0] + nm[1] * nm[1] + nm[2] * nm[2]);
	nm[0] /= dst; nm[1] /= dst; nm[2] /= dst;
	double dir[3] = { center[0] - hcp[lev][pid].coor[0],center[1] - hcp[lev][pid].coor[1], center[2] - hcp[lev][pid].coor[2] };
	dst = dir[0] * nm[0] + dir[1] * nm[1] + dir[2] * nm[2];
	dir[0] = dir[0] - dst*nm[0];
	dir[1] = dir[1] - dst*nm[1];
	dir[2] = dir[2] - dst*nm[2];

	//dst = dir[0] * nm[0] + dir[1] * nm[1] + dir[2] * nm[2];
	//cout << "inner product: " << dst << "\n";
	//getchar();

	double pold[3] = { hcp[lev][pid].coor[0],hcp[lev][pid].coor[1],hcp[lev][pid].coor[2] };

	//before smoothing
	double minJacob_0(1.e5), minJacob_1(1.e5), tmpJacob, GaussPos[3];
	for (uint i = 0; i < hcp[lev][pid].hex.size(); i++)
	{
		GetHexMinJacob_Adapt(lev, hcp[lev][pid].hex[i], tmpJacob, GaussPos);
		if (tmpJacob < minJacob_0) minJacob_0 = tmpJacob;
	}
	//double edlen = AverageEdgeLength(pid);
	//double stepSize_scale(stepSize*edlen);
	double stepSize_scale(stepSize);
	//after smoothing
	hcp[lev][pid].coor[0] = hcp[lev][pid].coor[0] + stepSize_scale*dir[0];
	hcp[lev][pid].coor[1] = hcp[lev][pid].coor[1] + stepSize_scale*dir[1];
	hcp[lev][pid].coor[2] = hcp[lev][pid].coor[2] + stepSize_scale*dir[2];
	for (uint i = 0; i < hcp[lev][pid].hex.size(); i++)
	{
		GetHexMinJacob_Adapt(lev, hcp[lev][pid].hex[i], tmpJacob, GaussPos);
		if (tmpJacob < minJacob_1) minJacob_1 = tmpJacob;
	}
	if (minJacob_1 < minJacob_0)
	{
		hcp[lev][pid].coor[0] = pold[0];
		hcp[lev][pid].coor[1] = pold[1];
		hcp[lev][pid].coor[2] = pold[2];
		return 0;
	}
	return 1;
}

void HexQuality::GetHexMinJacob_Adapt(int lev, int eid, double& minJacob_ele, double GaussPos[3])
{
	vector<double> Gpt(2);
	Gpt[0] = 0.; Gpt[1] = 1.;
	double detJ[2];
	minJacob_ele = 1.e5;
	GaussPos[0] = 0.; GaussPos[1] = 0.; GaussPos[2] = 0.;
	for (uint i = 0; i < Gpt.size(); i++)
	{
		for (uint j = 0; j < Gpt.size(); j++)
		{
			for (uint k = 0; k < Gpt.size(); k++)
			{
				JacobEval_Scale_Adapt(lev, eid, Gpt[i], Gpt[j], Gpt[k], detJ);
				if (detJ[0] < minJacob_ele)
				{
					minJacob_ele = detJ[0];
					GaussPos[0] = Gpt[i];
					GaussPos[1] = Gpt[j];
					GaussPos[2] = Gpt[k];
				}
			}
		}
	}
}

void HexQuality::GetHexVolAndCenter_Adapt(int lev, int eid, double& vol, double center[3])
{
	vol = 0.;
	center[0] = 0.; center[1] = 0.; center[2] = 0.;
	for (int i = 0; i < 8; i++)
	{
		center[0] += hcp[lev][hmesh[lev][eid].cnct[i]].coor[0];
		center[1] += hcp[lev][hmesh[lev][eid].cnct[i]].coor[1];
		center[2] += hcp[lev][hmesh[lev][eid].cnct[i]].coor[2];
	}
	center[0] /= 8.; center[1] /= 8.; center[2] /= 8.;

	vector<double> Gpt, wght;
	GetGaussPoint(2, Gpt, wght);
	double detJ[2];
	for (uint i = 0; i < Gpt.size(); i++)
	{
		for (uint j = 0; j < Gpt.size(); j++)
		{
			for (uint k = 0; k < Gpt.size(); k++)
			{
				JacobEval_Scale_Adapt(lev, eid, Gpt[i], Gpt[j], Gpt[k], detJ);
				vol += wght[i] * wght[j] * wght[k] * detJ[1];
			}
		}
	}
}

void HexQuality::GetQuadInfo_Adapt(int lev, int fcid, int ploc, double& area, double center[3], double nm[3])
{
	//face is counter-clock-wise oriented
	area = 0.;
	center[0] = 0.; center[1] = 0.; center[2] = 0.;
	nm[0] = 0.; nm[1] = 0.; nm[2] = 0.;
	for (int i = 0; i < 4; i++)
	{
		center[0] += hcp[lev][hface[lev][fcid].cnct[i]].coor[0];
		center[1] += hcp[lev][hface[lev][fcid].cnct[i]].coor[1];
		center[2] += hcp[lev][hface[lev][fcid].cnct[i]].coor[2];
	}
	center[0] /= 4.; center[1] /= 4.; center[2] /= 4.;

	int tp[2][3] = { { hface[lev][fcid].cnct[ploc],hface[lev][fcid].cnct[(ploc + 1) % 4],hface[lev][fcid].cnct[(ploc + 3) % 4] },
	{ hface[lev][fcid].cnct[(ploc + 2) % 4],hface[lev][fcid].cnct[(ploc + 3) % 4],hface[lev][fcid].cnct[(ploc + 1) % 4] } };
	double v0[2][3] = { { hcp[lev][tp[0][1]].coor[0] - hcp[lev][tp[0][0]].coor[0],hcp[lev][tp[0][1]].coor[1] - hcp[lev][tp[0][0]].coor[1],hcp[lev][tp[0][1]].coor[2] - hcp[lev][tp[0][0]].coor[2] },
	{ hcp[lev][tp[0][2]].coor[0] - hcp[lev][tp[0][0]].coor[0],hcp[lev][tp[0][2]].coor[1] - hcp[lev][tp[0][0]].coor[1],hcp[lev][tp[0][2]].coor[2] - hcp[lev][tp[0][0]].coor[2] } };
	double v1[2][3] = { { hcp[lev][tp[1][1]].coor[0] - hcp[lev][tp[1][0]].coor[0],hcp[lev][tp[1][1]].coor[1] - hcp[lev][tp[1][0]].coor[1],hcp[lev][tp[1][1]].coor[2] - hcp[lev][tp[1][0]].coor[2] },
	{ hcp[lev][tp[1][2]].coor[0] - hcp[lev][tp[1][0]].coor[0],hcp[lev][tp[1][2]].coor[1] - hcp[lev][tp[1][0]].coor[1],hcp[lev][tp[1][2]].coor[2] - hcp[lev][tp[1][0]].coor[2] } };
	double vc0[3] = { v0[0][1] * v0[1][2] - v0[0][2] * v0[1][1],v0[0][2] * v0[1][0] - v0[0][0] * v0[1][2],v0[0][0] * v0[1][1] - v0[0][1] * v0[1][0] };
	double vc1[3] = { v1[0][1] * v1[1][2] - v1[0][2] * v1[1][1],v1[0][2] * v1[1][0] - v1[0][0] * v1[1][2],v1[0][0] * v1[1][1] - v1[0][1] * v1[1][0] };
	double dst0 = sqrt(vc0[0] * vc0[0] + vc0[1] * vc0[1] + vc0[2] * vc0[2]);
	double dst1 = sqrt(vc1[0] * vc1[0] + vc1[1] * vc1[1] + vc1[2] * vc1[2]);

	area = (dst0 + dst1) / 2.;
	nm[0] = vc0[0] / dst0; nm[1] = vc0[1] / dst0; nm[2] = vc0[2] / dst0;
}

void HexQuality::JacobEval_Scale_Adapt(int lev, int eid, double u, double v, double w, double detJ[2])
{
	double Nu[2] = { 1. - u,u };
	double Nv[2] = { 1. - v,v };
	double Nw[2] = { 1. - w,w };
	double dNdu[2] = { -1.,1. };
	double dNdv[2] = { -1.,1. };
	double dNdw[2] = { -1.,1. };
	double dNdt[8][3];
	int ploc[8] = { 0,1,3,2,4,5,7,6 };
	int loc(0);
	for (int i = 0; i < 2; i++)
	{
		for (int j = 0; j < 2; j++)
		{
			for (int k = 0; k < 2; k++)
			{
				dNdt[ploc[loc]][0] = dNdu[k] * Nv[j] * Nw[i];
				dNdt[ploc[loc]][1] = Nu[k] * dNdv[j] * Nw[i];
				dNdt[ploc[loc]][2] = Nu[k] * Nv[j] * dNdw[i];
				loc++;
			}
		}
	}
	double dxdt[3][3] = { { 0.,0.,0. },{ 0.,0.,0. },{ 0.,0.,0. } };
	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			for (int k = 0; k < 8; k++)
			{
				dxdt[i][j] += hcp[lev][hmesh[lev][eid].cnct[k]].coor[i] * dNdt[k][j];
			}
		}
	}
	detJ[0] = dxdt[0][0] * (dxdt[1][1] * dxdt[2][2] - dxdt[1][2] * dxdt[2][1]) -
		dxdt[0][1] * (dxdt[1][0] * dxdt[2][2] - dxdt[1][2] * dxdt[2][0]) +
		dxdt[0][2] * (dxdt[1][0] * dxdt[2][1] - dxdt[1][1] * dxdt[2][0]);
	detJ[1] = detJ[0];
	double nm[3] = { sqrt(dxdt[0][0] * dxdt[0][0] + dxdt[1][0] * dxdt[1][0] + dxdt[2][0] * dxdt[2][0]),
		sqrt(dxdt[0][1] * dxdt[0][1] + dxdt[1][1] * dxdt[1][1] + dxdt[2][1] * dxdt[2][1]),
		sqrt(dxdt[0][2] * dxdt[0][2] + dxdt[1][2] * dxdt[1][2] + dxdt[2][2] * dxdt[2][2]) };
	double tmp(nm[0] * nm[1] * nm[2]);
	if (tmp > 1.e-12)
	{
		detJ[0] /= tmp;
	}
}

void HexQuality::GlobalMinJacob_Adapt(double& minJacob_glb, array<int, 2>& min_pos, vector<array<int,2>>& BadEle)
{
	minJacob_glb = 1.e5;
	min_pos[0] = -1; min_pos[0] = -1;
	BadEle.clear();
	const double eps(1.e-10);
	for (uint i = 0; i < hmesh.size(); i++)
	{
		for (uint j = 0; j < hmesh[i].size(); j++)
		{
			if (hmesh[i][j].act == 1)
			{
				double minJacob_ele, GaussPos[3];
				GetHexMinJacob_Adapt(i, j, minJacob_ele, GaussPos);
				if (minJacob_ele < eps)
				{
					array<int, 2> itmp = { i,j };
					BadEle.push_back(itmp);
				}
				if (minJacob_ele < minJacob_glb)
				{
					minJacob_glb = minJacob_ele;
					min_pos[0] = i; min_pos[1] = j;
				}
			}
		}
	}
	cout << "minJacob eid nBad: " << minJacob_glb << " " << BadEle.size() << "\n";
}

void HexQuality::BadElementFlag_Adapt(const vector<array<int, 2>>& rfid)
{
	for (uint i = 0; i < hmesh.size(); i++)
	{
		for (uint j = 0; j < hmesh[i].size(); j++)
		{
			hmesh[i][j].jacobFlag = 0;
		}
	}
	for (uint i = 0; i < rfid.size(); i++)
	{
		hmesh[rfid[i][0]][rfid[i][1]].jacobFlag = 1;
	}
}







////////////////////////////////////////////////////////////////////////////////


void HexQuality::Optimizing_1(int nSize, int nitrMax, double stepSize)
{
	//double edlen_max = MaxEdgeLength();
	//double eta(1.e-3);
	//double tol(eta*edlen_max);

	double tol(1.e-4);

	cout << "\nOptimizing...\n";
	for (int it = 0; it < nSize; it++)
	{
		int flag(0);
		double minJacob_glb;
		int min_pos;
		vector<int> BadEle;
		GlobalMinJacob(minJacob_glb, min_pos, BadEle);

		vector<int> pbad(cp.size(), 0);
		for (uint i = 0; i < BadEle.size(); i++)
		{
			for (int j = 0; j < 8; j++)
			{
				pbad[tmesh[BadEle[i]].cnct[j]] = 1;
			}
		}

		double dispMax(0.), disp;
		for (uint i = 0; i < cp.size(); i++)
		{
			if (cp[i].type != 1 /*&& pbad[i] == 1 && i==8*/)
			{
				disp = OptimizePoint(i, nitrMax, stepSize);
				if (disp > dispMax) dispMax = disp;
				//if (disp > 1.e14)
				//{
				//	getchar();
				//}
			}
			//else
			//{
			//	SmoothingPointBoundary(i, stepSize);
			//}
		}
		//cout << "Max disp: " << dispMax << "\n";
		//getchar();

		double minJacob_glb1;
		int min_pos1;
		vector<int> BadEle1;
		GlobalMinJacob(minJacob_glb1, min_pos1, BadEle1);

		//if (dispMax < tol)
		//{
		//	cout << "Converged after " << it << " steps!\n";
		//	break;
		//}

		for (uint i = 0; i < tmesh.size(); i++)
		{
			tmesh[i].trun = 0;
		}
		for (uint i = 0; i < BadEle1.size(); i++)
		{
			tmesh[BadEle1[i]].trun = 1;
		}

		//if (BadEle1.size() < 80)
		//{
		//	break;
		//}

		//if (flag == 0 /*|| fabs(minJacob_glb - minJacob_glb1) < 1.e-6*/)
		//{
		//	cout << "After " << it << "steps, no interior point moved!\n";
		//	break;
		//}

	}
	cout << "Done optimizing!\n";
	//OutputCM("../io/hex_input/optimize/cube_coarse4");
	//OutputCM("../io/hex_input/optimize/rockerarm2");
	//OutputCM("../io/hex_input/optimize/fertility2");
	//OutputCM("../io/hex_input/optimize/honda1");
}

double HexQuality::MaxEdgeLength()
{
	double lmax(0.), vtmp[3], dist;
	for (uint i = 0; i < tmedge.size(); i++)
	{
		vtmp[0] = cp[tmedge[i].pt[1]].coor[0] - cp[tmedge[i].pt[0]].coor[0];
		vtmp[1] = cp[tmedge[i].pt[1]].coor[1] - cp[tmedge[i].pt[0]].coor[1];
		vtmp[2] = cp[tmedge[i].pt[1]].coor[2] - cp[tmedge[i].pt[0]].coor[2];
		dist = sqrt(vtmp[0] * vtmp[0] + vtmp[1] * vtmp[1] + vtmp[2] * vtmp[2]);
		if (dist > lmax) lmax = dist;
	}
	return lmax;
}

double HexQuality::OptimizePoint(int pid, int nitrMax, double stepSize)
{
	//first translate and scale
	double trs[3];
	vector<double> scl;
	//TranslateScale(pid, trs, scl);
	double grad[3], dir[3], vdisp[3], stepScale, disp, keta;
	for (int it = 0; it < nitrMax; it++)
	{		
		//cout << "\n\n\npid: " << pid << "\n";
		GetAdvanceDirection(pid, grad, dir, keta);
		//stepScale = GetStepSize(pid, grad, dir, keta);
		stepScale = AverageEdgeLength(pid);
		stepScale *= stepSize;
		//stepScale = stepSize;
		double vdisp[3] = { stepScale*dir[0],stepScale*dir[1], stepScale*dir[2] };
		disp = sqrt(vdisp[0] * vdisp[0] + vdisp[1] * vdisp[1] + vdisp[2] * vdisp[2]);
		cp[pid].coor[0] -= vdisp[0];
		cp[pid].coor[1] -= vdisp[1];
		cp[pid].coor[2] -= vdisp[2];
		//cout << pid << ": " << keta << ", " << dir[0] << " " << dir[1] << " " << dir[2] << "\n";
		//getchar();
		//if (disp > 0.1)
		//{
		//	cout << pid << ": " << keta << ", " << dir[0] << " " << dir[1] << " " << dir[2] << "\n";
		//	getchar();
		//}
		//if (disp < 1.e-4)
		//{
		//	break;
		//}
	}
	//translate and scale back
	//TranslateScale_Reverse(pid, trs, scl);

	return disp;
}

void HexQuality::TranslateScale(int pid, double trs[3], vector<double>& scl)
{
	scl.clear();
	double tol(1.e-8);
	trs[0] = cp[pid].coor[0];
	trs[1] = cp[pid].coor[1];
	trs[2] = cp[pid].coor[2];
	cp[pid].coor[0] = 0.;
	cp[pid].coor[1] = 0.;
	cp[pid].coor[2] = 0.;
	double dmax(0.);
	for (uint i = 0; i < cp[pid].edge.size(); i++)
	{
		int edid(cp[pid].edge[i]);
		int edpt(tmedge[edid].pt[0]);
		if (edpt == pid) edpt = tmedge[edid].pt[1];
		//cp[edpt].coor[0] -= trs[0];
		//cp[edpt].coor[1] -= trs[1];
		//cp[edpt].coor[2] -= trs[2];
		double dst = sqrt(cp[edpt].coor[0]* cp[edpt].coor[0]+ cp[edpt].coor[1] * cp[edpt].coor[1] + cp[edpt].coor[2] * cp[edpt].coor[2]);
		if (dst > dmax)
		{
			dmax = dst;
		}
		//if (dst > tol)
		//{
		//	cp[edpt].coor[0] /= dst;
		//	cp[edpt].coor[1] /= dst;
		//	cp[edpt].coor[2] /= dst;
		//	scl.push_back(dst);
		//}
		//else
		//{
		//	//cout << "repeated points!\n";
		//	scl.push_back(1.);
		//	//cp[edpt].coor[0] /= tol;
		//	//cp[edpt].coor[1] /= tol;
		//	//cp[edpt].coor[2] /= tol;
		//	//scl.push_back(tol);
		//}
	}
	//for (uint i = 0; i < cp[pid].edge.size(); i++)
	//{
	//	int edid(cp[pid].edge[i]);
	//	int edpt(tmedge[edid].pt[0]);
	//	if (edpt == pid) edpt = tmedge[edid].pt[1];
	//	cp[edpt].coor[0] /= dmax;
	//	cp[edpt].coor[1] /= dmax;
	//	cp[edpt].coor[2] /= dmax;
	//	scl.push_back(dmax);
	//}

	for (uint i = 0; i < cp[pid].hex.size(); i++)
	{
		int hxid(cp[pid].hex[i]);
		for (int j = 0; j < 8; j++)
		{
			cp[tmesh[hxid].cnct[j]].update = 0;
		}
		scl.push_back(dmax);
	}
	cp[pid].update = 1;
	for (uint i = 0; i < cp[pid].hex.size(); i++)
	{
		int hxid(cp[pid].hex[i]);
		for (int j = 0; j < 8; j++)
		{
			if (cp[tmesh[hxid].cnct[j]].update == 0)
			{
				cp[tmesh[hxid].cnct[j]].update = 1;
				cp[tmesh[hxid].cnct[j]].coor[0] = (cp[tmesh[hxid].cnct[j]].coor[0] - trs[0]) / dmax;
				cp[tmesh[hxid].cnct[j]].coor[1] = (cp[tmesh[hxid].cnct[j]].coor[1] - trs[1]) / dmax;
				cp[tmesh[hxid].cnct[j]].coor[2] = (cp[tmesh[hxid].cnct[j]].coor[2] - trs[2]) / dmax;
			}
		}
	}
}

void HexQuality::TranslateScale_Reverse(int pid, double trs[3], const vector<double>& scl)
{
	cp[pid].coor[0] += trs[0];
	cp[pid].coor[1] += trs[1];
	cp[pid].coor[2] += trs[2];

	//for (uint i = 0; i < cp[pid].edge.size(); i++)
	//{
	//	int edid(cp[pid].edge[i]);
	//	int edpt(tmedge[edid].pt[0]);
	//	if (edpt == pid) edpt = tmedge[edid].pt[1];
	//	cp[edpt].coor[0] = scl[i] * cp[edpt].coor[0] + trs[0];
	//	cp[edpt].coor[1] = scl[i] * cp[edpt].coor[1] + trs[1];
	//	cp[edpt].coor[2] = scl[i] * cp[edpt].coor[2] + trs[2];
	//}

	for (uint i = 0; i < cp[pid].hex.size(); i++)
	{
		int hxid(cp[pid].hex[i]);
		for (int j = 0; j < 8; j++)
		{
			cp[tmesh[hxid].cnct[j]].update = 0;
		}
	}
	cp[pid].update = 1;
	for (uint i = 0; i < cp[pid].hex.size(); i++)
	{
		int hxid(cp[pid].hex[i]);
		for (int j = 0; j < 8; j++)
		{
			if (cp[tmesh[hxid].cnct[j]].update == 0)
			{
				cp[tmesh[hxid].cnct[j]].update = 1;
				cp[tmesh[hxid].cnct[j]].coor[0] = scl[i] * cp[tmesh[hxid].cnct[j]].coor[0] + trs[0];
				cp[tmesh[hxid].cnct[j]].coor[1] = scl[i] * cp[tmesh[hxid].cnct[j]].coor[1] + trs[1];
				cp[tmesh[hxid].cnct[j]].coor[2] = scl[i] * cp[tmesh[hxid].cnct[j]].coor[2] + trs[2];
			}
		}
	}
}

void HexQuality::GetAdvanceDirection(int pid, double dir[3], double dir1[3], double& keta)
{
	double alpha(1.e-3);
	double eps(1.e-4);
	double detJmin(1.e5), detJtmp, postmp[3];
	for (uint i = 0; i < cp[pid].hex.size(); i++)
	{
		GetHexMinJacob(cp[pid].hex[i], detJtmp, postmp);
		if (detJtmp < detJmin)
		{
			detJmin = detJtmp;
		}
	}
	double delta(-1.);
	if (detJmin < eps)
	{
		double beta(1.e-2);
		double to = alpha*(fabs(detJmin));
		delta = sqrt(to*to + to*fabs(detJmin));
	}
	//cout << "delta: " << delta << "\n";

	dir[0] = 0.; dir[1] = 0.; dir[2] = 0.;
	dir1[0] = 0.; dir1[1] = 0.; dir1[2] = 0.;
	keta = 0.;
	//vector<double> eeta(cp[pid].hex.size());
	MatrixXd hsmat = MatrixXd::Zero(3, 3);
	for (uint i = 0; i < cp[pid].hex.size(); i++)
	{
		int eid(cp[pid].hex[i]);
		//cout << "eid: " << eid << "\n";
		int* it = find(tmesh[eid].cnct, tmesh[eid].cnct + 8, pid);
		int loc(it - tmesh[eid].cnct);
		double obj, grad[3], grad2[3][3];
		objGrad(eid, loc, delta, obj, grad, grad2);
		dir[0] += obj*grad[0];
		dir[1] += obj*grad[1];
		dir[2] += obj*grad[2];
		keta += (obj*obj);
		//eeta[i] = obj;
		for (int k1 = 0; k1 < 3; k1++)
		{
			for (int k2 = 0; k2 < 3; k2++)
			{
				hsmat(k1, k2) += (grad[k1] * grad[k2] + obj*grad2[k1][k2]);
			}
		}
	}
	//double tmp(2./double(cp[pid].hex.size()));
	//dir[0] *= tmp;
	//dir[1] *= tmp;
	//dir[2] *= tmp;
	//keta /= double(cp[pid].hex.size());

	//for (int k1 = 0; k1 < 3; k1++)
	//{
	//	for (int k2 = 0; k2 < 3; k2++)
	//	{
	//		hsmat(k1, k2) *= tmp;
	//	}
	//}
	//double hsdet = hsmat.determinant();
	////cout << hsdet << "\n"; getchar();
	//MatrixXd hsinv = hsmat.inverse();
	////cout << hsinv << "\n"; getchar();
	//for (int k1 = 0; k1 < 3; k1++)
	//{
	//	for (int k2 = 0; k2 < 3; k2++)
	//	{
	//		dir1[k1] += (hsinv(k1, k2) * dir[k2]);
	//	}
	//	dir1[k1] = dir[k1];
	//}

	double dst = sqrt(dir[0] * dir[0] + dir[1] * dir[1] + dir[2] * dir[2]);
	dir[0] /= dst;
	dir[1] /= dst;
	dir[2] /= dst;
	dir1[0] = dir[0];
	dir1[1] = dir[1];
	dir1[2] = dir[2];
	//if (pid == 0)
	//{
	//	cout << "pid obj: " << pid << " " << keta << "\n";
	//	//cout << "element eta: ";
	//	//for (uint i = 0; i < eeta.size(); i++)
	//	//{
	//	//	cout << eeta[i] << " ";
	//	//}
	//	//cout << "\n";
	//	getchar();
	//}
	//cout << dir[0] << " " << dir[1] << " " << dir[2] << "\n";
	//cout << "pid obj: " << pid << " " << keta << "\n";
	//getchar();
}

double HexQuality::GetStepSize(int pid, double grad[3], double dir[3], double keta)
{
	double stepSize(1.);
	double c1(1.e-4);
	double c2(0.9);
	double eps(1.e-16);
	double rho(0.5);
	double coor[3] = { cp[pid].coor[0],cp[pid].coor[1], cp[pid].coor[2] };

	double norm_pk = grad[0] * dir[0] + grad[1] * dir[1] + grad[2] * dir[2];
	double Fleft, Fright;
	double gradtmp[3], dirtmp[3], normtmp;
	while (1)
	{
		Fright = keta - c1*stepSize*norm_pk;
		cp[pid].coor[0] = coor[0] - stepSize*dir[0];
		cp[pid].coor[1] = coor[1] - stepSize*dir[1];
		cp[pid].coor[2] = coor[2] - stepSize*dir[2];
		GetAdvanceDirection(pid, gradtmp, dirtmp, Fleft);
		normtmp = gradtmp[0] * dir[0] + gradtmp[1] * dir[1] + gradtmp[2] * dir[2];
		//cout << Fleft << " " << Fright << "\n";
		//getchar();
		if (Fleft <= Fright && fabs(normtmp) < fabs(c2*norm_pk))
		{
			break;
		}
		stepSize *= rho;
		if (stepSize < eps)
		{
			stepSize = 0.;
			break;
		}
	}
	cp[pid].coor[0] = coor[0];
	cp[pid].coor[1] = coor[1];
	cp[pid].coor[2] = coor[2];

	return stepSize;
}

void HexQuality::objGrad(int eid, int iloc, double delta, double& obj, double grad[3], double grad2[3][3])
{
	//const double eps(1.e-4);
	//double delta(1.e-3);
	double Jmat[3][3], DJmat[3][3][3], detJ, detJgrad[3];
	GetJacobMat(eid, iloc, Jmat, DJmat, detJ, detJgrad);
	double Fnorm2(0.), DFnorm2[3] = {0.,0.,0.};
	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			Fnorm2 += Jmat[i][j] * Jmat[i][j];
			DFnorm2[0] += DJmat[0][i][j] * Jmat[i][j];
			DFnorm2[1] += DJmat[1][i][j] * Jmat[i][j];
			DFnorm2[2] += DJmat[2][i][j] * Jmat[i][j];
		}
	}
	double h, coef, tmp;
	//if (detJ < eps)
	//if (delta > 0.)
	//{
	//	tmp = sqrt(detJ*detJ + 4.*delta*delta);
	//	h = .5*(detJ + tmp);
	//	coef = 1. / tmp;
	//}
	//else
	{
		h = detJ;
		coef = 1. / detJ;
	}
	obj = Fnorm2 / (3.*pow(h, 2. / 3.));
	grad[0] = 2.*obj*(DFnorm2[0] / Fnorm2 - detJgrad[0] * coef / 3.);
	grad[1] = 2.*obj*(DFnorm2[1] / Fnorm2 - detJgrad[1] * coef / 3.);
	grad[2] = 2.*obj*(DFnorm2[2] / Fnorm2 - detJgrad[2] * coef / 3.);
	//cout << "DFnorm2: " << DFnorm2[0] << " " << DFnorm2[1] << " " << DFnorm2[2] << "\n";
	//cout << "detGrad: " << detJgrad[0] << " " << detJgrad[1] << " " << detJgrad[2] << "\n";
	//cout << "Fnorm2 h: " << Fnorm2 << " " << h << " " << detJ << " " << tmp << "\n";
	//cout << "obj: " << obj << "\n";
	//cout << "grad: " << grad[0] << " " << grad[1] << " " << grad[2] << "\n";
	//getchar();

	double Fnorm4(Fnorm2*Fnorm2), coef3(pow(coef, 3.));
	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			grad2[i][j] = 0;
			double DF2(0.);
			for (int k1 = 0; k1 < 3; k1++)
			{
				for (int k2 = 0; k2 < 3; k2++)
				{
					DF2 += (DJmat[i][k1][k2] * DJmat[j][k1][k2]);
				}
			}
			grad2[i][j] = grad[i] * grad[j] / 3. + 2.*obj*(DF2 / Fnorm2 - 2.*DF2 / Fnorm4 + detJ*detJgrad[i] * detJgrad[j] * coef3 / 3.);
		}
	}
}

void HexQuality::GetJacobMat(int eid, int iloc, double Jmat[3][3], double DJmat[3][3][3], double& detJ, double grad[3])
{
	double u(0.), v(0.), w(0.);
	if (iloc == 1 || iloc == 2 || iloc == 5 || iloc == 6) u = 1.;
	if (iloc == 2 || iloc == 3 || iloc == 6 || iloc == 7) v = 1.;
	if (iloc == 4 || iloc == 5 || iloc == 6 || iloc == 7) w = 1.;
	double Nu[2] = { 1. - u,u };
	double Nv[2] = { 1. - v,v };
	double Nw[2] = { 1. - w,w };
	double dNdu[2] = { -1.,1. };
	double dNdv[2] = { -1.,1. };
	double dNdw[2] = { -1.,1. };
	double dNdt[8][3];
	int ploc[8] = { 0,1,3,2,4,5,7,6 };
	int loc(0);
	for (int i = 0; i < 2; i++)
	{
		for (int j = 0; j < 2; j++)
		{
			for (int k = 0; k < 2; k++)
			{
				dNdt[ploc[loc]][0] = dNdu[k] * Nv[j] * Nw[i];
				dNdt[ploc[loc]][1] = Nu[k] * dNdv[j] * Nw[i];
				dNdt[ploc[loc]][2] = Nu[k] * Nv[j] * dNdw[i];
				loc++;
			}
		}
	}
	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			Jmat[i][j] = 0.;
			for (int k = 0; k < 8; k++)
			{
				Jmat[i][j] += cp[tmesh[eid].cnct[k]].coor[i] * dNdt[k][j];
			}
		}
	}
	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			DJmat[0][i][j] = 0.;
			DJmat[1][i][j] = 0.;
			DJmat[2][i][j] = 0.;
		}
	}
	DJmat[0][0][0] = dNdt[iloc][0]; DJmat[0][0][1] = dNdt[iloc][1]; DJmat[0][0][2] = dNdt[iloc][2];
	DJmat[1][1][0] = dNdt[iloc][0]; DJmat[1][1][1] = dNdt[iloc][1]; DJmat[1][1][2] = dNdt[iloc][2];
	DJmat[2][2][0] = dNdt[iloc][0]; DJmat[2][2][1] = dNdt[iloc][1]; DJmat[2][2][2] = dNdt[iloc][2];
	double cng[3][3] = { {Jmat[1][1] * Jmat[2][2] - Jmat[1][2] * Jmat[2][1],Jmat[1][0] * Jmat[2][2] - Jmat[1][2] * Jmat[2][0], Jmat[1][0] * Jmat[2][1] - Jmat[1][1] * Jmat[2][0] },
	{ Jmat[0][1] * Jmat[2][2] - Jmat[0][2] * Jmat[2][1],Jmat[0][0] * Jmat[2][2] - Jmat[0][2] * Jmat[2][0], Jmat[0][0] * Jmat[2][1] - Jmat[0][1] * Jmat[2][0] }, 
	{ Jmat[0][1] * Jmat[1][2] - Jmat[0][2] * Jmat[1][1],Jmat[0][0] * Jmat[1][2] - Jmat[0][2] * Jmat[1][0], Jmat[0][0] * Jmat[1][1] - Jmat[0][1] * Jmat[1][0] } };
	
	detJ = Jmat[0][0] * cng[0][0] - Jmat[0][1] * cng[0][1] + Jmat[0][2] * cng[0][2];
	grad[0] = dNdt[iloc][0] * cng[0][0] - dNdt[iloc][1] * cng[0][1] + dNdt[iloc][2] * cng[0][2];
	grad[1] = -(dNdt[iloc][0] * cng[1][0] - dNdt[iloc][1] * cng[1][1] + dNdt[iloc][2] * cng[1][2]);
	grad[2] = dNdt[iloc][0] * cng[2][0] - dNdt[iloc][1] * cng[2][1] + dNdt[iloc][2] * cng[2][2];
}

void HexQuality::Optimizing_glb(int nStep, double stepSize)
{
	double beta(0.2);
}







/////////////////////////////////////////////////////////////////////////////////////////////

void HexQuality::LaplaceSmoothing(int nstep)
{
	cout << "\nLaplace Smoothing...\n";

	double tol(1.e-3);
	for (int it = 0; it < nstep; it++)
	{
		cout << "istep: " << it << "\n";
		double minJacob_glb0;
		int min_pos0(0);
		vector<int> BadEle0;
		GlobalMinJacob(minJacob_glb0, min_pos0, BadEle0);

		for (uint i = 0; i < cp.size(); i++)
		{
			if (cp[i].type != 1)
			{
				LaplaceSmooth_Interior(i);
			}
			//else
			//{
			//	if (cp[i].sharp == 0)
			//	{
			//		LaplaceSmooth_Boundary_NonSharp(i);
			//	}
			//	else if (cp[i].sharp == 1)
			//	{
			//		LaplaceSmooth_Boundary_Sharp(i);
			//	}
			//}
		}

		double minJacob_glb;
		int min_pos;
		vector<int> BadEle;
		GlobalMinJacob(minJacob_glb, min_pos, BadEle);

		//if (fabs(minJacob_glb - minJacob_glb0) < tol && BadEle0.size() == BadEle.size())
		//{
		//	cout << "Converged!\n"; 
		//	break;
		//}
		//else if (minJacob_glb < minJacob_glb0)
		//{

		//}

		for (uint i = 0; i < tmesh.size(); i++) tmesh[i].trun = 0;
		for (uint i = 0; i < BadEle.size(); i++) tmesh[BadEle[i]].trun = 1;
	}

	cout << "Done smoothing!\n";
}

void HexQuality::LaplaceSmooth_Interior(int pid)
{
	//before smoothing
	double minJacob_0(1.e5), minJacob_1(1.e5), tmpJacob, GaussPos[3];
	for (uint i = 0; i < cp[pid].hex.size(); i++)
	{
		GetHexMinJacob(cp[pid].hex[i], tmpJacob, GaussPos);
		if (tmpJacob < minJacob_0) minJacob_0 = tmpJacob;
	}

	//smoothing
	double ptmp[3] = { 0.,0.,0., };
	double pold[3] = { cp[pid].coor[0],cp[pid].coor[1], cp[pid].coor[2] };
	for (uint i = 0; i < cp[pid].edge.size(); i++)
	{
		int edid(cp[pid].edge[i]);
		int pnb(tmedge[edid].pt[0]);
		if (pnb == pid) pnb = tmedge[edid].pt[1];
		ptmp[0] += cp[pnb].coor[0];
		ptmp[1] += cp[pnb].coor[1];
		ptmp[2] += cp[pnb].coor[2];
	}
	cp[pid].coor[0] = ptmp[0] / cp[pid].edge.size();
	cp[pid].coor[1] = ptmp[1] / cp[pid].edge.size();
	cp[pid].coor[2] = ptmp[2] / cp[pid].edge.size();

	//after smoothing
	for (uint i = 0; i < cp[pid].hex.size(); i++)
	{
		GetHexMinJacob(cp[pid].hex[i], tmpJacob, GaussPos);
		if (tmpJacob < minJacob_1) minJacob_1 = tmpJacob;
	}

	//if (minJacob_1 < minJacob_0)
	//{
	//	cp[pid].coor[0] = pold[0];
	//	cp[pid].coor[1] = pold[1];
	//	cp[pid].coor[2] = pold[2];
	//}
}

void HexQuality::LaplaceSmooth_Boundary_NonSharp(int pid)
{
	//before smoothing
	double minJacob_0(1.e5), minJacob_1(1.e5), tmpJacob, GaussPos[3];
	for (uint i = 0; i < cp[pid].hex.size(); i++)
	{
		GetHexMinJacob(cp[pid].hex[i], tmpJacob, GaussPos);
		if (tmpJacob < minJacob_0) minJacob_0 = tmpJacob;
	}

	//smoothing
	double ptmp[3] = { 0.,0.,0., };
	double pold[3] = { cp[pid].coor[0],cp[pid].coor[1], cp[pid].coor[2] };
	int nedb(0);
	for (uint i = 0; i < cp[pid].edge.size(); i++)
	{
		if (tmedge[cp[pid].edge[i]].type == 1)
		{
			int edid(cp[pid].edge[i]);
			int pnb(tmedge[edid].pt[0]);
			if (pnb == pid) pnb = tmedge[edid].pt[1];
			ptmp[0] += cp[pnb].coor[0];
			ptmp[1] += cp[pnb].coor[1];
			ptmp[2] += cp[pnb].coor[2];
			nedb++;
		}
	}
	cp[pid].coor[0] = ptmp[0] / double(nedb);
	cp[pid].coor[1] = ptmp[1] / double(nedb);
	cp[pid].coor[2] = ptmp[2] / double(nedb);

	//after smoothing
	for (uint i = 0; i < cp[pid].hex.size(); i++)
	{
		GetHexMinJacob(cp[pid].hex[i], tmpJacob, GaussPos);
		if (tmpJacob < minJacob_1) minJacob_1 = tmpJacob;
	}

	//if (minJacob_1 < minJacob_0)
	//{
	//	cp[pid].coor[0] = pold[0];
	//	cp[pid].coor[1] = pold[1];
	//	cp[pid].coor[2] = pold[2];
	//}
}

void HexQuality::LaplaceSmooth_Boundary_Sharp(int pid)
{
	//before smoothing
	double minJacob_0(1.e5), minJacob_1(1.e5), tmpJacob, GaussPos[3];
	for (uint i = 0; i < cp[pid].hex.size(); i++)
	{
		GetHexMinJacob(cp[pid].hex[i], tmpJacob, GaussPos);
		if (tmpJacob < minJacob_0) minJacob_0 = tmpJacob;
	}

	//smoothing
	double ptmp[3] = { 0.,0.,0., };
	double pold[3] = { cp[pid].coor[0],cp[pid].coor[1], cp[pid].coor[2] };
	int nedb(0);
	for (uint i = 0; i < cp[pid].edge.size(); i++)
	{
		if (tmedge[cp[pid].edge[i]].sharp == 1)
		{
			int edid(cp[pid].edge[i]);
			int pnb(tmedge[edid].pt[0]);
			if (pnb == pid) pnb = tmedge[edid].pt[1];
			ptmp[0] += cp[pnb].coor[0];
			ptmp[1] += cp[pnb].coor[1];
			ptmp[2] += cp[pnb].coor[2];
			nedb++;
		}
	}
	cp[pid].coor[0] = ptmp[0] / double(nedb);
	cp[pid].coor[1] = ptmp[1] / double(nedb);
	cp[pid].coor[2] = ptmp[2] / double(nedb);

	//after smoothing
	for (uint i = 0; i < cp[pid].hex.size(); i++)
	{
		GetHexMinJacob(cp[pid].hex[i], tmpJacob, GaussPos);
		if (tmpJacob < minJacob_1) minJacob_1 = tmpJacob;
	}

	if (minJacob_1 < minJacob_0)
	{
		cp[pid].coor[0] = pold[0];
		cp[pid].coor[1] = pold[1];
		cp[pid].coor[2] = pold[2];
	}
}



bool HexQuality::ComputeElementNormal_Quad(int cnct[4], array<double, 3>& normal_quad)
{

	double vector_1[3], vector_2[3];


	vector_1[0] = hcp[0][cnct[1]].coor[0] - hcp[0][cnct[0]].coor[0];
	vector_1[1] = hcp[0][cnct[1]].coor[1] - hcp[0][cnct[0]].coor[1];
	vector_1[2] = hcp[0][cnct[1]].coor[2] - hcp[0][cnct[0]].coor[2];

	vector_2[0] = hcp[0][cnct[3]].coor[0] - hcp[0][cnct[0]].coor[0];
	vector_2[1] = hcp[0][cnct[3]].coor[1] - hcp[0][cnct[0]].coor[1];
	vector_2[2] = hcp[0][cnct[3]].coor[2] - hcp[0][cnct[0]].coor[2];
	//cout <<"Test"<< vector_1[0] << endl;
	CrossProduct(vector_1, vector_2, normal_quad);
	//cout << "Test" << normal_quad[0] << endl;


	return true;
}

bool HexQuality::CrossProduct(double vector_1[3], double vector_2[3], array<double, 3>& result)
{
	result[0] = vector_1[1] * vector_2[2] - vector_1[2] * vector_2[1];
	result[1] = vector_1[2] * vector_2[0] - vector_1[0] * vector_2[2];
	result[2] = vector_1[0] * vector_2[1] - vector_1[1] * vector_2[0];
	/*cout << "Test_result" << result[0] << endl;*/
	return true;
}

bool HexQuality::DotProduct(vector<array<double, 3>> check_normal_surface_from_edge, double &angle)
{
	double x1 = check_normal_surface_from_edge[0][0];
	double y1 = check_normal_surface_from_edge[0][1];
	double z1 = check_normal_surface_from_edge[0][2];
	double x2 = check_normal_surface_from_edge[1][0];
	double y2 = check_normal_surface_from_edge[1][1];
	double z2 = check_normal_surface_from_edge[1][2];


	/*cout << "test" << endl;
	for (uint k = 0; k < 3; k++)
	{
		cout << check_normal_surface_from_edge[0][k] << "  " ;
	}
	cout << endl;

	cout << "test" << endl;
	for (uint k = 0; k < 3; k++)
	{
		cout << check_normal_surface_from_edge[1][k] << "  ";
	}
	cout << endl;

	getchar();*/

	double dot = x1 * x2 + y1 * y2 + z1 * z2;
	double lenSq1 = x1 * x1 + y1 * y1 + z1 * z1;
	double lenSq2 = x2 * x2 + y2 * y2 + z2 * z2;
	angle = acos(dot / sqrt(lenSq1 * lenSq2));

	/*cout << "test2" << endl;
	cout<<dot<<" "<< lenSq1 <<" "<< lenSq2<<endl;


	cout << angle << endl;

	getchar();*/

	return true;
}

