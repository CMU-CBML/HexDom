//
// Author: Xinghua Liang
//
#pragma once
#ifndef _GeometricMethod_Liang_
#define _GeometricMethod_Liang_

#include <math.h>

class GeometricMethod
{
public:
	GeometricMethod(): ERR(1.0E-10), ERR2(1.0E-20), ERR3(1.0E-30), PI(3.14159265359),
					   E(2.71828183), A2R(0.01745329252), R2A(57.295779513), MAX(1.0E+100), MIN(-1.0E+100) {}
	~GeometricMethod() {}

	//Angle of P1-P0-P2
	double	angle(double p0[3], double p1[3], double p2[3]);	
	double	angle(double xo, double yo, double zo, double x1, double y1, double z1, double x2, double y2, double z2);
	float	angle(float p0[3], float p1[3], float p2[3])
			{return (float) angle((double) p0[0], (double) p0[1], (double) p0[2],
								  (double) p1[0], (double) p1[1], (double) p1[2],
								  (double) p2[0], (double) p2[1], (double) p2[2]);}
	double	angle(double rad)
	{
		return (rad*R2A);
	}
	int		angle(double rad, int error);

	double	angle3D(double p0[3], double p1[3], double p2[3]);


	//cp: center point; fp: former point; np: next point; newp: new point
	void	angularBisector(double cp[3], double fp[3], double np[3], double newp[3], double dist = 1.0);
	void	angularBisector(double cx, double cy, double fx, double fy, double nx, double ny, double &xx, double &yy, double dist = 1.0);

	void	boundingBox(double p1[3], double p2[3], double boxLow[3], double boxUp[3]);
	void	boundingBox(double p1[3], double p2[3], double p3[3], double boxLow[3], double boxUp[3]);

	void	centerPoint(double p1[3], double p2[3], double center[3])
	{
		center[0] = (p1[0] + p2[0])/2;
		center[1] = (p1[1] + p2[1])/2;
		center[2] = (p1[2] + p2[2])/2;
	}

	void	crossProduct(double p1[3], double p2[3], double result[3])
	{
		result[0] = p1[1]*p2[2] - p1[2]*p2[1];
		result[1] = p1[2]*p2[0] - p1[0]*p2[2];
		result[2] = p1[0]*p2[1] - p1[1]*p2[0];
	}
	void	crossProduct(float p1[3], float p2[3], float result[3])
	{
		result[0] = p1[1]*p2[2] - p1[2]*p2[1];
		result[1] = p1[2]*p2[0] - p1[0]*p2[2];
		result[2] = p1[0]*p2[1] - p1[1]*p2[0];
	}
	void	crossProduct(double x1, double y1, double z1, double x2, double y2, double z2,  double &nx, double &ny, double &nz)
	{
		nx = y1 * z2 - y2 * z1;	ny = z1 * x2 - z2 * x1; nz = x1 * y2 - x2 * y1;
	}


	void	curveMiddle(double p1[3], double p2[3], double p3[3], double pMid[3])
	{
		double s1, sMid = 0.5;

		s1 = distance(p1, p2);
		s1 = s1/(s1 + distance(p2, p3));
		for (int i=0; i<3; i++)
		{
			if (s1 < 0.1)
				pMid[i] = (p2[i] + p3[i])/2;
			else if (s1 > 0.9)
				pMid[i] = (p1[i] + p2[i])/2;
			else
				pMid[i] = (p1[i]*(1-s1)*(s1-sMid)*(sMid-1) + sMid*(p2[i]*(sMid-1) + p3[i]*s1*(s1-sMid)))/(s1*s1 - s1);
		}
	}


	double	distance(double p1[3], double p2[3])
	{
		return sqrt((p1[0] - p2[0])*(p1[0] - p2[0]) + (p1[1] - p2[1])*(p1[1] - p2[1]) + (p1[2] - p2[2])*(p1[2] - p2[2]));
	}
	float	distance(float p1[3], float p2[3])
	{
		return sqrt((p1[0] - p2[0])*(p1[0] - p2[0]) + (p1[1] - p2[1])*(p1[1] - p2[1]) + (p1[2] - p2[2])*(p1[2] - p2[2]));
	}
	double	distance(double x0, double y0, double x1, double y1)
	{
		return sqrt((x0 - x1)*(x0 - x1) + (y0 - y1)*(y0 - y1));
	}
	double	distance(double x0, double y0, double z0, double x1, double y1, double z1)
	{
		return sqrt((x0 - x1)*(x0 - x1) + (y0 - y1)*(y0 - y1) + (z0 - z1)*(z0 - z1));
	}

	//rad angle between two plane
	double dihedralAngle(double plane1_Center[3], double plane1_Prev[3], double plane1_Next[3],
						 double plane2_Center[3], double plane2_Prev[3], double plane2_Next[3])
	{
		double normal_1[3], normal_2[3];

		planeNormal(plane1_Center, plane1_Prev, plane1_Next, normal_1);
		planeNormal(plane2_Center, plane2_Prev, plane2_Next, normal_2);

		double coss = (normal_1[0] * normal_2[0] + normal_1[1] * normal_2[1] + normal_1[2] * normal_2[2])
					  /(norm(normal_1)*norm(normal_2));
		coss = coss>1.0 ? 1.0 : coss;
		coss = coss<-1.0 ? -1.0 : coss;
		return acos(coss);
	}

	double	dotProduct(float p1[3], float p2[3])
	{
		return (p1[0]*p2[0] + p1[1]*p2[1] + p1[2]*p2[2]);
	}
	double	dotProduct(double p1[3], double p2[3])
	{
		return (p1[0]*p2[0] + p1[1]*p2[1] + p1[2]*p2[2]);
	}
	double	dotProduct(double x1, double y1, double z1, double x2, double y2, double z2)
	{
		return (x1*x2 + y1*y2 + z1*z2);
	}

	bool	isTwoBoundingBoxNotIntersected(double box1_Low[3], double box1_Up[3], double box2_Low[3], double box2_Up[3]);

	void	intersectingPointOfThreePlane(double coef1[4], double coef2[4], double coef3[4], double point[3]);

	//For 2D
	//0: not inersected; 1: line 1 touch line 2; 2: intersected.
	int		lineSegmentIntersectingLineSegment(double x11, double y11, double x12, double y12, double x21, double y21, double x22, double y22);

	//For 3D
	//0: not inersected; 1: line touch plane 2; 2: intersected.
	int		lineSegmentIntersectingPlane(double line1[3], double line2[3], double pPrev[3], double pCenter[3], double pNext[3], double pos[3]);

	//For 3D
	//0: not inersected; 1: line touch point; 2: line intersect line; 3: intersected.
	int		lineSegmentIntersectingTriangle(double line1[3], double line2[3], double pPrev[3], double pCenter[3], double pNext[3], double pos[3]);

	//For 2D
	bool	lineIntersectingLine(double x11, double y11, double x12, double y12, double x21, double y21, double x22, double y22, double &xx, double &yy);
	//bool lineIntersectingLine(double x11, double y11, double x12, double y12, double x21, double y21, double x22, double y22, double &xx, double &yy);


	void	move(double movep[3], double refp[3], double dist)
	{
		dist /= distance(movep, refp);
		movep[0] = refp[0] + (movep[0] - refp[0])*dist;
		movep[1] = refp[1] + (movep[1] - refp[1])*dist;
		movep[2] = refp[2] + (movep[2] - refp[2])*dist;
		return;
	}

	//Normalize
	void	normalize(double *vector)
	{
		double sum = sqrt(vector[0]*vector[0] + vector[1]*vector[1] + vector[2]*vector[2]);
		if (sum == 0.0)
			return;
		vector[0] /= sum; vector[1] /= sum; vector[2] /= sum;
		return;
	}
	void	normalize(float *vector)
	{
		float sum = sqrt(vector[0]*vector[0] + vector[1]*vector[1] + vector[2]*vector[2]);
		if (sum == 0.0)
			return;
		vector[0] /= sum; vector[1] /= sum; vector[2] /= sum;
		return;
	}

	//Angle within 0 - 360
	double	normalizeAngle(double angle)
	{
		while (angle >= 360.0)
		{
			angle -= 360.0;
			return angle;
		}
		while (angle < 0.0)
			angle += 360.0;
		return angle;
	}

	//Norm of a vector
	double norm(double vector[3])
	{
		return sqrt(vector[0]*vector[0] + vector[1]*vector[1] + vector[2]*vector[2]);
	}

	//project for 2D
	void	project(double xp, double yp, double x1, double y1, double x2, double y2, double &xx, double &yy)
	{
		double kk;

		if (fabs(x1 - x2) < ERR)
		{
			xx = x1;
			yy = yp;
		}
		else if (fabs(y1 - y2) < ERR)
		{
			xx = xp;
			yy = y1;
		}
		else
		{
			kk = (y1 - y2)/(x1 - x2);
			yy = (kk*xp + kk*kk*yp - kk*x1 + y1)/(1+ kk*kk);
			xx = x1 - (y1 - yy)/kk;
		}
	}

	double	project2Line(double p[3], double lineP1[3], double lineP2[3], double projectP[3])
	{
		double v1[3], v2[3], len, len12;
		int i;

		for (i=0; i<3; i++)
		{
			v1[i] = p[i] - lineP1[i];
			v2[i] = lineP2[i] - lineP1[i];
		}

		len12 = norm(v2);
		normalize(v2);
		len = dotProduct(v1, v2);
		for (i=0; i<3; i++)
			projectP[i] = lineP1[i] + len/len12*(lineP2[i] - lineP1[i]);

		if (len < 0.0 || len > len12)
			return -sqrt(norm(v1)*norm(v1) - len*len);
		return sqrt(norm(v1)*norm(v1) - len*len);
	}
	float project2Line(float p[3], float lineP1[3], float lineP2[3], float projectP[3])
	{
		double fp[3], fLineP1[3], fLineP2[3], dProjectP[3], result;

		fp[0] = (double) p[0]; fp[1] = (double) p[1]; fp[2] = (double) p[2];
		fLineP1[0] = (double) lineP1[0]; fp[1] = (double) lineP1[1]; fp[2] = (double) lineP1[2];
		fLineP2[0] = (double) lineP2[0]; fp[1] = (double) lineP2[1]; fp[2] = (double) lineP2[2];

		result = project2Line(fp, fLineP1, fLineP2, dProjectP);

		projectP[0] = (float) dProjectP[0]; projectP[1] = (float) dProjectP[1]; projectP[2] = (float) dProjectP[2];
		return (float) result;
	}

	//For 3D
	//Return the distance from the point to the plane. Will be negtive if 
	double	projectToPlane(double point[3], double pPrev[3], double pCenter[3], double pNext[3], double pProject[3]);

	//For 3D
	//0: outside the triangle; 1: project to vertex; 2: project to line; 3: project to interior
	int		projectToTriangle(double point[3], double pPrev[3], double pCenter[3], double pNext[3], double pProject[3], double &pDist);
	int		projectToTriangle(double point[3], double pPrev[3], double pCenter[3], double pNext[3], double pProject[3])
	{
		double dist;
		return projectToTriangle(point, pPrev, pCenter, pNext, pProject, dist);
	}

	void	parallel(double xp, double yp, double x1, double y1, double x2, double y2, double &xx, double &yy)
	{
		double kk;

		if (fabs(x1 - x2) < ERR)	//vertical
		{
			xx = xp;
			yy = yp + y2 - y1;
		}
		else if (fabs(y1 - y2) < ERR)	//horizontal
		{
			xx = xp + x2 - x1;
			yy = yp;
		}
		else
		{
			kk = (y1 - y2)/(x1 - x2);
			xx = xp + x2 - x1;
			yy = (xx - xp)*kk + yp;
		}
	}

	void	planeBisector(double p1Center[3], double p1Prev[3], double p1Next[3], double p2Center[3], double p2Prev[3], double p2Next[3], double coef[4]);
	void	planeEquation(double pCenter[3], double pPrevious[3], double pNext[3], double pCoef[4])
	{
		//Plane equation: coef[0]*x + coef[1]*y + coef[2]*z + coef[3] = 0

		planeNormal(pCenter, pPrevious, pNext, pCoef);
		pCoef[3] = -(pCoef[0]*pCenter[0] + pCoef[1]*pCenter[1] + pCoef[2]*pCenter[2]);
	}
	void	planeNormal(double pCenter[3], double pPrevious[3], double pNext[3], double pNormal[3])
	{
		crossProduct(pNext[0] - pCenter[0], pNext[1] - pCenter[1], pNext[2] - pCenter[2],
			pPrevious[0] - pCenter[0], pPrevious[1] - pCenter[1], pPrevious[2] - pCenter[2], pNormal[0], pNormal[1], pNormal[2]);
	}

	void	perpendicular(double x1, double y1, double x2, double y2, double &xx, double &yy);

	void	QuickSort(double* &x, int* &I, int left, int right);
	void	QuickSort(int* &I1, int* &I2, int* &I3, int* &I4, int* &II, int left, int right);
	void	QuickSort(int* &I1, int* &I2, int* &II, int left, int right);

	//angle to rad
	double	rad(double angle)
		{return (angle*A2R);}

	void	rotate(double x1, double y1, double x2, double y2, double angle, double &xx, double &yy);

	void	sort(double * &matrix, int const dim);

	//triangle Area
	double	triangleArea(double p1[3], double p2[3], double p3[3])
	{
		double nx, ny, nz;
		crossProduct(p2[0] - p1[0], p2[1] - p1[1], p2[2] - p1[2], p3[0] - p1[0], p3[1] - p1[1], p3[2] - p1[2], nx, ny, nz);
		return 0.5 * sqrt(nx * nx + ny * ny + nz * nz);
	}
	double	triangleArea(double x1, double y1, double z1, double x2, double y2, double z2, double x3, double y3, double z3)
	{
		double nx, ny, nz;
		crossProduct(x2 - x1, y2 - y1, z2 - z1, x3 - x1, y3 - y1, z3 - z1, nx, ny, nz);
		return 0.5 * sqrt(nx * nx + ny * ny + nz * nz);
	}

	void	zero(int vector[], int size)
	{
		for (int i=0; i<size; ++i)
			vector[i] = 0;
		return;
	}
	void	zero(double vector[], int size)
	{
		for (int i=0; i<size; ++i)
			vector[i] = 0.0;
		return;
	}


	double const PI;
	double const ERR;
	double const ERR2;
	double const ERR3;
	double const E;
	double const A2R;
	double const R2A;
	double const MAX;
	double const MIN;
};

#endif _GeometricMethod_Liang_