#include "GeometricMethod.h"


double GeometricMethod::angle(double xo, double yo, double zo, 
						 double x1, double y1, double z1, 
						 double x2, double y2, double z2)
{
	double ax = x1 - xo; double ay = y1 - yo; double az = z1 - zo;
    double bx = x2 - xo; double by = y2 - yo; double bz = z2 - zo;

    double a = distance(xo, yo, zo, x1, y1, z1);
    double b = distance(xo, yo, zo, x2, y2, z2);
    double cos = (ax * bx + ay * by + az * bz) / (a * b);

	double nx, ny, nz;
	nx = x1 - x2; ny = y1 - y2; nz = z1 - z2;
	crossProduct(-bx, -by, -bz, nx, ny, nz, nx, ny, nz);
	//crossProduct(ax, ay, az, bx, by, bz, nx, ny, nz);

	double dot = nx + ny + nz;
	
    double theta = acos(fabs(cos)>1.0? 1.0 : fabs(cos));
    if (cos < 0 && dot >= 0) // 2nd quarter
    {
		theta = PI - theta;
    }
    else if (cos <= 0 && dot <= 0) // 3rd quarter
    {
		theta = PI + theta;
    }
    else if (cos > 0 && dot < 0) // 4th quarter
    {
		theta = 2 * PI - theta;
    }
    return theta;
}

double GeometricMethod::angle3D(double p0[3], double p1[3], double p2[3])
{
	double ax = p1[0] - p0[0]; double ay = p1[1] - p0[1]; double az = p1[2] - p0[2];
    double bx = p2[0] - p0[0]; double by = p2[1] - p0[1]; double bz = p2[2] - p0[2];

    double coss = (ax * bx + ay * by + az * bz) /(distance(p0, p1)*distance(p0, p2));
	
	coss = coss>1.0 ? 1.0 : coss;
	coss = coss<-1.0 ? -1.0 : coss;
    return acos(coss);
}

void GeometricMethod::planeBisector(double p1Center[3], double p1Prev[3], double p1Next[3], double p2Center[3], double p2Prev[3], double p2Next[3], double coef[4])
{
	int i;
	double v1[3], v2[3];

	//Bisector plane equation: coef[0]*x + coef[1]*y + coef[2]*z + coef[3] = 0

	planeNormal(p1Center, p1Prev, p1Next, v1);
	planeNormal(p2Center, p2Prev, p2Next, v2);
	normalize(v1);
	normalize(v2);

	//For any point P(x, y, z) on the bisector plane, dist(P, plane_1) = dist(P, plane2)
	//Give: v1[0]*x + v1[1]*y + v1[2]*z + v1[3] = v2[0]*x + v2[1]*y + v2[2]*z + v2[3] 
	//Resulting: (v1[0] - v2[0])*x + (v1[1] - v2[1])*y + (v1[2] - v2[2])*z + d = 0
	//Where d = coef[3];

	coef[3] = 0;
	for (i=0; i<3; i++)
	{
		coef[i] = v1[i] - v2[i];
		coef[3] += p2Center[i]*v2[i] - p1Center[i]*v1[i];
	}

	if (norm(coef) == 0.0)
	{
		coef[3] = 0;
		for (i=0; i<3; i++)
		{
			coef[i] = v1[i];
			coef[3] += -(p2Center[i] + p1Center[i])*v1[i]/2.0;
		}
	}
}


double GeometricMethod::angle(double p0[3], double p1[3], double p2[3])
{
	double ax = p1[0] - p0[0]; double ay = p1[1] - p0[1]; double az = p1[2] - p0[2];
    double bx = p2[0] - p0[0]; double by = p2[1] - p0[1]; double bz = p2[2] - p0[2];

    double a = distance(p0, p1);
    double b = distance(p0, p2);
    double coss = (ax * bx + ay * by + az * bz) / (a * b);

	double nx, ny, nz;
	nx = p1[0] - p2[0]; ny = p1[1] - p2[1]; nz = p1[2] - p2[2];
	crossProduct(-bx, -by, -bz, nx, ny, nz, nx, ny, nz);

	double dot = nx + ny + nz;
	
    double theta = acos(fabs(coss)>1.0 ? 1.0 : fabs(coss));
    if (coss < 0 && dot >= 0) // 2nd quarter
    {
		theta = PI - theta;
    }
    else if (coss <= 0 && dot <= 0) // 3rd quarter
    {
		theta = PI + theta;
    }
    else if (coss > 0 && dot < 0) // 4th quarter
    {
		theta = 2 * PI - theta;
    }
    return theta;
}

int GeometricMethod::angle(double rad, int error)
{
	int ang = (int) (rad/PI*180+0.5);
	
	if		(abs(180 - ang) <= error)
		return 180;
	else if	(abs(90 - ang) <= error)
		return 90;
	else if	(abs(270 - ang) <= error)
		return 270;
	else if	(abs(135 - ang) <= error)
		return 135;
	else if	(abs(225 - ang) <= error)
		return 225;
	else if	(abs(45 - ang) <= error)
		return 45;
	else if	(abs(0 - ang) <= error)
		return 0;
	//*/
	
	return ang;
}

void GeometricMethod::angularBisector(double cp[3], double fp[3], double np[3], double newp[3], double dist)
{
	double dDist, fx, fy, nx, ny;

	dDist = distance(cp[0], cp[1], fp[0], fp[1]);
	fx = (fp[0] - cp[0])/dDist; fy = (fp[1] - cp[1])/dDist;
	dDist = distance(cp[0], cp[1], np[0], np[1]);
	nx = (np[0] - cp[0])/dDist; ny = (np[1] - cp[1])/dDist;
	dDist = sqrt((fx - nx)*(fx - nx) + (fy - ny)*(fy - ny));
	newp[0] = cp[0] + dist*(ny - fy) / dDist;
	newp[1] = cp[1] - dist*(nx - fx) / dDist;
	newp[2] = cp[2];

	return;
}

void GeometricMethod::angularBisector(double cx, double cy, double fx, double fy, double nx, double ny, double &xx, double &yy, double dist)
{
	double dDist;

	dDist = distance(cx, cy, fx, fy);
	fx = (fx - cx)/dDist; fy = (fy - cy)/dDist;
	dDist = distance(cx, cy, nx, ny);
	nx = (nx - cx)/dDist; ny = (ny - cy)/dDist;
	dDist = sqrt((fx - nx)*(fx - nx) + (fy - ny)*(fy - ny));
	xx = cx + dist*(ny - fy) / dDist;
	yy = cy - dist*(nx - fx) / dDist;

	return;
}

int GeometricMethod::lineSegmentIntersectingLineSegment(double x11, double y11, double x12, double y12, double x21, double y21, double x22, double y22)
{
	#pragma region Check whether two lines intersect to each other.
	double minX1, minX2, minY1, minY2, maxX1, maxX2, maxY1, maxY2;

	minX1 = x11; maxX1 = x12;
	if (x11 > x12)
		{minX1 = x12; maxX1 = x11;}
	minX2 = x21; maxX2 = x22;
	if (x21 > x22)
		{minX2 = x22; maxX2 = x21;}
	if (maxX1 <  minX2 || minX1 > maxX2)
		return 0;

	minY1 = y11; maxY1 = y12;
	if (y11 > y12)
		{minY1 = y12; maxY1 = y11;}
	minY2 = y21; maxY2 = y22;
	if (y21 > y22)
		{minY2 = y22; maxY2 = y21;}
	if (maxY1 <  minY2 || minY1 > maxY2)
		return 0;

	double a1, a2, a3, a4;

	a1 = (x11 - x21)*(y12 - y21) - (x12 - x21)*(y11 - y21);
	if (fabs(a1) < ERR2)
	{
		a3 = (x21 - x11)*(y22 - y11) - (x22 - x11)*(y21 - y11);
		if (fabs(a3) < ERR2)
			return 1;
		a4 = (x21 - x12)*(y22 - y12) - (x22 - x12)*(y21 - y12);
		if (a3*a4 < ERR2)
			return 1;
		return 0;
	}

	a2 = (x11 - x22)*(y12 - y22) - (x12 - x22)*(y11 - y22);
	if (fabs(a2) < ERR2)
	{
		a3 = (x21 - x11)*(y22 - y11) - (x22 - x11)*(y21 - y11);
		if (fabs(a3) < ERR2)
			return 1;
		a4 = (x21 - x12)*(y22 - y12) - (x22 - x12)*(y21 - y12);
		if (a3*a4 < ERR2)
			return 1;
		return 0;
	}

	if (a1*a2 > ERR2)
		return 0;

	a3 = (x21 - x11)*(y22 - y11) - (x22 - x11)*(y21 - y11);
	if (fabs(a3) < ERR2)
		return 2;
	a4 = (x21 - x12)*(y22 - y12) - (x22 - x12)*(y21 - y12);
	if (a3*a4 <ERR2)
		return 2;

	return 0;

	#pragma endregion
}


void GeometricMethod::intersectingPointOfThreePlane(double coef1[4], double coef2[4], double coef3[4], double point[3])
{
	int i;

	for (i=0; i<4; i++)
	{
		coef1[i] /=coef1[0];
		coef2[i] /=coef2[0];
		coef3[i] /=coef3[0];
	}
	for (i=0; i<4; i++)
	{
		coef2[i] -=coef1[i];
		coef3[i] -=coef1[i];
	}
	for (i=1; i<4; i++)
	{
		coef2[i] /=coef2[1];
		coef3[i] /=coef3[1];
	}
	for (i=1; i<4; i++)
		coef3[i] -=coef2[1];

	point[2] = coef3[3]/coef3[2];
	point[1] = coef2[3] - coef2[2]*point[2];
	point[0] = coef1[3] - coef1[2]*point[2] - coef1[1]*point[1];
}

bool GeometricMethod::lineIntersectingLine(double x11, double y11, double x12, double y12, double x21, double y21, double x22, double y22, double &xx, double &yy)
{
	double k1, k2;

	xx = 100; yy = 100;

	if (fabs(x11 - x12) < 0.01*fabs(y11 - y12))
	{
		if (fabs(x21 - x22) < 0.01*fabs(y21 - y22))
			return false;
		else
		{
			xx = x11;
			yy = (y21 - y22)/(x21 - x22)*(xx - x21) + y21;
			return true;
		}
	}

	if (fabs(x21 - x22) < 0.01*fabs(y21 - y22))
	{
		if (fabs(x11 - x12) < 0.01*fabs(y11 - y12))
			return false;
		else
		{
			xx = x21;
			yy = (y11 - y12)/(x11 - x12)*(xx - x11) + y11;
			return true;
		}
	}
	
	k1 = (y11 - y12)/(x11 - x12);
	k2 = (y21 - y22)/(x21 - x22);

	if (fabs(k1 - k2) < 0.005)
		return false;

	xx = (k2*x21 - y21 - k1*x11 + y11)/(k2 - k1);
	yy = k1*(xx - x11) + y11;
	return true;
}

void GeometricMethod::boundingBox(double p1[3], double p2[3], double boxLow[3], double boxUp[3])
{
	if (p1[0] < p2[0])
		{boxLow[0] = p1[0]; boxUp[0] = p2[0];}
	else
		{boxLow[0] = p2[0]; boxUp[0] = p1[0];}

	if (p1[1] < p2[1])
		{boxLow[1] = p1[1]; boxUp[1] = p2[1];}
	else
		{boxLow[1] = p2[1]; boxUp[1] = p1[1];}

	if (p1[2] < p2[2])
		{boxLow[2] = p1[2]; boxUp[2] = p2[2];}
	else
		{boxLow[2] = p2[2]; boxUp[2] = p1[2];}

	return;
}

void GeometricMethod::boundingBox(double p1[3], double p2[3], double p3[3], double boxLow[3], double boxUp[3])
{
	boundingBox(p1, p2, boxLow, boxUp);

	if		(boxLow[0] > p3[0])
		boxLow[0] = p3[0];
	else if (boxUp[0] < p3[0])
		boxUp[0] = p3[0];

	if		(boxLow[1] > p3[1])
		boxLow[1] = p3[1];
	else if (boxUp[1] < p3[1])
		boxUp[1] = p3[1];

	if		(boxLow[2] > p3[2])
		boxLow[2] = p3[2];
	else if (boxUp[2] < p3[2])
		boxUp[2] = p3[2];

	return;
}

bool GeometricMethod::isTwoBoundingBoxNotIntersected(double box1_Low[3], double box1_Up[3], double box2_Low[3], double box2_Up[3])
{
	if	(box2_Low[0] > box1_Up[0] || box2_Up[0] < box1_Low[0])
		return true;
	if	(box2_Low[1] > box1_Up[1] || box2_Up[1] < box1_Low[1])
		return true;
	if	(box2_Low[2] > box1_Up[2] || box2_Up[2] < box1_Low[2])
		return true;

	return false;

	/*if		(box2_Low[0] > box1_Low[0] && box2_Low[0] < box1_Up[0])
	{
		if		(box2_Low[1] > box1_Low[1] && box2_Low[1] < box1_Up[1])
		{
			if		(box2_Low[2] > box1_Low[2] && box2_Low[2] < box1_Up[2])
				return true;
			else if (box2_Up[2] > box1_Low[2] && box2_Up[2] < box1_Up[2])
				return true;
		}
		else if (box2_Up[1] > box1_Low[1] && box2_Up[1] < box1_Up[1])
		{
			if		(box2_Low[2] > box1_Low[2] && box2_Low[2] < box1_Up[2])
				return true;
			else if (box2_Up[2] > box1_Low[2] && box2_Up[2] < box1_Up[2])
				return true;
		}
	}
	else if (box2_Up[0] > box1_Low[0] && box2_Up[0] < box1_Up[0])
	{
		if		(box2_Low[1] > box1_Low[1] && box2_Low[1] < box1_Up[1])
		{
			if		(box2_Low[2] > box1_Low[2] && box2_Low[2] < box1_Up[2])
				return true;
			else if (box2_Up[2] > box1_Low[2] && box2_Up[2] < box1_Up[2])
				return true;
		}
		else if (box2_Up[1] > box1_Low[1] && box2_Up[1] < box1_Up[1])
		{
			if		(box2_Low[2] > box1_Low[2] && box2_Low[2] < box1_Up[2])
				return true;
			else if (box2_Up[2] > box1_Low[2] && box2_Up[2] < box1_Up[2])
				return true;
		}
	}*/

	return false;
}

int GeometricMethod::lineSegmentIntersectingPlane(double line1[3], double line2[3], double pPrev[3], double pCenter[3], double pNext[3], double pos[3])
{
	int i;
	double pNormal[3], dVec1[3], dVec2[3], dDist, dDist1, dDist2;

	for (i=0; i<3; ++i)
	{
		dVec1[i] = line1[i] - pCenter[i];
		dVec2[i] = line2[i] - pCenter[i];
	}
	planeNormal(pCenter, pPrev, pNext, pNormal);
	normalize(pNormal);
	dDist1 = dotProduct(dVec1, pNormal);
	dDist2 = dotProduct(dVec2, pNormal);

	dDist = dDist1*dDist2;	
	if (fabs(dDist) < ERR)
	{
		if (fabs(dDist1) < ERR)
			{pos[0] = line1[0]; pos[1] = line1[1]; pos[2] = line1[2];}
		else
			{pos[0] = line2[0]; pos[1] = line2[1]; pos[2] = line2[2];}

		return 1;
	}
	else if (dDist > 0)
	{
		dDist = fabs(dDist1) - fabs(dDist2);
		dDist2 = fabs(dDist2)*distance(line1, line2)/dDist;
		pos[0] = line2[0]; pos[1] = line2[1]; pos[2] = line2[2];
		move(pos, line1, dDist2);

		return 0;
	}
	else
	{
		dDist = fabs(dDist1) + fabs(dDist2);
		dDist1 = fabs(dDist1)*distance(line1, line2)/dDist;
		pos[0] = line1[0]; pos[1] = line1[1]; pos[2] = line1[2];
		move(pos, line2, dDist1);

		return 2;
	}
}

int GeometricMethod::lineSegmentIntersectingTriangle(double line1[3], double line2[3], double pPrev[3], double pCenter[3], double pNext[3], double pos[3])
{	
	if (lineSegmentIntersectingPlane(line1, line2, pPrev, pCenter, pNext, pos) == 0)
		return 0;

	double pProject[3]={pos[0], pos[1], pos[2]};
	//Method 2: using angle
	//*
	double ang1, ang2, tAng;
	bool inLine=false;

	if (distance(pProject, pPrev) < ERR || distance(pProject, pCenter) < ERR || distance(pProject, pNext) < ERR)
		return 1; //project to vertex of triangle

	double TOL_1 = 0.035;
	double TOL_2 = TOL_1/2.0;
	tAng = angle3D(pPrev, pCenter, pNext);
	ang1 = angle3D(pPrev, pCenter, pProject);
	ang2 = angle3D(pPrev, pProject, pNext);
	if (fabs(tAng - ang1 - ang2) > TOL_1)
		return 0;
	if (fabs(tAng - ang1) <= TOL_2 || fabs(tAng - ang2) <= TOL_2)
		inLine = true;

	tAng = angle3D(pCenter, pPrev, pNext);
	ang1 = angle3D(pCenter, pPrev, pProject);
	ang2 = angle3D(pCenter, pProject ,pNext);
	if (fabs(tAng - ang1 - ang2) > TOL_1)
		return 0;
	if (fabs(tAng - ang1) <= TOL_2 || fabs(tAng - ang2) <= TOL_2)
		inLine = true;

	return inLine? 2 : 3;	//2: project to line; 3: project inside
	//*/
	//End of Method 2
}

inline double GeometricMethod::projectToPlane(double point[3], double pPrev[3], double pCenter[3], double pNext[3], double pProject[3])
{
	int i;
	double dNormal[3], dVec[3], dDist;

	for (i=0; i<3; ++i)
		{dVec[i] = point[i] - pCenter[i];}
	planeNormal(pCenter, pPrev, pNext, dNormal);
	normalize(dNormal);
	dDist = dotProduct(dVec, dNormal);
	for (i=0; i<3; ++i)
		pProject[i] = point[i] - dNormal[i] * dDist;

	return dDist;
}

int GeometricMethod::projectToTriangle(double point[3], double pPrev[3], double pCenter[3], double pNext[3], double pProject[3], double &pDist)
{
	pDist = projectToPlane(point, pPrev, pCenter, pNext, pProject);
	
	//Method 2: using angle
	//*
	double ang1, ang2, tAng;
	bool inLine=false;

	if (distance(pProject, pPrev) < ERR || distance(pProject, pCenter) < ERR || distance(pProject, pNext) < ERR)
		return 1; //project to vertex of triangle

	tAng = angle3D(pPrev, pCenter, pNext);
	ang1 = angle3D(pPrev, pCenter, pProject);
	ang2 = angle3D(pPrev, pNext, pProject);
	if (fabs(tAng - ang1 - ang2) > 1.0e-4)
		return 0;
	if (fabs(tAng - ang1) <= 1.0e-4 || fabs(tAng - ang2) <= 1.0e-4)
		inLine = true;

	tAng = angle3D(pCenter, pPrev, pNext);
	ang1 = angle3D(pCenter, pPrev, pProject);
	ang2 = angle3D(pCenter, pNext, pProject);
	if (fabs(tAng - ang1 - ang2) > 1.0e-4)
		return 0;
	if (fabs(tAng - ang1) <= 1.0e-4 || fabs(tAng - ang2) <= 1.0e-4)
		inLine = true;

	return inLine? 2 : 3;	//2: project to line; 3: project inside
	//*/
	//End of Method 2
}

void GeometricMethod::sort(double * &matrix, int const dim)
{
	int i, j;
	double dTmp;
	for (i=0; i<dim; i++)
		for (j=i+1; j<dim; j++)
			if (matrix[i] > matrix[j])
			{
				dTmp = matrix[i];
				matrix[i] = matrix[j];
				matrix[j] = dTmp;
			}
	return;
}

void GeometricMethod::QuickSort(double* &x, int* &I, int left, int right)
{
	#pragma region Quick Sort - Arrange the indices in I according to the values in x:
	int pivot, l_hold, r_hold;
	l_hold = left;
	r_hold = right;
	double xpivot = x[left];
	int Ipivot = I[left];
	while(left < right)
	{
		while(x[right] >= xpivot && left < right)
			right--;

		if (left != right)
		{
			x[left] = x[right];
			I[left] = I[right];
			left++;
		}

		while (x[left] <= xpivot && left < right)
			left++;

		if (left != right)
		{
			x[right] = x[left];
			I[right] = I[left];
			right--;
		}
	}
	x[left] = xpivot;
	I[left] = Ipivot;

	pivot = left;
	left = l_hold;
	right = r_hold;
	if (left < pivot)
		QuickSort(x, I, left, pivot - 1);
	if (right > pivot)
		QuickSort(x, I, pivot + 1, right);	
	#pragma endregion
}

void GeometricMethod::QuickSort(int* &I1, int* &I2, int* &I3, int* &I4, int* &II, int left, int right)
{
	#pragma region Quick Sort - Arrange the indices in I1, I2, I3, I4, II according to the values in I1:
	int pivot, l_hold, r_hold;
	l_hold = left;
	r_hold = right;
	int I1pivot = I1[left];
	int I2pivot = I2[left];
	int I3pivot = I3[left];
	int I4pivot = I4[left];
	int IIpivot = II[left];	
	while(left < right)
	{
		while(I1[right] >= I1pivot && left < right)
			right--;

		if (left != right)
		{
			I1[left] = I1[right];
			I2[left] = I2[right];
			I3[left] = I3[right];
			I4[left] = I4[right];
			II[left] = II[right];			
			left++;
		}

		while (I1[left] <= I1pivot && left < right)
			left++;

		if (left != right)
		{
			I1[right] = I1[left];
			I2[right] = I2[left];
			I3[right] = I3[left];
			I4[right] = I4[left];
			II[right] = II[left];
			right--;
		}
	}
	I1[left] = I1pivot;
	I2[left] = I2pivot;
	I3[left] = I3pivot;
	I4[left] = I4pivot;
	II[left] = IIpivot;	

	pivot = left;
	left = l_hold;
	right = r_hold;
	if (left < pivot)
		QuickSort(I1, I2, I3, I4, II, left, pivot - 1);
	if (right > pivot)
		QuickSort(I1, I2, I3, I4, II, pivot + 1, right);	
	#pragma endregion
}

void GeometricMethod::QuickSort(int* &I1, int* &I2, int* &II, int left, int right)
{
	#pragma region Quick Sort - Arrange the indices in I1, I2, II according to the values in I1:
	int pivot, l_hold, r_hold;
	l_hold = left;
	r_hold = right;
	int I1pivot = I1[left];
	int I2pivot = I2[left];	
	int IIpivot = II[left];	
	while(left < right)
	{
		while(I1[right] >= I1pivot && left < right)
			right--;

		if (left != right)
		{
			I1[left] = I1[right];
			I2[left] = I2[right];			
			II[left] = II[right];			
			left++;
		}

		while (I1[left] <= I1pivot && left < right)
			left++;

		if (left != right)
		{
			I1[right] = I1[left];
			I2[right] = I2[left];			
			II[right] = II[left];
			right--;
		}
	}
	I1[left] = I1pivot;
	I2[left] = I2pivot;	
	II[left] = IIpivot;	

	pivot = left;
	left = l_hold;
	right = r_hold;
	if (left < pivot)
		QuickSort(I1, I2, II, left, pivot - 1);
	if (right > pivot)
		QuickSort(I1, I2, II, pivot + 1, right);	
	#pragma endregion
}


void GeometricMethod::perpendicular(double x1, double y1, double x2, double y2, double &xx, double &yy)
{
	if (fabs(x1 - x2) <= ERR)
	{
		if (y2 > y1)
		{
			xx = x1 -1;
			yy = y1;
		}
		else
		{
			xx = x1 + 1;
			yy = y1;
		}
	}
	else if (fabs(y1 - y2) <= ERR)
	{
		if (x2 > x1)
		{
			xx = x1;
			yy = y1 + 1;
		}
		else
		{
			xx = x1;
			yy = y1 -1;
		}
	}
	else
	{
		double kk;
		kk = (y2 - y1)/(x2 - x1);
		if (kk > 0)
		{
			if (x2 > x1)
			{
				xx = x1 - 1.0;
				yy = y1 + 1.0/kk;
			}
			else
			{
				xx = x1 + 1.0;
				yy = y1 - 1.0/kk;
			}
		}
		else
		{
			if (x2 > x1)
			{
				xx = x1 + 1.0;
				yy = y1 - 1.0/kk;
			}
			else
			{
				xx = x1 - 1.0;
				yy = y1 + 1.0/kk;
			}
		}
	}

	return;
}


void GeometricMethod::rotate(double x1, double y1, double x2, double y2, double angle, double &xx, double &yy)
{
	//reduce angle within (0, 360)
	while (angle >= 360.0)
		angle -= 360.0;
	while (angle < 0.0)
		angle += 360.0;

	//*
	angle *= A2R;

	x2 = x2 - x1;
	y2 = y2 - y1;

	xx = x2*cos(angle) - y2*sin(angle);
	yy = x2*sin(angle) + y2*cos(angle);

	xx += x1;
	yy += y1;

	return;
	//*/
	//End of Method 1


	//Method 2: use algebra
	//Not recommended!
	/*
	//Angle = 0
	if (fabs(angle) < 0.001)
	{
		xx = x2;
		yy = y2;
		return;
	}

	//Angle = 180
	if (fabs(angle - 180.0) < 0.001)
	{
		xx = x1 - (x2 - x1);
		yy = y1 - (y2 - y1);
		return;
	}

	//Angle = 90, 270
	if (fabs(angle - 90.0) < 0.001)
	{
		perpendicular(x1, y1, x2, y2, xx, yy);
		return;
	}
	if (fabs(angle - 270.0) < 0.001)
	{
		perpendicular(x1, y1, x2, y2, xx, yy);
		xx = x1 - (xx - x1);
		yy = y1 - (yy - y1);
		return;
	}

	//else
	double k1, kk;

	if (fabs(x1 - x2) <= ERR)
	{
		if (y2 > y1)
		{
			kk = tan((angle + 90)*A2R);
			xx = x1 - 1;
			yy = -kk + y1;
		}
		else
		{
			kk = tan((angle - 90)*A2R);
			xx = x1 + 1;
			yy = kk + y1;
		}
	}
	else
	{
		k1 = (y1 - y2)/(x1 - x2);
		if ((y2 - y1)*k1 > 0.0)
			angle = normalizeAngle(atan(k1) + angle*A2R);
		else
			angle = normalizeAngle(atan(k1) + (angle + 180)*A2R);
		//if (useRadianInAngle)
		//	angle *= R2A;

		if (fabs(angle - 90) < 0.001)
		{
			xx = x1;
			yy = y1 + 1;
		}
		else if (fabs(angle - 270) < 0.001)
		{
			xx = x1;
			yy = y1 - 1;
		}
		else if (angle < 90 || angle > 270)
		{
			xx = x1 + 1;
			yy = tan(angle*A2R) + y1;
		}
		else
		{
			xx = x1 - 1;
			yy = y1 - tan(angle*A2R);
		}
	}
	//*/
}
