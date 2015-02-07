/*
 * typedefs.c
 * Helper functions for working with polygons and points
 *
 * Copyright (c) 2009, Paul Chote
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice, this
 * list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 * this list of conditions and the following disclaimer in the documentation
 * and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
 * ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#include <math.h>
#include <stdio.h>
#include "typedefs.h"
#include <cpgplot.h>

extern boolean debugMode;

/*
 * Creates a point with given parameters.
 */
point makePoint(double x, double y)
{
	point p;
	p.x = x;
	p.y = y;
	return p;
}

/*
 * returns TRUE if p is in the area a, false if it is not.
 */
boolean pointInArea(point p, searchArea a)
{
	return (p.x >= a.x && p.x <= (a.x + a.size) && p.y >= a.y && p.y <= (a.y + a.size)) ? TRUE : FALSE;
}

/*
 * Returns the point (ratio) distance along the line between two given points
 * ratio is a double between [0,1]. 0 = startPoint, 1 = endPoint
 */
point interpolatePosition(point startPoint, point endPoint, double ratio)
{
	point newPoint;
	newPoint.x = startPoint.x + (endPoint.x - startPoint.x)*ratio;
	newPoint.y = startPoint.y + (endPoint.y - startPoint.y)*ratio;
	return newPoint;
}

/*
 * Returns a point at the position of the requested corner.
 */
point areaCorner(searchArea a, corner c)
{
	switch (c) {
		default:
		case BOTTOM_LEFT:
			return makePoint(a.x, a.y);
		case TOP_LEFT:
			return makePoint(a.x, a.y+a.size);
		case TOP_RIGHT:
			return makePoint(a.x+a.size, a.y+a.size);
		case BOTTOM_RIGHT:
			return makePoint(a.x+a.size, a.y);
	}
}

/*
 * Creates a searchArea with given parameters.
 */
searchArea makeSearchArea(double x, double y, double size)
{
	searchArea s;
	s.x = x;
	s.y = y;
	s.size = size;
	return s;
}

/*
 * Creates a searchGrid with given parameters.
 */
searchGrid makeSearchGrid(searchArea a, source *source, event *event, boolean checkLenses, boolean checkCriticalCurve, int level)
{
	searchGrid s;
	s.searchArea = a;
	s.source = source;
	s.event = event;
	s.checkLenses = checkLenses;
	s.checkCriticalCurve = checkCriticalCurve;
	s.level = level;
	return s;
}

/*
 * Creates a source with given parameters.
 */
source makeSource(point origin, double radius)
{
	source s;
	s.origin = origin;
	s.radius = radius;
	return s;
}

/*
 * Creates a lens with given parameters.
 */
lens makeLens(point origin, double mass)
{
	lens l;
	l.origin = origin;
	l.mass = mass;
	return l;
}

/*
 * Creates an event with given parameters.
 */
event makeEvent(int numLenses, lens *lenses, double resolution)
{
	event e;
	e.numLenses = numLenses;
	e.lenses = lenses;
	e.resolution = resolution;
	return e;
}

/*
 * Find the contribution of a given lens to the jacobian determinant
 * of the lens equation at a point.
 */
double lensJacobianContribution(lens l, point p, int type)
{
	double dx = p.x - l.origin.x;
	double dy = p.y - l.origin.y;
	double dsq = dx * dx + dy * dy;
	
	switch (type) {
		default:
		case 0:
			return -l.mass / dsq + 2 * l.mass * dx * dx / (dsq * dsq);
		case 1:
			return 2 * l.mass * dx * dy / (dsq * dsq);
		case 2:
			return -l.mass / dsq + 2 * l.mass * dy * dy / (dsq * dsq);
	}
}

/*
 * Determines whether the source is inside, outside, or near the boundary of a given polygon
 */
intersectionType newTestPolygonAgainstSource(point *v, int vC, source *source)
{
	boolean totallyInsideSource = TRUE;
	boolean partiallyInsideSource = FALSE;
	int windingNumber = 0;
	
	int i;
	for(i=0;i<vC;i++)
	{
		/*
		 * New Algorithm:
		 * Check each vertex for three things:
		 *	1) Does it lie inside the source?
		 *		If yes, and at least one other vertex lies outside the source, 
		 *		then the polygon overlaps the source disk and we return OVERLAP
		 *
		 *	2) Is the vertex "close" to the source?
		 *		By close, do we mean "is there any chance that the source disk 
		 *		intersects the polygon on an edge between two vertices?".
		 *		Take a circular region around each vertex of radius equal to the 
		 *		distance to the next vertex. If the source disk intersects this 
		 *		boundary, we should subdivide the polygon and check again.
		 *
		 *	3) Does the source lie totally inside the polygon?
		 *		We check this by calculating the "winding number" of the point wrt
		 *		the polygon. For each edge of the polygon that lies on the right 
		 *		of the centre of the source, increment the winding number if the
		 *		edge points upwards, decrement it if it points downward.
		 *		If the winding number is non-zero, then the source is inside the polygon.
		 */
		
		point tv = v[i]; // This vertex
		point nv = v[(i+1)%vC]; // Next vertex
		point sp = source->origin; // Source position
		
		if (lineOnRightOfPoint( tv, nv, sp ))
		{
			if (tv.y > sp.y && nv.y <= sp.y) // Downward crossing
				windingNumber--;
			else if (tv.y <= sp.y && nv.y > sp.y) // Upward crossing
				windingNumber++;
		}
		
		// Find distance between vertex and source
		double vd = hypot(v[i].x-source->origin.x,v[i].y-source->origin.y);

		// Is vertex inside source?
		if (vd <= source->radius)
		{
			partiallyInsideSource = TRUE;
			if (!totallyInsideSource)
			{
				return OVERLAP;
			}
			continue;
		}
		else
		{
			totallyInsideSource = FALSE;
			if (partiallyInsideSource)
			{
				return OVERLAP;
			}
		}
		
		// Is source near the edge of the polygon?
	
		// Distance to the next vertex
		double nvd = hypot(v[i].x - v[(i+1)%vC].x,v[i].y - v[(i+1)%vC].y);

		// Does the source disk intersect this disk?
		if (vd <= source->radius+nvd)
		{
			//printf("nearby\n");
			return OVERLAP;
		}
	}
	
	// If all vertices are inside source, then so is polygon
	if (totallyInsideSource)			{
		return INSIDE_SOURCE;
	}

	// Is source inside or outside the polygon?
	if (windingNumber != 0)
	{
		return ENCLOSES_SOURCE;
	}	
	return NO_OVERLAP;
}

/*
 * Check whether a line (defined by v0 and v1) lies to the left or right of a point 
 * If the line is horizontal, define the point to be left if it is below or on the edge, right if it is above
 */
boolean lineOnRightOfPoint(point v0, point v1, point p)
{
	double dx = (p.y-v0.y)*(v1.x-v0.x)/(v1.y-v0.y)+(v0.x-p.x);
	
	if (dx == 0)
		return (p.y <= v0.y);
	else
		return (dx > 0);
}

/*
 * Find where a given source disk lies in relation to a given polygon
 */
intersectionType testPolygonAgainstSource(point *v, int vp, source *source)
{
	
	// for each vertex, check if it is inside the source.	
	boolean polygonEnclosedInSource = TRUE;
	boolean hit = FALSE;
	int edgeHits = 0;
	
	// for each vertex, check if it is inside the source.
	int i;
	for(i=0;i<vp;i++)
	{
		// calculate distance between vertex and centre of source object
		double vd =  hypot(v[i].x-source->origin.x,v[i].y-source->origin.y);
		
		// is vertex inside source object?
		if (vd <= source->radius)
			hit = TRUE;
		else
			polygonEnclosedInSource = FALSE;
	}
	
	// If all vertices are inside source, then so is polygon
	if (polygonEnclosedInSource) return INSIDE_SOURCE;
	// otherwise, if at lest one is inside, the grid partially overlaps the source
	if (hit) return OVERLAP;
	
	// Need to check if source is inside grid, or crosses an edge
	// group the vertices into line segments - each vertex is a part of 2 lines
	for (i=0;i<vp;i++)
	{
		// coordinates for the line end-points
		double u1 = v[i].x;
		double v1 = v[i].y;
		
		// loop the index back to zero once it reaches the end of the array so that we 
		//  can check the final line segment between the last point and the first point
		double u2 = v[(i+1) % vp].x;
		double v2 = v[(i+1) % vp].y;
		
		// commonly used expressions
		double du = u2 - u1;
		double dv = v2 - v1;
		double dux = u1 - source->origin.x;
		double dvy = v1 - source->origin.y;
		
		// check if ray from (x, +inf) to (x, y) hits this edge (the strict less than in x1 is to ensure that if a vertex is exactly source, it only gets counted once)
		double p = -dux/du;
		if (0 < p && p <= 1 && (v1 + p*dv >= source->origin.y))
			edgeHits++;
		
		// check if edge intersects source
		
		// coefficients of quadratic equation in p (the parameter describing where on the line segment the line crosses the source
		double a = du*du+dv*dv;
		double b = 2*(du*dux+dv*dvy);
		double c = dux*dux+dvy*dvy-source->radius*source->radius;
		
		// discriminant
		double d = b*b-4*a*c;
		
		// no solution if discriminant is negative
		if (d < 0)
			continue;
		
		// find solutions for p and check if it is on the line segment (0 <= p <=1)
		double p1 = (-b + sqrt(d))/(2*a);
		double p2 = (-b - sqrt(d))/(2*a);
		
		// If an edge is hit, there is an overlap; no need to continue
		if (( 0 <= p1 && p1 <= 1) || (0 <= p2 && p2 <= 1))
			return OVERLAP;
	}
	
	// e encloses source if there is no overlap and the test ray intersects an odd number of edges
	if ((edgeHits % 2) && !hit)
	{
		return ENCLOSES_SOURCE;
	}
	//if (!hit && pointInPolygon(source->origin, v, vp))
	//	return ENCLOSES_SOURCE;
	
	if (hit) return OVERLAP;
	
	return NO_OVERLAP;
}
