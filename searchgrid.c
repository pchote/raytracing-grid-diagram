/*
 * searchgrid.c
 * Helper functions for working with the search grid.
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

#include <stdio.h>
#include <math.h>
#include "typedefs.h"
#include "searchgrid.h"
#include <cpgplot.h>


extern boolean debugMode;
extern double totalArea;
extern double eliminated[];
extern int calculations[];

/*
 * Finds and draws to screen the images in a given area
 */
void search(searchGrid grid)
{	
	if (debugMode)
	{
		// Draw Image Plane Grid
		cpgsci(1); //yellow
		cpgsfs(2); //outline
		cpgrect((float)grid.searchArea.x, (float)(grid.searchArea.x+grid.searchArea.size), (float)grid.searchArea.y, (float)(grid.searchArea.y+grid.searchArea.size));
		cpgsfs(1); // fill	
	}

	/*
	 * Check for divergences
	 */
	if (grid.checkLenses)
	{
		int i;
		for (i=0; i < grid.event->numLenses; i++)
		{
			if (pointInArea(grid.event->lenses[i].origin, grid.searchArea))
			{
				if (grid.searchArea.size > grid.event->resolution)
					divideAndConquer(grid);
				return;
			}
		}
		grid.checkLenses = FALSE;
	}
	
	if (grid.checkCriticalCurve)
	{
		// Check if the grid straddles a critical curve
		// The sign of the jacobian determinant changes as you cross a critical curve

		int jacobianSign = jacobianSignAtPoint(areaCorner(grid.searchArea,BOTTOM_LEFT), grid);
		if (jacobianSign != jacobianSignAtPoint(areaCorner(grid.searchArea,TOP_LEFT), grid) ||
			jacobianSign != jacobianSignAtPoint(areaCorner(grid.searchArea,TOP_RIGHT), grid) ||
			jacobianSign != jacobianSignAtPoint(areaCorner(grid.searchArea,BOTTOM_RIGHT), grid))
		{
			if (grid.searchArea.size > grid.event->resolution)
				divideAndConquer(grid);
			
			return;
		}
		
		grid.checkCriticalCurve = FALSE;
	}
	
	// check for hit
	intersectionType hit = mapsToSource(grid);
	calculations[grid.level] += 40;
	
	if (hit == NO_OVERLAP)
	{
		cpgsci(2); // white
		eliminated[grid.level] += grid.searchArea.size*grid.searchArea.size;

		if (debugMode && grid.level < 6)
		{
			char buf[2];
			sprintf(buf, "%d",grid.level);
			cpgtext(grid.searchArea.x+grid.searchArea.size/2-0.04,grid.searchArea.y+grid.searchArea.size/2-0.04, buf);
		}
		return;
	}
	
	if (hit == INSIDE_SOURCE || grid.searchArea.size <= grid.event->resolution)
	{
		cpgsci(1); // white
		eliminated[grid.level] += grid.searchArea.size*grid.searchArea.size;
		cpgrect(grid.searchArea.x, grid.searchArea.x + grid.searchArea.size, grid.searchArea.y, grid.searchArea.y + grid.searchArea.size);
		return;
	}	
	divideAndConquer(grid);
}

/*
 * Split the area into quadrants and continue searching
 */
void divideAndConquer(searchGrid grid)
{		
	searchArea pArea = grid.searchArea;
		
	// Divide search area into 4
	double newSize = pArea.size/2;

	searchArea cRects[4];
	cRects[0] = makeSearchArea(pArea.x, pArea.y, newSize);
	cRects[1] = makeSearchArea(pArea.x+newSize, pArea.y, newSize);
	cRects[2] = makeSearchArea(pArea.x, pArea.y+newSize, newSize);
	cRects[3] = makeSearchArea(pArea.x+newSize, pArea.y+newSize, newSize);
	
	int i;
	for (i=0; i < 4; i++)
	{
		search(makeSearchGrid(cRects[i], grid.source, grid.event, grid.checkLenses, grid.checkCriticalCurve, grid.level+1));
	}
}

/*
 * Returns +/- 1 depending on the sign of the lens equation jacobian at a given point
 */
int jacobianSignAtPoint(point p, searchGrid grid)
{
	double dFxx = 1;
	double dFyy = 1;
	double dFxy = 0;
	
	int i;
	for (i=0; i < grid.event->numLenses; i++)
	{
		dFxx += lensJacobianContribution(grid.event->lenses[i], p, 0);
		dFxy += lensJacobianContribution(grid.event->lenses[i], p, 1);
		dFyy += lensJacobianContribution(grid.event->lenses[i], p, 2);
	}
	
	double jacobian = dFxx*dFyy-dFxy*dFxy;
	return (jacobian > 0) ? 1 : -1;
}

/*
 * Transforms the search area into the source plane and finds how it intersects the source
 */
intersectionType mapsToSource(searchGrid grid)
{
	// Generate a list of points around the edge of the grid to be transformed
	int minPoints = 10;
	int pointsPerSide = ((int)(grid.searchArea.size/grid.event->resolution) > minPoints) ? (int)(grid.searchArea.size/grid.event->resolution) : minPoints;
	int vC = pointsPerSide*4;
	int curV = 0;
	double du = grid.searchArea.size/pointsPerSide;
	int i,j;

	
	// Create an array of points around the edge of the search area
	point v[vC];
	
	// left
	for (i=0;i<pointsPerSide;i++)
		v[curV++] = makePoint(grid.searchArea.x, grid.searchArea.y + i*du);
	
	// top
	for (i=0;i<pointsPerSide;i++)
		v[curV++] = makePoint(grid.searchArea.x + i*du, grid.searchArea.y + grid.searchArea.size);
	
	// right
	for (i=0;i<pointsPerSide;i++)
		v[curV++] = makePoint(grid.searchArea.x + grid.searchArea.size, grid.searchArea.y + grid.searchArea.size - i*du);
	
	// bottom
	for (i=0;i<pointsPerSide;i++)
		v[curV++] = makePoint(grid.searchArea.x + grid.searchArea.size - i*du, grid.searchArea.y);
	
	// Transform the points into the source plane
	point vt[vC];
	for(i=0;i<vC;i++) // for each point
	{
		vt[i] = v[i];
		
		double dx, dy, lensDsq;
		for (j=0; j < grid.event->numLenses; j++) // for each lens
		{
			lens l = grid.event->lenses[j];
			dx = (v[i].x-l.origin.x);
			dy = (v[i].y-l.origin.y);
			lensDsq = dx*dx + dy*dy;
			
			vt[i].x -= l.mass*dx/lensDsq;
			vt[i].y -= l.mass*dy/lensDsq;
		}
	}
	
	intersectionType hit = testPolygonAgainstSource(vt, vC, grid.source); // Test transformed area against source
	
	return hit;
}
