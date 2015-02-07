/*
 * microlensing.c
 * Demonstrates a simple inverse-raytracing microlensing calculation technique.
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
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <cpgplot.h>

#include "typedefs.h"
#include "searchgrid.h"

#define MAX_LIGHTCURVE_POINTS 3000
#include <gsl/gsl_poly.h>
#include <float.h>

boolean debugMode = FALSE;

#define LEVELS 1000
double eliminated[1000];
int calculations[1000];

int main(int argc, char **argv)
{	
	/*
	 * Define event parameters
	 */
	const char *eventName = "Test";
	int numLenses = 2;
	lens lenses[numLenses];
	lenses[0] = makeLens(makePoint(0, 0), 1.0/1.5);
	lenses[1] = makeLens(makePoint(2,0), 0.5/1.5);
	
	double startTime = 5700;//4400;
	double endTime = 6000;//7000;
	double crossingTime = 800;
	double peakTime = 4500;
	double impactRadius = -0.17f;
	double peakMagnification = 40;
	double sourceRadius = 0.05;
	double searchResolution = 1e-2;
	double limbcoefficient = 0;
	double residuals = 0;
	
	int animationFrames = 100;
	double windowX = -1;
	double windowY = -2;
	double windowW = 4;
	
	searchArea a = makeSearchArea(windowX,windowY,windowW);
	event e = makeEvent(numLenses, lenses, searchResolution);
	point startPoint = makePoint((startTime - peakTime)/crossingTime, impactRadius);
	point endPoint = makePoint((endTime - peakTime)/crossingTime, impactRadius);
	source s = makeSource(startPoint, sourceRadius);


	/*
	 * Load caustic and critical curve data from Gravlens
	 */
		char dataLine[100];
	
	#define MAX_CAUSTIC_POINTS 10000
	float causticX[MAX_CAUSTIC_POINTS][2];
	float causticY[MAX_CAUSTIC_POINTS][2];
	float criticalX[MAX_CAUSTIC_POINTS][2];
	float criticalY[MAX_CAUSTIC_POINTS][2];
	int totalCausticPoints = 0;
	
	FILE *curveData = fopen("gravlens.curves", "r");
	
	if(curveData == NULL) {
		printf("Error: cannot open gravlens.curves\n");
		return EXIT_FAILURE;
	}
	
	while(fgets(dataLine, 100, curveData) != NULL)
	{
		if (totalCausticPoints >= MAX_CAUSTIC_POINTS)
		{
			printf("ERROR: Caustic data overruns point arrays");
			break;
		}
		sscanf(dataLine, "%f %f %f %f %f %f %f %f\n", &criticalX[totalCausticPoints][0], &criticalY[totalCausticPoints][0], &causticX[totalCausticPoints][0], &causticY[totalCausticPoints][0], &criticalX[totalCausticPoints][1], &criticalY[totalCausticPoints][1], &causticX[totalCausticPoints][1], &causticY[totalCausticPoints][1]);
		totalCausticPoints++;
	}
	fclose(curveData);
	
	/*
	 * Image plane window setup
	 */
	//int IPWindow = cpgopen("foo.ps/cps");
	int IPWindow = cpgopen("9/xs");
	if (IPWindow <= 0) return EXIT_FAILURE;
	
	cpgslw(4);
	cpgask(FALSE);
	cpgsvp(0.2, 0.8, 0.2, 0.8);
	
	// Set the window area to the search area above
	cpgwnad((float)a.x, (float)(a.x+a.size), (float)a.y, (float)(a.y+a.size));

	int i = 0, j;
	float x,y;
	char c;
	cpgslct(IPWindow);
	
	int jj;
	for (jj = 0; jj < LEVELS; jj++)
	{
		eliminated[jj] == 0;
	}
	
	do
	{
		switch (c)
		{
			case 'f':
				if (i <= animationFrames) i++;
			break;
			case 'b':
				if (i > 0) i--;
			break;
			case 'q':
				goto endloop;
			break;
			case 'g':
				debugMode = !debugMode; 
				break;
			default:
				break;
		}
		
		/*
		 * Image Plane
		 */		
		cpgbbuf(); // Start buffering output
		
		cpgsci(0); // Black
		cpgrect((float)a.x, (float)(a.x+a.size), (float)a.y, (float)(a.y+a.size)); // Erase display
		
		cpgsci(4); // Red
		cpgsfs(2); //outline
		
		if (!debugMode) {
			int ii;
			// cpgcirc(0, 0, 1); // Draw Einstein ring
			for (ii = 0; ii < totalCausticPoints; ii++) {
				cpgline(2, criticalX[ii], criticalY[ii]);
			}

			cpgsci(9);
			for (ii = 0; ii < totalCausticPoints; ii++) {
				cpgline(2, causticX[ii], causticY[ii]);
			}
		}

		cpgsfs(1); // fill
		
		// Place source, and find images (drawing is done inside search();)
		if (animationFrames > 0)
			s.origin = interpolatePosition(startPoint, endPoint, i/(double)animationFrames);
		
		cpgsci(2); // Red
		cpgcirc((float)s.origin.x, (float)s.origin.y, (float)s.radius); // Draw Source disk
		
		clock_t startT = clock();
		clock_t analyticT;
		clock_t numericT;
		
		search(makeSearchGrid(a, &s, &e, TRUE, TRUE, 1));
		numericT = clock()-startT;
		startT = clock();

		// Draw lenses
		cpgsci(8); // Yellow
		for (j=0; j < e.numLenses; j++) // for each lens
		{
			lens l = e.lenses[j];
		//	cpgcirc(l.origin.x, l.origin.y, windowW/200);
		}

		cpgsci(1);
		cpgbox("bcn", 0.0, 0, "bcvn", 0.0, 0); // Plot axes
		cpglab("x (R\\dE\\u)", "y (R\\dE\\u)", "");
		cpgebuf(); // Draw buffer to screen

		cpgslct(IPWindow);
		
		double totalArea = 0;
		int totalCalc = 0;
		for (jj = 1; jj < 20; jj++)
		{
			totalArea += eliminated[jj]/16.0;
			totalCalc += calculations[jj];
			//printf("%d %.6f %d\n", jj, totalArea, totalCalc);
		}
	} while (cpgcurs(&x,&y,&c));
endloop:;
	
	
	//time(&end);
	cpgend();
	//printf("runTime:%f",difftime(end,start));
	
	return EXIT_SUCCESS;
}
