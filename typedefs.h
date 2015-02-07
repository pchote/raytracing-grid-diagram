/*
 * typedefs.h
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

#ifndef TYPEDEFS_HEADER
#define TYPEDEFS_HEADER

typedef char boolean;
#define TRUE 1
#define FALSE 0
#define PI 3.14159265

typedef enum intersectionType {
	NO_OVERLAP = 3,
	INSIDE_SOURCE = 4,
	ENCLOSES_SOURCE = 5, 
	OVERLAP = 6
} intersectionType;

typedef enum corner {
	BOTTOM_LEFT = 0,
	TOP_LEFT = 1,
	TOP_RIGHT = 2,
	BOTTOM_RIGHT = 3
} corner;

typedef struct point {
	double x;
	double y;
} point;

typedef struct searchArea {
	double x;
	double y;
	double size;
} searchArea;

typedef struct lens {
	point origin;
	double mass;
} lens;

typedef struct source {
	point origin;
	double radius;
} source;

typedef struct event {
	int numLenses;
	lens *lenses;
	double resolution;
} event;

typedef struct searchGrid {
	event *event;
	source *source;
	searchArea searchArea;
	boolean checkLenses;
	boolean checkCriticalCurve;
	int level;
} searchGrid;

point makePoint(double x, double y);
boolean pointInArea(point p, searchArea a);
point interpolatePosition(point startPoint, point endPoint, double ratio);
point areaCorner(searchArea a, corner c);
searchArea makeSearchArea(double x, double y, double size);
searchGrid makeSearchGrid(searchArea a, source *source, event *event, boolean checkLenses, boolean checkCriticalCurve, int level);
source makeSource(point origin, double radius);
lens makeLens(point origin, double mass);
event makeEvent(int numLenses, lens *lenses, double resolution);
double lensJacobianContribution(lens l, point p, int type);
intersectionType testPolygonAgainstSource(point *v, int vp, source *source);
intersectionType testPolygonAgainstSourceWithLogging(point *v, int vp, source *source);
boolean pointInPolygon(point p, point *v, int vC);
boolean lineOnRightOfPoint( point tv, point nv, point sv );
#endif