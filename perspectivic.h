#ifndef __PERS__
#define __PERS__

#include "GUI.h"
#include "vecmat.h"

#include <iostream>
#include <vector>

using namespace std;

class Perspectivic {
private:
	struct Point {
		float x, y, z;
		Point(float _x, float _y, float _z) : x(_x), y(_y), z(_z) {}
	};

	struct Color {
		int R, G, B;
		Color(int _R, int _G, int _B) : R(_R), G(_G), B(_B) {}
	};

	struct Segment {
		Point begin, end;
		Color color;
		Segment(Point _begin, Point _end, Color _color) : begin(_begin), end(_end), color(_color) {}
	};

	vector<Perspectivic::Segment> data;

	Matrix4 t;
	Matrix4 t1;

	double maxY = -1.0e30;
	double minY = 1.0e30;

public:
	double countFunction(double x, double y);
	void RecountFunctionIntoData(vector<vector<double>> funValues);

	Matrix4 XRotation(double alpha);
	Matrix4 YRotation(double alpha);
	Matrix4 ZRotation(double alpha);

	void GenerateTransformMatrix(wxPanel* drawingPanel);

	void getMinYMaxY();

	void Repaint(wxPanel* drawingPanel);

};

#endif