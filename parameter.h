#pragma once
#include "h3api.h"

class Parameter
{
public:
	Parameter() {};
	Parameter(H3Index index, double comprehensive = 0, double grad = 0, int grid_num = 1) :
		index(index), total_grad(grad), total_comprehensive(comprehensive), grid_num(grid_num) {};
	~Parameter() {};
	H3Index index = 0;
	double total_grad = 0;
	double total_comprehensive = 0;
	double distance = 0;
	double avg_grad = 0;
	double avg_comprehensive = 0;
	int grid_num = 0;
	double g1;
	double g2;

	H3Index getIndex() {
		return index;
	}

	void setIndex(H3Index index) {
		this->index = index;
	}
	// ±È½Ïº¯Êý
	bool operator==(Parameter& t) {
		if (index == t.index)
			return true;
		else
			return false;
	}
};
 