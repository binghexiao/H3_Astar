#pragma once
#include "h3api.h"
#include <unordered_map>

class Vinoroi
{
public:
	// 当前格网的index
	H3Index index = 0;
	// 当前格网是不是站点
	bool is_site = false;

	// 距离当前格网最近站点格网的H3Index
	H3Index index_near_site = 0;

	// 距离当前格网最近的站点格网的distance
	double near_distance;
	Vinoroi() {};
	Vinoroi(H3Index index, bool is_site = false, H3Index index_near_site=1, double near_distance=1) :
		index(index), is_site(is_site), index_near_site(index_near_site), near_distance(near_distance) {};

	
};

