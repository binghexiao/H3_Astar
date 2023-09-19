#pragma once
#include "h3api.h"
#include <unordered_map>

class Vinoroi
{
public:
	// ��ǰ������index
	H3Index index = 0;
	// ��ǰ�����ǲ���վ��
	bool is_site = false;

	// ���뵱ǰ�������վ�������H3Index
	H3Index index_near_site = 0;

	// ���뵱ǰ���������վ�������distance
	double near_distance;
	Vinoroi() {};
	Vinoroi(H3Index index, bool is_site = false, H3Index index_near_site=1, double near_distance=1) :
		index(index), is_site(is_site), index_near_site(index_near_site), near_distance(near_distance) {};

	
};

