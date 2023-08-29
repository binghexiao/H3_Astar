#pragma once
#include <iostream>
#include <fstream>
#include <queue>
#include <vector>
#include <stack>
#include <list>
#include <algorithm>
#include <map>
#include <cmath>
#include "Point.h"
#include <pareto/front.h>
#include "type_def.h"
using namespace std;
static std::initializer_list<bool> direction_comparison;


class Multi_Astar {
public:
	Multi_Astar(H3_P&, std::initializer_list<bool>, H3_D&, H3_D&, int);
	~Multi_Astar();
	front_type search(H3Index, H3Index, int& maxSize, int& Size);
private:
	// 总路径代价集合，是非支配的
	front_type COSTS;
	// g,F,含当前节点的Parameter
	multimap<array_type, pair<array_type, Parameter>> OPEN; //待遍历的节点	
	// 格网地图
	H3_P& NodeMap;

	// 高程和通行能力
	H3_D& DEM;
	H3_D& Comprehensive;
	// baseLevel
	int baseLevel;



	Point& OpenPop(array_type&, Parameter& parameter);
	void NextStep(Point& current, array_type current_g, H3Index goal, int& Size, Parameter& parameter);
	array_type calG(array_type& parent_g, H3Index parent, H3Index current, bool& flag, Parameter& parameter);
	array_type calH(H3Index goal, H3Index current, H3Index parent);
	array_type calH(H3Index goal, H3Index current);
	vector<pair<list<H3Index>, array_type>> printPathByNotRecursion(H3Index start, H3Index end);
	void printPath(Point&, Point&, Point&, list<H3Index>&, array_type&, vector<pair<list<H3Index>, array_type>>&);
};