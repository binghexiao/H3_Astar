#pragma once
#ifndef ASTAR_H
#define ASTAR_H
#include <fstream>
#include <queue>
#include "Util.h"
#include "type_def.h"

using namespace std;
struct cmp //重写仿函数
{
	bool operator() (P_D_H3 a, P_D_H3 b)
	{
		return a.first > b.first; //小顶堆
	}
};

class Astar {
public:
	Astar(H3_N&, H3_D&, H3_D&);
	~Astar();
	vector<H3Index> search(H3Index, H3Index, int&, int&);

private:
	H3_N PointMap;
	H3_D DEM;
	H3_D Comprehensive;
	H3_B Openlist;
	H3_B Closelist;
	priority_queue<P_D_H3, vector<P_D_H3>, cmp> OPEN; //F最小优先队列
	void NextStep(H3Index current, H3Index endPos, int& size);
	double calcH(H3Index current, H3Index end);
	double calcH(H3Index current, H3Index end, H3Index parent);
	double calcG(H3Index current, H3Index start);
	bool check(H3Index p1, H3Index p2);
};
#endif