#pragma once
#include "h3api.h"
#include <list>
#include <inttypes.h>
#include <pareto/front.h>
#include <set>
#include "parameter.h"
const static int VectorDimension = 2;
using front_type = pareto::front<double, VectorDimension, Parameter>;
using point_type = pareto::point<double, VectorDimension>;
using array_type = std::array<double, VectorDimension>;

class Point
{
public:
	H3Index index;
	// g_m, 含父节点的Parameter
	front_type op;
	// g_m, 含父节点的Parameter
	front_type cl;
	// 邻接节点                                                                        
	// index,cost
	std::set<H3Index> neighbours;

	Point() {
	}

	Point(H3Index index) {
		//op = front_type({ true, true, true});
		//cl = front_type({ true, true, true});
		op = front_type({ true, true});
		cl = front_type({ true, true});
		this->index = index;
	}

	void add(H3Index index) {
		neighbours.insert(index);
	}

	void remove(H3Index index) {
		neighbours.erase(index);
	}

	H3Index getIndex() {
		return index;
	}

	// 比较函数
	bool operator==(Point& t)
	{
		if (index == t.getIndex())
			return true;
		else
			return false;
	}
};