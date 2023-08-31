#include "Astar.h"
// 注释
Astar::Astar(H3_N& PointMap, H3_D& DEM, H3_D& Comprehensive)
{
	this->PointMap = PointMap;
	this->DEM = DEM;
	this->Comprehensive = Comprehensive;
}


Astar::~Astar()
{

}


vector<H3Index> Astar::search(H3Index startPos, H3Index endPos, int& maxSize, int& size)
{
	vector<H3Index> path;
	// 判断起始点和终止点是否存在
	if (PointMap.count(startPos) == 0 || PointMap.count(endPos) == 0)
		return path;
	
	// 存入S点
	Node& start = PointMap[startPos];
	start.G = calcG(startPos, startPos);
	start.H = calcH(startPos, endPos);
	start.F = start.G + start.H;
	OPEN.push(make_pair(start.F, startPos));
	Openlist.insert(make_pair(startPos, true));
	while (!OPEN.empty())
	{
		//取F值最小的点
		H3Index current = OPEN.top().second;
		maxSize = max(maxSize, (int)OPEN.size());
		OPEN.pop();
		Openlist.erase(current);
		Closelist.insert({current, true});
		if (current == endPos)
		{
			// 输出路径
			while (current != 0) {
				path.push_back(current);
				Node& node = PointMap[current];
				current = node.parent;
			}
			return path;
		}
		else
			NextStep(current, endPos, size);
		//cout << "current Point " << current << " =============QueueSize:" << OPEN.size() << endl;
	}
}

void Astar::NextStep(H3Index current, H3Index endPos, int& size)
{
	auto neighbours = PointMap[current].neighbours;
	// 遍历邻接点
	for (auto it = neighbours.begin(); it != neighbours.end(); it++) {
		H3Index child = *it;
		// 判断通行能力是否达标、坡度是否能通行、是否已经遍历过（因为通行能力可能为0，所以可能会出现循环）
		if (Comprehensive[child] <= 0 || !check(current, child) || Closelist.count(child) != 0)
			continue;
		Node& node = PointMap[child];
		// 如果当前点在OPEN列表，那么视情况更新
		if (Openlist.count(child) != 0) {
			double tempG = calcG(child, current);
			if (tempG < node.G) {
				node.parent = current;
				node.G = tempG;
				node.F = node.G + node.H;
			}
		}
		// 如果是新点
		else {
			node.parent = current;
			node.G = calcG(child, current);
			node.H = calcH(child, endPos, current);
			node.F = node.G + node.H;
			OPEN.push(make_pair(node.F, child));
			size++;
			Openlist.insert(make_pair(child, true));
		}
	}
}

double Astar::calcH(H3Index current, H3Index end) {
	GeoPoint g1, g2;
	cellToPoint(current, &g1);
	cellToPoint(end, &g2);
	double distance = Util::calcdistance(g1, g2);
	double edge = getHexagonEdgeLengthAvgM(getResolution(end));
	int grid = distance / (1.732050 * edge) + 1;
	double h = fabs(DEM[current] - DEM[PointMap[current].parent]) / 1000;
	double s = h / distance;
	double grad = atan(s) * 180.0 / 3.1415926;
	grad = -0.00004061 * pow(grad, 3) + 0.002644 * pow(grad, 2) - 0.07607 * grad + 0.9977;
	if (grad > 1)
		grad = 1;
	if (grad < 0)
		grad = 0;
	//return (4 - grad - Comprehensive[current]) / 2 * distance;
	//return (2 - Comprehensive[current]) * distance;
	return Util::round(distance - Comprehensive[current] * 100, 0);
	//return  grid * 100 - Util::round(grad, 2);
}

double Astar::calcH(H3Index current, H3Index end, H3Index parent) {
	GeoPoint g1, g2, g3;
	cellToPoint(current, &g1);
	cellToPoint(end, &g2);
	cellToPoint(parent, &g3);
	double distance = Util::round(Util::calcdistance(g1, g2), 3);
	double edge = getHexagonEdgeLengthAvgM(getResolution(end));
	int grid = distance / (1.732050 * edge) + 1;
	double h = fabs(DEM[current] - DEM[PointMap[current].parent]);
	double s = h / distance;
	double grad = atan(s) * 180.0 / 3.1415926;
	grad = -0.00004061 * pow(grad, 3) + 0.002644 * pow(grad, 2) - 0.07607 * grad + 0.9977;
	if (grad > 1)
		grad = 1;
	if (grad < 0)
		grad = 0;
	// double weight = pow(cos(Util::angleOflocation(g2, g3, g1)) + 1, 2);
	// double weight = cos(degsToRads(Util::angleOflocation(g2, g3, g1))) + 1;
	// return Util::round(distance - (weight * Comprehensive[current]) * 100, 0);
	// return Util::round(distance - (1 - Comprehensive[current]) * 100, 0);
	// return Util::round(distance - (weight * grad) * 100, 0);
	return Util::round(distance + (1 - grad) * 100, 0);
}

double Astar::calcG(H3Index current, H3Index parent) {
	if (current == parent) 
		return 0;
	else {
	    GeoPoint g1, g2;
		cellToPoint(current, &g1);
		cellToPoint(parent, &g2);
		return PointMap[parent].G + Util::round(Util::calcdistance(g1, g2), 3);
		//return PointMap[parent].G + Util::calGridNum(parent, current) * 100;
	/*	return PointMap[parent].G + (1 - Comprehensive[current]);*/
	}
}

bool Astar::check(H3Index p1, H3Index p2) {
	GeoPoint g1, g2;
	cellToPoint(p1, &g1);
	cellToPoint(p2, &g2);
	double distance = Util::round(Util::calcdistance(g1, g2), 3);
	double h = fabs(DEM[p1] - DEM[p2]);
	double s = h / distance;
	double grad = atan(s) * 180.0 / 3.1415926;
	grad = -0.00004061 * pow(grad, 3) + 0.002644 * pow(grad, 2) - 0.07607 * grad + 0.9977;
	if (grad > 1 || grad < 0)
		return false;
	return true;
}

