#pragma once
#ifndef TYPE_DEF_H
#define TYPE_DEF_H
#include <vector>
#include <unordered_map>
#include "h3api.h"
#include "Node.h"
#include "Point.h"
typedef std::unordered_map<H3Index, std::vector<double>> H3_V;
typedef std::unordered_map<H3Index, Node> H3_N;
typedef std::unordered_map<H3Index, double> H3_D;
typedef std::unordered_map<H3Index, bool> H3_B;
typedef std::pair<H3Index, double> P_H3_D;
typedef std::pair<H3Index, Node> P_H3_N;
typedef std::pair<H3Index, std::vector<double>> P_H3_V;
typedef std::pair<double, H3Index> P_D_H3;
/////
typedef std::unordered_map<H3Index, Point> H3_P;

#endif

