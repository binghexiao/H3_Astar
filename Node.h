#pragma once
#include "h3api.h"
#include <unordered_set>
#include <list>
using namespace std;

class Node {
public:
	Node() {

	};

	Node(H3Index index) {
		this->index = index;
	};

	void add(H3Index index) {
		neighbours.insert(index);
	}

	void remove(H3Index index) {
		neighbours.erase(index);
	}

	H3Index index = 0;
	H3Index parent = 0;
	// F = G+H
	double F = 0;
	double H = 0;
	double G = 0;
	// ͨ������
	double speed = 0;
	// �ڽӽڵ㣬������һ��ڵ�
	unordered_set<H3Index> neighbours;
};