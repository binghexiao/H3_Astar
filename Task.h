#pragma once
#include <thread>
using namespace std;

class Task {
public:
	Task(int val = 0) :
		mParam(val)
	{

	}
	virtual ~Task() {

	}
	void setParam(int val) {
		mParam = val;
	}
	void run() {
		thread::id threadId = std::this_thread::get_id();
		//cout << "param = " << mParam << " ThreadId = " << threadId << endl;
		//string data = "param = " + std::to_string(mParam) + " ThreadId = " + getThreadIdOfString(threadId) + "\n";
		//myLogger.logData(data);
	}
private:
	int mParam;
};