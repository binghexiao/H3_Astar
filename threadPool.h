#pragma once
#include <vector>
#include <queue>
#include <thread>
#include <mutex>
#include "Task.h"
using namespace std;

class threadPool {
public:
	threadPool(int maxThreadsNum = 1) :
		mMaxThreadsNum(maxThreadsNum) {
		mReceStopOrder = false;
		//creat threads:
		for (int i = 0; i < mMaxThreadsNum; ++i) {
			shared_ptr<thread> th = std::make_shared<thread>(&threadPool::threadFunc, this);//Ok
			//shared_ptr<thread> th = std::make_shared<thread>(std::bind(&threadPool::threadFunc, this)); //Ok
			//shared_ptr<thread> th = std::make_shared<thread>(std::bind(&threadPool::threadFunc, *this));//build error

			//newThread.detach();
			threadCache.push_back(th);
		}
	}

	void joinThreads() {
		for (auto& it : threadCache) {
			it->join();
		}
	}
	virtual ~threadPool() {
		joinThreads();//等待所有线程结束
	}

	void pushTask(Task task) {

		std::lock_guard<std::mutex> lk(mMux);
		mTasksQueue.push(task);
		mCondVar.notify_one();
	}

	//构造函数中开辟的多个线程的入口
	void threadFunc() {
		while (true)
		{
			//终止函数,停止线程
			if (mReceStopOrder) {
				//myLogger.logData("break---");
				break;
			}
			std::unique_lock<std::mutex> lk(mMux);
			while (mTasksQueue.empty() && !mReceStopOrder)
			{
				mCondVar.wait(lk);
			}
			if (mReceStopOrder) {
				//myLogger.logData("break");
				break;
			}
			//get a task and execute it:
			Task task = mTasksQueue.front();
			mTasksQueue.pop();
			lk.unlock();//显式地解锁.
			task.run();
		}

		//thread::id threadId = std::this_thread::get_id();
		//string data = "ThreadId " + getThreadIdOfString(threadId) + " exit \n";
		//myLogger.logData(data);
	}
	void stopAllThreads() {
		//myLogger.logData("stopAllThreads was called!");
		mReceStopOrder = true;
		mCondVar.notify_all();//可能有的线程处于wait状态
	}

private:
	unsigned int mMaxThreadsNum;
	std::queue<Task> mTasksQueue; //缓存各个任务
	vector<shared_ptr<thread>> threadCache; //保存开辟的各个线程
	std::mutex mMux;
	std::condition_variable mCondVar;//当任务队列mTasksQueue不为空时，唤醒一个线程，从队列头取出一个task执行。
	std::atomic<bool> mReceStopOrder;//用于控制线程的停止
};