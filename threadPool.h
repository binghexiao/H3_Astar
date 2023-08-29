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
		joinThreads();//�ȴ������߳̽���
	}

	void pushTask(Task task) {

		std::lock_guard<std::mutex> lk(mMux);
		mTasksQueue.push(task);
		mCondVar.notify_one();
	}

	//���캯���п��ٵĶ���̵߳����
	void threadFunc() {
		while (true)
		{
			//��ֹ����,ֹͣ�߳�
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
			lk.unlock();//��ʽ�ؽ���.
			task.run();
		}

		//thread::id threadId = std::this_thread::get_id();
		//string data = "ThreadId " + getThreadIdOfString(threadId) + " exit \n";
		//myLogger.logData(data);
	}
	void stopAllThreads() {
		//myLogger.logData("stopAllThreads was called!");
		mReceStopOrder = true;
		mCondVar.notify_all();//�����е��̴߳���wait״̬
	}

private:
	unsigned int mMaxThreadsNum;
	std::queue<Task> mTasksQueue; //�����������
	vector<shared_ptr<thread>> threadCache; //���濪�ٵĸ����߳�
	std::mutex mMux;
	std::condition_variable mCondVar;//���������mTasksQueue��Ϊ��ʱ������һ���̣߳��Ӷ���ͷȡ��һ��taskִ�С�
	std::atomic<bool> mReceStopOrder;//���ڿ����̵߳�ֹͣ
};