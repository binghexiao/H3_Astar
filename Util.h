#pragma once
#include <iostream>
#include <map>
#include <unordered_map>
#include <fstream>
#include <sstream>
#include <string>
#include <set>
#include <ctime>
#include "h3api.h"
#include "Node.h"
#include "type_def.h"
#include <iomanip>
#include "C:/Python311/include/Python.h"
#include <chrono>

// 测试 嘻嘻
class Util {
public:
    static double PI;
    static double EARTH_RADIUS;
    // 生成路径Json
    static string getGeoJson(vector<H3Index> result, string name);
    static string getGeoJson(vector<pair<vector<H3Index>, point_type>> result, string name);
    static double calcdistance(GeoPoint g1, GeoPoint g2);
    // 构建邻接关系
    static void buildAijk(H3_N& PointMap);
    static void buildAijk(H3_P& PointMap);
    // 读入数据
    static void loadingData(vector<H3_N>& PointMap, vector<H3_D>& DEM, vector<H3_D>& Comprehensive, vector<H3_V>& Factor,
        int level, string filepath);
    // 生成多层次格网通行模型
    static void MultiHierarchy(vector<H3_N>& PointMap, vector<H3_D>& DEM, vector<H3_D>& Comprehensive, vector<H3_V>& Factor,
        int level, string filepath);
    // 计算地表覆盖因子占比
    static void ComputeFactor(vector<H3Index>& result, H3_V& Factor, vector<double>& factor);
    static void ComputeFactor(vector<pair<vector<H3Index>, point_type>>& result, H3_V& Factor, vector<vector<double>>& factor);
    // 求平均值 
    static double average(vector<double>& x);
    // 求方差
    static double variance(vector<double>& x);
    // 求标准差
    static double standardDev(vector<double>& x);
    // 保留小数
    static double round(double number, unsigned int bits) {
        /* stringstream ss;
         ss << fixed << setprecision(bits) << number;
         ss >> number;*/

        return (int)(number * pow(10, bits));
    };
    // 计算坡度平均值与标准差
    static void ComputeGrad(vector<H3Index>& result, H3_D& DEM, double& mean, double& dev);
    static void ComputeGrad(vector<pair<vector<H3Index>, point_type>>& result, H3_D& DEM, vector<double>& mean, vector<double>& dev);
    // 多层次路径规划
    static vector<H3Index> MultiHierarchySearch(vector<H3_N>& PointMap, vector<H3_D>& DEM, vector<H3_D>& Comprehensive, int targetLevel, int level, H3Index start, H3Index end, ofstream& log, int& maxSize, int& Size);
    // 普通路径规划
    static vector<H3Index> OrdinarySearch(H3_N& PointMap, H3_D& DEM, H3_D& Comprehensive, H3Index start, H3Index end, ofstream& log, int& maxSize, int& Size);
    // 统计数据    
    static void AnalysisData(H3_P& PointMap, H3_D& DEM, H3_D& Comprehensive, vector<pair<vector<H3Index>, point_type>>& paths, vector<double> & Length_of_Path);
    static void AnalysisData(H3_N& PointMap, H3_D& DEM, H3_D& Comprehensive, vector<H3Index>& path, double& Length_of_Path);
    // 并行优化
    static void ParallelSearch(H3_N PointMap, H3_D DEM, H3_D Comprehensive, vector<H3Index>& tmp, H3Index start, H3Index end);
    
    // 多目标路径规划
    static vector<pair<vector<H3Index>, point_type>> MultiObjectMultiHierarchySearch(vector<H3_P>& PointMap, vector<H3_D>& DEM, vector<H3_D>& Comprehensive, int targetLevel, int level,
        H3Index start, H3Index end, ofstream& log, int& maxSize, int& Size, double& time);

    // 多目标多层次格网通行模型生成
    static void loadingData(vector<H3_P>& PointMap, vector<H3_D>& DEM, vector<H3_D>& Comprehensive, vector<H3_V>& Factor,
        int level, string filepath);

    static void MultiObjectMultiHierarchy(vector<H3_P>& PointMap, vector<H3_D>& DEM, vector<H3_D>& Comprehensive, vector<H3_V>& Factor,
        int level, string filepath);

    // 计算路径
    static vector<pair<vector<H3Index>, point_type>> getPathList(H3_P& NodeMap, H3_D& DEM, H3_D& Comprehensive, H3Index start, H3Index end);
    static point_type calEdgeCost(H3_D& DEM, H3_D& Comprehensive, Parameter parent, H3Index current);

    static int calGridNum(H3Index p1, H3Index p2);

    static inline vector<double> ball2xyz(double lat, double lng);

    static inline vector<double> geo2xyz(GeoPoint& point);

    static double angleOflocation(GeoPoint l1, GeoPoint l2, GeoPoint l3);

    static void loading_site_data(string filePath, H3_D& cellList);

    static void SetNumber(H3_D &sitelist);
    static void callPython(vector<double> parList);
    
    template<class T>
    // T* arr的问题：
    // 1. T* 代表是一个指针 arr代表首地址
    //  (1) 传递普通数组 数组名就是首地址 没有问题
    //  (2) 传递vector的时候 写成 &a[0] 
    static string arr_to_string_list(T* arr, int N) {
        string s = "[";
        for (int i = 0; i < N; ++i) {
            s += to_string(arr[i]);
            if (i != N - 1) s += ",";
        }
        s += "]";
        return s;
    }

    template<class T, class V = int>
    static void plot(T* x, int N1, V* y = NULL, bool equal = false, string marker = "^") {
        PyRun_SimpleString("import matplotlib.pyplot as plt");
        /*if (equal) {
            PyRun_SimpleString("plt.axis(\"equal\")");
        }*/

        string cmd = "plt.plot(";
        string s1 = arr_to_string_list(x, N1);
        if (y != NULL) {
            string s2 = arr_to_string_list(y, N1);
            cmd += (s1 + "," + s2 + "," + "marker = \"^\"" + ")");
            PyRun_SimpleString(cmd.c_str());
        }
        else {
            cmd += (s1 + ")");
            PyRun_SimpleString(cmd.c_str());
        }
        PyRun_SimpleString("plt.show()");
    }

    template<class T, class V>
    static void scatter(T* x, int N1, V* y = NULL, bool equal = false)
    {
        PyRun_SimpleString("import matplotlib.pyplot as plt");
        if (equal) {
            PyRun_SimpleString("plt.axis(\"equal\")");
        }

        string cmd = "plt.scatter(";
        string s1 = arr_to_string_list(x, N1);
        if (y != NULL) {
            string s2 = arr_to_string_list(y, N1);
            cmd += (s1 + "," + s2 + ")");
            PyRun_SimpleString(cmd.c_str());
        }
        PyRun_SimpleString("plt.show()");
    }

    template<class T, class V>
    static void scatter_marker(T* x, int N1, V* y = NULL, bool equal = false)
    {
        PyRun_SimpleString("import matplotlib.pyplot as plt");
        if (equal) {
            PyRun_SimpleString("plt.axis(\"equal\")");
        }

        string cmd = "plt.scatter(";
        string s1 = arr_to_string_list(x, N1);
        if (y != NULL) {
            string s2 = arr_to_string_list(y, N1);
            cmd += (s1 + "," + s2 + "marker = \"^\"" + ")");
            PyRun_SimpleString(cmd.c_str());
        }
        PyRun_SimpleString("plt.show()");
    }

    static void pythonInitial() {
        Py_Initialize(); /*初始化python解释器,告诉编译器要用的python编译器*/
        string path = ".";
        // 要的结果是sys.path.append(".")
        // .前后的"需要加转义字符\

        string chdir_cmd = string("sys.path.append(\"") + path + "\")";

        const char* cstr_cmd = chdir_cmd.c_str();
        PyRun_SimpleString("import sys");
        PyRun_SimpleString(cstr_cmd);
    }
    static int cnt;

    static string PrintCurrentTime();
};

