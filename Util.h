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

class Util {
public:
    static double PI;
    static double EARTH_RADIUS;
    // ����·��Json
    static string getGeoJson(vector<H3Index> result, string name);
    static string getGeoJson(vector<pair<vector<H3Index>, point_type>> result, string name);
    static double calcdistance(GeoPoint g1, GeoPoint g2);
    // �����ڽӹ�ϵ
    static void buildAijk(H3_N& PointMap);
    static void buildAijk(H3_P& PointMap);
    // ��������
    static void loadingData(vector<H3_N>& PointMap, vector<H3_D>& DEM, vector<H3_D>& Comprehensive, vector<H3_V>& Factor,
        int level, string filepath);
    // ���ɶ��θ���ͨ��ģ��
    static void MultiHierarchy(vector<H3_N>& PointMap, vector<H3_D>& DEM, vector<H3_D>& Comprehensive, vector<H3_V>& Factor,
        int level, string filepath);
    // ����ر�������ռ��
    static void ComputeFactor(vector<H3Index>& result, H3_V& Factor, vector<double>& factor);
    static void ComputeFactor(vector<pair<vector<H3Index>, point_type>>& result, H3_V& Factor, vector<vector<double>>& factor);
    // ��ƽ��ֵ 
    static double average(vector<double>& x);
    // �󷽲�
    static double variance(vector<double>& x);
    // ���׼��
    static double standardDev(vector<double>& x);
    // ����С��
    static double round(double number, unsigned int bits) {
        /* stringstream ss;
         ss << fixed << setprecision(bits) << number;
         ss >> number;*/

        return (int)(number * pow(10, bits));
    };
    // �����¶�ƽ��ֵ���׼��
    static void ComputeGrad(vector<H3Index>& result, H3_D& DEM, double& mean, double& dev);
    static void ComputeGrad(vector<pair<vector<H3Index>, point_type>>& result, H3_D& DEM, vector<double>& mean, vector<double>& dev);
    // ����·���滮
    static vector<H3Index> MultiHierarchySearch(vector<H3_N>& PointMap, vector<H3_D>& DEM, vector<H3_D>& Comprehensive, int targetLevel, int level, H3Index start, H3Index end, ofstream& log, int& maxSize, int& Size);
    // ��ͨ·���滮
    static vector<H3Index> OrdinarySearch(H3_N& PointMap, H3_D& DEM, H3_D& Comprehensive, H3Index start, H3Index end, ofstream& log, int& maxSize, int& Size);
    // ͳ������    
    static void AnalysisData(H3_P& PointMap, H3_D& DEM, H3_D& Comprehensive, vector<pair<vector<H3Index>, point_type>>& paths, vector<double> & Length_of_Path);
    static void AnalysisData(H3_N& PointMap, H3_D& DEM, H3_D& Comprehensive, vector<H3Index>& path, double& Length_of_Path);
    // �����Ż�
    static void ParallelSearch(H3_N PointMap, H3_D DEM, H3_D Comprehensive, vector<H3Index>& tmp, H3Index start, H3Index end);
    
    // ��Ŀ��·���滮
    static vector<pair<vector<H3Index>, point_type>> MultiObjectMultiHierarchySearch(vector<H3_P>& PointMap, vector<H3_D>& DEM, vector<H3_D>& Comprehensive, int targetLevel, int level,
        H3Index start, H3Index end, ofstream& log, int& maxSize, int& Size, double& time);

    // ��Ŀ����θ���ͨ��ģ������
    static void loadingData(vector<H3_P>& PointMap, vector<H3_D>& DEM, vector<H3_D>& Comprehensive, vector<H3_V>& Factor,
        int level, string filepath);

    static void MultiObjectMultiHierarchy(vector<H3_P>& PointMap, vector<H3_D>& DEM, vector<H3_D>& Comprehensive, vector<H3_V>& Factor,
        int level, string filepath);

    // ����·��
    static vector<pair<vector<H3Index>, point_type>> getPathList(H3_P& NodeMap, H3_D& DEM, H3_D& Comprehensive, H3Index start, H3Index end);
    static point_type calEdgeCost(H3_D& DEM, H3_D& Comprehensive, Parameter parent, H3Index current);

    static int calGridNum(H3Index p1, H3Index p2);

    static inline vector<double> ball2xyz(double lat, double lng);

    static inline vector<double> geo2xyz(GeoPoint& point);

    static double angleOflocation(GeoPoint l1, GeoPoint l2, GeoPoint l3);
};