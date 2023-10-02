#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <set>
#include <ctime>
#include <unordered_map>
#include <fstream>
#include "h3api.h"
#include "Util.h"
#include "include/pareto/spatial_map.h"
#include "Multi_Astar.h"
#pragma warning(disable:4996)
using namespace std;

// beanch123
// 用来测试git 嘻嘻
void MultihierarchyVsOrdinary() {
    //Util::cnt = 0;
    // 0-15对应着H3的0-15层索引；16为以自定义层次为基准的多层次格网地图索引
    
    vector<H3_V> Factor(17, H3_V());
    vector<H3_D> DEM(17, H3_D());
    vector<H3_D> Comprehensive(17, H3_D());
    vector<H3_N> PointMap(17, H3_N());
    //用于统计参数数据
    ofstream log("D:/桌面/log.txt");
    //基准格网层级
    int level = 10;
    //多层次格网层级
    int target_level = 16;
    //基准格网数据路径
    string filepath = "D:/桌面/路径规划/已有资料/改进算法后实验结果/实验数据_10.csv";
    string m_filepath = "D:/桌面/路径规划/已有资料/改进算法后实验结果/other_compact_h3_res=10_算法用.csv";

    cout << "加载数据..." << endl;
    // 读取基准格网数据
    Util::loadingData(PointMap, DEM, Comprehensive, Factor, level, filepath);
    // 读取多层次格网数据，并生成多层次格网地图
    Util::MultiHierarchy(PointMap, DEM, Comprehensive, Factor, level, m_filepath);
    // 多粒度模型
    H3Index h3_start = 0x8a40088f1107fff;
    H3Index h3_end = 0x8a400881aa67fff;
    vector<H3Index> result;
    cout << "寻路启动..." << endl;
    cout << "多粒度：" << endl;
    log << "多粒度：" << endl;
    // 实验次数
    double Turn = 1;
    //用于统计时间数据
    clock_t start, finish;
    double totaltime = 0;
    int Max_of_the_Number_of_Open_Grid = 0;
    int Traverses_the_Number_of_Grid = 0;
    double Number_of_Path_Grid = 0;
    double Length_of_Path = 0;
    double mean = 0;
    double dev = 0;
    vector<double> factor(7, 0);
    for (int i = 0; i < Turn; i++) {
        start = clock();
        result = Util::MultiHierarchySearch(PointMap, DEM, Comprehensive, target_level, level, h3_start, h3_end, log, Max_of_the_Number_of_Open_Grid, Traverses_the_Number_of_Grid);
        finish = clock();
        totaltime += (double)(finish - start) / CLOCKS_PER_SEC;
        Number_of_Path_Grid += result.size();
        Util::AnalysisData(PointMap[target_level], DEM[target_level], Comprehensive[target_level], result, Length_of_Path);
        Util::ComputeFactor(result, Factor[target_level], factor);
        Util::ComputeGrad(result, DEM[target_level], mean, dev);
    }
    //时间，格网数，路径格网，路径长度
    log << "时间：" << totaltime / Turn << "  Open队列最大深度：" << Max_of_the_Number_of_Open_Grid / Turn << " 格网遍历数：" << Traverses_the_Number_of_Grid / Turn << "; 路径格网数：" << Number_of_Path_Grid / Turn << "; 路径长度：" << Length_of_Path / Turn << ";";
    // 因子
    log << " 因子：";
    for (int i = 0; i < factor.size(); i++) {
        log << factor[i] / Turn << " ";
    }
    // 坡度均值与标准差
    log << "; 均值：" << mean / Turn << "; 标准差：" << dev / Turn << endl;
    ofstream outfile("C:\\Users\\69025\\Desktop\\路径规划\\论文实验\\精细化数据\\多粒度.geojson");
    outfile << Util::getGeoJson(result, "多粒度");
    outfile.close();
    //普通模型
 /*   h3_start = 0x8a40088c20d7fff;
    h3_end = 0x8a40088e1717fff;*/
    totaltime = 0;
    Max_of_the_Number_of_Open_Grid = 0;
    Traverses_the_Number_of_Grid = 0;
    Number_of_Path_Grid = 0;
    Length_of_Path = 0;
    mean = 0;
    dev = 0;
    factor = vector<double>(factor.size(), 0);
    cout << "普通：" << endl;
    log << "普通：" << endl;
    for (int i = 0; i < Turn; i++) {
        start = clock();
        result = Util::OrdinarySearch(PointMap[level], DEM[level], Comprehensive[level], h3_start, h3_end, log, Max_of_the_Number_of_Open_Grid, Traverses_the_Number_of_Grid);
        finish = clock();
        totaltime += (double)(finish - start) / CLOCKS_PER_SEC;
        Number_of_Path_Grid += result.size();
        Util::AnalysisData(PointMap[level], DEM[level], Comprehensive[level], result, Length_of_Path);
        Util::ComputeFactor(result, Factor[level], factor);
        Util::ComputeGrad(result, DEM[level], mean, dev);
    }
    // 时间，格网数，路径格网，路径长度
    log << "时间：" << totaltime / Turn << "  Open队列最大深度：" << Max_of_the_Number_of_Open_Grid / Turn << " 格网遍历数：" << Traverses_the_Number_of_Grid / Turn << "; 路径格网数：" << Number_of_Path_Grid / Turn << "; 路径长度：" << Length_of_Path / Turn << ";";
    // 因子
    log << " 因子：";
    for (int i = 0; i < factor.size(); i++) {
        log << factor[i] / Turn << " ";
    }
    // 坡度均值与标准差
    log << "; 均值：" << mean / Turn << "; 标准差：" << dev / Turn << endl;
    cout << "寻路结束..." << endl;
    ofstream outfile1("C:\\Users\\69025\\Desktop\\路径规划\\论文实验\\精细化数据\\普通.geojson");
    outfile1 << Util::getGeoJson(result, "普通");
    outfile1.close();
    // 执行Python脚本
    string showInKepler = "python D:\\PythonProject\\Test\\new\\start.py -n 多粒度 -m 普通";
    system(showInKepler.c_str());
    log.close();
    Py_Finalize();
}

// 多目标
void MultiObjectVsSingleObject() {
    // 0-15对应着H3的0-15层索引；
    // by XGG：16为以自定义层次为基准的多层次格网地图索引  自定义层次是第十层  加入了影响因子  第17层（16）
    vector<H3_V> Factor(17, H3_V());//17个H3_V对象
    vector<H3_D> DEM(17, H3_D());
    vector<H3_D> Comprehensive(17, H3_D());
    vector<H3_P> PointMap(17, H3_P());
    //用于统计参数数据
    ofstream log("D:/桌面/log1001.txt",ios::app);
    //基准格网层级
    int level = 10;
    //多层次格网层级
    int target_level = 16;
    //基准格网数据路径
    // string filepath = "C:\\Users\\69025\\Desktop\\路径规划\\论文实验\\英文论文实验\\改进算法后实验结果\\实验数据_10.csv";
    //string filepath = "D:/桌面/路径规划/已有资料/改进算法后实验结果/实验数据_10.csv";
    string filepath = "D:/桌面/966数据_10.csv";
    // string m_filepath = "C:\\Users\\69025\\Desktop\\路径规划\\论文实验\\英文论文实验\\改进算法后实验结果\\other_compact_h3_res=10_算法用.csv";
    //string m_filepath = "D:\\桌面\\改进算法后实验结果\\other_compact_h3_res=10_算法用.csv";
    string m_filepath = "D:/桌面/路径规划/已有资料/改进算法后实验结果/other_compact_h3_res=10_算法用.csv";
    cout << "加载数据..." << endl;
    // 读取基准格网数据
    Util::loadingData(PointMap, DEM, Comprehensive, Factor, level, filepath);   
    // 读取多层次格网数据，并生成多层次格网地图
    //Util::MultiObjectMultiHierarchy(PointMap, DEM, Comprehensive, Factor, level, m_filepath);
    // 多粒度模型
    // 
     //部分数据的开始点和结束点
    H3Index h3_start = 0x8a4008b9d06ffff;
    H3Index h3_end = 0x8a4008b9eadffff;


    //// 完整数据的起点和终点
    /*H3Index h3_start = 0x8a400889c157fff;
    H3Index h3_end = 0x8a400881b487fff;*/
    
    vector<pair<vector<H3Index>, point_type>> result;
    cout << "寻路启动..." << endl;
    cout << "多粒度：" << endl;
    //auto time = Util::PrintCurrentTime();
    log << "多粒度：" << endl;
    //log << time << endl;
    //
    // 实验次数
    int Turn = 1;
    //用于统计时间数据
    double totaltime = 0;
    int Max_of_the_Number_of_Open_Grid = 0;
    int Traverses_the_Number_of_Grid = 0;
    double Number_of_Path_Grid = 0;
    vector<double> Length_of_Path;//路径长度
    vector<double> mean;//坡度均值
    vector<double> dev;//坡度标准差
    vector<vector<double>> factor;
    //for (int i = 0; i < Turn; i++) {
    //    double time;
    //    result = Util::MultiObjectMultiHierarchySearch(PointMap, DEM, Comprehensive, target_level, level, h3_start, h3_end, log, Max_of_the_Number_of_Open_Grid, Traverses_the_Number_of_Grid, time);
    //    if (i == 0) {
    //        Util::AnalysisData(PointMap[target_level], DEM[target_level], Comprehensive[target_level], result, Length_of_Path);
    //       Util::ComputeFactor(result, Factor[target_level], factor);
    //        Util::ComputeGrad(result, DEM[target_level], mean, dev);
    //    }
    //    totaltime += time;
    //}
    //if (result.size() == 0)
    //    cout << "未找到路径" << endl;
    //log << "时间：" << totaltime / Turn << "  Open队列最大深度：" << Max_of_the_Number_of_Open_Grid / Turn << " 格网遍历数：" << Traverses_the_Number_of_Grid / Turn << ";" << endl;
    //cout << "时间：" << totaltime / Turn << "  Open队列最大深度：" << Max_of_the_Number_of_Open_Grid / Turn << " 格网遍历数：" << Traverses_the_Number_of_Grid / Turn << ";" << endl;
    //// 分路径
    //for (int i = 0; i < result.size(); i++) {
    //    log << "Path {" << i << "} :";
    //    cout << "Path {" << i << "} :";
    //    log << "路径格网数：" << result[i].first.size() << "; 路径长度：" << Length_of_Path[i] << "; 因子：";
    //    cout << "路径格网数：" << result[i].first.size() << "; 路径长度：" << Length_of_Path[i] << "; 因子：";
    //    // 因子
    //    for (int j = 0; j < factor[i].size(); j++)
    //    {
    //        log << factor[i][j] * 100 << " ";
    //        cout << factor[i][j] * 100 << " ";
    //    }
    //        
    //    // 坡度均值与标准差
    //    log << "; 均值：" << mean[i] << "; 标准差：" << dev[i] << endl;
    //    cout << "; 均值：" << mean[i] << "; 标准差：" << dev[i] << endl;
    //}
    ////ofstream outfile("D:\\桌面\\改进算法后实验结果\\多粒度.geojson");
    //ofstream outfile("D:/桌面/多粒度.geojson");
    //outfile << Util::getGeoJson(result, "多粒度");
    //outfile.close();
    cout << "普通：" << endl;
    auto time = Util::PrintCurrentTime();
    log << time << endl;
    log << "普通：" << endl;
  //  //totaltime = 0;
  ///*  Max_of_the_Number_of_Open_Grid = 0;
  //  Traverses_the_Number_of_Grid = 0;
  //  Length_of_Path = vector<double>();
  //  mean = vector<double>();
  //  dev = vector<double>();
  //  factor = vector<vector<double>>();*/
    for (int i = 0; i < Turn; i++) {
        double time;
        //result = Util::MultiObjectMultiHierarchySearch(PointMap, DEM, Comprehensive, level, level, h3_start, h3_end, log, Max_of_the_Number_of_Open_Grid, Traverses_the_Number_of_Grid, time);
        result = Util::MultiObjectMultiHierarchySearch(PointMap, DEM, Comprehensive, level, level, h3_start, h3_end, log, Max_of_the_Number_of_Open_Grid, Traverses_the_Number_of_Grid, time);
        if (i == 0) {
            //Util::AnalysisData(PointMap[target_level], DEM[target_level], Comprehensive[target_level], result, Length_of_Path);
            Util::AnalysisData(PointMap[level], DEM[level], Comprehensive[level], result, Length_of_Path);
            cout << "analysisData运行完毕" << endl;
            Util::ComputeFactor(result, Factor[level], factor);
            Util::ComputeGrad(result, DEM[level], mean, dev);
        }
        totaltime += time;
    }
    if (result.size() == 0) {
        cout << "未找到路径" << endl;
        return;
    }
    log << "时间：" << totaltime / Turn << "  Open队列最大深度：" << Max_of_the_Number_of_Open_Grid / Turn << " 格网遍历数：" << Traverses_the_Number_of_Grid / Turn << ";" << endl;
    cout << "时间：" << totaltime / Turn << "  Open队列最大深度：" << Max_of_the_Number_of_Open_Grid / Turn << " 格网遍历数：" << Traverses_the_Number_of_Grid / Turn << ";" << endl;
    // 分路径
    for (int i = 0; i < result.size(); i++) {
        log << "Path {" << i << "} :";
        cout << "Path {" << i << "} :";
        log << "路径格网数：" << result[i].first.size() << "; 路径长度：" << Length_of_Path[i] << "; 因子：";
        cout << "路径格网数：" << result[i].first.size() << "; 路径长度：" << Length_of_Path[i] << "; 因子：";
        // 因子
        for (int j = 0; j < factor[i].size(); j++)
        {
            log << factor[i][j] * 100 << " ";
            cout << factor[i][j] * 100 << " ";
        }
            
        // 坡度均值与标准差
        log << "; 均值：" << mean[i] << "; 标准差：" << dev[i] << endl;
        cout << "; 均值：" << mean[i] << "; 标准差：" << dev[i] << endl;
        
    }
    ofstream outfile1("D:/桌面/普通.geojson");
    //ofstream outfile1("D:\\桌面\\改进算法后实验结果\\普通.geojson");
    outfile1 << Util::getGeoJson(result, "普通");
    outfile1.close();
    log.close();
    //log_size.close();

    
    
//    string showInKepler = "python D:/工作空间/pyworkspace/Test/newAlgorithm/start.py -n 多粒度 -m 普通";
//    //string showInKepler = "python D:\\PythonProject\\Test\\newAlgorithm\\start.py -n 多粒度 -m 普通";
//    system(showInKepler.c_str());
}

void Voronoi()
{
    // 1. 构建distance obst
    // 2. 找GVD边

    // 对应类的near_distance
    // 对应 is_edge
    // 对应is_site
    H3_D Distance;
    H3_B Edge;
    H3_B site;
    vector<H3Index> data;
    H3_D cellList;

    string filepath1 = "D:/桌面/H3index_site.csv";
    Util::loading_site_data(filepath1, cellList);
    cout << cellList.size() << endl;// 966条数据

    // 1.给站点编号
    // 2.计算网格到各个站点的距离 按照距离从小到大排序 得到distance 和 obst 矩阵
    




}
bool com(const double& a, const double& b)
{
    return a < b;
}
void test() {

    ofstream log_test;
    pareto::front<double, 2, string> m;
    pareto::spatial_map<double, 2, string> open;
    vector<double> x;
    m(1, 5) = "A";
    m(4, 1) = "B";
    m(4, 2) = "C";
    m(5, 2) = "D";
    m(4.5, 1.5) = "E";
    m(2, 3) = "F";
    open(1, 5) = "A";
    open(4, 1) = "B";
    open(4, 2) = "C";
    open(5, 2) = "D";
    open(4.5, 1.5) = "E";
    open(2, 3) = "F";
    /*for (auto it = m.begin(); it != m.end(); it++)
    {
        cout << m.crowding_distance(it) << endl;
    }*/
    //for (auto& [key, value] : m)
    //{
    //    //cout << m.crowding_distance(key) << endl;
    //    double cd = m.crowding_distance(key);
    //    x.push_back(cd);
    //}
    //sort(x.begin(), x.end(), com);

    //for (auto it = x.begin(); it != x.end(); it++)
    //{
    //    cout << *it << endl;
    //}
    std::vector<bool> is_minimization = { true, true };
    for (auto& [key, value] : m)
    {
        double cnt = 0;
        cout << value << "支配几个" << endl;
        for (auto& [key1, value1] : open)
        {
            if (key.dominates(key1, is_minimization))
            {
                cnt++;
            }
        }
        cout << cnt << endl;
    }

    /*cout << "Uniformity: " << m.uniformity() << std::endl;
    std::cout << "ACrowding distance: " << m.crowding_distance(m.begin()->first)<< std::endl;
    std::cout << "ACrowding distance: " << m.crowding_distance({1,5}) << std::endl;
    std::cout << "BCrowding distance: " << m.crowding_distance({4,1}) << std::endl;
    std::cout << "CCrowding distance: " << m.crowding_distance({4,2}) << std::endl;
    std::cout << "DCrowding distance: " << m.crowding_distance({5,2}) << std::endl;
    std::cout << "ECrowding distance: " << m.crowding_distance({4.5,1.5}) << std::endl;
    std::cout << "FCrowding distance: " << m.crowding_distance({2,3}) << std::endl;
    std::cout << "avg: " << m.average_crowding_distance() << std::endl;*/
    
    //auto it = m.find_nearest({4, 1});
    //std::cout << it->first << " -> " << it->second << std::endl;
}

void draw()
{
    double a[] = { 1.2, 2.0, 3, 4, 5, 6, 3, 2, 3, 2, 31, 2, 1, 4 };
    double b[] = { 2, 7, 8, 9, 5, 3, 2, 7, 3, 7, 1, 2, 5, 4 };
    /*vector<double> a = { 1,2,3,4,5 };
    vector<double> b = { 2,4,6,2,1 };*/
    Util::pythonInitial();
    //plot(a, 14); // 14是数组长度
    //Util::scatter(&a[0], 5, &b[0]);
    Util::plot(a,14, b);
    Py_Finalize(); /*结束python解释器，释放资源*/
    system("pause");
}
int main()
{
    //test();
    //Voronoi();
    //MultihierarchyVsOrdinary();
    //draw();
    Py_Initialize();
    string chdir_cmd = string("sys.path.append(\"") + "." + "\")";

    const char* cstr_cmd = chdir_cmd.c_str();
    PyRun_SimpleString("import sys");
    PyRun_SimpleString(cstr_cmd);
    // 语句1 
    PyRun_SimpleString("import matplotlib.pyplot as plt");

    MultiObjectVsSingleObject();
    // by xgg
    //array_type g_m({ 1, 2 });
    ////cout << g_m[0] << endl;
    ////cout << g_m[1] << endl;
    //point_type point_g_m(g_m.begin(), g_m.end());
    //cout << typeid(point_g_m).name() << endl;
    //cout << point_g_m[0] << endl;
    //cout << point_g_m[1] << endl;
    //cout << point_g_m[2] << endl;

    // Point-point dominance
  /*  using point_type = pareto::front<double, 2, unsigned>::key_type;
    point_type p1({ 163, 259 });
    point_type p2({ 163, 269 });
    std::vector<bool> is_minimization = { true, true};
    std::cout << (p1.dominates(p2, is_minimization) ? "p1 dominates p2" : "p1 does not dominate p2") << std::endl;
    std::cout << (p1.strongly_dominates(p2, is_minimization) ? "p1 strongly dominates p2" : "p1 does not strongly dominate p2") << std::endl;
    std::cout << (p1.non_dominates(p2, is_minimization) ? "p1 non-dominates p2" : "p1 does not non-dominate p2") << std::endl;*/

    //Front-point dominance
   /* pareto::front<double, 2, unsigned> pf;
    pf.insert(make_pair(point_type{ 153,279 }, 1));
    pf.insert(make_pair(point_type{ 163,259 }, 1));
    point_type p3({ 163,259 });

    for (auto it = pf.find_dominated(p3); it != pf.end(); ++it) {
        cout << it->first << endl;
    }

    std::cout << (pf.dominates(p3) ? "pf dominates p2" : "pf does not dominate p2") << std::endl;
    std::cout << (pf.strongly_dominates(p3) ? "pf strongly dominates p2" : "pf does not strongly dominate p2") << std::endl;
    std::cout << (pf.non_dominates(p3) ? "pf non-dominates p2" : "pf does not non-dominate p2") << std::endl;
    std::cout << (pf.is_partially_dominated_by(p3) ? "pf is partially dominated by p2" : "pf is not is partially dominated by p2") << std::endl;
    std::cout << (pf.is_completely_dominated_by(p3) ? "pf is completely dominated by p2" : "pf is not is completely dominated by p2") << std::endl;*/

    //GeoPoint p1;
    //GeoPoint p2;
    //GeoPoint p3;

    //p1.lat = degsToRads(40.823978);
    //p2.lat = degsToRads(40.71986);
    //p3.lat = degsToRads(40.484791);
    //p1.lon = degsToRads(124.639313);
    //p2.lon = degsToRads(124.78175);
    //p3.lon = degsToRads(124.81823);

    //cout << Util::angleOflocation(p1, p2, p3);
}