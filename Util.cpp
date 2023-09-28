#pragma once
#include "Util.h"
#include "Astar.h"
#include "Multi_Astar.h"
#include <thread>
#include "h3api.h"
double Util::PI = 3.1415926;
double Util::EARTH_RADIUS = 6371;


string Util::getGeoJson(vector<H3Index> result, string name) {
    GeoPoint geoHQ1;
    stringstream outfile;
    outfile.precision(6);
    outfile.setf(std::ios::fixed);
    outfile << "{\"type\":\"FeatureCollection\",\"features\":[{\"type\":\"Feature\",\"properties\":{\"" + name + "\":\"line\"},\"geometry\":{\"type\":\"LineString\",\"coordinates\":[";
    for (int i = 0; i < result.size(); i++) {
        H3Index index = result[i];
        cellToPoint(index, &geoHQ1);
        outfile << "[" << radsToDegs(geoHQ1.lon) << "," << radsToDegs(geoHQ1.lat) << "]";
        if (i != result.size() - 1)
            outfile << ",";
    }
    outfile << "]}}, {\"type\": \"Feature\", \"properties\" : { \"NAME\": \"point\"}, \"geometry\" : { \"type\":\"MultiPoint\",  \"coordinates\" : [";
    cellToPoint(result[0], &geoHQ1);
    outfile << "[" << radsToDegs(geoHQ1.lon) << "," << radsToDegs(geoHQ1.lat) << "],";
    cellToPoint(result[result.size() - 1], &geoHQ1);
    outfile << "[" << radsToDegs(geoHQ1.lon) << "," << radsToDegs(geoHQ1.lat) << "]";
    outfile << "]}}]}";
    return outfile.str();
}

string Util::getGeoJson(vector<pair<vector<H3Index>, point_type>> result, string name) {
    stringstream outfile;
    GeoPoint geoHQ1;
    outfile.precision(6);
    outfile.setf(std::ios::fixed);
    outfile << "{\"type\":\"FeatureCollection\",\"features\":[";
    for (int i = 0; i < result.size(); i++) {
        auto vec = result[i].first;
        auto cost = result[i].second;
        if (i != 0)
            outfile << ",";
        outfile << "{\"type\":\"Feature\",\"properties\":{\"" + name + "\":\"line\"},\"geometry\":{\"type\":\"LineString\",\"coordinates\":[";
        for (int j = 0; j < vec.size(); j++) {
            H3Index index = vec[j];
            cellToPoint(index, &geoHQ1);
            double lon = radsToDegs(geoHQ1.lon);
            double lat = radsToDegs(geoHQ1.lat);
            outfile << "[" << radsToDegs(geoHQ1.lon) << "," << radsToDegs(geoHQ1.lat) << "]";
            if (j != vec.size() - 1)
                outfile << ",";
        }
        outfile << "]}}, {\"type\": \"Feature\", \"properties\" : { \"NAME\": \"point\"}, \"geometry\" : { \"type\":\"MultiPoint\",  \"coordinates\" : [";
        cellToPoint(vec[0], &geoHQ1);
        outfile << "[" << radsToDegs(geoHQ1.lon) << "," << radsToDegs(geoHQ1.lat) << "],";
        cellToPoint(vec[vec.size() - 1], &geoHQ1);
        outfile << "[" << radsToDegs(geoHQ1.lon) << "," << radsToDegs(geoHQ1.lat) << "]";
        outfile << "]}}";
    }
    outfile << "]}";
    return outfile.str();
}

double Util::calcdistance(GeoPoint g1, GeoPoint g2) {

    double a = g1.lat - g2.lat;
    double b = g1.lon - g2.lon;

    double s = 2 * asin(sqrt(pow(sin(a / 2), 2) +
        cos(g1.lat) * cos(g2.lat) * pow(sin(b / 2), 2)));
    s = s * EARTH_RADIUS;
    return s;
}

// 构建邻接关系
void Util::buildAijk(H3_N& PointMap) {
    for (auto iter = PointMap.begin(); iter != PointMap.end(); iter++) {
        int k = 1;
        int maxNeighboring = maxGridDiskSize(k);
        H3Index* neighboring = new H3Index[maxNeighboring];
        gridDisk(iter->first, k, neighboring);
        for (int i = 1; i < maxNeighboring; i++) {
            // 如果不为0则是可以点存在
            if (neighboring[i] != 0) {
                H3Index point = neighboring[i];
                if (PointMap.count(point) != 0) {
                    iter->second.add(point);
                }
            }
        }
        delete[] neighboring;
    }
}

// 构建邻接关系
void Util::buildAijk(H3_P& PointMap) {
    for (auto iter = PointMap.begin(); iter != PointMap.end(); iter++) {
        int k = 1;
        int maxNeighboring = maxGridDiskSize(k);
        H3Index* neighboring = new H3Index[maxNeighboring];
        gridDisk(iter->first, k, neighboring);
        for (int i = 1; i < maxNeighboring; i++) {
            // 如果不为0则是可以点存在
            if (neighboring[i] != 0) {
                H3Index point = neighboring[i];
                if (PointMap.count(point) != 0) {
                    iter->second.add(point);
                }
            }
        }
        delete[] neighboring;
    }
}

// 读入数据
void Util::loadingData(vector<H3_N>& PointMap, vector<H3_D>& DEM, vector<H3_D>& Comprehensive, vector<H3_V>& Factor,
    int level, string filepath) {
    // 读文件
    ifstream fp(filepath);
    string line;
    auto& m_PointMap = PointMap[level];
    auto& m_DEM = DEM[level];
    auto& m_Comprehensive = Comprehensive[level];
    auto& m_Factor = Factor[level];

    getline(fp, line); //跳过列名，第一行不做处理
    while (getline(fp, line)) { //循环读取每行数据
        string data;
        istringstream readstr(line); //string数据流化
        //将一行数据按'，'分割
        getline(readstr, data, ','); //读取数据
        H3Index index = stringToH3(data.c_str());
        getline(readstr, data, ','); //读取数据
        double water = atof(data.c_str());
        getline(readstr, data, ','); //读取数据
        double building = atof(data.c_str());
        getline(readstr, data, ','); //读取数据
        double brush = atof(data.c_str());
        getline(readstr, data, ','); //读取数据
        double forest = atof(data.c_str());
        getline(readstr, data, ','); //读取数据
        double grass = atof(data.c_str());
        getline(readstr, data, ','); //读取数据
        double soft = atof(data.c_str());
        getline(readstr, data, ','); //读取数据
        double dem = atof(data.c_str());
        getline(readstr, data, ','); //读取数据
        double comprehensive = atof(data.c_str());
        getline(readstr, data, ','); //读取数据
        double hard = atof(data.c_str());
        m_PointMap.emplace(index, Node(index));
        m_DEM.emplace(index, dem);
        m_Comprehensive.emplace(index, comprehensive);
        m_Factor.emplace(index, vector<double>{grass, water, forest, building, soft, hard, brush});
    }
    // 构建邻接关系
    buildAijk(m_PointMap);
}

// 多目标
void Util::loadingData(vector<H3_P>& PointMap, vector<H3_D>& DEM, vector<H3_D>& Comprehensive, vector<H3_V>& Factor,
    int level, string filepath) {
    // 读文件
    ifstream fp(filepath);
    string line;
    auto& m_PointMap = PointMap[level];
    auto& m_DEM = DEM[level];
    auto& m_Comprehensive = Comprehensive[level];
    auto& m_Factor = Factor[level];

    // H3Index	water	building	brush	forest	grass	soft	dem	comprehensive	TotalWay4
    getline(fp, line); //跳过列名，第一行不做处理
    while (getline(fp, line)) { //循环读取每行数据
        string data;
        istringstream readstr(line); //string数据流化
        //将一行数据按'，'分割
        getline(readstr, data, ','); //读取数据
        H3Index index = stringToH3(data.c_str());
        getline(readstr, data, ','); //读取数据
        double water = atof(data.c_str());
        getline(readstr, data, ','); //读取数据
        double building = atof(data.c_str());
        getline(readstr, data, ','); //读取数据
        double brush = atof(data.c_str());
        getline(readstr, data, ','); //读取数据
        double forest = atof(data.c_str());
        getline(readstr, data, ','); //读取数据
        double grass = atof(data.c_str());
        getline(readstr, data, ','); //读取数据
        double soft = atof(data.c_str());
        getline(readstr, data, ','); //读取数据
        double dem = atof(data.c_str());
        getline(readstr, data, ','); //读取数据
        double comprehensive = atof(data.c_str());
        getline(readstr, data, ','); //读取数据
        double hard = atof(data.c_str());
        m_PointMap.emplace(index, Point(index));
        m_DEM.emplace(index, dem);
        m_Comprehensive.emplace(index, comprehensive);
        m_Factor.emplace(index, vector<double>{grass, water, forest, building, soft, hard, brush});
    }
    // 构建邻接关系
    buildAijk(m_PointMap);
}

// 多目标
void Util::MultiObjectMultiHierarchy(vector<H3_P>& PointMap, vector<H3_D>& DEM, vector<H3_D>& Comprehensive, vector<H3_V>& Factor, int level, string filepath)
{
    ifstream fp(filepath);
    string line;
    int target_level = 16;
    // 构建基础多层次格网地图
    PointMap[target_level] = PointMap[level];
    DEM[target_level] = DEM[level];
    Comprehensive[target_level] = Comprehensive[level];
    Factor[target_level] = Factor[level];
    // 构建多层次格网地图
    level--;
    // "H3Index", "dem", "comprehensive", "grass", "water", "forest", "", "soft", "hard", "brush"
    while (getline(fp, line)) {
        int num = atoi(line.c_str());
        for (int j = 0; j < num; j++) {
            //循环读取每行数据
            getline(fp, line);
            string data;
            //string数据流化
            istringstream readstr(line);
            //将一行数据按'，'分割
            getline(readstr, data, ','); //读取数据
            H3Index index = stringToH3(data.c_str());
            getline(readstr, data, ','); //读取数据
            double water = atof(data.c_str());
            getline(readstr, data, ','); //读取数据
            double building = atof(data.c_str());
            getline(readstr, data, ','); //读取数据
            double brush = atof(data.c_str());
            getline(readstr, data, ','); //读取数据
            double forest = atof(data.c_str());
            getline(readstr, data, ','); //读取数据
            double grass = atof(data.c_str());
            getline(readstr, data, ','); //读取数据
            double soft = atof(data.c_str());
            getline(readstr, data, ','); //读取数据
            double dem = atof(data.c_str());
            getline(readstr, data, ','); //读取数据
            double comprehensive = atof(data.c_str());
            getline(readstr, data, ','); //读取数据
            double hard = atof(data.c_str());
            // 存储多层次信息
            PointMap[target_level].emplace(index, Point(index));
            DEM[target_level].emplace(index, dem);
            Comprehensive[target_level].emplace(index, comprehensive);
            Factor[target_level].emplace(index, vector<double>{grass, water, forest, building, soft, hard, brush});
            // 找到下一个层级的六角网格
            int childrenNum = cellToChildrenSize(index, level + 1);
            H3Index* neighboring = new H3Index[childrenNum];
            cellToChildren(index, level + 1, neighboring);
            // 将下一级所有的邻接全部放入邻接集合
            unordered_set<H3Index> PointSet;
            for (int k = 0; k < childrenNum; k++) {
                Point& child = PointMap[target_level][neighboring[k]];
                auto& neighbours = child.neighbours;
                PointSet.insert(neighbours.begin(), neighbours.end());
            }
            // 去掉已被替换的格网 
            // by XGG：七个小的被一个大的替代
            for (int k = 0; k < 7; k++) {
                PointSet.erase(neighboring[k]);
                PointMap[target_level].erase(neighboring[k]);
            }

            Point& nowPoint = PointMap[target_level][index];
            // 构建多层次邻接关系
            for (auto it = PointSet.begin(); it != PointSet.end(); it++) {
                Point& point = PointMap[target_level][*it];

                for (int k = 0; k < 7; k++) {
                    if (point.neighbours.count(neighboring[k]) != 0) {
                        point.neighbours.erase(neighboring[k]);// by XGG：取消孩子和邻居的相邻关系
                    }
                }
                // 小的添加大的
                point.neighbours.emplace(index);
                // 大的添加小的
                nowPoint.neighbours.emplace(*it);
            }
        }
        level--;
    }
}


// 生成多层次格网通行模型
void Util::MultiHierarchy(vector<H3_N>& PointMap, vector<H3_D>& DEM, vector<H3_D>& Comprehensive, vector<H3_V>& Factor,
    int level, string filepath) {
    ifstream fp(filepath);
    string line;
    int target_level = 16;
    // 构建基础多层次格网地图
    PointMap[target_level] = PointMap[level];
    DEM[target_level] = DEM[level];
    Comprehensive[target_level] = Comprehensive[level];
    Factor[target_level] = Factor[level];
    // 构建多层次格网地图
    level--;
    while (getline(fp, line)) {
        int num = atoi(line.c_str());
        for (int j = 0; j < num; j++) {
            //循环读取每行数据
            getline(fp, line);
            string data;
            //string数据流化
            istringstream readstr(line);
            //将一行数据按'，'分割
            getline(readstr, data, ','); //读取数据
            H3Index index = stringToH3(data.c_str());
            getline(readstr, data, ','); //读取数据
            double water = atof(data.c_str());
            getline(readstr, data, ','); //读取数据
            double building = atof(data.c_str());
            getline(readstr, data, ','); //读取数据
            double brush = atof(data.c_str());
            getline(readstr, data, ','); //读取数据
            double forest = atof(data.c_str());
            getline(readstr, data, ','); //读取数据
            double grass = atof(data.c_str());
            getline(readstr, data, ','); //读取数据
            double soft = atof(data.c_str());
            getline(readstr, data, ','); //读取数据
            double dem = atof(data.c_str());
            getline(readstr, data, ','); //读取数据
            double comprehensive = atof(data.c_str());
            getline(readstr, data, ','); //读取数据
            double hard = atof(data.c_str());
            // 单层次单独存储
          /*  PointMap[level].emplace(index, Node(index));
            buildAijk(PointMap[level]);
            DEM[level].emplace(index, dem);
            Comprehensive[level].insert(make_pair(index, speed));
            Factor[level].insert(make_pair(index, vector<double>{grass, water, forest, building, soft, hard}));*/
            // 存储多层次信息
            PointMap[target_level].emplace(index, Node(index));
            DEM[target_level].emplace(index, dem);
            Comprehensive[target_level].emplace(index, comprehensive);
            Factor[target_level].emplace(index, vector<double>{grass, water, forest, building, soft, hard, brush});
            // 找到下一个层级的六角网格
            int childrenNum = cellToChildrenSize(index, level + 1);
            H3Index* neighboring = new H3Index[childrenNum];
            cellToChildren(index, level + 1, neighboring);
            // 将下一级所有的邻接全部放入邻接集合
            set<H3Index> PointSet;
            for (int k = 0; k < childrenNum; k++) {
                Node child = PointMap[target_level][neighboring[k]];
                auto neighbours = child.neighbours;
                PointSet.insert(neighbours.begin(), neighbours.end());
            }
            // 去掉已被替换的格网
            for (int k = 0; k < 7; k++) {
                PointSet.erase(neighboring[k]);
                PointMap[target_level].erase(neighboring[k]);
                /*DEM[target_level].erase(neighboring[k]);
                Comprehensive[target_level].erase(neighboring[k]);
                Factor[target_level].erase(neighboring[k]);*/
            }

            Node& nowPoint = PointMap[target_level][index];
            // 构建多层次邻接关系
            for (auto it = PointSet.begin(); it != PointSet.end(); it++) {
                Node& point = PointMap[target_level][*it];

                for (int k = 0; k < 7; k++) {
                    if (point.neighbours.count(neighboring[k]) != 0) {
                        point.neighbours.erase(neighboring[k]);
                    }
                }
                // 小的添加大的
                point.neighbours.emplace(index);
                // 大的添加小的
                nowPoint.neighbours.emplace(*it);
            }
        }
        level--;
    }
}

// 计算地表覆盖因子占比
void Util::ComputeFactor(vector<H3Index>& result, H3_V& Factor, vector<double>& factor) {
    // vector<double>{grass, water, forest, building, soft, hard, brush};
    vector<double> temp(factor.size(), 0);
    double total = 0;
    for (int i = 0; i < result.size(); i++) {
        for (int j = 0; j < temp.size(); j++) {
            temp[j] += Factor[result[i]][j];
            total += Factor[result[i]][j];
        }
    }
    for (int i = 0; i < temp.size(); i++) {
        factor[i] += temp[i] / total * 100;
    }
}

// 计算地表覆盖因子占比
void Util::ComputeFactor(vector<pair<vector<H3Index>, point_type>>& result, H3_V& Factor, vector<vector<double>>& factor) {
    
    cout << "运行computerFactor" << endl;
    // vector<double>{grass, water, forest, building, soft, hard, brush};
    for (int i = 0; i < result.size(); i++) {
        // 针对result的每个路径 在factor后面加一个七位向量
        factor.push_back(vector<double>(7, 0));
        double total = 0;
        for (int j = 0; j < result[i].first.size(); j++) {
            // 遍历这条路径的所有格网
            // 统计这条路径的所有 factor 的和
            // 1. result[i] 代表一条路径 pair<vector<H3index>, point_type>
            // 2. result[i].first 代表 vector<H3index>
            for (int k = 0; k < 7; k++) {
                total += Factor[result[i].first[j]][k];
                factor[i][k] += Factor[result[i].first[j]][k];
            }
        }
        for (int j = 0; j < 7; j++)
            factor[i][j] /= total;
    }
}

// 求平均值 
double Util::average(vector<double>& x)
{
    int len = x.size();
    double sum = 0;
    for (int i = 0; i < len; i++) // 求和
        sum += x[i];
    return sum / len; // 得到平均值
}

// 求方差
double Util::variance(vector<double>& x)
{
    int len = x.size();
    double ave = average(x);
    double sum = 0;
    for (int i = 0; i < len; i++) // 求和
        sum += pow(x[i] - ave, 2);
    return sum / len; // 得到平均值
}

// 求标准差
double Util::standardDev(vector<double>& x)
{
    double var = variance(x);
    return sqrt(var); // 得到标准差
}

// 计算坡度平均值与标准差
void Util::ComputeGrad(vector<H3Index>& result, H3_D& DEM, double& mean, double& dev) {
    double total = 0;
    vector<double> gradList;
    for (int i = 1; i < result.size(); i++) {
        GeoPoint g1, g2;
        cellToPoint(result[i - 1], &g1);
        cellToPoint(result[i], &g2);
        double distance = calcdistance(g1, g2);
        double h = fabs(DEM[result[i]] - DEM[result[i - 1]]) / 1000;
        double s = h / distance;
        double grad = atan(s) * 180.0 / 3.1415926;
        total += grad;
        gradList.emplace_back(grad);
    }

    // 输出平均值
    mean += total / gradList.size();
    // 输出标准差
    dev += standardDev(gradList);
}

// 计算坡度平均值与标准差
void Util::ComputeGrad(vector<pair<vector<H3Index>, point_type>>& result, H3_D& DEM, vector<double>& mean, vector<double>& dev) {
    cout << "运行cmputerGrad" << endl;
    for (int j = 0; j < result.size(); j++) {
        auto& res = result[j].first;
        double total = 0;
        vector<double> gradList;
        for (int i = 1; i < res.size(); i++) {
            GeoPoint g1, g2;
            cellToPoint(res[i - 1], &g1);
            cellToPoint(res[i], &g2);
            double distance = calcdistance(g1, g2);
            double h = fabs(DEM[res[i]] - DEM[res[i - 1]]) / 1000;
            double s = h / distance;
            double grad = atan(s) * 180.0 / 3.1415926;
            total += grad;
            gradList.emplace_back(grad);
        }

        // 输出平均值
        mean.push_back(total / gradList.size());
        // 输出标准差
        dev.push_back(standardDev(gradList));
    }
}

// 多层次路径规划
vector<H3Index> Util::MultiHierarchySearch(vector<H3_N>& PointMap, vector<H3_D>& DEM, vector<H3_D>& Comprehensive, int targetLevel, int level,
    H3Index start, H3Index end, ofstream& log, int& maxSize, int& Size)
{
    Astar alg(PointMap[targetLevel], DEM[targetLevel], Comprehensive[targetLevel]);
    vector<H3Index> result = alg.search(start, end, maxSize, Size);
    //vector<vector<H3Index>> temp;
    //temp.reserve(result.size());
    //vector<thread> threadpool;
    //for (int i = 0; i < result.size();) {
    //    int res = getResolution(result[i]);
    //    if (res != level) {
    //        int index = i;
    //        H3_N m_PointMap;
    //        H3_D m_DEM;
    //        H3_D m_Comprehensive;
    //        while (index < result.size() && getResolution(result[index]) != level) {
    //            int childrenNum = cellToChildrenSize(result[index], level);
    //            H3Index* neighboring = new H3Index[childrenNum];
    //            cellToChildren(result[index], level, neighboring);
    //            for (int i = 0; i < childrenNum; i++) {
    //                m_PointMap.emplace(neighboring[i], Node(neighboring[i]));
    //                m_DEM.emplace(neighboring[i], DEM[level][neighboring[i]]);
    //                m_Comprehensive.emplace(neighboring[i], Comprehensive[level][neighboring[i]]);
    //            }
    //            index++;
    //        }
    //        // 添加起点终点
    //        H3Index start = result[i - 1];
    //        H3Index end = result[index];
    //        m_PointMap.emplace(start, Node(start));
    //        m_DEM.emplace(start, DEM[level][start]);
    //        m_Comprehensive.emplace(start, Comprehensive[level][start]);
    //        m_PointMap.emplace(end, Node(end));
    //        m_DEM.emplace(end, DEM[level][end]);
    //        m_Comprehensive.emplace(end, Comprehensive[level][end]);
    //        // 构建邻接关系
    //        buildAijk(m_PointMap);
    //        temp.push_back({});
    //        threadpool.emplace_back(thread(ParallelSearch, m_PointMap, m_DEM, m_Comprehensive, ref(temp.back()), start, end));
    //        //ParallelSearch(m_PointMap, m_DEM, m_Comprehensive, temp.back(), start, end);
    //        i = index;
    //    }
    //    else {
    //        temp.push_back({ result[i] });
    //        i++;
    //    }
    //}
    //for (int i = 0; i < threadpool.size(); i++)
    //    threadpool[i].join();
    //vector<H3Index> res;
    //for (int i = 0; i < temp.size(); i++) {
    //    if (temp[i].size() <= 1)
    //        res.emplace_back(temp[i][0]);
    //    else
    //        res.insert(res.end(), temp[i].rbegin() + 1, temp[i].rend() - 1);
    //}
    //log << "遍历格网数：" << maxSize.size() << "; ";
    return result;
}

// 并行优化
void Util::ParallelSearch(H3_N PointMap, H3_D DEM, H3_D Comprehensive, vector<H3Index>& tmp, H3Index start, H3Index end) {
    Astar AstarTmp(PointMap, DEM, Comprehensive);
    int maxSize = 0;
    int size = 0;
    vector<H3Index> vectorTmp = AstarTmp.search(start, end, maxSize, size);
    tmp = vectorTmp;
}

// 普通路径规划
vector<H3Index> Util::OrdinarySearch(H3_N& PointMap, H3_D& DEM, H3_D& Comprehensive,
    H3Index start, H3Index end, ofstream& log, int& maxSize, int& Size)
{
    Astar alg(PointMap, DEM, Comprehensive);
    vector<H3Index> result = alg.search(start, end, maxSize, Size);
    return result;
}

// 多目标路径规划
vector<pair<vector<H3Index>, point_type>> Util::MultiObjectMultiHierarchySearch(vector<H3_P>& PointMap, vector<H3_D>& DEM, vector<H3_D>& Comprehensive, int targetLevel, int level, H3Index start, H3Index end,
    ofstream& log, int& maxSize, int& Size, double& time) {
    ofstream log_cost;

    cout << "正在运行MultiObjectMultiHierarchySearch" << endl;
    initializer_list<bool> direction_comparison = { true, true};
    Multi_Astar alg(PointMap[targetLevel], direction_comparison, DEM[targetLevel], Comprehensive[targetLevel], level);
    clock_t st, fin;
    st = clock();
    // maxsize Max_of_the_Number_of_Open_Grid open队列的最大深度
    // size Traverses_the_Number_of_Grid 遍历的格网数
    auto result = alg.search(start, end, maxSize, Size);
    fin = clock();
    time = (double)(fin - st) / CLOCKS_PER_SEC;
    cout << "大小" << result.size() << endl;
    for (auto it = result.begin(); it != result.end(); it++)
    {
        log_cost.open("D:/桌面/costs.txt", ios::app);
        log_cost << "路径 目标1:" << *it->first.begin() << endl;
        log_cost << "路径 目标2:" << *(it->first.end() -1) << endl;
    }
    // 注意起点和终点颠倒，因为需要倒序得到路径
    auto rnt = getPathList(PointMap[targetLevel], DEM[targetLevel], Comprehensive[targetLevel], end, start);
    return rnt;
}

void Util::AnalysisData(H3_N& PointMap, H3_D& DEM, H3_D& Comprehensive, vector<H3Index>& path, double& Length_of_Path) {
    // 计算个数和距离
    double Length = 0;
    for (int i = 1; i < path.size(); i++) {
        H3Index start = path[i];
        H3Index end = path[i - 1];
        GeoPoint g1, g2;
        cellToPoint(start, &g1);
        cellToPoint(end, &g2);
        double distance = Util::calcdistance(g1, g2);
        Length += distance;
    }
    // 输出路径实际长度
    Length_of_Path += Length * 1000;
}

void Util::AnalysisData(H3_P& PointMap, H3_D& DEM, H3_D& Comprehensive, vector<pair<vector<H3Index>, point_type>>& paths, vector<double> & Length_of_Path) {
    cout << "正在运行AnalysisData" << endl;
    // 计算个数和距离
    for (int i = 0; i < paths.size(); i++) {
        auto path = paths[i].first;
        double Length = 0;
        for (int j = 1; j < path.size(); j++) {
            H3Index start = path[j];
            H3Index end = path[j - 1];
            GeoPoint g1, g2;
            cellToPoint(start, &g1);
            cellToPoint(end, &g2);
            double distance = Util::calcdistance(g1, g2);
            Length += distance;
        }
        Length *= 1000;
        Length_of_Path.push_back(Length);
    }
}

//倒序得到路径，此时的start为初始终点，end为初始起点
vector<pair<vector<H3Index>, point_type>> Util::getPathList(H3_P& NodeMap, H3_D& DEM, H3_D& Comprehensive, H3Index start, H3Index end) {
    vector<pair<vector<H3Index>, point_type>> res;
    auto paths = NodeMap[start].cl;
    // para g_m，包含父节点的Parameter
    for (auto& [k, v] : paths) {
        vector<H3Index> vec;
        // current里保存的是父节点的index
        auto current = v;
        // 当前节点的index用next表示
        H3Index next = start;
        while (current.getIndex() != 0) {
            vec.push_back(next);
            auto g_m = calEdgeCost(DEM, Comprehensive, current, next);
            auto nextNode = NodeMap[current.getIndex()].cl;
            next = current.getIndex();
            current = nextNode[g_m];
        }
        vec.push_back(next);
        res.emplace_back(vec, k);
    }
    return res;
}

point_type Util::calEdgeCost(H3_D& DEM, H3_D& Comprehensive, Parameter parent, H3Index current) {
    point_type res;
    // 第一个分量：距离
    GeoPoint g1, g2;
    cellToPoint(current, &g1);
    cellToPoint(parent.getIndex(), &g2);
    double distance = Util::round(Util::calcdistance(g1, g2), 3);
    double h = fabs(DEM[current] - DEM[parent.getIndex()]);
    // 计算实际距离
    // int real_distance = Util::round(sqrt(pow(h, 2) + pow(distance, 2)), 0);
    // 第二个分量：坡度
    double s = h / distance;
    double grad = atan(s) * 180.0 / 3.1415926;
    // 这里使用1-grad，越小越好
    grad = -0.00004061 * pow(grad, 3) + 0.002644 * pow(grad, 2) - 0.07607 * grad + 0.9977;

    // res[0] = Util::round((1 - Comprehensive[current]), 2);
    // res[1] = parent_g[1] + Util::round((1 - Comprehensive[current]), 2);
    //res[0] = distance;
    // res[1] = parent_g[1] + Util::round((1 - grad), 2);
    //res[1] = calGridNum(parent, current) * 100;
    //[1] = distance;
    //res[1] = 0;
    parent.distance -= distance;
    parent.total_comprehensive -= Comprehensive[current];
    parent.total_grad -= grad;
    parent.grid_num--;
    parent.avg_grad = parent.total_grad / parent.grid_num;
    parent.avg_comprehensive = parent.total_comprehensive / parent.grid_num;

    res[0] = Util::round(parent.distance / parent.avg_comprehensive, 0);
    res[1] = Util::round(parent.distance / parent.avg_grad, 0);
    // 返回值
    return res;
}

int Util::calGridNum(H3Index p1, H3Index p2)
{
    int gridDistance = abs(getResolution(p1) - getResolution(p2));
    if (gridDistance == 0)
        return 1;
    else
        return 2 + (gridDistance - 1) * 3;
}

// lat,lng为弧度表示的经纬度，r为地球半径，由于是算夹角，r是多少不重要
// 返回一个x，y，z三元组
inline vector<double> Util::ball2xyz(double lat, double lng) {
    return { EARTH_RADIUS * cos(lat) * cos(lng), EARTH_RADIUS * cos(lat) * sin(lng), EARTH_RADIUS * sin(lat) };
}

// 将地理经纬度转换成笛卡尔坐标系
// 为角度表示的经纬度点
inline vector<double> Util::geo2xyz(GeoPoint& point) {
    return ball2xyz(point.lat, point.lon);
}

// 计算3个地理坐标点之间的夹角
double Util::angleOflocation(GeoPoint l1, GeoPoint l2, GeoPoint l3) {
    auto p1 = geo2xyz(l1);
    auto p2 = geo2xyz(l2);
    auto p3 = geo2xyz(l3);

    double x1 = p1[0];
    double x2 = p2[0];
    double x3 = p3[0];

    double y1 = p1[1];
    double y2 = p2[1];
    double y3 = p3[1];

    double z1 = p1[2];
    double z2 = p2[2];
    double z3 = p3[2];

    // 计算向量 P2P1 和 P2P3 的夹角
    auto _P1P2 = sqrt(pow((x2 - x1), 2) + pow((y2 - y1), 2) + pow((z2 - z1), 2));
    auto _P2P3 = sqrt(pow((x3 - x2), 2) + pow((y3 - y2), 2) + pow((z3 - z2), 2));

    auto P = (x1 - x2) * (x3 - x2) + (y1 - y2) * (y3 - y2) + (z1 - z2) * (z3 - z2); //P2P1*P2P3

    return (acos(P / (_P1P2 * _P2P3)) / PI) * 180;
}

void Util::loading_site_data(string filepath, H3_D &cellList)
{
    // 1. 参数1： 文件路径
    // 2. 参数2： data 存放所有站点的H3index
    // 3. 参数3： cellList 存放所有站点的键值对
    // 把站点的那个文件读到 site中去
    ifstream fp(filepath);
    string line;
    getline(fp, line); //跳过列名，第一行不做处理
    while (getline(fp, line)) { //循环读取每行数据
        string data;
        istringstream readstr(line); //string数据流化
        getline(readstr, data, ',');
        getline(readstr, data, ',');
        H3Index index = stringToH3(data.c_str());
        getline(readstr, data, ',');
        double is_site = atof(data.c_str());
        cellList.emplace(index, is_site);
        ////将一行数据按'，'分割
        //getline(readstr, data, ','); //读取数据
        //H3Index index = stringToH3(data.c_str());
        //getline(readstr, data, ','); //读取数据
        //double water = atof(data.c_str());
        //getline(readstr, data, ','); //读取数据
        //double building = atof(data.c_str());
        //getline(readstr, data, ','); //读取数据
        //double brush = atof(data.c_str());
        //getline(readstr, data, ','); //读取数据
        //double forest = atof(data.c_str());
        //getline(readstr, data, ','); //读取数据
        //double grass = atof(data.c_str());
        //getline(readstr, data, ','); //读取数据
        //double soft = atof(data.c_str());
        //getline(readstr, data, ','); //读取数据
        //double dem = atof(data.c_str());
        //getline(readstr, data, ','); //读取数据
        //double comprehensive = atof(data.c_str());
        //getline(readstr, data, ','); //读取数据
        //double hard = atof(data.c_str());
        /*m_PointMap.emplace(index, Point(index));
        m_DEM.emplace(index, dem);
        m_Comprehensive.emplace(index, comprehensive);
        m_Factor.emplace(index, vector<double>{grass, water, forest, building, soft, hard, brush});*/
    }
}

void Util::SetNumber(H3_D& siteList)
{
    // 1. sitelist的形式 h3_D 是一个unordered_map
    // h3index -- 1 或者 0


}

void Util::callPython(vector<double> parList)
{

}

string Util::PrintCurrentTime()
{
    auto now = chrono::system_clock::now();
    // 转换为 std::time_t 类型
    auto time = chrono::system_clock::to_time_t(now);
    // 打印时间
    return ctime(&time);
}

int Util::cnt = 0;