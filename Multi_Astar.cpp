#include "Multi_Astar.h"
#include "Util.h"

Multi_Astar::Multi_Astar(H3_P& m, std::initializer_list<bool> list, H3_D& dem, H3_D& comprehensive, int level) :NodeMap(m), DEM(dem), Comprehensive(comprehensive) {
    direction_comparison = list;
    baseLevel = level;
    COSTS = front_type(direction_comparison);
}

Multi_Astar::~Multi_Astar()
{

}

bool CompWithoutDirectionComparison(const pair<array_type, pair<array_type, Parameter>>& a,
    const pair<array_type, pair<array_type, Parameter>>& b) {
    // 这里设置优先级 F_m
    // < 升序 > 降序
    array_type arr_a = a.first;
    array_type arr_b = b.first;
    // 距离为第一个
    if (arr_a[0] == arr_b[0])
        return arr_a[1] < arr_b[1];// by xgg 第一个相同，就按第二个排
    else
        return arr_a[0] < arr_b[0];// by xgg 第一个不相同，就按第一个排
    //通行能力第一个
    //if (arr_a[1] == arr_b[1])
    //    return arr_a[0] < arr_b[0];
    //else
    //    return arr_a[1] < arr_b[1];
 /*   if (arr_a[0] == arr_b[0]) {
        if(arr_a[1] == arr_b[1])
            return arr_a[2] < arr_b[2];
        else
            return arr_a[1] < arr_b[1];
    } else
        return arr_a[0] < arr_b[0];*/
}

front_type Multi_Astar::search(H3Index startPos, H3Index endPos, int& maxSize, int& Size)
{
    list<H3Index> path;
    front_type pathlist;
    // 判断起始点和终止点是否存在
    if (NodeMap.count(startPos) == 0 || NodeMap.count(endPos) == 0)
        return pathlist;
    // 存入S点
    //array_type g{ 0, 0, 0};
    array_type g{ 0, 0};
    array_type h = calH(endPos, startPos);
    NodeMap[startPos].op.insert(make_pair(point_type(g.begin(), g.end()), Parameter{ 0 , Comprehensive[startPos]}));
    // 计算F_m:
    array_type F;
    for (int i = 0; i < VectorDimension; i++)
        F[i] = g[i] + h[i];
    // G,F,Parameter
    OPEN.insert(make_pair(g, make_pair(F, Parameter{ startPos, Comprehensive[startPos]})));
    while (!OPEN.empty())
    {
        // 获得点和它的G_M
        //array_type g_m({ 0, 0, 0});
        array_type g_m({ 0, 0 });
        // 从OPEN中取最优的非主宰向量
        maxSize = max(maxSize, (int)OPEN.size());
        // 格网数
        int grid_num = 0;
        // 弹出最佳F和其参数
        Parameter para;
        Point& current = OpenPop(g_m, para);
        //cout << "current Point " << current.getIndex() << "(" << g_m[0] << "," << g_m[1] << ") ;;;;;;;;;;;PriorityQueueSize:" << OPEN.size() << ";;;;;;;;result : " << COSTS.size() << endl;
        //cout << "; PriorityQueueSize:" << OPEN.size() << ";;;;;;;; result: " << COSTS.size() << endl;
        if (current.getIndex() == endPos)
        {
            // 因为终点的H为0，所以终点的F_m就等于G_m
            point_type point_g_m(g_m.begin(), g_m.end());
            // 记录路径
            // 
            COSTS.insert(make_pair(point_g_m, para));
            cout << "路径数：" << COSTS.size() << endl;
            // 清除在OPEN中被支配的数据
            decltype(OPEN) temp;
            for (auto& [p, v] : OPEN) {
                point_type point_f_m(v.first.begin(), v.first.end());
                if (point_g_m.dominates(point_f_m))
                    continue;
                temp.insert(make_pair(p, v));
            }
            OPEN = move(temp);
        }
        else
            NextStep(current, g_m, endPos, Size, para);
        //cout << "current Point " << current.getIndex() << "(" << g_m[0] << "," << g_m[1] << ");;;;;;;;;;;PriorityQueueSize:" << OPEN.size() << ";;;;;;;;result : " << COSTS.size() << endl;
    }
    //cout << "Total Result : " << COSTS.size() << endl;
    return COSTS;
}

// 弹出在非支配向量里面最优的一个
// 先将Open里的三元组全部放入非支配集合
// 再对非支配集合进行排序
// 再选取最优的一个
// 参数：第一个参数为G，第二个参数为从起点到当前的格网个数，主要是为了求出平均通行能力
Point& Multi_Astar::OpenPop(array_type& g_m, Parameter& parameter) {
    // OPEN g,F,Parameter
    // F,g,Parameter
    pareto::front<double, VectorDimension, pair<array_type, Parameter>> nondominantSet;
    // 遍历OPEN列表,并且插入到非支配向量集。
    for (auto it = OPEN.begin(); it != OPEN.end(); it++)
        nondominantSet.insert(make_pair(point_type(it->second.first.begin(), it->second.first.end()), make_pair(it->first, it->second.second)));
    // 按照不同分量的优先级排序
    // F,g,index
    vector<pair<array_type, pair<array_type, Parameter>>> nondominantVector;
    for (auto& [key, value] : nondominantSet)
        nondominantVector.push_back(make_pair(key.values(), value));
    // 按F分量来排序，升序
    sort(nondominantVector.begin(), nondominantVector.end(), CompWithoutDirectionComparison);
    // 取最优（最小）的代价向量，也就是第一个，这里是按g来定位
    auto target = nondominantVector[0];
    auto end = OPEN.upper_bound(target.second.first);
    auto it = OPEN.find(target.second.first);
    // g,F,Parameter
    while (it != end) {
        if (it->first == target.second.first && it->second.first == target.first && it->second.second == target.second.second)
            break;
        it++;
    }
    // 依据H3index 找到那个点
    assert(it != end);
    // 输出 index F[0] F[1]
    //cout << "0x" << std::hex << it->second.second << std::dec << "(" << (it->second.first)[0] << "," << (it->second.first)[1] << ")";
    // 找到这个格网
    Point& p = NodeMap[it->second.second.getIndex()];
    //将这个点的g_m从Gop移动到Gcl
    point_type point_g_m(it->first.begin(), it->first.end());
    // g_m, 父节点;
    assert(p.op.count(point_g_m) != 0);
    p.cl.insert(make_pair(point_g_m, p.op[point_g_m]));
    p.op.erase(point_g_m);
    g_m = it->first;
    parameter = it->second.second;
    //在Open列表中删除被取的三元组
    OPEN.erase(it);
    return p;
}

void Multi_Astar::NextStep(Point& current, array_type current_g, H3Index goal, int& Size, Parameter& parameter)
{
    // 遍历邻接点
    // it == H3Index;
    for (auto it = current.neighbours.begin(); it != current.neighbours.end(); it++) {
        // 获得邻接点
        Point& m = NodeMap[*it];
        if (current.getIndex() == *it)
            continue;
        // ============================================================================
        // 计算g_m:
        bool flag = true;
        auto para = parameter;
        array_type g_m = calG(current_g, current.getIndex(), *it, flag, para);
        // 需要排查邻接点是否可通行
        if (!flag)
            continue;
        Size++;
        // 计算h_m:
        array_type h_m = calH(goal, m.getIndex(), current.getIndex());
        // 计算F_m:
        array_type F_m;
        for (int i = 0; i < VectorDimension; i++)
            F_m[i] = g_m[i] + h_m[i];
        // ============================================================================
        point_type point_F_m(F_m.begin(), F_m.end());
        point_type point_g_m(g_m.begin(), g_m.end());
        // 如果m是新点
        if (m.op.empty() && m.cl.empty()) {
            // 如果F_m被COSTS支配，那么跳过它
            if (COSTS.dominates(point_F_m))
                continue;
            // 放入OPEN列表
            para.setIndex(m.getIndex());
            OPEN.insert(make_pair(g_m, make_pair(point_F_m.values(), para)));
            // 添加Gop_m,虽然插入的时候会判断是否主宰，但是这个时候op为空，不影响
            para.setIndex(current.getIndex());
            m.op.insert(make_pair(point_g_m, para));
        }
        // 如果不是新点
        // 如果g_m与op与cl中已存在的代价向量相等
        else if (m.op.contains(point_g_m)) {
            para.setIndex(current.getIndex());
            m.op.find(point_g_m)->second = para;
        }
        else if (m.cl.contains(point_g_m)) {
            para.setIndex(current.getIndex());
            m.cl.find(point_g_m)->second = para;
        }
        // 如果不是新点也不相等，那么需要满足
        // g_m不被op与cl中其他向量支配
        else if (!m.op.dominates(point_g_m) && !m.cl.dominates(point_g_m)) {
            //从op中删除被g_m支配的g'_m，并一同删除OPEN中的
            front_type temp;
            for (auto it = m.op.find_dominated(point_g_m); it != m.op.end(); ++it) {
                OPEN.erase(it->first.values());
                temp.insert(make_pair(it->first, it->second));
            }
            for (auto& [p, v] : temp) {
                m.op.erase(p);
            }
            if (m.getIndex() == 0x08a4008c40a0ffff) {
                cout << " ";
            }
            //从m的cl中删除被g_m支配的g'_m
            temp.clear();
            for (auto& [p, v] : m.cl) {
                if (point_g_m.dominates(p))
                    continue;
                temp.insert(make_pair(p, v));
            }
            m.cl = temp;
            //检查F_m是否被支配
            if (COSTS.dominates(point_F_m))
                continue;
            else {
                // 如果不支配，将当前点放入OPEN列表和op列表
                // g_m,F_m,index为当前节点的Parameter
                para.setIndex(m.getIndex());
                OPEN.insert(make_pair(g_m, make_pair(F_m, para)));
                // g_m,index为父节点的Parameter
                para.setIndex(current.getIndex());
                m.op.insert(make_pair(point_g_m, para));
            }
        }
    }
}

array_type Multi_Astar::calH(H3Index goal, H3Index current) {
    if (current == goal)
        return array_type({ 0,0 });// by Xgg 强制转换？
    else {
        array_type res;
        // 第一个分量：距离
        GeoPoint g1, g2;
        cellToPoint(current, &g1);
        cellToPoint(goal, &g2);
        double distance = Util::round(Util::calcdistance(g1, g2), 3);
        double h = fabs(DEM[current] - DEM[goal]);
        // 计算当前点到终点格网个数
        double edge = getHexagonEdgeLengthAvgM(getResolution(goal));
        int grid = distance / (1.732050 * edge) + 1;
        // 第二个分量：坡度与通行能力
        double s = h / distance;
        double grad = atan(s) * 180.0 / 3.1415926;
        // 这里使用1-grad，坡度越小量化的值越小
        grad = -0.00004061 * pow(grad, 3) + 0.002644 * pow(grad, 2) - 0.07607 * grad + 0.9977;
        if (grad > 1)
            grad = 1;
        if (grad < 0)
            grad = 0;
        // 问题在这，可能会出现反方向通行能力强的格网，所以这里得加上方向性
        // res[0] = distance / 100;
        // res[1] = Util::round((2 - grad - Comprehensive[current]) * distance, 2);
        // res[1] = 0;
        // res[1] = distance;
        // res[1] = int((1 - Comprehensive[current]) * grid);
        // res[1] = grid * 100 - Util::round(grad, 2);
        res[0] = Util::round(distance / 10, 1);
        res[1] = Util::round(distance / 10, 1);
        // 返回值
        return res;
    }
}

// 正在使用
array_type Multi_Astar::calH(H3Index goal, H3Index current, H3Index parent) {
    if (current == goal)
        return array_type({ 0,0});
    else {
        array_type res;
        // 第一个分量：距离
        GeoPoint g1, g2, g3;
        cellToPoint(current, &g1);
        cellToPoint(goal, &g2);
        cellToPoint(parent, &g3);
        double distance = Util::round(Util::calcdistance(g1, g2), 3);
        double h = fabs(DEM[current] - DEM[goal]);
        // 计算当前点到终点格网个数
        double edge = getHexagonEdgeLengthAvgM(getResolution(goal));
        int grid = distance / (1.732050 * edge) + 1;
        // 第二个分量：坡度与通行能力
        double s = h / distance;
        double grad = atan(s) * 180.0 / 3.1415926;
        // 这里使用1-grad，坡度越小量化的值越小
        grad = -0.00004061 * pow(grad, 3) + 0.002644 * pow(grad, 2) - 0.07607 * grad + 0.9977;
        if (grad > 1)
            grad = 1;
        if (grad < 0)
            grad = 0;
        // 问题在这，可能会出现反方向通行能力强的格网，所以这里得加上方向性
        //double weight = cos(degsToRads(Util::angleOflocation(g2, g3, g1))) + 1;
        //res[0] = Util::round(distance - (weight * Comprehensive[current]) * 100, 0);
        //res[1] = Util::round(distance - (weight * grad) * 100, 0);
        //res[0] = Util::round(distance - (1 - Comprehensive[current]) * 100, 0);
        //res[1] = Util::round(distance - (1 - grad) * 100, 0);
        //double h1 = distance - (weight + Comprehensive[current]) * 100;
        //double h2 = distance - (weight + grad) * 100;
        //res[0] = Util::round(h1 < 0 ? 0 : h1, 0);
        //res[1] = Util::round(h2 < 0 ? 0 : h2, 0);
        
        
        // res[0] = Util::round(1 - weight * Comprehensive[current], 2);
        // res[1] = Util::round((2 - grad - Comprehensive[current]) * distance, 2);
        // res[1] = 0;
        // res[1] = int((1 - Comprehensive[current]) * grid);
        res[0] = distance;
        res[1] = distance;
        // 返回值
        return  res;
    }
}

// flag为判断是否可通行，parameter为携带的参数
array_type Multi_Astar::calG(array_type& parent_g, H3Index parent, H3Index current, bool& flag, Parameter& parameter) {
    array_type res = { 0,0 };
    if (current == parent)
        return res;
    else {
        // 第一个分量：距离
        GeoPoint g1, g2;
        cellToPoint(current, &g1);
        cellToPoint(parent, &g2);
        double distance = Util::round(Util::calcdistance(g1, g2), 3);
        double h = fabs(DEM[current] - DEM[parent]);
        // 计算实际距离
        // int real_distance = Util::round(sqrt(pow(h, 2) + pow(distance, 2)), 0);
        // 第二个分量：坡度
        double s = h / distance;
        double grad = atan(s) * 180.0 / 3.1415926;
        // 这里使用1-grad，越小越好
        grad = -0.00004061 * pow(grad, 3) + 0.002644 * pow(grad, 2) - 0.07607 * grad + 0.9977;
        // 判断是否可通行
        if (Comprehensive[current] <= 0) {
            flag = false;
            return res;
        }
        // 判断坡度是否可通行
        if (grad > 1 || grad < 0) {
            flag = false;
            return res;
        }
        parameter.grid_num++;
        parameter.total_comprehensive += Comprehensive[current];
        parameter.total_grad += grad;
        parameter.avg_grad = parameter.total_grad / parameter.grid_num;
        parameter.avg_comprehensive = parameter.total_comprehensive / parameter.grid_num;
        parameter.distance += distance;
        //res[0] = parent_g[0] + distance;
        // res[1] = parent_g[1] + Util::calGridNum(parent, current) * 100;
        //res[1] = parent_g[1] + distance;
        res[0] = Util::round(parameter.distance / parameter.avg_comprehensive, 0);
        res[1] = Util::round(parameter.distance / parameter.avg_grad, 0);
        // 返回值
        return res;
    }
}

vector<pair<list<H3Index>, array_type>> Multi_Astar::printPathByNotRecursion(H3Index end, H3Index start) {

    ////主栈，用于计算处理路径
    //stack<H3Index> majorStack;
    ////副栈，用于存储待处理节点
    //stack<H3Index> minorStack;
    ////用于存储一个节点的子节点的弹出个数
    //map<H3Index, int> popCount;
    ////根节点入副栈
    //minorStack.push(start);
    //// 储存结果
    //map<list<H3Index>, array_type> pathList;
    //list<H3Index> currentList;
    //array_type cost;
    //cost.fill(0);
    //while (!minorStack.empty()) {
    //    // 出副栈，入主栈
    //    // 要处理的节点
    //    H3Index minLast = minorStack.top();
    //    minorStack.pop();
    //    majorStack.push(minLast);
    //    // 添加路径代价向量
    //    if (!currentList.empty()) {
    //        Point& minLastPoint = NodeMap[minLast];
    //        for (int i = 0; i < Multi_Astar::length; i++) {
    //            cost[i] += minLastPoint.children[currentList.back()][i];
    //        }
    //    }
    //    // 记录当前节点
    //    currentList.push_back(minLast);
    //    if (minLast == stringToH3("8931852C5A7FFFF")) {
    //        cout << "  " << endl;
    //    }
    //    // 将该节点的子节点入副栈
    //    Point& currentPoint = NodeMap[minLast];
    //    if (currentPoint.cl.size() != 0 && minLast != end) {
    //        for (auto it = currentPoint.cl.begin(); it != currentPoint.cl.end(); it++) {
    //            minorStack.push(it->second);
    //            // 如果重新添加进栈，则清空标记
    //            // popCount.erase(it->second);
    //        }
    //    }
    //    //栈顶元素
    //    H3Index majLast = majorStack.top();
    //    //循环条件：栈顶为叶子节点 或 栈顶节点孩子节点遍历完了
    //    while (majLast == end || NodeMap[majLast].cl.size() == 0 ||
    //        (popCount.count(majLast) != 0 && popCount[majLast] == NodeMap[majLast].cl.size())) {
    //        H3Index last = majorStack.top();
    //        majorStack.pop();
    //        //此时主栈为空，运算完毕
    //        if (majorStack.empty()) {
    //            return pathList;
    //        }
    //        majLast = majorStack.top();
    //        //第一次弹出孩子节点，弹出次数设为1
    //        if (popCount.count(majLast) == 0) {
    //            popCount.insert(make_pair(majLast, 1));
    //        }
    //        //非第一次弹出孩子节点，在原有基础上加1
    //        else {
    //            popCount[majLast] += 1;
    //        }
    //        //如果是叶子节点才将结果加入路径集中
    //        if (last == end) {
    //            pathList.insert(make_pair(currentList, cost));
    //        }
    //        // 减去代价值
    //        Point& lastPoint = NodeMap[last];
    //        for (int i = 0; i < Multi_Astar::length; i++) {
    //            cost[i] -= lastPoint.children[majLast][i];
    //        }
    //        // 弹出点
    //        currentList.pop_back();
    //    }
    //}
    //return pathList;
    return vector<pair<list<H3Index>, array_type>>();
}

void Multi_Astar::printPath(Point& current, Point& start, Point& end, list<H3Index>& path, array_type& cost, vector<pair<list<H3Index>, array_type>>& pathlist)
{
    /*   path.push_back(current.getIndex());
       if (current == start)
       {
           if (pathlist.size() > 10)
               return;
           pathlist.insert(make_pair(path, cost));
       }
       else {
           for (auto it = current.cl.begin(); it != current.cl.end(); it++) {
               Point& nextPoint = NodeMap[it->second];
               for (int i = 0; i < Multi_Astar::length; i++) {
                   cost[i] += nextPoint.children[current.getIndex()][i];
               }
               printPath(nextPoint, start, end, path, cost, pathlist);
               for (int i = 0; i < Multi_Astar::length; i++) {
                   cost[i] -= nextPoint.children[current.getIndex()][i];
               }
           }
       }
       path.pop_back();*/
}