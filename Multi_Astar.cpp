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

bool compareUseCD(const pair<array_type, pair<pair<array_type, Parameter>, double>>& a,
    pair<array_type, pair<pair<array_type, Parameter>, double>>& b) {
    //comp(a,b)函数的返回值是一个bool值，当返回值为true时不改变元素顺序，反之则需要调换元素。
    //可以把其中的a看作序列中前一个位置的元素，b看作后一个位置的元素:
    //如果a < b的时候comp(a, b) = true，那么a就会被放在b前面，排序呈升序。

    double dou_a = a.second.second;
    double dou_b = b.second.second;

    return dou_a < dou_b;
}
bool compareUseCD2(const pair<double, pair<array_type,pair<array_type, Parameter>>>& a,
    const pair<double, pair<array_type, pair<array_type, Parameter>>>& b) {
    //comp(a,b)函数的返回值是一个bool值，当返回值为true时不改变元素顺序，反之则需要调换元素。
    //可以把其中的a看作序列中前一个位置的元素，b看作后一个位置的元素:
    //如果a < b的时候comp(a, b) = true，那么a就会被放在b前面，排序呈升序。

    double dou_a = a.first;
    double dou_b = b.first;

    return dou_a > dou_b;
}



bool CompareUseDr(const pair<array_type, pair<array_type, Parameter>>& a,
    const pair<array_type, pair<array_type, Parameter>>& b) {
    // 这里设置优先级 F_m
    // < 升序 > 降序
    array_type arr_a = a.second.first;
    array_type arr_b = b.second.first;
    // 距离为第一个
    if (arr_a[0] == arr_b[0])
        return arr_a[1] < arr_b[1];// by xgg 第一个相同，就按第二个排
    else
        return arr_a[0] < arr_b[0];// by xgg 第一个不相同，就按第一个排

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

front_type Multi_Astar::search2(H3Index startPos, H3Index endPos, int& maxSize, int& Size)
{
    // maxsize Max_of_the_Number_of_Open_Grid open队列的最大深度
    // size Traverses_the_Number_of_Grid 遍历的格网数


    // 1. 迭代开始之前初始化startpos，将其加入到open列表中
    // 2. 开始执行迭代，每次从open中弹出来一个三元组：
    //      2.1 如果这个点是终点，把这条路的成本加入到COSTS中，并在OPEN中删除被这个成本支配的三元组
    //      2.2 如果这个点不是终点 执行nextstep
    int count = 0;
    //double g_2_min = DBL_MAX;
    list<H3Index> path;
    front_type pathlist;
    // 判断起始点和终止点是否存在
    if (NodeMap.count(startPos) == 0 || NodeMap.count(endPos) == 0)
        return pathlist;
    // 存入S点
    //array_type g{ 0, 0, 0};
    array_type g{ 0, 0 };
    array_type h = calH(endPos, startPos);
    NodeMap[startPos].op.insert(make_pair(point_type(g.begin(), g.end()), Parameter{ 0 , Comprehensive[startPos] }));
    // 计算F_m:
    array_type F;
    for (int i = 0; i < VectorDimension; i++)
        F[i] = g[i] + h[i];
    // G,F,Parameter
    OPEN.insert(make_pair(g, make_pair(F, Parameter{ startPos, Comprehensive[startPos] })));

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
        //Point& current = OpenPop(g_m, para);
        Point& current = OpenPop2(g_m, para, count);
        //Point& current = SelectSPEA2(g_m, para);
        //Point& current = SelectCrowdDisOPEN(g_m, para);
        //Point& current = SelectSingleSort(g_m, para);
        //cout << "current Point " << current.getIndex() << "(" << g_m[0] << "," << g_m[1] << ") ;;;;;;;;;;;PriorityQueueSize:" << OPEN.size() << ";;;;;;;;result : " << COSTS.size() << endl;
        //cout << "; PriorityQueueSize:" << OPEN.size() << ";;;;;;;; result: " << COSTS.size() << endl;

        // current g,
        // 获取这个点 Point& p = NodeMap[endPos]
        if (NodeMap[endPos].g_2_min <= g_m[1])
            continue;

        // 如果弹出的格网是终点
        if (current.getIndex() == endPos)
        {
            NodeMap[endPos].g_2_min = g_m[1];
            // 这个文件用来统计遍历的格网数
            ofstream log_size;
            // 因为终点的H为0，所以终点的F_m就等于G_m
            point_type point_g_m(g_m.begin(), g_m.end());
            // 记录路径
            // 
            COSTS.insert(make_pair(point_g_m, para));
            if (COSTS.size() == 1) {
                log_size.open("D:/桌面/size_and_maxSize0927.txt", ios::app);
                /*cout << "找到第一条路径之前遍历的格网数" << Util::cnt << endl;
                cout << "找到第一条路径之前遍历的格网数大" << Size << endl;
                log_size << "找到第一条路径之前遍历的格网数" << Util::cnt << endl;
                cout << "找到第一条路径之前最大OPEN队列深度" << maxSize << endl;;
                log_size << "找到第一条路径之前最大OPEN队列深度" << maxSize << endl;*/

                /*log_size.close();
                cout << "cnt" << Util::cnt << endl;
                cout << "count" << count << endl;*/

                cout << "wbb遍历格网数Size = " << Size << endl;
                cout << "改进后的cnt = " << Util::cnt << endl;
                cout << "找到第一条路的最大OPEN队列深度" << maxSize << endl;

                log_size << "时间戳： " << Util::PrintCurrentTime();
                log_size << "wbb遍历格网数Size = " << Size << endl;
                log_size << "改进后的cnt = " << Util::cnt << endl;
                log_size << "找到第一条路的最大OPEN队列深度" << maxSize << endl;
                log_size.close();
            }
            cout << "路径数：" << COSTS.size() << endl;
            continue;
        }


#pragma region 弹出的点不是终点，则需要遍历邻居
        else
            // 参数说明：
            // 1. current 当前弹出的格网
            // 2. 当前的g
            // 3. endpos终点
            // 4. para 弹出格网的para
            NodeMap[current.getIndex()].g_2_min = g_m[1];
            NextStep2(current, g_m, endPos, Size, para);
        //cout << "current Point " << current.getIndex() << "(" << g_m[0] << "," << g_m[1] << ");;;;;;;;;;;PriorityQueueSize:" << OPEN.size() << ";;;;;;;;result : " << COSTS.size() << endl;
    }
#pragma endregion 
    //cout << "Total Result : " << COSTS.size() << endl;
    return COSTS;
}

front_type Multi_Astar::search(H3Index startPos, H3Index endPos, int& maxSize, int& Size)
{
    // maxsize Max_of_the_Number_of_Open_Grid open队列的最大深度
    // size Traverses_the_Number_of_Grid 遍历的格网数
    

    // 1. 迭代开始之前初始化startpos，将其加入到open列表中
    // 2. 开始执行迭代，每次从open中弹出来一个三元组：
    //      2.1 如果这个点是终点，把这条路的成本加入到COSTS中，并在OPEN中删除被这个成本支配的三元组
    //      2.2 如果这个点不是终点 执行nextstep
    int count = 0;
    //double g_2_min = DBL_MAX;
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
        //Point& current = OpenPop(g_m, para);
        Point& current = OpenPop(g_m, para,count);
        //Point& current = SelectSPEA2(g_m, para);
        //Point& current = SelectCrowdDisOPEN(g_m, para);
        //Point& current = SelectSingleSort(g_m, para);
        //cout << "current Point " << current.getIndex() << "(" << g_m[0] << "," << g_m[1] << ") ;;;;;;;;;;;PriorityQueueSize:" << OPEN.size() << ";;;;;;;;result : " << COSTS.size() << endl;
        //cout << "; PriorityQueueSize:" << OPEN.size() << ";;;;;;;; result: " << COSTS.size() << endl;
//#pragma region g_2_min(Sgoal) <= g_current_2 把s消除
//        // current g,
//        if (g_2_min <= g_m[1])
//            continue;
//#pragma endregion

#pragma region 弹出的s等于终点的情况，在OPEN中删除被支配的列表
        // 如果弹出的格网是终点
        if (current.getIndex() == endPos)
        {
            // 这个文件用来统计遍历的格网数
            ofstream log_size;
            // 因为终点的H为0，所以终点的F_m就等于G_m
            point_type point_g_m(g_m.begin(), g_m.end());
            // 记录路径
            // 
            COSTS.insert(make_pair(point_g_m, para));
            if (COSTS.size() == 1) {
                log_size.open("D:/桌面/size_and_maxSize.txt", ios::app);
                /*cout << "找到第一条路径之前遍历的格网数" << Util::cnt << endl;
                cout << "找到第一条路径之前遍历的格网数大" << Size << endl;
                log_size << "找到第一条路径之前遍历的格网数" << Util::cnt << endl;
                cout << "找到第一条路径之前最大OPEN队列深度" << maxSize << endl;;
                log_size << "找到第一条路径之前最大OPEN队列深度" << maxSize << endl;*/
                
                /*log_size.close();
                cout << "cnt" << Util::cnt << endl;
                cout << "count" << count << endl;*/

                cout << "wbb遍历格网数Size = " << Size << endl;
                cout << "改进后的cnt = " << Util::cnt << endl;
                cout << "找到第一条路的最大OPEN队列深度" << maxSize << endl;

                log_size << "时间戳： " << Util::PrintCurrentTime();
                log_size << "wbb遍历格网数Size = " << Size << endl;
                log_size << "改进后的cnt = " << Util::cnt << endl;
                log_size << "找到第一条路的最大OPEN队列深度" << maxSize << endl;
                log_size.close();
            }
            cout << "路径数：" << COSTS.size() << endl;
            log_size << "路径数：" << COSTS.size() << endl;
            //cout << COSTS
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
#pragma endregion 

#pragma region 弹出的点不是终点，则需要遍历邻居
        else
            // 参数说明：
            // 1. current 当前弹出的格网
            // 2. 当前的g
            // 3. endpos终点
            // 4. para 弹出格网的para
            NextStep(current, g_m, endPos, Size, para);
        //cout << "current Point " << current.getIndex() << "(" << g_m[0] << "," << g_m[1] << ");;;;;;;;;;;PriorityQueueSize:" << OPEN.size() << ";;;;;;;;result : " << COSTS.size() << endl;
    }
#pragma endregion 
    //cout << "Total Result : " << COSTS.size() << endl;
    return COSTS;
}

Point& Multi_Astar::SelectSPEA2(array_type& g_m, Parameter& parameter) {
    Util::cnt++;
    pareto::front<double, VectorDimension, pair<array_type, Parameter>> nondominantSet;
    pareto::spatial_map<double, VectorDimension, pair<array_type, Parameter> > OPEN_SP;
    for (auto it = OPEN.begin(); it != OPEN.end(); it++)
    {

        OPEN_SP.insert(make_pair(point_type(it->second.first.begin(), it->second.first.end()),
            make_pair(it->first, it->second.second)));
        nondominantSet.insert(make_pair(point_type(it->second.first.begin(), it->second.first.end()),
            make_pair(it->first, it->second.second)));
    }
    vector<pair<double, pair<array_type, pair<array_type, Parameter>>>> nondominantVector;

    for (auto& [key, value] : nondominantSet)
    {
        double cnt = 0;
        //double crowd_distance = nondominantSet.crowding_distance(key);
        //nondominantVector.push_back(make_pair(crowd_distance, make_pair(key.values(), value)));
        for (auto& [k, v] : OPEN_SP)
        {
            if (key.dominates(k))
            {
                cnt++;
            }
        }
        nondominantVector.push_back(make_pair(cnt, make_pair(key.values(), value)));
    }
    sort(nondominantVector.begin(), nondominantVector.end(), compareUseCD2);
    auto target = nondominantVector[0].second;

    // 检查singlesort的排序是不是对
    ofstream log_test;
    log_test.open("D:/桌面/test.txt", ios::app);
    log_test << Util::PrintCurrentTime();

    for (auto it = nondominantVector.begin(); it != nondominantVector.end(); it++)
    {
        log_test << it->first << endl;
    }
    log_test.close();
//#pragma region 可视化
//
//
//    ****************************python可视化****************************
//    vector<double> show_first;
//    vector<double> show_second;
//    // g,F P
//    for (auto it = OPEN.begin(); it != OPEN.end(); it++)
//    {
//        double x = *(it->second.first.begin());
//        show_first.push_back(x);
//        auto temp1 = it->second.first.end() - 1;
//        show_second.push_back(*temp1);
//    }
//
//    double first[] = { *target.first.begin() };
//    auto temp2 = target.first.end() - 1;
//    double second[] = { *temp2 };
//    int N = show_first.size();
//    if (Util::cnt < 20)
//    {
//        string cmd = "plt.scatter(";
//        string s1 = Util::arr_to_string_list(&show_first[0], show_first.size());
//        string s2 = Util::arr_to_string_list(&show_second[0], show_second.size());
//        cmd = cmd + s1 + "," + s2 + ")";
//        PyRun_SimpleString(cmd.c_str());
//        // plt.scatter(x, y, marker = "^")
//        // 
//        // 语句2 
//        string cmd_target = "plt.scatter(" + to_string(first[0]) + "," + to_string(second[0]) + ","
//            + "marker = \"^\"" + ")";
//        PyRun_SimpleString(cmd_target.c_str());
//
//        PyRun_SimpleString("plt.ylim(4000,6000)");
//
//        //plt.savefig("D:/桌面/test.png")
//        string path_ = "D:/open_select_spea2/iter" + to_string(Util::cnt) + ".png";
//        // plt.savefig(" ")
//        string cmd2 = "plt.savefig(\"" + path_ + "\")";
//        PyRun_SimpleString(cmd2.c_str());
//        PyRun_SimpleString("plt.close()");
//
//        //Py_Finalize();
//        //system("pause");
//    }
//#pragma endregion**************************************pyhton 可视化****************************

    // target F, g P
    // upper_bound用于在指定范围内查找大于目标值的第一个元素
    // target F, g P
    auto end = OPEN.upper_bound(target.second.first);
    auto it = OPEN.find(target.second.first);   
    // g,F,Parameter
    while (it != end) {
        if (it->first == target.second.first && it->second.first == target.first && it->second.second == target.second.second)
            break;
        it++;
    }
    // 依据H3index 找到那个点
    // assert使用方法
    // 1. assert(expression) 
    // 如果expression的值为假的  即为0   将终止程序运行
    // 否则 无任何作用
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
    // erase可以传入要删除的迭代器也可以删除要删除键值对的键
    p.op.erase(point_g_m);
    g_m = it->first;
    parameter = it->second.second;
    //在Open列表中删除被取的三元组
    OPEN.erase(it);
    return p;
}
Point& Multi_Astar::SelectSingleSort(array_type& g_m, Parameter& parameter) {
    // 1. 将OPEN放到两个map中，map<<g, <F, P>>, double> 分别按照两个目标进行排序
    // 2. 将这两个OPEN加一起求平均值，放到vector中 排序 
    // F G P

    // 1. 插入到front里
    // 2. 计算cd
    // 3. 排序
    Util::cnt++;
    pareto::front<double, VectorDimension, pair<array_type, Parameter>> nondominantSet;
    // 遍历OPEN列表,并且插入到非支配向量集。
    for (auto it = OPEN.begin(); it != OPEN.end(); it++)
    {
        nondominantSet.insert(make_pair(point_type(it->second.first.begin(), it->second.first.end()), make_pair(it->first, it->second.second)));

    }
    vector<pair<double, pair<array_type, pair<array_type, Parameter>>>> nondominantVector;

    for (auto& [key, value] : nondominantSet)
    {
        double crowd_distance = nondominantSet.crowding_distance(key);
        nondominantVector.push_back(make_pair(crowd_distance, make_pair(key.values(), value)));                                                                                                                                                                                        
    }
    sort(nondominantVector.begin(), nondominantVector.end(), compareUseCD2);
    auto target = nondominantVector[0].second;

    // 检查singlesort的排序是不是对
    ofstream log_test;
    log_test.open("D:/桌面/test.txt", ios::app);
    log_test << Util::PrintCurrentTime();
    
    for (auto it = nondominantVector.begin(); it != nondominantVector.end(); it++)
    {
        log_test << it->first << endl;
    }
    log_test.close();
    //****************************python可视化****************************
    vector<double> show_first;
    vector<double> show_second;
    // g,F P
    for (auto it = OPEN.begin(); it != OPEN.end(); it++)
    {
        double x = *(it->second.first.begin());
        show_first.push_back(x);
        auto temp1 = it->second.first.end() - 1;
        show_second.push_back(*temp1);
    }

    double first[] = { *target.first.begin() };
    auto temp2 = target.first.end() - 1;
    double second[] = { *temp2 };
    int N = show_first.size();
    if (Util::cnt < 20)
    {
        string cmd = "plt.scatter(";
        string s1 = Util::arr_to_string_list(&show_first[0], show_first.size());
        string s2 = Util::arr_to_string_list(&show_second[0], show_second.size());
        cmd = cmd + s1 + "," + s2 + ")";
        PyRun_SimpleString(cmd.c_str());
        // plt.scatter(x, y, marker = "^")
        // 
        // 语句2 
        string cmd_target = "plt.scatter(" + to_string(first[0]) + "," + to_string(second[0]) + ","
            + "marker = \"^\"" + ")";
        PyRun_SimpleString(cmd_target.c_str());

        PyRun_SimpleString("plt.ylim(4000,6000)");

        //plt.savefig("D:/桌面/test.png")
        string path_ = "D:/open_select_single/iter" + to_string(Util::cnt) + ".png";
        // plt.savefig(" ")
        string cmd2 = "plt.savefig(\"" + path_ + "\")";
        PyRun_SimpleString(cmd2.c_str());
        PyRun_SimpleString("plt.close()");

        //Py_Finalize();
        //system("pause");
    }
    //**************************************pyhton 可视化****************************
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
Point& Multi_Astar::SelectCrowdDisOPEN(array_type& g_m, Parameter& parameter) {
    
    // 1. 将multimap转为spatialmap
    // OPEN  <g, <F, P>>
    // OPEN_SP  <F, <g, P>>
    // Front_with_CD <F, < <g, P>, CD > >
    Util::cnt++;
    pareto::spatial_map<double, VectorDimension, pair<array_type, Parameter> > OPEN_SP;
    pareto::front<double, VectorDimension, pair<pair<array_type, Parameter>, double>> Front_with_CD;

    for (auto it = OPEN.begin(); it != OPEN.end(); it++)
    {

        OPEN_SP.insert(make_pair(point_type(it->second.first.begin(), it->second.first.end()),
            make_pair(it->first, it->second.second)));
    }
    // 2. 计算距离  放到front<F, < <g, P>, cd >>
    for (auto it = OPEN_SP.begin(); it != OPEN_SP.end(); it++) {
        
        // 1. 使用 * 解引用，才能获取到it指向的元素
        // 2. it->first.end()指向的不是最后一个元素 
        //  2.1 iter = it->first.end() - 1
        //  2.2 *iter
        double cd = 0;
        double x_it = *(it->first.begin());
        auto temp = it->first.end() - 1;
        double y_it = *temp;

        for (auto iter = OPEN_SP.find_nearest({x_it, y_it}, 3); iter != OPEN_SP.end(); ++iter) {
            //std::cout << iter->first << " -> " << iter->second << std::endl;
            // it-> first 是 begin是x end是y
            double x_iter = *(iter->first.begin());
            auto temp2 = iter->first.end() - 1;
            double y_iter = *temp2;

            cd = cd + fabs(x_it - x_iter) + fabs(y_it - y_iter);
        }
        Front_with_CD.insert(make_pair(point_type(it->first.begin(), it->first.end()), make_pair(make_pair(it->second.first, it->second.second), cd)));
    }
    // 3. 把front中的三元组放到vector中，并进行排序
    // pair<array_type, pair<pair<array_type, Parameter>>>
    vector<pair<array_type, pair<pair<array_type, Parameter>, double>>> nondominantVector;

    for (auto& [key, value] : Front_with_CD) {
        nondominantVector.push_back(make_pair(key.values(), value));
    }
    sort(nondominantVector.begin(), nondominantVector.end(), compareUseCD);
    auto target = make_pair(nondominantVector[0].first, nondominantVector[0].second.first);

    // 检查selectcrowd的排序是不是对
    ofstream log_test;
    log_test.open("D:/桌面/test.txt", ios::app);
    log_test << Util::PrintCurrentTime();
    for (auto it = nondominantVector.begin(); it != nondominantVector.end(); it++)
    {
        log_test << it->second.second << endl;
    }
    
    log_test.close();
    // *******************************python可视化***************************************************
    // 1. 定义数组
    //      (1) OPEN的F, 以及选中的那个点
    // 2. py初始化
    // 3. 调用函数
    // 4. Py Finalize
    vector<double> show_first;
    vector<double> show_second;
    // g,F P
    for (auto it = OPEN.begin(); it != OPEN.end(); it++)
    {
        double x = *(it->second.first.begin());
        show_first.push_back(x);
        auto temp1 = it->second.first.end() - 1;
        show_second.push_back(*temp1);
    }
   
    double first[] = { *target.first.begin() };
    auto temp2 = target.first.end() - 1;
    double second[] = { *temp2 };
    int N = show_first.size();

    //  1. import matplotlib.pyplot as plt
    //  2. plt.scatter(x, y) 
    //  3. plt.scatter(x, y, marker = "^") 用来展示弹出来的点

    
    //Py_Initialize(); /*初始化python解释器,告诉编译器要用的python编译器*/
    //string path = ".";
    // 要的结果是sys.path.append(".")
    // .前后的"需要加转义字符\

    //string chdir_cmd = string("sys.path.append(\"") + "." + "\")";

    //const char* cstr_cmd = chdir_cmd.c_str();
    //PyRun_SimpleString("import sys");
    //PyRun_SimpleString(cstr_cmd);
    //// 语句1 
    //PyRun_SimpleString("import matplotlib.pyplot as plt");
    
    if (Util::cnt < 20)
    {
        string cmd = "plt.scatter(";
        string s1 = Util::arr_to_string_list(&show_first[0], show_first.size());
        string s2 = Util::arr_to_string_list(&show_second[0], show_second.size());
        cmd = cmd + s1 + "," + s2 + ")";
        PyRun_SimpleString(cmd.c_str());
        // plt.scatter(x, y, marker = "^")
        // 
        // 语句2 
        string cmd_target = "plt.scatter(" + to_string(first[0]) + "," + to_string(second[0]) + ","
            + "marker = \"^\"" + ")";
        PyRun_SimpleString(cmd_target.c_str());

        PyRun_SimpleString("plt.ylim(4000,6000)");

        //plt.savefig("D:/桌面/test.png")
        string path_ = "D:/open_select_crowd/iter" + to_string(Util::cnt) + ".png";
        // plt.savefig(" ")
        string cmd2 = "plt.savefig(\"" + path_ + "\")";
        PyRun_SimpleString(cmd2.c_str());
        PyRun_SimpleString("plt.close()");

        //Py_Finalize();
        //system("pause");
    }
    



    // *******************************python可视化***********************************************************

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
// 弹出在非支配向量里面最优的一个
// 先将Open里的三元组全部放入非支配集合
// 再对非支配集合进行排序
// 再选取最优的一个
// 参数：第一个参数为G，第二个参数为从起点到当前的格网个数，主要是为了求出平均通行能力
Point& Multi_Astar::OpenPop(array_type& g_m, Parameter& parameter, int& count) {
    Util::cnt++;
    count++;
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
    #pragma region python可视化
    //****************************python可视化****************************
    //vector<double> show_first;
    //vector<double> show_second;
    //// g,F P
    //for (auto it = OPEN.begin(); it != OPEN.end(); it++)
    //{
    //    double x = *(it->second.first.begin());
    //    show_first.push_back(x);
    //    auto temp1 = it->second.first.end() - 1;
    //    show_second.push_back(*temp1);
    //}

    //double first[] = { *target.first.begin() };
    //auto temp2 = target.first.end() - 1;
    //double second[] = { *temp2 };
    //int N = show_first.size();
    //if (Util::cnt < 20)
    //{
    //    string cmd = "plt.scatter(";
    //    string s1 = Util::arr_to_string_list(&show_first[0], show_first.size());
    //    string s2 = Util::arr_to_string_list(&show_second[0], show_second.size());
    //    cmd = cmd + s1 + "," + s2 + ")";
    //    PyRun_SimpleString(cmd.c_str());
    //    // plt.scatter(x, y, marker = "^")
    //    // 
    //    // 语句2 
    //    string cmd_target = "plt.scatter(" + to_string(first[0]) + "," + to_string(second[0]) + ","
    //        + "marker = \"^\"" + ")";
    //    PyRun_SimpleString(cmd_target.c_str());

    //    PyRun_SimpleString("plt.ylim(4000,6000)");

    //    //plt.savefig("D:/桌面/test.png")
    //    string path_ = "D:/open_select_pop/iter" + to_string(Util::cnt) + ".png";
    //    // plt.savefig(" ")
    //    string cmd2 = "plt.savefig(\"" + path_ + "\")";
    //    PyRun_SimpleString(cmd2.c_str());
    //    PyRun_SimpleString("plt.close()");

    //    //Py_Finalize();
    //    //system("pause");
    //}

    //**************************************pyhton 可视化****************************
#pragma endregion 
    ofstream log_test;
    log_test.open("D:/桌面/test0927.txt", ios::app);
    log_test << Util::PrintCurrentTime();
    log_test.close();
    #pragma region 找到弹出的那个点
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
#pragma endregion 找到弹出的那个点
    #pragma region 从op 移到 cl中
    //将这个点的g_m从Gop移动到Gcl
    point_type point_g_m(it->first.begin(), it->first.end());
    // g_m, 父节点;
    assert(p.op.count(point_g_m) != 0);
    p.cl.insert(make_pair(point_g_m, p.op[point_g_m]));
    p.op.erase(point_g_m);
#pragma endregion 
    // it 是open的迭代器 g, <F, P>
    g_m = it->first;
    parameter = it->second.second;
    //在Open列表中删除被取的三元组
    OPEN.erase(it);
    return p;
}

Point& Multi_Astar::OpenPop2(array_type& g_m, Parameter& parameter, int& count) {
    Util::cnt++;
    count++;
    // f g p
    vector<pair<array_type, pair<array_type, Parameter>>> nondominantVector;
    for (auto& [key, value] : OPEN)
        nondominantVector.push_back(make_pair(value.first, make_pair(key, value.second)));
    sort(nondominantVector.begin(), nondominantVector.end(), CompWithoutDirectionComparison);
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
    // 将这个点的g_m从Gop移动到Gcl
        point_type point_g_m(it->first.begin(), it->first.end());
    // g_m, 父节点;
    assert(p.op.count(point_g_m) != 0);
    p.cl.insert(make_pair(point_g_m, p.op[point_g_m]));
    p.op.erase(point_g_m);
    // it 是open的迭代器 g, <F, P>
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


#pragma region case1: 新点
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
#pragma endregion

#pragma region case2: m已经存在在Op和Cl中
        // 如果不是新点
        // 如果g_m与op与cl中已存在的代价向量相等
        else if (m.op.contains(point_g_m)) {
            // 这句的作用是设置父节点
            // 把current设为这个邻居的父亲
            para.setIndex(current.getIndex());
            m.op.find(point_g_m)->second = para;
        }
        else if (m.cl.contains(point_g_m)) {
            para.setIndex(current.getIndex());
            m.cl.find(point_g_m)->second = para;
        }
#pragma endregion

#pragma region case3：新的g
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
#pragma endregion
}

void Multi_Astar::NextStep2(Point& current, array_type current_g, H3Index goal, int& Size, Parameter& parameter)
{
    // 遍历邻接点
    // it == H3Index;
    //current.g_2_min = current_g[1];

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
        point_type point_g_m(g_m.begin(), g_m.end());

        //  情况1：g_m存在于op和cl中
        if (m.op.contains(point_g_m)) {
            para.setIndex(current.getIndex());
            m.op.find(point_g_m)->second = para;
            continue;
        }
        else if (m.cl.contains(point_g_m)) {
            para.setIndex(current.getIndex());
            m.cl.find(point_g_m)->second = para; para.setIndex(current.getIndex());
            m.cl.find(point_g_m)->second = para;
            continue;
        }
        //  情况2：被g支配
        if ( NodeMap[m.getIndex()].g_2_min <= g_m[1] || m.op.dominates(point_g_m))
            continue;
       
        // 计算h_m:
        array_type h_m = calH(goal, m.getIndex(), current.getIndex());
        // 计算F_m:
        array_type F_m;
        for (int i = 0; i < VectorDimension; i++)
            F_m[i] = g_m[i] + h_m[i];
        // ============================================================================
        point_type point_F_m(F_m.begin(), F_m.end());

        if (NodeMap[goal].g_2_min <= F_m[1])
            continue;

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

        // open 代表自己的index\
        // point中代表父亲的index
        para.setIndex(m.getIndex());
        OPEN.insert(make_pair(g_m, make_pair(F_m, para)));
        // g_m,index为父节点的Parameter
        para.setIndex(current.getIndex());
        m.op.insert(make_pair(point_g_m, para));
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