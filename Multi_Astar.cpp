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
    //comp(a,b)�����ķ���ֵ��һ��boolֵ��������ֵΪtrueʱ���ı�Ԫ��˳�򣬷�֮����Ҫ����Ԫ�ء�
    //���԰����е�a����������ǰһ��λ�õ�Ԫ�أ�b������һ��λ�õ�Ԫ��:
    //���a < b��ʱ��comp(a, b) = true����ôa�ͻᱻ����bǰ�棬���������

    double dou_a = a.second.second;
    double dou_b = b.second.second;

    return dou_a < dou_b;
}
bool compareUseCD2(const pair<double, pair<array_type,pair<array_type, Parameter>>>& a,
    const pair<double, pair<array_type, pair<array_type, Parameter>>>& b) {
    //comp(a,b)�����ķ���ֵ��һ��boolֵ��������ֵΪtrueʱ���ı�Ԫ��˳�򣬷�֮����Ҫ����Ԫ�ء�
    //���԰����е�a����������ǰһ��λ�õ�Ԫ�أ�b������һ��λ�õ�Ԫ��:
    //���a < b��ʱ��comp(a, b) = true����ôa�ͻᱻ����bǰ�棬���������

    double dou_a = a.first;
    double dou_b = b.first;

    return dou_a > dou_b;
}



bool CompareUseDr(const pair<array_type, pair<array_type, Parameter>>& a,
    const pair<array_type, pair<array_type, Parameter>>& b) {
    // �����������ȼ� F_m
    // < ���� > ����
    array_type arr_a = a.second.first;
    array_type arr_b = b.second.first;
    // ����Ϊ��һ��
    if (arr_a[0] == arr_b[0])
        return arr_a[1] < arr_b[1];// by xgg ��һ����ͬ���Ͱ��ڶ�����
    else
        return arr_a[0] < arr_b[0];// by xgg ��һ������ͬ���Ͱ���һ����

}
bool CompWithoutDirectionComparison(const pair<array_type, pair<array_type, Parameter>>& a,
    const pair<array_type, pair<array_type, Parameter>>& b) {
    // �����������ȼ� F_m
    // < ���� > ����
    array_type arr_a = a.first;
    array_type arr_b = b.first;
    // ����Ϊ��һ��
    if (arr_a[0] == arr_b[0])
        return arr_a[1] < arr_b[1];// by xgg ��һ����ͬ���Ͱ��ڶ�����
    else
        return arr_a[0] < arr_b[0];// by xgg ��һ������ͬ���Ͱ���һ����
    //ͨ��������һ��
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
    // maxsize Max_of_the_Number_of_Open_Grid open���е�������
    // size Traverses_the_Number_of_Grid �����ĸ�����


    // 1. ������ʼ֮ǰ��ʼ��startpos��������뵽open�б���
    // 2. ��ʼִ�е�����ÿ�δ�open�е�����һ����Ԫ�飺
    //      2.1 �����������յ㣬������·�ĳɱ����뵽COSTS�У�����OPEN��ɾ��������ɱ�֧�����Ԫ��
    //      2.2 �������㲻���յ� ִ��nextstep
    int count = 0;
    //double g_2_min = DBL_MAX;
    list<H3Index> path;
    front_type pathlist;
    // �ж���ʼ�����ֹ���Ƿ����
    if (NodeMap.count(startPos) == 0 || NodeMap.count(endPos) == 0)
        return pathlist;
    // ����S��
    //array_type g{ 0, 0, 0};
    array_type g{ 0, 0 };
    array_type h = calH(endPos, startPos);
    NodeMap[startPos].op.insert(make_pair(point_type(g.begin(), g.end()), Parameter{ 0 , Comprehensive[startPos] }));
    // ����F_m:
    array_type F;
    for (int i = 0; i < VectorDimension; i++)
        F[i] = g[i] + h[i];
    // G,F,Parameter
    OPEN.insert(make_pair(g, make_pair(F, Parameter{ startPos, Comprehensive[startPos] })));

    while (!OPEN.empty())
    {
        // ��õ������G_M
        //array_type g_m({ 0, 0, 0});
        array_type g_m({ 0, 0 });
        // ��OPEN��ȡ���ŵķ���������
        maxSize = max(maxSize, (int)OPEN.size());
        // ������
        int grid_num = 0;
        // �������F�������
        Parameter para;
        //Point& current = OpenPop(g_m, para);
        Point& current = OpenPop2(g_m, para, count);
        //Point& current = SelectSPEA2(g_m, para);
        //Point& current = SelectCrowdDisOPEN(g_m, para);
        //Point& current = SelectSingleSort(g_m, para);
        //cout << "current Point " << current.getIndex() << "(" << g_m[0] << "," << g_m[1] << ") ;;;;;;;;;;;PriorityQueueSize:" << OPEN.size() << ";;;;;;;;result : " << COSTS.size() << endl;
        //cout << "; PriorityQueueSize:" << OPEN.size() << ";;;;;;;; result: " << COSTS.size() << endl;

        // current g,
        // ��ȡ����� Point& p = NodeMap[endPos]
        if (NodeMap[endPos].g_2_min <= g_m[1])
            continue;

        // ��������ĸ������յ�
        if (current.getIndex() == endPos)
        {
            NodeMap[endPos].g_2_min = g_m[1];
            // ����ļ�����ͳ�Ʊ����ĸ�����
            ofstream log_size;
            // ��Ϊ�յ��HΪ0�������յ��F_m�͵���G_m
            point_type point_g_m(g_m.begin(), g_m.end());
            // ��¼·��
            // 
            COSTS.insert(make_pair(point_g_m, para));
            if (COSTS.size() == 1) {
                log_size.open("D:/����/size_and_maxSize0927.txt", ios::app);
                /*cout << "�ҵ���һ��·��֮ǰ�����ĸ�����" << Util::cnt << endl;
                cout << "�ҵ���һ��·��֮ǰ�����ĸ�������" << Size << endl;
                log_size << "�ҵ���һ��·��֮ǰ�����ĸ�����" << Util::cnt << endl;
                cout << "�ҵ���һ��·��֮ǰ���OPEN�������" << maxSize << endl;;
                log_size << "�ҵ���һ��·��֮ǰ���OPEN�������" << maxSize << endl;*/

                /*log_size.close();
                cout << "cnt" << Util::cnt << endl;
                cout << "count" << count << endl;*/

                cout << "wbb����������Size = " << Size << endl;
                cout << "�Ľ����cnt = " << Util::cnt << endl;
                cout << "�ҵ���һ��·�����OPEN�������" << maxSize << endl;

                log_size << "ʱ����� " << Util::PrintCurrentTime();
                log_size << "wbb����������Size = " << Size << endl;
                log_size << "�Ľ����cnt = " << Util::cnt << endl;
                log_size << "�ҵ���һ��·�����OPEN�������" << maxSize << endl;
                log_size.close();
            }
            cout << "·������" << COSTS.size() << endl;
            continue;
        }


#pragma region �����ĵ㲻���յ㣬����Ҫ�����ھ�
        else
            // ����˵����
            // 1. current ��ǰ�����ĸ���
            // 2. ��ǰ��g
            // 3. endpos�յ�
            // 4. para ����������para
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
    // maxsize Max_of_the_Number_of_Open_Grid open���е�������
    // size Traverses_the_Number_of_Grid �����ĸ�����
    

    // 1. ������ʼ֮ǰ��ʼ��startpos��������뵽open�б���
    // 2. ��ʼִ�е�����ÿ�δ�open�е�����һ����Ԫ�飺
    //      2.1 �����������յ㣬������·�ĳɱ����뵽COSTS�У�����OPEN��ɾ��������ɱ�֧�����Ԫ��
    //      2.2 �������㲻���յ� ִ��nextstep
    int count = 0;
    //double g_2_min = DBL_MAX;
    list<H3Index> path;
    front_type pathlist;
    // �ж���ʼ�����ֹ���Ƿ����
    if (NodeMap.count(startPos) == 0 || NodeMap.count(endPos) == 0)
        return pathlist;
    // ����S��
    //array_type g{ 0, 0, 0};
    array_type g{ 0, 0};
    array_type h = calH(endPos, startPos);
    NodeMap[startPos].op.insert(make_pair(point_type(g.begin(), g.end()), Parameter{ 0 , Comprehensive[startPos]}));
    // ����F_m:
    array_type F;
    for (int i = 0; i < VectorDimension; i++)
        F[i] = g[i] + h[i];
    // G,F,Parameter
    OPEN.insert(make_pair(g, make_pair(F, Parameter{ startPos, Comprehensive[startPos]})));

    while (!OPEN.empty())
    {
        // ��õ������G_M
        //array_type g_m({ 0, 0, 0});
        array_type g_m({ 0, 0 });
        // ��OPEN��ȡ���ŵķ���������
        maxSize = max(maxSize, (int)OPEN.size());
        // ������
        int grid_num = 0;
        // �������F�������
        Parameter para;
        //Point& current = OpenPop(g_m, para);
        Point& current = OpenPop(g_m, para,count);
        //Point& current = SelectSPEA2(g_m, para);
        //Point& current = SelectCrowdDisOPEN(g_m, para);
        //Point& current = SelectSingleSort(g_m, para);
        //cout << "current Point " << current.getIndex() << "(" << g_m[0] << "," << g_m[1] << ") ;;;;;;;;;;;PriorityQueueSize:" << OPEN.size() << ";;;;;;;;result : " << COSTS.size() << endl;
        //cout << "; PriorityQueueSize:" << OPEN.size() << ";;;;;;;; result: " << COSTS.size() << endl;
//#pragma region g_2_min(Sgoal) <= g_current_2 ��s����
//        // current g,
//        if (g_2_min <= g_m[1])
//            continue;
//#pragma endregion

#pragma region ������s�����յ���������OPEN��ɾ����֧����б�
        // ��������ĸ������յ�
        if (current.getIndex() == endPos)
        {
            // ����ļ�����ͳ�Ʊ����ĸ�����
            ofstream log_size;
            // ��Ϊ�յ��HΪ0�������յ��F_m�͵���G_m
            point_type point_g_m(g_m.begin(), g_m.end());
            // ��¼·��
            // 
            COSTS.insert(make_pair(point_g_m, para));
            if (COSTS.size() == 1) {
                log_size.open("D:/����/size_and_maxSize.txt", ios::app);
                /*cout << "�ҵ���һ��·��֮ǰ�����ĸ�����" << Util::cnt << endl;
                cout << "�ҵ���һ��·��֮ǰ�����ĸ�������" << Size << endl;
                log_size << "�ҵ���һ��·��֮ǰ�����ĸ�����" << Util::cnt << endl;
                cout << "�ҵ���һ��·��֮ǰ���OPEN�������" << maxSize << endl;;
                log_size << "�ҵ���һ��·��֮ǰ���OPEN�������" << maxSize << endl;*/
                
                /*log_size.close();
                cout << "cnt" << Util::cnt << endl;
                cout << "count" << count << endl;*/

                cout << "wbb����������Size = " << Size << endl;
                cout << "�Ľ����cnt = " << Util::cnt << endl;
                cout << "�ҵ���һ��·�����OPEN�������" << maxSize << endl;

                log_size << "ʱ����� " << Util::PrintCurrentTime();
                log_size << "wbb����������Size = " << Size << endl;
                log_size << "�Ľ����cnt = " << Util::cnt << endl;
                log_size << "�ҵ���һ��·�����OPEN�������" << maxSize << endl;
                log_size.close();
            }
            cout << "·������" << COSTS.size() << endl;
            log_size << "·������" << COSTS.size() << endl;
            //cout << COSTS
            // �����OPEN�б�֧�������
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

#pragma region �����ĵ㲻���յ㣬����Ҫ�����ھ�
        else
            // ����˵����
            // 1. current ��ǰ�����ĸ���
            // 2. ��ǰ��g
            // 3. endpos�յ�
            // 4. para ����������para
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

    // ���singlesort�������ǲ��Ƕ�
    ofstream log_test;
    log_test.open("D:/����/test.txt", ios::app);
    log_test << Util::PrintCurrentTime();

    for (auto it = nondominantVector.begin(); it != nondominantVector.end(); it++)
    {
        log_test << it->first << endl;
    }
    log_test.close();
//#pragma region ���ӻ�
//
//
//    ****************************python���ӻ�****************************
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
//        // ���2 
//        string cmd_target = "plt.scatter(" + to_string(first[0]) + "," + to_string(second[0]) + ","
//            + "marker = \"^\"" + ")";
//        PyRun_SimpleString(cmd_target.c_str());
//
//        PyRun_SimpleString("plt.ylim(4000,6000)");
//
//        //plt.savefig("D:/����/test.png")
//        string path_ = "D:/open_select_spea2/iter" + to_string(Util::cnt) + ".png";
//        // plt.savefig(" ")
//        string cmd2 = "plt.savefig(\"" + path_ + "\")";
//        PyRun_SimpleString(cmd2.c_str());
//        PyRun_SimpleString("plt.close()");
//
//        //Py_Finalize();
//        //system("pause");
//    }
//#pragma endregion**************************************pyhton ���ӻ�****************************

    // target F, g P
    // upper_bound������ָ����Χ�ڲ��Ҵ���Ŀ��ֵ�ĵ�һ��Ԫ��
    // target F, g P
    auto end = OPEN.upper_bound(target.second.first);
    auto it = OPEN.find(target.second.first);   
    // g,F,Parameter
    while (it != end) {
        if (it->first == target.second.first && it->second.first == target.first && it->second.second == target.second.second)
            break;
        it++;
    }
    // ����H3index �ҵ��Ǹ���
    // assertʹ�÷���
    // 1. assert(expression) 
    // ���expression��ֵΪ�ٵ�  ��Ϊ0   ����ֹ��������
    // ���� ���κ�����
    assert(it != end);
    // ��� index F[0] F[1]
    //cout << "0x" << std::hex << it->second.second << std::dec << "(" << (it->second.first)[0] << "," << (it->second.first)[1] << ")";
    // �ҵ��������
    Point& p = NodeMap[it->second.second.getIndex()];
    //��������g_m��Gop�ƶ���Gcl
    point_type point_g_m(it->first.begin(), it->first.end());
    // g_m, ���ڵ�;
    assert(p.op.count(point_g_m) != 0);
    p.cl.insert(make_pair(point_g_m, p.op[point_g_m]));
    // erase���Դ���Ҫɾ���ĵ�����Ҳ����ɾ��Ҫɾ����ֵ�Եļ�
    p.op.erase(point_g_m);
    g_m = it->first;
    parameter = it->second.second;
    //��Open�б���ɾ����ȡ����Ԫ��
    OPEN.erase(it);
    return p;
}
Point& Multi_Astar::SelectSingleSort(array_type& g_m, Parameter& parameter) {
    // 1. ��OPEN�ŵ�����map�У�map<<g, <F, P>>, double> �ֱ�������Ŀ���������
    // 2. ��������OPEN��һ����ƽ��ֵ���ŵ�vector�� ���� 
    // F G P

    // 1. ���뵽front��
    // 2. ����cd
    // 3. ����
    Util::cnt++;
    pareto::front<double, VectorDimension, pair<array_type, Parameter>> nondominantSet;
    // ����OPEN�б�,���Ҳ��뵽��֧����������
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

    // ���singlesort�������ǲ��Ƕ�
    ofstream log_test;
    log_test.open("D:/����/test.txt", ios::app);
    log_test << Util::PrintCurrentTime();
    
    for (auto it = nondominantVector.begin(); it != nondominantVector.end(); it++)
    {
        log_test << it->first << endl;
    }
    log_test.close();
    //****************************python���ӻ�****************************
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
        // ���2 
        string cmd_target = "plt.scatter(" + to_string(first[0]) + "," + to_string(second[0]) + ","
            + "marker = \"^\"" + ")";
        PyRun_SimpleString(cmd_target.c_str());

        PyRun_SimpleString("plt.ylim(4000,6000)");

        //plt.savefig("D:/����/test.png")
        string path_ = "D:/open_select_single/iter" + to_string(Util::cnt) + ".png";
        // plt.savefig(" ")
        string cmd2 = "plt.savefig(\"" + path_ + "\")";
        PyRun_SimpleString(cmd2.c_str());
        PyRun_SimpleString("plt.close()");

        //Py_Finalize();
        //system("pause");
    }
    //**************************************pyhton ���ӻ�****************************
    auto end = OPEN.upper_bound(target.second.first);
    auto it = OPEN.find(target.second.first);
    // g,F,Parameter
    while (it != end) {
        if (it->first == target.second.first && it->second.first == target.first && it->second.second == target.second.second)
            break;
        it++;
    }
    // ����H3index �ҵ��Ǹ���
    assert(it != end);
    // ��� index F[0] F[1]
    //cout << "0x" << std::hex << it->second.second << std::dec << "(" << (it->second.first)[0] << "," << (it->second.first)[1] << ")";
    // �ҵ��������
    Point& p = NodeMap[it->second.second.getIndex()];
    //��������g_m��Gop�ƶ���Gcl
    point_type point_g_m(it->first.begin(), it->first.end());
    // g_m, ���ڵ�;
    assert(p.op.count(point_g_m) != 0);
    p.cl.insert(make_pair(point_g_m, p.op[point_g_m]));
    p.op.erase(point_g_m);
    g_m = it->first;
    parameter = it->second.second;
    //��Open�б���ɾ����ȡ����Ԫ��
    OPEN.erase(it);
    return p;
}
Point& Multi_Astar::SelectCrowdDisOPEN(array_type& g_m, Parameter& parameter) {
    
    // 1. ��multimapתΪspatialmap
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
    // 2. �������  �ŵ�front<F, < <g, P>, cd >>
    for (auto it = OPEN_SP.begin(); it != OPEN_SP.end(); it++) {
        
        // 1. ʹ�� * �����ã����ܻ�ȡ��itָ���Ԫ��
        // 2. it->first.end()ָ��Ĳ������һ��Ԫ�� 
        //  2.1 iter = it->first.end() - 1
        //  2.2 *iter
        double cd = 0;
        double x_it = *(it->first.begin());
        auto temp = it->first.end() - 1;
        double y_it = *temp;

        for (auto iter = OPEN_SP.find_nearest({x_it, y_it}, 3); iter != OPEN_SP.end(); ++iter) {
            //std::cout << iter->first << " -> " << iter->second << std::endl;
            // it-> first �� begin��x end��y
            double x_iter = *(iter->first.begin());
            auto temp2 = iter->first.end() - 1;
            double y_iter = *temp2;

            cd = cd + fabs(x_it - x_iter) + fabs(y_it - y_iter);
        }
        Front_with_CD.insert(make_pair(point_type(it->first.begin(), it->first.end()), make_pair(make_pair(it->second.first, it->second.second), cd)));
    }
    // 3. ��front�е���Ԫ��ŵ�vector�У�����������
    // pair<array_type, pair<pair<array_type, Parameter>>>
    vector<pair<array_type, pair<pair<array_type, Parameter>, double>>> nondominantVector;

    for (auto& [key, value] : Front_with_CD) {
        nondominantVector.push_back(make_pair(key.values(), value));
    }
    sort(nondominantVector.begin(), nondominantVector.end(), compareUseCD);
    auto target = make_pair(nondominantVector[0].first, nondominantVector[0].second.first);

    // ���selectcrowd�������ǲ��Ƕ�
    ofstream log_test;
    log_test.open("D:/����/test.txt", ios::app);
    log_test << Util::PrintCurrentTime();
    for (auto it = nondominantVector.begin(); it != nondominantVector.end(); it++)
    {
        log_test << it->second.second << endl;
    }
    
    log_test.close();
    // *******************************python���ӻ�***************************************************
    // 1. ��������
    //      (1) OPEN��F, �Լ�ѡ�е��Ǹ���
    // 2. py��ʼ��
    // 3. ���ú���
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
    //  3. plt.scatter(x, y, marker = "^") ����չʾ�������ĵ�

    
    //Py_Initialize(); /*��ʼ��python������,���߱�����Ҫ�õ�python������*/
    //string path = ".";
    // Ҫ�Ľ����sys.path.append(".")
    // .ǰ���"��Ҫ��ת���ַ�\

    //string chdir_cmd = string("sys.path.append(\"") + "." + "\")";

    //const char* cstr_cmd = chdir_cmd.c_str();
    //PyRun_SimpleString("import sys");
    //PyRun_SimpleString(cstr_cmd);
    //// ���1 
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
        // ���2 
        string cmd_target = "plt.scatter(" + to_string(first[0]) + "," + to_string(second[0]) + ","
            + "marker = \"^\"" + ")";
        PyRun_SimpleString(cmd_target.c_str());

        PyRun_SimpleString("plt.ylim(4000,6000)");

        //plt.savefig("D:/����/test.png")
        string path_ = "D:/open_select_crowd/iter" + to_string(Util::cnt) + ".png";
        // plt.savefig(" ")
        string cmd2 = "plt.savefig(\"" + path_ + "\")";
        PyRun_SimpleString(cmd2.c_str());
        PyRun_SimpleString("plt.close()");

        //Py_Finalize();
        //system("pause");
    }
    



    // *******************************python���ӻ�***********************************************************

    auto end = OPEN.upper_bound(target.second.first);
    auto it = OPEN.find(target.second.first);
    // g,F,Parameter
    while (it != end) {
        if (it->first == target.second.first && it->second.first == target.first && it->second.second == target.second.second)
            break;
        it++;
    }
    // ����H3index �ҵ��Ǹ���
    assert(it != end);
    // ��� index F[0] F[1]
    //cout << "0x" << std::hex << it->second.second << std::dec << "(" << (it->second.first)[0] << "," << (it->second.first)[1] << ")";
    // �ҵ��������
    Point& p = NodeMap[it->second.second.getIndex()];
    //��������g_m��Gop�ƶ���Gcl
    point_type point_g_m(it->first.begin(), it->first.end());
    // g_m, ���ڵ�;
    assert(p.op.count(point_g_m) != 0);
    p.cl.insert(make_pair(point_g_m, p.op[point_g_m]));
    p.op.erase(point_g_m);
    g_m = it->first;
    parameter = it->second.second;
    //��Open�б���ɾ����ȡ����Ԫ��
    OPEN.erase(it);
    return p;

}
// �����ڷ�֧�������������ŵ�һ��
// �Ƚ�Open�����Ԫ��ȫ�������֧�伯��
// �ٶԷ�֧�伯�Ͻ�������
// ��ѡȡ���ŵ�һ��
// ��������һ������ΪG���ڶ�������Ϊ����㵽��ǰ�ĸ�����������Ҫ��Ϊ�����ƽ��ͨ������
Point& Multi_Astar::OpenPop(array_type& g_m, Parameter& parameter, int& count) {
    Util::cnt++;
    count++;
    // OPEN g,F,Parameter
    // F,g,Parameter
    pareto::front<double, VectorDimension, pair<array_type, Parameter>> nondominantSet;
    // ����OPEN�б�,���Ҳ��뵽��֧����������
    for (auto it = OPEN.begin(); it != OPEN.end(); it++)
        nondominantSet.insert(make_pair(point_type(it->second.first.begin(), it->second.first.end()), make_pair(it->first, it->second.second)));
    // ���ղ�ͬ���������ȼ�����
    // F,g,index
    vector<pair<array_type, pair<array_type, Parameter>>> nondominantVector;
    for (auto& [key, value] : nondominantSet)
        nondominantVector.push_back(make_pair(key.values(), value));
    // ��F��������������
    sort(nondominantVector.begin(), nondominantVector.end(), CompWithoutDirectionComparison);
    // ȡ���ţ���С���Ĵ���������Ҳ���ǵ�һ���������ǰ�g����λ
    auto target = nondominantVector[0];
    #pragma region python���ӻ�
    //****************************python���ӻ�****************************
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
    //    // ���2 
    //    string cmd_target = "plt.scatter(" + to_string(first[0]) + "," + to_string(second[0]) + ","
    //        + "marker = \"^\"" + ")";
    //    PyRun_SimpleString(cmd_target.c_str());

    //    PyRun_SimpleString("plt.ylim(4000,6000)");

    //    //plt.savefig("D:/����/test.png")
    //    string path_ = "D:/open_select_pop/iter" + to_string(Util::cnt) + ".png";
    //    // plt.savefig(" ")
    //    string cmd2 = "plt.savefig(\"" + path_ + "\")";
    //    PyRun_SimpleString(cmd2.c_str());
    //    PyRun_SimpleString("plt.close()");

    //    //Py_Finalize();
    //    //system("pause");
    //}

    //**************************************pyhton ���ӻ�****************************
#pragma endregion 
    ofstream log_test;
    log_test.open("D:/����/test0927.txt", ios::app);
    log_test << Util::PrintCurrentTime();
    log_test.close();
    #pragma region �ҵ��������Ǹ���
    auto end = OPEN.upper_bound(target.second.first);
    auto it = OPEN.find(target.second.first);
    // g,F,Parameter
    while (it != end) {
        if (it->first == target.second.first && it->second.first == target.first && it->second.second == target.second.second)
            break;
        it++;
    }
    // ����H3index �ҵ��Ǹ���
    assert(it != end);
    // ��� index F[0] F[1]
    //cout << "0x" << std::hex << it->second.second << std::dec << "(" << (it->second.first)[0] << "," << (it->second.first)[1] << ")";
    // �ҵ��������
    Point& p = NodeMap[it->second.second.getIndex()];
#pragma endregion �ҵ��������Ǹ���
    #pragma region ��op �Ƶ� cl��
    //��������g_m��Gop�ƶ���Gcl
    point_type point_g_m(it->first.begin(), it->first.end());
    // g_m, ���ڵ�;
    assert(p.op.count(point_g_m) != 0);
    p.cl.insert(make_pair(point_g_m, p.op[point_g_m]));
    p.op.erase(point_g_m);
#pragma endregion 
    // it ��open�ĵ����� g, <F, P>
    g_m = it->first;
    parameter = it->second.second;
    //��Open�б���ɾ����ȡ����Ԫ��
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
    // ����H3index �ҵ��Ǹ���
    assert(it != end);
    // ��� index F[0] F[1]
    //cout << "0x" << std::hex << it->second.second << std::dec << "(" << (it->second.first)[0] << "," << (it->second.first)[1] << ")";
    // �ҵ��������
    Point& p = NodeMap[it->second.second.getIndex()];
    // ��������g_m��Gop�ƶ���Gcl
        point_type point_g_m(it->first.begin(), it->first.end());
    // g_m, ���ڵ�;
    assert(p.op.count(point_g_m) != 0);
    p.cl.insert(make_pair(point_g_m, p.op[point_g_m]));
    p.op.erase(point_g_m);
    // it ��open�ĵ����� g, <F, P>
    g_m = it->first;
    parameter = it->second.second;
    //��Open�б���ɾ����ȡ����Ԫ��
    OPEN.erase(it);
    return p;

}

void Multi_Astar::NextStep(Point& current, array_type current_g, H3Index goal, int& Size, Parameter& parameter)
{

    // �����ڽӵ�
    // it == H3Index;
    for (auto it = current.neighbours.begin(); it != current.neighbours.end(); it++) {
        // ����ڽӵ�
        Point& m = NodeMap[*it];
        if (current.getIndex() == *it)
            continue;
        // ============================================================================
        // ����g_m:
        bool flag = true;
        auto para = parameter;
        array_type g_m = calG(current_g, current.getIndex(), *it, flag, para);
        // ��Ҫ�Ų��ڽӵ��Ƿ��ͨ��
        if (!flag)
            continue;
        Size++;
        // ����h_m:
        array_type h_m = calH(goal, m.getIndex(), current.getIndex());
        // ����F_m:
        array_type F_m;
        for (int i = 0; i < VectorDimension; i++)
            F_m[i] = g_m[i] + h_m[i];
        // ============================================================================
        point_type point_F_m(F_m.begin(), F_m.end());
        point_type point_g_m(g_m.begin(), g_m.end());


#pragma region case1: �µ�
        // ���m���µ�
        if (m.op.empty() && m.cl.empty()) {
            // ���F_m��COSTS֧�䣬��ô������
            if (COSTS.dominates(point_F_m))
                continue;
            // ����OPEN�б�
            para.setIndex(m.getIndex());
            OPEN.insert(make_pair(g_m, make_pair(point_F_m.values(), para)));
            // ���Gop_m,��Ȼ�����ʱ����ж��Ƿ����ף��������ʱ��opΪ�գ���Ӱ��
            para.setIndex(current.getIndex());
            m.op.insert(make_pair(point_g_m, para));
        }
#pragma endregion

#pragma region case2: m�Ѿ�������Op��Cl��
        // ��������µ�
        // ���g_m��op��cl���Ѵ��ڵĴ����������
        else if (m.op.contains(point_g_m)) {
            // �������������ø��ڵ�
            // ��current��Ϊ����ھӵĸ���
            para.setIndex(current.getIndex());
            m.op.find(point_g_m)->second = para;
        }
        else if (m.cl.contains(point_g_m)) {
            para.setIndex(current.getIndex());
            m.cl.find(point_g_m)->second = para;
        }
#pragma endregion

#pragma region case3���µ�g
        // ��������µ�Ҳ����ȣ���ô��Ҫ����
        // g_m����op��cl����������֧��
        else if (!m.op.dominates(point_g_m) && !m.cl.dominates(point_g_m)) {
            //��op��ɾ����g_m֧���g'_m����һͬɾ��OPEN�е�
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
            //��m��cl��ɾ����g_m֧���g'_m
            temp.clear();
            for (auto& [p, v] : m.cl) {
                if (point_g_m.dominates(p))
                    continue;
                temp.insert(make_pair(p, v));
            }
            m.cl = temp;
            //���F_m�Ƿ�֧��
            if (COSTS.dominates(point_F_m))
                continue;
            else {
                // �����֧�䣬����ǰ�����OPEN�б��op�б�
                // g_m,F_m,indexΪ��ǰ�ڵ��Parameter
                para.setIndex(m.getIndex());
                OPEN.insert(make_pair(g_m, make_pair(F_m, para)));
                // g_m,indexΪ���ڵ��Parameter
                para.setIndex(current.getIndex());
                m.op.insert(make_pair(point_g_m, para));
            }
        }

    }
#pragma endregion
}

void Multi_Astar::NextStep2(Point& current, array_type current_g, H3Index goal, int& Size, Parameter& parameter)
{
    // �����ڽӵ�
    // it == H3Index;
    //current.g_2_min = current_g[1];

    for (auto it = current.neighbours.begin(); it != current.neighbours.end(); it++) {
        
        // ����ڽӵ�
        Point& m = NodeMap[*it];
        if (current.getIndex() == *it)
            continue;
        // ============================================================================
        // ����g_m:
        bool flag = true;
        auto para = parameter;
        array_type g_m = calG(current_g, current.getIndex(), *it, flag, para);
        // ��Ҫ�Ų��ڽӵ��Ƿ��ͨ��
        if (!flag)
            continue;
        Size++;
        point_type point_g_m(g_m.begin(), g_m.end());

        //  ���1��g_m������op��cl��
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
        //  ���2����g֧��
        if ( NodeMap[m.getIndex()].g_2_min <= g_m[1] || m.op.dominates(point_g_m))
            continue;
       
        // ����h_m:
        array_type h_m = calH(goal, m.getIndex(), current.getIndex());
        // ����F_m:
        array_type F_m;
        for (int i = 0; i < VectorDimension; i++)
            F_m[i] = g_m[i] + h_m[i];
        // ============================================================================
        point_type point_F_m(F_m.begin(), F_m.end());

        if (NodeMap[goal].g_2_min <= F_m[1])
            continue;

        //��op��ɾ����g_m֧���g'_m����һͬɾ��OPEN�е�
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

        // open �����Լ���index\
        // point�д����׵�index
        para.setIndex(m.getIndex());
        OPEN.insert(make_pair(g_m, make_pair(F_m, para)));
        // g_m,indexΪ���ڵ��Parameter
        para.setIndex(current.getIndex());
        m.op.insert(make_pair(point_g_m, para));
    }

}
array_type Multi_Astar::calH(H3Index goal, H3Index current) {
    if (current == goal)
        return array_type({ 0,0 });// by Xgg ǿ��ת����
    else {
        array_type res;
        // ��һ������������
        GeoPoint g1, g2;
        cellToPoint(current, &g1);
        cellToPoint(goal, &g2);
        double distance = Util::round(Util::calcdistance(g1, g2), 3);
        double h = fabs(DEM[current] - DEM[goal]);
        // ���㵱ǰ�㵽�յ��������
        double edge = getHexagonEdgeLengthAvgM(getResolution(goal));
        int grid = distance / (1.732050 * edge) + 1;
        // �ڶ����������¶���ͨ������
        double s = h / distance;
        double grad = atan(s) * 180.0 / 3.1415926;
        // ����ʹ��1-grad���¶�ԽС������ֵԽС
        grad = -0.00004061 * pow(grad, 3) + 0.002644 * pow(grad, 2) - 0.07607 * grad + 0.9977;
        if (grad > 1)
            grad = 1;
        if (grad < 0)
            grad = 0;
        // �������⣬���ܻ���ַ�����ͨ������ǿ�ĸ�������������ü��Ϸ�����
        // res[0] = distance / 100;
        // res[1] = Util::round((2 - grad - Comprehensive[current]) * distance, 2);
        // res[1] = 0;
        // res[1] = distance;
        // res[1] = int((1 - Comprehensive[current]) * grid);
        // res[1] = grid * 100 - Util::round(grad, 2);
        res[0] = Util::round(distance / 10, 1);
        res[1] = Util::round(distance / 10, 1);
        // ����ֵ
        return res;
    }
}

// ����ʹ��
array_type Multi_Astar::calH(H3Index goal, H3Index current, H3Index parent) {
    if (current == goal)
        return array_type({ 0,0});
    else {
        array_type res;
        // ��һ������������
        GeoPoint g1, g2, g3;
        cellToPoint(current, &g1);
        cellToPoint(goal, &g2);
        cellToPoint(parent, &g3);
        double distance = Util::round(Util::calcdistance(g1, g2), 3);
        double h = fabs(DEM[current] - DEM[goal]);
        // ���㵱ǰ�㵽�յ��������
        double edge = getHexagonEdgeLengthAvgM(getResolution(goal));
        int grid = distance / (1.732050 * edge) + 1;
        // �ڶ����������¶���ͨ������
        double s = h / distance;
        double grad = atan(s) * 180.0 / 3.1415926;
        // ����ʹ��1-grad���¶�ԽС������ֵԽС
        grad = -0.00004061 * pow(grad, 3) + 0.002644 * pow(grad, 2) - 0.07607 * grad + 0.9977;
        if (grad > 1)
            grad = 1;
        if (grad < 0)
            grad = 0;
        // �������⣬���ܻ���ַ�����ͨ������ǿ�ĸ�������������ü��Ϸ�����
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
        // ����ֵ
        return  res;
    }
}

// flagΪ�ж��Ƿ��ͨ�У�parameterΪЯ���Ĳ���
array_type Multi_Astar::calG(array_type& parent_g, H3Index parent, H3Index current, bool& flag, Parameter& parameter) {
    array_type res = { 0,0 };
    if (current == parent)
        return res;
    else {
        // ��һ������������
        GeoPoint g1, g2;
        cellToPoint(current, &g1);
        cellToPoint(parent, &g2);
        double distance = Util::round(Util::calcdistance(g1, g2), 3);
        double h = fabs(DEM[current] - DEM[parent]);
        // ����ʵ�ʾ���
        // int real_distance = Util::round(sqrt(pow(h, 2) + pow(distance, 2)), 0);
        // �ڶ����������¶�
        double s = h / distance;
        double grad = atan(s) * 180.0 / 3.1415926;
        // ����ʹ��1-grad��ԽСԽ��
        grad = -0.00004061 * pow(grad, 3) + 0.002644 * pow(grad, 2) - 0.07607 * grad + 0.9977;
        // �ж��Ƿ��ͨ��
        if (Comprehensive[current] <= 0) {
            flag = false;
            return res;
        }
        // �ж��¶��Ƿ��ͨ��
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
        // ����ֵ
        return res;
    }
}

vector<pair<list<H3Index>, array_type>> Multi_Astar::printPathByNotRecursion(H3Index end, H3Index start) {

    ////��ջ�����ڼ��㴦��·��
    //stack<H3Index> majorStack;
    ////��ջ�����ڴ洢������ڵ�
    //stack<H3Index> minorStack;
    ////���ڴ洢һ���ڵ���ӽڵ�ĵ�������
    //map<H3Index, int> popCount;
    ////���ڵ��븱ջ
    //minorStack.push(start);
    //// ������
    //map<list<H3Index>, array_type> pathList;
    //list<H3Index> currentList;
    //array_type cost;
    //cost.fill(0);
    //while (!minorStack.empty()) {
    //    // ����ջ������ջ
    //    // Ҫ����Ľڵ�
    //    H3Index minLast = minorStack.top();
    //    minorStack.pop();
    //    majorStack.push(minLast);
    //    // ���·����������
    //    if (!currentList.empty()) {
    //        Point& minLastPoint = NodeMap[minLast];
    //        for (int i = 0; i < Multi_Astar::length; i++) {
    //            cost[i] += minLastPoint.children[currentList.back()][i];
    //        }
    //    }
    //    // ��¼��ǰ�ڵ�
    //    currentList.push_back(minLast);
    //    if (minLast == stringToH3("8931852C5A7FFFF")) {
    //        cout << "  " << endl;
    //    }
    //    // ���ýڵ���ӽڵ��븱ջ
    //    Point& currentPoint = NodeMap[minLast];
    //    if (currentPoint.cl.size() != 0 && minLast != end) {
    //        for (auto it = currentPoint.cl.begin(); it != currentPoint.cl.end(); it++) {
    //            minorStack.push(it->second);
    //            // ���������ӽ�ջ������ձ��
    //            // popCount.erase(it->second);
    //        }
    //    }
    //    //ջ��Ԫ��
    //    H3Index majLast = majorStack.top();
    //    //ѭ��������ջ��ΪҶ�ӽڵ� �� ջ���ڵ㺢�ӽڵ��������
    //    while (majLast == end || NodeMap[majLast].cl.size() == 0 ||
    //        (popCount.count(majLast) != 0 && popCount[majLast] == NodeMap[majLast].cl.size())) {
    //        H3Index last = majorStack.top();
    //        majorStack.pop();
    //        //��ʱ��ջΪ�գ��������
    //        if (majorStack.empty()) {
    //            return pathList;
    //        }
    //        majLast = majorStack.top();
    //        //��һ�ε������ӽڵ㣬����������Ϊ1
    //        if (popCount.count(majLast) == 0) {
    //            popCount.insert(make_pair(majLast, 1));
    //        }
    //        //�ǵ�һ�ε������ӽڵ㣬��ԭ�л����ϼ�1
    //        else {
    //            popCount[majLast] += 1;
    //        }
    //        //�����Ҷ�ӽڵ�Ž��������·������
    //        if (last == end) {
    //            pathList.insert(make_pair(currentList, cost));
    //        }
    //        // ��ȥ����ֵ
    //        Point& lastPoint = NodeMap[last];
    //        for (int i = 0; i < Multi_Astar::length; i++) {
    //            cost[i] -= lastPoint.children[majLast][i];
    //        }
    //        // ������
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