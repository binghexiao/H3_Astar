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

front_type Multi_Astar::search(H3Index startPos, H3Index endPos, int& maxSize, int& Size)
{
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
        Point& current = OpenPop(g_m, para);
        //cout << "current Point " << current.getIndex() << "(" << g_m[0] << "," << g_m[1] << ") ;;;;;;;;;;;PriorityQueueSize:" << OPEN.size() << ";;;;;;;;result : " << COSTS.size() << endl;
        //cout << "; PriorityQueueSize:" << OPEN.size() << ";;;;;;;; result: " << COSTS.size() << endl;
        if (current.getIndex() == endPos)
        {
            // ��Ϊ�յ��HΪ0�������յ��F_m�͵���G_m
            point_type point_g_m(g_m.begin(), g_m.end());
            // ��¼·��
            // 
            COSTS.insert(make_pair(point_g_m, para));
            cout << "·������" << COSTS.size() << endl;
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
        else
            NextStep(current, g_m, endPos, Size, para);
        //cout << "current Point " << current.getIndex() << "(" << g_m[0] << "," << g_m[1] << ");;;;;;;;;;;PriorityQueueSize:" << OPEN.size() << ";;;;;;;;result : " << COSTS.size() << endl;
    }
    //cout << "Total Result : " << COSTS.size() << endl;
    return COSTS;
}

// �����ڷ�֧�������������ŵ�һ��
// �Ƚ�Open�����Ԫ��ȫ�������֧�伯��
// �ٶԷ�֧�伯�Ͻ�������
// ��ѡȡ���ŵ�һ��
// ��������һ������ΪG���ڶ�������Ϊ����㵽��ǰ�ĸ�����������Ҫ��Ϊ�����ƽ��ͨ������
Point& Multi_Astar::OpenPop(array_type& g_m, Parameter& parameter) {
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
        // ��������µ�
        // ���g_m��op��cl���Ѵ��ڵĴ����������
        else if (m.op.contains(point_g_m)) {
            para.setIndex(current.getIndex());
            m.op.find(point_g_m)->second = para;
        }
        else if (m.cl.contains(point_g_m)) {
            para.setIndex(current.getIndex());
            m.cl.find(point_g_m)->second = para;
        }
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