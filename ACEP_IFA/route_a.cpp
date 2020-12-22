#include"route_a.h"

//密度计算 num表示需要增加的人群数量  num为-1时 表示有人走出这条路径
void density(myType a, myType b, int num = 1) {
	//vector<int>::iterator it;
	//it = find(sites.begin(), sites.end(), a);
	//if (it == sites.end()) {
	//	cerr << "density 更新失败" << endl;
	//	return;
	//}
	int pos = get_pos(a);
	for (auto& s : density_mtx[pos]) {
		if (s.first == b) {
			s.second.second += num;
			return;
		}
	}
	cerr << "density error" << endl;
}

bool isDensityOk(myType a, myType b) {
	int pos = get_pos(a);
	for (auto& s : density_mtx[pos]) {
		if (s.first == b) {
			if (s.second.second >= s.second.first)
				return false;
			else return true;
		}
	}
	cerr << "isDensityOk error:" << a << "," << b << endl;
	return false;
}

void updateHeur(myType a, myType b) {
	int pos = get_pos(a);
	for (int i = 0; i < heur_mtx.at(pos).size(); ++i) {
		if (heur_mtx.at(pos).at(i).first == b) {
			if (density_mtx.at(pos).at(i).second.second == 0)
				heur_mtx.at(pos).at(i).second = 1 / g.getweight(a, b);
			else
				heur_mtx.at(pos).at(i).second = 1 / (g.getweight(a, b) * density_mtx.at(pos).at(i).second.second);
			return;
		}
	}
	cerr << "updateHeur error : " << a << b << endl;
}

void updateTrans(bool tag) {
	if (tag == 0) {//tag==0 表示为初代
		vector<pair<myType, float> > tmp;
		pair<myType, float> p;
		for (auto s : pher_mtx) {
			//统计每一行信息素总和 较繁琐 待简化
			float sum = 0;
			for (auto j : s) {
				sum += j.second;
			}
			for (auto j : s) {
				p.first = j.first;
				p.second = j.second / sum;
				tmp.push_back(p);
			}
			trans_mtx.push_back(tmp);
			tmp.clear();
		}
	}
	else {//不是第一代 无需pushback操作
		//int cnt = 0;
		//for (auto &s : pher_mtx) {
		//	//统计每一行信息素总和 较繁琐 待简化
		//	float sum = 0;
		//	for (auto j : s) {
		//		sum += j.second;
		//	}
		//	for (int i = 0; i < s.size(); ++i) {
		//		trans_mtx[cnt][i].second = s[i].second / sum;
		//	}
		//	cnt++;
		//}
		for (auto s : sites) {
			updateNodeTrans(s);
		}
	}
}
void updateNodeTrans(myType a) {
	int pos = get_pos(a);
	float sum = 0;
	int tag;
	for (int i = 0; i < trans_mtx.at(pos).size(); ++i) {
		sum += pher_mtx.at(pos).at(i).second * pow(heur_mtx.at(pos).at(i).second, beta);
	}
	for (int j = 0; j < trans_mtx.at(pos).size(); ++j) {
		trans_mtx.at(pos).at(j).second =
			pher_mtx.at(pos).at(j).second * pow(heur_mtx.at(pos).at(j).second, beta) / sum;
	}
	return;
}
//初始化
void init() {
	int k, site;
	cout << "Please input safe sites'num" << endl;
	cin >> k;
	while (k--) {
		cin >> site;
		safe_site.push_back(site);
	}
	//获取全部节点 
	sites = g.getAllvertex();

	history_best.first = -1;
	//初始化信息素矩阵
	vector<myType> v;
	pair<myType, float> p;
	vector<pair<myType, float> > tmp;
	vector<pair<myType, pair<int, int> > > den_tmp;
	pair<myType, pair<int, int> > den_p_tmp;
	pair<int, int> den_p_p_tmp;

	for (auto s : sites) {
		v = g.Adj(s);
		for (auto j : v) {
			p.first = j;
			p.second = pher_start;//信息素初始值 pher_start
			tmp.push_back(p);
			//密度初始化

			den_p_p_tmp.first = cap;
			den_p_p_tmp.second = 0;
			den_p_tmp.first = j;//外层
			den_p_tmp.second = den_p_p_tmp;
			den_tmp.push_back(den_p_tmp);
		}
		pher_mtx.push_back(tmp);
		density_mtx.push_back(den_tmp);
		den_tmp.clear();
		tmp.clear();
		v.clear();
	}
	//flow_mtx = pher_mtx;//初始 随机更新 flow
	updateTrans(0);
	//初始化密度矩阵

	return;

}
void new_init() {
	int site;
	//获取全部节点 
	sites = g.getAllvertex();
	float choice;
	//目的节点集合 safe_site
	for (int i = 0; i < m; ++i) {
		do {
			choice = rand_n(sites.size());
		} while (sites[choice] == 48);// 521
		safe_site.push_back(sites[choice]);
	}
	cout << "safe_sites finish with size:" << safe_site.size() << endl;
	// 随机选定源节点

	//for (int i = 0; i < m; ++i) {
	//	do {
	//		choice = rand_n(sites.size());
	//	} while (sites[choice] == safe_site[i]);// 521
	//	//safe_site.push_back(sites[choice]);
	//	source[i] = sites[choice];
	//}
	history_best.first = -1;
	//初始化信息素矩阵
	vector<myType> v;
	pair<myType, float> p;
	vector<pair<myType, float> > tmp;
	vector<pair<myType, float> > tmp_heur;
	vector<pair<myType, pair<int, int> > > den_tmp;
	pair<myType, pair<int, int> > den_p_tmp;
	pair<int, int> den_p_p_tmp;
	float weight_road;
	for (auto s : sites) {
		v = g.Adj(s);
		for (auto j : v) {
			p.first = j;
			p.second = pher_start;//信息素初始值 pher_start
			tmp.push_back(p);
			//启发式初始化
			weight_road = g.getweight(s, j);
			p.second = 1 / weight_road;
			tmp_heur.push_back(p);
			//密度初始化
			den_p_p_tmp.first = cap;
			den_p_p_tmp.second = 0;
			den_p_tmp.first = j;//外层
			den_p_tmp.second = den_p_p_tmp;
			den_tmp.push_back(den_p_tmp);
		}
		heur_mtx.push_back(tmp_heur);
		pher_mtx.push_back(tmp);
		density_mtx.push_back(den_tmp);
		tmp_heur.clear();
		den_tmp.clear();
		tmp.clear();
		v.clear();
	}
	//flow_mtx = pher_mtx;//初始 随机更新 flow
	updateTrans(0);
	//初始化密度矩阵

	return;
}

//构建路径*************->
void routeConstricution() {
	if (G == 0) return;
	traffic_time = 1;
	vector<bool>::iterator it = find(ant_dest.begin(), ant_dest.end(), false);
	if (it == ant_dest.end()) return;//是否全到达终点
	ant_solution.clear();
	ant_solution.resize(m);
	myType current_pos;
	myType before_pos;
	//vector<int> v;//存储每只蚂蚁路径
	//for (int i = 0; i < m; ++i) {//i表示蚂蚁编号
	//	if (ant_dest[i] != true) {
	//		int current_pos = source[i];//当前节点
	//		v.push_back(current_pos);
	//		while (!isSafe(current_pos)) {
	//			float tmp = rand_p();//轮盘赌
	//			int choice = 0;//要选择的节点
	//			float sump = 0.0;
	//			while (sump <= tmp) {
	//				sump += trans_mtx[current_pos][choice++].second;
	//			}//生成0时候默认选择第一节点
	//			current_pos = trans_mtx[current_pos][--choice].first;//选择节点
	//			v.push_back(current_pos);
	//		}
	//		ant_solution.push_back(v);//存储每只蚂蚁的路径
	//		ant_dest[i] = true;
	//		v.clear();
	//	}
	//	else {
	//		continue;
	//	}
	//}
	fit_mtx.clear();
	fit_mtx.resize(m);
	for (int j = 0; j < m; ++j) {
		current_pos = source[j];//当前节点
		ant_solution[j].push_back(current_pos);
	}
	for (; find(ant_dest.begin(), ant_dest.end(), false) != ant_dest.end(); ++traffic_time) {//是否全部到达终点
		for (int i = 0; i < m; ++i) {
			if (ant_dest[i] != true) {
				current_pos = ant_solution[i].back();
				before_pos = current_pos;
				//轮盘赌选择节点
				int choice_test = 0;
				do {
					if (choice_test++ > 50) {
						cerr << choice_test << " ::choice_test erro: " << before_pos << "****" << current_pos << endl;
						current_pos = before_pos;
						break;
					}
					current_pos = before_pos;
					float tmp = rand_p();//轮盘赌
					int choice = 0;//要选择的节点
					float sump = 0.0;
					if (tmp == 0) {
						current_pos = trans_mtx[get_pos(current_pos)][0].first;
					}
					else {
						while (sump < tmp) {
							sump += trans_mtx[get_pos(current_pos)][choice++].second;
						}//生成0时候默认选择第一节点
						current_pos = trans_mtx[get_pos(current_pos)][--choice].first;//选择节点
					}
				} while (current_pos == before_pos || !isDensityOk(before_pos, current_pos) || current_pos == target);

				ant_solution[i].push_back(current_pos);

				if (isEnd(current_pos, i)) {
					ant_dest[i] = true;
				}
				else {
					density(before_pos, current_pos);
				}//走入新边  最后一段还未处理 安全路径之前  安全之路不计算拥挤度

				if (traffic_time != 1) {
					density(ant_solution[i].at(ant_solution[i].size() - 3), before_pos, -1);//上一条边 走出
					updateHeur(ant_solution[i].at(ant_solution[i].size() - 3), before_pos);// 更新启发式  上一次经过的路径
					updateNodeTrans(ant_solution[i].at(ant_solution[i].size() - 3));// 更新上上 节点的 转移概率
				}

				float vel = velocityCompute(before_pos, current_pos);
				if (before_pos == current_pos) fit_mtx[i] += 1;
				else fit_mtx[i] += g.getweight(before_pos, current_pos) / vel;
				updateLocalPher(before_pos, current_pos);
				updateHeur(before_pos, current_pos);
				updateNodeTrans(before_pos);
				//判断是否到达传送点
				int cur_tmp_safe = get_pos(safe_site[i]);
				if (current_pos == transmit_p[cur_tmp_safe][0]) {
					auto k = find(short_p[cur_tmp_safe].begin(), short_p[cur_tmp_safe].end(), current_pos);
					int cur_tmp_pos = *(k++ + 1);//迭代器为end的时候不能增加
					while (cur_tmp_pos != transmit_p[cur_tmp_safe][1]) {
						ant_solution[i].push_back(cur_tmp_pos);
						cur_tmp_pos = *(++k);
					}
					ant_solution[i].push_back(transmit_p[cur_tmp_safe][1]);
					//判断是否到达safe place
					if (isEnd(transmit_p[cur_tmp_safe][1], i)) {
						ant_dest[i] = true;
					}
				}
			}//end if ant_dest true
			else continue;
		}// end i
		//输出本代选择的节点
		//for (int i = 0; i < m; ++i) {
		//	cout << ant_solution[i].back() << endl;
		//}
	}

	ant_dest.clear();
	ant_dest.resize(m);
	G--;//代数增加
	return;
}

//适应度计算  只计算路径长度
vector<float> fitCompute() {
	vector<float> fittness;
	float sum;
	for (int i = 0; i < m; ++i) {
		float sum = 0.0;
		for (vector<int>::iterator it = ant_solution[i].begin(); it != ant_solution[i].end() - 1; ++it) {
			sum += g.getweight(*it, *(it + 1));
		}
		fittness.push_back(sum);
	}
	return fittness;
}

//更新till now best solution
void updateHistoryBeSou(float fittness, int serial) {
	if (history_best.first != -1) {
		if (history_best.first > fittness) {
			history_best.first = fittness;
			history_best.second = ant_solution[serial];
		}
	}
	else {
		history_best.first = fittness;
		history_best.second = ant_solution[serial];
	}
}

//计算每代flow  flow_mtx
/******** flow 与 density 名称exchange **********/
void flowCompute() {
	// 清零 flow_mtx
	flow_mtx.clear();
	flow_mtx.resize(sites.size());
	unordered_map<myType, float> tmp;
	for (auto s : ant_solution) {
		for (int j = 0; j < s.size() - 1; ++j) {
			if (flow_mtx[get_pos(s[j])].count(s[j + 1]) == 0)
				flow_mtx[get_pos(s[j])][s[j + 1]] = 1;
			else
				flow_mtx[get_pos(s[j])][s[j + 1]] += 1;
		}
	}


}

void updateGlobalPher() {
	vector<float> v{ fit_mtx };
	sort(v.begin(), v.end());
	//更新flow 每代只更新一次
	cur_global_route.clear();
	flowCompute();
	int dis;
	vector<float>::iterator it;
	vector<myType>::iterator it_t;
	const float weight_pher = 5 * pow(10, 3);
	for (int i = 0; i < w; ++i) {// 找到fittness 前几的解
		it = find(fit_mtx.begin(), fit_mtx.end(), v[i]);
		dis = distance(fit_mtx.begin(), it);// 排名第i的解 在fit_mtx中的编号
		if (i == 0) {
			updateHistoryBeSou(v[i], dis);//更新 迄今为止 最好的 solution,只需与本代最优解比较
		}
		for (int j = 0; j < ant_solution[dis].size() - 1; ++j) {
			if (cur_global_route.count({ ant_solution[dis].at(j),ant_solution[dis].at(j + 1) }))
				continue;
			it_t = find(sites.begin(), sites.end(), ant_solution[dis].at(j));
			int dist = distance(sites.begin(), it_t);//找到节点编号
			for (auto& s : pher_mtx[dist]) {
				if (s.first == ant_solution[dis].at(j + 1)) {// 要更新路径的 目的地节点
					cur_global_route.insert({ ant_solution[dis].at(j),ant_solution[dis].at(j + 1) });
					s.second = s.second * (1 - rho) + rho * ((w - 1 - i) * flow_mtx[dist].at(s.first) / (m * fit_mtx[dis] * weight_pher));
					break;
				}
			}
		}
	}
	// 历史最优更新
	for (int j = 0; j < history_best.second.size() - 1; ++j) {
		it_t = find(sites.begin(), sites.end(), history_best.second[j]);
		dis = distance(sites.begin(), it_t);//
		for (auto& s : pher_mtx[dis]) {
			if (s.first == history_best.second.at(j + 1)) { // 要更新路径的 目的地节点
				s.second = s.second * (1 - rho) + rho * (w * flow_mtx[dis].at(s.first) / (m * history_best.first));
				return;
			}
		}
	}
}

void updateLocalPher(myType a, myType b) {
	//vector<myType>::iterator it;
	//it = find(sites.begin(), sites.end(), a);
	int pos = get_pos(a);
	for (auto& s : pher_mtx[pos]) {
		if (s.first == b) {
			s.second = s.second * (1 - xik) + xik * pher_start;// 对信息素进行调整
			return;
		}
	}
	cerr << "updateLocalPher error" << endl;
}

//从a->b 的速度 T 时刻
float velocityCompute(myType a, myType b) {
	//vector<myType>::iterator it;
	//it = find(sites.begin(), sites.end(), a);
	//if (it == sites.end()) {
	//	cerr << "velocity 计算失败" << endl;
	//	return -1;
	//}
	int pos = get_pos(a);
	for (auto s : density_mtx[pos]) {
		if (s.first == b) {
			float weight = g.getweight(a, b);
			return v_start * exp(1 / static_cast<float>(s.second.second) * weight * (float)0.05 * (float)traffic_time);//计算公式 尽量使用static_cast转换或者不转换
		}
	}
	cerr << "velocityCompute error" << endl;
	return -1;
}

//从nyc bicycle dataset 读取数据 生成图
void readFile(string direction) {
	const int sites = 600;
	map<pair<int, int>, int > p;
	map<pair<int, int>, int> site_map;//哈希表 存出度
	ifstream fin(direction);
	string line;
	vector<vector<int> > Intact(sites, vector<int>(sites, 0));//sites  初始化
	const int ReadLine = 10000;
	for (int i = 0; i < ReadLine; ++i) {
		getline(fin, line);//读入一行数据
		//cout << "原始字符串：" << line << endl;
		istringstream sin(line);
		vector<int> res;
		string tmp;
		int t;
		for (int i = 0; i < 3; ++i) //将字符串流sin中的字符读入到field字符串中，以逗号为分隔符
		{
			getline(sin, tmp, ',');//读取三次
			t = stoi(tmp);
			res.push_back(t);
		}
		if (site_map.count({ res[1],res[2] }) == 0) site_map[{res[1], res[2]}] = 1; //出度
		else site_map.at({ res[1],res[2] })++;
		if (res[1] <= sites && res[2] < sites) {
			Intact[res[1]][res[2]] = res[0];
		}

		//cout << "处理之后的字符串：" << res[0] << ',' << res[1] << ',' << res[2] << endl;
		res.clear();
	}
	int cnt = 0, points = 0;
	int tmp;
	struct tmp1 //重写仿函数
	{
		bool operator() (vector<int> a, vector<int> b)
		{
			return a[1] < b[1]; //大顶堆
		}
	};
	priority_queue<vector<int>, vector<vector<int> >, tmp1>  min4;//大顶堆 存储四个当前最小值
	FILE* stream1;
	freopen_s(&stream1, "sites.txt", "w", stdout);//输出重定向
	cout << endl;
	for (int i = 0; i < sites; ++i) {
		bool tag = false;
		for (int j = i + 1; j < sites; ++j) {
			int number = 0;
			if (Intact[i][j] || Intact[j][i]) {
				if (!tag) {
					points++, tag = !tag;
				}
				if (site_map.count({ i,j }) != 0) number += site_map[{i, j}];
				if (site_map.count({ j,i }) != 0) number += site_map[{j, i}];
				if (number == 0) cerr << "line72 number is 0" << endl;
				tmp = Intact[i][j] + Intact[j][i] / number;
				if (min4.size() < 3) min4.push({ j,tmp });//if (min4.size() < 2) min4.push({ j,tmp });
				else if (min4.top()[1] > tmp) {
					min4.pop();
					min4.push({ j,tmp });
				}
			}// end if
		}// end j 
		while (!min4.empty())
		{
			cnt++;
			cout << i << " " << min4.top()[0] << " " << (float)min4.top()[1] / (float)50 << endl;
			min4.pop();
		}
	}// end i
	freopen_s(&stream1, "CON", "w", stdout);//输出重定向回 控制台
	ofstream fout("sites.txt", ios::out || ios::app);
	if (!fout.is_open()) cerr << "Can't open file" << endl;
	long location = 0;
	fout.seekp(location);
	//cout << cnt << endl << points << endl;
	fout << points << " " << cnt << " " << 1 << " " << 1 << " " << 2;
	fout.close();

	cout << "write finish" << "\a" << endl;
	return;
}

//v表示源节点,p表示v到各顶点的最短路径矩阵   
//res为二维矩阵  res[i][0]表示传送起始节点 res[i][1]表示传送终点
bool graphPreprocess(int v, int num, vector<vector<myType> >& p, vector<vector<myType> >& res) {
	for (int i = 0; i < num; ++i) {
		vector<int> tmb(2, -1);
		int n = p[i].size();
		if (n > 4) {
			float tmp, tmp_t;
			do {
				tmp = rand_n(n);
			} while (tmp == 0 || tmp == n - 1);
			tmb[0] = p[i][tmp];
			int k = n - tmp - 1;
			tmp_t = rand_n(k);
			tmb[1] = p[i][tmp + tmp_t + 1];
		}
		res.push_back(tmb);
	}
	return true;
}