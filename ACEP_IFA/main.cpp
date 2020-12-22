#include"path.h"
#include "Graph.h"
#include"route_a.h"
using namespace std;
UndiGraph<myType, float> g(310);

const int m = 100;//ants number
const float beta = 2.0;//heuristic information weight 
int G = 50;//迭代次数
const float rho = 0.3;//pheromone violation factor 
const float xik = 0.1;//
const float pher_start = 0.025;
const float v_start = 100;
const int w = 3;//number of ants allowed to release pheromone 
const int cap = 20;//capacity of each road
int traffic_time = 1; // 用以计算速度的时间
const float heur_strat = 0;

const int target = 48;//设定初始固定源节点  48//521
vector<myType> source(m, target);//每只蚂蚁的源节点  521
vector<myType> safe_site;//安全的节点  或 目的节点集合
vector<myType> sites;//全部节点的集合
vector<bool> ant_dest(m, false);//是否到达终点标识
vector<vector<myType> > ant_solution(m);
vector<unordered_map<myType, float> > flow_mtx;//flow hash function

vector<vector<pair<myType, float> > > heur_mtx;//heuristic information 
vector<vector<pair<myType, float> > > trans_mtx;//Transition matrix
vector<vector<pair<myType, float> > > pher_mtx;//Pheromone matrix

set<pair<int, int> > cur_global_route;
vector<vector<pair<myType, pair<int, int> > > > density_mtx;//density matrix  与flow 同时实现 最外层为节点编号
vector<float> fit_mtx(m, 0);//适应度矩阵
pair<float, vector<myType>> history_best;
vector<vector<myType> > transmit_p;// 传送矩阵
vector<vector<myType> > short_p;//short path

int main() {
	//readFile("F:/Edge_Download/Data/Data/NYC/newtest.csv");
	//SetConsoleTitle("ACEP");
	g.InputfileGraph("Dj.txt");//   C:/Users/shijian/Desktop/   sites.txt(2020/7/15之前) Dj.txt
	g.Display();
	new_init();
	auto k = g.getedgematrix();
	int num = g.getnumVertexes();
	vector<int> d(num);
	vector<int> p(num, -1);
	dijkstra(num, target, p, d, k);//构建最小生成树
	short_p = shortpath_construct(target, num, p);//最短路径的完整形式  
	//vector<vector<myType> > res;
	if (graphPreprocess(target, num, short_p, transmit_p))  cout << "graphPreprocess success\n";
	else cout << "graphPreprocess failed\n";
	FILE* stream2;
	freopen_s(&stream2, "res4.txt", "w", stdout);//输出重定向
	//g.Display();
	for (int cnts = 0;G;cnts++) {
		//freopen_s(&stream2, "CON", "w", stdout);//输出重定向回 控制台
		routeConstricution();
		updateGlobalPher();
		updateTrans();
		float sum = meanCompute();
		//cout << '\a' << endl;
		//freopen_s(&stream2, "res2.txt", "w", stdout);//输出重定向
		cout << sum<< endl;
		//cout << cnts<<":"<<meanCompute() << endl;
	}
	freopen_s(&stream2, "CON", "w", stdout);//输出重定向回 控制台
	cout << '\a';
	return 0;
}

