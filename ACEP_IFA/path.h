#ifndef SHORTPATH_H
#define SHORTPATH_H
#include"Graph.h"
#include<vector>
#include<string>
#include<stack>
using namespace std;
void dijkstra(int num, int v, vector<int>& p, vector<int>& d, const vector<vector<int>>& g);
vector<vector<int> > shortpath_construct(int v,int num, vector<int>& p);
#endif