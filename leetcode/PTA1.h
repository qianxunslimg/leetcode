 


#include <iostream>
#include <cstdio>
#include "string.h"
#include <string>
#include <vector>
#include <algorithm>
#include "math.h"
#include <cstring>
#include <stdio.h>
#include <stdlib.h>
#include "unordered_map"
#include "map"
#include "queue"
#include "sstream"
#include "unordered_set"
#include "windows.h"
#include "stack"
#include "numeric"
#include "set"

#define     O_CREAT                         0x0100
#define		O_RDONLY						0x0000
#define		O_WRONLY						0x0001
#define		O_RDWR						    0x0002
#define     O_BINARY                        0x8000  // file mode is binary (untranslated)
#define		O_TRUNC                         0x0200

#include "dirent.h"
#include "io.h"

using namespace std;

struct ListNode {
	int val;
	ListNode *next;
	ListNode() : val(0), next(nullptr) {}
	ListNode(int x) : val(x), next(nullptr) {}
	ListNode(int x, ListNode *next) : val(x), next(next) {}
 };

struct TreeNode {
	int val;
	TreeNode *left;
	TreeNode *right;
	TreeNode() : val(0), left(nullptr), right(nullptr) {}
	TreeNode(int x) : val(x), left(nullptr), right(nullptr) {}
	TreeNode(int x, TreeNode *left, TreeNode *right) : val(x), left(left), right(right) {}	
};

//class Node {
//public:
//	int val;
//	vector<Node*> children;
//
//	Node() {}
//
//	Node(int _val) {
//		val = _val;
//	}
//
//	Node(int _val, vector<Node*> _children) {
//		val = _val;
//		children = _children;
//	}
//};

//class Node {
//public:
//	int val;
//	Node* left;
//	Node* right;
//	Node* next;
//
//	Node() : val(0), left(NULL), right(NULL), next(NULL) {}
//
//	Node(int _val) : val(_val), left(NULL), right(NULL), next(NULL) {}
//
//	Node(int _val, Node* _left, Node* _right, Node* _next)
//		: val(_val), left(_left), right(_right), next(_next) {}
//};
//
//class Node {
//public:
//	int val;
//	Node* next;
//	Node* random;
//
//	Node(int _val) {
//		val = _val;
//		next = NULL;
//		random = NULL;
//	}
//};

class Node {
public:
	int val;
	Node* prev;
	Node* next;
	Node* child;

	Node(int _val) {
		val = _val;
		prev = NULL;
		next = NULL;
		child = NULL;
	}

	Node(int _val, Node* pre, Node* nex, Node* chil) {
		val = _val;
		prev = pre;
		next = nex;
		child = chil;
	}
};

TreeNode* createTree(vector<int> list, int start) {
	if (list[start] == '#') {
		return NULL;
	}

	TreeNode* root = new TreeNode(list[start]);

	int lnode = 2 * start + 1;
	int rnode = 2 * start + 2;
	if (lnode > list.size() - 1) {
		root->left = NULL;
	}
	else {
		root->left = createTree(list, lnode);
	}

	if (rnode > list.size() - 1) {
		root->right = NULL;
	}
	else {
		root->right = createTree(list, rnode);
	}

	return root;
}

ListNode* createList(vector<int> arr)
{
	ListNode* dummy = new ListNode(0);
	ListNode* cur = dummy;
	for (int i = 0; i < arr.size(); i++) {
		ListNode* tmp = new ListNode(arr[i]);
		cur->next = tmp;
		cur = cur->next;
	}
	return dummy->next;
}


Node* createNodeList(vector<int> arr) 
{
	Node* dummy = new Node(0);
	Node* cur = dummy;
	for (int i = 0; i < arr.size(); i++) {
		Node* tmp = new Node(arr[i]);
		cur->next = tmp;
		tmp->prev = cur;
		cur->child = nullptr;
		cur = cur->next;
	}
	return dummy->next;
}

void splitWithStringStream(const string& str, char delim, vector<string>& ret)
{
	stringstream ss(str);
	string s;
	while (getline(ss, s, delim))
	{
		ret.push_back(s);
	}
}

vector<string> split(string str)
{
	vector<string> res;
	stringstream ss(str);
	while (ss >> str) {
		res.push_back(str);
	}
	return res;
}


#define  MAX_DAQ_NUM 10
/*HRFK文件转浮点数*/
int decode_hrfkfile()
{
	char buf[16];
	char name[96];

	FILE *fin = fopen("F:/test_hv/N84.disp", "rb");

	if (fin == NULL) {
		//LOG_E("can not open file %s", filepath);
		return -1;
	}

	sprintf(name, "F:/test_hv/N84_1.disp");   /*删除前缀spark*/

	FILE *fout = fopen(name, "wt+");
	if (fout == NULL) {
		//LOG_E("can not open file %s", name);
		return -1;
	}

	int len;

	while (len = fread(buf, 1, 16, fin) > 0) {
		fprintf(fout, "%.6f %.6f\n", *(double *)(buf), *(double *)(buf + 8));
	}

	fclose(fin);
	fclose(fout);

	return 0;
}

/*拆分解码HV文件,ID号为当前连接设备*/
int decode_split_hvfile()
{
	int cnt = 0;    //有效连接设备数目
	char index[MAX_DAQ_NUM + 1] = { 0 };

	for (int i = 1; i <= 6; i++) {
		index[i] = 1;
		cnt++;
	}

	struct stat status;
	stat("F:/test_hv/N84.hv", &status);
	long size = status.st_size;

	//for (int i = 0; i < MAX_DAQ_NUM; i++) {
	//	unsigned char dev = dev_arr[i];
	//	if (dev == 0) break;
	//	index[dev] = 1;
	//	cnt++;
	//}

	/*HV文件拆分*/
	FILE *fin = fopen("F:/test_hv/N84.hv", "rb");

	if (fin == NULL) {
		//LOG_E("can not open file %s", mergedhv_path);
		return -1;
	}

	char name[96] = { 0 };
	char *buf = (char *)malloc(1024);
	int i = 1;

	while (i <= MAX_DAQ_NUM) {
		int sum = 0;
		if (index[i] == 0) {
			i++;
			continue;
		}

		sprintf(name, "F:/test_hv/N84_%d.hv", i);
		i++;

		FILE *fout = fopen(name, "wt+");
		if (fout == NULL) {
			//LOG_E("can not open file %s", name);
			return -1;
		}

		unsigned short integer = size / cnt / 1024;  //分包数
		unsigned short remainder = (size / cnt) % 1024;
		int k = 0;
		int len = 1024;

		/*每次只读取一个分hv文件大小的数据*/
		while (sum < size / cnt) {
			if (k == integer) len = remainder;    //不一定是1024倍数

			fread(buf, 1, len, fin);
			sum += len;
			k++;

			/*直接转成4列Hv数据*/
			for (int j = 0; j < len / 16; j++) {
				fprintf(fout, "%.6f %.6f %.6f %.6f\n",
					*(float *)(buf + 16 * j),
					*(float *)(buf + 16 * j + 4),
					*(float *)(buf + 16 * j + 8),
					*(float *)(buf + 16 * j + 12));
			}
		}
		fclose(fout);
	}

	fclose(fin);
	free(buf);

	return 0;
}