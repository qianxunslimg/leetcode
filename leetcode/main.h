/*!
 * \file main.h
 *
 * \author qianxunslimg
 * \date 十一月 2021
 *
 *
 */
#pragma once
#include "map"
#include "math.h"
#include "numeric"
#include "queue"
#include "set"
#include "sstream"
#include "stack"
#include "string.h"
#include "unordered_map"
#include "unordered_set"
#include "windows.h"
#include <algorithm>
#include <bitset>
#include <cstdio>
#include <cstring>
#include <functional>
#include <iostream>
#include <map>
#include <math.h>
#include <memory>
#include <queue>
#include <set>
#include <sstream>
#include <stack>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <tuple>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#define O_CREAT 0x0100
#define O_RDONLY 0x0000
#define O_WRONLY 0x0001
#define O_RDWR 0x0002
#define O_BINARY 0x8000 // file mode is binary (untranslated)
#define O_TRUNC 0x0200

//#include "io.h"

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
  TreeNode(int x, TreeNode *left, TreeNode *right)
      : val(x), left(left), right(right) {}
};

// class Node {
// public:
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

// class Node {
// public:
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
// class Node {
// public:
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
  Node *prev;
  Node *next;
  Node *child;

  Node(int _val) {
    val = _val;
    prev = NULL;
    next = NULL;
    child = NULL;
  }

  Node(int _val, Node *pre, Node *nex, Node *chil) {
    val = _val;
    prev = pre;
    next = nex;
    child = chil;
  }
};

TreeNode *createTree(vector<int> list, int start) {
  if (list[start] == '#') {
    return NULL;
  }

  TreeNode *root = new TreeNode(list[start]);

  int lnode = 2 * start + 1;
  int rnode = 2 * start + 2;
  if (lnode > list.size() - 1) {
    root->left = NULL;
  } else {
    root->left = createTree(list, lnode);
  }

  if (rnode > list.size() - 1) {
    root->right = NULL;
  } else {
    root->right = createTree(list, rnode);
  }

  return root;
}

ListNode *createList(vector<int> arr) {
  ListNode *dummy = new ListNode(0);
  ListNode *cur = dummy;
  for (int i = 0; i < arr.size(); i++) {
    ListNode *tmp = new ListNode(arr[i]);
    cur->next = tmp;
    cur = cur->next;
  }
  return dummy->next;
}

Node *createNodeList(vector<int> arr) {
  Node *dummy = new Node(0);
  Node *cur = dummy;
  for (int i = 0; i < arr.size(); i++) {
    Node *tmp = new Node(arr[i]);
    cur->next = tmp;
    tmp->prev = cur;
    cur->child = nullptr;
    cur = cur->next;
  }
  return dummy->next;
}
