#include "main.h"
#include <algorithm>
#include <bitset>
#include <iostream>
#include <queue>
#include <stack>
#include <unordered_map>
#include <unordered_set>
#include <vector>
using namespace std;

// Definition for singly - linked list.
struct ListNode {
  int val;
  ListNode *next;
  ListNode() : val(0), next(nullptr) {}
  ListNode(int x) : val(x), next(nullptr) {}
  ListNode(int x, ListNode *next) : val(x), next(next) {}
};

class Solution {
public:
  bool canConstruct(string ransomNote, string magazine) {
    int hash[26] = {0};
    for (int i = 0; magazine[i]; ++i) {
      hash[magazine[i] - 'a'] += 1;
    }
    for (int i = 0; ransomNote[i]; ++i) {
      hash[ransomNote[i] - 'a'] -= 1;
      if (hash[ransomNote[i] - 'a'] < 0)
        return false;
    }
    return true;
  }

  vector<int> intersection(vector<int> &nums1, vector<int> &nums2) {
    unordered_set<int> result_set;                           // 存放结果
    unordered_set<int> nums_set(nums1.begin(), nums1.end()); //初始化牛逼
    for (int num : nums2) {
      // 发现nums2的元素 在nums_set里又出现过
      if (nums_set.find(num) !=
          nums_set
              .end()) { //`find(key);`
                        //查找key是否存在,若存在，返回该键的元素的迭代器；若不存在，返回set.end();
        result_set.insert(num);
      }
    }
    return vector<int>(result_set.begin(), result_set.end());
  }

  vector<int> intersect(vector<int> &nums1, vector<int> &nums2) {
    if (nums1.size() > nums2.size()) {
      return intersect(nums2, nums1);
    }                          //保证nums1 size更小
    unordered_map<int, int> m; //记录出现的出现的数字和次数
    for (int num : nums1) {
      ++m[num];
    }
    vector<int> intersection;
    for (int num : nums2) {
      if (m.count(num)) {
        intersection.push_back(num);
        --m[num];
        if (m[num] == 0) {
          m.erase(num);
        }
      }
    }
    return intersection;
  }

  bool isHappy(int n) {
    unordered_set<int> set;
    while (1) {
      int sum = getSum(n);
      if (sum == 1)
        return 1;
      if (set.find(sum) != set.end()) // set中存在过
        return false;
      else
        set.insert(sum);
    }
  }

  int getSum(int n) {
    int sum = 0;
    while (n) {
      sum += (n % 10) * (n % 10);
      n /= 10;
    }
    return sum;
  }

  vector<int> twoSum(vector<int> &nums, int target) {
    std::unordered_map<int, int> map;
    for (int i = 0; i < nums.size(); i++) {
      auto iter = map.find(target - nums[i]);
      if (iter != map.end()) {
        return {iter->second, i};
      }
      map.insert(pair<int, int>(nums[i], i));
    }
    return {};
  }

  int fourSumCount(vector<int> &A, vector<int> &B, vector<int> &C,
                   vector<int> &D) {
    unordered_map<int, int> umap;
    // key:a+b的数值，value:a+b数值出现的次数
    // 遍历大A和大B数组，统计两个数组元素之和，和出现的次数，放到map中
    for (int a : A) {
      for (int b : B) {
        umap[a + b]++;
      }
    }
    int count =
        0; // 统计a+b+c+d = 0 出现的次数
           // 在遍历大C和大D数组，找到如果 0-(c+d)
           // 在map中出现过的话，就把map中key对应的value也就是出现次数统计出来。
    for (int c : C) {
      for (int d : D) {
        if (umap.find(0 - (c + d)) != umap.end()) {
          count += umap[0 - (c + d)];
        }
      }
    }
    return count;
  }

  vector<vector<int>> threeSum(vector<int> &nums) {
    vector<vector<int>> ans;
    if (nums.size() < 3 || nums.empty())
      return ans; // 特判
    int n = nums.size();

    sort(nums.begin(), nums.end()); //排序
    for (int i = 0; i < n; i++)     // 枚举最小值
    {
      if (nums[i] > 0)
        return ans;
      if (i > 0 && nums[i] == nums[i - 1])
        continue; // 最小元素去重！
      int l = i + 1;
      int r = n - 1;
      while (l < r) // 枚举中间值和最大值
      {
        int x = nums[l] + nums[r] + nums[i];
        if (x == 0) { // 符合条件，存储，并且去重，双端都移到下一个位置
          ans.push_back({nums[i], nums[l], nums[r]});
          while (l < r && nums[l] == nums[l + 1])
            l++;
          l++;
          while (l < r && nums[r] == nums[r - 1])
            r--;
          r--;
        } else if (x > 0) // 大了就让右边最大值变小
          r--;
        else // 小了就让左边中间值变大
          l++;
      }
    }
    return ans;
  }

  vector<vector<int>> fourSum(vector<int> &nums, int target) {}

  int findIntegers(int n) {
    int count = n + 1;
    for (int i = 0; i <= n; i++) {
      if (i & (i >> 1)) {
        count--;
      }
    }
    return count;
  }

  bool isValid(string s) {
    int n = s.size();
    if (n % 2 == 1) {
      return false;
    }

    unordered_map<char, char> pairs = {{')', '('}, {']', '['}, {'}', '{'}};
    stack<char> stk;
    for (char ch : s) {
      if (pairs.count(
              ch)) { // count函数用以统计key值在unordered_map中出现的次数。实际上，c++
                     // unordered_map不允许有重复的key。因此，如果key存在，则count返回1，如果不存在，则count返回0.
        if (stk.empty() || stk.top() != pairs[ch]) {
          return false;
        }
        stk.pop();
      } else {
        stk.push(ch);
      }
    }
    return stk.empty();
  }

  bool checkValidString(string s) {
    int n = s.length();
    if (n % 2 == 1) {
      if (s[n / 2] != '*') {
        return false;
      }
    }
    /*unordered_map<char, char> pairs = {{'(', ')'}};*/
    for (int i = 0; i < n / 2; i++) {
      if (s[i] == ')')
        return false;
      if ((s[i] == '(' && (s[n - 1 - i] == '*' || s[n - 1 - i] == ')')) ||
          (s[i] == '*' && (s[n - 1 - i] == '*' || s[n - 1 - i] == ')'))) {
        continue;
      } else {
        return false;
      }
    }
    return true;
  }

  //求取double（距离）的最大份数
  vector<int> getPices(vector<double> nums) {
    vector<int> res;
    double sum = 0;
    for (int i = 0; i < nums.size(); i++) {
      nums[i] *= 10;
      sum += nums[i];
    }
    double summ = sum;
    while (sum > 500) {
      sum /= 2;
    }
    for (int i = 0; i < nums.size(); i++) {
      int re = round(nums[i] * sum / summ);
      res.push_back(re);
    }
    return res;
  }

  //   vector<int> getPieces(vector<double> nums) {
  // 	  vector<int> res;
  // 	  for ()
  // 	  {
  // 	  }
  // 	  return res;
  //   }

  double chu(double a, double b) {
    if (a >= b) {
      return a / b;
    } else {
      return b / a;
    }
  }

  //   int numberOfBoomerangs(vector<vector<int>> &points) {
  //     int n = points.size();
  //     for (int i = 0; i < n; i++) {
  // 		for (int j = i; j<n ; j++)
  // 		{if ()
  // 		{
  // 		}
  // 		}
  //     }
  //   }
  //   int square_length(vector<int> point1, vector<int> point2) {
  //     return abs(point1[0] - point2[0]) * abs(point1[1] - point2[1]);
  //   }

  bool biggerstring(const string &a, const string &b) {
    return a.length() > b.length();
  }
  string findLongestWord(string s, vector<string> &dictionary) {
    unordered_map<int, int> s_map;
    int ns = s.length();
    int n = dictionary.size();
    for (int i = 0; i < ns; i++) {
      s_map[s[i]]++;
    }
    //     sort(dictionary.begin(), dictionary.end());
    //     reverse(dictionary.begin(), dictionary.end());
    vector<string> tmp = dictionary;
    // sort(tmp.begin(), tmp.end(), biggerstring);
    vector<string> cansize;
    for (int i = 0; i < n; i++) {
      for (int j = 0; i < dictionary[i].size(); j++) {
        s_map[dictionary[i][j]]--;
        if (s_map[dictionary[i][j]] < 0) {
          break;
        }
        if (j == dictionary[i].size() - 1) {
          return dictionary[i];
        }
      }
    }
  }

  string reverseStr(string s, int k) {
    int n = s.length();
    int round = n / 2 / k;
    int left = n % (2 * k);
    for (int i = 1; i <= round; i++) {
      reverse(s.begin() + k * i - k, s.begin() + k * i);
    }
    reverse(s.begin() + 2 * k, s.begin() + 3 * k);
    cout << s;
    return s;
  }

  string replaceSpace(string s) {
    int n = s.size();
    for (int i = 0; i < n; i++) {
      if (s[i] == ' ') {
        s.replace(i, 1, "%20");
      }
    }
    return s;
  }

  string reverseWords(string s) {
    vector<string> res;
    int n = s.length();
    int left, right;
    bool flag = false;
    for (int i = 0; i < n; i++) {
      if (s[i] == ' ') {
        continue;
      }
      if (s[i - 1] == ' ') {
        left = i;
      } else if (s[i + 1] == ' ') {
        right = i + 1;
      }
    }
  }

  int strStr(string haystack, string needle) {
    int nee = needle.length();
    int hay = haystack.size();
    int samesize = 0;
    for (int i = 0; i < nee; i++) {
      for (int j = 0; j < hay; j++) {
        if (haystack[j] == needle[i]) {
          samesize++;
          if (samesize == nee)
            return j - 1;
        } else
          samesize = 0;
      }
    }
    return -1;
  }

  string removeDuplicates(string s) {
    bool flag = false;
    stack<char> temp;
    char tempch;
    vector<char> res;
    string resu = {};
    for (int i = 0; i < s.size(); i++) {
      if (i == 0) {
        temp.push(s[i]);
        tempch = s[i];
        continue;
      }
      if (tempch == s[i] && !temp.empty()) {
        temp.pop();
        flag = true;
      } else {
        tempch = s[i];
        temp.push(s[i]);
      }
    }
    while (!temp.empty()) {
      res.push_back(temp.top());
      temp.pop();
    }
    for (int i = res.size() - 1; i >= 0; i--) {
      resu += res[i];
    }
    if (flag == true) {
      return removeDuplicates(resu);
    } else
      return resu;
  }

  vector<int> maxSlidingWindow(vector<int> &nums, int k) {
    int n = nums.size();
    vector<int> res;
    queue<int> temp;
    for (int i = 0; i < n - k + 1; i++) {
      for (int j = 0; j < k; j++) {
        temp.push(nums[i + j]);
      }
      res.push_back(getMaxFromQueue(temp));
      while (!temp.empty()) {
        temp.pop();
      }
    }
    return res;
  }

  int getMaxFromQueue(queue<int> que) {
    vector<int> res;
    while (!que.empty()) {
      res.push_back(que.front());
      que.pop();
    }
    return *std::max_element(res.begin(), res.end());
  }

  vector<int> maxSlidingWindow2(vector<int> &nums, int k) {
    int n = nums.size();
    vector<int> res;
    queue<int> temp;
    for (int i = 0; i < n - k + 1; i++) {
      if (i == 0) {
        for (int j = 0; j < k; j++) {
          temp.push(nums[i + j]);
        }
      } else {
        temp.pop();
        temp.push(nums[i + k - 1]);
      }
      res.push_back(getMaxFromQueue(temp));
    }
    return res;
  }

  vector<int> maxSlidingWindow3(vector<int> &nums, int k) {

    vector<int> res;
    deque<int> deque;

    for (int i = 0; i < nums.size(); i++) {

      if (!deque.empty() && deque.front() == i - k)
        deque.pop_front();

      while (!deque.empty() && nums[i] > nums[deque.back()])
        deque.pop_back();

      deque.push_back(i);

      if (i >= k - 1)
        res.push_back(nums[deque.front()]);
    }

    return res;
  }

  //给你一个整数数组 nums 和一个整数 k ，请你返回其中出现频率前 k
  //高的元素。你可以按 任意顺序 返回答案
  //输入: nums = [1,1,1,2,2,3], k = 2  输出: [1, 2]
  vector<int> topKFrequent(vector<int> &nums, int k) {
    vector<int> res;
    int n = nums.size();
    unordered_map<int, int> iMap;
    for (int i = 0; i < n; i++) {
      iMap[nums[i]]++;
    }
    vector<pair<int, int>> vtMap;
    for (auto it = iMap.begin(); it != iMap.end(); it++)
      vtMap.push_back(make_pair(it->first, it->second));

    sort(vtMap.begin(), vtMap.end(),
         [](const pair<int, int> &x, const pair<int, int> &y) -> int {
           return x.second > y.second;
         });
    for (int i = 0; i < k; i++) {
      res.push_back(vtMap[i].first);
    }
    return res;
  }
};

void solve() {
  Solution mysolution;
  vector<int> nums1 = {1, 3, -1, -3, 5, 3, 6, 7};
  int target = 3;
  string a = "aaaaaaa";
  string b = "aa";
  vector<double> test = {1.11, 2.5, 3.94, 10, 11.5, 11.111, 12.53, 10};
  vector<string> dictionary = {"ale", "apple", "monkey", "plea"};
  string s = "the sky is blue";
  string s1 = "aaaaa";
  string s2 = "abbcd";
  vector<int> str = mysolution.maxSlidingWindow3(nums1, target);
}
int main() {
  solve();
  return 0;
}
