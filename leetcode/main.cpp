#include "main.h"
#define M_PI 3.1415926535
/*!
 * \class Solution
 *
 * \brief
 *
 * \author qianxunslimg
 * \date 十一月 2021
 */
class Solution {
public:
  struct TreeNode {
    int val;
    TreeNode *left;
    TreeNode *right;
    TreeNode() : val(0), left(nullptr), right(nullptr) {}
    TreeNode(int x) : val(x), left(nullptr), right(nullptr) {}
    TreeNode(int x, TreeNode *left, TreeNode *right)
        : val(x), left(left), right(right) {}
  };

  vector<string> findRelativeRanks2(vector<int> &score) {
    int n = score.size();
    map<int, int, greater<int>> num2index;
    for (int i = 0; i < n; i++)
      num2index[score[i]] = i;
    vector<string> ans(n);
    int i = 0;
    for (auto mPair : num2index) {
      int index = mPair.second;
      if (i == 0) {
        ans[index] = "Gold Medal";
      } else if (i == 1) {
        ans[index] = "Silver Medal";
      } else if (i == 2) {
        ans[index] = "Bronze Medal";
      } else
        ans[index] += to_string(i + 1);
      i++;
    }
    return ans;
  }

  static bool sortPair(pair<int, int> a, pair<int, int> b) {
    return a.second > b.second;
  }
  vector<string> findRelativeRanks(vector<int> &score) {
    vector<string> res(score.size());

    vector<pair<int, int>> pairr;

    for (int i = 0; i < score.size(); i++) {
      pairr.push_back(pair<int, int>(i, score[i]));
    }

    sort(pairr.begin(), pairr.end(), sortPair);
    for (int i = 0; i < pairr.size(); i++) {
      switch (i) {
      case 0:
        res[pairr[i].first] = "Gold Medal";
        break;
      case 1:
        res[pairr[i].first] = "Silver Medal";
        break;
      case 2:
        res[pairr[i].first] = "Bronze Medal";
        break;
      default:
        res[pairr[i].first] = to_string(i + 1);
        break;
      }
    }
    return res;
  }

  //************************************
  // Method:    findTheDifference
  // FullName:  Solution::findTheDifference
  // Access:    public
  // Returns:   char
  // Qualifier:
  // Parameter: string s
  // Parameter: string t
  //************************************
  char findTheDifference(string s, string t) {
    sort(s.begin(), s.end());
    sort(t.begin(), t.end());
    for (int i = 0; i < s.size(); i++) {
      if (s[i] != t[i])
        return t[i];
    }
    return t[t.size() - 1];
  }

  //计数
  char findTheDifference2(string s, string t) {
    vector<int> cnt(26, 0);
    for (int i = 0; i < s.size(); i++) {
      cnt[s[i] - 'a']++;
    }
    for (int i = 0; i < t.size(); i++) {
      cnt[t[i] - 'a']--;
      if (cnt[t[i] - 'a'] < 0) {
        return t[i];
      }
    }
  }

  //求和
  char findTheDifference3(string s, string t) {
    int ac = 0;
    for (char ch : s) {
      ac -= ch;
    }
    for (char ch : t) {
      ac += ch;
    }
    return ac;
  }

  // runtime error
  int integerReplacement(int n) {
    int res = 0;
    while (n != 1) {
      if (n == 3) {
        res += 2;
        break;
      }
      if (n % 2) {
        if (!((n + 1) % 4)) {
          n = n + 1;
        } else {
          n = n - 1;
        }
      } else {
        n /= 2;
      }
      res++;
    }
    return res;
  }

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

  int arrangeCoins(int n) {
    return (int)((sqrt((long long)8 * n + 1) - 1) / 2);
  }

  int lengthOfLongestSubstring(string s) {
    if (s.size() == 0) {
      return 0;
    }
    vector<string> ss;
    for (int i = 0; i < s.size(); i++) {
      unordered_map<char, int> temp_map;
      for (int j = i; j < s.size(); j++) {
        ++temp_map[s[j]];
        if (temp_map[s[j]] > 1) {
          string sss = s.substr(i, j - i);
          ss.push_back(sss);
          break;
        }
        if (j == s.size() - 1) {
          string sss = s.substr(i, s.size() - i);
          ss.push_back(sss);
        }
      }
    }
    if (ss.size() == 0) {
      return s.size();
    }
    sort(ss.begin(), ss.end(),
         [](string &a, string &b) { return a.size() > b.size(); });

    return ss[0].size();
  }

  int lengthOfLongestSubstring2(string s) {
    // 哈希集合，记录每个字符是否出现过
    unordered_set<char> occ;
    int n = s.size();
    // 右指针，初始值为 -1，相当于我们在字符串的左边界的左侧，还没有开始移动
    int rk = -1, ans = 0;
    // 枚举左指针的位置，初始值隐性地表示为 -1
    for (int i = 0; i < n; ++i) {
      if (i != 0) {
        // 左指针向右移动一格，移除一个字符
        occ.erase(s[i - 1]);
      }
      while (rk + 1 < n && !occ.count(s[rk + 1])) {
        // 不断地移动右指针
        occ.insert(s[rk + 1]);
        ++rk;
      }
      // 第 i 到 rk 个字符是一个极长的无重复字符子串
      ans = max(ans, rk - i + 1);
    }
    return ans;
  }

  int lengthOfLongestSubstring3(string s) {
    if (s.size() == 0)
      return 0;
    unordered_set<char> lookup;
    int maxStr = 0;
    int left = 0;
    for (int i = 0; i < s.size(); i++) {
      while (lookup.find(s[i]) != lookup.end()) {
        lookup.erase(s[left]);
        left++;
      }
      maxStr = max(maxStr, i - left + 1);
      lookup.insert(s[i]);
    }
    return maxStr;
  }

  double findMedianSortedArrays(vector<int> &nums1, vector<int> &nums2) {
    vector<int> res = nums1;
    res.insert(res.end(), nums2.begin(), nums2.end());
    sort(res.begin(), res.end());
    int n = res.size();
    if (n % 2) {
      return (double)(res[n / 2] + res[n / 2 - 1]) / 2;
    } else {
      return res[n / 2];
    }
  }

  int countSubstrings(string s) {
    int n = s.size(), ans = 0;
    for (int i = 0; i < 2 * n - 1; ++i) {
      int l = i / 2, r = i / 2 + i % 2;
      while (l >= 0 && r < n && s[l] == s[r]) {
        --l;
        ++r;
        ++ans;
      }
    }
    return ans;
  }

  string longestPalindrome(string s) {
    int n = s.size();
    vector<string> ans;
    for (int i = 0; i < 2 * n - 1; ++i) {
      int l = i / 2, r = i / 2 + i % 2;
      while (l >= 0 && r < n && s[l] == s[r]) {
        --l;
        ++r;
      }
      if (r - l == 2) {
        ans.push_back(s.substr(l + 1, 1));
      } else
        ans.push_back(s.substr(l + 1, r - l - 1));
    }
    sort(ans.begin(), ans.end(),
         [](string &a, string &b) { return a.size() > b.size(); });
    return ans[0];
  }

  //输入：nums = [0, 2, 3, 4, 6, 8, 9]
  // 输出：["0", "2->4", "6", "8->9"]
  // 解释：区间范围是：
  // [0, 0] -- > "0"
  // [2, 4] -- > "2->4"
  // [6, 6] -- > "6"
  // [8, 9] -- > "8->9"
  vector<string> summaryRanges(vector<int> &nums) {
    string jiantou = "->";
    int size = nums.size();

    vector<string> res;
    string temp;
    for (int i = 0; i < size; i++) {
      //其实判断
      if (i == 0 || (i > 0 && (nums[i] - 1 != nums[i - 1]))) {
        temp = "";
        temp += to_string(nums[i]);
      }
      //终止判断
      if (nums[i] + 1 != nums[i + 1]) {
        res.push_back(temp);
      }
      if (i > 0 && nums[i] + 1 == nums[i + 1] && nums[i] - 1 == nums[i - 1]) {
      }
    }
    return res;
  }

  int lengthOfLastWord(string s) {
    int n = s.size();
    removeSpace(s);
    int count = 0;
    while (s[s.size() - 1] != ' ') {
      count++;
      s.pop_back();
      if (s.size() == 0)
        return count;
    }
    return count;
  }

  void removeSpace(string &s) {
    while (s[s.size() - 1] == ' ') {
      s.pop_back();
      if (s.size() == 0) {
        return;
      }
      return;
    }
  }

  int majorityElement(vector<int> &nums) {
    int n = nums.size();
    if (n == 0) {
      return 0;
    }
    unordered_map<int, int> res;
    for (int i = 0; i < n; i++) {
      res[nums[i]]++;
      if (res[nums[i]++] > n / 2) {
        return res[nums[i]++];
      }
    }
    return 0;
  }
  //   int cal1plus2n(int n) {
  //     bool a[n][n + 1];
  //     return sizeof(a) >> 1;
  //   }

  vector<int> plusOne(vector<int> &digits) {
    int n = digits.size();
    for (int i = n - 1; i >= 0; --i) {
      if (digits[i] == 9) {
        digits[i] = 0;
      } else {
        digits[i] += 1;
        break;
      }
    }
    if (digits[0] == 0) {
      // digits.insert(digits.begin(), 1);
      digits[0] = 1;
      digits.push_back(0);
    }
    return digits;
  }

  vector<int> nextGreaterElement(vector<int> &nums1, vector<int> &nums2) {
    unordered_map<int, int> hashmap;
    stack<int> st;
    for (int i = nums2.size() - 1; i >= 0; --i) {
      int num = nums2[i];
      while (!st.empty() && num >= st.top()) {
        st.pop();
      }
      hashmap[num] = st.empty() ? -1 : st.top();
      st.push(num);
    }
    vector<int> res(nums1.size());
    for (int i = 0; i < nums1.size(); ++i) {
      res[i] = hashmap[nums1[i]];
    }
    return res;
  }

  vector<string> findWords(vector<string> &words) {
    string s1 = "qwertyuiop";
    string s2 = "asdfghjkl";
    string s3 = "zxcvbnm";
    unordered_map<int, int> map;
    for (int i = 0; i < s1.size(); i++) {
      map[s1[i]] = 1;
    }
    for (int i = 0; i < s2.size(); i++) {
      map[s2[i]] = 2;
    }
    for (int i = 0; i < s3.size(); i++) {
      map[s3[i]] = 3;
    }

    vector<string> res;
    for (int i = 0; i < words.size(); i++) {
      string tmp = words[i];
      set<int> temp;
      transform(words[i].begin(), words[i].end(), words[i].begin(), ::tolower);
      for (int j = 0; j < words[i].size(); j++) {
        temp.insert(map[words[i][j]]);
        cout << words[i][j];
        cout << map[words[i][j]];
        cout << temp.size() << endl;
        if (temp.size() != 1) {
          break;
        }
      }
      if (temp.size() == 1) {
        res.push_back(tmp);
      }
    }

    return res;
  }

  int longestSubsequence(vector<int> &arr, int difference) {
    int ans = 0;
    unordered_map<int, int> dp;
    for (int v : arr) {
      dp[v] = dp[v - difference] + 1;
      ans = max(ans, dp[v]);
    }
    return ans;
  }

  int missingNumber(vector<int> &nums) {
    int n = nums.size();
    if (n == 0) {
      return 1;
    }
    sort(nums.begin(), nums.end());
    for (int i = 0; i < nums.size(); i++) {
      if (nums[i] != i) {
        return i;
      }
    }
    return nums[n - 1] + 1;
  }

  /*方法二：哈希集合
使用哈希集合，可以将时间复杂度降低到
   * O(n)O(n)。

首先遍历数组
   * \textit{nums}nums，将数组中的每个元素加入哈希集合，然后依次检查从 00 到 nn
   * 的每个整数是否在哈希集合中，不在哈希集合中的数字即为丢失的数字。由于哈希集合的每次添加元素和查找元素的时间复杂度都是
   * O(1)O(1)，因此总时间复杂度是
   * O(n)O(n)。*/
  int missingNumber2(vector<int> &nums) {
    int n = nums.size();
    unordered_set<int> set;
    for (int i = 0; i < n; i++) {
      set.insert(nums[i]);
    }
    for (int i = 0; i <= n; i++) {
      if (set.find(i) == set.end()) {
        return i;
      }
    }
    return 0;
  }

  int missingNumber3(vector<int> &nums) {
    int res = 0;
    int n = nums.size();
    for (int i = 0; i < n; i++) {
      res ^= nums[i];
    }
    for (int i = 0; i <= n; i++) {
      res ^= i;
    }
    return res;
  }

  vector<int> findDuplicates(vector<int> &nums) {
    /*int n = nums.size();
    vector<int> res;
    unordered_map<int,int> map;
    for(int i = 0; i<n; i++){
    map[nums[i]]++;
    }
    for(int i = 0; i<=n; i++){
    if (map[i] == 2){
    res.push_back(i);
    }
    }
    return res;*/

    int n = nums.size();
    vector<int> res;
    sort(nums.begin(), nums.end());
    auto it = unique(nums.begin(), nums.end());
    for (; it != nums.end(); it++) {
      res.push_back(*it);
    }
    return res;
  }

  int maxCount(int m, int n, vector<vector<int>> &ops) {
    unordered_map<int, int> map1;
    unordered_map<int, int> map2;
    int nn = ops.size();
    for (int i = 0; i < nn; i++) {
      while (ops[i][0]) {
        map1[ops[i][0]]++;
        ops[i][0]--;
      }
      while (ops[i][1]) {
        map2[ops[i][1]]++;
        ops[i][1]--;
      }
    }
    auto x = std::max_element(
        map1.begin(), map1.end(),
        [](const pair<int, int> &p1, const pair<int, int> &p2) {
          return p1.second > p2.second;
        });
    auto x2 = std::max_element(
        map2.begin(), map2.end(),
        [](const pair<int, int> &p1, const pair<int, int> &p2) {
          return p1.second > p2.second;
        });

    cout << x->first << " " << x2->first;
    return ((x->first - 1) > m ? m : (x->first - 1)) *
           ((x2->first - 1) > n ? n : (x2->first - 1));
  }
  // 495. 提莫攻击
  int findPoisonedDuration(vector<int> &timeSeries, int duration) {
    int n = timeSeries.size();
    int alltimes = 0;
    for (int i = 0; i < n; i++) {
      alltimes += duration;
      if (i < n - 1) {
        alltimes -= (timeSeries[i + 1] - timeSeries[i]) < duration
                        ? (duration - (timeSeries[i + 1] - timeSeries[i]))
                        : 0;
      }
    }
    return alltimes;
  }
  // 649. Dota2 参议院
  string predictPartyVictory(string senate) {
    int n = senate.size();
    unordered_map<char, int> mapp;
    for (int i = 0; i < n; i++) {
      mapp[senate[i]]++;
    }
    if (mapp['R'] != mapp['D']) {
      return mapp['R'] > mapp['D'] ? "Radiant" : "Dire";
    } else {
      return senate[0] == 'R' ? "Radiant" : "Dire";
    }
  }
  bool isUnique(string s) {
    set<char> sett;
    for (int i = 0; i < s.size(); i++) {
      sett.insert(s[i]);
    }
    return sett.size() == 1 ? 1 : 0;
  }

  // 45. 跳跃游戏
  int jump(vector<int> &nums) {
    if (nums.size() == 1)
      return 0;
    int reach = 0;
    int nextreach = nums[0];
    int step = 0;
    for (int i = 0; i < nums.size(); i++) {
      nextreach = max(i + nums[i], nextreach);
      if (nextreach >= nums.size() - 1)
        return (step + 1);
      if (i == reach) {
        step++;
        reach = nextreach;
      }
    }
    return step;
  }
  int jump2(vector<int> &nums) {
    int ans = 0;
    int end = 0;
    int maxPos = 0;
    for (int i = 0; i < nums.size() - 1; i++) {
      maxPos = max(nums[i] + i, maxPos);
      if (i == end) {
        end = maxPos;
        ans++;
      }
    }
    return ans;
  }

  int getMoneyAmount(int n) {
    vector<vector<int>> f(n + 1, vector<int>(n + 1));
    for (int i = n - 1; i >= 1; i--) {
      for (int j = i + 1; j <= n; j++) {
        int minCost = INT_MAX;
        for (int k = i; k < j; k++) {
          int cost = k + max(f[i][k - 1], f[k + 1][j]);
          minCost = min(minCost, cost);
        }
        f[i][j] = minCost;
      }
    }
    return f[1][n];
  }

  int game(vector<int> &guess, vector<int> &answer) {
    int res = 0;
    for (int i = 0; i < 3; i++) {
      if (guess[i] == answer[i]) {
        res++;
      }
    }
    return res;
  }

  //最小堆操作
  int headTest() {
    vector<int> ss{0, 1, 2, 3, 4, 5, 6};
    pop_heap(ss.begin(), ss.end(), greater<int>());
    int tmp = ss.back();
    ss.pop_back();

    ss.push_back(10);
    push_heap(ss.begin(), ss.end(), greater<int>());
    return 0;
  }

  //教训：少分情况讨论 多使用min函数
  int minTimeToType(string word) {
    int res = 0;
    char pre = 'a';
    for (char w : word) {
      int a, b;
      if (pre < w) {
        a = w - pre;
        b = pre - w + 26;
      } else {
        a = pre - w;
        b = w - pre + 26;
      }
      res += min(a, b) + 1;
      pre = w;
    }
    return res;
  }

  vector<vector<int>> flipAndInvertImage(vector<vector<int>> &image) {
    int n = image.size();
    int nn = image[0].size();
    for (int i = 0; i < n; i++) {
      reverse(image[i].begin(), image[i].end());
      for (int j = 0; j < nn; j++) {
        reverse10(image[i][j]);
      }
    }
    return image;
  }
  void reverse10(int &a) {
    if (a == 0) {
      a = 1;
    } else {
      a = 0;
    }
  }

  //   static bool judge(const vector<int> a, const vector<int> b) {
  //     return a[1] < b[1];
  //   }
  //     bool needMoreThanZero(const vector<int> &a) {
  //       if (a[0] == a[1])
  //         return false;
  //       if (a[0] <= a[1] - 1)
  //         return true;
  //       else
  //         return false;
  //     }
  //     double maxAverageRatio(vector<vector<int>> &classes, int extraStudents)
  //     {
  //       double res = 0;
  //       sort(classes.begin(), classes.end(), judge);
  //       for (int i = 0; i < classes.size(); i++) {
  //         if (extraStudents > 0) {
  //           if (needMoreThanZero(classes[i])) {
  //             classes[i][0]++;
  //             classes[i][1]++;
  //             extraStudents--;
  //           }
  //         }
  //         res += (double)classes[i][0] / classes[i][1];
  //       }
  //       return res / classes.size();
  //     }

  double maxAverageRatio(vector<vector<int>> &classes, int extraStudents) {
    priority_queue<tuple<double, int, int>> q;

    auto diff = [](int x, int y) -> double {
      return (double)(x + 1) / (y + 1) - (double)x / y;
    };

    double ans = 0.;
    for (const auto &c : classes) {
      int x = c[0], y = c[1];
      ans += (double)x / y;
      q.emplace(diff(x, y), x, y);
    }
    for (int _ = 0; _ < extraStudents; ++_) {
      // auto {d, x, y} = q.top();
      tuple<double, int, int> tmp = q.top();
      auto d = get<0>(tmp);
      auto x = get<1>(tmp);
      auto y = get<2>(tmp);
      q.pop();
      ans += d;
      q.emplace(diff(x + 1, y + 1), x + 1, y + 1);
    }
    return ans / classes.size();
  }

  bool isRectangleCover(vector<vector<int>> &rectangles) {
    int left = INT_MAX;
    int right = INT_MIN;
    int top = INT_MIN;
    int bottom = INT_MAX;

    int n = rectangles.size();
    int sumArea = 0;
    set<pair<int, int>> seet;
    for (int i = 0; i < n; i++) {
      left = min(left, rectangles[i][0]);
      bottom = min(bottom, rectangles[i][1]);
      right = max(right, rectangles[i][2]);
      top = max(top, rectangles[i][3]);

      //计算最小矩形总面积
      sumArea += (rectangles[i][3] - rectangles[i][1]) *
                 (rectangles[i][2] - rectangles[i][0]);

      if (seet.find(pair<int, int>(rectangles[i][0], rectangles[i][3])) !=
          seet.end())
        seet.erase(pair<int, int>(rectangles[i][0], rectangles[i][3]));
      else
        seet.insert(pair<int, int>(rectangles[i][0], rectangles[i][3]));
      if (seet.find(pair<int, int>(rectangles[i][0], rectangles[i][1])) !=
          seet.end())
        seet.erase(pair<int, int>(rectangles[i][0], rectangles[i][1]));
      else
        seet.insert(pair<int, int>(rectangles[i][0], rectangles[i][1]));
      if (seet.find(pair<int, int>(rectangles[i][2], rectangles[i][3])) !=
          seet.end())
        seet.erase(pair<int, int>(rectangles[i][2], rectangles[i][3]));
      else
        seet.insert(pair<int, int>(rectangles[i][2], rectangles[i][3]));
      if (seet.find(pair<int, int>(rectangles[i][2], rectangles[i][1])) !=
          seet.end())
        seet.erase(pair<int, int>(rectangles[i][2], rectangles[i][1]));
      else
        seet.insert(pair<int, int>(rectangles[i][2], rectangles[i][1]));
    }
    if (seet.size() == 4 &&
        seet.find(pair<int, int>(left, top)) != seet.end() &&
        seet.find(pair<int, int>(left, bottom)) != seet.end() &&
        seet.find(pair<int, int>(right, bottom)) != seet.end() &&
        seet.find(pair<int, int>(right, top)) != seet.end()) {
      return sumArea == (right - left) * (top - bottom);
    } else {
      return false;
    }
  }

  static bool judgeString(char s, char ss) { return (int)s > (int)ss; }

  int nextGreaterElement(int n) {
    long ress;
    string s = to_string(n);
    sort(s.begin(), s.end(), judgeString);
    ress = std::stoi(s);
    if (ress > INT_MAX || ress < INT_MIN) {
      return -1;
    }
    return (int)ress;
  }

  bool checkPowersOfThree(int n) {
    while (n % 3 == 1 || n % 3 == 0) {
      n = n / 3;
      if (n == 1 || n == 0) {
        return 1;
      }
    }
    return 0;
  }

  bool isPalindrome(string s) {
    int n = s.size();
    if (n == 0) {
      return 1;
    }
  }

  int arrayPairSum(vector<int> &nums) {
    int res = 0;
    sort(nums.begin(), nums.end());
    for (int i = 0; i < nums.size(); i += 2)
      res += nums[i];
    return res;
  }

  int maximumProduct(vector<int> &nums) {
    sort(nums.begin(), nums.end(),
         [](int a, int b) -> bool { return a > b; }); // Lambda表达式
    return max(nums[0] * nums[1] * nums[2],
               nums[nums.size() - 1] * nums[nums.size() - 2] * nums[0]);
  }

  vector<int> findErrorNums(vector<int> &nums) {
    vector<int> res(2);
    sort(nums.begin(), nums.end());
    unordered_map<int, int> mapp;
    bool flag = 0;
    for (int i = 0; i < nums.size(); i++) {
      mapp[nums[i]]++;
    }
    for (int i = 0; i < nums.size(); i++) {
      if (mapp[i + 1] == 0) {
        res[1] = i + 1;
      }
      if (mapp[i + 1] == 2) {
        res[0] = i + 1;
      }
    }
    return res;
  }

  string addBinary(string a, string b) {
    string res;
    int n = min(a.size(), b.size());
    bool flag = false;
    while (n) {
      if (a[n - 1] == '1' && b[n - 1] == '1') {
        res += to_string(0 + flag);
        flag = 1;
      } else {
        res += to_string(a[n - 1] + b[n - 1] + flag);
        flag = 0;
      }
      n--;
    }
    int ss = a.size() - n;
    while (ss) {
      if (a[ss - 1] == '1' && flag) {
        res += to_string(0);
        flag = 1;
      } else {
        res += to_string(a[n - 1] + flag);
        flag = 0;
      }
      ss--;
    }
    reverse(res.begin(), res.end());
    return res;
  }

  //************************************
  // Method:    longestPalindrome2
  // FullName:  Solution::longestPalindrome2
  // Access:    public
  // Returns:   int
  // Qualifier:
  // Parameter: string s
  //************************************
  //************************************
  // 函数名: longestPalindrome2
  // 完整名:
  // 权 限:  public
  // 返回值: int
  // 参 数:  string s
  //************************************
  int longestPalindrome2(string s) {
    unordered_map<char, int> mapp;
    for (char ch : s)
      mapp[ch]++;
    int max_odd = 0;
    int sum_even = 0;
    for (auto it = mapp.begin(); it != mapp.end(); it++) {
      if (!(it->second % 2))
        sum_even += it->second;
      else
        max_odd = max(max_odd, it->second);
    }
    return sum_even + max_odd;
  }

  //************************************
  // 函数名: findNthDigit1
  // 完整名: Solution::findNthDigit1
  // 权 限:  public
  // 返回值: int
  // 参 数:  int n
  // 备 注:
  //************************************
  int findNthDigit1(int n) {
    int i = 1;
    bool a;
    if (n >= 10)
      a = 1;
    while (n) {
      string s = to_string(i);
      cout << s << endl;
      n -= s.size();
      cout << "n = " << n << endl;
      if (n <= 0)
        return s[s.size() - 1 + n] - 48;

      i++;
    }
    return 0;
  }
  //************************************
  // 函数名: findNthDigit
  // 完整名: Solution::findNthDigit
  // 权 限:  public
  // 返回值: int
  // 参 数:  int n
  //************************************
  int findNthDigit2(int n) {
    int i = 1;
    int pren = n;
    vector<int> res;
    while (n > 0) {
      pren = n; // 1 10 100之后的位数
      int temp = pow(10, i - 1) * 9 * i;
      n = n - pow(10, i - 1) * 9 * i;
      i++;
    }
    int nown = pow(10, i - 2);
    for (int j = 0; j < pren; j++, nown++) {
      string s = to_string(nown);
      for (int jj = 0; jj < s.size(); jj++) {
        res.push_back(s[jj] - '0');
        if (res.size() == pren) {
          cout << res.back();
          return res.back();
        }
      }
    }
    return 0;
  }

  int findNthDigit(int n) {
    int d = 1, count = 9;
    while (n > (long)d * count) {
      n -= d * count;
      d++;
      count *= 10;
    }
    int index = n - 1;
    int start = (int)pow(10, d - 1);
    int num = start + index / d;
    int digitIndex = index % d;
    int digit = (num / (int)(pow(10, d - digitIndex - 1))) % 10;
    return digit;
  }

  int largestPerimeter(vector<int> &nums) {
    sort(nums.begin(), nums.end(), [](int a, int b) -> bool { return a > b; });
    for (int i = 0; i < nums.size(); i++) {
      for (int j = i + 1; j < nums.size(); j++) {
        for (int k = j + 1; k < nums.size(); k++) {
          if (nums[j] + nums[k] > nums[i]) {
            return nums[i] + nums[j] + nums[k];
          }
        }
      }
    }
    return 0;
  }

  int findMaxConsecutiveOnes(vector<int> &nums) {
    int max_res = INT_MIN;
    int pre = nums[0] == 1 ? 1 : 0;
    int temp = nums[0] == 1 ? 1 : 0;

    max_res = max(max_res, temp);
    for (int i = 1; i < nums.size(); i++) {
      if (nums[i] == 1) {
        temp++;
        max_res = max(max_res, temp);
      } else {
        pre = 0;
        temp = 0;
      }
    }
    return max_res;
  }
  int findMaxConsecutiveOnes2(vector<int> &nums) {
    int max_res = 0;
    int temp = 0;
    for (int i = 0; i < nums.size(); i++) {
      if (nums[i] == 1) {
        temp++;
        max_res = max(max_res, temp);
      } else {
        temp = 0;
      }
    }
    return max_res;
  }

  int findShortestSubArray(vector<int> &nums) {
    // 定义统计哈希表和最大频数
    unordered_map<int, int> mp;
    int max_fre = 0;
    // 统计频数，并计算最大频数
    for (int i = 0; i < nums.size(); i++) {
      mp[nums[i]]++;
      max_fre = max(max_fre, mp[nums[i]]);
    }
    // 重构map
    mp.erase(mp.begin(), mp.end());
    // 定义窗口和满足频数的最短长度
    int ans = nums.size();
    int left = 0, right = 0;
    while (right < nums.size()) {
      mp[nums[right]]++;
      // 频数达到要求后，移动左边界
      while (mp[nums[right]] == max_fre) {
        ans = min(ans, right - left + 1);
        mp[nums[left++]]--;
      }
      right++;
    }
    return ans;
  }

  int countSegments(string s) {
    int res = 0;
    for (int i = 0; i < s.size(); i++) {
      if ((i == 0 || s[i - 1] == ' ') && s[i] != ' ') {
        res++;
      }
    }
    return res;
  }

  //************************************
  // 函数名: largestSumAfterKNegations
  // 完整名: Solution::largestSumAfterKNegations
  // 权 限:  public
  // 返回值: int
  // 参 数:  vector<int> & nums
  // 参 数:  int k
  // 备 注:
  //************************************
  int largestSumAfterKNegations(vector<int> &nums, int k) {
    sort(nums.begin(), nums.end());
    int sum = 0;
    int n = nums.size();
    for (int i = 0; i < n; i++) {
      if (nums[i] < 0 && k > 0) {
        nums[i] = -1 * nums[i];
        k--;
      }
      sum += nums[i];
    }
    sort(nums.begin(), nums.end());
    return sum - (k % 2 == 0 ? 0 : 2 * nums[0]);
  }

  //************************************
  // 函数名: makeGood
  // 完整名: Solution::makeGood
  // 权 限:  public
  // 返回值: std::string
  // 参 数:  string s
  // 备 注:
  //************************************
  string makeGood(string s) {
    for (int i = 0; i < s.size() - 1; i++) {
      if (s.size() == 0) {
        return "";
      }
      if (islower(s[i])) {
        if (isupper(s[i + 1])) {
          if (s[i] == tolower(s[i + 1])) {
            s.erase(i, 2);
            i = -1;
          }
        }
      } else {
        if (islower(s[i + 1])) {
          if (s[i] == toupper(s[i + 1])) {
            s.erase(i, 2);
            i = -1;
          }
        }
      }
    }
    return s;
  }

  //************************************
  // 函数名: countStudents
  // 完整名: Solution::countStudents
  // 权 限:  public
  // 返回值: int
  // 参 数:  vector<int> & students
  // 参 数:  vector<int> & sandwiches
  // 备 注:  1700. 无法吃午餐的学生数量
  //************************************
  int countStudents(vector<int> &students, vector<int> &sandwiches) {
    queue<int> stu;
    queue<int> sand;
    for (int i = 0; i < students.size(); i++) {
      stu.push(students[i]);
      sand.push(sandwiches[i]);
    }
    while (!noStuSand(stu, sand)) {
      if (stu.front() == sand.front()) {
        stu.pop();
        sand.pop();
      } else {
        int temp = stu.front();
        stu.pop();
        stu.push(temp);
      }
    }
    return stu.size();
  }
  bool noStuSand(queue<int> stu, queue<int> sand) {
    while (!stu.empty()) {
      if (stu.front() != sand.front()) {
        stu.pop();
      } else {
        return 0;
      }
    }
    return true;
  }

  void inorder(TreeNode *node, vector<int> &res) {
    if (node == nullptr) {
      return;
    }
    inorder(node->left, res);
    res.push_back(node->val);
    inorder(node->right, res);
  }

  TreeNode *increasingBST(TreeNode *root) {
    vector<int> res;
    inorder(root, res);

    TreeNode *dummyNode = new TreeNode(-1);
    TreeNode *currNode = dummyNode;
    for (int value : res) {
      currNode->right = new TreeNode(value);
      currNode = currNode->right;
    }
    return dummyNode->right;
  }

  bool canConstruct2(string ransomNote, string magazine) {
    vector<int> cnt(26);
    for (int i = 0; i < magazine.size(); i++) {
      cnt[magazine[i] - 'a']++;
    }
    for (int i = 0; i < ransomNote.size(); i++) {
      cnt[ransomNote[i] - 'a']--;
      if (cnt[ransomNote[i] - 'a'] < 0) {
        return false;
      }
    }
    return 1;
  }

  double myPow(double x, int n) {
    double res = 1.0;
    for (int i = n; i != 0; i /= 2) {
      if (i % 2 != 0) {
        res *= x;
      }
      x *= x;
    }
    return n < 0 ? 1 / res : res;
  }

  int superPow(int a, vector<int> &b) {
    long bb = 0;
    for (int i = 0; i < b.size(); i++) {
      bb += b[i] * pow(10, b.size() - 1 - i);
    }
    return mypow(a, bb) % 1337;
  }
  long mypow(int a, long b) {
    if (b == 0)
      return 1;
    long tmp_res = mypow(a, b / 2);
    return b % 2 ? tmp_res * tmp_res * a : tmp_res * tmp_res;
  }

  int findKthPositive(vector<int> &arr, int k) {
    int res;
    int count = 0;
    int temp = 0;
    for (int i = 1;; i++) {
      if (arr[temp] != i) {
        count++;
        if (count == 5) {
          return i;
        }
        continue;
      } else {
        temp++;
      }
    }
    return 0;
  }

  int numIslands(vector<vector<char>> &grid) {
    int nr = grid.size();
    if (!nr)
      return 0;
    int nc = grid[0].size();

    int num_islands = 0;
    for (int r = 0; r < nr; ++r) {
      for (int c = 0; c < nc; ++c) {
        if (grid[r][c] == '1') {
          ++num_islands;
          grid[r][c] = '0';
          queue<pair<int, int>> neighbors;
          neighbors.push({r, c});
          while (!neighbors.empty()) {
            auto rc = neighbors.front();
            neighbors.pop();
            int row = rc.first, col = rc.second;
            if (row - 1 >= 0 && grid[row - 1][col] == '1') {
              neighbors.push({row - 1, col});
              grid[row - 1][col] = '0';
            }
            if (row + 1 < nr && grid[row + 1][col] == '1') {
              neighbors.push({row + 1, col});
              grid[row + 1][col] = '0';
            }
            if (col - 1 >= 0 && grid[row][col - 1] == '1') {
              neighbors.push({row, col - 1});
              grid[row][col - 1] = '0';
            }
            if (col + 1 < nc && grid[row][col + 1] == '1') {
              neighbors.push({row, col + 1});
              grid[row][col + 1] = '0';
            }
          }
        }
      }
    }

    return num_islands;
  }

  int islandPerimeter(vector<vector<int>> &grid) {
    int res = 0;
    for (int i = 0; i < grid.size(); i++) {
      for (int j = 0; j < grid[0].size(); j++) {
        if (grid[i][j] == 1) {
          res += dfs(grid, i, j);
        }
      }
    }
    return res;
  }

  int dfs(vector<vector<int>> grid, int i, int j) {
    int res = 4;
    if (i - 1 >= 0 && grid[i - 1][j] == 1)
      res--;
    if (j - 1 >= 0 && grid[i][j - 1] == 1)
      res--;
    if (i + 1 < grid.size() && grid[i + 1][j] == 1)
      res--;
    if (j + 1 < grid[0].size() && grid[i][j + 1] == 1)
      res--;
    return res;
  }

  bool validTicTacToe(vector<string> &board) {
    int x_count = 0;
    int o_count = 0;
    bool flag1 = 0;
    bool flag2 = 0;

    for (int i = 0; i < board.size(); i++) {
      unordered_map<int, int> mapp;
      unordered_map<int, int> mappp;
      for (int j = 0; j < board[i].size(); j++) {
        if (board[i][j] == 'X') {
          x_count++;
        } else if (board[i][j] == 'O') {
          o_count++;
        }
        mapp[board[i][j]]++;
        mappp[board[j][i]]++;
      }
      if (mapp.size() == 1 && board[i][0] != ' ') {
        if (board[i][0] == 'X') {
          flag1 = 1;
        } else {
          flag2 = 1;
        }
      }
      if (mappp.size() == 1 && board[0][i] != ' ') {
        if (board[0][i] == 'X') {
          flag1 = 1;
        } else {
          flag2 = 1;
        }
      }
    }

    if (board[0][0] == board[1][1] && board[2][2] == board[1][1]) {
      if (board[0][0] == 'X') {
        flag1 = 1;
      } else if (board[0][0] == 'O') {
        flag2 = 1;
      }
    }

    if (board[0][2] == board[1][1] && board[2][0] == board[1][1]) {
      if (board[0][2] == 'X') {
        flag1 = 1;
      } else if (board[0][2] == 'O') {
        flag2 = 1;
      }
    }

    if (flag1 == 1 && flag2 == 1) {
      return false;
    } else if (flag1 == 1) {
      return x_count == o_count + 1;
    } else if (flag2 == 1) {
      return x_count == o_count;
    } else {
      return x_count - o_count <= 1 && x_count >= o_count;
    }
  }

  int checkRecord(int n) {
    int all_situ = pow(n, 3);
    int aa = allAA(n);
    int lll = allLLL(n);
    int s2 = subTwice(n);
    return all_situ - allAA(n) - allLLL(n) + subTwice(n);
  }
  int allAA(int n) { return pow(n - 2, 3) * (n - 1); }
  int allLLL(int n) { return pow(n - 3, 3) * (n - 2); }
  int subTwice(int n) {
    int cc;
    cc = C(n - 5, 2);
    return pow(n - 5, 3) * cc;
  }

  long long C(int N, int M) {
    long long sum = 1;
    for (int i = 1; i <= M; i++) {
      sum = sum * (N - M + i) / i;
    }
    return sum;
  }

  vector<int> getAverages(vector<int> &nums, int k) {
    if (k == 0)
      return nums;
    if (2 * k + 1 > nums.size())
      return vector<int>(nums.size(), -1);
    vector<int> res;
    int win_sum = 0;
    for (int i = 0; i < k; i++) {
      res.push_back(-1);
    }
    for (int i = k; i < nums.size() - k; i++) {
      for (int j = i - k; j < i + k + 1; j++) {
        win_sum += nums[j];
      }
      res.push_back(win_sum / (2 * k + 1));
      win_sum = 0;
    }
    for (int i = nums.size() - k; i < nums.size(); i++) {
      res.push_back(-1);
    }
    return res;
  }

  int uniqueLetterString(string s) {
    vector<string> all = getAllSubstrings(s);
    int res = 0;

    for (int i = 0; i < all.size(); i++) {
      res += countUniqueChars(all[i]);
    }
    return res % 1000000007;
  }

  int countUniqueChars(string s) {
    int res = 0;
    unordered_map<int, int> mapp;
    for (int i = 0; i < s.size(); i++) {
      mapp[s[i]]++;
    }
    for (auto it = mapp.begin(); it != mapp.end(); it++) {
      if (it->second == 1)
        res++;
    }
    return res;
  }

  vector<string> getAllSubstrings(string str) {
    vector<string> res;
    if (str.length() == 0) {
      return res;
    } else {
      for (int i = 0; i < str.length(); i++) {
        for (int j = 1; j < str.length() - i + 1; j++) {
          std::cout << str.substr(i, j) << std::endl;
          res.push_back(str.substr(i, j));
        }
      }
    }
    return res;
  }

  int numJewelsInStones(string j, string s) {
    unordered_set<char> us(begin(j), end(j));
    return count_if(begin(s), end(s), [&](char c) { return us.count(c); });
  }

  string longestWord1(vector<string> &words) {
    string ans = "";
    unordered_set<string> se(begin(words), end(words));
    // for(string s:words)
    //     se.insert(s);
    for (string s : words) {
      if (s.size() > ans.size() ||
          (s.size() == ans.size() && s.compare(ans) < 0)) {
        bool f = true;
        for (int k = 1; k <= s.size(); k++) //判断前缀
        {
          if (!se.count(s.substr(0, k))) {
            f = false;
            break;
          }
        }
        if (f)
          ans = s;
      }
    }
    return ans;
  }

  string longestWord(vector<string> &words) {
    sort(words.begin(), words.end()); //从小到大排序
    unordered_set<string> mySet;
    string ans = ""; //从空开始记录
    mySet.insert(ans);
    for (auto &word : words) {
      string tmp = string(word.begin(), word.end() - 1);
      if (mySet.find(string(word.begin(), word.end() - 1)) != mySet.end()) {
        if (word.size() > ans.size())
          ans = word;       //如果长度更长，更新答案
        mySet.insert(word); //记录这个单词
      }
    }
    return ans;
  }

  string licenseKeyFormatting(string s, int k) {
    string res;
    int count = 0;
    for (int i = s.size() - 1; i >= 0; i--) {
      if (s[i] != '-') {
        res += toupper(s[i]);
        count++;
      }
      if (count == k) {
        res += '-';
        count = 0;
      }
    }

    reverse(res.begin(), res.end());
    if (res[0] == '-') {
      res.erase(0, 1);
    }
    return res;
  }

  bool isSubsequence(string s, string t) {
    unordered_set<int> sett(t.begin(), t.end());
    for (char ch : s) {
      if (sett.find(ch) == sett.end()) {
        return 0;
      }
    }
    return 1;
  }

  int visiblePoints(vector<vector<int>> &points, int angle,
                    vector<int> &location) {
    int sameCnt = 0;
    vector<double> polarDegrees;
    for (auto &point : points) {
      if (point[0] == location[0] && point[1] == location[1]) {
        sameCnt++;
        continue;
      }
      double degree = atan2(point[1] - location[1], point[0] - location[0]);
      polarDegrees.emplace_back(degree);
    }
    sort(polarDegrees.begin(), polarDegrees.end());

    int m = polarDegrees.size();
    for (int i = 0; i < m; ++i) {
      polarDegrees.emplace_back(polarDegrees[i] + 2 * M_PI);
    }

    int maxCnt = 0;
    double degree = angle * M_PI / 180;
    for (int i = 0; i < m; ++i) {
      auto it = upper_bound(polarDegrees.begin() + i, polarDegrees.end(),
                            polarDegrees[i] + degree);
      int curr = it - polarDegrees.begin() - i;
      maxCnt = max(maxCnt, curr);
    }
    return maxCnt + sameCnt;
  }

  //超时！
  int findRadius(vector<int> &houses, vector<int> &heaters) {
    int houses_size = houses.size();
    int heaters_size = heaters.size();
    int length = houses[houses_size - 1] - houses[0] + 1;
    int beginR = (length / heaters_size - 1) / 2;
    bool is_covered = false;
    while (!is_covered) {
      vector<int> temp = houses;
      for (int i = 0; i < heaters_size; i++) {
        heatR(temp, heaters[i], beginR);
      }
      if (temp == vector<int>(temp.size(), 0)) {
        is_covered = 1;
      } else
        beginR++;
    }
    return beginR;
  }

  void heatR(vector<int> &houses, int heatpos, int r) {
    for (int i = 0; i < houses.size(); i++) {
      if (abs(houses[i] - heatpos) <= r) {
        houses[i] = 0;
      }
    }
  }

  int dayOfYear(string date) {
    int ans = 0;
    bool isRun;
    int year = (date[0] - '0') * 1000 + (date[1] - '0') * 100 +
               (date[2] - '0') * 10 + (date[3] - '0') * 1;
    int month = (date[5] - '0') * 10 + (date[6] - '0') * 1;
    int day = (date[8] - '0') * 10 + (date[9] - '0') * 1;

    if (year % 4 == 0) {
      if (year % 100 == 0) {
        // // 这里如果被 400 整除是闰年
        if (year % 400 == 0)
          isRun = 1;
        else
          isRun = 0;
      } else
        isRun = 1;
    } else
      isRun = 0;
    int lastday;
    for (int i = 1; i <= month; i++) {
      if (i == 1 || i == 3 || i == 5 || i == 7 || i == 8 || i == 10 ||
          i == 12) {
        ans += 31;
        lastday = 31;
      } else if (i == 2) {
        if (isRun) {
          ans += 29;
          lastday = 29;
        } else {
          ans += 28;
          lastday = 28;
        }
      } else {
        ans += 30;
        lastday = 30;
      }
    }
    ans -= (lastday - day);
    return ans;
  }

  typedef pair<int, int> pii;
  int eatenApples(vector<int> &apples, vector<int> &days) {
    int ans = 0;
    priority_queue<pii, vector<pii>, greater<pii>> pq;
    int n = apples.size();
    int i = 0;
    while (i < n) {
      while (!pq.empty() && pq.top().first <= i) {
        pq.pop();
      }
      int rottenDay = i + days[i]; //腐烂天数
      int count = apples[i];
      if (count > 0) {
        pq.emplace(rottenDay, count);
      }
      if (!pq.empty()) {
        pii curr = pq.top();
        pq.pop();
        curr.second--;
        if (curr.second != 0) {
          pq.emplace(curr.first, curr.second);
        }
        ans++;
      }
      i++;
    }
    while (!pq.empty()) {
      while (!pq.empty() && pq.top().first <= i) {
        pq.pop();
      }
      if (pq.empty()) {
        break;
      }
      pii curr = pq.top();
      pq.pop();
      int num = min(curr.first - i, curr.second);
      ans += num;
      i += num;
    }
    return ans;
  }

  bool isNStraightHand(vector<int> &hand, int groupSize) {
    if (hand.size() % groupSize)
      return 0;
    sort(hand.begin(), hand.end());
    unordered_map<int, int> mapp;
    for (int i : hand)
      mapp[i]++;
    bool flag = 0;
    int packsize = 0;
    int pack = hand.size() / groupSize;
    int begin = hand[0];
    int maxx = hand[hand.size() - 1];
    while (pack > packsize) {
      flag = 0;
      int temp = begin;
      for (int j = temp; j < temp + groupSize; j++) {
        mapp[j]--;
        if (mapp[j] < 0)
          return 0;
        else if (mapp[j] > 0 && !flag) {
          begin = j;
          flag = 1;
        } else if (j == temp + groupSize - 1 && !flag) {
          if (j < maxx)
            begin = *upper_bound(hand.begin(), hand.end(), j);
        }
      }
      packsize++;
    }
    return 1;
  }

  int lastRemaining(int n) {
    vector<int> all;
    for (int i = 1; i <= n; i++) {
      all.push_back(i);
    }
    bool inverse = 0;
    while (all.size() != 1) {
      vector<int> temp;
      if (inverse) {
        for (int i = all.size() - 2; i >= 0; i -= 2) {
          temp.push_back(all[i]);
        }
        reverse(temp.begin(), temp.end());
      } else {
        for (int i = 1; i < all.size(); i += 2) {
          temp.push_back(all[i]);
        }
      }

      inverse = inverse == 0 ? 1 : 0;
      all = temp;
    }
    return all[0];
  }

  string modifyString(string s) {
    char pre = '0';
    char next = '0';
    for (int i = 0; i < s.size(); i++) {
      if (s[i] == '?') {
        pre = i > 0 ? s[i - 1] : '0';
        next = i < s.size() - 1 ? s[i + 1] : '0';
        char replace = 'a';
        while (replace == pre || replace == next) {
          replace++;
        }
        s[i] = replace;
      }
    }
    return s;
  }

  //   char slowestKey(vector<int> &releaseTimes, string keysPressed) {
  //     vector<int> count(26);
  //     for (int i = 0; i < releaseTimes.size(); i++) {
  //       if (i == 0) {
  //         count[keysPressed[i] - 'a'] += releaseTimes[i];
  //       } else {
  //         count[keysPressed[i] - 'a'] += releaseTimes[i] - releaseTimes[i -
  //         1];
  //       }
  //     }
  //     int maxx = *std::max_element(count.begin(), count.end());
  //     for (int i = 0; i < keysPressed.size(); i++)
  //       if (count[keysPressed[i] - 'a'] == maxx)
  //         return keysPressed[i];
  //     return 'a';
  //   }

  char slowestKey(vector<int> &releaseTimes, string keysPressed) {
    int maxdur = 0;
    vector<int> cnt(26);
    char res;
    for (int i = 0; i < releaseTimes.size(); i++) {
      if (i == 0) {
        if (maxdur < releaseTimes[i]) {
          maxdur = releaseTimes[i];
          cnt[keysPressed[i] - 'a'] = maxdur;
        }
      } else {
        if (maxdur < releaseTimes[i] - releaseTimes[i - 1]) {
          maxdur = releaseTimes[i] - releaseTimes[i - 1];
          cnt[keysPressed[i] - 'a'] = maxdur;
        }
      }
    }
    for (int i = 0; i < releaseTimes.size(); i++) {
      if (cnt[keysPressed[i] - 'a'] == maxdur) {
        res = keysPressed[i];
        break;
      }
    }
    return res;
  }

  //   bool isAdditiveNumber(string num) {
  //     int len = num.size();
  //     int numcount = 3;
  //     while (numcount < len) {
  // 		int num1 = num.substr()
  //     }
  //     for (int i = 0; i < len; i++) {
  //     }
  //   }

  vector<vector<int>> res;
  vector<int> path;

  vector<vector<int>> combine(int n, int k) {
    res.clear();
    path.clear();
    backTrackCombine(n, k, 1);
    return res;
  }
  void backTrackCombine(int n, int k, int startIndex) {
    if (path.size() == k) {
      res.push_back(path);
      return;
    }
    for (int i = startIndex; i <= n - (k - path.size()) + 1; i++) {
      path.push_back(i);
      backTrackCombine(n, k, i + 1);
      path.pop_back();
    }
  }

  bool increasingTriplet3(vector<int> &nums) {
    int size = nums.size();
    vector<pair<int, int>> all; //
    for (int i = 0; i < size; i++)
      all.push_back(make_pair(nums[i], i));
    sort(all.begin(), all.end(),
         [](pair<int, int> a, pair<int, int> b) -> bool {
           return a.first < b.first;
         });
    for (int i = size - 1; i >= 3; i--) {
      if (all[i].second >= 3) {
        return 1;
      }
    }
    return 0;
  }

  bool increasingTriplet(vector<int> &nums) {
    int len = nums.size();
    if (len < 3)
      return false;
    int small = INT_MAX, mid = INT_MAX;
    for (auto num : nums) {
      if (num <= small) {
        small = num;
      } else if (num <= mid) {
        mid = num;
      } else if (num > mid) {
        return true;
      }
    }
    return false;
  }

  bool increasingTriplet2(vector<int> &nums) {
    int size = nums.size();
    for (int i = 0; i < size; i++) {
      for (int j = i + 1; j < size; j++) {
        if (nums[i] >= nums[j]) {
          break;
        }
        for (int jj = j + 1; jj < size; jj++) {
          if (nums[j] < nums[jj]) {
            return 1;
          }
        }
      }
    }
    return 0;
  }
};

void solve() {
  Solution mysolution;
  //   int aa = mysolution.findRadius(
  //       vector<int>{282475249, 622650073, 984943658, 144108930, 470211272,
  //                   101027544, 457850878, 458777923},
  //       vector<int>{823564440, 115438165, 784484492, 74243042, 114807987,
  //                   137522503, 441282327, 16531729, 823378840, 143542612});

  vector<int> a{3, 0, 0, 0, 0, 2};
  vector<int> b{3, 0, 0, 0, 0, 2};
  bool aa = mysolution.increasingTriplet(vector<int>{0, 4, 2, 1, 0, -1, -3});
  cout << "pause!" << endl;
}

int main() {
  solve();
  vector<double> test(20, 5.0);
  cout << "program has finished..." << endl;

  //   string str1 = "hi,test,hello";
  //   string str2 = "hi,test";
  //   //字符串比较
  //   if (str1.compare(str2) > 0)
  //     printf("str1>str2\n");
  //   else if (str1.compare(str2) < 0)
  //     printf("str1<str2\n");
  //   else
  //     printf("str1==str2\n");
  //
  //   // str1的子串（从索引3开始，包含4个字符）与str2进行比较
  //   if (str1.compare(3, 4, str2) == 0)
  //     printf("str1的指定子串等于str2\n");
  //   else
  //     printf("str1的指定子串不等于str2\n");
  //
  //   // str1指定子串与str2的指定子串进行比较
  //   if (str1.compare(3, 4, str2, 3, 4) == 0)
  //     printf("str1的指定子串等于str2的指定子串\n");
  //   else
  //     printf("str1的指定子串不等于str2的指定子串\n");
  //
  //   // str1指定子串与字符串的前n个字	符进行比较
  //   if (str1.compare(0, 2, "hi,hello", 2) == 0)
  //     printf("str1的指定子串等于指定字符串的前2个字符组成的子串\n");
  //   else
  //     printf("str1的指定子串不等于指定字符串的前2个字符组成的子串\n");
  //   return 0;

  //   //八股取士
  //   // 定义简单的lambda表达式
  //   auto basicLambda = [] { cout << "Hello, world!" << endl; };
  //   basicLambda(); // 输出：Hello, world!
  //
  //   // 指明返回类型，托尾返回类型
  //   auto add = [](int a, int b) -> int { return a + b; };
  //   // ⾃动推断返回类型
  //   auto multiply = [](int a, int b) { return a * b; };
  //   int sum = add(2, 5);          // 输出：7
  //   int product = multiply(2, 5); // 输出：10
  //
  //   int x = 10;
  //   auto add_x = [x](int a) {
  //     return a + x;
  //   }; // 复制捕捉x,lambda  表达式⽆法修改此变ᰁ
  //   auto multiply_x = [&x](int a) {
  //     return a * x;
  //   }; // 引⽤捕捉x，lambda  表达式可以修改此变ᰁ
  //   cout << add_x(10) << " " << multiply_x(10) << endl;
  //   // 输出：20 100

  //   int val = 3;
  //   vector<int> v{1, 8, 5, 3, 6, 10};
  //   int count =
  //       std::count_if(v.begin(), v.end(), [val](int x) { return x > val;
  //       });
  //   // v中⼤于3的元素数ᰁ

  //// C++11的泛化常数
  // constexpr int N = 5; // N 变成了⼀个只读的值
  // int arr[N];          // OK

  //   //初始化列表
  //   class A {
  //   public:
  //     A(std::initializer_list<int> list);
  //   };
  //   A a = {1, 2, 3};

  //   vector<int> v{2, 3, 4};
  //   for (auto &x : v)
  //     x = 1;

  //   struct Base1 final {};
  //   struct Derived1 : Base1 {}; // 编译错：Base1不允许被继承
  //   struct Base2 {
  //     virtual void f1() final;
  //     virtual void f2();
  //   };
  //   struct Derived2 : Base2 {
  //     virtual void f1();             // 编译错：f1不允许᯿写
  //     virtual void f2(int) override; // 编译错：⽗类中没有 void f2(int)
  //   };

  //   //   classA 中定义了 classA(T value)
  //   //   构造函数，因此编译器不会默认⽣成⼀个⽆参数
  //   // 	  的构造函数了， 如果我们需要可以⼿动声明，或者直接 = default
  //   struct classA {
  //     classA() = default; // 声明⼀个⾃动⽣成的函数
  //     classA(int value);
  //     void *operator new(size_t) = delete; // 禁⽌⽣成new运算符
  //   };

  //   std::shared_ptr<double> p_first(new double);
  //   {
  //     std::shared_ptr<double> p_copy = p_first;
  //     *p_copy = 21.2;
  //   } // p_copy 被销毁，⾥⾯的 double 还有⼀个引⽤因此仍然保持

  //   typedef std::tuple<int, double, string> tuple_1;
  //   tuple_1 t1;
  //   typedef std::tuple<char, short, const char *> tuple_2;
  //   tuple_2 t2('X', 2, "Hola!");
  //   t1 = t2; // 隐式类型转换

  //////////////////////////////////////////////////////////////////////////
  //数据结构和算法
  //////////////////////////////////////////////////////////////////////////

  //   auto BubbleSort = [](std::vector<int> &nums, int n) {
  //     if (n <= 1)
  //       return;
  //     bool is_swap;
  //     for (int i = 1; i < n; ++i) {
  //       is_swap = false;
  //       //设定⼀个标记，若为false，则表示此次循环没有进⾏交换，也就是待排序列已经有序，排序已经完成。
  //       for (int j = 1; j < n - i + 1; ++j) {
  //         if (nums[j] < nums[j - 1]) {
  //           std::swap(nums[j], nums[j - 1]);
  //           is_swap = true; //表示有数据交换
  //         }
  //       }
  //       if (!is_swap)
  //         break; //没有数据交集，提前退出
  //     }
  //   };
  //   vector<int> a{34, 66, 2, 5, 95, 4, 46, 27};
  //   BubbleSort(a, a.size()); // cout => 2 4 5 27 34 46 66 95

  auto InsertSort = [](std::vector<int> &nums, int n) {
    if (n <= 1)
      return;
    for (int i = 0; i < n; ++i) {
      for (int j = i; j > 0 && nums[j] < nums[j - 1]; --j) {
        std::swap(nums[j], nums[j - 1]);
      }
    }
  };
  std::vector<int> nums = {4, 6, 5, 3, 2, 1};
  InsertSort(nums, 6); // cout => 1,2,3,4,5,6

  return 0;
}
