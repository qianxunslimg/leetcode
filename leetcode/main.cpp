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

  //   int findNthDigit(int n) {
  //     int d = 1, count = 9;
  //     while (n > (long)d * count) {
  //       n -= d * count;
  //       d++;
  //       count *= 10;
  //     }
  //     int index = n - 1;
  //     int start = (int)pow(10, d - 1);
  //     int num = start + index / d;
  //     int digitIndex = index % d;
  //     int digit = (num / (int)(pow(10, d - digitIndex - 1))) % 10;
  //     return digit;
  //   }

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
    int smalsl = INT_MAX, mid = INT_MAX;
    for (auto num : nums) {
      if (num <= smalsl) {
        smalsl = num;
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

  // vector<int> path;
  bool backTrace(int startIndex, vector<int> &nums) {
    if (path.size() == 3) {
      return true;
    }

    for (int i = startIndex; i < nums.size(); i++) {
      if (!path.empty() && nums[i] <= path.back())
        continue;
      path.push_back(nums[i]);
      if (backTrace(i + 1, nums))
        return true;
      path.pop_back();
    }
    return false;
  }

  bool increasingTriplet4(vector<int> &nums) {
    unordered_set<int> set(nums.begin(), nums.end());
    if (set.size() < 3)
      return false;
    return backTrace(0, nums);
  }

  int dominantIndex(vector<int> &nums) {
    int m1 = -1, m2 = -1;
    int index = -1;
    for (int i = 0; i < nums.size(); i++) {
      if (nums[i] > m1) {
        m2 = m1;
        m1 = nums[i];
        index = i;
      } else if (nums[i] > m2) {
        m2 = nums[i];
      }
    }
    return m1 >= m2 * 2 ? index : -1;
  }

  vector<vector<int>> kSmallestPairs(vector<int> &nums1, vector<int> &nums2,
                                     int k) {
    int pos1 = 0;
    int pos2 = 0;
    vector<vector<int>> ans;
    while (ans.size() < k) {
      ans.push_back(vector<int>{nums1[pos1], nums2[pos2]});
      if (pos1 == nums1.size() - 1 || pos2 == nums2.size() - 1) {
        if (pos1 == nums1.size() - 1) {
          pos2++;
        } else
          pos1++;
      } else if (nums1[pos1 + 1] < nums2[pos2 + 1])
        pos1++;
      else
        pos2++;
    }
    return ans;
  }

  vector<vector<int>> ans;
  // vector<int> path;

  void backtracking(vector<int> cand, int target, int sum) {
    if (sum > target)
      return;
    if (sum == target) {
      sort(path.begin(), path.end());
      ans.push_back(path);
      return;
    }
    for (int i = 0; i < cand.size(); i++) {
      sum += cand[i];
      path.push_back(cand[i]);
      backtracking(cand, target, sum);
      sum -= cand[i];
      path.pop_back();
    }
  }
  vector<vector<int>> combinationSum(vector<int> &candidates, int target) {
    backtracking(candidates, target, 0);
    sort(ans.begin(), ans.end());
    ans.erase(unique(ans.begin(), ans.end()), ans.end());
    return ans;
  }

  vector<int> decrypt(vector<int> &code, int k) {
    vector<int> ans(code.size(), 0);
    if (k == 0)
      return ans;
    int size = code.size();
    int sum = 0;
    if (k > 0) {
      for (int i = 0; i < code.size(); i++) {
        if (i < k) {
          sum += code[(i + 1) % size];
        } else if (i == k)
          ans[i - k] = sum;
        else {
          sum -= code[(i - k + 1) % size];
          sum += code[(i + 1) % size];
          ans[i - k] = sum;
        }
      }
    }
    return ans;
  }

  string reorganizeString(string s) {
    unordered_map<char, int> mapp;
    string ans;
    int maxCount = INT_MIN;
    for (auto ch : s) {
      mapp[ch]++;
      maxCount = maxCount > mapp[ch] ? maxCount : mapp[ch];
    }
    if (maxCount > ((s.size() + 1) / 2))
      return ans;
    vector<pair<int, int>> vec(mapp.begin(), mapp.end());
    sort(vec.begin(), vec.end(),
         [](pair<int, int> a, pair<int, int> b) -> bool {
           return a.second > b.second;
         }); // 给频率排个序
    ans = s;
    int index = 0; // 先按奇数位散开
    for (int i = 0; i < vec.size(); i++) {
      while (vec[i].second--) {
        ans[index] = vec[i].first;
        index += 2;
        if (index >= s.size())
          index = 1; // 奇数位插满了插偶数位
      }
    }
    return ans;
  }

  //   string longestNiceSubstring(string s) {
  //     string ans;
  //     int maxLen = INT_MIN;
  //     for (int i = 0; i < s.size(); i++) {
  //       for (int j = 1; j + i <= s.size(); j++) {
  //         if (j > maxLen) {
  //           string ss = s.substr(i, j);
  //           if (isNiceString(s.substr(i, j))) {
  //             maxLen = j;
  //             ans = s.substr(i, j);
  //             cout << ans << endl;
  //           }
  //         }
  //       }
  //     }
  //     return ans;
  //   }
  //   bool isNiceString(string s) {
  //     set<char> set1(s.begin(), s.end());
  //     set<char> set2;
  //     for (auto ch : s)
  //       set2.insert(tolower(ch));
  //     return set1.size() == (set2.size() * 2);
  //   }
  string longestNiceSubstring(string s) {
    int n = s.size();
    int maxPos = 0;
    int maxLen = 0;
    for (int i = 0; i < n; ++i) {
      int lower = 0;
      int upper = 0;
      for (int j = i; j < n; ++j) {
        if (islower(s[j])) {
          int l = s[j] - 'a';
          int ll = 1 << (s[j] - 'a');
          lower |= 1 << (s[j] - 'a');
        } else {
          upper |= 1 << (s[j] - 'A');
        }
        if (lower == upper && j - i + 1 > maxLen) {
          maxPos = i;
          maxLen = j - i + 1;
        }
      }
    }
    return s.substr(maxPos, maxLen);
  }

  string longestDiverseString(int a, int b, int c) {
    string res;
    vector<pair<int, char>> arr = {{a, 'a'}, {b, 'b'}, {c, 'c'}};

    while (true) {
      sort(arr.begin(), arr.end(),
           [](const pair<int, char> &p1, const pair<int, char> &p2) {
             return p1.first > p2.first;
           });
      bool hasNext = false;
      for (auto &pairr : arr) {
        int freq = pairr.first;
        char ch = pairr.second;
        int m = res.size();
        if (freq <= 0) {
          break;
        }
        if (m >= 2 && res[m - 2] == ch && res[m - 1] == ch) {
          continue;
        }
        hasNext = true;
        res.push_back(ch);
        freq--;
        break;
      }
      if (!hasNext) {
        break;
      }
    }

    return res;
  }

  //暴力
  int threeSumClosest(vector<int> &nums, int target) {
    int ans = INT_MAX;
    int min_dis = INT_MAX;
    int n = nums.size();
    for (int i = 0; i < n; i++) {
      for (int j = i + 1; j < n; j++) {
        for (int z = j + 1; z < n; z++) {
          int sum = nums[i] + nums[j] + nums[z]; //和
          int dis = abs(target - sum);
          if (dis < min_dis) {
            ans = sum;
            min_dis = dis;
          }
        }
      }
    }
    return ans;
  }

  vector<int> findSubstring(string s, vector<string> &words) {
    vector<int> res;
    map<string, int> keyMapCount, keyMapCountTemp;
    int len = words[0].size(), totalLen = len * words.size();
    for (auto it : words)
      keyMapCount[it]++;
    for (int i = 0; i < s.size() - totalLen + 1; i++) {
      int j = i;
      keyMapCountTemp = keyMapCount;
      for (; j < i + totalLen; j += len) {
        string temp = s.substr(j, len);
        if (keyMapCountTemp[temp] == 0)
          break;
        keyMapCountTemp[temp]--;
      }
      if (j == i + totalLen)
        res.push_back(i);
    }
    return res;
  }

  string convert(string s, int numRows) {
    if (numRows == 1)
      return s;

    vector<string> rows(min(numRows, int(s.size())));
    int curRow = 0;
    bool goingDown = false;

    for (char c : s) {
      rows[curRow] += c;
      if (curRow == 0 || curRow == numRows - 1)
        goingDown = !goingDown;
      curRow += goingDown ? 1 : -1;
    }

    string ret;
    for (string row : rows)
      ret += row;
    return ret;
  }

  //最优解， 反转一半的数字
  bool isPalindrome(int x) {
    // 特殊情况：
    // 如上所述，当 x < 0 时，x 不是回文数。
    // 同样地，如果数字的最后一位是 0，为了使该数字为回文，
    // 则其第一位数字也应该是 0
    // 只有 0 满足这一属性
    if (x < 0 || (x % 10 == 0 && x != 0)) {
      return false;
    }
    int revertedNumber = 0;
    while (x > revertedNumber) {
      revertedNumber = revertedNumber * 10 + x % 10;
      x /= 10;
    }
    // 当数字长度为奇数时，我们可以通过 revertedNumber/10 去除处于中位的数字。
    // 例如，当输入为 12321 时，在 while 循环的末尾我们可以得到 x =
    // 12，revertedNumber = 123，
    // 由于处于中位的数字不影响回文（它总是与自己相等），所以我们可以简单地将其去除。
    return x == revertedNumber || x == revertedNumber / 10;
  }

  bool isMatch(string s, string p) {
    bool dp[21][30] = {};
    dp[0][0] = true;
    int szs = s.size(), szp = p.size();
    for (int i = 0; i <= szs; ++i) {
      for (int j = 1; j <= szp; ++j) {
        if (i) {
          dp[i][j] |=
              (dp[i - 1][j - 1] &&
               (s[i - 1] == p[j - 1] || p[j - 1] == '.')); // s[i]、s[j]两两匹配
          if (j > 1)
            dp[i][j] |=
                (dp[i - 1][j] && (s[i - 1] == p[j - 2] || p[j - 2] == '.') &&
                 p[j - 1] == '*'); //*匹配多个：aaaa a*  aaaa .*
        }
        if (j > 1) {
          dp[i][j] |=
              (dp[i][j - 2] && p[j - 1] == '*'); //*将前一个字符跳过： ab c*abc*
        }
        dp[i][j] |=
            (dp[i][j - 1] &&
             p[j] == '*'); //*不用： a a*
                           // cout << i << ' ' << j << ' ' << dp[i][j] << endl;
      }
    }
    return dp[szs][szp];
  }

  string longestCommonPrefix(vector<string> &strs) {
    if (strs.empty())
      return "";
    const auto p = minmax_element(strs.begin(), strs.end());
    for (int i = 0; i < p.first->size(); ++i) {
      if (p.first->at(i) != p.second->at(i))
        return p.first->substr(0, i);
    }
    return *p.first;
  }

  vector<int> pancakeSort(vector<int> &arr) {
    vector<int> ans;
    int n = arr.size();
    auto isSorted = [](vector<int> arry) -> bool {
      vector<int> copy = arry;
      sort(arry.begin(), arry.end());
      return copy == arry;
    };
    while (!isSorted(arr)) {
      //       int maxx = *max_element(arr.begin(), arr.begin() + n);
      //       int pos = 0;
      //       for (int i = 0; i < n; i++) {
      //         if (maxx == arr[i]) {
      //           pos = i;
      //           break;
      //         }
      //       }
      int pos = max_element(arr.begin(), arr.begin() + n) - arr.begin();
      reverse(arr.begin(), arr.begin() + pos + 1);
      reverse(arr.begin(), arr.begin() + n);
      ans.push_back(pos + 1);
      ans.push_back(n);
      n--;
    }
    return ans;
  }

  bool isOneBitCharacter(vector<int> &bits) {
    int n = bits.size();
    bool ans = 0;
    int i = 0;
    for (; i < n - 1; i++) {
      if (bits[i] == 1) {
        i++;
      }
    }
    return ans;
  }
  string pushDominoes(string dominoes) {
    int n = dominoes.size();
    queue<int> q;
    vector<int> time(n, -1);
    vector<string> force(n);
    for (int i = 0; i < n; i++) {
      if (dominoes[i] != '.') {
        q.emplace(i);
        time[i] = 0;
        force[i].push_back(dominoes[i]);
      }
    }

    string res(n, '.');
    while (!q.empty()) {
      int i = q.front();
      q.pop();
      if (force[i].size() == 1) {
        char f = force[i][0];
        res[i] = f;
        int ni = (f == 'L') ? (i - 1) : (i + 1);
        if (ni >= 0 && ni < n) {
          int t = time[i];
          if (time[ni] == -1) {
            q.emplace(ni);
            time[ni] = t + 1;
            force[ni].push_back(f);
          } else if (time[ni] == t + 1) {
            force[ni].push_back(f);
          }
        }
      }
    }
    return res;
  }

  string reverseOnlyLetters(string s) {
    int left = 0, right = s.size() - 1;
    while (left < right) {
      if (isalpha(s[left]) && isalpha(s[right])) {
        swap(s[left], s[right]);
        left++;
        right--;
      } else {
        left += isalpha(s[left]) ? 0 : 1;
        right -= isalpha(s[right]) ? 0 : 1;
      }
    }
    return s;
  }

  vector<int> searchRange(vector<int> &nums, int target) {
    int l = 0, r = nums.size();
    int pos1 = 0;
    int pos2 = nums.size() - 1;
    while (l < r) {
      int mid = l + (r - l) / 2;
      if (nums[mid] > target) {
        r = mid;
      } else if (nums[mid] < target) {
        l = mid + 1;
      } else { //找到
        pos1 = pos2 = mid;
        while (pos1 >= 0 && nums[pos1] == target) {
          pos1--;
        }
        while (pos2 < nums.size() && nums[pos2] == target) {
          pos2++;
        }
        return {pos1 + 1, pos2 - 1};
      }
    }
    return {-1, -1};
  }

  // int search(vector<int> &nums, int target) {
  //  int l = 0, r = nums.size() - 1;
  //  bool findArea = 0;
  //  while (l <= r) {
  //    int mid = l + (r - l) / 2;
  //    if (!findArea) {
  //      if (nums[l] < nums[mid] && nums[mid] < nums[r]) {
  //        findArea = 1;
  //        continue;
  //      }
  //      if (nums[l] < target)
  //        r = mid - 1;
  //      else if (nums[l] > target)
  //        l = mid;
  //      else
  //        return l;
  //    }
  //    if (findArea) {
  //      if (target < nums[mid]) {
  //        r = mid - 1;
  //      }
  //      if (nums[mid] < target) {
  //        l = mid + 1;
  //      } else
  //        return mid;
  //    }
  //  }
  //  return -1;
  //}

  int search(vector<int> &nums, int target) {
    int l = 0, r = nums.size() - 1;
    bool findArea = 0;
    while (l <= r) {
      int mid = l + (r - l) / 2;

      //查找转折点
      if (!findArea) {
        if (nums[l] <= nums[mid] && nums[mid] <= nums[r]) {
          findArea = 1;
          continue;
        }
        if (nums[mid] > nums[l]) {
          l = mid + 1;
        } else if (nums[mid] < nums[r]) {
          r = mid - 1;
        }
      }
      if (findArea) {
        if (target < nums[mid]) {
          r = mid - 1;
        }
        if (nums[mid] < target) {
          l = mid + 1;
        } else
          return mid;
      }
    }
    return -1;
  }

  vector<vector<int>> dir{{1, 0}, {0, 1}, {-1, 0}, {0, -1}}; //顺时针方向
  vector<int> anss;
  int count = 0;
  vector<int> spiralOrder(vector<vector<int>> &matrix) {
    anss.clear();
    dfs(matrix, 0, 0, 0);
    return anss;
  }

  void dfs(vector<vector<int>> &matrix, int x, int y, int index) {
    if (x < 0 || y < 0 || x >= matrix.size() || y >= matrix[0].size() ||
        matrix[x][y] == -1) {
      index++;
      if (anss.size() == matrix.size() * (matrix[0].size()))
        return;
      dfs(matrix, x - dir[(index - 1) % 4][0] + dir[index % 4][0],
          y - dir[(index - 1) % 4][1] + dir[index % 4][1], index);
    }
    anss.push_back(matrix[x][y]);

    matrix[x][y] = -1;
    dfs(matrix, x + dir[index % 4][0], y + dir[index % 4][1], index);
  }

  vector<int> ansdd;
  vector<int> findBall(vector<vector<int>> &grid) {
    ansdd.clear();
    dfss(grid, 0, 1, 0);
    return ansdd;
  }

  void dfss(vector<vector<int>> &grid, int x, int y, int from) {
    if (y < 0 || y >= grid[0].size()) {
      ansdd.push_back(-1);
      return;
    }
    int temp = from + grid[x][y];
    if (temp == 0) { // 1 + -1 = 0 死角
      ansdd.push_back(-1);
      return;
    }
    if (x == grid.size() - 1 && from != 0) { //到底
      ansdd.push_back(y);
      return;
    }

    if (from == 0) { //从上面来 向左右运动
      dfss(grid, x, y + grid[x][y], grid[x][y]);
    }
    //方向相同 向下
    if (from == grid[x][y])
      dfss(grid, x + 1, y, 0);
  }

  //
  double myPow(double x, int n) {
    double res = 1.0;
    for (int i = n; i != 0; i /= 2) {
      if (i % 2 != 0) { //奇数额外乘上 同时能保证结束为1 乘到res上
        res *= x;
      }
      x *= x;
    }
    return n < 0 ? 1 / res : res;
  }

  double myyPow(double x, int n) {
    double res = 1.0;
    int i = n;
    while (i != 0) { //判断条件 可能是负数
      if (i % 2 != 0) {
        res *= x;
      }
      i /= 2;
      x *= x;
    }
    return n < 0 ? 1 / res : res;
  }

  // bool isNumber(string s) {
  //  int n = s.size();
  //  int l = 0;
  //  while (s[l] == ' ')
  //    l++;
  //  int r = n - 1;
  //  while (s[r - 1] == ' ')
  //    r--;

  //  int s1 = 0; // e 或者.的个数
  //  int epos = 0;
  //  for (int i = l; i <= r; i++) {
  //    if (isalpha(s[i])) {
  //      if (s[i] != 'e')
  //        return 0; //错误 非e字母
  //    }
  //    if (s1 > 1)
  //      return 0; // e . 超个数
  //    //一个符号字符（'+' 或 '-'）
  //    if (s[i] == '+' || s[i] == '-') {
  //      if (i == l)
  //        continue;
  //      else
  //        return 0; //则正负只能在首位
  //    }
  //    if (epos == 0) {
  //      if (s[i] == '.' || s[i] == 'e') {
  //        epos = i;
  //        s1++;
  //      }
  //      if (epos == r)
  //        return 0;
  //    } else {
  //      if (s[i] == '.' || s[i] == 'e') {
  //        return 0;
  //      }
  //    }
  //  }
  //  return 1;
  //}

  //   bool isNumber(string s) {
  //     int n = s.size();
  //     int l = 0;
  //     if (s[0] == 'e')
  //       return 0;
  //     while (s[l] == ' ')
  //       l++;
  //     if (l == n) {
  //       return 0;
  //     }
  //     int r = n - 1;
  //     while (s[r] == ' ')
  //       r--;
  //
  //     int s1 = 0; // e 或者.的个数
  //     int epos = 0;
  //     for (int i = l; i <= r; i++) {
  //       if (isalpha(s[i])) {
  //         if (s[i] != 'e' && s[i] != 'E')
  //           return 0; //错误 非e字母
  //       }
  //       if (s1 > 1)
  //         return 0; // e . 超个数
  //                   //一个符号字符（'+' 或 '-'）
  //       if (s[i] == '+' || s[i] == '-') {
  //         if (i == l)
  //           continue;
  //         else if (s[i - 1] == 'e' || s[i - 1] == 'E') {
  //           continue;
  //         } else
  //           return 0; //则正负只能在首位
  //       }
  //       if (s[i] == '.' || s[i] == 'e' || s[i] == 'E') {
  //         if (i == r) {
  //           return 0;
  //         }
  //         s1++;
  //       }
  //     }
  //     return 1;
  //   }

  int i = 0;
  //扫描空格
  void getSpace(string &s) { //移动前后space
    while (i < s.size() && s[i] == ' ')
      i++;
  }
  //扫描有符号整数
  bool getInt(string &s) { //匹配 +212313  -14556 1231
    if (i < s.size() && (s[i] == '+' || s[i] == '-'))
      i++;
    return getUint(s);
  }
  //扫描无符号整数
  bool getUint(string &s) { //匹配 12312
    int tmp = i;
    while (i < s.size() && isdigit(s[i]))
      i++;
    return i > tmp;
  }
  bool isNumber(string s) {
    if (s == "")
      return false;
    getSpace(s);
    bool flag = getInt(s);
    if (i < s.size() && s[i] == '.') {
      i++;
      //当为.时,后面必须是无符号整数，并且.的前后只要有一个为true就行
      //而且必须把getUint(s)放在前面，不然由于||的短视特征，小数点后面可能不会被扫描到
      flag = getUint(s) || flag;
    }
    if (i < s.size() && (s[i] == 'e' || s[i] == 'E')) {
      i++;
      flag = flag && getInt(s); // e的前后都必须为true
    }
    getSpace(s);
    return i == s.size() && flag;
  }

  //二分法细节 查找左边界
  //<写法
  int left_bound(vector<int> nums, int target) {
    if (nums.size() == 0)
      return -1;
    int left = 0;
    int right = nums.size(); // 注意

    while (left < right) { // 注意
      int mid = (left + right) / 2;
      if (nums[mid] == target) {
        right = mid;
      } else if (nums[mid] < target) {
        left = mid + 1;
      } else if (nums[mid] > target) {
        right = mid; // 注意
      }
    }
    // return left;    //返回>=target的左边界位置 [0,nums.size()]

    // 返回第一个target的位置 没有则返回-1；
    {
      if (left == nums.size())
        return -1; //[1,2,2,4]搜索8返回left 4，越界
                   // 类似之前算法的处理方式
      return nums[left] == target ? left : -1;
    }
  }

  //<= 写法 完全一致
  int left_bound2(vector<int> nums, int target) {
    if (nums.size() == 0)
      return -1;
    int left = 0;
    int right = nums.size() - 1; // 注意

    while (left <= right) { // 注意
      int mid = (left + right) / 2;
      if (nums[mid] == target) {
        right = mid - 1;
      } else if (nums[mid] < target) {
        left = mid + 1;
      } else if (nums[mid] > target) {
        right = mid - 1; // 注意
      }
    }
    // return left; //返回>=target的左边界位置 [0,nums.size()]

    // 返回第一个target的位置 没有则返回-1；
    {
      if (left == nums.size())
        return -1; //[1,2,2,4]搜索8返回left 4，越界
                   // 类似之前算法的处理方式
      return nums[left] == target ? left : -1;
    }
  }

  //二分法细节 查找右边界
  //<写法
  int right_bound(vector<int> nums, int target) {
    if (nums.size() == 0)
      return -1;
    int left = 0, right = nums.size();

    while (left < right) {
      int mid = (left + right) / 2;
      if (nums[mid] == target) {
        left = mid + 1; // 注意
      } else if (nums[mid] < target) {
        left = mid + 1;
      } else if (nums[mid] > target) {
        right = mid;
      }
    }
    // return left - 1; //返回>=target的右边界位置 [0,nums.size()]

    // 返回最后一个target的位置 没有则返回-1；
    {
      if (left == 0)
        return -1; //这个例子搜索0 就是返回left 0
      return nums[left - 1] == target ? (left - 1) : -1;
    }
  }

  //<=写法
  int right_bound2(vector<int> nums, int target) {
    if (nums.size() == 0)
      return -1;
    int left = 0, right = nums.size() - 1;

    while (left <= right) {
      int mid = (left + right) / 2;
      if (nums[mid] == target) {
        left = mid + 1; // 注意
      } else if (nums[mid] < target) {
        left = mid + 1;
      } else if (nums[mid] > target) {
        right = mid - 1;
      }
    }
    // return left - 1; //返回>=target的右边界位置 [0,nums.size()]

    // 返回最后一个target的位置 没有则返回-1；
    {
      if (left == 0)
        return -1; //这个例子搜索0 就是返回left 0
      return nums[left - 1] == target ? (left - 1) : -1;
    }
  }

  int kthSmallest(vector<vector<int>> &matrix, int k) {
    struct point {
      int val, x, y;
      point(int val, int x, int y) : val(val), x(x), y(y) {}
      bool operator>(const point &a) const { return this->val > a.val; }
    };
    priority_queue<point, vector<point>, greater<point>> que;
    int n = matrix.size();
    for (int i = 0; i < n; i++) {
      que.emplace(matrix[i][0], i, 0); //每行首元素压进
    }
    for (int i = 0; i < k - 1; i++) {
      point now = que.top();
      que.pop();
      if (now.y != n - 1) { //一行到头 会自动跳到下一行
        que.emplace(matrix[now.x][now.y + 1], now.x, now.y + 1);
      }
    }
    return que.top().val;
  }

  bool check(vector<vector<int>> &matrix, int mid, int k, int n) {
    int i = n - 1;
    int j = 0;
    int num = 0;
    //每次对于「猜测」的答案 midmid，计算矩阵中有多少数不大于 mid
    //如果数量不少于 k，那么说明最终答案 x 不大于 mid；
    //如果数量少于 k，那么说明最终答案 x 大于 mid。
    while (i >= 0 && j < n) {
      if (matrix[i][j] <= mid) {
        num += i + 1;
        j++;
      } else {
        i--;
      }
    }
    return num >= k;
  }

  int kthSmallest2(vector<vector<int>> &matrix, int k) {
    int n = matrix.size();
    int left = matrix[0][0];
    int right = matrix[n - 1][n - 1];
    while (left < right) {
      int mid = left + ((right - left) >> 1);
      if (check(matrix, mid, k, n)) { //<=mid的个数>=k 找左边界
        right = mid;                  //向左上角收缩
      } else {
        left = mid + 1; //向右下角扩大
      }
    }
    return left;
  }

  string minWindow(string s, string t) {
    unordered_map<char, int> need, window;
    for (char c : t)
      need[c]++;
    int left = 0, right = 0;
    int valid = 0;
    //记录最小覆盖字串的其实索引和长度
    int start = 0, len = INT_MAX;
    while (right < s.size()) {
      // c是移入窗口的字符
      char c = s[right];
      right++;
      // 进行窗口内数据的一系列更新
      if (need.count(c)) {
        window[c]++;
        if (window[c] == need[c])
          valid++;
      }

      //判断左窗口是否需要收缩
      while (valid == need.size()) { //窗口满足条件
                                     // 在这里更新最小覆盖子串
        if (right - left < len) {
          start = left;
          len = right - left;
        }
        // d 是将移出窗口的字符
        char d = s[left];
        // 左移窗口
        left++;
        // 进行窗口内数据的一系列更新
        if (need.count(d)) {
          if (window[d] == need[d]) {
            valid--;
          }
          window[d]--;
        }
      }
    }
    return len == INT_MAX ? "" : s.substr(start, len);
  }

  int lengthOfLongestSubstring22(string s) {
    int ans = 0;
    int left = 0, right = 0;
    unordered_map<char, int> window;
    while (right < s.size()) {
      char c = s[right];
      right++;
      window[c]++;
      while (window[c] > 1) { //有重复就要从left++ 知道消除当前重复
        char d = s[left];
        left++;
        window[d]--;
      }
      ans = max(ans, right - left);
    }
    return ans;
  }

  int lengthOfLIS(vector<int> &nums) {
    int n = nums.size();
    vector<int> dp(n + 1, 1);
    dp[1] = 1;
    int ans = 1;
    for (int i = 1; i < n; i++) {
      dp[i] = nums[i] > nums[i - 1] ? dp[i - 1] + 1 : dp[i - 1];
      ans = max(ans, dp[i]);
    }
    return ans;
  }

  long long minimalKSum(vector<int> &nums, int k) {
    sort(nums.begin(), nums.end());
    long long ans = (long long)k * (k + 1) / 2;
    long long last = k;
    long long pre = -1;
    for (int num : nums) {
      if (num == pre)
        continue;
      if (num <= last) {
        ans += last + 1 - num;
        last++;
        pre = num;
      } else
        return ans;
    }
    return ans;
  }

  int findNthDigit(long k) {
    for (int i = 1;; i++) {
      if (i * pow(10, i) > k) {
        string s = to_string(k / i);
        return to_string(k / i)[k % i] - '0';
      }
      k += pow(10, i);
    }
  }

  bool canPartition(vector<int> &nums) {
    int allSum = accumulate(begin(nums), end(nums), 0);
    if (allSum % 1)
      return 0;
    int target = allSum / 2;

    // dp[i]中的i表示背包内总和
    // 题目中说：每个数组中的元素不会超过 100，数组的大小不会超过 200
    // 总和不会大于20000，背包最大只需要其中一半，所以10001大小就可以
    vector<int> dp(10001, 0);
    // begin 0/1
    for (int i = 0; i < nums.size(); i++) {
      for (int j = target; j >= nums[i]; j--) {
        dp[j] = max(dp[j], dp[j - nums[i]] + nums[i]);
      }
    }
    return dp[target] == target;
  }

  vector<int> addToArrayForm(vector<int> &num, int k) {
    vector<int> res;
    int n = num.size();
    for (int i = n - 1; i >= 0; --i) {
      int sum = num[i] + k % 10;
      k /= 10;
      if (sum >= 10) {
        k++; //相当于carry进位
        sum -= 10;
      }
      res.push_back(sum);
    }
    for (; k > 0; k /= 10) {
      res.push_back(k % 10);
    }
    reverse(res.begin(), res.end());
    return res;
  }

  int nthUglyNumber(int n) {
    vector<int> factors = {2, 3, 5};
    unordered_set<long> seen;
    priority_queue<long, vector<long>, greater<long>> heap;
    seen.insert(1L);
    heap.push(1L);
    int ugly = 0;
    for (int i = 0; i < n; i++) {
      long curr = heap.top();
      heap.pop();
      ugly = (int)curr;
      for (int factor : factors) {
        long next = curr * factor;
        if (!seen.count(next)) {
          seen.insert(next);
          heap.push(next);
        }
      }
    }
    return ugly;
  }

  int missingNumber4(vector<int> &nums) {
    //二分查找
    int n = nums.size();
    int left = 0, right = n; //左闭右闭
    while (left < right) {
      int mid = left + (right - left) / 2;
      if (nums[mid] >= mid) {
        right--;
      } else
        left++;
    }
    return left;
  }

  // 普通回溯过不了，需要精准剪枝到第k个叶节点
  void dsssfs(int n, int k, unordered_set<int> &used, string &tmp,
              vector<int> &factorial) {
    if (tmp.size() == n) {
      return;
    }
    int ind =
        0; // 用来标记当前是第几次循环,直接用i的话有问题，比如说i是3，但是只是第一次循环，那就错了
    for (int i = 1; i <= n; ++i) {
      if (used.find(i) != used.end())
        continue;
      ++ind;
      // 需要看当前层切分后每个子节点包含的叶节点个数，所以要减一
      int size = factorial[n - used.size() - 1];

      if (k > (ind - 1) * size && k <= ind * size) {
        tmp.push_back(i + '0');
        used.insert(i);
        dsssfs(n, k - size * (ind - 1), used, tmp, factorial);
        // 无需回溯，因为从dfs出来后就已经是结果了
      }
    }
  }

  string getPermutation(int n, int k) {
    unordered_set<int> used;
    string tmp;
    // 提前吧阶乘算出来
    vector<int> factorial = {1, 1, 2, 6, 24, 120, 720, 5040, 40320, 362880};
    dsssfs(n, k, used, tmp, factorial);
    return tmp;
  }
  void nextPermuwwwtation(vector<int> &nums) {
    int i = nums.size() - 2;
    while (i >= 0 && nums[i] >= nums[i + 1]) {
      i--;
    }
    if (i >= 0) {
      int j = nums.size() - 1;
      while (j >= 0 && nums[i] >= nums[j]) {
        j--;
      }
      swap(nums[i], nums[j]);
    }
    reverse(nums.begin() + i + 1, nums.end());
  }

  void setZeroes(vector<vector<int>> &matrix) {
    int m = matrix.size();
    int n = matrix[0].size();
    int flag_col0 = false, flag_row0 = false;
    for (int i = 0; i < m; i++) {
      if (!matrix[i][0]) {
        flag_col0 = true;
      }
    }
    for (int j = 0; j < n; j++) {
      if (!matrix[0][j]) {
        flag_row0 = true;
      }
    }
    for (int i = 1; i < m; i++) {
      for (int j = 1; j < n; j++) {
        if (!matrix[i][j]) {
          matrix[i][0] = matrix[0][j] = 0;
        }
      }
    }
    for (int i = 1; i < m; i++) {
      for (int j = 1; j < n; j++) {
        if (!matrix[i][0] || !matrix[0][j]) {
          matrix[i][j] = 0;
        }
      }
    }
    if (flag_col0) {
      for (int i = 0; i < m; i++) {
        matrix[i][0] = 0;
      }
    }
    if (flag_row0) {
      for (int j = 0; j < n; j++) {
        matrix[0][j] = 0;
      }
    }
  }
};

void solve() {
  Solution mysolution;
  vector<int> a{3, 0, 0, 0, 0, 2};
  vector<int> b{3, 0, 0, 0, 0, 2};
  mysolution.setZeroes(
      vector<vector<int>>{{0, 1, 2, 0}, {3, 4, 5, 2}, {1, 3, 1, 5}});
  cout << "pause!" << endl;
}

int main() {
  solve();
  cout << "program has finished..." << endl;
  return 0;
}

// class B;
// class A {
// public:
//   shared_ptr<B> pb_;
//   ~A() { cout << "A delete\n"; }
// };
//
// class B {
// public:
//   shared_ptr<A> pa_;
//   ~B() { cout << "B delete\n"; }
// };

// void fun() {
//   shared_ptr<B> pb(new B());
//   shared_ptr<A> pa(new A());
//   pb->pa_ = pa;
//   pa->pb_ = pb;
//   cout << pb.use_count() << endl;
//   cout << pa.use_count() << endl;
// }

// int main() {
//   fun();
//   int a = 0;
//   return 0;
// }

//⼿写实现智能指针类
template <typename T> class SharedPtr {
private:
  size_t *m_count_;
  T *m_ptr_;

public:
  //构造函数
  SharedPtr() : m_ptr_(nullptr), m_count_(new size_t) {}
  SharedPtr(T *ptr) : m_ptr_(ptr), m_count_(new size_t) { m_count_ = 1; }
  //析构函数
  ~SharedPtr() {
    --(*m_count_);
    if (*m_count_ == 0) {
      delete m_ptr_;
      delete m_count_;
      m_ptr_ = nullptr;
      m_count_ = nullptr;
    }
  }
  //拷⻉构造函数
  SharedPtr(const SharedPtr &ptr) {
    m_count_ = ptr.m_count_;
    m_ptr_ = ptr.m_ptr_;
    ++(*m_count_);
  }
  //拷⻉赋值运算
  void operator=(const SharedPtr &ptr) { SharedPtr(std::move(ptr)); }
  //移动构造函数
  SharedPtr(SharedPtr &&ptr) : m_ptr_(ptr.m_ptr_), m_count_(ptr.m_count_) {
    ++(*m_count_);
  }
  //移动赋值运算
  void operator=(SharedPtr &&ptr) { SharedPtr(std::move(ptr)); }
  //解引⽤
  T &operator*() { return *m_ptr_; }
  //箭头运算
  T *operator->() { return m_ptr_; }
  //᯿载bool操作符
  operator bool() { return m_ptr_ == nullptr; }
  T *get() { return m_ptr_; }
  size_t use_count() { return *m_count_; }
  bool unique() { return *m_count_ == 1; }
  void swap(SharedPtr &ptr) { std::swap(*this, ptr); }
};
#include <assert.h>
//⼿写字符串函数 strcat，strcpy，strncpy，memset，memcpy实现
//把 src 所指向的字符串复制到 dest，注意：dest定义的空间应该⽐src⼤。
char *strcpy(char *dest, const char *src) {
  char *ret = dest;
  assert(dest != NULL); //优化点1：检查输⼊参数
  assert(src != NULL);
  while (*src != '\0')
    *(dest++) = *(src++);
  *dest = '\0'; //优化点2：⼿动地将最后的'\0'补上
  return ret;
}
//考虑内存᯿叠的字符串拷⻉函数 优化的点
char *strcpy(char *dest, char *src) {
  char *ret = dest;
  assert(dest != NULL);
  assert(src != NULL);
  memmove(dest, src, strlen(src) + 1);
  return ret;
}
//把 src 所指向的字符串追加到 dest 所指向的字符串的结尾。
char *strcat(char *dest, const char *src) {
  // 1. 将⽬的字符串的起始位置先保存，最后要返回它
  // 2. 先找到dest的结束位置,再把src拷⻉到dest中，记得在最后要加上'\0'
  char *ret = dest;
  assert(dest != NULL);
  assert(src != NULL);
  while (*dest != '\0')
    dest++;
  while (*src != '\0')
    *(dest++) = *(src++);
  *dest = '\0';
  return ret;
}
//把 str1 所指向的字符串和 str2 所指向的字符串进⾏⽐较。
//该函数返回值如下：
//如果返回值 < 0，则表示 str1 ⼩于 str2。
//如果返回值 > 0，则表示 str1 ⼤于 str2。
//如果返回值 = 0，则表示 str1 等于 str2。
int strcmp(const char *s1, const char *s2) {
  assert(s1 != NULL);
  assert(s2 != NULL);
  while (*s1 != '\0' && *s2 != '\0') {
    if (*s1 > *s2)
      return 1;
    else if (*s1 < *s2)
      return -1;
    else {
      s1++, s2++;
    }
  }
  //当有⼀个字符串已经⾛到结尾
  if (*s1 > *s2)
    return 1;
  else if (*s1 < *s2)
    return -1;
  else
    return 0;
}
//在字符串 str1 中查找第⼀次出现字符串 str2 的位置，不包含终⽌符 '\0'。
char *strstr(char *str1, char *str2) {
  char *s = str1;
  assert(str1 != '\0');
  assert(str2 != '\0');
  if (*str2 == '\0')
    return NULL;       //若str2为空，则直接返回空
  while (*s != '\0') { //若不为空，则进⾏查询
    char *s1 = s;
    char *s2 = str2;
    while (*s1 != '\0' && *s2 != '\0' && *s1 == *s2)
      s1++, s2++;
    if (*s2 == '\0')
      return str2; //若s2先结束
    if (*s2 != '\0' && *s1 == '\0')
      return NULL; //若s1先结束⽽s2还没结束，则返回空
    s++;
  }
  return NULL;
}
//模拟实现memcpy函数 从存储区 str2 复制 n 个字符到存储区 dst。
void *memcpy(void *dest, void *src, size_t num) {
  void *ret = dest;
  size_t i = 0;
  assert(dest != NULL);
  assert(src != NULL);
  for (i = 0; i < num; i++) {
    //因为void* 不能直接解引⽤，所以需要强转成char*再解引⽤
    //此处的void*实现了泛型编程
    *(char *)dest = *(char *)src;
    dest = (char *)dest + 1;
    src = (char *)src + 1;
  }
  return ret;
}

//考虑内存chong叠的memcpy函数 优化的点
void *memmove(void *dest, void *src, size_t num) {
  char *p1 = (char *)dest;
  char *p2 = (char *)src;
  if (p1 < p2) { // p1低地址p2⾼地址
    for (size_t i = 0; i != num; ++i)
      *(p1++) = *(p2++);
  } else {
    //从后往前赋值
    p1 += num - 1;
    p2 += num - 1;
    for (size_t i = 0; i != num; ++i)
      *(p1--) = *(p2--);
  }
  return dest;
}

// 7、说⼀下 ++i 和 i++ 的区别
// 菜鸡⼩贺： 这题我会！++i （前置加加）先⾃增 1再返回，i++ （后置加加）先返回 i
// 再⾃增 1。
// 前置加加不会产⽣临时对象，后置加加必须产⽣临时对象，临时对象会导致效率降低
// ++i 实现：
// int& int::operator++ () {
// 	*this += 1；
// 		return *this；
// }
// i++ 实现：
// const int int::operator（int）{
// 	int oldValue = *this；
// 	++（*this）；
// 	return oldValue； }

// 9、讲讲⼤端⼩端，如何检测
// ⼤端模式：是指数据的⾼字节保存在内存的低地址中，⽽数据的低字节保存在内存的⾼地址
// 端。
// ⼩端模式，是指数据的⾼字节保存在内存的⾼地址中，低位字节保存在在内存的低地址端。
// 直接读取存放在内存中的⼗六进制数值，取低位进⾏值判断
// ⽤共同体来进⾏判断
// union
// 共同体所有数据成员是共享⼀段内存的，后写⼊的成员数据将覆盖之前的成员数据，成
// 	员数据都有相同的⾸地址。Union 的⼤⼩为最⼤数据成员的⼤⼩。
// 	union 的成员数据共⽤内存，并且⾸地址都是低地址⾸字节。Int i =
// 1时：⼤端存储1放在最⾼ 	位，⼩端存储1放在最低位。当读取char
// ch时，是最低地址⾸字节，⼤⼩端会显示不同的值。
// int a = 0x12345678;
// int *c = &a;
// c[0] == 0x12 ⼤端模式
// c[0] == 0x78 ⼩段模式

void manyzhishidian() {
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

  //   auto InsertSort = [](std::vector<int> &nums, int n) {
  //     if (n <= 1)
  //       return;
  //     for (int i = 0; i < n; ++i) {
  //       for (int j = i; j > 0 && nums[j] < nums[j - 1]; --j) {
  //         std::swap(nums[j], nums[j - 1]);
  //       }
  //     }
  //   };
  //   std::vector<int> nums = {4, 6, 5, 3, 2, 1};
  //   InsertSort(nums, 6); // cout => 1,2,3,4,5,6

  // 对于32位的整数，大端机器会在内存的低地址存储高位，在高地址存储低位。
  // 小端机器恰好相反，内存的低地址存储低位，在高地址存储高位
  int a = 0x12345678;
  char *p;
  p = (char *)&a;
  if (*p == 0x78) {
    cout << *p << endl; // insert 小端
  } else if (*p == 0x12) {
    cout << *p << endl;
  }
}

// #include <algorithm>
// #include <iostream>
// using namespace std;
// class SingleInstance {
// public:
//   static SingleInstance *GetInstance() {
//     static SingleInstance ins;
//     return &ins;
//   }
//   ~SingleInstance(){};
//
// private:
//   //涉及到创建对象的函数都设置为private
//   SingleInstance() { std::cout << "SingleInstance() 饿汉" << std::endl; }
//   SingleInstance(const SingleInstance &other){};
//   SingleInstance &operator=(const SingleInstance &other) { return *this; }
// };
// int main() {
//   //因为不能创建对象所以通过静态成员函数的⽅法返回静态成员变ᰁ
//   SingleInstance *ins = SingleInstance::GetInstance();
//   return 0;
// }
// //输出 SingleInstance() 饿汉
