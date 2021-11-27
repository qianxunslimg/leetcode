#include "main.h"

class Solution {
public:
  //两种排序下标的方法
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
};

void solve() {
  Solution mysolution;
  vector<int> nums1 = {4, 1, 2};
  vector<int> nums2 = {3, 2, 3, 4, 6, 5};
  int target = 1;
  string a = "bza";
  string b = "Hello World";
  vector<int> test = {4, 3, 2, 7, 8, 2, 3, 1};
  vector<vector<int>> ve = {{2, 2}, {3, 3}};
  vector<string> dictionary = {"ale", "apple", "monkey", "plea"};
  vector<string> ress = {"Hello", "Alaska", "Dad", "Peace"};
  string s = "the sky is blue";
  string s1 = "aaaaa";
  string s2 = "abbcd";
  vector<vector<int>> img{
      {1, 1, 3, 3}, {3, 1, 4, 2}, {3, 2, 4, 4}, {1, 3, 2, 4}, {2, 3, 3, 4}};
  // int aa = mysolution.minTimeToType(a);

  int aa = mysolution.longestPalindrome2(
      "civilwartestingwhetherthatnaptionoranynartionsoconceivedandsodedicatedca"
      "nlongendureWeareqmetonagreatbattlefiemldoftzhatwarWehavecometodedicpatea"
      "portionofthatfieldasafinalrestingplaceforthosewhoheregavetheirlivesthatt"
      "hatnationmightliveItisaltogetherfangandproperthatweshoulddothisButinalar"
      "gersensewecannotdedicatewecannotconsecratewecannothallowthisgroundThebra"
      "velmenlivinganddeadwhostruggledherehaveconsecrateditfaraboveourpoorponwe"
      "rtoaddordetractTgheworldadswfilllittlenotlenorlongrememberwhatwesayhereb"
      "utitcanneverforgetwhattheydidhereItisforusthelivingrathertobededicatedhe"
      "retotheulnfinishedworkwhichtheywhofoughtherehavethusfarsonoblyadvancedIt"
      "isratherforustobeherededicatedtothegreattdafskremainingbeforeusthatfromt"
      "hesehonoreddeadwetakeincreaseddevotiontothatcauseforwhichtheygavethelast"
      "pfullmeasureofdevotionthatweherehighlyresolvethatthesedeadshallnothavedi"
      "edinvainthatthisnationunsderGodshallhaveanewbirthoffreedomandthatgovernm"
      "entofthepeoplebythepeopleforthepeopleshallnotperishfromtheearth");
}
int main() {
  solve();
  cout << "program has finished..." << endl;
  return 0;
}
