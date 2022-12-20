#include <iostream>
#include <iterator>
#include <string>
#include <regex>
using namespace std;

// Word may contain alphabet or integer, anything else is treated as delimeter
vector<string> split_string_into_words(string s){
	vector<string> ans;
    regex word_regex("(\\w+)");
    auto words_begin = sregex_iterator(s.begin(), s.end(), word_regex);
    auto words_end = sregex_iterator();

    for (sregex_iterator i = words_begin; i != words_end; ++i) {
        smatch match = *i;
        string match_str = match.str();
        ans.push_back(match_str);
    }
    return ans;
}

int main()
{
	string s="{0,1} {2,3} {3,4} {3,5}";
	// string s="n=6,m=4m00";
	// string s="ishan d";
    vector<string>words=split_string_into_words(s);
    for(auto word : words)
        cout << word<<"\n";
}
