#include <iostream>
#include <map>
#include <unordered_map>
#include <string>

using namespace std;

int main()
{
    map<string,unsigned long> employer;
    unordered_map<string,unsigned long> salary;
    //map<string,unsigned long> salary;
    string keyfeature;
    employer["Celine Dion"] = 1;
    employer["Whitney Houston"] = 1;
    
    for(auto e: employer)
        cout << "name: " << e.first
        << "\t id: "  << e.second << endl;
    
    unsigned total_payroll = 0;
    
    salary["-3AGG"] = 1;
    salary["+4AACC"] = 1;
    salary["-2GG"] = 1;
    keyfeature = "+4AACC";
    auto it = salary.find(keyfeature);
    if(it != salary.end())
        it->second = salary[keyfeature]+1 ;
    else
        salary[keyfeature] = 1;
        
    
    for(auto s: salary)
        cout  << s.first << ":"<<s.second << endl;
        //total_payroll += s.second;
    
    cout << "total counts " << total_payroll << endl;
    
    
    if (salary.count("+4AACC")>0)
        std::cout << "mymap has " << "+4AACC " << salary["+4AACC"]<< std::endl;
    else
        std::cout << "mymap has no " << "+4AACC" << std::endl;
    
    
    

    
    
    return 0;
}