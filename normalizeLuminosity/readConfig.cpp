// taken from https://github.com/nlohmann/json
#include "jsonLibrary.hpp"
#include <bits/stdc++.h>

using json = nlohmann::json;
using namespace std;

json readConfigFile(string configFilePath)
{
    cout << "-----------------------------------" << endl;
    cout << "Reading configuration file..." << endl;

    ifstream ifs(configFilePath);
    string content((istreambuf_iterator<char>(ifs)),
                   (istreambuf_iterator<char>()));

    cout << "Parsing to json..." << endl;
    json parsedStuff = json::parse(content);

    cout << "config file read: " << parsedStuff.dump() << endl;
    cout << "-----------------------------------" << endl
         << endl;

    return parsedStuff;
}

int main(int argc, char const *argv[])
{

    json parsedStuff = readConfigFile("config_xsec.json");

    for (json::iterator it = parsedStuff.begin(); it != parsedStuff.end(); ++it)
    {
        cout << it.key() << " : " << it.value() << "\n";
    }

    // cout << parsedStuff.dump() << endl;
    // cout << parsedStuff["zz"] << endl;
    // cout << parsedStuff["zz"]["xsec"] << endl;
    // cout << parsedStuff["zz"]["error"] << endl;

    float a = parsedStuff["xsections"]["zz"]["xsec"];
    cout << a*2 <<endl;

    return 0;
}
