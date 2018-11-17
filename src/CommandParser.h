//
// Created by ssunah on 11/28/16.
//

#ifndef FSC_COMMANDPARSER_H
#define FSC_COMMANDPARSER_H

#include <string>
#include <algorithm>
using namespace std;

class CommandParser {
private:
    vector<string> tokens_;

public:
    CommandParser(const int &argc, char **argv) {
        for (int i = 1; i < argc; ++i)
            tokens_.push_back(string(argv[i]));
    }

    const string GetCommandOption(const string &option) const {

        vector<string>::const_iterator itr;
        itr = find(tokens_.begin(), tokens_.end(), option);
        if (itr != tokens_.end() && ++itr != tokens_.end()) {
            return *itr;
        }
        return "";
    }

    bool CommandOptionExists(const string &option) const {
        return find(tokens_.begin(), tokens_.end(), option) != tokens_.end();
    }
};

#endif //FSC_COMMANDPARSER_H
