//
// Created by ssunah on 11/28/16.
//

#ifndef FSC_FSMCOMMAND_H
#define FSC_FSMCOMMAND_H

#include "CommandParser.h"
#include "Config.h"
#include "Types.h"
#include <map>
using namespace std;

enum OptionKeyword {
    QueryGraphFile = 1,     // -q, The query graph file path, compulsive parameter
    DataGraphFile = 2,      // -d, The data graph file path, compulsive parameter
    QueryGraphFormat = 3,   // -qf, The query graph file format, compulsive parameter
    DataGraphFormat = 4     // -df, The data graph file format, compulsive parameter
};

class FSCCommand : public CommandParser {
private:
    map<OptionKeyword, string> options_key;
    map<OptionKeyword, string> options_value;
    map<string, GraphFormat> format_key;
    map<OptionKeyword, GraphFormat > format_value;
private:
    void ProcessOptions() {

        for (map<OptionKeyword, string>::iterator iter = options_key.begin(); iter != options_key.end(); iter++) {
            if (!CommandOptionExists(iter->second)) {
                printf("Please specify %s option.\n", iter->second.c_str());
                exit(-1);
            }
        }

        // Query graph file path
        string value = GetCommandOption(options_key[OptionKeyword::QueryGraphFile]);
        if (value == "") {
            printf("The query graph file path can not be empty.\n");
            exit(-1);
        }
        options_value[OptionKeyword::QueryGraphFile] = value;

        // Data graph file path
        value = GetCommandOption(options_key[OptionKeyword::DataGraphFile]);
        if (value == "") {
            printf("The data graph file path can not be empty.\n");
            exit(-1);
        }
        options_value[OptionKeyword::DataGraphFile] = value;

        // Query graph format
        value = GetCommandOption(options_key[OptionKeyword::QueryGraphFormat]);
        if (format_key.find(value) == format_key.end()) {
            printf("The query graph format cannot be supported.\n");
            exit(-1);
        }
        format_value[OptionKeyword::QueryGraphFormat] = format_key[value];

        // Data graph format
        value = GetCommandOption(options_key[OptionKeyword::DataGraphFormat]);
        if (format_key.find(value) == format_key.end()) {
            printf("The data graph format cannot be supported.\n");
            exit(-1);
        }
        format_value[OptionKeyword::DataGraphFormat] = format_key[value];
    }

public:
    FSCCommand(const int &argc, char **argv) : CommandParser(argc, argv) {
        // Initialize options value
        options_key[OptionKeyword::QueryGraphFile] = "-q";
        options_key[OptionKeyword::DataGraphFile] = "-d";
        options_key[OptionKeyword::QueryGraphFormat] = "-qf";
        options_key[OptionKeyword::DataGraphFormat] = "-df";
        format_key["FA"] = GraphFormat::FA;
        format_key["FB"] = GraphFormat::FB;
        format_key["FC"] = GraphFormat::FC;
        ProcessOptions();
    }

    const string& GetDataGraphFilePath() {
        return options_value[OptionKeyword::DataGraphFile];
    }

    const string& GetQueryGraphFilePath() {
        return options_value[OptionKeyword::QueryGraphFile];
    }

    const GraphFormat GetQueryGraphFormat() {
        return format_value[OptionKeyword::QueryGraphFormat];
    }

    const GraphFormat GetDataGraphFormat() {
        return format_value[OptionKeyword::DataGraphFormat];
    }
};


#endif //FSC_FSMCOMMAND_H
