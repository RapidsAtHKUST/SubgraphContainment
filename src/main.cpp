#include <iostream>
#include <fstream>
#include "FSCCommand.h"
#include "GraphDB.h"
#include "FSC.h"
using namespace std;

int main(int argc, char** argv) {
    FSCCommand command(argc, argv);
    const string data_graph_path = command.GetDataGraphFilePath();
    const string query_graph_path = command.GetQueryGraphFilePath();
    const GraphFormat data_graph_format = command.GetDataGraphFormat();
    const GraphFormat query_graph_format = command.GetQueryGraphFormat();

    GraphDB graphDB;
    graphDB.BuildDB(data_graph_path, data_graph_format);
    vector<Graph*>* data_graphs = graphDB.GetDataGraphs();

    vector<Graph*>* query_graphs = graphDB.LoadQueryGraphs(query_graph_path, query_graph_format);

    std::cout << "The number of data graphs: " << data_graphs->size() << " ." << endl;
    std::cout << "The number of query graphs: " << query_graphs->size() << " ." << endl;
    FSC fsc(&graphDB);
    int* results = new int[data_graphs->size()];
    int count = 0;
    // Process the query graphs.
    for (size_t i = 0; i < query_graphs->size(); ++i) {
        cout << "-------------------------------------------------------------------" << endl;
        cout << "Process Query Graph: " << i << " ." << endl;
        fsc.Query((*query_graphs)[i], results, count);
    }
    cout << "-------------------------------------------------------------------" << endl;
    cout << "The average number of candidates : "  << (double)fsc.total_candidate_num_ / query_graphs->size() << " ." << endl;
    cout << "The average number of answers : " << (double)fsc.total_answer_num_ / query_graphs->size() << " ." << endl;
    cout << "The average processing time: " << fsc.total_time_ / query_graphs->size() << " us." << endl;
    cout << "The average filtering time: " << fsc.total_filtering_time_ / query_graphs->size() << " us." << endl;
    cout << "The average ordering time: " << fsc.total_ordering_time_ / query_graphs->size() << " us." << endl;
    cout << "The average enumeration time: " << fsc.total_enumeration_time_ / query_graphs->size() << " us." << endl;
    cout << "The average verification time: " << fsc.total_verification_time_ / query_graphs->size() << " us." << endl;
    cout << "The average per verification time: " << fsc.total_per_verification_time_ / query_graphs->size() << " us." << endl;
    cout << "The average verification time due to false prediction: " << fsc.total_save_time_ / query_graphs->size() << " us." << endl;
    cout << "The peak memory cost: " << fsc.peak_memory_cost_ << " MB." << endl;
    cout << "The average filtering precision: " << fsc.total_filtering_precision / query_graphs->size() << " ." << endl;
    delete[] results;
    return 0;
}