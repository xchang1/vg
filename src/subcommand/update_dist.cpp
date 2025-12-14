// update_dist.cpp: define the "vg update_dist" subcommand, which takes an old distance index and updates it to the current version
#include "subcommand.hpp"

#include "../vg.hpp"
#include "xg.hpp"
#include <vg/io/stream.hpp>
#include <vg/io/vpkg.hpp>
#include "../io/save_handle_graph.hpp"
#include "../snarl_distance_index.hpp"
#include "../gbwtgraph_helper.hpp"

using namespace std;
using namespace vg;
using namespace vg::subcommand;


void help_update_dist(char** argv) {
    cerr << "Re-write a v2 distance index to be a v3 distance index" << endl 
         << "usage: " << argv[0] << " update_dist [graph.old.dist] [graph.new.dist]" << endl;
}

int main_update_dist(int argc, char** argv) {

    if (argc == 2) {
        help_update_dist(argv);
        return 1;
    }

    std::string in_dist_name = argv[3];
    std::string out_dist_name = argv[4];

    cerr << "Re-writing old distance index " << in_dist_name << " to new file " << out_dist_name << endl;


    // Load old distance index
    ifstream instream;
    instream.open(in_dist_name);

    SnarlDistanceIndex distance_index;

    // Load the contents of the distance index as a vector, using the same prefix
    bdsg::yomo::UniqueMappedPointer<bdsg::MappedIntVector> saved_vector;
    saved_vector.construct(distance_index.get_prefix());
    saved_vector.load(instream, distance_index.get_prefix());
    instream.close();

    // The only difference is that the new vector must have the version as its second item
    bdsg::yomo::UniqueMappedPointer<bdsg::MappedIntVector> new_vector;
    new_vector.construct(distance_index.get_prefix());
    new_vector->width(saved_vector->width());
    new_vector->resize(saved_vector->size()+1);
    size_t old_first = saved_vector->at(std::size_t(0));
    new_vector->at(std::size_t(0)) = old_first;
    // COpy these from the distance index since they're private
    size_t version_number_sentinel = (1 << 10) - 1;
    size_t current_version_number = 3;
    new_vector->at(std::size_t(1)) = version_number_sentinel ^ current_version_number;  
    for (size_t i = 1 ; i < saved_vector->size( ); i++) {
        size_t old_val = saved_vector->at(i);
        new_vector->at(i+1) = old_val;
    }

    // Save the distance index in the new format, with the prefix
    ofstream outstream;
    outstream.open(out_dist_name);
    new_vector.save(outstream);
    outstream.close();


    return 0;
}

// Register subcommand
static Subcommand vg_update_dist("update_dist", "update distance index to new version", main_update_dist);
