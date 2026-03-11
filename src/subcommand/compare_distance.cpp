/**
 * \file compare_distance.cpp: Simulate seeds along a random walk, build a zipcode tree, and print the real distance along the walk and the distance from zipcodes to stdout
 */

#include <omp.h>
#include <unistd.h>
#include <getopt.h>
#include <iostream>
#include <cassert>
#include <vector>
#include <random>

#include "subcommand.hpp"

#include "../zip_code.hpp"
#include "../zip_code_tree.hpp"
#include "../snarl_seed_clusterer.hpp"
#include "../algorithms/chain_items.hpp"
#include <vg/io/vpkg.hpp>
#include <vg/io/stream.hpp>

//#define USE_CALLGRIND

#ifdef USE_CALLGRIND
#include <valgrind/callgrind.h>
#endif

using namespace std;
using namespace vg;
using namespace vg::subcommand;


void help_testzip(char** argv) {
    cerr
    << "usage: " << argv[0] << " testzip -x [graph] -d [dist] > distances.tsv" << endl 
    << "test distances found by zipcode trees by simulating reads and seeds along a random walk in the graph. Writes tsv of \"real_distance\tzipcode_distance\" to stdout" << endl
    << endl
    << "basic options:" << endl
    << "  -h, --help                    print this help message to stderr and exit" << endl
    << "  -x, --xg-name FILE            use this xg index or graph (required)" << endl
    << "  -d, --dist-name FILE          use this distance index (required)" << endl
    << "  -c, --read-count INT          simulate this many reads [1000]" << endl
    << "computational parameters:" << endl
    << "  -t, --threads INT             number of compute threads to use" << endl;
}

int main_testzip(int argc, char** argv) {

    if (argc == 2) {
        help_testzip(argv);
        return 1;
    }

    // initialize parameters with their default options
    string xg_name;
    string distance_name;
    size_t read_count = 1000;
    
    int c;
    optind = 2; // force optind past command positional argument
    while (true) {
        static struct option long_options[] =
        {
            {"help", no_argument, 0, 'h'},
            {"xg-name", required_argument, 0, 'x'},
            {"dist-name", required_argument, 0, 'd'},
            {"read-count", required_argument, 0, 'c'},
            {"threads", required_argument, 0, 't'},
            {0, 0, 0, 0}
        };

        int option_index = 0;
        c = getopt_long (argc, argv, "h?x:d:c:t:",
                         long_options, &option_index);


        // Detect the end of the options.
        if (c == -1)
            break;

        switch (c)
        {
            case 'x':
                xg_name = optarg;
                if (xg_name.empty()) {
                    cerr << "error:[vg testzip] Must provide XG file with -x." << endl;
                    exit(1);
                }
                break;
                               
            case 'd':
                distance_name = optarg;
                if (distance_name.empty()) {
                    cerr << "error:[vg testzip] Must provide distance index file with -d." << endl;
                    exit(1);
                }
                break;
            
            case 'c':
                read_count = parse<size_t>(optarg);
                break;
                
            case 't':
            {
                int num_threads = parse<int>(optarg);
                if (num_threads <= 0) {
                    cerr << "error:[vg testzip] Thread count (-t) set to " << num_threads 
                         << ", must set to a positive integer." << endl;
                    exit(1);
                }
                omp_set_num_threads(num_threads);
            }
                break;
                
            case 'h':
            case '?':
            default:
                help_testzip(argv);
                exit(1);
                break;
        }
    }
    
    
    if (xg_name.empty()) {
        cerr << "error:[vg testzip] Testing the zip code tree distances requires an XG index, must provide XG file (-x)" << endl;
        exit(1);
    }
    
    if (distance_name.empty()) {
        cerr << "error:[vg testzip] Testing the zip code tree distances requires a distance index, must provide distance index file (-d)" << endl;
        exit(1);
    }
    
    // create in-memory objects
    unique_ptr<PathHandleGraph> path_handle_graph = vg::io::VPKG::load_one<PathHandleGraph>(xg_name);
    PathHandleGraph* graph = path_handle_graph.get();

    unique_ptr<SnarlDistanceIndex> distance_index = vg::io::VPKG::load_one<SnarlDistanceIndex>(distance_name);
    distance_index->preload(true);

    // For now, just use 1500 as the length of the read
    size_t read_length = 1500;

    // Prepare random number generators for read start and read length 
    //Copied from https://en.cppreference.com/w/cpp/numeric/random/uniform_int_distribution.html
    std::random_device rd;  // a seed source for the random number engine
    std::mt19937 gen(rd()); // mersenne_twister_engine seeded with rd()
    std::uniform_int_distribution<> node_id_distr(graph->min_node_id(), graph->max_node_id());
    // Rough distribution of distances between seeds from real hifi reads
    std::normal_distribution<> seed_gap_distr{130, 123};

    std::cout << "truth_distance\tziptree_distance\tdiff" << endl;
    #pragma omp parallel for
    for (size_t i = 0 ; i < read_count ; i++) {

        // Start random walk from a random node
        handlegraph::id_t node_id = node_id_distr(gen);
        if (!graph->has_node(node_id)) {
            continue;
        }
        pos_t current_position = make_pos_t(node_id, false, 0);
        handle_t current_handle = graph->get_handle(node_id, false);


        std::vector<SnarlDistanceIndexClusterer::Seed> seeds;
        std::vector<fake_minimizer_t> minimizers;
        std::vector<vg::algorithms::Anchor> anchors;

        size_t seed_offset_in_path = 0;
        while (seed_offset_in_path < read_length) {

            // Get the next start of a seed
            size_t distance_to_traverse = seed_gap_distr(gen);
            if (distance_to_traverse > (read_length - seed_offset_in_path)) {
                break;
            }

            // Update the position in the "read"
            seed_offset_in_path += distance_to_traverse; 

            // If we hit a tip in the walk, break out of the outer loop without adding a new seed
            bool hit_tip = false;

            // Randomly walk through the graph distance_to_traverse bases
            while (distance_to_traverse > 0) {
                size_t current_node_length = graph->get_length(current_handle);
                assert(get_offset(current_position) <= current_node_length);
                size_t distance_to_end_of_node = current_node_length - get_offset(current_position); 
                if (distance_to_traverse < distance_to_end_of_node) {
                    // If we end the traversal in this node, put the position at the end of the traversal and stop
                    pos_t new_pos = make_pos_t(get_id(current_position), get_is_rev(current_position), get_offset(current_position) + distance_to_traverse);
                    current_position = new_pos;
                    distance_to_traverse = 0;
                } else {
                    // If we keep going, pick a random next node and reset the position to the start of this node


                    // Pick a random edge to follow. To do this, find how many edges there are and pick a random one
                    size_t next_step_count = 0;
                    graph->follow_edges(current_handle, false, [&](const handle_t& next) {
                        next_step_count++;
                        return true;
                    });
                    if (next_step_count == 0) {
                        // If there is nothing left to traverse, break out of the outer loop
                        hit_tip = true;
                        break;
                    }
                    std::uniform_int_distribution<> edge_distr(0, next_step_count-1);
                    size_t next_edge_num = edge_distr(gen);
                    size_t current_edge = 0;

                    bool found_next_node = graph->follow_edges(current_handle, false, [&](const handle_t& next_handle) {
                        if (next_edge_num == current_edge) {
                            // Reserve false for no tips
                            current_position = make_pos_t(graph->get_id(next_handle), graph->get_is_reverse(next_handle), 0);
                            current_handle = next_handle;
                            return true;
                        } else {
                            ++current_edge;
                            return true;
                        }
                    });

                    assert(distance_to_end_of_node <= distance_to_traverse);
                    distance_to_traverse -= distance_to_end_of_node;
                }
            }
            if (hit_tip) {
                break;
            }

            // Make the zipcode
            ZipCode zipcode;
            zipcode.fill_in_zipcode(*distance_index, current_position);

            //Make the seed
            seeds.emplace_back(current_position, minimizers.size(), zipcode); 

            //Make the minimizer
            fake_minimizer_t minimizer;
            minimizer.value.offset = seed_offset_in_path;
            minimizer.value.is_reverse = false;
            minimizers.emplace_back(std::move(minimizer));

            anchors.emplace_back(seed_offset_in_path, current_position, 1, 10, 10, 10, seeds.size()-1);

        }
        // Make the vector view of minimizers
        std::vector<size_t> minimizer_order(minimizers.size(), 0);
        for (size_t i = 0 ; i < minimizer_order.size() ; i++) {
            minimizer_order[i]=i;
        }

        VectorView<fake_minimizer_t> minimizer_vector{minimizers, minimizer_order};

        // Make the vector view of anchors
        std::vector<size_t> anchor_order(anchors.size(), 0);
        for (size_t i = 0 ; i < anchor_order.size() ; i++) {
            anchor_order[i]=i;
        }

        VectorView<vg::algorithms::Anchor> anchor_vector{anchors, anchor_order};


        //Make the zip code trees for these seeds
        ZipCodeForest forest;
        forest.fill_in_forest(seeds, minimizer_vector, *distance_index, std::numeric_limits<size_t>::max());

        for (const ZipCodeTree& ziptree : forest.trees) {

            vg::algorithms::transition_iterator for_each_transition = vg::algorithms::zip_tree_transition_iterator(seeds, ziptree, 3000, std::numeric_limits<size_t>::max());

            for_each_transition(anchor_vector, *distance_index, *graph, 3000, [&](size_t from_anchor, size_t to_anchor, size_t read_distance, size_t graph_distance) {
                #pragma omp critical (cout)
                {
                    std::cout << read_distance << "\t" << graph_distance << "\t" << ((int)read_distance -(int)graph_distance) << std::endl;
                }

            });
        }

    }

    return 0;
}

// Register subcommand
static Subcommand vg_testzip("testzip", "find distances between simulated seeds using zipcode trees", DEVELOPMENT, main_testzip);


