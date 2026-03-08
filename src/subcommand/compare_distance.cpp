/**
 * \file compare_distance.cpp: Simulate seeds along a path, build a zipcode tree, and print the real distance along the path and the distance from zipcodes to stdout
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

#include <bdsg/overlays/overlay_helper.hpp>

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
    << "test distances found by zipcode trees by simulating reads and seeds along a path in the graph. Writes tsv of \"real_distance\tzipcode_distance\" to stdout" << endl
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

    // Get a list of paths to include in the path position overlay
    std::unordered_set<std::string> paths_set;
    
    // go through all paths in the pangenome and save them
    path_handle_graph->for_each_path_matching(nullptr, nullptr, nullptr, [&] (handlegraph::path_handle_t path) {
        paths_set.emplace(path_handle_graph->get_path_name(path));
        return true;
    });

    bdsg::PathPositionOverlayHelper overlay_helper;
    PathPositionHandleGraph* graph = overlay_helper.apply(path_handle_graph.get(), paths_set);

    unique_ptr<SnarlDistanceIndex> distance_index = vg::io::VPKG::load_one<SnarlDistanceIndex>(distance_name);
    distance_index->preload(true);


    // Get all paths
    std::vector<path_handle_t> paths;
    graph->for_each_path_matching(nullptr, nullptr, nullptr, [&] (handlegraph::path_handle_t path_handle) {
        paths.emplace_back(path_handle);
        return true;
    });

    // For now, just use 1500 as the length of the read
    size_t read_length = 1500;

    // Prepare random number generators for read start and read length 
    //Copied from https://en.cppreference.com/w/cpp/numeric/random/uniform_int_distribution.html
    std::random_device rd;  // a seed source for the random number engine
    std::mt19937 gen(rd()); // mersenne_twister_engine seeded with rd()
    std::uniform_int_distribution<> path_distr(0, paths.size()-1);
    // Rough distribution of distances between seeds from real hifi reads
    std::normal_distribution<> seed_gap_distr{130, 123};

    std::cout << "truth_distance\tziptree_distance\tdiff" << endl;
    #pragma omp parallel for
    for (size_t i = 0 ; i < read_count ; i++) {


        const path_handle_t& path = paths.at(path_distr(gen));
        size_t path_length = graph->get_path_length(path);
        std::uniform_int_distribution<> read_start_distr(0, path_length);
        size_t read_start_offset = read_start_distr(gen);

        std::vector<SnarlDistanceIndexClusterer::Seed> seeds;
        std::vector<fake_minimizer_t> minimizers;
        std::vector<vg::algorithms::Anchor> anchors;

        size_t seed_offset = read_start_offset;
        while (seed_offset < read_start_offset + read_length && seed_offset < path_length) {

            step_handle_t step = graph->get_step_at_position(path, seed_offset);
            handle_t handle = graph->get_handle_of_step(step);

            // Get the offset of the start of the node on the path
            size_t node_start_offset = graph->get_position_of_step(step); 

            assert(node_start_offset <= seed_offset);
            assert((seed_offset - node_start_offset) < graph->get_length(handle));

            pos_t pos = make_pos_t(graph->get_id(handle), graph->get_is_reverse(handle), seed_offset - node_start_offset);

            // Make the zipcode
            ZipCode zipcode;
            zipcode.fill_in_zipcode(*distance_index, pos);

            //Make the seed
            seeds.emplace_back(pos, minimizers.size(), zipcode); 

            //Make the minimizer
            fake_minimizer_t minimizer;
            minimizer.value.offset = seed_offset;
            minimizer.value.is_reverse = false;
            minimizers.emplace_back(std::move(minimizer));

            anchors.emplace_back(seed_offset, pos, 1, 10, 10, 10, seeds.size()-1);

            // Get the next start of a seed
            seed_offset += seed_gap_distr(gen); 
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


