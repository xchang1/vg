/**
 * \file cluster_main.cpp: experimental snarl clustering test harness
 */

#include <omp.h>
#include <unistd.h>
#include <getopt.h>
#include <iostream>
#include <cassert>
#include <vector>
#include <unordered_set>
#include <chrono>

#include "subcommand.hpp"

#include "../seed_clusterer.hpp"
#include "../mapper.hpp"
#include "../annotation.hpp"
#include "../xg.hpp"
#include <vg/io/vpkg.hpp>
#include <vg/io/stream.hpp>
#include <vg/io/protobuf_emitter.hpp>
#include "cluster.hpp"
#include "min_distance.hpp"
#include "../algorithms/dijkstra.hpp"


#include <gbwtgraph/minimizer.h>

//#define USE_CALLGRIND

#ifdef USE_CALLGRIND
#include <valgrind/callgrind.h>
#endif

using namespace std;
using namespace vg;
using namespace vg::subcommand;

void help_plot(char** argv) {
    cerr
    << "usage: " << argv[0] << " cluster [options] input.gam > output.gam" << endl
    << "Find and cluster mapping seeds." << endl
    << endl
    << "basic options:" << endl
    << "  -x, --xg-name FILE            use this xg index (required)" << endl
    << "  -g, --gcsa-name FILE          use this GCSA2/LCP index pair (both FILE and FILE.lcp)" << endl
    << "  -m, --minimizer-name FILE     use this minimizer index" << endl
    << "  -d, --dist-name FILE          cluster using this distance index (required)" << endl
    << "  -c, --hit-cap INT             ignore minimizers with more than this many locations [10]" << endl
    << "computational parameters:" << endl
    << "  -t, --threads INT             number of compute threads to use" << endl;
}

int main_plot(int argc, char** argv) {

    if (argc == 2) {
        help_plot(argv);
        return 1;
    }

    // initialize parameters with their default options
    string xg_name;
    string gcsa_name;
    string minimizer_name;
    string distance_name;
    // How close should two hits be to be in the same cluster?
    size_t distance_limit = 1000;
    size_t hit_cap = 10;
    
    int c;
    optind = 2; // force optind past command positional argument
    while (true) {
        static struct option long_options[] =
        {
            {"help", no_argument, 0, 'h'},
            {"xg-name", required_argument, 0, 'x'},
            {"gcsa-name", required_argument, 0, 'g'},
            {"minimizer-name", required_argument, 0, 'm'},
            {"dist-name", required_argument, 0, 'd'},
            {"hit-cap", required_argument, 0, 'c'},
            {"threads", required_argument, 0, 't'},
            {0, 0, 0, 0}
        };

        int option_index = 0;
        c = getopt_long (argc, argv, "hx:g:m:d:c:t:",
                         long_options, &option_index);


        // Detect the end of the options.
        if (c == -1)
            break;

        switch (c)
        {
            case 'x':
                xg_name = optarg;
                if (xg_name.empty()) {
                    cerr << "error:[vg cluster] Must provide XG file with -x." << endl;
                    exit(1);
                }
                break;
                
            case 'g':
                gcsa_name = optarg;
                if (gcsa_name.empty()) {
                    cerr << "error:[vg cluster] Must provide GCSA file with -g." << endl;
                    exit(1);
                }
                break;
            
            case 'm':
                minimizer_name = optarg;
                if (minimizer_name.empty()) {
                    cerr << "error:[vg cluster] Must provide minimizer file with -m." << endl;
                    exit(1);
                }
                break;
                
                
            case 'd':
                distance_name = optarg;
                if (distance_name.empty()) {
                    cerr << "error:[vg cluster] Must provide distance index file with -d." << endl;
                    exit(1);
                }
                break;
            
            case 'c':
                hit_cap = parse<size_t>(optarg);
                break;
                
            case 't':
            {
                int num_threads = parse<int>(optarg);
                if (num_threads <= 0) {
                    cerr << "error:[vg cluster] Thread count (-t) set to " << num_threads << ", must set to a positive integer." << endl;
                    exit(1);
                }
                omp_set_num_threads(num_threads);
            }
                break;
                
            case 'h':
            case '?':
            default:
                help_plot(argv);
                exit(1);
                break;
        }
    }
    
    
    if (xg_name.empty()) {
        cerr << "error:[vg cluster] Finding clusters requires an XG index, must provide XG file (-x)" << endl;
        exit(1);
    }
    
    if (gcsa_name.empty() && minimizer_name.empty()) {
        cerr << "error:[vg cluster] Finding clusters requires a GCSA2 index or minimizer index (-g, -m)" << endl;
        exit(1);
    }
    
    
    if (distance_name.empty()) {
        cerr << "error:[vg cluster] Finding clusters requires a distance index, must provide distance index file (-d)" << endl;
        exit(1);
    }
    
    // create in-memory objects
    unique_ptr<PathPositionHandleGraph> xg_index = vg::io::VPKG::load_one<PathPositionHandleGraph>(xg_name);
    unique_ptr<xg::XG> xg = vg::io::VPKG::load_one<xg::XG>(xg_name);
    unique_ptr<gcsa::GCSA> gcsa_index;
    unique_ptr<gcsa::LCPArray> lcp_index;
    if (!gcsa_name.empty()) {
        gcsa_index = vg::io::VPKG::load_one<gcsa::GCSA>(gcsa_name);
        lcp_index = vg::io::VPKG::load_one<gcsa::LCPArray>(gcsa_name + ".lcp");
    }
    unique_ptr<gbwtgraph::DefaultMinimizerIndex> minimizer_index;
    if (!minimizer_name.empty()) {
        minimizer_index = vg::io::VPKG::load_one<gbwtgraph::DefaultMinimizerIndex>(minimizer_name);
    }
    unique_ptr<MinimumDistanceIndex> distance_index = vg::io::VPKG::load_one<MinimumDistanceIndex>(distance_name);
    TargetValueSearch tvs(*xg_index, new TipAnchoredMaxDistance(*distance_index),
                              new SnarlMinDistance(*distance_index));
    random_device seed_source;
    default_random_engine generator(seed_source());

    
    // Make the clusterer
    SnarlSeedClusterer clusterer(*distance_index);
    //distance_index->printSnarlStats();
    
    // Make a Mapper to look up MEM seeds
    unique_ptr<Mapper> mapper;
    if (gcsa_index) {
        // We will find MEMs using a Mapper
        mapper = make_unique<Mapper>(xg_index.get(), gcsa_index.get(), lcp_index.get());
    }
    // Otherwise we will find minimizers using the minimizer_index
    
    vector<pos_t> all_seeds;
    get_input_file(optind, argc, argv, [&](istream& in) {
        // Open up the input GAM
        
        // Make the output emitter
        vg::io::ProtobufEmitter<Alignment> emitter(cout);
        
#ifdef USE_CALLGRIND
        // We want to profile the clustering and the code around it.
        CALLGRIND_START_INSTRUMENTATION;
#endif
        
        vg::io::for_each_parallel<Alignment>(in, [&](Alignment& aln) {
            // For each input alignment
             
           size_t i1 = uniform_int_distribution<int>(0,147)(generator);
           size_t i2 = uniform_int_distribution<int>(0,147)(generator);


           if (i1 + i2 >= 148) {
               if (i2 > i1) {
                   i2 = 148 - i2;
               } else {
                   i1 = 148 - i2;
               }
           }

           vector<Alignment> alns = alignment_ends(aln, i1, i2);


           Position p1 = alignment_end(alns[0]);
           Position p2 = alignment_start(alns[1]);

           
           pos_t pos1 = make_pos_t(p1);
           pos_t pos2 = make_pos_t(p2);
           if (get_id(pos1) != 0 && get_id(pos2) != 0) {

           //Find min distance
           std::chrono::time_point<std::chrono::system_clock> start1 = std::chrono::system_clock::now();
           int64_t minDist = distance_index->minDistance(pos1, pos2);
           std::chrono::time_point<std::chrono::system_clock> end1 = std::chrono::system_clock::now();
           std::chrono::duration<double> min_time = end1-start1;
           

            //Dijkstra distance
            std::chrono::time_point<std::chrono::system_clock> start2 = std::chrono::system_clock::now();
            algorithms::dijkstra(&(*xg), xg_index->get_handle(get_id(pos1), is_rev(pos1)), 
                   [&](const handle_t& curr, size_t d) {
                       if (curr == xg_index->get_handle(get_id(pos2), is_rev(pos2))) {
                           return false;
                       } else {
                           return true;
                       }});
            std::chrono::time_point<std::chrono::system_clock> end2 = std::chrono::system_clock::now();
            std::chrono::duration<double> dtra_time = end2-start2;
           
           //Find old distance
           std::chrono::time_point<std::chrono::system_clock> start3 = std::chrono::system_clock::now();
           int64_t pathDist =abs(xg->min_approx_path_distance(
                                       get_id(pos1), get_id(pos2)));
           std::chrono::time_point<std::chrono::system_clock> end3 = std::chrono::system_clock::now();
           std::chrono::duration<double> path_time = end3-start3;
           
           
           
           
           
           cerr << min_time.count() << "\t" << dtra_time.count() << "\t" << path_time.count() << endl;
           }
           
        });
    });
        
   
    return 0;
}

// Register subcommand
static Subcommand vg_cluster("plot", "get stats for plot", DEVELOPMENT, main_plot);


