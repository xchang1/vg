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
    
    // Make a Mapper to look up MEM seeds
    unique_ptr<Mapper> mapper;
    if (gcsa_index) {
        // We will find MEMs using a Mapper
        mapper = make_unique<Mapper>(xg_index.get(), gcsa_index.get(), lcp_index.get());
    }
    // Otherwise we will find minimizers using the minimizer_index
    
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
            
            // We will find all the seed hits
            vector<pos_t> seeds;
            
            // If working with MEMs, this will hold all the MEMs
            vector<MaximalExactMatch> mems;
            // If working with minimizers, this will hold all the minimizers in the query
            vector<gbwtgraph::DefaultMinimizerIndex::minimizer_type> minimizers;
            // And either way this will map from seed to MEM or minimizer that generated it
            vector<size_t> seed_to_source;
            
            if (mapper) {
                // Find MEMs
                double lcp_avg, fraction_filtered;
                mems = mapper->find_mems_deep(aln.sequence().begin(), aln.sequence().end(), lcp_avg, fraction_filtered);
                
                // Convert to position seeds
                for (size_t i = 0; i < mems.size(); i++) {
                    auto& mem = mems[i];
                    for (gcsa::node_type n : mem.nodes) {
                        // Convert from GCSA node_type packing to a pos_t
                        seeds.push_back(make_pos_t(n));
                        // And remember which MEM the seed came from.
                        seed_to_source.push_back(i);
                    }
                }
            } else {
                // Find minimizers
                assert(minimizer_index);
                
                // Find minimizers in the query
                minimizers = minimizer_index->minimizers(aln.sequence());
                
                for (size_t i = 0; i < minimizers.size(); i++) {
                    // For each minimizer
                    if (hit_cap != 0 && minimizer_index->count(minimizers[i]) <= hit_cap) {
                        // The minimizer is infrequent enough to be informative, so feed it into clustering
                        
                        // Locate it in the graph. We do not have to reverse the hits for a
                        // reverse minimizers, as the clusterer only cares about node ids.
                        for (auto& hit : minimizer_index->find(minimizers[i])) {
                            // For each position, remember it and what minimizer it came from
                            seeds.push_back(hit);
                            seed_to_source.push_back(i);
                        }
                    }
                }
                
            }

            vector<size_t> correct_seed_indices;
            
            //Remove incorrect seeds
            if (aln.refpos_size() != 0) {
                // Take the first refpos as the true position.
                auto& true_pos = aln.refpos(0);
            
                for (size_t i = 0; i < seeds.size(); i++) {
                    // Find every seed's reference positions. This maps from path name to pairs of offset and orientation.
                    auto offsets = algorithms::nearest_offsets_in_paths(&(*xg_index), seeds[i], 100);
                    for (auto& hit_pos : offsets[xg_index->get_path_handle(true_pos.name())]) {
                        // Look at all the ones on the path the read's true position is on.
                        if (abs((int64_t)hit_pos.first - (int64_t) true_pos.offset()) < 200) {
                            // Call this seed hit close enough to be correct
                            correct_seed_indices.push_back(i);
                        }
                    }
                }
            }


            //Get distances between seeds
            for (size_t i1 = 0 ; i1 < correct_seed_indices.size() ; i1++) {
                pos_t pos1 = seeds[correct_seed_indices[i1]];
                for (size_t i2 = 0 ; i2 < correct_seed_indices.size() ; i2++ ) {
                    pos_t pos2 = seeds [correct_seed_indices[i2]];
                    if (minimizers[seed_to_source[correct_seed_indices[i2]]].offset <= minimizers[seed_to_source[correct_seed_indices[i1]]].offset) {
                        continue;
                    }
                    int64_t read_dist = minimizers[seed_to_source[correct_seed_indices[i2]]].offset - minimizers[seed_to_source[correct_seed_indices[i1]]].offset;

                    //Find min distance
                    std::chrono::time_point<std::chrono::system_clock> start1 = std::chrono::system_clock::now();
                    int64_t minDist = distance_index->minDistance(pos1, pos2);
                    std::chrono::time_point<std::chrono::system_clock> end1 = std::chrono::system_clock::now();
                    std::chrono::duration<double> min_time = end1-start1;
                    
                    //Find max distance
                    std::chrono::time_point<std::chrono::system_clock> start2 = std::chrono::system_clock::now();
                    int64_t maxDist = distance_index->maxDistance(pos1, pos2);
                    std::chrono::time_point<std::chrono::system_clock> end2 = std::chrono::system_clock::now();
                    std::chrono::duration<double> max_time = end2-start2;
                    
                    //Find old distance
                    std::chrono::time_point<std::chrono::system_clock> start3 = std::chrono::system_clock::now();
                    int64_t oldDist =abs(xg->min_approx_path_distance(
                                                get_id(pos1), get_id(pos2)));
                    std::chrono::time_point<std::chrono::system_clock> end3 = std::chrono::system_clock::now();
                    std::chrono::duration<double> old_time = end3-start3;
                    
                    
                    //Find max distance
                    std::chrono::time_point<std::chrono::system_clock> start4 = std::chrono::system_clock::now();
                    int64_t tvsDist = tvs.tv_path_length(pos1, pos2, read_dist, 20);
                    std::chrono::time_point<std::chrono::system_clock> end4 = std::chrono::system_clock::now();
                    std::chrono::duration<double> tvs_time = end4-start4;
                    
                    

                    if (minDist != -1) {
                    
                        cerr << read_dist << "\t" << minDist << "\t" << maxDist << "\t" << oldDist << "\t" << tvsDist << "\t" << min_time.count() << "\t" << max_time.count() << "\t" << old_time.count() << "\t" << tvs_time.count() << endl;

                    }
                }
            }

                    
            // Cluster the seeds. Get sets of input seed indexes that go together.
            // Make sure to time it.
            std::chrono::time_point<std::chrono::system_clock> start = std::chrono::system_clock::now();
            tuple<vector<vector<size_t>>,vector<vector<size_t>>> paired_clusters = clusterer.cluster_seeds(seeds, distance_limit);
            vector<vector<size_t>> clusters = std::move(std::get<0>(paired_clusters));
            std::chrono::time_point<std::chrono::system_clock> end = std::chrono::system_clock::now();
            std::chrono::duration<double> elapsed_seconds = end-start;
            cout << elapsed_seconds.count() << endl;
           
        });
    });
    
    return 0;
}

// Register subcommand
static Subcommand vg_cluster("plot", "get stats for plot", DEVELOPMENT, main_plot);


