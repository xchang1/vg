/**
 * \file gaffe_main.cpp: GFA (Graph Alignment Format) Fast Emitter: a new mapper that will be *extremely* fast once we actually write it
 */

#include <omp.h>
#include <unistd.h>
#include <getopt.h>
#include <iostream>
#include <cassert>
#include <vector>
#include <unordered_set>
#include <chrono>
#include <mutex>

#include "subcommand.hpp"

#include "../seed_clusterer.hpp"
#include "../mapper.hpp"
#include "../annotation.hpp"
#include "../minimizer.hpp"
#include <vg/io/vpkg.hpp>
#include <vg/io/stream.hpp>
#include "../alignment_emitter.hpp"
#include "../gapless_extender.hpp"
#include "../minimizer_mapper.hpp"

//#define USE_CALLGRIND

#ifdef USE_CALLGRIND
#include <valgrind/callgrind.h>
#endif

using namespace std;
using namespace vg;
using namespace vg::subcommand;

void help_giraffe(char** argv) {
    cerr
    << "usage: " << argv[0] << " gaffe [options] > output.gam" << endl
    << "Map unpaired reads using minimizers and gapless extension." << endl
    << endl
    << "basic options:" << endl
    << "  -x, --xg-name FILE            use this xg index (required)" << endl
    << "  -H, --gbwt-name FILE          use this GBWT index (required)" << endl
    << "  -m, --minimizer-name FILE     use this minimizer index (required)" << endl
    << "  -d, --dist-name FILE          cluster using this distance index (required)" << endl
    << "  -p, --progress                show progress" << endl
    << "input options:" << endl
    << "  -G, --gam-in FILE             read and realign GAM-format reads from FILE (may repeat)" << endl
    << "  -f, --fastq-in FILE           read and align FASTQ-format reads from FILE (may repeat)" << endl
    << "output options:" << endl
    << "  -M, --max-multimaps INT       produce up to INT alignments for each read [1]"
    << "  -N, --sample NAME             add this sample name" << endl
    << "  -R, --read-group NAME         add this read group" << endl
    << "  -n, --discard                 discard all output alignments (for profiling)" << endl
    << "computational parameters:" << endl
    << "  -F, --score-fraction FLOAT    select minimizers between hit caps until score is FLOAT of total [0.6]" << endl
    << "  -O, --no-chaining             disable seed chaining and all gapped alignment" << endl
    << "  -X, --xdrop                   use xdrop alignment for tails" << endl
    << "  -t, --threads INT             number of compute threads to use" << endl;
}

int main_giraffe(int argc, char** argv) {

    if (argc == 2) {
        help_giraffe(argv);
        return 1;
    }

    // initialize parameters with their default options
    string xg_name;
    string gbwt_name;
    string minimizer_name;
    string distance_name;
    // How close should two hits be to be in the same cluster?
    size_t distance_limit = 1000;
    bool progress = false;
    // Should we try chaining or just give up if we can't find a full length gapless alignment?
    bool do_chaining = true;
    // Whould we use the xdrop aligner for aligning tails?
    bool use_xdrop_for_tails = false;
    // What GAMs should we realign?
    vector<string> gam_filenames;
    // What FASTQs should we align.
    // Note: multiple FASTQs are not interpreted as paired.
    vector<string> fastq_filenames;
    // How many mappings per read can we emit?
    size_t max_multimaps = 1;
    // What sample name if any should we apply?
    string sample_name;
    // What read group if any should we apply?
    string read_group;
    // Should we throw out our alignments instead of outputting them?
    bool discard_alignments = false;
    
    vector<size_t> threads_to_run;
    
    int c;
    optind = 2; // force optind past command positional argument
    while (true) {
        static struct option long_options[] =
        {
            {"help", no_argument, 0, 'h'},
            {"xg-name", required_argument, 0, 'x'},
            {"gbwt-name", required_argument, 0, 'H'},
            {"minimizer-name", required_argument, 0, 'm'},
            {"dist-name", required_argument, 0, 'd'},
            {"progress", no_argument, 0, 'p'},
            {"gam-in", required_argument, 0, 'G'},
            {"fastq-in", required_argument, 0, 'f'},
            {"max-multimaps", required_argument, 0, 'M'},
            {"sample", required_argument, 0, 'N'},
            {"read-group", required_argument, 0, 'R'},
            {"discard", no_argument, 0, 'n'},
            {"no-chaining", no_argument, 0, 'O'},
            {"xdrop", no_argument, 0, 'X'},
            {"threads", required_argument, 0, 't'},
            {0, 0, 0, 0}
        };

        int option_index = 0;
        c = getopt_long (argc, argv, "hx:H:m:s:d:pG:f:M:N:R:nc:C:OXt:",
                         long_options, &option_index);


        // Detect the end of the options.
        if (c == -1)
            break;

        switch (c)
        {
            case 'x':
                xg_name = optarg;
                if (xg_name.empty()) {
                    cerr << "error:[vg gaffe] Must provide XG file with -x." << endl;
                    exit(1);
                }
                break;
                
            case 'H':
                gbwt_name = optarg;
                if (gbwt_name.empty()) {
                    cerr << "error:[vg gaffe] Must provide GBWT file with -H." << endl;
                    exit(1);
                }
                break;
                
            case 'm':
                minimizer_name = optarg;
                if (minimizer_name.empty()) {
                    cerr << "error:[vg gaffe] Must provide minimizer file with -m." << endl;
                    exit(1);
                }
                break;
                
            case 'd':
                distance_name = optarg;
                if (distance_name.empty()) {
                    cerr << "error:[vg gaffe] Must provide distance index file with -d." << endl;
                    exit(1);
                }
                break;

            case 'p':
                progress = true;
                break;
                
            case 'G':
                gam_filenames.push_back(optarg);
                break;
            
            case 'f':
                fastq_filenames.push_back(optarg);
                break;
                
            case 'M':
                max_multimaps = parse<size_t>(optarg);
                break;
            
            case 'N':
                sample_name = optarg;
                break;
                
            case 'R':
                read_group = optarg;
                break;
                
            case 'n':
                discard_alignments = true;
                break;

            case 'O':
                do_chaining = false;
                break;
                
            case 'X':
                use_xdrop_for_tails = true;
                break;
                
            case 't':
            {
                int num_threads = parse<int>(optarg);
                if (num_threads <= 0) {
                    cerr << "error:[vg gaffe] Thread count (-t) set to " << num_threads << ", must set to a positive integer." << endl;
                    exit(1);
                }
                threads_to_run.push_back(num_threads);
            }
                break;
                
            case 'h':
            case '?':
            default:
                help_giraffe(argv);
                exit(1);
                break;
        }
    }
    
    
    if (xg_name.empty()) {
        cerr << "error:[vg gaffe] Mapping requires an XG index (-x)" << endl;
        exit(1);
    }
    
    if (gbwt_name.empty()) {
        cerr << "error:[vg gaffe] Mapping requires a GBWT index (-H)" << endl;
        exit(1);
    }
    
    if (minimizer_name.empty()) {
        cerr << "error:[vg gaffe] Mapping requires a minimizer index (-m)" << endl;
        exit(1);
    }
    
    
    if (distance_name.empty()) {
        cerr << "error:[vg gaffe] Mapping requires a distance index (-d)" << endl;
        exit(1);
    }
    
    // create in-memory objects
    if (progress) {
        cerr << "Loading XG index " << xg_name << endl;
    }
    unique_ptr<XG> xg_index = vg::io::VPKG::load_one<XG>(xg_name);

    if (progress) {
        cerr << "Loading GBWT index " << gbwt_name << endl;
    }
    unique_ptr<gbwt::GBWT> gbwt_index = vg::io::VPKG::load_one<gbwt::GBWT>(gbwt_name);

    if (progress) {
        cerr << "Loading minimizer index " << minimizer_name << endl;
    }
    unique_ptr<MinimizerIndex> minimizer_index = vg::io::VPKG::load_one<MinimizerIndex>(minimizer_name);

    if (progress) {
        cerr << "Loading distance index " << distance_name << endl;
    }
    MinimumDistanceIndex distance_index;
    ifstream dist_in (distance_name);
    distance_index.load(dist_in);
    //unique_ptr<MinimumDistanceIndex> distance_index = vg::io::VPKG::load_one<MinimumDistanceIndex>(distance_name);
    

    //File to write the times to
    ofstream out ("mapping_times");
    out << "max_extensions\tmax_alignments\tcluster_score_threshold\tcluster_coverage_threshold\textension_set_threshold\thit_cap\thard_hit_cap\tminimizer_score_fraction\truntime(reads/thread/second)\tcorrect" << endl;

    // Set up the mapper
    if (progress) {
        cerr << "Initializing MinimizerMapper" << endl;
    }
    MinimizerMapper minimizer_mapper(xg_index.get(), gbwt_index.get(), minimizer_index.get(), &distance_index);
    minimizer_mapper.do_chaining = do_chaining;
 
    minimizer_mapper.max_multimaps = max_multimaps;
 
    minimizer_mapper.distance_limit = distance_limit;
 
    minimizer_mapper.use_xdrop_for_tails = use_xdrop_for_tails;

    size_t x = 0;
    // How many clusters should we extend?
    for (size_t max_extensions : {48})  {
        // How many extended clusters should we align, max?
        for (size_t max_alignments : { 8}) {
            //Throw away cluster with scores that are this amount below the best
            for (double cluster_score : {25, 75, 100, 0} ){
                //Throw away clusters with coverage this amount below the best 
                for (double cluster_coverage : {0.1, 0.3, 0.5} ) {
                    //Throw away extensions with scores that are this amount below the best
                    for (double extension_set : {5, 20, 70}) {
                        for (size_t hit_cap : { 10}) {
                            for (size_t hard_hit_cap : {300, 600}) {
                                for (double minimizer_score_fraction : {0.6}){
                                    if (progress) {
                                        cerr << "Aligning " << x++ << " of 80" << endl;
                                    }
                                    out << max_extensions << "\t" << max_alignments << "\t" << cluster_score << "\t" 
                                         << cluster_coverage << "\t" << extension_set << "\t" << hit_cap << "\t"
                                         << hard_hit_cap << "\t" << minimizer_score_fraction;

                                    minimizer_mapper.hit_cap = hit_cap;

                                    minimizer_mapper.hard_hit_cap = hard_hit_cap;

                                    minimizer_mapper.minimizer_score_fraction = minimizer_score_fraction;

                                    minimizer_mapper.max_extensions = max_extensions;

                                    minimizer_mapper.max_alignments = max_alignments;

                                    minimizer_mapper.cluster_score_threshold = cluster_score;

                                    minimizer_mapper.cluster_coverage_threshold = cluster_coverage;

                                    minimizer_mapper.extension_set_threshold = extension_set;


                                    minimizer_mapper.sample_name = sample_name;

                                    minimizer_mapper.read_group = read_group;
                                    
                                    if (threads_to_run.empty()) {
                                        // Use 0 to represent the default.
                                        threads_to_run.push_back(0);
                                    }

                                    //Map the real reads to get the times
                                    for (size_t thread_count_current : threads_to_run) {
                                    
                                        if (thread_count_current != 0) {
                                            omp_set_num_threads(thread_count_current);
                                        }
                                        
                                        // Work out the number of threads we will have
                                        size_t thread_count = omp_get_max_threads();

                                        // Set up counters per-thread for total reads mapped
                                        vector<size_t> reads_mapped_by_thread(thread_count, 0);
                                        
                                        // Have a place to log start time
                                        std::chrono::time_point<std::chrono::system_clock> start;
                                        
                                        {
                                            // Set up output to an emitter that will handle serialization
                                            // Discard alignments if asked to. Otherwise spit them out in GAM format.
                                            unique_ptr<AlignmentEmitter> alignment_emitter = discard_alignments ?
                                                make_unique<NullAlignmentEmitter>() :
                                                get_alignment_emitter("-", "GAM", {}, thread_count);

#ifdef USE_CALLGRIND
                                            // We want to profile the alignment, not the loading.
                                            CALLGRIND_START_INSTRUMENTATION;
#endif

                                            // Start timing overall mapping time now that indexes are loaded.
                                            start = std::chrono::system_clock::now();
                                            
                                            // Define how to align and output a read, in a thread.
                                            auto map_read = [&](Alignment& aln) {
                                                // Map the read with the MinimizerMapper.
                                                minimizer_mapper.map(aln, *alignment_emitter);
                                                // Record that we mapped a read.
                                                reads_mapped_by_thread.at(omp_get_thread_num())++;
                                            };
                                            
                                            for (auto& fastq_name : fastq_filenames) {
                                                // For every FASTQ file to map, map all its reads in parallel.
                                                fastq_unpaired_for_each_parallel(fastq_name, map_read);
                                            }
                                        
                                        } // Make sure alignment emitter is destroyed and all alignments are on disk.
                                        
                                        // Now mapping is done
                                        std::chrono::time_point<std::chrono::system_clock> end = std::chrono::system_clock::now();
                                        std::chrono::duration<double> elapsed_seconds = end-start;
                                        
                                        // How many reads did we map?
                                        size_t total_reads_mapped = 0;
                                        for (auto& reads_mapped : reads_mapped_by_thread) {
                                            total_reads_mapped += reads_mapped;
                                        }
                                        
                                        out  << "\t" << ((total_reads_mapped / elapsed_seconds.count()) / thread_count);
                                    }
                                    //Map the simulated reads to get accuracy
                                    size_t correct_count = 0;
                                    threads_to_run[0] = 0;
                                    for (size_t thread_count_current : threads_to_run) {
                                    
                                        if (thread_count_current != 0) {
                                            omp_set_num_threads(thread_count_current);
                                        }
                                        
                                        // Work out the number of threads we will have
                                        size_t thread_count = omp_get_max_threads();
                                        
                                        {
                                            // Set up output to an emitter that will handle serialization
                                            // Discard alignments if asked to. Otherwise spit them out in GAM format.
                                            unique_ptr<AlignmentEmitter> alignment_emitter = discard_alignments ?
                                                make_unique<NullAlignmentEmitter>() :
                                                get_alignment_emitter("-", "GAM", {}, thread_count);

#ifdef USE_CALLGRIND
                                            // We want to profile the alignment, not the loading.
                                            CALLGRIND_START_INSTRUMENTATION;
#endif

                                            // Define how to align and output a read, in a thread.
                                            auto map_read = [&](Alignment& aln) {
                                                // Map the read with the MinimizerMapper.
                                                bool correct = minimizer_mapper.map(aln, *alignment_emitter);
                                                if (correct) {
                                                    correct_count ++;
                                                }
                                            };
                                                
                                            for (auto& gam_name : gam_filenames) {
                                                // For every GAM file to remap
                                                get_input_file(gam_name, [&](istream& in) {
                                                    // Open it and map all the reads in parallel.
                                                    vg::io::for_each_parallel<Alignment>(in, map_read);
                                                });
                                            }
                                        
                                        } // Make sure alignment emitter is destroyed and all alignments are on disk.
                                    }
                                    out << "\t" << correct_count << endl;
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    out.close();
        
    return 0;
}

// Register subcommand
static Subcommand vg_giraffe("giraffe-test", "Graph (i r ) Alignment Format Fast Emitter", DEVELOPMENT, main_giraffe);


