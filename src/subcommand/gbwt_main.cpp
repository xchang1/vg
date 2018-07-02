/** \file gbwt_main.cpp
 *
 * Defines the "vg gbwt" subcommand, which wraps up access for commands we'd otherwise find
 * in the gbwt submodule.  */

#include <omp.h>
#include <unistd.h>
#include <getopt.h>

#include <list>
#include <fstream>

#include "subcommand.hpp"

#include "../xg.hpp"

using namespace std;
using namespace vg;
using namespace vg::subcommand;

#include <unistd.h>

#include <gbwt/dynamic_gbwt.h>


void help_gbwt(char** argv) {
    cerr << "usage: " << argv[0] << " [options] [args]" << endl
         << "Manipulate GBWTs." << endl
         << "merging:" << endl
         << "    -m, --merge            merge the GBWT files from the input args and write to output" << endl
         << "    -o, --output X         write output GBWT to X" << endl
         << "    -b, --batches N        use batches of N sequences for merging (default: " << gbwt::DynamicGBWT::MERGE_BATCH_SIZE << ")" << endl
         << "    -f, --fast             fast merging algorithm (node ids must not overlap; implies -m)" << endl
         << "    -p, --progress         show progress and statistics" << endl
         << "threads:" << endl
         << "    -c, --count-threads    print the number of threads" << endl
         << "    -e, --extract FILE     extract threads in SDSL format to FILE" << endl
         << endl;
}

int main_gbwt(int argc, char** argv)
{
    if (argc == 2) {
        help_gbwt(argv);
        return 1;
    }

    bool merge = false;
    size_t batch_size = gbwt::DynamicGBWT::MERGE_BATCH_SIZE;
    bool fast_merging = false;
    bool show_progress = false;
    bool count_threads = false;
    string gbwt_output, thread_output;

    int c;
    optind = 2; // force optind past command positional argument    
    while (true) {
        static struct option long_options[] =
            {
                // Merging
                {"merge", no_argument, 0, 'm'},
                {"output", required_argument, 0, 'o'},
                {"batches", required_argument, 0, 'b'},
                {"fast", no_argument, 0, 'f'},
                {"progress",  no_argument, 0, 'p'},

                // Threads
                {"count-threads", no_argument, 0, 'c'},
                {"extract", required_argument, 0, 'e'},

                {"help", no_argument, 0, 'h'},
                {0, 0, 0, 0}
            };

        int option_index = 0;
        c = getopt_long(argc, argv, "mo:b:fpce:h?", long_options, &option_index);

        /* Detect the end of the options. */
        if (c == -1)
            break;

        switch (c)
        {
        // Merging
        case 'm':
            merge = true;
            break;
        case 'o':
            gbwt_output = optarg;
            break;
        case 'b':
            batch_size = stoul(optarg);
            break;
        case 'f':
            fast_merging = true;
            merge = true;
            break;
        case 'p':
            show_progress = true;
            break;

        case 'c':
            count_threads = true;
            break;
        case 'e':
            thread_output = optarg;
            break;

        case 'h':
        case '?':
            /* getopt_long already printed an error message. */
            help_gbwt(argv);
            exit(1);
            break;

        default:
            abort ();
        }
    }

    gbwt::Verbosity::set(gbwt::Verbosity::SILENT);

    if (merge) {

        // Ugly hack here. GBWT prints to stdout, and we want to direct it to stderr.
        std::streambuf* cout_buf = cout.rdbuf();
        cout.rdbuf(cerr.rdbuf());

        size_t input_files = argc - optind;
        size_t total_inserted = 0;
        if (input_files <= 1) {
            cerr << "[vg gbwt] error: at least two input gbwt files required to merge" << endl;
            return 1;
        }
        if (gbwt_output.empty()) {
            cerr << "[vg gbwt] error: output file must be specified with -o" << endl;
        }
        if (show_progress) {
            gbwt::printHeader("Algorithm"); cout << (fast_merging ? "fast" : "insert") << endl;
            gbwt::printHeader("Input files"); cout << input_files << endl;
            gbwt::printHeader("Output name"); cout << gbwt_output << endl;
            if(!fast_merging) { gbwt::printHeader("Batch size"); cout << batch_size << endl; }
            cout << endl;
        }

        double start = gbwt::readTimer();

        if(fast_merging)
        {
            vector<gbwt::GBWT> indexes(argc - optind);
            for(int i = optind; i < argc; i++)
            {
                string input_name = argv[i];
                sdsl::load_from_file(indexes[i - optind], input_name);
                if (show_progress) {
                    gbwt::printStatistics(indexes[i - optind], input_name);
                }
                total_inserted += indexes[i - optind].size();
            }
            gbwt::GBWT merged(indexes);
            sdsl::store_to_file(merged, gbwt_output);
            if (show_progress) {
                gbwt::printStatistics(merged, gbwt_output);
            }
        }
        else
        {
            gbwt::DynamicGBWT index;
            {
                string input_name = argv[optind];
                sdsl::load_from_file(index, input_name);
                if (show_progress) {
                    gbwt::printStatistics(index, input_name);
                }
            }
            for (int curr = optind + 1; curr < argc; curr++)
            {
                string input_name = argv[curr];
                gbwt::GBWT next;
                sdsl::load_from_file(next, input_name);
                if (show_progress) {
                    gbwt::printStatistics(next, input_name);
                }
                index.merge(next, batch_size);
                total_inserted += next.size();
            }
            sdsl::store_to_file(index, gbwt_output);
            if (show_progress) { 
                gbwt::printStatistics(index, gbwt_output);
            }
        }
        
        double seconds = gbwt::readTimer() - start;

        if (show_progress) {

            cout << "Inserted " << total_inserted << " nodes in " << seconds << " seconds ("
                      << (total_inserted / seconds) << " nodes/second)" << endl;
            cout << "Memory usage " << gbwt::inGigabytes(gbwt::memoryUsage()) << " GB" << endl;
            cout << endl;
        }

        // Revert the hack.
        cout.rdbuf(cout_buf);
    }

    if (count_threads) {
        if (optind >= argc) {
            cerr << "[vg gbwt] error: no input files given" << endl;
            return 1;
        }

        size_t total_threads = 0;
        for (int i = optind; i < argc; i++) {
            gbwt::GBWT index;
            sdsl::load_from_file(index, argv[i]);
            total_threads += index.sequences() / 2; // Ignore reverse complements.
        }
        cout << total_threads << endl;
    }

    if (!thread_output.empty()) {
        if (optind + 1 != argc) {
            cerr << "[vg gbwt] error: option -e requires one input file" << endl;
            return 1;
        }

        gbwt::GBWT index;
        sdsl::load_from_file(index, argv[optind]);
        gbwt::size_type node_width = gbwt::bit_length(index.sigma() - 1);
        gbwt::text_buffer_type out(thread_output, std::ios::out, gbwt::MEGABYTE, node_width);
        for (gbwt::size_type id = 0; id < index.sequences(); id += 2) { // Ignore reverse complements.
            gbwt::vector_type sequence = index.extract(id);
            for (auto node : sequence) {
                out.push_back(node);
            }
            out.push_back(gbwt::ENDMARKER);
        }
        out.close();
    }

    return 0;
}


// Register subcommand
static Subcommand vg_gbwt("gbwt", "Manipuate GBWTs", main_gbwt);

