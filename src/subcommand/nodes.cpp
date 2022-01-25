#include "subcommand.hpp"
#include <vg/io/stream.hpp>
#include <vg/io/vpkg.hpp>
#include"../snarls.hpp"

#include <unistd.h>
#include <getopt.h>

using namespace vg;
using namespace vg::subcommand;
using namespace vg::io;

void help_nodes(char** argv) {
    cerr << "usage: " << argv[0] << " nodes [options] >nodes.txt" << endl
         << "options:" << endl
         << "    -g, --graph-name  [graph]      The graph file" << endl
         << "    -s, --snarls-name [snarls]    The snarl file" << endl;

}

int main_nodes(int argc, char** argv) {

    if (argc == 2) {
        help_nodes(argv);
        return 1;
    }

    string graph_name;
    string snarls_name;
    int c;
    optind = 2; // force optind past command positional argument
    while (true) {
        static struct option long_options[] =
            {
                //{"verbose", no_argument,       &verbose_flag, 1},
                {"graph-name", required_argument, 0, 'g'},
                {"snarls-name", required_argument, 0, 's'},
                {0, 0, 0, 0}
            };

        int option_index = 0;
        c = getopt_long (argc, argv, "s:g:",
                         long_options, &option_index);

        // Detect the end of the options.
        if (c == -1)
            break;

        switch (c)
        {

        case 'g':
            graph_name = optarg;
            break;
        case 's':
            snarls_name = optarg;
            break;

        default:
            abort ();
        }
    }
    if (optind < argc) {
        cerr << "[vg nodes] nodes does not accept positional arguments" << endl;
        return 1;
    }

    if (graph_name.empty()) {
        cerr << "[vg nodes] nodes requires -g" << endl;
        return 1;
    }
    if (snarls_name.empty()) {
        cerr << "[vg nodes] nodes requires -s" << endl;
        return 1;
    }

    //Load snarls
    ifstream snarl_stream(snarls_name);
    if (!snarl_stream) {
        cerr << "error: [vg index] cannot open Snarls file" << endl;
        exit(1);
    }
    SnarlManager* snarl_manager = new SnarlManager(snarl_stream);
    snarl_stream.close();

    //Load graph
    auto graph = vg::io::VPKG::load_one<HandleGraph>(graph_name);

    size_t max_snarl_size = 0;
    vg::id_t start_id = 0;
    bool start_orientation = false;
    //Go through the snarls and find the biggest one
    snarl_manager->for_each_snarl_preorder([&] (const Snarl* snarl)->void {
        size_t size = snarl_manager->shallow_contents(snarl, *graph, false).first.size();
        if (size > max_snarl_size) {
            start_id = snarl->start().node_id();
            start_orientation = snarl->start().backward();
            max_snarl_size = size; 
        }
    });

    //Print the id of every node in the largest snarl
    const Snarl* biggest_snarl = snarl_manager->into_which_snarl(start_id, start_orientation);
    auto snarl_contents = snarl_manager->deep_contents(biggest_snarl, *graph, true);
    for (const vg::id_t id: snarl_contents.first) {
        cout << id << endl;
    }
   
    return 0;

}

static Subcommand vg_msga("nodes", "print the nodes in the largest snarl to stdout", DEVELOPMENT, main_nodes);
