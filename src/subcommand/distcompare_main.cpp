/** \file crash_main.cpp
 *
 * Defines the "vg crash" subcommand, which throws errors to test the backtrace system.
 */
#include <omp.h>
#include <unistd.h>
#include <getopt.h>
#include <stdlib.h>
#include <iostream>
#include <cassert>
#include <csignal>
#include <random>
#include <time.h>
#include "subcommand.hpp"
#include "min_distance.hpp"
#include "cluster.hpp"
#include "../benchmark.hpp"
#include "../utility.hpp"
#include <vg/io/stream.hpp>
#include <vg/io/vpkg.hpp>

using namespace std;
using namespace vg;
using namespace vg::subcommand;
void help_distcompare(char** argv){
    cerr << "usage: " << argv[0] << " distcompare [options]" << endl
         << "Output info about distance calculations" << endl
         << endl
         << "options: " << endl
         << "    -x, --xg-name, FILE       use this xg index (required)" << endl
         << "    -v, --vg-name, FILE       use this vg index (required)" << endl
         << "    -d, --dist-name, FILE     use this distance index (required)" << endl
         << endl;
}
int main_distcompare(int argc, char** argv){
         
    string xg_name;
    string dist_name;
    int c;
    optind = 2; // force optind past command positional argument
    while (true) {
        static struct option long_options[] =
            {
                /* These options set a flag. */
                {"help", no_argument, 0, 'h'},
                {"xg-name", required_argument, 0, 'x'},
                {"dist-name", required_argument, 0, 'd'},
                {0, 0, 0, 0}
            };
        int option_index = 0;
        c = getopt_long (argc, argv, "hx:v:d:",
                         long_options, &option_index);
        /* Detect the end of the options. */
        if (c == -1)
            break;
        switch (c)
        {
        case 'x':
            xg_name = optarg;
            if (xg_name.empty()) {
                cerr << "error: [vg distcompare] Must provide XG file with -x." << endl;
                exit(1);
            }
            break;
        case 'd':
            dist_name = optarg;
            if (distance_name.empty()) {
                cerr << "error: [vg distcompare] Must provide distance index file with -d." << endl;
                exit(1);
            }
            break;
        case 'h':
        case '?':
        default:
            help_distcompare(argv);
            exit(1);
            break;
        }
    }
    //Check that all required arguments are given
    if (xg_name.empty()) {
        cerr << "error: [vg distcompare] Must provide XG file with -x." << endl;
        exit(1);
    }
    if (dist_name.empty()) {
        cerr << "error: [vg distcompare] Must provide distance index file with -d." << endl;
        exit(1);
    }
    
    //Get file streams for xg and vg
    

    xg::XG xg_index = vg::io::VPKG::load_one<PathPositionHandleGraph>(xg_name);
    MinimumDistanceIndex di = vg::io::VPKG::load_one<DistanceIndex>(dist_name);

    TargetValueSearch tvs(xg_index, new TipAnchoredMaxDistance(di), 
                              new SnarlMinDistance(di));
    t = clock() - t;
    
    cerr << "ticks per second " << CLOCKS_PER_SEC << endl;
    cout << "Time to create distance index: " << t << endl;
    cout << endl;
 
    random_device seed_source;
    default_random_engine generator(seed_source());
    vector<int64_t> minDists;
    vector<int64_t> maxDists;
    vector<int64_t> oldDists;
    vector<int64_t> tvsDists;


    vector<clock_t> minTimes;
    vector<clock_t> maxTimes;
    vector<clock_t> oldTimes;
    vector<clock_t> tvsTimes;


/*
        pos_t p1 = make_pos_t(2610, false, 0);
        pos_t p2 = make_pos_t(2634, false, 0);
cerr << di.minDistance(p1, p2) << " " <<  abs(xg_index.closest_shared_path_oriented_distance(
                    2610, 0, false, 2634, 0, false)) << endl ;
assert(di.minDistance(p1, p2) <= abs(xg_index.closest_shared_path_oriented_distance(2610, 0, false, 2634, 0, false)));
*/
    for (size_t i = 0 ; i < 100000 ; i ++) {
        //Min distances
        //Find distances between random positions in the graph
        //Outputs: my distance /t old distance /t time for my calculation /t 
        //           time for old calculation /t
        
    
        //Generate random positions by choosing snarls then nodes, then position
        size_t maxPos = xg_index.seq_length;
        size_t offset1 = uniform_int_distribution<int>(1, maxPos)(generator);
        size_t offset2 = uniform_int_distribution<int>(1, maxPos)(generator);
        vg::id_t nodeID1 = xg_index.node_at_seq_pos(offset1);
        vg::id_t nodeID2 = xg_index.node_at_seq_pos(offset2);

 
        pos_t pos1 = make_pos_t(nodeID1, false, 0);
        pos_t pos2 = make_pos_t(nodeID2, false, 0);
       
        //Find min distance 
        clock_t start1 = clock();
        int64_t minDist = di.minDistance(pos1, pos2);
        clock_t end1 = clock();
        clock_t t1 = end1 - start1;
        minDists.push_back(minDist);
        minTimes.push_back(t1);

        //Find max distance 
        clock_t start2 = clock();
        int64_t maxDist = di.maxDistance(pos1, pos2);
        clock_t end2 = clock();

        clock_t t2 = end2 - start2;
        maxDists.push_back(maxDist);
        maxTimes.push_back(t2);

        //Find old distance
        clock_t start3 = clock();
        int64_t oldDist = abs(xg_index.closest_shared_path_oriented_distance(nodeID1, 0, false, nodeID2, 0, false));
        clock_t end3 = clock();
        clock_t t3 = end3 - start3;
        oldDists.push_back(oldDist);
        oldTimes.push_back(t3);


        //Find max distance 
        clock_t start4 = clock();
        int64_t tvsDist = tvs.tv_path_length(pos1, pos2, oldDist, 20);
        clock_t end4 = clock();

        clock_t t4 = end4 - start4;
        tvsDists.push_back(tvsDist);
        tvsTimes.push_back(t4);

if (oldDist != INT64_MAX) {

cout << minDist << "\t" << maxDist << "\t" << oldDist << "\t" << tvsDist << "\t" << t1 << "\t" << t2 << "\t" << t3 << "\t" << t4 << "\t" << nodeID1 << "\t" << nodeID2 << endl;
}
    }                 
    
    return 0;
}
// Register subcommand
static Subcommand vg_crash("distcompare", "compare distance calculations", DEVELOPMENT, main_distcompare);
