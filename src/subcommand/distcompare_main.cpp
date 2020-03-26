/** \file crash_main.cpp
 *
 * Defines the "vg crash" subcommand, which throws errors to test the backtrace system.
 */
#include <vg/io/vpkg.hpp>
#include<chrono>

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
#include "snarls.hpp"
#include "vg.hpp"
#include "xg.hpp"
#include <bdsg/overlay_helper.hpp>
#include "../algorithms/dijkstra.hpp"


using namespace std;
using namespace vg;
using namespace vg::subcommand;

/*
int64_t dijkstra(xg::XG& graph, pos_t pos1, pos_t pos2){
    //Distance using djikstras algorithm

    auto cmp = [] (pair<pair<vg::id_t, bool>, int64_t> x,
                   pair<pair<vg::id_t, bool>, int64_t> y ) {
        return (x.second > y.second);
    };

    int64_t shortestDistance = -1;
    if (get_id(pos1) == get_id(pos2) && is_rev(pos1) == is_rev(pos2)) { //if positions are on the same node

        handle_t handle = graph.get_handle(get_id(pos1), 0);
        int64_t nodeSize = graph.get_length(handle);
        int64_t offset1 = get_offset(pos1);
        int64_t offset2 = get_offset(pos2);

        if (offset1 <= offset2) {
            shortestDistance = offset2-offset1+1; //+1 to be consistent
        }

    }


    priority_queue< pair<pair<vg::id_t, bool> , int64_t>,
                    vector<pair<pair<vg::id_t, bool>, int64_t>>,
                          decltype(cmp)> reachable(cmp);
    handle_t currHandle = graph.get_handle(get_id(pos1), is_rev(pos1));

    int64_t dist = graph.get_length(currHandle) - get_offset(pos1);

    auto addFirst = [&](const handle_t& h) -> bool {
        pair<vg::id_t, bool> node = make_pair(graph.get_id(h),
                                          graph.get_is_reverse(h));
        reachable.push(make_pair(node, dist));
        return true;
    };

    graph.follow_edges(currHandle, false, addFirst);
    unordered_set<pair<vg::id_t, bool>> seen;
    seen.insert(make_pair(get_id(pos1), is_rev(pos1)));
    while (reachable.size() > 0) {
        pair<pair<vg::id_t, bool>, int64_t> next = reachable.top();
        reachable.pop();
        pair<vg::id_t, bool> currID = next.first;
        dist = next.second;
        if (seen.count(currID) == 0) {

            seen.insert(currID);
            currHandle = graph.get_handle(currID.first, currID.second);
            int64_t currDist = graph.get_length(currHandle);

            auto addNext = [&](const handle_t& h) -> bool {
                pair<vg::id_t, bool> node = make_pair(graph.get_id(h),
                                                graph.get_is_reverse(h));
                reachable.push(make_pair(node, currDist + dist));
                return true;
            };
            graph.follow_edges(currHandle, false, addNext);

        }

        if (currID.first == get_id(pos2) && currID.second == is_rev(pos2)){
        //Dist is distance to beginning or end of node containing pos2

            if (is_rev(pos2) == currID.second) {
                dist = dist + get_offset(pos2) + 1;
            } else {
                handle_t handle = graph.get_handle(get_id(pos2), 0);
                dist = dist +  graph.get_length(handle); -
                            get_offset(pos2);
            }
            if (shortestDistance == -1) {shortestDistance = dist;}
            else {shortestDistance = min(dist, shortestDistance);}
        }

    }
    return shortestDistance == -1 ? -1 : shortestDistance-1;
};
*/

void help_distcompare(char** argv){
    cerr << "usage: " << argv[0] << " distcompare [options]" << endl
         << "Output info about distance calculations" << endl
         << endl
         << "options: " << endl
         << "    -x, --xg-name, FILE       use this xg index (required)" << endl
         << "    -s, --snarl-name, FILE       use this snarl index (required)" << endl
         << "    -d, --dist-name, FILE       use this distance index (required)" << endl
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
        c = getopt_long (argc, argv, "hx:d:",
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
            if (dist_name.empty()) {
                cerr << "error: [vg distcompare] Must provide distance file with -d." << endl;
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
        cerr << "error: [vg distcompare] Must provide dist file with -d." << endl;
        exit(1);
    }

    unique_ptr<xg::XG> xg_index = vg::io::VPKG::load_one<xg::XG>(xg_name);


    unique_ptr<MinimumDistanceIndex> distance_index = vg::io::VPKG::load_one<MinimumDistanceIndex>(dist_name);

 
    PathPositionHandleGraph* handle_graph = nullptr;
    unique_ptr<PathHandleGraph> path_handle_graph;
    bdsg::PathPositionOverlayHelper overlay_helper;
    if (!xg_name.empty()) {
        path_handle_graph = vg::io::VPKG::load_one<PathHandleGraph>(xg_name);
        handle_graph = overlay_helper.apply(path_handle_graph.get());
    }
    vg::PathOrientedDistanceMeasurer measurer(handle_graph);
    cerr << "Loaded xg again as a handle graph" << endl;
    random_device seed_source;
    default_random_engine generator(seed_source());



    for (size_t i = 0 ; i < 100000 ; i ++) {
        //Min distances
        //Find distances between random positions in the graph
        //
        //TODO: For 1kg whole_genome_graph specifically
        std::uniform_int_distribution<int64_t> randNodeIndex(1,306009791 );

        vg::id_t nodeID1 = randNodeIndex(generator);
        vg::id_t nodeID2 = randNodeIndex(generator);
        size_t len1 = xg_index->get_length(xg_index->get_handle(nodeID1, false));
        size_t len2 = xg_index->get_length(xg_index->get_handle(nodeID2, false));

        vg::off_t offset1 = std::uniform_int_distribution<int>(0,len1 - 1)(generator);
        vg::off_t offset2 = std::uniform_int_distribution<int>(0,len2 - 1)(generator);
        bool rev1 = false;
        bool rev2 = false;
        if (distance_index->minDistance(make_pos_t(nodeID1, false, offset1 ), make_pos_t(nodeID2, false, offset2 ))){
            rev1 = false;
            rev2 = false;
        } else if (distance_index->minDistance(make_pos_t(nodeID1, true, offset1 ), make_pos_t(nodeID2, true, offset2 ))){
            rev1 = true;
            rev2 = true;
        } else if (distance_index->minDistance(make_pos_t(nodeID1, false, offset1 ), make_pos_t(nodeID2, true, offset2 ))){
            rev1 = false;
            rev2 = true;
        } else {
            rev1 = true;
            rev2 = false;
        }

        pos_t pos1 = make_pos_t(nodeID1, rev1, offset1 );
        pos_t pos2 = make_pos_t(nodeID2, rev2, offset2 );


    
/*
        //Generate random positions by choosing snarls then nodes, then position
        size_t maxPos = xg_index->seq_length;
        size_t offset1 = std::uniform_int_distribution<int>(1, maxPos)(generator);
        size_t offset2 = std::uniform_int_distribution<int>(1, maxPos)(generator);
        vg::id_t nodeID1 = xg_index->node_at_seq_pos(offset1);
        vg::id_t nodeID2 = xg_index->node_at_seq_pos(offset2);
        bool rev1 = std::uniform_int_distribution<int>(0, 1)(generator);
        bool rev2 = std::uniform_int_distribution<int>(0, 1)(generator);


 
        pos_t pos1 = make_pos_t(nodeID1, rev1, 0);
        pos_t pos2 = make_pos_t(nodeID2, rev2, 0);
*/
       
        //Find min distance 
        std::chrono::steady_clock::time_point start1 = std::chrono::steady_clock::now();
        int64_t minDist = distance_index->minDistance(pos1, pos2);
        std::chrono::steady_clock::time_point end1 = std::chrono::steady_clock::now();
        auto t1 =std::chrono::duration_cast<std::chrono::nanoseconds>( end1 - start1).count()/ 1000000000.0;


        //Find old distance
        std::chrono::steady_clock::time_point start3 = std::chrono::steady_clock::now();
        int64_t oldDist = abs(measurer.oriented_distance(pos1, pos2));
        std::chrono::steady_clock::time_point end3 = std::chrono::steady_clock::now();
        auto t3 = std::chrono::duration_cast<std::chrono::nanoseconds>(end3 - start3).count()/ 1000000000.0;


        //Find dijkstra distance 
        std::chrono::steady_clock::time_point start4 = std::chrono::steady_clock::now();
        int64_t dijkstraDist = -1;
        algorithms::dijkstra(handle_graph, xg_index->get_handle(get_id(pos1), is_rev(pos1)),
                 [&] (const handle_t& handle, size_t dist) {
                     if (handle == xg_index->get_handle(get_id(pos2), is_rev(pos2))){
                         dijkstraDist = dist;
                         return false;
                     } else {
                         return true;
                     }
                 });
        std::chrono::steady_clock::time_point end4 = std::chrono::steady_clock::now();

        auto t4 = std::chrono::duration_cast<std::chrono::nanoseconds>(end4 - start4).count() / 1000000000.0;

if (true) {

//min dist time, dijkstra time, path-based time
cout <<  t1 << "\t" <<  t4 << "\t" <<  t3 << "\t" <<  minDist << "\t" << dijkstraDist << "\t" << oldDist << endl;
}

    }                 
    
    return 0;
}
// Register subcommand
static Subcommand vg_crash("distcompare", "compare distance calculations", DEVELOPMENT, main_distcompare);
