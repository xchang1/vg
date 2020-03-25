/** \file crash_main.cpp
 *
 * Defines the "vg crash" subcommand, which throws errors to test the backtrace system.
 */
#include <vg/io/vpkg.hpp>

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
static pair<unordered_set<Node*>, unordered_set<Edge*> > pb_contents(
    VG& graph, const pair<unordered_set<vg::id_t>, unordered_set<edge_t> >& contents) {
    pair<unordered_set<Node*>, unordered_set<Edge*> > ret;
    for (vg::id_t node_id : contents.first) {
        ret.first.insert(graph.get_node(node_id));
    }
    for (const edge_t& edge_handle : contents.second) {
        Edge* edge = graph.get_edge(NodeTraversal(graph.get_node(graph.get_id(edge_handle.first)),
                                                  graph.get_is_reverse(edge_handle.first)),
                                    NodeTraversal(graph.get_node(graph.get_id(edge_handle.second)),
                                                  graph.get_is_reverse(edge_handle.second)));
        ret.second.insert(edge);
    }
    return ret;
}

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
         
    string vg_name;
    string xg_name;
    string snarl_name;
    string dist_name;
    int c;
    optind = 2; // force optind past command positional argument
    while (true) {
        static struct option long_options[] =
            {
                /* These options set a flag. */
                {"help", no_argument, 0, 'h'},
                {"xg-name", required_argument, 0, 'x'},
                {"vg-name", required_argument, 0, 'v'},
                {"dist-name", required_argument, 0, 'd'},
                {"snarl-name", required_argument, 0, 's'},
                {0, 0, 0, 0}
            };
        int option_index = 0;
        c = getopt_long (argc, argv, "hx:d:s:v:",
                         long_options, &option_index);
        /* Detect the end of the options. */
        if (c == -1)
            break;
        switch (c)
        {
        case 'v':
            vg_name = optarg;
            if (vg_name.empty()) {
                cerr << "error: [vg distcompare] Must provide VG file with -v." << endl;
                exit(1);
            }
            break;
        case 'x':
            xg_name = optarg;
            if (xg_name.empty()) {
                cerr << "error: [vg distcompare] Must provide XG file with -x." << endl;
                exit(1);
            }
            break;
        case 's':
            snarl_name = optarg;
            if (snarl_name.empty()) {
                cerr << "error: [vg distcompare] Must provide snarl file with -s." << endl;
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
    if (vg_name.empty()) {
        cerr << "error: [vg distcompare] Must provide VG file with -v." << endl;
        exit(1);
    }
    if (xg_name.empty()) {
        cerr << "error: [vg distcompare] Must provide XG file with -x." << endl;
        exit(1);
    }
    if (snarl_name.empty()) {
        cerr << "error: [vg distcompare] Must provide snarl file with -s." << endl;
        exit(1);
    }
    if (dist_name.empty()) {
        cerr << "error: [vg distcompare] Must provide dist file with -d." << endl;
        exit(1);
    }

    unique_ptr<xg::XG> xg_index = vg::io::VPKG::load_one<xg::XG>(xg_name);
    unique_ptr<vg::VG> vg_index = vg::io::VPKG::load_one<vg::VG>(vg_name);


    unique_ptr<MinimumDistanceIndex> distance_index = vg::io::VPKG::load_one<MinimumDistanceIndex>(dist_name);

    ifstream snarl_stream(snarl_name);
    if (!snarl_stream) {
        cerr << "error:[vg distcompare] Cannot open snarl file " << snarl_name << endl;
        exit(1);
    }
    SnarlManager* snarl_manager = new SnarlManager(snarl_stream) ;
    snarl_stream.close();
 
    random_device seed_source;
    default_random_engine generator(seed_source());
    vector<clock_t> minTimes;
    vector<clock_t> oldTimes;
    vector<clock_t> dijkstraTimes;

          vector<const Snarl*> allSnarls;
          auto addSnarl = [&] (const Snarl* s) {
              allSnarls.push_back(s);
          };
          snarl_manager->for_each_snarl_preorder(addSnarl);
          std::uniform_int_distribution<int> randSnarlIndex(0, allSnarls.size()-1);

    for (size_t i = 0 ; i < 100000 ; i ++) {
        //Min distances
        //Find distances between random positions in the graph
        //Outputs: my distance /t old distance /t time for my calculation /t 
        //           time for old calculation /t
        
                const Snarl* snarl1 = allSnarls[randSnarlIndex(generator)];
                const Snarl* snarl2 = allSnarls[randSnarlIndex(generator)];

                pair<unordered_set<Node*>, unordered_set<Edge*>> contents1 =
                    pb_contents(*vg_index, snarl_manager->shallow_contents(snarl1, *vg_index, true));
                pair<unordered_set<Node*>, unordered_set<Edge*>> contents2 =
                    pb_contents(*vg_index, snarl_manager->shallow_contents(snarl2, *vg_index, true));

                vector<Node*> nodes1 (contents1.first.begin(), contents1.first.end());
                vector<Node*> nodes2 (contents2.first.begin(), contents2.first.end());


                std::uniform_int_distribution<int> randNodeIndex2(0,nodes2.size()-1);
                std::uniform_int_distribution<int> randNodeIndex1(0,nodes1.size()-1);

                Node* node1 = nodes1[randNodeIndex1(generator)];
                Node* node2 = nodes2[randNodeIndex2(generator)];
                vg::id_t nodeID1 = node1->id();
                vg::id_t nodeID2 = node2->id();

                vg::off_t offset1 = std::uniform_int_distribution<int>(0,node1->sequence().size() - 1)(generator);
                vg::off_t offset2 = std::uniform_int_distribution<int>(0,node2->sequence().size() - 1)(generator);

                pos_t pos1 = make_pos_t(nodeID1,
                  std::uniform_int_distribution<int>(0,1)(generator) == 0,offset1 );
                pos_t pos2 = make_pos_t(nodeID2,
                  std::uniform_int_distribution<int>(0,1)(generator) == 0, offset2 );

    
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
        clock_t start1 = clock();
        int64_t minDist = distance_index->minDistance(pos1, pos2);
        clock_t end1 = clock();
        clock_t t1 = end1 - start1;
        minTimes.push_back(t1);


        //Find old distance
        clock_t start3 = clock();
    PathPositionHandleGraph* handle_graph = nullptr;
    unique_ptr<PathHandleGraph> path_handle_graph;
    bdsg::PathPositionOverlayHelper overlay_helper;
    if (!xg_name.empty()) {
        path_handle_graph = vg::io::VPKG::load_one<PathHandleGraph>(xg_name);
        handle_graph = overlay_helper.apply(path_handle_graph.get());
    }

        vg::PathOrientedDistanceMeasurer measurer(handle_graph);
        int64_t oldDist = abs(measurer.oriented_distance(pos1, pos2));
        clock_t end3 = clock();
        clock_t t3 = end3 - start3;
        oldTimes.push_back(t3);


        //Find dijkstra distance 
        clock_t start4 = clock();
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
        clock_t end4 = clock();

        clock_t t4 = end4 - start4;
        dijkstraTimes.push_back(t4);

if (true) {

//min dist time, dijkstra time, path-based time
cout << t1 << "\t" << t4 << "\t" << t3 <<  endl;
cerr << minDist << "\t" << dijkstraDist << "\t" << oldDist << endl;
}
    }                 
    
    return 0;
}
// Register subcommand
static Subcommand vg_crash("distcompare", "compare distance calculations", DEVELOPMENT, main_distcompare);
