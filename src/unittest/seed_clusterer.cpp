#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <set>
#include "../json2pb.h"
#include "../vg.hpp"
#include "catch.hpp"
#include "../snarls.hpp"
#include "../cactus_snarl_finder.hpp"
#include "../genotypekit.hpp"
#include "random_graph.hpp"
#include "../seed_clusterer.hpp"
#include <random>
#include <time.h>
#include <structures/union_find.hpp>

//#define print

namespace vg {
namespace unittest {
    TEST_CASE( "looping chain of nested unary snarls",
                   "[cluster]" ) {
        VG graph;

        Node* n1 = graph.create_node("GCA");
        Node* n2 = graph.create_node("T");
        Node* n3 = graph.create_node("G");
        Node* n4 = graph.create_node("CTGA");
        Node* n5 = graph.create_node("GCA");
        Node* n6 = graph.create_node("T");

        Edge* e1 = graph.create_edge(n1, n2);
        Edge* e2 = graph.create_edge(n1, n3);
        Edge* e3 = graph.create_edge(n2, n4);
        Edge* e4 = graph.create_edge(n3, n4);
        Edge* e5 = graph.create_edge(n4, n5);
        Edge* e6 = graph.create_edge(n4, n6);
        Edge* e7 = graph.create_edge(n5, n6);
        Edge* e8 = graph.create_edge(n6, n6, false, true);
        Edge* e9 = graph.create_edge(n1, n1, true, false);

        CactusSnarlFinder bubble_finder(graph);
        SnarlManager snarl_manager = bubble_finder.find_snarls();

        MinimumDistanceIndex dist_index (&graph, &snarl_manager);
        SnarlSeedClusterer clusterer(dist_index);
        
        //graph.to_dot(cerr);

        SECTION( "One cluster taking loop" ) {
 
            id_t seed_nodes[] = {1, 4};
            //all are in the same cluster
            vector<SnarlSeedClusterer::Seed> seeds;
            for (id_t n : seed_nodes) {
                pos_t pos = make_pos_t(n, false, 0);
                std::tuple<bool, size_t, size_t, bool, size_t, size_t, size_t, size_t, bool> chain_info = dist_index.get_minimizer_distances(pos);
                seeds.push_back({ pos, 0, std::get<0>(chain_info), std::get<1>(chain_info), std::get<2>(chain_info),
                   std::get<3>(chain_info), std::get<4>(chain_info), std::get<5>(chain_info), std::get<6>(chain_info), std::get<7>(chain_info), std::get<8>(chain_info)});
            }

            vector<SnarlSeedClusterer::Cluster> clusters = clusterer.cluster_seeds(seeds, 6); 
            REQUIRE(clusters.size() == 1); 

        }
        SECTION( "One cluster on boundary" ) {
 
            id_t seed_nodes[] = {2, 4};
            //all are in the same cluster
            vector<SnarlSeedClusterer::Seed> seeds;
            for (id_t n : seed_nodes) {
                pos_t pos = make_pos_t(n, false, 0);
                std::tuple<bool, size_t, size_t, bool, size_t, size_t, size_t, size_t, bool> chain_info = dist_index.get_minimizer_distances(pos);
                seeds.push_back({ pos, 0, std::get<0>(chain_info), std::get<1>(chain_info), std::get<2>(chain_info),
                   std::get<3>(chain_info), std::get<4>(chain_info), std::get<5>(chain_info), std::get<6>(chain_info), std::get<7>(chain_info), std::get<8>(chain_info)});
            }

            vector<SnarlSeedClusterer::Cluster> clusters = clusterer.cluster_seeds(seeds, 3); 
            REQUIRE(clusters.size() == 1); 

        }
        SECTION( "One fragment cluster on boundary" ) {
 
            id_t seed_nodes[] = {2, 4};
            //all are in the same cluster
            vector<vector<SnarlSeedClusterer::Seed>> seeds (2);

            pos_t pos = make_pos_t(2, false, 0);
            std::tuple<bool, size_t, size_t, bool, size_t, size_t, size_t, size_t, bool> chain_info = dist_index.get_minimizer_distances(pos);
            seeds[0].push_back({ pos, 0, std::get<0>(chain_info), std::get<1>(chain_info), std::get<2>(chain_info),
               std::get<3>(chain_info), std::get<4>(chain_info), std::get<5>(chain_info), std::get<6>(chain_info), std::get<7>(chain_info), std::get<8>(chain_info)});

            pos = make_pos_t(4, false, 0);
            chain_info = dist_index.get_minimizer_distances(pos);
            seeds[1].push_back({ pos, 0, std::get<0>(chain_info), std::get<1>(chain_info), std::get<2>(chain_info),
               std::get<3>(chain_info), std::get<4>(chain_info), std::get<5>(chain_info), std::get<6>(chain_info), std::get<7>(chain_info), std::get<8>(chain_info)});

            vector<vector<SnarlSeedClusterer::Cluster>> clusters = clusterer.cluster_seeds(seeds, 3, 3); 
            REQUIRE(clusters.size() == 2); 
            REQUIRE(clusters[0][0].fragment == clusters[1][0].fragment);

        }
        SECTION( "One cluster on boundary" ) {
 
            id_t seed_nodes[] = {3, 4};
            //all are in the same cluster
            vector<SnarlSeedClusterer::Seed> seeds;
            for (id_t n : seed_nodes) {
                pos_t pos = make_pos_t(n, false, 0);
                std::tuple<bool, size_t, size_t, bool, size_t, size_t, size_t, size_t, bool> chain_info = dist_index.get_minimizer_distances(pos);
                seeds.push_back({ pos, 0, std::get<0>(chain_info), std::get<1>(chain_info), std::get<2>(chain_info),
                   std::get<3>(chain_info), std::get<4>(chain_info), std::get<5>(chain_info), std::get<6>(chain_info), std::get<7>(chain_info), std::get<8>(chain_info)});
            }

            vector<SnarlSeedClusterer::Cluster> clusters = clusterer.cluster_seeds(seeds, 3); 
            REQUIRE(clusters.size() == 1); 

        }
    }
    TEST_CASE( "chain with loop",
                   "[cluster]" ) {
        VG graph;

        Node* n1 = graph.create_node("GCA");
        Node* n2 = graph.create_node("T");
        Node* n3 = graph.create_node("G");
        Node* n4 = graph.create_node("CTGA");
        Node* n5 = graph.create_node("GCA");
        Node* n6 = graph.create_node("T");
        Node* n7 = graph.create_node("G");
        Node* n8 = graph.create_node("CTGA");

        Edge* e1 = graph.create_edge(n1, n2);
        Edge* e2 = graph.create_edge(n1, n3);
        Edge* e3 = graph.create_edge(n2, n3);
        Edge* e4 = graph.create_edge(n2, n4);
        Edge* e5 = graph.create_edge(n3, n4);
        Edge* e6 = graph.create_edge(n3, n5);
        Edge* e7 = graph.create_edge(n4, n5);
        Edge* e8 = graph.create_edge(n4, n6);
        Edge* e9 = graph.create_edge(n5, n6);
        Edge* e10 = graph.create_edge(n6, n7);
        Edge* e11 = graph.create_edge(n6, n7, false, true);
        Edge* e12 = graph.create_edge(n6, n8);
        Edge* e13 = graph.create_edge(n7, n8);

        CactusSnarlFinder bubble_finder(graph);
        SnarlManager snarl_manager = bubble_finder.find_snarls();

        MinimumDistanceIndex dist_index (&graph, &snarl_manager);
        SnarlSeedClusterer clusterer(dist_index);
        
        //graph.to_dot(cerr);

        SECTION( "One cluster taking loop" ) {
 
            id_t seed_nodes[] = {4, 5};
            //all are in the same cluster
            vector<SnarlSeedClusterer::Seed> seeds;
            for (id_t n : seed_nodes) {
                pos_t pos = make_pos_t(n, false, 0);
                std::tuple<bool, size_t, size_t, bool, size_t, size_t, size_t, size_t, bool> chain_info = dist_index.get_minimizer_distances(pos);
                seeds.push_back({ pos, 0, std::get<0>(chain_info), std::get<1>(chain_info), std::get<2>(chain_info),
                   std::get<3>(chain_info), std::get<4>(chain_info), std::get<5>(chain_info), std::get<6>(chain_info), std::get<7>(chain_info), std::get<8>(chain_info)});
            }

            vector<SnarlSeedClusterer::Cluster> clusters = clusterer.cluster_seeds(seeds, 11); 
            REQUIRE(clusters.size() == 1); 

        }
        SECTION( "One cluster not taking loop" ) {
 
            id_t seed_nodes[] = {4, 5, 3};
            //all are in the same cluster
            vector<SnarlSeedClusterer::Seed> seeds;
            for (id_t n : seed_nodes) {
                pos_t pos = make_pos_t(n, false, 0);
                std::tuple<bool, size_t, size_t, bool, size_t, size_t, size_t, size_t, bool> chain_info = dist_index.get_minimizer_distances(pos);
                seeds.push_back({ pos, 0, std::get<0>(chain_info), std::get<1>(chain_info), std::get<2>(chain_info),
                   std::get<3>(chain_info), std::get<4>(chain_info), std::get<5>(chain_info), std::get<6>(chain_info), std::get<7>(chain_info), std::get<8>(chain_info)});
            }

            vector<SnarlSeedClusterer::Cluster> clusters = clusterer.cluster_seeds(seeds, 3); 
            REQUIRE(clusters.size() == 1); 

        }
        SECTION( "One cluster not taking loop" ) {
 
            id_t seed_nodes[] = {4, 5, 6};
            //all are in the same cluster
            vector<SnarlSeedClusterer::Seed> seeds;
            for (id_t n : seed_nodes) {
                pos_t pos = make_pos_t(n, false, 0);
                std::tuple<bool, size_t, size_t, bool, size_t, size_t, size_t, size_t, bool> chain_info = dist_index.get_minimizer_distances(pos);
                seeds.push_back({ pos, 0, std::get<0>(chain_info), std::get<1>(chain_info), std::get<2>(chain_info),
                   std::get<3>(chain_info), std::get<4>(chain_info), std::get<5>(chain_info), std::get<6>(chain_info), std::get<7>(chain_info), std::get<8>(chain_info)});
            }

            vector<SnarlSeedClusterer::Cluster> clusters = clusterer.cluster_seeds(seeds, 8); 
            REQUIRE(clusters.size() == 1); 

        }
        SECTION( "Two clusters" ) {
 
            id_t seed_nodes[] = {4, 5, 1};
            //all are in the same cluster
            vector<SnarlSeedClusterer::Seed> seeds;
            for (id_t n : seed_nodes) {
                pos_t pos = make_pos_t(n, false, 0);
                std::tuple<bool, size_t, size_t, bool, size_t, size_t, size_t, size_t, bool> chain_info = dist_index.get_minimizer_distances(pos);
                seeds.push_back({ pos, 0, std::get<0>(chain_info), std::get<1>(chain_info), std::get<2>(chain_info),
                   std::get<3>(chain_info), std::get<4>(chain_info), std::get<5>(chain_info), std::get<6>(chain_info), std::get<7>(chain_info), std::get<8>(chain_info)});
            }

            vector<SnarlSeedClusterer::Cluster> clusters = clusterer.cluster_seeds(seeds, 3); 
            REQUIRE(clusters.size() == 3); 

        }
    }
    TEST_CASE( "multiple clusters in a chain",
                   "[cluster]" ) {
        VG graph;

        Node* n1 = graph.create_node("GCA");
        Node* n2 = graph.create_node("T");
        Node* n3 = graph.create_node("G");
        Node* n4 = graph.create_node("CTGA");
        Node* n5 = graph.create_node("GCA");
        Node* n6 = graph.create_node("T");
        Node* n7 = graph.create_node("G");
        Node* n8 = graph.create_node("CTGA");
        Node* n9 = graph.create_node("GCA");
        Node* n10 = graph.create_node("T");
        Node* n11 = graph.create_node("G");
        Node* n12 = graph.create_node("CTGA");
        Node* n13 = graph.create_node("GCA");

        Edge* e1 = graph.create_edge(n1, n2);
        Edge* e2 = graph.create_edge(n1, n9);
        Edge* e3 = graph.create_edge(n2, n3);
        Edge* e4 = graph.create_edge(n2, n4);
        Edge* e5 = graph.create_edge(n3, n4);
        Edge* e6 = graph.create_edge(n4, n5);
        Edge* e7 = graph.create_edge(n4, n5, false, true);
        Edge* e8 = graph.create_edge(n5, n6);
        Edge* e9 = graph.create_edge(n5, n6, true, false);
        Edge* e10 = graph.create_edge(n6, n7);
        Edge* e11 = graph.create_edge(n6, n8);
        Edge* e12 = graph.create_edge(n7, n8);
        Edge* e13 = graph.create_edge(n8, n10);
        Edge* e14 = graph.create_edge(n9, n10);
        Edge* e15 = graph.create_edge(n10, n11);
        Edge* e16 = graph.create_edge(n10, n12);
        Edge* e17 = graph.create_edge(n11, n13);
        Edge* e18 = graph.create_edge(n12, n13);

        CactusSnarlFinder bubble_finder(graph);
        SnarlManager snarl_manager = bubble_finder.find_snarls();

        const Snarl* snarl1 = snarl_manager.into_which_snarl(1, false);
        const Snarl* snarl2 = snarl_manager.into_which_snarl(2, false);
        const Snarl* snarl3 = snarl_manager.into_which_snarl(3, false);

        MinimumDistanceIndex dist_index (&graph, &snarl_manager);
        SnarlSeedClusterer clusterer(dist_index);
        
        //graph.to_dot(cerr);

        SECTION( "One cluster with seed struct" ) {
 
            id_t seed_nodes[] = {2, 3, 4, 7, 8, 9, 11};
            //all are in the same cluster
            vector<SnarlSeedClusterer::Seed> seeds;
            for (id_t n : seed_nodes) {
                pos_t pos = make_pos_t(n, false, 0);
                std::tuple<bool, size_t, size_t, bool, size_t, size_t, size_t, size_t, bool> chain_info = dist_index.get_minimizer_distances(pos);
                seeds.push_back({ pos, 0, std::get<0>(chain_info), std::get<1>(chain_info), std::get<2>(chain_info),
                   std::get<3>(chain_info), std::get<4>(chain_info), std::get<5>(chain_info), std::get<6>(chain_info), std::get<7>(chain_info), std::get<8>(chain_info)});
            }

            vector<SnarlSeedClusterer::Cluster> clusters = clusterer.cluster_seeds(seeds, 10); 
            REQUIRE(clusters.size() == 1); 

        }
        SECTION( "Check distance index assignments") {

            REQUIRE (std::get<0>(dist_index.get_minimizer_distances(make_pos_t(1 , false, 0))));
            REQUIRE (!std::get<0>(dist_index.get_minimizer_distances(make_pos_t(2 , false, 0))));
            REQUIRE (!std::get<0>(dist_index.get_minimizer_distances(make_pos_t(3 , false, 0))));
            REQUIRE (!std::get<0>(dist_index.get_minimizer_distances(make_pos_t(4 , false, 0))));
            REQUIRE (!std::get<0>(dist_index.get_minimizer_distances(make_pos_t(5 , false, 0))));
            REQUIRE (!std::get<0>(dist_index.get_minimizer_distances(make_pos_t(6 , false, 0))));
            REQUIRE (!std::get<0>(dist_index.get_minimizer_distances(make_pos_t(7 , false, 0))));
            REQUIRE (!std::get<0>(dist_index.get_minimizer_distances(make_pos_t(8 , false, 0))));
            REQUIRE (!std::get<0>(dist_index.get_minimizer_distances(make_pos_t(9 , false, 0))));
            REQUIRE (std::get<0>(dist_index.get_minimizer_distances(make_pos_t(10 , false, 0))));
            REQUIRE (!std::get<0>(dist_index.get_minimizer_distances(make_pos_t(11 , false, 0))));
            REQUIRE (!std::get<0>(dist_index.get_minimizer_distances(make_pos_t(12 , false, 0))));
            REQUIRE (std::get<0>(dist_index.get_minimizer_distances(make_pos_t(13 , false, 0))));

            REQUIRE (std::get<3>(dist_index.get_minimizer_distances(make_pos_t(11 , false, 0))));
            REQUIRE (std::get<3>(dist_index.get_minimizer_distances(make_pos_t(12 , false, 0))));
            REQUIRE (std::get<4>(dist_index.get_minimizer_distances(make_pos_t(11 , false, 0))) == 1);
            REQUIRE (std::get<5>(dist_index.get_minimizer_distances(make_pos_t(11 , false, 0))) == 1);
            REQUIRE (std::get<6>(dist_index.get_minimizer_distances(make_pos_t(11 , false, 0))) == 3);
            REQUIRE (std::get<7>(dist_index.get_minimizer_distances(make_pos_t(11 , false, 0))) == 1);
        }
        SECTION( "Two clusters" ) {
 
            vector<id_t> seed_nodes( {2, 3, 4, 7, 8, 10, 11});
            //Clusters should be {2, 3, 4}, {7, 8, 10, 11}
            //Distance from pos on 4 to pos on 7 is 8, including one position
            vector<SnarlSeedClusterer::Seed> seeds;
            for (id_t n : seed_nodes) {
                pos_t pos = make_pos_t(n, false, 0);
                std::tuple<bool, size_t, size_t, bool, size_t, size_t, size_t, size_t, bool> chain_info = dist_index.get_minimizer_distances(pos);
                seeds.push_back({pos, 0, std::get<0>(chain_info), std::get<1>(chain_info), std::get<2>(chain_info),
                   std::get<3>(chain_info), std::get<4>(chain_info), std::get<5>(chain_info), std::get<6>(chain_info), std::get<7>(chain_info), std::get<8>(chain_info)});
            }


            vector<SnarlSeedClusterer::Cluster> clusters = clusterer.cluster_seeds(seeds, 7); 
            vector<hash_set<size_t>> cluster_sets;
            for (auto& c : clusters) {
                hash_set<size_t> h;
                for (size_t s : c.seeds) {
                    h.insert(s);
                }
                cluster_sets.push_back(h);
            }
            REQUIRE( clusters.size() == 2);
            REQUIRE (( (cluster_sets[0].count(0) == 1 &&
                       cluster_sets[0].count(1) == 1 &&
                       cluster_sets[0].count(2) == 1 &&
                       cluster_sets[1].count(3) == 1 &&
                       cluster_sets[1].count(4) == 1 &&
                       cluster_sets[1].count(5) == 1 &&
                       cluster_sets[1].count(6) == 1  ) ||

                     ( cluster_sets[1].count(0) == 1 &&
                       cluster_sets[1].count(1) == 1 &&
                       cluster_sets[1].count(2) == 1 &&
                       cluster_sets[0].count(3) == 1 &&
                       cluster_sets[0].count(4) == 1 &&
                       cluster_sets[0].count(5) == 1 &&
                       cluster_sets[0].count(6) == 1  )));

        }
        SECTION( "One fragment cluster of the same node" ) {
 
            vector<id_t> seed_nodes( {2, 3});
            vector<id_t> seed_nodes1({2, 7, 8, 10, 11});
            //Clusters should be {2, 3, 4}, {2}, {7, 8, 10, 11}
            //One fragment cluster
            //Distance from pos on 4 to pos on 7 is 8, including one position
            vector<SnarlSeedClusterer::Seed> seeds ;
            for (id_t n : seed_nodes) {
                pos_t pos = make_pos_t(n, false, 0);
                std::tuple<bool, size_t, size_t, bool, size_t, size_t, size_t, size_t, bool> chain_info = dist_index.get_minimizer_distances(pos);
                seeds.push_back({ pos, 0, std::get<0>(chain_info), std::get<1>(chain_info), std::get<2>(chain_info),
                   std::get<3>(chain_info), std::get<4>(chain_info), std::get<5>(chain_info), std::get<6>(chain_info), std::get<7>(chain_info), std::get<8>(chain_info)});
            }
            vector<SnarlSeedClusterer::Seed> seeds1;
            for (id_t n : seed_nodes1) {
                pos_t pos = make_pos_t(n, false, 0);
                std::tuple<bool, size_t, size_t, bool, size_t, size_t, size_t, size_t, bool> chain_info = dist_index.get_minimizer_distances(pos);
                seeds1.push_back({ pos, 0, std::get<0>(chain_info), std::get<1>(chain_info), std::get<2>(chain_info),
                   std::get<3>(chain_info), std::get<4>(chain_info), std::get<5>(chain_info), std::get<6>(chain_info), std::get<7>(chain_info), std::get<8>(chain_info)});
            }
            vector<vector<SnarlSeedClusterer::Seed>> all_seeds;
            all_seeds.push_back(seeds);
            all_seeds.push_back(seeds1);


            vector<vector<SnarlSeedClusterer::Cluster>> paired_clusters = clusterer.cluster_seeds(all_seeds, 7, 15); 
            //Should be [[<[0,1,2], 0>],[<[3,4,5,6], 0>]] 
            REQUIRE( paired_clusters.size() == 2);
            REQUIRE( paired_clusters[0].size() == 1);
            REQUIRE( paired_clusters[1].size() == 2);
            REQUIRE( paired_clusters[0][0].fragment == paired_clusters[1][0].fragment);
            REQUIRE( paired_clusters[1][0].fragment == paired_clusters[1][1].fragment);
        }
        SECTION( "One fragment cluster" ) {
 
            vector<id_t> seed_nodes( {2, 3, 4});
            vector<id_t> seed_nodes1({7, 8, 10, 11});
            //Clusters should be {2, 3, 4}, {7, 8, 10, 11}
            //One fragment cluster
            //Distance from pos on 4 to pos on 7 is 8, including one position
            vector<SnarlSeedClusterer::Seed> seeds ;
            for (id_t n : seed_nodes) {
                pos_t pos = make_pos_t(n, false, 0);
                std::tuple<bool, size_t, size_t, bool, size_t, size_t, size_t, size_t, bool> chain_info = dist_index.get_minimizer_distances(pos);
                seeds.push_back({ pos, 0, std::get<0>(chain_info), std::get<1>(chain_info), std::get<2>(chain_info),
                   std::get<3>(chain_info), std::get<4>(chain_info), std::get<5>(chain_info), std::get<6>(chain_info), std::get<7>(chain_info), std::get<8>(chain_info)});
            }
            vector<SnarlSeedClusterer::Seed> seeds1;
            for (id_t n : seed_nodes1) {
                pos_t pos = make_pos_t(n, false, 0);
                std::tuple<bool, size_t, size_t, bool, size_t, size_t, size_t, size_t, bool> chain_info = dist_index.get_minimizer_distances(pos);
                seeds1.push_back({ pos, 0, std::get<0>(chain_info), std::get<1>(chain_info), std::get<2>(chain_info),
                   std::get<3>(chain_info), std::get<4>(chain_info), std::get<5>(chain_info), std::get<6>(chain_info), std::get<7>(chain_info), std::get<8>(chain_info)});
            }
            vector<vector<SnarlSeedClusterer::Seed>> all_seeds;
            all_seeds.push_back(seeds);
            all_seeds.push_back(seeds1);


            vector<vector<SnarlSeedClusterer::Cluster>> paired_clusters = clusterer.cluster_seeds(all_seeds, 7, 15); 
            //Should be [[<[0,1,2], 0>],[<[3,4,5,6], 0>]] 
            REQUIRE( paired_clusters.size() == 2);
            REQUIRE( paired_clusters[0].size() == 1);
            REQUIRE( paired_clusters[1].size() == 1);
            REQUIRE( paired_clusters[0][0].seeds.size() == 3);
            REQUIRE( paired_clusters[1][0].seeds.size() == 4);
            REQUIRE( paired_clusters[0][0].fragment == paired_clusters[1][0].fragment);
        }
        SECTION( "Two fragment clusters with seed structs" ) {
 
            vector<id_t> seed_nodes( {2, 3, 4});
            vector<id_t> seed_nodes1({7, 8, 10, 11});
            //Fragment clusters should be {2, 3, 4}, {7, 8, 10, 11}
            //Distance from pos on 4 to pos on 7 is 8, including one position
            vector<SnarlSeedClusterer::Seed> seeds;
            for (id_t n : seed_nodes) {
                pos_t pos = make_pos_t(n, false, 0);
                std::tuple<bool, size_t, size_t, bool, size_t, size_t, size_t, size_t, bool> chain_info = dist_index.get_minimizer_distances(pos);
                seeds.push_back({ pos, 0, std::get<0>(chain_info), std::get<1>(chain_info), std::get<2>(chain_info),
                   std::get<3>(chain_info), std::get<4>(chain_info), std::get<5>(chain_info), std::get<6>(chain_info), std::get<7>(chain_info), std::get<8>(chain_info)});
            }
            vector<SnarlSeedClusterer::Seed> seeds1;
            for (id_t n : seed_nodes1) {
                pos_t pos = make_pos_t(n, false, 0);
                std::tuple<bool, size_t, size_t, bool, size_t, size_t, size_t, size_t, bool> chain_info = dist_index.get_minimizer_distances(pos);
                seeds1.push_back({ pos, 0, std::get<0>(chain_info), std::get<1>(chain_info), std::get<2>(chain_info),
                   std::get<3>(chain_info), std::get<4>(chain_info), std::get<5>(chain_info), std::get<6>(chain_info), std::get<7>(chain_info), std::get<8>(chain_info)});
            }
            vector<vector<SnarlSeedClusterer::Seed>> all_seeds;
            all_seeds.push_back(seeds);
            all_seeds.push_back(seeds1);


            vector<vector<SnarlSeedClusterer::Cluster>> paired_clusters = clusterer.cluster_seeds(all_seeds, 2, 7); 
            // read_clusters = [ [[0,1,2]],[[3,4],[5,6]] ]
            // fragment_clusters = [ [0,1,2], [3,4,5,6] ]
            REQUIRE( paired_clusters.size() == 2) ;
            REQUIRE( paired_clusters[0].size() == 1);
            REQUIRE( paired_clusters[1].size() == 2);
            REQUIRE( paired_clusters[0][0].fragment != paired_clusters[1][0].fragment);
            REQUIRE( paired_clusters[0][0].fragment != paired_clusters[1][1].fragment);
            REQUIRE( paired_clusters[1][0].fragment == paired_clusters[1][1].fragment);

        }
        SECTION( "Two fragment clusters" ) {
 
            vector<id_t> seed_nodes( {2, 3, 4});
            vector<id_t> seed_nodes1({7, 8, 10, 11});
            //Fragment clusters should be {2, 3, 4}, {7, 8, 10, 11}
            //Distance from pos on 4 to pos on 7 is 8, including one position
            vector<SnarlSeedClusterer::Seed> seeds ;
            for (id_t n : seed_nodes) {
                pos_t pos = make_pos_t(n, false, 0);
                std::tuple<bool, size_t, size_t, bool, size_t, size_t, size_t, size_t, bool> chain_info = dist_index.get_minimizer_distances(pos);
                seeds.push_back({ pos, 0, std::get<0>(chain_info), std::get<1>(chain_info), std::get<2>(chain_info),
                   std::get<3>(chain_info), std::get<4>(chain_info), std::get<5>(chain_info), std::get<6>(chain_info), std::get<7>(chain_info), std::get<8>(chain_info)});
            }
            vector<SnarlSeedClusterer::Seed> seeds1;
            for (id_t n : seed_nodes1) {
                pos_t pos = make_pos_t(n, false, 0);
                std::tuple<bool, size_t, size_t, bool, size_t, size_t, size_t, size_t, bool> chain_info = dist_index.get_minimizer_distances(pos);
                seeds1.push_back({ pos, 0, std::get<0>(chain_info), std::get<1>(chain_info), std::get<2>(chain_info),
                   std::get<3>(chain_info), std::get<4>(chain_info), std::get<5>(chain_info), std::get<6>(chain_info), std::get<7>(chain_info), std::get<8>(chain_info)});
            }
            vector<vector<SnarlSeedClusterer::Seed>> all_seeds;
            all_seeds.push_back(seeds);
            all_seeds.push_back(seeds1);


            vector<vector<SnarlSeedClusterer::Cluster>> paired_clusters = clusterer.cluster_seeds(all_seeds, 2, 7); 
            // read_clusters = [ [[0,1,2]],[[3,4],[5,6]] ]
            // fragment_clusters = [ [0,1,2], [3,4,5,6] ]
            REQUIRE( paired_clusters.size() == 2) ;
            REQUIRE( paired_clusters[0].size() == 1);
            REQUIRE( paired_clusters[1].size() == 2);
            REQUIRE( paired_clusters[0][0].fragment != paired_clusters[1][0].fragment);
            REQUIRE( paired_clusters[0][0].fragment != paired_clusters[1][1].fragment);
            REQUIRE( paired_clusters[1][0].fragment == paired_clusters[1][1].fragment);

        }
    }//End test case

    TEST_CASE( "Revese in chain right","[cluster]" ) {
        VG graph;

        Node* n1 = graph.create_node("GCA");
        Node* n2 = graph.create_node("T");
        Node* n3 = graph.create_node("G");
        Node* n4 = graph.create_node("CTGA");
        Node* n5 = graph.create_node("G");
        Node* n6 = graph.create_node("T");
        Node* n7 = graph.create_node("G");
        Node* n8 = graph.create_node("G");
        Node* n9 = graph.create_node("AA");
        Node* n10 = graph.create_node("G");
        Node* n11 = graph.create_node("GGGGGGGGGG");//10

        Edge* e1 = graph.create_edge(n1, n2);
        Edge* e2 = graph.create_edge(n1, n10);
        Edge* e3 = graph.create_edge(n2, n3);
        Edge* e4 = graph.create_edge(n2, n4);
        Edge* e5 = graph.create_edge(n3, n5);
        Edge* e6 = graph.create_edge(n4, n5);
        Edge* e7 = graph.create_edge(n5, n6);
        Edge* e8 = graph.create_edge(n5, n11);
        Edge* e9 = graph.create_edge(n11, n7);
        Edge* e10 = graph.create_edge(n6, n7);
        Edge* e11 = graph.create_edge(n8, n8, false, true);
        Edge* e12 = graph.create_edge(n7, n8);
        Edge* e13 = graph.create_edge(n7, n9);
        Edge* e14 = graph.create_edge(n8, n9);
        Edge* e15 = graph.create_edge(n9, n10);

        CactusSnarlFinder bubble_finder(graph);
        SnarlManager snarl_manager = bubble_finder.find_snarls();
        MinimumDistanceIndex dist_index (&graph, &snarl_manager);

        SnarlSeedClusterer clusterer(dist_index);

        SECTION( "Same snarl" ) {
            vector<id_t> seed_nodes ({3, 4});
            vector<SnarlSeedClusterer::Seed> seeds;
            for (id_t n : seed_nodes) {
                pos_t pos = make_pos_t(n, false, 0);
                std::tuple<bool, size_t, size_t, bool, size_t, size_t, size_t, size_t, bool> chain_info = dist_index.get_minimizer_distances(pos);
                seeds.push_back({ pos, 0, std::get<0>(chain_info), std::get<1>(chain_info), std::get<2>(chain_info),
                   std::get<3>(chain_info), std::get<4>(chain_info), std::get<5>(chain_info), std::get<6>(chain_info), std::get<7>(chain_info), std::get<8>(chain_info)});
            }


            vector<SnarlSeedClusterer::Cluster> clusters = clusterer.cluster_seeds(seeds, 13); 

            REQUIRE( clusters.size() == 1);
        }
        SECTION( "Different snarl" ) {
            vector<SnarlSeedClusterer::Seed> seeds;

            vector<pos_t> pos_ts;
            pos_ts.emplace_back(3, false, 0);
            pos_ts.emplace_back(11, false, 9);
            for (pos_t pos : pos_ts) {
                std::tuple<bool, size_t, size_t, bool, size_t, size_t, size_t, size_t, bool> chain_info = dist_index.get_minimizer_distances(pos);
                seeds.push_back({ pos, 0, std::get<0>(chain_info), std::get<1>(chain_info), std::get<2>(chain_info),
                   std::get<3>(chain_info), std::get<4>(chain_info), std::get<5>(chain_info), std::get<6>(chain_info), std::get<7>(chain_info), std::get<8>(chain_info)});
            }



            vector<SnarlSeedClusterer::Cluster> clusters = clusterer.cluster_seeds(seeds, 8); 


            REQUIRE( clusters.size() == 1);
        }
    }//end test case
    TEST_CASE( "Reverse in chain left","[cluster]" ) {
        VG graph;

        Node* n1 = graph.create_node("GCA");
        Node* n2 = graph.create_node("T");
        Node* n3 = graph.create_node("G");
        Node* n4 = graph.create_node("CTGA");
        Node* n5 = graph.create_node("G");
        Node* n6 = graph.create_node("TGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG");
        Node* n7 = graph.create_node("G");
        Node* n8 = graph.create_node("G");
        Node* n9 = graph.create_node("AA");
        Node* n10 = graph.create_node("G");
        Node* n11 = graph.create_node("G");

        Edge* e1 = graph.create_edge(n1, n2);
        Edge* e2 = graph.create_edge(n1, n10);
        Edge* e3 = graph.create_edge(n2, n3);
        Edge* e4 = graph.create_edge(n2, n4);
        Edge* e5 = graph.create_edge(n3, n5);
        Edge* e6 = graph.create_edge(n4, n5);
        Edge* e7 = graph.create_edge(n5, n6);
        Edge* e8 = graph.create_edge(n5, n7);
        Edge* e9 = graph.create_edge(n6, n7);
        Edge* e10 = graph.create_edge(n7, n8);
        Edge* e11 = graph.create_edge(n7, n9);
        Edge* e12 = graph.create_edge(n8, n9);
        Edge* e13 = graph.create_edge(n9, n10);
        Edge* e14 = graph.create_edge(n11, n5);
        Edge* e15 = graph.create_edge(n11, n5, true, false);

        CactusSnarlFinder bubble_finder(graph);
        SnarlManager snarl_manager = bubble_finder.find_snarls();
        MinimumDistanceIndex dist_index (&graph, &snarl_manager);

        SnarlSeedClusterer clusterer(dist_index);

        SECTION( "One cluster" ) {
            vector<id_t> seed_nodes ({7, 7, 6});
            vector<SnarlSeedClusterer::Seed> seeds;
            for (id_t n : seed_nodes) {
                pos_t pos = make_pos_t(n, false, 0);
                std::tuple<bool, size_t, size_t, bool, size_t, size_t, size_t, size_t, bool> chain_info = dist_index.get_minimizer_distances(pos);
                seeds.push_back({ pos, 0, std::get<0>(chain_info), std::get<1>(chain_info), std::get<2>(chain_info),
                   std::get<3>(chain_info), std::get<4>(chain_info), std::get<5>(chain_info), std::get<6>(chain_info), std::get<7>(chain_info), std::get<8>(chain_info)});
            }

            vector<SnarlSeedClusterer::Cluster> clusters =  clusterer.cluster_seeds(seeds, 20); 


            REQUIRE( clusters.size() == 1);
        }
        SECTION( "two clusters" ) {
            vector<id_t> seed_nodes ({2, 6});
            vector<SnarlSeedClusterer::Seed> seeds;
            for (id_t n : seed_nodes) {
                pos_t pos = make_pos_t(n, false, 0);
                std::tuple<bool, size_t, size_t, bool, size_t, size_t, size_t, size_t, bool> chain_info = dist_index.get_minimizer_distances(pos);
                seeds.push_back({ pos, 0, std::get<0>(chain_info), std::get<1>(chain_info), std::get<2>(chain_info),
                   std::get<3>(chain_info), std::get<4>(chain_info), std::get<5>(chain_info), std::get<6>(chain_info), std::get<7>(chain_info), std::get<8>(chain_info)});
            }

            vector<SnarlSeedClusterer::Cluster> clusters =  clusterer.cluster_seeds(seeds, 20); 


        }
        SECTION( "different snarl" ) {
            vector<id_t> seed_nodes ({8, 6});
            vector<SnarlSeedClusterer::Seed> seeds;
            for (id_t n : seed_nodes) {
                pos_t pos = make_pos_t(n, false, 0);
                std::tuple<bool, size_t, size_t, bool, size_t, size_t, size_t, size_t, bool> chain_info = dist_index.get_minimizer_distances(pos);
                seeds.push_back({ pos, 0, std::get<0>(chain_info), std::get<1>(chain_info), std::get<2>(chain_info),
                   std::get<3>(chain_info), std::get<4>(chain_info), std::get<5>(chain_info), std::get<6>(chain_info), std::get<7>(chain_info), std::get<8>(chain_info)});
            }

            vector<SnarlSeedClusterer::Cluster> clusters =  clusterer.cluster_seeds(seeds, 20); 


            REQUIRE( clusters.size() == 1);
        }
    }//end test case


    TEST_CASE( "Loop on node","[cluster]" ) {
        VG graph;

        Node* n1 = graph.create_node("GCA");
        Node* n2 = graph.create_node("T");
        Node* n3 = graph.create_node("G");
        Node* n4 = graph.create_node("CTGA");
        Node* n5 = graph.create_node("GGGGGGGGGGGG");//12 Gs
        Node* n6 = graph.create_node("T");

        Edge* e1 = graph.create_edge(n1, n2);
        Edge* e2 = graph.create_edge(n1, n3);
        Edge* e3 = graph.create_edge(n2, n4);
        Edge* e4 = graph.create_edge(n3, n4);
        Edge* e27 = graph.create_edge(n4, n5);
        Edge* e5 = graph.create_edge(n4, n6);
        Edge* e6 = graph.create_edge(n5, n6);
        Edge* e7 = graph.create_edge(n5, n5);

        CactusSnarlFinder bubble_finder(graph);
        SnarlManager snarl_manager = bubble_finder.find_snarls();
        MinimumDistanceIndex dist_index (&graph, &snarl_manager);

        SnarlSeedClusterer clusterer(dist_index);

        SECTION( "One cluster taking node loop" ) {
            vector<SnarlSeedClusterer::Seed> seeds;
            vector<pos_t> pos_ts;
            pos_ts.emplace_back(5, false, 0);
            pos_ts.emplace_back(5, true, 0);

            for (pos_t pos : pos_ts){

                std::tuple<bool, size_t, size_t, bool, size_t, size_t, size_t, size_t, bool> chain_info = dist_index.get_minimizer_distances(pos);
                seeds.push_back({ pos, 0, std::get<0>(chain_info), std::get<1>(chain_info), std::get<2>(chain_info),
                   std::get<3>(chain_info), std::get<4>(chain_info), std::get<5>(chain_info), std::get<6>(chain_info), std::get<7>(chain_info), std::get<8>(chain_info)});
            }
            vector<SnarlSeedClusterer::Cluster> clusters =  clusterer.cluster_seeds(seeds, 3); 

            REQUIRE( clusters.size() == 1);
        }
    }
    TEST_CASE( "Clusters in snarl","[cluster]" ) {
        VG graph;

        Node* n1 = graph.create_node("GCA");
        Node* n2 = graph.create_node("T");
        Node* n3 = graph.create_node("G");
        Node* n4 = graph.create_node("CTGA");
        Node* n5 = graph.create_node("GGGGGGGGGGGG");//12 Gs
        Node* n6 = graph.create_node("T");
        Node* n7 = graph.create_node("G");
        Node* n8 = graph.create_node("G");
        Node* n9 = graph.create_node("AA");
        Node* n10 = graph.create_node("G");
        Node* n11 = graph.create_node("G");
        Node* n12 = graph.create_node("G");
        Node* n13 = graph.create_node("GA");
        Node* n14 = graph.create_node("G");
        Node* n15 = graph.create_node("G");
        Node* n16 = graph.create_node("G");

        Edge* e1 = graph.create_edge(n1, n2);
        Edge* e2 = graph.create_edge(n1, n13);
        Edge* e3 = graph.create_edge(n2, n3);
        Edge* e4 = graph.create_edge(n2, n16);
        Edge* e27 = graph.create_edge(n16, n9);
        Edge* e5 = graph.create_edge(n3, n4);
        Edge* e6 = graph.create_edge(n3, n5);
        Edge* e7 = graph.create_edge(n4, n6);
        Edge* e8 = graph.create_edge(n5, n6);
        Edge* e9 = graph.create_edge(n6, n7);
        Edge* e10 = graph.create_edge(n6, n8);
        Edge* e11 = graph.create_edge(n7, n8);
        Edge* e12 = graph.create_edge(n8, n9);
        Edge* e13 = graph.create_edge(n9, n10);
        Edge* e14 = graph.create_edge(n9, n11);
        Edge* e15 = graph.create_edge(n10, n11);
        Edge* e16 = graph.create_edge(n11, n12);
        Edge* e17 = graph.create_edge(n11, n2);
        Edge* e18 = graph.create_edge(n12, n1);
        Edge* e19 = graph.create_edge(n13, n14);
        Edge* e20 = graph.create_edge(n13, n15);
        Edge* e21 = graph.create_edge(n14, n15);
        Edge* e22 = graph.create_edge(n15, n12);
        Edge* e23 = graph.create_edge(n2, n2, true, false);
        Edge* e24 = graph.create_edge(n11, n11, false, true);
        Edge* e25 = graph.create_edge(n1, n1, true, false);
        Edge* e26 = graph.create_edge(n12, n12, false, true);

        CactusSnarlFinder bubble_finder(graph);
        SnarlManager snarl_manager = bubble_finder.find_snarls();
        MinimumDistanceIndex dist_index (&graph, &snarl_manager);

        SnarlSeedClusterer clusterer(dist_index);

        SECTION( "Two clusters in a chain and loop of snarl boundary" ) {
            vector<SnarlSeedClusterer::Seed> seeds;
            vector<pos_t> pos_ts;
            pos_ts.emplace_back(2, false, 0);
            pos_ts.emplace_back(3, false, 0);
            pos_ts.emplace_back(5, false, 0);
            pos_ts.emplace_back(16, false, 0);
            //New cluster
            pos_ts.emplace_back(5, false, 10);
            pos_ts.emplace_back(6, false, 0);
            pos_ts.emplace_back(8, false, 0);

            for (pos_t pos : pos_ts){

                std::tuple<bool, size_t, size_t, bool, size_t, size_t, size_t, size_t, bool> chain_info = dist_index.get_minimizer_distances(pos);
                seeds.push_back({ pos, 0, std::get<0>(chain_info), std::get<1>(chain_info), std::get<2>(chain_info),
                   std::get<3>(chain_info), std::get<4>(chain_info), std::get<5>(chain_info), std::get<6>(chain_info), std::get<7>(chain_info), std::get<8>(chain_info)});
            }
            vector<SnarlSeedClusterer::Cluster> clusters =  clusterer.cluster_seeds(seeds, 3); 

            REQUIRE( clusters.size() == 2);

            vector<hash_set<size_t>> cluster_sets;
            for (auto& c : clusters) {
                hash_set<size_t> h;
                for (size_t s : c.seeds) {
                    h.insert(s);
                }
                cluster_sets.push_back(h);
            }
            REQUIRE (( (cluster_sets[0].count(0) == 1 &&
                       cluster_sets[0].count(1) == 1 &&
                       cluster_sets[0].count(2) == 1 &&
                       cluster_sets[0].count(3) == 1 &&
                       cluster_sets[1].count(4) == 1 &&
                       cluster_sets[1].count(5) == 1 &&
                       cluster_sets[1].count(6) == 1)  ||

                     ( cluster_sets[1].count(0) == 1 &&
                       cluster_sets[1].count(1) == 1 &&
                       cluster_sets[1].count(2) == 1 &&
                       cluster_sets[0].count(3) == 1 &&
                       cluster_sets[0].count(4) == 1 &&
                       cluster_sets[0].count(5) == 1 &&
                       cluster_sets[0].count(6) == 1  )));
        }
        SECTION( "Four clusters" ) {
            vector<SnarlSeedClusterer::Seed> seeds;
            vector<pos_t> pos_ts;
            pos_ts.emplace_back(3, false, 0);
            pos_ts.emplace_back(5, false, 0);
            pos_ts.emplace_back(16, false, 0);
            //New cluster
            pos_ts.emplace_back(5, false, 8);
            //new_cluster
            pos_ts.emplace_back(6, false, 0);
            pos_ts.emplace_back(8, false, 0);
            //New_cluster
            pos_ts.emplace_back(13, false, 1);
            pos_ts.emplace_back(14, false, 0);
            pos_ts.emplace_back(15, false, 0);

            for (pos_t pos : pos_ts){
                std::tuple<bool, size_t, size_t, bool, size_t, size_t, size_t, size_t, bool> chain_info = dist_index.get_minimizer_distances(pos);
                seeds.push_back({ pos, 0, std::get<0>(chain_info), std::get<1>(chain_info), std::get<2>(chain_info),
                   std::get<3>(chain_info), std::get<4>(chain_info), std::get<5>(chain_info), std::get<6>(chain_info), std::get<7>(chain_info), std::get<8>(chain_info)});
            }
            vector<SnarlSeedClusterer::Cluster> clusters =  clusterer.cluster_seeds(seeds, 3); 

            REQUIRE( clusters.size() == 4);

            vector<vector<SnarlSeedClusterer::Seed>> all_seeds;

            all_seeds.push_back(seeds);
            vector<vector<SnarlSeedClusterer::Cluster>> paired_clusters = clusterer.cluster_seeds(all_seeds, 3, 3); 

            REQUIRE( paired_clusters.size() == 1);
            REQUIRE( paired_clusters[0].size() == 4);
            REQUIRE( paired_clusters[0][0].fragment != paired_clusters[0][1].fragment);
            REQUIRE( paired_clusters[0][0].fragment != paired_clusters[0][2].fragment);
            REQUIRE( paired_clusters[0][0].fragment != paired_clusters[0][3].fragment);
            REQUIRE( paired_clusters[0][1].fragment != paired_clusters[0][2].fragment);
            REQUIRE( paired_clusters[0][1].fragment != paired_clusters[0][3].fragment);
            REQUIRE( paired_clusters[0][2].fragment != paired_clusters[0][3].fragment);

            //New fragment clusters
        } SECTION ("Four fragment clusters") {
            vector<vector<SnarlSeedClusterer::Seed>> all_seeds;
            vector<SnarlSeedClusterer::Seed> seeds;
            vector<pos_t>pos_ts;
            pos_ts.emplace_back(3, false, 0);
            pos_ts.emplace_back(5, false, 0);
            pos_ts.emplace_back(16, false, 0);
            //New cluster
            pos_ts.emplace_back(6, false, 0);
            pos_ts.emplace_back(8, false, 0);
            for (pos_t pos : pos_ts){
                std::tuple<bool, size_t, size_t, bool, size_t, size_t, size_t, size_t, bool> chain_info = dist_index.get_minimizer_distances(pos);
                seeds.push_back({ pos, 0, std::get<0>(chain_info), std::get<1>(chain_info), std::get<2>(chain_info),
                   std::get<3>(chain_info), std::get<4>(chain_info), std::get<5>(chain_info), std::get<6>(chain_info), std::get<7>(chain_info), std::get<8>(chain_info)});
            }
            all_seeds.push_back(seeds);
            seeds.clear();
            pos_ts.clear();
            //New cluster
            pos_ts.emplace_back(5, false, 8);
            //New cluster
            pos_ts.emplace_back(13, false, 1);
            pos_ts.emplace_back(14, false, 0);
            pos_ts.emplace_back(15, false, 0);
            for (pos_t pos : pos_ts){
                std::tuple<bool, size_t, size_t, bool, size_t, size_t, size_t, size_t, bool> chain_info = dist_index.get_minimizer_distances(pos);
                seeds.push_back({ pos, 0, std::get<0>(chain_info), std::get<1>(chain_info), std::get<2>(chain_info),
                   std::get<3>(chain_info), std::get<4>(chain_info), std::get<5>(chain_info), std::get<6>(chain_info), std::get<7>(chain_info), std::get<8>(chain_info)});
            }
            all_seeds.push_back(seeds);

            vector<vector<SnarlSeedClusterer::Cluster>> paired_clusters = clusterer.cluster_seeds(all_seeds, 3, 3);

            REQUIRE( paired_clusters.size() == 2);
            REQUIRE( paired_clusters[0].size() == 2);
            REQUIRE( paired_clusters[1].size() == 2);
            REQUIRE( paired_clusters[0][0].fragment != paired_clusters[0][1].fragment);
            REQUIRE( paired_clusters[0][0].fragment != paired_clusters[1][0].fragment);
            REQUIRE( paired_clusters[0][0].fragment != paired_clusters[1][1].fragment);
            REQUIRE( paired_clusters[0][1].fragment != paired_clusters[1][0].fragment);
            REQUIRE( paired_clusters[0][1].fragment != paired_clusters[1][1].fragment);

            //New fragment clusters

            paired_clusters = clusterer.cluster_seeds(all_seeds, 3, 5);

            REQUIRE( paired_clusters.size() == 2);
            REQUIRE( paired_clusters[0].size() == 2);
            REQUIRE( paired_clusters[1].size() == 2);
            REQUIRE( paired_clusters[0][0].fragment == paired_clusters[0][1].fragment);
            REQUIRE( paired_clusters[0][0].fragment == paired_clusters[1][0].fragment);
            REQUIRE( paired_clusters[0][0].fragment != paired_clusters[1][1].fragment);
            REQUIRE( paired_clusters[0][1].fragment == paired_clusters[1][0].fragment);
            REQUIRE( paired_clusters[0][1].fragment != paired_clusters[1][1].fragment);
            REQUIRE( paired_clusters[1][0].fragment != paired_clusters[1][1].fragment);
        }
        SECTION( "Same node, same cluster" ) {
            vector<SnarlSeedClusterer::Seed> seeds;
            vector<pos_t> pos_ts;
            pos_ts.emplace_back(5, false, 0);
            pos_ts.emplace_back(5, false, 11);
            pos_ts.emplace_back(5, false, 5);

            for (pos_t pos : pos_ts){
                std::tuple<bool, size_t, size_t, bool, size_t, size_t, size_t, size_t, bool> chain_info = dist_index.get_minimizer_distances(pos);
                seeds.push_back({ pos, 0, std::get<0>(chain_info), std::get<1>(chain_info), std::get<2>(chain_info),
                   std::get<3>(chain_info), std::get<4>(chain_info), std::get<5>(chain_info), std::get<6>(chain_info), std::get<7>(chain_info), std::get<8>(chain_info)});
            }
            vector<SnarlSeedClusterer::Cluster> clusters = clusterer.cluster_seeds(seeds, 7); 


            REQUIRE( clusters.size() == 1);
        }
    }//end test case
    TEST_CASE("Top level unary snarl", "[cluster]") {
        VG graph;

        Node* n1 = graph.create_node("GCA");
        Node* n2 = graph.create_node("T");
        Node* n3 = graph.create_node("G");
        Node* n4 = graph.create_node("CTGA");
        Node* n5 = graph.create_node("GCA");
        Node* n6 = graph.create_node("T");
        Node* n7 = graph.create_node("G");

        Edge* e1 = graph.create_edge(n1, n2);
        Edge* e2 = graph.create_edge(n1, n3);
        Edge* e3 = graph.create_edge(n2, n7);
        Edge* e4 = graph.create_edge(n3, n4);
        Edge* e5 = graph.create_edge(n3, n5);
        Edge* e6 = graph.create_edge(n4, n6);
        Edge* e7 = graph.create_edge(n5, n6);
        Edge* e8 = graph.create_edge(n6, n7);
        Edge* e9 = graph.create_edge(n1, n1, true, false);

        CactusSnarlFinder bubble_finder(graph);
        SnarlManager snarl_manager = bubble_finder.find_snarls();
        MinimumDistanceIndex dist_index (&graph, &snarl_manager);
        SnarlSeedClusterer clusterer(dist_index);



        // We end up with a big unary snarl of 7 rev -> 7 rev
        // Inside that we have a chain of two normal snarls 2 rev -> 3 fwd,      and 3 fwd -> 6 fwd
        // And inside 2 rev -> 3 fwd, we get 1 rev -> 1 rev as another unar     y snarl.

        // We name the snarls for the distance index by their start nodes.

        SECTION("Top level cluster") {
            vector<id_t> ids({1, 2, 7});
            vector<SnarlSeedClusterer::Seed> seeds;
            for (id_t n : ids) {
                pos_t pos = make_pos_t(n, false, 0);
                std::tuple<bool, size_t, size_t, bool, size_t, size_t, size_t, size_t, bool> chain_info = dist_index.get_minimizer_distances(pos);
                seeds.push_back({ pos, 0, std::get<0>(chain_info), std::get<1>(chain_info), std::get<2>(chain_info),
                   std::get<3>(chain_info), std::get<4>(chain_info), std::get<5>(chain_info), std::get<6>(chain_info), std::get<7>(chain_info), std::get<8>(chain_info)});
            }

            vector<SnarlSeedClusterer::Cluster> clusters= clusterer.cluster_seeds(seeds, 10); 


            REQUIRE( clusters.size() == 1);
        }
        SECTION("One cluster") {
            vector<SnarlSeedClusterer::Seed> seeds;
            vector<pos_t> pos_ts;
            pos_ts.emplace_back(1, false, 0);
            pos_ts.emplace_back(2, false, 0);
            pos_ts.emplace_back(7, false, 0);
            pos_ts.emplace_back(4, false, 0);

            for (pos_t pos : pos_ts){
                std::tuple<bool, size_t, size_t, bool, size_t, size_t, size_t, size_t, bool> chain_info = dist_index.get_minimizer_distances(pos);
                seeds.push_back({ pos, 0, std::get<0>(chain_info), std::get<1>(chain_info), std::get<2>(chain_info),
                   std::get<3>(chain_info), std::get<4>(chain_info), std::get<5>(chain_info), std::get<6>(chain_info), std::get<7>(chain_info), std::get<8>(chain_info)});
            }
            vector<SnarlSeedClusterer::Cluster> clusters = clusterer.cluster_seeds(seeds, 10); 


            REQUIRE( clusters.size() == 1);
        }
        SECTION("One cluster") {
            vector<SnarlSeedClusterer::Seed> seeds;
            vector<pos_t> pos_ts;
            pos_ts.emplace_back(2, false, 0);
            pos_ts.emplace_back(4, false, 0);

            for (pos_t pos : pos_ts){
                std::tuple<bool, size_t, size_t, bool, size_t, size_t, size_t, size_t, bool> chain_info = dist_index.get_minimizer_distances(pos);
                seeds.push_back({ pos, 0, std::get<0>(chain_info), std::get<1>(chain_info), std::get<2>(chain_info),
                   std::get<3>(chain_info), std::get<4>(chain_info), std::get<5>(chain_info), std::get<6>(chain_info), std::get<7>(chain_info), std::get<8>(chain_info)});
            }
            vector<SnarlSeedClusterer::Cluster> clusters =  clusterer.cluster_seeds(seeds, 10); 



            REQUIRE( clusters.size() == 1);
        }
        SECTION("Two clusters") {
            vector<SnarlSeedClusterer::Seed> seeds;
            vector<pos_t> pos_ts;
            pos_ts.emplace_back(2, false, 0);
            pos_ts.emplace_back(4, false, 1);
            pos_ts.emplace_back(6, false, 0);

            for (pos_t pos : pos_ts){
                std::tuple<bool, size_t, size_t, bool, size_t, size_t, size_t, size_t, bool> chain_info = dist_index.get_minimizer_distances(pos);
                seeds.push_back({ pos, 0, std::get<0>(chain_info), std::get<1>(chain_info), std::get<2>(chain_info),
                   std::get<3>(chain_info), std::get<4>(chain_info), std::get<5>(chain_info), std::get<6>(chain_info), std::get<7>(chain_info), std::get<8>(chain_info)});
            }
            vector<SnarlSeedClusterer::Cluster> clusters = clusterer.cluster_seeds(seeds, 5); 


            REQUIRE( clusters.size() == 2);
        }
        SECTION("No clusters") {
            vector<SnarlSeedClusterer::Seed> seeds;

            vector<SnarlSeedClusterer::Cluster> clusters =  clusterer.cluster_seeds(seeds, 5); 


            REQUIRE( clusters.size() == 0);
        }
    }
    TEST_CASE( "Long chain",
                   "[cluster]" ) {
        VG graph;

        Node* n1 = graph.create_node("GCA");
        Node* n2 = graph.create_node("T");
        Node* n3 = graph.create_node("G");
        Node* n4 = graph.create_node("CTGA");
        Node* n5 = graph.create_node("GCA");
        Node* n6 = graph.create_node("T");
        Node* n7 = graph.create_node("G");
        Node* n8 = graph.create_node("CTGA");
        Node* n9 = graph.create_node("TTA");
        Node* n10 = graph.create_node("G");
        Node* n11 = graph.create_node("CTGA");
        Node* n12 = graph.create_node("G");
        Node* n13 = graph.create_node("CTGA");
        Node* n14 = graph.create_node("CTGA");

        Edge* e1 = graph.create_edge(n1, n2);
        Edge* e2 = graph.create_edge(n1, n6);
        Edge* e3 = graph.create_edge(n2, n3);
        Edge* e4 = graph.create_edge(n2, n4);
        Edge* e5 = graph.create_edge(n3, n5);
        Edge* e6 = graph.create_edge(n4, n5);
        Edge* e7 = graph.create_edge(n5, n6);
        Edge* e8 = graph.create_edge(n6, n7);
        Edge* e9 = graph.create_edge(n6, n8);
        Edge* e10 = graph.create_edge(n7, n8);
        Edge* e11 = graph.create_edge(n8, n9);
        Edge* e12 = graph.create_edge(n8, n12);
        Edge* e13 = graph.create_edge(n9, n10);
        Edge* e14 = graph.create_edge(n9, n11);
        Edge* e15 = graph.create_edge(n10, n11);
        Edge* e16 = graph.create_edge(n11, n12);
        Edge* e17 = graph.create_edge(n12, n13);
        Edge* e18 = graph.create_edge(n12, n14);
        Edge* e19 = graph.create_edge(n13, n14);

        CactusSnarlFinder bubble_finder(graph);
        SnarlManager snarl_manager = bubble_finder.find_snarls();
        MinimumDistanceIndex dist_index (&graph, &snarl_manager);
        SnarlSeedClusterer clusterer(dist_index);

        SECTION("Snarl then seed") {

            vector<id_t> ids({3, 5, 6, 11});
            vector<SnarlSeedClusterer::Seed> seeds;
            for (id_t n : ids) {
                pos_t pos = make_pos_t(n, false, 0);
                std::tuple<bool, size_t, size_t, bool, size_t, size_t, size_t, size_t, bool> chain_info = dist_index.get_minimizer_distances(pos);
                seeds.push_back({ pos, 0, std::get<0>(chain_info), std::get<1>(chain_info), std::get<2>(chain_info),
                   std::get<3>(chain_info), std::get<4>(chain_info), std::get<5>(chain_info), std::get<6>(chain_info), std::get<7>(chain_info), std::get<8>(chain_info)});
            }

            vector<SnarlSeedClusterer::Cluster> clusters =  clusterer.cluster_seeds(seeds, 5); 


            REQUIRE( clusters.size() == 2);

        }
        SECTION("Seed then snarl") {

            vector<id_t> ids({1, 2, 3, 5, 6, 11, 10});
            vector<SnarlSeedClusterer::Seed> seeds;
            for (id_t n : ids) {
                pos_t pos = make_pos_t(n, false, 0);
                std::tuple<bool, size_t, size_t, bool, size_t, size_t, size_t, size_t, bool> chain_info = dist_index.get_minimizer_distances(pos);
                seeds.push_back({ pos, 0, std::get<0>(chain_info), std::get<1>(chain_info), std::get<2>(chain_info),
                   std::get<3>(chain_info), std::get<4>(chain_info), std::get<5>(chain_info), std::get<6>(chain_info), std::get<7>(chain_info), std::get<8>(chain_info)});
            }

            vector<SnarlSeedClusterer::Cluster> clusters =  clusterer.cluster_seeds(seeds, 5); 


            REQUIRE( clusters.size() == 2);

        }
        SECTION("Only seeds") {

            vector<id_t> ids({1, 6, 14});
            vector<SnarlSeedClusterer::Seed> seeds;
            for (id_t n : ids) {
                pos_t pos = make_pos_t(n, false, 0);
                std::tuple<bool, size_t, size_t, bool, size_t, size_t, size_t, size_t, bool> chain_info = dist_index.get_minimizer_distances(pos);
                seeds.push_back({ pos, 0, std::get<0>(chain_info), std::get<1>(chain_info), std::get<2>(chain_info),
                   std::get<3>(chain_info), std::get<4>(chain_info), std::get<5>(chain_info), std::get<6>(chain_info), std::get<7>(chain_info), std::get<8>(chain_info)});
            }

            vector<SnarlSeedClusterer::Cluster> clusters =  clusterer.cluster_seeds(seeds, 4); 


            REQUIRE( clusters.size() == 2);

        }
        SECTION("Only seeds two reads") {

            vector<id_t> ids({1, 6, 14});
            vector<SnarlSeedClusterer::Seed> seeds;
            for (id_t n : ids) {
                pos_t pos = make_pos_t(n, false, 0);
                std::tuple<bool, size_t, size_t, bool, size_t, size_t, size_t, size_t, bool> chain_info = dist_index.get_minimizer_distances(pos);
                seeds.push_back({ pos, 0, std::get<0>(chain_info), std::get<1>(chain_info), std::get<2>(chain_info),
                   std::get<3>(chain_info), std::get<4>(chain_info), std::get<5>(chain_info), std::get<6>(chain_info), std::get<7>(chain_info), std::get<8>(chain_info)});
            }
            vector<id_t> ids1({8, 12});
            vector<SnarlSeedClusterer::Seed> seeds1;
            for (id_t n : ids1) {
                pos_t pos = make_pos_t(n, false, 0);
                std::tuple<bool, size_t, size_t, bool, size_t, size_t, size_t, size_t, bool> chain_info = dist_index.get_minimizer_distances(pos);
                seeds1.push_back({ pos, 0, std::get<0>(chain_info), std::get<1>(chain_info), std::get<2>(chain_info),
                   std::get<3>(chain_info), std::get<4>(chain_info), std::get<5>(chain_info), std::get<6>(chain_info), std::get<7>(chain_info), std::get<8>(chain_info)});
            }
            vector<vector<SnarlSeedClusterer::Seed>> all_seeds;
            all_seeds.emplace_back(seeds);
            all_seeds.emplace_back(seeds1);

            vector<vector<SnarlSeedClusterer::Cluster>> clusters =  clusterer.cluster_seeds(all_seeds, 4, 5); 


            REQUIRE( clusters.size() == 2);
            REQUIRE( clusters[0].size() == 2);
            REQUIRE( clusters[1].size() == 1);
            REQUIRE( clusters[0][0].fragment == clusters[0][1].fragment);
            REQUIRE( clusters[0][0].fragment == clusters[1][0].fragment);

        }
        SECTION("Only snarls") {

            vector<id_t> ids({4, 5, 9});
            vector<SnarlSeedClusterer::Seed> seeds;
            for (id_t n : ids) {
                pos_t pos = make_pos_t(n, false, 0);
                std::tuple<bool, size_t, size_t, bool, size_t, size_t, size_t, size_t, bool> chain_info = dist_index.get_minimizer_distances(pos);
                seeds.push_back({ pos, 0, std::get<0>(chain_info), std::get<1>(chain_info), std::get<2>(chain_info),
                   std::get<3>(chain_info), std::get<4>(chain_info), std::get<5>(chain_info), std::get<6>(chain_info), std::get<7>(chain_info), std::get<8>(chain_info)});
            }

            vector<SnarlSeedClusterer::Cluster> clusters =  clusterer.cluster_seeds(seeds, 9); 


            REQUIRE( clusters.size() == 1);

        }
        SECTION("Skip snarl") {

            vector<id_t> ids({7, 10, 13});
            vector<SnarlSeedClusterer::Seed> seeds;
            for (id_t n : ids) {
                pos_t pos = make_pos_t(n, false, 0);
                std::tuple<bool, size_t, size_t, bool, size_t, size_t, size_t, size_t, bool> chain_info = dist_index.get_minimizer_distances(pos);
                seeds.push_back({ pos, 0, std::get<0>(chain_info), std::get<1>(chain_info), std::get<2>(chain_info),
                   std::get<3>(chain_info), std::get<4>(chain_info), std::get<5>(chain_info), std::get<6>(chain_info), std::get<7>(chain_info), std::get<8>(chain_info)});
            }

            vector<SnarlSeedClusterer::Cluster> clusters =  clusterer.cluster_seeds(seeds, 6); 

            REQUIRE( clusters.size() == 1);
        }
    }

    TEST_CASE( "Disconnected graph",
                   "[cluster]" ) {
        VG graph;

        Node* n1 = graph.create_node("GCA");
        Node* n2 = graph.create_node("T");
        Node* n3 = graph.create_node("G");
        Node* n4 = graph.create_node("CTGA");
        Node* n5 = graph.create_node("GCA");
        Node* n6 = graph.create_node("T");
        Node* n7 = graph.create_node("G");
        Node* n8 = graph.create_node("CTGA");
//Disconnected
        Node* n9 = graph.create_node("T");
        Node* n10 = graph.create_node("G");
        Node* n11 = graph.create_node("CTGA");
        Node* n12 = graph.create_node("G");
        Node* n13 = graph.create_node("CTGA");

        Edge* e1 = graph.create_edge(n1, n2);
        Edge* e2 = graph.create_edge(n1, n3);
        Edge* e3 = graph.create_edge(n2, n3);
        Edge* e4 = graph.create_edge(n3, n4);
        Edge* e5 = graph.create_edge(n3, n5);
        Edge* e6 = graph.create_edge(n4, n5);
        Edge* e7 = graph.create_edge(n5, n6);
        Edge* e8 = graph.create_edge(n5, n7);
        Edge* e9 = graph.create_edge(n6, n8);
        Edge* e10 = graph.create_edge(n7, n8);

        Edge* e11 = graph.create_edge(n9, n10);
        Edge* e12 = graph.create_edge(n9, n11);
        Edge* e13 = graph.create_edge(n10, n11);
        Edge* e14 = graph.create_edge(n11, n12);
        Edge* e15 = graph.create_edge(n11, n13);
        Edge* e16 = graph.create_edge(n12, n13);

        CactusSnarlFinder bubble_finder(graph);
        SnarlManager snarl_manager = bubble_finder.find_snarls();
        MinimumDistanceIndex dist_index (&graph, &snarl_manager);
        SnarlSeedClusterer clusterer(dist_index);

        SECTION("Two clusters") {
            vector<SnarlSeedClusterer::Seed> seeds;
            vector<pos_t> pos_ts;
            pos_ts.emplace_back(2, false, 0);
            pos_ts.emplace_back(3, false, 0);
            pos_ts.emplace_back(9, false, 0);

            for (pos_t pos : pos_ts){
                std::tuple<bool, size_t, size_t, bool, size_t, size_t, size_t, size_t, bool> chain_info = dist_index.get_minimizer_distances(pos);
                seeds.push_back({ pos, 0, std::get<0>(chain_info), std::get<1>(chain_info), std::get<2>(chain_info),
                   std::get<3>(chain_info), std::get<4>(chain_info), std::get<5>(chain_info), std::get<6>(chain_info), std::get<7>(chain_info), std::get<8>(chain_info)});
            }
            vector<SnarlSeedClusterer::Cluster> clusters =  clusterer.cluster_seeds(seeds, 5); 


            REQUIRE( clusters.size() == 2);

        }
        SECTION("Two clusters with seed structs") {

            vector<id_t> ids({2, 3, 9});
            vector<SnarlSeedClusterer::Seed> seeds;
            for (id_t n : ids) {
                pos_t pos = make_pos_t(n, false, 0);
                std::tuple<bool, size_t, size_t, bool, size_t, size_t, size_t, size_t, bool> chain_info = dist_index.get_minimizer_distances(pos);
                seeds.push_back({ pos, 0, std::get<0>(chain_info), std::get<1>(chain_info), std::get<2>(chain_info),
                   std::get<3>(chain_info), std::get<4>(chain_info), std::get<5>(chain_info), std::get<6>(chain_info), std::get<7>(chain_info), std::get<8>(chain_info)});
            }

            vector<SnarlSeedClusterer::Cluster> clusters =  clusterer.cluster_seeds(seeds, 5); 


            REQUIRE( clusters.size() == 2);

        }
        SECTION("Two clusters with seed structs") {

            vector<id_t> ids({2, 3, 5, 9, 10});
            vector<SnarlSeedClusterer::Seed> seeds;
            for (id_t n : ids) {
                pos_t pos = make_pos_t(n, false, 0);
                std::tuple<bool, size_t, size_t, bool, size_t, size_t, size_t, size_t, bool> chain_info = dist_index.get_minimizer_distances(pos);
                seeds.push_back({ pos, 0, std::get<0>(chain_info), std::get<1>(chain_info), std::get<2>(chain_info),
                   std::get<3>(chain_info), std::get<4>(chain_info), std::get<5>(chain_info), std::get<6>(chain_info), std::get<7>(chain_info), std::get<8>(chain_info)});
            }

            vector<SnarlSeedClusterer::Cluster> clusters =  clusterer.cluster_seeds(seeds, 5); 


            REQUIRE( clusters.size() == 2);

        }
        SECTION("Two top level clusters") {

            vector<id_t> ids({1, 3, 11});
            vector<SnarlSeedClusterer::Seed> seeds;
            for (id_t n : ids) {
                pos_t pos = make_pos_t(n, false, 0);
                std::tuple<bool, size_t, size_t, bool, size_t, size_t, size_t, size_t, bool> chain_info = dist_index.get_minimizer_distances(pos);
                seeds.push_back({ pos, 0, std::get<0>(chain_info), std::get<1>(chain_info), std::get<2>(chain_info),
                   std::get<3>(chain_info), std::get<4>(chain_info), std::get<5>(chain_info), std::get<6>(chain_info), std::get<7>(chain_info), std::get<8>(chain_info)});
            }
            vector<id_t> ids1({5, 13});
            vector<SnarlSeedClusterer::Seed> seeds1;
            for (id_t n : ids1) {
                pos_t pos = make_pos_t(n, false, 0);
                std::tuple<bool, size_t, size_t, bool, size_t, size_t, size_t, size_t, bool> chain_info = dist_index.get_minimizer_distances(pos);
                seeds1.push_back({ pos, 0, std::get<0>(chain_info), std::get<1>(chain_info), std::get<2>(chain_info),
                   std::get<3>(chain_info), std::get<4>(chain_info), std::get<5>(chain_info), std::get<6>(chain_info), std::get<7>(chain_info), std::get<8>(chain_info)});
            }
            //Clusters are 
            //Read 1: {1, 3} in a fragment cluster with Read 2: {5}
            //Read 1: {11} in a fragment cluster with Read 2: {13}
            vector<vector<SnarlSeedClusterer::Seed>> all_seeds;
            all_seeds.emplace_back(seeds);
            all_seeds.emplace_back(seeds1);


            vector<vector<SnarlSeedClusterer::Cluster>> clusters =  clusterer.cluster_seeds(all_seeds, 5, 10); 


            REQUIRE( clusters.size() == 2);
            REQUIRE( clusters[0].size() == 2);
            REQUIRE( clusters[1].size() == 2);
            REQUIRE( clusters[0][0].fragment != clusters[0][1].fragment);
            REQUIRE( clusters[1][0].fragment != clusters[1][1].fragment);

            REQUIRE(( clusters[0][0].fragment == clusters[1][0].fragment || clusters[0][0].fragment == clusters[1][1].fragment));
            REQUIRE(( clusters[0][1].fragment == clusters[1][0].fragment || clusters[0][1].fragment == clusters[1][1].fragment));


        }
    }
    TEST_CASE("Top level loop creates looping chain", "[cluster]") {
        VG graph;
 
        Node* n1 = graph.create_node("G");
        Node* n2 = graph.create_node("T");
        Node* n3 = graph.create_node("G");
        Node* n4 = graph.create_node("CTGAAAAAAAAAAAA"); //15
        Node* n5 = graph.create_node("GCAA");
        Node* n6 = graph.create_node("T");
        Node* n7 = graph.create_node("G");
        Node* n8 = graph.create_node("A");
        Node* n9 = graph.create_node("T");
        Node* n10 = graph.create_node("G");
        Node* n11 = graph.create_node("GGGGG");
 
        Edge* e1 = graph.create_edge(n9, n1);
        Edge* e2 = graph.create_edge(n9, n11);
        Edge* e3 = graph.create_edge(n1, n2);
        Edge* e4 = graph.create_edge(n1, n8);
        Edge* e5 = graph.create_edge(n2, n3);
        Edge* e6 = graph.create_edge(n2, n4);
        Edge* e7 = graph.create_edge(n3, n5);
        Edge* e8 = graph.create_edge(n4, n5);
        Edge* e9 = graph.create_edge(n5, n6);
        Edge* e10 = graph.create_edge(n5, n7);
        Edge* e11 = graph.create_edge(n6, n7);
        Edge* e12 = graph.create_edge(n7, n8);
        Edge* e13 = graph.create_edge(n8, n10);
        Edge* e16 = graph.create_edge(n10, n9);
        Edge* e17 = graph.create_edge(n2, n2, true, false);
        Edge* e18 = graph.create_edge(n11, n10);

        CactusSnarlFinder bubble_finder(graph);
        SnarlManager snarl_manager = bubble_finder.find_snarls();
        MinimumDistanceIndex dist_index (&graph, &snarl_manager);
        SnarlSeedClusterer clusterer(dist_index);

        SECTION("Two clusters") {
            vector<SnarlSeedClusterer::Seed> seeds;
            vector<pos_t> pos_ts;
            pos_ts.emplace_back(2, false, 0);
            pos_ts.emplace_back(3, false, 0);
            pos_ts.emplace_back(8, false, 0);

            for (pos_t pos : pos_ts){
                std::tuple<bool, size_t, size_t, bool, size_t, size_t, size_t, size_t, bool> chain_info = dist_index.get_minimizer_distances(pos);
                seeds.push_back({ pos, 0, std::get<0>(chain_info), std::get<1>(chain_info), std::get<2>(chain_info),
                   std::get<3>(chain_info), std::get<4>(chain_info), std::get<5>(chain_info), std::get<6>(chain_info), std::get<7>(chain_info), std::get<8>(chain_info)});
            }
            vector<SnarlSeedClusterer::Cluster> clusters =  clusterer.cluster_seeds(seeds, 3); 


            REQUIRE( clusters.size() == 2);

        }
        SECTION("One cluster") {
            vector<SnarlSeedClusterer::Seed> seeds;
            vector<pos_t> pos_ts;
            pos_ts.emplace_back(1, false, 0);
            pos_ts.emplace_back(2, false, 0);
            pos_ts.emplace_back(7, false, 0);

            for (pos_t pos : pos_ts){
                std::tuple<bool, size_t, size_t, bool, size_t, size_t, size_t, size_t, bool> chain_info = dist_index.get_minimizer_distances(pos);
                seeds.push_back({ pos, 0, std::get<0>(chain_info), std::get<1>(chain_info), std::get<2>(chain_info),
                   std::get<3>(chain_info), std::get<4>(chain_info), std::get<5>(chain_info), std::get<6>(chain_info), std::get<7>(chain_info), std::get<8>(chain_info)});
            }
            vector<SnarlSeedClusterer::Cluster> clusters =  clusterer.cluster_seeds(seeds, 6); 


            REQUIRE( clusters.size() == 1);

        }
        SECTION("One cluster taking chain loop") {
            vector<id_t> ids({8, 9, 10});
            vector<SnarlSeedClusterer::Seed> seeds;
            for (id_t n : ids) {
                pos_t pos = make_pos_t(n, false, 0);
                std::tuple<bool, size_t, size_t, bool, size_t, size_t, size_t, size_t, bool> chain_info = dist_index.get_minimizer_distances(pos);
                seeds.push_back({ pos, 0, std::get<0>(chain_info), std::get<1>(chain_info), std::get<2>(chain_info),
                   std::get<3>(chain_info), std::get<4>(chain_info), std::get<5>(chain_info), std::get<6>(chain_info), std::get<7>(chain_info), std::get<8>(chain_info)});
            }

            vector<SnarlSeedClusterer::Cluster> clusters =  clusterer.cluster_seeds(seeds, 3); 


            REQUIRE( clusters.size() == 1);

        }
    }//End test case


    TEST_CASE( "Nested unary snarls","[cluster]" ) {
        VG graph;

        Node* n1 = graph.create_node("GCA");
        Node* n2 = graph.create_node("T");
        Node* n3 = graph.create_node("G");
        Node* n4 = graph.create_node("CTGA");
        Node* n5 = graph.create_node("G");
        Node* n6 = graph.create_node("T");
        Node* n7 = graph.create_node("G");
        Node* n8 = graph.create_node("G");

        Edge* e1 = graph.create_edge(n1, n2);
        Edge* e2 = graph.create_edge(n1, n3);
        Edge* e3 = graph.create_edge(n2, n4);
        Edge* e4 = graph.create_edge(n3, n4);
        Edge* e5 = graph.create_edge(n4, n5);
        Edge* e6 = graph.create_edge(n4, n6);
        Edge* e7 = graph.create_edge(n5, n6);
        Edge* e8 = graph.create_edge(n6, n7);
        Edge* e9 = graph.create_edge(n6, n8);
        Edge* e10 = graph.create_edge(n7, n8);
        Edge* e11 = graph.create_edge(n8, n8, false, true);

        CactusSnarlFinder bubble_finder(graph);
        SnarlManager snarl_manager = bubble_finder.find_snarls();
        MinimumDistanceIndex dist_index (&graph, &snarl_manager);

        SnarlSeedClusterer clusterer(dist_index);
        //Unary snarl at 8 nested in unary snarl at 6 nested in 
        //unary snarl at  4 nested in regular snarl at 2 (ending at 3)
        //nested in unary snarl at 1

        SECTION( "One cluster" ) {
            vector<id_t> ids({4, 3});
            vector<SnarlSeedClusterer::Seed> seeds;
            for (id_t n : ids) {
                pos_t pos = make_pos_t(n, false, 0);
                std::tuple<bool, size_t, size_t, bool, size_t, size_t, size_t, size_t, bool> chain_info = dist_index.get_minimizer_distances(pos);
                seeds.push_back({ pos, 0, std::get<0>(chain_info), std::get<1>(chain_info), std::get<2>(chain_info),
                   std::get<3>(chain_info), std::get<4>(chain_info), std::get<5>(chain_info), std::get<6>(chain_info), std::get<7>(chain_info), std::get<8>(chain_info)});
            }

            vector<SnarlSeedClusterer::Cluster> clusters =  clusterer.cluster_seeds(seeds, 10); 


            REQUIRE( clusters.size() == 1);
        }
        SECTION( "One cluster nested" ) {
            vector<id_t> ids({5, 3});
            vector<SnarlSeedClusterer::Seed> seeds;
            for (id_t n : ids) {
                pos_t pos = make_pos_t(n, false, 0);
                std::tuple<bool, size_t, size_t, bool, size_t, size_t, size_t, size_t, bool> chain_info = dist_index.get_minimizer_distances(pos);
                seeds.push_back({ pos, 0, std::get<0>(chain_info), std::get<1>(chain_info), std::get<2>(chain_info), std::get<3>(chain_info), std::get<4>(chain_info), std::get<5>(chain_info), std::get<6>(chain_info), std::get<7>(chain_info), std::get<8>(chain_info)});
            }

            vector<SnarlSeedClusterer::Cluster> clusters =  clusterer.cluster_seeds(seeds, 10); 


            REQUIRE( clusters.size() == 1);
        }
        SECTION( "Three clusters" ) {
            vector<id_t> ids({2, 3, 8});
            vector<SnarlSeedClusterer::Seed> seeds;
            for (id_t n : ids) {
                pos_t pos = make_pos_t(n, false, 0);
                std::tuple<bool, size_t, size_t, bool, size_t, size_t, size_t, size_t, bool> chain_info = dist_index.get_minimizer_distances(pos);
                seeds.push_back({ pos, 0, std::get<0>(chain_info), std::get<1>(chain_info), std::get<2>(chain_info),
                   std::get<3>(chain_info), std::get<4>(chain_info), std::get<5>(chain_info), std::get<6>(chain_info), std::get<7>(chain_info), std::get<8>(chain_info)});
            }

            vector<SnarlSeedClusterer::Cluster> clusters =  clusterer.cluster_seeds(seeds, 3); 



            REQUIRE( clusters.size() == 3);
        }
        SECTION( "One cluster taking loop" ) {
            vector<id_t> ids({2, 3});
            vector<SnarlSeedClusterer::Seed> seeds;
            for (id_t n : ids) {
                pos_t pos = make_pos_t(n, false, 0);
                std::tuple<bool, size_t, size_t, bool, size_t, size_t, size_t, size_t, bool> chain_info = dist_index.get_minimizer_distances(pos);
                seeds.push_back({ pos, 0, std::get<0>(chain_info), std::get<1>(chain_info), std::get<2>(chain_info),
                   std::get<3>(chain_info), std::get<4>(chain_info), std::get<5>(chain_info), std::get<6>(chain_info), std::get<7>(chain_info), std::get<8>(chain_info)});
            }

            vector<SnarlSeedClusterer::Cluster> clusters =  clusterer.cluster_seeds(seeds, 15); 


            REQUIRE( clusters.size() == 1);
        }
    }//end test case
    TEST_CASE( "Top level snarl","[cluster]" ) {
        VG graph;

        Node* n1 = graph.create_node("GCA");
        Node* n2 = graph.create_node("T");
        Node* n3 = graph.create_node("G");
        Node* n4 = graph.create_node("CTGA");
        Node* n5 = graph.create_node("GCA");

        Edge* e1 = graph.create_edge(n1, n2);
        Edge* e2 = graph.create_edge(n1, n3);
        Edge* e3 = graph.create_edge(n2, n5);
        Edge* e4 = graph.create_edge(n3, n4);
        Edge* e5 = graph.create_edge(n3, n5);
        Edge* e6 = graph.create_edge(n4, n5);

        CactusSnarlFinder bubble_finder(graph);
        SnarlManager snarl_manager = bubble_finder.find_snarls();
        MinimumDistanceIndex dist_index (&graph, &snarl_manager);

        SnarlSeedClusterer clusterer(dist_index);

        SECTION( "Top level seeds" ) {
            vector<id_t> ids({1, 2, 4});
            vector<SnarlSeedClusterer::Seed> seeds;
            for (id_t n : ids) {
                pos_t pos = make_pos_t(n, false, 0);
                std::tuple<bool, size_t, size_t, bool, size_t, size_t, size_t, size_t, bool> chain_info = dist_index.get_minimizer_distances(pos);
                seeds.push_back({ pos, 0, std::get<0>(chain_info), std::get<1>(chain_info), std::get<2>(chain_info),
                   std::get<3>(chain_info), std::get<4>(chain_info), std::get<5>(chain_info), std::get<6>(chain_info), std::get<7>(chain_info), std::get<8>(chain_info)});
            }

            vector<SnarlSeedClusterer::Cluster> clusters =  clusterer.cluster_seeds(seeds, 3); 


            REQUIRE( clusters.size() == 2);
        }
    }
    TEST_CASE( "Two tip right","[cluster]" ) {
        VG graph;

        Node* n1 = graph.create_node("GCA");
        Node* n2 = graph.create_node("T");
        Node* n3 = graph.create_node("GACCT");
        Node* n4 = graph.create_node("CTGA");
        Node* n5 = graph.create_node("G");
        Node* n6 = graph.create_node("CTGA");
        Node* n7 = graph.create_node("G");

        Edge* e1 = graph.create_edge(n1, n2);
        Edge* e2 = graph.create_edge(n1, n3);
        Edge* e3 = graph.create_edge(n2, n3);
        Edge* e4 = graph.create_edge(n3, n4);
        Edge* e5 = graph.create_edge(n3, n5);
        Edge* e6 = graph.create_edge(n6, n1);
        Edge* e7 = graph.create_edge(n6, n7);
        Edge* e8 = graph.create_edge(n7, n1);
        Edge* e9 = graph.create_edge(n1, n1, true, false);

        CactusSnarlFinder bubble_finder(graph);
        SnarlManager snarl_manager = bubble_finder.find_snarls();
        MinimumDistanceIndex dist_index (&graph, &snarl_manager);

        SnarlSeedClusterer clusterer(dist_index);

        SECTION( "Two cluster" ) {
            vector<id_t> ids({4, 5});
            vector<SnarlSeedClusterer::Seed> seeds;
            for (id_t n : ids) {
                pos_t pos = make_pos_t(n, false, 0);
                std::tuple<bool, size_t, size_t, bool, size_t, size_t, size_t, size_t, bool> chain_info = dist_index.get_minimizer_distances(pos);
                seeds.push_back({ pos, 0, std::get<0>(chain_info), std::get<1>(chain_info), std::get<2>(chain_info),
                   std::get<3>(chain_info), std::get<4>(chain_info), std::get<5>(chain_info), std::get<6>(chain_info), std::get<7>(chain_info), std::get<8>(chain_info)});
            }

            vector<SnarlSeedClusterer::Cluster> clusters =  clusterer.cluster_seeds(seeds, 5); 


            REQUIRE( clusters.size() == 2);
        }

        SECTION( "One clusters" ) {
            vector<id_t> ids({4, 5, 3});
            vector<SnarlSeedClusterer::Seed> seeds;
            for (id_t n : ids) {
                pos_t pos = make_pos_t(n, false, 0);
                std::tuple<bool, size_t, size_t, bool, size_t, size_t, size_t, size_t, bool> chain_info = dist_index.get_minimizer_distances(pos);
                seeds.push_back({ pos, 0, std::get<0>(chain_info), std::get<1>(chain_info), std::get<2>(chain_info),
                   std::get<3>(chain_info), std::get<4>(chain_info), std::get<5>(chain_info), std::get<6>(chain_info), std::get<7>(chain_info), std::get<8>(chain_info)});
            }

            vector<SnarlSeedClusterer::Cluster> clusters =  clusterer.cluster_seeds(seeds, 10); 


            REQUIRE( clusters.size() == 1);
        }

        SECTION( "One cluster loop" ) {
            vector<id_t> ids({4, 5});
            vector<SnarlSeedClusterer::Seed> seeds;
            for (id_t n : ids) {
                pos_t pos = make_pos_t(n, false, 0);
                std::tuple<bool, size_t, size_t, bool, size_t, size_t, size_t, size_t, bool> chain_info = dist_index.get_minimizer_distances(pos);
                seeds.push_back({ pos, 0, std::get<0>(chain_info), std::get<1>(chain_info), std::get<2>(chain_info),
                   std::get<3>(chain_info), std::get<4>(chain_info), std::get<5>(chain_info), std::get<6>(chain_info), std::get<7>(chain_info), std::get<8>(chain_info)});
            }

            vector<SnarlSeedClusterer::Cluster> clusters =  clusterer.cluster_seeds(seeds, 18); 


            REQUIRE( clusters.size() == 1);
        }
    }
    TEST_CASE( "Two tip left","[cluster]" ) {
        VG graph;

        Node* n1 = graph.create_node("GCA");
        Node* n2 = graph.create_node("T");
        Node* n3 = graph.create_node("G");
        Node* n4 = graph.create_node("CTGA");
        Node* n5 = graph.create_node("G");
        Node* n6 = graph.create_node("G");
        Node* n7 = graph.create_node("CTGA");

        Edge* e1 = graph.create_edge(n1, n3);
        Edge* e2 = graph.create_edge(n2, n3);
        Edge* e3 = graph.create_edge(n3, n4);
        Edge* e4 = graph.create_edge(n3, n5);
        Edge* e5 = graph.create_edge(n4, n5);
        Edge* e6 = graph.create_edge(n5, n6);
        Edge* e7 = graph.create_edge(n5, n7);
        Edge* e8 = graph.create_edge(n6, n7);
        Edge* e9 = graph.create_edge(n5, n5, false, true);

        CactusSnarlFinder bubble_finder(graph);
        SnarlManager snarl_manager = bubble_finder.find_snarls();
        MinimumDistanceIndex dist_index (&graph, &snarl_manager);

        SnarlSeedClusterer clusterer(dist_index);

        SECTION( "One cluster" ) {
            vector<id_t> ids({1, 2});
            vector<SnarlSeedClusterer::Seed> seeds;
            for (id_t n : ids) {
                pos_t pos = make_pos_t(n, false, 0);
                std::tuple<bool, size_t, size_t, bool, size_t, size_t, size_t, size_t, bool> chain_info = dist_index.get_minimizer_distances(pos);
                seeds.push_back({  pos, 0, std::get<0>(chain_info), std::get<1>(chain_info), std::get<2>(chain_info),
                   std::get<3>(chain_info), std::get<4>(chain_info), std::get<5>(chain_info), std::get<6>(chain_info), std::get<7>(chain_info), std::get<8>(chain_info)});
            }

            vector<SnarlSeedClusterer::Cluster> clusters =  clusterer.cluster_seeds(seeds, 5); 


            REQUIRE( clusters.size() == 2);
        }

        SECTION( "Two clusters" ) {
            vector<id_t> ids({1, 2, 3});
            vector<SnarlSeedClusterer::Seed> seeds;
            for (id_t n : ids) {
                pos_t pos = make_pos_t(n, false, 0);
                std::tuple<bool, size_t, size_t, bool, size_t, size_t, size_t, size_t, bool> chain_info = dist_index.get_minimizer_distances(pos);
                seeds.push_back({ pos, 0, std::get<0>(chain_info), std::get<1>(chain_info), std::get<2>(chain_info),
                   std::get<3>(chain_info), std::get<4>(chain_info), std::get<5>(chain_info), std::get<6>(chain_info), std::get<7>(chain_info), std::get<8>(chain_info)});
            }

            vector<SnarlSeedClusterer::Cluster> clusters =  clusterer.cluster_seeds(seeds, 5); 


            REQUIRE( clusters.size() == 1);
        }
        SECTION( "Two clusters with snarl" ) {
            vector<id_t> ids({1, 2, 4});
            vector<SnarlSeedClusterer::Seed> seeds;
            for (id_t n : ids) {
                pos_t pos = make_pos_t(n, false, 0);
                std::tuple<bool, size_t, size_t, bool, size_t, size_t, size_t, size_t, bool> chain_info = dist_index.get_minimizer_distances(pos);
                seeds.push_back({ pos, 0, std::get<0>(chain_info), std::get<1>(chain_info), std::get<2>(chain_info),
                   std::get<3>(chain_info), std::get<4>(chain_info), std::get<5>(chain_info), std::get<6>(chain_info), std::get<7>(chain_info), std::get<8>(chain_info)});
            }

            vector<SnarlSeedClusterer::Cluster> clusters =  clusterer.cluster_seeds(seeds, 5); 


            REQUIRE( clusters.size() == 1);
        }
        SECTION( "One cluster with loop" ) {
            vector<id_t> ids({1, 2});
            vector<SnarlSeedClusterer::Seed> seeds;
            for (id_t n : ids) {
                pos_t pos = make_pos_t(n, false, 0);
                std::tuple<bool, size_t, size_t, bool, size_t, size_t, size_t, size_t, bool> chain_info = dist_index.get_minimizer_distances(pos);
                seeds.push_back({ pos, 0, std::get<0>(chain_info), std::get<1>(chain_info), std::get<2>(chain_info),
                   std::get<3>(chain_info), std::get<4>(chain_info), std::get<5>(chain_info), std::get<6>(chain_info), std::get<7>(chain_info), std::get<8>(chain_info)});
            }

            vector<SnarlSeedClusterer::Cluster> clusters =  clusterer.cluster_seeds(seeds, 7); 


            REQUIRE( clusters.size() == 1);
        }
    }



/*
    TEST_CASE("Load graph", "[cluster]"){

        ifstream vg_stream("testGraph");
        VG vg(vg_stream);
        vg_stream.close();
        CactusSnarlFinder bubble_finder(vg);
        SnarlManager snarl_manager = bubble_finder.find_snarls();

        MinimumDistanceIndex dist_index (&vg, &snarl_manager);
        SnarlSeedClusterer clusterer(dist_index);

        int64_t read_lim = 20;// Distance between read clusters
        int64_t fragment_lim = 30;// Distance between fragment clusters



        vector<SnarlSeedClusterer::Seed> seeds;
        vector<pos_t> pos_ts;
        pos_ts.emplace_back(31, false, 75);
        pos_ts.emplace_back(30, false, 2);

        for (pos_t pos : pos_ts) {
            std::tuple<bool, size_t, size_t, bool, size_t, size_t, size_t, size_t, bool> chain_info = dist_index.get_minimizer_distances(pos);
            seeds.push_back({ pos, 0, std::get<0>(chain_info), std::get<1>(chain_info), std::get<2>(chain_info),
               std::get<3>(chain_info), std::get<4>(chain_info), std::get<5>(chain_info), std::get<6>(chain_info), std::get<7>(chain_info), std::get<8>(chain_info)});
        }
        vector<SnarlSeedClusterer::Cluster> clusters =  clusterer.cluster_seeds(seeds, 30); 
        REQUIRE(clusters.size() == 1);
    }//end test case
    */

    TEST_CASE("Random graphs", "[cluster]"){

/*
        ifstream vg_stream("testGraph");
        VG vg(vg_stream);
        vg_stream.close();
        CactusSnarlFinder bubble_finder(vg);
        SnarlManager snarl_manager = bubble_finder.find_snarls();

        MinimumDistanceIndex dist_index (&vg, &snarl_manager);
        dist_index.print_self();
        SnarlSeedClusterer clusterer(dist_index);



        vector<SnarlSeedClusterer::Seed> seeds;
        vector<pos_t> pos_ts;
        pos_ts.emplace_back(70, false, 5);
        pos_ts.emplace_back(9, false, 0);
        


        for (pos_t pos : pos_ts) {
            std::tuple<bool, size_t, size_t, bool, size_t, size_t, size_t, size_t, bool> chain_info = dist_index.get_minimizer_distances(pos);
            seeds.push_back({ pos, 0, std::get<0>(chain_info), std::get<1>(chain_info), std::get<2>(chain_info),
               std::get<3>(chain_info), std::get<4>(chain_info), std::get<5>(chain_info), std::get<6>(chain_info), std::get<7>(chain_info), std::get<8>(chain_info)});
        }
        vector<SnarlSeedClusterer::Cluster> clusters =  clusterer.cluster_seeds(seeds, 30); 

        assert(clusters.size() == 1);
        REQUIRE(false);
        */




        for (int i = 0; i < 0; i++) {
            // For each random graph
            
            default_random_engine generator(time(NULL));

            uniform_int_distribution<int> variant_count(5, 200);

            VG graph;
            random_graph(1000, 30, variant_count(generator), &graph);


            CactusSnarlFinder bubble_finder(graph);

            SnarlManager snarl_manager = bubble_finder.find_snarls();
            MinimumDistanceIndex dist_index (&graph, &snarl_manager);

            SnarlSeedClusterer clusterer(dist_index);


            vector<id_t> all_nodes;
            graph.for_each_handle([&](const handle_t& h)->bool{
                id_t id = graph.get_id(h);
                all_nodes.push_back(id);
                return true;
            });


            uniform_int_distribution<int> randPosIndex(0, all_nodes.size()-1);

            for (size_t k = 0; k < 20 ; k++) {

                vector<vector<SnarlSeedClusterer::Seed>> all_seeds;
                all_seeds.emplace_back();
                all_seeds.emplace_back();
                int64_t read_lim = 30;// Distance between read clusters
                int64_t fragment_lim = 50;// Distance between fragment clusters
                for (size_t read = 0 ; read < 2 ; read ++) {
                    uniform_int_distribution<int> randPosCount(1, 50);
                    for (int j = 0; j < randPosCount(generator); j++) {
                        //Check clusters of j random positions 
 
                        id_t nodeID1 = all_nodes[randPosIndex(generator)];
                        handle_t node1 = graph.get_handle(nodeID1);
 
                        off_t offset1 = uniform_int_distribution<int>(0,graph.get_length(node1) - 1)(generator);

                        pos_t pos = make_pos_t(nodeID1,
                            uniform_int_distribution<int>(0,1)(generator) == 0,offset1 );
                        std::tuple<bool, size_t, size_t, bool, size_t, size_t, size_t, size_t, bool> chain_info = dist_index.get_minimizer_distances(pos);
                        std::tuple<bool, size_t, size_t, bool, size_t, size_t, size_t, size_t, bool> encoded_chain_info = MIPayload::decode(MIPayload::encode(dist_index.get_minimizer_distances(pos)));

                        assert (chain_info == encoded_chain_info);

                        all_seeds[read].push_back({ pos, 0, std::get<0>(chain_info), std::get<1>(chain_info), std::get<2>(chain_info),
                           std::get<3>(chain_info), std::get<4>(chain_info), std::get<5>(chain_info), std::get<6>(chain_info), std::get<7>(chain_info), std::get<8>(chain_info)});

                    }
                }
                vector<vector<SnarlSeedClusterer::Cluster>> paired_clusters = clusterer.cluster_seeds(all_seeds, read_lim, fragment_lim); 
               
                vector<vector<pos_t>> fragment_clusters;

                for (size_t read_num = 0 ; read_num < 2 ; read_num ++) {
                    auto& one_read_clusters = paired_clusters[read_num];
                    if (one_read_clusters.size() > 0) {
                        for (size_t a = 0; a < one_read_clusters.size(); a++) {
                            // For each cluster -cluster this cluster to ensure that 
                            // there is only one
                            vector<size_t> clust = one_read_clusters[a].seeds;
                            size_t fragment_cluster = one_read_clusters[a].fragment;
                            if (fragment_cluster >= fragment_clusters.size()) {
                                fragment_clusters.resize(fragment_cluster+1);
                            }
                            
                            structures::UnionFind new_clusters (clust.size(), false);

                            for (size_t i1 = 0 ; i1 < clust.size() ; i1++) {
                                pos_t pos1 = all_seeds[read_num][clust[i1]].pos;
                                fragment_clusters[fragment_cluster].emplace_back(pos1);
                                size_t len1 = dist_index.node_length(get_id(pos1));
                                pos_t rev1 = make_pos_t(get_id(pos1), !is_rev(pos1),len1 - get_offset(pos1)-1); 

                                for (size_t b = 0 ; b < one_read_clusters.size() ; b++) {
                                    if (b != a) {
                                        //For each other cluster
                                        vector<size_t> clust2 = one_read_clusters[b].seeds;
                                        for (size_t i2 = 0 ; i2 < clust2.size() ; i2++) {
                                            //And each position in each other cluster,
                                            //make sure that this position is far away from i1
                                            pos_t pos2 = all_seeds[read_num][clust2[i2]].pos;
                                            size_t len2 = dist_index.node_length(get_id(pos2));
                                            pos_t rev2 = make_pos_t(get_id(pos2), 
                                                             !is_rev(pos2),
                                                             len2 - get_offset(pos2)-1); 

                                            int64_t dist1 = dist_index.min_distance(pos1, pos2);
                                            int64_t dist2 = dist_index.min_distance(pos1, rev2);
                                            int64_t dist3 = dist_index.min_distance(rev1, pos2);
                                            int64_t dist4 = dist_index.min_distance(rev1, rev2);
                                            int64_t dist = MinimumDistanceIndex::min_pos({dist1, 
                                                               dist2, dist3, dist4});
                                            if ( dist != -1 && dist <= read_lim) {
                                                dist_index.print_self();
                                                graph.serialize_to_file("testGraph");
                                                cerr << "These should have been in the same read cluster: " ;
                                                cerr << pos1 << " and " << pos2 << endl;
                                                cerr << dist1 << " " << dist2 << " " << dist3 << " " << dist4 << endl;
                                                REQUIRE(false);
                                            }
                                            
                                        }
                                    }
                                }
                                for (size_t i2 = 0 ; i2 < clust.size() ; i2++) {
                                    //For each position in the same cluster
                                    pos_t pos2 = all_seeds[read_num][clust[i2]].pos;
                                    size_t len2 = graph.get_length(graph.get_handle(get_id(pos2), false));
                                    pos_t rev2 = make_pos_t(get_id(pos2), 
                                                         !is_rev(pos2),
                                                         len2 - get_offset(pos2)-1); 
                                    int64_t dist1 = dist_index.min_distance(pos1, pos2);
                                    int64_t dist2 = dist_index.min_distance(pos1, rev2);
                                    int64_t dist3 = dist_index.min_distance(rev1, pos2);
                                    int64_t dist4 = dist_index.min_distance(rev1, rev2);
                                    int64_t dist = MinimumDistanceIndex::min_pos({dist1, 
                                                       dist2, dist3, dist4});
                                    if ( dist != -1 && dist <= read_lim) {
                                        new_clusters.union_groups(i1, i2);
                                    }

                                }
                            }
                            auto actual_clusters = new_clusters.all_groups();
                            if (actual_clusters.size() != 1) {
                                                dist_index.print_self();
                                graph.serialize_to_file("testGraph");
                                cerr << "These should be different read clusters: " << endl;
                                for (auto c : actual_clusters) {
                                    cerr << "cluster: " ; 
                                    for (size_t i1 : c) {
                                        cerr << all_seeds[read_num][clust[i1]].pos << " ";
                                    }
                                    cerr << endl;
                                }
                            }
                            REQUIRE(actual_clusters.size() == 1);
                        }
                    }
                }
                for (size_t a = 0; a < fragment_clusters.size(); a++) {
                    // For each cluster -cluster this cluster to ensure that 
                    // there is only one
                    vector<pos_t> clust = fragment_clusters[a];
                    
                    structures::UnionFind new_clusters (clust.size(), false);

                    for (size_t i1 = 0 ; i1 < clust.size() ; i1++) {
                        pos_t pos1 = clust[i1];
                        size_t len1 = graph.get_length(graph.get_handle(get_id(pos1), false));
                        pos_t rev1 = make_pos_t(get_id(pos1), 
                                            !is_rev(pos1),
                                            len1 - get_offset(pos1)-1); 

                        for (size_t b = 0 ; b < fragment_clusters.size() ; b++) {
                            if (b != a) {
                                //For each other cluster
                                vector<pos_t> clust2 = fragment_clusters[b];
                                for (size_t i2 = 0 ; i2 < clust2.size() ; i2++) {
                                    //And each position in each other cluster,
                                    //make sure that this position is far away from i1
                                    pos_t pos2 = clust2[i2];
                                    size_t len2 = graph.get_length(graph.get_handle(get_id(pos2), false));
                                    pos_t rev2 = make_pos_t(get_id(pos2), 
                                                     !is_rev(pos2),
                                                     len2 - get_offset(pos2)-1); 

                                    int64_t dist1 = dist_index.min_distance(pos1, pos2);
                                    int64_t dist2 = dist_index.min_distance(pos1, rev2);
                                    int64_t dist3 = dist_index.min_distance(rev1, pos2);
                                    int64_t dist4 = dist_index.min_distance(rev1, rev2);
                                    int64_t dist = MinimumDistanceIndex::min_pos({dist1, 
                                                       dist2, dist3, dist4});
                                    if ( dist != -1 && dist <= fragment_lim) {
                                        dist_index.print_self();
                                        graph.serialize_to_file("testGraph");
                                        cerr << "These should have been in the same fragment cluster: " ;
                                        cerr << pos1 << " and " << pos2 << endl;
                                        cerr << dist1 << " " << dist2 << " " << dist3 << " " << dist4 << endl;
                                        REQUIRE(false);
                                    }
                                    
                                }
                            }
                        }
                        for (size_t i2 = 0 ; i2 < clust.size() ; i2++) {
                            //For each position in the same cluster
                            pos_t pos2 = clust[i2];
                            size_t len2 = graph.get_length(graph.get_handle(get_id(pos2), false));
                            pos_t rev2 = make_pos_t(get_id(pos2), 
                                                 !is_rev(pos2),
                                                 len2 - get_offset(pos2)-1); 
                            int64_t dist1 = dist_index.min_distance(pos1, pos2);
                            int64_t dist2 = dist_index.min_distance(pos1, rev2);
                            int64_t dist3 = dist_index.min_distance(rev1, pos2);
                            int64_t dist4 = dist_index.min_distance(rev1, rev2);
                            int64_t dist = MinimumDistanceIndex::min_pos({dist1, 
                                               dist2, dist3, dist4});
                            if ( dist != -1 && dist <= fragment_lim) {
                                new_clusters.union_groups(i1, i2);
                            }

                        }
                    }
                    auto actual_clusters = new_clusters.all_groups();
                    if (actual_clusters.size() != 1) {
                                        dist_index.print_self();
                        graph.serialize_to_file("testGraph");
                        cerr << "These should be different fragment clusters: " << endl;
                        for (auto c : actual_clusters) {
                            cerr << "cluster: " ; 
                            for (size_t i1 : c) {
                                cerr << clust[i1] << " ";
                            }
                            cerr << endl;
                        }
                    }
                    REQUIRE(actual_clusters.size() == 1);
                }
            }
        }
    } //end test case
}
}
