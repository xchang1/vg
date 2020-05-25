/// \file multipath_alignment.cpp
///  
/// unit tests for multipath alignment construction and utility functions
///

#include <stdio.h>
#include <iostream>

#include <gbwt/dynamic_gbwt.h>

#include "../haplotypes.hpp"
#include "../json2pb.h"
#include <vg/vg.pb.h>
#include "../vg.hpp"
#include "../multipath_alignment.hpp"
#include "../utility.hpp"

#include "catch.hpp"

namespace vg {
    namespace unittest {
        
        TEST_CASE( "Multipath alignments correctly identify source subpaths",
                  "[alignment][multipath][mapping]" ) {
            
            SECTION( "Multipath alignment can identify source subpath in a linear subpath structure") {
                
                VG graph;
                
                Node* n1 = graph.create_node("GCA");
                Node* n2 = graph.create_node("T");
                Node* n3 = graph.create_node("G");
                Node* n4 = graph.create_node("CTGA");
                
                graph.create_edge(n1, n2);
                graph.create_edge(n1, n3);
                graph.create_edge(n2, n4);
                graph.create_edge(n3, n4);
                
                string read = string("GCATCTGA");
                multipath_alignment_t multipath_aln;
                multipath_aln.set_sequence(read);
                
                // add subpaths
                multipath_aln.add_subpath();
                multipath_aln.add_subpath();
                subpath_t* subpath0 = multipath_aln.mutable_subpath(0);
                subpath_t* subpath1 = multipath_aln.mutable_subpath(1);
                
                
                // set edges between subpaths
                subpath0->add_next(1);
                
                // designate mappings
                path_mapping_t* mapping0 = subpath0->mutable_path()->add_mapping();
                mapping0->mutable_position()->set_node_id(1);
                
                path_mapping_t* mapping1 = subpath0->mutable_path()->add_mapping();
                mapping1->mutable_position()->set_node_id(2);
                
                path_mapping_t* mapping2 = subpath1->mutable_path()->add_mapping();
                mapping2->mutable_position()->set_node_id(4);
                identify_start_subpaths(multipath_aln);
                
                REQUIRE(multipath_aln.start_size() == 1);
                REQUIRE(multipath_aln.start(0) == 0);
            }
            
            SECTION( "Multipath alignment can identify source subpath in a forked subpath structure") {
                
                VG graph;
                
                Node* n1 = graph.create_node("GCA");
                Node* n2 = graph.create_node("T");
                Node* n3 = graph.create_node("G");
                Node* n4 = graph.create_node("CTGA");
                
                graph.create_edge(n1, n2);
                graph.create_edge(n1, n3);
                graph.create_edge(n2, n4);
                graph.create_edge(n3, n4);
                
                string read = string("GCATCTGA");
                multipath_alignment_t multipath_aln;
                multipath_aln.set_sequence(read);
                
                // add subpaths
                multipath_aln.add_subpath();
                multipath_aln.add_subpath();
                multipath_aln.add_subpath();
                multipath_aln.add_subpath();
                subpath_t* subpath0 = multipath_aln.mutable_subpath(0);
                subpath_t* subpath1 = multipath_aln.mutable_subpath(1);
                subpath_t* subpath2 = multipath_aln.mutable_subpath(2);
                subpath_t* subpath3 = multipath_aln.mutable_subpath(3);
                
                // set edges between subpaths
                subpath0->add_next(1);
                subpath1->add_next(3);
                subpath2->add_next(3);
                
                // designate mappings
                path_mapping_t* mapping0 = subpath0->mutable_path()->add_mapping();
                mapping0->mutable_position()->set_node_id(1);
                
                path_mapping_t* mapping1 = subpath1->mutable_path()->add_mapping();
                mapping1->mutable_position()->set_node_id(2);
                
                path_mapping_t* mapping2 = subpath2->mutable_path()->add_mapping();
                mapping2->mutable_position()->set_node_id(3);
                
                path_mapping_t* mapping3 = subpath3->mutable_path()->add_mapping();
                mapping3->mutable_position()->set_node_id(4);
                
                
                identify_start_subpaths(multipath_aln);
                
                REQUIRE(multipath_aln.start_size() == 2);
                bool found_1 = false;
                bool found_2 = false;
                
                
                REQUIRE(multipath_aln.start_size() == 2);
                
                for (int i = 0; i < multipath_aln.start_size(); i++) {
                    if (multipath_aln.start(i) == 0) {
                        found_1 = true;
                    }
                    else if (multipath_aln.start(i) == 2) {
                        found_2 = true;
                    }
                    else {
                        REQUIRE(false);
                    }
                }
                REQUIRE( found_1 );
                REQUIRE( found_2 );
            }
        }
        
        TEST_CASE( "Multipath alignments correctly identifies optimal alignment within subpath DAG",
                  "[alignment][multipath][mapping]" ) {
            
            SECTION( "Multipath alignment can identify optimal alignment between two disjoint paths" ) {
                
                VG graph;
                
                Node* n1 = graph.create_node("GCA");
                Node* n2 = graph.create_node("T");
                Node* n3 = graph.create_node("G");
                Node* n4 = graph.create_node("CTGA");
                
                graph.create_edge(n1, n2);
                graph.create_edge(n1, n3);
                graph.create_edge(n2, n4);
                graph.create_edge(n3, n4);
                
                string read = string("T");
                multipath_alignment_t multipath_aln;
                multipath_aln.set_sequence(read);
                
                // add subpaths
                multipath_aln.add_subpath();
                multipath_aln.add_subpath();
                subpath_t* subpath0 = multipath_aln.mutable_subpath(0);
                subpath_t* subpath1 = multipath_aln.mutable_subpath(1);
                
                // set edges between subpaths
                
                // set scores
                subpath0->set_score(1);
                subpath1->set_score(0);
                
                // designate mappings
                path_mapping_t* mapping0 = subpath0->mutable_path()->add_mapping();
                mapping0->mutable_position()->set_node_id(2);
                
                path_mapping_t* mapping1 = subpath1->mutable_path()->add_mapping();
                mapping1->mutable_position()->set_node_id(3);
                
                // get optimal alignment
                identify_start_subpaths(multipath_aln);
                Alignment aln;
                optimal_alignment(multipath_aln, aln);
                
                REQUIRE(aln.path().mapping_size() == multipath_aln.subpath(0).path().mapping_size());
                for (int i = 0; i < aln.path().mapping_size(); i++) {
                    REQUIRE(aln.path().mapping(i).position().node_id() == multipath_aln.subpath(0).path().mapping(i).position().node_id());
                }
                REQUIRE(aln.score() == multipath_aln.subpath(0).score());
                
            }
            
            SECTION( "Multipath alignment can identify optimal alignment among paths that intersect" ) {
                
                VG graph;
                
                Node* n1 = graph.create_node("GCA");
                Node* n2 = graph.create_node("T");
                Node* n3 = graph.create_node("G");
                Node* n4 = graph.create_node("A");
                Node* n5 = graph.create_node("CTGA");
                
                graph.create_edge(n1, n3);
                graph.create_edge(n2, n3);
                graph.create_edge(n3, n4);
                graph.create_edge(n3, n5);
                
                string read = string("GCAGCTGA");
                multipath_alignment_t multipath_aln;
                multipath_aln.set_sequence(read);
                
                // add subpaths
                multipath_aln.add_subpath();
                multipath_aln.add_subpath();
                multipath_aln.add_subpath();
                multipath_aln.add_subpath();
                multipath_aln.add_subpath();
                subpath_t* subpath0 = multipath_aln.mutable_subpath(0);
                subpath_t* subpath1 = multipath_aln.mutable_subpath(1);
                subpath_t* subpath2 = multipath_aln.mutable_subpath(2);
                subpath_t* subpath3 = multipath_aln.mutable_subpath(3);
                subpath_t* subpath4 = multipath_aln.mutable_subpath(4);
                
                // set edges between subpaths
                subpath0->add_next(2);
                subpath1->add_next(2);
                subpath2->add_next(3);
                subpath2->add_next(4);
                
                // set scores
                subpath0->set_score(3);
                subpath1->set_score(0);
                subpath2->set_score(1);
                subpath3->set_score(0);
                subpath4->set_score(4);
                
                // designate mappings
                path_mapping_t* mapping0 = subpath0->mutable_path()->add_mapping();
                mapping0->mutable_position()->set_node_id(1);
                
                path_mapping_t* mapping1 = subpath1->mutable_path()->add_mapping();
                mapping1->mutable_position()->set_node_id(2);
                
                path_mapping_t* mapping2 = subpath2->mutable_path()->add_mapping();
                mapping2->mutable_position()->set_node_id(3);
                
                path_mapping_t* mapping3 = subpath3->mutable_path()->add_mapping();
                mapping3->mutable_position()->set_node_id(4);
                
                path_mapping_t* mapping4 = subpath4->mutable_path()->add_mapping();
                mapping4->mutable_position()->set_node_id(5);
                
                // get optimal alignment
                identify_start_subpaths(multipath_aln);
                Alignment aln;
                optimal_alignment(multipath_aln, aln);
                
                // follows correct path
                REQUIRE(aln.path().mapping_size() == 3);
                REQUIRE(aln.path().mapping(0).position().node_id() == 1);
                REQUIRE(aln.path().mapping(1).position().node_id() == 3);
                REQUIRE(aln.path().mapping(2).position().node_id() == 5);
                
                // has correct score
                REQUIRE(aln.score() == 8);
                
            }
            
            SECTION( "Multipath alignment correctly merge Mappings while finding optimal alignment" ) {
                
                VG graph;
                
                Node* n1 = graph.create_node("GCA");
                Node* n2 = graph.create_node("T");
                Node* n3 = graph.create_node("GTTGA");
                Node* n4 = graph.create_node("A");
                Node* n5 = graph.create_node("CTGA");
                
                graph.create_edge(n1, n3);
                graph.create_edge(n2, n3);
                graph.create_edge(n3, n4);
                graph.create_edge(n3, n5);
                
                string read = string("GCAGTGACTGA");
                multipath_alignment_t multipath_aln;
                multipath_aln.set_sequence(read);
                
                // add subpaths
                multipath_aln.add_subpath();
                multipath_aln.add_subpath();
                multipath_aln.add_subpath();
                multipath_aln.add_subpath();
                multipath_aln.add_subpath();
                multipath_aln.add_subpath();
                subpath_t* subpath0 = multipath_aln.mutable_subpath(0);
                subpath_t* subpath1 = multipath_aln.mutable_subpath(1);
                subpath_t* subpath2 = multipath_aln.mutable_subpath(2);
                subpath_t* subpath3 = multipath_aln.mutable_subpath(3);
                subpath_t* subpath4 = multipath_aln.mutable_subpath(4);
                subpath_t* subpath5 = multipath_aln.mutable_subpath(5);
                
                // set edges between subpaths
                subpath0->add_next(2);
                subpath1->add_next(2);
                subpath2->add_next(3);
                subpath3->add_next(4);
                subpath3->add_next(5);
                
                // set scores
                subpath0->set_score(3);
                subpath1->set_score(0);
                subpath2->set_score(3);
                subpath3->set_score(2);
                subpath4->set_score(0);
                subpath5->set_score(4);
                
                // designate mappings
                path_mapping_t* mapping0 = subpath0->mutable_path()->add_mapping();
                mapping0->mutable_position()->set_node_id(1);
                
                path_mapping_t* mapping1 = subpath1->mutable_path()->add_mapping();
                mapping1->mutable_position()->set_node_id(2);
                
                path_mapping_t* mapping2 = subpath2->mutable_path()->add_mapping();
                mapping2->mutable_position()->set_node_id(3);
                edit_t* edit2 = mapping2->add_edit();
                edit2->set_from_length(2);
                edit2->set_to_length(2);
                
                path_mapping_t* mapping3 = subpath3->mutable_path()->add_mapping();
                mapping3->mutable_position()->set_node_id(3);
                mapping3->mutable_position()->set_offset(2);
                edit_t* edit3 = mapping3->add_edit();
                edit3->set_from_length(3);
                edit3->set_to_length(3);
                
                path_mapping_t* mapping4 = subpath4->mutable_path()->add_mapping();
                mapping4->mutable_position()->set_node_id(4);
                
                path_mapping_t* mapping5 = subpath5->mutable_path()->add_mapping();
                mapping5->mutable_position()->set_node_id(5);
                
                // get optimal alignment
                identify_start_subpaths(multipath_aln);
                Alignment aln;
                optimal_alignment(multipath_aln, aln);
                
                // follows correct path
                REQUIRE(aln.path().mapping_size() == 3);
                REQUIRE(aln.path().mapping(0).position().node_id() == 1);
                REQUIRE(aln.path().mapping(1).position().node_id() == 3);
                REQUIRE(aln.path().mapping(2).position().node_id() == 5);
                
                // has correct ranks
                REQUIRE(aln.path().mapping(0).rank() == 1);
                REQUIRE(aln.path().mapping(1).rank() == 2);
                REQUIRE(aln.path().mapping(2).rank() == 3);
                
                // has correct score
                REQUIRE(aln.score() == 12);
                
                // has correct edits
                REQUIRE(aln.path().mapping(1).position().offset() == 0);
                REQUIRE(aln.path().mapping(1).edit(0).from_length() == 5);
                REQUIRE(aln.path().mapping(1).edit(0).to_length() == 5);
            }
            
            SECTION( "The optimal alignment can be forced to take low-scoring intervening subpaths" ) {
                
                string read = "GCAGTG";
                multipath_alignment_t multipath_aln;
                multipath_aln.set_sequence(read);
                
                // add subpaths
                multipath_aln.add_subpath();
                multipath_aln.add_subpath();
                multipath_aln.add_subpath();
                subpath_t* subpath0 = multipath_aln.mutable_subpath(0);
                subpath_t* subpath1 = multipath_aln.mutable_subpath(1);
                subpath_t* subpath2 = multipath_aln.mutable_subpath(2);
                
                // set edges between subpaths
                subpath0->add_next(1);
                subpath1->add_next(2);
                
                // set scores
                subpath0->set_score(3);
                subpath1->set_score(-4);
                subpath2->set_score(2);
                
                // designate mappings
                path_mapping_t* mapping0 = subpath0->mutable_path()->add_mapping();
                mapping0->mutable_position()->set_node_id(1);
                edit_t* edit0 = mapping0->add_edit();
                edit0->set_from_length(3);
                edit0->set_to_length(3);
                
                path_mapping_t* mapping1 = subpath1->mutable_path()->add_mapping();
                mapping1->mutable_position()->set_node_id(2);
                edit_t* edit1 = mapping0->add_edit();
                edit1->set_from_length(1);
                edit1->set_to_length(1);
                edit1->set_sequence("T");
                
                path_mapping_t* mapping2 = subpath2->mutable_path()->add_mapping();
                mapping2->mutable_position()->set_node_id(3);
                edit_t* edit2 = mapping2->add_edit();
                edit2->set_from_length(2);
                edit2->set_to_length(2);
                
                // get optimal alignment
                identify_start_subpaths(multipath_aln);
                Alignment aln;
                optimal_alignment(multipath_aln, aln, true);
                
                // follows correct path
                REQUIRE(aln.path().mapping_size() == 3);
                REQUIRE(aln.path().mapping(0).position().node_id() == 1);
                REQUIRE(aln.path().mapping(1).position().node_id() == 2);
                REQUIRE(aln.path().mapping(2).position().node_id() == 3);
                
                // has correct ranks
                REQUIRE(aln.path().mapping(0).rank() == 1);
                REQUIRE(aln.path().mapping(1).rank() == 2);
                REQUIRE(aln.path().mapping(2).rank() == 3);
                
                // has correct score
                REQUIRE(aln.score() == 1);
            }
        }
        
        TEST_CASE("Multipath alignment correctly identifies suboptimal alignments") {
        
            SECTION( "Multipath alignment can identify optimal and suboptimal alignment between two disjoint paths" ) {
                
                VG graph;
                
                Node* n1 = graph.create_node("GCA");
                Node* n2 = graph.create_node("T");
                Node* n3 = graph.create_node("G");
                Node* n4 = graph.create_node("CTGA");
                
                graph.create_edge(n1, n2);
                graph.create_edge(n1, n3);
                graph.create_edge(n2, n4);
                graph.create_edge(n3, n4);
                
                string read = string("T");
                multipath_alignment_t multipath_aln;
                multipath_aln.set_sequence(read);
                
                // add subpaths
                multipath_aln.add_subpath();
                multipath_aln.add_subpath();
                subpath_t* subpath0 = multipath_aln.mutable_subpath(0);
                subpath_t* subpath1 = multipath_aln.mutable_subpath(1);
                
                // set edges between subpaths
                
                // set scores
                subpath0->set_score(1);
                subpath1->set_score(0);
                
                // designate mappings
                path_mapping_t* mapping0 = subpath0->mutable_path()->add_mapping();
                mapping0->mutable_position()->set_node_id(2);
                
                path_mapping_t* mapping1 = subpath1->mutable_path()->add_mapping();
                mapping1->mutable_position()->set_node_id(3);
                
                // get top 10 alignments
                auto top10 = optimal_alignments(multipath_aln, 10);
                
                // There should be two alignments
                REQUIRE(top10.size() == 2);
                
                // The best alignment should be the same one we got for this example before as optimal
                REQUIRE(top10[0].path().mapping_size() == multipath_aln.subpath(0).path().mapping_size());
                for (int i = 0; i < top10[0].path().mapping_size(); i++) {
                    REQUIRE(top10[0].path().mapping(i).position().node_id() == multipath_aln.subpath(0).path().mapping(i).position().node_id());
                }
                REQUIRE(top10[0].score() == multipath_aln.subpath(0).score());
                
                // The second best alignment should be that 0-score subpath
                REQUIRE(top10[1].path().mapping_size() == multipath_aln.subpath(1).path().mapping_size());
                for (int i = 0; i < top10[1].path().mapping_size(); i++) {
                    REQUIRE(top10[1].path().mapping(i).position().node_id() == multipath_aln.subpath(1).path().mapping(i).position().node_id());
                }
                REQUIRE(top10[1].score() == multipath_aln.subpath(1).score());
                
            }
            
            SECTION( "Multipath alignment correctly merge Mappings while finding optimal and suboptimal alignments" ) {
                
                VG graph;
                
                Node* n1 = graph.create_node("GCA");
                Node* n2 = graph.create_node("T");
                Node* n3 = graph.create_node("GTTGA");
                Node* n4 = graph.create_node("A");
                Node* n5 = graph.create_node("CTGA");
                
                graph.create_edge(n1, n3);
                graph.create_edge(n2, n3);
                graph.create_edge(n3, n4);
                graph.create_edge(n3, n5);
                
                string read = string("GCAGTGACTGA");
                multipath_alignment_t multipath_aln;
                multipath_aln.set_sequence(read);
                
                // add subpaths
                multipath_aln.add_subpath();
                multipath_aln.add_subpath();
                multipath_aln.add_subpath();
                multipath_aln.add_subpath();
                multipath_aln.add_subpath();
                multipath_aln.add_subpath();
                subpath_t* subpath0 = multipath_aln.mutable_subpath(0);
                subpath_t* subpath1 = multipath_aln.mutable_subpath(1);
                subpath_t* subpath2 = multipath_aln.mutable_subpath(2);
                subpath_t* subpath3 = multipath_aln.mutable_subpath(3);
                subpath_t* subpath4 = multipath_aln.mutable_subpath(4);
                subpath_t* subpath5 = multipath_aln.mutable_subpath(5);
                
                // set edges between subpaths
                subpath0->add_next(2);
                subpath1->add_next(2);
                subpath2->add_next(3);
                subpath3->add_next(4);
                subpath3->add_next(5);
                
                // set scores
                subpath0->set_score(3);
                subpath1->set_score(0);
                subpath2->set_score(3);
                subpath3->set_score(2);
                subpath4->set_score(0);
                subpath5->set_score(4);
                
                // designate mappings
                // TODO: Note that these edits aren't quite realistic
                path_mapping_t* mapping0 = subpath0->mutable_path()->add_mapping();
                mapping0->mutable_position()->set_node_id(1);
                
                path_mapping_t* mapping1 = subpath1->mutable_path()->add_mapping();
                mapping1->mutable_position()->set_node_id(2);
                
                path_mapping_t* mapping2 = subpath2->mutable_path()->add_mapping();
                mapping2->mutable_position()->set_node_id(3);
                edit_t* edit2 = mapping2->add_edit();
                edit2->set_from_length(2);
                edit2->set_to_length(2);
                
                path_mapping_t* mapping3 = subpath3->mutable_path()->add_mapping();
                mapping3->mutable_position()->set_node_id(3);
                mapping3->mutable_position()->set_offset(2);
                edit_t* edit3 = mapping3->add_edit();
                edit3->set_from_length(3);
                edit3->set_to_length(3);
                
                path_mapping_t* mapping4 = subpath4->mutable_path()->add_mapping();
                mapping4->mutable_position()->set_node_id(4);
                
                path_mapping_t* mapping5 = subpath5->mutable_path()->add_mapping();
                mapping5->mutable_position()->set_node_id(5);
                
                // get top 10 alignments
                identify_start_subpaths(multipath_aln);
                vector<Alignment> top10 = optimal_alignments(multipath_aln, 10);
                
                SECTION("Top alignment is correct") {
                    // Exists
                    REQUIRE(top10.size() > 0);
                    
                    // follows correct path
                    REQUIRE(top10[0].path().mapping_size() == 3);
                    REQUIRE(top10[0].path().mapping(0).position().node_id() == 1);
                    REQUIRE(top10[0].path().mapping(1).position().node_id() == 3);
                    REQUIRE(top10[0].path().mapping(2).position().node_id() == 5);
                    
                    // has correct score
                    REQUIRE(top10[0].score() == 12);
                    
                    // has correct edits
                    REQUIRE(top10[0].path().mapping(1).position().offset() == 0);
                    REQUIRE(top10[0].path().mapping(1).edit(0).from_length() == 5);
                    REQUIRE(top10[0].path().mapping(1).edit(0).to_length() == 5);
                }
                
                SECTION("Secondary alignment is correct") {
                    
                    // Exists
                    REQUIRE(top10.size() > 1);
                    
                    // follows correct path
                    REQUIRE(top10[1].path().mapping_size() == 3);
                    REQUIRE(top10[1].path().mapping(0).position().node_id() == 2);
                    REQUIRE(top10[1].path().mapping(1).position().node_id() == 3);
                    REQUIRE(top10[1].path().mapping(2).position().node_id() == 5);
                    
                    // has correct score
                    REQUIRE(top10[1].score() == 9);
                    
                    // has correct edits
                    REQUIRE(top10[1].path().mapping(1).position().offset() == 0);
                    REQUIRE(top10[1].path().mapping(1).edit(0).from_length() == 5);
                    REQUIRE(top10[1].path().mapping(1).edit(0).to_length() == 5);
                
                }
                
                SECTION("Tertiary alignment is correct") {
                    
                    // Exists
                    REQUIRE(top10.size() > 2);
                    
                    // follows correct path
                    REQUIRE(top10[2].path().mapping_size() == 3);
                    REQUIRE(top10[2].path().mapping(0).position().node_id() == 1);
                    REQUIRE(top10[2].path().mapping(1).position().node_id() == 3);
                    REQUIRE(top10[2].path().mapping(2).position().node_id() == 4);
                    
                    // has correct score
                    REQUIRE(top10[2].score() == 8);
                    
                    // has correct edits
                    REQUIRE(top10[2].path().mapping(1).position().offset() == 0);
                    REQUIRE(top10[2].path().mapping(1).edit(0).from_length() == 5);
                    REQUIRE(top10[2].path().mapping(1).edit(0).to_length() == 5);
                
                }
                
                SECTION("Quaternary alignment is correct") {
                    
                    // Exists
                    REQUIRE(top10.size() > 3);
                    
                    // follows correct path
                    REQUIRE(top10[3].path().mapping_size() == 3);
                    REQUIRE(top10[3].path().mapping(0).position().node_id() == 2);
                    REQUIRE(top10[3].path().mapping(1).position().node_id() == 3);
                    REQUIRE(top10[3].path().mapping(2).position().node_id() == 4);
                    
                    // has correct score
                    REQUIRE(top10[3].score() == 5);
                    
                    // has correct edits
                    REQUIRE(top10[3].path().mapping(1).position().offset() == 0);
                    REQUIRE(top10[3].path().mapping(1).edit(0).from_length() == 5);
                    REQUIRE(top10[3].path().mapping(1).edit(0).to_length() == 5);
                
                }
                
                SECTION("Quinary alignment does not exist") {
                    REQUIRE(top10.size() < 5);
                }
            }
        
        }
        
        TEST_CASE("Multipath alignment correctly identifies haplotype consistent alignments") {
        
            // Make a wrapper to make the compiler shut up about narrowing
            auto visit = [](id_t id, bool orientation) {
                return (gbwt::vector_type::value_type) gbwt::Node::encode(id, orientation);
            };
            auto end = (gbwt::vector_type::value_type) gbwt::ENDMARKER;
        
            SECTION( "Multipath alignment can identify a single haplotype consistent alignment" ) {
                
                VG graph;
                
                Node* n1 = graph.create_node("GCA");
                Node* n2 = graph.create_node("T");
                Node* n3 = graph.create_node("G");
                Node* n4 = graph.create_node("CTGA");
                
                graph.create_edge(n1, n2);
                graph.create_edge(n1, n3);
                graph.create_edge(n2, n4);
                graph.create_edge(n3, n4);
                
                string read = string("GCATCTGA");
                multipath_alignment_t multipath_aln;
                multipath_aln.set_sequence(read);
                
                // add subpaths
                multipath_aln.add_subpath();
                multipath_aln.add_subpath();
                multipath_aln.add_subpath();
                multipath_aln.add_subpath();
                subpath_t* subpath0 = multipath_aln.mutable_subpath(0);
                subpath_t* subpath1 = multipath_aln.mutable_subpath(1);
                subpath_t* subpath2 = multipath_aln.mutable_subpath(2);
                subpath_t* subpath3 = multipath_aln.mutable_subpath(3);
                
                // set edges between subpaths
                subpath0->add_next(1);
                subpath0->add_next(2);
                subpath1->add_next(3);
                subpath2->add_next(3);
                
                // set scores
                subpath0->set_score(3);
                subpath1->set_score(1);
                subpath2->set_score(-1);
                subpath3->set_score(4);
                
                // designate mappings
                path_mapping_t* mapping0 = subpath0->mutable_path()->add_mapping();
                mapping0->mutable_position()->set_node_id(1);
                edit_t* edit0 = mapping0->add_edit();
                edit0->set_from_length(3);
                edit0->set_to_length(3);
                
                path_mapping_t* mapping1 = subpath1->mutable_path()->add_mapping();
                mapping1->mutable_position()->set_node_id(2);
                edit_t* edit1 = mapping1->add_edit();
                edit1->set_from_length(1);
                edit1->set_to_length(1);
                
                path_mapping_t* mapping2 = subpath2->mutable_path()->add_mapping();
                mapping2->mutable_position()->set_node_id(3);
                edit_t* edit2 = mapping2->add_edit();
                edit2->set_from_length(1);
                edit2->set_to_length(1);
                edit2->set_sequence("T");
                
                path_mapping_t* mapping3 = subpath3->mutable_path()->add_mapping();
                mapping3->mutable_position()->set_node_id(4);
                edit_t* edit3 = mapping3->add_edit();
                edit3->set_from_length(4);
                edit3->set_to_length(4);
                
                // Create haplotype database
                gbwt::Verbosity::set(gbwt::Verbosity::SILENT);
                gbwt::DynamicGBWT gbwt_index;
             
                // Make and insert a haplotype
                gbwt::vector_type hap1fwd = {
                    visit(1, false),
                    visit(3, false),
                    visit(4, false),
                    end
                };
                gbwt_index.insert(hap1fwd);
                // Do both orienatations
                gbwt::vector_type hap1rev = {
                    visit(4, true),
                    visit(3, true),
                    visit(1, true),
                    end
                };
                gbwt_index.insert(hap1rev);
                
                // Make a ScoreProvider
                haplo::GBWTScoreProvider<gbwt::DynamicGBWT> provider(gbwt_index);
                
                SECTION( "We can find the consistent suboptimal alignment" ) {

                    // get haplotype consistent alignments
                    auto consistent = haplotype_consistent_alignments(multipath_aln, provider, 0, 0);
                    
                    // There should be just one alignment
                    REQUIRE(consistent.size() == 1);
                    
                    // The alignment should have score 3 + 4 - 1 = 6
                    REQUIRE(consistent[0].score() == 6);
                    
                    // The alignment should be nodes 1, 3, 4, all forward
                    REQUIRE(consistent[0].path().mapping_size() == 3);
                    REQUIRE(consistent[0].path().mapping(0).position().node_id() == 1);
                    REQUIRE(consistent[0].path().mapping(0).position().is_reverse() == false);
                    REQUIRE(consistent[0].path().mapping(1).position().node_id() == 3);
                    REQUIRE(consistent[0].path().mapping(1).position().is_reverse() == false);
                    REQUIRE(consistent[0].path().mapping(2).position().node_id() == 4);
                    REQUIRE(consistent[0].path().mapping(2).position().is_reverse() == false);
                    
                }
                
                SECTION( "We can find the optimal alignment when asked for it" ) {
                    auto optimal_and_consistent = haplotype_consistent_alignments(multipath_aln, provider, 0, 0, true);
                    
                    // There should be two alignments
                    REQUIRE(optimal_and_consistent.size() == 2);
                    
                    // The first should have the optimal score of 3 + 4 + 1 = 8
                    REQUIRE(optimal_and_consistent[0].score() == 8);
                    // The second should be the haplotype-consistent alignment with score 6.
                    REQUIRE(optimal_and_consistent[1].score() == 6);
                }
            }
            
            SECTION( "Multipath alignment can identify multiple haplotype consistent alignments" ) {
                
                VG graph;
                
                Node* n1 = graph.create_node("GCA");
                Node* n2 = graph.create_node("T");
                Node* n3 = graph.create_node("G");
                Node* n4 = graph.create_node("CTGA");
                Node* n5 = graph.create_node("ACAC");
                Node* n6 = graph.create_node("CC");
                
                graph.create_edge(n1, n2);
                graph.create_edge(n1, n3);
                graph.create_edge(n2, n4);
                graph.create_edge(n3, n4);
                graph.create_edge(n4, n5);
                graph.create_edge(n4, n6);
                
                string read = string("GCATCTGAAC");
                multipath_alignment_t multipath_aln;
                multipath_aln.set_sequence(read);
                
                // add subpaths
                multipath_aln.add_subpath();
                multipath_aln.add_subpath();
                multipath_aln.add_subpath();
                multipath_aln.add_subpath();
                multipath_aln.add_subpath();
                multipath_aln.add_subpath();
                subpath_t* subpath0 = multipath_aln.mutable_subpath(0);
                subpath_t* subpath1 = multipath_aln.mutable_subpath(1);
                subpath_t* subpath2 = multipath_aln.mutable_subpath(2);
                subpath_t* subpath3 = multipath_aln.mutable_subpath(3);
                subpath_t* subpath4 = multipath_aln.mutable_subpath(4);
                subpath_t* subpath5 = multipath_aln.mutable_subpath(5);
                
                // set edges between subpaths
                subpath0->add_next(1);
                subpath0->add_next(2);
                subpath1->add_next(3);
                subpath2->add_next(3);
                subpath3->add_next(4);
                subpath3->add_next(5);
                
                // set scores
                subpath0->set_score(3);
                subpath1->set_score(1);
                subpath2->set_score(-1);
                subpath3->set_score(4);
                subpath4->set_score(2);
                subpath5->set_score(0);
                
                // designate mappings
                path_mapping_t* mapping0 = subpath0->mutable_path()->add_mapping();
                mapping0->mutable_position()->set_node_id(1);
                edit_t* edit0 = mapping0->add_edit();
                edit0->set_from_length(3);
                edit0->set_to_length(3);
                
                path_mapping_t* mapping1 = subpath1->mutable_path()->add_mapping();
                mapping1->mutable_position()->set_node_id(2);
                edit_t* edit1 = mapping1->add_edit();
                edit1->set_from_length(1);
                edit1->set_to_length(1);
                
                path_mapping_t* mapping2 = subpath2->mutable_path()->add_mapping();
                mapping2->mutable_position()->set_node_id(3);
                edit_t* edit2 = mapping2->add_edit();
                edit2->set_from_length(1);
                edit2->set_to_length(1);
                edit2->set_sequence("T");
                
                path_mapping_t* mapping3 = subpath3->mutable_path()->add_mapping();
                mapping3->mutable_position()->set_node_id(4);
                edit_t* edit3 = mapping3->add_edit();
                edit3->set_from_length(4);
                edit3->set_to_length(4);
                
                path_mapping_t* mapping4 = subpath4->mutable_path()->add_mapping();
                mapping4->mutable_position()->set_node_id(5);
                edit_t* edit4 = mapping4->add_edit();
                edit4->set_from_length(2);
                edit4->set_to_length(2);
                
                path_mapping_t* mapping5 = subpath5->mutable_path()->add_mapping();
                mapping5->mutable_position()->set_node_id(6);
                mapping5->add_edit();
                mapping5->add_edit();
                edit_t* edit5a = mapping5->mutable_edit(0);
                edit_t* edit5b = mapping5->mutable_edit(1);
                edit5a->set_from_length(1);
                edit5a->set_to_length(1);
                edit5a->set_sequence("A");
                edit5b->set_from_length(1);
                edit5b->set_to_length(1);
                
                // Create haplotype database
                gbwt::Verbosity::set(gbwt::Verbosity::SILENT);
                gbwt::DynamicGBWT gbwt_index;
             
                // Make and insert a haplotype
                gbwt::vector_type hap1fwd = {
                    visit(1, false),
                    visit(3, false),
                    visit(4, false),
                    visit(6, false),
                    end
                };
                gbwt_index.insert(hap1fwd);
                // Do both orienatations
                gbwt::vector_type hap1rev = {
                    visit(6, true),
                    visit(4, true),
                    visit(3, true),
                    visit(1, true),
                    end
                };
                gbwt_index.insert(hap1rev);
                
                // More haplotypes. All except 2, 6
                gbwt::vector_type hap2fwd = {
                    visit(1, false),
                    visit(3, false),
                    visit(4, false),
                    visit(5, false),
                    end
                };
                gbwt_index.insert(hap2fwd);
                gbwt::vector_type hap2rev = {
                    visit(5, true),
                    visit(4, true),
                    visit(3, true),
                    visit(1, true),
                    end
                };
                gbwt_index.insert(hap2rev);
                gbwt::vector_type hap3fwd = {
                    visit(1, false),
                    visit(2, false),
                    visit(4, false),
                    visit(5, false),
                    end
                };
                gbwt_index.insert(hap3fwd);
                gbwt::vector_type hap3rev = {
                    visit(5, true),
                    visit(4, true),
                    visit(2, true),
                    visit(1, true),
                    end
                };
                gbwt_index.insert(hap3rev);
                
                // Make a ScoreProvider
                haplo::GBWTScoreProvider<gbwt::DynamicGBWT> provider(gbwt_index);

                // get haplotype consistent alignments
                auto consistent = haplotype_consistent_alignments(multipath_aln, provider, 0, 0);
                
                // There should be 3 alignments
                REQUIRE(consistent.size() == 3);
                
                for (auto& aln : consistent) {
                    // Each should have 4 mappings
                    REQUIRE(aln.path().mapping_size() == 4);
                    // None of them should go to nodes 2 and 6
                    REQUIRE(!(aln.path().mapping(1).position().node_id() == 2 && aln.path().mapping(3).position().node_id() == 6));
                }
            }
        }
        
        TEST_CASE( "Reverse complementing multipath alignments works correctly",
                  "[alignment][multipath][mapping]" ) {
            
            VG graph;
            
            Node* n1 = graph.create_node("GCA");
            Node* n2 = graph.create_node("T");
            Node* n3 = graph.create_node("GCC");
            Node* n4 = graph.create_node("A");
            Node* n5 = graph.create_node("CTGA");
            
            graph.create_edge(n1, n3);
            graph.create_edge(n2, n3);
            graph.create_edge(n3, n4);
            graph.create_edge(n3, n5);
            
            string read = string("CACCCTGA");
            multipath_alignment_t multipath_aln;
            multipath_aln.set_sequence(read);
            
            // add subpaths
            multipath_aln.add_subpath();
            multipath_aln.add_subpath();
            multipath_aln.add_subpath();
            multipath_aln.add_subpath();
            multipath_aln.add_subpath();
            multipath_aln.add_subpath();
            subpath_t* subpath0 = multipath_aln.mutable_subpath(0);
            subpath_t* subpath1 = multipath_aln.mutable_subpath(1);
            subpath_t* subpath2 = multipath_aln.mutable_subpath(2);
            subpath_t* subpath3 = multipath_aln.mutable_subpath(3);
            subpath_t* subpath4 = multipath_aln.mutable_subpath(4);
            subpath_t* subpath5 = multipath_aln.mutable_subpath(5);
            
            // set edges between subpaths
            subpath0->add_next(2);
            subpath0->add_next(5);
            subpath1->add_next(2);
            subpath2->add_next(3);
            subpath2->add_next(4);
            
            // set scores
            subpath0->set_score(3);
            subpath1->set_score(0);
            subpath2->set_score(1);
            subpath3->set_score(0);
            subpath4->set_score(4);
            subpath5->set_score(0);
            
            // designate mappings
            path_mapping_t* mapping0 = subpath0->mutable_path()->add_mapping();
            mapping0->mutable_position()->set_node_id(1);
            mapping0->mutable_position()->set_offset(1);
            edit_t* edit00 = mapping0->add_edit();
            edit00->set_from_length(2);
            edit00->set_to_length(2);
            
            path_mapping_t* mapping1 = subpath1->mutable_path()->add_mapping();
            mapping1->mutable_position()->set_node_id(2);
            mapping1->mutable_position()->set_offset(1);
            edit_t* edit10 = mapping1->add_edit();
            edit10->set_from_length(0);
            edit10->set_to_length(2);
            edit10->set_sequence("CA");
            
            path_mapping_t* mapping2 = subpath2->mutable_path()->add_mapping();
            mapping2->mutable_position()->set_node_id(3);
            mapping2->add_edit();
            mapping2->add_edit();
            edit_t* edit20 = mapping2->mutable_edit(0);
            edit_t* edit21 = mapping2->mutable_edit(1);
            edit20->set_from_length(1);
            edit20->set_to_length(0);
            edit21->set_from_length(2);
            edit21->set_to_length(2);
            
            path_mapping_t* mapping3 = subpath3->mutable_path()->add_mapping();
            mapping3->mutable_position()->set_node_id(4);
            edit_t* edit30 = mapping3->add_edit();
            edit30->set_from_length(0);
            edit30->set_to_length(4);
            edit30->set_sequence("CTGA");
            
            path_mapping_t* mapping4 = subpath4->mutable_path()->add_mapping();
            mapping4->mutable_position()->set_node_id(5);
            edit_t* edit40 = mapping4->add_edit();
            edit40->set_from_length(4);
            edit40->set_to_length(4);
            
            path_mapping_t* mapping5 = subpath5->mutable_path()->add_mapping();
            mapping5->mutable_position()->set_node_id(3);
            mapping5->mutable_position()->set_offset(3);
            edit_t* edit50 = mapping5->add_edit();
            edit50->set_from_length(0);
            edit50->set_to_length(6);
            edit50->set_sequence("CCCTGA");
            
            // identify starts
            multipath_aln.add_start(0);
            multipath_aln.add_start(1);
            
            multipath_alignment_t rc_multipath_aln;
            
            auto node_length = [&graph](int64_t node_id) { return graph.get_node(node_id)->sequence().length(); };
            rev_comp_multipath_alignment(multipath_aln, node_length, rc_multipath_aln);
            
            // all subpaths preserved
            REQUIRE(rc_multipath_aln.subpath_size() == multipath_aln.subpath_size());
            
            // subpaths are in reversed order
            for (int i = 0; i < multipath_aln.subpath_size(); i++) {
                REQUIRE(rc_multipath_aln.subpath(i).score() == multipath_aln.subpath(multipath_aln.subpath_size() - i - 1).score());
                REQUIRE(rc_multipath_aln.subpath(i).path().mapping(0).position().node_id() == multipath_aln.subpath(multipath_aln.subpath_size() - i - 1).path().mapping(0).position().node_id());
            }
            
            // positions are on reverse
            for (int i = 0; i < rc_multipath_aln.subpath_size(); i++) {
                for (int j = 0; j < rc_multipath_aln.subpath(i).path().mapping_size(); j++) {
                    REQUIRE(rc_multipath_aln.subpath(i).path().mapping(j).position().is_reverse());
                }
            }
            
            // offsets set correctly
            REQUIRE(rc_multipath_aln.subpath(0).path().mapping(0).position().offset() == 0);
            REQUIRE(rc_multipath_aln.subpath(1).path().mapping(0).position().offset() == 0);
            REQUIRE(rc_multipath_aln.subpath(2).path().mapping(0).position().offset() == 1);
            REQUIRE(rc_multipath_aln.subpath(3).path().mapping(0).position().offset() == 0);
            REQUIRE(rc_multipath_aln.subpath(4).path().mapping(0).position().offset() == 0);
            REQUIRE(rc_multipath_aln.subpath(5).path().mapping(0).position().offset() == 0);
            
            // edits are correct
            REQUIRE(rc_multipath_aln.subpath(0).path().mapping(0).edit(0).from_length() == 0);
            REQUIRE(rc_multipath_aln.subpath(0).path().mapping(0).edit(0).to_length() == 6);
            REQUIRE(rc_multipath_aln.subpath(0).path().mapping(0).edit(0).sequence() == "TCAGGG");
            
            REQUIRE(rc_multipath_aln.subpath(1).path().mapping(0).edit(0).from_length() == 4);
            REQUIRE(rc_multipath_aln.subpath(1).path().mapping(0).edit(0).to_length() == 4);
            
            REQUIRE(rc_multipath_aln.subpath(2).path().mapping(0).edit(0).from_length() == 0);
            REQUIRE(rc_multipath_aln.subpath(2).path().mapping(0).edit(0).to_length() == 4);
            REQUIRE(rc_multipath_aln.subpath(2).path().mapping(0).edit(0).sequence() == "TCAG");
            
            REQUIRE(rc_multipath_aln.subpath(3).path().mapping(0).edit(0).from_length() == 2);
            REQUIRE(rc_multipath_aln.subpath(3).path().mapping(0).edit(0).to_length() == 2);
            
            REQUIRE(rc_multipath_aln.subpath(3).path().mapping(0).edit(1).from_length() == 1);
            REQUIRE(rc_multipath_aln.subpath(3).path().mapping(0).edit(1).to_length() == 0);
            
            REQUIRE(rc_multipath_aln.subpath(4).path().mapping(0).edit(0).from_length() == 0);
            REQUIRE(rc_multipath_aln.subpath(4).path().mapping(0).edit(0).to_length() == 2);
            REQUIRE(rc_multipath_aln.subpath(4).path().mapping(0).edit(0).sequence() == "TG");
            
            REQUIRE(rc_multipath_aln.subpath(5).path().mapping(0).edit(0).from_length() == 2);
            REQUIRE(rc_multipath_aln.subpath(5).path().mapping(0).edit(0).to_length() == 2);
            
            // starts are switched correctly
            for (int i : {0, 1, 2}) {
                bool found = false;
                for (int64_t j = 0; j < rc_multipath_aln.start_size(); j++) {
                    found = found || rc_multipath_aln.start(j) == i;
                }
                REQUIRE(found);
            }
            
            // now check that the in place reverse complement works correctly
            // (it should match the already verified reverse complement)
            rev_comp_multipath_alignment_in_place(&multipath_aln, node_length);
            
            REQUIRE(rc_multipath_aln.subpath_size() == multipath_aln.subpath_size());
            
            for (int64_t i = 0; i < multipath_aln.subpath_size(); i++) {
                const subpath_t& subpath_1 = multipath_aln.subpath(i);
                const subpath_t& subpath_2 = rc_multipath_aln.subpath(i);
                REQUIRE(subpath_1.score() == subpath_2.score());
                const path_t& path_1 = subpath_1.path();
                const path_t& path_2 = subpath_2.path();
                REQUIRE(path_1.mapping_size() == path_2.mapping_size());
                for (int64_t j = 0; j < path_1.mapping_size(); j++) {
                    const path_mapping_t& mapping_1 = path_1.mapping(j);
                    const path_mapping_t& mapping_2 = path_2.mapping(j);
                    REQUIRE(mapping_1.edit_size() == mapping_2.edit_size());
                    for (int64_t k = 0; k < mapping_1.edit_size(); k++) {
                        const edit_t& edit_1 = mapping_1.edit(k);
                        const edit_t& edit_2 = mapping_2.edit(k);
                        REQUIRE(edit_1.from_length() == edit_2.from_length());
                        REQUIRE(edit_1.to_length() == edit_2.to_length());
                        REQUIRE(edit_1.sequence() == edit_2.sequence());
                    }
                }
                REQUIRE(subpath_1.next_size() == subpath_2.next_size());
                // these might be out of order
                for (int64_t j = 0; j < subpath_1.next_size(); j++) {
                    bool found = false;
                    for (int64_t k = 0; k < subpath_2.next_size(); k++) {
                        found = found || subpath_1.next(j) == subpath_2.next(k);
                    }
                    REQUIRE(found);
                }
            }
            
            // starts are switched correctly
            REQUIRE(rc_multipath_aln.start_size() == multipath_aln.start_size());
            for (int64_t i = 0; i < rc_multipath_aln.start_size(); i++) {
                bool found = false;
                for (int64_t j = 0; j < multipath_aln.start_size(); j++) {
                    found = found || rc_multipath_aln.start(i) == multipath_aln.start(j);
                }
                REQUIRE(found);
            }
        }
        
        TEST_CASE( "Algorithm returns correct connected components for multipath_alignment_ts",
                  "[alignment][multipath]" ){
            
            SECTION("Works for a single component multipath_alignment_t") {
                
                
                string read = "CACCCTGA";
                multipath_alignment_t multipath_aln;
                multipath_aln.set_sequence(read);
                
                // add subpaths
                multipath_aln.add_subpath();
                multipath_aln.add_subpath();
                multipath_aln.add_subpath();
                subpath_t* subpath0 = multipath_aln.mutable_subpath(0);
                subpath_t* subpath1 = multipath_aln.mutable_subpath(1);
                subpath_t* subpath2 = multipath_aln.mutable_subpath(2);
                
                // set edges between subpaths
                subpath0->add_next(1);
                subpath1->add_next(2);
                
                // set scores
                subpath0->set_score(3);
                subpath1->set_score(4);
                subpath2->set_score(7);
                
                auto comps = connected_components(multipath_aln);
                
                REQUIRE(comps.size() == 1);
                
                bool found_1 = false;
                for (size_t i = 0; i < comps.size(); i++) {
                    bool this_is_1 = true;
                    for (int j : {0, 1, 2}) {
                        bool found_this_item = false;
                        for (size_t k = 0; k < comps[i].size(); k++) {
                            if (comps[i][k] == j) {
                                found_this_item = true;
                                break;
                            }
                        }
                        this_is_1 = this_is_1 && found_this_item;
                    }
                    found_1 = this_is_1 || found_1;
                }
                
                REQUIRE(found_1);
            }
            
            SECTION("Works for a two component multipath_alignment_t") {
                
                string read = "CACCCTGA";
                multipath_alignment_t multipath_aln;
                multipath_aln.set_sequence(read);
                
                // add subpaths
                multipath_aln.add_subpath();
                multipath_aln.add_subpath();
                multipath_aln.add_subpath();
                subpath_t* subpath0 = multipath_aln.mutable_subpath(0);
                subpath_t* subpath1 = multipath_aln.mutable_subpath(1);
                subpath_t* subpath2 = multipath_aln.mutable_subpath(2);
                
                // set edges between subpaths
                subpath0->add_next(1);
                
                // set scores
                subpath0->set_score(3);
                subpath1->set_score(4);
                subpath2->set_score(7);
                
                auto comps = connected_components(multipath_aln);
                
                REQUIRE(comps.size() == 2);
                
                bool found_1 = false;
                bool found_2 = false;
                for (size_t i = 0; i < comps.size(); i++) {
                    bool is_this_comp;
                    is_this_comp = true;
                    for (int j : {0, 1}) {
                        bool found_this_item = false;
                        for (size_t k = 0; k < comps[i].size(); k++) {
                            if (comps[i][k] == j) {
                                found_this_item = true;
                                break;
                            }
                        }
                        is_this_comp = is_this_comp && found_this_item;
                    }
                    found_1 = is_this_comp || found_1;
                    
                    is_this_comp = true;
                    for (int j : {2}) {
                        bool found_this_item = false;
                        for (size_t k = 0; k < comps[i].size(); k++) {
                            if (comps[i][k] == j) {
                                found_this_item = true;
                                break;
                            }
                        }
                        is_this_comp = is_this_comp && found_this_item;
                    }
                    found_2 = is_this_comp || found_2;
                }
                
                REQUIRE(found_1);
                REQUIRE(found_2);
            }
            
            SECTION("Works for a multi-component multipath_alignment_t with complicated edge structure") {
                
                multipath_alignment_t multipath_aln;
                
                // add subpaths
                multipath_aln.add_subpath();
                multipath_aln.add_subpath();
                multipath_aln.add_subpath();
                multipath_aln.add_subpath();
                multipath_aln.add_subpath();
                multipath_aln.add_subpath();
                multipath_aln.add_subpath();
                subpath_t* subpath0 = multipath_aln.mutable_subpath(0);
                subpath_t* subpath1 = multipath_aln.mutable_subpath(1);
                subpath_t* subpath2 = multipath_aln.mutable_subpath(2);
                subpath_t* subpath3 = multipath_aln.mutable_subpath(3);
                subpath_t* subpath4 = multipath_aln.mutable_subpath(4);
                subpath_t* subpath5 = multipath_aln.mutable_subpath(5);
                subpath_t* subpath6 = multipath_aln.mutable_subpath(6);
                
                // set edges between subpaths
                subpath0->add_next(1);
                subpath1->add_next(5);
                subpath2->add_next(4);
                subpath3->add_next(5);
                
                // set scores
                subpath0->set_score(3);
                subpath1->set_score(0);
                subpath2->set_score(1);
                subpath3->set_score(0);
                subpath4->set_score(4);
                subpath5->set_score(0);
                subpath5->set_score(9);
                
                auto comps = connected_components(multipath_aln);
                
                REQUIRE(comps.size() == 3);
                
                bool found_1 = false;
                bool found_2 = false;
                bool found_3 = false;
                for (size_t i = 0; i < comps.size(); i++) {
                    bool is_this_comp;
                    is_this_comp = true;
                    for (int j : {0, 1, 3, 5}) {
                        bool found_this_item = false;
                        for (size_t k = 0; k < comps[i].size(); k++) {
                            if (comps[i][k] == j) {
                                found_this_item = true;
                                break;
                            }
                        }
                        is_this_comp = is_this_comp && found_this_item;
                    }
                    found_1 = is_this_comp || found_1;
                    
                    is_this_comp = true;
                    for (int j : {2, 4}) {
                        bool found_this_item = false;
                        for (size_t k = 0; k < comps[i].size(); k++) {
                            if (comps[i][k] == j) {
                                found_this_item = true;
                                break;
                            }
                        }
                        is_this_comp = is_this_comp && found_this_item;
                    }
                    found_2 = is_this_comp || found_2;
                    
                    is_this_comp = true;
                    for (int j : {6}) {
                        bool found_this_item = false;
                        for (size_t k = 0; k < comps[i].size(); k++) {
                            if (comps[i][k] == j) {
                                found_this_item = true;
                                break;
                            }
                        }
                        is_this_comp = is_this_comp && found_this_item;
                    }
                    found_3 = is_this_comp || found_3;
                }
                
                REQUIRE(found_1);
                REQUIRE(found_2);
                REQUIRE(found_3);
            }
        }
        
        TEST_CASE("We can extract subgraphs of a multipath_alignment_t",
                  "[alignment][multipath]" ) {
            
            SECTION("Works correctly when splitting by connected components") {
                
                multipath_alignment_t multipath_aln;
                
                // add subpaths
                multipath_aln.add_subpath();
                multipath_aln.add_subpath();
                multipath_aln.add_subpath();
                multipath_aln.add_subpath();
                multipath_aln.add_subpath();
                multipath_aln.add_subpath();
                multipath_aln.add_subpath();
                subpath_t* subpath0 = multipath_aln.mutable_subpath(0);
                subpath_t* subpath1 = multipath_aln.mutable_subpath(1);
                subpath_t* subpath2 = multipath_aln.mutable_subpath(2);
                subpath_t* subpath3 = multipath_aln.mutable_subpath(3);
                subpath_t* subpath4 = multipath_aln.mutable_subpath(4);
                subpath_t* subpath5 = multipath_aln.mutable_subpath(5);
                subpath_t* subpath6 = multipath_aln.mutable_subpath(6);
                
                // set edges between subpaths
                subpath0->add_next(1);
                subpath1->add_next(5);
                subpath2->add_next(4);
                subpath3->add_next(5);
                
                // set scores (hijacking them here to label subpaths with their original index)
                subpath0->set_score(0);
                subpath1->set_score(1);
                subpath2->set_score(2);
                subpath3->set_score(3);
                subpath4->set_score(4);
                subpath5->set_score(5);
                subpath5->set_score(6);
                
                auto comps = connected_components(multipath_aln);
                
                REQUIRE(comps.size() == 3);
                
                for (size_t i = 0; i < comps.size(); i++) {
                    
                    multipath_alignment_t sub;
                    extract_sub_multipath_alignment(multipath_aln, comps[i], sub);
                                        
                    bool is_this_comp;
                    bool finds_all_subpaths;
                    
                    is_this_comp = true;
                    for (int j : {0, 1, 3, 5}) {
                        bool found_this_item = false;
                        for (size_t k = 0; k < comps[i].size(); k++) {
                            if (comps[i][k] == j) {
                                found_this_item = true;
                                break;
                            }
                        }
                        is_this_comp = is_this_comp && found_this_item;
                    }
                    
                    if (is_this_comp) {
                        
                        REQUIRE(sub.subpath_size() == 4);
                        
                        finds_all_subpaths = true;
                        for (int j : {0, 1, 3, 5}) {
                            bool found_this_item = false;
                            for (size_t k = 0; k < sub.subpath_size(); k++) {
                                if (sub.subpath(k).score() == j) {
                                    REQUIRE(sub.subpath(k).next_size() == multipath_aln.subpath(j).next_size());
                                    for (size_t l = 0; l < sub.subpath(k).next_size(); l++) {
                                        bool found_this_next = false;
                                        for (size_t m = 0; m < multipath_aln.subpath(j).next_size(); m++) {
                                            if (sub.subpath(sub.subpath(k).next(l)).score() == multipath_aln.subpath(j).next(m)) {
                                                found_this_next = true;
                                            }
                                        }
                                        found_this_item = found_this_item && found_this_next;
                                    }
                                }
                            }
                            finds_all_subpaths = finds_all_subpaths && found_this_item;
                        }
                    }
                    
                    is_this_comp = true;
                    for (int j : {2, 4}) {
                        bool found_this_item = false;
                        for (size_t k = 0; k < comps[i].size(); k++) {
                            if (comps[i][k] == j) {
                                found_this_item = true;
                                break;
                            }
                        }
                        is_this_comp = is_this_comp && found_this_item;
                    }
                    
                    if (is_this_comp) {
                        
                        REQUIRE(sub.subpath_size() == 2);
                        
                        finds_all_subpaths = true;
                        for (int j : {2, 4}) {
                            bool found_this_item = false;
                            for (size_t k = 0; k < sub.subpath_size(); k++) {
                                if (sub.subpath(k).score() == j) {
                                    REQUIRE(sub.subpath(k).next_size() == multipath_aln.subpath(j).next_size());
                                    for (size_t l = 0; l < sub.subpath(k).next_size(); l++) {
                                        bool found_this_next = false;
                                        for (size_t m = 0; m < multipath_aln.subpath(j).next_size(); m++) {
                                            if (sub.subpath(sub.subpath(k).next(l)).score() == multipath_aln.subpath(j).next(m)) {
                                                found_this_next = true;
                                            }
                                        }
                                        found_this_item = found_this_item && found_this_next;
                                    }
                                }
                            }
                            finds_all_subpaths = finds_all_subpaths && found_this_item;
                        }
                    }
                    
                    is_this_comp = true;
                    for (int j : {6}) {
                        bool found_this_item = false;
                        for (size_t k = 0; k < comps[i].size(); k++) {
                            if (comps[i][k] == j) {
                                found_this_item = true;
                                break;
                            }
                        }
                        is_this_comp = is_this_comp && found_this_item;
                    }
                    
                    if (is_this_comp) {
                        
                        REQUIRE(sub.subpath_size() == 1);
                        
                        finds_all_subpaths = true;
                        for (int j : {6}) {
                            bool found_this_item = false;
                            for (size_t k = 0; k < sub.subpath_size(); k++) {
                                if (sub.subpath(k).score() == j) {
                                    REQUIRE(sub.subpath(k).next_size() == multipath_aln.subpath(j).next_size());
                                    for (size_t l = 0; l < sub.subpath(k).next_size(); l++) {
                                        bool found_this_next = false;
                                        for (size_t m = 0; m < multipath_aln.subpath(j).next_size(); m++) {
                                            if (sub.subpath(sub.subpath(k).next(l)).score() == multipath_aln.subpath(j).next(m)) {
                                                found_this_next = true;
                                            }
                                        }
                                        found_this_item = found_this_item && found_this_next;
                                    }
                                }
                            }
                            finds_all_subpaths = finds_all_subpaths && found_this_item;
                        }
                    }
                }
            }
        }
        
        TEST_CASE( "Non-branching paths in a multipath alignment can be merged", "[alignment][multipath]") {
            SECTION("Non-branching paths can be merged across an edge") {
                
                multipath_alignment_t mpaln;
                
                mpaln.add_subpath();
                mpaln.add_subpath();
                subpath_t* sp1 = mpaln.mutable_subpath(0);
                subpath_t* sp2 = mpaln.mutable_subpath(1);
                
                path_mapping_t* m11 = sp1->mutable_path()->add_mapping();
                position_t* p11 = m11->mutable_position();
                p11->set_node_id(1);
                
                edit_t* e111 = m11->add_edit();
                e111->set_from_length(1);
                e111->set_to_length(1);
                                
                path_mapping_t* m21 = sp2->mutable_path()->add_mapping();
                position_t* p21 = m21->mutable_position();
                p21->set_node_id(2);
                
                edit_t* e211 = m21->add_edit();
                e211->set_from_length(2);
                e211->set_to_length(2);
                
                sp1->add_next(1);
                
                merge_non_branching_subpaths(mpaln);
                
                REQUIRE(mpaln.subpath_size() == 1);
                
                REQUIRE(mpaln.subpath(0).path().mapping_size() == 2);
                
                REQUIRE(mpaln.subpath(0).path().mapping(0).position().node_id() == 1);
                REQUIRE(mpaln.subpath(0).path().mapping(0).position().is_reverse() == false);
                REQUIRE(mpaln.subpath(0).path().mapping(0).position().offset() == 0);
                
                REQUIRE(mpaln.subpath(0).path().mapping(0).edit_size() == 1);
                REQUIRE(mpaln.subpath(0).path().mapping(0).edit(0).from_length() == 1);
                REQUIRE(mpaln.subpath(0).path().mapping(0).edit(0).to_length() == 1);
                REQUIRE(mpaln.subpath(0).path().mapping(0).edit(0).sequence() == "");
                
                REQUIRE(mpaln.subpath(0).path().mapping(1).position().node_id() == 2);
                REQUIRE(mpaln.subpath(0).path().mapping(1).position().is_reverse() == false);
                REQUIRE(mpaln.subpath(0).path().mapping(1).position().offset() == 0);
                
                REQUIRE(mpaln.subpath(0).path().mapping(1).edit_size() == 1);
                REQUIRE(mpaln.subpath(0).path().mapping(1).edit(0).from_length() == 2);
                REQUIRE(mpaln.subpath(0).path().mapping(1).edit(0).to_length() == 2);
                REQUIRE(mpaln.subpath(0).path().mapping(1).edit(0).sequence() == "");
                
                REQUIRE(mpaln.subpath(0).next_size() == 0);
            }
            
            SECTION("Non-branching paths can be merged within a node") {
                
                multipath_alignment_t mpaln;
                
                mpaln.add_subpath();
                mpaln.add_subpath();
                subpath_t* sp1 = mpaln.mutable_subpath(0);
                subpath_t* sp2 = mpaln.mutable_subpath(1);
                
                path_mapping_t* m11 = sp1->mutable_path()->add_mapping();
                position_t* p11 = m11->mutable_position();
                p11->set_node_id(1);
                p11->set_is_reverse(true);
                
                edit_t* e111 = m11->add_edit();
                e111->set_from_length(1);
                e111->set_to_length(1);
                                
                path_mapping_t* m21 = sp2->mutable_path()->add_mapping();
                position_t* p21 = m21->mutable_position();
                p21->set_node_id(1);
                p21->set_offset(1);
                p21->set_is_reverse(true);
                
                edit_t* e211 = m21->add_edit();
                e211->set_from_length(0);
                e211->set_to_length(2);
                
                sp1->add_next(1);
                
                merge_non_branching_subpaths(mpaln);
                
                REQUIRE(mpaln.subpath_size() == 1);
                
                REQUIRE(mpaln.subpath(0).path().mapping_size() == 1);
                
                REQUIRE(mpaln.subpath(0).path().mapping(0).position().node_id() == 1);
                REQUIRE(mpaln.subpath(0).path().mapping(0).position().is_reverse() == true);
                REQUIRE(mpaln.subpath(0).path().mapping(0).position().offset() == 0);
                
                REQUIRE(mpaln.subpath(0).path().mapping(0).edit_size() == 2);
                
                REQUIRE(mpaln.subpath(0).path().mapping(0).edit(0).from_length() == 1);
                REQUIRE(mpaln.subpath(0).path().mapping(0).edit(0).to_length() == 1);
                REQUIRE(mpaln.subpath(0).path().mapping(0).edit(0).sequence() == "");
                
                REQUIRE(mpaln.subpath(0).path().mapping(0).edit(1).from_length() == 0);
                REQUIRE(mpaln.subpath(0).path().mapping(0).edit(1).to_length() == 2);
                REQUIRE(mpaln.subpath(0).path().mapping(0).edit(1).sequence() == "");
                
                REQUIRE(mpaln.subpath(0).next_size() == 0);
            }
            
            SECTION("Non-branching paths can be merged within an edit") {
                
                multipath_alignment_t mpaln;
                
                mpaln.add_subpath();
                mpaln.add_subpath();
                subpath_t* sp1 = mpaln.mutable_subpath(0);
                subpath_t* sp2 = mpaln.mutable_subpath(1);
                
                path_mapping_t* m11 = sp1->mutable_path()->add_mapping();
                position_t* p11 = m11->mutable_position();
                p11->set_node_id(1);
                p11->set_is_reverse(true);
                
                edit_t* e111 = m11->add_edit();
                e111->set_from_length(1);
                e111->set_to_length(1);
                                
                path_mapping_t* m21 = sp2->mutable_path()->add_mapping();
                position_t* p21 = m21->mutable_position();
                p21->set_node_id(1);
                p21->set_offset(1);
                p21->set_is_reverse(true);
                
                edit_t* e211 = m21->add_edit();
                e211->set_from_length(2);
                e211->set_to_length(2);
                
                sp1->add_next(1);
                
                merge_non_branching_subpaths(mpaln);
                
                REQUIRE(mpaln.subpath_size() == 1);
                
                REQUIRE(mpaln.subpath(0).path().mapping_size() == 1);
                
                REQUIRE(mpaln.subpath(0).path().mapping(0).position().node_id() == 1);
                REQUIRE(mpaln.subpath(0).path().mapping(0).position().is_reverse() == true);
                REQUIRE(mpaln.subpath(0).path().mapping(0).position().offset() == 0);
                
                REQUIRE(mpaln.subpath(0).path().mapping(0).edit_size() == 1);
                REQUIRE(mpaln.subpath(0).path().mapping(0).edit(0).from_length() == 3);
                REQUIRE(mpaln.subpath(0).path().mapping(0).edit(0).to_length() == 3);
                REQUIRE(mpaln.subpath(0).path().mapping(0).edit(0).sequence() == "");
                
                REQUIRE(mpaln.subpath(0).next_size() == 0);
            }
            
            SECTION("Non-branching paths can be distinguished from branching paths") {
                
                multipath_alignment_t mpaln;
                
                mpaln.add_subpath();
                mpaln.add_subpath();
                mpaln.add_subpath();
                mpaln.add_subpath();
                subpath_t* sp1 = mpaln.mutable_subpath(0);
                subpath_t* sp2 = mpaln.mutable_subpath(1);
                subpath_t* sp3 = mpaln.mutable_subpath(2);
                subpath_t* sp4 = mpaln.mutable_subpath(3);
                
                path_mapping_t* m11 = sp1->mutable_path()->add_mapping();
                position_t* p11 = m11->mutable_position();
                p11->set_node_id(1);
                p11->set_is_reverse(true);
                
                edit_t* e111 = m11->add_edit();
                e111->set_from_length(1);
                e111->set_to_length(1);
                                
                path_mapping_t* m21 = sp2->mutable_path()->add_mapping();
                position_t* p21 = m21->mutable_position();
                p21->set_node_id(1);
                p21->set_offset(1);
                p21->set_is_reverse(true);
                
                m21->add_edit();
                m21->add_edit();
                edit_t* e211 = m21->mutable_edit(0);
                edit_t* e212 = m21->mutable_edit(1);
                
                e211->set_from_length(2);
                e211->set_to_length(2);
                
                e212->set_from_length(2);
                e212->set_to_length(0);
                                
                path_mapping_t* m31 = sp3->mutable_path()->add_mapping();
                position_t* p31 = m31->mutable_position();
                p31->set_node_id(2);
                p31->set_offset(0);
                
                edit_t* e311 = m31->add_edit();
                e311->set_from_length(1);
                e311->set_to_length(1);
                                
                path_mapping_t* m41 = sp4->mutable_path()->add_mapping();
                position_t* p41 = m41->mutable_position();
                p41->set_node_id(3);
                p41->set_offset(0);
                
                edit_t* e411 = m41->add_edit();
                e411->set_from_length(1);
                e411->set_to_length(1);
                
                sp1->add_next(1);
                sp2->add_next(2);
                sp2->add_next(3);
                
                merge_non_branching_subpaths(mpaln);
                                
                REQUIRE(mpaln.subpath_size() == 3);
                
                REQUIRE(mpaln.subpath(0).path().mapping_size() == 1);
                
                REQUIRE(mpaln.subpath(0).path().mapping(0).position().node_id() == 1);
                REQUIRE(mpaln.subpath(0).path().mapping(0).position().is_reverse() == true);
                REQUIRE(mpaln.subpath(0).path().mapping(0).position().offset() == 0);
                
                REQUIRE(mpaln.subpath(0).path().mapping(0).edit_size() == 2);
                
                REQUIRE(mpaln.subpath(0).path().mapping(0).edit(0).from_length() == 3);
                REQUIRE(mpaln.subpath(0).path().mapping(0).edit(0).to_length() == 3);
                REQUIRE(mpaln.subpath(0).path().mapping(0).edit(0).sequence() == "");
                
                REQUIRE(mpaln.subpath(0).path().mapping(0).edit(1).from_length() == 2);
                REQUIRE(mpaln.subpath(0).path().mapping(0).edit(1).to_length() == 0);
                REQUIRE(mpaln.subpath(0).path().mapping(0).edit(1).sequence() == "");
                
                REQUIRE(mpaln.subpath(0).next_size() == 2);
                
                REQUIRE(mpaln.subpath(0).next(0) == 1);
                REQUIRE(mpaln.subpath(0).next(1) == 2);
                
                REQUIRE(mpaln.subpath(1).path().mapping_size() == 1);
                
                REQUIRE(mpaln.subpath(1).path().mapping(0).position().node_id() == 2);
                REQUIRE(mpaln.subpath(1).path().mapping(0).position().is_reverse() == false);
                REQUIRE(mpaln.subpath(1).path().mapping(0).position().offset() == 0);
                
                REQUIRE(mpaln.subpath(1).path().mapping(0).edit_size() == 1);
                
                REQUIRE(mpaln.subpath(1).path().mapping(0).edit(0).from_length() == 1);
                REQUIRE(mpaln.subpath(1).path().mapping(0).edit(0).to_length() == 1);
                REQUIRE(mpaln.subpath(1).path().mapping(0).edit(0).sequence() == "");
                
                REQUIRE(mpaln.subpath(1).next_size() == 0);
                
                REQUIRE(mpaln.subpath(2).path().mapping_size() == 1);
                
                REQUIRE(mpaln.subpath(2).path().mapping(0).position().node_id() == 3);
                REQUIRE(mpaln.subpath(2).path().mapping(0).position().is_reverse() == false);
                REQUIRE(mpaln.subpath(2).path().mapping(0).position().offset() == 0);
                
                REQUIRE(mpaln.subpath(2).path().mapping(0).edit_size() == 1);
                
                REQUIRE(mpaln.subpath(2).path().mapping(0).edit(0).from_length() == 1);
                REQUIRE(mpaln.subpath(2).path().mapping(0).edit(0).to_length() == 1);
                REQUIRE(mpaln.subpath(2).path().mapping(0).edit(0).sequence() == "");
                
                REQUIRE(mpaln.subpath(2).next_size() == 0);
            }
        }
        
        TEST_CASE( "Single path alignments with disjoint subpaths can be found", "[alignment][multipath]") {
        
            string multipath_json = R"(
            {"sequence":"TGAGGTGGTCTCAGCTGGAGGTGAGGAACTTGTTGGGAACTAGAGTGAAAGTCATTCTTGCTATGTTAGGAAAAACTGGGTGCATTTTGCTCCTGCCCTAGAGATCCGTGGAACACTGAACTTGAGAAAGATGATTTACCGTATCTGG","quality":"ISEjIyMjIyMiIh0iHh0fIiIjIyMjIyMhIiAkJCMjIyMjIyMjHyIgIyMjIyMiHx0dIiIjIyMjIyMjIyIiICIhIiApKCkpKCcnJSQlJSkoKSkoKCYnJCMiIyMpKCgoKCgnJyQlISclKSkpKSkpKCgmKSgpKCkpKSkpKSkpKCkoKSkpKSkpKSkoJyUnJiclJSUlJR8fHw==","name":"seed_90_fragment_287_1","subpath":[{"path":{"mapping":[{"position":{"node_id":"1724335","offset":"12","is_reverse":true},"edit":[{"from_length":14,"to_length":14}],"rank":"1"}]},"next":[2,1],"score":19},{"path":{"mapping":[{"position":{"node_id":"1724335","offset":"26","is_reverse":true},"edit":[{"from_length":1},{"to_length":1,"sequence":"C"}],"rank":"1"}]},"next":[3],"score":-12},{"path":{"mapping":[{"position":{"node_id":"1724335","offset":"26","is_reverse":true},"edit":[{"from_length":1,"to_length":1,"sequence":"C"}],"rank":"1"}]},"next":[3],"score":-4},{"path":{"mapping":[{"position":{"node_id":"1724335","offset":"27","is_reverse":true},"edit":[{"from_length":5,"to_length":5}],"rank":"1"}]},"next":[7,6,5,4],"score":5},{"path":{"mapping":[{"position":{"node_id":"1724333","is_reverse":true},"edit":[{"from_length":1},{"to_length":1,"sequence":"G"}],"rank":"1"}]},"next":[8],"score":-12},{"path":{"mapping":[{"position":{"node_id":"1724334","is_reverse":true},"edit":[{"from_length":1},{"to_length":1,"sequence":"G"}],"rank":"1"}]},"next":[8],"score":-12},{"path":{"mapping":[{"position":{"node_id":"1724333","is_reverse":true},"edit":[{"from_length":1,"to_length":1,"sequence":"G"}],"rank":"1"}]},"next":[8],"score":-4},{"path":{"mapping":[{"position":{"node_id":"1724334","is_reverse":true},"edit":[{"from_length":1,"to_length":1}],"rank":"1"}]},"next":[8],"score":1},{"path":{"mapping":[{"position":{"node_id":"1724332","is_reverse":true},"edit":[{"from_length":4,"to_length":4}],"rank":"1"}]},"next":[12,11,10,9],"score":4},{"path":{"mapping":[{"position":{"node_id":"1724331","is_reverse":true},"edit":[{"from_length":1},{"to_length":1,"sequence":"G"}],"rank":"1"}]},"next":[13],"score":-12},{"path":{"mapping":[{"position":{"node_id":"1724330","is_reverse":true},"edit":[{"from_length":1},{"to_length":1,"sequence":"G"}],"rank":"1"}]},"next":[13],"score":-12},{"path":{"mapping":[{"position":{"node_id":"1724330","is_reverse":true},"edit":[{"from_length":1,"to_length":1,"sequence":"G"}],"rank":"1"}]},"next":[13],"score":-4},{"path":{"mapping":[{"position":{"node_id":"1724331","is_reverse":true},"edit":[{"from_length":1,"to_length":1}],"rank":"1"}]},"next":[13],"score":1},{"path":{"mapping":[{"position":{"node_id":"1724329","is_reverse":true},"edit":[{"from_length":15,"to_length":15}],"rank":"1"}]},"next":[17,16,15,14],"score":15},{"path":{"mapping":[{"position":{"node_id":"1724327","is_reverse":true},"edit":[{"from_length":1},{"to_length":1,"sequence":"A"}],"rank":"1"}]},"next":[18],"score":-12},{"path":{"mapping":[{"position":{"node_id":"1724328","is_reverse":true},"edit":[{"from_length":1},{"to_length":1,"sequence":"A"}],"rank":"1"}]},"next":[18],"score":-12},{"path":{"mapping":[{"position":{"node_id":"1724327","is_reverse":true},"edit":[{"from_length":1,"to_length":1,"sequence":"A"}],"rank":"1"}]},"next":[18],"score":-4},{"path":{"mapping":[{"position":{"node_id":"1724328","is_reverse":true},"edit":[{"from_length":1,"to_length":1}],"rank":"1"}]},"next":[18],"score":1},{"path":{"mapping":[{"position":{"node_id":"1724326","is_reverse":true},"edit":[{"from_length":21,"to_length":21}],"rank":"1"}]},"next":[22,21,20,19],"score":21},{"path":{"mapping":[{"position":{"node_id":"1724325","is_reverse":true},"edit":[{"from_length":1},{"to_length":1,"sequence":"T"}],"rank":"1"}]},"next":[23],"score":-12},{"path":{"mapping":[{"position":{"node_id":"1724324","is_reverse":true},"edit":[{"from_length":1},{"to_length":1,"sequence":"T"}],"rank":"1"}]},"next":[23],"score":-12},{"path":{"mapping":[{"position":{"node_id":"1724324","is_reverse":true},"edit":[{"from_length":1,"to_length":1,"sequence":"T"}],"rank":"1"}]},"next":[23],"score":-4},{"path":{"mapping":[{"position":{"node_id":"1724325","is_reverse":true},"edit":[{"from_length":1,"to_length":1}],"rank":"1"}]},"next":[23],"score":1},{"path":{"mapping":[{"position":{"node_id":"1724323","is_reverse":true},"edit":[{"from_length":1,"to_length":1}],"rank":"1"}]},"next":[32,24],"score":1},{"path":{"mapping":[{"position":{"node_id":"1724323","offset":"1","is_reverse":true},"edit":[{"from_length":2}],"rank":"1"}]},"next":[25],"score":-7},{"path":{"mapping":[{"position":{"node_id":"1724323","offset":"3","is_reverse":true},"edit":[{"from_length":1,"to_length":1}],"rank":"1"}]},"next":[31,30,29,28,27,26],"score":1},{"path":{"mapping":[{"position":{"node_id":"1724322","is_reverse":true},"edit":[{"from_length":1},{"to_length":1,"sequence":"T"}],"rank":"1"}]},"next":[35],"score":-12},{"path":{"mapping":[{"position":{"node_id":"1724323","offset":"4","is_reverse":true},"edit":[{"to_length":1,"sequence":"T"}],"rank":"1"}]},"next":[35],"score":-6},{"path":{"mapping":[{"position":{"node_id":"1724322","is_reverse":true},"edit":[{"from_length":1,"to_length":1}],"rank":"1"}]},"next":[35],"score":1},{"path":{"mapping":[{"position":{"node_id":"1724322","is_reverse":true},"edit":[{"from_length":1},{"to_length":1,"sequence":"T"}],"rank":"1"}]},"next":[45],"score":-12},{"path":{"mapping":[{"position":{"node_id":"1724323","offset":"4","is_reverse":true},"edit":[{"to_length":1,"sequence":"T"}],"rank":"1"}]},"next":[45],"score":-6},{"path":{"mapping":[{"position":{"node_id":"1724322","is_reverse":true},"edit":[{"from_length":1,"to_length":1}],"rank":"1"}]},"next":[45],"score":1},{"path":{"mapping":[{"position":{"node_id":"1724323","offset":"1","is_reverse":true},"edit":[{"from_length":1}],"rank":"1"}]},"next":[33],"score":-6},{"path":{"mapping":[{"position":{"node_id":"1724323","offset":"2","is_reverse":true},"edit":[{"from_length":2,"to_length":2}],"rank":"1"}]},"next":[45,44,35,34],"score":2},{"path":{"mapping":[{"position":{"node_id":"1724322","is_reverse":true},"edit":[{"from_length":1}],"rank":"1"}]},"next":[35],"score":-6},{"path":{"mapping":[{"position":{"node_id":"1724321","is_reverse":true},"edit":[{"from_length":2,"to_length":2}],"rank":"1"}]},"next":[43,42,41,40,39,38,37,36],"score":2},{"path":{"mapping":[{"position":{"node_id":"1724319","is_reverse":true},"edit":[{"from_length":1},{"to_length":1,"sequence":"G"}],"rank":"1"}]},"next":[50],"score":-12},{"path":{"mapping":[{"position":{"node_id":"1724320","is_reverse":true},"edit":[{"from_length":1},{"to_length":1,"sequence":"G"}],"rank":"1"}]},"next":[50],"score":-12},{"path":{"mapping":[{"position":{"node_id":"1724319","is_reverse":true},"edit":[{"from_length":1,"to_length":1,"sequence":"G"}],"rank":"1"}]},"next":[50],"score":-4},{"path":{"mapping":[{"position":{"node_id":"1724320","is_reverse":true},"edit":[{"from_length":1,"to_length":1}],"rank":"1"}]},"next":[50],"score":1},{"path":{"mapping":[{"position":{"node_id":"1724319","is_reverse":true},"edit":[{"from_length":1},{"to_length":1,"sequence":"G"}],"rank":"1"}]},"next":[63],"score":-12},{"path":{"mapping":[{"position":{"node_id":"1724320","is_reverse":true},"edit":[{"from_length":1},{"to_length":1,"sequence":"G"}],"rank":"1"}]},"next":[63],"score":-12},{"path":{"mapping":[{"position":{"node_id":"1724319","is_reverse":true},"edit":[{"from_length":1,"to_length":1,"sequence":"G"}],"rank":"1"}]},"next":[63],"score":-4},{"path":{"mapping":[{"position":{"node_id":"1724320","is_reverse":true},"edit":[{"from_length":1,"to_length":1}],"rank":"1"}]},"next":[63],"score":1},{"path":{"mapping":[{"position":{"node_id":"1724322","is_reverse":true},"edit":[{"from_length":1}],"rank":"1"}]},"next":[45],"score":-6},{"path":{"mapping":[{"position":{"node_id":"1724321","is_reverse":true},"edit":[{"from_length":2,"to_length":2}],"rank":"1"}]},"next":[62,61,60,59,49,48,47,46],"score":2},{"path":{"mapping":[{"position":{"node_id":"1724319","is_reverse":true},"edit":[{"from_length":1},{"to_length":1,"sequence":"G"}],"rank":"1"}]},"next":[50],"score":-12},{"path":{"mapping":[{"position":{"node_id":"1724320","is_reverse":true},"edit":[{"from_length":1},{"to_length":1,"sequence":"G"}],"rank":"1"}]},"next":[50],"score":-12},{"path":{"mapping":[{"position":{"node_id":"1724319","is_reverse":true},"edit":[{"from_length":1,"to_length":1,"sequence":"G"}],"rank":"1"}]},"next":[50],"score":-4},{"path":{"mapping":[{"position":{"node_id":"1724320","is_reverse":true},"edit":[{"from_length":1,"to_length":1}],"rank":"1"}]},"next":[50],"score":1},{"path":{"mapping":[{"position":{"node_id":"1724318","is_reverse":true},"edit":[{"from_length":8,"to_length":8}],"rank":"1"}]},"next":[58,57,56,55,54,53,52,51],"score":8},{"path":{"mapping":[{"position":{"node_id":"1724317","is_reverse":true},"edit":[{"from_length":1},{"to_length":1,"sequence":"G"}],"rank":"1"}]},"next":[68],"score":-12},{"path":{"mapping":[{"position":{"node_id":"1724316","is_reverse":true},"edit":[{"from_length":1},{"to_length":1,"sequence":"G"}],"rank":"1"}]},"next":[68],"score":-12},{"path":{"mapping":[{"position":{"node_id":"1724316","is_reverse":true},"edit":[{"from_length":1,"to_length":1,"sequence":"G"}],"rank":"1"}]},"next":[68],"score":-4},{"path":{"mapping":[{"position":{"node_id":"1724317","is_reverse":true},"edit":[{"from_length":1,"to_length":1}],"rank":"1"}]},"next":[68],"score":1},{"path":{"mapping":[{"position":{"node_id":"1724317","is_reverse":true},"edit":[{"from_length":1},{"to_length":1,"sequence":"G"}],"rank":"1"}]},"next":[81],"score":-12},{"path":{"mapping":[{"position":{"node_id":"1724316","is_reverse":true},"edit":[{"from_length":1},{"to_length":1,"sequence":"G"}],"rank":"1"}]},"next":[81],"score":-12},{"path":{"mapping":[{"position":{"node_id":"1724316","is_reverse":true},"edit":[{"from_length":1,"to_length":1,"sequence":"G"}],"rank":"1"}]},"next":[81],"score":-4},{"path":{"mapping":[{"position":{"node_id":"1724317","is_reverse":true},"edit":[{"from_length":1,"to_length":1}],"rank":"1"}]},"next":[81],"score":1},{"path":{"mapping":[{"position":{"node_id":"1724319","is_reverse":true},"edit":[{"from_length":1},{"to_length":1,"sequence":"G"}],"rank":"1"}]},"next":[63],"score":-12},{"path":{"mapping":[{"position":{"node_id":"1724320","is_reverse":true},"edit":[{"from_length":1},{"to_length":1,"sequence":"G"}],"rank":"1"}]},"next":[63],"score":-12},{"path":{"mapping":[{"position":{"node_id":"1724319","is_reverse":true},"edit":[{"from_length":1,"to_length":1,"sequence":"G"}],"rank":"1"}]},"next":[63],"score":-4},{"path":{"mapping":[{"position":{"node_id":"1724320","is_reverse":true},"edit":[{"from_length":1,"to_length":1}],"rank":"1"}]},"next":[63],"score":1},{"path":{"mapping":[{"position":{"node_id":"1724318","is_reverse":true},"edit":[{"from_length":8,"to_length":8}],"rank":"1"}]},"next":[80,79,78,77,67,66,65,64],"score":8},{"path":{"mapping":[{"position":{"node_id":"1724317","is_reverse":true},"edit":[{"from_length":1},{"to_length":1,"sequence":"G"}],"rank":"1"}]},"next":[68],"score":-12},{"path":{"mapping":[{"position":{"node_id":"1724316","is_reverse":true},"edit":[{"from_length":1},{"to_length":1,"sequence":"G"}],"rank":"1"}]},"next":[68],"score":-12},{"path":{"mapping":[{"position":{"node_id":"1724316","is_reverse":true},"edit":[{"from_length":1,"to_length":1,"sequence":"G"}],"rank":"1"}]},"next":[68],"score":-4},{"path":{"mapping":[{"position":{"node_id":"1724317","is_reverse":true},"edit":[{"from_length":1,"to_length":1}],"rank":"1"}]},"next":[68],"score":1},{"path":{"mapping":[{"position":{"node_id":"1724315","is_reverse":true},"edit":[{"from_length":4,"to_length":4}],"rank":"1"}]},"next":[76,75,74,73,72,71,70,69],"score":4},{"path":{"mapping":[{"position":{"node_id":"1724313","is_reverse":true},"edit":[{"from_length":1},{"to_length":1,"sequence":"A"}],"rank":"1"}]},"next":[86],"score":-12},{"path":{"mapping":[{"position":{"node_id":"1724314","is_reverse":true},"edit":[{"from_length":1},{"to_length":1,"sequence":"A"}],"rank":"1"}]},"next":[86],"score":-12},{"path":{"mapping":[{"position":{"node_id":"1724313","is_reverse":true},"edit":[{"from_length":1,"to_length":1,"sequence":"A"}],"rank":"1"}]},"next":[86],"score":-4},{"path":{"mapping":[{"position":{"node_id":"1724314","is_reverse":true},"edit":[{"from_length":1,"to_length":1}],"rank":"1"}]},"next":[86],"score":1},{"path":{"mapping":[{"position":{"node_id":"1724313","is_reverse":true},"edit":[{"from_length":1},{"to_length":1,"sequence":"A"}],"rank":"1"}]},"next":[99],"score":-12},{"path":{"mapping":[{"position":{"node_id":"1724314","is_reverse":true},"edit":[{"from_length":1},{"to_length":1,"sequence":"A"}],"rank":"1"}]},"next":[99],"score":-12},{"path":{"mapping":[{"position":{"node_id":"1724313","is_reverse":true},"edit":[{"from_length":1,"to_length":1,"sequence":"A"}],"rank":"1"}]},"next":[99],"score":-4},{"path":{"mapping":[{"position":{"node_id":"1724314","is_reverse":true},"edit":[{"from_length":1,"to_length":1}],"rank":"1"}]},"next":[99],"score":1},{"path":{"mapping":[{"position":{"node_id":"1724317","is_reverse":true},"edit":[{"from_length":1},{"to_length":1,"sequence":"G"}],"rank":"1"}]},"next":[81],"score":-12},{"path":{"mapping":[{"position":{"node_id":"1724316","is_reverse":true},"edit":[{"from_length":1},{"to_length":1,"sequence":"G"}],"rank":"1"}]},"next":[81],"score":-12},{"path":{"mapping":[{"position":{"node_id":"1724316","is_reverse":true},"edit":[{"from_length":1,"to_length":1,"sequence":"G"}],"rank":"1"}]},"next":[81],"score":-4},{"path":{"mapping":[{"position":{"node_id":"1724317","is_reverse":true},"edit":[{"from_length":1,"to_length":1}],"rank":"1"}]},"next":[81],"score":1},{"path":{"mapping":[{"position":{"node_id":"1724315","is_reverse":true},"edit":[{"from_length":4,"to_length":4}],"rank":"1"}]},"next":[98,97,96,95,85,84,83,82],"score":4},{"path":{"mapping":[{"position":{"node_id":"1724313","is_reverse":true},"edit":[{"from_length":1},{"to_length":1,"sequence":"A"}],"rank":"1"}]},"next":[86],"score":-12},{"path":{"mapping":[{"position":{"node_id":"1724314","is_reverse":true},"edit":[{"from_length":1},{"to_length":1,"sequence":"A"}],"rank":"1"}]},"next":[86],"score":-12},{"path":{"mapping":[{"position":{"node_id":"1724313","is_reverse":true},"edit":[{"from_length":1,"to_length":1,"sequence":"A"}],"rank":"1"}]},"next":[86],"score":-4},{"path":{"mapping":[{"position":{"node_id":"1724314","is_reverse":true},"edit":[{"from_length":1,"to_length":1}],"rank":"1"}]},"next":[86],"score":1},{"path":{"mapping":[{"position":{"node_id":"1724312","is_reverse":true},"edit":[{"from_length":6,"to_length":6}],"rank":"1"}]},"next":[94,93,92,91,90,89,88,87],"score":6},{"path":{"mapping":[{"position":{"node_id":"1724310","is_reverse":true},"edit":[{"from_length":1},{"to_length":1,"sequence":"T"}],"rank":"1"}]},"next":[104],"score":-12},{"path":{"mapping":[{"position":{"node_id":"1724311","is_reverse":true},"edit":[{"from_length":1},{"to_length":1,"sequence":"T"}],"rank":"1"}]},"next":[104],"score":-12},{"path":{"mapping":[{"position":{"node_id":"1724310","is_reverse":true},"edit":[{"from_length":1,"to_length":1,"sequence":"T"}],"rank":"1"}]},"next":[104],"score":-4},{"path":{"mapping":[{"position":{"node_id":"1724311","is_reverse":true},"edit":[{"from_length":1,"to_length":1}],"rank":"1"}]},"next":[104],"score":1},{"path":{"mapping":[{"position":{"node_id":"1724310","is_reverse":true},"edit":[{"from_length":1},{"to_length":1,"sequence":"T"}],"rank":"1"}]},"next":[117],"score":-12},{"path":{"mapping":[{"position":{"node_id":"1724311","is_reverse":true},"edit":[{"from_length":1},{"to_length":1,"sequence":"T"}],"rank":"1"}]},"next":[117],"score":-12},{"path":{"mapping":[{"position":{"node_id":"1724310","is_reverse":true},"edit":[{"from_length":1,"to_length":1,"sequence":"T"}],"rank":"1"}]},"next":[117],"score":-4},{"path":{"mapping":[{"position":{"node_id":"1724311","is_reverse":true},"edit":[{"from_length":1,"to_length":1}],"rank":"1"}]},"next":[117],"score":1},{"path":{"mapping":[{"position":{"node_id":"1724313","is_reverse":true},"edit":[{"from_length":1},{"to_length":1,"sequence":"A"}],"rank":"1"}]},"next":[99],"score":-12},{"path":{"mapping":[{"position":{"node_id":"1724314","is_reverse":true},"edit":[{"from_length":1},{"to_length":1,"sequence":"A"}],"rank":"1"}]},"next":[99],"score":-12},{"path":{"mapping":[{"position":{"node_id":"1724313","is_reverse":true},"edit":[{"from_length":1,"to_length":1,"sequence":"A"}],"rank":"1"}]},"next":[99],"score":-4},{"path":{"mapping":[{"position":{"node_id":"1724314","is_reverse":true},"edit":[{"from_length":1,"to_length":1}],"rank":"1"}]},"next":[99],"score":1},{"path":{"mapping":[{"position":{"node_id":"1724312","is_reverse":true},"edit":[{"from_length":6,"to_length":6}],"rank":"1"}]},"next":[116,115,114,113,103,102,101,100],"score":6},{"path":{"mapping":[{"position":{"node_id":"1724310","is_reverse":true},"edit":[{"from_length":1},{"to_length":1,"sequence":"T"}],"rank":"1"}]},"next":[104],"score":-12},{"path":{"mapping":[{"position":{"node_id":"1724311","is_reverse":true},"edit":[{"from_length":1},{"to_length":1,"sequence":"T"}],"rank":"1"}]},"next":[104],"score":-12},{"path":{"mapping":[{"position":{"node_id":"1724310","is_reverse":true},"edit":[{"from_length":1,"to_length":1,"sequence":"T"}],"rank":"1"}]},"next":[104],"score":-4},{"path":{"mapping":[{"position":{"node_id":"1724311","is_reverse":true},"edit":[{"from_length":1,"to_length":1}],"rank":"1"}]},"next":[104],"score":1},{"path":{"mapping":[{"position":{"node_id":"1724309","is_reverse":true},"edit":[{"from_length":15,"to_length":15}],"rank":"1"}]},"next":[112,111,110,109,108,107,106,105],"score":15},{"path":{"mapping":[{"position":{"node_id":"1724308","is_reverse":true},"edit":[{"from_length":1},{"to_length":1,"sequence":"C"}],"rank":"1"}]},"next":[122],"score":-12},{"path":{"mapping":[{"position":{"node_id":"1724307","is_reverse":true},"edit":[{"from_length":1},{"to_length":1,"sequence":"C"}],"rank":"1"}]},"next":[122],"score":-12},{"path":{"mapping":[{"position":{"node_id":"1724307","is_reverse":true},"edit":[{"from_length":1,"to_length":1,"sequence":"C"}],"rank":"1"}]},"next":[122],"score":-4},{"path":{"mapping":[{"position":{"node_id":"1724308","is_reverse":true},"edit":[{"from_length":1,"to_length":1}],"rank":"1"}]},"next":[122],"score":1},{"path":{"mapping":[{"position":{"node_id":"1724308","is_reverse":true},"edit":[{"from_length":1},{"to_length":1,"sequence":"C"}],"rank":"1"}]},"next":[133],"score":-12},{"path":{"mapping":[{"position":{"node_id":"1724307","is_reverse":true},"edit":[{"from_length":1},{"to_length":1,"sequence":"C"}],"rank":"1"}]},"next":[133],"score":-12},{"path":{"mapping":[{"position":{"node_id":"1724307","is_reverse":true},"edit":[{"from_length":1,"to_length":1,"sequence":"C"}],"rank":"1"}]},"next":[133],"score":-4},{"path":{"mapping":[{"position":{"node_id":"1724308","is_reverse":true},"edit":[{"from_length":1,"to_length":1}],"rank":"1"}]},"next":[133],"score":1},{"path":{"mapping":[{"position":{"node_id":"1724310","is_reverse":true},"edit":[{"from_length":1},{"to_length":1,"sequence":"T"}],"rank":"1"}]},"next":[117],"score":-12},{"path":{"mapping":[{"position":{"node_id":"1724311","is_reverse":true},"edit":[{"from_length":1},{"to_length":1,"sequence":"T"}],"rank":"1"}]},"next":[117],"score":-12},{"path":{"mapping":[{"position":{"node_id":"1724310","is_reverse":true},"edit":[{"from_length":1,"to_length":1,"sequence":"T"}],"rank":"1"}]},"next":[117],"score":-4},{"path":{"mapping":[{"position":{"node_id":"1724311","is_reverse":true},"edit":[{"from_length":1,"to_length":1}],"rank":"1"}]},"next":[117],"score":1},{"path":{"mapping":[{"position":{"node_id":"1724309","is_reverse":true},"edit":[{"from_length":15,"to_length":15}],"rank":"1"}]},"next":[132,131,130,129,121,120,119,118],"score":15},{"path":{"mapping":[{"position":{"node_id":"1724308","is_reverse":true},"edit":[{"from_length":1},{"to_length":1,"sequence":"C"}],"rank":"1"}]},"next":[122],"score":-12},{"path":{"mapping":[{"position":{"node_id":"1724307","is_reverse":true},"edit":[{"from_length":1},{"to_length":1,"sequence":"C"}],"rank":"1"}]},"next":[122],"score":-12},{"path":{"mapping":[{"position":{"node_id":"1724307","is_reverse":true},"edit":[{"from_length":1,"to_length":1,"sequence":"C"}],"rank":"1"}]},"next":[122],"score":-4},{"path":{"mapping":[{"position":{"node_id":"1724308","is_reverse":true},"edit":[{"from_length":1,"to_length":1}],"rank":"1"}]},"next":[122],"score":1},{"path":{"mapping":[{"position":{"node_id":"1724306","is_reverse":true},"edit":[{"from_length":26,"to_length":26}],"rank":"1"}]},"next":[128,127,126,125,124,123],"score":26},{"path":{"mapping":[{"position":{"node_id":"1724305","is_reverse":true},"edit":[{"from_length":1},{"to_length":1,"sequence":"A"}],"rank":"1"}]},"next":[137],"score":-12},{"path":{"mapping":[{"position":{"node_id":"1724306","offset":"26","is_reverse":true},"edit":[{"to_length":1,"sequence":"A"}],"rank":"1"}]},"next":[137],"score":-6},{"path":{"mapping":[{"position":{"node_id":"1724305","is_reverse":true},"edit":[{"from_length":1,"to_length":1}],"rank":"1"}]},"next":[137],"score":1},{"path":{"mapping":[{"position":{"node_id":"1724305","is_reverse":true},"edit":[{"from_length":1},{"to_length":1,"sequence":"A"}],"rank":"1"}]},"next":[145],"score":-12},{"path":{"mapping":[{"position":{"node_id":"1724306","offset":"26","is_reverse":true},"edit":[{"to_length":1,"sequence":"A"}],"rank":"1"}]},"next":[145],"score":-6},{"path":{"mapping":[{"position":{"node_id":"1724305","is_reverse":true},"edit":[{"from_length":1,"to_length":1}],"rank":"1"}]},"next":[145],"score":1},{"path":{"mapping":[{"position":{"node_id":"1724308","is_reverse":true},"edit":[{"from_length":1},{"to_length":1,"sequence":"C"}],"rank":"1"}]},"next":[133],"score":-12},{"path":{"mapping":[{"position":{"node_id":"1724307","is_reverse":true},"edit":[{"from_length":1},{"to_length":1,"sequence":"C"}],"rank":"1"}]},"next":[133],"score":-12},{"path":{"mapping":[{"position":{"node_id":"1724307","is_reverse":true},"edit":[{"from_length":1,"to_length":1,"sequence":"C"}],"rank":"1"}]},"next":[133],"score":-4},{"path":{"mapping":[{"position":{"node_id":"1724308","is_reverse":true},"edit":[{"from_length":1,"to_length":1}],"rank":"1"}]},"next":[133],"score":1},{"path":{"mapping":[{"position":{"node_id":"1724306","is_reverse":true},"edit":[{"from_length":26,"to_length":26}],"rank":"1"}]},"next":[144,143,142,136,135,134],"score":26},{"path":{"mapping":[{"position":{"node_id":"1724305","is_reverse":true},"edit":[{"from_length":1},{"to_length":1,"sequence":"A"}],"rank":"1"}]},"next":[137],"score":-12},{"path":{"mapping":[{"position":{"node_id":"1724306","offset":"26","is_reverse":true},"edit":[{"to_length":1,"sequence":"A"}],"rank":"1"}]},"next":[137],"score":-6},{"path":{"mapping":[{"position":{"node_id":"1724305","is_reverse":true},"edit":[{"from_length":1,"to_length":1}],"rank":"1"}]},"next":[137],"score":1},{"path":{"mapping":[{"position":{"node_id":"1724304","is_reverse":true},"edit":[{"from_length":4,"to_length":4}],"rank":"1"}]},"next":[141,140,139,138],"score":4},{"path":{"mapping":[{"position":{"node_id":"1724303","is_reverse":true},"edit":[{"from_length":1,"to_length":1}],"rank":"1"},{"position":{"node_id":"1724301","is_reverse":true},"edit":[{"from_length":1,"to_length":1,"sequence":"C"}],"rank":"2"},{"position":{"node_id":"1724299","is_reverse":true},"edit":[{"from_length":3,"to_length":3}],"rank":"3"},{"position":{"node_id":"1724297","is_reverse":true},"edit":[{"from_length":1,"to_length":1,"sequence":"T"}],"rank":"4"},{"position":{"node_id":"1724296","is_reverse":true},"edit":[{"from_length":4,"to_length":4}],"rank":"5"}]},"score":5},{"path":{"mapping":[{"position":{"node_id":"1724302","is_reverse":true},"edit":[{"from_length":1,"to_length":1,"sequence":"C"}],"rank":"1"},{"position":{"node_id":"1724300","is_reverse":true},"edit":[{"from_length":1,"to_length":1,"sequence":"C"}],"rank":"2"},{"position":{"node_id":"1724299","is_reverse":true},"edit":[{"from_length":3,"to_length":3}],"rank":"3"},{"position":{"node_id":"1724298","is_reverse":true},"edit":[{"from_length":1,"to_length":1}],"rank":"4"},{"position":{"node_id":"1724296","is_reverse":true},"edit":[{"from_length":4,"to_length":4}],"rank":"5"}]},"score":5},{"path":{"mapping":[{"position":{"node_id":"1724303","is_reverse":true},"edit":[{"from_length":1,"to_length":1}],"rank":"1"},{"position":{"node_id":"1724300","is_reverse":true},"edit":[{"from_length":1,"to_length":1,"sequence":"C"}],"rank":"2"},{"position":{"node_id":"1724299","is_reverse":true},"edit":[{"from_length":3,"to_length":3}],"rank":"3"},{"position":{"node_id":"1724298","is_reverse":true},"edit":[{"from_length":1,"to_length":1}],"rank":"4"},{"position":{"node_id":"1724296","is_reverse":true},"edit":[{"from_length":4,"to_length":4}],"rank":"5"}]},"score":10},{"path":{"mapping":[{"position":{"node_id":"1724303","is_reverse":true},"edit":[{"from_length":1,"to_length":1}],"rank":"1"},{"position":{"node_id":"1724301","is_reverse":true},"edit":[{"from_length":1,"to_length":1,"sequence":"C"}],"rank":"2"},{"position":{"node_id":"1724299","is_reverse":true},"edit":[{"from_length":3,"to_length":3}],"rank":"3"},{"position":{"node_id":"1724298","is_reverse":true},"edit":[{"from_length":1,"to_length":1}],"rank":"4"},{"position":{"node_id":"1724296","is_reverse":true},"edit":[{"from_length":4,"to_length":4}],"rank":"5"}]},"score":10},{"path":{"mapping":[{"position":{"node_id":"1724305","is_reverse":true},"edit":[{"from_length":1},{"to_length":1,"sequence":"A"}],"rank":"1"}]},"next":[145],"score":-12},{"path":{"mapping":[{"position":{"node_id":"1724306","offset":"26","is_reverse":true},"edit":[{"to_length":1,"sequence":"A"}],"rank":"1"}]},"next":[145],"score":-6},{"path":{"mapping":[{"position":{"node_id":"1724305","is_reverse":true},"edit":[{"from_length":1,"to_length":1}],"rank":"1"}]},"next":[145],"score":1},{"path":{"mapping":[{"position":{"node_id":"1724304","is_reverse":true},"edit":[{"from_length":4,"to_length":4}],"rank":"1"}]},"next":[149,148,147,146],"score":4},{"path":{"mapping":[{"position":{"node_id":"1724303","is_reverse":true},"edit":[{"from_length":1,"to_length":1}],"rank":"1"},{"position":{"node_id":"1724301","is_reverse":true},"edit":[{"from_length":1,"to_length":1,"sequence":"C"}],"rank":"2"},{"position":{"node_id":"1724299","is_reverse":true},"edit":[{"from_length":3,"to_length":3}],"rank":"3"},{"position":{"node_id":"1724297","is_reverse":true},"edit":[{"from_length":1,"to_length":1,"sequence":"T"}],"rank":"4"},{"position":{"node_id":"1724296","is_reverse":true},"edit":[{"from_length":4,"to_length":4}],"rank":"5"}]},"score":5},{"path":{"mapping":[{"position":{"node_id":"1724302","is_reverse":true},"edit":[{"from_length":1,"to_length":1,"sequence":"C"}],"rank":"1"},{"position":{"node_id":"1724300","is_reverse":true},"edit":[{"from_length":1,"to_length":1,"sequence":"C"}],"rank":"2"},{"position":{"node_id":"1724299","is_reverse":true},"edit":[{"from_length":3,"to_length":3}],"rank":"3"},{"position":{"node_id":"1724298","is_reverse":true},"edit":[{"from_length":1,"to_length":1}],"rank":"4"},{"position":{"node_id":"1724296","is_reverse":true},"edit":[{"from_length":4,"to_length":4}],"rank":"5"}]},"score":5},{"path":{"mapping":[{"position":{"node_id":"1724303","is_reverse":true},"edit":[{"from_length":1,"to_length":1}],"rank":"1"},{"position":{"node_id":"1724300","is_reverse":true},"edit":[{"from_length":1,"to_length":1,"sequence":"C"}],"rank":"2"},{"position":{"node_id":"1724299","is_reverse":true},"edit":[{"from_length":3,"to_length":3}],"rank":"3"},{"position":{"node_id":"1724298","is_reverse":true},"edit":[{"from_length":1,"to_length":1}],"rank":"4"},{"position":{"node_id":"1724296","is_reverse":true},"edit":[{"from_length":4,"to_length":4}],"rank":"5"}]},"score":10},{"path":{"mapping":[{"position":{"node_id":"1724303","is_reverse":true},"edit":[{"from_length":1,"to_length":1}],"rank":"1"},{"position":{"node_id":"1724301","is_reverse":true},"edit":[{"from_length":1,"to_length":1,"sequence":"C"}],"rank":"2"},{"position":{"node_id":"1724299","is_reverse":true},"edit":[{"from_length":3,"to_length":3}],"rank":"3"},{"position":{"node_id":"1724298","is_reverse":true},"edit":[{"from_length":1,"to_length":1}],"rank":"4"},{"position":{"node_id":"1724296","is_reverse":true},"edit":[{"from_length":4,"to_length":4}],"rank":"5"}]},"score":10}],"mapping_quality":60,"start":[0],"paired_read_name":"seed_90_fragment_287_2"}
            
            
            )"; // vim syntax highlighting gives up unless I put "
            
            MultipathAlignment mpaln_pb;
            json2pb(mpaln_pb, multipath_json.c_str(), multipath_json.size());
            multipath_alignment_t mpaln;
            from_proto_multipath_alignment(mpaln_pb, mpaln);
            
            auto alns = optimal_alignments_with_disjoint_subpaths(mpaln, 5);
            
            REQUIRE(alns.size() >= 1);
            REQUIRE(alns.size() <= 5);
            
            // We also should get the right best alignment.
            Alignment best;
            optimal_alignment(mpaln, best);
            
            REQUIRE(best.score() == alns[0].score());
            // TODO: there's no guarantee that we will pick the same path among
            // equally good alignments (i.e. resolve mismatches to SNPs the
            // same way) 
        
        }
        
    
    
        TEST_CASE( "Multipath alignment linearization doesn't produce large fake softclips", "[alignment][multipath][mapping]" ) {
                
            string multipath_json = R"(
{"sequence":"GATTACAA","subpath":[
    {"path":{"mapping":[{"position":{"node_id":"1"},"edit":[{"from_length":4,"to_length":4}],"rank":"1"}]},"next":[5,6,7,8],"score":9},
    {"path":{"mapping":[{"position":{"node_id":"4"},"edit":[{"from_length":1,"to_length":1}],"rank":"1"}]},"next":[9,10],"score":1},
    {"path":{"mapping":[{"position":{"node_id":"4"},"edit":[{"from_length":1,"to_length":1}],"rank":"1"}]},"next":[11,12],"score":1},
    {"path":{"mapping":[{"position":{"node_id":"7"},"edit":[{"from_length":1,"to_length":1}],"rank":"1"}]},"score":6},
    {"path":{"mapping":[{"position":{"node_id":"7"},"edit":[{"from_length":1,"to_length":1}],"rank":"1"}]},"score":6},
    {"path":{"mapping":[{"position":{"node_id":"2"},"edit":[{"from_length":1,"to_length":1}],"rank":"1"}]},"next":[1],"score":1},
    {"path":{"mapping":[{"position":{"node_id":"3"},"edit":[{"from_length":1,"to_length":1,"sequence":"A"}],"rank":"1"}]},"next":[1],"score":-4},
    {"path":{"mapping":[{"position":{"node_id":"2"},"edit":[{"from_length":1,"to_length":1}],"rank":"1"}]},"next":[2],"score":1},
    {"path":{"mapping":[{"position":{"node_id":"3"},"edit":[{"from_length":1,"to_length":1,"sequence":"A"}],"rank":"1"}]},"next":[2],"score":-4},
    {"path":{"mapping":[{"position":{"node_id":"5"},"edit":[{"from_length":1,"to_length":1}],"rank":"1"}]},"next":[3],"score":1},
    {"path":{"mapping":[{"position":{"node_id":"6"},"edit":[{"from_length":1,"to_length":1,"sequence":"A"}],"rank":"1"}]},"next":[3],"score":-4},
    {"path":{"mapping":[{"position":{"node_id":"5"},"edit":[{"from_length":1,"to_length":1}],"rank":"1"}]},"next":[4],"score":1},
    {"path":{"mapping":[{"position":{"node_id":"6"},"edit":[{"from_length":1,"to_length":1,"sequence":"A"}],"rank":"1"}]},"next":[4],"score":-4}
],"start":[0]}
            )";
            
            MultipathAlignment mpaln_pb;
            json2pb(mpaln_pb, multipath_json.c_str(), multipath_json.size());
            multipath_alignment_t mpaln;
            from_proto_multipath_alignment(mpaln_pb, mpaln);
            
            // Topologically sort the multipath_alignment_t so we can linearize it.
            topologically_order_subpaths(mpaln);
            
            // Generate the best linearization with optimal_alignments
            auto alns = optimal_alignments(mpaln, 1);
            REQUIRE(alns.size() == 1);
            
            // Also generate it with just optimal_alignment;
            Alignment opt;
            optimal_alignment(mpaln, opt);
            
            // Make sure they match
            REQUIRE(pb2json(alns[0]) == pb2json(opt));
            
            // Make sure they are all perfect matches, which they should be,
            // because the best linearization is a perfect match.
            REQUIRE(alns[0].path().mapping_size() == 5);
            REQUIRE(alns[0].path().mapping(0).position().node_id() == 1);
            REQUIRE(mapping_is_match(alns[0].path().mapping(0)));
            REQUIRE(alns[0].path().mapping(1).position().node_id() == 2);
            REQUIRE(mapping_is_match(alns[0].path().mapping(1)));
            REQUIRE(alns[0].path().mapping(2).position().node_id() == 4);
            REQUIRE(mapping_is_match(alns[0].path().mapping(2)));
            REQUIRE(alns[0].path().mapping(3).position().node_id() == 5);
            REQUIRE(mapping_is_match(alns[0].path().mapping(3)));
            REQUIRE(alns[0].path().mapping(4).position().node_id() == 7);
            REQUIRE(mapping_is_match(alns[0].path().mapping(4)));
            
            
        }
    }
}



