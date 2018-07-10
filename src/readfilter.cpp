#include "readfilter.hpp"
#include "IntervalTree.h"

#include <fstream>
#include <sstream>

namespace vg {

using namespace std;

bool ReadFilter::trim_ambiguous_ends(xg::XG* index, Alignment& alignment, int k) {
    assert(index != nullptr);

    // Define a way to get node length, for flipping alignments
    function<int64_t(id_t)> get_node_length = [&index](id_t node) {
        return index->node_length(node);
    };

    // Because we need to flip the alignment, make sure it is new-style and
    // doesn't have any Mappings with no Edits.
    for(size_t i = 0; i < alignment.path().mapping_size(); i++) {
        if(alignment.path().mapping(i).edit_size() == 0) {
            // Complain!
            throw runtime_error("Found mapping wit no edits in " + pb2json(alignment));
        }
    }

    // TODO: we're going to flip the alignment twice! This is a waste of time!
    // Define some kind of oriented view or something, or just two duplicated
    // trimming functions, so we can just trim once without flipping.

    // Trim the end
    bool end_changed = trim_ambiguous_end(index, alignment, k);
    // Flip and trim the start
    
    Alignment flipped = reverse_complement_alignment(alignment, get_node_length);
    
    if(trim_ambiguous_end(index, flipped, k)) {
        // The start needed trimming
        
        // Flip the trimmed flipped alignment back
        alignment = reverse_complement_alignment(flipped, get_node_length);
        // We definitely changed something
        return true;
    }
    
    // We maybe changed something
    return end_changed;
    
}

bool ReadFilter::trim_ambiguous_end(xg::XG* index, Alignment& alignment, int k) {
    // What mapping in the alignment is the leftmost one starting in the last k
    // bases? (Except when that would be the first mapping, we use the second.)
    // Start out with it set to the past-the-end value.
    size_t trim_start_mapping = alignment.path().mapping_size();
    
    // How many real non-softclip bases have we seen reading in from the end of
    // the read?
    size_t real_base_count = 0;
    // How many softclip bases have we seen in from the end of the read?
    size_t softclip_base_count = 0;
    for(size_t i = alignment.path().mapping_size() - 1; i != -1 && i != 0; i--) {
        // Scan in from the end of the read.
        
        auto* mapping = alignment.mutable_path()->mutable_mapping(i);
        
        // We should always have edits in our mappings.
        assert(mapping->edit_size() > 0);
        
        for(int j = mapping->edit_size() - 1; j != -1; j--) {
            // Visit every edit in the mapping
            auto& edit = mapping->edit(j);
            
            
            if(real_base_count == 0 && edit.from_length() == 0) {
                // This is a trailing insert. Put it as a softclip
                softclip_base_count += edit.to_length();
            } else {
                // This is some other kind of thing. Record it as real bases.
                real_base_count += edit.to_length();
            }
        }
        
        if(real_base_count <= k) {
            // This mapping starts fewer than k non-softclipped alignment
            // bases from the end of the read.
            trim_start_mapping = i;
        } else {
            // This mapping starts more than k in from the end. So the
            // previous one, if we had one, must be the right one.
            break;
        }
    }
    
    if(trim_start_mapping == alignment.path().mapping_size()) {
        // No mapping was found that starts within the last k non-softclipped
        // bases. So there's nothing to do.
        return false;
    }
    
    if(real_base_count == 0) {
        // We have an anchoring mapping, but all the mappings we could trim are
        // softclips, so there's no point. TODO: will we ever get softclips
        // placed as the only thing on a node?
        return false;
    }
    
    // Which is the last assumed-non-ambiguous mapping from which we can anchor
    // our search?
    size_t root_mapping = trim_start_mapping - 1;
    
    // What's the sequence, including that root node, that we are looking for?
    // We need the sequence of the nodes, rather than the read's sequence,
    // because you can still be ambiguous even if you have a SNP on top of the
    // ambiguous thing.
    
    // We need to ignore all the offsets and from_lengths, except for the from
    // length on the last node to let us know if we end early. It's sort of
    // nonsense to have offsets and non-full from_lengths on internal mappings,
    // and everything is easiest if we use the full length sequence of the root
    // node.
    stringstream target_sequence_stream;
    for(size_t i = root_mapping; i < alignment.path().mapping_size(); i++) {
        // Collect the appropriately oriented from sequence from each mapping
        auto& mapping = alignment.path().mapping(i);
        string sequence = index->node_sequence(mapping.position().node_id());
        if(mapping.position().is_reverse()) {
            // Have it in the right orientation
            sequence = reverse_complement(sequence);
        }
        
        if(i == root_mapping) {
            // Use the full length of the node and ignore any offset
            target_sequence_stream << sequence;
        } else {
            // Use the offset plus the total from_length of all the
            // edits (in case we're the last node and ending early). We made
            // sure all non-root nodes had edits earlier.
            
            size_t from_length = mapping.position().offset();
            for(size_t j = 0; j < mapping.edit_size(); j++) {
                from_length += mapping.edit(j).from_length();
            }
            
            // Put in the sequence that the mapping visits
            target_sequence_stream << sequence.substr(0, from_length);
        }
    }
    string target_sequence = target_sequence_stream.str();
        
#ifdef debug
    #pragma omp critical(cerr)
    cerr << "Need to look for " << target_sequence << " right of mapping " << root_mapping << endl;
#endif
    
    // We're not going to recurse hundreds of nodes deep, so we can use the real
    // stack and a real recursive function.
    
    // Do the DFS into the given node, after already having matched the given
    // number of bases of the target sequence. See if you can match any more
    // bases of the target sequence.
    
    // Return the total number of leaves in all subtrees that match the full
    // target sequence, and the depth in bases of the shallowest point at which
    // multiple subtrees with full lenght matches are unified.

    // We keep a maximum number of visited nodes here, just to prevent recursion
    // from going on forever in worst-case-type graphs
    size_t dfs_visit_count = 0;
    function<pair<size_t, size_t>(id_t, bool, size_t)> do_dfs = 
        [&](id_t node_id, bool is_reverse, size_t matched) -> pair<size_t, size_t> {

        ++dfs_visit_count;
      
        // Grab the node sequence and match more of the target sequence.
        string node_sequence = index->node_sequence(node_id);
        if(is_reverse) {
            node_sequence = reverse_complement(node_sequence);
        }
        
        
#ifdef debug
        #pragma omp critical(cerr)
        cerr << "Node " << node_id <<  " " << (is_reverse ? "rev" : "fwd") << ": "
            << node_sequence << " at offset " << matched << " in " << target_sequence << endl;
#endif
        
        // Now count up the new matches between this node and the target sequence.
        size_t new_matches;
        for(
            // Start with no matches
            new_matches = 0;
            // Keep going as long as we're inside both strings and haven't had a mismatch
            new_matches < node_sequence.size() && 
            matched + new_matches < target_sequence.size() && 
            node_sequence[new_matches] == target_sequence[matched + new_matches];
            // Count up all the matches we find
            new_matches++
        );

        if(matched + new_matches == target_sequence.size()) {
            // We found a tail end of a complete match of the target sequence
            // on this node.
            
#ifdef debug
            #pragma omp critical(cerr)
            cerr << "Node " << node_id << " is a matching leaf" << endl;
#endif
            
            // Return one match and unification at full length (i.e. nothing can
            // be discarded).
            return make_pair(1, target_sequence.size());
        }
        
        if(new_matches < node_sequence.size()) {
            // We didn't make it all the way through this node, nor did we
            // finish the target sequence; there's a mismatch between the node
            // and the target sequence.
            
#ifdef debug
            #pragma omp critical(cerr)
            cerr << "Node " << node_id << " has a mismatch" << endl;
#endif
            
            // If we mismatch, return 0 matches and unification at full length.
            return make_pair(0, target_sequence.size());
        }        
        
        // If we get through the whole node sequence without mismatching or
        // running out of target sequence, keep going.
        
#ifdef debug
        #pragma omp critical(cerr)
        cerr << "Node " << node_id << " has " << new_matches << " internal new matches" << endl;
#endif
    
        // Get all the edges we can take off of the right side of this oriented
        // node.
        auto edges = is_reverse ? index->edges_on_start(node_id) : index->edges_on_end(node_id);
        
#ifdef debug
        #pragma omp critical(cerr)
        cerr << "Recurse into " << edges.size() << " children" << endl;
#endif
        
        // We're going to call all the children and collect the results, and
        // then aggregate them. It might be slightly faster to aggregate while
        // calling, but that might be less clear.
        vector<pair<size_t, size_t>> child_results;
        
        for(auto& edge : edges) {
            // check the user-supplied visit count before recursing any more
            if (dfs_visit_count < defray_count) {
                if(edge.from() == node_id && edge.from_start() == is_reverse) {
                    // The end we are leaving matches this edge's from, so we can
                    // just go to its to end and recurse on it.
                    child_results.push_back(do_dfs(edge.to(), edge.to_end(), matched + node_sequence.size()));
                } else if(edge.to() == node_id && edge.to_end() == !is_reverse) {
                    // The end we are leaving matches this edge's to, so we can just
                    // recurse on its from end.
                    child_results.push_back(do_dfs(edge.from(), !edge.from_start(), matched + node_sequence.size()));
                } else {
                    // XG is feeding us nonsense up with which we should not put.
                    throw runtime_error("Edge " + pb2json(edge) + " does not attach to " +
                                        to_string(node_id) + (is_reverse ? " start" : " end"));
                }
            }
#ifdef debug
            else {
                #pragma omp critical(cerr)
                cerr << "Aborting read filter DFS at node " << node_id << " after " << dfs_visit_count << " visited" << endl;
            }
#endif
            
        }
        
        // Sum up the total leaf matches, which will be our leaf match count.
        size_t total_leaf_matches = 0;
        // If we don't find multiple children with leaf matches, report
        // unification at the min unification depth of any subtree (and there
        // will only be one that isn't at full length).
        size_t children_with_leaf_matches = 0;
        size_t unification_depth = target_sequence.size();
        
        for(auto& result : child_results) {
            total_leaf_matches += result.first;
            if(result.first > 0) {
                children_with_leaf_matches++;
            }
            unification_depth = min(unification_depth, result.second);
        }
        if(children_with_leaf_matches > 1) {
            // If multiple children have nonzero leaf match counts, report
            // unification at the end of this node.
            unification_depth = matched + node_sequence.size();
        }
        
        return make_pair(total_leaf_matches, unification_depth);
    };
    
    // Search from the root mapping's node looking right in its orientation in
    // the mapping
    auto result = do_dfs(alignment.path().mapping(root_mapping).position().node_id(),
        alignment.path().mapping(root_mapping).position().is_reverse(), 0);
    
#ifdef debug
    #pragma omp critical(cerr)
    cerr << "Found " << result.first << " matching leaves with closest unification at " << result.second << endl;
#endif
    
    // We keep this much of the target sequence.
    size_t target_sequence_to_keep = result.second;
    
    if(target_sequence_to_keep == target_sequence.size()) {
        // Nothing to trim!
        return false;
    }
    
    // Figure out how many mappings we need to keep from the root in order to
    // get that much sequence; we know it is either full length or at a mapping
    // boundary. We handle the root special because it's always full length and
    // we have to cut after its end.
    size_t kept_sequence_accounted_for = index->node_length(alignment.path().mapping(root_mapping).position().node_id());
    size_t first_mapping_to_drop;
    for(first_mapping_to_drop = root_mapping + 1;
        first_mapping_to_drop < alignment.path().mapping_size();
        first_mapping_to_drop++) {
        // Consider starting dropping at each mapping after the root.
        if(kept_sequence_accounted_for == target_sequence_to_keep) {
            // OK this mapping really is the first one to drop.
            break;
        } else {
            // Keep going. Account for the sequence from this mapping.
            auto& mapping = alignment.path().mapping(first_mapping_to_drop);
            
            // We know it's not the root mapping, and it can't be the non-full-
            // length end mapping (because we would have kept the full length
            // target sequence and not had to cut). So assume full node is used.
            kept_sequence_accounted_for += index->node_length(mapping.position().node_id());
        }
    }
    
    // OK we know the first mapping to drop. We need to work out the to_size,
    // including all softclips, from there to the end, so we know how much to
    // trim off of the sequence and quality.
    size_t to_length = 0;
    for(size_t i = first_mapping_to_drop; i < alignment.path().mapping_size(); i++) {
        // Go through all the mappings
        auto& mapping = alignment.path().mapping(i);
        for(size_t j = 0; j < mapping.edit_size(); j++) {
            // Add up the to_length of all the edits
            to_length += mapping.edit(j).to_length();
        }
    }
    
    // Trim sequence
    alignment.set_sequence(alignment.sequence().substr(0, alignment.sequence().size() - to_length));
    
    // Trim quality
    if(!alignment.quality().empty()) {
        // If we have a quality, it always ought to have been the same length as the sequence.
        assert(alignment.quality().size() > to_length);
        alignment.set_quality(alignment.quality().substr(0, alignment.quality().size() - to_length));
    }
    
    // Now we can discard the extra mappings
    size_t to_delete = alignment.path().mapping_size() - first_mapping_to_drop;
    alignment.mutable_path()->mutable_mapping()->DeleteSubrange(first_mapping_to_drop, to_delete);
    
    // Now the alignment is fixed!
    return true;
}

// quick and dirty filter to see if removing reads that can slip around
// and still map perfectly helps vg call.  returns true if at either
// end of read sequence, at least k bases are repetitive, checking repeats
// of up to size 2k
bool ReadFilter::has_repeat(Alignment& aln, int k) {
    if (k == 0) {
        return false;
    }
    const string& s = aln.sequence();
    for (int i = 1; i <= 2 * k; ++i) {
        int covered = 0;
        bool ffound = true;
        bool bfound = true;
        for (int j = 1; (ffound || bfound) && (j + 1) * i < s.length(); ++j) {
            ffound = ffound && s.substr(0, i) == s.substr(j * i, i);
            bfound = bfound && s.substr(s.length() - i, i) == s.substr(s.length() - i - j * i, i);
            if (ffound || bfound) {
                covered += i;
            }
        }
        if (covered >= k) {
            return true;
        }
    }
    return false;
}

bool ReadFilter::is_split(xg::XG* index, Alignment& alignment) {
    if(index == nullptr) {
        // Can't tell if the read is split.
        throw runtime_error("XG index required to check for split reads");
    }
    
    
    for(size_t i = 0; i + 1 < alignment.path().mapping_size(); i++) {
        // For each mapping and the one after it
        auto& pos1 = alignment.path().mapping(i).position();
        auto& pos2 = alignment.path().mapping(i + 1).position();
        
        id_t from_id = pos1.node_id();
        bool from_start = pos1.is_reverse();
        id_t to_id = pos2.node_id();
        bool to_end = pos2.is_reverse();
        
        // Can we find the same articulation of the edge as the alignment uses
        bool found = index->has_edge(from_id, from_start, to_id, to_end);   
        
        if(!found) {
            // Check the other articulation of the edge
            std::swap(from_id, to_id);
            std::swap(from_start, to_end);
            from_start = !from_start;
            to_end = !to_end;
            
            
            found = index->has_edge(from_id, from_start, to_id, to_end);
        }
        
        if(!found) {
            // We found a skip!
            if(verbose) {
                cerr << "Warning: read " << alignment.name() << " has an unknown edge "
                    << from_id << " " << from_start << " " << to_id << " " << to_end
                    << ". Removing!" << endl;
            }
            return true;
        } 
    }
    
    // No wandering jumps between nodes found
    return false;
}


bool ReadFilter::sample_bool(double probability) {
    // Have a per-thread RNG
    static thread_local minstd_rand rng;
    
    // Have a distribution with the given success probability
    bernoulli_distribution dist(probability);
    
    // Sample from the distribution, backed by the RNG
    return dist(rng);
}


int ReadFilter::filter(istream* alignment_stream, xg::XG* xindex) {

    // name helper for output
    function<string(int)> chunk_name = [this](int num) -> string {
        stringstream ss;
        ss << outbase << "-" << num << ".gam";
        return ss.str();
    };

    // index regions by their inclusive ranges
    vector<Interval<int, int64_t> > interval_list;
    vector<Region> regions;
    // use strings instead of ofstreams because worried about too many handles
    vector<string> chunk_names;
    vector<bool> chunk_new; // flag if write or append

    // parse a bed, for now this is only way to do regions.  note
    // this operation converts from 0-based BED to 1-based inclusive VCF
    if (!regions_file.empty()) {
        if (outbase.empty()) {
            cerr << "-B option required with -R" << endl;
            return 1;
        }
        parse_bed_regions(regions_file, regions);
        if (regions.empty()) {
            cerr << "No regions read from BED file, doing whole graph" << endl;
        }
    }

    if(defray_length > 0 && xindex == nullptr) {
        cerr << "xg index required for end de-fraying" << endl;
        return 1;
    }
    
    if (regions.empty()) {
        // empty region, do everything
        // we handle empty intervals as special case when looking up, otherwise,
        // need to insert giant interval here.
        chunk_names.push_back(outbase.empty() ? "-" : chunk_name(0));
    } else {
        // otherwise, need to extract regions with xg
        if (xindex == nullptr) {
            cerr << "xg index required for -R option" << endl;
            return 1;
        }
    
        // fill in the map using xg index
        // relies entirely on the assumption that are path chunks
        // are perfectly spanned by an id range
        for (int i = 0; i < regions.size(); ++i) {
            Graph graph;
            int rank = xindex->path_rank(regions[i].seq);
            int64_t path_size = rank == 0 ? 0 : xindex->path_length(regions[i].seq);

            if (regions[i].start >= path_size) {
                cerr << "Unable to find region in index: " << regions[i].seq << ":" << regions[i].start
                     << "-" << regions[i].end << endl;
            } else {
                // clip region on end of path
                regions[i].end = min(path_size - 1, regions[i].end);
                // do path node query
                xindex->get_path_range(regions[i].seq, regions[i].start, regions[i].end, graph);
                if (context_size > 0) {
                    xindex->expand_context(graph, context_size);
                }
            }
            // find node range of graph, without bothering to build vg indices..
            int64_t min_id = numeric_limits<int64_t>::max();
            int64_t max_id = 0;
            for (int j = 0; j < graph.node_size(); ++j) {
                min_id = min(min_id, (int64_t)graph.node(j).id());
                max_id = max(max_id, (int64_t)graph.node(j).id());
            }
            // map the chunk id to a name
            chunk_names.push_back(chunk_name(i));

            // map the node range to the chunk id.
            if (graph.node_size() > 0) {
                interval_list.push_back(Interval<int, int64_t>(min_id, max_id, i));
                assert(chunk_names.size() == i + 1);
            }
        }
    }

    // index chunk regions
    IntervalTree<int, int64_t> region_map(interval_list);

    // which chunk(s) does a gam belong to?
    function<void(Alignment&, vector<int>&)> get_chunks = [&region_map, &regions](Alignment& aln, vector<int>& chunks) {
        // speed up case where no chunking
        if (regions.empty()) {
            chunks.push_back(0);
        } else {
            int64_t min_aln_id = numeric_limits<int64_t>::max();
            int64_t max_aln_id = -1;
            for (int i = 0; i < aln.path().mapping_size(); ++i) {
                const Mapping& mapping = aln.path().mapping(i);
                min_aln_id = min(min_aln_id, (int64_t)mapping.position().node_id());
                max_aln_id = max(max_aln_id, (int64_t)mapping.position().node_id());
            }
            vector<Interval<int, int64_t> > found_ranges;
            region_map.findOverlapping(min_aln_id, max_aln_id, found_ranges);
            for (auto& interval : found_ranges) {
                chunks.push_back(interval.value);
            }
        }
    };

    // buffered output (one buffer per chunk), one set of chunks per thread
    // buffer[THREAD][CHUNK] = vector<Alignment>
    vector<vector<vector<Alignment> > > buffer(threads);
    for (int i = 0; i < buffer.size(); ++i) {
        buffer[i].resize(chunk_names.size());
    }
    
    static const int buffer_size = 1000; // we let this be off by 1

    // remember if write or append
    vector<bool> chunk_append(chunk_names.size(), append_regions);

    // flush a buffer specified by cur_buffer to target in chunk_names, and clear it
    function<void(int, int)> flush_buffer = [&buffer, &chunk_names, &chunk_append](int tid, int cur_buffer) {
        ofstream outfile;
        auto& outbuf = chunk_names[cur_buffer] == "-" ? cout : outfile;
        if (chunk_names[cur_buffer] != "-") {
            outfile.open(chunk_names[cur_buffer], chunk_append[cur_buffer] ? ios::app : ios_base::out);
            chunk_append[cur_buffer] = true;
        }
        function<Alignment&(uint64_t)> write_buffer = [&buffer, &tid, &cur_buffer](uint64_t i) -> Alignment& {
            return buffer[tid][cur_buffer][i];
        };
        stream::write(outbuf, buffer[tid][cur_buffer].size(), write_buffer);
        buffer[tid][cur_buffer].clear();
    };

    // add alignment to all appropriate buffers, flushing as necessary
    function<void(int, Alignment&, const vector<int>&)> update_buffers = [
        &buffer, &region_map, &get_chunks, &flush_buffer](int tid, Alignment& aln,
                                                          const vector<int>& aln_chunks) {
        for (auto chunk : aln_chunks) {
            buffer[tid][chunk].push_back(aln);
            if (buffer[tid][chunk].size() >= buffer_size) {
                // flush buffer (could get fancier and allow parallel writes to different
                // files, but unlikely to be worth effort as we're mostly trying to
                // speed up defray and not write IO)
#pragma omp critical (ReadFilter_flush_buffer)
                {
                    flush_buffer(tid, chunk);
                }
            }
        }
    };

    // keep counts of what's filtered to report (in verbose mode)
    vector<Counts> counts_vec(threads);
            
    // we assume that every primary alignment has 0 or 1 secondary alignment
    // immediately following in the stream
    function<void(Alignment&)> lambda = [&](Alignment& aln) {
        int tid = omp_get_thread_num();        
        Counts& counts = counts_vec[tid];
        double score = (double)aln.score();
        double denom = aln.sequence().length();
        // toggle substitution score
        if (sub_score == true) {
            // hack in ident to replace old counting logic.
            score = aln.identity() * aln.sequence().length();
            assert(score <= denom);
        } else if (rescore == true) {
            // We need to recalculate the score with the base aligner always
            const static Aligner unadjusted;
            BaseAligner* aligner = (BaseAligner*)&unadjusted;
            
            // Rescore and assign the score
            aln.set_score(aligner->score_ungapped_alignment(aln));
            // Also use the score
            score = aln.score();
        }

        // toggle absolute or fractional score
        if (frac_score == true) {
            if (denom > 0.) {
                score /= denom;
            }
            else {
                assert(score == 0.);
            }
        }
        // compute overhang
        int overhang = 0;
        if (aln.path().mapping_size() > 0) {
            const auto& left_mapping = aln.path().mapping(0);
            if (left_mapping.edit_size() > 0) {
                overhang = left_mapping.edit(0).to_length() - left_mapping.edit(0).from_length();
            }
            const auto& right_mapping = aln.path().mapping(aln.path().mapping_size() - 1);
            if (right_mapping.edit_size() > 0) {
                const auto& edit = right_mapping.edit(right_mapping.edit_size() - 1);
                overhang = max(overhang, edit.to_length() - edit.from_length());
            }
        } else {
            overhang = aln.sequence().length();
        }
        // compute end matches.
        int end_matches = 0;
        // from the left
        for (int i = 0; i < aln.path().mapping_size() && end_matches < min_end_matches; ++i) {
            for (int j = 0; j < aln.path().mapping(i).edit_size() && end_matches < min_end_matches; ++j) {
                const Edit& edit = aln.path().mapping(i).edit(j);
                if (edit.from_length() == edit.to_length() && edit.sequence().empty()) {
                    end_matches += edit.to_length();
                } else {
                    i = aln.path().mapping_size();
                    break;
                }
            }
        }
        if (end_matches >= min_end_matches) {
            end_matches = 0;
            // from the right
            for (int i = aln.path().mapping_size() - 1; i >= 0 && end_matches < min_end_matches; --i) {
                for (int j = aln.path().mapping(i).edit_size() - 1; j >= 0 && end_matches < min_end_matches; --j) {
                    const Edit& edit = aln.path().mapping(i).edit(j);
                    if (edit.from_length() == edit.to_length() && edit.sequence().empty()) {
                        end_matches += edit.to_length();
                    } else {
                        i = -1;
                        break;
                    }
                }
            }
        }                

        // offset in count tuples
        int co = aln.is_secondary() ? 1 : 0;
        
        ++counts.read[co];
        bool keep = true;
        // filter (current) alignment
        if (!name_prefix.empty() && !std::equal(name_prefix.begin(), name_prefix.end(), aln.name().begin())) {
            // There's a prefix and a mismatch against it
            ++counts.wrong_name[co];
            keep = false;    
        }
        if ((keep || verbose) && ((aln.is_secondary() && score < min_secondary) ||
            (!aln.is_secondary() && score < min_primary))) {
            ++counts.min_score[co];
            keep = false;
        }
        if ((keep || verbose) && overhang > max_overhang) {
            ++counts.max_overhang[co];
            keep = false;
        }
        if ((keep || verbose) && end_matches < min_end_matches) {
            ++counts.min_end_matches[co];
            keep = false;
        }
        if ((keep || verbose) && aln.mapping_quality() < min_mapq) {
            ++counts.min_mapq[co];
            keep = false;
        }

        // do region check before heavier filters
        vector<int> aln_chunks;
        if (keep || verbose) {
            get_chunks(aln, aln_chunks);
            if (aln_chunks.empty()) {
                keep = false;
            }
        }
        
        if ((keep || verbose) && drop_split && is_split(xindex, aln)) {
            ++counts.split[co];
            keep = false;
        }
        if ((keep || verbose) && has_repeat(aln, repeat_size)) {
            ++counts.repeat[co];
            keep = false;
        }
        if ((keep || verbose) && defray_length && trim_ambiguous_ends(xindex, aln, defray_length)) {
            ++counts.defray[co];
            // We keep these, because the alignments get modified.
        }
        if ((keep || verbose) && downsample_probability != 1.0 && !sample_bool(downsample_probability)) {
            ++counts.random[co];
            keep = false;
        }
        if (!keep) {
            ++counts.filtered[co];
        }

        // add to write buffer
        if (keep) {
            update_buffers(tid, aln, aln_chunks);
        }
    };
    stream::for_each_parallel(*alignment_stream, lambda);

    for (int tid = 0; tid < buffer.size(); ++tid) {
        for (int chunk = 0; chunk < buffer[tid].size(); ++chunk) {
            if (buffer[tid][chunk].size() > 0 || 
                // we deliberately write empty gams at this point for empty chunks:
                (chunk_append[chunk] == false && chunk_names[chunk] != "-" )) {
                flush_buffer(tid, chunk);
            }
        }
    }

    if (verbose) {
        Counts& counts = counts_vec[0];
        for (int i = 1; i < counts_vec.size(); ++i) {
            counts += counts_vec[i];
        }
        size_t tot_reads = counts.read[0] + counts.read[1];
        size_t tot_filtered = counts.filtered[0] + counts.filtered[1];
        cerr << "Total Filtered (primary):          " << counts.filtered[0] << " / "
             << counts.read[0] << endl
             << "Total Filtered (secondary):        " << counts.filtered[1] << " / "
             << counts.read[1] << endl
             << "Read Name Filter (primary):        " << counts.wrong_name[0] << endl
             << "Read Name Filter (secondary):      " << counts.wrong_name[1] << endl
             << "Min Identity Filter (primary):     " << counts.min_score[0] << endl
             << "Min Identity Filter (secondary):   " << counts.min_score[1] << endl
             << "Max Overhang Filter (primary):     " << counts.max_overhang[0] << endl
             << "Max Overhang Filter (secondary):   " << counts.max_overhang[1] << endl
             << "Min End Match Filter (primary):    " << counts.min_end_matches[0] << endl
             << "Min End Match Filter (secondary):  " << counts.min_end_matches[1] << endl            
             << "Split Read Filter (primary):       " << counts.split[0] << endl
             << "Split Read Filter (secondary):     " << counts.split[1] << endl
             << "Repeat Ends Filter (primary):      " << counts.repeat[0] << endl
             << "Repeat Ends Filter (secondary):    " << counts.repeat[1] << endl
             << "Min Quality Filter (primary):      " << counts.min_mapq[0] << endl
             << "Min Quality Filter (secondary):   " << counts.min_mapq[1] << endl
             << "Random Filter (primary):      " << counts.random[0] << endl
             << "Random Filter (secondary):   " << counts.random[1] << endl
                        
            
             << endl;
    }
    
    return 0;

}

}
