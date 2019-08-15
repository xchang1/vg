#include "sampler.hpp"

//#define debug_ngs_sim

namespace vg {


void Sampler::set_source_paths(const vector<string>& source_paths,
                               const vector<pair<string, double>>& transcript_expressions) {
    if (!source_paths.empty() && !transcript_expressions.empty()) {
        cerr << "error:[Sampler] cannot sample from list of paths and from list of transcripts simultaneously" << endl;
        exit(1);
    }
    else if (!transcript_expressions.empty()) {
        this->source_paths.clear();
        vector<double> expression_values;
        for (const pair<string, double>& transcript_expression : transcript_expressions) {
            this->source_paths.push_back(transcript_expression.first);
            expression_values.push_back(transcript_expression.second);
        }
        path_sampler = vg::discrete_distribution<>(expression_values.begin(), expression_values.end());
    }
    else if (!source_paths.empty()) {
        this->source_paths = source_paths;
        vector<size_t> path_lengths;
        for (auto& source_path : source_paths) {
            path_lengths.push_back(xgidx->get_path_length(xgidx->get_path_handle(source_path)));
        }
        path_sampler = vg::discrete_distribution<>(path_lengths.begin(), path_lengths.end());
    }
    else {
        path_sampler = vg::discrete_distribution<>();
    }
}
    

/// We have a helper function to convert path positions and orientations to
/// pos_t values.
pos_t position_at(PathPositionHandleGraph* xgidx, const string& path_name, const size_t& path_offset, bool is_reverse) {
    path_handle_t path_handle = xgidx->get_path_handle(path_name);
    step_handle_t step = xgidx->get_step_at_position(path_handle, path_offset);
    handle_t handle = xgidx->get_handle_of_step(step);
    
    // Work out where in that mapping we should be.
    size_t node_offset = path_offset - xgidx->get_position_of_step(step);

    if (is_reverse) {
        // Flip the node offset around to be from the end and not the start
        node_offset = xgidx->get_length(handle) - node_offset - 1;
    }

    // Make a pos_t for where we are, on the appropriate strand
    pos_t pos = make_pos_t(xgidx->get_id(handle), xgidx->get_is_reverse(handle) != is_reverse, node_offset);
    
    return pos;
}

pos_t Sampler::position(void) {
    // We sample from the entire graph sequence, 1-based.
    vg::uniform_int_distribution<size_t> xdist(1, total_seq_length);
    size_t offset = xdist(rng);
    id_t id = dynamic_cast<VectorizableHandleGraph*>(xgidx)->node_at_vector_offset(offset);
    vg::uniform_int_distribution<size_t> flip(0, 1);
    bool rev = forward_only ? false : flip(rng);
    // 1-0 base conversion
    size_t node_offset = offset - dynamic_cast<VectorizableHandleGraph*>(xgidx)->node_vector_offset(id) - 1;
    // Ignore flipping the node offset because we're uniform over both strands
    return make_pos_t(id, rev, node_offset);
}

string Sampler::sequence(size_t length) {
    pos_t pos = position();
    cerr << pos << endl;
    string seq;
    while (seq.size() < length) {
        auto nextc = next_pos_chars(pos);
        if (nextc.empty()) break;
        vector<pos_t> nextp;
        for (auto& n : nextc) nextp.push_back(n.first);
        // pick one at random
        vg::uniform_int_distribution<int> next_dist(0, nextc.size()-1);
        // update our position
        pos = nextp.at(next_dist(rng));
        // append to our sequence
        seq += nextc[pos];
    }
    return seq;
}


vector<Edit> Sampler::mutate_edit(const Edit& edit,
                                  const pos_t& position,
                                  double base_error,
                                  double indel_error,
                                  const string& bases,
                                  vg::uniform_real_distribution<double>& rprob,
                                  vg::uniform_int_distribution<int>& rbase) {

    // we will build up a mapping representing the modified edit
    Mapping new_mapping;
    //*new_mapping.mutable_position() = make_position(position);
    // determine to-length of edit
    size_t to_length = edit.to_length();
    // we will keep track of the current base using this
    pos_t curr_pos = position;
    
    // We punt if we aren't a well-defined kind of edit; from and to lengths
    // must be equal, or from length must be 0.
    if (edit_is_match(edit) || edit_is_sub(edit)
        || edit_is_insertion(edit)) {
        
#ifdef debug
        cerr << "Handle edit " << pb2json(edit) << endl;
#endif
        
        // distribute mutations across the to_length
        for (size_t k = 0; k < to_length;) {
            // Until we've consumed all the characters produced by the original edit
        
            // Get the character that the original edit has here
            char c;
            if (edit_is_match(edit)) {
                // It's not stored in the edit, so find it in the reference
                c = pos_char(curr_pos);
            } else {
                // It must be stored in the edit itself
                c = edit.sequence().at(k);
            }
            
#ifdef debug
            cerr << "At to_length position " << k << " of " << edit.to_length() << endl;
#endif
            
            // This is the edit we're going to create, to deside the fate of
            // this character from the old edit.
            Edit* e = nullptr;
            
            if (rprob(rng) <= base_error) {
                // We should do a substitution relative to the old edit.
                
                // pick another base than what c is
                char n;
                do {
                    n = bases[rbase(rng)];
                } while (n == c);
                // make the edit for the sub
                e = new_mapping.add_edit();
                string s(1, n);
                e->set_sequence(s);
                if (!edit_is_insertion(edit)) {
                    // We should stay aligned against whatever we were aligned
                    // against before.
                    e->set_from_length(1);
                }
                e->set_to_length(1);
                
#ifdef debug
                cerr << "Produced relative substitution " << pb2json(*e) << endl;
#endif
            } else if (rprob(rng) <= indel_error) {
                // We have an indel.
                // Note that we're using a simple geometric indel dsitribution here
                if (rprob(rng) < 0.5) {
                    // This should be an insertion relative to the original edit.
                    char n = bases[rbase(rng)];
                    e = new_mapping.add_edit();
                    string s(1, n);
                    e->set_sequence(s);
                    e->set_to_length(1);
                    
#ifdef debug
                    cerr << "Produced relative insertion " << pb2json(*e) << endl;
#endif
                    
                    // None of the graph sequence is consumed. But we are
                    // inserting relative to the old edit, so we want to give
                    // the base we just inserted before another shot. We'll
                    // continue now, before advancing our position in the edit
                    // we're modifying.
                    continue;
                } else {
                    // This should be a deletion relative to the edit we're
                    // modifying.
                    
                    // The old edit base isn't going to come out, but we need to
                    // consume it anyway.
                    k++;
                    
                    if (edit_is_insertion(edit)) {
                        // We just need to consume this base of the old edit and
                        // not emit any new edit.
#ifdef debug
                        cerr << "Skipped base for relative deletion" << endl;
#endif
                        continue;
                    } else {
                        // We have to delete the base that was substituted or
                        // matched against.
                        e = new_mapping.add_edit();
                        e->set_from_length(1);
#ifdef debug
                        cerr << "Produced relative deletion " << pb2json(*e) << endl;
#endif
                    }
                }
            } else {
                // make the edit for the 1bp match relative to the old edit
                // (which may actually be an insertion edit or a substitution
                // edit relative to the graph)
                e = new_mapping.add_edit();
                if (!edit_is_match(edit)) {
                    // We're not a match, so we need the sequence set.
                    string s(1, c);
                    e->set_sequence(s);
                }
                if (!edit_is_insertion(edit)) {
                    // We should stay aligned against whatever we were aligned
                    // against before.
                    e->set_from_length(1);
                }
                e->set_to_length(1);
                
#ifdef debug
                cerr << "Produced relative match " << pb2json(*e) << endl;
#endif
            }
            
            // Now advance in the old edit by the number of old edit bases used
            // in this new edit.
            k += e->to_length();
            // And in the graph by the number of graph bases consumed
            get_offset(curr_pos) += e->from_length();

        }
    } else if (edit_is_deletion(edit)) {
        // special case: 0 (deletion)
        // Just copy over the deletion edit
        *(new_mapping.add_edit()) = edit;
    }
    
#ifdef debug
    cerr << "Before merging adjacent edits: " << pb2json(new_mapping) << endl;
#endif

    // Merge adjacent edits, but don't get rid of leading or trailing deletions
    // (as with simplify), because we want a path that reflects the real
    // simulated history and because we don't output any modifications to the
    // position out of this function.
    new_mapping = merge_adjacent_edits(new_mapping);
    
#ifdef debug
    cerr << "Replacing " << pb2json(edit) << " with " << pb2json(new_mapping) << endl;
#endif
    
    assert(mapping_from_length(new_mapping) == edit.from_length());
    
    // copy the new edits
    vector<Edit> new_edits;
    for (size_t i = 0; i < new_mapping.edit_size(); ++i) {
        new_edits.push_back(new_mapping.edit(i));
    }
    // and send them back
    return new_edits;
}

Alignment Sampler::mutate(const Alignment& aln,
                          double base_error,
                          double indel_error) {

    if (base_error == 0 && indel_error == 0) return aln;

    string bases = "ATGC";
    vg::uniform_real_distribution<double> rprob(0, 1);
    vg::uniform_int_distribution<int> rbase(0, 3);

    Alignment mutaln;
    for (size_t i = 0; i < aln.path().mapping_size(); ++i) {
        auto& orig_mapping = aln.path().mapping(i);
        Mapping* new_mapping = mutaln.mutable_path()->add_mapping();
        *new_mapping->mutable_position() = orig_mapping.position();
        // for each edit in the mapping
        for (size_t j = 0; j < orig_mapping.edit_size(); ++j) {
            auto& orig_edit = orig_mapping.edit(j);
            auto new_edits = mutate_edit(orig_edit, make_pos_t(orig_mapping.position()),
                                         base_error, indel_error,
                                         bases, rprob, rbase);
            for (auto& edit : new_edits) {
                *new_mapping->add_edit() = edit;
            }
        }
    }
    
    // Don't simplify the alignment, because it's nice to see the deletions as
    // originally generated. Also, simplified alignments won't validate.
    
    // re-derive the alignment's sequence.
    mutaln.set_sequence(alignment_seq(mutaln));
    mutaln.set_name(aln.name());
    mutaln.clear_refpos();
    algorithms::annotate_with_initial_path_positions(*xgidx, mutaln);
    return mutaln;
}

string Sampler::alignment_seq(const Alignment& aln) {
    // get the graph corresponding to the alignment path
    VG g;
    for (int i = 0; i < aln.path().mapping_size(); ++ i) {
        auto& m = aln.path().mapping(i);
        if (m.has_position() && m.position().node_id()) {
            auto id = aln.path().mapping(i).position().node_id();
            algorithms::extract_id_range(*xgidx, id, id, g);
        }
    }
    algorithms::expand_subgraph_by_steps(*xgidx, g, 2);
    return g.path_string(aln.path());
}

vector<Alignment> Sampler::alignment_pair(size_t read_length, size_t fragment_length, double fragment_std_dev, double base_error, double indel_error) {
    // simulate forward/reverse pair by first simulating a long read
    vg::normal_distribution<> norm_dist(fragment_length, fragment_std_dev);
    // bound at read length so we always get enough sequence
    int frag_len = max((int)read_length, (int)round(norm_dist(rng)));
    auto fragment = alignment_with_error(frag_len, base_error, indel_error);
    // then taking the ends
    auto fragments = alignment_ends(fragment, read_length, read_length);
    auto& aln1 = fragments.front();
    auto& aln2 = fragments.back();
    { // name the alignments
        string data;
        aln1.SerializeToString(&data);
        aln2.SerializeToString(&data);
        int n;
#pragma omp critical(nonce)
        n = nonce++;
        data += std::to_string(n);
        const string hash = sha1head(data, 16);
        aln1.set_name(hash + "_1");
        aln2.set_name(hash + "_2");
    }
    // set the appropriate flags for pairing
    aln1.mutable_fragment_next()->set_name(aln2.name());
    aln2.mutable_fragment_prev()->set_name(aln1.name());
    // reverse complement the back fragment
    fragments.back() = reverse_complement_alignment(fragments.back(),
                                                  (function<int64_t(int64_t)>) ([&](int64_t id) {
                                                          return (int64_t)node_length(id);
                                                      }));
    return fragments;
}

// generates a perfect alignment from the graph or the selected source path if one exists
Alignment Sampler::alignment(size_t length) {
    if (source_paths.empty()) {
        return alignment_to_graph(length);
    } else {
        return alignment_to_path(source_paths[path_sampler(rng)], length);
    }
}

// generates a perfect alignment from the graph
Alignment Sampler::alignment_to_path(const string& source_path, size_t length) {

    // Pick a starting point along the path and an orientation
    path_handle_t path_handle = xgidx->get_path_handle(source_path);
    uint64_t path_length = xgidx->get_path_length(path_handle);
    vg::uniform_int_distribution<size_t> xdist(0, path_length - 1);
    size_t path_offset = xdist(rng);
    vg::uniform_int_distribution<size_t> flip(0, 1);
    bool rev = forward_only ? false : flip(rng);
    
    // We will fill in this string
    string seq;
    // And this Alignment
    Alignment aln;
    Path* path = aln.mutable_path();
    
    while (seq.size() < length) {

        // Make a pos_t for where we are, on the appropriate strand
        pos_t pos = position_at(xgidx, source_path, path_offset, rev);

        // Add that character to the sequence
        seq.push_back(pos_char(pos));
        
        // Add a perfect match edit for 1 base
        Mapping* mapping = path->add_mapping();
        *mapping->mutable_position() = make_position(pos);
        Edit* edit = mapping->add_edit();
        edit->set_from_length(1);
        edit->set_to_length(1);
        
        // Advance along the path in the appropriate direction
        if (rev) {
            if (path_offset == 0) {
                // Out of path!
                break;
            }
            path_offset--;
        } else {
            if (path_offset == path_length - 1) {
                // Out of path!
                break;
            }
            path_offset++;
        }
    }
    
    // save our sequence in the alignment
    aln.set_sequence(seq);
    // Simplify the alignment to merge redundant mappings. There are no deletions to get removed.
    aln = simplify(aln); 
    
    { // name the alignment
        string data;
        aln.SerializeToString(&data);
        int n;
#pragma omp critical(nonce)
        n = nonce++;
        data += std::to_string(n);
        const string hash = sha1head(data, 16);
        aln.set_name(hash);
    }
    // And set its identity
    aln.set_identity(identity(aln.path()));
    aln.clear_refpos();
    algorithms::annotate_with_initial_path_positions(*xgidx, aln);
    return aln;
}

// generates a perfect alignment from the graph
Alignment Sampler::alignment_to_graph(size_t length) {
    string seq;
    Alignment aln;
    Path* path = aln.mutable_path();
    pos_t pos = position();
    char c = pos_char(pos);
    // we do something wildly inefficient but conceptually clean
    // for each position in the mapping we add a mapping
    do {
        // add in the char for the current position
        seq += c;
        Mapping* mapping = path->add_mapping();
        *mapping->mutable_position() = make_position(pos);
        Edit* edit = mapping->add_edit();
        edit->set_from_length(1);
        edit->set_to_length(1);
        // decide the next position
        auto nextc = next_pos_chars(pos);
        // no new positions mean we are done; we've reached the end of the graph
        if (nextc.empty()) break;
        // what positions do we go to next?
        vector<pos_t> nextp;
        for (auto& n : nextc) nextp.push_back(n.first);
        // pick one at random
        vg::uniform_int_distribution<int> next_dist(0, nextc.size()-1);
        // update our position
        pos = nextp.at(next_dist(rng));
        // update our char
        c = nextc[pos];
    } while (seq.size() < length);
    // save our sequence in the alignment
    aln.set_sequence(seq);
    // Simplify the alignment to merge redundant mappings. There are no deletions to get removed.
    aln = simplify(aln); 
    
    { // name the alignment
        string data;
        aln.SerializeToString(&data);
        int n;
#pragma omp critical(nonce)
        n = nonce++;
        data += std::to_string(n);
        const string hash = sha1head(data, 16);
        aln.set_name(hash);
    }
    // And set its identity
    aln.set_identity(identity(aln.path()));
    algorithms::annotate_with_initial_path_positions(*xgidx, aln);
    return aln;
}

Alignment Sampler::alignment_with_error(size_t length,
                                        double base_error,
                                        double indel_error) {
    size_t maxiter = 100;
    Alignment aln;
    size_t iter = 0;
    if (base_error > 0 || indel_error > 0) {
        // sample a longer-than necessary alignment, then trim
        while (iter++ < maxiter) {
            aln = mutate(
                alignment(length + 2 * ((double) length * indel_error)),
                base_error, indel_error);
            if (!(no_Ns && aln.sequence().find('N') != string::npos)) {
                if (aln.sequence().size() == length) {
                    break;
                } else if (aln.sequence().size() > length) {
                    aln = strip_from_end(aln, aln.sequence().size() - length);
                    break;
                }
            }
        }
    } else {
        size_t iter = 0;
        while (iter++ < maxiter) {
            aln = alignment(length);
            if (aln.sequence().size() == length
                && !(no_Ns && aln.sequence().find('N') != string::npos)) {
                break;
            }
        }
    }
    if (iter == maxiter) {
        cerr << "[vg::Sampler] Warning: could not generate alignment of sufficient length. "
             << "Graph may be too small, or indel rate too high." << endl;
    }
    aln.set_identity(identity(aln.path()));
    
    // Check the alignment to make sure we didn't mess it up
    assert(is_valid(aln));
    algorithms::annotate_with_initial_path_positions(*xgidx, aln);
    return aln;
}

size_t Sampler::node_length(id_t id) {
    return xgidx->get_length(xgidx->get_handle(id));
}

char Sampler::pos_char(pos_t pos) {
    return xgidx->get_base(xgidx->get_handle(id(pos), is_rev(pos)), offset(pos));
}

map<pos_t, char> Sampler::next_pos_chars(pos_t pos) {
    return algorithms::next_pos_chars(*xgidx, pos);
}

bool Sampler::is_valid(const Alignment& aln) {
    for (auto i = 0; i + 1 < aln.path().mapping_size(); i++) {
        // For each mapping except the very last (which might not use its whole
        // node)
        auto& mapping = aln.path().mapping(i);
        
        // What's the number of bases it consumes?
        auto observed_from = mapping_from_length(mapping);
        
        // How many bases are accounted for?
        auto accounted_bases = observed_from + mapping.position().offset();
        
        // How many bases need to be accounted for?
        auto expected_bases = xgidx->get_length(xgidx->get_handle(mapping.position().node_id()));
        
        if (accounted_bases != expected_bases) {
            cerr << "[vg::Sampler] Warning: alignment mapping " << i << " accounts for "
                << accounted_bases << " bases of graph sequence, but needs to account for "
                << expected_bases << endl;
            cerr << pb2json(aln) << endl;
            return false;
        }
    }
    
    // For now, we just say an alignment is valid if it accounts for all the
    // bases on its source nodes.
    return true;
}


const string NGSSimulator::alphabet = "ACGT";

NGSSimulator::NGSSimulator(PathPositionHandleGraph& graph,
                           const string& ngs_fastq_file,
                           bool interleaved_fastq,
                           const vector<string>& source_paths_input,
                           const vector<pair<string, double>>& transcript_expressions,
                           double substition_polymorphism_rate,
                           double indel_polymorphism_rate,
                           double indel_error_proportion,
                           double insert_length_mean,
                           double insert_length_stdev,
                           double error_multiplier,
                           bool retry_on_Ns,
                           size_t seed) :
      graph(graph)
    , node_cache(100)
    , edge_cache(100)
    , sub_poly_rate(substition_polymorphism_rate)
    , indel_poly_rate(indel_polymorphism_rate)
    , indel_error_prop(indel_error_proportion)
    , insert_mean(insert_length_mean)
    , insert_sd(insert_length_stdev)
    , retry_on_Ns(retry_on_Ns)
    , prng(seed ? seed : random_device()())
    , strand_sampler(0, 1)
    , background_sampler(0, alphabet.size() - 1)
    , mut_sampler(0, alphabet.size() - 2)
    , prob_sampler(0.0, 1.0)
    , insert_sampler(insert_length_mean, insert_length_stdev)
    , seed(seed)
    , source_paths(source_paths_input)
    , joint_initial_distr(seed - 1)
{
    graph.for_each_handle([&](const handle_t& handle) {
            total_seq_length += graph.get_length(handle);
        });
    
    if (!source_paths.empty() && !transcript_expressions.empty()) {
        cerr << "error:[NGSSimulator] cannot simultaneously limit sampling to paths and match an expresssion profile" << endl;
        exit(1);
    }
    
    if (substition_polymorphism_rate < 0.0 || substition_polymorphism_rate > 1.0
        || indel_polymorphism_rate < 0.0 || indel_polymorphism_rate > 1.0
        || indel_error_proportion < 0.0 || indel_error_proportion > 1.0) {
        cerr << "error:[NGSSimulator] All proportions must be between 0.0 and 1.0" << endl;
        exit(1);
    }
    
    if (substition_polymorphism_rate + indel_polymorphism_rate > 1.0) {
        cerr << "error:[NGSSimulator] Indel polymorphism rate and substitution polymorphism rate cannot sum to greater than 1.0" << endl;
        exit(1);
    }
    
    if (insert_length_mean <= 0.0) {
        cerr << "error:[NGSSimulator] Mean insert length must be positive" << endl;
        exit(1);
    }
    
    if (insert_length_stdev < 0.0) {
        cerr << "error:[NGSSimulator] Insert length standard deviation must be positive" << endl;
        exit(1);
    }
    
    if (insert_length_mean < 5.0 * insert_length_stdev) {
        cerr << "warning:[NGSSimulator] Recommended that insert length mean (" << insert_length_mean
            << ") > 5 * insert length standard deviation (" << insert_length_stdev << ")" << endl;
    }
    
    if (source_paths.empty() && transcript_expressions.empty()) {
        start_pos_samplers.emplace_back(1, total_seq_length);
    }
    else if (!source_paths.empty()) {
        vector<size_t> path_sizes;
        for (const auto& source_path : source_paths) {
            path_sizes.push_back(graph.get_path_length(graph.get_path_handle(source_path)));
            start_pos_samplers.emplace_back(0, path_sizes.back() - 1);
        }
        path_sampler = vg::discrete_distribution<>(path_sizes.begin(), path_sizes.end());
    }
    else {
        vector<double> expression_values;
        for (const pair<string, double>& transcript_expression : transcript_expressions) {
            source_paths.push_back(transcript_expression.first);
            start_pos_samplers.emplace_back(0, graph.get_path_length(graph.get_path_handle(transcript_expression.first)) - 1);
            expression_values.push_back(transcript_expression.second);
        }
        path_sampler = vg::discrete_distribution<>(expression_values.begin(), expression_values.end());
    }
    
    // memoize phred conversions
    phred_prob.resize(256);
    for (int i = 1; i < phred_prob.size(); i++) {
        phred_prob[i] = error_multiplier * phred_to_prob(i);
    }
    
    for (size_t i = 0; i < alphabet.size(); i++) {
        mutation_alphabets[alphabet[i]] = string();
        for (size_t j = 0; j < alphabet.size(); j++) {
            if (j == i) {
                continue;
            }
            mutation_alphabets[alphabet[i]].push_back(alphabet[j]);
        }
    }
    
    unordered_map<size_t, size_t> length_count;
    if (interleaved_fastq) {
        fastq_paired_interleaved_for_each(ngs_fastq_file, [&](const Alignment& aln_1, const Alignment& aln_2) {
            length_count[aln_1.quality().size()]++;
            length_count[aln_2.quality().size()]++;
            record_read_pair_quality(aln_1, aln_2);
        });
    }
    else {
        fastq_unpaired_for_each(ngs_fastq_file, [&](const Alignment& aln) {
            length_count[aln.quality().size()]++;
            record_read_quality(aln);
        });
    }
    
    size_t modal_length = 0;
    size_t modal_length_count = 0;
    size_t total_reads = 0;
    for (const pair<size_t, size_t>& length_record : length_count) {
        if (length_record.second > modal_length_count) {
            modal_length_count = length_record.second;
            modal_length = length_record.first;
        }
        total_reads += length_record.second;
    }
    
    if (((double) modal_length_count) / total_reads < 0.5) {
        cerr << "warning:[NGSSimulator] Auto-detected read length of " << modal_length << " encompasses less than half of training reads, NGSSimulator is optimized for training data in which most reads are the same length" << endl;
    }
    
    if (modal_length > insert_length_mean - 2.0 * insert_length_stdev) {
        cerr << "warning:[NGSSimulator] Auto-detected read length of " << modal_length << " is long compared to mean insert length " << insert_length_mean << " and standard deviation " << insert_length_stdev << ", sampling may take additional time and statistical properties of insert length distribution may not reflect input parameters" << endl;
    }
    
    while (transition_distrs_1.size() > modal_length) {
        transition_distrs_1.pop_back();
    }
    while (transition_distrs_2.size() > modal_length) {
        transition_distrs_2.pop_back();
    }
    
    if (transition_distrs_1.size() != transition_distrs_2.size() && transition_distrs_2.size() > 0) {
        cerr << "error:[NGSSimulator] One fragment end in training data has no reads at the modal length, cannot produce joint samples" << endl;
        exit(1);
    }
    
    finalize();
    
#ifdef debug_ngs_sim
    cerr << "finished initializing simulator" << endl;
#endif
}

Alignment NGSSimulator::sample_read() {
    
    
    Alignment aln;
    // sample a quality string based on the trained distribution
    pair<string, vector<bool>> qual_and_masks = sample_read_quality();
    
#ifdef debug_ngs_sim
    cerr << "sampled qualities and N-mask:" << endl;
    cerr << string_quality_short_to_char(qual_and_masks.first) << endl;
    for (bool mask : qual_and_masks.second) {
        cerr << (mask ? "1" : "0");
    }
    cerr << endl;
#endif
    
    assert(qual_and_masks.first.size() == qual_and_masks.second.size());
    
    aln.set_quality(qual_and_masks.first);
    
    // attempt samples until we get one that succeeds without walking
    // off the end of the graph
    while (!aln.has_path()) {
        // This is our offset along the source path, if in use
        size_t offset;
        // And our direction to go along the source path, if in use
        bool is_reverse;
        // And our position in the graph, which we use whether there's a source path or not.
        pos_t pos;
        // And our path (if dealing with source_paths)
        string source_path;
        // Populate them
        sample_start_pos(offset, is_reverse, pos, source_path);
        // align the first end at this position on the source path or graph
        sample_read_internal(aln, offset, is_reverse, pos, source_path);
        
        // make sure we didn't sample sequence from
        if (retry_on_Ns) {
            if (aln.sequence().find('N') != string::npos) {
                aln.clear_path();
                aln.clear_sequence();
            }
        }
    }
    
    // mask out any of the sequence that we sampled to be an 'N'
    apply_N_mask(*aln.mutable_sequence(), qual_and_masks.second);
    
    aln.set_name(get_read_name());
    algorithms::annotate_with_initial_path_positions(graph, aln);
    return aln;
}

pair<Alignment, Alignment> NGSSimulator::sample_read_pair() {
    pair<Alignment, Alignment> aln_pair;
    pair<pair<string, vector<bool>>, pair<string, vector<bool>>> qual_and_mask_pair = sample_read_quality_pair();
    
#ifdef debug_ngs_sim
    cerr << "sampled qualities and N-masks:" << endl;
    cerr << string_quality_short_to_char(qual_and_mask_pair.first.first) << endl;
    for (bool mask : qual_and_mask_pair.first.second) {
        cerr << (mask ? "1" : "0");
    }
    cerr << endl;
    cerr << string_quality_short_to_char(qual_and_mask_pair.second.first) << endl;
    for (bool mask : qual_and_mask_pair.second.second) {
        cerr << (mask ? "1" : "0");
    }
    cerr << endl;
#endif

    
    assert(qual_and_mask_pair.first.first.size() == qual_and_mask_pair.first.second.size());
    assert(qual_and_mask_pair.second.first.size() == qual_and_mask_pair.second.second.size());
    
    aln_pair.first.set_quality(qual_and_mask_pair.first.first);
    aln_pair.second.set_quality(qual_and_mask_pair.second.first);
    
    // reverse the quality string so that it acts like it's reading from the opposite end
    // when we walk forward from the beginning of the first read
    std::reverse(aln_pair.second.mutable_quality()->begin(),
                 aln_pair.second.mutable_quality()->end());
    
    
    while (!aln_pair.first.has_path() || !aln_pair.second.has_path()) {
        int64_t insert_length = (int64_t) round(insert_sampler(prng));
        if (insert_length < (int64_t) transition_distrs_1.size()) {
            // don't make reads where the insert length is shorter than one end of the read
            continue;
        }
        
        // This is our offset along the source path, if in use
        size_t offset;
        // And our direction to go along the source path, if in use
        bool is_reverse;
        // And our position in the graph, which we use whether there's a source path or not.
        pos_t pos;
        // And our path (if dealing with source_paths)
        string source_path;
        // Populate them
        sample_start_pos(offset, is_reverse, pos, source_path);
        // align the first end at this position on the source path or graph
        sample_read_internal(aln_pair.first, offset, is_reverse, pos, source_path);
        
        if (retry_on_Ns) {
            if (aln_pair.first.sequence().find('N') != string::npos) {
                aln_pair.first.clear_path();
                aln_pair.first.clear_sequence();
            }
        }
        
        if (!aln_pair.first.has_path()) {
            continue;
        }
        
        // walk out the unsequenced part of the insert in the graph
        int64_t remaining_length = insert_length - 2 * transition_distrs_1.size();
        if (remaining_length >= 0) {
            // we need to move forward from the end of the first read
            if (advance_by_distance(offset, is_reverse, pos, remaining_length, source_path)) {
                // we hit the end of the graph trying to walk
                continue;
            }
        } else {
            // we need to walk backwards from the end of the first read
            pos = walk_backwards(aln_pair.first.path(), -remaining_length);
            // Make sure to update the offset along the path as well. If offset
            // and is_reverse aren't being used (becasue we aren't in path
            // mode), the result won't be used either.
            offset += is_reverse ? -remaining_length : remaining_length;
        }
        // guard against running off the end of nodes
        // XXX this should not be happening
        // it seems to occur in some graphs due to the behavior of advance_by_distance
        if (vg::offset(pos) >= graph.get_length(graph.get_handle(id(pos)))) continue;

        // align the second end starting at the walked position
        sample_read_internal(aln_pair.second, offset, is_reverse, pos, source_path); 
        
        if (retry_on_Ns) {
            if (aln_pair.second.sequence().find('N') != string::npos) {
                aln_pair.second.clear_path();
                aln_pair.second.clear_sequence();
            }
        }
    }
    
    // unreverse the second read in the pair
    aln_pair.second = reverse_complement_alignment(aln_pair.second, [&](id_t node_id) {
            return graph.get_length(graph.get_handle(node_id));
    });
    
    // mask out any of the sequence that we sampled to be an 'N'
    apply_N_mask(*aln_pair.first.mutable_sequence(), qual_and_mask_pair.first.second);
    apply_N_mask(*aln_pair.second.mutable_sequence(), qual_and_mask_pair.second.second);
    
    string name = get_read_name();
    aln_pair.first.set_name(name + "_1");
    aln_pair.second.set_name(name + "_2");
    algorithms::annotate_with_initial_path_positions(graph, aln_pair.first);
    algorithms::annotate_with_initial_path_positions(graph, aln_pair.second);
    return aln_pair;
}

void NGSSimulator::sample_read_internal(Alignment& aln, size_t& offset, bool& is_reverse, pos_t& curr_pos,
                                        const string& source_path) {
   
    // Make sure we are starting inside the node
    // XXX this is broken
    auto first_node_length = graph.get_length(graph.get_handle(id(curr_pos)));
    if (vg::offset(curr_pos) >= first_node_length) {
        cerr << "somithgn wrong << " << vg::offset(curr_pos) << " " << first_node_length << endl;
        cerr << vg::id(curr_pos) << ":" << vg::is_rev(curr_pos) << ":" << vg::offset(curr_pos) << endl;
    }
    assert(vg::offset(curr_pos) < first_node_length);
   
    aln.clear_path();
    aln.clear_sequence();
    
    char graph_char = graph.get_base(graph.get_handle(id(curr_pos), is_rev(curr_pos)), vg::offset(curr_pos));
    bool hit_end = false;
    
    // walk a path and generate a read sequence at the same time
    while (aln.sequence().size() < aln.quality().size() && !hit_end) {
        // sample insertion in the true graph path
        while (aln.sequence().size() < aln.quality().size() && prob_sampler(prng) < indel_poly_rate * 0.5) {
            // TODO: no allowance for indel errors on inserted sequence
            
#ifdef debug_ngs_sim
            cerr << "insertion polymorphism at read idx " << aln.sequence().size() << ", graph pos " << curr_pos << endl;
#endif
            
            apply_insertion(aln, curr_pos);
        }
        if (aln.sequence().size() >= aln.quality().size() || hit_end) {
            break;
        }
        
        // sample errors
        double err_sample = prob_sampler(prng);
        double err_prob = phred_prob[aln.quality()[aln.sequence().size()]];
        while (err_sample < err_prob * indel_error_prop && !hit_end) {
            // indel errors
            if (prob_sampler(prng) < 0.5) {
#ifdef debug_ngs_sim
                cerr << "insertion error at read idx " << aln.sequence().size() << ", graph pos " << curr_pos << endl;
#endif
                // insert error
                apply_insertion(aln, curr_pos);
                
                if (aln.sequence().size() >= aln.quality().size() || hit_end) {
                    break;
                }
            }
            else {
#ifdef debug_ngs_sim
                cerr << "deletion error at read idx " << aln.sequence().size() << ", graph pos " << curr_pos << endl;
#endif
                // deletion error
                apply_deletion(aln, curr_pos);
                hit_end = advance(offset, is_reverse, curr_pos, graph_char, source_path);
            }
            
            err_sample = prob_sampler(prng);
            err_prob = phred_prob[aln.quality()[aln.sequence().size()]];
        }
        if (aln.sequence().size() >= aln.quality().size() || hit_end) {
            break;
        }
        
        // get the true graph char, possibly with a substitution polymorphism
        char poly_graph_char = graph_char;
        if (prob_sampler(prng) < sub_poly_rate) {
            poly_graph_char = mutation_alphabets[poly_graph_char != 'N' ? poly_graph_char : alphabet[background_sampler(prng)]][mut_sampler(prng)];
        }
        
        // by default the read matches the true graph char
        char read_char = poly_graph_char;
        
        // sample substitution errors with the remaining err sample
        if (err_sample < err_prob) {
            // substitution error
            read_char = mutation_alphabets[read_char != 'N' ? read_char : alphabet[background_sampler(prng)]][mut_sampler(prng)];
        }
        
#ifdef debug_ngs_sim
        cerr << "aligned base at read idx " << aln.sequence().size() << ", graph pos " << curr_pos << endl;
#endif
        
        // add an aligned base (allowing errors to mask polymorphisms)
        apply_aligned_base(aln, curr_pos, graph_char, read_char);
        hit_end = advance(offset, is_reverse, curr_pos, graph_char, source_path);
        
        if (aln.sequence().size() >= aln.quality().size() || hit_end) {
            break;
        }
        
        // sample deletions in the true graph path
        while (prob_sampler(prng) < indel_poly_rate * 0.5 && !hit_end) {
#ifdef debug_ngs_sim
            cerr << "deletion polymorphism at read idx " << aln.sequence().size() << ", graph pos " << curr_pos << endl;
#endif
            
            apply_deletion(aln, curr_pos);
            hit_end = advance(offset, is_reverse, curr_pos, graph_char, source_path);
        }
    }
    
    // remove the sequence and path if we hit the end the graph before finishing
    // the alignment
    if (aln.sequence().size() != aln.quality().size()) {
        aln.clear_path();
        aln.clear_sequence();
    }
    
#ifdef debug_ngs_sim
    cerr << "completed read: " << pb2json(aln) << endl;
#endif
}

bool NGSSimulator::advance(size_t& offset, bool& is_reverse, pos_t& pos, char& graph_char, const string& source_path) {
    if (source_path.empty()) {
        return advance_on_graph(pos, graph_char);
    } else {
        return advance_on_path(offset, is_reverse, pos, graph_char, source_path);
    }
}

bool NGSSimulator::advance_on_graph(pos_t& pos, char& graph_char) {
    
    // choose a next position at random
    map<pos_t, char> next_pos_chars = algorithms::next_pos_chars(graph, pos);
    if (next_pos_chars.empty()) {
        return true;
    }
    
    vg::uniform_int_distribution<size_t> pos_distr(0, next_pos_chars.size() - 1);
    size_t next = pos_distr(prng);
    auto iter = next_pos_chars.begin();
    for (size_t i = 0; i != next; i++) {
        iter++;
    }
    pos = iter->first;
    graph_char = iter->second;
    
    return false;
}

bool NGSSimulator::advance_on_path(size_t& offset, bool& is_reverse, pos_t& pos, char& graph_char, const string& source_path) {
    size_t path_length = graph.get_path_length(graph.get_path_handle(source_path));
    if (is_reverse) {
        // Go left on the path
        if (offset == 0) {
            // We hit the end
            return true;
        }
        offset--;
    } else {
        // Go right on the path
        if (offset == path_length - 1) {
            // We hit the end
            return true;
        }
        offset++;
    }
    
    // Set position according to position on path
    pos = position_at(&graph, source_path, offset, is_reverse);
    
    // And look up the character
    graph_char = graph.get_base(graph.get_handle(id(pos), is_rev(pos)), vg::offset(pos));
    
    return false;
}

bool NGSSimulator::advance_by_distance(size_t& offset, bool& is_reverse, pos_t& pos, size_t distance,
                                       const string& source_path) {
    if (source_path.empty()) {
        return advance_on_graph_by_distance(pos, distance);
    } else {
        return advance_on_path_by_distance(offset, is_reverse, pos, distance, source_path);
    }
}


bool NGSSimulator::advance_on_graph_by_distance(pos_t& pos, size_t distance) {
    int64_t remaining = distance;
    handle_t handle = graph.get_handle(id(pos), is_rev(pos));
    int64_t node_length = graph.get_length(handle) - offset(pos);
    while (remaining >= node_length) {
        remaining -= node_length;
        vector<handle_t> nexts;
        graph.follow_edges(handle, false, [&](const handle_t& next) {
                nexts.push_back(next);
            });
        if (nexts.empty()) {
            return true;
        }
        size_t choice = vg::uniform_int_distribution<size_t>(0, nexts.size() - 1)(prng);
        handle = nexts[choice];
        node_length = graph.get_length(handle);
    }

    get_id(pos) = graph.get_id(handle);
    get_is_rev(pos) = graph.get_is_reverse(handle);
    get_offset(pos) += remaining;
    
    return false;
}

bool NGSSimulator::advance_on_path_by_distance(size_t& offset, bool& is_reverse, pos_t& pos, size_t distance,
                                               const string& source_path) {
    
    size_t path_length = graph.get_path_length(graph.get_path_handle(source_path));
    if (is_reverse) {
        // Go left on the path
        if (offset < distance) {
            // We hit the end
            return true;
        }
        offset -= distance;
    } else {
        // Go right on the path
        if (offset + distance >= path_length) {
            // We hit the end
            return true;
        }
        offset += distance;
    }
    
    // Set position according to position on path
    pos = position_at(&graph, source_path, offset, is_reverse);
    
    return false;
}

pos_t NGSSimulator::walk_backwards(const Path& path, size_t distance) {
    // Starting at the past-the-end of the path, walk back to the given nonzero distance.
    // Walking back the whole path length puts you at the start of the path.
    
    if (distance > path_to_length(path)) {
        throw runtime_error("Cannot walk back " + to_string(distance) + " on path of length " + to_string(path_to_length(path)));
    }
    assert(distance > 0);

    // walk backwards until we find the mapping it's on
    int64_t remaining = distance;
    int64_t mapping_idx = path.mapping_size() - 1;
    int64_t mapping_length = mapping_to_length(path.mapping(mapping_idx));
    while (remaining > mapping_length) {
        remaining -= mapping_length;
        mapping_idx--;
        mapping_length = mapping_to_length(path.mapping(mapping_idx));
    }
    // Now we know the position we want is inside this mapping.
    const Mapping& mapping = path.mapping(mapping_idx);
    const Position& mapping_pos = mapping.position();
    // walk forward from the beginning of the mapping, edit by edit, until we've passed where it is on the read
    int64_t remaining_flipped = mapping_length - remaining;
    int64_t edit_idx = 0;
    int64_t prefix_from_length = 0;
    int64_t prefix_to_length = 0;
    while (prefix_to_length <= remaining_flipped) {
        prefix_from_length += mapping.edit(edit_idx).from_length();
        prefix_to_length += mapping.edit(edit_idx).to_length();
        edit_idx++;
    }
    // Go back one edit to the edit that covers the position we want to be at
    --edit_idx;
    // use the mapping's position and the distance we traveled on the graph to get the offset of
    // the beginning of this edit
    int64_t offset = mapping_pos.offset() + prefix_from_length - mapping.edit(edit_idx).from_length();
    if (mapping.edit(edit_idx).from_length() == mapping.edit(edit_idx).to_length()) {
        // if the edit is a match/mismatch, we can walk part of it
        
        // How many extra bases did this mapping have when we got past where we wanted to be?
        auto extra_bases = prefix_to_length - remaining_flipped;
        
        // Take all the non-extra bases of the mapping.
        offset += mapping.edit(edit_idx).to_length() - extra_bases;
    } else {
        // Otherwise it's an insert, so just land at the spot it is inserted at.
        
        // But we will have a problem if the insert is the last thing in its mapping
        if (prefix_from_length == mapping_from_length(mapping)) {
            // We are trying to put our starting position past the end of this
            // mapping, because we are landing inside an insert that is the
            // last thing in its mapping.
            
            // Note that this can happen at the end of the path, if the whole
            // path ends in an insert. So we can't just go to the next
            // position.
            
            // There's not really a quite correct thing to do here, but since
            // we can't go right we have to go left.
            assert(offset > 0);
            offset--;
        }
    }
    
    // Get the length of the node we landed on
    auto node_length = graph.get_length(graph.get_handle(mapping_pos.node_id()));
    // The position we pick should not be past the end of the node.
    if (offset >= node_length) {
        cerr << pb2json(path) << endl;
        cerr << "Covering mapping: " << pb2json(mapping) << endl;
        cerr << "Covering edit: " << pb2json(mapping.edit(edit_idx)) << endl;
        throw runtime_error("Could not go back " + to_string(distance) + " in path of length " +
            to_string(path_to_length(path)) + "; hit node " + to_string(mapping_pos.node_id()) + " length " +
            to_string(node_length) + " end at offset " + to_string(offset));
    }

    return make_pos_t(mapping_pos.node_id(), mapping_pos.is_reverse(), offset);
}

void NGSSimulator::apply_aligned_base(Alignment& aln, const pos_t& pos, char graph_char,
                                      char read_char) {
    Path* path = aln.mutable_path();
    aln.mutable_sequence()->push_back(read_char);
    bool is_match = (graph_char == read_char);
    
    if (path->mapping_size() == 0) {
        Mapping* new_mapping = path->add_mapping();
        new_mapping->set_rank(1);
        
        Position* mapping_pos = new_mapping->mutable_position();
        mapping_pos->set_node_id(id(pos));
        mapping_pos->set_is_reverse(is_rev(pos));
        mapping_pos->set_offset(offset(pos));
        
        Edit* new_edit = new_mapping->add_edit();
        new_edit->set_from_length(1);
        new_edit->set_to_length(1);
        if (!is_match) {
            new_edit->mutable_sequence()->push_back(read_char);
        }
    }
    else {
        Mapping* last_mapping = path->mutable_mapping(path->mapping_size() - 1);
        if (last_mapping->position().node_id() == id(pos) &&
            last_mapping->position().is_reverse() == is_rev(pos)) {
            
            Edit* last_edit = last_mapping->mutable_edit(last_mapping->edit_size() - 1);
            if (last_edit->from_length() > 0 && last_edit->to_length() > 0) {
                if (last_edit->sequence().size() > 0 && !is_match) {
                    last_edit->set_from_length(last_edit->from_length() + 1);
                    last_edit->set_to_length(last_edit->to_length() + 1);
                    last_edit->mutable_sequence()->push_back(read_char);
                }
                else if (last_edit->sequence().size() == 0 && is_match) {
                    last_edit->set_from_length(last_edit->from_length() + 1);
                    last_edit->set_to_length(last_edit->to_length() + 1);
                }
                else {
                    Edit* new_edit = last_mapping->add_edit();
                    new_edit->set_from_length(1);
                    new_edit->set_to_length(1);
                    if (!is_match) {
                        new_edit->mutable_sequence()->push_back(read_char);
                    }
                }
            }
            else {
                Edit* new_edit = last_mapping->add_edit();
                new_edit->set_from_length(1);
                new_edit->set_to_length(1);
                if (!is_match) {
                    new_edit->mutable_sequence()->push_back(read_char);
                }
            }
        }
        else {
            Mapping* new_mapping = path->add_mapping();
            new_mapping->set_rank(last_mapping->rank() + 1);
            
            Position* mapping_pos = new_mapping->mutable_position();
            mapping_pos->set_node_id(id(pos));
            mapping_pos->set_is_reverse(is_rev(pos));
            
            Edit* new_edit = new_mapping->add_edit();
            new_edit->set_from_length(1);
            new_edit->set_to_length(1);
            if (!is_match) {
                new_edit->mutable_sequence()->push_back(read_char);
            }
        }
    }
}

void NGSSimulator::apply_deletion(Alignment& aln, const pos_t& pos) {
    Path* path = aln.mutable_path();
    if (path->mapping_size() == 0) {
        // don't introduce a deletion at the beginning of a read
        // TODO: check for deletions at the end of a read?
        return;
    }
    
    Mapping* last_mapping = path->mutable_mapping(path->mapping_size() - 1);
    if (last_mapping->position().node_id() == id(pos) &&
        last_mapping->position().is_reverse() == is_rev(pos)) {
        
        Edit* last_edit = last_mapping->mutable_edit(last_mapping->edit_size() - 1);
        if (last_edit->from_length() > 0 && last_edit->to_length() == 0) {
            last_edit->set_from_length(last_edit->from_length() + 1);
        }
        else {
            Edit* new_edit = last_mapping->add_edit();
            new_edit->set_from_length(1);
        }
    }
    else {
        Mapping* new_mapping = path->add_mapping();
        new_mapping->set_rank(last_mapping->rank() + 1);
        
        Position* mapping_pos = new_mapping->mutable_position();
        mapping_pos->set_node_id(id(pos));
        mapping_pos->set_is_reverse(is_rev(pos));
        
        Edit* new_edit = new_mapping->add_edit();
        new_edit->set_from_length(1);
    }
}

void NGSSimulator::apply_insertion(Alignment& aln, const pos_t& pos) {
    Path* path = aln.mutable_path();
    char insert_char = alphabet[background_sampler(prng)];
    aln.mutable_sequence()->push_back(insert_char);
    
    if (path->mapping_size() == 0) {
        Mapping* new_mapping = path->add_mapping();
        new_mapping->set_rank(1);
        
        Position* mapping_pos = new_mapping->mutable_position();
        mapping_pos->set_node_id(id(pos));
        mapping_pos->set_is_reverse(is_rev(pos));
        mapping_pos->set_offset(offset(pos));
        
        Edit* new_edit = new_mapping->add_edit();
        new_edit->set_to_length(1);
        new_edit->set_sequence(string(1, insert_char));
    }
    else {
        Mapping* last_mapping = path->mutable_mapping(path->mapping_size() - 1);
        Edit* last_edit = last_mapping->mutable_edit(last_mapping->edit_size() - 1);
        if (last_edit->from_length() == 0 && last_edit->to_length() > 0) {
            last_edit->set_to_length(last_edit->to_length() + 1);
            last_edit->mutable_sequence()->push_back(insert_char);
        }
        else {
            Edit* new_edit = last_mapping->add_edit();
            new_edit->set_to_length(1);
            new_edit->set_sequence(string(1, insert_char));
        }
    }
}

void NGSSimulator::sample_start_pos(size_t& offset, bool& is_reverse, pos_t& pos, string& source_path) {
    if (source_paths.empty()) {
        pos = sample_start_graph_pos();
        offset = 0;
        is_reverse = false;
        source_path = "";
    } else {
        tie(offset, is_reverse, pos, source_path) = sample_start_path_pos();
    }
}

pos_t NGSSimulator::sample_start_graph_pos() {
    // The start pos sampler has been set up in graph space, 1-based
    assert(start_pos_samplers.size() == 1);
    size_t idx = start_pos_samplers[0](prng);
    
    id_t id = dynamic_cast<VectorizableHandleGraph&>(graph).node_at_vector_offset(idx);
    bool rev = strand_sampler(prng);
    size_t node_offset = idx - dynamic_cast<VectorizableHandleGraph&>(graph).node_vector_offset(id) - 1;
    
    return make_pos_t(id, rev, node_offset);
}

tuple<size_t, bool, pos_t, string> NGSSimulator::sample_start_path_pos() {
    // choose a path
    size_t source_path_idx = path_sampler(prng);
    string source_path = source_paths[source_path_idx];
    // The start pos sampler hasd been set up in path space, 0-based
    size_t offset = start_pos_samplers[source_path_idx](prng);
    bool rev = strand_sampler(prng);
    pos_t pos = position_at(&graph, source_path, offset, rev);
    
    return make_tuple(offset, rev, pos, source_path);
}

string NGSSimulator::get_read_name() {
    stringstream sstrm;
    sstrm << "seed_" << seed << "_fragment_" << sample_counter;
    sample_counter++;
    return sstrm.str();
}

void NGSSimulator::record_read_quality(const Alignment& aln, bool read_2) {
    const string& quality = aln.quality();
    const string& sequence = aln.sequence();
    auto& transition_distrs = read_2 ? transition_distrs_2 : transition_distrs_1;
    if (quality.empty()) {
        return;
    }
    while (transition_distrs.size() < quality.size()) {
        transition_distrs.emplace_back(seed ? seed + transition_distrs.size() + 1 : random_device()());
    }
    // record the initial quality and N-mask
    transition_distrs[0].record_transition(pair<uint8_t, bool>(0, false),
                                           pair<uint8_t, bool>(quality[0], sequence[0] == 'N'));
    // record the subsequent quality and N-mask transitions
    for (size_t i = 1; i < transition_distrs.size(); i++) {
        transition_distrs[i].record_transition(pair<uint8_t, bool>(quality[i - 1], sequence[i - 1] == 'N'),
                                               pair<uint8_t, bool>(quality[i], sequence[i] == 'N'));
    }
}
    
void NGSSimulator::record_read_pair_quality(const Alignment& aln_1, const Alignment& aln_2) {
    // record the transitions within the reads separates
    record_read_quality(aln_1, false);
    record_read_quality(aln_2, true);
    // record the joint distribution of the first quality and N-mask
    if (!aln_1.quality().empty() && !aln_2.quality().empty()) {
        joint_initial_distr.record_transition(pair<uint8_t, bool>(0, false),
                                              make_pair(pair<uint8_t, bool>(aln_1.quality()[0], aln_1.sequence()[0] == 'N'),
                                                        pair<uint8_t, bool>(aln_2.quality()[0], aln_2.sequence()[0] == 'N')));
    }
}

void NGSSimulator::finalize() {
    for (auto& markov_distr : transition_distrs_1) {
        markov_distr.finalize();
    }
    for (auto& markov_distr : transition_distrs_2) {
        markov_distr.finalize();
    }
    joint_initial_distr.finalize();
}

pair<string, vector<bool>> NGSSimulator::sample_read_quality() {
    // only use the first trained distribution (on the assumption that it better reflects the properties of
    // single-ended sequencing)
    return sample_read_quality_internal(transition_distrs_1[0].sample_transition(pair<uint8_t, bool>(0, false)),
                                        true);
}
    
pair<pair<string, vector<bool>>, pair<string, vector<bool>>> NGSSimulator::sample_read_quality_pair() {
    if (transition_distrs_2.empty()) {
        // no paired training data, sample qual strings independently
        return make_pair(sample_read_quality(), sample_read_quality());
    }
    else {
        // paired training data, sample the start quality jointly
        auto first_quals_and_masks = joint_initial_distr.sample_transition(pair<uint8_t, bool>(0, false));
        return make_pair(sample_read_quality_internal(first_quals_and_masks.first, true),
                         sample_read_quality_internal(first_quals_and_masks.second, false));
    }
}
    
                                
pair<string, vector<bool>> NGSSimulator::sample_read_quality_internal(pair<uint8_t, bool> first,
                                                                      bool transitions_1) {
    
    auto& transition_distrs = transitions_1 ? transition_distrs_1 : transition_distrs_2;
    string quality(transition_distrs.size(), first.first);
    vector<bool> n_masks(transition_distrs.size(), first.second);
    pair<uint8_t, bool> at = first;
    for (size_t i = 1; i < transition_distrs.size(); i++) {
        at = transition_distrs[i].sample_transition(at);
        quality[i] = at.first;
        n_masks[i] = at.second;
    }
    return make_pair(quality, n_masks);
}
                                              
void NGSSimulator::apply_N_mask(string& sequence, const vector<bool>& n_mask) {
    assert(sequence.size() == n_mask.size());
    for (size_t i = 0; i < n_mask.size(); i++) {
        if (n_mask[i]) {
            sequence[i] = 'N';
        }
    }
}
    
template<class From, class To>
NGSSimulator::MarkovDistribution<From, To>::MarkovDistribution(size_t seed) : prng(seed) {
    // nothing to do
}

template<class From, class To>
void NGSSimulator::MarkovDistribution<From, To>::record_transition(From from, To to) {
    if (!cond_distrs.count(from)) {
        cond_distrs[from] = vector<size_t>(value_at.size(), 0);
    }
    
    if (!column_of.count(to)) {
        column_of[to] = value_at.size();
        value_at.push_back(to);
        for (pair<const From, vector<size_t>>& cond_distr : cond_distrs) {
            cond_distr.second.push_back(0);
        }
    }
    
    cond_distrs[from][column_of[to]]++;
}

template<class From, class To>
void NGSSimulator::MarkovDistribution<From, To>::finalize() {
    for (pair<const From, vector<size_t>>& cond_distr : cond_distrs) {
        for (size_t i = 1; i < cond_distr.second.size(); i++) {
            cond_distr.second[i] += cond_distr.second[i - 1];
        }
        
        samplers[cond_distr.first] = vg::uniform_int_distribution<size_t>(1, cond_distr.second.back());
    }
}

template<class From, class To>
To NGSSimulator::MarkovDistribution<From, To>::sample_transition(From from) {
    // return randomly if a transition has never been observed
    if (!cond_distrs.count(from)) {
        return value_at[vg::uniform_int_distribution<size_t>(0, value_at.size() - 1)(prng)];
    }
    
    size_t sample_val = samplers[from](prng);
    vector<size_t>& cdf = cond_distrs[from];
    
    if (sample_val <= cdf[0]) {
        return value_at[0];
    }
    
    size_t low = 0;
    size_t hi = cdf.size() - 1;
    while (hi > low + 1) {
        int64_t mid = (hi + low) / 2;
        
        if (sample_val <= cdf[mid]) {
            hi = mid;
        }
        else {
            low = mid;
        }
    }
    return value_at[hi];
}

}
