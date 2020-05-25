#include "aligner.hpp"

#include "hash_map.hpp"

//#define debug_print_score_matrices

using namespace vg;
using namespace std;


static const double quality_scale_factor = 10.0 / log(10.0);
static const double exp_overflow_limit = log(std::numeric_limits<double>::max());

GSSWAligner::~GSSWAligner(void) {
    free(nt_table);
    free(score_matrix);
}

gssw_graph* GSSWAligner::create_gssw_graph(const HandleGraph& g) const {
    
    vector<handle_t> topological_order = algorithms::lazier_topological_order(&g);
    
    gssw_graph* graph = gssw_graph_create(g.get_node_count());
    unordered_map<int64_t, gssw_node*> nodes;
    
    // compute the topological order
    for (const handle_t& handle : topological_order) {
        auto cleaned_seq = nonATGCNtoN(g.get_sequence(handle));
        gssw_node* node = gssw_node_create(nullptr,       // TODO: the ID should be enough, don't need Node* too
                                           g.get_id(handle),
                                           cleaned_seq.c_str(),
                                           nt_table,
                                           score_matrix); // TODO: this arg isn't used, could edit
                                                          // in gssw
        nodes[node->id] = node;
        gssw_graph_add_node(graph, node);
    }
    
    g.for_each_edge([&](const edge_t& edge) {
        if(!g.get_is_reverse(edge.first) && !g.get_is_reverse(edge.second)) {
            // This is a normal end to start edge.
            gssw_nodes_add_edge(nodes[g.get_id(edge.first)], nodes[g.get_id(edge.second)]);
        }
        else if (g.get_is_reverse(edge.first) && g.get_is_reverse(edge.second)) {
            // This is a start to end edge, but isn't reversing and can be converted to a normal end to start edge.
            
            // Flip the start and end
            gssw_nodes_add_edge(nodes[g.get_id(edge.second)], nodes[g.get_id(edge.first)]);
        }
        else {
            // TODO: It's a reversing edge, which gssw doesn't support yet. What
            // we should really do is do a topological sort to break cycles, and
            // then flip everything at the lower-rank end of this edge around,
            // so we don't have to deal with its reversing-ness. But for now we
            // just die so we don't get nonsense into gssw.
#pragma omp critical
            {
                // We need the critical section so we don't throw uncaught
                // exceptions in multiple threads at once, leading to C++ trying
                // to run termiante in parallel. This doesn't make it safe, just
                // slightly safer.
                cerr << "Can't gssw over reversing edge " << g.get_id(edge.first) << (g.get_is_reverse(edge.first) ? "-" : "+") << " -> " << g.get_id(edge.second) << (g.get_is_reverse(edge.second) ? "-" : "+") << endl;
                // TODO: there's no safe way to kill the program without a way
                // to signal the master to do it, via a shared variable in the
                // clause that made us parallel.
            }
            exit(1);
        }
        return true;
    });
    
    return graph;
    
}

unordered_set<vg::id_t> GSSWAligner::identify_pinning_points(const HandleGraph& graph) const {
    
    unordered_set<vg::id_t> return_val;
    
    // start at the sink nodes
    vector<handle_t> sinks = algorithms::tail_nodes(&graph);
    
    // walk backwards to find non-empty nodes if necessary
    for (const handle_t& handle : sinks) {
        vector<handle_t> stack(1, handle);
        while (!stack.empty()) {
            handle_t here =  stack.back();
            stack.pop_back();
            
            if (graph.get_length(here) > 0) {
                return_val.insert(graph.get_id(here));
            }
            else {
                graph.follow_edges(here, true, [&](const handle_t& prev) {
                    // TODO: technically this won't filter out all redundant walks, but it should
                    // handle all cases we're practically interested in and it doesn't require a
                    // second set object
                    if (!return_val.count(graph.get_id(prev))) {
                        stack.push_back(prev);
                    }
                });
            }
        }
    }
    
    return return_val;
}

void GSSWAligner::gssw_mapping_to_alignment(gssw_graph* graph,
                                            gssw_graph_mapping* gm,
                                            Alignment& alignment,
                                            bool pinned,
                                            bool pin_left) const {
    alignment.clear_path();
    alignment.set_score(gm->score);
    alignment.set_query_position(0);
    Path* path = alignment.mutable_path();
    //alignment.set_cigar(graph_cigar(gm));
    
    gssw_graph_cigar* gc = &gm->cigar;
    gssw_node_cigar* ncs = gc->elements;
    //cerr << "gm->position " << gm->position << endl;
    string& to_seq = *alignment.mutable_sequence();
    //cerr << "-------------" << endl;
    
#ifdef debug_print_score_matrices
    gssw_graph_print_score_matrices(graph, to_seq.c_str(), to_seq.size(), stderr);
#endif
    
    int to_pos = 0;
    int from_pos = gm->position;
    
    for (int i = 0; i < gc->length; ++i) {
        // check that the current alignment has a non-zero length
        gssw_cigar* c = ncs[i].cigar;
        int l = c->length;
        if (l == 0) continue;
        gssw_cigar_element* e = c->elements;
        
        gssw_node* node = ncs[i].node;
        Mapping* mapping = path->add_mapping();
        
        if (i > 0) {
            // reset for each node after the first
            from_pos = 0;
        }
        
        mapping->mutable_position()->set_node_id(node->id);
        mapping->mutable_position()->set_offset(from_pos);
        mapping->set_rank(path->mapping_size());
        
        //cerr << node->id << ":" << endl;
        
        for (int j=0; j < l; ++j, ++e) {
            int32_t length = e->length;
            //cerr << e->length << e->type << endl;
            
            Edit* edit;
            switch (e->type) {
                case 'M':
                case 'X':
                case 'N': {
                    //cerr << "j = " << j << ", type = " << e->type << endl;
                    // do the sequences match?
                    // emit a stream of "SNPs" and matches
                    int h = from_pos;
                    int last_start = from_pos;
                    int k = to_pos;
                    for ( ; h < from_pos + length; ++h, ++k) {
                        //cerr << h << ":" << k << " " << node->seq[h] << " " << to_seq[k] << endl;
                        if (node->seq[h] != to_seq[k]) {
                            // emit the last "match" region
                            if (h - last_start > 0) {
                                edit = mapping->add_edit();
                                edit->set_from_length(h-last_start);
                                edit->set_to_length(h-last_start);
                            }
                            // set up the SNP
                            edit = mapping->add_edit();
                            edit->set_from_length(1);
                            edit->set_to_length(1);
                            edit->set_sequence(to_seq.substr(k,1));
                            last_start = h+1;
                        }
                    }
                    // handles the match at the end or the case of no SNP
                    if (h - last_start > 0) {
                        edit = mapping->add_edit();
                        edit->set_from_length(h-last_start);
                        edit->set_to_length(h-last_start);
                    }
                    to_pos += length;
                    from_pos += length;
                } break;
                case 'D':
                    edit = mapping->add_edit();
                    edit->set_from_length(length);
                    edit->set_to_length(0);
                    from_pos += length;
                    break;
                case 'I':
                    edit = mapping->add_edit();
                    edit->set_from_length(0);
                    edit->set_to_length(length);
                    edit->set_sequence(to_seq.substr(to_pos, length));
                    to_pos += length;
                    break;
                case 'S':
                    // note that soft clips and insertions are semantically equivalent
                    // and can only be differentiated by their position in the read
                    // with soft clips coming at the start or end
                    edit = mapping->add_edit();
                    edit->set_from_length(0);
                    edit->set_to_length(length);
                    edit->set_sequence(to_seq.substr(to_pos, length));
                    to_pos += length;
                    break;
                default:
                    cerr << "error:[Aligner::gssw_mapping_to_alignment] "
                    << "unsupported cigar op type " << e->type << endl;
                    exit(1);
                    break;
                    
            }
        }
    }
    
    // compute and set identity
    alignment.set_identity(identity(alignment.path()));
}

void GSSWAligner::unreverse_graph(gssw_graph* graph) const {
    // this is only for getting correct reference-relative edits, so we can get away with only
    // reversing the sequences and not paying attention to the edges
    
    for (size_t i = 0; i < graph->size; i++) {
        gssw_node* node = graph->nodes[i];
        for (int j = 0, stop = node->len / 2; j < stop; j++) {
            std::swap(node->seq[j], node->seq[node->len - j - 1]);
        }
    }
}

void GSSWAligner::unreverse_graph_mapping(gssw_graph_mapping* gm) const {
    
    gssw_graph_cigar* graph_cigar = &(gm->cigar);
    gssw_node_cigar* node_cigars = graph_cigar->elements;
    
    // reverse the order of the node cigars
    int32_t num_switching_nodes = graph_cigar->length / 2;
    int32_t last_idx = graph_cigar->length - 1;
    for (int32_t i = 0; i < num_switching_nodes; i++) {
        std::swap(node_cigars[i], node_cigars[last_idx - i]);
    }
    
    // reverse the actual cigar string for each node cigar
    for (int32_t i = 0; i < graph_cigar->length; i++) {
        gssw_cigar* node_cigar = node_cigars[i].cigar;
        gssw_cigar_element* elements = node_cigar->elements;
        
        int32_t num_switching_elements = node_cigar->length / 2;
        last_idx = node_cigar->length - 1;
        for (int32_t j = 0; j < num_switching_elements; j++) {
            std::swap(elements[j], elements[last_idx - j]);
        }
    }
    
    // compute the position in the first node
    if (graph_cigar->length > 0) {
        gssw_cigar_element* first_node_elements = node_cigars[0].cigar->elements;
        int32_t num_first_node_elements = node_cigars[0].cigar->length;
        uint32_t num_ref_aligned = 0; // the number of characters on the node sequence that are aligned
        for (int32_t i = 0; i < num_first_node_elements; i++) {
            switch (first_node_elements[i].type) {
                case 'M':
                case 'X':
                case 'N':
                case 'D':
                    num_ref_aligned += first_node_elements[i].length;
                    break;
                    
            }
        }
        gm->position = node_cigars[0].node->len - num_ref_aligned - (graph_cigar->length == 1 ? gm->position : 0);
    }
    else {
        gm->position = 0;
    }
}

string GSSWAligner::graph_cigar(gssw_graph_mapping* gm) const {

    stringstream s;
    gssw_graph_cigar* gc = &gm->cigar;
    gssw_node_cigar* nc = gc->elements;
    int to_pos = 0;
    int from_pos = gm->position;
    //string& to_seq = *alignment.mutable_sequence();
    s << from_pos << '@';
    for (int i = 0; i < gc->length; ++i, ++nc) {
        if (i > 0) from_pos = 0; // reset for each node after the first
        Node* from_node = (Node*) nc->node->data;
        s << from_node->id() << ':';
        gssw_cigar* c = nc->cigar;
        int l = c->length;
        gssw_cigar_element* e = c->elements;
        for (int j=0; j < l; ++j, ++e) {
            s << e->length << e->type;
        }
        if (i + 1 < gc->length) {
            s << ",";
        }
    }
    return s.str();
}

void GSSWAligner::init_mapping_quality(const int8_t* score_matrix, double gc_content) {
    
    // TODO: repetitive with QualAdjAligner constructor
    
    // convert gc content into base-wise frequencies
    double* nt_freqs = (double*) malloc(sizeof(double) * 4);
    nt_freqs[0] = 0.5 * (1 - gc_content);
    nt_freqs[1] = 0.5 * gc_content;
    nt_freqs[2] = 0.5 * gc_content;
    nt_freqs[3] = 0.5 * (1 - gc_content);
    
    log_base = gssw_recover_log_base(score_matrix, nt_freqs, 4, 1e-12);
    
    free(nt_freqs);
}

int32_t GSSWAligner::score_gap(size_t gap_length) const {
    return gap_length ? -gap_open - (gap_length - 1) * gap_extension : 0;
}

double GSSWAligner::maximum_mapping_quality_exact(const vector<double>& scaled_scores, size_t* max_idx_out,
                                                  const vector<double>* multiplicities) {
    
    // work in log transformed values to avoid risk of overflow
    double log_sum_exp = numeric_limits<double>::lowest();
    double max_score = numeric_limits<double>::lowest();
    
    // go in reverse order because this has fewer numerical problems when the scores are sorted (as usual)
    for (int64_t i = scaled_scores.size() - 1; i >= 0; i--) {
        // get the value of one copy of the score and check if it's the max
        double score = scaled_scores.at(i);
        if (score >= max_score) {
            // Since we are going in reverse order, make sure to break ties in favor of the earlier item.
            *max_idx_out = i;
            max_score = score;
        }
        
        // add all copies of the score
        if (multiplicities && multiplicities->at(i) > 1.0) {
            score += log(multiplicities->at(i));
        }
        
        // accumulate the sum of all score
        log_sum_exp = add_log(log_sum_exp, score);
    }
    
    // if necessary, assume a null alignment of 0.0 for comparison since this is local
    if (scaled_scores.size() == 1) {
        if (multiplicities && multiplicities->at(0) <= 1.0) {
            log_sum_exp = add_log(log_sum_exp, 0.0);
        }
        else if (!multiplicities) {
            log_sum_exp = add_log(log_sum_exp, 0.0);
        }
    }
    
    double direct_mapq = -quality_scale_factor * subtract_log(0.0, max_score - log_sum_exp);
    return std::isinf(direct_mapq) ? (double) numeric_limits<int32_t>::max() : direct_mapq;
}

// TODO: this algorithm has numerical problems that would be difficult to solve without increasing the
// time complexity: adding the probability of the maximum likelihood tends to erase the contribution
// of the other terms so that when you subtract them off you get scores of 0 or infinity

//vector<double> Aligner::all_mapping_qualities_exact(vector<double>& scaled_scores) {
//
//    double max_score = *max_element(scaled_scores.begin(), scaled_scores.end());
//    size_t size = scaled_scores.size();
//
//    vector<double> mapping_qualities(size);
//
//    if (max_score * size < exp_overflow_limit) {
//        // no risk of double overflow, sum exp directly (half as many transcendental function evals)
//        vector<double> exp_scaled_scores(size);
//        for (size_t i = 0; i < size; i++) {
//            exp_scaled_scores[i] = exp(scaled_scores[i]);
//        }
//        double denom = std::accumulate(exp_scaled_scores.begin(), exp_scaled_scores.end(), 0.0);
//        for (size_t i = 0; i < size; i++) {
//            mapping_qualities[i] = -10.0 * log10((denom - exp_scaled_scores[i]) / denom);
//        }
//    }
//    else {
//        // work in log transformed valued to avoid risk of overflow
//        double log_sum_exp = scaled_scores[0];
//        for (size_t i = 1; i < size; i++) {
//            log_sum_exp = add_log(log_sum_exp, scaled_scores[i]);
//        }
//        for (size_t i = 0; i < size; i++) {
//            mapping_qualities[i] = -10.0 * log10(1.0 - exp(scaled_scores[i] - log_sum_exp));
//        }
//    }
//    return mapping_qualities;
//}

double GSSWAligner::maximum_mapping_quality_approx(const vector<double>& scaled_scores, size_t* max_idx_out,
                                                   const vector<double>* multiplicities) {
    assert(!scaled_scores.empty());
    
    // determine the maximum score and the count of the next highest score
    double max_score = scaled_scores.at(0);
    size_t max_idx = 0;
    
    // we start with the possibility of a null score of 0.0
    double next_score = 0.0;
    double next_count = 1.0;
    
    if (multiplicities) {
        if (multiplicities->at(0) > 1.0) {
            // there are extra copies of this one, so we'll init with those
            next_score = max_score;
            next_count = multiplicities->at(0) - 1.0;
        }
    }
    
    for (int32_t i = 1; i < scaled_scores.size(); ++i) {
        double score = scaled_scores.at(i);
        if (score > max_score) {
            if (multiplicities && multiplicities->at(i) > 1.0) {
                // there are extra counts of the new highest score due to multiplicity
                next_score = score;
                next_count = multiplicities->at(i) - 1.0;
            }
            else if (next_score == max_score) {
                // the next highest was the same score as the old max, so we can
                // add its count back in
                next_count += 1.0;
            }
            else {
                // the old max score is now the second highest
                next_score = max_score;
                next_count = multiplicities ? multiplicities->at(max_idx) : 1.0;
            }
            max_score = score;
            max_idx = i;
        }
        else if (score > next_score) {
            // the new score is the second highest
            next_score = score;
            next_count = multiplicities ? multiplicities->at(i) : 1.0;
        }
        else if (score == next_score) {
            // the new score ties the second highest, so we combine their counts
            next_count += multiplicities ? multiplicities->at(i) : 1.0;
        }
    }
   
    // record the index of the highest score
    *max_idx_out = max_idx;

    return max(0.0, quality_scale_factor * (max_score - next_score - (next_count > 1.0 ? log(next_count) : 0.0)));
}

double GSSWAligner::group_mapping_quality_exact(const vector<double>& scaled_scores, const vector<size_t>& group,
                                                const vector<double>* multiplicities) const {
    
    // work in log transformed values to avoid risk of overflow
    double total_log_sum_exp = numeric_limits<double>::lowest();
    double non_group_log_sum_exp = numeric_limits<double>::lowest();
    
    // go in reverse order because this has fewer numerical problems when the scores are sorted (as usual)
    int64_t group_idx = group.size() - 1;
    for (int64_t i = scaled_scores.size() - 1; i >= 0; i--) {
        
        // the score of one alignment
        double score = scaled_scores.at(i);
        
        // the score all the multiples of this score combined
        double multiple_score = score;
        if (multiplicities && multiplicities->at(i) > 1.0) {
            multiple_score += log(multiplicities->at(i));
        }
        
        total_log_sum_exp = add_log(total_log_sum_exp, multiple_score);
        
        if (group_idx >= 0 ? i == group[group_idx] : false) {
            // this is the next index in the group
            group_idx--;
            if (multiplicities && multiplicities->at(i) > 1.0) {
                // there's some remaining multiples of this score that don't get added into the group
                 non_group_log_sum_exp = add_log(non_group_log_sum_exp,
                                                 score + log(multiplicities->at(i) - 1.0));
            }
        }
        else {
            // this index is not part of the group
            non_group_log_sum_exp = add_log(non_group_log_sum_exp, multiple_score);
        }
    }
    
    if (scaled_scores.size() == 1) {
        if (multiplicities && multiplicities->at(0) <= 1.0) {
            // assume a null alignment of 0.0 for comparison since this is local
            non_group_log_sum_exp = add_log(non_group_log_sum_exp, 0.0);
            total_log_sum_exp = add_log(total_log_sum_exp, 0.0);
        }
        else if (!multiplicities) {
            //TODO: repetitive, do I need to be this careful to not deref a null?
            // assume a null alignment of 0.0 for comparison since this is local
            non_group_log_sum_exp = add_log(non_group_log_sum_exp, 0.0);
            total_log_sum_exp = add_log(total_log_sum_exp, 0.0);
        }
    }
    
    double direct_mapq = quality_scale_factor * (total_log_sum_exp - non_group_log_sum_exp);
    return (std::isinf(direct_mapq) || direct_mapq > numeric_limits<int32_t>::max()) ?
           (double) numeric_limits<int32_t>::max() : direct_mapq;
}

void GSSWAligner::compute_mapping_quality(vector<Alignment>& alignments,
                                          int max_mapping_quality,
                                          bool fast_approximation,
                                          double cluster_mq,
                                          bool use_cluster_mq,
                                          int overlap_count,
                                          double mq_estimate,
                                          double maybe_mq_threshold,
                                          double identity_weight) const {
    
    if (log_base <= 0.0) {
        cerr << "error:[Aligner] must call init_mapping_quality before computing mapping qualities" << endl;
        exit(EXIT_FAILURE);
    }

    if (alignments.empty()) {
        return;
    }
    
    vector<double> scaled_scores(alignments.size());
    for (size_t i = 0; i < alignments.size(); i++) {
        scaled_scores[i] = log_base * alignments[i].score();
    }

    double mapping_quality;
    size_t max_idx;
    if (!fast_approximation) {
        mapping_quality = maximum_mapping_quality_exact(scaled_scores, &max_idx);
    }
    else {
        mapping_quality = maximum_mapping_quality_approx(scaled_scores, &max_idx);
    }

    if (use_cluster_mq) {
        mapping_quality = prob_to_phred(sqrt(phred_to_prob(cluster_mq + mapping_quality)));
    }

    if (overlap_count) {
        mapping_quality -= quality_scale_factor * log(overlap_count);
    }

    auto& max_aln = alignments.at(max_idx);
    int l = max(alignment_to_length(max_aln), alignment_from_length(max_aln));
    double identity = 1. - (double)(l * match - max_aln.score()) / (match + mismatch) / l;

    mapping_quality /= 2;

    mapping_quality *= pow(identity, identity_weight);

    if (mq_estimate < maybe_mq_threshold && mq_estimate < mapping_quality) {
        mapping_quality = prob_to_phred(sqrt(phred_to_prob(mq_estimate + mapping_quality)));
    }

    if (mapping_quality > max_mapping_quality) {
        mapping_quality = max_mapping_quality;
    }

    if (alignments[max_idx].score() == 0) {
        mapping_quality = 0;
    }

    alignments[max_idx].set_mapping_quality(max(0, (int32_t) round(mapping_quality)));
    for (int i = 1; i < alignments.size(); ++i) {
        alignments[0].add_secondary_score(alignments[i].score());
    }
}

int32_t GSSWAligner::compute_mapping_quality(const vector<double>& scores, bool fast_approximation,
                                             const vector<double>* multiplicities) const {
    
    vector<double> scaled_scores(scores.size());
    for (size_t i = 0; i < scores.size(); i++) {
        scaled_scores[i] = log_base * scores[i];
    }
    size_t idx;
    return (int32_t) (fast_approximation ? maximum_mapping_quality_approx(scaled_scores, &idx, multiplicities)
                                         : maximum_mapping_quality_exact(scaled_scores, &idx, multiplicities));
}

int32_t GSSWAligner::compute_group_mapping_quality(const vector<double>& scores, const vector<size_t>& group,
                                                   const vector<double>* multiplicities) const {
    
    // make a non-const local version in case we need to sort it
    vector<size_t> non_const_group;
    const vector<size_t>* grp_ptr = &group;
    
    // ensure that group is in sorted order as following function expects
    if (!is_sorted(group.begin(), group.end())) {
        non_const_group = group;
        sort(non_const_group.begin(), non_const_group.end());
        grp_ptr = &non_const_group;
    }
    
    vector<double> scaled_scores(scores.size(), 0.0);
    for (size_t i = 0; i < scores.size(); i++) {
        scaled_scores[i] = log_base * scores[i];
    }
    return group_mapping_quality_exact(scaled_scores, *grp_ptr, multiplicities);
}

void GSSWAligner::compute_paired_mapping_quality(pair<vector<Alignment>, vector<Alignment>>& alignment_pairs,
                                                 const vector<double>& frag_weights,
                                                 int max_mapping_quality1,
                                                 int max_mapping_quality2,
                                                 bool fast_approximation,
                                                 double cluster_mq,
                                                 bool use_cluster_mq,
                                                 int overlap_count1,
                                                 int overlap_count2,
                                                 double mq_estimate1,
                                                 double mq_estimate2,
                                                 double maybe_mq_threshold,
                                                 double identity_weight) const {
    
    if (log_base <= 0.0) {
        cerr << "error:[Aligner] must call init_mapping_quality before computing mapping qualities" << endl;
        exit(EXIT_FAILURE);
    }
    
    size_t size = min(alignment_pairs.first.size(),
                      alignment_pairs.second.size());
    
    if (size == 0) {
        return;
    }
    
    vector<double> scaled_scores(size);
    
    for (size_t i = 0; i < size; i++) {
        auto& aln1 = alignment_pairs.first[i];
        auto& aln2 = alignment_pairs.second[i];
        scaled_scores[i] = log_base * (aln1.score() + aln2.score());
        // + frag_weights[i]);
        // ^^^ we could also incorporate the fragment weights, but this does not seem to help performance in the current form
    }

    size_t max_idx;
    double mapping_quality;
    if (!fast_approximation) {
        mapping_quality = maximum_mapping_quality_exact(scaled_scores, &max_idx);
    }
    else {
        mapping_quality = maximum_mapping_quality_approx(scaled_scores, &max_idx);
    }
    
    if (use_cluster_mq) {
        mapping_quality = prob_to_phred(sqrt(phred_to_prob(cluster_mq + mapping_quality)));
    }

    double mapping_quality1 = mapping_quality;
    double mapping_quality2 = mapping_quality;

    if (overlap_count1) {
        mapping_quality1 -= quality_scale_factor * log(overlap_count1);
    }
    if (overlap_count2) {
        mapping_quality2 -= quality_scale_factor * log(overlap_count2);
    }

    auto& max_aln1 = alignment_pairs.first.at(max_idx);
    int len1 = max(alignment_to_length(max_aln1), alignment_from_length(max_aln1));
    double identity1 = 1. - (double)(len1 * match - max_aln1.score()) / (match + mismatch) / len1;
    auto& max_aln2 = alignment_pairs.second.at(max_idx);
    int len2 = max(alignment_to_length(max_aln2), alignment_from_length(max_aln2));
    double identity2 = 1. - (double)(len2 * match - max_aln2.score()) / (match + mismatch) / len2;

    mapping_quality1 /= 2;
    mapping_quality2 /= 2;

    mapping_quality1 *= pow(identity1, identity_weight);
    mapping_quality2 *= pow(identity2, identity_weight);

    double mq_estimate = min(mq_estimate1, mq_estimate2);
    if (mq_estimate < maybe_mq_threshold && mq_estimate < mapping_quality1) {
        mapping_quality1 = prob_to_phred(sqrt(phred_to_prob(mq_estimate + mapping_quality1)));
    }
    if (mq_estimate < maybe_mq_threshold && mq_estimate < mapping_quality2) {
        mapping_quality2 = prob_to_phred(sqrt(phred_to_prob(mq_estimate + mapping_quality2)));
    }

    if (mapping_quality1 > max_mapping_quality1) {
        mapping_quality1 = max_mapping_quality1;
    }
    if (mapping_quality2 > max_mapping_quality2) {
        mapping_quality2 = max_mapping_quality2;
    }

    if (alignment_pairs.first[max_idx].score() == 0) {
        mapping_quality1 = 0;
    }
    if (alignment_pairs.second[max_idx].score() == 0) {
        mapping_quality2 = 0;
    }

    mapping_quality = max(0, (int32_t)round(min(mapping_quality1, mapping_quality2)));

    alignment_pairs.first[max_idx].set_mapping_quality(mapping_quality);
    alignment_pairs.second[max_idx].set_mapping_quality(mapping_quality);

    for (int i = 1; i < alignment_pairs.first.size(); ++i) {
        alignment_pairs.first[0].add_secondary_score(alignment_pairs.first[i].score());
    }
    for (int i = 1; i < alignment_pairs.second.size(); ++i) {
        alignment_pairs.second[0].add_secondary_score(alignment_pairs.second[i].score());
    }

}

double GSSWAligner::mapping_quality_score_diff(double mapping_quality) const {
    return mapping_quality / (quality_scale_factor * log_base);
}

double GSSWAligner::estimate_next_best_score(int length, double min_diffs) const {
    return ((length - min_diffs) * match - min_diffs * mismatch);
}

double GSSWAligner::max_possible_mapping_quality(int length) const {
    double max_score = log_base * length * match;
    vector<double> v = { max_score };
    size_t max_idx;
    return maximum_mapping_quality_approx(v, &max_idx);
}

double GSSWAligner::estimate_max_possible_mapping_quality(int length, double min_diffs, double next_min_diffs) const {
    double max_score = log_base * ((length - min_diffs) * match - min_diffs * mismatch);
    double next_max_score = log_base * ((length - next_min_diffs) * match - next_min_diffs * mismatch);
    vector<double> v = { max_score, next_max_score };
    size_t max_idx;
    return maximum_mapping_quality_approx(v, &max_idx);
}

double GSSWAligner::score_to_unnormalized_likelihood_ln(double score) const {
    // Log base needs to be set, or this can't work.
    assert(log_base != 0);
    // Likelihood is proportional to e^(lambda * score), so ln is just the exponent.
    return log_base * score;
}

size_t GSSWAligner::longest_detectable_gap(const Alignment& alignment, const string::const_iterator& read_pos) const {
    return longest_detectable_gap(alignment.sequence().size(), read_pos - alignment.sequence().begin());
}

size_t GSSWAligner::longest_detectable_gap(size_t read_length, size_t read_pos) const {
    // algebraic solution for when score is > 0 assuming perfect match other than gap
    int64_t overhang_length = min(read_pos, read_length - read_pos);
    int64_t numer = match * overhang_length + full_length_bonus;
    int64_t gap_length = (numer - gap_open) / gap_extension + 1;
    return gap_length >= 0 && overhang_length > 0 ? gap_length : 0;
}

size_t GSSWAligner::longest_detectable_gap(const Alignment& alignment) const {
    // longest detectable gap across entire read is in the middle
    return longest_detectable_gap(alignment.sequence().size(), alignment.sequence().size() / 2);
}

int32_t GSSWAligner::score_gappy_alignment(const Alignment& aln, const function<size_t(pos_t, pos_t, size_t)>& estimate_distance,
    bool strip_bonuses) const {
    
    int score = 0;
    int read_offset = 0;
    auto& path = aln.path();

    // We keep track of whether the last edit was a deletion for coalescing
    // adjacent deletions across node boundaries
    bool last_was_deletion = false;

    for (int i = 0; i < path.mapping_size(); ++i) {
        // For each mapping
        auto& mapping = path.mapping(i);
        for (int j = 0; j < mapping.edit_size(); ++j) {
            // For each edit in the mapping
            auto& edit = mapping.edit(j);
            
            // Score the edit according to its type
            if (edit_is_match(edit)) {
                score += score_exact_match(aln, read_offset, edit.to_length());
                last_was_deletion = false;
            } else if (edit_is_sub(edit)) {
                score -= mismatch * edit.sequence().size();
                last_was_deletion = false;
            } else if (edit_is_deletion(edit)) {
                if (last_was_deletion) {
                    // No need to charge a gap open
                    score -= edit.from_length() * gap_extension;
                } else {
                    // We need a gap open
                    score -= edit.from_length() ? gap_open + (edit.from_length() - 1) * gap_extension : 0;
                }

                if (edit.from_length()) {
                    // We already charged a gap open
                    last_was_deletion = true;
                }
                // If there's a 0-length deletion, leave the last_was_deletion flag unchanged.
            } else if (edit_is_insertion(edit) && !((i == 0 && j == 0) ||
                                                    (i == path.mapping_size()-1 && j == mapping.edit_size()-1))) {
                // todo how do we score this qual adjusted?
                score -= edit.to_length() ? gap_open + (edit.to_length() - 1) * gap_extension : 0;
                last_was_deletion = false;
                // No need to track if the last edit was an insertion because
                // insertions will be all together in a single edit at a point.
            } else {
                // Edit has no score effect. Probably a softclip.
                last_was_deletion = false;
            }
            read_offset += edit.to_length();
        }
        // score any intervening gaps in mappings using approximate distances
        if (i+1 < path.mapping_size()) {
            // what is the distance between the last position of this mapping
            // and the first of the next
            Position last_pos = mapping.position();
            last_pos.set_offset(last_pos.offset() + mapping_from_length(mapping));
            Position next_pos = path.mapping(i+1).position();
            // Estimate the distance
            int dist = estimate_distance(make_pos_t(last_pos), make_pos_t(next_pos), aln.sequence().size());
            if (dist > 0) {
                // If it's nonzero, score it as a deletion gap
                score -= gap_open + (dist - 1) * gap_extension;
            }
        }
    }

    if (!strip_bonuses) {
        // We should report any bonuses used in the DP in the final score
        if (!softclip_start(aln)) {
            score += full_length_bonus;
        }
        if (!softclip_end(aln)) {
            score += full_length_bonus;
        }
    }
    
    return score;
}

int32_t GSSWAligner::score_ungapped_alignment(const Alignment& aln, bool strip_bonuses) const {
    return score_gappy_alignment(aln, [](pos_t, pos_t, size_t){return (size_t) 0;}, strip_bonuses);
}

int32_t GSSWAligner::remove_bonuses(const Alignment& aln, bool pinned, bool pin_left) const {
    int32_t score = aln.score();
    if (softclip_start(aln) == 0 && !(pinned && pin_left)) {
        // No softclip at the start, and a left end bonus was applied.
        score -= full_length_bonus;
    }
    if (softclip_end(aln) == 0 && !(pinned && !pin_left)) {
        // No softclip at the end, and a right end bonus was applied.
        score -= full_length_bonus;
    }
    return score;
}

Aligner::Aligner(const int8_t* _score_matrix,
                 int8_t _gap_open,
                 int8_t _gap_extension,
                 int8_t _full_length_bonus,
                 double _gc_content)
{
    // TODO: now that everything is in terms of score matrices, having match/mismatch is a bit
    // misleading, but a fair amount of code depends on them
    match = _score_matrix[0];
    mismatch = -_score_matrix[1];
    gap_open = _gap_open;
    gap_extension = _gap_extension;
    full_length_bonus = _full_length_bonus;

    // table to translate chars to their integer value
    nt_table = gssw_create_nt_table();
    
    // calculate the scale of the scores so we can use it in mapping quality calculations
    init_mapping_quality(_score_matrix, _gc_content);
    
    // add in the 5th row and column of 0s for N matches like GSSW wants
    score_matrix = (int8_t*) malloc(sizeof(int8_t) * 25);
    for (size_t i = 0, j = 0; i < 25; ++i) {
        if (i % 5 == 4 || i / 5 == 4) {
            score_matrix[i] = 0;
        }
        else {
            score_matrix[i] = _score_matrix[j];
            ++j;
        }
    }
    
    // make an XdropAligner for each thread
    int num_threads = get_thread_count();
    xdrops.reserve(num_threads);
    for (size_t i = 0; i < num_threads; ++i) {
        xdrops.emplace_back(_score_matrix, _gap_open, _gap_extension, _full_length_bonus);
    }
}

void Aligner::align_internal(Alignment& alignment, vector<Alignment>* multi_alignments, const HandleGraph& g,
                             bool pinned, bool pin_left,int32_t max_alt_alns, bool traceback_aln) const {
    // bench_start(bench);
    // check input integrity
    if (pin_left && !pinned) {
        cerr << "error:[Aligner] cannot choose pinned end in non-pinned alignment" << endl;
        exit(EXIT_FAILURE);
    }
    if (multi_alignments && !pinned) {
        cerr << "error:[Aligner] multiple traceback is not implemented in local alignment, only pinned and global" << endl;
        exit(EXIT_FAILURE);
    }
    if (!multi_alignments && max_alt_alns != 1) {
        cerr << "error:[Aligner] cannot specify maximum number of alignments in single alignment" << endl;
        exit(EXIT_FAILURE);
    }
    if (max_alt_alns <= 0) {
        cerr << "error:[Aligner] cannot do less than 1 alignment" << endl;
        exit(EXIT_FAILURE);
    }

    // alignment pinning algorithm is based on pinning in bottom right corner, if pinning in top
    // left we need to reverse all the sequences first and translate the alignment back later
    
    // make a place to reverse the graph and sequence if necessary
    ReverseGraph reversed_graph(&g, false);
    string reversed_sequence;

    // choose forward or reversed objects
    const HandleGraph* oriented_graph = &g;
    const string* align_sequence = &alignment.sequence();
    if (pin_left) {
        // choose the reversed graph
        oriented_graph = &reversed_graph;
        
        // make and assign the reversed sequence
        reversed_sequence.resize(align_sequence->size());
        reverse_copy(align_sequence->begin(), align_sequence->end(), reversed_sequence.begin());
        align_sequence = &reversed_sequence;
    }
    
    // to save compute, we won't make these unless we're doing pinning
    unordered_set<vg::id_t> pinning_ids;
    NullMaskingGraph* null_masked_graph = nullptr;
    const HandleGraph* align_graph = oriented_graph;
    if (pinned) {
        pinning_ids = identify_pinning_points(*oriented_graph);
        null_masked_graph = new NullMaskingGraph(oriented_graph);
        align_graph = null_masked_graph;
    }
    
    // convert into gssw graph
    gssw_graph* graph = create_gssw_graph(*align_graph);
    
    // perform dynamic programming
    gssw_graph_fill_pinned(graph, align_sequence->c_str(),
                           nt_table, score_matrix,
                           gap_open, gap_extension, full_length_bonus,
                           pinned ? 0 : full_length_bonus, 15, 2, traceback_aln);

    // traceback either from pinned position or optimal local alignment
    if (traceback_aln) {
        if (pinned) {
            // we can only run gssw's DP on non-empty graphs, but we may have masked the entire graph
            // if it consists of only empty nodes, so don't both with the DP in that case
            gssw_graph_mapping** gms = nullptr;
            if (align_graph->get_node_count() > 0) {
                gssw_node** pinning_nodes = (gssw_node**) malloc(pinning_ids.size() * sizeof(gssw_node*));
                size_t j = 0;
                for (size_t i = 0; i < graph->size; i++) {
                    gssw_node* node = graph->nodes[i];
                    if (pinning_ids.count(node->id)) {
                        pinning_nodes[j] = node;
                        j++;
                    }
                }
                
                // trace back pinned alignment
                gms = gssw_graph_trace_back_pinned_multi (graph,
                                                          max_alt_alns,
                                                          true,
                                                          align_sequence->c_str(),
                                                          align_sequence->size(),
                                                          pinning_nodes,
                                                          pinning_ids.size(),
                                                          nt_table,
                                                          score_matrix,
                                                          gap_open,
                                                          gap_extension,
                                                          full_length_bonus,
                                                          0);
                
                free(pinning_nodes);
            }
            
            // did we both 1) do DP (i.e. the graph is non-empty), and 2) find a traceback with positive score?
            if (gms ? gms[0]->score > 0 : false) {
                
                if (pin_left) {
                    // translate nodes and mappings into original sequence so that the cigars come out right
                    unreverse_graph(graph);
                    for (int32_t i = 0; i < max_alt_alns; i++) {
                        unreverse_graph_mapping(gms[i]);
                    }
                }
                
                // have a mapping, can just convert normally
                gssw_mapping_to_alignment(graph, gms[0], alignment, pinned, pin_left);
                
                if (multi_alignments) {
                    // determine how many non-null alignments were returned
                    int32_t num_non_null = max_alt_alns;
                    for (int32_t i = 1; i < max_alt_alns; i++) {
                        if (gms[i]->score <= 0) {
                            num_non_null = i;
                            break;
                        }
                    }
                    
                    // reserve to avoid illegal access errors that occur when the vector reallocates
                    multi_alignments->reserve(num_non_null);
                    
                    // copy the primary alignment
                    multi_alignments->emplace_back(alignment);
                    
                    // convert the alternate alignments and store them at the back of the vector (this will not
                    // execute if we are doing single alignment)
                    for (int32_t i = 1; i < num_non_null; i++) {
                        // make new alignment object
                        multi_alignments->emplace_back();
                        Alignment& next_alignment = multi_alignments->back();
                        
                        // copy over sequence information from the primary alignment
                        next_alignment.set_sequence(alignment.sequence());
                        next_alignment.set_quality(alignment.quality());
                        
                        // get path of the alternate alignment
                        gssw_mapping_to_alignment(graph, gms[i], next_alignment, pinned, pin_left);
                    }
                }
            }
            else if (g.get_node_count() > 0) {
                // we didn't get any alignments either because the graph was empty and we couldn't run
                // gssw DP or because they had score 0 and gssw didn't want to do traceback. however,
                // we can infer the location of softclips based on the pinning nodes, so we'll just make
                // those manually
                
                // find the sink nodes of the oriented graph, which may be empty
                auto pinning_points = algorithms::tail_nodes(oriented_graph);
                // impose a consistent ordering for machine independent behavior
                sort(pinning_points.begin(), pinning_points.end(), [&](const handle_t& h1, const handle_t& h2) {
                    return oriented_graph->get_id(h1) < oriented_graph->get_id(h2);
                });
                
                for (size_t i = 0; i < max_alt_alns && i < pinning_points.size(); i++) {
                    // make a record in the multi alignments if we're using them
                    if (multi_alignments) {
                        multi_alignments->emplace_back();
                    }
                    // choose an alignment object to construct the path in
                    Alignment& softclip_alignment = i == 0 ? alignment : multi_alignments->back();
                    
                    handle_t& pinning_point = pinning_points[i];
                    
                    Mapping* mapping = alignment.mutable_path()->add_mapping();
                    mapping->set_rank(1);
                    
                    // locate at the beginning or end of the node
                    Position* position = mapping->mutable_position();
                    position->set_node_id(oriented_graph->get_id(pinning_point));
                    position->set_offset(pin_left ? 0 : oriented_graph->get_length(pinning_point));
                    
                    // soft clip
                    Edit* edit = mapping->add_edit();
                    edit->set_to_length(alignment.sequence().length());
                    edit->set_sequence(alignment.sequence());
                    
                    // we want to also have the first alignment in the multi-alignment vector
                    if (i == 0 && multi_alignments) {
                        multi_alignments->back() = alignment;
                    }
                }
            }
            if (gms) {
                for (int32_t i = 0; i < max_alt_alns; i++) {
                    gssw_graph_mapping_destroy(gms[i]);
                }
                free(gms);
            }
        }
        else {
            // trace back local alignment
            gssw_graph_mapping* gm = gssw_graph_trace_back (graph,
                                                            align_sequence->c_str(),
                                                            align_sequence->size(),
                                                            nt_table,
                                                            score_matrix,
                                                            gap_open,
                                                            gap_extension,
                                                            full_length_bonus,
                                                            full_length_bonus);
        
            gssw_mapping_to_alignment(graph, gm, alignment, pinned, pin_left);
            gssw_graph_mapping_destroy(gm);
        }
    } else {
        // get the alignment position and score
        alignment.set_score(graph->max_node->alignment->score1);
        Mapping* m = alignment.mutable_path()->add_mapping();
        Position* p = m->mutable_position();
        p->set_node_id(graph->max_node->id);
        p->set_offset(graph->max_node->alignment->ref_end1); // mark end position; for de-duplication
    }
        
    // this might be null if we're not doing pinned alignment, but delete doesn't care
    delete null_masked_graph;
    
    gssw_graph_destroy(graph);
    // bench_end(bench);
}

void Aligner::align(Alignment& alignment, const HandleGraph& g, bool traceback_aln) const {
    
    align_internal(alignment, nullptr, g, false, false, 1, traceback_aln);
}

void Aligner::align(Alignment& alignment, const HandleGraph& g,
                    const std::vector<handle_t>& topological_order) const {

    // Create a gssw_graph and a mapping from handles to nodes.
    gssw_graph* graph = gssw_graph_create(topological_order.size());
    hash_map<handle_t, gssw_node*> nodes;
    nodes.reserve(topological_order.size());

    // Create the nodes. Use offsets in the topological order as node ids.
    for (size_t i = 0; i < topological_order.size(); i++) {
        handle_t handle = topological_order[i];
        auto cleaned_seq = nonATGCNtoN(g.get_sequence(handle));
        gssw_node* node = gssw_node_create(nullptr,
                                           i,
                                           cleaned_seq.c_str(),
                                           nt_table,
                                           score_matrix);
        nodes[handle] = node;
        gssw_graph_add_node(graph, node);
    }

    // Create the edges.
    for (const handle_t& from : topological_order) {
        gssw_node* from_node = nodes[from];
        g.follow_edges(from, false, [&](const handle_t& to) {
            auto iter = nodes.find(to);
            if (iter != nodes.end()) {
                gssw_nodes_add_edge(from_node, iter->second);
            }
        });
    }

    // Align the read to the subgraph.
    gssw_graph_fill_pinned(graph, alignment.sequence().c_str(),
                           nt_table, score_matrix,
                           gap_open, gap_extension, full_length_bonus, full_length_bonus,
                           15, 2, true);
    gssw_graph_mapping* gm = gssw_graph_trace_back(graph,
                                                   alignment.sequence().c_str(), alignment.sequence().length(),
                                                   nt_table, score_matrix,
                                                   gap_open, gap_extension, full_length_bonus, full_length_bonus);

    // Convert the mapping to Alignment.
    this->gssw_mapping_to_alignment(graph, gm, alignment, false, false);
    Path& path = *(alignment.mutable_path());
    for (size_t i = 0; i < path.mapping_size(); i++) {
        Position& pos = *(path.mutable_mapping(i)->mutable_position());
        handle_t handle = topological_order[pos.node_id()];
        pos.set_node_id(g.get_id(handle));
        pos.set_is_reverse(g.get_is_reverse(handle));
    }

    // Destroy the temporary objects.
    gssw_graph_mapping_destroy(gm);
    gssw_graph_destroy(graph);
}

void Aligner::align_pinned(Alignment& alignment, const HandleGraph& g, bool pin_left, bool xdrop,
                           uint16_t xdrop_max_gap_length) const {
    
    if (xdrop) {
        // XdropAligner manages its own stack, so it can never be threadsafe without be recreated
        // for every alignment, which meshes poorly with its stack implementation. We achieve
        // thread-safety by having one per thread, which makes this method const-ish.
        XdropAligner& xdrop = const_cast<XdropAligner&>(xdrops[omp_get_thread_num()]);
        
        // wrap the graph so that empty pinning points are handled correctly
        DozeuPinningOverlay overlay(&g, !pin_left);
        
        if (overlay.get_node_count() == 0 && g.get_node_count() > 0) {
            // the only nodes in the graph are empty nodes for pinning, which got masked.
            // we can still infer a pinned alignment based purely on the pinning point but
            // dozeu won't handle this correctly
            g.for_each_handle([&](const handle_t& handle) {
                bool can_pin = g.follow_edges(handle, pin_left, [&](const handle_t& next) {return false;});
                if (can_pin) {
                    // manually make the softclip
                    Mapping* mapping = alignment.mutable_path()->add_mapping();
                    Position* pos = mapping->mutable_position();
                    pos->set_node_id(g.get_id(handle));
                    pos->set_is_reverse(false);
                    pos->set_offset(pin_left ? 0 : g.get_length(handle));
                    
                    mapping->set_rank(1);
                    
                    Edit* edit = mapping->add_edit();
                    edit->set_from_length(0);
                    edit->set_to_length(alignment.sequence().size());
                    edit->set_sequence(alignment.sequence());
                    alignment.set_score(0);
                    return false;
                }
                return true;
            });
        }
        else {
            // do the alignment
            xdrop.align_pinned(alignment, overlay, pin_left, xdrop_max_gap_length);
            
            if (overlay.performed_duplications()) {
                // the overlay is not a strict subset of the underlying graph, so we may
                // need to translate some node IDs
                translate_oriented_node_ids(*alignment.mutable_path(), [&](id_t node_id) {
                    handle_t under = overlay.get_underlying_handle(overlay.get_handle(node_id));
                    return make_pair(g.get_id(under), g.get_is_reverse(under));
                });
            }
        }
    }
    else {
        align_internal(alignment, nullptr, g, true, pin_left, 1, true);
    }
}

void Aligner::align_pinned_multi(Alignment& alignment, vector<Alignment>& alt_alignments, const HandleGraph& g,
                                 bool pin_left, int32_t max_alt_alns) const {
    
    if (alt_alignments.size() != 0) {
        cerr << "error:[Aligner::align_pinned_multi] output vector must be empty for pinned multi-aligning" << endl;
        exit(EXIT_FAILURE);
    }
    
    align_internal(alignment, &alt_alignments, g, true, pin_left, max_alt_alns, true);
}

void Aligner::align_global_banded(Alignment& alignment, const HandleGraph& g,
                                  int32_t band_padding, bool permissive_banding) const {
    
    // We need to figure out what size ints we need to use.
    // Get upper and lower bounds on the scores. TODO: if these overflow int64 we're out of luck
    int64_t best_score = alignment.sequence().size() * match;
    size_t total_bases = 0;
    g.for_each_handle([&](const handle_t& handle) {
        total_bases += g.get_length(handle);
    });
    int64_t worst_score = (alignment.sequence().size() + total_bases) * -max(max(mismatch, gap_open), gap_extension);
    
    // TODO: put this all into another template somehow?
    
    if (best_score <= numeric_limits<int8_t>::max() && worst_score >= numeric_limits<int8_t>::min()) {
        // We'll fit in int8
        BandedGlobalAligner<int8_t> band_graph(alignment,
                                               g,
                                               band_padding,
                                               permissive_banding,
                                               false);
        
        band_graph.align(score_matrix, nt_table, gap_open, gap_extension);
    } else if (best_score <= numeric_limits<int16_t>::max() && worst_score >= numeric_limits<int16_t>::min()) {
        // We'll fit in int16
        BandedGlobalAligner<int16_t> band_graph(alignment,
                                                g,
                                                band_padding,
                                                permissive_banding,
                                                false);
        
        band_graph.align(score_matrix, nt_table, gap_open, gap_extension);
    } else if (best_score <= numeric_limits<int32_t>::max() && worst_score >= numeric_limits<int32_t>::min()) {
        // We'll fit in int32
        BandedGlobalAligner<int32_t> band_graph(alignment,
                                                g,
                                                band_padding,
                                                permissive_banding,
                                                false);
        
        band_graph.align(score_matrix, nt_table, gap_open, gap_extension);
    } else {
        // Fall back to int64
        BandedGlobalAligner<int64_t> band_graph(alignment,
                                                g,
                                                band_padding,
                                                permissive_banding,
                                                false);
        
        band_graph.align(score_matrix, nt_table, gap_open, gap_extension);
    }
}

void Aligner::align_global_banded_multi(Alignment& alignment, vector<Alignment>& alt_alignments, const HandleGraph& g,
                                        int32_t max_alt_alns, int32_t band_padding, bool permissive_banding) const {
                                        
    // We need to figure out what size ints we need to use.
    // Get upper and lower bounds on the scores. TODO: if these overflow int64 we're out of luck
    int64_t best_score = alignment.sequence().size() * match;
    size_t total_bases = 0;
    g.for_each_handle([&](const handle_t& handle) {
        total_bases += g.get_length(handle);
    });
    int64_t worst_score = (alignment.sequence().size() + total_bases) * -max(max(mismatch, gap_open), gap_extension);
    
    if (best_score <= numeric_limits<int8_t>::max() && worst_score >= numeric_limits<int8_t>::min()) {
        // We'll fit in int8
        BandedGlobalAligner<int8_t> band_graph(alignment,
                                               g,
                                               alt_alignments,
                                               max_alt_alns,
                                               band_padding,
                                               permissive_banding,
                                               false);
        
        band_graph.align(score_matrix, nt_table, gap_open, gap_extension);
    } else if (best_score <= numeric_limits<int16_t>::max() && worst_score >= numeric_limits<int16_t>::min()) {
        // We'll fit in int16
        BandedGlobalAligner<int16_t> band_graph(alignment,
                                                g,
                                                alt_alignments,
                                                max_alt_alns,
                                                band_padding,
                                                permissive_banding,
                                                false);
        
        band_graph.align(score_matrix, nt_table, gap_open, gap_extension);
    } else if (best_score <= numeric_limits<int32_t>::max() && worst_score >= numeric_limits<int32_t>::min()) {
        // We'll fit in int32
        BandedGlobalAligner<int32_t> band_graph(alignment,
                                                g,
                                                alt_alignments,
                                                max_alt_alns,
                                                band_padding,
                                                permissive_banding,
                                                false);
        
        band_graph.align(score_matrix, nt_table, gap_open, gap_extension);
    } else {
        // Fall back to int64
        BandedGlobalAligner<int64_t> band_graph(alignment,
                                                g,
                                                alt_alignments,
                                                max_alt_alns,
                                                band_padding,
                                                permissive_banding,
                                                false);
        
        band_graph.align(score_matrix, nt_table, gap_open, gap_extension);
    }
}

void Aligner::align_xdrop(Alignment& alignment, const HandleGraph& g, const vector<MaximalExactMatch>& mems,
                          bool reverse_complemented, uint16_t max_gap_length) const
{
    align_xdrop(alignment, g, algorithms::lazier_topological_order(&g), mems, reverse_complemented, max_gap_length);
}

void Aligner::align_xdrop(Alignment& alignment, const HandleGraph& g, const vector<handle_t>& order,
                          const vector<MaximalExactMatch>& mems, bool reverse_complemented, uint16_t max_gap_length) const
{
    // XdropAligner manages its own stack, so it can never be threadsafe without be recreated
    // for every alignment, which meshes poorly with its stack implementation. We achieve
    // thread-safety by having one per thread, which makes this method const-ish.
    XdropAligner& xdrop = const_cast<XdropAligner&>(xdrops[omp_get_thread_num()]);
    xdrop.align(alignment, g, order, mems, reverse_complemented, max_gap_length);
    if (!alignment.has_path() && mems.empty()) {
        // dozeu couldn't find an alignment, probably because it's seeding heuristic failed
        // we'll just fall back on GSSW
        // TODO: This is a bit inconsistent. GSSW gives a full-length bonus at both ends, while
        // dozeu only gives it once.
        align(alignment, g, order);
    }
}


// Scoring an exact match is very simple in an ordinary Aligner

int32_t Aligner::score_exact_match(const Alignment& aln, size_t read_offset, size_t length) const {
    return match * length;
}

int32_t Aligner::score_exact_match(const string& sequence) const {
    return match * sequence.length();
}

int32_t Aligner::score_exact_match(string::const_iterator seq_begin, string::const_iterator seq_end) const {
    return match * (seq_end - seq_begin);
}

int32_t Aligner::score_exact_match(const string& sequence, const string& base_quality) const {
    return score_exact_match(sequence);
}

int32_t Aligner::score_exact_match(string::const_iterator seq_begin, string::const_iterator seq_end,
                                   string::const_iterator base_qual_begin) const {
    return score_exact_match(seq_begin, seq_end);
}

int32_t Aligner::score_mismatch(size_t length) const {
    return -match * length;
}

int32_t Aligner::score_partial_alignment(const Alignment& alignment, const HandleGraph& graph, const Path& path,
                                         string::const_iterator seq_begin) const {
    
    int32_t score = 0;
    string::const_iterator read_pos = seq_begin;
    for (size_t i = 0; i < path.mapping_size(); i++) {
        const Mapping& mapping = path.mapping(i);
        
        for (size_t j = 0; j < mapping.edit_size(); j++) {
            const Edit& edit = mapping.edit(j);
            
            if (edit.from_length() > 0) {
                if (edit.to_length() > 0) {
                    if (edit.sequence().empty()) {
                        // match
                        score += match * edit.from_length();
                    }
                    else {
                        // mismatch
                        score -= mismatch * edit.from_length();
                    }
                    
                    // apply full length bonus
                    if (read_pos == alignment.sequence().begin()) {
                        score += full_length_bonus;
                    }
                    if (read_pos + edit.from_length() == alignment.sequence().end()) {
                        score += full_length_bonus;
                    }
                }
                else {
                    // deletion
                    score -= gap_open + (edit.from_length() - 1) * gap_extension;
                }
            }
            else if (edit.to_length() > 0) {
                // don't score soft clips
                if (read_pos != alignment.sequence().begin() &&
                    read_pos + edit.to_length() != alignment.sequence().end()) {
                    // insert
                    score -= gap_open + (edit.to_length() - 1) * gap_extension;
                }
            }
            
            read_pos += edit.to_length();
        }
    }
    return score;
}

QualAdjAligner::QualAdjAligner(const int8_t* _score_matrix,
                               int8_t _gap_open,
                               int8_t _gap_extension,
                               int8_t _full_length_bonus,
                               double _gc_content)
{
    
    // TODO: now that everything is in terms of score matrices, having match/mismatch is a bit
    // misleading, but a fair amount of code depends on them
    match = _score_matrix[0];
    mismatch = -_score_matrix[1];
    gap_open = _gap_open;
    gap_extension = _gap_extension;
    full_length_bonus = _full_length_bonus;
    
    // table to translate chars to their integer value
    nt_table = gssw_create_nt_table();
    
    // calculate the scale of the scores so we can use it in mapping quality calculations
    init_mapping_quality(_score_matrix, _gc_content);
    
    // TODO: this interface could really be improved in GSSW, oh well though
    
    // convert gc content into base-wise frequencies
    double* nt_freqs = (double*) malloc(sizeof(double) * 4);
    nt_freqs[0] = 0.5 * (1 - _gc_content);
    nt_freqs[1] = 0.5 * _gc_content;
    nt_freqs[2] = 0.5 * _gc_content;
    nt_freqs[3] = 0.5 * (1 - _gc_content);
    
    // find max score and set it to be the max so that scaling doesn't change score
    // TODO: hacky
    int8_t max_score = max(gap_open, gap_extension);
    for (size_t i = 0; i < 16; ++i) {
        max_score = max<int8_t>(max_score, abs(_score_matrix[i]));
    }
    
    uint8_t max_base_qual = 255;
    size_t num_nts = 4;
    
    // find the quality-adjusted, scaled scores
    int8_t* scaled_matrix = gssw_scaled_adjusted_qual_matrix(max_score, max_base_qual, &gap_open, &gap_extension,
                                                             _score_matrix, nt_freqs, num_nts, 1e-10);
    
        
    // finally, add in the 0s to the 5-th row and column for Ns
    score_matrix = gssw_add_ambiguous_char_to_adjusted_matrix(scaled_matrix, max_base_qual, num_nts);
    
    
    // make a QualAdjXdropAligner for each thread
    int num_threads = get_thread_count();
    xdrops.reserve(num_threads);
    for (size_t i = 0; i < num_threads; ++i) {
        xdrops.emplace_back(_score_matrix, scaled_matrix, _gap_open, _gap_extension, _full_length_bonus);
    }
    
    // free the temporary arrays we allocated
    free(scaled_matrix);
    free(nt_freqs);
}

void QualAdjAligner::align_internal(Alignment& alignment, vector<Alignment>* multi_alignments, const HandleGraph& g,
                                    bool pinned, bool pin_left, int32_t max_alt_alns, bool traceback_aln) const {
    
    // check input integrity
    if (pin_left && !pinned) {
        cerr << "error:[Aligner] cannot choose pinned end in non-pinned alignment" << endl;
        exit(EXIT_FAILURE);
    }
    if (multi_alignments && !pinned) {
        cerr << "error:[Aligner] multiple traceback is not implemented in local alignment, only pinned and global" << endl;
        exit(EXIT_FAILURE);
    }
    if (!multi_alignments && max_alt_alns != 1) {
        cerr << "error:[Aligner] cannot specify maximum number of alignments in single alignment" << endl;
        exit(EXIT_FAILURE);
    }
    if (max_alt_alns <= 0) {
        cerr << "error:[Aligner] cannot do less than 1 alignment" << endl;
        exit(EXIT_FAILURE);
    }
    
    // alignment pinning algorithm is based on pinning in bottom right corner, if pinning in top
    // left we need to reverse all the sequences first and translate the alignment back later
    
    // make a place to reverse the graph and sequence if necessary
    ReverseGraph reversed_graph(&g, false);
    string reversed_sequence;
    string reversed_quality;
    
    // choose forward or reversed objects
    const HandleGraph* oriented_graph = &g;
    const string* align_sequence = &alignment.sequence();
    const string* align_quality = &alignment.quality();
    if (pin_left) {
        // choose the reversed graph
        oriented_graph = &reversed_graph;
        
        // make and assign the reversed sequence
        reversed_sequence.resize(align_sequence->size());
        reverse_copy(align_sequence->begin(), align_sequence->end(), reversed_sequence.begin());
        align_sequence = &reversed_sequence;
        
        // make and assign the reversed quality
        reversed_quality.resize(align_quality->size());
        reverse_copy(align_quality->begin(), align_quality->end(), reversed_quality.begin());
        align_quality = &reversed_quality;
    }
    
    if (align_quality->size() != align_sequence->size()) {
        cerr << "error:[QualAdjAligner] Read " << alignment.name() << " has sequence and quality strings with different lengths. Cannot perform base quality adjusted alignment. Consider toggling off base quality adjusted alignment at the command line." << endl;
        exit(EXIT_FAILURE);
    }
    
    // to save compute, we won't make these unless we're doing pinning
    unordered_set<vg::id_t> pinning_ids;
    NullMaskingGraph* null_masked_graph = nullptr;
    const HandleGraph* align_graph = oriented_graph;
    if (pinned) {
        pinning_ids = identify_pinning_points(*oriented_graph);
        null_masked_graph = new NullMaskingGraph(oriented_graph);
        align_graph = null_masked_graph;
    }
    
    // convert into gssw graph
    gssw_graph* graph = create_gssw_graph(*align_graph);
    
    // perform dynamic programming
    // offer a full length bonus on each end, or only on the left if the right end is pinned.
    gssw_graph_fill_pinned_qual_adj(graph, align_sequence->c_str(), align_quality->c_str(),
                                    nt_table, score_matrix,
                                    gap_open, gap_extension,
                                    full_length_bonus, pinned ? 0 : full_length_bonus, 15, 2, traceback_aln);
    
    // traceback either from pinned position or optimal local alignment
    if (traceback_aln) {
        if (pinned) {
            gssw_graph_mapping** gms = nullptr;
            if (align_graph->get_node_count() > 0) {
                
                gssw_node** pinning_nodes = (gssw_node**) malloc(pinning_ids.size() * sizeof(gssw_node*));
                size_t j = 0;
                for (size_t i = 0; i < graph->size; i++) {
                    gssw_node* node = graph->nodes[i];
                    if (pinning_ids.count(node->id)) {
                        pinning_nodes[j] = node;
                        j++;
                    }
                }
                
                // trace back pinned alignment
                gms = gssw_graph_trace_back_pinned_qual_adj_multi (graph,
                                                                   max_alt_alns,
                                                                   true,
                                                                   align_sequence->c_str(),
                                                                   align_quality->c_str(),
                                                                   align_sequence->size(),
                                                                   pinning_nodes,
                                                                   pinning_ids.size(),
                                                                   nt_table,
                                                                   score_matrix,
                                                                   gap_open,
                                                                   gap_extension,
                                                                   full_length_bonus,
                                                                   0);
                
                free(pinning_nodes);
            }
            
            // did we both 1) do DP (i.e. the graph is non-empty), and 2) find a traceback with positive score?
            if (gms && gms[0]->score > 0) {
                
                if (pin_left) {
                    // translate graph and mappings into original node space
                    unreverse_graph(graph);
                    for (int32_t i = 0; i < max_alt_alns; i++) {
                        unreverse_graph_mapping(gms[i]);
                    }
                }
                
                // have a mapping, can just convert normally
                gssw_mapping_to_alignment(graph, gms[0], alignment, pinned, pin_left);
                
                if (multi_alignments) {
                    // determine how many non-null alignments were returned
                    int32_t num_non_null = max_alt_alns;
                    for (int32_t i = 1; i < max_alt_alns; i++) {
                        if (gms[i]->score <= 0) {
                            num_non_null = i;
                            break;
                        }
                    }
                    
                    // reserve to avoid illegal access errors that occur when the vector reallocates
                    multi_alignments->reserve(num_non_null);
                    
                    // copy the primary alignment
                    multi_alignments->emplace_back(alignment);
                    
                    // convert the alternate alignments and store them at the back of the vector (this will not
                    // execute if we are doing single alignment)
                    for (int32_t i = 1; i < num_non_null; i++) {
                        // make new alignment object
                        multi_alignments->emplace_back();
                        Alignment& next_alignment = multi_alignments->back();
                        
                        // copy over sequence information from the primary alignment
                        next_alignment.set_sequence(alignment.sequence());
                        next_alignment.set_quality(alignment.quality());
                        
                        // get path of the alternate alignment
                        gssw_mapping_to_alignment(graph, gms[i], next_alignment, pinned, pin_left);
                    }
                }
            }
            else if (g.get_node_count() > 0) {
                /// we didn't get any alignments either because the graph was empty and we couldn't run
                // gssw DP or because they had score 0 and gssw didn't want to do traceback. however,
                // we can infer the location of softclips based on the pinning nodes, so we'll just make
                // those manually
                
                // find the sink nodes of the oriented graph, which may be empty
                auto pinning_points = algorithms::tail_nodes(oriented_graph);
                // impose a consistent ordering for machine independent behavior
                sort(pinning_points.begin(), pinning_points.end(), [&](const handle_t& h1, const handle_t& h2) {
                    return oriented_graph->get_id(h1) < oriented_graph->get_id(h2);
                });
                
                for (size_t i = 0; i < max_alt_alns && i < pinning_points.size(); i++) {
                    // make a record in the multi alignments if we're using them
                    if (multi_alignments) {
                        multi_alignments->emplace_back();
                    }
                    // choose an alignment object to construct the path in
                    Alignment& softclip_alignment = i == 0 ? alignment : multi_alignments->back();
                    
                    handle_t& pinning_point = pinning_points[i];
                    
                    Mapping* mapping = alignment.mutable_path()->add_mapping();
                    mapping->set_rank(1);
                    
                    // locate at the beginning or end of the node
                    Position* position = mapping->mutable_position();
                    position->set_node_id(oriented_graph->get_id(pinning_point));
                    position->set_offset(pin_left ? 0 : oriented_graph->get_length(pinning_point));
                    
                    // soft clip
                    Edit* edit = mapping->add_edit();
                    edit->set_to_length(alignment.sequence().length());
                    edit->set_sequence(alignment.sequence());
                    
                    // we want to also have the first alignment in the multi-alignment vector
                    if (i == 0 && multi_alignments) {
                        multi_alignments->back() = alignment;
                    }
                }
            }
            
            if (gms) {
                for (int32_t i = 0; i < max_alt_alns; i++) {
                    gssw_graph_mapping_destroy(gms[i]);
                }
                free(gms);
            }
        }
        else {
            // trace back local alignment
            gssw_graph_mapping* gm = gssw_graph_trace_back_qual_adj (graph,
                                                                     align_sequence->c_str(),
                                                                     align_quality->c_str(),
                                                                     align_sequence->size(),
                                                                     nt_table,
                                                                     score_matrix,
                                                                     gap_open,
                                                                     gap_extension,
                                                                     full_length_bonus,
                                                                     full_length_bonus);
        
            gssw_mapping_to_alignment(graph, gm, alignment, pinned, pin_left);
            gssw_graph_mapping_destroy(gm);
        }
    } else {
        // get the alignment position and score
        alignment.set_score(graph->max_node->alignment->score1);
        Mapping* m = alignment.mutable_path()->add_mapping();
        Position* p = m->mutable_position();
        p->set_node_id(graph->max_node->id);
        p->set_offset(graph->max_node->alignment->ref_end1); // mark end position; for de-duplication
    }
        
    // this might be null if we're not doing pinned alignment, but delete doesn't care
    delete null_masked_graph;
    
    gssw_graph_destroy(graph);
    
}

void QualAdjAligner::align(Alignment& alignment, const HandleGraph& g, bool traceback_aln) const {
    
    align_internal(alignment, nullptr, g, false, false, 1, traceback_aln);
}

void QualAdjAligner::align_pinned(Alignment& alignment, const HandleGraph& g, bool pin_left, bool xdrop,
                                  uint16_t xdrop_max_gap_length) const {
    if (xdrop) {
        // QualAdjXdropAligner manages its own stack, so it can never be threadsafe without be recreated
        // for every alignment, which meshes poorly with its stack implementation. We achieve
        // thread-safety by having one per thread, which makes this method const-ish.
        QualAdjXdropAligner& xdrop = const_cast<QualAdjXdropAligner&>(xdrops[omp_get_thread_num()]);
        
        // wrap the graph so that empty pinning points are handled correctly
        DozeuPinningOverlay overlay(&g, !pin_left);
        if (overlay.get_node_count() == 0 && g.get_node_count() > 0) {
            // the only nodes in the graph are empty nodes for pinning, which got masked.
            // we can still infer a pinned alignment based purely on the pinning point but
            // dozeu won't handle this correctly
            g.for_each_handle([&](const handle_t& handle) {
                bool can_pin = g.follow_edges(handle, pin_left, [&](const handle_t& next) {return false;});
                if (can_pin) {
                    // manually make the softclip
                    Mapping* mapping = alignment.mutable_path()->add_mapping();
                    Position* pos = mapping->mutable_position();
                    pos->set_node_id(g.get_id(handle));
                    pos->set_is_reverse(false);
                    pos->set_offset(pin_left ? 0 : g.get_length(handle));
                    
                    mapping->set_rank(1);
                    
                    Edit* edit = mapping->add_edit();
                    edit->set_from_length(0);
                    edit->set_to_length(alignment.sequence().size());
                    edit->set_sequence(alignment.sequence());
                    alignment.set_score(0);
                    return false;
                }
                return true;
            });
        }
        else {
            //
            
            xdrop.align_pinned(alignment, overlay, pin_left, xdrop_max_gap_length);
            
            if (overlay.performed_duplications()) {
                // the overlay is not a strict subset of the underlying graph, so we may
                // need to translate some node IDs
                translate_oriented_node_ids(*alignment.mutable_path(), [&](id_t node_id) {
                    handle_t under = overlay.get_underlying_handle(overlay.get_handle(node_id));
                    return make_pair(g.get_id(under), g.get_is_reverse(under));
                });
            }
        }
    }
    else {
        align_internal(alignment, nullptr, g, true, pin_left, 1, true);
    }
}

void QualAdjAligner::align_pinned_multi(Alignment& alignment, vector<Alignment>& alt_alignments, const HandleGraph& g,
                                        bool pin_left, int32_t max_alt_alns) const {
    align_internal(alignment, &alt_alignments, g, true, pin_left, max_alt_alns, true);
}

void QualAdjAligner::align_global_banded(Alignment& alignment, const HandleGraph& g,
                                         int32_t band_padding, bool permissive_banding) const {
    
    int64_t best_score = alignment.sequence().size() * match;
    size_t total_bases = 0;
    g.for_each_handle([&](const handle_t& handle) {
        total_bases += g.get_length(handle);
    });
    int64_t worst_score = (alignment.sequence().size() + total_bases) * -max(max(mismatch, gap_open), gap_extension);
    
    // TODO: put this all into another template somehow?
    
    if (best_score <= numeric_limits<int8_t>::max() && worst_score >= numeric_limits<int8_t>::min()) {
        // We'll fit in int8
        BandedGlobalAligner<int8_t> band_graph(alignment,
                                               g,
                                               band_padding,
                                               permissive_banding,
                                               true);
        
        band_graph.align(score_matrix, nt_table, gap_open, gap_extension);
    } else if (best_score <= numeric_limits<int16_t>::max() && worst_score >= numeric_limits<int16_t>::min()) {
        // We'll fit in int16
        BandedGlobalAligner<int16_t> band_graph(alignment,
                                                g,
                                                band_padding,
                                                permissive_banding,
                                                true);
        
        band_graph.align(score_matrix, nt_table, gap_open, gap_extension);
    } else if (best_score <= numeric_limits<int32_t>::max() && worst_score >= numeric_limits<int32_t>::min()) {
        // We'll fit in int32
        BandedGlobalAligner<int32_t> band_graph(alignment,
                                                g,
                                                band_padding,
                                                permissive_banding,
                                                true);
        
        band_graph.align(score_matrix, nt_table, gap_open, gap_extension);
    } else {
        // Fall back to int64
        BandedGlobalAligner<int64_t> band_graph(alignment,
                                                g,
                                                band_padding,
                                                permissive_banding,
                                                true);
        
        band_graph.align(score_matrix, nt_table, gap_open, gap_extension);
    }
}

void QualAdjAligner::align_global_banded_multi(Alignment& alignment, vector<Alignment>& alt_alignments, const HandleGraph& g,
                                               int32_t max_alt_alns, int32_t band_padding, bool permissive_banding) const {
    
    // We need to figure out what size ints we need to use.
    // Get upper and lower bounds on the scores. TODO: if these overflow int64 we're out of luck
    int64_t best_score = alignment.sequence().size() * match;
    size_t total_bases = 0;
    g.for_each_handle([&](const handle_t& handle) {
        total_bases += g.get_length(handle);
    });
    int64_t worst_score = (alignment.sequence().size() + total_bases) * -max(max(mismatch, gap_open), gap_extension);
    
    if (best_score <= numeric_limits<int8_t>::max() && worst_score >= numeric_limits<int8_t>::min()) {
        // We'll fit in int8
        BandedGlobalAligner<int8_t> band_graph(alignment,
                                               g,
                                               alt_alignments,
                                               max_alt_alns,
                                               band_padding,
                                               permissive_banding,
                                               true);
        
        band_graph.align(score_matrix, nt_table, gap_open, gap_extension);
    } else if (best_score <= numeric_limits<int16_t>::max() && worst_score >= numeric_limits<int16_t>::min()) {
        // We'll fit in int16
        BandedGlobalAligner<int16_t> band_graph(alignment,
                                                g,
                                                alt_alignments,
                                                max_alt_alns,
                                                band_padding,
                                                permissive_banding,
                                                true);
        
        band_graph.align(score_matrix, nt_table, gap_open, gap_extension);
    } else if (best_score <= numeric_limits<int32_t>::max() && worst_score >= numeric_limits<int32_t>::min()) {
        // We'll fit in int32
        BandedGlobalAligner<int32_t> band_graph(alignment,
                                                g,
                                                alt_alignments,
                                                max_alt_alns,
                                                band_padding,
                                                permissive_banding,
                                                true);
        
        band_graph.align(score_matrix, nt_table, gap_open, gap_extension);
    } else {
        // Fall back to int64
        BandedGlobalAligner<int64_t> band_graph(alignment,
                                                g,
                                                alt_alignments,
                                                max_alt_alns,
                                                band_padding,
                                                permissive_banding,
                                                true);
        
        band_graph.align(score_matrix, nt_table, gap_open, gap_extension);
    }
}

void QualAdjAligner::align_xdrop(Alignment& alignment, const HandleGraph& g, const vector<MaximalExactMatch>& mems,
                                 bool reverse_complemented, uint16_t max_gap_length) const
{
    align_xdrop(alignment, g, algorithms::lazier_topological_order(&g), mems, reverse_complemented, max_gap_length);
}

void QualAdjAligner::align_xdrop(Alignment& alignment, const HandleGraph& g, const vector<handle_t>& order,
                                 const vector<MaximalExactMatch>& mems, bool reverse_complemented,
                                 uint16_t max_gap_length) const
{
    // QualAdjXdropAligner manages its own stack, so it can never be threadsafe without being recreated
    // for every alignment, which meshes poorly with its stack implementation. We achieve
    // thread-safety by having one per thread, which makes this method const-ish.
    QualAdjXdropAligner& xdrop = const_cast<QualAdjXdropAligner&>(xdrops[omp_get_thread_num()]);
    xdrop.align(alignment, g, order, mems, reverse_complemented, max_gap_length);
    if (!alignment.has_path() && mems.empty()) {
        // dozeu couldn't find an alignment, probably because it's seeding heuristic failed
        // we'll just fall back on GSSW
        // TODO: This is a bit inconsistent. GSSW gives a full-length bonus at both ends, while
        // dozeu only gives it once.
        align(alignment, g, true);
    }
}

int32_t QualAdjAligner::score_exact_match(const Alignment& aln, size_t read_offset, size_t length) const {
    auto& sequence = aln.sequence();
    auto& base_quality = aln.quality();
    int32_t score = 0;
    for (int32_t i = 0; i < length; i++) {
        // index 5 x 5 score matrices (ACGTN)
        // always have match so that row and column index are same and can combine algebraically
        score += score_matrix[25 * base_quality[read_offset + i] + 6 * nt_table[sequence[read_offset + i]]];
    }
    return score;
}

int32_t QualAdjAligner::score_exact_match(const string& sequence, const string& base_quality) const {
    int32_t score = 0;
    for (int32_t i = 0; i < sequence.length(); i++) {
        // index 5 x 5 score matrices (ACGTN)
        // always have match so that row and column index are same and can combine algebraically
        score += score_matrix[25 * base_quality[i] + 6 * nt_table[sequence[i]]];
    }
    return score;
}


int32_t QualAdjAligner::score_exact_match(string::const_iterator seq_begin, string::const_iterator seq_end,
                                          string::const_iterator base_qual_begin) const {
    int32_t score = 0;
    for (auto seq_iter = seq_begin, qual_iter = base_qual_begin; seq_iter != seq_end; seq_iter++) {
        // index 5 x 5 score matrices (ACGTN)
        // always have match so that row and column index are same and can combine algebraically
        score += score_matrix[25 * (*qual_iter) + 6 * nt_table[*seq_iter]];
        qual_iter++;
    }
    return score;
}

int32_t QualAdjAligner::score_partial_alignment(const Alignment& alignment, const HandleGraph& graph, const Path& path,
                                                string::const_iterator seq_begin) const {
    
    int32_t score = 0;
    string::const_iterator read_pos = seq_begin;
    string::const_iterator qual_pos = alignment.quality().begin() + (seq_begin - alignment.sequence().begin());
    
    for (size_t i = 0; i < path.mapping_size(); i++) {
        const Mapping& mapping = path.mapping(i);
        
        // get the sequence of this node on the proper strand
        string node_seq = graph.get_sequence(graph.get_handle(mapping.position().node_id(),
                                                              mapping.position().is_reverse()));
        
        string::const_iterator ref_pos = node_seq.begin() + mapping.position().offset();
        
        for (size_t j = 0; j < mapping.edit_size(); j++) {
            const Edit& edit = mapping.edit(j);
            
            if (edit.from_length() > 0) {
                if (edit.to_length() > 0) {
                    
                    for (auto siter = read_pos, riter = ref_pos, qiter = qual_pos;
                         siter != read_pos + edit.from_length(); siter++, qiter++, riter++) {
                        score += score_matrix[25 * (*qiter) + 5 * nt_table[*riter] + nt_table[*siter]];
                    }
                    
                    // apply full length bonus
                    if (read_pos == alignment.sequence().begin()) {
                        score += full_length_bonus;
                    }
                    if (read_pos + edit.from_length() == alignment.sequence().end()) {
                        score += full_length_bonus;
                    }
                }
                else {
                    // deletion
                    score -= gap_open + (edit.from_length() - 1) * gap_extension;
                }
            }
            else if (edit.to_length() > 0) {
                // don't score soft clips
                if (read_pos != alignment.sequence().begin() &&
                    read_pos + edit.to_length() != alignment.sequence().end()) {
                    // insert
                    score -= gap_open + (edit.to_length() - 1) * gap_extension;
                }
            }
            
            read_pos += edit.to_length();
            qual_pos += edit.to_length();
            ref_pos += edit.from_length();
        }
    }
    return score;
}

AlignerClient::AlignerClient(double gc_content_estimate) : gc_content_estimate(gc_content_estimate) {
    
    // Adopt the default scoring parameters and make the aligners
    set_alignment_scores(default_score_matrix,
                         default_gap_open, default_gap_extension,
                         default_full_length_bonus);
}

const GSSWAligner* AlignerClient::get_aligner(bool have_qualities) const {
    return (have_qualities && adjust_alignments_for_base_quality) ?
        (GSSWAligner*) get_qual_adj_aligner() :
        (GSSWAligner*) get_regular_aligner();
}

const QualAdjAligner* AlignerClient::get_qual_adj_aligner() const {
    assert(qual_adj_aligner.get() != nullptr);
    return qual_adj_aligner.get();
}

const Aligner* AlignerClient::get_regular_aligner() const {
    assert(regular_aligner.get() != nullptr);
    return regular_aligner.get();
}

int8_t* AlignerClient::parse_matrix(istream& matrix_stream) {
    int8_t* matrix = (int8_t*) malloc(16 * sizeof(int8_t));
    for (size_t i = 0; i < 16; i++) {
        if (!matrix_stream.good()) {
            std::cerr << "error: vg Aligner::parse_matrix requires a 4x4 whitespace separated integer matrix\n";
            throw "";
        }
        int score;
        matrix_stream >> score;
        if (score > 127 || score < -127) {
            std::cerr << "error: vg Aligner::parse_matrix requires values in the range [-127,127]\n";
            throw "";
        }
        matrix[i] = score;
    }
    return matrix;
}

void AlignerClient::set_alignment_scores(int8_t match, int8_t mismatch, int8_t gap_open, int8_t gap_extend, 
                                         int8_t full_length_bonus) {
    
    int8_t* matrix = (int8_t*) malloc(sizeof(int8_t) * 16);
    for (size_t i = 0; i < 16; ++i) {
        if (i % 5 == 0) {
            matrix[i] = match;
        }
        else {
            matrix[i] = -mismatch;
        }
    }
    
    qual_adj_aligner = unique_ptr<QualAdjAligner>(new QualAdjAligner(matrix, gap_open, gap_extend,
                                                                     full_length_bonus, gc_content_estimate));
    regular_aligner = unique_ptr<Aligner>(new Aligner(matrix, gap_open, gap_extend,
                                                      full_length_bonus, gc_content_estimate));
                  
    free(matrix);
}


void AlignerClient::set_alignment_scores(const int8_t* score_matrix, int8_t gap_open, int8_t gap_extend, int8_t full_length_bonus) {
    
    qual_adj_aligner = unique_ptr<QualAdjAligner>(new QualAdjAligner(score_matrix, gap_open, gap_extend,
                                                                     full_length_bonus, gc_content_estimate));
    regular_aligner = unique_ptr<Aligner>(new Aligner(score_matrix, gap_open, gap_extend,
                                                      full_length_bonus, gc_content_estimate));
    
}

void AlignerClient::set_alignment_scores(std::istream& matrix_stream, int8_t gap_open, int8_t gap_extend, int8_t full_length_bonus) {
    int8_t* score_matrix = parse_matrix(matrix_stream);
    qual_adj_aligner = unique_ptr<QualAdjAligner>(new QualAdjAligner(score_matrix, gap_open, gap_extend,
                                                                     full_length_bonus, gc_content_estimate));
    regular_aligner = unique_ptr<Aligner>(new Aligner(score_matrix, gap_open, gap_extend,
                                                      full_length_bonus, gc_content_estimate));
    free(score_matrix);
}
