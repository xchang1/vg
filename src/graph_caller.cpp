#include "graph_caller.hpp"
#include "algorithms/expand_context.hpp"

namespace vg {

GraphCaller::GraphCaller(SnarlCaller& snarl_caller,
                         SnarlManager& snarl_manager,
                         ostream& out_stream) :
    snarl_caller(snarl_caller), snarl_manager(snarl_manager), out_stream(out_stream) {
}

GraphCaller::~GraphCaller() {
}

void GraphCaller::call_top_level_snarls(bool recurse_on_fail) {

    // Used to recurse on children of parents that can't be called
    vector<const Snarl*> snarl_queue;

    // Run the snarl caller on a snarl, and queue up the children if it fails
    auto process_snarl = [&](const Snarl* snarl) {
        bool was_called = call_snarl(*snarl);
        if (!was_called && recurse_on_fail) {
            const vector<const Snarl*>& children = snarl_manager.children_of(snarl);
#pragma omp critical (snarl_queue)
            {
                snarl_queue.insert(snarl_queue.end(), children.begin(), children.end());
            }
        }
    };

    // Start with the top level snarls
    snarl_manager.for_each_top_level_snarl_parallel(process_snarl);

    // Then recurse on any children the snarl caller failed to handle
    while (!snarl_queue.empty()) {
        vector<const Snarl*> cur_queue;
        std::swap(snarl_queue, cur_queue);
#pragma omp parallel for
        for (int i = 0; i < cur_queue.size(); ++i) {
            process_snarl(cur_queue[i]);
        }
    }
  
}

VCFOutputCaller::VCFOutputCaller(const string& sample_name) : sample_name(sample_name) {
}

VCFOutputCaller::~VCFOutputCaller() {
}

string VCFOutputCaller::vcf_header(const PathHandleGraph& graph, const vector<string>& contigs,
                                   const vector<size_t>& contig_length_overrides) const {
    stringstream ss;
    ss << "##fileformat=VCFv4.2" << endl;
    for (int i = 0; i < contigs.size(); ++i) {
        const string& contig = contigs[i];
        path_handle_t path_handle = graph.get_path_handle(contig);
        string path_name = graph.get_path_name(path_handle);
        size_t length;
        if (i < contig_length_overrides.size()) {
            // length override provided
            length = contig_length_overrides[i];
        } else {
            length = 0;
            for (handle_t handle : graph.scan_path(graph.get_path_handle(contig))) {
                length += graph.get_length(handle);
            }
        }
        ss << "##contig=<ID=" << contig << ",length=" << length << ">" << endl;
    }
    return ss.str();
}

void VCFOutputCaller::add_variant(vcflib::Variant& var) const {
#pragma omp critical(add_variant)
    {
        output_variants.push_back(var);
    }
}

void VCFOutputCaller::write_variants(ostream& out_stream) const {
    std::sort(output_variants.begin(), output_variants.end(), [](const vcflib::Variant& v1, const vcflib::Variant& v2) {
            return v1.sequenceName < v2.sequenceName || (v1.sequenceName == v2.sequenceName && v1.position < v2.position);
        });
    for (auto v : output_variants) {
        v.setVariantCallFile(output_vcf);
        out_stream << v << endl;
    }
}

VCFGenotyper::VCFGenotyper(const PathHandleGraph& graph,
                           SnarlCaller& snarl_caller,
                           SnarlManager& snarl_manager,
                           vcflib::VariantCallFile& variant_file,
                           const string& sample_name,
                           const vector<string>& ref_paths,
                           FastaReference* ref_fasta,
                           FastaReference* ins_fasta,
                           ostream& out_stream) :
    GraphCaller(snarl_caller, snarl_manager, out_stream),
    VCFOutputCaller(sample_name),
    graph(graph),
    input_vcf(variant_file),
    traversal_finder(graph, snarl_manager, variant_file, ref_paths, ref_fasta, ins_fasta, snarl_caller.get_skip_allele_fn()) {

    scan_contig_lengths();    
}

VCFGenotyper::~VCFGenotyper() {

}

bool VCFGenotyper::call_snarl(const Snarl& snarl) {

    // could be that our graph is a subgraph of the graph the snarls were computed from
    // so bypass snarls we can't process
    if (!graph.has_node(snarl.start().node_id()) || !graph.has_node(snarl.end().node_id())) {
        return false;
    }
    
    // get our traversals out of the finder
    vector<pair<SnarlTraversal, vector<int>>> alleles;
    vector<vcflib::Variant*> variants;
    std::tie(alleles, variants) = traversal_finder.find_allele_traversals(snarl);

    if (!alleles.empty()) {

        // hmm, maybe find a way not to copy?
        vector<SnarlTraversal> travs;
        travs.reserve(alleles.size());
        for (const auto& ta : alleles) {
            travs.push_back(ta.first);
        }

        // find the reference traversal
        // todo: is it the reference always first?
        int ref_trav_idx = -1;
        for (int i = 0; i < alleles.size() && ref_trav_idx < 0; ++i) {
            if (std::all_of(alleles[i].second.begin(), alleles[i].second.end(), [](int x) {return x == 0;})) {
                ref_trav_idx = i;
            }
        }

        // use our support caller to choose our genotype (int traversal coordinates)
        vector<int> trav_genotype = snarl_caller.genotype(snarl, travs, ref_trav_idx, 2);

        assert(trav_genotype.empty() || trav_genotype.size() == 2);

        // map our genotype back to the vcf
        for (int i = 0; i < variants.size(); ++i) {
            vector<int> vcf_alleles;
            set<int> used_vcf_alleles;
            string vcf_genotype;
            vector<SnarlTraversal> vcf_traversals(variants[i]->alleles.size());            
            if (trav_genotype.empty()) {
                vcf_genotype = "./.";
            } else {
                // map our traversal genotype to a vcf variant genotype
                // using the information out of the traversal finder
                for (int j = 0; j < trav_genotype.size(); ++j) {
                    int trav_allele = trav_genotype[j];
                    int vcf_allele = alleles[trav_allele].second[i];
                    vcf_genotype += std::to_string(vcf_allele);
                    if (j < trav_genotype.size() - 1) {
                        vcf_genotype += "/";
                    }
                    vcf_alleles.push_back(vcf_allele);
                    used_vcf_alleles.insert(vcf_allele);
                    vcf_traversals[vcf_allele] = travs[trav_allele];
                }
                // add traversals that correspond to vcf genotypes that are not
                // present in the traversal_genotypes
                for (int j = 0; j < travs.size(); ++j) {
                    int vcf_allele = alleles[j].second[i];
                    if (!used_vcf_alleles.count(vcf_allele)) {
                        vcf_traversals[vcf_allele] = travs[j];
                        used_vcf_alleles.insert(vcf_allele);
                    }
                }
            }
            // create an output variant from the input one
            vcflib::Variant out_variant;
            out_variant.sequenceName = variants[i]->sequenceName;
            out_variant.position = variants[i]->position;
            out_variant.id = variants[i]->id;
            out_variant.ref = variants[i]->ref;
            out_variant.alt = variants[i]->alt;
            out_variant.alleles = variants[i]->alleles;
            out_variant.filter = "PASS";
            out_variant.updateAlleleIndexes();

            // add the genotype
            out_variant.format.push_back("GT");
            auto& genotype_vector = out_variant.samples[sample_name]["GT"];
            genotype_vector.push_back(vcf_genotype);

            // add some info
            snarl_caller.update_vcf_info(snarl, vcf_traversals, vcf_alleles, sample_name, out_variant);

            // print the variant
            add_variant(out_variant);
        }
        return true;
    }
    
    return false;

}

string VCFGenotyper::vcf_header(const PathHandleGraph& graph, const vector<string>& ref_paths,
                                const vector<size_t>& contig_length_overrides) const {
    assert(contig_length_overrides.empty()); // using this override makes no sense

    // get the contig length overrides from the VCF
    vector<size_t> vcf_contig_lengths;
    auto length_map = scan_contig_lengths();
    for (int i = 0; i < ref_paths.size(); ++i) {
        vcf_contig_lengths.push_back(length_map[ref_paths[i]]);
    }
    
    string header = VCFOutputCaller::vcf_header(graph, ref_paths, vcf_contig_lengths);
    header += "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n";
    snarl_caller.update_vcf_header(header);
    header += "##FILTER=<ID=PASS,Description=\"All filters passed\">\n";
    header += "##SAMPLE=<ID=" + sample_name + ">\n";
    header += "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + sample_name;
    assert(output_vcf.openForOutput(header));
    header += "\n";
    return header;
}

unordered_map<string, size_t> VCFGenotyper::scan_contig_lengths() const {

    unordered_map<string, size_t> ref_lengths;
    
    // copied from dumpContigsFromHeader.cpp in vcflib
    vector<string> headerLines = split(input_vcf.header, "\n");
    for(vector<string>::iterator it = headerLines.begin(); it != headerLines.end(); it++) {
        if((*it).substr(0,8) == "##contig"){
            string contigInfo = (*it).substr(10, (*it).length() -11);
            vector<string> info = split(contigInfo, ",");
            string id;
            int64_t length = -1;
            for(vector<string>::iterator sub = info.begin(); sub != info.end(); sub++) {
                vector<string> subfield = split((*sub), "=");
                if(subfield[0] == "ID"){
                    id = subfield[1];
                }
                if(subfield[0] == "length"){
                    length = parse<int>(subfield[1]);
                }
            }
            if (!id.empty() && length >= 0) {
                ref_lengths[id] = length;
            }
        }
    }

    return ref_lengths;
}


LegacyCaller::LegacyCaller(const PathPositionHandleGraph& graph,
                           SupportBasedSnarlCaller& snarl_caller,
                           SnarlManager& snarl_manager,
                           const string& sample_name,
                           const vector<string>& ref_paths,
                           const vector<size_t>& ref_path_offsets) :
    GraphCaller(snarl_caller, snarl_manager, out_stream),
    VCFOutputCaller(sample_name),
    graph(graph),
    ref_paths(ref_paths) {

    for (int i = 0; i < ref_paths.size(); ++i) {
        ref_offsets[ref_paths[i]] = i < ref_path_offsets.size() ? ref_path_offsets[i] : 0;
    }
    
    is_vg = dynamic_cast<const VG*>(&graph) != nullptr;
    if (is_vg) {
        // our graph is in vg format.  we index the paths and make a traversal finder just
        // like in the old call code
        for (auto ref_path : ref_paths) {
            path_indexes.push_back(new PathIndex(graph, ref_path));
        }
        // map snarl to the first reference path that spans it
        function<PathIndex*(const Snarl&)> get_path_index = [&](const Snarl& site) -> PathIndex* {
            return find_index(site, path_indexes).second;
        };
        // initialize our traversal finder
        traversal_finder = new RepresentativeTraversalFinder(graph, snarl_manager,
                                                             max_search_depth,
                                                             max_search_width,
                                                             max_bubble_paths,
                                                             0,
                                                             0,
                                                             get_path_index,
                                                             [&](id_t id) { return snarl_caller.get_min_node_support(id);},
                                                             [&](edge_t edge) { return snarl_caller.get_edge_support(edge);});

    } else {
        // our graph is not in vg format.  we will make graphs for each site as needed and work with those
        traversal_finder = nullptr;
    }
}

LegacyCaller::~LegacyCaller() {
    delete traversal_finder;
    for (PathIndex* path_index : path_indexes) {
        delete path_index;
    }
}

bool LegacyCaller::call_snarl(const Snarl& snarl) {

    // if we can't handle the snarl, then the GraphCaller framework will recurse on its children
    if (!is_traversable(snarl)) {
        return false;
    }
           
    RepresentativeTraversalFinder* rep_trav_finder;
    vector<PathIndex*> site_path_indexes;
    function<PathIndex*(const Snarl&)> get_path_index;
    VG vg_graph;
    SupportBasedSnarlCaller& support_caller = dynamic_cast<SupportBasedSnarlCaller&>(snarl_caller);
    
    if (is_vg) {
        // our graph is in VG format, so we've sorted this out in the constructor
        rep_trav_finder = traversal_finder;
        get_path_index = [&](const Snarl& site) {
            return find_index(site, path_indexes).second;
        };
        
    } else {
        // our graph isn't in VG format.  we are using a (hopefully temporary) workaround
        // of converting the subgraph into VG.
        pair<unordered_set<id_t>, unordered_set<edge_t> > contents = snarl_manager.deep_contents(&snarl, graph, true);
        size_t total_snarl_length = 0;
        for (auto node_id : contents.first) {
            handle_t new_handle = vg_graph.create_handle(graph.get_sequence(graph.get_handle(node_id)), node_id);
            if (node_id != snarl.start().node_id() && node_id != snarl.end().node_id()) {
                total_snarl_length += vg_graph.get_length(new_handle);
            }
        }
        for (auto edge : contents.second) {
            vg_graph.create_edge(vg_graph.get_handle(graph.get_id(edge.first), vg_graph.get_is_reverse(edge.first)),
                                 vg_graph.get_handle(graph.get_id(edge.second), vg_graph.get_is_reverse(edge.second)));
            total_snarl_length += 1;
        }
        // add the paths to the subgraph
        algorithms::expand_context_with_paths(&graph, &vg_graph, 1);
        // and index them
        for (auto& ref_path : ref_paths) {
            if (vg_graph.has_path(ref_path)) {
                site_path_indexes.push_back(new PathIndex(vg_graph, ref_path));
            } else {
                site_path_indexes.push_back(nullptr);
            }
        }
        get_path_index = [&](const Snarl& site) -> PathIndex* {
            return find_index(site, site_path_indexes).second;
        };
        // determine the support threshold for the traversal finder.  if we're using average
        // support, then we don't use any (set to 0), other wise, use the minimum support for a call
        SupportBasedSnarlCaller& support_caller = dynamic_cast<SupportBasedSnarlCaller&>(snarl_caller);
        size_t threshold = support_caller.get_average_traversal_support_switch_threshold();
        double support_cutoff = total_snarl_length <= threshold ? support_caller.get_min_total_support_for_call() : 0;
        rep_trav_finder = new RepresentativeTraversalFinder(vg_graph, snarl_manager,
                                                            max_search_depth,
                                                            max_search_width,
                                                            max_bubble_paths,
                                                            support_cutoff,
                                                            support_cutoff,
                                                            get_path_index,
                                                            [&](id_t id) { return support_caller.get_min_node_support(id);},
                                                            // note: because our traversal finder and support caller have
                                                            // different graphs, they can't share edge handles
                                                            [&](edge_t edge) { return support_caller.get_edge_support(
                                                                    vg_graph.get_id(edge.first), vg_graph.get_is_reverse(edge.first),
                                                                    vg_graph.get_id(edge.second), vg_graph.get_is_reverse(edge.second));});
                                                            
    }

    PathIndex* path_index = get_path_index(snarl);
    if (path_index != nullptr) {
        string path_name = find_index(snarl, is_vg ? path_indexes : site_path_indexes).first;

        // orient the snarl along the reference path
        if (get_ref_position(snarl, path_name).second == true) {
            snarl_manager.flip(&snarl);
        }

        // recursively genotype the site beginning here at the top level snarl
        vector<SnarlTraversal> called_traversals;
        // these integers map the called traversals to their positions in the list of all traversals
        // of the top level snarl.  
        vector<int> genotype;
        std::tie(called_traversals, genotype) = top_down_genotype(snarl, *rep_trav_finder, 2);
    
        if (!called_traversals.empty()) {
            // regenotype our top-level traversals now that we know they aren't nested, and we have a
            // good idea of all the sizes
            std::tie(called_traversals, genotype) = re_genotype(snarl, *rep_trav_finder, called_traversals, genotype, 2);

            // emit our vcf variant
            emit_variant(snarl, *rep_trav_finder, called_traversals, genotype, path_name);
        }
    }        
    if (!is_vg) {
        // delete the temporary vg subgraph and traversal finder we created for this snarl
        delete rep_trav_finder;
        for (PathIndex* path_index : site_path_indexes) {
            delete path_index;
        }
    }

    return true;
}

string LegacyCaller::vcf_header(const PathHandleGraph& graph, const vector<string>& ref_paths,
                                const vector<size_t>& contig_length_overrides) const {
    string header = VCFOutputCaller::vcf_header(graph, ref_paths, contig_length_overrides);
    header += "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n";
    snarl_caller.update_vcf_header(header);
    header += "##FILTER=<ID=PASS,Description=\"All filters passed\">\n";
    header += "##SAMPLE=<ID=" + sample_name + ">\n";
    header += "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + sample_name;
    assert(output_vcf.openForOutput(header));
    header += "\n";
    return header;
}

pair<vector<SnarlTraversal>, vector<int>> LegacyCaller::top_down_genotype(const Snarl& snarl, TraversalFinder& trav_finder, int ploidy) const {

    // get the traversals through the site
    vector<SnarlTraversal> traversals = trav_finder.find_traversals(snarl);

    // use our support caller to choose our genotype
    vector<int> trav_genotype = snarl_caller.genotype(snarl, traversals, 0, ploidy);
    if (trav_genotype.empty()) {
        return make_pair(vector<SnarlTraversal>(), vector<int>());
    }

    assert(trav_genotype.size() == ploidy);

    vector<SnarlTraversal> called_travs(ploidy);

    // do we have two paths going through a given traversal?  This is handled
    // as a special case below
    bool hom = trav_genotype.size() == 2 && trav_genotype[0] == trav_genotype[1];
    
    for (int i = 0; i < trav_genotype.size() && (!hom || i < 1); ++i) {
        int allele = trav_genotype[i];
        const SnarlTraversal& traversal = traversals[allele];
        Visit prev_end;
        for (int j = 0; j < traversal.visit_size(); ++j) {
            if (traversal.visit(j).node_id() > 0) {
                *called_travs[i].add_visit() = traversal.visit(j);
                if (hom && i == 0) {
                    *called_travs[1].add_visit() = traversal.visit(j);
                }
            } else {
                // recursively determine the traversal
                const Snarl* into_snarl = snarl_manager.into_which_snarl(traversal.visit(j));
                bool flipped = traversal.visit(j).backward();
                if (flipped) {
                    // we're always processing our snarl from start to end, so make sure
                    // it lines up with the parent (note that we've oriented the root along the ref path)
                    snarl_manager.flip(into_snarl);
                }
                vector<SnarlTraversal> child_genotype = top_down_genotype(*into_snarl,
                                                                          trav_finder, hom ? 2: 1).first;                
                if (child_genotype.empty()) {
                    return make_pair(vector<SnarlTraversal>(), vector<int>());
                }
                bool back_to_back = j > 0 && traversal.visit(j - 1).node_id() == 0 && prev_end == into_snarl->start();

                for (int k = back_to_back ? 1 : 0; k < child_genotype[0].visit_size(); ++k) {
                    *called_travs[i].add_visit() = child_genotype[0].visit(k);
                }
                if (hom) {
                    assert(child_genotype.size() == 2 && i == 0);
                    for (int k = back_to_back ? 1 : 0; k < child_genotype[1].visit_size(); ++k) {
                        *called_travs[1].add_visit() = child_genotype[1].visit(k);
                    }
                }
                prev_end = into_snarl->end();
                if (flipped) {
                    // leave our snarl like we found it
                    snarl_manager.flip(into_snarl);
                }
            }
        }
    }

    return make_pair(called_travs, trav_genotype);
}

SnarlTraversal LegacyCaller::get_reference_traversal(const Snarl& snarl, TraversalFinder& trav_finder) const {

    // get the ref traversal through the site
    // todo: don't avoid so many traversal recomputations
    SnarlTraversal traversal = trav_finder.find_traversals(snarl)[0];
    SnarlTraversal out_traversal;

    Visit prev_end;
    for (int i = 0; i < traversal.visit_size(); ++i) {
        const Visit& visit = traversal.visit(i);
        if (visit.node_id() != 0) {
            *out_traversal.add_visit() = visit;
        } else {
            const Snarl* into_snarl = snarl_manager.into_which_snarl(visit);
            if (visit.backward()) {
                snarl_manager.flip(into_snarl);
            }
            bool back_to_back = i > 0 && traversal.visit(i - 1).node_id() == 0 && prev_end == into_snarl->start();

            SnarlTraversal child_ref = get_reference_traversal(*into_snarl, trav_finder);
            for (int j = back_to_back ? 1 : 0; j < child_ref.visit_size(); ++j) {
                *out_traversal.add_visit() = child_ref.visit(j);
            }
            prev_end = into_snarl->end();
            if (visit.backward()) {
                // leave our snarl like we found it
                snarl_manager.flip(into_snarl);
            }
        }
    }
    return out_traversal;    
}

pair<vector<SnarlTraversal>, vector<int>> LegacyCaller::re_genotype(const Snarl& snarl, TraversalFinder& trav_finder,
                                                                    const vector<SnarlTraversal>& in_traversals,
                                                                    const vector<int>& in_genotype,
                                                                    int ploidy) const {
    assert(in_traversals.size() == in_genotype.size());

    // create a set of unique traversal candidates that must include the reference first
    vector<SnarlTraversal> rg_traversals;
    // add our reference traversal to the front
    for (int i = 0; i < in_traversals.size() && !rg_traversals.empty(); ++i) {
        if (in_genotype[i] == 0) {
            rg_traversals.push_back(in_traversals[i]);
        }
    }
    if (rg_traversals.empty()) {
        rg_traversals.push_back(get_reference_traversal(snarl, trav_finder));
    }
    set<int> gt_set = {0};
    for (int i = 0; i < in_traversals.size(); ++i) {
        if (!gt_set.count(in_genotype[i])) {
            rg_traversals.push_back(in_traversals[i]);
            gt_set.insert(in_genotype[i]);
        }
    }

    // re-genotype the candidates
    vector<int> rg_genotype = snarl_caller.genotype(snarl, rg_traversals, 0, ploidy);

    // convert our output to something that emit_variant() will understand
    // todo: this is needlessly inefficient and should be streamlined to operate
    // directly on rg_traversals/rg_genotype once testing done
    vector<SnarlTraversal> out_traversals;
    vector<int> out_genotype;

    for (int i = 0; i < rg_genotype.size(); ++i) {
        out_traversals.push_back(rg_traversals[rg_genotype[i]]);
        out_genotype.push_back(rg_genotype[i]);
    }

    return make_pair(out_traversals, out_genotype);
}

void LegacyCaller::emit_variant(const Snarl& snarl, TraversalFinder& trav_finder, const vector<SnarlTraversal>& called_traversals,
                                const vector<int>& genotype, const string& ref_path_name) const {
    
    // convert traversal to string
    function<string(const SnarlTraversal&)> trav_string = [&](const SnarlTraversal& trav) {
        string seq;
        for (int i = 0; i < trav.visit_size(); ++i) {
            seq += graph.get_sequence(graph.get_handle(trav.visit(i).node_id(), trav.visit(i).backward()));
        }
        return seq;
    };

    vcflib::Variant out_variant;

    // when calling alt/alt, the reference traversal doesn't end up in called_traversals.
    // this should get changed, but in the meantime we add it back here (as we need it for
    // the VCF output)
    // udpate: the reference traversal will be there when re-genotyping, but we can leave this logic
    // in case we want to ever add an option to toggle this.
    vector<SnarlTraversal> site_traversals;
    vector<int> site_genotype;
    for (int i = 0; i < genotype.size(); ++i) {
        if (genotype[i] == 0) {
            site_traversals.push_back(called_traversals[i]);
            break;
        }
    }
    if (site_traversals.empty()) {
        site_traversals.push_back(get_reference_traversal(snarl, trav_finder));
    }
    out_variant.ref = trav_string(site_traversals[0]);
    
    // deduplicate alleles and compute the site traversals and genotype
    map<string, int> allele_to_gt;    
    allele_to_gt[out_variant.ref] = 0;    
    for (int i = 0; i < genotype.size(); ++i) {
        if (genotype[i] == 0) {
            site_genotype.push_back(0);
        } else {
            string allele_string = trav_string(called_traversals[i]);
            if (allele_to_gt.count(allele_string)) {
                site_genotype.push_back(allele_to_gt[allele_string]);
            } else {
                site_traversals.push_back(called_traversals[i]);
                site_genotype.push_back(allele_to_gt.size());
                allele_to_gt[allele_string] = site_genotype.back();
            }
        }
    }

    out_variant.alt.resize(allele_to_gt.size() - 1);
    out_variant.alleles.resize(allele_to_gt.size());
    for (auto& allele_gt : allele_to_gt) {
        if (allele_gt.second > 0) {
            out_variant.alt[allele_gt.second - 1] = allele_gt.first;
        }
        out_variant.alleles[allele_gt.second] = allele_gt.first;
    }

    // fill out the rest of the variant
    out_variant.sequenceName = ref_path_name;
    // +1 to convert to 1-based VCF
    out_variant.position = get_ref_position(snarl, ref_path_name).first + ref_offsets.find(ref_path_name)->second + 1; 
    out_variant.id = std::to_string(snarl.start().node_id()) + "_" + std::to_string(snarl.end().node_id());
    out_variant.filter = "PASS";
    out_variant.updateAlleleIndexes();

    // add the genotype
    out_variant.format.push_back("GT");
    auto& genotype_vector = out_variant.samples[sample_name]["GT"];
    
    stringstream vcf_gt;
    for (int i = 0; i < site_genotype.size(); ++i) {
        vcf_gt << site_genotype[i];
        if (i != site_genotype.size() - 1) {
            vcf_gt << "/";
        }
    }
    genotype_vector.push_back(vcf_gt.str());

    // clean up the alleles to not have so man common prefixes
    flatten_common_allele_ends(out_variant, true);
    flatten_common_allele_ends(out_variant, false);

    // add some support info
    snarl_caller.update_vcf_info(snarl, site_traversals, site_genotype, sample_name, out_variant);

    if (!out_variant.alt.empty()) {
        add_variant(out_variant);
    }
}

bool LegacyCaller::is_traversable(const Snarl& snarl) {
    // we need this to be true all the way down to use the RepresentativeTraversalFinder on our snarl.
    bool ret = snarl.start_end_reachable() && snarl.directed_acyclic_net_graph() &&
       graph.has_node(snarl.start().node_id()) && graph.has_node(snarl.end().node_id());
    if (ret == true) {
        const vector<const Snarl*>& children = snarl_manager.children_of(&snarl);
        for (int i = 0; i < children.size() && ret; ++i) {
            ret = is_traversable(*children[i]);
        }
    }
    return ret;
}

pair<string, PathIndex*> LegacyCaller::find_index(const Snarl& snarl, const vector<PathIndex*> path_indexes) const {
    assert(path_indexes.size() == ref_paths.size());
    for (int i = 0; i < path_indexes.size(); ++i) {
        PathIndex* path_index = path_indexes[i];
        if (path_index != nullptr &&
            path_index->by_id.count(snarl.start().node_id()) &&
            path_index->by_id.count(snarl.end().node_id())) {
            // This path threads through this site
            return make_pair(ref_paths[i], path_index);
        }
    }
    return make_pair("", nullptr);
}

pair<size_t, bool> LegacyCaller::get_ref_position(const Snarl& snarl, const string& ref_path_name) const {
    path_handle_t path_handle = graph.get_path_handle(ref_path_name);

    handle_t start_handle = graph.get_handle(snarl.start().node_id(), snarl.start().backward());
    map<size_t, step_handle_t> start_steps;
    graph.for_each_step_on_handle(start_handle, [&](step_handle_t step) {
            if (graph.get_path_handle_of_step(step) == path_handle) {
                start_steps[graph.get_position_of_step(step)] = step;
            }
        });

    handle_t end_handle = graph.get_handle(snarl.end().node_id(), snarl.end().backward());
    map<size_t, step_handle_t> end_steps;
    graph.for_each_step_on_handle(end_handle, [&](step_handle_t step) {
            if (graph.get_path_handle_of_step(step) == path_handle) {
                end_steps[graph.get_position_of_step(step)] = step;
            }
        });

    assert(start_steps.size() > 0 && end_steps.size() > 0);
    step_handle_t start_step = start_steps.begin()->second;
    step_handle_t end_step = end_steps.begin()->second;
    bool scan_backward = graph.get_is_reverse(graph.get_handle_of_step(start_step));

    // if we're on a cycle, we keep our start step and find the end step by scanning the path
    if (start_steps.size() > 1 || end_steps.size() > 1) {
        bool found_end = false;
        if (scan_backward) {
            for (step_handle_t cur_step = start_step; graph.has_previous_step(end_step) && !found_end;
                 cur_step = graph.get_previous_step(cur_step)) {
                if (graph.get_handle_of_step(cur_step) == end_handle) {
                    end_step = cur_step;
                    found_end = true;
                }
            }
            assert(found_end);
        } else {
            for (step_handle_t cur_step = start_step; graph.has_next_step(end_step) && !found_end;
                 cur_step = graph.get_next_step(cur_step)) {
                if (graph.get_handle_of_step(cur_step) == end_handle) {
                    end_step = cur_step;
                    found_end = true;
                }
            }
            assert(found_end);
        }
    }
    
    size_t start_position = start_steps.begin()->first;
    size_t end_position = end_step == end_steps.begin()->second ? end_steps.begin()->first : graph.get_position_of_step(end_step);
    bool backward = end_position < start_position;

    return make_pair(backward ? end_position : start_position, backward);
}

void LegacyCaller::flatten_common_allele_ends(vcflib::Variant& variant, bool backward) const {
    if (variant.alt.size() == 0) {
        return;
    }
    size_t min_len = variant.alleles[0].length();
    for (int i = 1; i < variant.alleles.size(); ++i) {
        min_len = std::min(min_len, variant.alleles[i].length());
    }
    // want to leave at least one in the reference position
    if (min_len > 0) {
        --min_len;
    }

    bool match = true;
    int shared_prefix_len = 0;
    for (int i = 0; i < min_len && match; ++i) {
        char c1 = std::toupper(variant.alleles[0][!backward ? i : variant.alleles[0].length() - 1 - i]);
        for (int j = 1; j < variant.alleles.size() && match; ++j) {
            char c2 = std::toupper(variant.alleles[j][!backward ? i : variant.alleles[j].length() - 1 - i]);
            match = c1 == c2;
        }
        if (match) {
            ++shared_prefix_len;
        }
    }

    if (!backward) {
        variant.position += shared_prefix_len;
    }
    for (int i = 0; i < variant.alleles.size(); ++i) {
        if (!backward) {
            variant.alleles[i] = variant.alleles[i].substr(shared_prefix_len);
        } else {
            variant.alleles[i] = variant.alleles[i].substr(0, variant.alleles[i].length() - shared_prefix_len);
        }
        if (i == 0) {
            variant.ref = variant.alleles[i];
        } else {
            variant.alt[i - 1] = variant.alleles[i];
        }
    }
}

}

