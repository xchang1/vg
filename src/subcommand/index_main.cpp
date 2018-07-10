// index.cpp: define the "vg index" subcommand, which makes xg, GCSA2, GBWT, and RocksDB indexes

#include <omp.h>
#include <unistd.h>
#include <getopt.h>

#include <random>
#include <string>
#include <vector>

#include "subcommand.hpp"

#include "../vg.hpp"
#include "../index.hpp"
#include "../stream.hpp"
#include "../vg_set.hpp"
#include "../utility.hpp"
#include "../region.hpp"

#include <gcsa/gcsa.h>
#include <gcsa/algorithms.h>
#include <gbwt/dynamic_gbwt.h>
#include <gbwt/variants.h>

using namespace std;
using namespace vg;
using namespace vg::subcommand;

void help_index(char** argv) {
    cerr << "usage: " << argv[0] << " index [options] <graph1.vg> [graph2.vg ...]" << endl
         << "Creates an index on the specified graph or graphs. All graphs indexed must " << endl
         << "already be in a joint ID space." << endl
         << "general options:" << endl
         << "    -b, --temp-dir DIR     use DIR for temporary files" << endl
         << "    -t, --threads N        number of threads to use" << endl
         << "    -p, --progress         show progress" << endl
         << "xg options:" << endl
         << "    -x, --xg-name FILE     use this file to store a succinct, queryable version of the graph(s)" << endl
         << "    -F, --thread-db FILE   read thread database from FILE (may repeat)" << endl
         << "gbwt options:" << endl
         << "    -v, --vcf-phasing FILE generate threads from the haplotypes in the VCF file FILE" << endl
         << "    -T, --store-threads    generate threads from the embedded paths" << endl
         << "    -M, --store-gam FILE   generate threads from the alignments in FILE (many allowed)" << endl
         << "    -G, --gbwt-name FILE   store the threads as GBWT in FILE" << endl
         << "    -H, --write-haps FILE  store the threads as sequences in FILE" << endl
         << "    -F, --thread-db FILE   write thread database to FILE" << endl
         << "    -P, --force-phasing    replace unphased genotypes with randomly phased ones" << endl
         << "    -B, --batch-size N     number of samples per batch (default 200)" << endl
         << "    -R, --range X..Y       process samples X to Y (inclusive)" << endl
         << "    -r, --rename V=P       rename contig V in the VCFs to path P in the graph (may repeat)" << endl
         << "    -I, --region C:S-E     operate on only the given 1-based region of the given VCF contig (may repeat)" << endl
         << "    -E, --exclude SAMPLE   exclude any samples with the given name from haplotype indexing" << endl
         << "gcsa options:" << endl
         << "    -g, --gcsa-out FILE    output a GCSA2 index instead of a rocksdb index" << endl
         << "    -i, --dbg-in FILE      use kmers from FILE instead of input VG (may repeat)" << endl
         << "    -f, --mapping FILE     use this node mapping in GCSA2 construction" << endl
         << "    -k, --kmer-size N      index kmers of size N in the graph (default " << gcsa::Key::MAX_LENGTH << ")" << endl
         << "    -X, --doubling-steps N use this number of doubling steps for GCSA2 construction (default " << gcsa::ConstructionParameters::DOUBLING_STEPS << ")" << endl
         << "    -Z, --size-limit N     limit temporary disk space usage to N gigabytes (default " << gcsa::ConstructionParameters::SIZE_LIMIT << ")" << endl
         << "    -V, --verify-index     validate the GCSA2 index using the input kmers (important for testing)" << endl
         << "rocksdb options:" << endl
         << "    -d, --db-name  <X>     store the RocksDB index in <X>" << endl
         << "    -m, --store-mappings   input is .gam format, store the mappings in alignments by node" << endl
         << "    -a, --store-alignments input is .gam format, store the alignments by node" << endl
         << "    -A, --dump-alignments  graph contains alignments, output them in sorted order" << endl
         << "    -N, --node-alignments  input is (ideally, sorted) .gam format," << endl
         << "                           cross reference nodes by alignment traversals" << endl
         << "    -D, --dump             print the contents of the db to stdout" << endl
         << "    -C, --compact          compact the index into a single level (improves performance)" << endl;
}

// Convert gbwt::node_type to ThreadMapping.
xg::XG::ThreadMapping gbwt_to_thread_mapping(gbwt::node_type node) {
    xg::XG::ThreadMapping thread_mapping = { (int64_t)(gbwt::Node::id(node)), gbwt::Node::is_reverse(node) };
    return thread_mapping;
}

// Convert Path to a GBWT path.
gbwt::vector_type path_to_gbwt(const Path& path) {
    gbwt::vector_type result(path.mapping_size());
    for (size_t i = 0; i < result.size(); i++) {
        result[i] = gbwt::Node::encode(path.mapping(i).position().node_id(), path.mapping(i).position().is_reverse());
    }
    return result;
}

// Find all predecessor nodes of the path, ignoring self-loops.
gbwt::vector_type predecessors(const xg::XG& xg_index, const Path& path) {
    gbwt::vector_type result;
    if (path.mapping_size() == 0) {
        return result;
    }

    vg::id_t first_node = path.mapping(0).position().node_id();
    bool is_reverse = path.mapping(0).position().is_reverse();
    auto pred_edges = (is_reverse ? xg_index.edges_on_end(first_node) : xg_index.edges_on_start(first_node));
    for (auto& edge : pred_edges) {
        if (edge.from() == edge.to()) {
            continue; // Self-loop.
        }
        if (edge.from() == first_node) {  // Reverse the edge if it is from this node.
            result.push_back(gbwt::Node::encode(edge.to(), !(edge.to_end())));
        } else {
            result.push_back(gbwt::Node::encode(edge.from(), edge.from_start()));
        }
    }

    return result;
}

// Thread database files written by vg index -G and read by vg index -x.
// These should probably be in thread_database.cpp or something like that.
void write_thread_db(const std::string& filename, const std::vector<std::string>& thread_names, size_t haplotype_count);
void read_thread_db(const std::vector<std::string>& filenames, std::vector<std::string>& thread_names, size_t& haplotype_count);

int main_index(int argc, char** argv) {

    if (argc == 2) {
        help_index(argv);
        return 1;
    }

    // Which indexes to build.
    bool build_xg = false, build_gbwt = false, write_threads = false, build_gpbwt = false, build_gcsa = false, build_rocksdb = false;

    // Files we should read.
    string vcf_name, mapping_name;
    vector<string> thread_db_names;
    vector<string> dbg_names;

    // Files we should write.
    string xg_name, gbwt_name, threads_name, gcsa_name, rocksdb_name;

    // General
    bool show_progress = false;

    // GBWT
    bool index_haplotypes = false, index_paths = false, index_gam = false;
    vector<string> gam_file_names;
    bool force_phasing = false;
    size_t samples_in_batch = 200; // Samples per batch.
    std::pair<size_t, size_t> sample_range(0, ~(size_t)0); // The semiopen range of samples to process.
    map<string, string> path_to_vcf; // Path name conversion from --rename.
    map<string, pair<size_t, size_t>> regions; // Region restrictions for contigs, in VCF name space, as 0-based exclusive-end ranges.
    unordered_set<string> excluded_samples; // Excluded sample names from --exclude.

    // GCSA
    gcsa::size_type kmer_size = gcsa::Key::MAX_LENGTH;
    gcsa::ConstructionParameters params;
    bool verify_gcsa = false;

    // RocksDB
    bool dump_index = false;
    bool store_alignments = false;
    bool store_node_alignments = false;
    bool store_mappings = false;
    bool dump_alignments = false;

    // Unused?
    bool compact = false;

    int c;
    optind = 2; // force optind past command positional argument
    while (true) {
        static struct option long_options[] =
        {
            // General
            {"temp-dir", required_argument, 0, 'b'},
            {"threads", required_argument, 0, 't'},
            {"progress",  no_argument, 0, 'p'},

            // XG
            {"xg-name", required_argument, 0, 'x'},
            {"thread-db", required_argument, 0, 'F'},

            // GBWT
            {"vcf-phasing", required_argument, 0, 'v'},
            {"store-threads", no_argument, 0, 'T'},
            {"store-gam", required_argument, 0, 'M'},
            {"gbwt-name", required_argument, 0, 'G'},
            {"write-haps", required_argument, 0, 'H'},
            {"force-phasing", no_argument, 0, 'P'},
            {"batch-size", required_argument, 0, 'B'},
            {"range", required_argument, 0, 'R'},
            {"rename", required_argument, 0, 'r'},
            {"region", required_argument, 0, 'I'},
            {"exclude", required_argument, 0, 'E'},

            // GCSA
            {"gcsa-name", required_argument, 0, 'g'},
            {"dbg-in", required_argument, 0, 'i'},
            {"mapping", required_argument, 0, 'f'},
            {"kmer-size", required_argument, 0, 'k'},
            {"doubling-steps", required_argument, 0, 'X'},
            {"size-limit", required_argument, 0, 'Z'},
            {"verify-index", no_argument, 0, 'V'},

            // RocksDB
            {"db-name", required_argument, 0, 'd'},
            {"store-mappings", no_argument, 0, 'm'},
            {"store-alignments", no_argument, 0, 'a'},
            {"dump-alignments", no_argument, 0, 'A'},
            {"node-alignments", no_argument, 0, 'N'},
            {"dump", no_argument, 0, 'D'},
            {"compact", no_argument, 0, 'C'},
            {0, 0, 0, 0}
        };

        int option_index = 0;
        c = getopt_long (argc, argv, "b:t:px:F:v:TG:H:PB:R:r:I:E:g:i:f:k:X:Z:Vd:maANDP:CM:h",
                long_options, &option_index);

        // Detect the end of the options.
        if (c == -1)
            break;

        switch (c)
        {
        // General
        case 'b':
            temp_file::set_dir(optarg);
            break;
        case 't':
            omp_set_num_threads(atoi(optarg));
            break;
        case 'p':
            show_progress = true;
            break;

        // XG
        case 'x':
            build_xg = true;
            xg_name = optarg;
            break;
        case 'F':
            thread_db_names.push_back(optarg);
            break;

        // GBWT
        case 'v':
            index_haplotypes = true;
            build_xg = true;
            vcf_name = optarg;
            break;
        case 'T':
            index_paths = true;
            build_xg = true;
            break;
        case 'M':
            index_gam = true;
            build_gbwt = true;
            gam_file_names.push_back(optarg);
            break;
        case 'G':
            build_gbwt = true;
            gbwt_name = optarg;
            break;
        case 'H':
            write_threads = true;
            threads_name = optarg;
            break;
        case 'P':
            force_phasing = true;
            break;
        case 'B':
            samples_in_batch = std::max(std::stoul(optarg), 1ul);
            break;
        case 'R':
            {
                // Parse first..last
                string temp(optarg);
                size_t found = temp.find("..");
                if(found == string::npos || found == 0 || found + 2 == temp.size()) {
                    cerr << "error: [vg index] could not parse range " << temp << endl;
                    exit(1);
                }
                sample_range.first = std::stoul(temp.substr(0, found));
                sample_range.second = std::stoul(temp.substr(found + 2)) + 1;
            }
            break;
        case 'r':
            {
                // Parse the rename old=new
                string key_value(optarg);
                auto found = key_value.find('=');
                if (found == string::npos || found == 0 || found + 1 == key_value.size()) {
                    cerr << "error: [vg index] could not parse rename " << key_value << endl;
                    exit(1);
                }
                // Parse out the two parts
                string vcf_contig = key_value.substr(0, found);
                string graph_contig = key_value.substr(found + 1);
                // Add the name mapping
                path_to_vcf[graph_contig] = vcf_contig;
            }
            break;
        case 'I':
            {
                // We want to parse this region specifier
                string region(optarg);
                
                Region parsed;
                parse_region(region, parsed);
                if (parsed.start <= 0 || parsed.end <= 0) {
                    // We need both range bounds, and we can't accept 0 since input is 1-based.
                    cerr << "error: [vg index] could not parse 1-based region " << optarg << endl;
                }
                
                // Make sure to correct the coordinates to 0-based exclusive-end, from 1-based inclusive-end
                regions[parsed.seq] = make_pair((size_t) (parsed.start - 1), (size_t) parsed.end);
            }
            break;
        case 'E':
            excluded_samples.insert(optarg);
            break;

        // GCSA
        case 'g':
            build_gcsa = true;
            gcsa_name = optarg;
            break;
        case 'i':
            dbg_names.push_back(optarg);
            break;
        case 'f':
            mapping_name = optarg;
            break;
        case 'k':
            kmer_size = std::max(std::stoul(optarg), 1ul);
            break;
        case 'X':
            params.setSteps(std::stoul(optarg));
            break;
        case 'Z':
            params.setLimit(std::stoul(optarg));
            break;
        case 'V':
            verify_gcsa = true;
            break;

        // RocksDB
        case 'd':
            build_rocksdb = true;
            rocksdb_name = optarg;
            break;
        case 'm':
            store_mappings = true;
            break;
        case 'a':
            store_alignments = true;
            break;
        case 'A':
            dump_alignments = true;
            break;
        case 'N':
            store_node_alignments = true;
            break;
        case 'D':
            dump_index = true;
            break;

        case 'C':
            compact = true;
            break;

        case 'h':
        case '?':
            help_index(argv);
            exit(1);
            break;
        default:
            abort ();
        }
    }

    vector<string> file_names;
    while (optind < argc) {
        string file_name = get_input_file_name(optind, argc, argv);
        file_names.push_back(file_name);
    }

    if (xg_name.empty() && gbwt_name.empty() && threads_name.empty() && gcsa_name.empty() && rocksdb_name.empty()) {
        cerr << "error: [vg index] index type not specified" << endl;
        return 1;
    }

    if ((build_gbwt || write_threads) && !(index_haplotypes || index_paths || index_gam)) {
        cerr << "error: [vg index] cannot build GBWT without threads" << endl;
        return 1;
    }

    if (file_names.size() <= 0 && dbg_names.empty()){
        //cerr << "No graph provided for indexing. Please provide a .vg file or GCSA2-format deBruijn graph to index." << endl;
        //return 1;
    }
    
    if (build_gcsa && kmer_size > gcsa::Key::MAX_LENGTH) {
        cerr << "error: [vg index] GCSA2 cannot index with kmer size greater than " << gcsa::Key::MAX_LENGTH << endl;
        return 1;
    }

    if ((build_gbwt || write_threads) && thread_db_names.size() > 1) {
        cerr << "error: [vg index] cannot use multiple thread database files with -G or -H" << endl;
        return 1;
    }

    // Build XG
    xg::XG* xg_index = new xg::XG();
    map<string, Path> alt_paths;
    if (build_xg) {
        if (file_names.empty()) {
            // VGset or something segfaults when we feed it no graphs.
            cerr << "error: [vg index] at least one graph is required to build an xg index" << endl;
            return 1;
        }
        VGset graphs(file_names);
        build_gpbwt = !build_gbwt & !write_threads;
        graphs.to_xg(*xg_index, index_paths & build_gpbwt, Paths::is_alt, index_haplotypes ? &alt_paths : nullptr);
        if (show_progress) {
            cerr << "Built base XG index" << endl;
        }
    }

    // Generate threads
    if (index_haplotypes || index_paths || index_gam) {

        if (!build_gbwt && !write_threads && !build_gpbwt) {
            cerr << "error: [vg index] No output format specified for the threads" << endl;
            return 1;
        }

        // Use the same temp directory as VG.
        gbwt::TempFile::setDirectory(temp_file::get_dir());

        // if we already made the xg index we can determine the
        size_t id_width;
        if (!index_gam) {
            id_width = gbwt::bit_length(gbwt::Node::encode(xg_index->get_max_id(), true));
        } else { // indexing a GAM
            if (show_progress) {
                cerr << "Finding maximum node id in GAM..." << endl;
            }
            vg::id_t max_id = 0;
            function<void(Alignment&)> lambda = [&](Alignment& aln) {
                gbwt::vector_type buffer;
                for (auto& m : aln.path().mapping()) {
                    max_id = max(m.position().node_id(), max_id);
                }
            };
            for (auto& file_name : gam_file_names) {
                get_input_file(file_name, [&](istream& in) {
                    stream::for_each_parallel(in, lambda);
                });
            }
            id_width = gbwt::bit_length(gbwt::Node::encode(max_id, true));
        }
        
        if (show_progress) {
            cerr << "Node id width: " << id_width << endl;
        }

        vector<string> thread_names;                // Store thread names in insertion order.
        vector<xg::XG::thread_t> all_phase_threads; // Store all threads if building gPBWT.
        size_t haplotype_count = 0;

        // Do we build GBWT?
        gbwt::GBWTBuilder* gbwt_builder = 0;
        if (build_gbwt) {
            if (show_progress) { cerr << "Building GBWT index" << endl; }
            gbwt::Verbosity::set(gbwt::Verbosity::SILENT);  // Make the construction thread silent.
            gbwt_builder = new gbwt::GBWTBuilder(id_width);
        }

        // Do we write threads?
        gbwt::text_buffer_type binary_file;
        if (write_threads) {
            if (show_progress) { cerr << "Writing the threads to " << threads_name << endl; }
            binary_file = gbwt::text_buffer_type(threads_name, std::ios::out, gbwt::MEGABYTE, id_width);
        }

        // Store a thread and its name.
        auto store_thread = [&](const gbwt::vector_type& to_save, const std::string& thread_name) {
            if (build_gbwt) {
                gbwt_builder->insert(to_save, true); // Insert in both orientations.
            }
            if (write_threads) {
                for (auto node : to_save) { binary_file.push_back(node); }
                binary_file.push_back(gbwt::ENDMARKER);
            }
            if (build_gpbwt) {
                xg::XG::thread_t temp;
                temp.reserve(to_save.size());
                for (auto node : to_save) { temp.push_back(gbwt_to_thread_mapping(node)); }
                all_phase_threads.push_back(temp);
            }
            thread_names.push_back(thread_name);
        };

        // Convert paths to threads
        if (index_paths & !build_gpbwt) {
            if (show_progress) {
                cerr << "Converting paths to threads..." << endl;
            }
            for (size_t path_rank = 1; path_rank <= xg_index->max_path_rank(); path_rank++) {
                const xg::XGPath& path = xg_index->get_path(xg_index->path_name(path_rank));
                if (path.ids.size() == 0) {
                    continue;
                }
                gbwt::vector_type buffer(path.ids.size());
                for (size_t i = 0; i < path.ids.size(); i++) {
                    buffer[i] = gbwt::Node::encode(path.node(i), path.is_reverse(i));
                }
                store_thread(buffer, xg_index->path_name(path_rank));
            }
            haplotype_count++; // We assume that the XG index contains the reference paths.
        }

        if (index_gam) {
            if (show_progress) {
                cerr << "Converting GAM to threads..." << endl;
            }
            function<void(Alignment&)> lambda = [&](Alignment& aln) {
                gbwt::vector_type buffer;
                for (auto& m : aln.path().mapping()) {
                    buffer.push_back(gbwt::Node::encode(m.position().node_id(), m.position().is_reverse()));
                }
                store_thread(buffer, aln.name());
                haplotype_count++;
            };
            for (auto& file_name : gam_file_names) {
                get_input_file(file_name, [&](istream& in) {
                    stream::for_each_parallel(in, lambda);
                });
            }
        }

        // Generate haplotypes
        if (index_haplotypes) {
            vcflib::VariantCallFile variant_file;
            variant_file.open(vcf_name);
            if (!variant_file.is_open()) {
                cerr << "error: [vg index] could not open " << vcf_name << endl;
                return 1;
            } else if (show_progress) {
                cerr << "Opened variant file " << vcf_name << endl;
            }
            std::mt19937 rng(0xDEADBEEF);
            std::uniform_int_distribution<std::mt19937::result_type> random_bit(0, 1);

            // How many samples are there?
            size_t num_samples = variant_file.sampleNames.size();
            if (num_samples == 0) {
                cerr << "error: [vg index] The variant file does not contain phasings" << endl;
                return 1;
            }

            // Remember the sample names
            const vector<string>& sample_names = variant_file.sampleNames;

            // Determine the range of samples.
            sample_range.second = std::min(sample_range.second, num_samples);
            haplotype_count += 2 * (sample_range.second - sample_range.first);  // Assuming a diploid genome
            if (show_progress) {
                cerr << "Processing samples " << sample_range.first << " to " << (sample_range.second - 1) << " with batch size " << samples_in_batch << endl;
            }

            // Process each VCF contig corresponding to an XG path.
            size_t max_path_rank = xg_index->max_path_rank();
            for (size_t path_rank = 1; path_rank <= max_path_rank; path_rank++) {
                string path_name = xg_index->path_name(path_rank);
                string vcf_contig_name = path_to_vcf.count(path_name) ? path_to_vcf[path_name] : path_name;
                if (show_progress) {
                    cerr << "Processing path " << path_name << " as VCF contig " << vcf_contig_name << endl;
                }

                // Structures to parse the VCF file into.
                const xg::XGPath& path = xg_index->get_path(path_name);
                gbwt::VariantPaths variants(path.ids.size());
                std::vector<gbwt::PhasingInformation> phasings;

                // Add the reference to VariantPaths.
                for (size_t i = 0; i < path.ids.size(); i++) {
                    variants.appendToReference(gbwt::Node::encode(path.node(i), path.is_reverse(i)));
                }
                variants.indexReference();

                // Create a PhasingInformation for each batch.
                for (size_t batch_start = sample_range.first; batch_start < sample_range.second; batch_start += samples_in_batch) {
                    phasings.emplace_back(batch_start, std::min(samples_in_batch, sample_range.second - batch_start));
                }

                // Set the VCF region or process the entire contig.
                if (regions.count(vcf_contig_name)) {
                    auto region = regions[vcf_contig_name];
                    if (show_progress) {
                        cerr << "- Setting region " << region.first << " to " << region.second << endl;
                    }
                    variant_file.setRegion(vcf_contig_name, region.first, region.second);
                } else {
                    variant_file.setRegion(vcf_contig_name);
                }

                // Parse the variants and the phasings.
                vcflib::Variant var(variant_file);
                size_t variants_processed = 0;
                std::vector<bool> was_diploid(sample_range.second, true); // Was the sample diploid at the previous site?
                while (variant_file.is_open() && variant_file.getNextVariant(var) && var.sequenceName == vcf_contig_name) {
                    // Skip variants with non-DNA sequence, as they are not included in the graph.
                    bool isDNA = allATGC(var.ref);
                    for (vector<string>::iterator a = var.alt.begin(); a != var.alt.end(); ++a) {
                         if (!allATGC(*a)) isDNA = false;
                    }
                    if (!isDNA) {
                        continue;
                    }

                    // Determine the reference nodes for the current variant and create a variant site.
                    // If the variant is not an insertion, there should be a path for the ref allele.
                    var.position--; // Use a 0-based position to get the correct var_name.
                    std::string var_name = make_variant_id(var);
                    std::string ref_path_name = "_alt_" + var_name + "_0";
                    auto ref_path_iter = alt_paths.find(ref_path_name);
                    gbwt::vector_type ref_path;
                    size_t ref_pos = variants.invalid_position();
                    if (ref_path_iter != alt_paths.end() && ref_path_iter->second.mapping_size() != 0) {
                        ref_path = path_to_gbwt(ref_path_iter->second);
                        ref_pos = variants.firstOccurrence(ref_path.front());
                        if (ref_pos == variants.invalid_position()) {
                            cerr << "warning: [vg index] Invalid ref path for " << var_name << " at "
                                 << var.sequenceName << ":" << var.position << endl;
                            continue;
                        }
                    } else { // Try using alt paths instead.
                        bool found = false;
                        for (size_t alt_index = 1; alt_index < var.alleles.size(); alt_index++) {
                            std::string alt_path_name = "_alt_" + var_name + "_" + to_string(alt_index);
                            auto alt_path_iter = alt_paths.find(alt_path_name);
                            if (alt_path_iter != alt_paths.end()) {
                                gbwt::vector_type pred_nodes = predecessors(*xg_index, alt_path_iter->second);
                                for (auto node : pred_nodes) {
                                    size_t pred_pos = variants.firstOccurrence(node);
                                    if (pred_pos != variants.invalid_position()) {
                                        ref_pos = pred_pos + 1;
                                        found = true;
                                        break;
                                    }
                                }
                                if (found) {
                                    break;
                                }
                            }
                        }
                        if (!found) {
                            cerr << "warning: [vg index] Alt and ref paths for " << var_name
                                 << " at " << var.sequenceName << ":" << var.position
                                 << " missing/empty! Was the variant skipped during construction?" << endl;
                            continue;
                        }
                    }
                    variants.addSite(ref_pos, ref_pos + ref_path.size());

                    // Add alternate alleles to the site.
                    for (size_t alt_index = 1; alt_index < var.alleles.size(); alt_index++) {
                        std::string alt_path_name = "_alt_" + var_name + "_" + to_string(alt_index);
                        auto alt_path_iter = alt_paths.find(alt_path_name);
                        if (alt_path_iter != alt_paths.end()) {
                            variants.addAllele(path_to_gbwt(alt_path_iter->second));
                        } else {
                            variants.addAllele(ref_path);
                        }
                    }

                    // Store the phasings in PhasingInformation structures.
                    for (size_t batch = 0; batch < phasings.size(); batch++) {
                        std::vector<gbwt::Phasing> current_phasings;
                        for (size_t sample = phasings[batch].offset(); sample < phasings[batch].limit(); sample++) {
                            string& sample_name = variant_file.sampleNames[sample];
                            current_phasings.emplace_back(var.getGenotype(sample_name), was_diploid[sample]);
                            was_diploid[sample] = current_phasings.back().diploid;
                            if(force_phasing) {
                                current_phasings.back().forcePhased([&]() {
                                   return random_bit(rng);
                                });
                            }
                        }
                        phasings[batch].append(current_phasings);
                    }
                    variants_processed++;
                } // End of variants.
                if (show_progress) {
                    cerr << "- Parsed " << variants_processed << " variants" << endl;
                    size_t phasing_bytes = 0;
                    for (size_t batch = 0; batch < phasings.size(); batch++) {
                        phasing_bytes += phasings[batch].bytes();
                    }
                    cerr << "- Phasing information: " << gbwt::inMegabytes(phasing_bytes) << " MB" << endl;
                }

                // Save memory:
                // - Delete the alt paths if we no longer need them.
                // - Delete the XG index if we no longer need it.
                // - Close the phasings files.
                if (path_rank == max_path_rank) {
                    alt_paths.clear();
                    if (xg_name.empty()) {
                        delete xg_index;
                        xg_index = nullptr;
                    }
                }
                for (size_t batch = 0; batch < phasings.size(); batch++) {
                    phasings[batch].close();
                }

                // Generate the haplotypes.
                for (size_t batch = 0; batch < phasings.size(); batch++) {
                    gbwt::generateHaplotypes(variants, phasings[batch],
                        [&](gbwt::size_type sample) -> bool {
                            return (excluded_samples.find(sample_names[sample]) == excluded_samples.end());
                        },
                        [&](const gbwt::Haplotype& haplotype) {
                            stringstream sn;
                            sn  << "_thread_" << sample_names[haplotype.sample]
                                << "_" << path_name
                                << "_" << haplotype.phase
                                << "_" << haplotype.count;
                            store_thread(haplotype.path, sn.str());
                        });
                    if (show_progress) {
                        cerr << "- Processed samples " << phasings[batch].offset() << " to " << (phasings[batch].offset() + phasings[batch].size() - 1) << endl;
                    }
                }
            } // End of contig.
        } // End of haplotypes.

        // Store the thread database. Write it to disk if a filename is given,
        // or store it in the XG index if building gPBWT or if the XG index
        // will be written to disk.
        alt_paths.clear();
        if (build_gbwt) {
            gbwt_builder->finish();
            if (show_progress) { cerr << "Saving GBWT to disk..." << endl; }
            sdsl::store_to_file(gbwt_builder->index, gbwt_name);
            delete gbwt_builder; gbwt_builder = nullptr;
        }
        if (write_threads) {
            binary_file.close();
        }
        if (build_gbwt || write_threads) {
            if (!thread_db_names.empty()) {
                write_thread_db(thread_db_names.front(), thread_names, haplotype_count);
            } else if (!xg_name.empty()) {
                if (show_progress) {
                    cerr << "Storing " << thread_names.size() << " thread names from "
                         << haplotype_count << " haplotypes in the XG index..." << endl;
                }
                xg_index->set_thread_names(thread_names);
                xg_index->set_haplotype_count(haplotype_count);
            }
        }
        if (build_gpbwt) {
            if (show_progress) {
                cerr << "Inserting all phase threads into DAG..." << endl;
            }
            xg_index->insert_threads_into_dag(all_phase_threads, thread_names);
            xg_index->set_haplotype_count(haplotype_count);            
        }
    }

    // Save XG
    if (!xg_name.empty()) {
        if (!thread_db_names.empty()) {
            vector<string> thread_names;
            size_t haplotype_count = 0;
            read_thread_db(thread_db_names, thread_names, haplotype_count);
            if (show_progress) {
                cerr << thread_names.size() << " threads for "
                     << haplotype_count << " haplotypes in "
                     << thread_db_names.size() << " file(s)" << endl;
            }
            xg_index->set_thread_names(thread_names);
            xg_index->set_haplotype_count(haplotype_count);
        }

        if (show_progress) {
            cerr << "Saving XG index to disk..." << endl;
        }
        ofstream db_out(xg_name);
        xg_index->serialize(db_out);
        db_out.close();
    }
    delete xg_index; xg_index = nullptr;

    // Build GCSA
    if (build_gcsa) {

        // Configure GCSA2 verbosity so it doesn't spit out loads of extra info
        if (!show_progress) {
            gcsa::Verbosity::set(gcsa::Verbosity::SILENT);
        }

        // Use the same temp directory as VG.
        gcsa::TempFile::setDirectory(temp_file::get_dir());

        double start = gcsa::readTimer();

        // Generate temporary kmer files
        bool delete_kmer_files = false;
        if (dbg_names.empty()) {
            if (show_progress) {
                cerr << "Generating kmer files..." << endl;
            }
            VGset graphs(file_names);
            graphs.show_progress = show_progress;
            size_t kmer_bytes = params.getLimitBytes();
            dbg_names = graphs.write_gcsa_kmers_binary(kmer_size, kmer_bytes);
            params.reduceLimit(kmer_bytes);
            delete_kmer_files = true;
        }

        // Build the index
        if (show_progress) {
            cerr << "Building the GCSA2 index..." << endl;
        }
        gcsa::InputGraph input_graph(dbg_names, true, gcsa::Alphabet(), mapping_name);
        gcsa::GCSA gcsa_index(input_graph, params);
        gcsa::LCPArray lcp_array(input_graph, params);
        if (show_progress) {
            double seconds = gcsa::readTimer() - start;
            cerr << "GCSA2 index built in " << seconds << " seconds, "
                 << gcsa::inGigabytes(gcsa::memoryUsage()) << " GB" << endl;
            cerr << "I/O volume: " << gcsa::inGigabytes(gcsa::readVolume()) << " GB read, "
                 << gcsa::inGigabytes(gcsa::writeVolume()) << " GB write" << endl;
        }

        // Save the indexes
        if (show_progress) {
            cerr << "Saving the index to disk..." << endl;
        }
        sdsl::store_to_file(gcsa_index, gcsa_name);
        sdsl::store_to_file(lcp_array, gcsa_name + ".lcp");

        // Verify the index
        if (verify_gcsa) {
            if (show_progress) {
                cerr << "Verifying the index..." << endl;
            }
            if (!gcsa::verifyIndex(gcsa_index, &lcp_array, input_graph)) {
                cerr << "warning: [vg index] GCSA2 index verification failed" << endl;
            }
        }

        // Delete the temporary kmer files
        if (delete_kmer_files) {
            for (auto& filename : dbg_names) {
                temp_file::remove(filename);
            }
        }
    }

    if (build_rocksdb) {

        Index index;

        if (compact) {
            index.open_for_write(rocksdb_name);
            index.compact();
            index.flush();
            index.close();
        }

        if (store_node_alignments && file_names.size() > 0) {
            index.open_for_bulk_load(rocksdb_name);
            int64_t aln_idx = 0;
            function<void(Alignment&)> lambda = [&index,&aln_idx](Alignment& aln) {
                index.cross_alignment(aln_idx++, aln);
            };
            for (auto& file_name : file_names) {
                get_input_file(file_name, [&](istream& in) {
                    stream::for_each_parallel(in, lambda);
                });
            }
            index.flush();
            index.close();
        }

        if (store_alignments && file_names.size() > 0) {
            index.open_for_bulk_load(rocksdb_name);
            function<void(Alignment&)> lambda = [&index](Alignment& aln) {
                index.put_alignment(aln);
            };
            for (auto& file_name : file_names) {
                get_input_file(file_name, [&](istream& in) {
                    stream::for_each_parallel(in, lambda);
                });
            }
            index.flush();
            index.close();
        }

        if (dump_alignments) {
            vector<Alignment> output_buf;
            index.open_read_only(rocksdb_name);
            auto lambda = [&output_buf](const Alignment& aln) {
                output_buf.push_back(aln);
                stream::write_buffered(cout, output_buf, 100);
            };
            index.for_each_alignment(lambda);
            stream::write_buffered(cout, output_buf, 0);
            index.close();
        }

        if (store_mappings && file_names.size() > 0) {
            index.open_for_bulk_load(rocksdb_name);
            function<void(Alignment&)> lambda = [&index](Alignment& aln) {
                const Path& path = aln.path();
                for (int i = 0; i < path.mapping_size(); ++i) {
                    index.put_mapping(path.mapping(i));
                }
            };
            for (auto& file_name : file_names) {
                get_input_file(file_name, [&](istream& in) {
                    stream::for_each_parallel(in, lambda);
                });
            }
            index.flush();
            index.close();
        }

        if (dump_index) {
            index.open_read_only(rocksdb_name);
            index.dump(cout);
            index.close();
        }

    }

    if (show_progress) {
        cerr << "Memory usage: " << gcsa::inGigabytes(gcsa::memoryUsage()) << " GB" << endl;
    }
    return 0;
}

void write_thread_db(const std::string& filename, const std::vector<std::string>& thread_names, size_t haplotype_count) {
    std::ofstream out(filename, std::ios_base::binary);
    if (!out) {
        std::cerr << "error: [vg index] cannot write thread database to " << filename << std::endl;
    }

    out.write(reinterpret_cast<const char*>(&haplotype_count), sizeof(haplotype_count));
    size_t thread_count = thread_names.size();
    out.write(reinterpret_cast<const char*>(&thread_count), sizeof(thread_count));
    for (const std::string& name : thread_names) {
        size_t name_length = name.length();
        out.write(reinterpret_cast<const char*>(&name_length), sizeof(name_length));
        out.write(name.data(), name_length);
    }
    out.close();
}

void read_thread_db(const std::vector<std::string>& filenames, std::vector<std::string>& thread_names, size_t& haplotype_count) {
    thread_names.clear();
    haplotype_count = 0;

    for (const std::string& filename : filenames) {
        std::ifstream in(filename, std::ios_base::binary);
        if (!in) {
            std::cerr << "error: [vg index] cannot read thread database from " << filename << std::endl;
            std::exit(EXIT_FAILURE);
        }

        size_t new_haplotype_count = 0;
        in.read(reinterpret_cast<char*>(&new_haplotype_count), sizeof(new_haplotype_count));
        haplotype_count = std::max(haplotype_count, new_haplotype_count);
        size_t threads_remaining = 0;
        in.read(reinterpret_cast<char*>(&threads_remaining), sizeof(threads_remaining));
        while (threads_remaining > 0) {
            size_t name_length = 0;
            in.read(reinterpret_cast<char*>(&name_length), sizeof(name_length));
            std::vector<char> buffer(name_length);
            in.read(buffer.data(), name_length);
            thread_names.emplace_back(buffer.begin(), buffer.end());
            threads_remaining--;
        }
        in.close();
    }
}

// Register subcommand
static Subcommand vg_construct("index", "index graphs or alignments for random access or mapping", PIPELINE, 2, main_index);
