/**
 * \file mpmap_main.cpp: multipath mapping of reads to a graph
 */

#include <omp.h>
#include <unistd.h>
#include <getopt.h>

#include "subcommand.hpp"

#include <vg/io/vpkg.hpp>
#include "../multipath_mapper.hpp"
#include "../multipath_alignment_emitter.hpp"
#include "../path.hpp"
#include "../watchdog.hpp"
#include "../watchdog.hpp"
#include <bdsg/overlays/overlay_helper.hpp>
#include <bdsg/odgi.hpp>
#include <bdsg/packed_graph.hpp>
#include <bdsg/hash_graph.hpp>
#include <xg.hpp>

//#define record_read_run_times

#ifdef record_read_run_times
#define READ_TIME_FILE "_read_times.tsv"
#include <ctime>
#include <iostream>
#endif

#ifdef mpmap_instrument_mem_statitics
#define MEM_STATS_FILE "_mem_statistics.tsv"
#endif

using namespace std;
using namespace vg;
using namespace vg::subcommand;

void help_mpmap(char** argv) {
    cerr
    << "usage: " << argv[0] << " mpmap [options] -x graph.[xg|pg|og|hg] -g index.gcsa [-f reads1.fq [-f reads2.fq] | -G reads.gam] > aln.gamp" << endl
    << "Multipath align reads to a graph." << endl
    << endl
    << "basic options:" << endl
    << "graph/index:" << endl
    << "  -x, --graph-name FILE         graph (required; XG format recommended but other formats are valid, see `vg convert`) " << endl
    << "  -g, --gcsa-name FILE          use this GCSA2/LCP index pair for MEMs (required; both FILE and FILE.lcp)" << endl
    //<< "  -H, --gbwt-name FILE          use this GBWT haplotype index for population-based MAPQs" << endl
    << "  -d, --dist-name FILE          use this snarl distance index for clustering" << endl
    //<< "      --linear-index FILE       use this sublinear Li and Stephens index file for population-based MAPQs" << endl
    //<< "      --linear-path PATH        use the given path name as the path that the linear index is against" << endl
    << "  -s, --snarls FILE             align to alternate paths in these snarls" << endl
    << "input:" << endl
    << "  -f, --fastq FILE              input FASTQ (possibly gzipped), can be given twice for paired ends (for stdin use -)" << endl
    << "  -i, --interleaved             input contains interleaved paired ends" << endl
    << "algorithm presets:" << endl
    << "  -n, --nt-type TYPE            sequence type preset: 'dna' for genomic data, 'rna' for transcriptomic data [dna]" << endl
    << "  -l, --read-length TYPE        read length preset: 'very-short', 'short', or 'long' (approx. <50bp, 50-500bp, and >500bp) [short]" << endl
    << "  -e, --error-rate TYPE         error rate preset: 'low' or 'high' (approx. PHRED >20 and <20) [low]" << endl
    << "output:" << endl
    << "  -S, --single-path-mode        output single-path alignments (GAM) instead of multipath alignments (GAMP)" << endl
    << "  -N, --sample NAME             add this sample name to output" << endl
    << "  -R, --read-group NAME         add this read group to output" << endl
    << "  -p, --suppress-progress       do not report progress to stderr (slightly reduces thread contention)" << endl
    //<< "algorithm:" << endl
    //<< "       --min-dist-cluster        use the minimum distance based clusterer (requires a distance index from -d)" << endl
//    << "scoring:" << endl
//    << "  -E, --long-read-scoring       set alignment scores to long-read defaults: -q1 -z1 -o1 -y1 -L0 (can be overridden)" << endl
    << "computational parameters:" << endl
    << "  -t, --threads INT             number of compute threads to use [all available]" << endl
    << endl
    << "advanced options:" << endl
    << "algorithm:" << endl
    //<< "  -v, --tvs-clusterer           use the target value search-based clusterer (requires a distance index from -d)" << endl
    //<< "  -X, --snarl-max-cut INT       do not align to extra paths in a snarl if there is an exact match this long (0 for no limit) [5]" << endl
    //<< "  -a, --alt-paths INT           align to (up to) this many alternate paths in snarls [10]" << endl
    //<< "      --suppress-tail-anchors   don't produce extra anchors when aligning to alternate paths in snarls" << endl
    << "  -T, --same-strand             read pairs are from the same strand of the DNA/RNA molecule" << endl
    << "  -M, --max-multimaps INT       report (up to) this many mappings per read [1]" << endl
    << "  -Q, --mq-max INT              cap mapping quality estimates at this much [60]" << endl
    << "  -b, --frag-sample INT         look for this many unambiguous mappings to estimate the fragment length distribution [1000]" << endl
    << "  -I, --frag-mean FLOAT         mean for a pre-determined fragment length distribution (also requires -D)" << endl
    << "  -D, --frag-stddev FLOAT       standard deviation for a pre-determined fragment length distribution (also requires -I)" << endl
    << "  -B, --no-calibrate            do not auto-calibrate mismapping dectection" << endl
    << "  -G, --gam-input FILE          input GAM (for stdin, use -)" << endl
    //<< "  -P, --max-p-val FLOAT         background model p-value must be less than this to avoid mismapping detection [0.00001]" << endl
    << "  -U, --report-group-mapq       add an annotation for the collective mapping quality of all reported alignments" << endl
    //<< "      --padding-mult FLOAT      pad dynamic programming bands in inter-MEM alignment FLOAT * sqrt(read length) [1.0]" << endl
    << "  -u, --map-attempts INT        perform (up to) this many mappings per read (0 for no limit) [24 paired / 64 unpaired]" << endl
    //<< "      --max-paths INT           consider (up to) this many paths per alignment for population consistency scoring, 0 to disable [10]" << endl
    //<< "      --top-tracebacks          consider paths for each alignment based only on alignment score and not based on haplotypes" << endl
    //<< "  -r, --reseed-length INT       reseed SMEMs for internal MEMs if they are at least this long (0 for no reseeding) [28]" << endl
    //<< "  -W, --reseed-diff FLOAT       require internal MEMs to have length within this much of the SMEM's length [0.45]" << endl
    //<< "  -K, --clust-length INT        minimum MEM length used in clustering [automatic]" << endl
    //<< "  -F, --stripped-match          use stripped match algorithm instead of MEMs" << endl
    << "  -c, --hit-max INT             use at most this many hits for any match seeds (0 for no limit) [1024 DNA / 100 RNA]" << endl
    //<< "  --approx-exp FLOAT            let the approximate likelihood miscalculate likelihood ratios by this power [10.0 DNA / 5.0 RNA]" << endl
    //<< "  --recombination-penalty FLOAT use this log recombination penalty for GBWT haplotype scoring [20.7]" << endl
    //<< "  --always-check-population     always try to population-score reads, even if there is only a single mapping" << endl
    //<< "  --delay-population            do not apply population scoring at intermediate stages of the mapping algorithm" << endl
    //<< "  --force-haplotype-count INT   assume that INT haplotypes ought to run through each fixed part of the graph, if nonzero [0]" << endl
    //<< "  -C, --drop-subgraph FLOAT     drop alignment subgraphs whose MEMs cover this fraction less of the read than the best subgraph [0.2]" << endl
    //<< "  --prune-exp FLOAT             prune MEM anchors if their approximate likelihood is this root less than the optimal anchors [1.25]" << endl
    << "scoring:" << endl
    << "  -A, --no-qual-adjust          do not perform base quality adjusted alignments even when base qualities are available" << endl
    << "  -q, --match INT               use this match score [1]" << endl
    << "  -z, --mismatch INT            use this mismatch penalty [4 low error, 1 high error]" << endl
    << "  -o, --gap-open INT            use this gap open penalty [6 low error, 1 high error]" << endl
    << "  -y, --gap-extend INT          use this gap extension penalty [1]" << endl
    << "  -L, --full-l-bonus INT        add this score to alignments that align each end of the read [mismatch+1 short, 0 long]" << endl
    << "  -w, --score-matrix FILE       read a 4x4 integer substitution scoring matrix from a file (in the order ACGT)" << endl
    << "  -m, --remove-bonuses          remove full length alignment bonuses in reported scores" << endl;
    //<< "computational parameters:" << endl
    //<< "  -Z, --buffer-size INT         buffer this many alignments together (per compute thread) before outputting to stdout [200]" << endl;    
}



int main_mpmap(int argc, char** argv) {

    if (argc == 2) {
        help_mpmap(argv);
        return 1;
    }

    // initialize parameters with their default options
    #define OPT_PRUNE_EXP 1000
    #define OPT_RECOMBINATION_PENALTY 1001
    #define OPT_ALWAYS_CHECK_POPULATION 1002
    #define OPT_FORCE_HAPLOTYPE_COUNT 1004
    #define OPT_SUPPRESS_TAIL_ANCHORS 1005
    #define OPT_TOP_TRACEBACKS 1006
    #define OPT_MIN_DIST_CLUSTER 1007
    #define OPT_APPROX_EXP 1008
    #define OPT_MAX_PATHS 1009
    #define OPT_GREEDY_MIN_DIST 1010
    #define OPT_COMPONENT_MIN_DIST 1011
    #define OPT_BAND_PADDING_MULTIPLIER 1012
    #define OPT_HARD_HIT_MAX_MULTIPLIER 1013
    #define OPT_MAX_RESCUE_ATTEMPTS 1014
    #define OPT_STRIP_LENGTH 1015
    #define OPT_STRIP_COUNT 1016
    #define OPT_SECONDARY_RESCUE_ATTEMPTS 1017
    #define OPT_SECONDARY_MAX_DIFF 1018
    #define OPT_NO_CLUSTER 1019
    #define OPT_NO_GREEDY_MEM_RESTARTS 1020
    #define OPT_GREEDY_MEM_RESTART_MAX_LCP 1021
    #define OPT_NO_OUTPUT 1022
    string matrix_file_name;
    string graph_name;
    string gcsa_name;
    string gbwt_name;
    string sublinearLS_name;
    string sublinearLS_ref_path;
    string snarls_name;
    string distance_index_name;
    string fastq_name_1;
    string fastq_name_2;
    string gam_file_name;
    int match_score = default_match;
    int mismatch_score = default_mismatch;
    int gap_open_score = default_gap_open;
    int gap_extension_score = default_gap_extension;
    int full_length_bonus = default_full_length_bonus;
    bool interleaved_input = false;
    int default_snarl_cut_size = 5;
    int snarl_cut_size = default_snarl_cut_size;
    int max_branch_trim_length = 5;
    bool synthesize_tail_anchors = false;
    int max_paired_end_map_attempts = 24;
    int max_single_end_map_attempts = 64;
    int max_single_end_mappings_for_rescue = max_single_end_map_attempts;
    int max_rescue_attempts = 10;
    int population_max_paths = 10;
    int population_paths_hard_cap = 1000;
    bool top_tracebacks = false;
    // How many distinct single path alignments should we look for in a multipath, for MAPQ?
    // TODO: create an option.
    int localization_max_paths = 5;
    int max_num_mappings = 1;
    int buffer_size = 200;
    int hit_max = 1024;
    int hit_max_arg = numeric_limits<int>::min();
    int hard_hit_max_muliplier = 3;
    int min_mem_length = 1;
    int min_clustering_mem_length = 0;
    int min_clustering_mem_length_arg = numeric_limits<int>::min();
    bool use_stripped_match_alg = false;
    int default_strip_length = 10;
    int stripped_match_alg_strip_length = default_strip_length;
    int stripped_match_alg_max_length = 0; // no maximum yet
    int default_strip_count = 10;
    int stripped_match_alg_target_count = default_strip_count;
    bool use_greedy_mem_restarts = true;
    // TODO: it would be best if these parameters responded to the size of the graph...
    int greedy_restart_min_length = 30;
    int greedy_restart_max_lcp = 25;
    int greedy_restart_max_count = 2;
    bool greedy_restart_assume_substitution = true;
    int reseed_length = 28;
    int reseed_length_arg = numeric_limits<int>::min();
    double reseed_diff = 0.45;
    double reseed_diff_arg = numeric_limits<double>::lowest();
    double reseed_exp = 0.065;
    bool use_adaptive_reseed = true;
    double cluster_ratio = 0.2;
    bool use_tvs_clusterer = false;
    bool use_min_dist_clusterer = false;
    bool greedy_min_dist = false;
    bool component_min_dist = true;
    bool no_clustering = false;
    bool qual_adjusted = true;
    bool strip_full_length_bonus = false;
    MappingQualityMethod mapq_method = Exact;
    bool report_group_mapq = false;
    double band_padding_multiplier = 1.0;
    int max_dist_error = 12;
    int default_num_alt_alns = 10;
    int num_alt_alns = default_num_alt_alns;
    double suboptimal_path_exponent = 1.25;
    double likelihood_approx_exp = 10.0;
    double likelihood_approx_exp_arg = numeric_limits<double>::lowest();
    double recombination_penalty = 20.7;
    bool always_check_population = false;
    size_t force_haplotype_count = 0;
    bool single_path_alignment_mode = false;
    int max_mapq = 60;
    size_t frag_length_sample_size = 1000;
    double frag_length_robustness_fraction = 0.95;
    double frag_length_mean = NAN;
    double frag_length_stddev = NAN;
    bool same_strand = false;
    bool auto_calibrate_mismapping_detection = true;
    double max_mapping_p_value = 0.00001;
    size_t num_calibration_simulations = 100;
    vector<size_t> calibration_read_lengths{50, 100, 150, 250, 450};
    bool use_weibull_calibration = false;
    size_t order_length_repeat_hit_max = 3000;
    size_t sub_mem_count_thinning = 4;
    size_t sub_mem_thinning_burn_in = 16;
    double secondary_rescue_score_diff = 0.8;
    size_t secondary_rescue_attempts = 4;
    size_t secondary_rescue_attempts_arg = numeric_limits<size_t>::max();
    size_t rescue_only_min = numeric_limits<size_t>::max(); // disabling this for now
    size_t rescue_only_anchor_max = 16;
    string sample_name = "";
    string read_group = "";
    bool prefilter_redundant_hits = true;
    bool precollapse_order_length_hits = true;
    int max_sub_mem_recursion_depth = 1;
    int max_map_attempts_arg = 0;
    int secondary_rescue_subopt_diff = 35;
    int min_median_mem_coverage_for_split = 0;
    bool suppress_cluster_merging = false;
    bool dynamic_max_alt_alns = true;
    bool simplify_topologies = true;
    int max_alignment_gap = 5000;
    double pessimistic_tail_gap_multiplier = 0.0; // i.e. none
    int match_score_arg = std::numeric_limits<int>::min();
    int mismatch_score_arg = std::numeric_limits<int>::min();
    int gap_open_score_arg = std::numeric_limits<int>::min();
    int gap_extension_score_arg = std::numeric_limits<int>::min();
    int full_length_bonus_arg = std::numeric_limits<int>::min();
    int reversing_walk_length = 1;
    bool no_output = false;

    // default presets
    string nt_type = "dna";
    string read_length = "short";
    string error_rate = "low";
    
    // logging and warning
    bool suppress_progress = false;
    int fragment_length_warning_factor = 50;
    
    int c;
    optind = 2; // force optind past command positional argument
    while (true) {
        static struct option long_options[] =
        {
            {"help", no_argument, 0, 'h'},
            {"graph-name", required_argument, 0, 'x'},
            {"gcsa-name", required_argument, 0, 'g'},
            {"gbwt-name", required_argument, 0, 'H'},
            {"dist-name", required_argument, 0, 'd'},
            {"linear-index", required_argument, 0, 1},
            {"linear-path", required_argument, 0, 2},
            {"fastq", required_argument, 0, 'f'},
            {"gam-input", required_argument, 0, 'G'},
            {"sample", required_argument, 0, 'N'},
            {"read-group", required_argument, 0, 'R'},
            {"interleaved", no_argument, 0, 'i'},
            {"same-strand", no_argument, 0, 'T'},
            {"single-path-mode", no_argument, 0, 'S'},
            {"snarls", required_argument, 0, 's'},
            {"synth-tail-anchors", no_argument, 0, OPT_SUPPRESS_TAIL_ANCHORS},
            {"tvs-clusterer", no_argument, 0, 'v'},
            {"snarl-max-cut", required_argument, 0, 'X'},
            {"alt-paths", required_argument, 0, 'a'},
            {"frag-sample", required_argument, 0, 'b'},
            {"frag-mean", required_argument, 0, 'I'},
            {"frag-stddev", required_argument, 0, 'D'},
            {"max-rescues", required_argument, 0, OPT_MAX_RESCUE_ATTEMPTS},
            {"max-secondary-rescues", required_argument, 0, OPT_SECONDARY_RESCUE_ATTEMPTS},
            {"secondary-diff", required_argument, 0, OPT_SECONDARY_MAX_DIFF},
            {"no-calibrate", no_argument, 0, 'B'},
            {"max-p-val", required_argument, 0, 'P'},
            {"mq-max", required_argument, 0, 'Q'},
            {"report-group-mapq", no_argument, 0, 'U'},
            {"padding-mult", required_argument, 0, OPT_BAND_PADDING_MULTIPLIER},
            {"map-attempts", required_argument, 0, 'u'},
            {"max-paths", required_argument, 0, OPT_MAX_PATHS},
            {"top-tracebacks", no_argument, 0, OPT_TOP_TRACEBACKS},
            {"max-multimaps", required_argument, 0, 'M'},
            {"reseed-length", required_argument, 0, 'r'},
            {"reseed-diff", required_argument, 0, 'W'},
            {"clustlength", required_argument, 0, 'K'},
            {"stripped-match", no_argument, 0, 'F'},
            {"strip-length", no_argument, 0, OPT_STRIP_LENGTH},
            {"strip-count", no_argument, 0, OPT_STRIP_COUNT},
            {"no-greedy-restart", no_argument, 0, OPT_NO_GREEDY_MEM_RESTARTS},
            {"greedy-max-lcp", required_argument, 0, OPT_GREEDY_MEM_RESTART_MAX_LCP},
            {"hit-max", required_argument, 0, 'c'},
            {"hard-hit-mult", required_argument, 0, OPT_HARD_HIT_MAX_MULTIPLIER},
            {"approx-exp", required_argument, 0, OPT_APPROX_EXP},
            {"recombination-penalty", required_argument, 0, OPT_RECOMBINATION_PENALTY},
            {"always-check-population", no_argument, 0, OPT_ALWAYS_CHECK_POPULATION},
            {"force-haplotype-count", required_argument, 0, OPT_FORCE_HAPLOTYPE_COUNT},
            {"min-dist-cluster", no_argument, 0, OPT_MIN_DIST_CLUSTER},
            {"greedy-min-dist", no_argument, 0, OPT_GREEDY_MIN_DIST},
            {"component-min-dist", no_argument, 0, OPT_COMPONENT_MIN_DIST},
            {"no-cluster", no_argument, 0, OPT_NO_CLUSTER},
            {"drop-subgraph", required_argument, 0, 'C'},
            {"prune-exp", required_argument, 0, OPT_PRUNE_EXP},
            {"long-read-scoring", no_argument, 0, 'E'},
            {"read-length", required_argument, 0, 'l'},
            {"nt-type", required_argument, 0, 'n'},
            {"error-rate", required_argument, 0, 'e'},
            {"match", required_argument, 0, 'q'},
            {"mismatch", required_argument, 0, 'z'},
            {"score-matrix", required_argument, 0, 'w'},
            {"gap-open", required_argument, 0, 'o'},
            {"gap-extend", required_argument, 0, 'y'},
            {"full-l-bonus", required_argument, 0, 'L'},
            {"remove-bonuses", no_argument, 0, 'm'},
            {"no-qual-adjust", no_argument, 0, 'A'},
            {"threads", required_argument, 0, 't'},
            {"buffer-size", required_argument, 0, 'Z'},
            {"no-output", no_argument, 0, OPT_NO_OUTPUT},
            {0, 0, 0, 0}
        };

        int option_index = 0;
        c = getopt_long (argc, argv, "hx:g:H:d:f:G:N:R:iSs:vX:u:a:b:I:D:BP:Q:UpM:r:W:K:Fc:C:R:En:l:e:q:z:w:o:y:L:mAt:Z:",
                         long_options, &option_index);


        // Detect the end of the options.
        if (c == -1)
            break;

        switch (c)
        {
            case 'x':
                graph_name = optarg;
                if (graph_name.empty()) {
                    cerr << "error:[vg mpmap] Must provide Graph file with -x." << endl;
                    exit(1);
                }
                break;
                
            case 'g':
                gcsa_name = optarg;
                if (gcsa_name.empty()) {
                    cerr << "error:[vg mpmap] Must provide GCSA file with -g." << endl;
                    exit(1);
                }
                break;
                
            case 'H':
                gbwt_name = optarg;
                if (gbwt_name.empty()) {
                    cerr << "error:[vg mpmap] Must provide GBWT index file with -H" << endl;
                    exit(1);
                }
                break;
                
            case 'd':
                distance_index_name = optarg;
                if (distance_index_name.empty()) {
                    cerr << "error:[vg mpmap] Must provide distance index file with -d" << endl;
                    exit(1);
                }
                if (!use_tvs_clusterer) {
                    use_min_dist_clusterer = true;
                }
                break;
                
            case OPT_MAX_RESCUE_ATTEMPTS:
                max_rescue_attempts = parse<int>(optarg);
                break;
                
            case OPT_SECONDARY_RESCUE_ATTEMPTS:
                secondary_rescue_attempts_arg = parse<int>(optarg);
                break;
                
            case OPT_SECONDARY_MAX_DIFF:
                secondary_rescue_score_diff = parse<double>(optarg);
                break;
                
            case 1: // --linear-index
                sublinearLS_name = optarg;
                break;
            
            case 2: // --linear-path
                sublinearLS_ref_path = optarg;
                break;
                
            case 'f':
                if (fastq_name_1.empty()) {
                    fastq_name_1 = optarg;
                    if (fastq_name_1.empty()) {
                        cerr << "error:[vg mpmap] Must provide FASTQ file with -f" << endl;
                        exit(1);
                    }
                }
                else if (fastq_name_2.empty()) {
                    fastq_name_2 = optarg;
                    if (fastq_name_2.empty()) {
                        cerr << "error:[vg mpmap] Must provide FASTQ file with -f" << endl;
                        exit(1);
                    }
                }
                else {
                    cerr << "error:[vg mpmap] Cannot specify more than two FASTQ files" << endl;
                    exit(1);
                }
                break;
                
            case 'G':
                gam_file_name = optarg;
                if (gam_file_name.empty()) {
                    cerr << "error:[vg mpmap] Must provide GAM file with -G." << endl;
                    exit(1);
                }
                break;
                
            case 'N':
                sample_name = optarg;
                if (sample_name.empty()) {
                    cerr << "error:[vg mpmap] Must provide sample name with -N." << endl;
                    exit(1);
                }
                break;
                
            case 'R':
                read_group = optarg;
                if (read_group.empty()) {
                    cerr << "error:[vg mpmap] Must provide read group with -R." << endl;
                    exit(1);
                }
                break;
                
            case 'i':
                interleaved_input = true;
                break;
                
            case 'T':
                same_strand = true;
                break;
                
            case 'p':
                suppress_progress = true;
                break;
                
            case 'S':
                single_path_alignment_mode = true;
                break;
                
            case 's':
                snarls_name = optarg;
                if (snarls_name.empty()) {
                    cerr << "error:[vg mpmap] Must provide snarl file with -s." << endl;
                    exit(1);
                }
                break;
                
            case OPT_SUPPRESS_TAIL_ANCHORS:
                synthesize_tail_anchors = true;
                break;
                
            case 'v':
                use_tvs_clusterer = true;
                use_min_dist_clusterer = false;
                break;
                
            case 'X':
                snarl_cut_size = parse<int>(optarg);
                break;
                
            case 'a':
                num_alt_alns = parse<int>(optarg);
                break;
                
            case 'b':
                frag_length_sample_size = parse<int>(optarg);
                break;
                
            case 'I':
                frag_length_mean = parse<double>(optarg);
                break;
                
            case 'D':
                frag_length_stddev = parse<double>(optarg);
                break;
                
            case 'B':
                auto_calibrate_mismapping_detection = false;
                break;
                
            case 'P':
                max_mapping_p_value = parse<double>(optarg);
                break;
                
            case 'Q':
                max_mapq = parse<int>(optarg);
                break;
                
            case 'U':
                report_group_mapq = true;
                break;
                
            case OPT_BAND_PADDING_MULTIPLIER:
                band_padding_multiplier = parse<double>(optarg);
                break;
                
            case 'u':
                max_map_attempts_arg = parse<int>(optarg);
                // let 0 be a sentinel for no limit and also a sentinel for not giving an arg
                if (max_map_attempts_arg == 0) {
                    max_map_attempts_arg = numeric_limits<int>::max();
                }
                break;
                
            case OPT_MAX_PATHS:
                population_max_paths = parse<int>(optarg);
                break;
                
            case OPT_TOP_TRACEBACKS:
                top_tracebacks = true;
                break;
                
            case 'M':
                max_num_mappings = parse<int>(optarg);
                break;
                
            case 'r':
                reseed_length_arg = parse<int>(optarg);
                break;
                
            case 'W':
                reseed_diff_arg = parse<double>(optarg);
                break;
                
            case 'K':
                min_clustering_mem_length_arg = parse<int>(optarg);
                break;
                
            case 'F':
                use_stripped_match_alg = true;
                break;
                
            case OPT_STRIP_LENGTH:
                stripped_match_alg_strip_length = parse<int>(optarg);
                break;
                
            case OPT_STRIP_COUNT:
                stripped_match_alg_target_count = parse<int>(optarg);
                break;
                
            case OPT_NO_GREEDY_MEM_RESTARTS:
                use_greedy_mem_restarts = false;
                break;
                
            case OPT_GREEDY_MEM_RESTART_MAX_LCP:
                greedy_restart_max_lcp = parse<int>(optarg);
                break;
                
            case 'c':
                hit_max_arg = parse<int>(optarg);
                break;
                
            case OPT_HARD_HIT_MAX_MULTIPLIER:
                hard_hit_max_muliplier = parse<int>(optarg);
                break;
                
            case OPT_APPROX_EXP:
                likelihood_approx_exp_arg = parse<double>(optarg);
                break;
                
            case OPT_RECOMBINATION_PENALTY:
                recombination_penalty = parse<double>(optarg);
                break;
                
            case OPT_ALWAYS_CHECK_POPULATION:
                always_check_population = true;
                break;
                
            case OPT_FORCE_HAPLOTYPE_COUNT:
                force_haplotype_count = parse<size_t>(optarg);
                break;
                
            case OPT_MIN_DIST_CLUSTER:
                // This the default behavior
                //use_min_dist_clusterer = true;
                break;
                
            case OPT_GREEDY_MIN_DIST:
                greedy_min_dist = true;
                component_min_dist = false;
                break;
                
            case OPT_COMPONENT_MIN_DIST:
                // the default now
                component_min_dist = true;
                break;
                
            case OPT_NO_CLUSTER:
                no_clustering = true;
                break;
                
            case 'C':
                cluster_ratio = parse<double>(optarg);
                break;
                
            case OPT_PRUNE_EXP:
                suboptimal_path_exponent = parse<double>(optarg);
                break;
                
            case 'E':
                cerr << "warning:[vg mpmap] Long read scoring option (--long-read-scoring) is deprecated. Instead, use read length preset (--read-length)." << endl;
                read_length = "long";
                break;
                
            case 'l':
                read_length = optarg;
                break;
                
            case 'e':
                error_rate = optarg;
                break;
                
            case 'n':
                nt_type = optarg;
                break;
                
            case 'q':
                match_score_arg = parse<int>(optarg);
                break;
                
            case 'z':
                mismatch_score_arg = parse<int>(optarg);
                break;
                
            case 'w':
                matrix_file_name = optarg;
                if (matrix_file_name.empty()) {
                    cerr << "error:[vg mpmap] Must provide matrix file with --matrix-file." << endl;
                    exit(1);
                }
                break;
                
            case 'o':
                gap_open_score_arg = parse<int>(optarg);
                break;
                
            case 'y':
                gap_extension_score_arg = parse<int>(optarg);
                break;
                
            case 'm':
                strip_full_length_bonus = true;
                break;
                
            case 'L':
                full_length_bonus_arg = parse<int>(optarg);
                break;
                
            case 'A':
                qual_adjusted = false;
                break;
                
            case 't':
            {
                int num_threads = parse<int>(optarg);
                if (num_threads <= 0) {
                    cerr << "error:[vg mpmap] Thread count (-t) set to " << num_threads << ", must set to a positive integer." << endl;
                    exit(1);
                }
                omp_set_num_threads(num_threads);
            }
                break;
                
            case 'Z':
                buffer_size = parse<int>(optarg);
                break;
                
            case OPT_NO_OUTPUT:
                no_output = true;
                break;
                
            case 'h':
            case '?':
            default:
                help_mpmap(argv);
                exit(1);
                break;
        }
    }
    
    // normalize capitalization on preset options
    if (read_length == "Long" || read_length == "LONG") {
        read_length = "long";
    }
    else if (read_length == "Very-Short" || read_length == "Very-short" || read_length == "VERY-SHORT") {
        read_length = "very-short";
    }
    else if (read_length != "Short" || read_length != "SHORT") {
        read_length = "short";
    }
    
    if (nt_type == "RNA") {
        nt_type = "rna";
    }
    else if (nt_type == "DNA") {
        nt_type = "dna";
    }
    
    if (error_rate == "Low" || error_rate == "LOW") {
        error_rate = "low";
    }
    else if (error_rate == "High" || error_rate == "HIGH") {
        error_rate = "high";
    }
        
    // set baseline parameters according to presets
    
    if (error_rate == "high") {
        // alignment scores that don't penalize gaps or mismatches as much
        mismatch_score = 1;
        gap_open_score = 1;
        // do less DP on tails (having a presumption that long tails with no seeds
        // will probably be soft-clipped)
        pessimistic_tail_gap_multiplier = 3.0;
        // we won't assume that errors are likely to be substitutions
        greedy_restart_assume_substitution = false;
    }
    
    if (read_length == "long") {
        // we don't care so much about soft-clips on long reads
        full_length_bonus = 0;
        // long read technologies (even low error ones) don't have the same bias toward
        // substitution variants as NGS
        greedy_restart_assume_substitution = false;
    }
    else if (read_length == "very-short") {
        // clustering is unlikely to improve accuracy in very short data
        // TODO: should i add a way to override this option?
        no_clustering = true;
        // we don't want to throw away short matches a priori in very short data
        min_clustering_mem_length = 1;
    }
    else if (read_length != "short") {
        // short is the default
        cerr << "error:[vg mpmap] Cannot identify read length preset (-l): " << read_length << endl;
        exit(1);
    }
    
    if (nt_type == "rna") {
        // RNA preset
        if (distance_index_name.empty()) {
            cerr << "warning:[vg mpmap] It is HIGHLY recommended to use a distance index (-d) for clustering on splice graphs. Both accuracy and speed will suffer without one." << endl;
        }
        // seed finding, cluster pruning, and rescue parameters tuned for a lower repeat content
        secondary_rescue_attempts = 1;
        max_single_end_mappings_for_rescue = 32;
        hit_max = 100;
        reseed_length = 28; // TODO: returned to the DNA value
        reseed_diff = 0.6;
        likelihood_approx_exp = 3.5;
    }
    else if (nt_type != "dna") {
        // DNA is the default
        cerr << "error:[vg mpmap] Cannot identify sequencing type preset (-n): " << nt_type << endl;
        exit(1);
    }
        
    if (single_path_alignment_mode && read_length != "long") {
        // we get better performance by splitting up clusters a bit more when we're forcing alignments to go to only one place
        // long reads on the other hand display underclustering behavior with some parameter settings
        min_median_mem_coverage_for_split = 2;
        suppress_cluster_merging = true;
    }
    
    if (read_length == "long" && !single_path_alignment_mode) {
        // we sometimes need to synthesize anchors for the tails to get good multipath alignments on long reads
        synthesize_tail_anchors = true;
    }
        
    if (single_path_alignment_mode) {
        // simplifying topologies is redundant work if we're just going to take the maximum weight path anyway
        simplify_topologies = false;
    }
    
    if (single_path_alignment_mode &&
        (population_max_paths == 0 || (sublinearLS_name.empty() && gbwt_name.empty()))) {
        // adjust parameters that produce irrelevant extra work single path mode
        if (!snarls_name.empty()) {
            cerr << "warning:[vg mpmap] Snarl file (-s) is ignored in single path mode (-S) without multipath population scoring (--max-paths)." << endl;
        }
        
        if (snarl_cut_size != default_snarl_cut_size) {
            cerr << "warning:[vg mpmap] Snarl cut limit (-X) is ignored in single path mode (-S) without multipath population scoring (--max-paths)." << endl;
        }
        
        if (num_alt_alns != default_num_alt_alns) {
            cerr << "warning:[vg mpmap] Number of alternate alignments (-a) is ignored in single path mode (-S) without multipath population scoring (--max-paths)." << endl;
        }
        
        // don't cut inside snarls or load the snarl manager
        snarl_cut_size = 0;
        snarls_name = "";
        
        // only get 1 traceback for an inter-MEM or tail alignment
        dynamic_max_alt_alns = false;
        num_alt_alns = 1;
    }
        
    // set the overrides to preset-controlled parameters
    if (hit_max_arg != numeric_limits<int>::min()) {
        hit_max = hit_max_arg;
    }
    if (reseed_length_arg != numeric_limits<int>::min()) {
        reseed_length = reseed_length_arg;
    }
    if (reseed_diff_arg != numeric_limits<double>::lowest()) {
        reseed_diff = reseed_diff_arg;
    }
    if (min_clustering_mem_length_arg != numeric_limits<int>::min()) {
        min_clustering_mem_length = min_clustering_mem_length_arg;
    }
    if (secondary_rescue_attempts_arg != numeric_limits<size_t>::max()) {
        secondary_rescue_attempts = secondary_rescue_attempts_arg;
    }
    if (likelihood_approx_exp_arg != numeric_limits<double>::lowest()) {
        likelihood_approx_exp = likelihood_approx_exp_arg;
    }
    if (match_score_arg != std::numeric_limits<int>::min()) {
        match_score = match_score_arg;
    }
    if (mismatch_score_arg != std::numeric_limits<int>::min()) {
        mismatch_score = mismatch_score_arg;
    }
    if (gap_open_score_arg != std::numeric_limits<int>::min()) {
        gap_open_score = gap_open_score_arg;
    }
    if (gap_extension_score_arg != std::numeric_limits<int>::min()) {
        gap_extension_score = gap_extension_score_arg;
    }
    if (full_length_bonus_arg != std::numeric_limits<int>::min()) {
        full_length_bonus = full_length_bonus_arg;
    }
    else if (read_length != "long") {
        // TODO: not so elegant
        // the full length bonus should override a mismatch unless we're in long read mode
        full_length_bonus = min<int>(mismatch_score + 1, std::numeric_limits<int8_t>::max());
    }
        
    // choose either the user supplied max or the default for paired/unpaired
    int max_map_attempts = max_map_attempts_arg ? max_map_attempts_arg : ((interleaved_input || !fastq_name_2.empty()) ?
                                                                          max_paired_end_map_attempts : max_single_end_map_attempts);
    if (max_map_attempts_arg) {
        max_single_end_mappings_for_rescue = max_map_attempts_arg;
    }
    
    // hits that are much more frequent than the number of hits we sample are unlikely to produce high MAPQs, so
    // we can usually ignore them
    int hard_hit_max = hard_hit_max_muliplier * hit_max;
    
    // check for valid parameters
    
    if (std::isnan(frag_length_mean) != std::isnan(frag_length_stddev)) {
        cerr << "error:[vg mpmap] Cannot specify only one of fragment length mean (-I) and standard deviation (-D)." << endl;
        exit(1);
    }
    
    if (!std::isnan(frag_length_mean) && frag_length_mean < 0) {
        cerr << "error:[vg mpmap] Fragment length mean (-I) must be nonnegative." << endl;
        exit(1);
    }
    
    if (!std::isnan(frag_length_stddev) && frag_length_stddev < 0) {
        cerr << "error:[vg mpmap] Fragment length standard deviation (-D) must be nonnegative." << endl;
        exit(1);
    }
    
    if (interleaved_input && !fastq_name_2.empty()) {
        cerr << "error:[vg mpmap] Cannot designate both interleaved paired ends (-i) and separate paired end file (-f)." << endl;
        exit(1);
    }
    
    
    if (!fastq_name_1.empty() && !gam_file_name.empty()) {
        cerr << "error:[vg mpmap] Cannot designate both FASTQ input (-f) and GAM input (-G) in same run." << endl;
        exit(1);
    }
    
    if (fastq_name_1.empty() && gam_file_name.empty()) {
        cerr << "error:[vg mpmap] Must designate reads to map from either FASTQ (-f) or GAM (-G) file." << endl;
        exit(1);
    }
    
    if (!interleaved_input && fastq_name_2.empty() && same_strand) {
        cerr << "warning:[vg mpmap] Ignoring same strand parameter (-e) because no paired end input provided." << endl;
    }
    
    if (num_alt_alns <= 0) {
        cerr << "error:[vg mpmap] Number of alternate snarl paths (-a) set to " << num_alt_alns << ", must set to a positive integer." << endl;
        exit(1);
    }
    
    if (frag_length_sample_size <= 0) {
        cerr << "error:[vg mpmap] Fragment length distribution sample size (-b) set to " << frag_length_sample_size << ", must set to a positive integer." << endl;
        exit(1);
    }
    
    if (snarl_cut_size < 0) {
        cerr << "error:[vg mpmap] Max snarl cut size (-X) set to " << snarl_cut_size << ", must set to a positive integer or 0 for no maximum." << endl;
        exit(1);
    }
    
    if (max_mapping_p_value <= 0.0) {
        cerr << "error:[vg mpmap] Max mapping p-value (-P) set to " << max_mapping_p_value << ", must set to a positive number." << endl;
        exit(1);
    }
    
    if (max_mapq <= 0 && mapq_method != None) {
        cerr << "error:[vg mpmap] Maximum mapping quality (-Q) set to " << max_mapq << ", must set to a positive integer." << endl;
        exit(1);
    }
    
    if (band_padding_multiplier < 0.0) {
        cerr << "error:[vg mpmap] Band padding (-p) set to " << band_padding_multiplier << ", must set to a nonnegative number." << endl;
        exit(1);
    }
    
    if (max_map_attempts_arg < 0) {
        cerr << "error:[vg mpmap] Maximum number of mapping attempts (-u) set to " << max_map_attempts_arg << ", must set to a positive integer or 0 for no maximum." << endl;
        exit(1);
    }
    
    if (population_max_paths < 0) {
        cerr << "error:[vg mpmap] Maximum number of paths per alignment for population scoring (--max-paths) set to " << population_max_paths << ", must set to a nonnegative integer." << endl;
        exit(1);
    }
    
    if (population_max_paths != 10 && population_max_paths != 0 && gbwt_name.empty() && sublinearLS_name.empty()) {
        // Don't allow anything but the default or the "disabled" setting without an index.
        // TODO: This restriction makes neat auto-generation of command line options for different conditions hard.
        cerr << "error:[vg mpmap] Maximum number of paths per alignment for population scoring (--max-paths) is specified but population database (-H or --linear-index) was not provided." << endl;
        exit(1);
    }
    
    if (always_check_population && gbwt_name.empty() && sublinearLS_name.empty()) {
        cerr << "error:[vg mpmap] Cannot --always-check-population if no population database (-H or --linear-index) is provided." << endl;
        exit(1);
    }
    
    if (force_haplotype_count != 0 && gbwt_name.empty() && sublinearLS_name.empty()) {
        cerr << "warning:[vg mpmap] Cannot --force-haplotype-count if no population database (-H or --linear-index) is provided. Ignoring option." << endl;
    }
    
    if (!sublinearLS_name.empty() && !gbwt_name.empty()) {
        cerr << "error:[vg mpmap] GBWT index (-H) and linear haplotype index (--linear-index) both specified. Only one can be used." << endl;
        exit(1);
    }
    
    if (!sublinearLS_name.empty() && sublinearLS_ref_path.empty()) {
        cerr << "error:[vg mpmap] Linear haplotype index (--linear-index) cannot be used without a single reference path (--linear-path)." << endl;
        exit(1);
    }
    
    if (sublinearLS_name.empty() && !sublinearLS_ref_path.empty()) {
        cerr << "error:[vg mpmap] Linear haplotype ref path (--linear-path) cannot be used without an index (--linear-index)." << endl;
        exit(1);
    }
    
    if (max_num_mappings > max_map_attempts && max_map_attempts != 0) {
        cerr << "warning:[vg mpmap] Reporting up to " << max_num_mappings << " mappings, but only computing up to " << max_map_attempts << " mappings." << endl;
    }
    
    if (max_rescue_attempts < 0) {
        cerr << "error:[vg mpmap] Maximum number of rescue attempts (--max-rescues) set to " << max_rescue_attempts << ", must set to a non-negative integer (0 for no rescue)." << endl;
        exit(1);
    }
    
    if (max_rescue_attempts > max_single_end_mappings_for_rescue) {
        cerr << "warning:[vg mpmap] Maximum number of rescue attempts (--max-rescues) of " << max_rescue_attempts << " is greater than number of mapping attempts for rescue " << max_single_end_mappings_for_rescue << endl;
    }
    
    if (secondary_rescue_attempts < 0) {
        cerr << "error:[vg mpmap] Maximum number of rescue attempts for secondary mappings (--max-secondary-rescues) set to " << secondary_rescue_attempts << ", must set to a non-negative integer (0 for no rescue)." << endl;
        exit(1);
    }
    
    if (secondary_rescue_score_diff < 0.0) {
        cerr << "error:[vg mpmap] Max score difference for candidates clusters for secondary rescue (--secondary-diff) set to " << secondary_rescue_score_diff << ", must set to a non-negative number." << endl;
        exit(1);
    }
    
    if (max_num_mappings <= 0) {
        cerr << "error:[vg mpmap] Maximum number of mappings per read (-M) set to " << max_num_mappings << ", must set to a positive integer." << endl;
        exit(1);
    }
    
    if (reseed_length < 0) {
        cerr << "error:[vg mpmap] Reseeding length (-r) set to " << reseed_length << ", must set to a positive integer or 0 for no reseeding." << endl;
        exit(1);
    }
    
    if ((reseed_diff <= 0 || reseed_diff >= 1.0) && reseed_length != 0) {
        cerr << "error:[vg mpmap] Reseeding length difference (-W) set to " << reseed_diff << ", must set to a number between 0.0 and 1.0." << endl;
        exit(1);
    }
    
    if (reseed_exp < 0 && reseed_length != 0 && use_adaptive_reseed) {
        cerr << "error:[vg mpmap] Reseeding exponent set to " << reseed_exp << ",  must set to a nonnegative number." << endl;
        exit(1);
    }
    
    if (hard_hit_max_muliplier < 0) {
        cerr << "error:[vg mpmap] Hard MEM hit max multipler (--hard-hit-mult) set to " << hard_hit_max_muliplier << ", must set to a positive integer or 0 for no maximum." << endl;
        exit(1);
    }
    
    if (hit_max < 0) {
        cerr << "error:[vg mpmap] MEM hit max (-c) set to " << hit_max << ", must set to a positive integer or 0 for no maximum." << endl;
        exit(1);
    }
    
    if (hard_hit_max < hit_max && hit_max && hard_hit_max) {
        cerr << "warning:[vg mpmap] MEM hit query limit (-c) set to " << hit_max << ", which is higher than the threshold to ignore a MEM (" << hard_hit_max << ")" << endl;
    }
    
    if (min_mem_length < 0) {
        cerr << "error:[vg mpmap] Minimum MEM length set to " << min_mem_length << ", must set to a positive integer or 0 for no maximum." << endl;
        exit(1);
    }
    
    if (stripped_match_alg_strip_length <= 0) {
        cerr << "error:[vg mpmap] Match strip length (--strip-length) set to " << stripped_match_alg_strip_length << ", must set to a positive integer or 0 for no maximum." << endl;
        exit(1);
    }
    
    if (stripped_match_alg_max_length < 0) {
        cerr << "error:[vg mpmap] Maximum seed match length set to " << stripped_match_alg_max_length << ", must set to a positive integer or 0 for no maximum." << endl;
        exit(1);
    }
    
    if (stripped_match_alg_target_count < 0) {
        cerr << "error:[vg mpmap] Target seed count (--strip-count) set to " << stripped_match_alg_target_count << ", must set to a positive integer or 0 for no maximum." << endl;
        exit(1);
    }
    
    if (stripped_match_alg_target_count != default_strip_count && !use_stripped_match_alg) {
        cerr << "warning:[vg mpmap] Target stripped match count (--strip-count) set to " << stripped_match_alg_target_count << ", but stripped algorithm (--stripped-match) was not selected. Ignoring strip count." << endl;
    }
    
    if (stripped_match_alg_strip_length != default_strip_length && !use_stripped_match_alg) {
        cerr << "warning:[vg mpmap] Strip length (--strip-length) set to " << stripped_match_alg_strip_length << ", but stripped algorithm (--stripped-match) was not selected. Ignoring strip length." << endl;
    }
    
    if (likelihood_approx_exp < 1.0) {
        cerr << "error:[vg mpmap] Likelihood approximation exponent (--approx-exp) set to " << likelihood_approx_exp << ", must set to at least 1.0." << endl;
        exit(1);
    }
    
    if (cluster_ratio < 0.0 || cluster_ratio >= 1.0) {
        cerr << "error:[vg mpmap] Cluster drop ratio (-C) set to " << cluster_ratio << ", must set to a number between 0.0 and 1.0." << endl;
        exit(1);
    }
    
    if (use_tvs_clusterer && distance_index_name.empty()) {
        cerr << "error:[vg mpmap] The Target Value Search clusterer (-v) requires a distance index (-d)." << endl;
        exit(1);
    }
    
    if (use_min_dist_clusterer && distance_index_name.empty()) {
        cerr << "error:[vg mpmap] The minimum distance clusterer (--min-dist-cluster) requires a distance index (-d)." << endl;
        exit(1);
    }
    
    if (use_min_dist_clusterer && use_tvs_clusterer) {
        cerr << "error:[vg mpmap] Cannot perform both minimum distance clustering (--min-dist-cluster) and target value clustering (-v)." << endl;
        exit(1);
    }
    
    if (greedy_min_dist && !use_min_dist_clusterer) {
        cerr << "warning:[vg mpmap] greedy minimum distance clustering (--greedy-min-dist) is ignored if not using minimum distance clustering (-d)" << endl;
    }
    
    if (greedy_min_dist && component_min_dist) {
        cerr << "error:[vg mpmap] cannot simultaneously use greedy (--greedy-min-dist) and component (--component-min-dist) clustering" << endl;
        exit(1);
    }
    
    if (no_clustering && !distance_index_name.empty()) {
        cerr << "warning:[vg mpmap] No clustering option (--no-cluster) causes distance index (-d) to be ignored. This option is activated by default for 'very-short' read lengths (-l)." << endl;
    }
    
    if (suboptimal_path_exponent < 1.0) {
        cerr << "error:[vg mpmap] Suboptimal path likelihood root (--prune-exp) set to " << suboptimal_path_exponent << ", must set to at least 1.0." << endl;
        exit(1);
    }
    
    if (max_alignment_gap < 0) {
        cerr << "error:[vg mpmap] Max alignment grap set to " << max_alignment_gap << ", must set to a non-negative integer." << endl;
        exit(1);
    }
        
    if (buffer_size <= 0) {
        cerr << "error:[vg mpmap] Buffer size (-Z) set to " << buffer_size << ", must set to a positive integer." << endl;
        exit(1);
    }
    
    
    if ((match_score_arg != std::numeric_limits<int>::min() || mismatch_score_arg != std::numeric_limits<int>::min()) && !matrix_file_name.empty())  {
        cerr << "error:[vg mpmap] Cannot choose custom scoring matrix (-w) and custom match/mismatch score (-q/-z) simultaneously." << endl;
        exit(1);
    }
    
    if (match_score > std::numeric_limits<int8_t>::max() || mismatch_score > std::numeric_limits<int8_t>::max()
        || gap_open_score > std::numeric_limits<int8_t>::max() || gap_extension_score > std::numeric_limits<int8_t>::max()
        || full_length_bonus > std::numeric_limits<int8_t>::max() || match_score < 0 || mismatch_score < 0
        || gap_open_score < 0 || gap_extension_score < 0 || full_length_bonus < 0) {
        cerr << "error:[vg mpmap] All alignment scoring parameters (-qzoyL) must be between 0 and " << (int) std::numeric_limits<int8_t>::max() << endl;
        exit(1);
    }
    
    // ensure required parameters are provided
    
    if (graph_name.empty()) {
        cerr << "error:[vg mpmap] Multipath mapping requires a graph (-x)" << endl;
        exit(1);
    }
    
    if (gcsa_name.empty()) {
        cerr << "error:[vg mpmap] Multipath mapping requires a GCSA2 index (-g)" << endl;
        exit(1);
    }
    
    
#ifdef mpmap_instrument_mem_statitics
    if (auto_calibrate_mismapping_detection) {
        cerr << "error:[vg mpmap] set calibration off when profiling MEM statistics" << endl;
        exit(1);
    }
#endif
    
    // create in-memory objects
    
    ifstream graph_stream(graph_name);
    if (!graph_stream) {
        cerr << "error:[vg mpmap] Cannot open graph file " << graph_name << endl;
        exit(1);
    }
    
    ifstream gcsa_stream(gcsa_name);
    if (!gcsa_stream) {
        cerr << "error:[vg mpmap] Cannot open GCSA2 file " << gcsa_name << endl;
        exit(1);
    }
    
    string lcp_name = gcsa_name + ".lcp";
    ifstream lcp_stream(lcp_name);
    if (!lcp_stream) {
        cerr << "error:[vg mpmap] Cannot open LCP file " << lcp_name << endl;
        exit(1);
    }

    ifstream matrix_stream;
    if (!matrix_file_name.empty()) {
      matrix_stream.open(matrix_file_name);
      if (!matrix_stream) {
          cerr << "error:[vg mpmap] Cannot open scoring matrix file " << matrix_file_name << endl;
          exit(1);
      }
    }
    
    ifstream distance_index_stream;
    if (!distance_index_name.empty() && !no_clustering) {
        distance_index_stream.open(distance_index_name);
        if (!distance_index_stream) {
            cerr << "error:[vg mpmap] Cannot open distance index file " << distance_index_name << endl;
            exit(1);
        }
    }
    
    ifstream snarl_stream;
    if (!snarls_name.empty()) {
        snarl_stream.open(snarls_name);
        if (!snarl_stream) {
            cerr << "error:[vg mpmap] Cannot open Snarls file " << snarls_name << endl;
            exit(1);
        }
    }
    
    ifstream gbwt_stream;
    ifstream ls_stream;
    if (!gbwt_name.empty()) {
        gbwt_stream.open(gbwt_name);
        if (!gbwt_stream) {
            cerr << "error:[vg mpmap] Cannot open GBWT file " << gbwt_name << endl;
            exit(1);
        }
    }
    else if (!sublinearLS_name.empty()) {
        // We want to use sublinear Li and Stephens as our haplotype scoring approach
        ls_stream.open(sublinearLS_name);
        if (!ls_stream) {
            cerr << "error:[vg mpmap] Cannot open sublinear Li & Stevens file " << sublinearLS_name << endl;
            exit(1);
        }
    }
    
    if (!suppress_progress) {
        cerr << "[vg mpmap] Executing command:";
        for (size_t i = 0; i < argc; ++i) {
            cerr << " " << argv[i];
        }
        cerr << endl;
    }
    
    // Configure GCSA2 verbosity so it doesn't spit out loads of extra info
    gcsa::Verbosity::set(gcsa::Verbosity::SILENT);
    
    // Configure its temp directory to the system temp directory
    gcsa::TempFile::setDirectory(temp_file::get_dir());
    
    // Load required indexes
    if (!suppress_progress) {
        cerr << "[vg mpmap] Loading graph from " << graph_name << endl;
    }
    unique_ptr<PathHandleGraph> path_handle_graph = vg::io::VPKG::load_one<PathHandleGraph>(graph_stream);
    
    if (!suppress_progress) {
        // let's be a friendly guide to selecting a graph
        
        // get the graphs magic number if it has one
        uint32_t magic_num = 0;
        {
            SerializableHandleGraph* serializable = dynamic_cast<SerializableHandleGraph*>(path_handle_graph.get());
            if (serializable) {
                magic_num = serializable->get_magic_number();
            }
        }
        
        // compare to known magic numbers
        string type;
        if (magic_num == xg::XG().get_magic_number()) {
            type = "XG";
        }
        else if (magic_num == bdsg::PackedGraph().get_magic_number()) {
            type = "PackedGraph";
        }
        else if (magic_num == bdsg::HashGraph().get_magic_number()) {
            type = "HashGraph";
        }
        else if (magic_num == bdsg::ODGI().get_magic_number()) {
            type = "ODGI";
        }
        
        if (!type.empty()) {
            // we found the type, give an appropriate message about it
            cerr << "[vg mpmap] Graph is in " << type << " format. ";
            
            if (type == "XG") {
                cerr << "XG is a good graph format for most mapping use cases. PackedGraph may be selected if memory usage is too high. ";
            }
            else if (type == "HashGraph") {
                cerr << "HashGraph can have high memory usage. ";
            }
            else if (type == "PackedGraph") {
                cerr << "PackedGraph is memory efficient, but has some slow queries. ";
            }
            else if (type == "ODGI") {
                cerr << "ODGI is fairly memory efficient, but can be slow on certain queries. ";
            }
        }
        else {
            // probably a VG graph
            cerr << "[vg mpmap] Graph is not in XG format. " << endl;
        }
        
        // are they using a graph combo that I don't recommend?
        if (type != "XG" && (!use_min_dist_clusterer || type == "HashGraph" || type.empty())) {
            // min dist clustering alleviates the issues with slow path queries because we don't need to do
            // so many, but I want to dissuade people from using HashGraph and VG for mapping regardless
            cerr << "XG format is recommended for most mapping tasks. ";
        }
        
        cerr << "See `vg convert` if you want to change graph formats." << endl;
    }
    
    if (path_handle_graph->get_path_count() == 0 && distance_index_name.empty()) {
        cerr << "warning:[vg mpmap] Using a distance index (-d) for clustering is highly recommended for graphs that lack embedded paths. Speed and accuracy are likely to suffer severely without one." << endl;
    }
    
    bdsg::PathPositionOverlayHelper overlay_helper;
    PathPositionHandleGraph* path_position_handle_graph = overlay_helper.apply(path_handle_graph.get());
    
    if (!suppress_progress) {
        cerr << "[vg mpmap] Loading GCSA2 from " << gcsa_name << endl;
    }
    unique_ptr<gcsa::GCSA> gcsa_index = vg::io::VPKG::load_one<gcsa::GCSA>(gcsa_stream);
    unique_ptr<gcsa::LCPArray> lcp_array;
    if (!use_stripped_match_alg) {
        // The stripped algorithm doesn't use the LCP, but we aren't doing it
        if (!suppress_progress) {
            cerr << "[vg mpmap] Loading LCP from " << lcp_name << endl;
        }
        lcp_array = vg::io::VPKG::load_one<gcsa::LCPArray>(lcp_stream);
    }
    
    // Load optional indexes
    
    unique_ptr<gbwt::GBWT> gbwt;
    haplo::linear_haplo_structure* sublinearLS = nullptr;
    haplo::ScoreProvider* haplo_score_provider = nullptr;
    if (!gbwt_name.empty()) {
        if (!suppress_progress) {
            cerr << "[vg mpmap] Loading GBWT from " << gbwt_name << endl;
        }
        // Load the GBWT from its container
        gbwt = vg::io::VPKG::load_one<gbwt::GBWT>(gbwt_stream);

        if (gbwt.get() == nullptr) {
          // Complain if we couldn't.
          cerr << "error:[vg mpmap] unable to load gbwt index file" << endl;
          exit(1);
        }
    
        // We have the GBWT available for scoring haplotypes
        haplo_score_provider = new haplo::GBWTScoreProvider<gbwt::GBWT>(*gbwt);
    }
    else if (!sublinearLS_name.empty()) {
        if (!suppress_progress) {
            cerr << "[vg mpmap] Loading LS index from " << sublinearLS_name << endl;
        }
        
        // TODO: we only support a single ref contig, and we use these
        // hardcoded mutation and recombination likelihoods
        
        sublinearLS = new linear_haplo_structure(ls_stream, -9 * 2.3, -6 * 2.3, *path_position_handle_graph,
                                                 path_position_handle_graph->get_path_handle(sublinearLS_ref_path));
        haplo_score_provider = new haplo::LinearScoreProvider(*sublinearLS);
    }
    
    unique_ptr<SnarlManager> snarl_manager;
    if (!snarls_name.empty()) {
        if (!suppress_progress) {
            cerr << "[vg mpmap] Loading snarls from " << snarls_name << endl;
        }
        snarl_manager = vg::io::VPKG::load_one<SnarlManager>(snarl_stream);
    }
    
    unique_ptr<MinimumDistanceIndex> distance_index;
    if (!distance_index_name.empty() && !no_clustering) {
        if (!suppress_progress) {
            cerr << "[vg mpmap] Loading distance index from " << distance_index_name << endl;
        }
        
        // Load the index
        distance_index = vg::io::VPKG::load_one<MinimumDistanceIndex>(distance_index_stream);
        
    }
    
    MultipathMapper multipath_mapper(path_position_handle_graph, gcsa_index.get(), lcp_array.get(), haplo_score_provider,
        snarl_manager.get(), distance_index.get());
    
    // set alignment parameters
    if (matrix_stream.is_open()) {
        multipath_mapper.set_alignment_scores(matrix_stream, gap_open_score, gap_extension_score, full_length_bonus);
    }
    else {
        multipath_mapper.set_alignment_scores(match_score, mismatch_score, gap_open_score, gap_extension_score, full_length_bonus);
    }
    multipath_mapper.adjust_alignments_for_base_quality = qual_adjusted;
    multipath_mapper.strip_bonuses = strip_full_length_bonus;
    multipath_mapper.band_padding_multiplier = band_padding_multiplier;
    multipath_mapper.init_band_padding_memo();
    
    // set mem finding parameters
    multipath_mapper.hit_max = hit_max;
    multipath_mapper.hard_hit_max = hard_hit_max;
    multipath_mapper.mem_reseed_length = reseed_length;
    multipath_mapper.fast_reseed = true;
    multipath_mapper.fast_reseed_length_diff = reseed_diff;
    multipath_mapper.sub_mem_count_thinning = sub_mem_count_thinning;
    multipath_mapper.sub_mem_thinning_burn_in = sub_mem_thinning_burn_in;
    multipath_mapper.order_length_repeat_hit_max = order_length_repeat_hit_max;
    multipath_mapper.min_mem_length = min_mem_length;
    multipath_mapper.stripped_match_alg_strip_length = stripped_match_alg_strip_length;
    multipath_mapper.stripped_match_alg_max_length = stripped_match_alg_max_length;
    multipath_mapper.stripped_match_alg_target_count = stripped_match_alg_target_count;
    multipath_mapper.use_greedy_mem_restarts = use_greedy_mem_restarts;
    multipath_mapper.greedy_restart_min_length = greedy_restart_min_length;
    multipath_mapper.greedy_restart_max_count = greedy_restart_max_count;
    multipath_mapper.greedy_restart_max_lcp = greedy_restart_max_lcp;
    multipath_mapper.greedy_restart_assume_substitution = greedy_restart_assume_substitution;
    multipath_mapper.use_stripped_match_alg = use_stripped_match_alg;
    multipath_mapper.adaptive_reseed_diff = use_adaptive_reseed;
    multipath_mapper.adaptive_diff_exponent = reseed_exp;
    multipath_mapper.use_approx_sub_mem_count = false;
    multipath_mapper.prefilter_redundant_hits = prefilter_redundant_hits;
    multipath_mapper.precollapse_order_length_hits = precollapse_order_length_hits;
    multipath_mapper.max_sub_mem_recursion_depth = max_sub_mem_recursion_depth;
    multipath_mapper.max_mapping_p_value = max_mapping_p_value;
    if (min_clustering_mem_length) {
        multipath_mapper.min_clustering_mem_length = min_clustering_mem_length;
    }
    else {
        multipath_mapper.set_automatic_min_clustering_length();
    }
    
    // set mapping quality parameters
    multipath_mapper.mapping_quality_method = mapq_method;
    multipath_mapper.max_mapping_quality = max_mapq;
    multipath_mapper.report_group_mapq = report_group_mapq;
    multipath_mapper.use_weibull_calibration = use_weibull_calibration;
    // Use population MAPQs when we have the right option combination to make that sensible.
    multipath_mapper.use_population_mapqs = (haplo_score_provider != nullptr && population_max_paths > 0);
    multipath_mapper.population_max_paths = population_max_paths;
    multipath_mapper.population_paths_hard_cap = population_paths_hard_cap;
    multipath_mapper.top_tracebacks = top_tracebacks;
    multipath_mapper.recombination_penalty = recombination_penalty;
    multipath_mapper.always_check_population = always_check_population;
    multipath_mapper.force_haplotype_count = force_haplotype_count;
    
    // set pruning and clustering parameters
    multipath_mapper.no_clustering = no_clustering;
    multipath_mapper.use_tvs_clusterer = use_tvs_clusterer;
    multipath_mapper.use_min_dist_clusterer = use_min_dist_clusterer;
    multipath_mapper.greedy_min_dist = greedy_min_dist;
    multipath_mapper.component_min_dist = component_min_dist;
    multipath_mapper.max_expected_dist_approx_error = max_dist_error;
    multipath_mapper.mem_coverage_min_ratio = cluster_ratio;
    multipath_mapper.log_likelihood_approx_factor = likelihood_approx_exp;
    multipath_mapper.num_mapping_attempts = max_map_attempts;
    multipath_mapper.min_median_mem_coverage_for_split = min_median_mem_coverage_for_split;
    multipath_mapper.suppress_cluster_merging = suppress_cluster_merging;
    multipath_mapper.use_tvs_clusterer = use_tvs_clusterer;
    multipath_mapper.reversing_walk_length = reversing_walk_length;
    multipath_mapper.max_alt_mappings = max_num_mappings;
    multipath_mapper.max_alignment_gap = max_alignment_gap;
    multipath_mapper.pessimistic_tail_gap_multiplier = pessimistic_tail_gap_multiplier;
    
    // set pair rescue parameters
    multipath_mapper.max_rescue_attempts = max_rescue_attempts;
    multipath_mapper.max_single_end_mappings_for_rescue = max(max_single_end_mappings_for_rescue, max_rescue_attempts);
    multipath_mapper.secondary_rescue_subopt_diff = secondary_rescue_subopt_diff;
    multipath_mapper.secondary_rescue_score_diff = secondary_rescue_score_diff;
    multipath_mapper.secondary_rescue_attempts = secondary_rescue_attempts;
    multipath_mapper.rescue_only_min = rescue_only_min;
    multipath_mapper.rescue_only_anchor_max = rescue_only_anchor_max;
    multipath_mapper.fragment_length_warning_factor = fragment_length_warning_factor;
    
    // set multipath alignment topology parameters
    multipath_mapper.max_snarl_cut_size = snarl_cut_size;
    multipath_mapper.max_branch_trim_length = max_branch_trim_length;
    multipath_mapper.suppress_tail_anchors = !synthesize_tail_anchors;
    multipath_mapper.num_alt_alns = num_alt_alns;
    multipath_mapper.dynamic_max_alt_alns = dynamic_max_alt_alns;
    multipath_mapper.simplify_topologies = simplify_topologies;
    multipath_mapper.max_suboptimal_path_score_ratio = suboptimal_path_exponent;

#ifdef mpmap_instrument_mem_statitics
    multipath_mapper._mem_stats.open(MEM_STATS_FILE);
#endif
    
    // if directed to, auto calibrate the mismapping detection to the graph
    if (auto_calibrate_mismapping_detection) {
        if (!suppress_progress) {
            cerr << "[vg mpmap] Building null model to calibrate mismapping detection (can take some time)." << endl;
        }
        multipath_mapper.calibrate_mismapping_detection(num_calibration_simulations, calibration_read_lengths);
    }
    
    
    // Count our threads 
    int thread_count = get_thread_count();
    
    // Establish a watchdog to find reads that take too long to map.
    // If we see any, we will issue a warning.
    unique_ptr<Watchdog> watchdog(new Watchdog(thread_count, chrono::minutes(read_length == "long" ? 40 : 5)));
    
    // are we doing paired ends?
    if (interleaved_input || !fastq_name_2.empty()) {
        // make sure buffer size is even (ensures that output will be interleaved)
        if (buffer_size % 2 == 1) {
            buffer_size++;
        }

        if (!std::isnan(frag_length_mean) && !std::isnan(frag_length_stddev)) {
            // Force a fragment length distribution
            multipath_mapper.force_fragment_length_distr(frag_length_mean, frag_length_stddev);
        }
        else {
            // choose the sample size and tail-fraction for estimating the fragment length distribution
            multipath_mapper.set_fragment_length_distr_params(frag_length_sample_size, frag_length_sample_size,
                                                              frag_length_robustness_fraction);
        }

    }
    
#ifdef record_read_run_times
    ofstream read_time_file(READ_TIME_FILE);
#endif
    
    // a probably over-engineered way to report progress across threads with minimal contention
    const uint64_t progress_frequency = read_length == "long" ? 250000 : 5000000;
    const uint64_t thread_progress_frequency = 1000;
    assert(progress_frequency % thread_progress_frequency == 0);
    uint64_t num_reads_mapped = 0;
    vector<uint64_t> thread_num_reads_mapped(thread_count, 0);
    
    function<void(int)> register_mapping = [&](int thread_num) {
        if (!suppress_progress) {
            uint64_t num_mapped = ++thread_num_reads_mapped[thread_num];
            if (num_mapped == thread_progress_frequency) {
                uint64_t n;
#pragma omp atomic capture
                n = num_reads_mapped += num_mapped;
                if (n % progress_frequency == 0) {
#pragma omp critical
                    {
                        cerr << "[vg mpmap] Mapped " << n << (!interleaved_input && fastq_name_2.empty() ? " reads" : " read pairs") << endl;
                    }
                }
                thread_num_reads_mapped[thread_num] = 0;
            }
        }
    };
    
    // init a writer for the output
    MultipathAlignmentEmitter* emitter = new MultipathAlignmentEmitter(cout, thread_count, single_path_alignment_mode);
    
    // a buffer to hold read pairs that can't be unambiguously mapped before the fragment length distribution
    // is estimated
    // note: sufficient to have only one buffer because multithreading code enforces single threaded mode
    // during distribution estimation
    vector<pair<Alignment, Alignment>> ambiguous_pair_buffer;
    
    // do unpaired multipath alignment and write to buffer
    function<void(Alignment&)> do_unpaired_alignments = [&](Alignment& alignment) {
#ifdef record_read_run_times
        clock_t start = clock();
#endif

        auto thread_num = omp_get_thread_num();

        if (watchdog) {
            watchdog->check_in(thread_num, alignment.name());
        }
        
        bool is_rna = uses_Us(alignment);
        if (is_rna) {
            convert_Us_to_Ts(alignment);
        }

        vector<multipath_alignment_t> mp_alns;
        multipath_mapper.multipath_map(alignment, mp_alns);
        
        if (is_rna) {
            for (multipath_alignment_t& mp_aln : mp_alns) {
                convert_Ts_to_Us(mp_aln);
            }
        }
        
        if (!no_output) {
            emitter->emit_singles(move(mp_alns));
        }
        
        if (watchdog) {
            watchdog->check_out(thread_num);
        }
        
        register_mapping(thread_num);
        
#ifdef record_read_run_times
        clock_t finish = clock();
#pragma omp critical
        read_time_file << alignment.name() << "\t" << alignment.sequence().size() << "\t" << double(finish - start) / CLOCKS_PER_SEC << endl;
#endif
    };
    
    // do paired multipath alignment and write to buffer
    function<void(Alignment&, Alignment&)> do_paired_alignments = [&](Alignment& alignment_1, Alignment& alignment_2) {
        // get reads on the same strand so that oriented distance estimation works correctly
        // but if we're clearing the ambiguous buffer we already RC'd these on the first pass

        auto thread_num = omp_get_thread_num();

#ifdef record_read_run_times
        clock_t start = clock();
#endif

        if (watchdog) {
            watchdog->check_in(thread_num, alignment_1.name());
        }
        
        bool is_rna = (uses_Us(alignment_1) || uses_Us(alignment_2));
        if (is_rna) {
            convert_Us_to_Ts(alignment_1);
            convert_Us_to_Ts(alignment_2);
        }
        
        if (!same_strand) {
            // remove the path so we won't try to RC it (the path may not refer to this graph)
            alignment_2.clear_path();
            reverse_complement_alignment_in_place(&alignment_2, [&](vg::id_t node_id) {
                return path_position_handle_graph->get_length(path_position_handle_graph->get_handle(node_id));
            });
        }
        
        size_t num_buffered = ambiguous_pair_buffer.size();
        
        vector<pair<multipath_alignment_t, multipath_alignment_t>> mp_aln_pairs;
        multipath_mapper.multipath_map_paired(alignment_1, alignment_2, mp_aln_pairs, ambiguous_pair_buffer);
        
        
        if (!same_strand) {
            for (auto& mp_aln_pair : mp_aln_pairs) {
                rev_comp_multipath_alignment_in_place(&mp_aln_pair.second, [&](vg::id_t node_id) { return path_position_handle_graph->get_length(path_position_handle_graph->get_handle(node_id));
                });
            }
        }
        
        if (is_rna) {
            for (pair<multipath_alignment_t, multipath_alignment_t>& mp_aln_pair : mp_aln_pairs) {
                convert_Ts_to_Us(mp_aln_pair.first);
                convert_Ts_to_Us(mp_aln_pair.second);
            }
        }
        
        if (!no_output) {
            emitter->emit_pairs(move(mp_aln_pairs));
        }
        
        if (watchdog) {
            watchdog->check_out(thread_num);
        }
        
        if (num_buffered == ambiguous_pair_buffer.size()) {
            // the read didn't get buffered during the frag length estimation phase
            register_mapping(thread_num);
        }
        
#ifdef record_read_run_times
        clock_t finish = clock();
#pragma omp critical
        read_time_file << alignment_1.name() << "\t" << alignment_2.name() << "\t" << alignment_1.sequence().size() << "\t" << alignment_2.sequence().size() << "\t" << double(finish - start) / CLOCKS_PER_SEC << endl;
#endif
    };
    
    // do unpaired, independent multipath alignment, and write to buffer as paired
    function<void(Alignment&, Alignment&)> do_independent_paired_alignments = [&](Alignment& alignment_1, Alignment& alignment_2) {
        // get reads on the same strand so that oriented distance estimation works correctly
        // but if we're clearing the ambiguous buffer we already RC'd these on the first pass

        auto thread_num = omp_get_thread_num();

#ifdef record_read_run_times
        clock_t start = clock();
#endif

        if (watchdog) {
            watchdog->check_in(thread_num, alignment_1.name());
        }
        
        bool is_rna = (uses_Us(alignment_1) || uses_Us(alignment_2));
        if (is_rna) {
            convert_Us_to_Ts(alignment_1);
            convert_Us_to_Ts(alignment_2);
        }
        
        // Align independently
        vector<multipath_alignment_t> mp_alns_1, mp_alns_2;
        multipath_mapper.multipath_map(alignment_1, mp_alns_1);
        multipath_mapper.multipath_map(alignment_2, mp_alns_2);
        
        if (is_rna) {
            for (multipath_alignment_t& mp_aln : mp_alns_1) {
                convert_Ts_to_Us(mp_aln);
            }
            for (multipath_alignment_t& mp_aln : mp_alns_2) {
                convert_Ts_to_Us(mp_aln);
            }
        }
            
        // keep an equal number to protect interleaving
        mp_alns_1.resize(min(mp_alns_1.size(), mp_alns_2.size()));
        mp_alns_2.resize(min(mp_alns_1.size(), mp_alns_2.size()));
        
        // combine into a single vector
        vector<multipath_alignment_t> mp_alns_combined;
        mp_alns_combined.reserve(mp_alns_1.size() + mp_alns_2.size());
        for (size_t i = 0; i < mp_alns_1.size(); ++i) {
            mp_alns_combined.emplace_back(std::move(mp_alns_1[i]));
            mp_alns_combined.emplace_back(std::move(mp_alns_2[i]));
        }
        
        if (!no_output) {
            emitter->emit_singles(move(mp_alns_combined));
        }
        
        if (watchdog) {
            watchdog->check_out(thread_num);
        }
        
        register_mapping(thread_num);
        
#ifdef record_read_run_times
        clock_t finish = clock();
#pragma omp critical
        read_time_file << alignment_1.name() << "\t" << alignment_2.name() << "\t" << alignment_1.sequence().size() << "\t" << alignment_2.sequence().size() << "\t" << double(finish - start) / CLOCKS_PER_SEC << endl;
#endif
    };
    
    // for streaming paired input, don't spawn parallel tasks unless this evalutes to true
    function<bool(void)> multi_threaded_condition = [&](void) {
        return multipath_mapper.has_fixed_fragment_length_distr();
    };
    
    // FASTQ input
    if (!fastq_name_1.empty()) {
        if (!suppress_progress) {
            cerr << "[vg mpmap] Mapping reads from " << (fastq_name_1 == "-" ? "STDIN" : fastq_name_1) << (fastq_name_2.empty() ? "" : " and " + (fastq_name_2 == "-" ? "STDIN" : fastq_name_2)) << " using " << thread_count << " threads" << endl;
        }
        
        if (interleaved_input) {
            fastq_paired_interleaved_for_each_parallel_after_wait(fastq_name_1, do_paired_alignments,
                                                                  multi_threaded_condition);
        }
        else if (fastq_name_2.empty()) {
            fastq_unpaired_for_each_parallel(fastq_name_1, do_unpaired_alignments);
        }
        else {
            fastq_paired_two_files_for_each_parallel_after_wait(fastq_name_1, fastq_name_2, do_paired_alignments,
                                                                multi_threaded_condition);
        }
    }
    
    // GAM input
    if (!gam_file_name.empty()) {
        if (!suppress_progress) {
            cerr << "[vg mpmap] Mapping reads from " << (gam_file_name == "-" ? "STDIN" : gam_file_name) << " using " << thread_count << " threads" << endl;
        }
        
        function<void(istream&)> execute = [&](istream& gam_in) {
            if (!gam_in) {
                cerr << "error:[vg mpmap] Cannot open GAM file " << gam_file_name << endl;
                exit(1);
            }
            
            if (interleaved_input) {
                vg::io::for_each_interleaved_pair_parallel_after_wait(gam_in, do_paired_alignments,
                                                                      multi_threaded_condition);
            }
            else {
                vg::io::for_each_parallel(gam_in, do_unpaired_alignments);
            }
        };
        get_input_file(gam_file_name, execute);
    }

    // take care of any read pairs that we couldn't map unambiguously before the fragment length distribution
    // had been estimated
    if (!ambiguous_pair_buffer.empty()) {
        if (multipath_mapper.has_fixed_fragment_length_distr()) {
#pragma omp parallel for
            for (size_t i = 0; i < ambiguous_pair_buffer.size(); i++) {
                pair<Alignment, Alignment>& aln_pair = ambiguous_pair_buffer[i];
                // we reverse complemented the alignment on the first pass, so switch back so we don't break
                // the alignment functions expectations
                // TODO: slightly wasteful, inelegant
                if (!same_strand) {
                    reverse_complement_alignment_in_place(&aln_pair.second,
                                                          [&](vg::id_t node_id) {
                        return path_position_handle_graph->get_length(path_position_handle_graph->get_handle(node_id));
                    });
                }
                do_paired_alignments(aln_pair.first, aln_pair.second);
            }
        }
        else {
            cerr << "warning:[vg mpmap] Could not find " << frag_length_sample_size << " (-b) unambiguous read pair mappings to estimate fragment length ditribution. This can happen due to data issues (e.g. unpaired reads being mapped as pairs) or because the sample size is too large for the read set. Mapping read pairs as independent single-ended reads." << endl;
            
#pragma omp parallel for
            for (size_t i = 0; i < ambiguous_pair_buffer.size(); i++) {
                pair<Alignment, Alignment>& aln_pair = ambiguous_pair_buffer[i];
                // we reverse complemented the alignment on the first pass, so switch back so we don't break
                // the alignment function's expectations
                // TODO: slightly wasteful, inelegant
                if (!same_strand) {
                    reverse_complement_alignment_in_place(&aln_pair.second,
                                                          [&](vg::id_t node_id) {
                        return path_position_handle_graph->get_length(path_position_handle_graph->get_handle(node_id));
                    });
                }
                do_independent_paired_alignments(aln_pair.first, aln_pair.second);
            }
        }
    }
    
    // flush output
    delete emitter;
    cout.flush();
    
    if (!suppress_progress) {
        for (auto uncounted_mappings : thread_num_reads_mapped) {
            num_reads_mapped += uncounted_mappings;
        }
        cerr << "[vg mpmap] Mapping finished. Mapped " << num_reads_mapped;
        if (fastq_name_2.empty() && !interleaved_input) {
            cerr << " reads.";
        }
        else {
            cerr << " read pairs.";
        }
        cerr << endl;
    }
    
#ifdef record_read_run_times
    read_time_file.close();
#endif
    
    //cerr << "MEM length filtering efficiency: " << ((double) OrientedDistanceClusterer::MEM_FILTER_COUNTER) / OrientedDistanceClusterer::MEM_TOTAL << " (" << OrientedDistanceClusterer::MEM_FILTER_COUNTER << "/" << OrientedDistanceClusterer::MEM_TOTAL << ")" << endl;
    //cerr << "MEM cluster filtering efficiency: " << ((double) OrientedDistanceClusterer::PRUNE_COUNTER) / OrientedDistanceClusterer::CLUSTER_TOTAL << " (" << OrientedDistanceClusterer::PRUNE_COUNTER << "/" << OrientedDistanceClusterer::CLUSTER_TOTAL << ")" << endl;
    //cerr << "subgraph filtering efficiency: " << ((double) MultipathMapper::PRUNE_COUNTER) / MultipathMapper::SUBGRAPH_TOTAL << " (" << MultipathMapper::PRUNE_COUNTER << "/" << MultipathMapper::SUBGRAPH_TOTAL << ")" << endl;
    //cerr << "attempted to split " << OrientedDistanceClusterer::SPLIT_ATTEMPT_COUNTER << " of " << OrientedDistanceClusterer::PRE_SPLIT_CLUSTER_COUNTER << " clusters with " << OrientedDistanceClusterer::SUCCESSFUL_SPLIT_ATTEMPT_COUNTER << " splits successful (" << 100.0 * double(OrientedDistanceClusterer::SUCCESSFUL_SPLIT_ATTEMPT_COUNTER) / OrientedDistanceClusterer::SPLIT_ATTEMPT_COUNTER << "%) resulting in " << OrientedDistanceClusterer::POST_SPLIT_CLUSTER_COUNTER << " total clusters (" << OrientedDistanceClusterer::POST_SPLIT_CLUSTER_COUNTER - OrientedDistanceClusterer::PRE_SPLIT_CLUSTER_COUNTER << " new)" << endl;
    //cerr << "entered secondary rescue " << MultipathMapper::SECONDARY_RESCUE_TOTAL << " times with " << MultipathMapper::SECONDARY_RESCUE_COUNT << " actually attempting rescues, totaling " << MultipathMapper::SECONDARY_RESCUE_ATTEMPT << " rescues (" << double(MultipathMapper::SECONDARY_RESCUE_ATTEMPT) / MultipathMapper::SECONDARY_RESCUE_COUNT << " average per attempt)" << endl;
    
    if (haplo_score_provider != nullptr) {
        delete haplo_score_provider;
    }
    
    if (sublinearLS != nullptr) {
        delete sublinearLS;
    }
   
    return 0;
}

// Register subcommand
static Subcommand vg_mpmap("mpmap", "multipath alignments of reads to a graph", main_mpmap);


