/** \file sim_main.cpp
 *
 * Defines the "vg sim" subcommand, which generates potential reads from a graph.
 */


#include <omp.h>
#include <unistd.h>
#include <getopt.h>

#include <list>
#include <fstream>

#include "subcommand.hpp"

#include "../vg.hpp"
#include "../xg.hpp"
#include "../aligner.hpp"
#include "../sampler.hpp"
#include <vg/io/protobuf_emitter.hpp>
#include <vg/io/vpkg.hpp>
#include <bdsg/overlay_helper.hpp>

using namespace std;
using namespace vg;
using namespace vg::subcommand;

// Gets the transcript IDs and TPM values from an RSEM output .tsv file
vector<pair<string, double>> parse_rsem_expression_file(istream& rsem_in) {
    vector<pair<string, double>> return_val;
    string line;
    // skip the header line
    getline(rsem_in, line);
    line.clear();
    while (getline(rsem_in, line)) {
        vector<string> tokens;
        stringstream strm(line);
        string token;
        while (getline(strm, token, '\t')) {
            tokens.push_back(move(token));
            token.clear();
        }
        if (tokens.size() != 8) {
            cerr << "[vg sim] error: Cannot parse transcription file. Expected 8-column TSV file as produced by RSEM, got " << tokens.size() << " columns." << endl;
            exit(1);
        }
        return_val.emplace_back(tokens[0], parse<double>(tokens[5]));
        line.clear();
    }
    return return_val;
}

void help_sim(char** argv) {
    cerr << "usage: " << argv[0] << " sim [options]" << endl
         << "Samples sequences from the xg-indexed graph." << endl
         << endl
         << "options:" << endl
         << "    -x, --xg-name FILE          use the xg index (or graph) in FILE" << endl
         << "    -F, --fastq FILE            superpose errors matching the error profile of NGS reads in FILE (ignores -l,-f)" << endl
         << "    -I, --interleaved           reads in FASTQ (-F) are interleaved read pairs" << endl
         << "    -P, --path PATH             simulate from the given names path (multiple allowed, cannot also give -T)" << endl
         << "    -T, --tx-expr-file FILE     simulate from an expression profile formatted as RSEM output (cannot also give -P)" << endl
         << "    -l, --read-length N         write reads of length N" << endl
         << "    -n, --num-reads N           simulate N reads or read pairs" << endl
         << "    -s, --random-seed N         use this specific seed for the PRNG" << endl
         << "    -e, --sub-rate FLOAT        base substitution rate (default 0.0)" << endl
         << "    -i, --indel-rate FLOAT      indel rate (default 0.0)" << endl
         << "    -d, --indel-err-prop FLOAT  proportion of trained errors from -F that are indels (default 0.0)" << endl
         << "    -S, --scale-err FLOAT       scale trained error probabilities from -F by this much (default 1.0)" << endl
         << "    -f, --forward-only          don't simulate from the reverse strand" << endl
         << "    -p, --frag-len N            make paired end reads with given fragment length N" << endl
         << "    -v, --frag-std-dev FLOAT    use this standard deviation for fragment length estimation" << endl
         << "    -N, --allow-Ns              allow reads to be sampled from the graph with Ns in them" << endl
         << "    -a, --align-out             generate true alignments on stdout rather than reads" << endl
         << "    -J, --json-out              write alignments in json" << endl;
}

int main_sim(int argc, char** argv) {

    if (argc == 2) {
        help_sim(argv);
        return 1;
    }

    int read_length = 100;
    int num_reads = 1;
    int seed_val = time(NULL);
    double base_error = 0;
    double indel_error = 0;
    bool forward_only = false;
    bool align_out = false;
    bool json_out = false;
    int fragment_length = 0;
    double fragment_std_dev = 0;
    bool reads_may_contain_Ns = false;
    string xg_name;
    bool strip_bonuses = false;
    bool interleaved = false;
    double indel_prop = 0.0;
    double error_scale_factor = 1.0;
    string fastq_name;
    // What path should we sample from? Empty string = the whole graph.
    vector<string> path_names;
    // Alternatively, which transcripts with how much expression?
    string rsem_file_name;
    vector<pair<string, double>> transcript_expressions;

    int c;
    optind = 2; // force optind past command positional argument
    while (true) {
        static struct option long_options[] =
        {
            {"help", no_argument, 0, 'h'},
            {"xg-name", required_argument, 0, 'x'},
            {"fastq", required_argument, 0, 'F'},
            {"interleaved", no_argument, 0, 'I'},
            {"path", required_argument, 0, 'P'},
            {"tx-expr-file", required_argument, 0, 'T'},
            {"read-length", required_argument, 0, 'l'},
            {"num-reads", required_argument, 0, 'n'},
            {"random-seed", required_argument, 0, 's'},
            {"forward-only", no_argument, 0, 'f'},
            {"align-out", no_argument, 0, 'a'},
            {"json-out", no_argument, 0, 'J'},
            {"allow-Ns", no_argument, 0, 'N'},
            {"sub-rate", required_argument, 0, 'e'},
            {"indel-rate", required_argument, 0, 'i'},
            {"indel-err-prop", required_argument, 0, 'd'},
            {"scale-err", required_argument, 0, 'S'},
            {"frag-len", required_argument, 0, 'p'},
            {"frag-std-dev", required_argument, 0, 'v'},
            {0, 0, 0, 0}
        };

        int option_index = 0;
        c = getopt_long (argc, argv, "hl:n:s:e:i:fax:Jp:v:Nd:F:P:T:S:I",
                long_options, &option_index);

        // Detect the end of the options.
        if (c == -1)
            break;

        switch (c)
        {

        case 'x':
            xg_name = optarg;
            break;
            
        case 'F':
            fastq_name = optarg;
            break;
            
        case 'I':
            interleaved = true;
            break;
            
        case 'P':
            path_names.push_back(optarg);
            break;
            
        case 'T':
            rsem_file_name = optarg;
            break;

        case 'l':
            read_length = parse<int>(optarg);
            break;

        case 'n':
            num_reads = parse<int>(optarg);
            break;

        case 's':
            seed_val = parse<int>(optarg);
            if (seed_val == 0) {
                // Don't let the user specify seed 0 as we will confuse it with no deterministic seed.
                cerr << "error[vg sim]: seed 0 cannot be used. Omit the seed option if you want nondeterministic results." << endl;
                exit(1);
            }
            break;

        case 'e':
            base_error = parse<double>(optarg);
            break;

        case 'i':
            indel_error = parse<double>(optarg);
            break;
            
        case 'd':
            indel_prop = parse<double>(optarg);
            break;
            
        case 'S':
            error_scale_factor = parse<double>(optarg);
            break;

        case 'f':
            forward_only = true;
            break;

        case 'a':
            align_out = true;
            break;

        case 'J':
            json_out = true;
            align_out = true;
            break;

        case 'N':
            reads_may_contain_Ns = true;
            break;

        case 'p':
            fragment_length = parse<int>(optarg);
            break;

        case 'v':
            fragment_std_dev = parse<double>(optarg);
            break;
            
        case 'h':
        case '?':
            help_sim(argv);
            exit(1);
            break;

        default:
            abort ();
        }
    }

    if (xg_name.empty()) {
        cerr << "[vg sim] error: we need an xg index to sample reads from" << endl;
        return 1;
    }
    
    if (!rsem_file_name.empty()) {
        ifstream rsem_in(rsem_file_name);
        if (!rsem_in) {
            cerr << "[vg sim] error: could not open transcription profile file " << rsem_file_name << endl;
            return 1;
        }
        transcript_expressions = parse_rsem_expression_file(rsem_in);
    }

    unique_ptr<PathHandleGraph> path_handle_graph;
    bdsg::PathPositionVectorizableOverlayHelper overlay_helper;
    PathPositionHandleGraph* xgidx = nullptr;

    if (!xg_name.empty()) {
        path_handle_graph = vg::io::VPKG::load_one<PathHandleGraph>(xg_name);
        xgidx = dynamic_cast<PathPositionHandleGraph*>(overlay_helper.apply(path_handle_graph.get()));
    }
    
    for (auto& path_name : path_names) {
        if (xgidx->has_path(path_name) == false) {
            cerr << "[vg sim] error: path \""<< path_name << "\" not found in index" << endl;
            return 1;
        }
    }
    
    for (auto& transcript_expression : transcript_expressions) {
        if (xgidx->has_path(transcript_expression.first) == false) {
            cerr << "[vg sim] error: transcript path for \""<< transcript_expression.first << "\" not found in index" << endl;
            return 1;
        }
    }
    
    unique_ptr<vg::io::ProtobufEmitter<Alignment>> aln_emitter;
    if (align_out && !json_out) {
        // Make an emitter to emit Alignments
        aln_emitter = unique_ptr<vg::io::ProtobufEmitter<Alignment>>(new vg::io::ProtobufEmitter<Alignment>(cout));
    }
    
    if (fastq_name.empty()) {
        // Use the fixed error rate sampler
        
        // Make a sample to sample reads with
        Sampler sampler(xgidx, seed_val, forward_only, reads_may_contain_Ns, path_names, transcript_expressions);
        
        Aligner rescorer(default_match, default_mismatch, default_gap_open, default_gap_extension, default_full_length_bonus);

        // We define a function to score a using the aligner
        auto rescore = [&] (Alignment& aln) {
            // Score using exact distance.
            aln.set_score(rescorer.score_ungapped_alignment(aln, strip_bonuses));
        };
        
        size_t max_iter = 1000;
        int nonce = 1;
        for (int i = 0; i < num_reads; ++i) {
            // For each read we are going to generate
            
            if (fragment_length) {
                // fragment_lenght is nonzero so make it two paired reads
                auto alns = sampler.alignment_pair(read_length, fragment_length, fragment_std_dev, base_error, indel_error);
                
                size_t iter = 0;
                while (iter++ < max_iter) {
                    // For up to max_iter iterations
                    if (alns.front().sequence().size() < read_length
                        || alns.back().sequence().size() < read_length) {
                        // If our read was too short, try again
                        alns = sampler.alignment_pair(read_length, fragment_length, fragment_std_dev, base_error, indel_error);
                    }
                }
                
                // write the alignment or its string
                if (align_out) {
                    // write it out as requested
                    
                    // We will need scores
                    rescore(alns.front());
                    rescore(alns.back());
                    
                    if (json_out) {
                        cout << pb2json(alns.front()) << endl;
                        cout << pb2json(alns.back()) << endl;
                    } else {
                        aln_emitter->write_copy(alns.front());
                        aln_emitter->write_copy(alns.back());
                    }
                } else {
                    cout << alns.front().sequence() << "\t" << alns.back().sequence() << endl;
                }
            } else {
                // Do single-end reads
                auto aln = sampler.alignment_with_error(read_length, base_error, indel_error);
                
                size_t iter = 0;
                while (iter++ < max_iter) {
                    // For up to max_iter iterations
                    if (aln.sequence().size() < read_length) {
                        // If our read is too short, try again
                        auto aln_prime = sampler.alignment_with_error(read_length, base_error, indel_error);
                        if (aln_prime.sequence().size() > aln.sequence().size()) {
                            // But only keep the new try if it is longer
                            aln = aln_prime;
                        }
                    }
                }
                
                // write the alignment or its string
                if (align_out) {
                    // write it out as requested
                    
                    // We will need scores
                    rescore(aln);
                    
                    if (json_out) {
                        cout << pb2json(aln) << endl;
                    } else {
                        aln_emitter->write_copy(aln);
                    }
                } else {
                    cout << aln.sequence() << endl;
                }
            }
        }
        
    }
    else {
        // Use the trained error rate
        
        Aligner aligner(default_match, default_mismatch, default_gap_open, default_gap_extension, 5);
        
        NGSSimulator sampler(*xgidx,
                             fastq_name,
                             interleaved,
                             path_names,
                             transcript_expressions,
                             base_error,
                             indel_error,
                             indel_prop,
                             fragment_length ? fragment_length : std::numeric_limits<double>::max(), // suppresses warnings about fragment length
                             fragment_std_dev ? fragment_std_dev : 0.000001, // eliminates errors from having 0 as stddev without substantial difference
                             error_scale_factor,
                             !reads_may_contain_Ns,
                             seed_val);
        
        if (fragment_length) {
            for (size_t i = 0; i < num_reads; i++) {
                pair<Alignment, Alignment> read_pair = sampler.sample_read_pair();
                read_pair.first.set_score(aligner.score_ungapped_alignment(read_pair.first, strip_bonuses));
                read_pair.second.set_score(aligner.score_ungapped_alignment(read_pair.second, strip_bonuses));
                
                if (align_out) {
                    if (json_out) {
                        cout << pb2json(read_pair.first) << endl;
                        cout << pb2json(read_pair.second) << endl;
                    }
                    else {
                        aln_emitter->write_copy(read_pair.first);
                        aln_emitter->write_copy(read_pair.second);
                    }
                }
                else {
                    cout << read_pair.first.sequence() << "\t" << read_pair.second.sequence() << endl;
                }
            }
        }
        else {
            for (size_t i = 0; i < num_reads; i++) {
                Alignment read = sampler.sample_read();
                read.set_score(aligner.score_ungapped_alignment(read, strip_bonuses));
                
                if (align_out) {
                    if (json_out) {
                        cout << pb2json(read) << endl;
                    }
                    else {
                        aln_emitter->write_copy(read);
                    }
                }
                else {
                    cout << read.sequence() << endl;
                }
            }
        }
    }
    
    return 0;
}

// Register subcommand
static Subcommand vg_sim("sim", "simulate reads from a graph", TOOLKIT, main_sim);

