// construct.cpp: define the "vg construct" subcommand.

#include <omp.h>
#include <unistd.h>
#include <getopt.h>
#include <memory>

#include "subcommand.hpp"

#include "../stream.hpp"
#include "../constructor.hpp"
#include "../msa_converter.hpp"
#include "../region.hpp"

using namespace std;
using namespace vg;
using namespace vg::subcommand;

void help_construct(char** argv) {
    cerr << "usage: " << argv[0] << " construct [options] >new.vg" << endl
         << "options:" << endl
         << "construct from a reference and variant calls:" << endl
         << "    -r, --reference FILE   input FASTA reference (may repeat)" << endl
         << "    -v, --vcf FILE         input VCF (may repeat)" << endl
         << "    -n, --rename V=F       rename contig V in the VCFs to contig F in the FASTAs (may repeat)" << endl
         << "    -a, --alt-paths        save paths for alts of variants by variant ID" << endl
         << "    -R, --region REGION    specify a particular chromosome or 1-based inclusive region" << endl
         << "    -C, --region-is-chrom  don't attempt to parse the region (use when the reference" << endl
         << "                           sequence name could be inadvertently parsed as a region)" << endl
         << "    -z, --region-size N    variants per region to parallelize (default: 1024)" << endl
         << "    -t, --threads N        use N threads to construct graph (defaults to numCPUs)" << endl
         << "    -S, --handle-sv        include structural variants in construction of graph." << endl
         << "    -I, --insertions FILE  a FASTA file containing insertion sequences "<< endl
         << "                           (referred to in VCF) to add to graph." << endl
         << "    -f, --flat-alts N      don't chop up alternate alleles from input VCF" << endl
         << "construct from a multiple sequence alignment:" << endl
         << "    -M, --msa FILE         input multiple sequence alignment" << endl
         << "    -F, --msa-format       format of the MSA file (options: fasta, maf, clustal; default fasta)" << endl
         << "    -d, --drop-msa-paths   don't add paths for the MSA sequences into the graph" << endl
         << "shared construction options:" << endl
         << "    -m, --node-max N       limit the maximum allowable node sequence size (defaults to 1000)" << endl
         << "                           nodes greater than this threshold will be divided" << endl
         << "                           Note: nodes larger than ~1024 bp can't be GCSA2-indexed" << endl
         << "    -p, --progress         show progress" << endl;

}

int main_construct(int argc, char** argv) {

    if (argc == 2) {
        help_construct(argv);
        return 1;
    }

    // Make a constructor to fill in
    Constructor constructor;

    // We also parse some arguments separately.
    vector<string> fasta_filenames;
    vector<string> vcf_filenames;
    vector<string> insertion_filenames;
    string region;
    bool region_is_chrom = false;
    string msa_filename;
    int max_node_size = 1000;
    bool keep_paths = true;
    string msa_format = "fasta";
    bool show_progress = false;

    int c;
    optind = 2; // force optind past command positional argument
    while (true) {
        static struct option long_options[] =
            {
                /* These options set a flag. */
                //{"verbose", no_argument,       &verbose_flag, 1},
                {"vcf", required_argument, 0, 'v'},
                {"reference", required_argument, 0, 'r'},
                {"msa", required_argument, 0, 'M'},
                {"msa-format", required_argument, 0, 'F'},
                {"drop-msa-paths", no_argument, 0, 'd'},
                {"rename", required_argument, 0, 'n'},
                {"alt-paths", no_argument, 0, 'a'},
                {"handle-sv", no_argument, 0, 'S'},
                {"insertions", required_argument, 0, 'I'},
                {"progress",  no_argument, 0, 'p'},
                {"region-size", required_argument, 0, 'z'},
                {"threads", required_argument, 0, 't'},
                {"region", required_argument, 0, 'R'},
                {"region-is-chrom", no_argument, 0, 'C'},
                {"node-max", required_argument, 0, 'm'},\
                {"flat-alts", no_argument, 0, 'f'},
                {0, 0, 0, 0}
            };

        int option_index = 0;
        c = getopt_long (argc, argv, "v:r:n:ph?z:t:R:m:as:CfSI:M:dF:",
                         long_options, &option_index);

        /* Detect the end of the options. */
        if (c == -1)
            break;

        switch (c)
        {
        case 'v':
            vcf_filenames.push_back(optarg);
            break;

        case 'M':
            msa_filename = optarg;
            break;
            
        case 'F':
            msa_format = optarg;
            break;
            
        case 'd':
            keep_paths = false;
            break;

        case 'r':
            fasta_filenames.push_back(optarg);
            break;

        case 'S':
            constructor.do_svs = true;
            break;

        case 'I':
            insertion_filenames.push_back(optarg);
            break;

            
        case 'n':
            {
                // Parse the rename old=new
                string key_value(optarg);
                auto found = key_value.find('=');
                if (found == string::npos || found == 0 || found + 1 == key_value.size()) {
                    cerr << "error:[vg construct] could not parse rename " << key_value << endl;
                    exit(1);
                }
                // Parse out the two parts
                string vcf_contig = key_value.substr(0, found);
                string fasta_contig = key_value.substr(found + 1);
                // Add the name mapping
                constructor.add_name_mapping(vcf_contig, fasta_contig);
            }
            break;

        case 'a':
            constructor.alt_paths = true;
            break;

        case 'p':
            show_progress = true;
            break;

        case 'z':
            constructor.vars_per_chunk = atoi(optarg);
            break;

        case 'R':
            region = optarg;
            break;

        case 'C':
            region_is_chrom = true;
            break;

        case 't':
            omp_set_num_threads(atoi(optarg));
            break;

        case 'm':
            max_node_size = atoi(optarg);
            break;

        case 'f':
            constructor.flat = true;
            break;

        case 'h':
        case '?':
            /* getopt_long already printed an error message. */
            help_construct(argv);
            exit(1);
            break;

        default:
            abort ();

        }
    }
    
    // We need a callback to handle pieces of graph as they are produced.
    auto callback = [&](Graph& big_chunk) {
        // Wrap the chunk in a vg object that can properly divide it into
        // reasonably sized serialized chunks.
        VG* g = new VG(big_chunk, false, true);
#pragma omp critical (cout)
        g->serialize_to_ostream(cout);
    };
    
    constructor.max_node_size = max_node_size;
    constructor.show_progress = show_progress;
    
    if (constructor.max_node_size == 0) {
        // Make sure we can actually make nodes
        cerr << "error:[vg construct] max node size cannot be 0" << endl;
        exit(1);
    }
    
    if (!msa_filename.empty() && !fasta_filenames.empty()) {
        cerr << "error:[vg construct] cannot construct from a reference/VCF and an MSA simultaneously" << endl;
        exit(1);
    }
    
    if (!fasta_filenames.empty()) {
        
        
        if (!region.empty()) {
            // We want to limit to a certain region
            if (!region_is_chrom) {
                // We are allowed to parse the region.
                // Break out sequence name and region bounds
                string seq_name;
                int64_t start_pos = -1, stop_pos = -1;
                parse_region(region,
                             seq_name,
                             start_pos,
                             stop_pos);
                
                if (start_pos > 0 && stop_pos > 0) {
                    // These are 0-based, so if both are nonzero we got a real set of coordinates
                    if (constructor.show_progress) {
                        cerr << "Restricting to " << seq_name << " from " << start_pos << " to " << stop_pos << endl;
                    }
                    constructor.allowed_vcf_names.insert(seq_name);
                    // Make sure to correct the coordinates to 0-based exclusive-end, from 1-based inclusive-end
                    constructor.allowed_vcf_regions[seq_name] = make_pair(start_pos - 1, stop_pos);
                } else if (start_pos < 0 && stop_pos < 0) {
                    // We just got a name
                    cerr << "Restricting to " << seq_name << " from 1 to end" << endl;
                    constructor.allowed_vcf_names.insert(seq_name);
                } else {
                    // This doesn't make sense. Does it have like one coordinate?
                    cerr << "error:[vg construct] could not parse " << region << endl;
                    exit(1);
                }
            } else {
                // We have been told not to parse the region
                cerr << "Restricting to " << region << " from 1 to end" << endl;
                constructor.allowed_vcf_names.insert(region);
            }
        }
        
        // This will own all the VCF files
        vector<unique_ptr<vcflib::VariantCallFile>> variant_files;
        for (auto& vcf_filename : vcf_filenames) {
            // Make sure each VCF file exists. Otherwise Tabix++ may exit with a non-
            // helpful message.
            
            // We can't invoke stat woithout a place for it to write. But all we
            // really want is its return value.
            struct stat temp;
            if(stat(vcf_filename.c_str(), &temp)) {
                cerr << "error:[vg construct] file \"" << vcf_filename << "\" not found" << endl;
                return 1;
            }
            vcflib::VariantCallFile* variant_file = new vcflib::VariantCallFile();
            variant_files.emplace_back(variant_file);
            variant_file->open(vcf_filename);
            if (!variant_file->is_open()) {
                cerr << "error:[vg construct] could not open" << vcf_filename << endl;
                return 1;
            }
        }
        
        if (fasta_filenames.empty()) {
            cerr << "error:[vg construct] a reference is required for graph construction" << endl;
            return 1;
        }
        vector<unique_ptr<FastaReference>> references;
        for (auto& fasta_filename : fasta_filenames) {
            // Open each FASTA file
            FastaReference* reference = new FastaReference();
            references.emplace_back(reference);
            reference->open(fasta_filename);
        }
        
        vector<unique_ptr<FastaReference> > insertions;
        for (auto& insertion_filename : insertion_filenames){
            // Open up those insertion files
            FastaReference* insertion = new FastaReference();
            insertions.emplace_back(insertion);
            insertion->open(insertion_filename);
        }
        
        // Make vectors of just bare pointers
        vector<vcflib::VariantCallFile*> vcf_pointers;
        for(auto& vcf : variant_files) {
            vcf_pointers.push_back(vcf.get());
        }
        vector<FastaReference*> fasta_pointers;
        for(auto& fasta : references) {
            fasta_pointers.push_back(fasta.get());
        }
        vector<FastaReference*> ins_pointers;
        for (auto& ins : insertions){
            ins_pointers.push_back(ins.get());
        }
        
        if (ins_pointers.size() > 1){
            cerr << "Error: only one insertion file may be provided." << endl;
            exit(1);
        }
        
        // Construct the graph.
        constructor.construct_graph(fasta_pointers, vcf_pointers,
                                    ins_pointers, callback);
        
        // NB: If you worry about "still reachable but possibly lost" warnings in valgrind,
        // this would free all the memory used by protobuf:
        //ShutdownProtobufLibrary();
    }
    else if (!msa_filename.empty()) {
        
        ifstream msa_file(msa_filename);
        if (!msa_file) {
            cerr << "error:[vg construct] could not open MSA file " << msa_filename << endl;
            exit(1);
        }
        
        MSAConverter msa_converter;
        msa_converter.show_progress = show_progress;
        
        msa_converter.load_alignments(msa_file, msa_format);
        VG msa_graph = msa_converter.make_graph(keep_paths, max_node_size);
        
        callback(msa_graph.graph);
    }
    else {
        cerr << "error:[vg construct] a reference or an MSA is required for construct" << endl;
        exit(1);
    }

    return 0;
}

// Register subcommand
static Subcommand vg_construct("construct", "graph construction", PIPELINE, 1, main_construct);

