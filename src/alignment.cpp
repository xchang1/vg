#include "alignment.hpp"
#include "gafkluge.hpp"

#include <sstream>
#include <regex>

namespace vg {

int hts_for_each(string& filename, function<void(Alignment&)> lambda, const PathPositionHandleGraph* graph) {

    samFile *in = hts_open(filename.c_str(), "r");
    if (in == NULL) return 0;
    bam_hdr_t *hdr = sam_hdr_read(in);
    map<string, string> rg_sample;
    parse_rg_sample_map(hdr->text, rg_sample);
    bam1_t *b = bam_init1();
    while (sam_read1(in, hdr, b) >= 0) {
        Alignment a = bam_to_alignment(b, rg_sample, hdr, graph);
        lambda(a);
    }
    bam_destroy1(b);
    bam_hdr_destroy(hdr);
    hts_close(in);
    return 1;

}

int hts_for_each(string& filename, function<void(Alignment&)> lambda) {
    return hts_for_each(filename, lambda, nullptr);
}

int hts_for_each_parallel(string& filename, function<void(Alignment&)> lambda,
                          const PathPositionHandleGraph* graph) {

    samFile *in = hts_open(filename.c_str(), "r");
    if (in == NULL) return 0;
    bam_hdr_t *hdr = sam_hdr_read(in);
    map<string, string> rg_sample;
    parse_rg_sample_map(hdr->text, rg_sample);

    int thread_count = get_thread_count();
    vector<bam1_t*> bs; bs.resize(thread_count);
    for (auto& b : bs) {
        b = bam_init1();
    }

    bool more_data = true;
#pragma omp parallel shared(in, hdr, more_data, rg_sample)
    {
        int tid = omp_get_thread_num();
        while (more_data) {
            bam1_t* b = bs[tid];
            // We need to track our own read operation's success separate from
            // the global flag, or someone else encountering EOF will cause us
            // to drop our read on the floor.
            bool got_read = false;
#pragma omp critical (hts_input)
            if (more_data) {
                got_read = sam_read1(in, hdr, b) >= 0;
                more_data &= got_read;
            }
            // Now we're outside the critical section so we can only rely on our own variables.
            if (got_read) {
                Alignment a = bam_to_alignment(b, rg_sample, hdr, graph);
                lambda(a);
            }
        }
    }

    for (auto& b : bs) bam_destroy1(b);
    bam_hdr_destroy(hdr);
    hts_close(in);
    return 1;

}

int hts_for_each_parallel(string& filename, function<void(Alignment&)> lambda) {
    return hts_for_each_parallel(filename, lambda, nullptr);
}

bam_hdr_t* hts_file_header(string& filename, string& header) {
    samFile *in = hts_open(filename.c_str(), "r");
    if (in == NULL) {
        cerr << "[vg::alignment] could not open " << filename << endl;
        exit(1);
    }
    bam_hdr_t *hdr = sam_hdr_read(in);
    header = hdr->text;
    bam_hdr_destroy(hdr);
    hts_close(in);
    return hdr;
}

bam_hdr_t* hts_string_header(string& header,
                             map<string, int64_t>& path_length,
                             map<string, string>& rg_sample) {
    stringstream hdr;
    hdr << "@HD\tVN:1.5\tSO:unknown\n";
    for (auto& p : path_length) {
        hdr << "@SQ\tSN:" << p.first << "\t" << "LN:" << p.second << "\n";
    }
    for (auto& s : rg_sample) {
        hdr << "@RG\tID:" << s.first << "\t" << "SM:" << s.second << "\n";
    }
    hdr << "@PG\tID:0\tPN:vg\n";
    header = hdr.str();
    string sam = "data:," + header;
    samFile *in = sam_open(sam.c_str(), "r");
    bam_hdr_t *h = sam_hdr_read(in);
    sam_close(in);
    return h;
}

bool get_next_alignment_from_fastq(gzFile fp, char* buffer, size_t len, Alignment& alignment) {

    alignment.Clear();
    bool is_fasta = false;
    // handle name
    if (0!=gzgets(fp,buffer,len)) {
        buffer[strlen(buffer)-1] = '\0';
        string name = buffer;
        if (name[0] == '@') {
            is_fasta = false;
        } else if (name[0] == '>') {
            is_fasta = true;
        } else {
            throw runtime_error("Found unexpected delimiter " + name.substr(0,1) + " in fastq/fasta input");
        }
        name = name.substr(1, name.find(' ') - 1); // trim off leading @ and things after the first whitespace
        // keep trailing /1 /2
        alignment.set_name(name);
    } else { return false; }
    // handle sequence
    if (0!=gzgets(fp,buffer,len)) {
        buffer[strlen(buffer)-1] = '\0';
        alignment.set_sequence(buffer);
    } else {
        cerr << "[vg::alignment.cpp] error: incomplete fastq record" << endl; exit(1);
    }
    // handle "+" sep
    if (!is_fasta) {
        if (0!=gzgets(fp,buffer,len)) {
        } else {
            cerr << "[vg::alignment.cpp] error: incomplete fastq record" << endl; exit(1);
        }
        // handle quality
        if (0!=gzgets(fp,buffer,len)) {
            buffer[strlen(buffer)-1] = '\0';
            string quality = string_quality_char_to_short(buffer);
            //cerr << string_quality_short_to_char(quality) << endl;
            alignment.set_quality(quality);
        } else {
            cerr << "[vg::alignment.cpp] error: incomplete fastq record" << endl; exit(1);
        }
    }

    return true;

}

bool get_next_interleaved_alignment_pair_from_fastq(gzFile fp, char* buffer, size_t len, Alignment& mate1, Alignment& mate2) {
    return get_next_alignment_from_fastq(fp, buffer, len, mate1) && get_next_alignment_from_fastq(fp, buffer, len, mate2);
}

bool get_next_alignment_pair_from_fastqs(gzFile fp1, gzFile fp2, char* buffer, size_t len, Alignment& mate1, Alignment& mate2) {
    return get_next_alignment_from_fastq(fp1, buffer, len, mate1) && get_next_alignment_from_fastq(fp2, buffer, len, mate2);
}

size_t unpaired_for_each_parallel(function<bool(Alignment&)> get_read_if_available,
                                  function<void(Alignment&)> lambda,
                                  uint64_t batch_size = DEFAULT_PARALLEL_BATCHSIZE) {
    assert(batch_size % 2 == 0);    
    size_t nLines = 0;
    vector<Alignment> *batch = nullptr;
    // number of batches currently being processed
    uint64_t batches_outstanding = 0;
#pragma omp parallel default(none) shared(batches_outstanding, batch, nLines, get_read_if_available, lambda, batch_size)
#pragma omp single
    {
        
        // max # of such batches to be holding in memory
        uint64_t max_batches_outstanding = batch_size;
        // max # we will ever increase the batch buffer to
        const uint64_t max_max_batches_outstanding = 1 << 13; // 8192
        
        // alignments to hold the incoming data
        Alignment aln;
        // did we find the end of the file yet?
        bool more_data = true;
        
        while (more_data) {
            // init a new batch
            batch = new std::vector<Alignment>();
            batch->reserve(batch_size);
            
            // load up to the batch-size number of reads
            for (int i = 0; i < batch_size; i++) {
                
                more_data = get_read_if_available(aln);
                
                if (more_data) {
                    batch->emplace_back(std::move(aln));
                    nLines++;
                }
                else {
                    break;
                }
            }
            
            // did we get a batch?
            if (batch->size()) {
                
                // how many batch tasks are outstanding currently, including this one?
                uint64_t current_batches_outstanding;
#pragma omp atomic capture
                current_batches_outstanding = ++batches_outstanding;
                
                if (current_batches_outstanding >= max_batches_outstanding) {
                    // do this batch in the current thread because we've spawned the maximum number of
                    // concurrent batch tasks
                    for (auto& aln : *batch) {
                        lambda(aln);
                    }
                    delete batch;
#pragma omp atomic capture
                    current_batches_outstanding = --batches_outstanding;
                    
                    if (4 * current_batches_outstanding / 3 < max_batches_outstanding
                        && max_batches_outstanding < max_max_batches_outstanding) {
                        // we went through at least 1/4 of the batch buffer while we were doing this thread's batch
                        // this looks risky, since we want the batch buffer to stay populated the entire time we're
                        // occupying this thread on compute, so let's increase the batch buffer size
                        
                        max_batches_outstanding *= 2;
                    }
                }
                else {
                    // spawn a new task to take care of this batch
#pragma omp task default(none) firstprivate(batch) shared(batches_outstanding, lambda)
                    {
                        for (auto& aln : *batch) {
                            lambda(aln);
                        }
                        delete batch;
#pragma omp atomic update
                        batches_outstanding--;
                    }
                }
            }
        }
    }
    return nLines;
}

size_t paired_for_each_parallel_after_wait(function<bool(Alignment&, Alignment&)> get_pair_if_available,
                                           function<void(Alignment&, Alignment&)> lambda,
                                           function<bool(void)> single_threaded_until_true,
                                           uint64_t batch_size = DEFAULT_PARALLEL_BATCHSIZE) {

    assert(batch_size % 2 == 0);
    size_t nLines = 0;
    vector<pair<Alignment, Alignment> > *batch = nullptr;
    // number of batches currently being processed
    uint64_t batches_outstanding = 0;
    
#pragma omp parallel default(none) shared(batches_outstanding, batch, nLines, get_pair_if_available, single_threaded_until_true, lambda, batch_size)
#pragma omp single
    {

        // max # of such batches to be holding in memory
        uint64_t max_batches_outstanding = batch_size;
        // max # we will ever increase the batch buffer to
        const uint64_t max_max_batches_outstanding = 1 << 13; // 8192
        
        // alignments to hold the incoming data
        Alignment mate1, mate2;
        // did we find the end of the file yet?
        bool more_data = true;
        
        while (more_data) {
            // init a new batch
            batch = new std::vector<pair<Alignment, Alignment>>();
            batch->reserve(batch_size);
            
            // load up to the batch-size number of pairs
            for (int i = 0; i < batch_size; i++) {
                
                more_data = get_pair_if_available(mate1, mate2);
                
                if (more_data) {
                    batch->emplace_back(std::move(mate1), std::move(mate2));
                    nLines++;
                }
                else {
                    break;
                }
            }
            
            // did we get a batch?
            if (batch->size()) {
                // how many batch tasks are outstanding currently, including this one?
                uint64_t current_batches_outstanding;
#pragma omp atomic capture
                current_batches_outstanding = ++batches_outstanding;
                
                bool do_single_threaded = !single_threaded_until_true();
                if (current_batches_outstanding >= max_batches_outstanding || do_single_threaded) {
                    // do this batch in the current thread because we've spawned the maximum number of
                    // concurrent batch tasks or because we are directed to work in a single thread
                    for (auto& p : *batch) {
                        lambda(p.first, p.second);
                    }
                    delete batch;
#pragma omp atomic capture
                    current_batches_outstanding = --batches_outstanding;
                    
                    if (4 * current_batches_outstanding / 3 < max_batches_outstanding
                        && max_batches_outstanding < max_max_batches_outstanding
                        && !do_single_threaded) {
                        // we went through at least 1/4 of the batch buffer while we were doing this thread's batch
                        // this looks risky, since we want the batch buffer to stay populated the entire time we're
                        // occupying this thread on compute, so let's increase the batch buffer size
                        // (skip this adjustment if you're in single-threaded mode and thus expect the buffer to be
                        // empty)
                        
                        max_batches_outstanding *= 2;
                    }
                }
                else {
                    // spawn a new task to take care of this batch
#pragma omp task default(none) firstprivate(batch) shared(batches_outstanding, lambda)
                    {
                        for (auto& p : *batch) {
                            lambda(p.first, p.second);
                        }
                        delete batch;
#pragma omp atomic update
                        batches_outstanding--;
                    }
                }
            }
        }
    }
    
    return nLines;
}

size_t fastq_unpaired_for_each_parallel(const string& filename, function<void(Alignment&)> lambda) {
    
    gzFile fp = (filename != "-") ? gzopen(filename.c_str(), "r") : gzdopen(fileno(stdin), "r");
    if (!fp) {
        cerr << "[vg::alignment.cpp] couldn't open " << filename << endl; exit(1);
    }
    
    size_t len = 2 << 22; // 4M
    char* buf = new char[len];
    
    function<bool(Alignment&)> get_read = [&](Alignment& aln) {
        return get_next_alignment_from_fastq(fp, buf, len, aln);;
    };
    
    
    size_t nLines = unpaired_for_each_parallel(get_read, lambda);
    
    delete[] buf;
    gzclose(fp);
    return nLines;
    
}

size_t fastq_paired_interleaved_for_each_parallel(const string& filename, function<void(Alignment&, Alignment&)> lambda) {
    return fastq_paired_interleaved_for_each_parallel_after_wait(filename, lambda, [](void) {return true;});
}
    
size_t fastq_paired_two_files_for_each_parallel(const string& file1, const string& file2, function<void(Alignment&, Alignment&)> lambda) {
    return fastq_paired_two_files_for_each_parallel_after_wait(file1, file2, lambda, [](void) {return true;});
}
    
size_t fastq_paired_interleaved_for_each_parallel_after_wait(const string& filename,
                                                             function<void(Alignment&, Alignment&)> lambda,
                                                             function<bool(void)> single_threaded_until_true) {
    
    gzFile fp = (filename != "-") ? gzopen(filename.c_str(), "r") : gzdopen(fileno(stdin), "r");
    if (!fp) {
        cerr << "[vg::alignment.cpp] couldn't open " << filename << endl; exit(1);
    }
    
    size_t len = 1 << 18; // 256k
    char* buf = new char[len];
    
    function<bool(Alignment&, Alignment&)> get_pair = [&](Alignment& mate1, Alignment& mate2) {
        return get_next_interleaved_alignment_pair_from_fastq(fp, buf, len, mate1, mate2);
    };
    
    size_t nLines = paired_for_each_parallel_after_wait(get_pair, lambda, single_threaded_until_true);
    
    delete[] buf;
    gzclose(fp);
    return nLines;
}
    
size_t fastq_paired_two_files_for_each_parallel_after_wait(const string& file1, const string& file2,
                                                           function<void(Alignment&, Alignment&)> lambda,
                                                           function<bool(void)> single_threaded_until_true) {
    
    gzFile fp1 = (file1 != "-") ? gzopen(file1.c_str(), "r") : gzdopen(fileno(stdin), "r");
    if (!fp1) {
        cerr << "[vg::alignment.cpp] couldn't open " << file1 << endl; exit(1);
    }
    gzFile fp2 = (file2 != "-") ? gzopen(file2.c_str(), "r") : gzdopen(fileno(stdin), "r");
    if (!fp2) {
        cerr << "[vg::alignment.cpp] couldn't open " << file2 << endl; exit(1);
    }
    
    size_t len = 1 << 18; // 256k
    char* buf = new char[len];
    
    function<bool(Alignment&, Alignment&)> get_pair = [&](Alignment& mate1, Alignment& mate2) {
        return get_next_alignment_pair_from_fastqs(fp1, fp2, buf, len, mate1, mate2);
    };
    
    size_t nLines = paired_for_each_parallel_after_wait(get_pair, lambda, single_threaded_until_true);
    
    delete[] buf;
    gzclose(fp1);
    gzclose(fp2);
    return nLines;
}

size_t fastq_unpaired_for_each(const string& filename, function<void(Alignment&)> lambda) {
    gzFile fp = (filename != "-") ? gzopen(filename.c_str(), "r") : gzdopen(fileno(stdin), "r");
    if (!fp) {
        cerr << "[vg::alignment.cpp] couldn't open " << filename << endl; exit(1);
    }
    size_t len = 2 << 22; // 4M
    size_t nLines = 0;
    char *buffer = new char[len];
    Alignment alignment;
    while(get_next_alignment_from_fastq(fp, buffer, len, alignment)) {
        lambda(alignment);
        nLines++;
    }
    gzclose(fp);
    delete[] buffer;
    return nLines;
}

size_t fastq_paired_interleaved_for_each(const string& filename, function<void(Alignment&, Alignment&)> lambda) {
    gzFile fp = (filename != "-") ? gzopen(filename.c_str(), "r") : gzdopen(fileno(stdin), "r");
    if (!fp) {
        cerr << "[vg::alignment.cpp] couldn't open " << filename << endl; exit(1);
    }
    size_t len = 2 << 18; // 256k
    size_t nLines = 0;
    char *buffer = new char[len];
    Alignment mate1, mate2;
    while(get_next_interleaved_alignment_pair_from_fastq(fp, buffer, len, mate1, mate2)) {
        lambda(mate1, mate2);
        nLines++;
    }
    gzclose(fp);
    delete[] buffer;
    return nLines;
}


size_t fastq_paired_two_files_for_each(const string& file1, const string& file2, function<void(Alignment&, Alignment&)> lambda) {
    gzFile fp1 = (file1 != "-") ? gzopen(file1.c_str(), "r") : gzdopen(fileno(stdin), "r");
    if (!fp1) {
        cerr << "[vg::alignment.cpp] couldn't open " << file1 << endl; exit(1);
    }
    gzFile fp2 = (file2 != "-") ? gzopen(file2.c_str(), "r") : gzdopen(fileno(stdin), "r");
    if (!fp2) {
        cerr << "[vg::alignment.cpp] couldn't open " << file2 << endl; exit(1);
    }
    size_t len = 2 << 18; // 256k
    size_t nLines = 0;
    char *buffer = new char[len];
    Alignment mate1, mate2;
    while(get_next_alignment_pair_from_fastqs(fp1, fp2, buffer, len, mate1, mate2)) {
        lambda(mate1, mate2);
        nLines++;
    }
    gzclose(fp1);
    gzclose(fp2);
    delete[] buffer;
    return nLines;

}

bool get_next_alignment_from_gaf(const HandleGraph& graph, htsFile* fp, kstring_t& s_buffer, gafkluge::GafRecord& g_buffer,
                                 Alignment& alignment) {
    if (hts_getline(fp, '\n', &s_buffer) <= 0) {
        return false;
    }

    gafkluge::parse_gaf_record(ks_str(&s_buffer), g_buffer);
    gaf_to_alignment(graph, g_buffer, alignment);
    return true;
}

bool get_next_interleaved_alignment_pair_from_gaf(const HandleGraph& graph, htsFile* fp, kstring_t& s_buffer,
                                                  gafkluge::GafRecord& g_buffer, Alignment& mate1, Alignment& mate2) {
    return get_next_alignment_from_gaf(graph, fp, s_buffer, g_buffer, mate1) &&
        get_next_alignment_from_gaf(graph, fp, s_buffer, g_buffer, mate2);
}

size_t gaf_unpaired_for_each(const HandleGraph& graph, const string& filename, function<void(Alignment&)> lambda) {

    htsFile* in = hts_open(filename.c_str(), "r");
    if (in == NULL) {
        cerr << "[vg::alignment.cpp] couldn't open " << filename << endl; exit(1);
    }
    
    kstring_t s_buffer = KS_INITIALIZE;
    Alignment aln;
    gafkluge::GafRecord gaf;
    size_t count = 0;

    while (get_next_alignment_from_gaf(graph, in, s_buffer, gaf, aln) == true) {
        lambda(aln);
        ++count;
    }
    
    hts_close(in);

    return count;
}

size_t gaf_paired_interleaved_for_each(const HandleGraph& graph, const string& filename,
                                       function<void(Alignment&, Alignment&)> lambda) {

    htsFile* in = hts_open(filename.c_str(), "r");
    if (in == NULL) {
        cerr << "[vg::alignment.cpp] couldn't open " << filename << endl; exit(1);
    }
    
    kstring_t s_buffer = KS_INITIALIZE;
    Alignment aln1, aln2;
    gafkluge::GafRecord gaf;
    size_t count = 0;

    while (get_next_interleaved_alignment_pair_from_gaf(graph, in, s_buffer, gaf, aln1, aln2) == true) {
        lambda(aln1, aln2);
        count += 2;
    }
    
    hts_close(in);

    return count;
}

size_t gaf_unpaired_for_each_parallel(const HandleGraph& graph, const string& filename,
                                      function<void(Alignment&)> lambda,
                                      uint64_t batch_size) {

    htsFile* in = hts_open(filename.c_str(), "r");
    if (in == NULL) {
        cerr << "[vg::alignment.cpp] couldn't open " << filename << endl; exit(1);
    }

    kstring_t s_buffer = KS_INITIALIZE;
    Alignment aln1, aln2;
    gafkluge::GafRecord gaf;
    
    function<bool(Alignment&)> get_read = [&](Alignment& aln) {
        return get_next_alignment_from_gaf(graph, in, s_buffer, gaf, aln);
    };
        
    size_t nLines = unpaired_for_each_parallel(get_read, lambda, batch_size);
    
    hts_close(in);
    return nLines;

}

size_t gaf_paired_interleaved_for_each_parallel(const HandleGraph& graph, const string& filename,
                                                function<void(Alignment&, Alignment&)> lambda,
                                                uint64_t batch_size) {
    return gaf_paired_interleaved_for_each_parallel_after_wait(graph, filename, lambda, [](void) {return true;}, batch_size);
}

size_t gaf_paired_interleaved_for_each_parallel_after_wait(const HandleGraph& graph, const string& filename,
                                                           function<void(Alignment&, Alignment&)> lambda,
                                                           function<bool(void)> single_threaded_until_true,
                                                           uint64_t batch_size) {
    
    htsFile* in = hts_open(filename.c_str(), "r");
    if (in == NULL) {
        cerr << "[vg::alignment.cpp] couldn't open " << filename << endl; exit(1);
    }

    kstring_t s_buffer = KS_INITIALIZE;
    Alignment aln1, aln2;
    gafkluge::GafRecord gaf;
    
    function<bool(Alignment&, Alignment&)> get_pair = [&](Alignment& mate1, Alignment& mate2) {
        return get_next_interleaved_alignment_pair_from_gaf(graph, in, s_buffer, gaf, mate1, mate2);
    };
    
    size_t nLines = paired_for_each_parallel_after_wait(get_pair, lambda, single_threaded_until_true, batch_size);

    hts_close(in);
    return nLines;    
}

gafkluge::GafRecord alignment_to_gaf(const HandleGraph& graph, const Alignment& aln, bool cs_cigar, bool base_quals) {

    gafkluge::GafRecord gaf;

    //1 string Query sequence name
    gaf.query_name = aln.name();

    //2 int Query sequence length
    gaf.query_length = aln.sequence().length();

    //12 int Mapping quality (0-255; 255 for missing)
    //Note: protobuf can't distinguish between 0 and missing so we just copy it through
    gaf.mapq = aln.mapping_quality();

    if (aln.has_path() && aln.path().mapping_size() > 0) {    
        //3 int Query start (0-based; closed)
        gaf.query_start = 0; //(aln.path().mapping_size() ? first_path_position(aln.path()).offset() : "*") << "\t"
        //4 int Query end (0-based; open)
        gaf.query_end = aln.sequence().length();
        //5 char Strand relative to the path: "+" or "-"
        gaf.strand = '+'; // always positive relative to the path
        //7 int Path length
        gaf.path_length = 0;
        //8 int Start position on the path (0-based)
        gaf.path_start = 0;
        //9 int End position on the path (0-based)
        uint64_t end_position = 0;
        //10 int Number of residue matches
        gaf.matches = 0;
        gaf.path.reserve(aln.path().mapping_size());
        string cs_cigar_str;
        size_t running_match_length = 0;
        size_t total_to_len = 0;
        for (size_t i = 0; i < aln.path().mapping_size(); ++i) {
            auto& mapping = aln.path().mapping(i);
            size_t offset = mapping.position().offset();
            string node_seq;
            handle_t handle = graph.get_handle(mapping.position().node_id(), mapping.position().is_reverse());
            if (i == 0) {
                // use path_start to store the offset of the first node
                gaf.path_start = offset;
            } else if (cs_cigar == true && offset > 0) {
                // to support split-mappings we gobble up the beginnings
                // of nodes using deletions since unlike GAM, we can only
                // set the offset of the first node
                if (node_seq.empty()) {
                    node_seq = graph.get_sequence(handle);
                }
                cs_cigar_str += "-" + node_seq.substr(0, offset);
            }
            for (size_t j = 0; j < mapping.edit_size(); ++j) {
                auto& edit = mapping.edit(j);
                if (edit_is_match(edit)) {
                    gaf.matches += edit.from_length();
                }
                if (cs_cigar == true) {
                    // CS-cigar string
                    if (edit_is_match(edit)) {
                        // Merge up matches that span edits/mappings
                        running_match_length += edit.from_length();
                    } else {
                        if (running_match_length > 0) {
                            // Matches are : followed by the match length
                            cs_cigar_str += ":" + std::to_string(running_match_length);
                            running_match_length = 0;
                        }
                        if (edit_is_sub(edit)) {
                            if (node_seq.empty()) {
                                node_seq = graph.get_sequence(handle);
                            }
                            // Substitions expressed one base at a time, preceded by *
                            for (size_t k = 0; k < edit.from_length(); ++k) {
                                cs_cigar_str += "*" + node_seq.substr(offset + k, 1) + edit.sequence().substr(k, 1); 
                            }
                        } else if (edit_is_deletion(edit)) {
                            if (node_seq.empty()) {
                                node_seq = graph.get_sequence(handle);
                            }
                            // Deletion is - followed by deleted sequence
                            assert(offset + edit.from_length() <= node_seq.length());
                            cs_cigar_str += "-" + node_seq.substr(offset, edit.from_length());
                        } else if (edit_is_insertion(edit)) {
                            // Insertion is "+" followed by inserted sequence
                            cs_cigar_str += "+" + edit.sequence();
                        }
                    }
                }
                offset += edit.from_length();
                total_to_len += edit.to_length();
            }

            bool skip_step = false;
            if (i < aln.path().mapping_size() - 1 && offset != graph.get_length(handle)) {
                if (mapping.position().node_id() != aln.path().mapping(i + 1).position().node_id() ||
                    mapping.position().is_reverse() != aln.path().mapping(i + 1).position().is_reverse()) {
                    // we are hopping off the middle of a node, need to gobble it up with a deletion
                    if (node_seq.empty()) {
                        node_seq = graph.get_sequence(handle);
                    }
                    if (running_match_length > 0) {
                        // Matches are : followed by the match length
                        cs_cigar_str += ":" + std::to_string(running_match_length);
                        running_match_length = 0;
                    }
                    cs_cigar_str += "-" + node_seq.substr(offset);
                } else {
                    // we have a duplicate node mapping.  vg map actually produces these sometimes
                    // where an insert gets its own mapping even though its from_length is 0
                    // the gaf cigar format assumes nodes are fully covered, so we squish it out.
                    skip_step = true;
                }
            }
            
            //6 string Path matching /([><][^\s><]+(:\d+-\d+)?)+|([^\s><]+)/
            if (!skip_step) {
                auto& position = mapping.position();
                gafkluge::GafStep step;
                step.name = std::to_string(position.node_id());
                step.is_stable = false;
                step.is_reverse = position.is_reverse();
                step.is_interval = false;
                uint64_t node_length = graph.get_length(graph.get_handle(position.node_id()));
                gaf.path_length += node_length;
                if (i == 0) {
                    gaf.path_start = position.offset();
                }
                if (i == aln.path().mapping_size()-1) {
                    gaf.path_end = node_length - position.offset() - mapping_from_length(aln.path().mapping(i));
                }
                gaf.path.push_back(std::move(step));
            }
        }
        if (cs_cigar && running_match_length > 0) {
            cs_cigar_str += ":" + std::to_string(running_match_length);
            running_match_length = 0;
        }

        // We can support gam alignments without sequences by inferring the sequence length from edits
        if (gaf.query_length == 0 && total_to_len > 0) {
            gaf.query_length = total_to_len;
            gaf.query_end = total_to_len;
        } 

        //11 int Alignment block length
        gaf.block_length = std::max(gaf.path_end - gaf.path_start, gaf.query_length);

        // optional cs-cigar string
        if (cs_cigar) {
            gaf.opt_fields["cs"] = make_pair("Z", std::move(cs_cigar_str));
        }

        // convert the identity into the dv divergence field
        // https://lh3.github.io/minimap2/minimap2.html#10
        if (aln.identity() > 0) {
            stringstream dv_str;
            dv_str << std::floor((1. - aln.identity()) * 10000. + 0.5) / 10000.;
            gaf.opt_fields["dv"] = make_pair("f", dv_str.str());
        }

        // convert the score into the AS field
        // https://lh3.github.io/minimap2/minimap2.html#10
        if (aln.score() > 0) {
            gaf.opt_fields["AS"] = make_pair("i", std::to_string(aln.score()));
        }

        // optional base qualities
        if (base_quals) { 
            gaf.opt_fields["bq"] = make_pair("Z", string_quality_short_to_char(aln.quality()));
        }   
                
    }

    return gaf;
    
}

void gaf_to_alignment(const HandleGraph& graph, const gafkluge::GafRecord& gaf, Alignment& aln) {

    aln.Clear();

    aln.set_name(gaf.query_name);

    for (size_t i = 0; i < gaf.path.size(); ++i) {
        const auto& gaf_step = gaf.path[i];
        // only support unstable gaf at this point
        assert(gaf_step.is_stable == false);
        assert(gaf_step.is_interval == false);
        Mapping* mapping = aln.mutable_path()->add_mapping();
        mapping->mutable_position()->set_node_id(std::stol(gaf_step.name));
        mapping->mutable_position()->set_is_reverse(gaf_step.is_reverse);
        if (i == 0) {
            mapping->mutable_position()->set_offset(gaf.path_start);
        }
        mapping->set_rank(i + 1);
    }

    if (gaf.mapq != 255) {
        // We let 255 be equivalent to 0, which isn't great
        aln.set_mapping_quality(gaf.mapq);
    }

    if (!gaf.path.empty()) {
        size_t cur_mapping = 0;
        int64_t cur_offset = gaf.path_start;
        handle_t cur_handle = graph.get_handle(aln.path().mapping(cur_mapping).position().node_id(),
                                               aln.path().mapping(cur_mapping).position().is_reverse());
        size_t cur_len = graph.get_length(cur_handle);
        string& sequence = *aln.mutable_sequence();
        // Use the CS cigar string to add Edits into our Path, as well as set the sequence
        gafkluge::for_each_cs(gaf, [&] (const string& cs_cigar) {
                assert(cur_offset < cur_len);

                if (cs_cigar[0] == ':') {
                    int64_t match_len = stol(cs_cigar.substr(1));
                    while (match_len > 0) {
                        int64_t current_match = std::min(match_len, (int64_t)graph.get_length(cur_handle) - cur_offset);
                        Edit* edit = aln.mutable_path()->mutable_mapping(cur_mapping)->add_edit();
                        edit->set_from_length(current_match);
                        edit->set_to_length(current_match);
                        sequence += graph.get_sequence(cur_handle).substr(cur_offset, current_match);
                        match_len -= current_match;
                        cur_offset += current_match;
                        if (match_len > 0) {
                            assert(cur_mapping < aln.path().mapping_size() - 1);
                            ++cur_mapping;
                            cur_offset = 0;
                            cur_handle = graph.get_handle(aln.path().mapping(cur_mapping).position().node_id(),
                                                          aln.path().mapping(cur_mapping).position().is_reverse());
                            cur_len = graph.get_length(cur_handle);
                        }
                    }
                } else if (cs_cigar[0] == '+') {
                    size_t tgt_mapping = cur_mapping;
                    // left-align insertions to try to be more consistent with vg
                    if (cur_offset == 0 && cur_mapping > 0 && (!aln.path().mapping(cur_mapping - 1).position().is_reverse()
                                                               || cur_mapping == aln.path().mapping_size())) {
                        --tgt_mapping;
                    }
                    Edit* edit = aln.mutable_path()->mutable_mapping(tgt_mapping)->add_edit();
                    edit->set_from_length(0);
                    edit->set_to_length(cs_cigar.length() - 1);
                    edit->set_sequence(cs_cigar.substr(1));
                    sequence += edit->sequence();
                } else if (cs_cigar[0] == '-') {
                    string del = cs_cigar.substr(1);
                    assert(del.length() <= graph.get_length(cur_handle) - cur_offset);
                    assert(del == graph.get_sequence(cur_handle).substr(cur_offset, del.length()));
                    Edit* edit = aln.mutable_path()->mutable_mapping(cur_mapping)->add_edit();
                    edit->set_to_length(0);
                    edit->set_from_length(del.length());
                    cur_offset += del.length();
                    // unlike matches, we don't allow deletions to span multiple nodes
                    assert(cur_offset <= graph.get_length(cur_handle));
                } else if (cs_cigar[0] == '*') {
                    assert(cs_cigar.length() == 3);
                    char from = cs_cigar[1];
                    char to = cs_cigar[2];
                    assert(graph.get_sequence(cur_handle)[cur_offset] == from);
                    Edit* edit = aln.mutable_path()->mutable_mapping(cur_mapping)->add_edit();
                    // todo: support multibase snps
                    edit->set_from_length(1);
                    edit->set_to_length(1);
                    edit->set_sequence(string(1, to));
                    sequence += edit->sequence();
                    ++cur_offset;
                }
            
                // advance to the next mapping if we've pushed the offset past the current node
                assert(cur_offset <= cur_len);
                if (cur_offset == cur_len) {
                    ++cur_mapping;
                    cur_offset = 0;
                    if (cur_mapping < aln.path().mapping_size()) {
                        cur_handle = graph.get_handle(aln.path().mapping(cur_mapping).position().node_id(),
                                                      aln.path().mapping(cur_mapping).position().is_reverse());
                        cur_len = graph.get_length(cur_handle);
                    }
                }
            });
    }

    for (auto opt_it : gaf.opt_fields) {
        if (opt_it.first == "dv") {
            // get the identity from the dv divergence field
            // https://lh3.github.io/minimap2/minimap2.html#10
            aln.set_identity(1. - std::stof(opt_it.second.second));
        } else if (opt_it.first == "AS") {
            // get the score from the AS field
            // https://lh3.github.io/minimap2/minimap2.html#10
            aln.set_score(std::stoi(opt_it.second.second));
        } else if (opt_it.first == "bq") {
            // get the quality from the bq field
            aln.set_quality(string_quality_char_to_short(opt_it.second.second));
        }
    }
}

void parse_rg_sample_map(char* hts_header, map<string, string>& rg_sample) {
    string header(hts_header);
    vector<string> header_lines = split_delims(header, "\n");

    for (auto& line : header_lines) {

        // get next line from header, skip if empty
        if ( line.empty() ) { continue; }

        // lines of the header look like:
        // "@RG     ID:-    SM:NA11832      CN:BCM  PL:454"
        //                     ^^^^^^^\ is our sample name
        if (line.find("@RG") == 0) {
            vector<string> rg_parts = split_delims(line, "\t ");
            string name;
            string rg_id;
            for (auto& part : rg_parts) {
                size_t colpos = part.find(":");
                if (colpos != string::npos) {
                    string fieldname = part.substr(0, colpos);
                    if (fieldname == "SM") {
                        name = part.substr(colpos+1);
                    } else if (fieldname == "ID") {
                        rg_id = part.substr(colpos+1);
                    }
                }
            }
            if (name.empty()) {
                cerr << "[vg::alignment] Error: could not find 'SM' in @RG line " << endl << line << endl;
                exit(1);
            }
            if (rg_id.empty()) {
                cerr << "[vg::alignment] Error: could not find 'ID' in @RG line " << endl << line << endl;
                exit(1);
            }
            map<string, string>::iterator s = rg_sample.find(rg_id);
            if (s != rg_sample.end()) {
                if (s->second != name) {
                    cerr << "[vg::alignment] Error: multiple samples (SM) map to the same read group (RG)" << endl
                          << endl
                          << "samples " << name << " and " << s->second << " map to " << rg_id << endl
                          << endl
                          << "It will not be possible to determine what sample an alignment belongs to" << endl
                          << "at runtime." << endl
                          << endl
                          << "To resolve the issue, ensure that RG ids are unique to one sample" << endl
                          << "across all the input files to freebayes." << endl
                          << endl
                          << "See bamaddrg (https://github.com/ekg/bamaddrg) for a method which can" << endl
                          << "add RG tags to alignments." << endl;
                    exit(1);
                }
            }
            // if it's the same sample name and RG combo, no worries
            rg_sample[rg_id] = name;
        }
    }
}

short quality_char_to_short(char c) {
    return static_cast<short>(c) - 33;
}

char quality_short_to_char(short i) {
    return static_cast<char>(i + 33);
}

void alignment_quality_short_to_char(Alignment& alignment) {
    alignment.set_quality(string_quality_short_to_char(alignment.quality()));
}

string string_quality_short_to_char(const string& quality) {
    string buffer; buffer.resize(quality.size());
    for (int i = 0; i < quality.size(); ++i) {
        buffer[i] = quality_short_to_char(quality[i]);
    }
    return buffer;
}

void alignment_quality_char_to_short(Alignment& alignment) {
    alignment.set_quality(string_quality_char_to_short(alignment.quality()));
}

string string_quality_char_to_short(const string& quality) {
    string buffer; buffer.resize(quality.size());
    for (int i = 0; i < quality.size(); ++i) {
        buffer[i] = quality_char_to_short(quality[i]);
    }
    return buffer;
}

// Internal conversion function for both paired and unpaired codepaths
string alignment_to_sam_internal(const Alignment& alignment,
                                 const string& refseq,
                                 const int32_t refpos,
                                 const bool refrev,
                                 const vector<pair<int, char>>& cigar,
                                 const string& mateseq,
                                 const int32_t matepos,
                                 bool materev,
                                 const int32_t tlen,
                                 bool paired,
                                 const int32_t tlen_max) {

    // Determine flags, using orientation, next/prev fragments, and pairing status.
    int32_t flags = determine_flag(alignment, refseq, refpos, refrev, mateseq, matepos, materev, tlen, paired, tlen_max);
    
    string alignment_name;
    if (paired) {
        // We need to strip the /1 and /2 or _1 and _2 from paired reads so the two ends have the same name.
        alignment_name = regex_replace(alignment.name(), regex("[/_][12]$"), "");
    } else {
        // Keep the alignment name as is because even if the name looks paired, the reads are semantically unpaired.
        alignment_name = alignment.name();
    }
    
    // Have One True Flag for whether the read is mapped (and should have its
    // mapping stuff set) or unmapped (and should have things *'d out).
    bool mapped = !(flags & BAM_FUNMAP);
        
    if (mapped) {
        // Make sure we have everything
        assert(!refseq.empty());
        assert(refpos != -1);
        assert(!cigar.empty());
        assert(alignment.has_path());
        assert(alignment.path().mapping_size() > 0);
    }

    // We apply the convention of unmapped reads getting their mate's coordinates
    // See section 2.4.1 https://samtools.github.io/hts-specs/SAMv1.pdf
    bool use_mate_loc = !mapped && paired && !mateseq.empty();
    
    stringstream sam;
    
    sam << (!alignment_name.empty() ? alignment_name : "*") << "\t"
        << flags << "\t"
        << (mapped ? refseq : use_mate_loc ? mateseq : "*") << "\t"
        << (use_mate_loc ? matepos + 1 : refpos + 1) << "\t"
        << (mapped ? alignment.mapping_quality() : 0) << "\t"
        << (mapped ? cigar_string(cigar) : "*") << "\t"
        << (mateseq == "" ? "*" : (mateseq == refseq ? "=" : mateseq)) << "\t"
        << matepos + 1 << "\t"
        << tlen << "\t"
        // Make sure sequence always comes out in reference forward orientation by looking at the flags.
        << (!alignment.sequence().empty() ? (refrev ? reverse_complement(alignment.sequence()) : alignment.sequence()) : "*") << "\t";
    if (!alignment.quality().empty()) {
        auto quality = alignment.quality();
        if (refrev) {
            // Quality also needs to be flipped
            std::reverse(quality.begin(), quality.end());
        }
        for (int i = 0; i < quality.size(); ++i) {
            sam << quality_short_to_char(quality[i]);
        }
    } else {
        sam << "*";
    }
    //<< (alignment.has_quality() ? string_quality_short_to_char(alignment.quality()) : string(alignment.sequence().size(), 'I'));
    if (!alignment.read_group().empty()) sam << "\tRG:Z:" << alignment.read_group();
    sam << "\n";
    return sam.str();
}

int32_t determine_flag(const Alignment& alignment,
                       const string& refseq,
                       const int32_t refpos,
                       const bool refrev,
                       const string& mateseq,
                       const int32_t matepos,
                       bool materev,
                       const int32_t tlen,
                       bool paired,
                       const int32_t tlen_max) {
    
    // Determine flags, using orientation, next/prev fragments, and pairing status.
    int32_t flags = sam_flag(alignment, refrev, paired);
    
    // We've observed some reads with the unmapped flag set and also a CIGAR string set, which shouldn't happen.
    // We will check for this. The CIGAR string will only be set in the output if the alignment has a path.
    assert((bool)(flags & BAM_FUNMAP) != (alignment.has_path() && alignment.path().mapping_size()));
    
    if (!((bool)(flags & BAM_FUNMAP)) && paired && !refseq.empty() && refseq == mateseq) {
        // Properly paired if both mates mapped to same sequence, in inward-facing orientations.
        // We know they're on the same sequence, so check orientation.
        
        // If we are first, mate needs to be reverse, and if mate is first, we need to be reverse.
        // If we are at the same position either way is fine.
        bool facing = ((refpos <= matepos) && !refrev && materev) || ((matepos <= refpos) && refrev && !materev);
        
        // We are close enough if there is not tlen limit, or if there is one and we do not exceed it
        bool close_enough = (tlen_max == 0) || abs(tlen) <= tlen_max;
        
        if (facing && close_enough) {
            // We can't find anything wrong with this pair; it's properly paired.
            flags |= BAM_FPROPER_PAIR;
        }
        
        // TODO: Support sequencing technologies where "proper" pairing may
        // have a different meaning or expected combination of orientations.
    }
    
    if (paired && mateseq.empty()) {
        // Set the flag for the mate being unmapped
        flags |= BAM_FMUNMAP;
    }
    
    if (paired && materev) {
        // Set the flag for the mate being reversed
        flags |= BAM_FMREVERSE;
    }
    
    return flags;
}

string alignment_to_sam(const Alignment& alignment,
                        const string& refseq,
                        const int32_t refpos,
                        const bool refrev,
                        const vector<pair<int, char>>& cigar,
                        const string& mateseq,
                        const int32_t matepos,
                        bool materev,
                        const int32_t tlen,
                        const int32_t tlen_max) {
    
    return alignment_to_sam_internal(alignment, refseq, refpos, refrev, cigar, mateseq, matepos, materev, tlen, true, tlen_max);

}

string alignment_to_sam(const Alignment& alignment,
                        const string& refseq,
                        const int32_t refpos,
                        const bool refrev,
                        const vector<pair<int, char>>& cigar) {
    
    return alignment_to_sam_internal(alignment, refseq, refpos, refrev, cigar, "", -1, false, 0, false, 0);

}

// Internal conversion function for both paired and unpaired codepaths
bam1_t* alignment_to_bam_internal(bam_hdr_t* header,
                                  const Alignment& alignment,
                                  const string& refseq,
                                  const int32_t refpos,
                                  const bool refrev,
                                  const vector<pair<int, char>>& cigar,
                                  const string& mateseq,
                                  const int32_t matepos,
                                  bool materev,
                                  const int32_t tlen,
                                  bool paired,
                                  const int32_t tlen_max) {
    
    // this table does seem to be reproduced in htslib publicly, so I'm copying
    // it from the CRAM conversion code
    static const char nt_encoding[256] = {
        15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,
        15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,
        15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,
        15,15,15,15,15,15,15,15,15,15,15,15,15, 0,15,15,
        15, 1,14, 2,13,15,15, 4,11,15,15,12,15, 3,15,15,
        15,15, 5, 6, 8,15, 7, 9,15,10,15,15,15,15,15,15,
        15, 1,14, 2,13,15,15, 4,11,15,15,12,15, 3,15,15,
        15,15, 5, 6, 8,15, 7, 9,15,10,15,15,15,15,15,15,
        15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,
        15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,
        15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,
        15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,
        15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,
        15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,
        15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,
        15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15
    };
    
    // init an empty BAM record
    bam1_t* bam = bam_init1();
    
    // strip the pair order identifiers
    string alignment_name;
    if (paired) {
        // We need to strip the /1 and /2 or _1 and _2 from paired reads so the two ends have the same name.
        alignment_name = regex_replace(alignment.name(), regex("[/_][12]$"), "");
    } else {
        // Keep the alignment name as is because even if the name looks paired, the reads are semantically unpaired.
        alignment_name = alignment.name();
    }
    
    // calculate the size in bytes of the variable length fields (which are all concatenated in memory)
    int qname_nulls = 4 - alignment_name.size() % 4;
    int qname_data_size = alignment_name.size() + qname_nulls;
    int cigar_data_size = 4 * cigar.size();
    int seq_data_size = (alignment.sequence().size() + 1) / 2; // round up
    int qual_data_size = alignment.sequence().size(); // we will allocate this even if quality doesn't exist
    
    // allocate the joint variable length fields
    int var_field_data_size = qname_data_size + cigar_data_size + seq_data_size + qual_data_size;
    bam->data = (uint8_t*) calloc(var_field_data_size, sizeof(uint8_t));
    
    // TODO: what ID is this? CRAM seems to ignore it, so maybe we can too...
    //bam->id = 0;
    bam->l_data = var_field_data_size; // current length of data
    bam->m_data = var_field_data_size; // max length of data
    
    bam1_core_t& core = bam->core;
    // mapping position
    core.pos = refpos;
    // ID of sequence mapped to
    core.tid = sam_hdr_name2tid(header, refseq.c_str());
    // MAPQ
    core.qual = alignment.mapping_quality();
    // number of nulls (above 1) used to pad read name string
    core.l_extranul = qname_nulls - 1;
    // bit flag
    core.flag = determine_flag(alignment, refseq, refpos, refrev, mateseq, matepos, materev, tlen, paired, tlen_max);
    // length of read name, including nulls
    core.l_qname = qname_data_size;
    // number of cigar operations
    core.n_cigar = cigar.size();
    // length of read
    core.l_qseq = alignment.sequence().size();
    // ID of sequence mate is mapped to
    core.mtid = sam_hdr_name2tid(header, mateseq.c_str()); // TODO: what if there is no mate
    // mapping position of mate
    core.mpos = matepos;
    // insert length of fragment
    core.isize = tlen;
    
    // all variable-length data, concatenated; structure: qname-cigar-seq-qual-aux
    
    // write query name, padded by nulls
    uint8_t* name_data = bam->data;
    for (size_t i = 0; i < alignment_name.size(); ++i) {
        name_data[i] = (uint8_t) alignment_name[i];
    }
    for (size_t i = 0; i < qname_nulls; ++i) {
        name_data[i + alignment_name.size()] = '\0';
    }
    
    // encode cigar and copy into data

    uint32_t* cigar_data = (uint32_t*) (name_data + qname_data_size);
    
    auto refend = core.pos;
    for (size_t i = 0; i < cigar.size(); ++i) {
        uint32_t op;
        switch (cigar[i].second) {
            case 'M':
            case 'm':
                op = BAM_CMATCH;
                refend += cigar[i].first;
                break;
            case 'I':
            case 'i':
                op = BAM_CINS;
                break;
            case 'D':
            case 'd':
                op = BAM_CDEL;
                refend += cigar[i].first;
                break;
            case 'N':
            case 'n':
                op = BAM_CREF_SKIP;
                refend += cigar[i].first;
                break;
            case 'S':
            case 's':
                op = BAM_CSOFT_CLIP;
                break;
            case 'H':
            case 'h':
                op = BAM_CHARD_CLIP;
                break;
            case 'P':
            case 'p':
                op = BAM_CPAD;
                break;
            case '=':
                op = BAM_CEQUAL;
                refend += cigar[i].first;
                break;
            case 'X':
            case 'x':
                op = BAM_CDIFF;
                refend += cigar[i].first;
                break;
            default:
                throw runtime_error("Invalid CIGAR operation " + string(1, cigar[i].second));
                break;
        }
        cigar_data[i] = bam_cigar_gen(cigar[i].first, op);
    }
    
    
    // now we know where it ends, we can compute the bin
    // copied from cram/cram_samtools.h
    core.bin = hts_reg2bin(refpos, refend - 1, 14, 5); // TODO: not sure if end is past-the-last
    
    // convert sequence to 4-bit (nibble) encoding
    uint8_t* seq_data = (uint8_t*) (cigar_data + cigar.size());
    const string* seq = &alignment.sequence();
    string rev_seq;
    if (refrev) {
        rev_seq = reverse_complement(*seq);
        seq = &rev_seq;
    }
    for (size_t i = 0; i < alignment.sequence().size(); i += 2) {
        if (i + 1 < alignment.sequence().size()) {
            seq_data[i / 2] = (nt_encoding[seq->at(i)] << 4) | nt_encoding[seq->at(i + 1)];
        }
        else {
            seq_data[i / 2] = nt_encoding[seq->at(i)] << 4;
        }
    }
    
    // write the quality directly (it should already have the +33 offset removed)
    uint8_t* qual_data = seq_data + seq_data_size;
    for (size_t i = 0; i < alignment.sequence().size(); ++i) {
        if (alignment.quality().empty()) {
            // hacky, but this seems to be what they do in CRAM anyway
            qual_data[i] = '\xff';
        }
        else {
            qual_data[i] = alignment.quality().at(i);
        }
    }
    
    if (!alignment.read_group().empty()) {
        bam_aux_append(bam, "RG", 'Z', alignment.read_group().size() + 1, (uint8_t*) alignment.read_group().c_str());
    }
    // TODO: this does not seem to be a standardized field (https://samtools.github.io/hts-specs/SAMtags.pdf)
//    if (!alignment.sample_name()) {
//
//    }
        
    return bam;
}

bam1_t* alignment_to_bam(bam_hdr_t* bam_header,
                         const Alignment& alignment,
                         const string& refseq,
                         const int32_t refpos,
                         const bool refrev,
                         const vector<pair<int, char>>& cigar,
                         const string& mateseq,
                         const int32_t matepos,
                         bool materev,
                         const int32_t tlen,
                         const int32_t tlen_max) {

    return alignment_to_bam_internal(bam_header, alignment, refseq, refpos, refrev, cigar, mateseq, matepos, materev, tlen, true, tlen_max);

}

bam1_t* alignment_to_bam(bam_hdr_t* bam_header,
                         const Alignment& alignment,
                         const string& refseq,
                         const int32_t refpos,
                         const bool refrev,
                         const vector<pair<int, char>>& cigar) {
    
    return alignment_to_bam_internal(bam_header, alignment, refseq, refpos, refrev, cigar, "", -1, false, 0, false, 0);

}

string cigar_string(const vector<pair<int, char> >& cigar) {
    vector<pair<int, char> > cigar_comp;
    pair<int, char> cur = make_pair(0, '\0');
    for (auto& e : cigar) {
        if (cur == make_pair(0, '\0')) {
            cur = e;
        } else {
            if (cur.second == e.second) {
                cur.first += e.first;
            } else {
                cigar_comp.push_back(cur);
                cur = e;
            }
        }
    }
    cigar_comp.push_back(cur);
    stringstream cigarss;
    for (auto& e : cigar_comp) {
        cigarss << e.first << e.second;
    }
    return cigarss.str();
}

string mapping_string(const string& source, const Mapping& mapping) {
    string result;
    int p = mapping.position().offset();
    for (const auto& edit : mapping.edit()) {
        // mismatch/sub state
// *matches* from_length == to_length, or from_length > 0 and offset unset
// *snps* from_length == to_length; sequence = alt
        // mismatch/sub state
        if (edit.from_length() == edit.to_length()) {
            if (!edit.sequence().empty()) {
                result += edit.sequence();
            } else {
                result += source.substr(p, edit.from_length());
            }
            p += edit.from_length();
        } else if (edit.from_length() == 0 && edit.sequence().empty()) {
// *skip* from_length == 0, to_length > 0; implies "soft clip" or sequence skip
            //cigar.push_back(make_pair(edit.to_length(), 'S'));
        } else if (edit.from_length() > edit.to_length()) {
// *deletions* from_length > to_length; sequence may be unset or empty
            result += edit.sequence();
            p += edit.from_length();
        } else if (edit.from_length() < edit.to_length()) {
// *insertions* from_length < to_length; sequence contains relative insertion
            result += edit.sequence();
            p += edit.from_length();
        }
    }
    return result;
}

inline void append_cigar_operation(const int length, const char operation, vector<pair<int, char>>& cigar) {
    if (cigar.empty() || operation != cigar.back().second) {
        cigar.emplace_back(length, operation);
    }
    else {
        cigar.back().first += length;
    }
}

void mapping_cigar(const Mapping& mapping, vector<pair<int, char>>& cigar) {
    for (const auto& edit : mapping.edit()) {
        if (edit.from_length() && edit.from_length() == edit.to_length()) {
// *matches* from_length == to_length, or from_length > 0 and offset unset
            // match state
            append_cigar_operation(edit.from_length(), 'M', cigar);
            //cerr << "match " << edit.from_length() << endl;
        } else {
            // mismatch/sub state
// *snps* from_length == to_length; sequence = alt
            if (edit.from_length() == edit.to_length()) {
                append_cigar_operation(edit.from_length(), 'M', cigar);
                //cerr << "match " << edit.from_length() << endl;
            } else if (edit.from_length() > edit.to_length()) {
// *deletions* from_length > to_length; sequence may be unset or empty
                int32_t del = edit.from_length() - edit.to_length();
                int32_t eq = edit.to_length();
                if (eq) append_cigar_operation(eq, 'M', cigar);
                append_cigar_operation(del, 'D', cigar);
                //cerr << "del " << edit.from_length() - edit.to_length() << endl;
            } else if (edit.from_length() < edit.to_length()) {
// *insertions* from_length < to_length; sequence contains relative insertion
                int32_t ins = edit.to_length() - edit.from_length();
                int32_t eq = edit.from_length();
                if (eq) append_cigar_operation(eq, 'M', cigar);
                append_cigar_operation(ins, 'I', cigar);
                //cerr << "ins " << edit.to_length() - edit.from_length() << endl;
            }
        }
    }
}

int64_t cigar_mapping(const bam1_t *b, Mapping* mapping) {
    int64_t ref_length = 0;
    int64_t query_length = 0;

    const auto cigar = bam_get_cigar(b);

    for (int k = 0; k < b->core.n_cigar; k++) {
        Edit* e = mapping->add_edit();
        const int op = bam_cigar_op(cigar[k]);
        const int ol = bam_cigar_oplen(cigar[k]);
        if (bam_cigar_type(cigar[k])&1) {
            // Consume query
            e->set_to_length(ol);
            string sequence; sequence.resize(ol);
            for (int i = 0; i < ol; i++ ) {
               sequence[i] = "=ACMGRSVTWYHKDBN"[bam_seqi(bam_get_seq(b), query_length + i)];
            }
            e->set_sequence(sequence);
            query_length += ol;
        } else {
            e->set_to_length(0);
        }
        if (bam_cigar_type(cigar[k])&2) {
            // Consume ref
            e->set_from_length(ol);
            ref_length += ol;
        } else {
            e->set_from_length(0);
        }
    }
    return ref_length;
}

void mapping_against_path(Alignment& alignment, const bam1_t *b, char* chr, const PathPositionHandleGraph* graph, bool on_reverse_strand) {

    if (b->core.pos == -1) return;

    Mapping mapping;

    int64_t length = cigar_mapping(b, &mapping);

    Alignment aln = target_alignment(graph, chr, b->core.pos, b->core.pos + length, "", on_reverse_strand, mapping);

    *alignment.mutable_path() = aln.path();

    Position* refpos = alignment.add_refpos();
    refpos->set_name(chr);
    refpos->set_offset(b->core.pos);
    refpos->set_is_reverse(on_reverse_strand);
}

vector<pair<int, char>> cigar_against_path(const Alignment& alignment, bool on_reverse_strand, int64_t& pos, size_t path_len, size_t softclip_suppress) {
    vector<pair<int, char> > cigar;
    
    if (!alignment.has_path() || alignment.path().mapping_size() == 0) return cigar;
    const Path& path = alignment.path();
    int l = 0;

    for (const auto& mapping : path.mapping()) {
        mapping_cigar(mapping, cigar);
    }
    
    if(on_reverse_strand) {
        // Flip CIGAR ops into forward strand ordering
        reverse(cigar.begin(), cigar.end());
    }

    // handle soft clips, which are just insertions at the start or end
    // back
    if (cigar.back().second == 'I') {
        // make sure we stay in the reference sequence when suppressing the softclips
        if (cigar.back().first <= softclip_suppress
            && pos + alignment_from_length(alignment) + cigar.back().first <= path_len) {
            cigar.back().second = 'M';
        } else {
            cigar.back().second = 'S';
        }
    }
    // front
    if (cigar.front().second == 'I') {
        // make sure we stay in the reference sequence when suppressing the softclips
        if (cigar.front().first <= softclip_suppress
            && pos - cigar.front().first >= 0) {
            cigar.front().second = 'M';
            pos -= cigar.front().first;
        } else {
            cigar.front().second = 'S';
        }
    }

    return cigar;
}

pair<int32_t, int32_t> compute_template_lengths(const int64_t& pos1, const vector<pair<int, char>>& cigar1,
    const int64_t& pos2, const vector<pair<int, char>>& cigar2) {

    // Compute signed distance from outermost matched/mismatched base of each
    // alignment to the outermost matched/mismatched base of the other.
    
    // We work with CIGARs because it's easier than reverse complementing
    // Alignment objects without node lengths.
    
    // Work out the low and high mapped bases for each side
    auto find_bounds = [](const int64_t& pos, const vector<pair<int, char>>& cigar) {
        // Initialize bounds to represent no mapped bases
        int64_t low = numeric_limits<int64_t>::max();
        int64_t high = numeric_limits<int64_t>::min();
        
        // Track position in the reference
        int64_t here = pos;
        for (auto& item : cigar) {
            // Trace along the cigar
            if (item.second == 'M') {
                // Bases are matched. Count them in the bounds and execute the operation
                low = min(low, here);
                here += item.first;
                high = max(high, here - 1);
            } else if (item.second == 'D') {
                // Only other way to advance in the reference
                here += item.first;
            }
        }
        
        return make_pair(low, high);
    };
    
    auto bounds1 = find_bounds(pos1, cigar1);
    auto bounds2 = find_bounds(pos2, cigar2);
    
    // Compute the separation
    int32_t dist = 0;
    if (bounds1.first < bounds2.second) {
        // The reads are in order
        dist = bounds2.second - bounds1.first;
    } else if (bounds2.first < bounds1.second) {
        // The reads are out of order so the other bounds apply
        dist = bounds1.second - bounds2.first;
    }
    
    if (pos1 < pos2) {
        // Count read 1 as the overall "leftmost", so its value will be positive
        return make_pair(dist, -dist);
    } else {
        // Count read 2 as the overall leftmost
        return make_pair(-dist, dist);
    }

}

int32_t sam_flag(const Alignment& alignment, bool on_reverse_strand, bool paired) {
    int16_t flag = 0;

    if (paired) {
        // Respect the alignment's internal crossreferences.
        // Allow for multiple-read-long fragments. 
        
        flag |= BAM_FPAIRED;
        if (!alignment.has_fragment_next()) {
            // This is the last read in a pair
            flag |= BAM_FREAD2;
        }
        if (!alignment.has_fragment_prev()) {
            // This is the first read in a pair
            flag |= BAM_FREAD1;
        }
        
        // Invalid paired GAM is caught, for surject, on GAM input
        // TODO: catch reads with pair partners when they shouldn't be paired?
    }

    if (!alignment.has_path() || alignment.path().mapping_size() == 0) {
        // unmapped
        flag |= BAM_FUNMAP;
    } 
    if (on_reverse_strand) {
        flag |= BAM_FREVERSE;
    }
    if (alignment.is_secondary()) {
        flag |= BAM_FSECONDARY;
    }
    
    
    
    return flag;
}

Alignment bam_to_alignment(const bam1_t *b, map<string, string>& rg_sample, const bam_hdr_t *bh,
                           const PathPositionHandleGraph* graph) {

    Alignment alignment;

    // get the sequence and qual
    int32_t lqseq = b->core.l_qseq;
    string sequence; sequence.resize(lqseq);

    uint8_t* qualptr = bam_get_qual(b);
    string quality;//(lqseq, 0);
    quality.assign((char*)qualptr, lqseq);

    // process the sequence into chars
    uint8_t* seqptr = bam_get_seq(b);
    for (int i = 0; i < lqseq; ++i) {
        sequence[i] = "=ACMGRSVTWYHKDBN"[bam_seqi(seqptr, i)];
    }

    // get the read group and sample name
    uint8_t *rgptr = bam_aux_get(b, "RG");
    char* rg = (char*) (rgptr+1);
    //if (!rg_sample
    string sname;
    if (!rg_sample.empty()) {
        sname = rg_sample[string(rg)];
    }

    // Now name the read after the scaffold
    string read_name = bam_get_qname(b);

    // Decide if we are a first read (/1) or second (last) read (/2)
    if(b->core.flag & BAM_FREAD1) {
        read_name += "/1";
    }
    if(b->core.flag & BAM_FREAD2) {
        read_name += "/2";
    }
    
    // If we are marked as both first and last we get /1/2, and if we are marked
    // as neither the scaffold name comes through unchanged as the read name.
    // TODO: produce correct names for intermediate reads on >2 read scaffolds.

    // add features to the alignment
    alignment.set_name(read_name);
    // was the sequence reverse complemented?
    if (b->core.flag & BAM_FREVERSE) {
        
        alignment.set_sequence(reverse_complement(sequence));
        
        string rev_quality;
        rev_quality.resize(quality.size());
        reverse_copy(quality.begin(), quality.end(), rev_quality.begin());
        alignment.set_quality(rev_quality);
    }
    else {
        
        alignment.set_sequence(sequence);
        alignment.set_quality(quality);
        
    }
    
    if (graph != nullptr && bh != nullptr) {
        alignment.set_mapping_quality(b->core.qual);
        mapping_against_path(alignment, b, bh->target_name[b->core.tid], graph, b->core.flag & BAM_FREVERSE);
    }
    
    // TODO: htslib doesn't wrap this flag for some reason.
    alignment.set_is_secondary(b->core.flag & BAM_FSECONDARY);
    if (sname.size()) {
        alignment.set_sample_name(sname);
        alignment.set_read_group(rg);
    }

    return alignment;
}

Alignment bam_to_alignment(const bam1_t *b, map<string, string>& rg_sample) {
    return bam_to_alignment(b, rg_sample, nullptr, nullptr);
}

int alignment_to_length(const Alignment& a) {
    int l = 0;
    for (const auto& m : a.path().mapping()) {
        l += to_length(m);
    }
    return l;
}

int alignment_from_length(const Alignment& a) {
    int l = 0;
    for (const auto& m : a.path().mapping()) {
        l += from_length(m);
    }
    return l;
}

Alignment strip_from_start(const Alignment& aln, size_t drop) {
    if (!drop) return aln;
    Alignment res;
    res.set_name(aln.name());
    res.set_score(aln.score());
    res.set_sequence(aln.sequence().substr(drop));
    if (!aln.has_path()) return res;
    *res.mutable_path() = cut_path(aln.path(), drop).second;
    assert(res.has_path());
    if (alignment_to_length(res) != res.sequence().size()) {
        cerr << "failed!!! drop from start 轰" << endl;
        cerr << "drop " << drop << " from start" << endl << pb2json(aln) << endl;
        cerr << "wanted " << aln.sequence().size() - drop << " got " << alignment_to_length(res) << endl;
        cerr << pb2json(res) << endl << endl;
        assert(false);
    }
    return res;
}

Alignment strip_from_end(const Alignment& aln, size_t drop) {
    if (!drop) return aln;
    Alignment res;
    res.set_name(aln.name());
    res.set_score(aln.score());
    //cerr << "drop " << drop << " from end" << endl;
    size_t cut_at = aln.sequence().size()-drop;
    //cerr << "Cut at " << cut_at << endl;
    res.set_sequence(aln.sequence().substr(0, cut_at));
    if (!aln.has_path()) return res;
    *res.mutable_path() = cut_path(aln.path(), cut_at).first;
    assert(res.has_path());
    if (alignment_to_length(res) != res.sequence().size()) {
        cerr << "failed!!! drop from end 轰" << endl;
        cerr << pb2json(res) << endl << endl;
        assert(false);
    }
    return res;
}

Alignment trim_alignment(const Alignment& aln, const Position& pos1, const Position& pos2) {
    // cut the alignment into 3 (possibly empty) pieces
    auto p = cut_path(aln.path(), pos1);
    auto path1 = p.first;
    p = cut_path(p.second, pos2);
    auto path2 = p.first;
    auto path3 = p.second;
    // measure the length of the left and right bits, and use this to trim the current alignment
    auto trimmed = aln;
    if (path1.mapping_size()) {
        trimmed = strip_from_start(trimmed, path_to_length(path1));
    }
    if (path3.mapping_size()) {
        trimmed = strip_from_end(trimmed, path_to_length(path3));
    }
    return trimmed;
}

vector<Alignment> alignment_ends(const Alignment& aln, size_t len1, size_t len2) {
    vector<Alignment> ends;
    ends.push_back(strip_from_end(aln, aln.sequence().size()-len1));
    ends.push_back(strip_from_start(aln, aln.sequence().size()-len2));
    return ends;
}

Alignment alignment_middle(const Alignment& aln, int len) {
    int trim = (aln.sequence().size() - len)/2;
    return strip_from_start(strip_from_end(aln, trim), trim);
}

vector<Alignment> reverse_complement_alignments(const vector<Alignment>& alns, const function<int64_t(int64_t)>& node_length) {
    vector<Alignment> revalns;
    for (auto& aln : alns) {
        revalns.push_back(reverse_complement_alignment(aln, node_length));
    }
    return revalns;
}

Alignment reverse_complement_alignment(const Alignment& aln,
                                       const function<int64_t(id_t)>& node_length) {
    // We're going to reverse the alignment and all its mappings.
    // TODO: should we/can we do this in place?
    
    Alignment reversed = aln;
    reversed.set_sequence(reverse_complement(aln.sequence()));
    string quality = aln.quality();
    std::reverse(quality.begin(), quality.end());
    reversed.set_quality(quality);

    if(aln.has_path()) {
        // Now invert the order of the mappings, and for each mapping, flip the
        // is_reverse flag, and adjust offsets to count from the other end. The
        // edits within mappings also get put in reverse order, and get their
        // sequences reverse complemented.
        *reversed.mutable_path() = reverse_complement_path(aln.path(), node_length);
    }
    
    return reversed;
}
    
void reverse_complement_alignment_in_place(Alignment* aln,
                                           const function<int64_t(id_t)>& node_length) {

    reverse_complement_in_place(*aln->mutable_sequence());
    string* quality = aln->mutable_quality();
    std::reverse(quality->begin(), quality->end());
    
    if (aln->has_path()) {
        reverse_complement_path_in_place(aln->mutable_path(), node_length);
    }
}

// merge that properly handles long indels
// assumes that alignments should line up end-to-end
Alignment merge_alignments(const vector<Alignment>& alns) {

    if (alns.size() == 0) {
        Alignment aln;
        return aln;
    } else if (alns.size() == 1) {
        return alns.front();
    }

    // execute a serial merge
    // buliding up the alignment
    Alignment merged;
    merged.set_name(alns.front().name());

    size_t len = 0;
    for (size_t i = 0; i < alns.size(); ++i) {
        len += alns[i].sequence().size();
    }
    merged.mutable_sequence()->reserve(len);
    if (alns.front().quality().size()) merged.mutable_quality()->reserve(len);

    // get the alignments ready for merge
    for (size_t i = 0; i < alns.size(); ++i) {
        Alignment aln = alns[i];
        if (!aln.has_path()) {
            Mapping m;
            Edit* e = m.add_edit();
            e->set_to_length(aln.sequence().size());
            e->set_sequence(aln.sequence());
            *aln.mutable_path()->add_mapping() = m;
        }
        if (i == 0) {
            merged = aln;
        } else {
            if (!merged.quality().empty()) merged.mutable_quality()->append(aln.quality());
            extend_path(*merged.mutable_path(), aln.path());
            merged.mutable_sequence()->append(aln.sequence());
        }
    }
    return merged;
}

Alignment& extend_alignment(Alignment& a1, const Alignment& a2, bool debug) {
    //if (debug) cerr << "extending alignment " << endl << pb2json(a1) << endl << pb2json(a2) << endl;
    a1.set_sequence(a1.sequence() + a2.sequence());
    if (!a1.quality().empty()) a1.set_quality(a1.quality() + a2.quality());
    extend_path(*a1.mutable_path(), a2.path());
    //if (debug) cerr << "extended alignments, result is " << endl << pb2json(a1) << endl;
    return a1;
}

// use a deep copy of the alignments, concatenating them
Alignment merge_alignments(const Alignment& a1, const Alignment& a2, bool debug) {
    //cerr << "overlap is " << overlap << endl;
    // if either doesn't have a path, then treat it like a massive softclip
    if (debug) cerr << "merging alignments " << endl << pb2json(a1) << endl << pb2json(a2) << endl;
    // concatenate them
    Alignment a3;
    a3.set_name(a1.name());
    a3.set_sequence(a1.sequence() + a2.sequence());
    *a3.mutable_path() = concat_paths(a1.path(), a2.path());
    if (debug) cerr << "merged alignments, result is " << endl << pb2json(a3) << endl;
    return a3;
}

void translate_nodes(Alignment& a, const unordered_map<id_t, pair<id_t, bool> >& ids, const std::function<size_t(int64_t)>& node_length) {
    Path* path = a.mutable_path();
    for(size_t i = 0; i < path->mapping_size(); i++) {
        // Grab each mapping (includes its position)
        Mapping* mapping = path->mutable_mapping(i);
        auto pos = mapping->position();
        auto oldp = ids.find(pos.node_id());
        if (oldp != ids.end()) {
            auto& old = oldp->second;
            mapping->mutable_position()->set_node_id(old.first);
            if (old.second) {
                mapping->mutable_position()->set_is_reverse(true);
            }
        }
    }
}

void flip_nodes(Alignment& a, const set<int64_t>& ids, const std::function<size_t(int64_t)>& node_length) {
    Path* path = a.mutable_path();
    for(size_t i = 0; i < path->mapping_size(); i++) {
        // Grab each mapping (includes its position)
        Mapping* mapping = path->mutable_mapping(i);
        if(ids.count(mapping->position().node_id())) {
            // We need to flip this mapping
            *mapping = reverse_complement_mapping(*mapping, node_length);
        } 
    }
}

int non_match_start(const Alignment& alignment) {
    int length = 0;
    auto& path = alignment.path();
    for (int i = 0; i < path.mapping_size(); ++i) {
        auto& mapping = path.mapping(i);
        for (int j = 0; j < mapping.edit_size(); ++j) {
            auto& edit = mapping.edit(j);
            if (edit_is_match(edit)) {
                return length;
            }
            length += edit.to_length();
        }
    }
    return length;
}

int non_match_end(const Alignment& alignment) {
    int length = 0;
    auto& path = alignment.path();
    for (int i = path.mapping_size()-1; i >= 0; --i) {
        auto& mapping = path.mapping(i);
        for (int j = mapping.edit_size()-1; j >= 0; --j) {
            auto& edit = mapping.edit(j);
            if (edit_is_match(edit)) {
                return length;
            }
            length += edit.to_length();
        }
    }
    return length;
}

int softclip_start(const Alignment& alignment) {
    if (alignment.path().mapping_size() > 0) {
        auto& path = alignment.path();
        auto& first_mapping = path.mapping(0);
        auto& first_edit = first_mapping.edit(0);
        if (first_edit.from_length() == 0 && first_edit.to_length() > 0) {
            return first_edit.to_length();
        }
    }
    return 0;
}

int softclip_end(const Alignment& alignment) {
    if (alignment.path().mapping_size() > 0) {
        auto& path = alignment.path();
        auto& last_mapping = path.mapping(path.mapping_size()-1);
        auto& last_edit = last_mapping.edit(last_mapping.edit_size()-1);
        if (last_edit.from_length() == 0 && last_edit.to_length() > 0) {
            return last_edit.to_length();
        }
    }
    return 0;
}

int softclip_trim(Alignment& alignment) {
    // Trim the softclips off of every read
    // Work out were to cut
    int cut_start = softclip_start(alignment);
    int cut_end = softclip_end(alignment);
    // Cut the sequence and quality
    alignment.set_sequence(alignment.sequence().substr(cut_start, alignment.sequence().size() - cut_start - cut_end));
    if (alignment.quality().size() != 0) {
        alignment.set_quality(alignment.quality().substr(cut_start, alignment.quality().size() - cut_start - cut_end));
    }
    // Trim the path
    *alignment.mutable_path() = trim_hanging_ends(alignment.path());
    return cut_start + cut_end;
}

int query_overlap(const Alignment& aln1, const Alignment& aln2) {
    if (!alignment_to_length(aln1) || !alignment_to_length(aln2)
        || !aln1.path().mapping_size() || !aln2.path().mapping_size()
        || aln1.sequence().size() != aln2.sequence().size()) {
        return 0;
    }
    int qb1 = softclip_start(aln1);
    int qe1 = softclip_end(aln1);
    int qb2 = softclip_start(aln2);
    int qe2 = softclip_end(aln2);
    int l = aln1.sequence().size();
    return l - ((qe1 > qe2 ? qe1 : qe2) + (qb1 > qb2 ? qb1 : qb2));
}

int edit_count(const Alignment& alignment) {
    int i = 0;
    auto& path = alignment.path();
    for (int j = path.mapping_size(); j < path.mapping_size(); ++j) {
        i += path.mapping(j).edit_size();
    }
    return i;
}

size_t to_length_after_pos(const Alignment& aln, const Position& pos) {
    return path_to_length(cut_path(aln.path(), pos).second);
}

size_t from_length_after_pos(const Alignment& aln, const Position& pos) {
    return path_from_length(cut_path(aln.path(), pos).second);
}

size_t to_length_before_pos(const Alignment& aln, const Position& pos) {
    return path_to_length(cut_path(aln.path(), pos).first);
}

size_t from_length_before_pos(const Alignment& aln, const Position& pos) {
    return path_from_length(cut_path(aln.path(), pos).first);
}

const string hash_alignment(const Alignment& aln) {
    string data;
    aln.SerializeToString(&data);
    return sha1sum(data);
}

Alignment simplify(const Alignment& a, bool trim_internal_deletions) {
    auto aln = a;
    *aln.mutable_path() = simplify(aln.path(), trim_internal_deletions);
    if (!aln.path().mapping_size()) {
        aln.clear_path();
    }
    return aln;
}
    
void normalize_alignment(Alignment& alignment) {
    
    enum edit_type_t {None, Match, Mismatch, Insert, Delete, N};
    
    size_t cumul_to_length = 0;
    
    // we only build the normalized path if we find things we need to normalize
    // (this makes the whole algorithm a little fucky, but it should be less overhead)
    bool doing_normalization = false;
    Path normalized;
    
    const Path& path = alignment.path();
    const string& seq = alignment.sequence();
    
    auto ensure_init_normalized_path = [&](size_t i, size_t j) {
        // we won't copy the already normalized prefix unless we have to
        if (!doing_normalization) {
            for (size_t k = 0; k < i; k++) {
                *normalized.add_mapping() = path.mapping(k);
            }
            Mapping* mapping = normalized.add_mapping();
            *mapping->mutable_position() = path.mapping(i).position();
            mapping->set_rank(path.mapping_size());
            for (size_t k = 0; k < j; k++) {
                *mapping->add_edit() = path.mapping(i).edit(k);
            }
            doing_normalization = true;
        }
    };
    
    edit_type_t prev = None;
    
    for (size_t i = 0; i < path.mapping_size(); ++i) {
        
        const Mapping& mapping = path.mapping(i);
        prev = None;
        
        if (doing_normalization) {
            // we're maintaining the normalized path, so we need to add mappings
            // as we go
            Mapping* norm_mapping = normalized.add_mapping();
            *norm_mapping->mutable_position() = mapping.position();
            norm_mapping->set_rank(normalized.mapping_size());
        }
        
        for (size_t j = 0; j < mapping.edit_size(); ++j) {
            
            const Edit& edit = mapping.edit(j);
            
            if (edit.from_length() > 0 && edit.to_length() == 0) {
                
                if (prev == Delete || doing_normalization) {
                    // we need to modify the normalized path this round
                    ensure_init_normalized_path(i, j);
                    Mapping* norm_mapping = normalized.mutable_mapping(normalized.mapping_size() - 1);
                    if (prev == Delete) {
                        // merge with the previous
                        Edit* norm_edit = norm_mapping->mutable_edit(norm_mapping->edit_size() - 1);
                        norm_edit->set_from_length(norm_edit->from_length() + edit.from_length());
                    }
                    else {
                        // just copy
                        *norm_mapping->add_edit() = edit;
                    }
                }
                
                prev = Delete;
            }
            else if (edit.from_length() == 0 && edit.to_length() > 0) {
                
                if (prev == Insert || doing_normalization) {
                    // we need to modify the normalized path this round
                    ensure_init_normalized_path(i, j);
                    Mapping* norm_mapping = normalized.mutable_mapping(normalized.mapping_size() - 1);
                    if (prev == Insert) {
                        // merge with the previous
                        Edit* norm_edit = norm_mapping->mutable_edit(norm_mapping->edit_size() - 1);
                        norm_edit->set_to_length(norm_edit->to_length() + edit.to_length());
                        norm_edit->mutable_sequence()->append(edit.sequence());
                    }
                    else {
                        // just copy
                        *norm_mapping->add_edit() = edit;
                    }
                }
                
                cumul_to_length += edit.to_length();
                prev = Insert;
            }
            else {
                auto begin = seq.begin() + cumul_to_length;
                auto end = begin + edit.to_length();
                
                auto first_N = find(begin, end, 'N');
                
                edit_type_t type =  edit.sequence().empty() ? Match : Mismatch;
                
                if (prev == type || first_N != end || doing_normalization) {
                    // we have to do some normalization here
                    ensure_init_normalized_path(i, j);
                    
                    Mapping* norm_mapping = normalized.mutable_mapping(normalized.mapping_size() - 1);
                    if (first_N == end && prev != type) {
                        // just need to copy, no fancy normalization
                        *norm_mapping->add_edit() = edit;
                        prev = type;
                    }
                    else if (first_N == end) {
                        // we need to extend the previous edit, but we don't need
                        // to worry about Ns
                        Edit* norm_edit = norm_mapping->mutable_edit(norm_mapping->edit_size() - 1);
                        norm_edit->set_from_length(norm_edit->from_length() + edit.from_length());
                        norm_edit->set_to_length(norm_edit->to_length() + edit.to_length());
                        if (type == Mismatch) {
                            norm_edit->mutable_sequence()->append(edit.sequence());
                        }
                    }
                    else {
                        bool on_Ns = first_N == begin;
                        auto next_pos = begin;
                        // iterate until we've handled the whole edit sequence
                        while (next_pos != end) {
                            // find the next place where we switch from N to non-N or the reverse
                            auto next_end = find_if(next_pos, end, [&](char c) {
                                return c == 'N' != on_Ns;
                            });
                            
                            if ((prev == N && on_Ns) || (prev == type && !on_Ns)) {
                                // we need to merge with the previous edit
                                Edit* norm_edit = norm_mapping->mutable_edit(norm_mapping->edit_size() - 1);
                                norm_edit->set_from_length(norm_edit->from_length() + edit.from_length());
                                norm_edit->set_to_length(norm_edit->to_length() + edit.to_length());
                                
                                // we copy sequence for Ns and for mismatches only
                                if ((prev == N && on_Ns) || (prev == type && !on_Ns && type == Mismatch)) {
                                    norm_edit->mutable_sequence()->append(next_pos, next_end);
                                }
                            }
                            else {
                                // we can just copy
                                Edit* norm_edit = norm_mapping->add_edit();
                                norm_edit->set_from_length(next_end - next_pos);
                                norm_edit->set_to_length(next_end - next_pos);
                                *norm_edit->mutable_sequence() = string(next_pos, next_end);
                            }
                            
                            next_pos = next_end;
                            prev = on_Ns ? N : type;
                            on_Ns = !on_Ns;
                        }
                    }
                }
                else {
                    // no normalization yet
                    prev = type;
                }
                
                cumul_to_length += edit.to_length();
            }
        }
    }
    
    if (doing_normalization) {
        // we found things we needed to normalize away, so we must have built the normalized
        // path, now replace the original with it
        *alignment.mutable_path() = move(normalized);
    }
}

bool uses_Us(const Alignment& alignment) {
    
    for (char nt : alignment.sequence()) {
        switch (nt) {
            case 'U':
                return true;
                break;
                
            case 'T':
                return false;
                break;
                
            default:
                break;
        }
    }
    return false;
}

void convert_alignment_char(Alignment& alignment, char from, char to) {
    auto& seq = *alignment.mutable_sequence();
    for (size_t i = 0; i < seq.size(); ++i) {
        if (seq[i] == from) {
            seq[i] = to;
        }
    }
    if (alignment.has_path()) {
        for (Mapping& mapping : *alignment.mutable_path()->mutable_mapping()) {
            for (Edit& edit : *mapping.mutable_edit()) {
                if (!edit.sequence().empty()) {
                    auto& eseq = *edit.mutable_sequence();
                    for (size_t i = 0; i < eseq.size(); ++i) {
                        if (eseq[i] == from) {
                            eseq[i] = to;
                        }
                    }
                }
            }
        }
    }
}

void convert_Us_to_Ts(Alignment& alignment) {
    convert_alignment_char(alignment, 'U', 'T');
}

void convert_Ts_to_Us(Alignment& alignment) {
    convert_alignment_char(alignment, 'T', 'U');
}

map<id_t, int> alignment_quality_per_node(const Alignment& aln) {
    map<id_t, int> quals;
    int to_pos = 0; // offset in quals
    for (size_t i = 0; i < aln.path().mapping_size(); ++i) {
        auto& mapping = aln.path().mapping(i);
        auto to_len = mapping_to_length(mapping);
        if (mapping.has_position()) {
            auto& q = quals[mapping.position().node_id()];
            for (size_t j = 0; j < to_len; ++j) {
                q += aln.quality()[to_pos + j];
            }
        }
        to_pos += mapping_to_length(mapping);
    }
    return quals;
}

string middle_signature(const Alignment& aln, int len) {
    return signature(alignment_middle(aln, len));
}

pair<string, string> middle_signature(const Alignment& aln1, const Alignment& aln2, int len) {
    return make_pair(middle_signature(aln1, len), middle_signature(aln1, len));
}

string signature(const Alignment& aln) {
    stringstream s;
    if (aln.has_path() && aln.path().mapping_size()) {
        auto& pos1 = aln.path().mapping(0).position();
        s << pos1.node_id();
        s << (pos1.is_reverse() ? "-" : "+");
        s << ":" << pos1.offset();
        s << "_";
        auto& last = aln.path().mapping(aln.path().mapping_size()-1);
        auto& pos2 = last.position();
        s << pos2.node_id();
        s << (pos2.is_reverse() ? "-" : "+");
        s << ":" << pos2.offset() + mapping_from_length(last);
    }
    return s.str();
}

pair<string, string> signature(const Alignment& aln1, const Alignment& aln2) {
    return make_pair(signature(aln1), signature(aln2));
}

void parse_bed_regions(istream& bedstream,
                       const PathPositionHandleGraph* graph,
                       vector<Alignment>* out_alignments) {
    out_alignments->clear();
    if (!bedstream) {
        cerr << "Unable to open bed file." << endl;
        return;
    }
    string row;
    string seq;
    // Record start position
    size_t sbuf;
    // Record end position
    size_t ebuf;
    string name;
    size_t score = 0;
    string strand;

    for (int line = 1; getline(bedstream, row); ++line) {
        if (row.size() < 2 || row[0] == '#') {
            continue;
        }
        istringstream ss(row);
        ss >> seq;
        
        if (!graph->has_path(seq)) {
            // This path doesn't exist, and we'll get a segfault or worse if
            // we go look for positions in it.
            cerr << "warning: path \"" << seq << "\" not found in index, skipping" << endl;
            continue;
        }
        
        path_handle_t path_handle = graph->get_path_handle(seq);
        
        ss >> sbuf;
        ss >> ebuf;

        if (ss.fail()) {
            // Skip lines that can't be parsed
            cerr << "Error parsing bed line " << line << ": " << row << endl;
            continue;
        } 
        
        if (sbuf >= ebuf && !graph->get_is_circular(path_handle)) {
            // The start of the region can be after the end of the region only if the underlying path is circular.
            // That's not the case, so complain and skip the region.
            cerr << "warning: path \"" << seq << "\" is not circular, skipping end-spanning region on line "
                << line << ": " << row << endl;
            continue;
        }
        
        // Try parsing the optional fields. If they fail, ignore the problem, because they're optional.
        ss >> name;
        ss >> score;
        ss >> strand;

        bool is_reverse = false;
        if(!ss.fail() && strand.compare("-") == 0) {
            is_reverse = true;
        }

        // Make the Alignment
        Alignment alignment = target_alignment(graph, seq, sbuf, ebuf, name, is_reverse);
        alignment.set_score(score);

        out_alignments->push_back(alignment);
    }
}

void parse_gff_regions(istream& gffstream,
                       const PathPositionHandleGraph* graph,
                       vector<Alignment>* out_alignments) {
    out_alignments->clear();
    if (!gffstream) {
        cerr << "Unable to open gff3/gtf file." << endl;
        return;
    }
    string row;
    string seq;
    string source;
    string type;
    string buf;
    size_t sbuf;
    size_t ebuf;
    string name = "";
    string score;
    string strand;
    string num;
    string annotations;

    for (int line = 1; getline(gffstream, row); ++line) {
        if (row.size() < 2 || row[0] == '#') {
            continue;
        }
        istringstream ss(row);
        getline(ss, seq, '\t');
        getline(ss, source, '\t');
        getline(ss, type, '\t');
        getline(ss, buf, '\t');
        // Convert to 0-based 
        sbuf = atoi(buf.c_str()) - 1;
        getline(ss, buf, '\t');
        // 1-based inclusive == 0-based exclusive
        ebuf = atoi(buf.c_str());

        if (ss.fail() || !(sbuf < ebuf)) {
            cerr << "Error parsing gtf/gff line " << line << ": " << row << endl;
        } else {
            getline(ss, score, '\t');
            getline(ss, strand, '\t');
            getline(ss, num, '\t');
            getline(ss, annotations, '\t');
            vector<string> vals = split(annotations, ";");

            string name = "";

            for (auto& s : vals) {
                if (s.find("Name=") == 0) {
                    name = s.substr(5);
                }
            }

            // Skips annotations where the name can not be parsed. Empty names can 
            // results in undefinable behavior downstream. 
            if (name.empty()) {
                cerr << "warning: could not parse annotation name (Name=), skipping line " << line << endl;  
                continue;              
            }

            bool is_reverse = false;
            if(!ss.fail() && strand.compare("-") == 0) {
                is_reverse = true;
            }

            if (!graph->has_path(seq)) {
                // This path doesn't exist, and we'll get a segfault or worse if
                // we go look for positions in it.
                cerr << "warning: path \"" << seq << "\" not found in index, skipping" << endl;
            } else {
                Alignment alignment = target_alignment(graph, seq, sbuf, ebuf, name, is_reverse);

                out_alignments->push_back(alignment);
            }
        }
    }
}

Position alignment_start(const Alignment& aln) {
    Position pos;
    if (aln.path().mapping_size()) {
        pos = aln.path().mapping(0).position();
    }
    return pos;
}

Position alignment_end(const Alignment& aln) {
    Position pos;
    if (aln.path().mapping_size()) {
        auto& last = aln.path().mapping(aln.path().mapping_size()-1);
        pos = last.position();
        pos.set_offset(pos.offset() + mapping_from_length(last));
    }
    return pos;
}

map<string ,vector<pair<size_t, bool> > > alignment_refpos_to_path_offsets(const Alignment& aln) {
    map<string, vector<pair<size_t, bool> > > offsets;
    for (auto& refpos : aln.refpos()) {
        offsets[refpos.name()].push_back(make_pair(refpos.offset(), refpos.is_reverse()));
    }
    return offsets;
}

void alignment_set_distance_to_correct(Alignment& aln, const Alignment& base) {
    auto base_offsets = alignment_refpos_to_path_offsets(base);
    return alignment_set_distance_to_correct(aln, base_offsets);
}

void alignment_set_distance_to_correct(Alignment& aln, const map<string ,vector<pair<size_t, bool> > >& base_offsets) {
    auto aln_offsets = alignment_refpos_to_path_offsets(aln);
    // bail out if we can't compare
    if (!(aln_offsets.size() && base_offsets.size())) return;
    // otherwise find the minimum distance and relative orientation
    Position result;
    size_t min_distance = std::numeric_limits<size_t>::max();
    for (auto& path : aln_offsets) {
        auto& name = path.first;
        auto& aln_positions = path.second;
        auto f = base_offsets.find(name);
        if (f == base_offsets.end()) continue;
        auto& base_positions = f->second;
        for (auto& p1 : aln_positions) {
            for (auto& p2 : base_positions) {
                // disable relative inversions
                if (p1.second != p2.second) continue;
                // are they in the same orientation?
                size_t dist = abs((int64_t)p1.first - (int64_t)p2.first);
                if (dist < min_distance) {
                    min_distance = dist;
                    result.set_name(name);
                    result.set_is_reverse(p1.second != p2.second);
                    result.set_offset(dist);
                }
            }
        }
    }
    // set the distance to correct if we got one
    if (min_distance < std::numeric_limits<size_t>::max()) {
        *aln.mutable_to_correct() = result;
    }
}

bool alignment_is_valid(Alignment& aln, const HandleGraph* hgraph) {
    for (size_t i = 0; i < aln.path().mapping_size(); ++i) {
        const Mapping& mapping = aln.path().mapping(i);
        if (!hgraph->has_node(mapping.position().node_id())) {
            cerr << "Invalid Alignment:\n" << pb2json(aln) <<"\nNode " << mapping.position().node_id()
                 << " not found in graph" << endl;
            return false;
        }
        size_t node_len = hgraph->get_length(hgraph->get_handle(mapping.position().node_id()));
        if (mapping_from_length(mapping) + mapping.position().offset() > node_len) {
            cerr << "Invalid Alignment:\n" << pb2json(aln) << "\nLength of node "
                 << mapping.position().node_id() << " (" << node_len << ") exceeded by Mapping with offset "
                 << mapping.position().offset() << " and from-length " << mapping_from_length(mapping) << ":\n"
                 << pb2json(mapping) << endl;
            return false;
        }
    }
    return true;
}

Alignment target_alignment(const PathPositionHandleGraph* graph, const string& name, size_t pos1, size_t pos2,
                           const string& feature, bool is_reverse, Mapping& cigar_mapping) {
    Alignment aln;
    
    path_handle_t path_handle = graph->get_path_handle(name);
    
    if (pos2 < pos1) {
        // Looks like we want to span the origin of a circular path
        if (!graph->get_is_circular(path_handle)) {
            // But the path isn't circular, which is a problem
            throw runtime_error("Cannot extract Alignment from " + to_string(pos1) +
                                " to " + to_string(pos2) + " across the junction of non-circular path " + name);
        }
        
        // How long is the path?
        auto path_len = graph->get_path_length(path_handle);
        
        if (pos1 >= path_len) {
            // We want to start off the end of the path, which is no good.
            throw runtime_error("Cannot extract Alignment starting at " + to_string(pos1) +
                                " which is past end " + to_string(path_len) + " of path " + name);
        }
        
        if (pos2 > path_len) {
            // We want to end off the end of the path, which is no good either.
            throw runtime_error("Cannot extract Alignment ending at " + to_string(pos2) +
                                " which is past end " + to_string(path_len) + " of path " + name);
        }
        
        // Split the proivided Mapping of edits at the path end/start junction
        auto part_mappings = cut_mapping_offset(cigar_mapping, path_len - pos1);
        
        // We extract from pos1 to the end
        Alignment aln1 = target_alignment(graph, name, pos1, path_len, feature, is_reverse, part_mappings.first);
        
        // And then from the start to pos2
        Alignment aln2 = target_alignment(graph, name, 0, pos2, feature, is_reverse, part_mappings.second);
        
        if (is_reverse) {
            // The alignments were flipped, so the second has to be first
            return merge_alignments(aln2, aln1);
        } else {
            // The alignments get merged in the same order
            return merge_alignments(aln1, aln2);
        }
    }
    
    // Otherwise, the base case is that we don't go over the circular path junction
    
    
    step_handle_t step = graph->get_step_at_position(path_handle, pos1);
    size_t step_start = graph->get_position_of_step(step);
    handle_t handle = graph->get_handle_of_step(step);
    
    int64_t trim_start = pos1 - step_start;
    {
        Mapping* first_mapping = aln.mutable_path()->add_mapping();
        first_mapping->mutable_position()->set_node_id(graph->get_id(handle));
        first_mapping->mutable_position()->set_is_reverse(graph->get_is_reverse(handle));
        first_mapping->mutable_position()->set_offset(trim_start);
        
        auto mappings = cut_mapping_offset(cigar_mapping, graph->get_length(handle)-trim_start);
        first_mapping->clear_edit();
        
        string from_seq = graph->get_sequence(handle);
        int from_pos = trim_start;
        for (size_t j = 0; j < mappings.first.edit_size(); ++j) {
            if (mappings.first.edit(j).to_length() == mappings.first.edit(j).from_length()) {// if (mappings.first.edit(j).sequence() != nullptr) {
                // do the sequences match?
                // emit a stream of "SNPs" and matches
                int last_start = from_pos;
                int k = 0;
                Edit* edit;
                for (int to_pos = 0 ; to_pos < mappings.first.edit(j).to_length() ; ++to_pos, ++from_pos) {
                    //cerr << h << ":" << k << " " << from_seq[h] << " " << to_seq[k] << endl;
                    if (from_seq[from_pos] != mappings.first.edit(j).sequence()[to_pos]) {
                        // emit the last "match" region
                        if (from_pos - last_start > 0) {
                            edit = first_mapping->add_edit();
                            edit->set_from_length(from_pos-last_start);
                            edit->set_to_length(from_pos-last_start);
                        }
                        // set up the SNP
                        edit = first_mapping->add_edit();
                        edit->set_from_length(1);
                        edit->set_to_length(1);
                        edit->set_sequence(from_seq.substr(to_pos,1));
                        last_start = from_pos+1;
                    }
                }
                // handles the match at the end or the case of no SNP
                if (from_pos - last_start > 0) {
                    edit = first_mapping->add_edit();
                    edit->set_from_length(from_pos-last_start);
                    edit->set_to_length(from_pos-last_start);
                }
                // to_pos += length;
                // from_pos += length;
            } else {
                // Edit* edit = first_mapping->add_edit();
                // *edit = mappings.first.edit(j);
                *first_mapping->add_edit() = mappings.first.edit(j);
                from_pos += mappings.first.edit(j).from_length();
            }
        }
        cigar_mapping = mappings.second;
    }
    // get p to point to the next step (or past it, if we're a feature on a single node)
    int64_t p = step_start + graph->get_length(handle);
    step = graph->get_next_step(step);
    while (p < pos2) {
        handle = graph->get_handle_of_step(step);
        
        auto mappings = cut_mapping_offset(cigar_mapping, graph->get_length(handle));
        
        Mapping m;
        m.mutable_position()->set_node_id(graph->get_id(handle));
        m.mutable_position()->set_is_reverse(graph->get_is_reverse(handle));
        
        string from_seq = graph->get_sequence(handle);
        int from_pos = 0;
        for (size_t j = 0 ; j < mappings.first.edit_size(); ++j) {
            if (mappings.first.edit(j).to_length() == mappings.first.edit(j).from_length()) {
                // do the sequences match?
                // emit a stream of "SNPs" and matches
                int last_start = from_pos;
                int k = 0;
                Edit* edit;
                for (int to_pos = 0 ; to_pos < mappings.first.edit(j).to_length() ; ++to_pos, ++from_pos) {
                    //cerr << h << ":" << k << " " << from_seq[h] << " " << to_seq[k] << endl;
                    if (from_seq[from_pos] != mappings.first.edit(j).sequence()[to_pos]) {
                        // emit the last "match" region
                        if (from_pos - last_start > 0) {
                            edit = m.add_edit();
                            edit->set_from_length(from_pos-last_start);
                            edit->set_to_length(from_pos-last_start);
                        }
                        // set up the SNP
                        edit = m.add_edit();
                        edit->set_from_length(1);
                        edit->set_to_length(1);
                        edit->set_sequence(from_seq.substr(to_pos,1));
                        last_start = from_pos+1;
                    }
                }
                // handles the match at the end or the case of no SNP
                if (from_pos - last_start > 0) {
                    edit = m.add_edit();
                    edit->set_from_length(from_pos-last_start);
                    edit->set_to_length(from_pos-last_start);
                }
                // to_pos += length;
                // from_pos += length;
            } else {
                *m.add_edit() = mappings.first.edit(j);
                from_pos += mappings.first.edit(j).from_length();
            }
        }
        cigar_mapping = mappings.second;
        *aln.mutable_path()->add_mapping() = m;
        p += mapping_from_length(aln.path().mapping(aln.path().mapping_size()-1));
        step = graph->get_next_step(step);
    }
    aln.set_name(feature);
    if (is_reverse) {
        reverse_complement_alignment_in_place(&aln, [&](vg::id_t node_id) { return graph->get_length(graph->get_handle(node_id)); });
    }
    return aln;
}

Alignment target_alignment(const PathPositionHandleGraph* graph, const string& name, size_t pos1, size_t pos2,
                           const string& feature, bool is_reverse) {
    Alignment aln;
    
    path_handle_t path_handle = graph->get_path_handle(name);
    
    if (pos2 < pos1) {
        // Looks like we want to span the origin of a circular path
        if (!graph->get_is_circular(path_handle)) {
            // But the path isn't circular, which is a problem
            throw runtime_error("Cannot extract Alignment from " + to_string(pos1) +
                                " to " + to_string(pos2) + " across the junction of non-circular path " + name);
        }
        
        // How long is the path?
        auto path_len = graph->get_path_length(path_handle);
        
        if (pos1 >= path_len) {
            // We want to start off the end of the path, which is no good.
            throw runtime_error("Cannot extract Alignment starting at " + to_string(pos1) +
                                " which is past end " + to_string(path_len) + " of path " + name);
        }
        
        if (pos2 > path_len) {
            // We want to end off the end of the path, which is no good either.
            throw runtime_error("Cannot extract Alignment ending at " + to_string(pos2) +
                                " which is past end " + to_string(path_len) + " of path " + name);
        }
        
        // We extract from pos1 to the end
        Alignment aln1 = target_alignment(graph, name, pos1, path_len, feature, is_reverse);
        
        // And then from the start to pos2
        Alignment aln2 = target_alignment(graph, name, 0, pos2, feature, is_reverse);
        
        if (is_reverse) {
            // The alignments were flipped, so the second has to be first
            return merge_alignments(aln2, aln1);
        } else {
            // The alignments get merged in the same order
            return merge_alignments(aln1, aln2);
        }
    }
    
    // If we get here, we do the normal non-circular path case.
    
    step_handle_t step = graph->get_step_at_position(path_handle, pos1);
    size_t step_start = graph->get_position_of_step(step);
    handle_t handle = graph->get_handle_of_step(step);
    
    int64_t trim_start = pos1 - step_start;
    {
        Mapping* first_mapping = aln.mutable_path()->add_mapping();
        first_mapping->mutable_position()->set_node_id(graph->get_id(handle));
        first_mapping->mutable_position()->set_is_reverse(graph->get_is_reverse(handle));
        first_mapping->mutable_position()->set_offset(trim_start);
        
        Edit* e = first_mapping->add_edit();
        size_t edit_len = min<size_t>(graph->get_length(handle) - trim_start, pos2 - pos1);
        e->set_from_length(edit_len);
        e->set_to_length(edit_len);
    }
    // get p to point to the next step (or past it, if we're a feature on a single node)
    int64_t p = step_start + graph->get_length(handle);
    step = graph->get_next_step(step);
    while (p < pos2) {
        handle = graph->get_handle_of_step(step);
        
        Mapping* m = aln.mutable_path()->add_mapping();
        m->mutable_position()->set_node_id(graph->get_id(handle));
        m->mutable_position()->set_is_reverse(graph->get_is_reverse(handle));
        
        Edit* e = m->add_edit();
        size_t edit_len = min<size_t>(graph->get_length(handle), pos2 - p);
        e->set_from_length(edit_len);
        e->set_to_length(edit_len);
        
        p += graph->get_length(handle);
        step = graph->get_next_step(step);
    }
    
    aln.set_name(feature);
    if (is_reverse) {
        reverse_complement_alignment_in_place(&aln, [&](vg::id_t node_id) { return graph->get_length(graph->get_handle(node_id)); });
    }
    return aln;
}
}
