#!/usr/bin/env bash

BASH_TAP_ROOT=../deps/bash-tap
. ../deps/bash-tap/bash-tap-bootstrap

PATH=../bin:$PATH # for vg

export LC_ALL="en_US.utf8" # force ekg's favorite sort order 

plan tests 51

# Single graph without haplotypes
vg construct -r small/x.fa -v small/x.vcf.gz > x.vg

vg index -x x.xg x.vg
is $? 0 "building an XG index of a graph"

vg index -g x.gcsa x.vg
is $? 0 "building a GCSA index of a graph"

vg index -x x2.xg -g x2.gcsa x.vg
is $? 0 "building both indexes at once"

cmp x.xg x2.xg && cmp x.gcsa x2.gcsa && cmp x.gcsa.lcp x2.gcsa.lcp
is $? 0 "the indexes are identical"

rm -f x.vg
rm -f x.xg x.gcsa x.gcsa.lcp
rm -f x2.xg x2.gcsa x2.gcsa.lcp


# Single graph with haplotypes
vg construct -r small/x.fa -v small/x.vcf.gz -a > x.vg

vg index -G x.gbwt -v small/x.vcf.gz -F x.threads x.vg
is $? 0 "building a GBWT index of a graph with haplotypes"

vg index -x x.xg -F x.threads x.vg
is $? 0 "building an XG index of a graph with haplotypes"

vg index -g x.gcsa x.vg
is $? 0 "building a GCSA index of a graph with haplotypes"

vg index -x x2.xg -G x2.gbwt -v small/x.vcf.gz -g x2.gcsa x.vg
is $? 0 "building all indexes at once"

cmp x.xg x2.xg && cmp x.gbwt x2.gbwt && cmp x.gcsa x2.gcsa && cmp x.gcsa.lcp x2.gcsa.lcp
is $? 0 "the indexes are identical"

rm -f x.vg
rm -f x.threads
rm -f x.xg x.gbwtx.gcsa x.gcsa.lcp
rm -f x2.xg x2.gbwt x2.gcsa x2.gcsa.lcp

# Subregion graph with haplotypes
vg construct -r small/x.fa -v small/x.vcf.gz -a --region x:100-200 > x.part.vg

vg index -x x.part.xg -G x.part.gbwt --region x:100-200 -v small/x.vcf.gz x.part.vg 2>log.txt
is $? 0 "building GBWT index for a regional graph"

is "$(cat log.txt | wc -c)" "0" "no warnings about missing variants produced"

rm -f x.part.vg x.part.xg x.part.gbwt log.txt

# Multiple graphs without haplotypes
vg construct -r small/xy.fa -v small/xy2.vcf.gz -R x -C > x.vg 2> /dev/null
vg construct -r small/xy.fa -v small/xy2.vcf.gz -R y -C > y.vg 2> /dev/null
vg ids -j x.vg y.vg

vg index -x xy.xg x.vg y.vg
is $? 0 "building an XG index of multiple graphs"

vg index -g xy.gcsa -k 2 x.vg y.vg
is $? 0 "building a GCSA index of multiple graphs"

vg index -x xy2.xg -g xy2.gcsa -k 2 x.vg y.vg
is $? 0 "building both indexes at once"

cmp xy.xg xy2.xg && cmp xy.gcsa xy2.gcsa && cmp xy.gcsa.lcp xy2.gcsa.lcp
is $? 0 "the indexes are identical"

rm -f x.vg y.vg
rm -f xy.xg xy.gcsa xy.gcsa.lcp
rm -f xy2.xg xy2.gcsa xy2.gcsa.lcp


# Multiple graphs with haplotypes
vg construct -r small/xy.fa -v small/xy2.vcf.gz -R x -C -a > x.vg 2> /dev/null
vg construct -r small/xy.fa -v small/xy2.vcf.gz -R y -C -a > y.vg 2> /dev/null
vg ids -j x.vg y.vg

vg index -G x.gbwt -v small/xy2.vcf.gz -F x.threads x.vg && vg index -G y.gbwt -v small/xy2.vcf.gz -F y.threads y.vg && vg gbwt -m -f -o xy.gbwt x.gbwt y.gbwt
is $? 0 "building a GBWT index of multiple graphs with haplotypes"

vg index -x xy.xg -F x.threads -F y.threads x.vg y.vg
is $? 0 "building an XG index of multiple graphs with haplotypes"

vg index -g xy.gcsa -k 2 x.vg y.vg
is $? 0 "building a GCSA index of multiple graphs with haplotypes"

vg index -x xy2.xg -g xy2.gcsa -k 2 -G xy2.gbwt -v small/xy2.vcf.gz x.vg y.vg
is $? 0 "building all three indexes at once"

cmp xy.xg xy2.xg && cmp xy.gcsa xy2.gcsa && cmp xy.gcsa.lcp xy2.gcsa.lcp && cmp xy.gbwt xy2.gbwt
is $? 0 "the indexes are identical"

rm -f x.vg y.vg
rm -f x.gbwt y.gbwt x.threads y.threads
rm -f xy.xg xy.gbwt xy.gcsa xy.gcsa.lcp
rm -f xy2.xg xy2.gbwt xy2.gcsa xy2.gcsa.lcp


# GBWT construction options
vg construct -r small/xy.fa -v small/xy2.vcf.gz -R x -C -a > x.vg 2> /dev/null

vg index -G x_ref.gbwt -T x.vg
is $? 0 "GBWT can be built for paths"

vg index -G x_both.gbwt -T -v small/xy2.vcf.gz x.vg
is $? 0 "GBWT can be built for both paths and haplotypes"

rm -f x_ref.gbwt x_both.gbwt

vg index -x x.xg x.vg
vg sim -s 2384259 -n 100 -l 100 -x x.xg -a >sim.gam
vg index -G x_gam.gbwt -M sim.gam -x x_gam.xg x.vg

is $(vg paths -g x_gam.gbwt -T -x x_gam.xg -V | vg view -c - | jq -cr '.path[].name'  | sort | md5sum | cut -f 1 -d\ ) $(vg view -a sim.gam | jq -r .name | sort | md5sum | cut -f 1 -d\ ) "we can build a GBWT from alignments and index it by name with xg thread naming"

rm -f x.vg x.xg sim.gam x_gam.gbwt

# Other tests
vg construct -r small/x.fa -v small/x.vcf.gz >x.vg
vg index -x x.xg x.vg bogus123.vg
is $? 134 "fail with nonexistent file"
rm -rf x.idx

vg kmers -k 16 -gB x.vg >x.graph
vg index -i x.graph -g x.gcsa
is $? 0 "a prebuilt deBruijn graph in GCSA2 format may be used"
rm -f x.gcsa x.gcsa.lcp x.graph

vg index -x x.xg -g x.gcsa -k 11 x.vg
vg map -T <(vg sim -s 1337 -n 100 -x x.xg) -d x > x1337.gam
vg index -a x1337.gam -d x.vg.aln
is $(vg index -D -d x.vg.aln | wc -l) 101 "index can store alignments"
is $(vg index -A -d x.vg.aln | vg view -a - | wc -l) 100 "index can dump alignments"

# repeat with an unmapped read (sequence from phiX)
rm -rf x.vg.aln
(vg map -s "CTGATGAGGCCGCCCCTAGTTTTGTTTCTGGTGCTATGGCTAAAGCTGGTAAAGGACTTC" -d x -P 0.9; vg map -s "CTGATGAGGCCGCCCCTAGTTTTGTTTCTGGTGCTATGGCTAAAGCTGGTAAAGGACTTC" -d x ) | vg index -a - -d x.vg.aln
is $(vg index -A -d x.vg.aln | vg view -a - | wc -l) 2 "alignment index stores unmapped reads"

# load the same data again & see that the records are duplicated.
# that's not really a wise use-case, but it tests that we don't
# overwrite anything in an existing alignment index, by virtue
# of keying each entry with a unqiue nonce.
rm -rf x.vg.aln
vg index -a x1337.gam -d x.vg.aln
vg index -a x1337.gam -d x.vg.aln
is $(vg index -D -d x.vg.aln | wc -l) 201 "alignment index can be loaded using sequential invocations; next_nonce persistence"

vg map -T <(vg sim -s 1337 -n 100 -x x.xg) -d x | vg index -m - -d x.vg.map
is $(vg index -D -d x.vg.map | wc -l) 1477 "index stores all mappings"

rm -rf x.idx x.vg.map x.vg.aln x1337.gam

vg construct -r small/x.fa -v small/x.vcf.gz -a >x.vg
vg index -x x.xg -v small/x.vcf.gz x.vg
is $? 0 "building an xg index containing a gPBWT"

vg find -t -x x.xg >part.vg
is "$(cat x.vg part.vg | vg view -j - | jq '.path[].name' | grep '_thread' | wc -l)" 2 "the gPBWT can be queried for two threads for each haplotype"

is $(vg find -x x.xg -q _thread_1_x_0 | vg paths -L -v - | wc -l) 1 "a specific thread may be pulled from the graph by name"

vg index -x x.xg -v small/x.vcf.gz x.vg --exclude 1
vg find -t -x x.xg >part.vg
is "$(cat x.vg part.vg | vg view -j - | jq '.path[].name' | grep '_thread' | wc -l)" 0 "samples can be excluded from haplotype indexing"

rm -f x.vg x.xg part.vg x.gcsa


vg construct -r small/xy.fa -v small/xy.vcf.gz -a >xy.vg
vg index -x xy.xg -v small/xy.vcf.gz xy.vg
is $(vg find -x xy.xg -t | vg paths -L -v - | wc -l) 4 "a thread is stored per haplotype, sample, and reference sequence"
is $(vg find -x xy.xg -q _thread_1_y | vg paths -L -v - | wc -l) 2 "we have the expected number of threads per chromosome"
rm -f xy.vg xy.xg

vg construct -r small/x.fa -v small/x.vcf.gz -a >x.vg
vg index -x x.xg -v small/x.vcf.gz -H haps.bin x.vg
is $(du -b haps.bin | cut -f 1) 345 "threads may be exported to binary for use in GBWT construction"

rm -f x.vg x.xg part.vg x.gcsa haps.bin x.gbwt

vg construct -r small/x.fa -v small/x.vcf.gz >x.vg
vg construct -r small/x.fa -v small/x.vcf.gz >y.vg
vg construct -r small/x.fa -v small/x.vcf.gz >z.vg

vg concat x.vg y.vg z.vg >q.vg

vg index -x q.xg q.vg

is $? 0 "storage of multiple graphs in an index succeeds"

rm x.vg y.vg z.vg q.vg
rm -rf x.idx q.xg

# Now test backward nodes
vg index -x r.xg reversing/reversing_x.vg
is $? 0 "can index backward nodes"

vg index -k 11 -g r.gcsa reversing/reversing_x.vg
is $? 0 "can index kmers for backward nodes"

vg map -T <(vg sim -s 1338 -n 100 -x r.xg) -d r | vg index -a - -d r.aln.idx
is $(vg index -D -d r.aln.idx | wc -l) 101 "index can store alignments to backward nodes"

rm -rf r.aln.idx r.xg r.gcsa

vg index -x c.xg -g c.gcsa -k 11 cyclic/all.vg
vg map -T <(vg sim -s 1337 -n 100 -x c.xg) -d c | vg index -a - -d all.vg.aln
is $(vg index -D -d all.vg.aln | wc -l) 101 "index can store alignments to cyclic graphs"

rm -rf all.vg.aln c.xg c.gcsa

is $(vg index -g x.gcsa -k 16 -V <(vg view -Fv cyclic/two_node.gfa) 2>&1 |  grep 'Index verification complete' | wc -l) 1 "GCSA2 index works on cyclic graphs with heads and tails"

is $(vg index -g x.gcsa -k 16 -V cyclic/no_heads.vg 2>&1 |  grep 'Index verification complete' | wc -l) 1 "GCSA2 index works on cyclic graphs with no heads or tails"

is $(vg index -g x.gcsa -k 16 -V cyclic/self_loops.vg 2>&1 |  grep 'Index verification complete' | wc -l) 1 "GCSA2 index works on cyclic graphs with self loops"

is $(vg index -g x.gcsa -k 16 -V cyclic/all.vg 2>&1 |  grep 'Index verification complete' | wc -l) 1 "GCSA2 index works on general cyclic graphs"

rm -f x.gcsa x.gcsa.lcp

is $(vg construct -r tiny/tiny.fa -v tiny/tiny.vcf.gz | vg index -g t.gcsa -k 16 -V - 2>&1 |  grep 'Index verification complete' | wc -l) 1 "GCSA2 indexing of a tiny graph works"

is $(vg construct -r tiny/tiny.fa | vg index -g t.gcsa -k 16 -V - 2>&1 | grep 'Index verification complete' | wc -l) 1 "GCSA2 indexing succeeds on a single-node graph"

is $(vg index -g t.gcsa reversing/cactus.vg -k 16 -V 2>&1 | grep 'Index verification complete' | wc -l) 1 "GCSA2 indexing succeeds on graph with heads but no tails"

vg construct -r ins_and_del/ins_and_del.fa -v ins_and_del/ins_and_del.vcf.gz -a >ins_and_del.vg
is $(vg index -x ins_and_del.vg.xg -v ins_and_del/ins_and_del.vcf.gz ins_and_del.vg 2>&1 | wc -l) 0 "indexing with allele paths handles combination insert-and-deletes"

vg construct -m 1025 -r 1mb1kgp/z.fa > big.vg

is $(vg index -g big.gcsa big.vg -k 16 2>&1 | head -n10 | grep 'Found kmer with offset' | wc -l) 1 "a useful error message is produced when nodes are too large"

rm -f big.vg

rm -f t.gcsa
rm -f x.vg

rm -f r.gcsa.lcp c.gcsa.lcp t.gcsa.lcp ins_and_del.vg ins_and_del.vg.xg
