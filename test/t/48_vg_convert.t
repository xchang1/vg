#!/usr/bin/env bash

BASH_TAP_ROOT=../deps/bash-tap
. ../deps/bash-tap/bash-tap-bootstrap

PATH=../bin:$PATH # for vg

plan tests 22

vg construct -r complex/c.fa -v complex/c.vcf.gz > c.vg
cat <(vg view c.vg | grep ^S | sort) <(vg view c.vg | grep L | uniq | wc -l) <(vg paths -v c.vg -E) > c.info

vg convert c.vg -x > c.xg
vg convert c.xg -v > c1.vg
cat <(vg view c1.vg | grep ^S | sort) <(vg view c1.vg | grep L | uniq | wc -l) <(vg paths -v c1.vg -E) > c1.info
diff c.info c1.info
is "$?" 0 "vg convert maintains same nodes throughout xg conversion"

rm -f c.xg c1.vg c1.info

vg convert c.vg -a > c.hg
vg convert c.hg -v > c1.vg
cat <(vg view c1.vg | grep ^S | sort) <(vg view c1.vg | grep L | uniq | wc -l) <(vg paths -v c1.vg -E) > c1.info
diff c.info c1.info
is "$?" 0 "vg convert maintains same nodes throughout hash-graph conversion"

rm -f c.hg c1.vg c1.info

vg convert c.vg -p > c.pg
vg convert c.pg -v > c1.vg
cat <(vg view c1.vg | grep ^S | sort) <(vg view c1.vg | grep L | uniq | wc -l) <(vg paths -v c1.vg -E) > c1.info
diff c.info c1.info
is "$?" 0 "vg convert maintains same nodes throughout packed-graph conversion"

rm -f c.pg c1.vg c1.info

vg convert c.vg -o > c.odgi
vg convert c.odgi -v > c1.vg
cat <(vg view c1.vg | grep ^S | sort) <(vg view c1.vg | grep L | uniq | wc -l) <(vg paths -v c1.vg -E) > c1.info
diff c.info c1.info
is "$?" 0 "vg convert maintains same nodes throughout ODGI conversion"

rm -f c.vg c.odgi c1.vg c.info c1.info

# some less rigorous tests I made without noticing that the earlier ones had already been written
vg construct -r small/x.fa -v small/x.vcf.gz > x.vg
vg view x.vg > x.gfa

is "$(vg convert -a x.vg | vg view - | wc -l)" "$(wc -l < x.gfa)" "hash graph conversion looks good"
is "$(vg convert -p x.vg | vg view - | wc -l)" "$(wc -l < x.gfa)" "packed graph conversion looks good"
is "$(vg convert -v x.vg | vg view - | wc -l)" "$(wc -l < x.gfa)" "vg conversion looks good"
is "$(vg convert -o x.vg | vg view - | wc -l)" "$(wc -l < x.gfa)" "odgi conversion looks good"
is "$(vg convert -x x.vg | vg find -n 1 -c 300 -x - | vg view - | wc -l)" "$(wc -l < x.gfa)" "xg conversion looks good"

is "$(vg convert -g -a x.gfa | vg view - | wc -l)" "$(wc -l < x.gfa)" "on disk gfa conversion looks good"
is "$(cat x.gfa | vg convert -g -a - | vg view - | wc -l)" "$(wc -l < x.gfa)" "streaming gfa conversion looks good"
is "$(vg convert -g -x x.gfa | vg find -n 1 -c 300 -x - | vg view - | wc -l)" "$(wc -l < x.gfa)" "gfa to xg conversion looks good"

rm x.vg x.gfa
rm -f c.vg c.pg c1.vg c.info c1.info

vg construct -r small/x.fa -v small/x.vcf.gz > x.vg
vg index x.vg -g x.gcsa
vg sim -x x.vg -n 10 -s 23 -a > sim.gam
vg map -x x.vg -g x.gcsa -G sim.gam > sim-rm.gam
vg convert x.vg -G sim-rm.gam -t 1 > sim-rm.gaf
vg convert x.vg -F sim-rm.gaf -t 1 | vg convert x.vg -G - -t 1 > sim-rm2.gaf
diff sim-rm.gaf sim-rm2.gaf
is "$?" 0 "vg convert gam -> gaf -> gam -> gaf makes same gaf twice"

vg convert x.vg -G sim-rm.gam | vg convert x.vg -F - | vg convert x.vg -G - | sort > sim-rm2-mt-sort.gaf
sort sim-rm2.gaf > sim-rm2-sort.gaf
diff sim-rm2-sort.gaf sim-rm2-mt-sort.gaf
is "$?" 0 "vg convert gam -> gaf -> gam -> gaf gives same result multithreaded as with -t 1"

vg convert x.vg -G sim-rm.gam | bgzip | vg convert x.vg -F - | vg convert x.vg -G - | sort > sim-rm2-mtbg-sort.gaf
diff sim-rm2-sort.gaf sim-rm2-mtbg-sort.gaf
is "$?" 0 "vg convert gam -> gaf.gz -> gam -> gaf gives same result multithreaded as with -t 1"

# some snps and indels
vg map -s "TAATGGATATGTTAAGCTTTTTTTTTCTTTGATTTATTTGAAAAGACGTTTGACAATCTATCGGGTAATGTGGGGAAA" -x x.vg -g x.gcsa > mut.gam
# reverse complement of above
vg map -s "TTTCCCCACATTACCCGATAGATTGTCAAACGTCTTTTCAAATAAATCAAAGAAAAAAAAAGCTTAACATATCCATTA" -x x.vg -g x.gcsa >> mut.gam

vg convert x.vg -G mut.gam -t 1 > mut.gaf
vg convert x.vg -F mut.gaf -t 1 > mut-back.gam
vg view -a mut.gam | jq .path > mut.path
vg view -a mut-back.gam | jq .path > mut-back.path
# Json comparison that is not order dependent: https://stackoverflow.com/a/31933234
is $(jq --argfile a mut.path --argfile b mut-back.path -n 'def post_recurse(f): def r: (f | select(. != null) | r), .; r; def post_recurse: post_recurse(.[]?); ($a | (post_recurse | arrays) |= sort) as $a | ($b | (post_recurse | arrays) |= sort) as $b | $a == $b') true "vg convert gam -> gaf -> gam produces same gam Paths with snps and indels"

vg view -a mut.gam | jq .sequence > mut.seq
vg view -a mut-back.gam | jq .sequence > mut-back.seq
diff mut.seq mut-back.seq
is "$?" 0 "vg convert gam -> gaf -> gam preserves sequence"

vg convert x.vg -G mut-back.gam -t 1 > mut-back.gaf
diff mut.gaf mut-back.gaf
is "$?" 0 "vg convert gam -> gaf -> gam -> gaf makes same gaf twice in presence of indels and snps"

rm -f x.vg x.gcsa sim.gam sim-rm.gam sim-rm.gaf sim-rm2.gaf sim-rm2-mt-sort.gaf sim-rm2-mtbg-sort.gaf sim-rm2-sort.gaf mut.gam mut-back.gam mut.gaf mut-back.gaf mut.path mut-back.path mut.seq mut-back.seq

vg construct -r 1mb1kgp/z.fa -v 1mb1kgp/z.vcf.gz > z.vg 2> /dev/null
vg sim -n 10000 -s 23 -a -x z.vg > sim.gam
vg construct -r 1mb1kgp/z.fa > zflat.vg
vg index zflat.vg -g zflat.gcsa
vg map -x zflat.vg -g zflat.gcsa -G sim.gam > sim-map.gam
vg convert zflat.vg -G sim-map.gam | bgzip > sim-map.gaf.gz
vg convert zflat.vg -F sim-map.gaf.gz > sim-map-back.gam
vg view -a sim-map.gam | jq .sequence | sort > sim-map.sequence
vg view -a sim-map-back.gam | jq .sequence | sort > sim-map-back.sequence
diff sim-map.sequence sim-map-back.sequence
is "$?" 0 "vg convert gam -> gaf -> gam preserves sequences of 1mb1kgp simulated reads"

vg convert zflat.vg -G sim-map-back.gam | sort > sim-map-back.gaf
bgzip -dc sim-map.gaf.gz | sort > sim-map.gaf
diff sim-map-back.gaf sim-map.gaf
is "$?" 0 "vg convert gam -> gaf -> gam ->gaf makes same gaf each time on 1mb1kgp simulated reads"

printf '{"name": "split", "path": {"mapping": [{"edit": [{"from_length": 13, "to_length": 13}], "position": {"node_id": "1", "offset": "10"}}, {"edit": [{"from_length": 2, "to_length": 2}], "position": {"node_id": "3", "offset": "5"}}]}}' | vg view -JaG - > split.gam
vg convert zflat.vg -G split.gam > split.gaf
is "$(awk '{print $14}' split.gaf)" "cs:Z::13-CCAGTGCTC-GCATC:2" "split alignment converted using deletions to represent internal offsets"
vg convert zflat.vg -F split.gaf | vg convert zflat.vg -G - > split-back.gaf
diff split.gaf split-back.gaf
is "$?" 0 "vg convert gam -> gaf ->gam -> gaf makes same gaf each time for split alignment"

rm -f z.vg zflat.vg sim.gam sim-map.gam sim-map-back.gam sim-map.gaf.gz sim-map.sequence sim-map-back.sequence sim-map-back.gaf sim-map.gaf split.gam split.gaf split-back.gaf

