#!/usr/bin/env bash

# giraffe-wrangler.sh: Run and profile vg gaffe and analyze the results.

set -e

usage() {
    # Print usage to stderr
    exec 1>&2
    printf "Usage: $0 [Options] FASTA XG_INDEX GCSA_INDEX GBWT_INDEX MINIMIZER_INDEX DISTANCE_INDEX SIM_GAM REAL_FASTQ\n"
    printf "\n"
    printf "Options:\n"
    printf "  -t N  Use N threads\n"
    printf "\n"
    exit 1
}

# Define the thread count for everyone. Can be changed with -t.
# Should fit on a NUMA node
THREAD_COUNT=24

while getopts ":t:" o; do
    case "${o}" in
        t)
            THREAD_COUNT=$OPTARG
            ;;
        ?)
            usage
            ;;
    esac
done
shift $((OPTIND-1))
if [[ "$#" -lt "8" ]]; then
    # Too few arguments
    usage
fi

echo "Using ${THREAD_COUNT} threads"

fetch_input() {
    # Download the specified file, if not a file already.
    # Dumps all files into the current directory as their basenames
    # Output the new filename
    if [[ "${1}" == s3://* ]] ; then
        aws s3 --quiet cp "${1}" "$(basename "${1}")"
        basename "${1}"
    else
        echo "${1}"
    fi
}

FASTA="$(fetch_input "${1}")"
for EXT in amb ann bwt fai pac sa ; do
    # Make sure we have all the indexes adjacent to the FASTA
    fetch_input "${1}.${EXT}" >/dev/null
done
shift
XG_INDEX="$(fetch_input "${1}")"
# Make sure we have the GBWTGraph pre-made
GBWT_GRAPH="$(fetch_input "${1%.xg}.gg")"
shift
GCSA_INDEX="$(fetch_input "${1}")"
LCP_INDEX="$(fetch_input "${1}.lcp")"
shift
GBWT_INDEX="$(fetch_input "${1}")"
shift
MINIMIZER_INDEX="$(fetch_input "${1}")"
shift
DISTANCE_INDEX="$(fetch_input "${1}")"
shift
SIM_GAM="$(fetch_input "${1}")"
shift
REAL_FASTQ="$(fetch_input "${1}")"
shift

if [ -f "$GBWT_GRAPH" ]; then
    GIRAFFE_GRAPH=(-g "${GBWT_GRAPH}")
else
    GIRAFFE_GRAPH=(-x "${XG_INDEX}")
fi

echo "Indexes:"
echo "${XG_INDEX}"
echo "${GBWT_GRAPH}"
echo "${GCSA_INDEX}"
echo "${LCP_INDEX}"
echo "${GBWT_INDEX}"
echo "${MINIMIZER_INDEX}"
echo "${DISTANCE_INDEX}"
echo "${SIM_GAM}"
echo "${REAL_FASTQ}"

# Define the Giraffe parameters
GIRAFFE_OPTS=(-s75 -u 0.1 -v 1 -w 5 -C 600)

# Define a work directory
# TODO: this requires GNU mptemp
WORK="$(mktemp -d)"

# Check for NUMA. If we have NUMA and no numactl results may be unreliable
NUMA_COUNT="$(lscpu | grep "NUMA node(s)" | cut -f3- -d' ' | tr -d ' ')"
NUMA_PREFIX=""
NUMA_WARNING=0

if [[ "${NUMA_COUNT}" -gt "1" ]] ; then
    if which numactl >/dev/null 2>&1 ; then
        # Run everything on one NUMA node
        NUMA_PREFIX="numactl --cpunodebind=0 --membind=0"
    else
        # We should warn in the report that NUMA may confound the results
        NUMA_WARNING=1
    fi
fi

if which perf >/dev/null 2>&1 ; then
    # Record profile.
    # Do this first because perf is likely to be misconfigured and we want to fail fast.
    
    # If we don't strip bin/vg to make it small, the addr2line calls that perf
    # script makes take forever because the binary is huge
    strip -d bin/vg
    
    ${NUMA_PREFIX} perf record -F 100 --call-graph dwarf -o "${WORK}/perf.data"  vg gaffe "${GIRAFFE_GRAPH[@]}" -m "${MINIMIZER_INDEX}" -H "${GBWT_INDEX}" -d "${DISTANCE_INDEX}" -f "${REAL_FASTQ}" -t "${THREAD_COUNT}" "${GIRAFFE_OPTS[@]}" >"${WORK}/perf.gam"
    perf script -i "${WORK}/perf.data" >"${WORK}/out.perf"
    deps/FlameGraph/stackcollapse-perf.pl "${WORK}/out.perf" >"${WORK}/out.folded"
    deps/FlameGraph/flamegraph.pl "${WORK}/out.folded" > "${WORK}/profile.svg"
fi

# Run simulated reads, with stats
${NUMA_PREFIX} vg gaffe --track-correctness -x "${XG_INDEX}" "${GIRAFFE_GRAPH[@]}" -m "${MINIMIZER_INDEX}" -H "${GBWT_INDEX}" -d "${DISTANCE_INDEX}" -G "${SIM_GAM}" -t "${THREAD_COUNT}" "${GIRAFFE_OPTS[@]}" >"${WORK}/mapped.gam"

# And map to compare with them
${NUMA_PREFIX} vg map -x "${XG_INDEX}" -g "${GCSA_INDEX}" -G "${SIM_GAM}" -t "${THREAD_COUNT}" >"${WORK}/mapped-map.gam"

# Annotate and compare against truth
vg annotate -p -x "${XG_INDEX}" -a "${WORK}/mapped.gam" >"${WORK}/annotated.gam"
vg annotate -p -x "${XG_INDEX}" -a "${WORK}/mapped-map.gam" >"${WORK}/annotated-map.gam"

# GAM compare against truth. Use gamcompare to count correct reads to save a JSON scan.
CORRECT_COUNT="$(vg gamcompare -r 100 "${WORK}/annotated.gam" "${SIM_GAM}" 2>&1 >/dev/null | sed 's/[^0-9]//g')"
CORRECT_COUNT_MAP="$(vg gamcompare -r 100 "${WORK}/annotated-map.gam" "${SIM_GAM}" 2>&1 >/dev/null | sed 's/[^0-9]//g')"

# Compute identity of mapped reads
MEAN_IDENTITY="$(vg view -aj "${WORK}/mapped.gam" | jq -c 'select(.path) | .identity' | awk '{x+=$1} END {print x/NR}')"
MEAN_IDENTITY_MAP="$(vg view -aj "${WORK}/mapped-map.gam" | jq -c 'select(.path) | .identity' | awk '{x+=$1} END {print x/NR}')"

# Compute loss stages
# Let giraffe facts errors out
vg view -aj "${WORK}/mapped.gam" | scripts/giraffe-facts.py "${WORK}/facts" >"${WORK}/facts.txt"

# Now do the real reads

# Get RPS
${NUMA_PREFIX} vg gaffe -p "${GIRAFFE_GRAPH[@]}" -m "${MINIMIZER_INDEX}" -H "${GBWT_INDEX}" -d "${DISTANCE_INDEX}" -f "${REAL_FASTQ}" -t "${THREAD_COUNT}" "${GIRAFFE_OPTS[@]}" >"${WORK}/real.gam" 2>"${WORK}/log.txt"

GIRAFFE_RPS="$(cat "${WORK}/log.txt" | grep "reads per second" | sed 's/[^0-9.]//g')"

# Get RPS for bwa-mem
REAL_READ_COUNT="$(cat "${REAL_FASTQ}" | wc -l)"
((REAL_READ_COUNT /= 4))

${NUMA_PREFIX} bwa mem -t "${THREAD_COUNT}" "${FASTA}" "${REAL_FASTQ}" >"${WORK}/mapped.bam" 2>"${WORK}/bwa-log.txt"

# Now we get all the batch times from BWA and use those to compute RPS values.
# This is optimistic but hopefully consistent.
BWA_RPS_ALL_THREADS="$(cat "${WORK}/bwa-log.txt" | grep "Processed" | sed 's/[^0-9]*\([0-9]*\) reads in .* CPU sec, \([0-9]*\.[0-9]*\) real sec/\1 \2/g' | tr ' ' '\t' | awk '{sum1+=$1; sum2+=$2} END {print sum1/sum2}')"

BWA_RPS="$(echo "${BWA_RPS_ALL_THREADS} / ${THREAD_COUNT}" | bc -l)"

echo "==== Giraffe Wrangler Report for vg $(vg version -s) ===="

if [[ "${NUMA_WARNING}" == "1" ]] ; then
    echo "WARNING! Unable to restrict to a single NUMA node! Results may have high variance!"
fi

if which perf >/dev/null 2>&1 ; then
    # Output perf stuff
    mv "${WORK}/perf.data" ./perf.data
    mv "${WORK}/profile.svg" ./profile.svg
    echo "Profiling information saved as ./perf.data"
    echo "Interactive flame graph (for browsers) saved as ./profile.svg"
fi

# Print the report
echo "Giraffe got ${CORRECT_COUNT} simulated reads correct with ${MEAN_IDENTITY} average identity per mapped base"
echo "Map got ${CORRECT_COUNT_MAP} simulated reads correct with ${MEAN_IDENTITY_MAP} average identity per mapped base"
echo "Giraffe aligned real reads at ${GIRAFFE_RPS} reads/second vs. bwa-mem's ${BWA_RPS} reads/second on ${THREAD_COUNT} threads"

cat "${WORK}/facts.txt"

rm -Rf "${WORK}"




