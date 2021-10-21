#!/bin/sh -e
fail() {
    echo "Error: $1"
    exit 1
}

notExists() {
	[ ! -f "$1" ]
}

#pre processing
[ -z "$MMSEQS" ] && echo "Please set the environment variable \$MMSEQS to your MMSEQS binary." && exit 1;
# check number of input variables
[ "$#" -ne 4 ] && echo "Please provide <queryDB> <targetDB> <outputDB> <tmpDir>" && exit 1;
# check if files exist
[ ! -f "$1.dbtype" ] && echo "$1.dbtype not found!" && exit 1;
[ ! -f "$2.dbtype" ] && echo "$2.dbtype not found!" && exit 1;
# TO DO??? add check if $3.dbtype already exists before entire workfolw ???

QUERY="$1"
TARGET="$2"
OUTPUT="$3"
TMP_PATH="$4"

#clustering before search?

if notExists "${TMP_PATH}/result.index"; then
    # shellcheck disable=SC2086
    "${MMSEQS}" search "${QUERY}" "${TARGET}" "${TMP_PATH}/result" "${TMP_PATH}/search" ${SEARCH_PAR} \
        || fail "search failed"
fi

#filter self hits? a lot of self hits could overwhlem the results produced in the search module
if notExists "${TMP_PATH}/result_filtered.index"; then
    # shellcheck disable=SC2086
    #TODO: parameterize cEval_thr
    "${MMSEQS}" filterdb "${TMP_PATH}/result" "${TMP_PATH}/result_filtered" --filter-column "3" --comparison-operator "le" --comparison-value "0.99" ${THREADS_PAR} \
        || fail "filterdb failed"
fi

if notExists "${TMP_PATH}/aggregate.index"; then
    # aggregation: take for each target set the best hit
    # shellcheck disable=SC2086
    "${MMSEQS}" besthitperset "${QUERY}" "${TARGET}" "${TMP_PATH}/result_filtered" "${TMP_PATH}/aggregate" ${BESTHITBYSET_PAR} \
        || fail "aggregate best hit failed"
fi

if notExists "${TMP_PATH}/aggregate_merged.index"; then
    # shellcheck disable=SC2086
    "${MMSEQS}" mergeresultsbyset "${QUERY}_set_to_member" "${TMP_PATH}/aggregate" "${TMP_PATH}/aggregate_merged" ${THREADS_PAR} \
        || fail "mergesetresults failed"
fi


if notExists "${TMP_PATH}/cEval.index"; then
    # shellcheck disable=SC2086
    "${MMSEQS}" combinepvalperset "${QUERY}" "${TARGET}" "${TMP_PATH}/aggregate_merged" "${TMP_PATH}/cEval" "${TMP_PATH}" ${COMBINEPVALPERSET_PAR} \
        || fail "combinepvalperset failed"
fi


if notExists "${TMP_PATH}/match.index"; then
    # shellcheck disable=SC2086
    #TODO: parameterize cEval_thr
    "${MMSEQS}" filterdb "${TMP_PATH}/cEval" "${TMP_PATH}/match" --filter-column "2" --comparison-operator "le" --comparison-value "0.01" ${THREADS_PAR} \
        || fail "filterdb failed"
fi

#TODO: check if comebinepvalperset can handle qid being in the first column, if so, add prefix to prior to aggregate_merged
if notExists "${TMP_PATH}/aggregate_prefixed.index"; then
    # shellcheck disable=SC2086
    "${MMSEQS}" prefixid  "${TMP_PATH}/aggregate" "${TMP_PATH}/aggregate_prefixed" ${THREADS_PAR} \
        || fail "prefixid failed"
fi

if notExists "${TMP_PATH}/aggregate_prefixed_merged.index"; then
    # shellcheck disable=SC2086
    "${MMSEQS}" mergeresultsbyset "${QUERY}_set_to_member" "${TMP_PATH}/aggregate_prefixed" "${TMP_PATH}/aggregate_prefixed_merged" ${THREADS_PAR} \
        || fail "mergesetresults failed"
fi

##provisional##
if notExists "${TMP_PATH}/matches.index"; then
    # shellcheck disable=SC2086
    "${MMSEQS}" filtermatches "${QUERY}" "${TARGET}" "${TMP_PATH}/aggregate_prefixed_merged" "${TMP_PATH}/match" "${TMP_PATH}/matches" ${THREADS_PAR} \
        || fail "filtermatches failed"
fi

if notExists "${TMP_PATH}/clusters.index"; then
    # shellcheck disable=SC2086
    "${MMSEQS}" clusterhits "${QUERY}" "${TARGET}" "${TMP_PATH}/matches" "${OUTPUT}" ${THREADS_PAR} \
        || fail "clusterhits failed"
fi

if [ -n "${REMOVE_TMP}" ]; then
    rmdir "${TMP_PATH}/search"
    # shellcheck disable=SC2086
    "$MMSEQS" rmdb "${TMP_PATH}/result" ${VERBOSITY}
    # shellcheck disable=SC2086
    "$MMSEQS" rmdb "${TMP_PATH}/aggregate" ${VERBOSITY}
    rm -f "${TMP_PATH}/multihitsearch.sh"
fi

