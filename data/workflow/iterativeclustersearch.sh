#!/bin/sh -e
# Iterative cluster search workflow script
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
[ "$#" -ne 4 ] && echo "Please provide <queryDB> <targetDB> <outDB> <tmp>" && exit 1;
# check if files exist
[ ! -f "$1.dbtype" ] && echo "$1.dbtype not found!" && exit 1;
[ ! -f "$2.dbtype" ] && echo "$2.dbtype not found!" && exit 1;
[   -f "$3.dbtype" ] && echo "$3.dbtype exists already!" && exit 1;
[ ! -d "$4" ] && echo "tmp directory $4 not found!" && mkdir -p "$4";

QUERY="$1"
#for accessing the metadata
QUERY_STEP_0="$1"
TARGET="$2"
TMP_PATH="$4"

STEP=0
# processing
[ -z "$NUM_IT" ] && NUM_IT=3;
while [ $STEP -lt $NUM_IT ]; do
    # call prefilter module
    if notExists "$TMP_PATH/pref_tmp_${STEP}.done"; then
        PARAM="PREFILTER_PAR_$STEP"
        eval TMP="\$$PARAM"
        if [ $STEP -eq 0 ]; then
            # shellcheck disable=SC2086
            $RUNNER "$MMSEQS" prefilter "${QUERY}" "${TARGET}" "$TMP_PATH/pref_$STEP" ${TMP} \
                || fail "Prefilter died"
        else
            # shellcheck disable=SC2086
            $RUNNER "$MMSEQS" prefilter "${QUERY}" "${TARGET}" "$TMP_PATH/pref_tmp_$STEP" ${TMP} \
                || fail "Prefilter died"
        fi
        touch "$TMP_PATH/pref_tmp_${STEP}.done"
    fi

    if [ $STEP -ge 1 ]; then
        if notExists "$TMP_PATH/pref_$STEP.done"; then
            STEPONE=$((STEP-1))
            # shellcheck disable=SC2086
            "$MMSEQS" subtractdbs "$TMP_PATH/pref_tmp_$STEP" "$TMP_PATH/aln_$STEPONE" "$TMP_PATH/pref_$STEP" $SUBSTRACT_PAR \
                || fail "Substract died"
            "$MMSEQS" rmdb "$TMP_PATH/pref_tmp_$STEP"
        fi
        touch "$TMP_PATH/pref_$STEP.done"
    fi

	# call alignment module
	if notExists "$TMP_PATH/aln_tmp_$STEP.done"; then
	    PARAM="ALIGNMENT_PAR_$STEP"
        eval TMP="\$$PARAM"
        if [ $STEP -eq 0 ]; then
            # shellcheck disable=SC2086
            $RUNNER "$MMSEQS" "${ALIGN_MODULE}" "${QUERY}" "${TARGET}" "$TMP_PATH/pref_$STEP" "$TMP_PATH/aln_$STEP" ${TMP} \
                || fail "Alignment died"
        else
            # shellcheck disable=SC2086
            $RUNNER "$MMSEQS" "${ALIGN_MODULE}" "${QUERY}" "${TARGET}" "$TMP_PATH/pref_$STEP" "$TMP_PATH/aln_tmp_$STEP" ${TMP} \
                || fail "Alignment died"
        fi
        touch "$TMP_PATH/aln_tmp_$STEP.done"
    fi



    if [ $STEP -gt 0 ]; then
        if notExists "$TMP_PATH/aln_$STEP.done"; then
            STEPONE=$((STEP-1))
            #TODO: check if the aln_tmp_step contains any new hits, if not, exit prematurely
            if [ -s "$TMP_PATH/aln_tmp_$STEP.index" ]; then 
                # shellcheck disable=SC2086
                "$MMSEQS" mergedbs "${QUERY_STEP_0}" "$TMP_PATH/aln_$STEP" "$TMP_PATH/aln_$STEPONE" "$TMP_PATH/aln_tmp_$STEP" \
                    || fail "mergedbs died"
                "$MMSEQS" rmdb "$TMP_PATH/aln_tmp_$STEP"
                "$MMSEQS" rmdb "$TMP_PATH/aln_$STEPONE"
                touch "$TMP_PATH/aln_$STEP.done"
            else
                echo "No new hits found in step $STEP, printing previous results"
                # shellcheck disable=SC2086
                "${MMSEQS}" tsv2db "${TMP_PATH}/cluster_sorted_$STEPONE" "$3" --output-dbtype 5 ${VERBOSITY} \
                || fail "tsv2db failed"
                break
            fi
        fi
    fi


    # clustersearch pipeline
    if notExists "$TMP_PATH/cluster_aln_$STEP.done"; then

        if notExists "${TMP_PATH}/aggregate_$STEP.index"; then
        # aggregation: take for each target set the best hit
        # shellcheck disable=SC2086
            "${MMSEQS}" besthitperset "${QUERY_STEP_0}" "${TARGET}" "${TMP_PATH}/aln_$STEP" "${TMP_PATH}/aggregate_$STEP" ${BESTHITBYSET_PAR} \
                || fail "aggregate best hit failed"
        fi

        if notExists "${TMP_PATH}/aggregate_merged_$STEP.index"; then
            # shellcheck disable=SC2086
            "${MMSEQS}" mergeresultsbyset "${QUERY_STEP_0}_set_to_member" "${TMP_PATH}/aggregate_$STEP" "${TMP_PATH}/aggregate_merged_$STEP" ${THREADS_PAR} \
                || fail "mergesetresults failed"
        fi


        if notExists "${TMP_PATH}/cEval_$STEP.index"; then
            # shellcheck disable=SC2086
            "${MMSEQS}" combinepvalperset "${QUERY_STEP_0}" "${TARGET}" "${TMP_PATH}/aggregate_merged_$STEP" "${TMP_PATH}/cEval_$STEP" "${TMP_PATH}" ${COMBINEPVALPERSET_PAR} \
                || fail "combinepvalperset failed"
        fi

        #TODO: check if this step is needed at all
        if notExists "${TMP_PATH}/match_$STEP.index"; then
            # shellcheck disable=SC2086
            #TODO: parameterize cEval_thr
            "${MMSEQS}" filterdb "${TMP_PATH}/cEval_$STEP" "${TMP_PATH}/match_$STEP" --filter-column "2" --comparison-operator "le" --comparison-value "0.01" ${THREADS_PAR} \
                || fail "filterdb failed"
        fi

        #TODO: check if comebinepvalperset can handle qid being in the first column, if so, add prefix to prior to aggregate_merged
        if notExists "${TMP_PATH}/aggregate_prefixed_$STEP.index"; then
            # shellcheck disable=SC2086
            "${MMSEQS}" prefixid  "${TMP_PATH}/aggregate_$STEP" "${TMP_PATH}/aggregate_prefixed_$STEP" ${THREADS_PAR} \
                || fail "prefixid failed"
        fi

        if notExists "${TMP_PATH}/aggregate_prefixed_merged_$STEP.index"; then
            # shellcheck disable=SC2086
            "${MMSEQS}" mergeresultsbyset "${QUERY_STEP_0}_set_to_member" "${TMP_PATH}/aggregate_prefixed_$STEP" "${TMP_PATH}/aggregate_prefixed_merged_$STEP" ${THREADS_PAR} \
                || fail "mergesetresults failed"
        fi

        if notExists "${TMP_PATH}/matches_$STEP.index"; then
            # shellcheck disable=SC2086
            "${MMSEQS}" filtermatches "${QUERY_STEP_0}" "${TARGET}" "${TMP_PATH}/aggregate_prefixed_merged_$STEP" "${TMP_PATH}/match_$STEP" "${TMP_PATH}/matches_$STEP" ${FILTERMATCHES_PAR} \
                || fail "filtermatches failed"
        fi

        #db-output set to false to not print the null bytes
        if notExists "${TMP_PATH}/cluster_$STEP.index"; then
            # shellcheck disable=SC2086
            "${MMSEQS}" clusterhits "${QUERY_STEP_0}" "${TARGET}" "${TMP_PATH}/matches_$STEP" "${TMP_PATH}/cluster_$STEP" ${CLUSTERHITS_PAR} \
                || fail "clusterhits failed"
        fi

        #TODO:sort by qid to group by qid?
        if notExists "${TMP_PATH}/cluster_sorted_$STEP.index"; then
            sort -k1,1n "${TMP_PATH}/cluster_$STEP" > "${TMP_PATH}/cluster_sorted_$STEP" \
            || fail "sort step $STEP died"
        fi

        if [ $STEP -ne $((NUM_IT  - 1)) ]; then
            if notExists "${TMP_PATH}/cluster_aln_$STEP.index"; then
                #TODO:--output-dbtype 5 which is the alignment dbtype
                # shellcheck disable=SC2086
                "${MMSEQS}" tsv2db "${TMP_PATH}/cluster_sorted_$STEP" "${TMP_PATH}/cluster_aln_$STEP" --output-dbtype 5 ${VERBOSITY} \
                    || fail "tsv2db failed"
            fi
        else
            # shellcheck disable=SC2086
            "${MMSEQS}" tsv2db "${TMP_PATH}/cluster_sorted_$STEP" "$3" --output-dbtype 5 ${VERBOSITY} \
                || fail "tsv2db failed"
        fi

        if [ -n "${REMOVE_TMP}" ]; then
            # shellcheck disable=SC2086
            "$MMSEQS" rmdb "${TMP_PATH}/aln_$STEP" ${VERBOSITY}
            # shellcheck disable=SC2086
            "$MMSEQS" rmdb "${TMP_PATH}/aggregate_$STEP" ${VERBOSITY}
            # shellcheck disable=SC2086
            "$MMSEQS" rmdb "${TMP_PATH}/aggregate_merged_$STEP" ${VERBOSITY}
            # shellcheck disable=SC2086
            "$MMSEQS" rmdb "${TMP_PATH}/cEval_$STEP" ${VERBOSITY}
            # shellcheck disable=SC2086
            "$MMSEQS" rmdb "${TMP_PATH}/match_$STEP" ${VERBOSITY}
            # shellcheck disable=SC2086
            "$MMSEQS" rmdb "${TMP_PATH}/aggregate_prefixed_$STEP" ${VERBOSITY}
            # shellcheck disable=SC2086
            "$MMSEQS" rmdb "${TMP_PATH}/aggregate_prefixed_merged_$STEP" ${VERBOSITY}
            # shellcheck disable=SC2086
            "$MMSEQS" rmdb "${TMP_PATH}/matches_$STEP" ${VERBOSITY}
            # shellcheck disable=SC2086
            "$MMSEQS" rmdb "${TMP_PATH}/cluster_$STEP" ${VERBOSITY}
        fi

    touch "$TMP_PATH/cluster_aln_$STEP.done"
    fi

# create profiles
    if [ $STEP -ne $((NUM_IT  - 1)) ]; then
        if notExists "$TMP_PATH/profile_$STEP.dbtype"; then
            PARAM="PROFILE_PAR_$STEP"
            eval TMP="\$$PARAM"
            # shellcheck disable=SC2086
            $RUNNER "$MMSEQS" result2profile "${QUERY}" "${TARGET}" "$TMP_PATH/cluster_aln_$STEP" "$TMP_PATH/profile_$STEP" ${TMP} \
                || fail "Create profile died"
        fi
    fi
	QUERY="$TMP_PATH/profile_$STEP"
	STEP=$((STEP+1))
done

if [ -n "$REMOVE_TMP" ]; then
    STEP=0
    while [ "$STEP" -lt "$NUM_IT" ]; do
        if [ $STEP -gt 0 ]; then
            rm -f -- "$TMP_PATH/aln_$STEP.done" "$TMP_PATH/pref_$STEP.done"
        fi
        rm -f -- "$TMP_PATH/aln_tmp_$STEP.done" "$TMP_PATH/pref_tmp_${STEP}.done"
        # shellcheck disable=SC2086
        "$MMSEQS" rmdb "${TMP_PATH}/pref_$STEP" ${VERBOSITY}
        # shellcheck disable=SC2086
        "$MMSEQS" rmdb "${TMP_PATH}/aln_$STEP" ${VERBOSITY}
        # shellcheck disable=SC2086
        "$MMSEQS" rmdb "${TMP_PATH}/profile_$STEP" ${VERBOSITY}
        # shellcheck disable=SC2086
        "$MMSEQS" rmdb "${TMP_PATH}/profile_${STEP}_h" ${VERBOSITY}
        STEP=$((STEP+1))
    done
    rm -f "$TMP_PATH/iterativevlustersearch.sh"
fi
