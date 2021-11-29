#!/bin/sh -e
fail() {
    echo "Error: $1"
    exit 1
}

notExists() {
	[ ! -f "$1" ]
}

abspath() {
    if [ -d "$1" ]; then
        (cd "$1"; pwd)
    elif [ -f "$1" ]; then
        if [ -z "${1##*/*}" ]; then
            echo "$(cd "${1%/*}"; pwd)/${1##*/}"
        else
            echo "$(pwd)/$1"
        fi
    elif [ -d "$(dirname "$1")" ]; then
        echo "$(cd "$(dirname "$1")"; pwd)/$(basename "$1")"
    fi
}

hasCommand () {
    command -v "$1" >/dev/null 2>&1 || { echo "Please make sure that $1 is in \$PATH."; exit 1; }
}

hasCommand awk
hasCommand sort

[ -z "$MMSEQS" ] && echo "Please set the environment variable \$MMSEQS to your MMSEQS binary." && exit 1;

export MMSEQS_FORCE_MERGE=1

OUTDB="$(abspath "${OUTDB}")"

if notExists "${OUTDB}"; then
    # shellcheck disable=SC2086
    "${MMSEQS}" createdb "$@" "${TMP_PATH}/seqDB" ${CREATEDB_PAR} \
        || fail "createdb failed"
fi

if notExists "${OUTDB}"; then
    # shellcheck disable=SC2086
    "${MMSEQS}" gff2db $(cat "${GFFDIR}") "${TMP_PATH}/seqDB" "${OUTDB}" ${GFF2DB_PAR} \
        || fail "gff2db failed"
fi


if [ "$("${MMSEQS}" dbtype "${OUTDB}")" = "Nucleotide" ]; then
    mv -f "${OUTDB}" "${OUTDB}_nucl"
    mv -f "${OUTDB}.index" "${OUTDB}_nucl.index"
#    mv -f "${OUTDB}.lookup" "${OUTDB}_nucl.lookup"
#    mv -f "${OUTDB}.source" "${OUTDB}_nucl.source"
    mv -f "${OUTDB}.dbtype" "${OUTDB}_nucl.dbtype"

#    mv -f "${OUTDB}_h" "${OUTDB}_nucl_h"
#    mv -f "${OUTDB}_h.index" "${OUTDB}_nucl_h.index"
#    mv -f "${OUTDB}_h.dbtype" "${OUTDB}_nucl_h.dbtype"

    if notExists "${OUTDB}_nucl_contig_to_set.index"; then
        awk '{ print $1"\t"$3; }' "${OUTDB}.lookup" | sort -k1,1n -k2,2n > "${OUTDB}_member_to_set.tsv"
        "${MMSEQS}" tsv2db "${OUTDB}_member_to_set.tsv" "${OUTDB}_member_to_set" --output-dbtype 5 \
            || fail "tsv2db failed"
    fi

    if notExists "${OUTDB}_nucl_set_to_contig.index"; then
        awk '{ print $3"\t"$1; }' "${OUTDB}.lookup" | sort -k1,1n -k2,2n > "${OUTDB}_set_to_member.tsv"
        "${MMSEQS}" tsv2db "${OUTDB}_set_to_member.tsv" "${OUTDB}_set_to_member" --output-dbtype 5 \
            || fail "tsv2db failed"
    fi

    if notExists "${OUTDB}.index"; then
        # shellcheck disable=SC2086
        "${MMSEQS}" translatenucs "${OUTDB}_nucl" "${OUTDB}" ${TRANSLATENUCS_PAR} \
            || fail "translatenucs failed"
    fi

    if notExists "${OUTDB}_set_size.index"; then
        # shellcheck disable=SC2086
        "${MMSEQS}" result2stats "${OUTDB}" "${OUTDB}" "${OUTDB}_set_to_member" "${OUTDB}_set_size" ${RESULT2STATS_PAR} \
            || fail "result2stats failed"
    fi

else
    fail "protein mode not implemented"
fi

if [ -n "${REMOVE_TMP}" ]; then
    echo "Remove temporary files"
    rm -f "${TMP_PATH}/multihitdb.sh"
fi
