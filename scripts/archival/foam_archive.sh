#!/bin/bash
# Archive completed OpenFOAM time directories to compressed tarballs.
# Works for serial, uncollated (processor*) and collated (processors*) output.
# Safe to run while the solver is writing: skips the newest KEEP times and
# anything modified in the last AGE minutes.

set -u

CASE=${CASE:-$PWD}
KEEP=${KEEP:-2}                 # newest time dirs left untouched for restart
AGE=${AGE:-20}                  # minutes of quiescence before a time is archived
ARCH=${ARCH:-$CASE/archive}
THREADS=${THREADS:-8}
LEVEL=${LEVEL:-3}               # zstd level, 3 fast, 12 slower and smaller
OFFLOAD=${OFFLOAD:-}            # optional target dir, e.g. $HOME/foam_archive
SCRATCH_PCT=${SCRATCH_PCT:-80}  # offload once scratch use exceeds this percent
DRYRUN=${DRYRUN:-0}

mkdir -p "$ARCH"
LOG=$ARCH/archive.log
LOCK=$ARCH/archive.lock

log() { echo "$(date '+%F %T') $*" >> "$LOG"; }

# single instance
exec 9>"$LOCK"
flock -n 9 || { log "another instance running, exit"; exit 0; }

# pick compressor
if command -v zstd >/dev/null 2>&1; then
    COMP="zstd -T$THREADS -$LEVEL -q"
    TEST="zstd -t -q"
    EXT=tar.zst
elif command -v pigz >/dev/null 2>&1; then
    COMP="pigz -p $THREADS -6"
    TEST="pigz -t"
    EXT=tar.gz
else
    COMP="gzip -6"
    TEST="gzip -t"
    EXT=tar.gz
fi

cd "$CASE" || { log "case dir missing"; exit 1; }

# detect output layout
if compgen -G "processors[0-9]*" >/dev/null; then
    MODE=collated
    SRC=$(ls -d processors[0-9]* | head -1)
elif [ -d processor0 ]; then
    MODE=uncollated
    SRC=processor0
else
    MODE=serial
    SRC=.
fi
log "mode $MODE source $SRC compressor $EXT"

# list numeric time dirs, ascending
times=$(find "$SRC" -maxdepth 1 -mindepth 1 -type d -printf '%f\n' \
        | grep -E '^[0-9]+(\.[0-9]+)?([eE][-+]?[0-9]+)?$' \
        | grep -vE '^0(\.0+)?$' \
        | sort -g)

n=$(echo "$times" | grep -c . )
if [ "$n" -le "$KEEP" ]; then
    log "nothing eligible, $n times present"
    exit 0
fi

# drop the newest KEEP
eligible=$(echo "$times" | head -n $((n - KEEP)))

for t in $eligible; do
    out=$ARCH/time_$t.$EXT
    [ -f "$out" ] && continue

    case $MODE in
        collated)   paths=$(ls -d processors[0-9]*/"$t" 2>/dev/null) ;;
        uncollated) paths=$(ls -d processor*/"$t" 2>/dev/null) ;;
        serial)     paths=$t ;;
    esac
    [ -z "$paths" ] && continue

    # still being written
    busy=0
    for p in $paths; do
        if [ -n "$(find "$p" -newermt "-$AGE minutes" -print -quit 2>/dev/null)" ]; then
            busy=1; break
        fi
    done
    [ "$busy" -eq 1 ] && { log "time $t still active, skip"; continue; }

    raw=$(du -sbc $paths 2>/dev/null | tail -1 | cut -f1)

    if [ "$DRYRUN" -eq 1 ]; then
        log "dryrun would archive time $t raw bytes $raw"
        continue
    fi

    if ! tar -cf - $paths 2>>"$LOG" | $COMP > "$out.part" 2>>"$LOG"; then
        log "tar failed for time $t, archive kept as part file"
        rm -f "$out.part"
        continue
    fi

    if ! $TEST "$out.part" 2>>"$LOG"; then
        log "verify failed for time $t, originals kept"
        rm -f "$out.part"
        continue
    fi

    mv "$out.part" "$out"
    comp=$(stat -c %s "$out")
    ratio=$(awk -v a="$raw" -v b="$comp" 'BEGIN{printf "%.2f", a/b}')
    rm -rf $paths
    log "archived time $t raw $raw compressed $comp ratio $ratio"
done

# offload when scratch is filling
if [ -n "$OFFLOAD" ]; then
    use=$(df --output=pcent "$CASE" | tail -1 | tr -dc '0-9')
    if [ "${use:-0}" -ge "$SCRATCH_PCT" ]; then
        mkdir -p "$OFFLOAD"
        moved=0
        for f in $(ls -tr "$ARCH"/time_*."$EXT" 2>/dev/null); do
            if rsync -a --remove-source-files "$f" "$OFFLOAD"/ 2>>"$LOG"; then
                moved=$((moved + 1))
                use=$(df --output=pcent "$CASE" | tail -1 | tr -dc '0-9')
                [ "${use:-0}" -lt "$SCRATCH_PCT" ] && break
            else
                log "offload failed for $(basename "$f")"
                break
            fi
        done
        log "offload moved $moved archives, scratch use now $use percent"
    fi
fi

log "pass complete"
