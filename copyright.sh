#!/bin/bash
# add copyright notice to all source code
# adapted from Wouter Deconinck's one-liner

function get_ext { echo $1 | sed 's;.*\.;;g'; }

function get_authors {
  git blame -w -M -C -C --line-porcelain "$1" |\
    grep -I '^author ' |\
    sort -f |\
    uniq -ic |\
    sort -n --reverse |\
    awk 'NR==1 {ref=$1; $1=$2="";print $0} {if($1>0.1*ref){$1=$2="";print $0}}' |\
    sed 's/^\s*//g' |\
    sed ':x {N;s/\n/, /g; bx}' |\
    sed 's/chris/Chris/g' |\
    sed 's/dilks/Dilks/g' |\
    sed 's/cpecar/Connor Pecar/g' |\
    sed 's/bspage912/Brian Page/g' |\
    sed 's/sanghwapark/Sanghwa Park/g'
}

git ls-files |\
  grep -v 'deps/adage' |\
  grep -v 'deps/delphes_EIC' |\
  grep -v 'deps/mstwpdf' |\
  grep -v '\.keep' |\
  grep -v '\.gitignore' |\
  grep -v '\.md' |\
  grep -v 'Makefile' |\
  grep -v 'config.mk' |\
  grep -v '\.config' |\
  grep -v 'LICENSE' |\
  grep -v '\.odg' |\
  grep -v '\.png' |\
  grep -vE '^grids/' |\
  grep -v '\.dat' |\
  grep -v '\.tmp' |\
  while read f; do

    case `get_ext $f` in
      C|cxx|h|ipp) cmt='//'; ;;
      *) cmt='#'; ;;
    esac

    authors=`get_authors $f`
    echo $f $authors

    copyright="$cmt SPDX-License-Identifier: LGPL-3.0-or-later\n$cmt Copyright (C) 2022 $authors\n"
    
    case `get_ext $f` in
      C|cxx|h|ipp|yml)
        echo -e $copyright | cat - $f > $f.tmp
        ;;
      *)
        shebang=`head -n1 $f`
        tail -n+2 $f > $f.code
        echo -e "$shebang\n\n$copyright" | cat - $f.code > $f.tmp
        rm $f.code
        ;;
    esac
    mv $f.tmp $f
  done
