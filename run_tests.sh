#!/bin/bash


let num_failed=0
function must_fail {
    cmd=$*
    STATUS=PASSED
    bash -c "$cmd" 2> /dev/null
    if [ 0 -eq $? ];  then	
	STATUS=FAILED
	let num_failed=num_failed+1
    fi
    echo $STATUS $cmd
}

function must_succeed {
    cmd=$*
    STATUS=PASSED
    bash -c "$cmd" 2> /dev/null
    if [ 0 -ne $? ];  then	
	STATUS=FAILED
	bash -c "$cmd"
	let num_failed=num_failed+1
    fi
    echo $STATUS $cmd
}

echo "Running a few basic tests..."

must_fail "./bin/sv_genefusions_overlap.R example/id1.sv.bedpe example/id1.sum.tsv.gz 100000000000 n "
must_fail "./bin/sv_genefusions_overlap.R example/id1.sv.bedpe example/id1.sum.tsv.gz 100000000000"
must_fail "./bin/sv_genefusions_overlap.R example/id1.sv.bedpe example/id1.sum.tsv.gz"
must_fail "./bin/sv_genefusions_overlap.R example/id1.sv.bedpe example/id1.sum.tsv.gz aaaa"
must_fail "./bin/sv_genefusions_overlap.R example/id1.sv.bedpe"
must_fail "./bin/sv_genefusions_overlap.R"

must_succeed "./bin/sv_genefusions_overlap.R example/id1.sv.bedpe example/id1.sum.tsv.gz 100000000000 n out_file1"
must_succeed "./bin/sv_genefusions_overlap.R example/id1.sv.bedpe example/id1.sum.tsv.gz 100000000000 y out_file1"

echo Failed tests: $num_failed
exit $num_failed
