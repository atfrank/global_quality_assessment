#!/bin/bash
rnas="1XHP 1YSV 1Z2J 1ZC5 28SR 2FDT 2JWV 2K66 2KOC 2L1V 2L3E 2LBJ 2LBL 2LDL 2LDT 2LHP 2LI4 2LK3 2LP9 2LPA 2LPS 2LQZ 2LU0 2LUB 2LUN 2LV0 2M12 2M21 2M22 2M24 2M4W 2M5U 2M8K 2MEQ 2MHI 2MNC 2MXL 2N2O 2N2P 2N4L 2N6S 2N6T 2N6W 2N6X 2NCI 2QH2 2QH4 2RVO 2Y95 4A4S 4A4T 4A4U 5A17 5A18 5IEM 5KH8 5KMZ 5KQE 5LSN 5LWJ 5N5C 5UF3 5UZT 5V16 5V17 5WQ1 6EZ0"

paste <(echo "rna") <(echo "recall0") <(echo "recall1") <(echo "f1_0") <(echo "f1_1") <(echo "acc") --delimiters ' ' >> nn_loo_acc.txt
for rna in $rnas
do
  # get recall 0
  recall0="$(tail -n7 log/nn_log_${rna}.txt | head -n1 | awk '{print $3}')"
  # get recall 1
  recall1="$(tail -n6 log/nn_log_${rna}.txt | head -n1 | awk '{print $3}')"
  # get f1 0
  f1_0="$(tail -n7 log/nn_log_${rna}.txt | head -n1 | awk '{print $4}')"
  # get f1 1
  f1_1="$(tail -n6 log/nn_log_${rna}.txt | head -n1 | awk '{print $4}')"
  # get accuracy
  acc="$(tail -n4 log/nn_log_${rna}.txt | head -n1 | awk '{print $5}')"
  # paste together
  paste <(echo "$rna") <(echo "$recall0") <(echo "$recall1") <(echo "$f1_0") <(echo "$f1_1") <(echo "$acc") --delimiters ' ' >> nn_loo_acc.txt
done
