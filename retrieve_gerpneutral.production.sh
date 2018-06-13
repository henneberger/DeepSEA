#retrive PhastCons, PhyloP, and GERP scores
sed 's/chr//g' $1 |grep -v NA|awk '{printf $1 "\t" $2-1 "\t" $2 "\n"}'|sort -u  -k1,1 -k2,2n|cut -f 1,3|tabix  ./resources/phastCons/primates_nohuman.tsv.gz  -R - >$1.evo1
sed 's/chr//g' $1 |grep -v NA|awk '{printf $1 "\t" $2-1 "\t" $2 "\n"}'|sort -u -k1,1 -k2,2n|cut -f 1,3|tabix  ./resources/phyloP/primates_nohuman.tsv.gz  -R - >$1.evo2
sed 's/chr//g' $1 |grep -v NA|awk '{printf $1 "\t" $2-1 "\t" $2 "\n"}'|sort -u -k1,1 -k2,2n|cut -f 1,3|tabix  ./resources/Gerp/gerp_scores.tsv.gz  -R - >$1.evo3
