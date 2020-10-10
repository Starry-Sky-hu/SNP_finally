#!/bin/bash
#PBS -N sh_SNPtoTable_s
#PBS -o /public/agis/huangsanwen_group/fengshuangshuang/huyong/source/CallSNP/zzz-qsub.history/qsub.sh_SNPtoTable_s.24754.log
#PBS -e /public/agis/huangsanwen_group/fengshuangshuang/huyong/source/CallSNP/zzz-qsub.history/qsub.sh_SNPtoTable_s.24754.err
#PBS -q low
# qsub parameter: "-q queue1 -l nodes=1:ppn=3"
uname -a
cd /public/agis/huangsanwen_group/fengshuangshuang/huyong/source/CallSNP
date
echo "job: sh_SNPtoTable_s"
echo $'sh SNPtoTable.sh'
sh SNPtoTable.sh
date
