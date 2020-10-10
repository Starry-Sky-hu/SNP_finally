#!/bin/bash
#PBS -N sh_ParentSNP_sh
#PBS -o /public/agis/huangsanwen_group/fengshuangshuang/huyong/source/CallSNP/zzz-qsub.history/qsub.sh_ParentSNP_sh.11304.log
#PBS -e /public/agis/huangsanwen_group/fengshuangshuang/huyong/source/CallSNP/zzz-qsub.history/qsub.sh_ParentSNP_sh.11304.err
#PBS -q low
# qsub parameter: "-q queue7 -l nodes=1:ppn=5"
uname -a
cd /public/agis/huangsanwen_group/fengshuangshuang/huyong/source/CallSNP
date
echo "job: sh_ParentSNP_sh"
echo $'sh ParentSNP.sh'
sh ParentSNP.sh
date
