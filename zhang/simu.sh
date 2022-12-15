for j in 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1
do
python zhang/data/MTX_synthetic/synth_mgx_mtx.py --spike-exp 0.1 --spike-exp-strength $j zhang/data/MTX_synthetic/hmp1-ii_metaphlan2-mtd-qcd.pcl Stool zhang/data/simu/true-exp_$j
python zhang/data/MTX_synthetic/synth_mgx_mtx.py --spike-dep --spike-exp 0.1 --spike-exp-strength $j zhang/data/MTX_synthetic/hmp1-ii_metaphlan2-mtd-qcd.pcl Stool zhang/data/simu/true-combo-dep-exp_$j
python zhang/data/MTX_synthetic/synth_mgx_mtx.py --spike-bug 0.5 --spike-exp 0.1 --spike-exp-strength $j zhang/data/MTX_synthetic/hmp1-ii_metaphlan2-mtd-qcd.pcl Stool zhang/data/simu/true-combo-bug-exp_$j
python zhang/data/MTX_synthetic/synth_mgx_mtx.py --spike-groups --spike-exp 0.1 --spike-exp-strength $j zhang/data/MTX_synthetic/hmp1-ii_metaphlan2-mtd-qcd.pcl Stool zhang/data/simu/group-true-exp_$j
done


