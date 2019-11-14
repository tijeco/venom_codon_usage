fastq-dump --defline-seq '@$sn[_$rn]/$ri' --split-files SRR1964035
mv SRR1964035_1.fq Conus_consors_venom_1.fq
mv SRR1964035_2.fq Conus_consors_venom_2.fq

fastq-dump --defline-seq '@$sn[_$rn]/$ri' --split-files SRR1958882
mv SRR1958882_1.fq Conus_consors_body_1.fq
mv SRR1958882_2.fq Conus_consors_body_2.fq

fastq-dump --defline-seq '@$sn[_$rn]/$ri' --split-files SRR5185496
mv SRR5185496_1.fq Eutolmus_rufibarbis_venom_1.fq
mv SRR5185496_2.fq Eutolmus_rufibarbis_venom_2.fq

fastq-dump --defline-seq '@$sn[_$rn]/$ri' --split-files SRR5185497
mv SRR5185497_1.fq Eutolmus_rufibarbis_body_1.fq
mv SRR5185497_2.fq Eutolmus_rufibarbis_body_2.fq

fastq-dump --defline-seq '@$sn[_$rn]/$ri' --split-files SRR5185498
mv SRR5185498_1.fq Machimus_arthriticus_venom_1.fq
mv SRR5185498_2.fq Machimus_arthriticus_venom_2.fq

fastq-dump --defline-seq '@$sn[_$rn]/$ri' --split-files SRR5185499
mv SRR5185499_1.fq Machimus_arthriticus_body_1.fq
mv SRR5185499_2.fq Machimus_arthriticus_body_2.fq

fastq-dump --defline-seq '@$sn[_$rn]/$ri' --split-files SRR2592960
fastq-dump --defline-seq '@$sn[_$rn]/$ri' --split-files SRR3061379
fastq-dump --defline-seq '@$sn[_$rn]/$ri' --split-files SRR3984663
cat SRR2592960_1.fq SRR3061379_1.fq SRR3984663_1.fq > Mesobuthus_Martensii_venom_1.fq
cat SRR2592960_2.fq SRR3061379_2.fq SRR3984663_2.fq > Mesobuthus_Martensii_venom_2.fq




fastq-dump --defline-seq '@$sn[_$rn]/$ri' --split-files SRR3984597
fastq-dump --defline-seq '@$sn[_$rn]/$ri' --split-files SRR2592319
fastq-dump --defline-seq '@$sn[_$rn]/$ri' --split-files SRR3056832
cat SRR3984597_1.fq SRR2592319_1.fq SRR3056832_1.fq > Mesobuthus_martensii_body_1.fq
cat SRR3984597_2.fq SRR2592319_2.fq SRR3056832_2.fq > Mesobuthus_martensii_body_2.fq


fastq-dump --defline-seq '@$sn[_$rn]/$ri' --split-files SRR5285122
mv SRR5285122_1.fq Latrodectus_Hesperus_venom_1.fq
mv SRR5285122_2.fq Latrodectus_Hesperus_venom_2.fq

fastq-dump --defline-seq '@$sn[_$rn]/$ri' --split-files SRR5285105
mv SRR5285105_1.fq Latrodectus_Hesperus_body_1.fq
mv SRR5285105_2.fq Latrodectus_Hesperus_body_2.fq

fastq-dump --defline-seq '@$sn[_$rn]/$ri' --split-files SRR5285141
mv SRR5285141_1.fq Steatoda_Grossa_venom_1.fq
mv SRR5285141_2.fq Steatoda_Grossa_venom_2.fq

fastq-dump --defline-seq '@$sn[_$rn]/$ri' --split-files SRR5285126
mv SRR5285126_1.fq Steatoda_Grossa_body_1.fq
mv SRR5285126_2.fq Steatoda_Grossa_body_2.fq
