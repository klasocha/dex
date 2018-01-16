#!/bin/bash
treat=(DEX CTRL)
gtf_dir="/home/seq/resources/mm10/mm10.tss.gtf"


#SalDex-ASTRO-invivo
for j in `seq 0 1`; do
	for k in `seq 1 5`; do
		input=/home/seq/bam/SalDex-ASTRO-invivo/A${k}_${treat[j]}.bam
		output=cufflinks/${treat[j]}-AST-${k}
		cufflinks $input --no-effective-length-correction -G $gtf_dir -o $output -p 20
	done
done



#SalDex-100nM-4h-NEU-invivo
treat_local=(DEX SAL)
for j in `seq 0 1`; do
	input=/home/seq/bam/SalDex-100nM-4h-NEU-invivo/N1${treat_local[j]}.bam
        output=cufflinks/${treat[j]}-NEU-1
	cufflinks $input --no-effective-length-correction -G $gtf_dir -o $output -p 20

        input=/home/seq/bam/SalDex-100nM-4h-NEU-invivo/N${treat_local[j]}.bam
        output=cufflinks/${treat[j]}-NEU-2
	cufflinks $input --no-effective-length-correction -G $gtf_dir -o $output -p 20

        input=/home/seq/bam/SalDex-100nM-4h-NEU-invivo/N3${treat_local[j]}.bam
        output=cufflinks/${treat[j]}-NEU-3
	cufflinks $input --no-effective-length-correction -G $gtf_dir -o $output -p 20

        input=/home/seq/bam/SalDex-100nM-4h-NEU-invivo/N4${treat_local[j]}.bam
        output=cufflinks/${treat[j]}-NEU-4
	cufflinks $input --no-effective-length-correction -G $gtf_dir -o $output -p 20
done


#Maestro-NAc-DEX-acute
treat_local=(dex ctrl)
for j in `seq 0 1`; do
        for k in `seq 1 4`; do
                input=/home/seq/bam/Maestro-NAc-DEX-acute/NAc-acute_${treat_local[j]}_${k}.bam
                output=cufflinks/${treat[j]}-NAC-${k}
		cufflinks $input --no-effective-length-correction -G $gtf_dir -o $output -p 20
        done
done


#slawgo-dex
treat_local=(DEX SAL)
for j in `seq 0 1`; do
        for k in `seq 1 5`; do
                input=/home/seq/bam/slawgo-dex/STR_${treat_local[j]}${k}.bam
                output=cufflinks/${treat[j]}-STR-${k}
		cufflinks $input --no-effective-length-correction -G $gtf_dir -o $output -p 20
        done
done

