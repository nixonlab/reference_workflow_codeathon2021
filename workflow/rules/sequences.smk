# -*- coding: utf-8 -*-

rule setup_ncbi_analysis_set:
    input:
        'databases/remotefiles/GCA_000001405.15_GRCh38_{aset}_analysis_set.fna.gz'
    output:
        'databases/sequences/GCA_000001405.15_GRCh38_{aset}_analysis_set/genome.fna.gz',
        'databases/sequences/GCA_000001405.15_GRCh38_{aset}_analysis_set/genome.fna.gz.gzi',
        'databases/sequences/GCA_000001405.15_GRCh38_{aset}_analysis_set/genome.fna.gz.fai',
        'databases/sequences/GCA_000001405.15_GRCh38_{aset}_analysis_set/genome.dict'
    conda:
        '../envs/utils.yaml'
    shell:
        '''
mkdir -p $(dirname {output[0]})        
gunzip -c {input} | bgzip > {output[0]}
samtools faidx {output[0]}
picard CreateSequenceDictionary -R {output[0]}
        '''


rule setup_gdc_analysis_set:
    input:
        'databases/remotefiles/GRCh38.d1.vd1.fa.tar.gz'
    output:
        'databases/sequences/GRCh38.d1.vd1/genome.fna.gz',
        'databases/sequences/GRCh38.d1.vd1/genome.fna.gz.gzi',
        'databases/sequences/GRCh38.d1.vd1/genome.fna.gz.fai',
        'databases/sequences/GRCh38.d1.vd1/genome.dict'
    conda:
        '../envs/utils.yaml'
    shell:
        '''
mkdir -p $(dirname {output[0]})
tar -Oxzf {input} | bgzip > {output[0]}
samtools faidx {output[0]}
picard CreateSequenceDictionary -R {output[0]}
        '''

rule setup_T2T_analysis_set:
    input:
        'databases/remotefiles/chm13v2.0{aset}.fa.gz'
    output:
        'databases/sequences/chm13v2.0{aset}/genome.fna.gz',
        'databases/sequences/chm13v2.0{aset}/genome.fna.gz.gzi',
        'databases/sequences/chm13v2.0{aset}/genome.fna.gz.fai',
        'databases/sequences/chm13v2.0{aset}/genome.dict'
    wildcard_constraints:
        aset = '(_maskedY|_maskedY_rCRS|_noY)*'
    conda:
        '../envs/utils.yaml'
    shell:
        '''
mkdir -p $(dirname {output[0]})
gunzip -c {input} | bgzip > {output[0]}
samtools faidx {output[0]}
picard CreateSequenceDictionary -R {output[0]}
        '''

localrules: genbank_fasta
rule genbank_fasta:
    output:
        'databases/sequences/genbank/{acc}.fna'
    input:
        rules.download_genbank.output[0]
    wildcard_constraints:
        acc = '([A-Z]\d{5}|[A-Z]{2}\d{6}|[A-Z]{2}\d{8}|[A-Z]{2}_\d{6,8})(\.\d{1,3})?'
    conda: '../envs/biopython.yaml'
    shell:
        '''
workflow/scripts/gbk2fna.py {input[0]} {output[0]}
        '''
