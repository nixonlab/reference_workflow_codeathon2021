# -*- coding: utf-8 -*-

rule extract_cellranger_ref:
    input:
        'databases/remotefiles/refdata-{refver}.tar.gz'
    output:
        'databases/indexes/cellranger/refdata-{refver}/reference.json',
        'databases/indexes/cellranger/refdata-{refver}/fasta/genome.fa',
        'databases/indexes/cellranger/refdata-{refver}/fasta/genome.fa.fai',
        'databases/indexes/cellranger/refdata-{refver}/genes/genes.gtf',
        'databases/indexes/cellranger/refdata-{refver}/star/SAindex'
    wildcard_constraints:
        refver = f'({"|".join(config["cellranger_versions"])})'
    shell:
        '''
mkdir -p databases/indexes/cellranger
tar -xzvf {input[0]} -C databases/indexes/cellranger
        '''


rule star_index_GRCh38_2020A_F4HSA_gencode32:
    input:
        genome_fna = 'databases/indexes/cellranger/refdata-gex-GRCh38-2020-A/fasta/genome.fa',
        hiv_fna = 'databases/localfiles/HIV1.H4-FSA.fna',
        genome_gtf = 'databases/indexes/cellranger/refdata-gex-GRCh38-2020-A/genes/genes.gtf'
    output:
        protected(directory("databases/indexes/STAR_GRCh38_2020A_F4HSA_gencode32")),
        expand('databases/indexes/STAR_GRCh38_2020A_F4HSA_gencode32/{basename}',
            basename = star_index_files
        )
    params:
        sjdbOverhang=100,
        genomeDir = lambda wildcards, output: os.path.dirname(output[1])
    conda: "../envs/star.yaml"
    threads: snakemake.utils.available_cpu_count()
    shell:
        '''
tdir=$(mktemp -d {config[local_tmp]}/{rule}.XXXXXX)
echo "temporary directory: $tdir"
cat {input.genome_fna} {input.hiv_fna} > $tdir/genome.fna
cat {input.genome_gtf} > $tdir/transcripts.gtf

mkdir -p $tdir/starindex
STAR\
 --runThreadN {threads}\
 --runMode genomeGenerate\
 --genomeDir $tdir/starindex\
 --genomeFastaFiles $tdir/genome.fna\
 --sjdbGTFfile $tdir/transcripts.gtf\
 --sjdbOverhang {params.sjdbOverhang}

mkdir -p {output[0]}
cp $tdir/starindex/* {output[0]}
        '''

rule cellranger_annotation_gtf:
    input:
        allgtf = 'databases/indexes/cellranger/refdata-{refver}/genes/genes.gtf'
    output:
        allgtf = 'databases/annotations/cellranger.{refver}/genes.gtf.gz',
        alltbi = 'databases/annotations/cellranger.{refver}/genes.gtf.gz.tbi',
        chroms = 'databases/annotations/cellranger.{refver}/chroms.txt'
    wildcard_constraints:
        refver = f'({"|".join(config["cellranger_versions"])})'
    conda: '../envs/utils.yaml'
    shell:
        '''
mkdir -p $(dirname {output.chroms})
cat {input.allgtf} | grep -v '^#' | cut -f1 | uniq > {output.chroms}

(
    set +o pipefail
    cat {input.allgtf} | head -n1000 | grep "^#"
    echo "## PG: bedtools sort -g {output.chroms} -i -"
    echo "## PG: bgzip"
    cat {input.allgtf} | grep -v "^#" | bedtools sort -g {output.chroms} -i -
) | bgzip -c > {output.allgtf}
tabix -p gff {output.allgtf}
        '''

rule cellranger_annotation_metadata:
    input:
        allgtf = rules.cellranger_annotation_gtf.output.allgtf
    output:
        expand('databases/annotations/cellranger.{{refver}}/metadata.{table}.txt.gz',
            table = ['gene_features', 'tx_features',
                     'gid_gname', 'gid_gtype', 'gid_hgnc', 'gid_tid',
                     'tid_tname', 'tid_ttype', 'tid_hgnc', 'tid_gid',
                    ]
        )
    wildcard_constraints:
        refver = f'({"|".join(config["cellranger_versions"])})'
    params:
        out_prefix = lambda wc: f'databases/annotations/cellranger.{wc.refver}/'
    shell:
        '''
workflow/scripts/gtf_metadata.py {input.allgtf} {params.out_prefix}
        '''

localrules: cellranger_annotation_rds

rule cellranger_annotation_rds:
    input:
        rules.cellranger_annotation_metadata.output
    output:
        expand('databases/annotations/cellranger.{{refver}}/metadata.{table}.rds',
            table=['gene_features', 'tx_features',
                   'gid_gname', 'gid_gtype', 'gid_hgnc', 'gid_tid',
                   'tid_tname', 'tid_ttype', 'tid_hgnc', 'tid_gid',
                   ]
        )
    wildcard_constraints:
        refver = f'({"|".join(config["cellranger_versions"])})'
    conda: '../envs/rbase.yaml'
    script:
        '../scripts/metadata_rds.R'


rule star_index_GRCh38_2020A_HIVADA_gencode32:
    input:
        genome_fna = [
            'databases/indexes/cellranger/refdata-gex-GRCh38-2020-A/fasta/genome.fa',
            'databases/sequences/genbank/AF004394.1.fna',
            'databases/sequences/genbank/AY245535.1.fna',
        ],
        genome_gtf = 'databases/indexes/cellranger/refdata-gex-GRCh38-2020-A/genes/genes.gtf'
    output:
        protected(directory("databases/indexes/STAR_GRCh38_2020A_HIVADA_gencode32")),
        expand('databases/indexes/STAR_GRCh38_2020A_HIVADA_gencode32/{basename}',
            basename = star_index_files
        )
    params:
        sjdbOverhang=100,
        genomeDir = lambda wildcards, output: os.path.dirname(output[1])
    conda: "../envs/star.yaml"
    threads: snakemake.utils.available_cpu_count()
    shell:
        '''
tdir=$(mktemp -d {config[local_tmp]}/{rule}.XXXXXX)
echo "temporary directory: $tdir"
cat {input.genome_fna} > $tdir/genome.fna
cat {input.genome_gtf} > $tdir/transcripts.gtf

mkdir -p $tdir/starindex
STAR\
 --runThreadN {threads}\
 --runMode genomeGenerate\
 --genomeDir $tdir/starindex\
 --genomeFastaFiles $tdir/genome.fna\
 --sjdbGTFfile $tdir/transcripts.gtf\
 --sjdbOverhang {params.sjdbOverhang}

mkdir -p {output[0]}
cp $tdir/starindex/* {output[0]}
        '''
