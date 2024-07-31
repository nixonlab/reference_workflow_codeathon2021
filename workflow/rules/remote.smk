import json
from snakemake import WorkflowError
# localrules: build_remotefile_db
# rule build_remotefile_db:
#     output:
#         'resources/remotefile_db.json'
#     input:
#         'config/config.yaml'
#     script:
#         '../scripts/build_remotefile_db.py'


def getparams_download_remote_db(wildcards, input, output):
    with open(input[0], 'r') as dbfh:
        db = json.load(dbfh)
    if wildcards.f in db:
        return dict(db[wildcards.f])
    else:
        raise WorkflowError(f"""
        Remote file '{wildcards.f}' not found in remotefile database. The database
        may need to be rebuilt. For example:
            $  workflow/scripts/build_remotefile_db.py config/remote.yaml > resources/remotefile_db.json
        Also ensure that the file is included in `config/remote.yaml`.
        """)

rule download_remotefile:
    """ Downloads a remote file and checks the md5sum.
        Filenames, URLs and md5 checksums are stored in 
        resources/remotefiles_db.txt
    """
    output:
        'databases/remotefiles/{f}'
    input:
        ancient(config['remotefile_db'])
    params:
        getparams_download_remote_db
    shell:
        '''
curl -L {params[0][url]} > {output[0]}
echo {params[0][md5]}  {output[0]} | md5sum -c -
        '''

localrules: download_genbank
rule download_genbank:
    output:
        'databases/remotefiles/genbank/{acc}.gbk'
    wildcard_constraints:
        acc = '([A-Z]\d{5}|[A-Z]{2}\d{6}|[A-Z]{2}\d{8}|[A-Z]{2}_\d{6,8})(\.\d{1,3})?'
    conda: '../envs/entrez.yaml'
    shell:
        '''
efetch -db nuccore -id {wildcards.acc} -format gb > {output[0]}
        '''
