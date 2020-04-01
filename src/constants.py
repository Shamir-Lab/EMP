import json
import os


SCORES_FILE_NANE="scores.tsv"
ENSEMBL_TO_GENE_SYMBOLS = "ensembl2gene_symbol.txt"
ENSEMBL_TO_ENTREZ = "ensembl2entrez.txt"

dir_path = os.path.join(os.path.dirname(os.path.realpath(__file__)),"../")
PATH_TO_CONF = "config/env_conf.json"
config_json = json.load(open(os.path.join(dir_path, PATH_TO_CONF)))

GO_OBO_URL = 'http://purl.obolibrary.org/obo/go/go-basic.obo'
GO_ASSOCIATION_GENE2GEO_URL = 'https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2go.gz'
GO_FILE_NAME = 'go_bp.obo' #'go-basic.obo'
GO_ASSOCIATION_FILE_NAME = "gene2go"