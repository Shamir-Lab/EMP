import json
import os


SCORES_FILE_NANE="scores.tsv"
ENSEMBL_TO_GENE_SYMBOLS = "ensembl2gene_symbol.txt"
ENSEMBL_TO_ENTREZ = "ensembl2entrez.txt"

dir_path = os.path.join(os.path.dirname(os.path.realpath(__file__)),"../")
PATH_TO_CONF = "config/conf.json"
config_json = json.load(open(os.path.join(dir_path, PATH_TO_CONF)))

GO_DIR=config_json["go_folder"]
GO_OBO_URL = config_json["go_obo_url"]
GO_ASSOCIATION_GENE2GEO_URL = config_json["go_association_gene2geo_url"]
GO_FILE_NAME = config_json["go_file_name"]
GO_ASSOCIATION_FILE_NAME = config_json["go_association_file_name"]
ROOT_GO_ID=config_json["root_go_id"]
EV_EXCLUDE=set(config_json["ev_exclude"].split(','))

ALGO_DIR=config_json["algo_dir"]

