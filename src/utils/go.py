import re
import gzip
import shutil


import os
import src.constants as constants

from goatools.obo_parser import GODag
from goatools.go_enrichment import GOEnrichmentStudy
from goatools.associations import read_ncbi_gene2go
from src.utils.ensembl2entrez import ensembl2entrez_convertor
from src.utils.ensembl_converter import load_gene_list
from src.utils.download_resources import download
from src.utils import go_hierarcies

HG_GO_ROOT = "GO root"
HG_GO_ID = "GO id"
HG_GO_NAME = "GO name"
HG_PVAL = "pval"
HG_QVAL = "qval"
HG_VALUE = "value"
HG_TABLE_HEADERS = [HG_GO_ROOT, HG_GO_ID, HG_GO_NAME, HG_VALUE, HG_PVAL, HG_QVAL]

assoc=None

dict_result, go2geneids, geneids2go, entrez2ensembl = go_hierarcies.build_hierarcy(
    roots=['GO:0008150'])
vertices = list(dict_result.values())[0]['vertices']

terms_to_genes = {}
ids_to_names = None


def check_enrichment(gene_list):
    ensembl_for_url = re.sub("\.\d{1,2},", ",", gene_list)
    url = "http://david.abcc.ncifcrf.gov/api.jsp?type=ENSEMBL_GENE_ID&ids={}&tool=chartReport&annot=GOTERM_BP_DIRECT,GOTERM_CC_DIRECT,GOTERM_MF_DIRECT,KEGG_PATHWAY".format(ensembl_for_url)
    return url



def check_group_enrichment(tested_gene_file_name, total_gene_file_name, go_folder, th=1):
    if len(tested_gene_file_name) == 0 or len(total_gene_file_name) == 0: return []

    if type(total_gene_file_name) == str:
        total_gene_list = load_gene_list(total_gene_file_name)
    else:
        total_gene_list = total_gene_file_name

    if type(tested_gene_file_name) == str:
        tested_gene_list = load_gene_list(tested_gene_file_name)
    else:
        tested_gene_list = tested_gene_file_name

    if not os.path.exists(os.path.join(go_folder, constants.GO_FILE_NAME)):
        download(constants.GO_OBO_URL, constants.GO_DIR)

    obo_dag = GODag(os.path.join(go_folder, constants.GO_FILE_NAME))

    if not os.path.exists(os.path.join(go_folder, constants.GO_ASSOCIATION_FILE_NAME)):
        download(constants.GO_ASSOCIATION_GENE2GEO_URL, constants.GO_DIR)
        with gzip.open(os.path.join(go_folder, os.path.basename(constants.GO_ASSOCIATION_GENE2GEO_URL)), 'rb') as f_in:
            with open(os.path.join(go_folder, constants.GO_ASSOCIATION_FILE_NAME),'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)

    global assoc
    if assoc == None:
        assoc = read_ncbi_gene2go(os.path.join(go_folder, constants.GO_ASSOCIATION_FILE_NAME), no_top=True)

    g = GOEnrichmentStudy([int(cur) for cur in ensembl2entrez_convertor(total_gene_list)],
                          assoc, obo_dag, log=None) # "bonferroni", "fdr_bh"
    g_res = g.run_study([int(cur) for cur in ensembl2entrez_convertor(tested_gene_list)])

    GO_results = [(cur.NS, cur.GO, cur.goterm.name, cur.pop_count, cur.p_uncorrected) for cur in g_res ] # , cur.p_fdr_bh    if cur.p_fdr_bh <= th


    hg_report = [{HG_GO_ROOT : cur[0], HG_GO_ID : cur[1], HG_GO_NAME : cur[2], HG_VALUE : cur[3], HG_PVAL : cur[4]} for cur in GO_results] # , HG_QVAL : cur[5]
    hg_report.sort(key=lambda x: x[HG_PVAL]) # HG_QVAL

    return hg_report


def get_all_genes_for_term(vertices, cur_root, term, in_subtree):
    if term in terms_to_genes:
        return terms_to_genes[cur_root]

    all_genes = set()
    if in_subtree:
        try:
            all_genes.update(go2geneids[cur_root])
        except Exception:
            # print(e)
            pass
            # print("cur_root : {} was not found".format(cur_root))

    for cur_child in vertices[cur_root]["obj"].children:
        all_genes.update(get_all_genes_for_term(vertices, cur_child.id, term, in_subtree))

    terms_to_genes[cur_root] = all_genes
    return all_genes


def get_go_names(GO_ids, go_folder):
    GO_names = []
    global ids_to_names
    if ids_to_names is None:
        ids_to_names = {}
        f = open(os.path.join(go_folder, 'go-basic.obo'))
        parsed = f.read().split("\n\n")
        for cur_obo in parsed[1:]:
            if cur_obo.split("\n")[0] != "[Term]": continue
            ids_to_names[cur_obo.split("\n")[1][4:]] = cur_obo.split("\n")[2][6:]

    for cur_id in GO_ids:
        if cur_id in ids_to_names:
            GO_names.append(ids_to_names[cur_id])
        else:
            GO_names.append(cur_id)

    # print "\n".join(GO_names)
    return GO_names



