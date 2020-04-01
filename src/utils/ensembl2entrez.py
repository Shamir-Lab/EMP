import os
from src import constants


g2e_dict = None
e2g_dict = None
ensembl2entrez_dict = None
entrez2ensembl_dict = None

def load_gene_dictionary(gene_list_file_name, gene_list_path=None):
    if gene_list_path == None:
        gene_list_path = os.path.join(constants.dir_path,"data",gene_list_file_name)
    f = open(gene_list_path,'r')
    lines = [l.strip() for l in f]
    f.close()
    return lines


def get_entrez2ensembl_dictionary():
    lines_dict = load_gene_dictionary(constants.ENSEMBL_TO_ENTREZ)

    entrez2ensembl = {}
    for cur in lines_dict:
        splited_line = cur.split()
        if len(splited_line) != 2: continue
        if splited_line[0].find('.') > 0:
            limit = splited_line[0].find('.')
        else:
            limit = len(splited_line[0])
            entrez2ensembl[splited_line[1]] = splited_line[0][:limit]
    return entrez2ensembl

def get_ensembl2entrez_dictionary():

    lines_dict = load_gene_dictionary(constants.ENSEMBL_TO_ENTREZ)
    ensembl2gene_symbols = {}
    for cur in lines_dict:
        splited_line = cur.split()
        if len(splited_line) !=2: continue
        if splited_line[0].find('.') > 0:
            limit = splited_line[0].find('.')
        else:
            limit = len(splited_line[0])
        ensembl2gene_symbols[splited_line[0][:limit]] = splited_line[1]
    return ensembl2gene_symbols


def ensembl2entrez_convertor(e_ids):
    global ensembl2entrez_dict
    if ensembl2entrez_dict is None:
        ensembl2entrez_dict = get_ensembl2entrez_dictionary()
    results = []
    for cur in e_ids:
        if cur.split(".")[0] in ensembl2entrez_dict:
            results.append(ensembl2entrez_dict[cur.split(".")[0]])
    return results


def entrez2ensembl_convertor(entrez_ids):
    global entrez2ensembl_dict
    if entrez2ensembl_dict is None:
        entrez2ensembl_dict = get_entrez2ensembl_dictionary()
    results = []
    for cur in entrez_ids:
        cur=str(cur)
        if entrez2ensembl_dict.has_key(cur):
            results.append(entrez2ensembl_dict[cur])
    return results