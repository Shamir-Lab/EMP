import pandas as pd
import numpy as np
import os
import src.constants as constants
import time
import shutil
import sys
import json
import pandas as pd
from functools import reduce

# from df_helpers import to_full_list
# from df_helpers import to_full_np
# from utils.scripts import format_script
from src.utils.ensembl2gene_symbol import e2g_convertor
import zipfile

from src.utils.go import check_group_enrichment
import src.utils.go as go
import multiprocessing
from src.utils.daemon_multiprocessing import func_star

SH_MODULE_NAME = "module"
SH_NUM_GENES = "#_genes"
SH_ENRICHED = "enriched_groups"
SH_DETAILS = "more_details"
SH_TABLE_HEADERS = [SH_MODULE_NAME, SH_NUM_GENES, SH_ENRICHED, SH_DETAILS]

MODULE_TH = 10

def zipdir(path_to_zip, zip_file_path):
    ziph = zipfile.ZipFile(zip_file_path, 'w', zipfile.ZIP_DEFLATED)
    for root, dirs, files in os.walk(path_to_zip):
        for file in files:
            ziph.write(os.path.join(root, file))

def get_network_genes(network_file):
    df_network = pd.read_csv(network_file, sep="\t")
    src = np.array(df_network.iloc[:,0])
    dst = np.array(df_network.iloc[:,2])
    vertices = list(set(np.append(src, dst)))
    return vertices

def remove_subgraph_self_loops(nodes_to_remove, network_file_name):
    if len(nodes_to_remove) == 0:
        return network_file_name
    network_df = pd.read_csv(network_file_name, sep="\t")
    filtered_network = network_df[network_df.iloc[:,0]!=network_df.iloc[:,2]]
    new_file_name = os.path.splitext(network_file_name) + "_no_loops" +".sif"
    filtered_network.to_csv(new_file_name, sep="\t", index=False)
    return filtered_network

def remove_subgraph_by_nodes(nodes_to_remove, network_file_name, ts=str(time.time())):
    if len(nodes_to_remove) == 0:
        return network_file_name
    network_df = pd.read_csv(network_file_name, sep="\t")
    filtered_network = network_df[~(network_df.iloc[:,0].isin(nodes_to_remove) | network_df.iloc[:,2].isin(nodes_to_remove))]
    new_file_name = os.path.splitext(network_file_name)[0] + ts +".sif"
    filtered_network.to_csv(new_file_name, sep="\t", index=False)
    return new_file_name



# def summary_intergrative_reports(all_hg_reports, modules_summary, total_hg_report, algo_name, module_genes, report_file_name):
#
#     general_algo_report(algo_name, all_hg_reports, module_genes, modules_summary, report_file_name, total_hg_report)
#
#     if constants.EMB_MODE:
#         module_enrichment_report(algo_name, report_file_name, "hg_samples", total_hg_report)


def module_enrichment_report(algo_name, report_file_name, hg_sample_file_name, hg_report, output_folder):
    samples = [{go.HG_GO_ID : cur_term[go.HG_GO_ID], go.HG_GO_NAME : cur_term[go.HG_GO_NAME], go.HG_GO_ROOT : cur_term[go.HG_GO_ROOT], go.HG_VALUE : cur_term[go.HG_VALUE], go.HG_PVAL : cur_term[go.HG_PVAL], go.HG_QVAL : -1 } for cur_term in hg_report] # cur_term[go.HG_QVAL]
    df_emb = pd.DataFrame(samples)
    df_emb.to_csv(os.path.join(output_folder, "{}_{}.tsv".format(report_file_name, hg_sample_file_name)),sep="\t", index=False)
    return df_emb


def general_algo_report(algo_name, all_hg_reports, module_genes, modules_summary, report_file_name, total_hg_report, dataset_name):
    data = {}
    if len(modules_summary) > 0 :
        df_summary = pd.DataFrame(modules_summary)
        data = {"num_of_modules": df_summary.index.size,
                "module_size_avg": df_summary[SH_NUM_GENES].mean(),
                "module_size_std": df_summary[SH_NUM_GENES].std(),
                "total_num_genes": len(module_genes)
                }


    if len(all_hg_reports) > 0:
        df_all_hg = [pd.DataFrame(x) for x in all_hg_reports]
        enrichment_dist = [x.index.size for x in df_all_hg]
        pval_dist = [np.array(x[go.HG_PVAL].astype(np.float).apply(lambda x: -np.log10(x))) if x.index.size > 0 else np.array([]) for x in df_all_hg]
        modules_enrichment_data = {"module_enriched_terms_avg": np.average(enrichment_dist),
                                   "module_enriched_terms_std": np.std(enrichment_dist),
                                   "module_enriched_terms_signal_avg_avg": np.average([np.average(x) if len(x) > 1 else 0
                                                                                        for x in pval_dist]),
                                   "module_enriched_terms_signal_avg_std": np.average([np.std(x) if len(x) > 1 else 0
                                        for x in pval_dist]),
                                   "module_enriched_terms_signal_std_avg": np.std([np.average(x) if len(x) > 1 else 0
                                        for x in pval_dist]),
                                   "module_enriched_terms_signal_std_std": np.std([np.std(x) if len(x) > 1 else 0
                                        for x in pval_dist])}


        data.update(modules_enrichment_data)
        data["module_enriched_terms_signal_score"] = \
            data['module_enriched_terms_signal_avg_avg'] / ((data[
                                                                 'module_enriched_terms_signal_avg_std'] +
                                                             data[
                                                                 'module_enriched_terms_signal_std_avg'] +
                                                             data[
                                                                 'module_enriched_terms_signal_std_std']) * \
                                                            (data[
                                                                 "total_num_genes"] /
                                                             data[
                                                                 "module_size_avg"] *
                                                             data[
                                                                 "num_of_modules"]))

    if len(total_hg_report) > 0:
        df_total_hg = pd.DataFrame(total_hg_report)
        all_enrichment_data = {
            "total_enriched_terms_avg" : df_total_hg[go.HG_PVAL].astype(np.float).apply(lambda x: -np.log10(x)).mean(),
            "total_enriched_terms_std" : df_total_hg[go.HG_PVAL].astype(np.float).apply(lambda x: -np.log10(x)).std(),
            "total_num_enriched_terms": len(total_hg_report)
        }
        data.update(all_enrichment_data)

    df = pd.DataFrame()
    if len(data) >0:
        df = pd.DataFrame([data])

    df.to_csv(
        os.path.join(constants.OUTPUT_GLOBAL_DIR, dataset_name, algo_name,
                     "{}_general.tsv".format(report_file_name)), sep="\t", index=False)



def output_modules(output_file_name, modules, score_file_name, output_base_dir=""):
    output_data = create_modules_output(modules, score_file_name)
    open(output_file_name, 'w+').write(output_base_dir + "\n")
    json.dump(output_data, open(output_file_name, 'a+'))
    sys.stdout.write(output_file_name)

def reduce_to_dict(x,y):
    if y["id"] in x:
        x[y["id"]]["modules"] = x[y["id"]]["modules"] + y["modules"]
    else:
        x[y["id"]]=y
    return x

def merge_two_dicts(x, y):

    z = x.copy()
    z.update(y)
    return z

def create_modules_output(modules, score_file_name):
    scores=None
    if score_file_name is not None:
       print("score_file_name: {}".format(score_file_name))
       print(pd.read_csv(score_file_name,sep="\t").columns)
       scores = pd.read_csv(score_file_name,sep="\t").set_index("id")

       if constants.IS_PVAL_SCORES:
            scores["score"] = scores["pval"].apply(lambda x: -np.log10(x))

    zero_scores = [ {"score" : 0, "id" : gene} for module in modules for gene in module if scores is None or gene not in scores.index]
    if len(zero_scores) !=0:
        zero_scores = pd.DataFrame(zero_scores).set_index("id")
        zero_scores=zero_scores[~zero_scores.index.duplicated(keep='first')]
        scores = pd.concat([scores, zero_scores],axis=0)
    return [merge_two_dicts({"id" : k}, v) for k,v in reduce(reduce_to_dict, [{"eid": gene, "modules": [i], "id": gene, "gene_symbol": e2g_convertor([gene])[0], "score" : scores.loc[gene,"score"]} for i, module in enumerate(modules) for gene in module],\
            {}).iteritems()]

def draw_network(modules, score_file_name, network_file_name, h_src="ID_interactor_A", h_dst="ID_interactor_B"):
    active_genes = [y for x in modules for y in x]
    output = [{"data" : x, "label" : x["eid"], "selected" : True } for x in create_modules_output(modules, score_file_name)]
    active_edges = [[x[h_src], x[h_dst]] for i, x in pd.read_csv(os.path.join(constants.NETWORKS_DIR, network_file_name), sep="\t").iterrows() if x[h_src] in active_genes and x[h_dst] in active_genes]
    additional_edges = [[x[h_src], x[h_dst]] for i, x in pd.read_csv(os.path.join(constants.NETWORKS_DIR, network_file_name), sep="\t").iterrows() if not (x[h_src] in active_genes and x[h_dst] in active_genes) and (x[h_src] in active_genes or x[h_dst] in active_genes)]
    additional_nodes = [y for x in (active_edges + additional_edges) for y in x if y if y not in active_genes]
    additional_nodes = list(set(additional_nodes))

    return output + [{"data" : {"id" : x, "eid" : x, "modules" : []}, "label" : ""} for x in additional_nodes] + [{"data": {"id" : x[0]+"_"+x[1], "source":x[0], "target":x[1]}, "label" : ""} for x in additional_edges] + [{"data": {"id" : x[0]+"_"+x[1], "source":x[0], "target":x[1]}, "label" : "-"} for x in active_edges]



def build_all_reports(algo_name, modules, all_bg_genes, go_folder, output_folder):


    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    manager=multiprocessing.Manager()
    all_hg_reports = manager.list()
    modules_summary = manager.list()

    params=[]
    # p=multiprocessing.Pool(3)
    for i, module in enumerate(modules):
        # params.append([module_report, [algo_name, i, module, all_bg_genes[i], score_file_name, network_file_name, dataset_name, all_hg_reports,
        #      modules_summary]])
        module_report(algo_name, i, module, all_bg_genes[i], go_folder, output_folder, modules_summary)

    # p.map(func_star, params)

    modules_summary=list(modules_summary)
    all_hg_reports=list(all_hg_reports)

    df_summary = pd.DataFrame([], columns=['#_genes'])
    df_summary.index.name="module"
    bg_genes = []
    if len(modules) > 0:
        df_summary=pd.DataFrame(modules_summary).set_index("module")
        bg_genes = all_bg_genes[0]

    df_summary.to_csv(
        os.path.join(output_folder, "modules_summary.tsv"), sep="\t")
    # generate_algo_report(algo_name, modules, bg_genes, all_hg_reports, modules_summary, "all_modules")


    return output_folder



# def generate_algo_report(algo_name, modules, bg_genes, all_hg_reports, modules_summary, report_name):
#     hg_report = []
#     module_genes = list(set([gene for module in modules for gene in module]))
#     if (constants.HG_MODE or constants.EMB_MODE) and constants.ALGO_HG_MODE:
#         hg_report = check_group_enrichment(module_genes, bg_genes, algo_name)
#
#     summary_intergrative_reports(all_hg_reports, modules_summary, hg_report, algo_name, module_genes, report_name)


def module_report(algo_name, module_index, module, bg_genes, go_folder, output_folder, modules_summary):
    print("summarize module {} for algo {}".format(module_index, algo_name))

    open(os.path.join(output_folder, "{}_module_genes_{}.txt".format(algo_name, module_index)), "w+").write(
        "\n".join(module)) 
    try:
        open(os.path.join(output_folder, "{}_bg_genes_{}.txt".format(algo_name, module_index)), "w+").write(
        "\n".join(bg_genes))
    except Exception as e:
        # print([a for a in bg_genes if type(a) != str])
        print(e)
        raise e
         
    modules_summary_row = {SH_MODULE_NAME: module_index, SH_NUM_GENES: len(module)}

    hg_report = check_group_enrichment(list(module), list(bg_genes), go_folder)
    modules_summary_row[SH_ENRICHED] = len(hg_report)
    module_enrichment_report(algo_name, "module_" + str(module_index), "separated_modules_hg_samples", hg_report, output_folder)

    if modules_summary is not None:
        modules_summary.append(modules_summary_row)

    return modules_summary_row



