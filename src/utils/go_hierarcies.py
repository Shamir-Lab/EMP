"""Print a GO term's lower-level hierarchy."""

__copyright__ = "Copyright (C) 2016-2018, DV Klopfenstein, H Tang. All rights reserved."
__author__ = "DV Klopfenstein"

import sys
sys.path.insert(0, '../')
import collections as cx
from goatools.gosubdag.go_paths import GoPaths
# from goatools.godag.consts import Consts
"""driver imports"""
import timeit
import datetime
from goatools.obo_parser import GODag
from goatools.gosubdag.gosubdag import GoSubDag
from goatools.gosubdag.rpt.write_hierarchy import WrHierGO
""""gene to go mapper imports"""
from  goatools.associations import read_ncbi_gene2go
from src import constants
from src.utils.ensembl2entrez import get_entrez2ensembl_dictionary
import wget, os
import shutil
import gzip

class WrHierGO(object):
    """Write hierarchy object."""

    kws_dct = set(['max_indent'])
    kws_set = set(['no_indent', 'concise'])

    def __init__(self, gosubdag, **kws):
        self.gosubdag = gosubdag  # GoSubDag arg, children=True, must be used
        self.usrdct = {k:v for k, v in kws.items() if k in kws}
        self.usrset = set([k for k, v in kws.items() if k in kws and v])
        # ' {NS} {dcnt:6,} L{level:02} D{depth:02} {D1:5} {GO_name}'

    def prt_hier_all(self, prt=sys.stdout):
        """Write hierarchy for all GO Terms in obo file."""
        # Print: [biological_process, molecular_function, and cellular_component]
        gos_printed = set()
        for goid in ['GO:0008150', 'GO:0003674', 'GO:0005575']:
            gos_printed.update(self.prt_hier_down(goid, prt))
        return gos_printed

    def prt_hier_down(self, goid, prt=sys.stdout):
        """Write hierarchy for all GO IDs below GO ID in arg, goid."""
        obj = _WrHierPrt(self, prt)
        obj.prt_hier_rec(goid)
        return obj.gos_printed

    def ext_hier_down(self, goid, prt=sys.stdout):
        """Write hierarchy for all GO IDs below GO ID in arg, goid."""

        obj = _WrHierPrt(self, prt)
        obj.ext_hier_rec(goid)

        return obj

    def prt_hier_up(self, goids, prt=sys.stdout):
        """Write hierarchy for all GO IDs below GO ID in arg, goid."""
        go2goterm_all = {go:self.gosubdag.go2obj[go] for go in goids}
        objp = GoPaths()
        gos_printed = set()
        for namespace, go2term_ns in self._get_namespace2go2term(go2goterm_all).items():
            go_root = self.consts.NAMESPACE2GO[namespace]
            goids_all = set()  # GO IDs from user-specfied GO to root
            for goid, goterm in go2term_ns.items():
                goids_all.add(goid)
                paths = objp.get_paths_from_to(goterm, goid_end=None, dn0_up1=True)
                goids_all.update(set(o.id for p in paths for o in p))
            # Only include GO IDs from user-specified GO to the root
            if 'include_only' not in self.usrdct:
                self.usrdct['include_only'] = set()
            self.usrdct['include_only'].update(goids_all)
            # Mark the user-specfied GO term
            if 'go_marks' not in self.usrdct:
                self.usrdct['go_marks'] = set()
            self.usrdct['go_marks'].update(go2term_ns.keys())
            obj = _WrHierPrt(self, prt)  # , goids_all, set(go2term_ns.keys()))
            gos_printed.update(obj.gos_printed)
            obj.prt_hier_rec(go_root)
        return gos_printed

    @staticmethod
    def _get_namespace2go2term(go2terms):
        """Group GO IDs by namespace."""
        namespace2go2term = cx.defaultdict(dict)
        for goid, goterm in go2terms.items():
            namespace2go2term[goterm.namespace][goid] = goterm
        return namespace2go2term


class _WrHierPrt(object):
    """Print GO hierarchy."""

    def __init__(self, obj, prt=sys.stdout):
        self.gosubdag = obj.gosubdag
        self.max_indent = obj.usrdct.get('max_indent')
        self.include_only = obj.usrdct['include_only'] if 'include_only' in obj.usrdct else None
        self.go_marks = obj.usrdct['go_marks'] if 'go_marks' in obj.usrdct else set()
        self.concise_prt = 'concise' in obj.usrset
        self.indent = 'no_indent' not in obj.usrset
        self.go2geneids = obj.usrdct.get('go2geneids')
        # vars
        self.prt = prt
        self.edges = {}
        self.vertices = {}
        self.gos_printed = set()
        self.prtfmt = self._init_prtfmt()
        self.dash_len = obj.usrdct.get('dash_len', 6) + 12

    def ext_hier_rec(self, goid, depth=1):
        """Write hierarchy for a GO Term record and all GO IDs down to the leaf level."""
        ntgo = self.gosubdag.go2nt[goid]
        ntobj = self.gosubdag.go2obj[goid]
        # Shortens hierarchy report by only printing the hierarchy
        # for the sub-set of user-specified GO terms which are connected.
        if self.include_only and goid not in self.include_only:
            return
        nrp = self.concise_prt and goid in self.gos_printed
        if goid in self.vertices:
            self.vertices[goid]["weight"] += 1
            self.vertices[goid]["depth"].append(depth)
        else:
            self.vertices[goid] = {"name": ntgo.GO_name,
                                   "weight": 0,
                                   "NS": ntgo.NS,
                                   "depth": [depth],
                                   "L" : ntgo.level,
                                   "D" : ntgo.depth,
				   "obj" : ntobj,
				   "n_children" : len(ntobj.children)}

        self.gos_printed.add(goid)

        depth += 1
        if self.max_indent is not None and depth > self.max_indent:
            return
        for child in ntobj.children:
            if child.id in self.go2geneids or True:
                if "{}={}".format(goid, child.id) in self.edges:
                    self.edges["{}={}".format(goid, child.id)]["weight"] += 1
                else:
                    self.edges["{}={}".format(goid, child.id)] = {"weight" : 0}
                self.ext_hier_rec(child.id, depth)

    def prt_hier_rec(self, goid, depth=1):
        """Write hierarchy for a GO Term record and all GO IDs down to the leaf level."""
        ntgo = self.gosubdag.go2nt[goid]
        ntobj = self.gosubdag.go2obj[goid]
        # Shortens hierarchy report by only printing the hierarchy
        # for the sub-set of user-specified GO terms which are connected.
        if self.include_only and goid not in self.include_only:
            return
        nrp = self.concise_prt and goid in self.gos_printed
        if self.go_marks:
            self.prt.write('{} '.format('>' if goid in self.go_marks else ' '))

        # '-' is default character indicating hierarchy level
        # '=' is used to indicate a hierarchical path printed in detail previously.
        dashgo = self._str_dashgoid(ntgo, depth, not nrp or not ntobj.children)
        self.prt.write('{DASHGO:{N}}'.format(DASHGO=dashgo, N=self.dash_len))

        self.prt.write("{GO_INFO}\n".format(GO_INFO=self.prtfmt.format(**ntgo._asdict())))
        self.gos_printed.add(goid)
        # Do not print hierarchy below this turn if it has already been printed
        if nrp:
            return
        depth += 1
        if self.max_indent is not None and depth > self.max_indent:
            return
        for child in ntobj.children:
            self.prt_hier_rec(child.id, depth)

    @staticmethod
    def _str_dash(depth, single_or_double):
        """Return a string containing dashes (optional) and GO ID."""
        # '-' is default character indicating hierarchy level
        # '=' is used to indicate a hierarchical path printed in detail previously.
        letter = '-' if single_or_double else '='
        return ''.join([letter]*depth)

    def _str_dashgoid(self, ntgo, depth, single_or_double):
        """Return a string containing dashes (optional) and GO ID."""
        dashes = self._str_dash(depth, single_or_double) if self.indent else ""
        return "{DASHES} {GO}{alt:1}".format(DASHES=dashes, GO=ntgo.GO, alt=ntgo.alt)

    def _init_prtfmt(self):
        """Initialize print format."""
        prtfmt = self.gosubdag.prt_attr['fmt']
        prtfmt = prtfmt.replace('{GO} # ', '')
        prtfmt = prtfmt.replace('{D1:5} ', '')
        return prtfmt


def write_hier_all(gosubdag, out, root_term):
    """write_hier.py: Prints the entire mini GO hierarchy, with counts of children."""
    out.write('\nTEST ALL: Print all hierarchies:\n')
    objwr = WrHierGO(gosubdag)
    gos_printed = objwr.prt_hier_down(root_term, out)
    print(len(gos_printed))
    # assert gos_printed == set(objwr.gosubdag.go2nt)


#################################################################
# Sub-routines to tests
#################################################################
def extract_hier_all(gosubdag, root_term, go2geneids):
    """write_hier.py: Prints the entire mini GO hierarchy, with counts of children."""
    # out.write('\nTEST EXTRACTION: Print all hierarchies:\n')
    objwr = WrHierGO(gosubdag,  go2geneids=go2geneids)
    obj = objwr.ext_hier_down(root_term)
    return (obj.vertices, obj.edges)



def write_hier_norep(gosubdag, out):
    """Shortens hierarchy report by only printing branches once.
         Prints the 'entire hierarchy' of GO:0000005 the 1st time seen:
           ---     1 GO:0000005    L-02    D-02
           ----     0 GO:0000010   L-03    D-04
         Prints just GO:0000005 (ommit child GO:10) the 2nd time seen:
           ===     1 GO:0000005    L-02    D-02
         '=' is used in hierarchy mark to indicate that the pathes
             below the marked term have already been printed.
    """
    # out.write('\nTEST ALL: Print branches just once:\n')
    objwr = WrHierGO(gosubdag, concise=True)
    gos_printed = objwr.prt_hier_down("GO:0000001", out)
    assert gos_printed == set(objwr.gosubdag.go2nt)


def write_hier_lim(gosubdag, out):
    """Limits hierarchy list to GO Terms specified by user."""
    go_omit = ['GO:0000005', 'GO:0000010']
    go_ids = [go_id for go_id in gosubdag.go2obj if go_id not in go_omit]
    # out.write('\nTEST OMIT: 05 and 10:\n')
    objwr = WrHierGO(gosubdag, include_only=go_ids)
    gos_printed = objwr.prt_hier_down("GO:0000001", out)
    assert not gos_printed.intersection(go_omit), "SHOULD NOT PRINT {GOs}".format(GOs=go_omit)


def write_hier_mrk(gosubdag, out):
    """Print all paths, but mark GO Terms of interest. """
    mark_lst = ['GO:0000001', 'GO:0000003', 'GO:0000006', 'GO:0000008', 'GO:0000009']
    objwr = WrHierGO(gosubdag, go_marks=mark_lst)
    objwr.prt_hier_down("GO:0000001", out)


def fetch_go_hierarcy(go_folder):
    obo_file_location = os.path.join(constants.GO_DIR, constants.GO_FILE_NAME)
    if not os.path.exists(os.path.join(constants.GO_DIR, constants.GO_FILE_NAME)):
        wget.download(constants.GO_OBO_URL, os.path.join(constants.GO_DIR, constants.GO_FILE_NAME))

    print("Downloading gene-GO associations")
    association_file_location = os.path.join(constants.GO_DIR, constants.GO_ASSOCIATION_FILE_NAME)
    # if not os.path.exists(association_file_location):
    #     wget.download(constants.GO_ASSOCIATION_GENE2GEO_URL,
    #                   os.path.join(constants.GO_DIR, constants.GO_ASSOCIATION_FILE_NAME))

    if not os.path.exists(os.path.join(go_folder, constants.GO_ASSOCIATION_FILE_NAME)):
        if not os.path.exists(os.path.join(go_folder, constants.GO_ASSOCIATION_FILE_NAME+".gz")):
            wget.download(constants.GO_ASSOCIATION_GENE2GEO_URL, os.path.join(constants.GO_DIR, os.path.split(constants.GO_ASSOCIATION_GENE2GEO_URL)[1]))
        with gzip.open(os.path.join(go_folder, os.path.basename(constants.GO_ASSOCIATION_GENE2GEO_URL)), 'rb') as f_in:
            with open(os.path.join(go_folder, constants.GO_ASSOCIATION_FILE_NAME),'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)

    print("Loading gene-GO associations")

    go2geneids = read_ncbi_gene2go(association_file_location, taxids=[9606], go2geneids=True)
    geneids2go = read_ncbi_gene2go(association_file_location, taxids=[9606])

    ## backward compatibility to goatools python 2.7##
    # all_go_ids=set().union(*list(geneids2go.values()))
    # for cur_id in all_go_ids:
    #     go2geneids[cur_id]=set()
    ############################


    return (go2geneids, geneids2go)

#################################################################
# Driver
#################################################################
def build_hierarcy(go_folder, roots=['GO:0008150']): #  0008150 0005575 0003674

    go2geneids, geneids2go = fetch_go_hierarcy(go_folder)

    """Run numerous tests for various reports."""
    dag_fin = os.path.join(constants.GO_DIR, constants.GO_FILE_NAME)
    tic = timeit.default_timer()
    godag = GODag(dag_fin, optional_attrs=['relationship'])
    gosubdag = GoSubDag(godag.keys(), godag)
    toc = timeit.default_timer()
    dict_result = {}
    for cur_term in roots:
        vertices, edges = extract_hier_all(gosubdag, cur_term ,go2geneids)

        # all_go_ids=set(vertices.keys())
        # for cur_id in all_go_ids:
        #     if not cur_id in go2geneids:
        #         go2geneids[cur_id]=set()

        msg = "Elapsed HMS: {}\n\n".format(str(datetime.timedelta(seconds=(toc-tic))))
        sys.stdout.write(msg)
        dict_result[cur_term] = {"vertices" : vertices, "edges": edges}
    return dict_result, go2geneids, geneids2go, get_entrez2ensembl_dictionary()



########################################### ######################
# main
#################################################################
if __name__ == '__main__':

    print(build_hierarcy())
