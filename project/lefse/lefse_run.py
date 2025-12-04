#!/usr/bin/env python3

import sys, argparse
from lefse import (init, load_data, get_class_means, test_kw_r, 
                   test_rep_wilcoxon_r, test_cohens_d, test_svm, save_res)


def read_params(args):
    parser = argparse.ArgumentParser(description='LEfSe 1.1.01 with Enhanced Debugging')
    parser.add_argument('input_file', metavar='INPUT_FILE', type=str, help="the input file")
    parser.add_argument('output_file', metavar='OUTPUT_FILE', type=str,
                help="the output file containing the data for the visualization module")
    parser.add_argument('-o',dest="out_text_file", metavar='str', type=str, default="",
                help="set the file for exporting the result (only concise textual form)")
    parser.add_argument('-a',dest="anova_alpha", metavar='float', type=float, default=0.05,
                help="set the alpha value for the Anova test (default 0.05)")
    parser.add_argument('-w',dest="wilcoxon_alpha", metavar='float', type=float, default=0.05,
                help="set the alpha value for the Wilcoxon test (default 0.05)")
    parser.add_argument('-l',dest="cohens_abs_th", metavar='float', type=float, default=2.0,
                help="set the threshold on the absolute value of the logarithmic COHENS score (default 2.0)")
    parser.add_argument('--nlogs',dest="nlogs", metavar='int', type=int, default=3,
        help="max log ingluence of COHENS D score")
    parser.add_argument('--verbose',dest="verbose", metavar='int', choices=[0,1], type=int, default=0,
        help="verbose execution (default 0)")
    parser.add_argument('--wilc',dest="wilc", metavar='int', choices=[0,1], type=int, default=1,
        help="wheter to perform the Wicoxon step (default 1)")
    parser.add_argument('-r',dest="rank_tec", metavar='str', choices=['cohens','svm'], type=str, default='cohens',
        help="select COHENS or SVM for effect size (default COHENS)")
    parser.add_argument('--svm_norm',dest="svm_norm", metavar='int', choices=[0,1], type=int, default=1,
        help="whether to normalize the data in [0,1] for SVM feature waiting (default 1 strongly suggested)")
    parser.add_argument('-b',dest="n_boots", metavar='int', type=int, default=30,
                help="set the number of bootstrap iteration for COHENS (default 30)")
    parser.add_argument('-e',dest="only_same_subcl", metavar='int', type=int, default=0,
                help="set whether perform the wilcoxon test only among the subclasses with the same name (default 0)")
    parser.add_argument('-c',dest="curv", metavar='int', type=int, default=0,
                help="set whether perform the wilcoxon test ing the Curtis's approach [BETA VERSION] (default 0)")
    parser.add_argument('-f',dest="f_boots", metavar='float', type=float, default=0.67,
                help="set the subsampling fraction value for each bootstrap iteration (default 0.66666)")
    parser.add_argument('-s',dest="strict", choices=[0,1,2], type=int, default=0,
                help="set the multiple testing correction options. 0 no correction (more strict, default), 1 correction for independent comparisons, 2 correction for dependent comparison")
    parser.add_argument('--min_c',dest="min_c", metavar='int', type=int, default=10,
                help="minimum number of samples per subclass for performing wilcoxon test (default 10)")
    parser.add_argument('-t',dest="title", metavar='str', type=str, default="",
                help="set the title of the analysis (default input file without extension)")
    parser.add_argument('-y',dest="multiclass_strat", choices=[0,1], type=int, default=0,
                help="(for multiclass tasks) set whether the test is performed in a one-against-one ( 1 - more strict!) or in a one-against-all setting ( 0 - less strict) (default 0)")
    args = parser.parse_args()

    params = vars(args)
    if params['title'] == "":
        params['title'] = params['input_file'].split("/")[-1].split('.')[0]

    return params


def lefse_run():
    """Main LEfSe execution with debugging and validation"""
    params = None
    logger = None

    init()
    params = read_params(sys.argv)
    
    feats, cls, class_sl, subclass_sl, class_hierarchy = load_data(params['input_file'])

    kord, cls_means = get_class_means(class_sl, feats)
    
    wilcoxon_res = {}
    kw_n_ok = 0
    nf = 0
    kw_rejected = 0
    wilc_rejected = 0
    
    total_features = len(feats)
    feat_items = list(feats.items())
    
    for feat_name, feat_values in feat_items:
        nf += 1
        
        try:
            kw_ok, pv = test_kw_r(cls, feat_values, params['anova_alpha'], 
                                    sorted(cls.keys()))
            
            if not kw_ok:
                del feats[feat_name]
                wilcoxon_res[feat_name] = "-"
                kw_rejected += 1
                continue
                
            kw_n_ok += 1
            
            # Wilcoxon test
            if not params['wilc']: 
                continue
            
            res_wilcoxon_rep = test_rep_wilcoxon_r(
                subclass_sl, class_hierarchy, feat_values,
                params['wilcoxon_alpha'], params['multiclass_strat'],
                params['strict'], feat_name, params['min_c'],
                params['only_same_subcl'], params['curv']
            )
            
            wilcoxon_res[feat_name] = str(pv) if res_wilcoxon_rep else "-"
            
            if not res_wilcoxon_rep:
                del feats[feat_name]
                wilc_rejected += 1
                
        except Exception as e:
            if feat_name in feats:
                del feats[feat_name]
            wilcoxon_res[feat_name] = "-"
            continue
    
    if len(feats) > 0:
        class_counts = {}
        for class_name, samples in cls.items():
            class_counts[class_name] = len(samples)
        
        
        if params['cohens_abs_th'] < 0.0:
            cohens_res = dict([(k, 0.0) for k, v in feats.items()])
            cohens_res_th = dict([(k, v) for k, v in feats.items()])
        else:
            if params['rank_tec'] == 'cohens':
                cohens_res, cohens_res_th = test_cohens_d(
                    cls, feats, class_sl, params['n_boots'],
                    params['f_boots'], params['cohens_abs_th'],
                    1e-7, params['nlogs']
                )
                
            elif params['rank_tec'] == 'svm':      
                cohens_res, cohens_res_th = test_svm(
                    cls, feats, class_sl, params['n_boots'],
                    params['f_boots'], params['cohens_abs_th'],
                    0.0, params['svm_norm']
                )
                
                if cohens_res is None or cohens_res_th is None:
                    raise ValueError("SVM computation failed - check data quality")
                
            else:
                cohens_res = dict([(k, 0.0) for k, v in feats.items()])
                cohens_res_th = dict([(k, v) for k, v in feats.items()])
        

    else:
        cohens_res, cohens_res_th = {}, {}

    outres = {}
    outres['cohens_res_th'] = cohens_res_th
    outres['cohens_res'] = cohens_res
    outres['cls_means'] = cls_means
    outres['cls_means_kord'] = kord
    outres['wilcox_res'] = wilcoxon_res
    
    save_res(outres, params["output_file"])


if __name__ == '__main__':
    lefse_run()