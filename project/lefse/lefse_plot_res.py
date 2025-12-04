#!/usr/bin/env python3

import os, sys
import matplotlib
matplotlib.use('Agg')
from pylab import *
from collections import defaultdict

from lefse import *
import argparse

colors = ['r','g','b','m','c','y','k','w']

def read_params(args):
    parser = argparse.ArgumentParser(description='Plot Cohen\'s D results')
    parser.add_argument('input_file', metavar='INPUT_FILE', type=str, help="tab delimited input file")
    parser.add_argument('output_file', metavar='OUTPUT_FILE', type=str, help="the file for the output image")
    parser.add_argument('--feature_font_size', dest="feature_font_size", type=int, default=7, help="font size for feature labels")
    parser.add_argument('--format', dest="format", choices=["png","svg","pdf"], default='png', type=str, help="output file format")
    parser.add_argument('--dpi', dest="dpi", type=int, default=150, help="resolution in dots per inch")
    parser.add_argument('--title', dest="title", type=str, default="Cohen's D Effect Size")
    parser.add_argument('--title_font_size', dest="title_font_size", type=int, default=14)
    parser.add_argument('--class_legend_font_size', dest="class_legend_font_size", type=int, default=10)
    parser.add_argument('--width', dest="width", type=float, default=8.0)
    parser.add_argument('--height', dest="height", type=float, default=6.0, help="height for vertical histograms")
    parser.add_argument('--left_space', dest="ls", type=float, default=0.25)
    parser.add_argument('--right_space', dest="rs", type=float, default=0.1)
    parser.add_argument('--orientation', dest="orientation", type=str, choices=["h","v"], default="h")
    parser.add_argument('--autoscale', dest="autoscale", type=int, choices=[0,1], default=1)
    parser.add_argument('--background_color', dest="back_color", type=str, choices=["k","w"], default="w", help="background color")
    parser.add_argument('--subclades', dest="n_scl", type=int, default=1, help="number of label levels to display (-1 for all)")
    parser.add_argument('--max_feature_len', dest="max_feature_len", type=int, default=60, help="maximum length of feature strings")
    parser.add_argument('--all_feats', dest="all_feats", type=str, default="")
    parser.add_argument('--otu_only', dest="otu_only", default=False, action='store_true', help="plot only species-resolved OTUs")
    parser.add_argument('--report_features', dest="report_features", default=False, action='store_true', help="report features to STDOUT")
    parser.add_argument('--effect_size_threshold', dest="threshold", type=float, default=None, help="Cohen's d threshold line (e.g., 0.5 for medium effect)")
    parser.add_argument('--show_effect_labels', dest="show_effect_labels", default=True, action='store_true', help="show effect size interpretation labels")
    args = parser.parse_args()
    return vars(args)

def get_effect_size_label(d_value):
    """Return Cohen's d effect size interpretation"""
    d = abs(d_value)
    if d < 0.2:
        return "negligible"
    elif d < 0.5:
        return "small"
    elif d < 0.8:
        return "medium"
    else:
        return "large"

def read_data(input_file, output_file, otu_only):
    with open(input_file, 'r') as inp:
        if not otu_only:
            rows = [line.strip().split()[:-1] for line in inp.readlines() if len(line.strip().split()) > 3]
        else:
            rows = [line.strip().split()[:-1] for line in inp.readlines() 
                   if len(line.strip().split()) > 3 and len(line.strip().split()[0].split('.')) == 8]
    
    classes = list(set([v[2] for v in rows if len(v) > 2]))
    if len(classes) < 1: 
        print("No differentially abundant features found in " + input_file)
        os.system("touch " + output_file)
        sys.exit()
    
    data = {}
    data['rows'] = rows
    data['cls'] = classes
    return data

def plot_histo_hor(path, params, data, bcl, report_features):
    cls2 = []
    if params['all_feats'] != "":
        cls2 = sorted(params['all_feats'].split(":"))
    cls = sorted(data['cls'])
    
    if bcl: 
        data['rows'].sort(key=lambda ab: fabs(float(ab[3])) * (cls.index(ab[2]) * 2 - 1))
    else: 
        mmax = max([fabs(float(a)) for a in list(zip(*list(data['rows'])))[3]])
        data['rows'].sort(key=lambda ab: fabs(float(ab[3])) / mmax + (cls.index(ab[2]) + 1))
    
    pos = arange(len(data['rows']))
    head = 0.75
    tail = 0.5
    ht = head + tail
    ints = max(len(pos) * 0.2, 1.5)
    fig = plt.figure(figsize=(params['width'], ints + ht), 
                     edgecolor=params['back_color'], facecolor=params['back_color'])
    ax = fig.add_subplot(111, frame_on=False, facecolor=params['back_color'])
    ls, rs = params['ls'], 1.0 - params['rs']
    plt.subplots_adjust(left=ls, right=rs, 
                       top=1 - head * (1.0 - ints / (ints + ht)), 
                       bottom=tail * (1.0 - ints / (ints + ht)))

    fig.canvas.manager.set_window_title('Cohen\'s D Effect Size Results')

    l_align = {'horizontalalignment': 'left', 'verticalalignment': 'baseline'}
    r_align = {'horizontalalignment': 'right', 'verticalalignment': 'baseline'}
    added = []
    m = 1 if data['rows'][0][2] == cls[0] else -1
    out_data = defaultdict(list)
    
    for i, v in enumerate(data['rows']):
        if report_features:
            otu = v[0].split('.')[7].replace('_', '.') if len(v[0].split('.')) > 7 else v[0]
            score = v[3]
            otu_class = v[2]
            effect_label = get_effect_size_label(float(score))
            out_data[otu] = [score, otu_class, effect_label]
        
        indcl = cls.index(v[2])
        lab = str(v[2]) if str(v[2]) not in added else ""
        added.append(str(v[2])) 
        col = colors[indcl % len(colors)] 
        if len(cls2) > 0: 
            col = colors[cls2.index(v[2]) % len(colors)]
        vv = fabs(float(v[3])) * (m * (indcl * 2 - 1)) if bcl else fabs(float(v[3]))
        ax.barh(pos[i], vv, align='center', color=col, label=lab, height=0.8, 
               edgecolor=params['fore_color'], alpha=0.8)
    
    mv = max([abs(float(v[3])) for v in data['rows']])
    
    if report_features:
        print('Feature\tCohen\'s_D\tClass\tEffect_Size')
        for otu in out_data:
            print('%s\t%s\t%s\t%s' % (otu, out_data[otu][0], out_data[otu][1], out_data[otu][2]))
    
    if params['threshold'] is not None:
        ax.axvline(x=params['threshold'], color='gray', linestyle='--', linewidth=1, alpha=0.7, label=f'd={params["threshold"]}')
        if bcl:
            ax.axvline(x=-params['threshold'], color='gray', linestyle='--', linewidth=1, alpha=0.7)
    
    if params['show_effect_labels'] and bcl:
        xlim = ax.get_xlim()
        ax.axvspan(0, 0.2, alpha=0.05, color='gray', label='negligible')
        ax.axvspan(0.2, 0.5, alpha=0.1, color='yellow', label='small')
        ax.axvspan(0.5, 0.8, alpha=0.1, color='orange', label='medium')
        ax.axvspan(0.8, xlim[1], alpha=0.1, color='red', label='large')
        ax.axvspan(0, -0.2, alpha=0.05, color='gray')
        ax.axvspan(-0.2, -0.5, alpha=0.1, color='yellow')
        ax.axvspan(-0.5, -0.8, alpha=0.1, color='orange')
        ax.axvspan(-0.8, xlim[0], alpha=0.1, color='red')
    
    for i, r in enumerate(data['rows']):
        indcl = cls.index(data['rows'][i][2])
        if params['n_scl'] < 0: 
            rr = r[0]
        else: 
            rr = ".".join(r[0].split(".")[-params['n_scl']:])
        if len(rr) > params['max_feature_len']: 
            rr = rr[:params['max_feature_len'] // 2 - 2] + " [..]" + rr[-params['max_feature_len'] // 2 + 2:]
        if m * (indcl * 2 - 1) < 0 and bcl: 
            ax.text(mv / 40.0, float(i) - 0.3, rr, l_align, 
                   size=params['feature_font_size'], color=params['fore_color'])
        else: 
            ax.text(-mv / 40.0, float(i) - 0.3, rr, r_align, 
                   size=params['feature_font_size'], color=params['fore_color'])
    
    ax.set_title(params['title'], size=params['title_font_size'], 
                y=1.0 + head * (1.0 - ints / (ints + ht)) * 0.8, 
                color=params['fore_color'], weight='bold')
    
    ax.set_yticks([])
    ax.set_xlabel("Cohen's D Effect Size", fontsize=12, weight='bold')
    ax.xaxis.grid(True, alpha=0.3)
    
    xlim = ax.get_xlim()
    if params['autoscale']: 
        ran = arange(0.0001, round(round((abs(xlim[0]) + abs(xlim[1])) / 10, 4) * 100, 0) / 100)
        if len(ran) > 1 and len(ran) < 100:
            ax.set_xticks(arange(xlim[0], xlim[1] + 0.0001, 
                                min(xlim[1] + 0.0001, round(round((abs(xlim[0]) + abs(xlim[1])) / 10, 4) * 100, 0) / 100)))
    
    ax.set_ylim((pos[0] - 1, pos[-1] + 1))
    leg = ax.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3, ncol=5, 
                   borderaxespad=0., frameon=False, 
                   prop={'size': params['class_legend_font_size']})

    def get_col_attr(x):
        return hasattr(x, 'set_color') and not hasattr(x, 'set_facecolor')
    
    for o in leg.findobj(get_col_attr):
        o.set_color(params['fore_color'])
    for o in ax.findobj(get_col_attr):
        o.set_color(params['fore_color'])
    
    plt.savefig(path, format=params['format'], facecolor=params['back_color'], 
               edgecolor=params['fore_color'], dpi=params['dpi'], bbox_inches='tight')
    plt.close()

def plot_histo_ver(path, params, data, report_features):
    cls = data['cls']
    mmax = max([fabs(float(a)) for a in zip(*data['rows'])[3]])
    data['rows'].sort(key=lambda ab: fabs(float(ab[3])) / mmax + (cls.index(ab[2]) + 1))
    pos = arange(len(data['rows'])) 
    
    if params['n_scl'] < 0: 
        nam = [d[0] for d in data['rows']]
    else: 
        nam = [d[0].split(".")[-min(d[0].count("."), params['n_scl'])] for d in data['rows']]
    
    fig = plt.figure(edgecolor=params['back_color'], facecolor=params['back_color'], 
                     figsize=(params['width'], params['height'])) 
    ax = fig.add_subplot(111, facecolor=params['back_color'])
    plt.subplots_adjust(top=0.9, left=params['ls'], right=params['rs'], bottom=0.3) 
    fig.canvas.manager.set_window_title('Cohen\'s D Effect Size Results')
    
    added = []
    out_data = defaultdict(list)
    
    for i, v in enumerate(data['rows']):
        if report_features:
            otu = v[0].split('.')[7].replace('_', '.') if len(v[0].split('.')) > 7 else v[0]
            score = v[3]
            otu_class = v[2]
            effect_label = get_effect_size_label(float(score))
            out_data[otu] = [score, otu_class, effect_label]
        
        indcl = data['cls'].index(v[2])
        lab = str(v[2]) if str(v[2]) not in added else ""
        added.append(str(v[2])) 
        col = colors[indcl % len(colors)]
        vv = fabs(float(v[3])) 
        ax.bar(pos[i], vv, align='center', color=col, label=lab, alpha=0.8)
    
    if params['threshold'] is not None:
        ax.axhline(y=params['threshold'], color='gray', linestyle='--', 
                  linewidth=1, alpha=0.7, label=f'd={params["threshold"]}')
    
    if report_features:
        print('Feature\tCohen\'s_D\tClass\tEffect_Size')
        for otu in out_data:
            print('%s\t%s\t%s\t%s' % (otu, out_data[otu][0], out_data[otu][1], out_data[otu][2]))
    
    xticks(pos, nam, rotation=-20, ha='left', size=params['feature_font_size'])  
    ax.set_title(params['title'], size=params['title_font_size'], weight='bold')
    ax.set_ylabel("Cohen's D Effect Size", fontsize=12, weight='bold')
    ax.yaxis.grid(True, alpha=0.3) 
    
    a, b = ax.get_xlim()
    dx = float(len(pos)) / float((b - a))
    ax.set_xlim((0 - dx, max(pos) + dx)) 
    
    plt.savefig(path, format=params['format'], facecolor=params['back_color'], 
               edgecolor=params['fore_color'], dpi=params['dpi'], bbox_inches='tight')
    plt.close() 

def plot_res():
    params = read_params(sys.argv)
    params['fore_color'] = 'w' if params['back_color'] == 'k' else 'k'
    data = read_data(params['input_file'], params['output_file'], params['otu_only'])
    
    if params['orientation'] == 'v': 
        plot_histo_ver(params['output_file'], params, data, params['report_features'])
    else: 
        plot_histo_hor(params['output_file'], params, data, len(data['cls']) == 2, 
                      params['report_features'])

if __name__ == '__main__':
    plot_res()