####################
# ZFC
# Author: Wolfson Liu, Wei Tang
# Email: wolfsonliu@live.com, tangwei@stu.pku.edu.cn
####################

import pandas as pd
import numpy as np
import scipy.stats as stats
from sklearn import linear_model
from .statsfunc import df_normalization
from .statsfunc import df_smallcount
from .statsfunc import ecdf
from .statsfunc import p_adjust
from .statsfunc import df_robust_rank_aggregation
from .statsfunc import df_mean_rank_aggregation
from .fit_by_lowess import lowess_regression,fit_qc_plot,cal_zfc


def zfoldchange(inputdata,
                outprefix,
                top_n_sgrna=None,
                top_n_gene=None,
                min_ctrl_counts=1,
                min_ctrl_and_exp_counts=1,
                iteration=100,
                use_small_count=False,
                normalization='median'):
    # The data DF should contain: [gene, guide, barcode, ctrl, exp]

    for x in ['gene', 'guide', 'barcode', 'ctrl', 'exp']:
        assert x in inputdata.columns, 'data should have column: {}'.format(x)
    inputdata.replace(np.nan, 0, inplace=True)
    # ------------------
    # Step 1: Normalization of raw counts

    # remove ctrl < min_ctrl_counts, by group
    guide_to_drop = inputdata.loc[inputdata['ctrl'] < min_ctrl_counts, 'guide']
    inputdata = inputdata.loc[~inputdata['guide'].isin(guide_to_drop), :]

    # remove ctrl+exp < min_ctrl_and_exp_counts, by group
    grouped = inputdata.groupby("guide").agg({'ctrl':'sum','exp':'sum'})
    grouped = grouped.reset_index()
    result = grouped.loc[grouped['ctrl'] + grouped['exp'] < min_ctrl_and_exp_counts, 'guide']

    inputdata = inputdata.loc[~inputdata['guide'].isin(result), :]
    inputdata['gene_guide'] = inputdata[['gene', 'guide']].apply('_'.join, axis=1)

    df1 = inputdata.pivot(index='gene_guide', columns='barcode', values='ctrl').add_prefix('ctrl_')
    df2 = inputdata.pivot(index='gene_guide', columns='barcode', values='exp').add_prefix('exp_')

    data_ = pd.merge(df1, df2, left_index=True, right_index=True)
    del df1, df2


    norm = df_normalization(
        data_,
        normalization
    )

    if use_small_count:

        for a in [0.05, 0.1, 0.15]:
            smallcount = df_smallcount(norm, a)
            if len(smallcount) > 0:
                break

        if len(smallcount) == 0:
            smallcount = np.array(
                [data['ctrl'][data['ctrl'] > 0].quantile(0.15)]
            )

        norm = norm + smallcount.mean()
    else:
        norm = norm + 1

    norm = norm.reset_index()
    norm = pd.melt(norm, id_vars='gene_guide', var_name='lfc', value_name='value') 
    norm[['design', 'barcode']] = norm['lfc'].str.split('_', expand=True)
    norm['gene_guide_barcode'] = norm['gene_guide'].str.cat(norm['barcode'], sep='_')
    norm = norm.pivot(index='gene_guide_barcode', columns='design', values='value')

    bar_df = pd.concat(
        [
            inputdata[['gene', 'guide', 'barcode']],
            norm
        ],axis=1, sort=False
    )

    # Step 2: Calculate fold change
    bar_df['fc'] = bar_df['exp'] / bar_df['ctrl']
    bar_df['lfc'] = np.log2(bar_df['fc'])

    del norm



    # ------------------
    # Step 3: Calculate fold change std

    # get valid data range
    lfc_q1 = bar_df['lfc'].quantile(0.025)
    lfc_q2 = bar_df['lfc'].quantile(0.975)
    ctrl_q1 = bar_df['ctrl'].quantile(0.025)
    ctrl_q2 = bar_df['ctrl'].quantile(0.975)

    # select valid data for the model
    model_data = bar_df.loc[
        (bar_df['lfc'] >= lfc_q1) & (bar_df['lfc'] <= lfc_q2) &
        (bar_df['ctrl'] >= ctrl_q1) & (bar_df['ctrl'] <= ctrl_q2)
    ]

    raw_model_data = model_data
    bar_df = bar_df.sort_values(by='ctrl')

    res1 = cal_zfc(raw_model_data,bar_df,'CACT',outprefix=outprefix, frac=0.4, it=10, delta=0.0)
    res2 = cal_zfc(raw_model_data,bar_df,'GCAG',outprefix=outprefix, frac=0.4, it=10, delta=0.0)
    res3 = cal_zfc(raw_model_data,bar_df,'AGCA',outprefix=outprefix, frac=0.4, it=10, delta=0.0)

    bar_df = pd.concat([res1, res2, res3])

    # ------------------
    # Step 4: Calculate raw zlfc using lfc and lfc_std

    bar_df['rank_down'] = bar_df['zlfc'].rank(
        method='average', ascending=True
    ) / len(bar_df['zlfc'])

    bar_df['rank_up'] = bar_df['zlfc'].rank(
        method='average', ascending=False
    ) / len(bar_df['zlfc'])

    bar_df['zlfc'] = bar_df['zlfc'].astype(float)
    bar_df['p'] = stats.norm.cdf(bar_df['zlfc'])

    # Calculate P value using normal distribution
    bar_df['p'] = stats.norm.cdf(bar_df['zlfc'])
    bar_df['p'] = bar_df['p'].map(
        lambda x: x if x <= 0.5 else 1 - x
    )
    bar_df = bar_df.set_index(['gene', 'guide', 'barcode'])

    # ------------------
    # Step 5: Calculate sgRNA mean zscore of fold change

    sg_df = pd.DataFrame(
        {
            'zlfc': bar_df.groupby(level=[0, 1])['zlfc'].mean(),
            'count': bar_df.groupby(level=[0, 1])['zlfc'].count()
        }
    )

    # generate null distribution of sgRNA zlfc
    sg_null_list = dict()
    sg_ecdf_list = dict()

    for i in range(iteration):
        sample_bar = bar_df[['zlfc']].copy()
        sample_bar.loc[:, 'zlfc'] = bar_df['zlfc'].sample(
            frac=1
        ).values
        sg_zlfc = sample_bar.groupby(
            level=[0, 1]
        ).mean()

        for c in sg_df['count'].unique():
            if c not in sg_null_list:
                sg_null_list[c] = list()
            sg_null_list[c].extend(
                sg_zlfc.loc[sg_df['count'] == c, 'zlfc'].values.tolist()
            )
        for c in sg_null_list:
            sg_ecdf_list[c] = ecdf(sg_null_list[c])

    # Calculate p value and FDR of sgRNA zlfc
    sg_df['tp'] = sg_df.apply(
        lambda a: sg_ecdf_list[a['count']](a['zlfc']),
        axis=1
    )
    sg_df['p'] = sg_df['tp'].map(
        lambda a: a if a <= 0.5 else 1 - a
    )
    sg_df['p_adj'] = p_adjust(sg_df['p'], 'BH')
    del sg_df['tp']
    del sg_null_list
    del sg_ecdf_list

    # -------------------
    # Step 6: Calculate Gene mean zscore of fold change

    g_df = pd.DataFrame(
        {
            'zlfc': sg_df.groupby(level=0)['zlfc'].mean(),
            'count': sg_df.groupby(level=0)['zlfc'].count(),
        }
    )

    # generate null distribution of gene zlfc
    g_null_list = dict()
    g_ecdf_list = dict()

    for i in range(iteration):
        sample_sg = sg_df[['zlfc']].copy()
        sample_sg.loc[:, 'zlfc'] = sg_df['zlfc'].sample(
            frac=1
        ).values
        g_zlfc = sample_sg.groupby(
            level=0
        ).mean()

        for c in g_df['count'].unique():
            if c not in g_null_list:
                g_null_list[c] = list()
            g_null_list[c].extend(
                g_zlfc.loc[g_df['count'] == c, 'zlfc'].values.tolist()
            )
        for c in g_null_list:
            g_ecdf_list[c] = ecdf(g_null_list[c])

    # Calculate p value and FDR of gene zlfc
    g_df['tp'] = g_df.apply(
        lambda a: g_ecdf_list[a['count']](a['zlfc']),
        axis=1
    )
    g_df['p'] = g_df['tp'].map(
        lambda a: a if a <= 0.5 else 1 - a
    )
    g_df['p_adj'] = p_adjust(g_df['p'], 'BH')
    del g_df['tp']
    del g_null_list
    del g_ecdf_list

    # ------------------
    # Step 7: Rank aggregation

    # sgRNA Rank
    if top_n_sgrna is None:
        top_n_sgrna = int(
            bar_df.groupby(level=[0, 1])['zlfc'].count().median()
        )
    sg_b_down = bar_df[['rank_down']].copy()
    sg_b_down.loc[:, 'groupid'] = sg_b_down.groupby(
        level=[0, 1]
    )['rank_down'].rank(
        method='first', ascending=True
    ).astype(int)
    sg_b_down.reset_index(drop=False, inplace=True)
    del sg_b_down['barcode']
    sg_b_down.set_index(['gene', 'guide', 'groupid'], inplace=True)
    sg_b_down = sg_b_down.unstack(level=2)
    sg_b_down.columns = sg_b_down.columns.levels[1]
    sg_b_down = sg_b_down[list(range(1, top_n_sgrna + 1))]

    sg_b_up = bar_df[['rank_up']].copy()
    sg_b_up.loc[:, 'groupid'] = sg_b_up.groupby(
        level=[0, 1]
    )['rank_up'].rank(
        method='first', ascending=True
    ).astype(int)
    sg_b_up.reset_index(drop=False, inplace=True)
    del sg_b_up['barcode']
    sg_b_up.set_index(['gene', 'guide', 'groupid'], inplace=True)
    sg_b_up = sg_b_up.unstack(level=2)
    sg_b_up.columns = sg_b_up.columns.levels[1]
    sg_b_up = sg_b_up[list(range(1, top_n_sgrna + 1))]

    sg_df.loc[:, 'RRA_Score_down'] = df_robust_rank_aggregation(sg_b_down)
    sg_df.loc[:, 'RRA_Score_down_adj'] = p_adjust(
        sg_df['RRA_Score_down'], 'BH'
    )
    sg_df.loc[:, 'RRA_Score_up'] = df_robust_rank_aggregation(sg_b_up)
    sg_df.loc[:, 'RRA_Score_up_adj'] = p_adjust(
        sg_df['RRA_Score_up'], 'BH'
    )

    sg_df.loc[:, 'Mean_Rank_down'] = df_mean_rank_aggregation(sg_b_down)
    sg_df.loc[:, 'Mean_Rank_up'] = df_mean_rank_aggregation(sg_b_up)

    # gene Rank
    if top_n_gene is None:
        top_n_gene = int(
            bar_df.groupby(level=0)['zlfc'].count().median()
        )
    g_b_down = bar_df[['rank_down']].copy()
    g_b_down.loc[:, 'groupid'] = g_b_down.groupby(
        level=0
    )['rank_down'].rank(
        method='first', ascending=True
    ).astype(int)
    g_b_down.reset_index(drop=False, inplace=True)
    del g_b_down['barcode']
    del g_b_down['guide']
    g_b_down.set_index(['gene', 'groupid'], inplace=True)
    g_b_down = g_b_down.unstack(level=1)
    g_b_down.columns = g_b_down.columns.levels[1]
    g_b_down = g_b_down[list(range(1, top_n_gene + 1))]

    g_b_up = bar_df[['rank_up']].copy()
    g_b_up.loc[:, 'groupid'] = g_b_up.groupby(
        level=0
    )['rank_up'].rank(
        method='first', ascending=True
    ).astype(int)
    g_b_up.reset_index(drop=False, inplace=True)
    del g_b_up['barcode']
    del g_b_up['guide']
    g_b_up.set_index(['gene', 'groupid'], inplace=True)
    g_b_up = g_b_up.unstack(level=1)
    g_b_up.columns = g_b_up.columns.levels[1]
    g_b_up = g_b_up[list(range(1, top_n_gene + 1))]

    g_df.loc[:, 'RRA_Score_down'] = df_robust_rank_aggregation(g_b_down)
    g_df.loc[:, 'RRA_Score_down_adj'] = p_adjust(
        g_df['RRA_Score_down'], 'BH'
    )
    g_df.loc[:, 'RRA_Score_up'] = df_robust_rank_aggregation(g_b_up)
    g_df.loc[:, 'RRA_Score_up_adj'] = p_adjust(
        g_df['RRA_Score_up'], 'BH'
    )

    g_df.loc[:, 'Mean_Rank_down'] = df_mean_rank_aggregation(g_b_down)
    g_df.loc[:, 'Mean_Rank_up'] = df_mean_rank_aggregation(g_b_up)

    return (bar_df, sg_df, g_df)
