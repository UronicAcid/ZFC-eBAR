import numpy as np
import statsmodels.api as sm
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
import pandas as pd
def lowess_regression(x, y, frac=0.4, it=10, delta=0.0):
    """
    Perform LOWESS regression and return interpolation function.
    :param x: Control mean values (independent variable)
    :param y: LFC standard deviation (dependent variable)
    :param frac: Smoothing factor for LOWESS
    :param it: Number of iterations for robust estimation
    :param delta: A distance over which to combine points
    """
    lowess = sm.nonparametric.lowess
    lowess_regression = lowess(y, x, frac=frac, it=it, delta=delta)
    
    x_ = lowess_regression[:, 0]
    yest = lowess_regression[:, 1]

    # Interpolation for points within the range
    f = interp1d(x_, yest)

    # Extrapolation slopes and intercepts
    slope1 = (yest[1] - yest[0]) / (x_[1] - x_[0])
    intercept1 = yest[0] - slope1 * x_[0]
    slope2 = (yest[-1] - yest[-2]) / (x_[-1] - x_[-2])
    intercept2 = yest[-1] - slope2 * x_[-1]

    return x_, yest, f, slope1, intercept1, slope2, intercept2

def fit_qc_plot(bar_df_split, train_data, x_, yest, barcode, display=True):
    """
    Plot LOWESS regression results, raw data, histograms, and save to PDF.
    :param bar_df_split: DataFrame with computed values for plotting
    :param train_data: DataFrame with control means and standard deviations
    :param x_: X-values from LOWESS regression
    :param yest: Estimated Y-values from LOWESS regression
    :param barcode: Barcode label for the plot
    :param display: Whether to display the plots in Jupyter-Lab
    """
    fig, axs = plt.subplots(2, 3, figsize=(15, 10))

    # Plot 1: Scatter plot of raw data
    axs[0, 0].scatter(train_data['ctrlmean'], train_data['lfcstd'], s=6, color='gray')
    axs[0, 0].plot(x_, yest, color='black')
    axs[0, 0].set_title(f'{barcode} - Raw Data with LOWESS Fit')
    axs[0, 0].set_xlabel('Control Mean')
    axs[0, 0].set_ylabel('LFC Std')

    # Plot 2: Histogram of log(ctrl)
    log_ctrl = np.log(bar_df_split['ctrl'])
    mean_log_ctrl = np.mean(log_ctrl)
    median_log_ctrl = np.median(log_ctrl)
    axs[0, 1].hist(log_ctrl, bins=50, color='lightgray')
    axs[0, 1].axvline(mean_log_ctrl, color='black', label=f'Mean: {mean_log_ctrl:.2f}')
    axs[0, 1].axvline(median_log_ctrl, color='red', label=f'Median: {median_log_ctrl:.2f}')
    axs[0, 1].set_xlabel('log(ctrl)')
    axs[0, 1].set_ylabel('Density')
    # Annotating mean and median values on the plot
    axs[0, 1].annotate(f'Mean: {mean_log_ctrl:.2f}', xy=(mean_log_ctrl, 0.02), xytext=(mean_log_ctrl + 0.5, 0.03))
    axs[0, 1].annotate(f'Median: {median_log_ctrl:.2f}', xy=(median_log_ctrl, 0.02), xytext=(median_log_ctrl - 1, 0.03))

    axs[0, 1].legend()

    # Plot 3: LFC vs log(ctrl)
    axs[0, 2].scatter(np.log(bar_df_split['ctrl']), bar_df_split['lfc'], alpha=0.1, s=1, color='gray')
    axs[0, 2].scatter(np.log(bar_df_split['ctrl']), 2 * bar_df_split['lfc_std'], s=1, c="red")
    axs[0, 2].scatter(np.log(bar_df_split['ctrl']), -2 * bar_df_split['lfc_std'], s=1, c="red")
    axs[0, 2].axhline(0, color='black', linestyle='--')
    axs[0, 2].set_xlabel('log(ctrl)')
    axs[0, 2].set_ylabel('lfc')

    # Plot 4: zlfc hist
    zlfc = bar_df_split['zlfc']
    mean_zlfc = np.mean(zlfc)
    median_zlfc = np.median(zlfc)
    axs[1, 0].hist(zlfc, bins=50, color='lightgray')
    axs[1, 0].axvline(mean_zlfc, color='black', label=f'Mean: {mean_zlfc:.2f}')
    axs[1, 0].axvline(median_zlfc, color='red', label=f'Median: {median_zlfc:.2f}')
    axs[1, 0].set_xlabel('zlfc')
    axs[1, 0].set_ylabel('Density')
    # Annotating mean and median values on the plot
    axs[1, 0].annotate(f'Mean: {mean_zlfc:.2f}', xy=(mean_zlfc, 0.02), xytext=(mean_zlfc + 0.5, 0.03))
    axs[1, 0].annotate(f'Median: {median_zlfc:.2f}', xy=(median_zlfc, 0.02), xytext=(median_zlfc - 1, 0.03))

    axs[1, 0].legend()
    
    
    
    # Plot 5: 
    axs[1, 1].scatter(np.log(bar_df_split.loc[:, 'ctrl']),bar_df_split.loc[:, 'lfc'],alpha=1,s=1,color='gray')
    axs[1, 1].axhline(y=0, color='r', linestyle='--')
    axs[1, 1].set_xlabel('log(ctrl)')
    axs[1, 1].set_ylabel('lfc')

    # Plot 6: LFC Std against control mean (for comparison)    
    axs[1, 2].scatter(np.log(bar_df_split.loc[:, 'ctrl']),bar_df_split.loc[:, 'zlfc'],alpha=1,s=1)
    axs[1, 2].axhline(y=0, color='r', linestyle='--')
    axs[1, 2].set_xlabel('log(ctrl)')
    axs[1, 2].set_ylabel('zlfc')

    plt.tight_layout()

    # Save figure
    #plt.savefig(os.path.join(os.getcwd(), f'{barcode}_lowess_results.png'))

    # Display if required
    if display:
        plt.show()

    return fig

def cal_zfc(model_data,bar_df, barcode, frac=0.4, it=10, delta=0.0):
    # Filter raw data for the specified barcode
    model_data = model_data[model_data['barcode'] == barcode]

    # Discretize the 'ctrl' data into 500 bins
    bins = pd.cut(model_data['ctrl'], 500).astype(str)

    # Create DataFrame for mean and standard deviation by bin
    train_data = pd.DataFrame({
        'lfcstd': model_data.groupby(bins)['lfc'].std(),
        'ctrlmean': model_data.groupby(bins)['ctrl'].mean()
    })

    # Remove rows with missing values in 'lfcstd'
    train_data = train_data.loc[~train_data['lfcstd'].isnull()]

    # Perform LOWESS regression
    x_, yest, f, slope1, intercept1, slope2, intercept2 = lowess_regression(train_data['ctrlmean'].values, train_data['lfcstd'].values, frac, it, delta)

    # Cutoff calculation
    bar_df_split = bar_df.loc[bar_df['barcode'] == barcode, :]
    cutoff = min(np.std(bar_df_split['lfc'][-100:]), yest[-1])

    # Interpolation function to estimate lfc_std for given ctrl values
    def interpolate(x_new, y_cutoff):
        if x_new < x_[0]:
            y_new = slope1 * x_new + intercept1
        elif x_new > x_[-1]:
            y_new = slope2 * x_new + intercept2
            if y_new < y_cutoff:
                y_new = y_cutoff
        else:
            y_new = f(x_new)
        return y_new

    # Apply interpolation to calculate 'lfc_std'
    bar_df_split.loc[:, 'lfc_std']= bar_df_split['ctrl'].apply(lambda x: interpolate(x_new=x, y_cutoff=cutoff))

    # Calculate zlfc
    bar_df_split.loc[:, 'zlfc'] = bar_df_split['lfc'] / bar_df_split['lfc_std']

    # Plot and save results
    fig = fit_qc_plot(bar_df_split, train_data, x_, yest, barcode, display=False)
    fig.savefig('_'.join([barcode, 'qc.png']))

    return bar_df_split

