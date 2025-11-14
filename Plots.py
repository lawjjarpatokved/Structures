import numpy as np
import matplotlib.pyplot as plt
import os
import seaborn as sns

def line_plot(
    x, y,
    xlabel='', ylabel='', title='',
    filename=None,
    ax=None, label=None, show=False,
    color=None, linestyle='--', marker='.', linewidth=2, markersize=2,
    grid=True, grid_style=':', grid_alpha=0.6,
    legend=True, legend_loc='best',
    figsize=(6, 5),
    title_fontsize=14, label_fontsize=12, tick_fontsize=10, legend_fontsize=10,
    tight_layout=True,
    facecolor='white',
    fontfamily='Times New Roman',
    padding_fraction=0.05, 
):

    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    m = np.isfinite(x) & np.isfinite(y)
    if not np.any(m):
        print(f"[warn] Nothing to plot for '{title}' (no finite points).")
        return ax

    # Create new figure/axes if not provided
    if ax is None:
        fig, ax = plt.subplots(figsize=figsize, facecolor=facecolor)

    if color is None:
        color = sns.color_palette("pastel", 1)[0]

    # --- Plot line ---
    ax.plot(
        x[m], y[m],
        linestyle=linestyle, color=color,
        marker=marker, linewidth=linewidth, markersize=markersize,
        label=label
    )

    # --- Labels and title ---
    ax.set_xlabel(xlabel, fontsize=label_fontsize, fontfamily=fontfamily)
    ax.set_ylabel(ylabel, fontsize=label_fontsize, fontfamily=fontfamily)
    ax.set_title(title, fontsize=title_fontsize, fontfamily=fontfamily, fontweight='normal')

    ax.relim()
    ax.autoscale_view()
    ax.margins(x=padding_fraction, y=padding_fraction)

    # --- Grid, legend, style ---
    if grid:
        ax.grid(True, which='both', linestyle=grid_style, alpha=grid_alpha)

    if legend and label is not None:
        ax.legend(loc=legend_loc, fontsize=legend_fontsize, frameon=False)

    # --- Tick settings ---
    ax.tick_params(axis='both', labelsize=tick_fontsize)
    for tick_label in ax.get_xticklabels() + ax.get_yticklabels():
        tick_label.set_fontfamily(fontfamily)

    if tight_layout:
        plt.tight_layout()

    if filename:
        plt.savefig(filename, dpi=300, bbox_inches='tight', facecolor=facecolor)

    if show:
        plt.show()

    plt.close()
    return ax


def plot_single_bar(
    data, filename,
    x_labels=None, title=None, xlabel=None, ylabel=None,
    colors=None, edgecolor='black', alpha=0.8,
    width=0.3, spacing=0.7,
    value_labels=True, value_format='{:.3f}', value_fontsize=10, value_color='black', value_offset=0.02,
    grid=False, grid_style='--', grid_alpha=0.5,
    figsize=(6, 6),
    title_fontsize=14, label_fontsize=12, tick_fontsize=10, legend_fontsize=10,
    show=False, tight_layout=True,
    fontfamily='Times New Roman'
                        ):


    n = len(data)
    x = np.arange(n) * (width + spacing)

    if colors is None:
        colors = sns.color_palette("tab10", n)

    fig, ax = plt.subplots(figsize=figsize, facecolor='white')

    # --- Bars ---
    bars = []
    for i, val in enumerate(data):
        b = ax.bar(
            x[i], val,
            width=width,
            color=colors[i % len(colors)],
            edgecolor=edgecolor,
            alpha=alpha,
            label=x_labels[i] if x_labels else f'Bar {i + 1}'
        )
        bars.append(b)

        # --- Value label ---
        if value_labels:
            ax.text(
                x[i],
                val + value_offset * max(data),
                value_format.format(val),
                ha='center', va='bottom',
                fontsize=value_fontsize,
                color=value_color,
                fontfamily=fontfamily
            )

    # --- Labels and title ---
    if xlabel:
        ax.set_xlabel(xlabel, fontsize=label_fontsize, fontfamily=fontfamily)
    if ylabel:
        ax.set_ylabel(ylabel, fontsize=label_fontsize, fontfamily=fontfamily)
    if title:
        ax.set_title(title, fontsize=title_fontsize, fontfamily=fontfamily, fontweight='normal')

    # --- Grid ---
    if grid:
        ax.grid(axis='y', linestyle=grid_style, alpha=grid_alpha)

    ymax = max(data) * 1.25
    ax.set_ylim(0, ymax)
    ax.set_xticks([])
    ax.set_xticklabels([])

    # --- Legend ---
    ax.legend(
        loc='upper right',
        fontsize=legend_fontsize,
        frameon=False
    )

    for spine in ['top', 'right']:
        ax.spines[spine].set_visible(False)

    # --- Tick font ---
    ax.tick_params(axis='both', labelsize=tick_fontsize)
    for tick_label in ax.get_yticklabels():
        tick_label.set_fontfamily(fontfamily)

    if tight_layout:
        plt.tight_layout()

    plt.savefig(filename, dpi=300, bbox_inches='tight', facecolor='white')
    if show:
        plt.show()
    plt.close(fig)


def plot_grouped_bars(*datasets, labels=None, x_labels=None, title=None, ylabel=None, width=0.15):

    """
    Plots a grouped bar chart for any number of input lists/arrays.

    Parameters
    ----------
    *datasets : list of lists or arrays
        Any number of iterable datasets of equal length.
        Example: plot_grouped_bars(data1, data2, data3)
    labels : list of str, optional
        Legend labels for each dataset.
    x_labels : list of str, optional
        Labels for the x-axis categories.
    title : str, optional
        Title of the plot.
    ylabel : str, optional
        Label for the y-axis.
    width : float, optional
        Bar width (default 0.15).
    """

    num_groups = len(datasets[0])
    num_datasets = len(datasets)

    x = np.arange(num_groups)
    offsets = np.linspace(-width * (num_datasets - 1) / 2,
                          width * (num_datasets - 1) / 2, num_datasets)

    plt.figure(figsize=(8, 5))

    for i, data in enumerate(datasets):
        label = labels[i] if labels else f'Data {i+1}'
        plt.bar(x + offsets[i], data, width=width, label=label)

    if x_labels:
        plt.xticks(x, x_labels)

    if ylabel:
        plt.ylabel(ylabel)
    if title:
        plt.title(title)

    plt.legend()
    plt.tight_layout()
    plt.show()



