
import matplotlib.pyplot as plt
import numpy as np
import math
import pandas as pd


def remove_duplicate_motifs(combined):
    to_remove = []
    for i, motif_considered in enumerate(combined):
        for j, other_motif in enumerate(combined):
            if motif_considered != other_motif and motif_considered.words == other_motif.words and i < j:
                motif_considered.alignments = motif_considered.alignments.union(other_motif.alignments)
                to_remove.append(other_motif)

    result = [c for c in combined if c not in to_remove]
    return result


def get_normalized_sequence_for_window(window):
    sequence = window.sequence[window.start_index:window.end_index+1]
    return [float(x) for x in window.sax.whiten(sequence)]

def get_normalized_sequences_for_motif(motif):
    result = []
    for alignment in list(motif.alignments):
        result.append(get_normalized_sequence_for_window(alignment.window_a))
        result.append(get_normalized_sequence_for_window(alignment.window_b))
    return result

def plot_motif_window(window_height_lower, window_height_upper, window, color, alpha):
    plt.axvspan(window.start_index, window.end_index, window_height_lower, window_height_upper,
                color=color, alpha=alpha)



def plot_project_sequence(project_sequence, title=None, figsize=(16, 7), title_font_size=20,
                     axis_font_size=18, axis_x_dates=False,
                     xlabel="week", ylabel="#defects",
                     xtick_interval=50, ytick_interval=10, between_coverage=""):

    fig = plt.figure()
    project_sequence.sequence.plot(title=title, figsize=figsize, grid=True, linewidth=2, color="black")



    #print(fig.axes[0].axes.title)
    fig.axes[0].axes.title.set_fontsize(title_font_size)
    #fig.axes[0].axes.label.set_fontsize(axis_font_size)
    #fig.axes[0].axes.label.set_fontsize(axis_font_size)

    if title is not None:
        title = title + "; Between cov. = {}".format(between_coverage)
        plt.title(title)
    if xlabel is not None:
        plt.xlabel(xlabel)
    if ylabel is not None:
        plt.ylabel(ylabel)

    if axis_x_dates:
        plt.xticks(rotation=90)
    else:
        plt.xticks([x for x in range(0, len(project_sequence.sequence), xtick_interval)],
                   [x for x in range(project_sequence.start_index(),
                                     project_sequence.start_index()+project_sequence.length(), xtick_interval)])

    max_value = int(np.max(project_sequence.sequence))

    fig.axes[0].set_ylim(ymin=0, ymax=max_value)
    plt.yticks([x for x in np.arange(0, max_value, ytick_interval)],
               [int(x) for x in np.arange(0, max_value, ytick_interval)])

    return fig



def plot_time_series(ts, title=None, figsize=(16, 7), title_font_size=20,
                     axis_font_size=18, axis_x_dates=False,
                     xlabel="week", ylabel="#defects",
                     xtick_interval=50):
    fig = plt.figure()
    ts.plot(figsize=figsize, title=title, grid=True, linewidth=2,
            xticks=range(0, len(ts), xtick_interval), color="black")

    fig.axes[0].title.set_fontsize(title_font_size)
    fig.axes[0].xaxis.label.set_fontsize(axis_font_size)
    fig.axes[0].yaxis.label.set_fontsize(axis_font_size)

    if title is not None:
        plt.title(title)
    if xlabel is not None:
        plt.xlabel(xlabel)
    if ylabel is not None:
        plt.ylabel(ylabel)

    if axis_x_dates:
        plt.xticks(rotation=90)
    else:
        plt.xticks([x for x in range(1, len(ts), xtick_interval)], [x for x in range(1, len(ts), xtick_interval)])

    return fig


def plot_window(window, color, alpha):
    plt.axvspan(window.start_index, window.end_index,
                    color=color, alpha=alpha)

def overlay_alignments(fig, ts, alignments, color, alpha=0.05):

    # determine if the sequence is a or b
    alignment = alignments[0]
    window_a_plot = False
    window_b_plot = False
    if alignment.window_a.sequence.equals(ts):
        window_a_plot = True
    if alignment.window_b.sequence.equals(ts):
        window_b_plot = True

    for alignment in alignments:
        if window_a_plot:
            plot_window(alignment.window_a, color, alpha)
        if window_b_plot:
            plot_window(alignment.window_b, color, alpha)

    return fig


def overlay_weeks(fig, weeks, color, alpha=0.2, placement=None):
    for week in weeks:
        if placement == 'b':
            plt.axvspan(week, week + 1, 0, 0.95, facecolor=color, alpha=alpha)
        elif placement == 't':
            plt.axvspan(week, week + 1, 0.05, 1, facecolor=color, alpha=alpha)
        else:
            plt.axvspan(week, week + 1, facecolor=color, alpha=alpha)
    return fig


def overlay_motifs(project_sequence, motifs, motifs_all=None, pos=0, color="gray", alpha=0.35):
    if motifs_all is None:
        motifs_all = motifs

    number_of_motifs = len(motifs_all)
    window_height = 1.0 / number_of_motifs
    window_height_value = max(project_sequence.sequence) / number_of_motifs

    for i, motif in enumerate(motifs_all):
        window_height_lower = i * window_height
        window_height_upper = window_height_lower + window_height
        if pos == 1:
            window_height_lower = i * window_height
            window_height_upper = window_height_lower + (window_height / 2)
        elif pos == 2:
            window_height_lower = i * window_height + (window_height / 2)
            window_height_upper = window_height_lower + (window_height / 2)

        window_height_value_lower = i * window_height_value
        window_height_value_upper = window_height_value_lower + window_height_value

        if motif in motifs:
            for align in list(motif.alignments):
                if align.window_a.sequence.equals(project_sequence.sequence):
                    plot_motif_window(window_height_lower, window_height_upper, align.window_a, color, alpha)
                else:
                    plot_motif_window(window_height_lower, window_height_upper, align.window_b, color, alpha)

        motif_name = str(motif)
        if len(motif_name) > 64:
            motif_name = motif_name[0:63]+"..."
        plt.text(4, (window_height_value_lower + window_height_value_upper) // 2, motif_name)


def plot_motif_window(window_height_lower, window_height_upper, window, color, alpha):
    plt.axvspan(window.start_index, window.end_index, window_height_lower, window_height_upper,
                color=color, alpha=alpha)


def plot_project_sequence_with_motifs(project_sequence, motifs_within=[], motifs_between=[],
                                      between_color='pink', within_color='lightgreen',
                                      figsize=(16, 7),
                                      title=None, xtick_interval=24, axis_x_dates=False,
                                      between_coverage="",  ylabel="#defects"):
    combined_motifs = remove_duplicate_motifs(motifs_within + motifs_between)

    number_of_motifs = len(combined_motifs)

    if number_of_motifs == 0:
        f = plot_project_sequence(project_sequence, between_coverage=between_coverage,
                                                            xtick_interval=xtick_interval,
                                                            axis_x_dates=axis_x_dates,
                                                            figsize=figsize,
                                                            title=title,  ylabel=ylabel)
        return

    window_height = 1.0 / number_of_motifs
    window_height_value = max(project_sequence.sequence) / number_of_motifs

    f = plot_project_sequence(project_sequence, between_coverage=between_coverage,
                                                        xtick_interval=xtick_interval,
                                                        ytick_interval=window_height_value,
                                                        axis_x_dates=axis_x_dates,
                                                        figsize=figsize,
                                                        title=title, ylabel=ylabel)

    motifs_per_height_window = 0
    motifs_per_height_window = motifs_per_height_window + 1 if len(motifs_within) > 0 else motifs_per_height_window
    motifs_per_height_window = motifs_per_height_window + 1 if len(motifs_between) > 0 else motifs_per_height_window
    pos = 0

    if len(motifs_within) > 0:
        if motifs_per_height_window > 1:
            pos = 1
        overlay_motifs(project_sequence, motifs_within, motifs_all=combined_motifs,
                       pos=pos, color=within_color)

    if len(motifs_between) > 0:
        if motifs_per_height_window > 1:
            pos = 2
        overlay_motifs(project_sequence, motifs_between, motifs_all=combined_motifs,
                       pos=pos, color=between_color)
    return f

def get_colors_for_sequences_in_motif(motif, color_a="Blues", color_b="Oranges"):
    result = []
    color_a_map = plt.get_cmap(color_a)
    color_b_map = plt.get_cmap(color_b)
    for i, alignment in enumerate(list(motif.alignments)):
        result.append(color_a_map(i*50+140))
        result.append(color_b_map(i*50+140))
    return result

def plot_motif_alignments(motif, title=None, figsize=(10, 5), xlabel='week', xtick_div_by=10, ax=None):
    normalized_sequences = get_normalized_sequences_for_motif(motif)
    end_index = len(normalized_sequences[0])
    xtick_step = math.ceil((end_index+1) / xtick_div_by)
    normalized_sequences_df = pd.DataFrame(normalized_sequences).T
    result = normalized_sequences_df.plot(legend=None, grid=True,
                                          xticks=range(0, end_index, xtick_step), title=title, figsize=figsize, ax=ax,
                                          color=get_colors_for_sequences_in_motif(motif))
    plt.ylabel(r'#defects (normalized $(y-\mu)$ / $\sigma$)')
    plt.xlabel(xlabel)


def plot_all_motif_alignments(motifs_within, motifs_between, columns=3, width=5, height=3):
    unique_motifs = remove_duplicate_motifs(motifs_within + motifs_between)

    if len(unique_motifs) == 0:
        return

    figsize = (width * columns, height * math.ceil(len(unique_motifs) / columns))

    #fig = plt.figure(figsize=figsize)
    fig, axs = plt.subplots(math.ceil(len(unique_motifs) / columns), columns, figsize=figsize)
    many_rows = type(axs[0]) is np.ndarray
    j = 0
    for i, motif in enumerate(unique_motifs):
        motif_name = str(motif)
        if len(motif_name) > 40:
            motif_name = motif_name[0:39]+"..."
        if many_rows:
            plot_motif_alignments(motif, figsize=None, ax=axs[i // columns][j], title=motif_name)
        else:
            plot_motif_alignments(motif, figsize=None, ax=axs[j], title=motif_name)
        j = (j + 1) % columns
    plt.tight_layout()


def plot_subplot_motif(subplot, motif, project_sequence, xtick_interval=4, 
                        axis_x_dates=True, ytick_interval=10,
                        title ="",  between_color = "lightgreen"):
    

    subplot.grid(color='lightgray', linestyle='-', linewidth=1)
    plt.title(f'Motif: {str(motif)}, {title}')
    plt.xlabel("week")
    plt.ylabel("#defects")

    if axis_x_dates:
        plt.xticks([x for x in range(0, len(project_sequence.sequence), xtick_interval)], 
                   [project_sequence.sequence.index[x] for x in range(project_sequence.start_index(),
                                     project_sequence.start_index()+project_sequence.length(), xtick_interval)],
                   rotation=90)
    else:
        plt.xticks([x for x in range(0, len(project_sequence.sequence), xtick_interval)],
                   [x for x in range(project_sequence.start_index(),
                                     project_sequence.start_index()+project_sequence.length(), xtick_interval)])

    max_value = int(np.max(project_sequence.sequence))
    subplot.set_ylim(ymin=0, ymax=max_value)

    plt.plot(project_sequence.sequence, linewidth=1.5, color="gray")
    
    for align in list(motif.alignments):
        if align.window_a.sequence.equals(project_sequence.sequence):
            plt.plot(align.window_a.sequence[align.window_a.start_index:align.window_a.end_index], 
                                             color='black', alpha=1, linewidth=4)
        else:
            plt.plot(align.window_b.sequence[align.window_b.start_index:align.window_b.end_index],
                                             color='black', alpha=1, linewidth=4)

    
def plot_comparisons_motif(figsize, project_a_seq, project_b_seq, title_a, title_b, motif):
    fig = plt.figure(figsize=figsize)

    # plot top
    subplot_a = plt.subplot(2, 1, 1)
    plot_subplot_motif(subplot_a, motif, project_a_seq, xtick_interval=4, 
                            axis_x_dates=True, ytick_interval=10,
                            title =title_a,  between_color = "lightgreen")


    subplot_b = plt.subplot(2, 1, 2)
    plot_subplot_motif(subplot_b, motif, project_b_seq, xtick_interval=4, 
                            axis_x_dates=True, ytick_interval=10,
                            title =title_b,  between_color = "lightgreen")
    return fig
    