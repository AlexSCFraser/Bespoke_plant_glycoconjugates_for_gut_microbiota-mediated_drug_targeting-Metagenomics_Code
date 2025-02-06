import logging
import os

import matplotlib
import seaborn
from matplotlib import pyplot as plt

from GraphText import adjusted_labels


def render_coverage_heatmap_xyloglucan_old(alignment_threshold, coverage_attributes, graphics_output_folder,
                                       dataframe, name_prefix, name_title_insert, pul_names, max_name=["Maximum\ncoverage"]):
    fig, ax = plt.subplots()
    create_coverage_heatmap_axis(alignment_threshold, coverage_attributes, graphics_output_folder, dataframe,
                                 name_prefix, name_title_insert, pul_names, ax)

    ax.set(ylabel="One row per patient")
    ax.figure.set_dpi(500)
    ax.figure.suptitle(f"{name_title_insert} Coverage Heatmap")
    matplotlib.pyplot.title(f"{alignment_threshold * 100: 2.0f}% percent identity threshold")
    plt.savefig(os.path.join(graphics_output_folder, f"{name_prefix}xygul_coverage_heatmap.pdf"), dpi=500)
    plt.savefig(os.path.join(graphics_output_folder, f"{name_prefix}xygul_coverage_heatmap.png"), dpi=500)
    matplotlib.pyplot.show()


def render_coverage_heatmap_xyloglucan(alignment_threshold, coverage_attributes, graphics_output_folder, dataframe,
                                       name_prefix, name_title_insert, pul_names, subtitle, title,
                                       percent_identity_threshold, xlabel=None, species_name=None,
                                       vertical_group_labels=False, no_data_colour="gainsboro",
                                       yticks=None, ytick_labels=None, row_label="Individual", xtick_fontsize=None,
                                       cell_dim: (float, float) = None, linewidth=None, cbar_fraction=0.1):

    marg_top = 0.7
    marg_bottom = 0.7
    marg_left = 0.7
    marg_right = 0.7

    # fig.figsize = (fig_width, fig_height)
    # fig = plt.figure(figsize=(fig_width, fig_height))
    if linewidth is None:
        linewidth = 1.5 if dataframe.shape[0] < 20 else 0.2 if dataframe.shape[0] < 100 else 0.2 if dataframe.shape[0] < 150 else 0
    if xtick_fontsize is None:
        xtick_fontsize = 5 if len(pul_names) > 8 else 6 if len(pul_names) > 5 else 7
    try:
        num_categories = len(dataframe.index.levels[0])
    except AttributeError:
        num_categories = 1
    if num_categories > 1:
        slices = []
        lengths = []
        for item in dataframe.index.levels[0]:
            try:
                slices.append(dataframe.index.get_loc(item))
                try:
                    lengths.append(len([pos for pos in dataframe.index.get_loc(item) if pos]))
                except TypeError:
                    lengths.append(slices[-1].stop - slices[-1].start)
            except KeyError:
                num_categories -= 1
                pass
        subtitle_height = sum(lengths)/30
        fig, axes = plt.subplots(num_categories+1, 1, height_ratios=[subtitle_height] + lengths, layout="compressed")
        # fig = seaborn.FacetGrid(dataframe_plot_presence, row='phenotype')
        # fig.map_dataframe(create_presence_heatmap_axis, count_name, groups, is_binary, linewidth,
        #                                               no_data_colour, pul_names, xlabel)
        fig.supylabel(row_label)
        axes[0].set_title(subtitle)
        axes[0].axis("off")
        for i, slice in enumerate(slices):
            ylabel = dataframe.index.levels[0][i]
            create_coverage_heatmap_axis(alignment_threshold, coverage_attributes, graphics_output_folder,
                                         dataframe[slice], name_prefix, name_title_insert, pul_names, axes[i + 1],
                                         ylabel=ylabel, no_data_colour=no_data_colour, cbar_fraction=cbar_fraction,
                                         linewidth=linewidth, vertical_ylabel=vertical_group_labels,
                                         tick_font_size=xtick_fontsize)

        for axis in axes[1:-1]:
            axis.collections[0].colorbar.remove()
            axis.tick_params(top=False, bottom=False, left=False, right=False,
                             labelleft=False, labelbottom=False)
            axis.set_xlabel(None)
        pass
    else:
        if cell_dim:
            cell_height, cell_width = cell_dim
            marg_top = 0.7
            marg_bottom = 1.25
            marg_left = 0.4
            marg_right = 0.4

            # determine figure dimensions
            ax_height = cell_height * dataframe.shape[0]
            ax_width = cell_width * dataframe.shape[1]

            fig_height = ax_height + marg_top + marg_bottom
            fig_width = ax_width + marg_left + marg_right
            l, b, w, h = marg_left / fig_width, marg_bottom / fig_height, ax_width / fig_width, ax_height / fig_height
            fig, axes = plt.subplots(figsize=(fig_width, fig_height), nrows=2)
            ax = axes[0]
            create_coverage_heatmap_axis(alignment_threshold, coverage_attributes, graphics_output_folder, dataframe,
                                         name_prefix, name_title_insert, pul_names, ax, ylabel=row_label,
                                         no_data_colour=no_data_colour, cbar=False, cbar_fraction=0.02, linewidth=linewidth,
                                         yticks=yticks, ytick_labels=ytick_labels, tick_font_size=xtick_fontsize,
                                         vertical_ylabel=vertical_group_labels)
            ax.set_position([l, b, w, h])
            cax = axes[1]
            cbar_bottom = marg_bottom / fig_height / 2
            cbar_height = marg_bottom / fig_height / 15
            cbar_left = marg_left / fig_width
            cbar_width = ax_width / fig_width

            cax.set_position([cbar_left, cbar_bottom, cbar_width, cbar_height])
            colorbar = fig.colorbar(ax.collections[0], cax=cax, orientation='horizontal',
                                    label=f"{name_title_insert} coverage", format="%1.1f")
            colorbar.outline.set_color('black')
            colorbar.outline.set_linewidth(0.75)
            colorbar.ax.xaxis.set_tick_params(width=0.75)
        else:
            fig, ax = plt.subplots()
            create_coverage_heatmap_axis(alignment_threshold, coverage_attributes, graphics_output_folder, dataframe,
                                         name_prefix, name_title_insert, pul_names, ax, ylabel=row_label,
                                         no_data_colour=no_data_colour, cbar_fraction=cbar_fraction, linewidth=linewidth,
                                         yticks=yticks, ytick_labels=ytick_labels, tick_font_size=xtick_fontsize,
                                         vertical_ylabel=vertical_group_labels)
        ax.set_title(subtitle)

    fig.set_dpi(500)
    full_title = f"{f'{name_prefix} ' if name_prefix else ''}{title}"
    fig.suptitle(full_title)

    plt.savefig(os.path.join(graphics_output_folder, f"{full_title} xygul coverage heatmap {percent_identity_threshold:.2f}_perc_ident.eps"), dpi=500)
    plt.savefig(os.path.join(graphics_output_folder, f"{full_title} xygul coverage heatmap {percent_identity_threshold:.2f}_perc_ident.png"), dpi=500)
    plt.savefig(os.path.join(graphics_output_folder, f"{full_title} xygul coverage heatmap {percent_identity_threshold:.2f}_perc_ident.svg"), dpi=500)
    try:
        matplotlib.pyplot.show()
    except AttributeError:
        logging.getLogger().exception("Error showing image to user")
    pass


def create_coverage_heatmap_axis(alignment_threshold, coverage_attributes, graphics_output_folder, dataframe,
                                 name_prefix, name_title_insert, column_names, ax, xlabel=None, ylabel=None,
                                 no_data_colour="gainsboro", cbar_fraction=0.05, linewidth=0, vertical_ylabel=False,
                                 yticks=None, ytick_labels=None, tick_font_size=7, cbar=True):

    seaborn.heatmap(data=dataframe, robust=True,
                    cmap="mako",
                    cbar=cbar,
                    cbar_kws={"location": "bottom", "fraction": cbar_fraction, "label": f"{name_title_insert} coverage",
                              "format": "%1.1f"},
                    yticklabels=[],
                    ax=ax,
                    linewidth=linewidth,
                    linecolor='black',
                    vmin=0,
                    vmax=1
                    )
    if cbar:
        ax.collections[0].colorbar.outline.set_color('black')
        ax.collections[0].colorbar.outline.set_linewidth(0.5)
        ax.collections[0].colorbar.ax.xaxis.set_tick_params(width=0.5)

    ax.set_xticklabels(column_names)
    pul_name_angle, x_labels = adjusted_labels(ax.get_xticklabels(), name_title_insert)
    ax.set_xticklabels(x_labels, rotation=pul_name_angle, fontsize=tick_font_size)


    if vertical_ylabel:
        # this was 0 before, seemed backwards?
        ylabel_rot = 90
    else:
        ylabel_rot = 0
    ax.set_facecolor(no_data_colour)

    if yticks is not None and ytick_labels is not None:
        ax.set_yticks(yticks, labels=ytick_labels, fontsize=tick_font_size, rotation=ylabel_rot)

    if ylabel:
        ax.set_ylabel(ylabel=ylabel, rotation=ylabel_rot, labelpad=0.0)
    if xlabel:
        ax.set(xlabel=xlabel)

    # Move xticklabels up slightly
    if pul_name_angle == 45:
        ax.tick_params(axis='x', which='both', pad=-4)

    return ax


def render_presence_heatmap_xyloglucan(dataframe_plot_presence, subtitle, title, percent_identity_threshold,
                                       coverage_threshold, name_prefix=None, graphics_output_folder=None,
                                       is_binary=False, xlabel=None, species_name=None, no_data_colour="gainsboro",
                                       row_annotations=None, cbar_label=None, yticks=None, ytick_labels=None,):
    # if metahit and trying to make control comparison figure
    # pul_names = ['B. uniformis\nATCC 8492\nxyGUL 1', 'B. uniformis\nATCC 8492\nxyGUL 2', 'B. ovatus \nATCC 8483\nxyGUL', 'B. cellulosilyticus\nDSM 14838\nxyGUL']
    # dataframe_plot_presence = dataframe_plot_presence.reindex(columns=(["xygul_count"] + [item.replace("ATCC 8492\n", "").replace(" \nATCC 8483\nxyGUL", "").replace("B. cellulosilyticus\nDSM 14838\nxyGUL", "B.\ncellulosilyticus") + "\npresent" for item in pul_names]))
    # name_prefix += " comparison"
    linewidth = 1.5 if dataframe_plot_presence.shape[0] < 20 else 0.2 if dataframe_plot_presence.shape[0] < 100 else 0.2 if dataframe_plot_presence.shape[0] < 150 else 0
    xtick_fontsize = 5 if dataframe_plot_presence.shape[1] > 8 else 7
    try:
        num_categories = len(dataframe_plot_presence.index.levels[0])
    except AttributeError:
        num_categories = 1
    if num_categories > 1:
        slices = []
        lengths = []
        for item in dataframe_plot_presence.index.levels[0]:
            try:
                slices.append(dataframe_plot_presence.index.get_loc(item))
                try:
                    lengths.append(len([pos for pos in dataframe_plot_presence.index.get_loc(item) if pos]))
                except TypeError:
                    lengths.append(slices[-1].stop - slices[-1].start)
            except KeyError:
                num_categories -= 1
                pass
        subtitle_height = sum(lengths)/30
        fig, axes = plt.subplots(num_categories+1, 1, height_ratios=[subtitle_height] + lengths, layout="compressed")
        fig.supylabel("Individual")
        axes[0].set_title(subtitle)
        axes[0].axis("off")
        for i, slice in enumerate(slices):
            ylabel = dataframe_plot_presence.index.levels[0][i]
            create_presence_heatmap_axis(dataframe_plot_presence[slice], is_binary, linewidth, no_data_colour, xlabel,
                                         ylabel, axes[i + 1], 0.1, species_name=species_name,
                                         row_annotations=row_annotations, cbar_label=cbar_label,
                                         tick_font_size=xtick_fontsize)
        for axis in axes[1:-1]:
            axis.collections[0].colorbar.remove()
            axis.tick_params(top=False, bottom=False, left=False,
                             labelleft=False, labelbottom=False)
            axis.set_xlabel(None)
    else:
        fig, ax = plt.subplots()
        create_presence_heatmap_axis(dataframe_plot_presence, is_binary, linewidth, no_data_colour, xlabel, "Individual",
                                     ax, row_annotations=row_annotations, yticks=yticks, ytick_labels=ytick_labels,
                                     tick_font_size=xtick_fontsize)
        ax.set_title(subtitle)
    fig.set_dpi(500)
    fulltitle = f"{f'{name_prefix} ' if name_prefix else ''}{title}"
    fig.suptitle(fulltitle)
    if graphics_output_folder:
        os.makedirs(graphics_output_folder, exist_ok=True)
        plt.savefig(os.path.join(graphics_output_folder, f"{fulltitle} {species_name + ' 'if species_name else ''}presence heatmap {percent_identity_threshold:.2f}_perc_ident {coverage_threshold:.2f}_coverage.eps"), dpi=500)
        plt.savefig(os.path.join(graphics_output_folder, f"{fulltitle} {species_name + ' 'if species_name else ''}presence heatmap {percent_identity_threshold:.2f}_perc_ident {coverage_threshold:.2f}_coverage.png"), dpi=500)
    try:
        matplotlib.pyplot.show()
    except AttributeError:
        logging.getLogger().exception("Error showing image to user")
    pass


def create_presence_heatmap_axis(data, is_binary, linewidth, no_data_colour, xlabel, ylabel, ax=None,
                                 cbar_fraction=0.05, species_name=None, row_annotations=None, cbar_label=None,
                                 yticks=None, ytick_labels=None, tick_font_size=7):
    if data.columns[0] == "phenotype":
        phenotypes = data["phenotype"].unique()
        data = data.drop("phenotype", 1)
    if is_binary:
        h_linewidth = 0.4
        v_linewidth = 0.4
        linewidth = 0
        cbar_label = cbar_label
        annot_ticklabels = list(map(lambda x: row_annotations[x] if x in row_annotations.keys() else '', data.index.get_level_values("patient_id").values))
        seaborn.heatmap(data=data, vmin=-1, vmax=1,
                        cmap=[no_data_colour, "#d95f02", "#1b9e77"],
                        # cmap="Dark2",
                        mask=data.isnull(),
                        cbar_kws={"location": "bottom", "fraction": cbar_fraction, "label": cbar_label,
                                  "format": "%1.0f", "ticks": [-0.67, 0, 0.67],
                                  "drawedges": "True" },
                        yticklabels=annot_ticklabels,
                        # xticklabels=count_name + column_names,
                        linewidth=linewidth,
                        linecolor='black', ax=ax)
        ax.tick_params(axis='y', labelrotation=0, labelright=True, labelleft=False)
        ax.yaxis.tick_right()
        non_empty_indices = [[i, label] for i, label in enumerate(annot_ticklabels) if label != '']
        tick_idxs = [i + 0.5 for i, label in non_empty_indices]
        y_tick_labels = [label for i, label in non_empty_indices]
        ax.set_yticks(tick_idxs, y_tick_labels)
        # for item in groups:
        #     ax2.axhline(item, color='blue', lw=1)
        for i in range(data.shape[0] + 1):
            ax.axhline(i, color='black', lw=h_linewidth)
        for i in [0, data.shape[1] - 1, data.shape[1]]:
            ax.axvline(i, color='black', lw=v_linewidth)
        ax.collections[0].colorbar.set_ticklabels(["No sample", "Not present", "Present"])
        ax.collections[0].colorbar.outline.set(visible=True, lw=1, edgecolor="black")
        ax.collections[0].colorbar.ax.invert_xaxis()
    else:
        seaborn.heatmap(data=data, vmin=0, vmax=4,
                        cmap=["black", "blue", "green", "yellow", "red"],
                        cbar_kws={"location": "bottom", "fraction": cbar_fraction, "label": "XyGUL count",
                                        "ticks": [0.42, 1.2, 2, 2.8, 3.6],
                                        "format": "%1.0f", "drawedges": "True"},
                        yticklabels=[],
                        # xticklabels=count_name + column_names,
                        linewidth=linewidth, linecolor='black',
                        ax=ax
                        )
        cbar_widths = 1
        ax.collections[0].colorbar.outline.set(visible=True, lw=cbar_widths*.75, edgecolor="black")
        # ax.collections[0].colorbar.dividers.set(visible=True, linewidth=0.5, color="black")

        ax.collections[0].colorbar.dividers.set_color('black')
        ax.collections[0].colorbar.dividers.set_linewidth(cbar_widths)

    if yticks is not None and ytick_labels is not None:
        ax.set_yticks(yticks, labels=ytick_labels, fontsize=tick_font_size)

    ax.set_facecolor(no_data_colour)
    pul_name_angle, x_labels = adjusted_labels(ax.get_xticklabels())
    ax.set_xticklabels(x_labels, rotation=pul_name_angle, fontsize=tick_font_size)
    ax.set(ylabel=ylabel)

    # Move xticklabels up slightly
    if pul_name_angle == 45:
        ax.tick_params(axis='x', which='both', pad=-4)

    if xlabel:
        ax.set(xlabel=xlabel)
    return ax
