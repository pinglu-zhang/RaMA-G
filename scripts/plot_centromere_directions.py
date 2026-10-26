"""Plot complete experimental direction evidence; never runs alignment.

Figure contract: show how raw matches, deduplicated collinear length, and
contextual event hypotheses differ. Quantitative three-panel grid, Python only,
183 x 220 mm; editable SVG/PDF plus PNG. Every input observation is retained;
repeated provenance rows for the SAME anchor ID are displayed once. Region
scores are not summed across overlapping windows. No biological truth inferred.
"""
import argparse
import csv
import json
from pathlib import Path


def rows(path):
    with Path(path).open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("anchors", "regions", "tasks", "candidates", "output"):
        parser.add_argument("--" + name, required=True, type=Path)
    parser.add_argument("--title", default="Experimental centromere direction evidence")
    args = parser.parse_args()
    # Fail explicitly when the selected Python plotting runtime is unavailable.
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.collections import LineCollection, PatchCollection
    from matplotlib.colors import LinearSegmentedColormap, Normalize
    from matplotlib.lines import Line2D
    from matplotlib.patches import Rectangle
    from matplotlib.ticker import FuncFormatter

    matplotlib.rcParams.update({
        "font.family": "sans-serif", "font.sans-serif": ["DejaVu Sans"],
        "font.size": 7, "axes.titlesize": 8, "axes.labelsize": 7,
        "xtick.labelsize": 6, "ytick.labelsize": 6, "legend.fontsize": 6,
        "pdf.fonttype": 42, "svg.fonttype": "none",
        "axes.spines.top": False, "axes.spines.right": False,
        "legend.frameon": False,
    })
    evidence, regions, tasks, candidates = (rows(p) for p in
        (args.anchors, args.regions, args.tasks, args.candidates))
    root = next(t for t in tasks if int(t["parent_id"]) == 0)
    rlen, qlen = int(root["reference_end"]), int(root["query_end"])
    unique, initial = {}, set()
    keys = ("reference_begin", "reference_end", "query_begin", "query_end", "strand")
    for row in evidence:
        key = int(row["anchor_id"])
        if key in unique and any(unique[key][k] != row[k] for k in keys):
            raise ValueError("Anchor ID refers to inconsistent coordinate mappings")
        unique[key] = row
        if int(row["depth"]) == 0:
            initial.add(key)
    blue, orange = "#14558A", "#B64B12"
    colors = {"+": blue, "-": orange}
    fig, axes = plt.subplots(3, 1, figsize=(183 / 25.4, 220 / 25.4),
                             layout="constrained")
    fig.suptitle(args.title, fontsize=9)

    def dotplot(ax, selected, neutral=False):
        for strand in ("+", "-"):
            segments = []
            for a in selected:
                if a["strand"] != strand:
                    continue
                rb, re, qb, qe = (int(a[k]) for k in keys[:4])
                if not (0 <= rb < re <= rlen and 0 <= qb < qe <= qlen and re - rb == qe - qb):
                    raise ValueError("Invalid anchor geometry")
                segments.append(((rb, qb if strand == "+" else qe),
                                 (re, qe if strand == "+" else qb)))
            ax.add_collection(LineCollection(segments, linewidths=0.55,
                colors="#A0A0A0" if neutral else colors[strand], alpha=0.35 if neutral else 0.7))
        ax.set(xlim=(0, rlen), ylim=(0, qlen), xlabel="Reference coordinate (bp)",
               ylabel="Original query coordinate (bp)")

    dotplot(axes[0], [unique[i] for i in sorted(initial)])
    axes[0].set_title(f"a  Initial bidirectional matches ({len(initial):,} distinct IDs)", loc="left", pad=24)
    axes[0].legend(handles=[Line2D([0], [0], color=blue, label="Forward match"),
                           Line2D([0], [0], color=orange, label="Reverse match")],
                   loc="lower left", bbox_to_anchor=(0, 1.0), ncols=2)

    f = [int(r["forward_support"]) for r in regions]
    r = [int(r["reverse_support"]) for r in regions]
    d = [float(row["advantage"]) for row in regions]
    cmap = LinearSegmentedColormap.from_list("direction", [orange, "#D8D8D8", blue])
    points = axes[1].scatter(f, r, c=d, cmap=cmap, norm=Normalize(-1, 1),
                             s=8, linewidths=0, alpha=0.65)
    largest = max([1] + f + r)
    axes[1].plot([0, largest], [0, largest], color="#686868", lw=0.7, ls="--")
    axes[1].set_xscale("symlog", linthresh=1)
    axes[1].set_yscale("symlog", linthresh=1)
    # Plain numeric labels avoid undersized logarithmic superscripts.
    axes[1].xaxis.set_major_formatter(FuncFormatter(lambda value, _: f"{value:g}"))
    axes[1].yaxis.set_major_formatter(FuncFormatter(lambda value, _: f"{value:g}"))
    axes[1].set(xlabel="Forward effective chain length (bp; symlog)",
                ylabel="Reverse effective chain length (bp; symlog)")
    axes[1].set_title(f"b  All local windows ({len(regions):,}; overlapping supports are not added)", loc="left")
    fig.colorbar(points, ax=axes[1], label="Direction advantage D (not probability)", pad=0.02)

    dotplot(axes[2], list(unique.values()), neutral=True)
    outlines, strong_outlines, strong_colors = [], [], []
    classes = {}
    for e in candidates:
        rb, re, qb, qe = (int(e[k]) for k in keys[:4])
        category = e["classification"]
        classes[category] = classes.get(category, 0) + 1
        rectangle = Rectangle((rb, qb), re - rb, qe - qb)
        strong = category == "bilateral_inversion_candidate"
        if strong:
            strong_outlines.append(rectangle)
            strong_colors.append(colors[e["strand"]])
        else:
            outlines.append(rectangle)
    axes[2].add_collection(PatchCollection(outlines, facecolors="none", edgecolors="#808080",
                                          linewidths=0.3, linestyles="dotted", alpha=0.07))
    axes[2].add_collection(PatchCollection(strong_outlines, facecolors="none", edgecolors=strong_colors,
                                          linewidths=0.8, linestyles="solid", alpha=0.9))
    axes[2].set_title("c  Recursive matches and all hypothesis bounds (not alignment spans)", loc="left", pad=24)
    axes[2].legend(handles=[Line2D([0], [0], color=orange, label="Two-flank candidate: solid outline"),
                           Line2D([0], [0], color="#808080", ls=":", label="Other / ambiguous: dotted outline")],
                   loc="lower left", bbox_to_anchor=(0, 1.0), ncols=2)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    targets = [args.output.with_suffix(ext) for ext in (".svg", ".pdf", ".png", ".qa.json")]
    if any(p.exists() for p in targets):
        raise FileExistsError("Refusing to replace an existing figure bundle")
    for path in targets[:3]:
        fig.savefig(path, dpi=600)
    plt.close(fig)
    targets[3].write_text(json.dumps({
        "input_provenance_rows": len(evidence), "distinct_anchor_ids": len(unique),
        "initial_anchor_ids": len(initial), "regions": len(regions),
        "hypotheses": len(candidates), "classification_counts": classes,
        "excluded_observations": 0,
        "deduplication": "Same anchor ID repeated for search provenance is drawn once",
        "support": "Length per same 2D region; no addition across windows",
        "uncertainty": "Hypothesis rectangles are NOT breakpoint confidence intervals; TSV contains bounds",
        "statistics": "Deterministic single sequence-pair evidence; no replicate CI or p-values",
        "raster_or_vector": "Vector match segments and bounds; PNG is a 600 dpi preview",
        "visual_inspection": "Required after rendering; not established by this script",
        "hypothesis_opacity": "All unresolved/ambiguous bounds at alpha 0.07; two-flank candidates at 0.9; neither opacity is probability",
    }, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
