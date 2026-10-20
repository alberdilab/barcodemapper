# The HTML report

Every run writes one self-contained HTML report (`<stem>.html`). It works offline in any modern browser: the charts, styles and logo are embedded, so you can email the file or archive it with the results.

The header shows the BarcodeMapper version, the number of samples and levels, the time the report was made, and a theme button (◐) that switches between the light and dark themes.

![Example report](../images/report.png)

## Layout

With several stringency levels the report has one tab per view, under the header:

- **Overview** compares all levels.
- One **level tab** (ultrastrict, strict, standard, permissive and, if used, custom) holds the detailed results for that level.

With a single level there are no tabs and the report opens directly on that level's detailed view.

## Overview tab

The Overview answers one question: how does stringency change what is retained and how far it is resolved? Cards at the top give the sample count and the headline numbers for each level, followed by five sections.

| Section | What it shows |
| --- | --- |
| 01 Summary | Per level: retained fragments, share of aligned fragments retained, classified share, and number of distinct taxa at each rank. |
| 02 Resolution by rank | Share of fragments resolved down to each rank, per level. |
| 03 Kingdom composition | Kingdom composition of the classified fragments, per level. |
| 04 Samples | Classified fragments for every sample at every level. Use *Search* to filter and click a column header to sort. |
| 05 Distribution | A dot plot comparing levels sample by sample. Sort by sample or level, and switch between fragment counts and % classified. |

Typical use: read the Summary to see how much the stricter levels discard, check *Resolution by rank* to see whether the loss buys deeper identification, then pick the level whose tab you will inspect.

## Level tabs

Each level tab has three sections.

### 01 Composition

Stacked bars of the taxonomic composition of every sample.

- **Sort** samples by name, total fragments or largest kingdom share.
- **Absolute / Share** switches between fragment counts and the fraction of marker-matching fragments.
- Click entries in the legend to hide or show categories.
- *Table view of composition data* opens the numbers behind the bars, with a CSV download.

### 02 Detection: taxonomy cladogram and detection grid

A taxon-by-sample grid of cumulative counts. Rows are the taxonomy tree; columns are samples. A taxon's count includes everything below it, so rows at different ranks overlap and should not be summed.

Navigating the tree:

- Click a row, or use **Right/Left** to expand or collapse a branch. **Up/Down** moves between rows, and **Enter** or **Space** also toggles a branch.
- **Depth** opens the tree down to a chosen rank; **Expand all** and **Collapse all** do the rest.
- **Find taxon** searches the whole tree and shows the matching branches. It does not change which branches you had expanded, and expanding or collapsing clears it.
- Drag the edge of the taxonomy column to resize it, or focus the handle and use Left/Right. The width is kept for the browser session.

Changing the display:

- **Bubbles / Heatmap** switches the cell style.
- **Fragments / Share of sample / Presence** switches between counts, the fraction of each sample's fragments, and detected or not.
- **Samples** orders the columns by name, total fragments, or the abundance of the taxon you have focused.

Taxon colours follow the phylum, and the same colours are used in the composition bars, the grid and the sunburst, so a colour means the same thing everywhere in the report.

Clicking a sample's column header selects it for the sunburst below; *View sunburst* jumps there.

### 03 Taxonomy

An interactive sunburst of one sample's taxonomy.

- Choose the sample in the **Samples** list beside it: filter by name, or step with the ↑ / ↓ buttons. The list appears only next to this section.
- Click a segment to zoom into that taxon; click the centre to go back up.
- The card shows the sample's total and classified fragments.
- Only taxa whose whole lineage is present are drawn.

## Downloading the data

Every table can be downloaded as CSV (**Download CSV**): the Overview tables, the folded composition table, and from the detection grid either **Download visible CSV** (the rows currently shown) or **Download all taxa CSV**. The per-sample table follows its search filter. Files carry a UTF-8 signature so spreadsheet software reads them correctly. Text cells starting with `=`, `+`, `-`, `@`, a double quote, tab, carriage return or line feed are prefixed with an apostrophe to stop spreadsheets from treating them as formulas; numeric cells stay numeric.

## Reading the numbers

- Counts are marker-matching fragments, not organism abundances (see [](limitations.md)).
- Relative bars use phylum-assigned plus unresolved and unclassified fragments as the denominator.
- Counts at a given rank are not necessarily subsets across levels: accepting different alignments can change the lowest common ancestor.
- Both mates of a pair count as one fragment.
