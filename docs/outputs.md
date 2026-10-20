# Outputs

For `-o results/taxonomy.tsv`, the default run writes:

| File | Contents |
| --- | --- |
| `results/taxonomy.ultrastrict.tsv` | Cumulative taxonomy counts at 0% edits and 120 aligned bases. |
| `results/taxonomy.strict.tsv` | Cumulative taxonomy counts at 2% edits and 100 aligned bases. |
| `results/taxonomy.standard.tsv` | Cumulative taxonomy counts at 3% edits and 80 aligned bases. |
| `results/taxonomy.permissive.tsv` | Cumulative taxonomy counts at 5% edits and 60 aligned bases. |
| `results/taxonomy.html` | Offline, self-contained report comparing levels and showing bar and sunburst plots for each. |
| `results/taxonomy.manifest.json` | Version, parameters, level settings and output paths, tool versions, database name, checksum and source releases, inputs, and fastp/assignment counts. |

When only one level is selected, its TSV is written exactly to `-o` (here, `results/taxonomy.tsv`). With several levels, each TSV uses `<stem>.<level><suffix>`; the HTML report and manifest still use the `-o` stem. Each TSV has a `Taxonomy` column and one integer count column per sample, including zero-hit samples. Rows are sorted by `Taxonomy`, with `Unclassified` first.

The working directory additionally holds fastp reports, per-sample TSVs under `extract/<level>/`, `taxonomy/<level>/` and `aggregate/<level>/`, and extraction stats in `extract/<sample>.stats.json`. Unless `--no-bam` is used, it also holds one BAM per sample in `bowtie/`.

## Reading the count table

Counts are propagated to every resolved ancestor. A fragment assigned to species therefore also appears in its genus, family, and higher-rank rows. **Do not sum all taxonomy rows** to obtain a total.

- A fragment with no shared named taxon is `Unclassified`.
- A sample with no accepted alignments has a zero column.
- Both mates of a pair share a read ID and receive one lowest-common-ancestor assignment.

## HTML report

The report compares the levels on an Overview tab and gives each level a composition view, a taxonomy detection grid and a per-sample sunburst, with CSV export for every table. See [](report.md) for a tour of its sections and controls.

## Manifest

The JSON manifest records everything needed to reproduce a run: BarcodeMapper version, parameters and ordered levels with their TSV paths, tool versions, the database name, build date, source releases and checksum, the input files, and fastp and per-level assignment counts. `samples.<sample>.assigned_read_ids` is a per-level object `{level: count}`. Keep it with the count tables.
