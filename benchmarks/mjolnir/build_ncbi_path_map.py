"""Build common taxonomy paths from an NCBI taxdump or Kraken2 ktaxonomy."""

import argparse
import csv
import hashlib
import json
from collections import Counter
from contextlib import ExitStack
from itertools import chain
from pathlib import Path
from tempfile import NamedTemporaryFile


RANKS = "kpcofgs"
NCBI_RANKS = {"kingdom": "k", "phylum": "p", "class": "c", "order": "o",
              "family": "f", "genus": "g", "species": "s"}
KRAKEN_RANKS = {rank.upper(): rank for rank in RANKS}
KRAKEN_RANKS["D"] = "domain"
DEFAULT_NAME_MAP = Path(__file__).with_name("metadata") / "05-default-databases" / "name-map.tsv"


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def dmp_rows(path, minimum):
    with path.open(encoding="utf-8") as handle:
        for number, line in enumerate(handle, 1):
            fields = [field.rstrip("\t| ") for field in line.rstrip("\r\n").split("\t|\t")]
            if len(fields) < minimum:
                raise ValueError(f"Invalid taxonomy row at {path}:{number}")
            yield fields


def name_map_rows(path):
    mapping = {}
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.reader(handle, delimiter="\t")
        if next(reader, None) != ["Rank", "NCBIName", "CommonName"]:
            raise ValueError(f"Invalid name-map header in {path}")
        for number, row in enumerate(reader, 2):
            if (len(row) != 3 or row[0] not in RANKS or not row[1] or not row[2]
                    or any(char in row[2] for char in ";\t\r\n")):
                raise ValueError(f"Invalid name-map row at {path}:{number}")
            key = row[0], row[1]
            if key in mapping and mapping[key] != row[2]:
                raise ValueError(f"Conflicting name map for {key!r}")
            mapping[key] = row[2]
    return mapping


def taxdump_data(directory, name_map):
    files = {name: directory / f"{name}.dmp" for name in
             ("nodes", "names", "merged", "delnodes")}
    nodes = {}
    input_rows = Counter({name: 0 for name in files})
    for fields in dmp_rows(files["nodes"], 3):
        input_rows["nodes"] += 1
        taxid, parent, rank = fields[:3]
        if taxid in nodes:
            raise ValueError(f"Duplicate taxid {taxid!r}")
        nodes[taxid] = parent, NCBI_RANKS.get(rank, rank if rank in ("superkingdom", "domain") else "")

    names = {}
    mapped = {}
    semicolon_names = 0
    for fields in dmp_rows(files["names"], 4):
        input_rows["names"] += 1
        taxid, name, _, name_class = fields[:4]
        rank = nodes.get(taxid, (None, ""))[1]
        if not rank:
            continue
        letter = "k" if rank in ("superkingdom", "domain") else rank
        replacement = name_map.get((letter, name))
        if replacement is not None:
            if taxid in mapped and mapped[taxid] != replacement:
                raise ValueError(f"Conflicting synonym mappings for taxid {taxid}")
            mapped[taxid] = replacement
        if name_class == "scientific name":
            if ";" in name:
                semicolon_names += 1
            names[taxid] = name.replace(";", ",").replace("\t", " ")

    merged = {}
    for fields in dmp_rows(files["merged"], 2):
        input_rows["merged"] += 1
        merged[fields[0]] = fields[1]
    deleted = set()
    for fields in dmp_rows(files["delnodes"], 1):
        input_rows["delnodes"] += 1
        deleted.add(fields[0])
    return nodes, names, mapped, merged, deleted, files, input_rows, semicolon_names


def kingdom_taxids_from_file(path):
    taxids = set()
    with path.open(encoding="utf-8") as handle:
        for number, line in enumerate(handle, 1):
            taxid = line.partition("#")[0].strip()
            if not taxid:
                continue
            if not taxid.isdecimal() or int(taxid) < 1:
                raise ValueError(f"Invalid kingdom taxid at {path}:{number}")
            taxids.add(taxid)
    return taxids, {"type": "taxid_file", "path": str(path), "sha256": sha256(path)}


def kingdom_taxids_from_taxdump(directory):
    nodes = directory / "nodes.dmp"
    taxids = {fields[0] for fields in dmp_rows(nodes, 3) if fields[2] == "kingdom"}
    return taxids, {"type": "taxdump", "directory": str(directory),
                    "nodes_dmp": str(nodes), "nodes_sha256": sha256(nodes)}


def ktaxonomy_data(path, name_map, kingdom_taxids):
    nodes = {}
    names = {}
    mapped = {}
    semicolon_names = 0
    input_rows = Counter({"ktaxonomy": 0})
    has_numeric_rank_code = False
    kingdom_taxids_found = 0
    for fields in dmp_rows(path, 5):
        input_rows["ktaxonomy"] += 1
        taxid, parent, code, _, name = fields[:5]
        if code[:1] in ("D", "K") and code[1:].isdigit():
            has_numeric_rank_code = True
        rank = KRAKEN_RANKS.get(code, "")
        if taxid in nodes:
            raise ValueError(f"Duplicate taxid {taxid!r}")
        if taxid in kingdom_taxids:
            kingdom_taxids_found += 1
            if rank != "domain":
                rank = "k"
        nodes[taxid] = parent, rank
        if rank:
            letter = "k" if rank == "domain" else rank
            replacement = name_map.get((letter, name))
            if replacement is not None:
                mapped[taxid] = replacement
            if ";" in name:
                semicolon_names += 1
            names[taxid] = name.replace(";", ",").replace("\t", " ")
    return (nodes, names, mapped, {}, set(), {"ktaxonomy": path}, input_rows,
            semicolon_names, has_numeric_rank_code, kingdom_taxids_found)


def build_paths(nodes, names, mapped):
    empty = ("",) * len(RANKS)
    cache = {}

    def use_kingdom_name(rank, domain_taxid):
        return rank in ("superkingdom", "domain") or rank == "k" and domain_taxid == "2759"

    def lineage(taxid):
        pending = []
        seen = set()
        current = taxid
        while current not in cache:
            if current not in nodes:
                raise ValueError(f"Missing parent taxid {current!r}")
            if current in seen:
                raise ValueError(f"Taxonomy cycle at taxid {current!r}")
            seen.add(current)
            pending.append(current)
            parent = nodes[current][0]
            if parent == current:
                break
            current = parent
        levels, domain_taxid = cache.get(current, (empty, None))
        for node in reversed(pending):
            rank = nodes[node][1]
            if rank:
                if node not in names:
                    raise ValueError(f"Missing scientific name for taxid {node}")
                if use_kingdom_name(rank, domain_taxid):
                    changed = list(levels)
                    changed[0] = mapped.get(node, names[node])
                    levels = tuple(changed)
                elif rank != "k" and rank not in ("superkingdom", "domain"):
                    changed = list(levels)
                    changed[RANKS.index(rank)] = mapped.get(node, names[node])
                    levels = tuple(changed)
                if rank in ("superkingdom", "domain"):
                    domain_taxid = node
            cache[node] = levels, domain_taxid
        return cache[taxid][0]

    return lineage


def build_map(taxdump, ktaxonomy, output, audit, name_map=DEFAULT_NAME_MAP, report=None,
              kingdom_taxids=None, kingdom_from_taxdump=None):
    outputs = [output, audit] + ([report] if report is not None else [])
    if len(set(outputs)) != len(outputs):
        raise ValueError("Output paths must be distinct")
    for path in outputs:
        if path.exists():
            raise ValueError(f"Output already exists: {path}")

    mapping = name_map_rows(name_map)
    if kingdom_taxids is not None and kingdom_from_taxdump is not None:
        raise ValueError("Choose only one kingdom taxid source")
    kingdom_source = None
    kingdom_taxids_found = 0
    if taxdump is not None:
        if kingdom_taxids is not None or kingdom_from_taxdump is not None:
            raise ValueError("Kingdom taxid options require --ktaxonomy")
        nodes, names, mapped, merged, deleted, files, input_rows, semicolon_names = taxdump_data(taxdump, mapping)
    else:
        selected = set()
        if kingdom_taxids is not None:
            selected, kingdom_source = kingdom_taxids_from_file(kingdom_taxids)
        elif kingdom_from_taxdump is not None:
            selected, kingdom_source = kingdom_taxids_from_taxdump(kingdom_from_taxdump)
        (nodes, names, mapped, merged, deleted, files, input_rows, semicolon_names,
         has_numeric_rank_code, kingdom_taxids_found) = ktaxonomy_data(ktaxonomy, mapping, selected)
        if has_numeric_rank_code and kingdom_source is None:
            raise ValueError("Numeric D/K rank codes in ktaxonomy require "
                             "--kingdom-taxids or --kingdom-from-taxdump")
        if kingdom_source is not None and not kingdom_taxids_found:
            raise ValueError("No supplied kingdom taxids found in ktaxonomy")
    lineage = build_paths(nodes, names, mapped)
    lost = Counter()
    resolved_depth_counts = Counter()
    missing_kingdom_count = 0
    rows = 0
    merged_rows = 0
    merged_missing_target = 0
    temporary = {}
    published = []
    try:
        with ExitStack() as stack:
            handle = stack.enter_context(NamedTemporaryFile(
                "w", dir=output.parent, prefix=f".{output.name}.", suffix=".tmp",
                newline="", encoding="utf-8", delete=False))
            temporary[output] = Path(handle.name)
            report_handle = None
            if report is not None:
                report_handle = stack.enter_context(NamedTemporaryFile(
                    "w", dir=report.parent, prefix=f".{report.name}.", suffix=".tmp",
                    encoding="utf-8", delete=False))
                temporary[report] = Path(report_handle.name)
                report_handle.write("TaxID\n")
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(("TaxID", "Taxonomy"))
            for taxid in chain((taxid for taxid in nodes if taxid not in merged), merged):
                if taxid in deleted:
                    continue
                target = taxid
                seen = set()
                while target in merged:
                    if target in seen:
                        raise ValueError(f"Merged taxid cycle at {taxid}")
                    seen.add(target)
                    target = merged[target]
                if target not in nodes or target in deleted:
                    merged_missing_target += 1
                    continue
                levels = lineage(target)
                for rank, level in zip(RANKS, levels):
                    if not level:
                        lost[rank] += 1
                depth = 0
                for level in levels:
                    if not level:
                        break
                    depth += 1
                resolved_depth_counts[depth] += 1
                if not levels[0]:
                    missing_kingdom_count += 1
                    if report_handle:
                        report_handle.write(f"{taxid}\n")
                    writer.writerow((taxid, ""))
                else:
                    writer.writerow((taxid, ";".join(f"{rank}__{level}" for rank, level in zip(RANKS, levels))))
                rows += 1
                merged_rows += taxid in merged

        audit_data = {
            "input_sha256": {name: sha256(path) for name, path in files.items()},
            "input_rows": dict(input_rows),
            "name_map_sha256": sha256(name_map),
            "name_map_rows": len(mapping),
            "kingdom_source": kingdom_source,
            "kingdom_taxids_found": kingdom_taxids_found,
            "node_count": len(nodes),
            "merged_count": len(merged),
            "merged_rows": merged_rows,
            "merged_missing_target": merged_missing_target,
            "deleted_count": len(deleted),
            "output_rows": rows,
            "missing_kingdom_count": missing_kingdom_count,
            "taxa_lost_by_rank": {rank: lost[rank] for rank in RANKS},
            "resolved_depth_counts": {str(depth): resolved_depth_counts[depth] for depth in range(8)},
            "name_map_changes": sum(mapped[taxid] != names[taxid] for taxid in mapped),
            "semicolon_names_replaced": semicolon_names,
        }
        with NamedTemporaryFile("w", dir=audit.parent, prefix=f".{audit.name}.",
                                suffix=".tmp", encoding="utf-8", delete=False) as handle:
            temporary[audit] = Path(handle.name)
            json.dump(audit_data, handle, indent=2)
            handle.write("\n")
        for path, staged in temporary.items():
            if path.exists():
                raise ValueError(f"Output already exists: {path}")
            staged.replace(path)
            published.append(path)
    except Exception:
        for staged in temporary.values():
            staged.unlink(missing_ok=True)
        for path in published:
            path.unlink(missing_ok=True)
        raise


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument("--taxdump", type=Path)
    source.add_argument("--ktaxonomy", type=Path)
    kingdoms = parser.add_mutually_exclusive_group()
    kingdoms.add_argument("--kingdom-taxids", type=Path)
    kingdoms.add_argument("--kingdom-from-taxdump", type=Path)
    parser.add_argument("--name-map", type=Path, default=DEFAULT_NAME_MAP)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--audit", type=Path)
    parser.add_argument("--report", type=Path)
    args = parser.parse_args(argv)
    audit = args.audit or args.output.with_suffix(args.output.suffix + ".audit.json")
    build_map(args.taxdump, args.ktaxonomy, args.output, audit, args.name_map, args.report,
              args.kingdom_taxids, args.kingdom_from_taxdump)


if __name__ == "__main__":
    main()
