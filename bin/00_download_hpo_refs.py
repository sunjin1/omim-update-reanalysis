#!/usr/bin/env python3
"""
Download reference files required by the OMIM update-driven reanalysis pipeline.

Default outputs:
  ref/hp.obo
  ref/genes_to_phenotype.txt

These files are intentionally not bundled with the repository so users can
retrieve the current HPO release at run time.
"""

import argparse
import os
import sys
import urllib.request
from datetime import datetime

HP_OBO_URL = "http://purl.obolibrary.org/obo/hp.obo"
GENES_TO_PHENOTYPE_URL = "https://purl.obolibrary.org/obo/hp/hpoa/genes_to_phenotype.txt"
PHENOTYPE_TO_GENES_URL = "https://purl.obolibrary.org/obo/hp/hpoa/phenotype_to_genes.txt"


def download(url: str, out_path: str, force: bool = False, timeout: int = 180) -> bool:
    os.makedirs(os.path.dirname(os.path.abspath(out_path)), exist_ok=True)

    if os.path.exists(out_path) and os.path.getsize(out_path) > 0 and not force:
        print(f"[SKIP] {out_path} already exists. Use --force to overwrite.")
        return False

    tmp_path = out_path + ".tmp"
    print(f"[DOWNLOAD] {url}")
    print(f"           -> {out_path}")

    req = urllib.request.Request(
        url,
        headers={"User-Agent": "OMIM-reanalysis-pipeline/1.0"},
    )
    with urllib.request.urlopen(req, timeout=timeout) as response, open(tmp_path, "wb") as fout:
        fout.write(response.read())

    if os.path.getsize(tmp_path) == 0:
        raise RuntimeError(f"Downloaded empty file: {url}")

    os.replace(tmp_path, out_path)
    return True


def write_manifest(path: str, records):
    manifest = os.path.join(path, "download_manifest.tsv")
    with open(manifest, "w", encoding="utf-8") as f:
        f.write("file\turl\tdownload_time\tsize_bytes\n")
        now = datetime.now().isoformat(timespec="seconds")
        for out_path, url in records:
            size = os.path.getsize(out_path) if os.path.exists(out_path) else 0
            f.write(f"{out_path}\t{url}\t{now}\t{size}\n")
    print(f"[DONE] manifest written to {manifest}")


def parse_args():
    p = argparse.ArgumentParser(description="Download HPO reference files for OMIM reanalysis.")
    p.add_argument("--out-dir", default="ref", help="Reference output directory. Default: ref")
    p.add_argument("--hp-obo", default=None, help="Output path for hp.obo. Default: <out-dir>/hp.obo")
    p.add_argument("--genes-to-phenotype", default=None, help="Output path for genes_to_phenotype.txt. Default: <out-dir>/genes_to_phenotype.txt")
    p.add_argument("--phenotype-to-genes", default=None, help="Optional output path for phenotype_to_genes.txt. If provided, this file is downloaded too.")
    p.add_argument("--force", action="store_true", help="Overwrite existing files.")
    return p.parse_args()


def main():
    args = parse_args()
    out_dir = args.out_dir
    os.makedirs(out_dir, exist_ok=True)

    hp_obo = args.hp_obo or os.path.join(out_dir, "hp.obo")
    g2p = args.genes_to_phenotype or os.path.join(out_dir, "genes_to_phenotype.txt")

    records = []
    download(HP_OBO_URL, hp_obo, force=args.force)
    records.append((hp_obo, HP_OBO_URL))

    download(GENES_TO_PHENOTYPE_URL, g2p, force=args.force)
    records.append((g2p, GENES_TO_PHENOTYPE_URL))

    if args.phenotype_to_genes:
        download(PHENOTYPE_TO_GENES_URL, args.phenotype_to_genes, force=args.force)
        records.append((args.phenotype_to_genes, PHENOTYPE_TO_GENES_URL))

    write_manifest(out_dir, records)
    print("[DONE] HPO reference download completed.")


if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        print(f"[ERROR] {e}", file=sys.stderr)
        sys.exit(1)
