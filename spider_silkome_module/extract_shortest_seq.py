import re
from pathlib import Path

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from loguru import logger
from tqdm import tqdm
import typer

from spider_silkome_module.config import EXTERNAL_DATA_DIR, INTERIM_DATA_DIR

app = typer.Typer()


def parse_clstr_file(clstr_path: Path) -> dict[int, list[tuple[int, str]]]:
    """
    Parse CD-HIT cluster file and return a dict of cluster_id -> [(length, seq_id), ...].
    """
    clusters: dict[int, list[tuple[int, str]]] = {}
    current_cluster = -1

    with open(clstr_path) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">Cluster"):
                current_cluster = int(line.split()[1])
                clusters[current_cluster] = []
            else:
                # Parse line like: 0	1854aa, >5826_MaSp2_CTD... *
                match = re.match(r"\d+\s+(\d+)aa,\s+>([^.]+)\.\.\.", line)
                if match:
                    length = int(match.group(1))
                    seq_id = match.group(2)
                    clusters[current_cluster].append((length, seq_id))

    return clusters


def get_shortest_seq_ids(clusters: dict[int, list[tuple[int, str]]]) -> set[str]:
    """
    Get the shortest sequence ID from each cluster.
    """
    shortest_ids = set()
    for cluster_id, members in clusters.items():
        if members:
            shortest = min(members, key=lambda x: x[0])
            shortest_ids.add(shortest[1])
            logger.debug(f"Cluster {cluster_id}: shortest = {shortest[1]} ({shortest[0]}aa)")
    return shortest_ids


@app.command()
def main(
    clstr_path: Path = INTERIM_DATA_DIR / "miniprot_mapping_20260127" / "cdhit_rep_seq.fa.clstr",
    fasta_path: Path = EXTERNAL_DATA_DIR / "spider-silkome-database.v1.prot.fixed.renamed.fasta",
    output_path: Path = INTERIM_DATA_DIR / "miniprot_mapping_20260127" / "shortest_seqs.fasta",
):
    """
    Extract the shortest sequence from each cluster in CD-HIT output.
    """
    logger.info(f"Parsing cluster file: {clstr_path}")
    clusters = parse_clstr_file(clstr_path)
    logger.info(f"Found {len(clusters)} clusters")

    shortest_ids = get_shortest_seq_ids(clusters)
    logger.info(f"Extracting {len(shortest_ids)} shortest sequences")

    logger.info(f"Reading FASTA file: {fasta_path}")
    seq_dict = SeqIO.to_dict(SeqIO.parse(fasta_path, "fasta"))

    records: list[SeqRecord] = []
    missing_ids = []
    for seq_id in tqdm(shortest_ids, desc="Extracting sequences"):
        if seq_id in seq_dict:
            record = seq_dict[seq_id]
            seq_len = len(record.seq)
            record.id = f"{seq_id}_{seq_len}aa"
            record.description = ""
            records.append(record)
        else:
            missing_ids.append(seq_id)

    if missing_ids:
        logger.warning(f"Missing {len(missing_ids)} sequences in FASTA file")
        for mid in missing_ids[:10]:
            logger.warning(f"  Missing: {mid}")

    records.sort(key=lambda r: int(r.id.split("_")[0]))

    logger.info(f"Writing {len(records)} sequences to {output_path}")
    SeqIO.write(records, output_path, "fasta")

    logger.success(f"Extracted {len(records)} shortest sequences from {len(clusters)} clusters")


if __name__ == "__main__":
    app()
