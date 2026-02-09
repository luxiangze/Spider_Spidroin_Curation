import json
import re
from pathlib import Path

from loguru import logger
import typer

from spider_silkome_module.config import INTERIM_DATA_DIR

app = typer.Typer()

BASE_URL = "http://192.168.101.184:9501/spider-genome"


def parse_tree_file(tree_path: Path) -> dict:
    """
    Parse the tree file and extract file structure.
    Returns a dict with species -> files mapping.
    """
    structure = {
        "ref_gff": {},
        "bgi_rna_bam": {},
        "bgi_rna_bigwig": {},
        "ont_rna_bam": {},
        "ont_rna_bigwig": {},
        "hmmer_search": {},
        "miniprot": {},
    }

    current_section = None
    current_species = None
    current_species_dir = None

    with open(tree_path) as f:
        for line in f:
            line = line.rstrip()
            if not line:
                continue

            # Detect section
            if "01.ref_gff" in line:
                current_section = "ref_gff"
            elif "BGI_RNA_align" in line:
                current_section = "bgi_rna"
            elif "ONT_RNA_align" in line:
                current_section = "ont_rna"
            elif "hmmer_search_output" in line:
                current_section = "hmmer_search"
            elif "miniprot_output" in line:
                current_section = "miniprot"
            elif "├─ bam" in line or "│  ├─ bam" in line:
                if current_section == "bgi_rna":
                    current_section = "bgi_rna_bam"
                elif current_section == "ont_rna":
                    current_section = "ont_rna_bam"
            elif "├─ bigwig" in line or "│  └─ bigwig" in line or "│  ├─ bigwig" in line:
                if current_section in ["bgi_rna", "bgi_rna_bam"]:
                    current_section = "bgi_rna_bigwig"
                elif current_section in ["ont_rna", "ont_rna_bam"]:
                    current_section = "ont_rna_bigwig"

            # Parse species directory in ref_gff section
            if current_section == "ref_gff":
                match = re.search(r"├─ (\d+)\.([A-Za-z_\.0-9]+)", line)
                if not match:
                    match = re.search(r"└─ (\d+)\.([A-Za-z_\.0-9]+)", line)
                if match:
                    species_num = match.group(1)
                    species_name = match.group(2)
                    current_species = species_name
                    current_species_dir = f"{species_num}.{species_name}"
                    structure["ref_gff"][species_name] = {
                        "dir": current_species_dir,
                        "files": [],
                    }
                elif current_species and ("├─" in line or "└─" in line):
                    file_match = re.search(r"[├└]─ ([^\s]+)", line)
                    if file_match:
                        filename = file_match.group(1)
                        structure["ref_gff"][current_species]["files"].append(filename)

            # Parse BAM files
            elif current_section == "bgi_rna_bam":
                match = re.search(r"[├└]─ ([A-Za-z_]+)\.(bam|bam\.csi)", line)
                if match:
                    species = match.group(1)
                    ext = match.group(2)
                    if species not in structure["bgi_rna_bam"]:
                        structure["bgi_rna_bam"][species] = []
                    structure["bgi_rna_bam"][species].append(f"{species}.{ext}")

            elif current_section == "ont_rna_bam":
                match = re.search(r"[├└]─ ([A-Za-z_]+)\.(sorted\.bam|sorted\.bam\.csi)", line)
                if match:
                    species = match.group(1)
                    ext = match.group(2)
                    if species not in structure["ont_rna_bam"]:
                        structure["ont_rna_bam"][species] = []
                    structure["ont_rna_bam"][species].append(f"{species}.{ext}")

            # Parse BigWig files
            elif current_section == "bgi_rna_bigwig":
                match = re.search(r"[├└]─ ([A-Za-z_]+)\.bw", line)
                if match:
                    species = match.group(1)
                    structure["bgi_rna_bigwig"][species] = f"{species}.bw"

            elif current_section == "ont_rna_bigwig":
                match = re.search(r"[├└]─ ([A-Za-z_]+)\.bw", line)
                if match:
                    species = match.group(1)
                    structure["ont_rna_bigwig"][species] = f"{species}.bw"

            # Parse hmmer_search GFF
            elif current_section == "hmmer_search":
                match = re.search(r"[├└]─ ([A-Za-z_]+)$", line)
                if match:
                    current_species = match.group(1)
                gff_match = re.search(r"[├└]─ ([A-Za-z_]+)\.gff", line)
                if gff_match and current_species:
                    structure["hmmer_search"][current_species] = f"{current_species}/{current_species}.gff"

            # Parse miniprot GFF
            elif current_section == "miniprot":
                match = re.search(r"[├└]─ ([A-Za-z_]+)$", line)
                if match:
                    current_species = match.group(1)
                gff_match = re.search(r"[├└]─ ([A-Za-z_]+)\.gff", line)
                if gff_match and current_species:
                    structure["miniprot"][current_species] = f"{current_species}/{current_species}.gff"

    return structure


def create_assembly(species_name: str, species_info: dict) -> dict:
    """
    Create assembly config for a species.
    """
    species_dir = species_info["dir"]
    files = species_info["files"]

    # Find genome file
    genome_file = None
    for f in files:
        if f.endswith(".genome.fa.gz") or (f.endswith(".fa.gz") and "genome" not in f):
            genome_file = f
            break
        if f.endswith(".fna.gz"):
            genome_file = f
            break

    if not genome_file:
        for f in files:
            if f.endswith(".fa.gz"):
                genome_file = f
                break

    if not genome_file:
        return None

    assembly = {
        "name": species_name,
        "sequence": {
            "type": "ReferenceSequenceTrack",
            "trackId": f"{species_name}-ReferenceSequenceTrack",
            "adapter": {
                "type": "BgzipFastaAdapter",
                "fastaLocation": {
                    "uri": f"{BASE_URL}/01.ref_gff/{species_dir}/{genome_file}",
                    "locationType": "UriLocation",
                },
                "faiLocation": {
                    "uri": f"{BASE_URL}/01.ref_gff/{species_dir}/{genome_file}.fai",
                    "locationType": "UriLocation",
                },
                "gziLocation": {
                    "uri": f"{BASE_URL}/01.ref_gff/{species_dir}/{genome_file}.gzi",
                    "locationType": "UriLocation",
                },
            },
        },
    }

    return assembly


def create_gff_track(species_name: str, species_info: dict) -> dict | None:
    """
    Create GFF annotation track for a species.
    """
    species_dir = species_info["dir"]
    files = species_info["files"]

    # Find GFF file (prefer .gff.gz with .tbi)
    gff_gz_file = None
    for f in files:
        if f.endswith(".gff.gz") and not f.endswith(".tbi"):
            gff_gz_file = f
        elif f.endswith(".gff3.gz") and not f.endswith(".tbi"):
            gff_gz_file = f

    if gff_gz_file:
        tbi_file = f"{gff_gz_file}.tbi"
        if tbi_file in files:
            return {
                "type": "FeatureTrack",
                "trackId": f"{species_name}-gene-annotation",
                "name": "Gene Annotation",
                "assemblyNames": [species_name],
                "category": ["Annotation"],
                "adapter": {
                    "type": "Gff3TabixAdapter",
                    "gffGzLocation": {
                        "uri": f"{BASE_URL}/01.ref_gff/{species_dir}/{gff_gz_file}",
                        "locationType": "UriLocation",
                    },
                    "index": {
                        "location": {
                            "uri": f"{BASE_URL}/01.ref_gff/{species_dir}/{tbi_file}",
                            "locationType": "UriLocation",
                        },
                    },
                },
            }

    return None


def create_bam_track(
    species_name: str,
    bam_file: str,
    track_name: str,
    track_id_prefix: str,
    base_path: str,
    category: list[str],
) -> dict:
    """
    Create BAM alignment track.
    """
    return {
        "type": "AlignmentsTrack",
        "trackId": f"{species_name}-{track_id_prefix}",
        "name": track_name,
        "assemblyNames": [species_name],
        "category": category,
        "adapter": {
            "type": "BamAdapter",
            "bamLocation": {
                "uri": f"{BASE_URL}/{base_path}/{bam_file}",
                "locationType": "UriLocation",
            },
            "index": {
                "location": {
                    "uri": f"{BASE_URL}/{base_path}/{bam_file}.csi",
                    "locationType": "UriLocation",
                },
                "indexType": "CSI",
            },
        },
    }


def create_bigwig_track(
    species_name: str,
    bw_file: str,
    track_name: str,
    track_id_prefix: str,
    base_path: str,
    category: list[str],
) -> dict:
    """
    Create BigWig coverage track.
    """
    return {
        "type": "QuantitativeTrack",
        "trackId": f"{species_name}-{track_id_prefix}",
        "name": track_name,
        "assemblyNames": [species_name],
        "category": category,
        "adapter": {
            "type": "BigWigAdapter",
            "bigWigLocation": {
                "uri": f"{BASE_URL}/{base_path}/{bw_file}",
                "locationType": "UriLocation",
            },
        },
    }


def create_gff_track_simple(
    species_name: str,
    gff_path: str,
    track_name: str,
    track_id_prefix: str,
    base_path: str,
    category: list[str],
) -> dict:
    """
    Create simple GFF track (non-indexed).
    """
    return {
        "type": "FeatureTrack",
        "trackId": f"{species_name}-{track_id_prefix}",
        "name": track_name,
        "assemblyNames": [species_name],
        "category": category,
        "adapter": {
            "type": "Gff3Adapter",
            "gffLocation": {
                "uri": f"{BASE_URL}/{base_path}/{gff_path}",
                "locationType": "UriLocation",
            },
        },
    }


def generate_config(structure: dict) -> dict:
    """
    Generate complete JBrowse2 config from parsed structure.
    """
    assemblies = []
    tracks = []

    # Create assemblies and annotation tracks
    for species_name, species_info in structure["ref_gff"].items():
        assembly = create_assembly(species_name, species_info)
        if assembly:
            assemblies.append(assembly)

            # Add gene annotation track if available
            gff_track = create_gff_track(species_name, species_info)
            if gff_track:
                tracks.append(gff_track)

    # Create RNA-seq tracks
    for species_name in structure["bgi_rna_bam"]:
        tracks.append(
            create_bam_track(
                species_name,
                f"{species_name}.bam",
                "BGI RNA-seq Alignment",
                "bgi-rna-bam",
                "BGI_RNA_align/bam",
                ["RNA-seq", "BGI"],
            )
        )

    for species_name, bw_file in structure["bgi_rna_bigwig"].items():
        tracks.append(
            create_bigwig_track(
                species_name,
                bw_file,
                "BGI RNA-seq Coverage",
                "bgi-rna-coverage",
                "BGI_RNA_align/bigwig",
                ["RNA-seq", "BGI"],
            )
        )

    for species_name in structure["ont_rna_bam"]:
        tracks.append(
            create_bam_track(
                species_name,
                f"{species_name}.sorted.bam",
                "ONT RNA-seq Alignment",
                "ont-rna-bam",
                "ONT_RNA_align/bam",
                ["RNA-seq", "ONT"],
            )
        )

    for species_name, bw_file in structure["ont_rna_bigwig"].items():
        tracks.append(
            create_bigwig_track(
                species_name,
                bw_file,
                "ONT RNA-seq Coverage",
                "ont-rna-coverage",
                "ONT_RNA_align/bigwig",
                ["RNA-seq", "ONT"],
            )
        )

    # Create HMMER search tracks
    for species_name, gff_path in structure["hmmer_search"].items():
        tracks.append(
            create_gff_track_simple(
                species_name,
                gff_path,
                "HMMER Search Results",
                "hmmer-search",
                "hmmer_search_output",
                ["Spidroin Annotation", "HMMER"],
            )
        )

    # Create miniprot tracks
    for species_name, gff_path in structure["miniprot"].items():
        tracks.append(
            create_gff_track_simple(
                species_name,
                gff_path,
                "Miniprot Mapping",
                "miniprot",
                "miniprot_output",
                ["Spidroin Annotation", "Miniprot"],
            )
        )

    return {"assemblies": assemblies, "tracks": tracks}


@app.command()
def main(
    tree_file: Path = INTERIM_DATA_DIR / "rustfs_spider-genome_tree.txt",
    output_file: Path = INTERIM_DATA_DIR / "jbrowse2_config.json",
    base_url: str = BASE_URL,
):
    """
    Generate JBrowse2 configuration from rustfs tree file.
    """
    global BASE_URL
    BASE_URL = base_url

    logger.info(f"Parsing tree file: {tree_file}")
    structure = parse_tree_file(tree_file)

    logger.info(f"Found {len(structure['ref_gff'])} species in ref_gff")
    logger.info(f"Found {len(structure['bgi_rna_bam'])} species with BGI RNA BAM")
    logger.info(f"Found {len(structure['ont_rna_bam'])} species with ONT RNA BAM")
    logger.info(f"Found {len(structure['hmmer_search'])} species with HMMER results")
    logger.info(f"Found {len(structure['miniprot'])} species with miniprot results")

    logger.info("Generating JBrowse2 config...")
    config = generate_config(structure)

    logger.info(f"Generated {len(config['assemblies'])} assemblies and {len(config['tracks'])} tracks")

    logger.info(f"Writing config to: {output_file}")
    with open(output_file, "w") as f:
        json.dump(config, f, indent=2)

    logger.success("JBrowse2 config generated successfully")


if __name__ == "__main__":
    app()
