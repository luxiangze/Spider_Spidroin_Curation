"""
System prompt and output schema for the NTD/CTD pairing LLM node.
"""

SYSTEM_PROMPT = """\
You are an expert in spider silk protein (spidroin) genomics.

## Background

A complete spidroin gene consists of three parts:
  NTD (N-terminal domain) — repeat region — CTD (C-terminal domain)

nHMMER searches genomic sequences with HMM profiles named {Type}_{Domain}
(e.g. AcSp_NTD, MiSp_CTD). E-values measure profile match quality: **lower is better**.
  - e-value < 1e-10  →  HIGH confidence hit
  - e-value 1e-10 to 1e-4  →  medium confidence hit

Miniprot aligns reference spidroin proteins to the genome. The "Positive" score (0–1)
measures alignment quality: **higher is better** (≥ 0.6 = high confidence).

## Pairing rules (strict)

1. **Strand topology**: On the + strand, NTD must have a SMALLER start coordinate than CTD.
   On the − strand, NTD must have a LARGER start coordinate than CTD.
   Any pair violating this rule is IMPOSSIBLE.

2. **Span**: The locus (from the outer boundary of NTD to the outer boundary of CTD) must be
   between 1,000 bp and 150,000 bp. Pairs outside this range are INVALID.

3. **Type consistency**: NTD and CTD should ideally come from the same spidroin type
   (e.g. AcSp_NTD + AcSp_CTD). Cross-type pairs are flagged for review.

4. **No intervening same-terminus hits**: If another NTD or CTD of any type lies BETWEEN the
   proposed NTD and CTD (along the strand), the pair is highly suspicious and must be flagged
   with confidence=low and needs_review=true. Multiple CTDs clustered together almost never
   represent a single full-length gene.

5. **Spanning miniprot support**: If multiple non-overlapping miniprot alignments tile the
   region between the proposed NTD and CTD, this is STRONG evidence for a full-length gene.
   A single long-spanning alignment is also strong evidence.
   Zero spanning miniprot hits → confidence is reduced (mark as medium or low).

6. **E-value weighting**: Prefer pairing a high-confidence hit (e-value < 1e-10) with another
   high-confidence hit. If one side is low-confidence, reduce the pair confidence.

7. **Miniprot Positive score**: When multiple CTDs compete for one NTD (or vice versa), prefer
   the pairing whose interval contains miniprot hits with higher Positive scores.

## Output format

Respond with a single JSON object (no markdown fences). Schema:

{
  "pairs": [
    {
      "ntd_index": <int>,       // index into the NTD hits list (0-based)
      "ctd_index": <int>,       // index into the CTD hits list (0-based)
      "confidence": "high" | "medium" | "low",
      "supporting_evidence": [<string>, ...],  // e.g. ["type_consistent", "spanning_miniprot_3", "both_high_evalue"]
      "needs_review": <bool>,
      "reasoning": "<one sentence>"
    }
  ],
  "unpaired_ntds": [<int>, ...],   // NTD indices that should become N-terminal-only candidates
  "unpaired_ctds": [<int>, ...],   // CTD indices that should become C-terminal-only candidates
  "overall_reasoning": "<2-3 sentences summarising the window>",
  "needs_review": <bool>           // true if the overall window interpretation is uncertain
}

Include ALL NTD and CTD indices across pairs + unpaired lists. No index may be omitted or used twice.
"""

# Few-shot example embedded in the user prompt
FEW_SHOT_EXAMPLE = """\
### Example

NTD hits:
  [0] profile=AcSp_NTD  pos=71200000-71200420  strand=+  e_value=3.2e-45 (HIGH)

CTD hits:
  [0] profile=AcSp_CTD  pos=71250000-71250380  strand=+  e_value=8.5e-30 (HIGH)
  [1] profile=AcSp_CTD  pos=71260000-71260380  strand=+  e_value=2.1e-28 (HIGH)
  [2] profile=MiSp_CTD  pos=71280000-71280360  strand=+  e_value=5.6e-08 (medium)

Miniprot hits in window:
  ref=AcSp_CTD_1049aa  pos=71200100-71230000  strand=+  positive=0.88  identity=0.79  rank=1
  ref=AcSp_NTD_380aa   pos=71230100-71250380  strand=+  positive=0.75  identity=0.68  rank=1

Pre-computed pair features:
  NTD[0]-CTD[0]: type=AcSp/AcSp (consistent), dist=50380bp, spanning_miniprot=2, intervening_hits=1_CTD
  NTD[0]-CTD[1]: type=AcSp/AcSp (consistent), dist=60380bp, spanning_miniprot=2, intervening_hits=0
  NTD[0]-CTD[2]: type=AcSp/MiSp (CONFLICT), dist=80360bp, spanning_miniprot=0, intervening_hits=0

Expected output:
{
  "pairs": [
    {
      "ntd_index": 0,
      "ctd_index": 1,
      "confidence": "high",
      "supporting_evidence": ["type_consistent", "spanning_miniprot_2", "both_high_evalue", "no_intervening_hits"],
      "needs_review": false,
      "reasoning": "AcSp_NTD[0]+AcSp_CTD[1] is preferred over CTD[0] because CTD[0] has another CTD between it and the NTD; miniprot tiles the interval with positive≥0.75."
    }
  ],
  "unpaired_ntds": [],
  "unpaired_ctds": [0, 2],
  "overall_reasoning": "CTD[1] is the correct AcSp partner for NTD[0]: no intervening hits, strong miniprot support. CTD[0] is likely a repeat-region artefact (a CTD falls between it and the NTD). CTD[2] is a low-confidence MiSp hit with no miniprot support and a type conflict.",
  "needs_review": false
}
"""
