"""
agents.pairing — LangGraph-based NTD/CTD intelligent pairing module.

Handles ambiguous locus assembly cases where the greedy algorithm is insufficient,
e.g., multiple adjacent CTDs, mixed-type clusters, or uncertain pairings.
"""

from agents.pairing.pairing_agent import run_pairing_agent

__all__ = ["run_pairing_agent"]
