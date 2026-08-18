"""NCBI remote BLAST client — sequence-bait homology search.

This is the classic route to an unnamed paralog: take a known member's
full protein sequence and blastp it against NCBI nr, which returns homologs
*regardless of gene naming* — including "uncharacterized protein LOC…"
entries and misannotated fragments invisible to symbol search.

Like Foldseek, this is a bait client: it needs `query.include_sequence=True`
plus a prior sequence-DB result. The orchestrator harvests the longest
sequence per gene symbol and hands it in via `bait_sequences`; one BLAST job
runs per gene. Jobs go through Biopython's `NCBIWWW.qblast` (the public
BLAST URL API: submit → poll → XML), so expect 2–10 min per bait on the
shared queue.

Hit gene symbols are inferred from RefSeq titles: an explicit known-symbol
match (\"itpr2\", \"…receptor type 2\") keeps the known name; a family-root
match without a number becomes \"<root>-like\"; otherwise the LOC id or
\"unnamed\" — anything not in the known-paralog list flows into the
discovery scorer as a candidate.
"""

from __future__ import annotations

import re
from io import StringIO
from typing import Optional

from ..core.models import ProteinVariant, SearchQuery
from .base import DatabaseClient


class BlastClient(DatabaseClient):
    name = "BLAST"

    def __init__(
        self,
        email: str = "",
        api_key: str = "",
        timeout_s: int = 60,
        bait_sequences: Optional[dict[str, str]] = None,
        program: str = "blastp",
        database: str = "nr",
        hitlist_size: int = 100,
        expect: float = 1e-5,
        fetch_top_sequences: int = 15,
    ) -> None:
        super().__init__(email=email, api_key=api_key, timeout_s=timeout_s)
        self.bait_sequences = bait_sequences or {}
        self.program = program
        self.database = database
        self.hitlist_size = hitlist_size
        self.expect = expect
        self.fetch_top_sequences = fetch_top_sequences

    def search(self, query: SearchQuery) -> list[ProteinVariant]:
        if not self.bait_sequences:
            return []
        from Bio.Blast import NCBIWWW, NCBIXML
        if self.email:
            NCBIWWW.email = self.email

        variants: list[ProteinVariant] = []
        seen: set[str] = set()
        per_gene_cap = max(1, query.max_results_per_source)
        for symbol in query.gene_symbols:
            bait = self.bait_sequences.get(symbol)
            if not bait:
                continue
            entrez_query = (
                f'"{query.species}"[Organism]' if query.species else "(none)"
            )
            try:
                handle = NCBIWWW.qblast(
                    self.program, self.database, bait,
                    hitlist_size=self.hitlist_size,
                    expect=self.expect,
                    entrez_query=entrez_query,
                )
                record = NCBIXML.read(handle)
            except Exception:
                continue
            n_kept = 0
            for aln in record.alignments:
                if n_kept >= per_gene_cap:
                    break
                v = self._alignment_to_variant(aln, symbol, query.gene_symbols)
                if v and v.key() not in seen:
                    seen.add(v.key())
                    variants.append(v)
                    n_kept += 1
        if query.include_sequence and variants:
            self._fetch_sequences(variants)
        return variants

    # ---- parsing ---------------------------------------------------------

    def _alignment_to_variant(
        self, aln, bait_symbol: str, known_symbols: list[str]
    ) -> Optional[ProteinVariant]:
        title = (aln.hit_def or "").split(" >")[0].strip()
        if not title or not aln.hsps:
            return None
        hsp = aln.hsps[0]
        identity_pct = 100.0 * hsp.identities / max(1, hsp.align_length)
        species = ""
        m = re.search(r"\[([^\[\]]+)\]\s*$", title)
        if m:
            species = m.group(1)
        symbol = self._symbol_for(title, known_symbols)
        desc = (f"{title} · {self.program} vs {bait_symbol} bait "
                f"({identity_pct:.0f}% id over {hsp.align_length} aa, "
                f"e={hsp.expect:.1e})")
        return ProteinVariant(
            source=self.name,
            accession=aln.accession or aln.hit_id,
            gene_symbol=symbol,
            species=species,
            length_aa=aln.length or None,
            description=desc,
            gene_id=aln.hit_id or "",
            url=f"https://www.ncbi.nlm.nih.gov/protein/{aln.accession}",
            raw={
                "identity_pct": f"{identity_pct:.1f}",
                "evalue": f"{hsp.expect:.2e}",
                "bit_score": f"{hsp.bits:.0f}",
                "align_length": str(hsp.align_length),
                "bait_symbol": bait_symbol,
            },
        )

    @staticmethod
    def _symbol_for(title: str, known_symbols: list[str]) -> str:
        """Best-effort gene symbol from a RefSeq/GenBank protein title."""
        tl = title.lower()
        for sym in known_symbols:
            s = sym.lower()
            m = re.match(r"([a-z]+)(\d[a-z0-9]*)$", s)
            if s in tl:
                return sym
            if m and m.group(1) in tl and f"component {m.group(2)}" in tl:
                return sym
        roots = {re.match(r"[a-z]+", s.lower()).group(0)
                 for s in known_symbols if re.match(r"[a-z]+", s.lower())}
        for root in roots:
            if root in tl:
                return f"{root}-like"
        m = re.search(r"\bLOC\d+\b", title)
        if m:
            return m.group(0)
        return "unnamed"

    # ---- sequence fetch --------------------------------------------------

    def _fetch_sequences(self, variants: list[ProteinVariant]) -> None:
        """Fetch full sequences for the top hits via Entrez (one batch)."""
        try:
            from Bio import Entrez, SeqIO
        except ImportError:
            return
        targets = [v for v in variants if not v.sequence][: self.fetch_top_sequences]
        if not targets:
            return
        if self.email:
            Entrez.email = self.email
        try:
            handle = Entrez.efetch(
                db="protein",
                id=",".join(v.accession for v in targets),
                rettype="fasta", retmode="text",
            )
            text = handle.read()
            handle.close()
        except Exception:
            return
        by_acc: dict[str, str] = {}
        for rec in SeqIO.parse(StringIO(text), "fasta"):
            acc = rec.id.split("|")[-1] if "|" in rec.id else rec.id
            by_acc[acc.split(".")[0]] = str(rec.seq)
        for v in targets:
            seq = by_acc.get(v.accession.split(".")[0])
            if seq:
                v.sequence = seq
                if not v.length_aa:
                    v.length_aa = len(seq)
