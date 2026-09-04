"""S3 — the committed seed manifests for `itpr.hmm` and `ryr.hmm`.

Two profiles, not one. D14 says ITPR vs RYR is a positive test at every
stage; at HMM scale that test is *best-profile assignment with a bit-score
margin*, and a margin needs two comparable profiles built to the same design.

Selection rules (applied by hand here, checked in code by s3_build_seed.py):

  R1  Every seed comes from the S2 census (`results/census_v2/census_v2.tsv`)
      and carries that census's call. The call is the architecture rule of
      D14b, which S2 audited against 6,191 gene symbols with 0 disagreements
      — so the seeds are labelled by a rule, not by a name, and the profile
      assignment they support stays an independent sequence-level test.
  R2  One seed per species per paralog. No species contributes more than
      three seeds to a profile, so no well-sequenced genome dominates it.
  R3  Complete architecture (5/5 signatures) unless the seed is there for
      phylogenetic breadth no complete record covers, in which case its
      `arch` column records what it actually carries and `note` says why.
  R4  Length within the family band unless noted. The band supports the
      choice; it never makes it (D14).
  R5  The two sets are disjoint by accession and by species-paralog, and
      the ITPR set is built to span the same taxonomic depth as the RyR set
      so neither profile is the broader one by construction.

The one deliberate exception to R3 is Q9NA13, *Dictyostelium* iplA — a
characterised IP3 receptor carrying 2 of the 5 signatures, which S2 leaves
`unassigned` and which a signature census cannot find (D21). It is in the
ITPR seed set precisely so the profile can find its relatives.
"""

# (accession, clade, arch_expected, note) — arch_expected is the census
# n_itpr_arch value this manifest was written against; s3_build_seed.py
# fails if the census disagrees, so a census re-run cannot silently change
# what the profile was built from.
ITPR_SEEDS: list[tuple[str, str, int, str]] = [
    # --- ITPR1 -----------------------------------------------------------
    ("Q14643", "ITPR1", 5, "human, Swiss-Prot, the reference"),
    ("P11881", "ITPR1", 5, "mouse, Swiss-Prot"),
    ("A0A8V0ZJT0", "ITPR1", 5, "chicken"),
    ("A0A8J0T3U5", "ITPR1", 5, "Xenopus tropicalis"),
    ("A0A8M2BDG0", "ITPR1", 5, "zebrafish itpr1a (3R duplicate)"),
    ("H3A1P5", "ITPR1", 5, "coelacanth — the sarcopterygian outgroup"),
    # --- ITPR2 -----------------------------------------------------------
    ("Q14571", "ITPR2", 5, "human, Swiss-Prot"),
    ("Q9Z329", "ITPR2", 5, "mouse, Swiss-Prot"),
    ("F1P1X4", "ITPR2", 5, "chicken"),
    ("A0A6I8RGU0", "ITPR2", 5, "Xenopus tropicalis"),
    ("A0A8M9PQQ3", "ITPR2", 5, "zebrafish"),
    # --- ITPR3 -----------------------------------------------------------
    ("Q14573", "ITPR3", 5, "human, Swiss-Prot"),
    ("P70227", "ITPR3", 5, "mouse, Swiss-Prot"),
    ("A0A8V1ACE6", "ITPR3", 5, "chicken"),
    ("A0A8J0SYK9", "ITPR3", 5, "Xenopus tropicalis"),
    ("A0A8M6Z266", "ITPR3", 5, "zebrafish"),
    # --- early vertebrates, unnamed loci ---------------------------------
    ("A0A4W3JHJ0", "vertebrate_basal", 5, "elephant shark — chondrichthyan"),
    ("A0AAJ7U4X1", "vertebrate_basal", 5, "sea lamprey — cyclostome, pre-2R"),
    ("W5MYT2", "vertebrate_basal", 5, "spotted gar — non-teleost actinopt."),
    # --- invertebrates, the single-Itpr grade ----------------------------
    ("P29993", "invertebrate", 5, "Drosophila Itpr, Swiss-Prot"),
    ("Q9Y0A1", "invertebrate", 5, "C. elegans itr-1, Swiss-Prot"),
    ("Q8WSR4", "invertebrate", 5, "starfish IP3R, Swiss-Prot"),
    ("A0A7M7MY00", "invertebrate", 5, "sea urchin — echinoderm"),
    ("A0A6P8IUM8", "invertebrate", 5, "sea anemone — cnidarian"),
    ("A0AAN0JIF7", "invertebrate", 5, "Amphimedon — sponge, basal metazoan"),
    ("A0A1S3H6V5", "invertebrate", 4, "Lingula — brachiopod; 4/5, lophotroch. breadth"),
    ("B3RXX7", "invertebrate", 4, "Trichoplax — placozoan; 4/5, basal breadth"),
    # --- non-metazoan grade ----------------------------------------------
    ("L8GF85", "non_metazoan", 5, "Acanthamoeba — Amoebozoa/Discosea"),
    ("Q9NA13", "non_metazoan", 2, "Dictyostelium iplA — characterised receptor, "
                                  "2/5 signatures; the R3 exception (D21)"),
    ("A0AA88KEU5", "non_metazoan", 5, "Naegleria — Discoba"),
    ("A0A058ZDJ2", "non_metazoan", 5, "Fonticula alba — Holomycota"),
    ("A0A8T1W1U0", "non_metazoan", 5, "Phytophthora — SAR/oomycete"),
    ("A0ABR2W6A1", "non_metazoan", 5, "Basidiobolus — early-diverging fungus"),
    ("A0AAE0LDZ1", "non_metazoan", 3, "Cymbomonas — Chlorophyta; 3/5, the green "
                                      "algal record the plant question rests on"),
]

RYR_SEEDS: list[tuple[str, str, int, str]] = [
    # --- RYR1 ------------------------------------------------------------
    ("P21817", "RYR1", 5, "human, Swiss-Prot"),
    ("E9PZQ0", "RYR1", 5, "mouse, Swiss-Prot"),
    ("A0A6I8Q6P7", "RYR1", 5, "Xenopus tropicalis"),
    ("A0A8M6YTI2", "RYR1", 5, "zebrafish ryr1b"),
    # --- RYR2 ------------------------------------------------------------
    ("Q92736", "RYR2", 5, "human, Swiss-Prot"),
    ("E9Q401", "RYR2", 5, "mouse, Swiss-Prot"),
    ("A0A8V0Y8K1", "RYR2", 5, "chicken"),
    ("A0AC58GZ77", "RYR2", 5, "zebrafish ryr2a"),
    # --- RYR3 ------------------------------------------------------------
    ("Q15413", "RYR3", 5, "human, Swiss-Prot"),
    ("A2AGL3", "RYR3", 5, "mouse, Swiss-Prot"),
    ("A0A8V0ZGY9", "RYR3", 5, "chicken"),
    ("A0A8J1INJ9", "RYR3", 5, "Xenopus tropicalis"),
    # --- early vertebrates -----------------------------------------------
    ("A0A4W3I1N0", "vertebrate_basal", 5, "elephant shark — chondrichthyan"),
    ("A0AAJ7XD62", "vertebrate_basal", 5, "sea lamprey — cyclostome"),
    # --- invertebrates ---------------------------------------------------
    ("Q24498", "invertebrate", 5, "Drosophila RyR, Swiss-Prot"),
    ("Q94279", "invertebrate", 5, "C. elegans unc-68, Swiss-Prot"),
    ("Q8TA74", "invertebrate", 5, "sea urchin — echinoderm"),
    ("X1Z8I7", "invertebrate", 5, "Capitella — annelid"),
    ("A0A1S3I436", "invertebrate", 5, "Lingula — brachiopod"),
    ("B3RK74", "invertebrate", 5, "Trichoplax — placozoan"),
    ("A0AAV7K8I1", "invertebrate", 4, "Oopsacas — glass sponge; 4/5 shared signatures but 3 RyR-specific ones, so the call is firm"),
    # --- non-metazoan ----------------------------------------------------
    ("A0A0D2WXI6", "non_metazoan", 5, "Capsaspora — filasterean, sister to Metazoa"),
]

# Both profiles are built from these sets and nothing else.
PROFILES = {
    "itpr": {"seeds": ITPR_SEEDS, "hmm_name": "ITPR_family",
             "expect_call": "ITPR"},
    "ryr": {"seeds": RYR_SEEDS, "hmm_name": "RYR_family",
            "expect_call": "RYR"},
}

# jackhmmer single-sequence seeds (brief step 3: human ITPR1, an
# invertebrate Itpr, a protist ITPR). Tag → (accession, why).
JACKHMMER_SEEDS: dict[str, tuple[str, str]] = {
    "itpr1_human": ("Q14643", "the reference paralog"),
    "itpr_fly": ("P29993", "the single-Itpr invertebrate grade"),
    "itpr_acanthamoeba": ("L8GF85", "the non-metazoan grade — a complete "
                                    "architecture outside Metazoa"),
}

MIN_SEED_LEN = 2000       # ITPR band floor; RyR seeds are far above it
MAX_SPECIES_PER_PROFILE = 3   # R2
