# S20 — the non-vertebrate sweep: the family's true range

_Rendered from the committed tables on 2026-09-05 10:49 by `scripts/s20_report.py` (D13)._

S3 swept the vertebrate reference proteomes and found where the three paralogs live. This task asks the opposite question — how far the family reaches — and the one it was set up to answer: the databases hold a few tens of plant and fungal records while *Arabidopsis* and *S. cerevisiae* hold none. Which of those is a fact about genomes and which about databases?

## The declared search space

| group | proteomes | proteins | residues | sampling | swept |
|---|---|---|---|---|---|
| metazoa_nonvert | 546 | 11,396,529 | 4.56 G | every reference proteome | yes |
| fungi | 1,527 | 17,339,036 | 7.79 G | every reference proteome | yes |
| viridiplantae | 432 | 15,114,263 | 5.81 G | every reference proteome | yes |
| protista_other | 252 | 3,636,524 | 1.78 G | every reference proteome | yes |
| archaea | 634 | 1,754,814 | 0.50 G | every reference proteome | yes |
| bacteria_genus | 3,537 | 13,903,732 | 4.48 G | largest_per_genus | yes |


**6,928 reference proteomes, 63,144,898 proteins, 24.93 G residues.** The four eukaryote groups partition Eukaryota with S3's `vertebrata`, so no proteome is swept twice and the two sweeps' denominators add.


Bacteria are sampled and archaea are not: one proteome per genus (first token of the organism name), the one with the most proteins — the most sensitive member of each genus, so the negative claim is made in the places most likely to break it. Archaea's reference set is small enough to take whole, so no sampling caveat attaches to it.


23 proteome(s) UniProt lists are not published in the current release FTP tree and 404 permanently; they are excluded from the denominator above and recorded in `proteome_unavailable_<group>.tsv` with what they took out of it.


## The instrument

The same two profiles S3 built and calibrated — `itpr.hmm` and `ryr.hmm` — with the same margin assignment: a target is called only when the winning profile beats the loser by more than 10% of its own score, clears 30 bits, and spans at least 200 match states (D22). Reusing S3's instrument rather than building a new one is what makes the vertebrate and non-vertebrate numbers comparable at all.


**One search, two sensitivities.** Every `hmmsearch` ran at `-E 10` and the primary call was taken by filtering the same output at E ≤ 1e-5. `-E` is a reporting threshold and does not touch hmmsearch's acceleration filters, so the strict set is exactly what a strict run would have produced and the relaxed set is a superset of it from the same search — not a second experiment that might have differed some other way.


The relaxed panel adds the family's four Pfam domain models:

| Pfam | name | match states | why it is in the panel |
|---|---|---|---|
| PF08709 | Ins145_P3_rec | 213 | IP3-binding core — the signature that names the family |
| PF02815 | MIR | 185 | MIR — shared with POMT1/2 |
| PF01365 | RYDR_ITPR | 201 | RYDR_ITPR — shared with the ryanodine receptors |
| PF08454 | RIH_assoc | 99 | RIH_assoc — shared with the ryanodine receptors |

| group | itpr.hmm targets | ryr.hmm targets | → ITPR | → RYR | → unassigned | search min |
|---|---|---|---|---|---|---|
| metazoa_nonvert | 2,513 | 8,125 | 1,194 | 1,156 | 6,257 | 18 |
| fungi | 44 | 2,942 | 35 | 20 | 2,897 | 22 |
| viridiplantae | 71 | 13,770 | 54 | 1 | 13,736 | 13 |
| protista_other | 798 | 2,580 | 729 | 121 | 1,871 | 9 |
| archaea | 0 | 7 | 0 | 0 | 7 | 1 |
| bacteria_genus | 0 | 59 | 0 | 0 | 59 | 6 |


## The family's range across the eukaryotic tree

6,928 reference proteomes were swept with both family profiles at the S3 protocol's threshold, and **662/6928 (10%) carry at least one ITPR call**.

| group | proteomes swept | with an ITPR call | ITPR records |
|---|---|---|---|
| Metazoa (non-vertebrate) | 546 | 509/546 (93%) | 1,194 |
| SAR | 167 | 65/167 (39%) | 602 |
| Fungi | 1,527 | 28/1527 (2%) | 35 |
| Discoba | 35 | 23/35 (66%) | 58 |
| Viridiplantae | 432 | 15/432 (3%) | 54 |
| Eukaryota (other) | 37 | 13/37 (35%) | 57 |
| Amoebozoa | 13 | 9/13 (69%) | 12 |
| Archaea | 634 | 0/634 (0%) | 0 |
| Bacteria | 3,537 | 0/3537 (0%) | 0 |


Across every group, **45 clades of the 135 swept carry an ITPR call**. Clade is the phylum where UniProt states one and the next-deepest named group where it does not — the choanoflagellates, filastereans, cryptophytes and apusozoans have no phylum rank, and those are precisely the lineages a range result about this family has to be able to name. Most-covered first:

| clade | proteomes | with an ITPR call |
|---|---|---|
| Arthropoda | 311 | 292/311 (94%) |
| Nematoda | 110 | 106/110 (96%) |
| Mollusca | 32 | 32/32 (100%) |
| Platyhelminthes | 40 | 26/40 (65%) |
| Oomycota | 40 | 26/40 (65%) |
| Euglenozoa | 30 | 20/30 (67%) |
| Mucoromycota | 34 | 18/34 (53%) |
| Chlorophyta | 48 | 15/48 (31%) |
| Ciliophora | 16 | 15/16 (94%) |
| Cnidaria | 14 | 14/14 (100%) |
| Evosea | 12 | 8/12 (67%) |
| Dinophyceae | 10 | 8/10 (80%) |
| Chordata | 8 | 8/8 (100%) |
| Rotifera | 8 | 8/8 (100%) |
| Chytridiomycota | 16 | 6/16 (38%) |
| Annelida | 6 | 6/6 (100%) |
| Bolidophyceae | 6 | 6/6 (100%) |
| Echinodermata | 5 | 5/5 (100%) |
| Pelagophyceae | 4 | 4/4 (100%) |
| Haptophyta | 4 | 4/4 (100%) |
| Heterolobosea | 4 | 3/4 (75%) |
| Porifera | 3 | 3/3 (100%) |
| Perkinsozoa | 3 | 3/3 (100%) |
| Tardigrada | 2 | 2/2 (100%) |
| Basidiobolomycota | 2 | 2/2 (100%) |


Census v5 holds **8,809 ITPR records across 1,402 taxa**.


## D14 outside the vertebrates

28,137 targets were scored by at least one profile. **2,769 were scored by both above the 30-bit floor** — the only ones where the two families can be said to compete at all — and of those **1/2769 (0%) fall inside D7's 10% no-call band**, where this instrument declines to choose.

| outcome | targets |
|---|---|
| declined: shared module only | 24,725 |
| called | 3,310 |
| declined: under the bit-score floor | 101 |
| declined: inside the no-call band | 1 |


The 200-state floor (D22) is doing most of the declining: a protein sharing one module with a full-length channel model scores against it, and outside the vertebrates that is the common case rather than the exception.


### What each call actually rests on

A call is only as good as how much of the model it spans, and the two profiles are very different sizes — `itpr.hmm` is 2,684 match states, `ryr.hmm` 4,930 — so 200 states is 7 % of one model and 4 % of the other. Breaking the calls out by evidence class is what stops that asymmetry from being read as biology:

| group | call | records | architecture (≥50 % of the model) | partial | median coverage |
|---|---|---|---|---|---|
| fungi | ITPR | 35 | 28 | 7 | 91% |
| fungi | RYR | 20 | 0 | 20 | 5% |
| metazoa_nonvert | ITPR | 1,194 | 617 | 577 | 54% |
| metazoa_nonvert | RYR | 1,156 | 441 | 715 | 27% |
| protista_other | ITPR | 729 | 466 | 263 | 75% |
| protista_other | RYR | 121 | 2 | 119 | 4% |
| viridiplantae | ITPR | 54 | 23 | 31 | 41% |
| viridiplantae | RYR | 1 | 0 | 1 | 5% |


**1134/2012 (56%) of the ITPR calls span at least half the model, against 443/1298 (34%) of the RYR ones.** The RYR calls out here are overwhelmingly module-level matches to a 4,930-state channel model by proteins a few hundred to two thousand residues long — the shared architecture D14 exists to see through, not ryanodine receptors. Read as gene counts they would put RyRs in fungi and green algae; read as coverage they do not.


The 443 that *are* architecture-level are worth naming, because where the sister family is tells you when the two families split:

| accession | species | length | ryr bits | coverage | margin | a profile seed? |
|---|---|---|---|---|---|---|
| A0AAD8DKG6 | *Mythimna separata* | 5857 | 7536.3 | 100% | 91% | no |
| A0A1J1IT56 | *Clunio marinus* | 6190 | 7484.5 | 96% | 90% | no |
| Q24498 | *Drosophila melanogaster* | 5127 | 7265.2 | 100% | 91% | yes — circular |
| A0A9P9YFB4 | *Drosophila gunungcola* | 5165 | 7257.9 | 100% | 91% | no |
| A0A6P8K2C8 | *Drosophila mauritiana* | 5135 | 7230.2 | 100% | 90% | no |
| A0A0Q9XB63 | *Drosophila mojavensis* | 5137 | 7222.4 | 100% | 91% | no |
| A0A6J2SWD8 | *Drosophila lebanonensis* | 5133 | 7221.3 | 100% | 90% | no |
| A0AB39ZTW8 | *Drosophila suzukii* | 5135 | 7220.7 | 100% | 90% | no |
| B4JVD7 | *Drosophila grimshawi* | 5174 | 7217.1 | 100% | 90% | no |
| A0A6P4FG33 | *Drosophila rhopaloa* | 5130 | 7216.7 | 100% | 90% | no |
| A0A0R1DMC8 | *Drosophila yakuba* | 5146 | 7215.9 | 100% | 90% | no |
| A0A3B0J281 | *Drosophila guanche* | 5142 | 7209.2 | 100% | 90% | no |
| A0ABM1PIB9 | *Drosophila arizonae* | 5149 | 7208.8 | 100% | 90% | no |
| A0AAD4K9G6 | *Drosophila rubida* | 5131 | 7207.4 | 100% | 90% | no |
| A0A6I8VTB1 | *Drosophila pseudoobscura pseudoobscura* | 5143 | 7196.2 | 100% | 90% | no |
| A0A9C6W5B0 | *Drosophila albomicans* | 5148 | 7195.9 | 100% | 90% | no |
| A0A1I8MDA8 | *Musca domestica* | 5131 | 7195.4 | 100% | 90% | no |
| A0ABM4G9K7 | *Drosophila kikkawai* | 5117 | 7195.1 | 100% | 90% | no |
| X1Z8I7 | *Capitella teleta* | 5038 | 7194.7 | 98% | 91% | yes — circular |
| A0A0Q5WNI3 | *Drosophila erecta* | 5134 | 7193.9 | 100% | 90% | no |
| A0A0Q9W2M0 | *Drosophila virilis* | 5143 | 7193.6 | 100% | 90% | no |
| A0A0L0BUZ9 | *Lucilia cuprina* | 5136 | 7192.3 | 100% | 90% | no |
| B3MFE1 | *Drosophila ananassae* | 5106 | 7186.7 | 100% | 90% | no |
| A0A1I8P9C9 | *Stomoxys calcitrans* | 5129 | 7184.7 | 100% | 90% | no |
| A0A0Q9X0Z3 | *Drosophila willistoni* | 5142 | 7179.8 | 100% | 90% | no |
| A0ABQ9EYP8 | *Tegillarca granosa* | 5194 | 7168.9 | 97% | 89% | no |
| A0ABM3JCZ7 | *Bactrocera dorsalis* | 5142 | 7167.9 | 100% | 90% | no |
| A0A6J1LN04 | *Drosophila hydei* | 5132 | 7164.6 | 100% | 90% | no |
| A0A484BQ74 | *Drosophila navojoa* | 5144 | 7164.3 | 100% | 90% | no |
| A0A7R8YZN5 | *Hermetia illucens* | 5141 | 7163.1 | 100% | 90% | no |
| A0AAU9G7M0 | *Drosophila madeirensis* | 5148 | 7147.6 | 100% | 90% | no |
| A0A8S3X2A4 | *Parnassius apollo* | 5166 | 7147.3 | 100% | 90% | no |
| A0A9C5Z257 | *Glossina fuscipes* | 5125 | 7146.7 | 100% | 90% | no |
| A0A9P0TSB5 | *Pieris brassicae* | 5162 | 7145.9 | 100% | 90% | no |
| A0A9N9R5I0 | *Diatraea saccharalis* | 5122 | 7139.3 | 100% | 90% | no |
| A0A8J9VTZ2 | *Brenthis ino* | 5139 | 7136.9 | 100% | 90% | no |
| A0A182PDB6 | *Anopheles epiroticus* | 5107 | 7132.9 | 100% | 90% | no |
| A0A921YRC0 | *Manduca sexta* | 5127 | 7129.5 | 100% | 90% | no |
| A0ACE0G7E8 | *Anthophora plagiata* | 5492 | 7127.8 | 100% | 90% | no |
| A0AAD9RL42 | *Odynerus spinipes* | 5110 | 7125.1 | 100% | 90% | no |
| A0ABN8B934 | *Chilo suppressalis* | 5128 | 7119.4 | 100% | 90% | no |
| A0A9P0A197 | *Bemisia tabaci* | 5135 | 7114.4 | 100% | 90% | no |
| A0AAN9YVG3 | *Gryllus longicercus* | 5148 | 7114.1 | 100% | 90% | no |
| A0A1S4H4U6 | *Anopheles gambiae* | 5111 | 7112.5 | 100% | 90% | no |
| A0A182FUC1 | *Anopheles albimanus* | 5119 | 7111.3 | 100% | 90% | no |
| A0ABP1N951 | *Xylocopa violacea* | 5113 | 7111.2 | 100% | 90% | no |
| A0A903Z5W0 | *Anopheles quadriannulatus* | 5106 | 7111.0 | 100% | 90% | no |
| A0A6E8W3D1 | *Anopheles coluzzii* | 5109 | 7110.0 | 100% | 90% | no |
| A0ABM3GP40 | *Neodiprion lecontei* | 5128 | 7109.0 | 100% | 90% | no |
| A0ABM3M9C1 | *Galleria mellonella* | 5164 | 7106.0 | 100% | 90% | no |
| A0A9P0IA62 | *Spodoptera littoralis* | 5121 | 7105.0 | 100% | 90% | no |
| A0A6P3XEL1 | *Dinoponera quadriceps* | 5107 | 7104.0 | 100% | 90% | no |
| A0A2A3ESZ8 | *Apis cerana cerana* | 5051 | 7103.7 | 100% | 90% | no |
| A0A7E5WRJ5 | *Trichoplusia ni* | 5114 | 7103.7 | 100% | 90% | no |
| A0AAG5DRS8 | *Anopheles atroparvus* | 5111 | 7103.4 | 100% | 90% | no |
| A0A6I8T8R3 | *Aedes aegypti* | 5124 | 7103.3 | 100% | 90% | no |
| A0A9P0ARV5 | *Brassicogethes aeneus* | 5109 | 7103.2 | 100% | 90% | no |
| A0A8R2M7J7 | *Bombyx mori* | 5167 | 7102.3 | 100% | 90% | no |
| A0A6J2KJN8 | *Bombyx mandarina* | 5161 | 7101.7 | 100% | 90% | no |
| A0ACC2Q884 | *Mythimna loreyi* | 5118 | 7101.5 | 100% | 90% | no |
| A0AAV2PB15 | *Lasius platythorax* | 5118 | 7100.6 | 100% | 90% | no |
| A0A182WFR7 | *Anopheles minimus* | 5113 | 7100.2 | 100% | 90% | no |
| A0A9P0MMS1 | *Acanthoscelides obtectus* | 5142 | 7097.8 | 100% | 90% | no |
| A0AAE1GY24 | *Frankliniella fusca* | 5165 | 7095.6 | 100% | 90% | no |
| A0A9N9T9N1 | *Phyllotreta striolata* | 5111 | 7093.6 | 100% | 90% | no |
| A0AAW1A662 | *Tetragonisca angustula* | 5132 | 7093.6 | 100% | 90% | no |
| A0AAV1LIQ1 | *Parnassius mnemosyne* | 5092 | 7093.3 | 100% | 90% | no |
| A0A9R0DY57 | *Spodoptera frugiperda* | 5165 | 7091.3 | 100% | 90% | no |
| E2BZJ9 | *Harpegnathos saltator* | 5080 | 7090.6 | 100% | 90% | no |
| A0A7M7GWS7 | *Apis mellifera* | 5142 | 7089.2 | 100% | 90% | no |
| A0ABQ7PQG9 | *Plutella xylostella* | 5164 | 7089.1 | 100% | 90% | no |
| A0AAJ7CFC2 | *Cephus cinctus* | 5147 | 7088.7 | 100% | 90% | no |
| A0AAW2FGZ5 | *Cardiocondyla obscurior* | 5136 | 7085.5 | 100% | 90% | no |
| A0ABM3LVF0 | *Bicyclus anynana* | 5150 | 7084.8 | 100% | 90% | no |
| A0A158P0A3 | *Atta cephalotes* | 5126 | 7083.3 | 100% | 90% | no |
| A0ABM4AVN2 | *Vanessa tameamea* | 5162 | 7082.4 | 100% | 90% | no |
| A0ABM1HYJ6 | *Polistes dominula* | 5161 | 7082.3 | 100% | 90% | no |
| A0AA40G5E2 | *Melipona bicolor* | 5069 | 7081.4 | 100% | 90% | no |
| A0A6P8L8T1 | *Bombus impatiens* | 5158 | 7081.0 | 100% | 90% | no |
| A0A9C6SJC6 | *Bombus terrestris* | 5158 | 7081.0 | 100% | 90% | no |
| A0A834M383 | *Rhynchophorus ferrugineus* | 5144 | 7080.2 | 100% | 90% | no |
| A0A8B8GT52 | *Sipha flava* | 5102 | 7079.7 | 100% | 90% | no |
| A0A9P0D3I1 | *Psylliodes chrysocephalus* | 5121 | 7079.7 | 100% | 90% | no |
| A0A182Y6U7 | *Anopheles stephensi* | 5085 | 7079.6 | 100% | 90% | no |
| A0A6P8LNL7 | *Bombus bifarius* | 5158 | 7079.4 | 100% | 90% | no |
| A0A9N9RUB9 | *Chironomus riparius* | 5105 | 7077.1 | 100% | 90% | no |
| A0A9N9TDL6 | *Diabrotica balteata* | 5121 | 7075.4 | 100% | 90% | no |
| A0A6J2XGA9 | *Sitophilus oryzae* | 5122 | 7073.8 | 100% | 90% | no |
| A0A8N1SB51 | *Pogonomyrmex barbatus* | 5130 | 7073.1 | 100% | 90% | no |
| A0A084WAS3 | *Anopheles sinensis* | 5077 | 7072.3 | 100% | 90% | no |
| A0A6J1Q3K4 | *Temnothorax curvispinosus* | 5140 | 7069.5 | 100% | 90% | no |
| A0A9J7EHB1 | *Spodoptera litura* | 5098 | 7068.0 | 99% | 90% | no |
| A0A195B655 | *Atta colombica* | 5074 | 7064.3 | 100% | 90% | no |
| A0A151IVU0 | *Trachymyrmex cornetzi* | 5093 | 7062.9 | 100% | 90% | no |
| A0A1A9W9U7 | *Glossina brevipalpis* | 5096 | 7062.0 | 100% | 90% | no |
| A0ABM5L3U1 | *Diabrotica virgifera virgifera* | 5130 | 7061.1 | 100% | 90% | no |
| A0A9P0GQZ3 | *Phaedon cochleariae* | 5115 | 7060.4 | 100% | 90% | no |
| A0ABD1FID6 | *Hypothenemus hampei* | 5104 | 7059.8 | 100% | 90% | no |
| A0A7F5QZ06 | *Agrilus planipennis* | 5123 | 7057.8 | 99% | 90% | no |
| A0ABD2P1J4 | *Cryptolaemus montrouzieri* | 5108 | 7056.6 | 100% | 90% | no |
| A0AAR5QBW0 | *Dendroctonus ponderosae* | 5120 | 7051.8 | 100% | 90% | no |
| E0VEK3 | *Pediculus humanus subsp. corporis* | 5058 | 7051.7 | 100% | 90% | no |
| A0A9W2ZVE5 | *Biomphalaria glabrata* | 5224 | 7050.5 | 100% | 90% | no |
| A0A6P8Z5B9 | *Thrips palmi* | 5163 | 7048.4 | 100% | 90% | no |
| A0A151K2J5 | *Cyphomyrmex costatus* | 5133 | 7048.3 | 100% | 90% | no |
| A0A026W7Q8 | *Ooceraea biroi* | 5061 | 7047.2 | 100% | 90% | no |
| A0A9P0J9A2 | *Aphis gossypii* | 5104 | 7046.9 | 100% | 90% | no |
| A0A6G0U418 | *Aphis glycines* | 5104 | 7046.8 | 100% | 90% | no |
| A0AAN9TKK1 | *Parthenolecanium corni* | 5107 | 7045.6 | 100% | 90% | no |
| A0AAE0TEQ7 | *Potamilus streckersoni* | 5165 | 7044.9 | 100% | 90% | no |
| A0A8S1C2V6 | *Cloeon dipterum* | 5116 | 7043.3 | 100% | 90% | no |
| A0A9N9MGQ0 | *Ceutorhynchus assimilis* | 5124 | 7042.6 | 100% | 90% | no |
| A0AAV7XDJ1 | *Megalurothrips usitatus* | 5138 | 7040.6 | 100% | 90% | no |
| A0ACF7J0R1 | *Uroleucon formosanum* | 5101 | 7040.0 | 100% | 90% | no |
| A0ABD3WLX5 | *Sinanodonta woodiana* | 5152 | 7038.9 | 100% | 90% | no |
| A0AAV2I6U5 | *Lymnaea stagnalis* | 5191 | 7038.2 | 100% | 90% | no |
| A0ACF8WPY9 | *Semiaphis heraclei* | 5101 | 7037.3 | 100% | 90% | no |
| A0A8R2H7N0 | *Acyrthosiphon pisum* | 5101 | 7037.2 | 100% | 90% | no |
| A0AAV0WCS0 | *Macrosiphum euphorbiae* | 5101 | 7036.7 | 100% | 90% | no |
| A0A195FCN9 | *Trachymyrmex septentrionalis* | 5023 | 7031.4 | 100% | 90% | no |
| A0A182ND92 | *Anopheles dirus* | 5078 | 7030.7 | 99% | 90% | no |
| A0AAN9GD88 | *Littorina saxatilis* | 5241 | 7028.4 | 100% | 90% | no |
| A0AAJ7N459 | *Ceratina calcarata* | 5153 | 7025.4 | 100% | 90% | no |
| A0ABD1DFM7 | *Culex pipiens pipiens* | 5051 | 7025.1 | 100% | 90% | no |
| A0A182I5V4 | *Anopheles arabiensis* | 5081 | 7024.1 | 99% | 90% | no |
| A0A7M7T6K6 | *Nasonia vitripennis* | 5176 | 7023.2 | 100% | 90% | no |
| A0AA36B9A0 | *Octopus vulgaris* | 5170 | 7016.1 | 100% | 90% | no |
| A0A9P0H1M9 | *Nezara viridula* | 5073 | 7015.7 | 100% | 90% | no |
| A0A232EYM9 | *Trichomalopsis sarcophagae* | 5113 | 7012.0 | 100% | 90% | no |
| A0A8I6SMR1 | *Cimex lectularius* | 5121 | 7011.4 | 100% | 90% | no |
| A0AAN7ZFM5 | *Pyrocoelia pectoralis* | 5103 | 7005.3 | 100% | 90% | no |
| A0A5N4A4Q8 | *Photinus pyralis* | 5100 | 7005.1 | 100% | 90% | no |
| A0A836FDY8 | *Acromyrmex charruanus* | 5128 | 7005.0 | 100% | 90% | no |
| A0AAE0Y8B9 | *Elysia crispata* | 5260 | 6996.1 | 100% | 89% | no |
| A0A834KKM1 | *Vespula vulgaris* | 5055 | 6995.7 | 98% | 90% | no |
| A0ACB9TXB3 | *Holotrichia oblita* | 5073 | 6994.9 | 100% | 90% | no |
| A0ACF7KUW5 | *Mactra antiquata* | 5161 | 6993.3 | 100% | 90% | no |
| A0A836ED53 | *Pseudoatta argentina* | 5138 | 6992.7 | 100% | 90% | no |
| A0A182MMK6 | *Anopheles culicifacies* | 5039 | 6991.9 | 99% | 90% | no |
| A0A836FFX1 | *Acromyrmex heyeri* | 5136 | 6986.1 | 100% | 90% | no |
| T1ICG6 | *Rhodnius prolixus* | 5107 | 6983.4 | 100% | 90% | no |
| A0A835CM71 | *Aphidius gifuensis* | 5145 | 6983.2 | 100% | 90% | no |
| A0AA88Y1D7 | *Pinctada imbricata* | 5203 | 6979.1 | 100% | 90% | no |
| A0A8B8DQA7 | *Crassostrea virginica* | 5210 | 6972.8 | 100% | 90% | no |
| A0ABN7B5W3 | *Nesidiocoris tenuis* | 5083 | 6972.3 | 98% | 90% | no |
| A0A835GAU2 | *Spodoptera exigua* | 5106 | 6971.9 | 100% | 90% | no |
| A0ABD2VZJ9 | *Trichogramma kaykai* | 5126 | 6967.7 | 100% | 89% | no |
| A0ACF9AK27 | *Megaselia abdita* | 5125 | 6961.9 | 100% | 90% | no |
| A0A9R1TY94 | *Fopius arisanus* | 5124 | 6949.9 | 100% | 90% | no |
| A0A1B0G1T7 | *Glossina morsitans morsitans* | 5092 | 6949.5 | 98% | 90% | no |
| A0ABM1MTS5 | *Nicrophorus vespilloides* | 5054 | 6940.9 | 100% | 91% | no |
| A0AAN7SKY1 | *Aquatica leii* | 5086 | 6939.6 | 100% | 90% | no |
| A0AAJ6YF89 | *Ceratosolen solmsi marchali* | 5090 | 6936.3 | 100% | 90% | no |
| V4AA46 | *Lottia gigantea* | 5045 | 6934.3 | 100% | 90% | no |
| A0AA39KH72 | *Microctonus hyperodae* | 5071 | 6927.9 | 99% | 90% | no |
| A0A182JT73 | *Anopheles christyi* | 4995 | 6927.8 | 98% | 90% | no |
| A0A8B7ZN74 | *Acanthaster planci* | 5283 | 6915.3 | 100% | 90% | no |
| D6WF72 | *Tribolium castaneum* | 5011 | 6914.2 | 98% | 90% | no |
| A0A4S2KV39 | *Temnothorax longispinosus* | 6896 | 6913.6 | 99% | 91% | no |
| A0A9J6CJ21 | *Polypedilum vanderplanki* | 5072 | 6912.5 | 100% | 90% | no |
| A0A182QKY2 | *Anopheles farauti* | 5013 | 6907.6 | 99% | 90% | no |
| W5JDV8 | *Anopheles darlingi* | 5004 | 6904.2 | 98% | 90% | no |
| A0A1A9VXJ3 | *Glossina austeni* | 5011 | 6904.1 | 98% | 90% | no |
| E9FTU9 | *Daphnia pulex* | 5119 | 6901.7 | 100% | 90% | no |
| A0AAD5L6S9 | *Daphnia sinensis* | 5164 | 6901.4 | 100% | 90% | no |
| A0ABD0YNJ7 | *Ranatra chinensis* | 5037 | 6899.6 | 99% | 90% | no |
| A0ABM1ZYA0 | *Aedes albopictus* | 5013 | 6892.0 | 98% | 90% | no |
| A0ABQ9ZR46 | *Daphnia magna* | 5172 | 6883.5 | 100% | 90% | no |
| A0A8S4S2J2 | *Pararge aegeria aegeria* | 5125 | 6874.7 | 96% | 90% | no |
| A0ACC2N2H0 | *Eretmocerus hayati* | 5019 | 6872.7 | 98% | 90% | no |
| A0AAW1DJ36 | *Rhynocoris fuscipes* | 4992 | 6864.3 | 98% | 90% | no |
| A0A0M8ZYB0 | *Melipona quadrifasciata* | 5082 | 6834.6 | 99% | 91% | no |
| A0AAV1ZCR3 | *Larinioides sclopetarius* | 5064 | 6827.8 | 100% | 89% | no |
| A0ACM7X729 | *Dermacentor variabilis* | 5213 | 6824.0 | 100% | 89% | no |
| A0A182UXB4 | *Anopheles merus* | 4955 | 6822.8 | 97% | 90% | no |
| A0A182TIU4 | *Anopheles melas* | 4934 | 6822.4 | 97% | 90% | no |
| A0ACM8F2K1 | *Ixodes scapularis* | 5156 | 6815.8 | 100% | 89% | no |
| A0AAV6UMA4 | *Oedothorax gibbosus* | 5046 | 6808.3 | 98% | 89% | no |
| A0A6H5I7I9 | *Trichogramma brassicae* | 5053 | 6804.3 | 99% | 89% | no |
| A0ACC1CG49 | *Dendrolimus kikuchii* | 4995 | 6796.2 | 98% | 90% | no |
| A0A2C9C3E8 | *Caenorhabditis elegans* | 5202 | 6795.4 | 100% | 90% | no |
| A0A2T7NWF8 | *Pomacea canaliculata* | 5035 | 6786.6 | 97% | 90% | no |
| A0A833VJX3 | *Frieseomelitta varia* | 4968 | 6777.4 | 97% | 91% | no |
| A0ACM8CH91 | *Saccoglossus kowalevskii* | 5221 | 6766.5 | 100% | 90% | no |
| A0ACF8R2W1 | *Solemya velum* | 5149 | 6765.7 | 92% | 90% | no |
| A0A8J2W1B1 | *Daphnia galeata* | 5083 | 6761.8 | 98% | 90% | no |
| A0A8J2HCS1 | *Cotesia congregata* | 5056 | 6744.2 | 98% | 90% | no |
| A0AAW0YJS6 | *Cherax quadricarinatus* | 5066 | 6741.4 | 100% | 90% | no |
| A0A164M377 | *Daphnia magna* | 5136 | 6723.9 | 99% | 90% | no |
| E3LPH3 | *Caenorhabditis remanei* | 5220 | 6713.5 | 100% | 90% | no |
| A0AAE9JLN6 | *Caenorhabditis briggsae* | 5277 | 6709.0 | 100% | 90% | no |
| A0ABP1QU76 | *Orchesella dallaii* | 5046 | 6706.7 | 100% | 90% | no |
| A0A8R1HP75 | *Caenorhabditis japonica* | 5076 | 6697.0 | 100% | 90% | no |
| A0A260YZI1 | *Caenorhabditis remanei* | 5269 | 6657.5 | 100% | 90% | no |
| A0A8S1E5U6 | *Caenorhabditis bovis* | 5259 | 6656.8 | 100% | 90% | no |
| A0AA39FIZ9 | *Microctonus aethiopoides* | 4909 | 6642.7 | 96% | 91% | no |
| A0A8B6HKR7 | *Mytilus galloprovincialis* | 5025 | 6627.6 | 94% | 90% | no |
| A0A1S3I436 | *Lingula anatina* | 4715 | 6618.9 | 88% | 90% | yes — circular |
| G0MY32 | *Caenorhabditis brenneri* | 5370 | 6618.6 | 100% | 90% | no |
| A0AAE1NS39 | *Petrolisthes manimaculis* | 5556 | 6607.0 | 98% | 90% | no |
| A0A8W8I0E5 | *Magallana gigas* | 5084 | 6598.6 | 93% | 90% | no |
| A0A261BHR6 | *Caenorhabditis latens* | 5209 | 6597.1 | 100% | 90% | no |
| A0ABM1SHV1 | *Limulus polyphemus* | 5186 | 6595.9 | 98% | 90% | no |
| A0A9Q0RJZ4 | *Blomia tropicalis* | 5152 | 6594.8 | 100% | 89% | no |
| A0A3S3QRZ4 | *Dinothrombium tinctorium* | 5145 | 6592.3 | 100% | 89% | no |
| A0A212F8V6 | *Danaus plexippus plexippus* | 4951 | 6579.8 | 96% | 93% | no |
| A0A9P1N8E4 | *Caenorhabditis angaria* | 5212 | 6573.9 | 100% | 90% | no |
| A0A443RD66 | *Dinothrombium tinctorium* | 5129 | 6564.0 | 99% | 89% | no |
| T1K6F3 | *Tetranychus urticae* | 5235 | 6558.6 | 100% | 90% | no |
| A0A7I4YT81 | *Haemonchus contortus* | 5055 | 6541.3 | 97% | 89% | no |
| A0ABQ8IR18 | *Dermatophagoides pteronyssinus* | 5324 | 6518.8 | 100% | 89% | no |
| A8XUM8 | *Caenorhabditis briggsae* | 6513 | 6518.0 | 100% | 89% | no |
| A0A7I4YSE1 | *Haemonchus contortus* | 5054 | 6517.3 | 97% | 89% | no |
| A0A154NXT6 | *Dufourea novaeangliae* | 4734 | 6507.5 | 93% | 91% | no |
| A0AA36DST8 | *Cylicocyclus nassatus* | 5216 | 6504.5 | 100% | 89% | no |
| A0A8T0EY51 | *Argiope bruennichi* | 5244 | 6500.7 | 95% | 89% | no |
| A0A7M7KDA0 | *Varroa destructor* | 5165 | 6498.1 | 100% | 90% | no |
| A0A8J5QYZ6 | *Cotesia typhae* | 4885 | 6486.9 | 95% | 93% | no |
| A0AAE1F6E0 | *Petrolisthes cinctipes* | 5490 | 6479.8 | 96% | 90% | no |
| A0AAJ7SIU4 | *Galendromus occidentalis* | 5151 | 6471.4 | 100% | 90% | no |
| A0ABP0GL45 | *Clavelina lepadiformis* | 5090 | 6450.9 | 100% | 89% | no |
| A0AA36D6A0 | *Mesorhabditis spiculigera* | 5069 | 6440.7 | 100% | 89% | no |
| A0A9C6U8W4 | *Frankliniella occidentalis* | 4775 | 6430.2 | 92% | 91% | no |
| A0A6L2Q4W7 | *Coptotermes formosanus* | 4766 | 6404.1 | 94% | 93% | no |
| A0A7I8XKB2 | *Bursaphelenchus xylophilus* | 5250 | 6395.9 | 100% | 89% | no |
| A0AAF3EW07 | *Mesorhabditis belari* | 5218 | 6395.4 | 100% | 89% | no |
| H2XSF3 | *Ciona intestinalis* | 5075 | 6388.8 | 100% | 89% | no |
| A0AAF3EW03 | *Mesorhabditis belari* | 5247 | 6378.6 | 100% | 89% | no |
| A0A811LFJ3 | *Bursaphelenchus okinawaensis* | 5318 | 6369.3 | 100% | 89% | no |
| A0A090MUA7 | *Strongyloides ratti* | 5303 | 6351.6 | 100% | 89% | no |
| A0AAF5DMV8 | *Strongyloides stercoralis* | 5295 | 6349.7 | 100% | 89% | no |
| A0AAF5DJD3 | *Strongyloides stercoralis* | 5273 | 6349.6 | 100% | 89% | no |
| A0AAN5DCL9 | *Pristionchus mayeri* | 5244 | 6349.3 | 98% | 89% | no |
| A0AAV7I734 | *Cotesia glomerata* | 4857 | 6341.6 | 94% | 93% | no |
| A0A0N4ZPH7 | *Parastrongyloides trichosuri* | 5293 | 6338.4 | 100% | 89% | no |
| A0A1P6C3Y0 | *Brugia malayi* | 5034 | 6332.0 | 100% | 90% | no |
| A0A0K0FNZ1 | *Strongyloides venezuelensis* | 5313 | 6331.0 | 100% | 89% | no |
| A0A8J9ZXK5 | *Branchiostoma lanceolatum* | 4975 | 6323.2 | 100% | 89% | no |
| A0A9J7MUA9 | *Branchiostoma floridae* | 4992 | 6310.8 | 100% | 89% | no |
| A0AA39HTB0 | *Steinernema hermaphroditum* | 5288 | 6310.0 | 100% | 89% | no |
| A0A158Q797 | *Elaeophora elaphi* | 5055 | 6307.0 | 100% | 89% | no |
| A0A836JHR0 | *Acromyrmex insinuator* | 4697 | 6300.8 | 91% | 91% | no |
| A0A158Q013 | *Brugia malayi* | 5033 | 6297.5 | 100% | 90% | no |
| A0A6P4XJ96 | *Branchiostoma belcheri* | 4960 | 6293.4 | 99% | 89% | no |
| A0A158PQQ2 | *Brugia pahangi* | 5139 | 6278.9 | 100% | 89% | no |
| A0A834VC47 | *Sarcoptes scabiei* | 5042 | 6271.4 | 98% | 91% | no |
| A0ABD2II26 | *Heterodera schachtii* | 5380 | 6266.0 | 100% | 89% | no |
| A0A8J2LX31 | *Cercopithifilaria johnstoni* | 5013 | 6252.8 | 100% | 90% | no |
| A0A0N5BIM7 | *Strongyloides papillosus* | 5177 | 6250.6 | 98% | 89% | no |
| A0A3P7G7P7 | *Wuchereria bancrofti* | 5126 | 6249.5 | 97% | 89% | no |
| A0A7M7T5D0 | *Strongylocentrotus purpuratus* | 4609 | 6230.3 | 86% | 93% | no |
| A0A813MD17 | *Brachionus calyciflorus* | 5152 | 6206.0 | 99% | 89% | no |
| A0A914H380 | *Globodera rostochiensis* | 5463 | 6142.7 | 100% | 89% | no |
| A0AAD4R815 | *Ditylenchus destructor* | 5178 | 6139.7 | 99% | 89% | no |
| A0ABR3H4X6 | *Loxostege sticticalis* | 4357 | 6133.8 | 85% | 93% | no |
| A0A914H250 | *Globodera rostochiensis* | 5247 | 6127.5 | 100% | 89% | no |
| A0A1I7VT37 | *Loa loa* | 4961 | 6118.9 | 98% | 89% | no |
| A0A914H5M0 | *Globodera rostochiensis* | 5428 | 6113.4 | 100% | 89% | no |
| A0A2A6CRA2 | *Pristionchus pacificus* | 5221 | 6112.8 | 99% | 89% | no |
| A0A8R1TZT6 | *Onchocerca volvulus* | 5043 | 6109.4 | 99% | 90% | no |
| B3RK74 | *Trichoplax adhaerens* | 4949 | 6104.0 | 96% | 89% | yes — circular |
| A0A813T754 | *Rotaria sordida* | 5245 | 6080.3 | 100% | 89% | no |
| A0A813R8P9 | *Rotaria sordida* | 5245 | 6078.7 | 100% | 89% | no |
| A0A0N5AR96 | *Syphacia muris* | 4873 | 6075.7 | 95% | 90% | no |
| A0A4Y2G0P6 | *Araneus ventricosus* | 4723 | 6070.4 | 92% | 93% | no |
| A0A4Y2G3N8 | *Araneus ventricosus* | 4701 | 6069.3 | 92% | 93% | no |
| H2Y884 | *Ciona savignyi* | 5011 | 6059.7 | 99% | 88% | no |
| A0A814WFU3 | *Adineta steineri* | 5242 | 6057.6 | 100% | 89% | no |
| A0A813NCM0 | *Adineta steineri* | 5242 | 6055.0 | 100% | 89% | no |
| A0A4Y2G666 | *Araneus ventricosus* | 4701 | 6053.5 | 92% | 93% | no |
| A0A085LRI9 | *Trichuris suis* | 5145 | 5967.3 | 100% | 88% | no |
| A0A814GWS8 | *Adineta steineri* | 5143 | 5963.3 | 100% | 89% | no |
| A0A813ZMY6 | *Adineta steineri* | 5143 | 5956.4 | 100% | 89% | no |
| A0A5S6QXR5 | *Trichuris muris* | 5075 | 5951.1 | 100% | 88% | no |
| A0A0V0SMI2 | *Trichinella nelsoni* | 5112 | 5919.9 | 100% | 88% | no |
| A0A0V1DD79 | *Trichinella britovi* | 5150 | 5918.1 | 100% | 88% | no |
| A0A0V1AFC9 | *Trichinella patagoniensis* | 5150 | 5917.4 | 100% | 88% | no |
| A0A0V1LTI1 | *Trichinella nativa* | 5150 | 5916.5 | 100% | 88% | no |
| A0A553P1H5 | *Tigriopus californicus* | 4758 | 5915.5 | 94% | 89% | no |
| A0ABR1DUS7 | *Necator americanus* | 4901 | 5915.5 | 90% | 93% | no |
| A0A0V0UGK5 | *Trichinella murrelli* | 5150 | 5915.3 | 100% | 88% | no |
| A0A183UDS7 | *Toxocara canis* | 4894 | 5905.8 | 93% | 90% | no |
| A0A077YYK0 | *Trichuris trichiura* | 5172 | 5887.8 | 100% | 89% | no |
| W2TAC6 | *Necator americanus* | 4907 | 5876.8 | 92% | 90% | no |
| A0A0V1FVZ4 | *Trichinella pseudospiralis* | 5098 | 5829.9 | 100% | 88% | no |
| A0A915B3X3 | *Parascaris univalens* | 4739 | 5823.0 | 92% | 93% | no |
| A0A813QCU3 | *Didymodactylos carnosus* | 5233 | 5822.6 | 96% | 89% | no |
| A0A0V1HV04 | *Trichinella zimbabwensis* | 5064 | 5809.1 | 99% | 88% | no |
| A0A915B6F5 | *Parascaris univalens* | 4738 | 5795.7 | 92% | 93% | no |
| A0ABN7S0R9 | *Oikopleura dioica* | 5011 | 5789.2 | 100% | 88% | no |
| E4XPI8 | *Oikopleura dioica* | 4968 | 5770.0 | 100% | 88% | no |
| A0A158P902 | *Angiostrongylus cantonensis* | 4798 | 5768.5 | 94% | 90% | no |
| A0AAV5V408 | *Pristionchus fissidentatus* | 4684 | 5761.6 | 90% | 90% | no |
| A0A0V1JF63 | *Trichinella pseudospiralis* | 5151 | 5701.5 | 100% | 89% | no |
| A0AA85FAW1 | *Schistosoma rodhaini* | 5653 | 5696.1 | 100% | 88% | no |
| A0A5K4EV29 | *Schistosoma mansoni* | 5663 | 5681.3 | 100% | 88% | no |
| A0A3Q0KS03 | *Schistosoma mansoni* | 5669 | 5680.1 | 100% | 88% | no |
| A0AA85FB23 | *Schistosoma rodhaini* | 5663 | 5677.2 | 100% | 88% | no |
| A0AA85FAP7 | *Schistosoma rodhaini* | 5669 | 5676.1 | 100% | 88% | no |
| A0A4Z2DNF4 | *Schistosoma japonicum* | 5688 | 5664.3 | 100% | 87% | no |
| A0AAE1ZFP2 | *Schistosoma mekongi* | 5688 | 5662.9 | 100% | 87% | no |
| A0A0V1C2E6 | *Trichinella spiralis* | 5115 | 5638.4 | 100% | 88% | no |
| A0A4S2LT93 | *Opisthorchis felineus* | 5609 | 5629.3 | 97% | 87% | no |
| A0ACE1DKF8 | *Calicophoron daubneyi* | 5674 | 5614.6 | 100% | 88% | no |
| A0ABQ7S901 | *Fragariocoptes setiger* | 6761 | 5613.9 | 100% | 88% | no |
| A0A922LUH1 | *Schistosoma haematobium* | 5605 | 5611.2 | 99% | 88% | no |
| A0A6A4WJ60 | *Amphibalanus amphitrite* | 4233 | 5598.6 | 82% | 93% | no |
| A0A813NLG5 | *Didymodactylos carnosus* | 5157 | 5589.6 | 97% | 89% | no |
| A0ABR4Q6E9 | *Taenia crassiceps* | 5781 | 5552.9 | 100% | 87% | no |
| A0AAV4NM55 | *Caerostris darwini* | 4289 | 5535.8 | 84% | 90% | no |
| A0AAV4NN14 | *Caerostris darwini* | 4267 | 5534.9 | 84% | 90% | no |
| A0A4E0RMD6 | *Fasciola hepatica* | 5599 | 5522.4 | 99% | 88% | no |
| A0AAV4NLR5 | *Caerostris darwini* | 4311 | 5519.3 | 84% | 90% | no |
| A0AAV4NN99 | *Caerostris darwini* | 4267 | 5518.5 | 84% | 90% | no |
| A0A1I8H4L8 | *Macrostomum lignano* | 4842 | 5517.4 | 92% | 90% | no |
| A0A504YGJ9 | *Fasciola gigantica* | 5566 | 5475.9 | 98% | 88% | no |
| A0A8J4THQ5 | *Paragonimus heterotremus* | 5627 | 5456.0 | 99% | 87% | no |
| A0A8T0D2K4 | *Paragonimus westermani* | 5548 | 5446.2 | 98% | 88% | no |
| A0AA85FAC6 | *Schistosoma rodhaini* | 5572 | 5412.0 | 98% | 90% | no |
| A0A8S9Z2A3 | *Paragonimus skrjabini miyazakii* | 5491 | 5398.1 | 97% | 87% | no |
| A0A8E0VKM4 | *Fasciolopsis buskii* | 5541 | 5390.6 | 98% | 88% | no |
| A0A0R3UJQ6 | *Mesocestoides corti* | 5743 | 5388.8 | 99% | 87% | no |
| A0A6J3KG89 | *Bombus vosnesenskii* | 4060 | 5343.1 | 78% | 94% | no |
| A0A8T1LXK5 | *Clonorchis sinensis* | 5698 | 5327.5 | 96% | 89% | no |
| A0A158PDE3 | *Angiostrongylus costaricensis* | 4328 | 5315.3 | 83% | 93% | no |
| A0AAD8B2N2 | *Biomphalaria pfeifferi* | 4080 | 5280.4 | 78% | 90% | no |
| A0A915EYX6 | *Echinococcus canadensis* | 6373 | 5260.5 | 97% | 87% | no |
| A0ACH3JT27 | *Kerria lacca* | 3977 | 5248.3 | 77% | 94% | no |
| A0A5J4NXI2 | *Paragonimus westermani* | 5606 | 5241.9 | 98% | 88% | no |
| A0ACE7F6K6 | *Sergentomyia squamirostris* | 3985 | 5207.8 | 76% | 90% | no |
| B4HSA0 | *Drosophila sechellia* | 3814 | 5200.1 | 71% | 88% | no |
| A0A811V8B2 | *Ceratitis capitata* | 3902 | 5179.3 | 75% | 95% | no |
| B4QG94 | *Drosophila simulans* | 3487 | 5137.9 | 70% | 88% | no |
| A0A5E4QI94 | *Leptidea sinapis* | 3842 | 5122.0 | 75% | 95% | no |
| A0A8J6H6D8 | *Tenebrio molitor* | 4154 | 5095.2 | 81% | 96% | no |
| A0A2A2JM86 | *Diploscapter pachys* | 4249 | 5079.2 | 78% | 94% | no |
| A0A267GJE7 | *Macrostomum lignano* | 4427 | 5048.0 | 86% | 92% | no |
| A0A8J5JM50 | *Homarus americanus* | 4088 | 5043.8 | 78% | 94% | no |
| A0A1S3I400 | *Lingula anatina* | 3518 | 5032.3 | 69% | 90% | no |
| A0A7E6F791 | *Octopus sinensis* | 4022 | 5018.0 | 75% | 95% | no |
| A0A0B2V2Q3 | *Toxocara canis* | 4668 | 5009.4 | 75% | 97% | no |
| A0ABD2CS54 | *Vespula maculifrons* | 3780 | 5006.9 | 73% | 94% | no |
| A0ABD2ABN5 | *Vespula squamosa* | 3780 | 5003.0 | 73% | 94% | no |
| A0A5S6PZZ0 | *Trichuris muris* | 4429 | 5001.7 | 87% | 90% | no |
| A0A9Q1CGM9 | *Holothuria leucospilota* | 4119 | 4998.5 | 77% | 94% | no |
| A0ABM1TA57 | *Limulus polyphemus* | 4082 | 4979.0 | 76% | 89% | no |
| A0A9P0BYQ2 | *Chrysodeixis includens* | 3766 | 4960.8 | 73% | 94% | no |
| A0AA38MLI7 | *Zophobas morio* | 3844 | 4959.8 | 73% | 95% | no |
| A0A5E4LY32 | *Cinara cedri* | 3669 | 4948.0 | 72% | 89% | no |
| T1J5U5 | *Strigamia maritima* | 4116 | 4897.5 | 76% | 91% | no |
| A0ABR1B0T0 | *Polyplax serrata* | 3872 | 4871.4 | 72% | 94% | no |
| A0AA85F8T4 | *Schistosoma rodhaini* | 5079 | 4832.6 | 87% | 91% | no |
| A0AA85F8G5 | *Schistosoma rodhaini* | 5082 | 4832.2 | 87% | 91% | no |
| A0AAD9N3M2 | *Paralvinella palmiformis* | 3804 | 4830.2 | 70% | 95% | no |
| A0ACC0KBC1 | *Choristoneura fumiferana* | 3625 | 4829.8 | 70% | 94% | no |
| A0AAW1M9L9 | *Popillia japonica* | 3675 | 4824.1 | 70% | 94% | no |
| A0A075A5M8 | *Opisthorchis viverrini* | 5727 | 4821.9 | 90% | 88% | no |
| A0AA85F8J4 | *Schistosoma rodhaini* | 4962 | 4821.6 | 85% | 91% | no |
| A0AA85JZA6 | *Trichobilharzia regenti* | 4994 | 4801.8 | 85% | 91% | no |
| A0A2G5TLD0 | *Caenorhabditis nigoni* | 3896 | 4793.6 | 72% | 95% | no |
| A0ABD0L2J9 | *Batillaria attramentaria* | 3798 | 4771.1 | 71% | 94% | no |
| A0AAV7K8I1 | *Oopsacas minuta* | 4959 | 4744.3 | 95% | 90% | yes — circular |
| A0A834U2V9 | *Vespula germanica* | 3737 | 4711.7 | 70% | 94% | no |
| A0A0D2WXI6 | *Capsaspora owczarzaki (strain ATCC 30864)* | 6625 | 4662.5 | 96% | 86% | yes — circular |
| A0AA85F9M4 | *Schistosoma rodhaini* | 4847 | 4636.6 | 83% | 91% | no |
| A0A7R8WZ42 | *Darwinula stevensoni* | 4163 | 4624.4 | 74% | 89% | no |
| A0A8S1I045 | *Caenorhabditis auriculariae* | 3871 | 4606.5 | 72% | 95% | no |
| A0A016U0Q9 | *Ancylostoma ceylanicum* | 3941 | 4602.1 | 72% | 94% | no |
| A0A564YS66 | *Hymenolepis diminuta* | 4964 | 4572.0 | 85% | 89% | no |
| A0A6G0Z6V5 | *Aphis craccivora* | 3427 | 4494.5 | 65% | 94% | no |
| A0A1B0A7F3 | *Glossina pallidipes* | 3171 | 4476.7 | 63% | 89% | no |
| A0A1P6C413 | *Brugia malayi* | 3732 | 4462.1 | 73% | 95% | no |
| A0A1B0BLM4 | *Glossina palpalis gambiensis* | 3168 | 4457.5 | 63% | 89% | no |
| A0A7E4VY38 | *Panagrellus redivivus* | 4039 | 4457.5 | 72% | 95% | no |
| A0A3P7SYB0 | *Dracunculus medinensis* | 4014 | 4444.7 | 71% | 89% | no |
| A0A815JFZ8 | *Adineta ricciae* | 4154 | 4439.6 | 80% | 95% | no |
| A0A194QMY4 | *Papilio xuthus* | 3275 | 4434.6 | 63% | 94% | no |
| A0A815JQP0 | *Adineta ricciae* | 4154 | 4433.1 | 80% | 95% | no |
| A0A3P6UQP9 | *Litomosoides sigmodontis* | 3662 | 4431.2 | 72% | 94% | no |
| A0A8S1BAV5 | *Arctia plantaginis* | 3315 | 4430.1 | 63% | 94% | no |
| A0AAU9UM41 | *Euphydryas editha* | 3253 | 4398.0 | 62% | 94% | no |
| A0A821VFV6 | *Pieris macdunnoughi* | 3242 | 4390.6 | 63% | 94% | no |
| A0A482X586 | *Laodelphax striatellus* | 3308 | 4383.5 | 63% | 94% | no |
| A0A8J2QY98 | *Danaus chrysippus* | 3237 | 4374.1 | 63% | 94% | no |
| A0A7R9BQ83 | *Notodromas monacha* | 3631 | 4368.7 | 67% | 89% | no |
| A0AAV2Q4J1 | *Meganyctiphanes norvegica* | 3455 | 4362.4 | 68% | 95% | no |
| A0A0R3WD36 | *Taenia asiatica* | 5267 | 4361.5 | 82% | 92% | no |
| A0ABN8ILP6 | *Iphiclides podalirius* | 3263 | 4360.9 | 63% | 94% | no |
| A0A183BIF4 | *Globodera pallida* | 3683 | 4348.0 | 69% | 93% | no |
| A0A8B8DQN4 | *Crassostrea virginica* | 3436 | 4339.2 | 63% | 94% | no |
| A0AAN8EZV3 | *Trichostrongylus colubriformis* | 3596 | 4335.1 | 68% | 88% | no |
| E2AGR1 | *Camponotus floridanus* | 3132 | 4326.4 | 62% | 88% | no |
| A0A3P7G6U0 | *Hydatigera taeniaeformis* | 5152 | 4299.4 | 84% | 89% | no |
| A0A068Y0W3 | *Echinococcus multilocularis* | 4760 | 4293.4 | 80% | 88% | no |
| A0AA85F8H2 | *Schistosoma rodhaini* | 4713 | 4268.8 | 79% | 92% | no |
| A0A4C1VKE3 | *Eumeta variegata* | 3532 | 4254.6 | 66% | 95% | no |
| A0A0R3QJC6 | *Brugia timori* | 3602 | 4194.9 | 67% | 95% | no |
| A0A922HHF9 | *Dermatophagoides farinae* | 3322 | 4132.9 | 63% | 93% | no |
| A0A194RNL9 | *Papilio machaon* | 2892 | 4110.3 | 57% | 88% | no |
| A0A8K0KDG3 | *Ladona fulva* | 3379 | 4109.2 | 63% | 95% | no |
| A0A0N5CZR3 | *Thelazia callipaeda* | 3363 | 4094.6 | 66% | 94% | no |
| A0AAV4CU45 | *Plakobranchus ocellatus* | 3331 | 4064.5 | 58% | 93% | no |
| A0A0V1N0Y7 | *Trichinella papuae* | 3473 | 4060.0 | 71% | 88% | no |
| A0AAN8ZTZ2 | *Halocaridina rubra* | 3266 | 4041.8 | 63% | 95% | no |
| A0A226ESW6 | *Folsomia candida* | 3178 | 3994.9 | 60% | 88% | no |
| A0A8X6LSL1 | *Trichonephila clavata* | 3380 | 3987.7 | 66% | 93% | no |
| A0A8X6LRW5 | *Trichonephila clavata* | 3395 | 3982.1 | 66% | 93% | no |
| A0ABM1DXK4 | *Priapulus caudatus* | 3274 | 3976.1 | 63% | 95% | no |
| A0A7R9PVA5 | *Medioppia subpectinata* | 3592 | 3936.2 | 67% | 90% | no |
| A0AAW0SKF4 | *Scylla paramamosain* | 2969 | 3853.1 | 58% | 93% | no |
| A0A1I7T9I4 | *Caenorhabditis tropicalis* | 2979 | 3837.4 | 56% | 87% | no |
| A0A1I7T9I5 | *Caenorhabditis tropicalis* | 2965 | 3822.4 | 56% | 87% | no |
| A0A0D6LLP6 | *Ancylostoma ceylanicum* | 3360 | 3791.3 | 61% | 95% | no |
| A0AA85F929 | *Schistosoma rodhaini* | 4169 | 3725.2 | 70% | 93% | no |
| A0A915Q2X8 | *Setaria digitata* | 3050 | 3712.9 | 58% | 87% | no |
| A0A815M856 | *Adineta ricciae* | 3298 | 3648.9 | 62% | 93% | no |
| G7YIT2 | *Clonorchis sinensis* | 4087 | 3633.7 | 69% | 93% | no |
| A0AA88HIW1 | *Artemia franciscana* | 2646 | 3589.4 | 51% | 88% | no |
| A0AAV8W0J6 | *Exocentrus adspersus* | 2680 | 3588.7 | 51% | 92% | no |
| A0A1S4EE60 | *Diaphorina citri* | 3137 | 3586.0 | 52% | 93% | no |
| A0A914BA11 | *Patiria miniata* | 2795 | 3576.9 | 52% | 92% | no |
| A0A1W0WNM8 | *Hypsibius exemplaris* | 3045 | 3548.7 | 57% | 94% | no |
| A0A1V9WZ08 | *Tropilaelaps mercedesae* | 2990 | 3492.1 | 57% | 95% | no |
| A0ABY7FBH2 | *Mya arenaria* | 2884 | 3447.3 | 55% | 94% | no |
| A0A1D1VU71 | *Ramazzottius varieornatus* | 2875 | 3421.8 | 57% | 85% | no |
| A0A267FNK0 | *Macrostomum lignano* | 2694 | 3366.7 | 51% | 87% | no |
| A0A1D2N817 | *Orchesella cincta* | 2712 | 3363.5 | 52% | 90% | no |
| A0A368H5X8 | *Ancylostoma caninum* | 2871 | 3363.5 | 50% | 93% | no |
| A0AAW1UQP2 | *Henosepilachna vigintioctopunctata* | 2763 | 3236.3 | 53% | 98% | no |
| A0A267DX28 | *Macrostomum lignano* | 2874 | 3213.1 | 57% | 93% | no |
| A0ABD2QM95 | *Cichlidogyrus casuarinus* | 3114 | 3153.4 | 56% | 88% | no |
| A0A3P6SJE8 | *Onchocerca ochengi* | 2445 | 3078.9 | 51% | 86% | no |
| A0A0L7LMR0 | *Operophtera brumata* | 3507 | 3069.6 | 59% | 93% | no |
| A0ACB0ZJ18 | *Meloidogyne enterolobii* | 3021 | 3032.4 | 51% | 97% | no |
| W6U4N9 | *Echinococcus granulosus* | 3636 | 3022.7 | 57% | 92% | no |
| A0A0N4V1U6 | *Enterobius vermicularis* | 2974 | 2922.1 | 53% | 86% | no |
| A0A3P8ALE8 | *Heligmosomoides polygyrus* | 2653 | 2740.3 | 53% | 98% | no |
| A0A820LXG1 | *Rotaria socialis* | 2627 | 2681.9 | 50% | 96% | no |
| A0AAV7JZ53 | *Oopsacas minuta* | 4658 | 2612.6 | 89% | 81% | no |
| A0AAN0JF00 | *Amphimedon queenslandica* | 4753 | 2475.8 | 86% | 82% | no |
| A0ABR3KI48 | *Trichinella spiralis* | 2447 | 2443.0 | 52% | 93% | no |
| A0AA35WZQ4 | *Geodia barretti* | 3421 | 2005.5 | 63% | 91% | no |
| F2UC37 | *Salpingoeca rosetta* | 5340 | 1387.2 | 59% | 65% | no |


A record that is itself one of `ryr.hmm`'s seeds scores well against a profile built partly from it, so its score is not independent evidence and the column says so. The rest are.


The contested targets — kept and reported, not resolved:

| accession | species | itpr bits | ryr bits | margin |
|---|---|---|---|---|
| A0AAN0JU79 | *Amphimedon queenslandica* | 44.4 | 45.1 | 1.6% |


## Land plants: a genome fact or a database fact?

S2 enumerated the InterPro-seeded space and found **Streptophyta — the land-plant lineage — with 0 ITPR calls from 15 records in 13 taxa**. That is a claim about a space a protein enters only by already carrying a family Pfam annotation. Sweeping the proteomes themselves asks the question without that filter.

**CONFIRMED** — 0/384 (0%) Streptophyta reference proteomes carry an ITPR call: searching the proteomes themselves finds no more than the seeded space did.

| clade | proteomes swept | with an ITPR call |
|---|---|---|
| Chlorophyta | 48 | 15/48 (31%) |
| Streptophyta | 384 | 0/384 (0%) |


Viridiplantae overall: 15/432 (3%) proteomes.


Every Viridiplantae record called ITPR by either instrument was chased individually (64 records):

| verdict | records | what it means |
|---|---|---|
| fragment | 42 | real but incomplete gene model |
| real_gene | 22 | survives every test |


The surviving records are not evenly spread: **Cymbomonas tetramitiformis** alone contributes 12 of them, against a median of 1 per species. Whether that is real copy-number expansion or a duplicated assembly is a question for an assembly search, not a proteome one *(pending: S23)*.


The 22 that survive, by lineage:

| accession | species | phylum | length | nearest outside its kingdom |
|---|---|---|---|---|
| A0A150GRF6 | *Gonium pectorale (Green alga)* | Chlorophyta | 2580 | 23.2 % (Vertebrata) |
| A0A2J8AAQ9 | *Tetrabaena socialis* | Chlorophyta | 2680 | 22.6 % (Vertebrata) |
| A0A2K3CTW4 | *Chlamydomonas reinhardtii (Chlamydomonas smithii)* | Chlorophyta | 3210 | 21.6 % (Metazoa (non-vertebrate)) |
| A0A835TB45 | *Chlamydomonas incerta* | Chlorophyta | 3182 | 20.9 % (Vertebrata) |
| A0A835Y1X9 | *Edaphochlamys debaryana* | Chlorophyta | 3268 | 23.2 % (Vertebrata) |
| A0A836B9Y3 | *Chlamydomonas schloesseri* | Chlorophyta | 3219 | 21.9 % (Vertebrata) |
| A0A8S1J4S2 | *Ostreobium quekettii* | Chlorophyta | 2446 | 20.5 % (Fungi) |
| A0A9W6C1N2 | *Pleodorina starrii* | Chlorophyta | 3372 | 22.6 % (Vertebrata) |
| A0AAE0BBD3 | *Cymbomonas tetramitiformis* | Chlorophyta | 2899 | 20.6 % (Vertebrata) |
| A0AAE0BCU2 | *Cymbomonas tetramitiformis* | Chlorophyta | 2953 | 19.9 % (Vertebrata) |
| A0AAE0ESH2 | *Cymbomonas tetramitiformis* | Chlorophyta | 3164 | 28.4 % (SAR) |
| A0AAE0ET13 | *Cymbomonas tetramitiformis* | Chlorophyta | 2660 | 20.4 % (Metazoa (non-vertebrate)) |
| A0AAE0ETT8 | *Cymbomonas tetramitiformis* | Chlorophyta | 3163 | 28.4 % (SAR) |
| A0AAE0EZ22 | *Cymbomonas tetramitiformis* | Chlorophyta | 2761 | 21.6 % (Metazoa (non-vertebrate)) |
| A0AAE0G204 | *Cymbomonas tetramitiformis* | Chlorophyta | 2900 | 39.9 % (SAR) |
| A0AAE0G465 | *Cymbomonas tetramitiformis* | Chlorophyta | 2922 | 22.4 % (Vertebrata) |
| A0AAE0GMP6 | *Cymbomonas tetramitiformis* | Chlorophyta | 2998 | 38.5 % (SAR) |
| A0AAE0GPV6 | *Cymbomonas tetramitiformis* | Chlorophyta | 3355 | 24.3 % (SAR) |
| A0AAE0GZ30 | *Cymbomonas tetramitiformis* | Chlorophyta | 2462 | 21.5 % (Metazoa (non-vertebrate)) |
| A0AAE0LDZ1 | *Cymbomonas tetramitiformis* | Chlorophyta | 3375 | 24.3 % (SAR) |


## Dikarya: the yeasts and moulds

S2 found **Dikarya contributing no records at all** to the seeded space — not zero calls from some records, zero records. Every fungal call sat in an early-diverging phylum (Mucoromycota 15, Chytridiomycota 6, Basidiobolomycota 3, Entomophthoromycota 1).

**CONFIRMED** — 0/1353 (0%) Ascomycota + Basidiomycota reference proteomes carry an ITPR call: searching the proteomes themselves finds no more than the seeded space did.

| clade | proteomes swept | with an ITPR call |
|---|---|---|
| Mucoromycota | 34 | 18/34 (53%) |
| Chytridiomycota | 16 | 6/16 (38%) |
| Basidiobolomycota | 2 | 2/2 (100%) |
| Zoopagomycota | 3 | 1/3 (33%) |
| Entomophthoromycota | 2 | 1/2 (50%) |
| Ascomycota | 1034 | 0/1034 (0%) |
| Basidiomycota | 319 | 0/319 (0%) |
| Kickxellomycota | 35 | 0/35 (0%) |
| Microsporidia | 29 | 0/29 (0%) |
| Glomeromycota | 27 | 0/27 (0%) |
| Mortierellomycota | 18 | 0/18 (0%) |
| Rozellomycota | 2 | 0/2 (0%) |
| Blastocladiomycota | 2 | 0/2 (0%) |
| Neocallimastigomycota | 2 | 0/2 (0%) |
| Monoblepharomycota | 1 | 0/1 (0%) |
| Olpidiomycota | 1 | 0/1 (0%) |


Fungi overall: 28/1527 (2%) proteomes.


Every Fungi record called ITPR by either instrument was chased individually (35 records):

| verdict | records | what it means |
|---|---|---|
| real_gene | 25 | survives every test |
| fragment | 10 | real but incomplete gene model |


The 25 that survive, by lineage:

| accession | species | phylum | length | nearest outside its kingdom |
|---|---|---|---|---|
| A0A1Y1XHY3 | *Basidiobolus meristosporus CBS 931.73* | Basidiobolomycota | 2495 | 24.2 % (Metazoa (non-vertebrate)) |
| A0A1Y1YQ92 | *Basidiobolus meristosporus CBS 931.73* | Basidiobolomycota | 2278 | 29.8 % (Amoebozoa) |
| A0ABR2W6A1 | *Basidiobolus ranarum* | Basidiobolomycota | 2502 | 24.6 % (Metazoa (non-vertebrate)) |
| A0A1Y2B255 | *Rhizoclosmatium globosum* | Chytridiomycota | 3246 | 33.2 % (Amoebozoa) |
| A0A507FK04 | *Chytriomyces confervae* | Chytridiomycota | 3329 | 33.5 % (Amoebozoa) |
| A0AAD5U9K5 | *Clydaea vesicula* | Chytridiomycota | 2759 | 30.5 % (Amoebozoa) |
| A0AAD5XWK9 | *Clydaea vesicula* | Chytridiomycota | 2327 | 20.5 % (Metazoa (non-vertebrate)) |
| A0ABR4NG49 | *Polyrhizophydium stewartii* | Chytridiomycota | 3104 | 23.6 % (Metazoa (non-vertebrate)) |
| A0ABR4NID5 | *Polyrhizophydium stewartii* | Chytridiomycota | 4446 | 30.5 % (Amoebozoa) |
| A0ACC2TPX2 | *Entomophthora muscae* | Entomophthoromycota | 2567 | 23.5 % (Metazoa (non-vertebrate)) |
| A0A0C9M4R7 | *Mucor ambiguus* | Mucoromycota | 2599 | 23.2 % (Metazoa (non-vertebrate)) |
| A0A162V9Q6 | *Phycomyces blakesleeanus (strain ATCC 8743b / DSM 1359 / FGSC 10004 / NBRC 33097 / NRRL 1555)* | Mucoromycota | 2551 | 23.5 % (Metazoa (non-vertebrate)) |
| A0A168KKP7 | *Mucor lusitanicus CBS 277.49* | Mucoromycota | 2620 | 33.1 % (Amoebozoa) |
| A0A168RJ50 | *Absidia glauca (Pin mould)* | Mucoromycota | 2540 | 22.7 % (Metazoa (non-vertebrate)) |
| A0A1C7NSC1 | *Choanephora cucurbitarum* | Mucoromycota | 2547 | 22.3 % (Metazoa (non-vertebrate)) |
| A0A1X2G6S3 | *Hesseltinella vesiculosa* | Mucoromycota | 2642 | 23.9 % (Metazoa (non-vertebrate)) |
| A0A1X2IXV9 | *Absidia repens* | Mucoromycota | 2651 | 24.2 % (Metazoa (non-vertebrate)) |
| A0A433QY11 | *Jimgerdemannia flammicorona* | Mucoromycota | 2599 | 23.6 % (Metazoa (non-vertebrate)) |
| A0A8H7BU48 | *Apophysomyces ossiformis* | Mucoromycota | 2391 | 22.4 % (Metazoa (non-vertebrate)) |
| A0A8H7VCW8 | *Mucor plumbeus* | Mucoromycota | 2604 | 22.5 % (Metazoa (non-vertebrate)) |


## The negative claims, at a stated sensitivity

Every claim was made at **E ≤ 1e-5** with the two full-length family profiles and re-made at **E ≤ 10** with those two *plus* the family's four Pfam domain models — a 213-state IP₃-core model can reach something a 2,684-state channel model cannot.


Reported at E ≤ 10 / of those, spanning at least half the model at E ≤ 1e-5:

| lineage | PF01365 | PF02815 | PF08454 | PF08709 | itpr | ryr |
|---|---|---|---|---|---|---|
| land plants | 357 / **0** | 849 / **633** | 137 / **1** | 448 / **0** | 744 / **0** | 32,423 / **0** |
| Chlorophyta (control) | 36 / **30** | 81 / **37** | 34 / **29** | 44 / **26** | 146 / **23** | 1,038 / **0** |
| Dikarya | 88 / **0** | 4,527 / **4,376** | 82 / **0** | 25 / **0** | 435 / **0** | 6,097 / **0** |
| Mucoromycota (control) | 12 / **3** | 351 / **329** | 18 / **16** | 100 / **16** | 120 / **16** | 709 / **0** |


**The controls fire and the claims do not.** PF08709 is the IP₃-binding core, the signature that names the family:

- **land plants**: 0 substantial PF08709 match(es), against 26 in *Chlorophyta*, the nearest lineage in the same kingdom where the family is present.

- **Dikarya**: 0 substantial PF08709 match(es), against 16 in *Mucoromycota*, the nearest lineage in the same kingdom where the family is present.


And the search is demonstrably sensitive in those very genomes: PF02815 (MIR, which the family shares with POMT1/2 and every eukaryote therefore carries) returns 633 substantial matches in land plants, 4,376 substantial matches in Dikarya. The promiscuous domain finds thousands where the family-defining one finds none, which is what an absence looks like when the instrument is working rather than blind.


Across every relaxed group, **51 target(s) span at least half a full-length family profile** at the relaxed threshold.

| group | profile | accession | species | length | E-value | model coverage |
|---|---|---|---|---|---|---|
| fungi | itpr | A0ABR2W6A1 | *Basidiobolus ranarum* | 2502 | 0.0 | 97% |
| fungi | itpr | A0A1Y1XHY3 | *Basidiobolus meristosporus CBS 931.73* | 2495 | 0.0 | 94% |
| fungi | itpr | A0A1Y1YQ92 | *Basidiobolus meristosporus CBS 931.73* | 2278 | 0.0 | 85% |
| fungi | itpr | A0A1Y2B255 | *Rhizoclosmatium globosum* | 3246 | 0.0 | 97% |
| fungi | itpr | A0A433QY11 | *Jimgerdemannia flammicorona* | 2599 | 0.0 | 90% |
| fungi | itpr | A0A507FK04 | *Chytriomyces confervae* | 3329 | 0.0 | 90% |
| fungi | itpr | A0ACC2TPX2 | *Entomophthora muscae* | 2567 | 0.0 | 88% |
| fungi | itpr | A0AAD5K078 | *Phascolomyces articulosus* | 2588 | 0.0 | 91% |
| fungi | itpr | A0AAD7UQL8 | *Lichtheimia ornata* | 2610 | 0.0 | 91% |
| fungi | itpr | A0ABR4NID5 | *Polyrhizophydium stewartii* | 4446 | 0.0 | 97% |
| fungi | itpr | A0ABR4NG49 | *Polyrhizophydium stewartii* | 3104 | 0.0 | 97% |
| fungi | itpr | A0A162V9Q6 | *Phycomyces blakesleeanus (strain ATCC 8743b / DSM 1359 / FGSC 10004 / NBRC 33097 / NRRL 1555)* | 2551 | 0.0 | 93% |
| fungi | itpr | A0ABR3B852 | *Phycomyces blakesleeanus* | 2551 | 0.0 | 93% |
| fungi | itpr | A0A1X2IXV9 | *Absidia repens* | 2651 | 0.0 | 97% |
| fungi | itpr | S2J437 | *Mucor circinelloides f. circinelloides (strain 1006PhL)* | 2592 | 0.0 | 97% |
| fungi | itpr | A0AAN7D8C4 | *Mucor velutinosus* | 2613 | 0.0 | 97% |
| fungi | itpr | A0A8H7VCW8 | *Mucor plumbeus* | 2604 | 0.0 | 97% |
| fungi | itpr | A0A1X2G6S3 | *Hesseltinella vesiculosa* | 2642 | 0.0 | 97% |
| fungi | itpr | A0A168KKP7 | *Mucor lusitanicus CBS 277.49* | 2620 | 0.0 | 97% |
| fungi | itpr | A0A0C9M4R7 | *Mucor ambiguus* | 2599 | 0.0 | 97% |


## Iterated search, per group (D10)

A single profile pass can only find what it is already close enough to. Each group was therefore searched again with `jackhmmer` from a seed native to that group — the highest-scoring full-length ITPR call in the group's own sweep — iterated to convergence under D10's coded kill criterion.

| group | seed | rounds | converged | D10 verdict | rule fired | targets in the accepted model |
|---|---|---|---|---|---|---|
| viridiplantae | A0AAE0LDZ1 | 5 | yes | clean | none | 70 |


### What each converged model is built from

The completeness statement the iterated search exists to make. A model seeded in one lineage and iterated to convergence either reaches a neighbouring lineage or it does not, and the targets **only iteration found** are the ones that matter: if they sit in the same lineage as the seed, the absence next door is not a sensitivity artefact.

| group | clade | profile call | targets | found by the single pass | iteration only |
|---|---|---|---|---|---|
| viridiplantae | Chlorophyta | ITPR | 54 | 54 | 0 |
| viridiplantae | Chlorophyta | unassigned | 10 | 10 | 0 |
| viridiplantae | Chlorophyta | not called by the single pass | 6 | 0 | 6 |


6 target(s) across all groups entered a converged model that the single profile pass never reported, and they fall in 1 clade(s) — the same ones the single pass already reached.


## Census v5

**17,882 records**, of which 785 are new from this sweep. Calls: ITPR 8,809, RYR 8,468, unassigned 605.


2 records changed their call against census v4, and 0 are `conflict` — kept and reported rather than resolved (D23).

| group | ITPR | RYR | conflict | unassigned | total |
|---|---|---|---|---|---|
| Amoebozoa | 12 | 1 | 0 | 0 | 13 |
| Bacteria | 0 | 0 | 0 | 1 | 1 |
| Discoba | 59 | 19 | 0 | 3 | 81 |
| Eukaryota (other) | 57 | 5 | 0 | 10 | 72 |
| Fungi | 35 | 20 | 0 | 14 | 69 |
| Metazoa (non-vertebrate) | 1860 | 2147 | 0 | 211 | 4218 |
| SAR | 610 | 99 | 0 | 20 | 729 |
| Vertebrata | 6112 | 6176 | 0 | 324 | 12612 |
| Viridiplantae | 64 | 1 | 0 | 22 | 87 |


## Figures

![jackhmmer_s20](figures/jackhmmer_s20.png)

![plant_fungal_chase](figures/plant_fungal_chase.png)

![profile_separation](figures/profile_separation.png)

![range_by_phylum](figures/range_by_phylum.png)


---

Groups: metazoa_nonvert, fungi, viridiplantae, protista_other (eukaryotes), archaea, bacteria_genus (prokaryotes).

