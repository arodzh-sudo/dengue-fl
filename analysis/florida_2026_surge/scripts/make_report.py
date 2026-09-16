#!/usr/bin/env python3
"""Write the Florida dengue report once, as HTML and as Word."""

import os
import struct
import zipfile
from xml.sax.saxutils import escape

import analyze_v2
from analyze_v2 import build_dir_argument


def report_dir():
    return os.path.join(analyze_v2.BASE, "report")


def figures_dir():
    return os.path.join(report_dir(), "figures")

TITLE = "Florida dengue 2026: what the genomes support"
SUBTITLE = ("A bioinformatics read of the 2026 Nextstrain build and of the interpretation "
            "returned by the Yale group")


def block(kind, **fields):
    fields["kind"] = kind
    return fields


def content():
    return [
        block("meta", rows=[
            ("Build", "2026-09-15"),
            ("Specimens sequenced", "117"),
            ("Genomes in the trees", "92"),
            ("Collection dates covered", "2025-09-01 to 2026-08-10"),
            ("Written", "2026-09-16"),
        ]),

        block("h2", text="Summary"),
        block("bullets", items=[
            "Of the 117 specimens sequenced, 92 produced a genome good enough to place on a "
            "tree: 54 DENV2, 22 DENV4, 14 DENV3 and 2 DENV1. Twenty of the 92 are recorded as "
            "acquired in Florida, one of them a mosquito pool.",
            "The Hillsborough and Pinellas outbreak is a single cluster of 13 closely related "
            "genomes, all locally acquired, 11 from Hillsborough and 2 from Pinellas. The "
            "Aedes aegypti pool MosquitoPool_K26-10948 belongs to the same cluster. The genomes "
            "differ from one another by 0 to 11 SNPs. A SNP here is a single-nucleotide difference "
            "between two genomes, counted along the tree, with insertions and deletions not "
            "counted.",
            "The pool genome is identical to the case genome TVU26000479 and differs by one "
            "SNP from TVU26000568. It is the only direct evidence in this dataset that "
            "the outbreak virus was present in local mosquitoes.",
            "The two Pinellas genomes are each other's closest relatives and sit inside the "
            "Hillsborough cluster, which is the pattern expected when the virus moved from one "
            "county to the other rather than arriving twice.",
            "The cluster descends from a part of the tree reconstructed as Cuba, where every "
            "Florida genome belongs to a traveler returning from Cuba. The closest of those "
            "travelers, TVU25001086, still differs by 19 SNPs, so the introduction "
            "itself was not sequenced.",
            "TVU26000552, the case whose origin was under investigation, differs by 3 "
            "SNPs from a Dade traveler returning from Cuba and by 36 to 42 from the "
            "outbreak cluster. It is not part of the outbreak and it is not related to the "
            "mosquito pool.",
            "Two of the conclusions we were sent rest on genomes that this build leaves out "
            "because too little of the genome assembled. TVU26000433 reached 36.6% of the "
            "genome and TVA26000028 reached 52.2%.",
        ]),

        block("h2", text="What changed since the build of 2026-09-09"),
        block("p", text=(
            "The numbers here differ from the earlier report because the analysis pipeline "
            "changed, not because new cases appeared. Four changes matter when the two are read "
            "side by side.")),
        block("bullets", items=[
            "Genomes that the quality checker flags for review now enter the build when at "
            "least 80% of the genome assembled. Fifteen of the 92 genomes arrive this way and "
            "four of them belong to the Hillsborough and Pinellas cluster, which is why that "
            "cluster holds 13 genomes here against 9 in the earlier report.",
            "A specimen sequenced more than once now enters the tree once, with the better run "
            "kept.",
            "A Florida genome already deposited in GenBank is matched to its public record "
            "instead of appearing twice. No genome in this run needed that treatment, since "
            "none of these samples has been deposited yet.",
            "The build now reads the sequencing run folders and the sample table directly, so "
            "every number in this report traces back to a file the pipeline read.",
        ]),
        block("p", text=(
            "This run covers specimens collected in 2025 and 2026. Florida genomes from 2022 to "
            "2024 remain available as public GenBank records, 53 of them for DENV2 and 446 for "
            "DENV3, so comparisons with earlier Florida cases are still possible.")),

        block("h2", text="From specimen to genome"),
        block("p", text=(
            "Sequencing turns a specimen into an assembled genome. Two separate checks then "
            "decide what that genome can be used for.")),
        block("p", text=(
            "The first check is the serotype call, which says which of the four dengue viruses "
            "the specimen carries. Four specimens never received one because too little virus "
            "was recovered to tell. Without a serotype there is no reference tree to place them "
            "on, so they go no further.")),
        block("p", text=(
            "The second check is VADR, the annotation tool that NCBI runs on viral genomes "
            "before they are accepted into GenBank. It compares the assembly against the "
            "reference gene models for that virus and reports whether the genes come out intact. "
            "The result is recorded as PASS, REVIEW or FAIL. REVIEW does not mean the sequence "
            "is wrong. It means the checker found something a person should look at, most often "
            "a gene that ends early or is interrupted, which is what missing sequence looks like "
            "to an annotation tool. The check describes how complete and well assembled a genome "
            "is, not which virus it is, so a genome can have a confident serotype call and still "
            "be marked REVIEW.")),
        block("p", text=(
            "The sequencing pipeline keeps every PASS assembly and, among the REVIEW assemblies, "
            "only those where at least 80% of the genome was reconstructed. This build uses "
            "both. The cutoff matters because a tree places a genome using the differences it "
            "can actually see. The less of a genome is present, the more of its position comes "
            "from inference rather than from observation. Below roughly 80% that position is "
            "not solid enough to read a transmission story from it. Eighteen REVIEW genomes fall "
            "below the cutoff and are left out.")),
        block("table",
              header=["Step", "Specimens", "Note"],
              rows=[
                  ["Sequenced", "117", "one per specimen"],
                  ["Serotype assigned", "113", "4 could not be typed"],
                  ["Quality check PASS", "77", "used"],
                  ["Quality check REVIEW", "33", "15 used, 18 below the 80% cutoff"],
                  ["Quality check FAIL", "7", "not used"],
                  ["Genomes in the trees", "92", "77 PASS and 15 REVIEW"],
              ]),
        block("figure", name="fig1_weeks", number=1, caption=(
            "Sequenced specimens collected in 2026, by week of collection. These are specimens "
            "that were sequenced, not reported cases, so the heights describe the sequencing "
            "workload rather than incidence. Source: results/sequenced_weeks_2026.tsv.")),
        block("p", text=(
            "One laboratory issue needs flagging. TVU25001126 was typed as DENV3 in one "
            "sequencing run and as DENV4 in its repeat, with 98.5% and 98.7% of the genome "
            "assembled against the respective references. Two runs this complete should not "
            "disagree about the serotype. The genome carried in the build under that identifier "
            "is the DENV4 one. The discrepancy is being checked in the laboratory.")),

        block("h2", text="DENV2 in Hillsborough and Pinellas"),
        block("p", text=(
            "Every locally acquired DENV2 case sequenced in Hillsborough and Pinellas belongs to "
            "one cluster of closely related genomes, together with the mosquito pool. Nothing else "
            "belongs to it. All 13 genomes in the cluster are Florida cases from those two "
            "counties. The branch leading to the cluster carries 15 SNPs, three of "
            "which change the protein sequence: NS5 T363I, M393T and S763N.")),
        block("figure", name="fig2_denv2_florida_tree", number=2, caption=(
            "All 54 Florida DENV2 genomes in the context of the global lineage that contains "
            "them, lineage 2II_F.1.1.2. Grey marks are public GenBank genomes. Filled circles "
            "are cases acquired in Florida, open circles are travel-associated cases and the "
            "diamond is the mosquito pool. Horizontal position is the date of collection. "
            "Source: dengue_denv2_genome.json.")),
        block("figure", name="fig3_cluster", number=3, caption=(
            "The Hillsborough and Pinellas cluster and its nearest travel-associated relative, on "
            "a time scale. Labels on the branches give the amino acid changes and, in brackets, "
            "the number of SNPs. The shaded band covers the range of dates the "
            "analysis gives for the common ancestor of the cluster. "
            "Source: results/denv2_cluster_structure.tsv.")),
        block("p", text=(
            "The structure inside the cluster is shallow, which is what recent local spread looks "
            "like. Twelve of the 13 genomes descend from a point dated to early March 2026. "
            "The smaller sub-groups below it date to May, June and July. The two Pinellas "
            "genomes, TVU26000565 and TVA26000031, are each other's closest relatives. They "
            "share the change E V484I on a branch dated 2026-06-26 that sits among Hillsborough "
            "genomes.")),
        block("p", text=(
            "The mosquito pool is identical to TVU26000479 across the genome and differs by one "
            "SNP from TVU26000568. Both of those genomes are REVIEW assemblies with 11.6% "
            "and 3.5% of their positions unresolved. Distances are counted along the tree, which "
            "fills unresolved positions with an inferred base, so identical here means "
            "indistinguishable at the positions that were actually read.")),
        block("figure", name="fig4_distances", number=4, caption=(
            "SNPs between every genome in the cluster and the three other cases "
            "of interest. Positions that were unresolved in an assembly carry an inferred base, "
            "so the bars on the right give the share of each assembly that was unresolved. "
            "Source: results/denv2_distances.tsv.")),

        block("h2", text="Why the two analyses differ on Pinellas"),
        block("p", text=(
            "Three Pinellas specimens were sequenced. Two of them, TVU26000565 and TVA26000031, "
            "are in this build and form the pair described above. The third, TVA26000028, "
            "assembled to 52.2% of the genome and was marked REVIEW, so it falls below the 80% "
            "cutoff and is absent here. The analysis carried out elsewhere includes it. That "
            "genome is what makes Pinellas appear to have two separate introductions.")),
        block("p", text=(
            "With roughly half of the genome missing, most of the positions that would fix its "
            "place on the tree are simply not there. The position it receives is driven by "
            "the fraction that remains. A genome that incomplete can land away from its true "
            "relatives, so a second introduction is one reading of it rather than the only one. "
            "On the genomes that can be placed with confidence, Pinellas looks like a single "
            "introduction out of Hillsborough. A second introduction is not ruled out. It cannot "
            "be settled by TVA26000028.")),
        block("p", text=(
            "The same caution applies to TVU26000433, which reached 36.6% of the genome. It is "
            "the genome used elsewhere to argue about the source of the Hillsborough cases. For "
            "the same reason this build cannot use it to answer that question.")),

        block("h2", text="Where the cluster came from"),
        block("p", text=(
            "The cluster descends from a point on the tree that the analysis reconstructs as Cuba "
            "with high confidence. The Florida genomes attached to that part of the tree, "
            "TVU25001086, TVU25001508, TVU26000429, TVU26000451 and TVA26000020, all belong to "
            "travelers returning from Cuba. That is the basis for the link to Cuba. It comes "
            "from travel histories recorded with Florida samples.")),
        block("p", text=(
            "Two things limit how far that can be taken. The reconstruction at the cluster's own "
            "ancestor reads USA with high confidence, but it has no choice: every genome below "
            "that point is a Florida case. And the closest sequenced traveler, TVU25001086, "
            "collected in October 2025, still differs by 19 SNPs, which is far too many "
            "for it to be the source. Whatever introduction started this outbreak was not among "
            "the specimens sequenced.")),
        block("figure", name="fig5a_exposure_tree", number=5, caption=(
            "The two clades that hold every Florida DENV2 genome of this build, on a "
            "time scale, with each genome colored by where the infection was acquired. The upper "
            "clade carries the Hillsborough and Pinellas outbreak, visible as the block of "
            "genomes acquired in Florida. Every genome in the lower one belongs to a traveler "
            "returning from Cuba. Source: results/denv2_sublineages.tsv.")),
        block("figure", name="fig5_lineage_context", number=6, caption=(
            "Genomes of lineage 2II_F.1.1.2 in this build, by country and year of collection. "
            "Public sequences are subsampled to 4000 per tree, so the counts describe where the "
            "lineage has been sequenced and shared, not where it is common. "
            "Source: results/lineage_2II_F_1_1_2.tsv.")),
        block("p", text=(
            "The wider history of this lineage, from South America into the Caribbean, comes "
            "from the published literature rather than from this build. We cannot check it "
            "here. What this build can say is that the lineage is well sampled in Brazil, "
            "Colombia, Costa Rica and Peru. Very few dengue genomes from Cuba are shared "
            "publicly, the most recent one here having been collected in 2022, so there is "
            "nothing recent from Cuba to compare the Florida outbreak against. The link to Cuba "
            "therefore rests on the travel histories reported with the Florida cases. The tree "
            "can neither confirm nor sharpen it.")),
        block("figure", name="fig7a_map_exposure", number=7, caption=(
            "The DENV2 tree drawn on a world map, colored by the country where infection was "
            "most likely acquired. Lines are inferred movements between countries, based on the "
            "travel histories recorded with each sample rather than on where it was collected."),
              screenshot=True),
        block("h2", text="Miami-Dade and the two clades behind it"),
        block("p", text=(
            "Dade contributes 32 of the 54 DENV2 genomes. Twenty-nine belong to travelers "
            "returning from Cuba, one to a traveler returning from Mexico and two were acquired "
            "in Florida. They do not form one group. Thirty of the 32 sit in two separate parts "
            "of the tree, far enough apart that they represent different arrivals rather than "
            "one.")),
        block("p", text=(
            "The first of the two carries the outbreak. Its common ancestor dates to May 2023 "
            "and the branch leading to it carries E S7A, NS2A D67G, NS2B I73V, NS5 I276V and "
            "H556Y. It holds 33 Florida genomes, 16 of them from Dade. It mixes travelers "
            "with locally acquired cases. Its closest public relatives are two Peruvian genomes "
            "from spring 2024. One step further out sit 34 genomes collected between 2022 "
            "and 2024, 17 from Colombia and 13 from Brazil. One public Florida genome, "
            "PQ155032 from March 2024, falls inside the clade itself, which means this virus was "
            "already reaching Florida two years before the current cases.")),
        block("figure", name="fig8_outbreak_clade", number=8, caption=(
            "The clade that carries the outbreak, with its nearest public relatives at the top. "
            "The Hillsborough and Pinellas cluster is drawn as a wedge so the Dade structure "
            "stays readable. TVU26000530 is the locally acquired Dade case. "
            "Source: results/denv2_clade_context.tsv.")),
        block("p", text=(
            "The second clade has only ever been seen in travelers. Its common ancestor dates to "
            "November 2023 and its branch carries NS2A I33L and NS4B N242S. It holds 19 Florida "
            "genomes, 14 of them from Dade. Every one belongs to a traveler returning from Cuba "
            "and no locally acquired case anywhere in Florida falls into it. Its neighbors are a "
            "different set: a Brazilian genome from August 2024, a Costa Rica genome from April "
            "2024 and two public Florida genomes from spring 2024, with Brazilian genomes from "
            "2023 below them.")),
        block("figure", name="fig9_traveler_clade", number=9, caption=(
            "The clade seen only in travelers, with its nearest public relatives at the top. "
            "TVU26000019 is the Dade case with travel to Mexico. "
            "Source: results/denv2_clade_context.tsv.")),
        block("p", text=(
            "The two therefore descend from different parts of the South American expansion of "
            "2II_F.1.1.2, one through a Peruvian and Colombian branch and one through a "
            "Brazilian branch. They are not two offshoots of a single Cuban virus. Whatever "
            "happened in Cuba, Florida has been receiving at least two distinct chains of "
            "arrival.")),
        block("p", text=(
            "The Mexico case belongs to neither. TVU26000019, collected in Dade on 2025-12-09, "
            "sits beside the traveler-only clade, 14 SNPs from that Costa Rica genome and 36 or "
            "more from everything else. It is a separate import from a separate journey. Dengue "
            "reaches Florida from more places than Cuba.")),
        block("p", text=(
            "Three pairs of Dade travelers are very close to one another. Both members of each "
            "pair reported travel to Cuba, so the simple explanation is a shared source there "
            "rather than transmission in Florida. A short unsampled local chain would look the "
            "same on a tree, which is why these pairs are worth naming. The collection dates and "
            "the case records will settle it faster than the genomes can.")),
        block("table",
              header=["Pair", "County", "Collected", "SNPs", "Clade"],
              rows=[
                  ["TVU26000471 and TVU26000504", "Dade and Dade", "2026-07-06, same day", "1",
                   "carries the outbreak"],
                  ["TVU26000426 and TVU26000448", "Dade and Dade", "2026-06-20 and 2026-06-26",
                   "1", "travelers only"],
                  ["TVU26000458 and TVU26000469", "Dade and Dade", "2026-07-08 and 2026-07-07",
                   "2", "travelers only"],
              ]),
        block("p", text=(
            "The two locally acquired Dade cases behave differently from each other. TVU26000530 "
            "sits inside the clade that carries the outbreak but 19 SNPs from anything else in "
            "it, so it represents its own introduction with no onward spread sequenced. "
            "TVU25001119, collected in October 2025, belongs to neither clade and is 38 SNPs "
            "from its closest relative in the build.")),
        block("p", text=(
            "Taken together, the Dade genomes show repeated importation into a county where two "
            "distinct groups are circulating among travelers. Nothing in them shows a sustained "
            "local chain of the kind seen in Hillsborough.")),

        block("h2", text="Locally acquired cases outside the cluster"),
        block("p", text=(
            "Four other locally acquired infections were sequenced. None of them belongs to the "
            "Hillsborough and Pinellas outbreak.")),
        block("table",
              header=["Sample", "Serotype", "County", "Collected", "Closest genome in the build"],
              rows=[
                  ["TVU26000530", "DENV2", "Dade", "2026-07-14",
                   "TVU26000545, a Dade traveler from Cuba, 19 SNPs"],
                  ["TVU25001119", "DENV2", "Dade", "2025-10-27",
                   "PQ155014, a public US genome from 2024, 38 SNPs"],
                  ["TVU26000521", "DENV3", "Dade", "2026-07-30",
                   "TVU26000373, a Dade traveler from Colombia, 41 SNPs"],
                  ["TVU26000522", "DENV3", "Dade", "2026-08-01",
                   "OQ821536, a Cuban genome from 2022, 37 SNPs"],
              ]),
        block("p", text=(
            "A distance of 19 to 41 SNPs to the closest available genome carries the same "
            "message in each case. The chain of transmission that produced these infections has "
            "not been sequenced, here or anywhere else the build can see. The two Dade DENV3 "
            "cases differ from each other by 80 SNPs, so they are unrelated to each other "
            "as well. They are recorded as locally acquired on the basis of the epidemiological "
            "investigation. The genomes neither support nor contradict that. They only show that "
            "no close relative has been sequenced.")),

        block("h2", text="Two cases with no origin recorded"),
        block("p", text=(
            "TVU26000552 and TVU25001298 reached the build without a recorded case origin, so "
            "they appear on the trees but carry no origin label. What the tree can say about "
            "each:")),
        block("bullets", items=[
            "TVU26000552, Hillsborough, DENV2, collected 2026-07-25, differs by 3 SNPs "
            "from TVU26000553, a Dade case with reported travel to Cuba collected two days "
            "later. It differs by 36 to 42 from every genome in the Hillsborough cluster. Its genome "
            "matches the same imported clade as TVU26000553 and does not match the local "
            "outbreak. The tree cannot say whether the infection was acquired abroad or in "
            "Florida from a shared source. Three SNPs is not close enough to call direct "
            "transmission between the two patients.",
            "TVU25001298, Orange, DENV3, collected 2025-11-11, differs by 20 SNPs from "
            "JVV25002080, a Palm Beach genome from December 2025. It differs by 30 or more from "
            "everything else. The tree points to no particular origin.",
        ]),

        block("h2", text="The other three serotypes"),
        block("p", text=(
            "DENV2 accounts for the surge, but it is not the only virus arriving. The other "
            "three serotypes add 38 genomes and a wider set of travel histories.")),
        block("figure", name="fig6_serotypes", number=10, caption=(
            "All 92 genomes by serotype and date of collection, colored by county. Filled marks "
            "are cases acquired in Florida and open marks are travel-associated cases. "
            "Source: results/tips.tsv.")),
        block("p", text=(
            "DENV4 is the second largest set, 22 genomes collected between September 2025 and "
            "January 2026, all but one of lineage 4II_B.1.3. Twenty belong to travelers "
            "returning from Cuba and two were acquired in Florida, both in Dade, in October and "
            "November 2025. Fifteen of the 22 were reported from Dade, the rest from Lee, "
            "Monroe, Pinellas, Sarasota and Palm Beach. The genomes are not a cluster. The "
            "closest pair differs by 11 SNPs and they scatter across a clade of about a hundred "
            "public genomes, so this is close to 20 separate introductions over five months "
            "rather than one chain of transmission. Figure 11 shows each of them beside the "
            "public genome it is closest to.")),
        block("figure", name="fig11_denv4_tree", number=11, caption=(
            "The 22 Florida DENV4 genomes with the public genomes nearest to each of them. "
            "Branches that carry no Florida genome are folded away, so a branch label counts "
            "every SNP along the path it replaces. Source: dengue_denv4_genome.json.")),
        block("p", text=(
            "DENV3 is the most varied set, 14 genomes collected between October 2025 and August "
            "2026. The travel histories point to Puerto Rico for four cases, Nicaragua for "
            "three, Colombia for one, Cuba for one and Indonesia for one. Three of the Puerto "
            "Rico travelers, reported from Broward, Dade and Palm Beach, sit together in one "
            "small part of the tree, which is what a shared source on the island looks like. "
            "Three cases were acquired in Florida: TVU25000991 in Brevard and the two unrelated "
            "Dade cases of July and August 2026 described earlier. Thirteen genomes belong to "
            "lineage 3III_B.3.2 and one to a different genotype altogether.")),
        block("figure", name="fig12_denv3_tree", number=12, caption=(
            "The 14 Florida DENV3 genomes with the public genomes nearest to each of them, drawn "
            "the same way as Figure 11. The travel histories line up with the tree: the Puerto "
            "Rico cases group together, as do the cases from Nicaragua. "
            "Source: dengue_denv3_genome.json.")),
        block("p", text=(
            "DENV1 appears twice, both travel-associated, both collected in March 2026. "
            "TVU26000280 was reported from Dade with travel to Nicaragua and TVU26000300 from "
            "Orange with travel to Colombia. They differ by 436 SNPs and belong to different "
            "genotypes, so they have nothing to do with each other beyond the serotype.")),

        block("h2", text="The conclusions we were sent, checked against this build"),
        block("table",
              header=["The statement", "What this build shows"],
              rows=[
                  ["Most local Florida cases fall into one cluster, with possible spread from "
                   "Hillsborough to Pinellas.",
                   "Supported. All 13 locally acquired Hillsborough and Pinellas genomes form "
                   "one cluster. The two Pinellas genomes pair together inside it, surrounded "
                   "by Hillsborough genomes."],
                  ["The Hillsborough cases most likely came from Cuba. The Pinellas cases came "
                   "from Hillsborough.",
                   "Partly supported. The parent part of the tree reconstructs as Cuba and every "
                   "Florida genome on it belongs to a Cuba traveler, but no sequenced import "
                   "comes within 19 SNPs of the cluster. County is never reconstructed in "
                   "this build, so the Pinellas direction rests on the shape of the tree and on "
                   "which specimens were sequenced."],
                  ["Lineage 2II_F.1.1.2 is expanding globally and reached the Caribbean around "
                   "2022 to 2023.",
                   "Outside what this build can test. It holds 258 genomes of the lineage, "
                   "concentrated in Brazil, Colombia, Costa Rica and Peru, with no recent "
                   "genomes from Cuba or the wider Caribbean available for comparison."],
                  ["The Hillsborough and Pinellas cases are a single cluster with spread between "
                   "counties. The Miami-Dade case is independent.",
                   "Supported. Differences within the cluster run from 0 to 11 SNPs while "
                   "the closest genome outside it differs by 19. TVU26000530 in Dade differs "
                   "from the cluster by 43 to 50."],
                  ["The cluster was introduced between August 2025 and January 2026, with local "
                   "transmission continuing through August.",
                   "Restated. On this set of genomes the common ancestor of the cluster dates to "
                   "2026-01-27, with a range of 2025-09-30 to 2026-03-22. That range describes "
                   "when the sampled genomes last shared an ancestor. It is an upper limit on "
                   "how early the cluster could have split rather than the date the virus reached "
                   "Florida. It moved when the four REVIEW genomes were added."],
                  ["TVU26000552 is unrelated to the other Hillsborough cases and to the mosquito "
                   "pool. The pool is identical to TVU26000479.",
                   "Supported on both points, with the caveat that 11.6% of TVU26000479 was "
                   "unresolved, so identical means indistinguishable at the positions read."],
                  ["TVU26000433 lies inside the local cluster, one SNP from TVA26000028, so "
                   "it cannot be the source.",
                   "Cannot be checked here. Neither genome is in this build. TVU26000433 reached "
                   "36.6% of the genome and TVA26000028 reached 52.2%, both below the cutoff "
                   "explained above. A position built on two genomes missing roughly two thirds "
                   "and one half of their sequence carries little information. Their closeness to "
                   "each other is partly a consequence of how much each is missing."],
                  ["Several travel-associated cases from Miami-Dade are genetically similar and "
                   "may share a source.",
                   "Confirmed and set out in the Miami-Dade section. Dade's genomes fall into "
                   "two separate clades and three pairs of travelers sit within 1 to 2 SNPs "
                   "of each other. Shared exposure in Cuba and an unsampled local chain look the "
                   "same on a tree, so the case records decide it rather than the genomes."],
              ]),

        block("h2", text="Limits on all of the above"),
        block("bullets", items=[
            "The 80% cutoff keeps 18 REVIEW genomes out of the build, including the two that "
            "several of the conclusions above were built on. Lowering it would place more "
            "specimens and place them less reliably.",
            "SNPs are counted along the tree, so positions that were unresolved in an "
            "assembly carry an inferred base. Any distance involving a REVIEW genome should be "
            "read together with its unresolved share, given in Figure 4.",
            "Public genomes are subsampled to 4000 per tree, clustered by year and region, while "
            "every Florida genome is kept. Counts of genomes describe sampling and deposition "
            "rather than how common a virus is.",
            "County and state appear on the map but are never reconstructed onto the tree, so "
            "the analysis draws no movement between counties. Any county to county statement is "
            "read off the shape of the tree by a person.",
            "The sample table covers specimens that were sequenced. It carries no case "
            "denominator, so nothing here can be expressed as a share of cases.",
            "Collection dates are given to the day and counties are named.",
        ]),

        block("h2", text="How the analysis was done"),
        block("p", text=(
            "Public context comes from GenBank. Every dengue genome in the database is "
            "downloaded, curated into a common format, split by serotype and assigned a lineage "
            "using the v-gen-lab Nextclade datasets, which is where names such as 2II_F.1.1.2 "
            "come from. That dataset is not version locked, so a later run can assign slightly "
            "different names.")),
        block("p", text=(
            "Florida genomes come from the Daytona_dengue sequencing runs. Our workflow reads "
            "each run's summary report, takes the assembly from the folder its quality grade "
            "places it in, joins the laboratory identifier to the sample table to pick up "
            "collection date, county, case origin and travel country, keeps a single run per "
            "specimen with PASS preferred over REVIEW. It also checks whether a genome has already "
            "been deposited in GenBank so that it enters the tree only once. Lineage names are "
            "carried across from the sequencing pipeline rather than recalculated.")),
        block("p", text=(
            "For each serotype the Florida genomes and the public ones are combined. Public "
            "genomes are subsampled to 4000 by year and region while every Florida genome is "
            "exempted from that subsampling. The sequences are aligned to the reference genome "
            "for that serotype with augur, a tree is built and then dated using a fixed "
            "evolutionary rate for that serotype, ancestral sequences and amino acid changes are "
            "reconstructed along the branches. Geography is reconstructed for region and "
            "country as well as for the travel-based exposure fields. County and state are "
            "carried for display only. The result is the set of files that Auspice opens.")),
        block("p", text=(
            "This report adds two scripts of its own. One reads the build and the run reports "
            "and writes every count, cluster membership, distance and date used here into a "
            "results folder. The other draws the figures from those files. The Auspice maps in "
            "Figures 6 and 7 are screenshots of the same build.")),

        block("footer", text=(
            "Sources: the DENV1 to DENV4 genome builds of 2026-09-15, together with the run "
            "summary report, the input and validation reports and the sample table from the same "
            "run. Every number was produced by analyze_v2.py and every drawn figure by "
            "make_figures.py, both in analysis/florida_2026_surge/scripts/. Counts of genomes "
            "are not counts of cases.")),
    ]


# ---------------------------------------------------------------- HTML
STYLE = """
:root { --surface: #ffffff; --ink: #0b0b0b; --soft: #52514e; --muted: #8a8a85;
        --rule: #e4e3df; --accent: #2a78d6; }
* { box-sizing: border-box; }
body { margin: 0; background: #f2f1ee; color: var(--ink);
       font-family: Inter, 'Helvetica Neue', Arial, sans-serif; line-height: 1.55; }
main { max-width: 980px; margin: 0 auto; padding: 48px 32px 96px; background: var(--surface);
       min-height: 100vh; }
h1 { font-size: 30px; line-height: 1.2; margin: 0 0 8px; }
.dek { font-size: 16px; color: var(--soft); margin: 0 0 24px; max-width: 70ch; }
h2 { font-size: 19px; margin: 44px 0 12px; padding-bottom: 6px; border-bottom: 2px solid var(--rule); }
p { max-width: 78ch; margin: 0 0 14px; }
ul { max-width: 78ch; padding-left: 20px; }
li { margin-bottom: 9px; }
table { border-collapse: collapse; width: 100%; margin: 8px 0 18px; font-size: 14px; }
th { text-align: left; font-weight: 600; color: var(--soft); border-bottom: 2px solid var(--rule);
     padding: 8px 10px; vertical-align: top; }
td { border-bottom: 1px solid var(--rule); padding: 8px 10px; vertical-align: top; }
td:first-child, th:first-child { padding-left: 0; }
figure { margin: 26px 0; }
figure img { width: 100%; height: auto; border: 1px solid var(--rule); border-radius: 6px;
             background: var(--surface); }
figcaption { font-size: 12.5px; color: var(--muted); margin-top: 8px; max-width: 86ch; }
figcaption b { color: var(--soft); }
.meta { display: flex; flex-wrap: wrap; gap: 26px; border-top: 1px solid var(--rule);
        border-bottom: 1px solid var(--rule); padding: 14px 0; margin-bottom: 8px; }
.meta div span { display: block; font-size: 11px; text-transform: uppercase;
                 letter-spacing: .04em; color: var(--muted); }
.meta div b { font-size: 15px; font-weight: 600; }
code { font-family: 'SF Mono', Menlo, Consolas, monospace; font-size: 13px; }
.footer { margin-top: 48px; padding-top: 14px; border-top: 1px solid var(--rule);
          font-size: 12.5px; color: var(--muted); max-width: 86ch; }
@media (max-width: 720px) { main { padding: 28px 16px 64px; } h1 { font-size: 24px; } }
"""


def render_html(blocks, path):
    out = ['<!DOCTYPE html><html lang="en"><head><meta charset="utf-8">',
           '<meta name="viewport" content="width=device-width, initial-scale=1">',
           "<title>{}</title><style>{}</style></head><body><main>".format(escape(TITLE), STYLE),
           "<h1>{}</h1>".format(escape(TITLE)),
           '<p class="dek">{}</p>'.format(escape(SUBTITLE))]
    for item in blocks:
        kind = item["kind"]
        if kind == "meta":
            out.append('<div class="meta">' + "".join(
                "<div><span>{}</span><b>{}</b></div>".format(escape(k), escape(v))
                for k, v in item["rows"]) + "</div>")
        elif kind == "h2":
            out.append("<h2>{}</h2>".format(escape(item["text"])))
        elif kind == "p":
            out.append("<p>{}</p>".format(mark(item["text"])))
        elif kind == "bullets":
            out.append("<ul>" + "".join("<li>{}</li>".format(mark(i)) for i in item["items"]) + "</ul>")
        elif kind == "table":
            head = "".join("<th>{}</th>".format(escape(c)) for c in item["header"])
            body = "".join("<tr>{}</tr>".format(
                "".join("<td>{}</td>".format(mark(c)) for c in row)) for row in item["rows"])
            out.append("<table><thead><tr>{}</tr></thead><tbody>{}</tbody></table>".format(head, body))
        elif kind == "figure":
            source = "figures/{}.{}".format(item["name"], "png" if item.get("screenshot") else "svg")
            out.append('<figure><img src="{}" alt="{}"><figcaption><b>Figure {}.</b> {}'
                       "</figcaption></figure>".format(source, escape(item["caption"][:120]),
                                                       item["number"], mark(item["caption"])))
        elif kind == "footer":
            out.append('<div class="footer">{}</div>'.format(mark(item["text"])))
    out.append("</main></body></html>")
    with open(path, "w", encoding="utf-8") as handle:
        handle.write("\n".join(out))
    print("wrote", os.path.basename(path))


def mark(text):
    """Escape, then set identifiers and file names in the monospace face."""
    escaped = escape(text)
    words = []
    for word in escaped.split(" "):
        bare = word.strip(".,;:()")
        if bare and (bare.startswith(("TVU", "TVA", "JVV", "MosquitoPool", "NODE_", "OQ", "PQ", "PP",
                                      "PV", "PX", "OR"))
                     or bare.endswith((".json", ".tsv", ".txt", ".py"))
                     or bare in ("case_origin", "metadata.txt", "augur", "2II_F.1.1.2",
                                 "3III_B.3.2", "4II_B.1.3")):
            word = word.replace(bare, "<code>{}</code>".format(bare))
        words.append(word)
    return " ".join(words)


# ---------------------------------------------------------------- Word
EMU_PER_PX = 9525
PAGE_WIDTH_EMU = 5943600


def png_size(path):
    with open(path, "rb") as handle:
        header = handle.read(24)
    return struct.unpack(">II", header[16:24])


def w_paragraph(runs, style=None, spacing=(0, 120)):
    properties = "<w:pPr>"
    if style:
        properties += '<w:pStyle w:val="{}"/>'.format(style)
    properties += '<w:spacing w:before="{}" w:after="{}"/></w:pPr>'.format(*spacing)
    return "<w:p>" + properties + runs + "</w:p>"


def w_run(text, bold=False, size=20, color="0B0B0B", mono=False):
    fonts = '<w:rFonts w:ascii="{0}" w:hAnsi="{0}"/>'.format("Consolas" if mono else "Calibri")
    return ('<w:r><w:rPr>{}{}<w:sz w:val="{}"/><w:color w:val="{}"/></w:rPr>'
            '<w:t xml:space="preserve">{}</w:t></w:r>'.format(
                fonts, "<w:b/>" if bold else "", size, color, escape(text)))


def w_image(index, width_emu, height_emu):
    return ('<w:p><w:pPr><w:spacing w:before="160" w:after="80"/></w:pPr><w:r><w:drawing>'
            '<wp:inline distT="0" distB="0" distL="0" distR="0">'
            '<wp:extent cx="{w}" cy="{h}"/><wp:docPr id="{i}" name="Figure {i}"/>'
            '<a:graphic xmlns:a="http://schemas.openxmlformats.org/drawingml/2006/main">'
            '<a:graphicData uri="http://schemas.openxmlformats.org/drawingml/2006/picture">'
            '<pic:pic xmlns:pic="http://schemas.openxmlformats.org/drawingml/2006/picture">'
            '<pic:nvPicPr><pic:cNvPr id="{i}" name="Figure {i}"/><pic:cNvPicPr/></pic:nvPicPr>'
            '<pic:blipFill><a:blip r:embed="rId{i}"/><a:stretch><a:fillRect/></a:stretch>'
            "</pic:blipFill>"
            '<pic:spPr><a:xfrm><a:off x="0" y="0"/><a:ext cx="{w}" cy="{h}"/></a:xfrm>'
            '<a:prstGeom prst="rect"><a:avLst/></a:prstGeom></pic:spPr></pic:pic>'
            "</a:graphicData></a:graphic></wp:inline></w:drawing></w:r></w:p>".format(
                w=width_emu, h=height_emu, i=index))


def w_table(header, rows):
    widths = [int(9000 / len(header))] * len(header)

    def cell(text, bold=False, shade=None):
        properties = '<w:tcPr><w:tcW w:w="{}" w:type="dxa"/>'.format(widths[0])
        if shade:
            properties += '<w:shd w:val="clear" w:fill="{}"/>'.format(shade)
        properties += "</w:tcPr>"
        return "<w:tc>" + properties + w_paragraph(w_run(text, bold=bold, size=18),
                                                   spacing=(40, 40)) + "</w:tc>"

    out = ['<w:tbl><w:tblPr><w:tblW w:w="9000" w:type="dxa"/>'
           '<w:tblBorders><w:top w:val="single" w:sz="4" w:color="E4E3DF"/>'
           '<w:bottom w:val="single" w:sz="4" w:color="E4E3DF"/>'
           '<w:insideH w:val="single" w:sz="4" w:color="E4E3DF"/></w:tblBorders></w:tblPr>']
    out.append("<w:tr>" + "".join(cell(c, bold=True, shade="F4F3F0") for c in header) + "</w:tr>")
    for row in rows:
        out.append("<w:tr>" + "".join(cell(c) for c in row) + "</w:tr>")
    out.append("</w:tbl>" + w_paragraph("", spacing=(0, 120)))
    return "".join(out)


def render_docx(blocks, path):
    missing = [item["name"] for item in blocks if item["kind"] == "figure"
               and not os.path.exists(os.path.join(figures_dir(), item["name"] + ".png"))]
    if missing:
        print("no PNG for {}, skipping the Word file".format(", ".join(missing)))
        print("run rasterize.sh once an SVG converter is available")
        return
    body, images = [], []
    body.append(w_paragraph(w_run(TITLE, bold=True, size=36)))
    body.append(w_paragraph(w_run(SUBTITLE, size=22, color="52514E")))

    for item in blocks:
        kind = item["kind"]
        if kind == "meta":
            body.append(w_table(["", ""], [[k, v] for k, v in item["rows"]]))
        elif kind == "h2":
            body.append(w_paragraph(w_run(item["text"], bold=True, size=26), spacing=(280, 100)))
        elif kind == "p":
            body.append(w_paragraph(w_run(item["text"])))
        elif kind == "bullets":
            for entry in item["items"]:
                body.append(w_paragraph(w_run("•  " + entry), spacing=(0, 80)))
        elif kind == "table":
            body.append(w_table(item["header"], item["rows"]))
        elif kind == "figure":
            name = item["name"] + ".png"
            source = os.path.join(figures_dir(), name)
            width, height = png_size(source)
            scale = min(PAGE_WIDTH_EMU / (width * EMU_PER_PX), 1.0)
            images.append(source)
            body.append(w_image(len(images), int(width * EMU_PER_PX * scale),
                                int(height * EMU_PER_PX * scale)))
            body.append(w_paragraph(
                w_run("Figure {}. ".format(item["number"]), bold=True, size=17, color="52514E")
                + w_run(item["caption"], size=17, color="8A8A85"), spacing=(0, 200)))
        elif kind == "footer":
            body.append(w_paragraph(w_run(item["text"], size=17, color="8A8A85"),
                                    spacing=(240, 0)))

    document = (
        '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
        '<w:document xmlns:w="http://schemas.openxmlformats.org/wordprocessingml/2006/main" '
        'xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships" '
        'xmlns:wp="http://schemas.openxmlformats.org/drawingml/2006/wordprocessingDrawing">'
        "<w:body>" + "".join(body) +
        '<w:sectPr><w:pgSz w:w="12240" w:h="15840"/>'
        '<w:pgMar w:top="1080" w:right="1080" w:bottom="1080" w:left="1080"/></w:sectPr>'
        "</w:body></w:document>")

    relationships = ['<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
                     '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">']
    overrides = []
    for index, source in enumerate(images, start=1):
        target = "media/{}".format(os.path.basename(source))
        relationships.append('<Relationship Id="rId{}" Type="http://schemas.openxmlformats.org/'
                             'officeDocument/2006/relationships/image" Target="{}"/>'.format(
                                 index, target))
        overrides.append('<Default Extension="png" ContentType="image/png"/>')
    relationships.append("</Relationships>")

    types = ('<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
             '<Types xmlns="http://schemas.openxmlformats.org/package/2006/content-types">'
             '<Default Extension="rels" ContentType="application/vnd.openxmlformats-package.'
             'relationships+xml"/><Default Extension="xml" ContentType="application/xml"/>'
             '<Default Extension="png" ContentType="image/png"/>'
             '<Override PartName="/word/document.xml" ContentType="application/vnd.'
             'openxmlformats-officedocument.wordprocessingml.document.main+xml"/></Types>')

    root_rels = ('<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
                 '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">'
                 '<Relationship Id="rId1" Type="http://schemas.openxmlformats.org/officeDocument/'
                 '2006/relationships/officeDocument" Target="word/document.xml"/></Relationships>')

    stem = path
    for attempt in range(1, 6):
        try:
            open(path, "ab").close()
            break
        except PermissionError:
            path = stem.replace(".docx", "-rebuilt{}.docx".format(attempt))
    if path != stem:
        print("the Word file is open, writing", os.path.basename(path), "instead")
    with zipfile.ZipFile(path, "w", zipfile.ZIP_DEFLATED) as archive:
        archive.writestr("[Content_Types].xml", types)
        archive.writestr("_rels/.rels", root_rels)
        archive.writestr("word/document.xml", document)
        archive.writestr("word/_rels/document.xml.rels", "".join(relationships))
        for source in images:
            archive.write(source, "word/media/" + os.path.basename(source))
    print("wrote", os.path.basename(path))


def main():
    blocks = content()
    render_html(blocks, os.path.join(report_dir(), "florida_dengue_2026_v3.html"))
    render_docx(blocks, os.path.join(report_dir(), "florida_dengue_2026_v3.docx"))


if __name__ == "__main__":
    build_dir_argument(__doc__)
    main()
