# Gene annotation table

Each probe can be located in multiple genomic feature classes: a probe can be in the promoter of gene A and an intron of gene B. The primary feature class of a probe is the feature class with the highest precedence, according to the following precedence order:

1. Promoter: 1500 bp upstream to 500 bp downstream of the TSS
2. UTRs (5' and 3')
3. Somewhere else in the gene body (exons or introns which are not located in the promoter or UTR regions)
4. Distant cis-regulatory domain (DCRD) around the TSS (+-100kb around the TSS)
5. Intergenic

All annotations for the feature class with highest precedence are shown in the table: if a probe is in the promoter region of two different genes, both genes are reported as primary hits, and added to the *gene_name* column as csv string.

Annotation was performed against the canoncical transcripts (based on the APPRIS score) in gencode M25.

The annotation table has the following columns:

Table 1: columns for the gene annotation tables.
| #Chromosome                | Chromosome; note the '#' to make this a BED header line                                                                                                   |
| Start                      | Probe interval start: 0-based, left-closed (BED format)                                                                                                   |
| End                        | Probe interval end: 0-based, right open (BED format)                                                                                                      |
| IlmnID                     | The Illumina probe ID (IlmnID field in manifest)                                                                                                          |
| Score                      | Empty ('.')                                                                                                                                               |
| Strand                     | Strand of the probed cytosine                                                                                                                             |
| Motif                      | One of: CG, CHH, CHG, D                                                                                                                                   |
| Genomic_region_class       | Highest precedence annotation class found for the probe, e.g. Promoter                                                                                    |
| Distance_to_genomic_region | Distance                                                                                                                                                  |
|                            | - from probe to TSS for promoter and DCRD annotations                                                                                                     |
|                            | - from probe to the interval center for exons, introns and UTRs                                                                                           |
|                            | (negative distances are upstream of the TSS)                                                                                                              |
| Gene_name                  | Gene symbol                                                                                                                                               |
| Gene_id                    | Ensemble gene id                                                                                                                                          |
| Transcript_id              | Ensemble transcript id                                                                                                                                    |
| Gene_strand                | Gene strand, if not intergenic                                                                                                                            |
| Cpg_island_region_class    | one of 'cpg_island', 'cpg_shore', 'cpg_shelve', 'open_sea', 'NA' (for un-annotatable probes)                                                              |
| Cpg_island_region_start    | 5'-coordinate of the closest CpG island                                                                                                                   |
| Cpg_island_region_end      | 3'-coordinate of the closest CpG island                                                                                                                   |
| Cpg_island_distance        | Distance of probe to nearest boundary of CpG island (0 if within CpG island), positive if upstream and negative if downstream of CpG island (on + strand) |


## Why provide multiple primary gene annotations for some probes?

It is relatively common that one probe lies in the promoter of multiple genes or in the intron of multiple genes etc. One could prioritize a single gene annotation based on the distance of the probe to the TSS. But biologically, the difference between a promoter location 500 bp away from the TSS or 1000 bp away from the TSS is minimal. We therefore prioritize only the feature class (e.g. only promoter annos, discard all lower precedence annos). We do not prioritize annotations within that feature class. The distance to the TSS or exon/interval center (see above) is still provided and could be used for filtering as required by the use case.

Table 2: Crosstabulation of the number of genes sharing the same highest-precedence feature class for the same probe. For example, 4704 probes are in the promoter of two genes.
|                       | number of annotated genes |       |      |      |      |      |     |
| Primary feature class | 1                         | 2     | 3    | 4    | 5    | 6    | ... |
| Promoter              | 56553                     | 4704  | 184  | 3    |      |      |     |
| 5'-UTR                | 1003                      | 4     |      |      |      |      |     |
| 3'-UTR                | 5284                      | 89    |      |      |      |      |     |
| exon                  | 26634                     | 186   | 3    |      |      |      |     |
| intron                | 90028                     | 1929  | 84   | 4    | 10   | 6    |     |
| DCRD                  | 21944                     | 13133 | 8339 | 5246 | 3632 | 2371 |     |
| intergenic            | 34506                     |       |      |      |      |      |     |
