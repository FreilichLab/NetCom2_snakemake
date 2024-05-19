<img src="dep_sign.png" width=120, height=120 align="left" />

# NetCom2 - snakemake version workflow description

A pipeline for generating predictions for the selective targeting of microbial groups based on one the processing of assembled and annotated metagenomics. A detailed description of this workflow is provided at Ofaim S, Ofek-Lalzar M, Sela N, Jinag J, Kashi Y, Minz D, Freilich S: **Analysis of Microbial Functions in the Rhizosphere Using a Metabolic-Network Based Framework for Metagenomics Interpretation.** Frontiers in microbiology 2017, **8**:1606.**https://www.frontiersin.org/articles/10.3389/fmicb.2017.01606/full**.
Based on https://github.com/ot483/NetCom2.git
This a supplementary to main README. 

Key steps (snakemke rules) are described to facilitate with the understating of Netcom2 workflow.

## NetCom2 snakemake workflow

Netcom2 workflow is composed three steps (final output file):

1. determine differential abundance (EC\_diff\_abundance.tsv)
see EC\_diff\_abundance\_file\_graph.pdf for more details

2. Netcom

3. removal network of community knockout (merged\_netcom2.txt)
see netcom2\_file\_graph.pdf for more details

## determine differential abundance (EC\_diff\_abundance.tsv)

For this step Netcom2 inputs are used as in configuration file:

1. Gene taxonomy annotation:  gene\_taxonomy\_annotation.

2. Gene functional annotation: gene\_functional\_annotation

3. Contig count table: contig\_count\_table

Using *Removal\_pipline* rule contig level taxonomy assignment and EC count table is generated to determine differential abundance using edgeR tool.
The resulting file is the input for Netcom analysis is preformed using online tool.

## Netcom

The resulting file is, EC\_diff\_abundance.tsv, the input for Netcom analysis
is preformed using online tool.

Netcom provides the following outputs (relevant for Netcom2, see Readme.txt in Netcom results):

1.  Lists of differentially abundant enzymes ({treatment}\_ECs.txt), and their pathway
association ({treatment}\_Enzymes\_pathway.csv).

2. Prediction of environmental resources (seeds) that are unique to each
treatment ({treatment}\_resources.txt), and their pathway association ({treatment}\_Resources\_pathway.csv).

3. Prediction of environmental compounds that are produced by
the microbial community (i.e., expansion method) ({treatment}\_compounds.txt).  And pathway association of compounds
that are treatment-specific ({treatment}\_Compounds\_pathway.csv).

## Removal network of community knockout (merged\_netcom2.txt)

The final output of Necom2 are treatment specific removal networks.
These networks highlight the effects of community knockouts on metabolic networks found in Netcom analysis. 

Community knockouts are communities when the dominant taxa are removed.
As part of community knockout process complementary compounds and enzymes (elements) are found.
Complementary elements are those which are *not* found in dominant taxa and are found in the full set.
The full set are Netcom's input (for enzymes) and output (for environmental compounds (point 2)).
These data are stored in files containing "complementary" in their names.

Dominant taxa are found using the calculate\_enzyme\_diverity\_score rule.
And the effects of community knockouts on metabolic networks and corresponding pathways
are found using the rules Commuinty\_knockouts\_simulations and prepare\_pathway\_compounds\_for\_visual.
Finally the visualisation are created using the rule visualize\_knockout\_network.
Since these visitations are treatment specific, functional count table is split by treatment prior to calculation of  enzyme diversity score.

The final outputs are located in
results/community\_konckouts\_simulation in pdf and png formats as

(tretment\_vaule)\_treatment\_removal\_network\_knockout\_png.png and  

(tretment\_vaule)\_treatment\_removal\_network\_knockout\_pdf.pdf.

In these networks, legends are ordered by amount of compounds in pathways.

## To generate  NetCom2 snakemake workflow file graphs

```shell

# activate snakemake enviroment #

conda activate snakemake

# change directory to NetCom2_snakemake (working direcotry) #

cd NetCom2_snakemake

snakemake --filegraph 'results/edgeR/EC_diff_abundance.tsv' | dot -Tpng > EC_diff_abundance_file_graph.png

snakemake --filegraph results/merged_netcom2.txt  | dot -Tpdf > netcom2_file_graph.pdf
```

<img src="Network.jpg" width=600, height=600 align="center" />


## Contributors

[Gon Carmi](https://www.freilich-lab.com/members) \

[Shiri Freilich](https://www.freilich-lab.com/) 

## References

Ofaim S, Zarecki R, Porob S, Gat D, Lahav T, Kashi Y, Aly R, Eizenberg H, Ronen Z, Freilich S: **Genome-scale reconstruction of Paenarthrobacter aurescens TC1 metabolic model towards the study of atrazine bioremediation**. Sci Rep 2020, 10(1):13019.

Tal O, Bartuv R, Vetcos M, Medina S, Jiang J, Freilich S: **NetCom: A Network-Based Tool for Predicting Metabolic Activities of Microbial Communities Based on Interpretation of Metagenomics Data**. Microorganisms 2021, 9(9):1838.

Ofaim S, Ofek-Lalzar M, Sela N, Jinag J, Kashi Y, Minz D, Freilich S: **Analysis of Microbial Functions in the Rhizosphere Using a Metabolic-Network Based Framework for Metagenomics Interpretation.** Frontiers in microbiology 2017, **8**:1606.**https://www.frontiersin.org/articles/10.3389/fmicb.2017.01606/full**.

## Funding

This work was funded by the United States - Israel Binational Agricultural Research and Development Fund (BARD) [grant number [US-5390-21]

