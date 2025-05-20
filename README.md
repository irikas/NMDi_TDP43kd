# Data Analysis & Visualization for Project Using NMDi to ID Direct Targets of TDP-43 
<p> The scripts included in this repository were used to analyze data regarding the composition of TDP-43-associated cryptic exons, differential gene expression as a result of TDP-43 KD, toxicity of single target knockdown, and survival of cells after over-expression of selected TDP-43 targets in shTDP43-treated cells. Some raw data files are also included and others can be available upon request. </p> 
  
## Publication
<p>Sinha I.R.<sup>(#)</sup>, Ye, Y<sup>(#)</sup>, Li Y, Sandal PS, Burns G.D., Cruz AF, Sun S, Wong PC, Ling JP. Expanded Detection of Cryptic Exons Through Nonsense-Mediated Decay Inhibition Reveals Broad Landscape of TDP-43 Targets. <i>bioRxiv</i>. doi: TBD </p>

<p>Irika R. Sinha and Yingzhi Ye are co-first authors (#)</p>

## Contents
### Scripts
<table style="width:100% min-width: fit-content" border="0">
  <tr>
    <th>File Name</th>
    <th>Purpose</th>
  </tr>
  <tr>
    <td>250329_CrypticPSIfill_InitialVis.Rmd*</td>
    <td>Initial analysis and visualization of cryptic exon data (ex. PSI calculations, CE type PI charts, NMD effect on PSI) </td>
  </tr>
  <tr>
    <td>250409_CrypticPSI_CE_ALS_Genes.Rmd*</td>
    <td>Comparison of known ALS-FTD risk genes with CE-including genes for overlap</td>
  </tr>
  <tr>
    <td>250410_bulkseq.Rmd*</td>
    <td>Differential gene expression analysis after TDP-43 and adjustment for function transcripts</td>
  </tr>
  <tr>
    <td>250411_UGrepeat.Rmd*</td>
    <td>Create function to analyze dinucleotide repeats in cryptic exon splice sites +/- 600bp. Includes UG and CA repeat analysis.</td>
  </tr>
  <tr>
    <td>250425_shRNA_imaging.Rmd</td>
    <td>Calculation of shRNA toxicity & visualization</td>
  </tr>
  <tr>
    <td>250506_YY_Rescue.Rmd</td>
    <td>Calculation of survival after shTDP43 + target gene overexpression & visualization</td>
  </tr>
</table>
<p>*R Markdown document included as html. Download to open properly.</p>

### Data
<table style="width:100% min-width: fit-content">
  <tr>
    <th>File Name</th>
    <th>Description</th>
  </tr>
  <tr>
    <td>250409_ALSgenes.csv</td>
    <td>ALS-associated risk genes as per alsod.ac.uk</td>
  </tr>
  <tr>
    <td><p>250423_CellProfiler_Settings_p1-5.csv</p><p>250423_CellProfiler_Settings_p6-8.csv</p></td>
    <td>Exported CellProfiler settings for shRNA image analysis + quantification of nuclei and PI specks</td>
  </tr>
  <tr>
    <td>250509_CE_table.csv</td>
    <td>Identified cryptic exon genes, coordinates, splice junctions, and calculated PSIs</td>
  </tr>
  <tr>
    <td>250509_IS_Rescue_YY_Quant_Pool.csv</td>
    <td>Manual counts of nuclei and PI specks after rescue experiments</td>
  </tr>
  <tr>
    <td>GRch38sequence_report.tsv</td>
    <td>Includes information on chromosome length</td>
  </tr>
  <tr>
    <td>250410_nfcore_RNAseq_pipeline_report.html</td>
    <td>nf-core/RNAseq pipeline info output</td>
  </tr>
</table>

### Output tables
<table style="width:100% min-width: fit-content">
  <tr>
    <th>File Name</th>
    <th>Description</th>
  </tr>
  <tr>
    <td>250412_cryptic_df_UG_freq_L.csv</td>
    <td>[UG]<sub>n</sub> near 5' end of CEs</td>
  </tr>
    <tr>
    <td>250412_cryptic_df_UG_freq_R.csv</td>
    <td>[UG]<sub>n</sub> near 3' end of CEs</td>
  </tr>
</table>

### Packages Used
<ul>
  <li>CRAN: BiocManager, gghighlight, ggplot2, ggpubr, ggrepel, ggstats, pheatmap, tidyverse, viridis</li>
  <li>Bioconductor: apeglm, biomaRt, Biostrings, DESeq2, clusterProfiler, org.Hs.eg.db</li>
</ul>

### Pipelines Used
<a href="https://nf-co.re/rnaseq/3.14.0">nf-core/RNAseq</a>

## Acknowledgements
This work was supported in part by the Alzheimer’s Association (to J.P.L.), the Institute for Data-Intensive Engineering and Science (to J.P.L.), , the Robert Packard Center for ALS Research at Johns Hopkins (to P.C.W.), the Target ALS Foundation (to P.C.W.), ALS Finding a Cure (to P.C.W.), the ALS Association (to P.C.W.), and the US Food and Drug Administration (no. 1U01FD008129 to P.C.W.). 

This material is based upon work supported by the National Science Foundation (NSF) Graduate Research Fellowship Program under Grant No. DGE2139757 (to I.R.S.). Any opinions, findings, and conclusions or recommendations expressed in this material are those of the authors and do not necessarily reflect the views of the NSF.

This work was supported by resources from the Advanced Research Computing at Hopkins (ARCH) core facility (rockfish.jhu.edu), which is supported by the NSF (no. OAC 1920103).

We thank Katherine E. Irwin and Anya A. Kim for their troubleshooting support and suggestions.


