# Bioinformatics Data Processing - Final assignment
## Overview
The typical metagenomics workflow involves several steps:
- DNA extraction
- Sequencing
- Assembly
- Genome binning
- Taxonomic and functional analysis

The assignment is based on data from two metagenomics assemblers **metaMDBG** and **myolasm** applied to nanopore sequencing data from a hot spring microbial community. The goal is to compare the performance of both assemblers using only the text-based metadata.
The metadata are from three sources:
- Contig FASTA headers
- **CheckM2** - genome completeness and contamination assessment
- **GTDBtk** - taxonomic assignment

The provided data are stored in the folder `results` in three subfolders: `checkm2`, `gtdbtk` and `contig_headers`.

For the analysis I will use Python with Jupyter Notebook and following libraries:
- pandas
- seaborn


```python
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
```

## FASTA headers
Each file contains the headers of the assembled contigs. On each line there is a single header.
### metaMDBG
First line of the file: `>ctg0 length=7710 coverage=2.25 circular=no`.

The [program manual](https://github.com/GaetanBenoitDev/metaMDBG) describes the headers:
- ctgID: the name of the contig
- length: the length of the contig in bps
- coverage: an estimated read coverage for the contig
- circular: whether the contig is circular or no

where each value variable is delimited by a space and the values are delimited by equals sign.


```python
meta_contigs = pd.read_csv('results/contig_headers/metamdbg_assembly_headers.txt', sep=' ', header=None,
                           names=['contig_id', 'length', 'coverage', 'circular'])
meta_contigs["contig_id"] = meta_contigs["contig_id"].str.replace(">", "")
meta_contigs["length"] = meta_contigs["length"].apply(lambda x: int(x.split('=')[1]))
meta_contigs["coverage"] = meta_contigs["coverage"].apply(lambda x: float(x.split('=')[1]))
meta_contigs["circular"] = meta_contigs["circular"].apply(lambda x: True if x.split('=')[1] == "yes" else False)
meta_contigs.head()
```




<div>

<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>contig_id</th>
      <th>length</th>
      <th>coverage</th>
      <th>circular</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>ctg0</td>
      <td>7710</td>
      <td>2.25</td>
      <td>False</td>
    </tr>
    <tr>
      <th>1</th>
      <td>ctg1</td>
      <td>8229</td>
      <td>1.96</td>
      <td>False</td>
    </tr>
    <tr>
      <th>2</th>
      <td>ctg2</td>
      <td>10850</td>
      <td>3.32</td>
      <td>False</td>
    </tr>
    <tr>
      <th>3</th>
      <td>ctg3</td>
      <td>9023</td>
      <td>4.85</td>
      <td>False</td>
    </tr>
    <tr>
      <th>4</th>
      <td>ctg4</td>
      <td>56052</td>
      <td>7.90</td>
      <td>False</td>
    </tr>
  </tbody>
</table>
</div>



There are 50109 contigs in the metaMDBG assembly.


```python
meta_contigs.nunique()
```




    contig_id    50109
    length       22387
    coverage      2056
    circular         2
    dtype: int64



### myolasm
First line of the file: `>u3840050ctg_len-10849_circular-no_depth-2-2-2_duplicated-no mult=1.00`.

The [program manual](https://myloasm-docs.github.io/output/) describes the headers:
- contig_id
- lenght: len-X - the length of the contig in nucleotides
- circularity: circular-X - whether the contig is circular or not
- depth: depth-X1-X2-X3 - the estimated depth of coverage:
    - X1: is the depth while allowing alignments of approximately > 99% true nucleotide similarity
    - X2: is the depth for > 99.75% similarity
    - X3: is the depth for 100% similarity
- duplicity: duplicated-X - whether the contig is duplicated or not
- mult: mult-X - is the estimated fraction of repetitiveness and useful for quality control. If mult is > 1.1, then duplicated is possibly. If > 1.5, then yes. Otherwise, it is no.

The variables are delimited by an underscore the values are delimited by a hyphen.


```python
mylo_contigs = pd.read_csv('results/contig_headers/myloasm_assembly_headers.txt', sep='_', header=None,
                           names=['contig_id', 'length', 'circular', 'coverage', 'duplicity_multiplicity'])
mylo_contigs["contig_id"] = mylo_contigs["contig_id"].str.replace(">", "")
mylo_contigs["length"] = mylo_contigs["length"].apply(lambda x: int(x.split('-')[1]))
mylo_contigs["circular"] = mylo_contigs["circular"].apply(lambda x: True if x.split('-')[1] == "yes" else False)
mylo_contigs["coverage"] = mylo_contigs["coverage"].apply(lambda x: int(x.split('-')[1]))
mylo_contigs["duplicity"] = mylo_contigs["duplicity_multiplicity"].apply(
    lambda x: True if x.split(' ')[0].split("-")[1] == "yes" else False)
mylo_contigs["multiplicity"] = mylo_contigs["duplicity_multiplicity"].apply(
    lambda x: float(x.split(' ')[1].split("=")[1]))
mylo_contigs.drop(columns=['duplicity_multiplicity'], inplace=True)
mylo_contigs.head()
```




<div>

<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>contig_id</th>
      <th>length</th>
      <th>circular</th>
      <th>coverage</th>
      <th>duplicity</th>
      <th>multiplicity</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>u3840050ctg</td>
      <td>10849</td>
      <td>False</td>
      <td>2</td>
      <td>False</td>
      <td>1.0</td>
    </tr>
    <tr>
      <th>1</th>
      <td>u418451ctg</td>
      <td>23116</td>
      <td>False</td>
      <td>2</td>
      <td>False</td>
      <td>1.0</td>
    </tr>
    <tr>
      <th>2</th>
      <td>u911112ctg</td>
      <td>78285</td>
      <td>False</td>
      <td>27</td>
      <td>False</td>
      <td>1.0</td>
    </tr>
    <tr>
      <th>3</th>
      <td>u2187011ctg</td>
      <td>16454</td>
      <td>False</td>
      <td>2</td>
      <td>False</td>
      <td>1.0</td>
    </tr>
    <tr>
      <th>4</th>
      <td>u3322747ctg</td>
      <td>14907</td>
      <td>False</td>
      <td>1</td>
      <td>False</td>
      <td>1.0</td>
    </tr>
  </tbody>
</table>
</div>



There are 53010 contigs in the myloasm assembly.


```python
mylo_contigs.nunique()
```




    contig_id       53010
    length          27047
    circular            2
    coverage          312
    duplicity           2
    multiplicity      113
    dtype: int64



The analysis of contig length distributions reveals a systematic difference in assembly contiguity: the metaMDBG assembly exhibits a distribution statistically skewed toward shorter lengths compared to the myolasm assembly. However, apart from this disparity in central tendency (mean/median length), the assemblers yield broadly comparable results across the measured range of contig sizes.


```python
circular_length_df = {"Circular": meta_contigs[meta_contigs["circular"]]["length"],
                      "Linear": meta_contigs[~meta_contigs["circular"]]["length"]}
sns.histplot(circular_length_df, log_scale=True)
plt.title("[metaMDBG] Contigs length distribution")
plt.xlabel("Contig length (bp)")
plt.yscale("log")
```


    
![png](BZD_files/BZD_11_0.png)
    



```python
circular_length_df = {"Circular": mylo_contigs[mylo_contigs["circular"]]["length"],
                      "Linear": mylo_contigs[~mylo_contigs["circular"]]["length"]}
sns.histplot(circular_length_df, log_scale=True)
plt.title("[myloasm] Contigs length distribution")
plt.xlabel("Contig length (bp)")
plt.yscale("log")
```


    
![png](BZD_files/BZD_12_0.png)
    


The scatterplots visualize the relationship between contig length and read coverage. For the metaMDBG assembly, a weak positive correlation is observed between contig length and coverage, strictly confined to the circular contigs. Conversely, the myolasm assembly demonstrates no discernible correlation between these variables in either circular or non-circular contigs.


```python
ax = sns.lmplot(x="length", y="coverage", data=meta_contigs, hue="circular", scatter_kws={'alpha': 0.6}, line_kws={"lw": 3})
plt.title("[metaMDBG] Contig length by coverage distribution")
plt.xlabel("Contig length (bp)")
plt.ylabel("Coverage")
plt.tight_layout()
```


    
![png](BZD_files/BZD_14_0.png)
    



```python
ax = sns.lmplot(x="length", y="coverage", data=mylo_contigs, hue="circular", scatter_kws={'alpha': 0.6}, line_kws={"lw": 3})
plt.title("[myloasm] Contig length by coverage distribution")
plt.xlabel("Contig length (bp)")
plt.ylabel("Coverage")
plt.tight_layout()
```


    
![png](BZD_files/BZD_15_0.png)
    


When focusing exclusively on long circular contigs (e.g., >500 kb), the normalized length histograms confirm that both assemblers produce contigs with a similar distribution profile. Nevertheless, the myolasm assembler successfully reconstructed nearly double the quantity of these target contigs compared to the metaMDBG assembler.


```python
meta_500 = meta_contigs[(meta_contigs["length"] > 500_000) & meta_contigs["circular"]].drop(columns=["circular"])
mylo_500 = mylo_contigs[(mylo_contigs["length"] > 500_000) & mylo_contigs["circular"]].drop(columns=["circular"])
print(f"Count of contigs:\n\tmetaMDBG: {meta_500.shape[0]},\n\tmyolasm: {mylo_500.shape[0]}")
sns.histplot({"metaMDBG": meta_500["length"], "myolasm": mylo_500["length"]}, stat='density', common_norm=False)
plt.title("Comparison of contig lengths for contigs longer than 500kb")
plt.xlabel("Contig length (bp)")
```

    Count of contigs:
    	metaMDBG: 36,
    	myolasm: 64
    




    Text(0.5, 0, 'Contig length (bp)')




    
![png](BZD_files/BZD_17_2.png)
    


Regarding read coverage for the long circular contigs, the myolasm assembly exhibits a coverage distribution that is distinctly skewed toward lower values. This suggests that the larger number of contigs reconstructed by myolasm (as noted previously) often corresponds to regions with lower sequence depth.


```python
sns.histplot({"metaMDBG": meta_500["coverage"], "myolasm": mylo_500["coverage"]}, stat='density', common_norm=False,
             log_scale=True)
plt.title("Comparison of contig coverages")
plt.xlabel("Contig coverage")
```




    Text(0.5, 0, 'Contig coverage')




    
![png](BZD_files/BZD_19_1.png)
    


## CheckM2
The next step after the assembly is genome binning. Metagenomics use genome bins as output, they are also known as MAGs (Metagenome-Assembled-Genomes). In the provided data, genome binning was not performed. The analysis is only performed to find (almost) complete MAGs, i.e. circular single contigs. CheckM2 evaluates genome/MAG completeness and contamination using set of markers. The columns of interest are **completeness** and **contaminations**.


```python
meta_m2 = pd.read_csv("results/checkm2/metamdbg/quality_report.tsv", delimiter="\t")[
    ["Name", "Completeness", "Contamination"]]
mylo_m2 = pd.read_csv("results/checkm2/myloasm/quality_report.tsv", delimiter="\t")[
    ["Name", "Completeness", "Contamination"]]
display(meta_m2.head())
display(mylo_m2.head())
print(f"Number of MAGs by assembler:\n\tmetaMDBG: {meta_m2.shape[0]},\n\tmyolasm: {mylo_m2.shape[0]}")
```


<div>

<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>Name</th>
      <th>Completeness</th>
      <th>Contamination</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>ctg10973</td>
      <td>100.00</td>
      <td>5.31</td>
    </tr>
    <tr>
      <th>1</th>
      <td>ctg11110</td>
      <td>95.23</td>
      <td>0.08</td>
    </tr>
    <tr>
      <th>2</th>
      <td>ctg11788</td>
      <td>96.85</td>
      <td>0.19</td>
    </tr>
    <tr>
      <th>3</th>
      <td>ctg12262</td>
      <td>7.77</td>
      <td>0.00</td>
    </tr>
    <tr>
      <th>4</th>
      <td>ctg12678</td>
      <td>96.32</td>
      <td>0.27</td>
    </tr>
  </tbody>
</table>
</div>



<div>

<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>Name</th>
      <th>Completeness</th>
      <th>Contamination</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>u1005671ctg</td>
      <td>50.78</td>
      <td>0.60</td>
    </tr>
    <tr>
      <th>1</th>
      <td>u1022233ctg</td>
      <td>95.41</td>
      <td>0.40</td>
    </tr>
    <tr>
      <th>2</th>
      <td>u1027066ctg</td>
      <td>89.94</td>
      <td>2.92</td>
    </tr>
    <tr>
      <th>3</th>
      <td>u1041087ctg</td>
      <td>11.27</td>
      <td>0.02</td>
    </tr>
    <tr>
      <th>4</th>
      <td>u1047122ctg</td>
      <td>91.01</td>
      <td>0.21</td>
    </tr>
  </tbody>
</table>
</div>


    Number of MAGs by assembler:
    	metaMDBG: 48,
    	myolasm: 103
    


```python
def quality(row):  # Helper function to define quality of the MAG
    if row["Completeness"] > 90 and row["Contamination"] < 5:
        return "High"
    elif row["Completeness"] > 50 and row["Contamination"] < 10:
        return "Medium"
    else:
        return "Low"


meta_m2["Quality"] = pd.Categorical(meta_m2.apply(quality, axis=1), categories=["Low", "Medium", "High"], ordered=True)
mylo_m2["Quality"] = pd.Categorical(mylo_m2.apply(quality, axis=1), categories=["Low", "Medium", "High"], ordered=True)
```


```python
g = sns.lmplot(meta_m2, x="Completeness", y="Contamination", hue="Quality", fit_reg=False)
plt.title("[metaMDBG] CheckM2 completeness and contamination")
g.ax.hlines(5, xmin=90, xmax=100, color="green", linestyle="--")
g.ax.vlines(90, ymin=0, ymax=5, color="green", linestyle="--")
g.ax.hlines(10, xmin=50, xmax=100, color="orange", linestyle="--")
g.ax.vlines(50, ymin=0, ymax=10, color="orange", linestyle="--")
print(meta_m2.Quality.value_counts().sort_index())
```

    Quality
    Low       12
    Medium     2
    High      34
    Name: count, dtype: int64
    


    
![png](BZD_files/BZD_23_1.png)
    



```python
g = sns.lmplot(mylo_m2, x="Completeness", y="Contamination", hue="Quality", fit_reg=False)
plt.title("[myloasm] CheckM2 completeness and contamination")
g.ax.hlines(5, xmin=90, xmax=100, color="green", linestyle="--")
g.ax.vlines(90, ymin=0, ymax=5, color="green", linestyle="--")
g.ax.hlines(10, xmin=50, xmax=100, color="orange", linestyle="--")
g.ax.vlines(50, ymin=0, ymax=10, color="orange", linestyle="--")
print(mylo_m2.Quality.value_counts().sort_index())
```

    Quality
    Low       29
    Medium    11
    High      63
    Name: count, dtype: int64
    


    
![png](BZD_files/BZD_24_1.png)
    



```python
completeness_df = {"metaMDBG": meta_m2["Completeness"], "myolasm": mylo_m2["Completeness"]}
sns.histplot(completeness_df, stat='density', common_norm=False)
plt.title("Comparison by completeness")
```




    Text(0.5, 1.0, 'Comparison by completeness')




    
![png](BZD_files/BZD_25_1.png)
    



```python
contamination_df = {"metaMDBG": meta_m2["Contamination"], "myolasm": mylo_m2["Contamination"]}
sns.histplot(contamination_df, stat='density', common_norm=False)
plt.title("Comparison by contamination")
```




    Text(0.5, 1.0, 'Comparison by contamination')




    
![png](BZD_files/BZD_26_1.png)
    



```python
meta_m2_qual = pd.DataFrame(meta_m2.Quality.value_counts()).reset_index()
meta_m2_qual["Name"] = "metaMDBG"
mylo_m2_qual = pd.DataFrame(mylo_m2.Quality.value_counts()).reset_index()
mylo_m2_qual["Name"] = "myloasm"
df_counts = pd.concat([meta_m2_qual, mylo_m2_qual])
sns.barplot(
    data=df_counts,
    x='Quality',
    y='count',
    hue='Name',
)
plt.title("Comparison by quality")
```




    Text(0.5, 1.0, 'Comparison by quality')




    
![png](BZD_files/BZD_27_1.png)
    


CheckM2 analysis reveals that contigs from the metaMDBG assembly are marginally more complete, though this finding is difficult to interpret definitively given the twofold greater number of contigs produced by the myolasm assembly. In terms of contamination, the myolasm assembly exhibits a slightly lower contamination rate, though the difference between the two assemblers is minimal.

## GTDB-Tk


```python
def prepare_dataset(location, contigs_df):
    df = pd.read_csv(location, delimiter="\t")[["user_genome", "classification"]]
    df = df[~df["classification"].str.contains("Unclassified")]
    df["classification"] = df["classification"].apply(lambda x: x.split(";")[1][3:])
    df = pd.merge(df, contigs_df, left_on=["user_genome"], right_on=["contig_id"], how="inner").drop(columns=["contig_id"])
    df = df[(df.length > 500_000) & (df.circular == True)]
    return df
```


```python
meta_gtdb_ar53 = prepare_dataset("results/gtdbtk/metamdbg/classify/gtdbtk.ar53.summary.tsv", meta_contigs)
meta_gtdb_bac120 = prepare_dataset("results/gtdbtk/metamdbg/classify/gtdbtk.bac120.summary.tsv", meta_contigs)
mylo_gtdb_ar53 = prepare_dataset("results/gtdbtk/myloasm/classify/gtdbtk.ar53.summary.tsv", mylo_contigs)
mylo_gtdb_bac120 = prepare_dataset("results/gtdbtk/myloasm/classify/gtdbtk.bac120.summary.tsv", mylo_contigs)
```

Taxonomic analysis demonstrates that both assemblers identified a highly similar number of archaeal MAGs, and the distribution of these identifications across phyla is consistent between the two assemblies. As figure below shows both assemblers found very similar number of Archaea genomes in the same phyla.


```python
meta_gtdb_ar53_counts = pd.DataFrame(meta_gtdb_ar53.classification.value_counts())
meta_gtdb_ar53_counts["Name"] = "metaMDBG"
mylo_gtdb_ar53_counts = pd.DataFrame(mylo_gtdb_ar53.classification.value_counts())
mylo_gtdb_ar53_counts["Name"] = "myloasm"
df_counts = pd.concat([meta_gtdb_ar53_counts, mylo_gtdb_ar53_counts])
sns.barplot(df_counts, x="classification", y="count", hue="Name")
plt.title("Comparison of number of large contigs per Archaea phylum between assemblers")
```




    Text(0.5, 1.0, 'Comparison of number of large contigs per Archaea phylum between assemblers')




    
![png](BZD_files/BZD_32_1.png)
    



```python
pd.concat([meta_gtdb_ar53_counts.drop(columns=["Name"]).rename({"count": "metaMDBG"}, axis=1),
           mylo_gtdb_ar53_counts.drop(columns=["Name"]).rename({"count": "myloasm"}, axis=1)], axis=1).fillna(0).astype(int)
```




<div>

<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>metaMDBG</th>
      <th>myloasm</th>
    </tr>
    <tr>
      <th>classification</th>
      <th></th>
      <th></th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>Aenigmatarchaeota</th>
      <td>4</td>
      <td>3</td>
    </tr>
    <tr>
      <th>Micrarchaeota</th>
      <td>3</td>
      <td>3</td>
    </tr>
    <tr>
      <th>Thermoproteota</th>
      <td>2</td>
      <td>4</td>
    </tr>
    <tr>
      <th>Nanobdellota</th>
      <td>2</td>
      <td>1</td>
    </tr>
  </tbody>
</table>
</div>



The figure below clearly illustrates that the myolasm assembler recovered a substantially greater number of bacterial MAGs compared to the metaMDBG assembler. Furthermore, while the taxonomic diversity was largely similar, the myolasm assembly contained MAGs from two unique phyla, whereas metaMDBG uniquely contributed one distinct phylum.


```python
meta_gtdb_bac120_counts = pd.DataFrame(meta_gtdb_bac120.classification.value_counts())
meta_gtdb_bac120_counts["Name"] = "metaMDBG"
mylo_gtdb_bac120_counts = pd.DataFrame(mylo_gtdb_bac120.classification.value_counts())
mylo_gtdb_bac120_counts["Name"] = "myloasm"
df_counts = pd.concat([meta_gtdb_bac120_counts, mylo_gtdb_bac120_counts])
ax = sns.barplot(df_counts, x="classification", y="count", hue="Name")
ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha='right')
plt.title("Comparison of number of large contigs per Bacteria phylum between assemblers")
plt.tight_layout()
```

    C:\Users\vilif\AppData\Local\Temp\ipykernel_15652\4155299355.py:7: UserWarning: set_ticklabels() should only be used with a fixed number of ticks, i.e. after set_ticks() or using a FixedLocator.
      ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha='right')
    


    
![png](BZD_files/BZD_35_1.png)
    



```python
pd.concat([meta_gtdb_bac120_counts.drop(columns=["Name"]).rename({"count": "metaMDBG"}, axis=1),
           mylo_gtdb_bac120_counts.drop(columns=["Name"]).rename({"count": "myloasm"}, axis=1)],
          axis=1).fillna(0).astype(int)
```




<div>

<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>metaMDBG</th>
      <th>myloasm</th>
    </tr>
    <tr>
      <th>classification</th>
      <th></th>
      <th></th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>Patescibacteriota</th>
      <td>16</td>
      <td>29</td>
    </tr>
    <tr>
      <th>Desulfobacterota</th>
      <td>2</td>
      <td>3</td>
    </tr>
    <tr>
      <th>Omnitrophota</th>
      <td>2</td>
      <td>8</td>
    </tr>
    <tr>
      <th>Firestonebacteria</th>
      <td>1</td>
      <td>0</td>
    </tr>
    <tr>
      <th>Planctomycetota</th>
      <td>1</td>
      <td>1</td>
    </tr>
    <tr>
      <th>Zixibacteria</th>
      <td>1</td>
      <td>1</td>
    </tr>
    <tr>
      <th>Chloroflexota</th>
      <td>1</td>
      <td>2</td>
    </tr>
    <tr>
      <th>Electryoneota</th>
      <td>1</td>
      <td>1</td>
    </tr>
    <tr>
      <th>Nitrospirota</th>
      <td>0</td>
      <td>2</td>
    </tr>
    <tr>
      <th>Bacillota</th>
      <td>0</td>
      <td>1</td>
    </tr>
    <tr>
      <th>Bipolaricaulota</th>
      <td>0</td>
      <td>1</td>
    </tr>
    <tr>
      <th>Caldisericota</th>
      <td>0</td>
      <td>1</td>
    </tr>
    <tr>
      <th>Acidobacteriota</th>
      <td>0</td>
      <td>1</td>
    </tr>
    <tr>
      <th>TA06</th>
      <td>0</td>
      <td>1</td>
    </tr>
    <tr>
      <th>UBA9089</th>
      <td>0</td>
      <td>1</td>
    </tr>
  </tbody>
</table>
</div>



### Summary
The myolasm assembly demonstrated superior contiguity and yield, producing a higher total number of contigs and nearly double the quantity of identifiable MAGs. These MAGs maintained a comparable quality distribution to those derived from the metaMDBG assembly.

In terms of taxonomic analysis (GTDB-Tk), both assemblers performed similarly in the Archaea domain. However, the myolasm assembly was significantly more successful in the Bacteria domain, recovering a greater absolute number of MAGs and resolving them into a more diverse set of phyla.
