# Count variants in each sample
This Python Jupyter notebook counts occurrences of each barcode in each sample from Illumina barcode sequencing, and adds these counts to the codon variant table.

## Set up analysis
### Import Python modules.
Use [plotnine](https://plotnine.readthedocs.io/en/stable/) for ggplot2-like plotting.

The analysis relies heavily on the Bloom lab's [dms_variants](https://jbloomlab.github.io/dms_variants) package:


```python
import itertools
import multiprocessing
import multiprocessing.pool
import os
import warnings

import alignparse
import alignparse.targets

import dms_variants.codonvarianttable
from dms_variants.constants import CBPALETTE
import dms_variants.illuminabarcodeparser
import dms_variants.utils
import dms_variants.plotnine_themes

from IPython.display import display, HTML

import pandas as pd

from plotnine import *

import yaml
```

Set [plotnine](https://plotnine.readthedocs.io/en/stable/) theme to the gray-grid one defined in `dms_variants`:


```python
theme_set(dms_variants.plotnine_themes.theme_graygrid())
```

Versions of key software:


```python
print(f"Using alignparse version {alignparse.__version__}")
print(f"Using dms_variants version {dms_variants.__version__}")
```

    Using alignparse version 0.6.0
    Using dms_variants version 1.4.3


Ignore warnings that clutter output:


```python
warnings.simplefilter('ignore')
```

### Parameters for notebook
Read the configuration file:


```python
with open('config.yaml') as f:
    config = yaml.safe_load(f)
```

Make output directory if needed:


```python
os.makedirs(config['counts_dir'], exist_ok=True)
os.makedirs(config['figs_dir'], exist_ok=True)
```

## Input variant tables
Initialize the table of barcode-variant pairs from the respective `process_ccs` notebooks for each background.


```python
variants = pd.read_csv(config['codon_variant_table_file_BA32'], na_filter=None)
#variants = variants.append(pd.read_csv(config['codon_variant_table_file_BA32'], na_filter=None))

variants = variants.reset_index(drop=True)

display(HTML(variants.head().to_html(index=False)))
```


<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th>target</th>
      <th>library</th>
      <th>barcode</th>
      <th>variant_call_support</th>
      <th>codon_substitutions</th>
      <th>aa_substitutions</th>
      <th>n_codon_substitutions</th>
      <th>n_aa_substitutions</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>BA32</td>
      <td>pool1</td>
      <td>AAAAAAAAAATTCAGT</td>
      <td>1</td>
      <td>TAC121ACT</td>
      <td>Y121T</td>
      <td>1</td>
      <td>1</td>
    </tr>
    <tr>
      <td>BA32</td>
      <td>pool1</td>
      <td>AAAAAAAAAGATTTAC</td>
      <td>9</td>
      <td>TAC35TGT</td>
      <td>Y35C</td>
      <td>1</td>
      <td>1</td>
    </tr>
    <tr>
      <td>BA32</td>
      <td>pool1</td>
      <td>AAAAAAAAAGTAGAAG</td>
      <td>18</td>
      <td>GTC11AAA</td>
      <td>V11K</td>
      <td>1</td>
      <td>1</td>
    </tr>
    <tr>
      <td>BA32</td>
      <td>pool1</td>
      <td>AAAAAAAAATTTCAAT</td>
      <td>2</td>
      <td>GCG33AGA</td>
      <td>A33R</td>
      <td>1</td>
      <td>1</td>
    </tr>
    <tr>
      <td>BA32</td>
      <td>pool1</td>
      <td>AAAAAAAACACCGACG</td>
      <td>1</td>
      <td></td>
      <td></td>
      <td>0</td>
      <td>0</td>
    </tr>
  </tbody>
</table>


Are there any barcodes in the same library that are shared across targets?
If so, we need to get rid of those as they will be confounded in barcode parsing:


```python
dup_barcodes = (
    variants
    .groupby(['library', 'barcode'])
    .size()
    .rename('duplicate_count')
    .reset_index()
    .query('duplicate_count > 1')
    )

print('Here are duplicated barcodes:')
display(HTML(dup_barcodes.head().to_html(index=False)))

print(f"\nRemoving the {len(dup_barcodes)} duplicated barcodes."
      f"Started with {len(variants)} barcodes:")
variants = (
    variants
    .merge(dup_barcodes, on=['library', 'barcode'], how='outer')
    .query('duplicate_count.isnull()', engine='python')
    )
print(f"After removing duplicates, there are {len(variants)} barcodes.")
```

    Here are duplicated barcodes:



<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th>library</th>
      <th>barcode</th>
      <th>duplicate_count</th>
    </tr>
  </thead>
  <tbody>
  </tbody>
</table>


    
    Removing the 0 duplicated barcodes.Started with 193880 barcodes:
    After removing duplicates, there are 193880 barcodes.


Pull out a target sequence for matching to the barcode and flanking sequence regions. Note, in this pipeline this is ok because our different backgrounds don't have differing flanks or other features within the actual N16 region covered in Illumina sequencing. If ever placing in-line barcodes here in the future, we would need to modify this.


```python
# get wildtype gene sequence for primary target
targets = alignparse.targets.Targets(seqsfile=config['amplicons_BA32'],
                                     feature_parse_specs=config['feature_parse_specs_BA32'])
```

## Setup to parse barcodes
Read data frame with list of all barcode runs.


```python
# barcode runs with R1 files by semicolon string split
barcode_runs = (pd.read_csv(config['barcode_runs'])
                .assign(R1=lambda x: x['R1'].str.split('; '))
                )

display(HTML(barcode_runs.to_html(index=False)))
```


<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th>library</th>
      <th>sample</th>
      <th>sample_type</th>
      <th>sort_bin</th>
      <th>concentration</th>
      <th>date</th>
      <th>number_cells</th>
      <th>R1</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>pool1</td>
      <td>TiteSeq_01_bin1</td>
      <td>TiteSeq</td>
      <td>1</td>
      <td>1.0</td>
      <td>260618</td>
      <td>1208755</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/ALT/2026/260723_illumina_SARS2-BA32/260618_s1-b1_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>TiteSeq_01_bin2</td>
      <td>TiteSeq</td>
      <td>2</td>
      <td>1.0</td>
      <td>260618</td>
      <td>680292</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/ALT/2026/260723_illumina_SARS2-BA32/260618_s1-b2_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>TiteSeq_01_bin3</td>
      <td>TiteSeq</td>
      <td>3</td>
      <td>1.0</td>
      <td>260618</td>
      <td>610637</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/ALT/2026/260723_illumina_SARS2-BA32/260618_s1-b3_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>TiteSeq_01_bin4</td>
      <td>TiteSeq</td>
      <td>4</td>
      <td>1.0</td>
      <td>260618</td>
      <td>2656753</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/ALT/2026/260723_illumina_SARS2-BA32/260618_s1-b4_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>TiteSeq_02_bin1</td>
      <td>TiteSeq</td>
      <td>1</td>
      <td>2.0</td>
      <td>260618</td>
      <td>1545467</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/ALT/2026/260723_illumina_SARS2-BA32/260618_s2-b1_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>TiteSeq_02_bin2</td>
      <td>TiteSeq</td>
      <td>2</td>
      <td>2.0</td>
      <td>260618</td>
      <td>872944</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/ALT/2026/260723_illumina_SARS2-BA32/260618_s2-b2_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>TiteSeq_02_bin3</td>
      <td>TiteSeq</td>
      <td>3</td>
      <td>2.0</td>
      <td>260618</td>
      <td>702650</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/ALT/2026/260723_illumina_SARS2-BA32/260618_s2-b3_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>TiteSeq_02_bin4</td>
      <td>TiteSeq</td>
      <td>4</td>
      <td>2.0</td>
      <td>260618</td>
      <td>2841176</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/ALT/2026/260723_illumina_SARS2-BA32/260618_s2-b4_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>TiteSeq_03_bin1</td>
      <td>TiteSeq</td>
      <td>1</td>
      <td>3.0</td>
      <td>260618</td>
      <td>1925289</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/ALT/2026/260723_illumina_SARS2-BA32/260618_s3-b1_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>TiteSeq_03_bin2</td>
      <td>TiteSeq</td>
      <td>2</td>
      <td>3.0</td>
      <td>260618</td>
      <td>847546</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/ALT/2026/260723_illumina_SARS2-BA32/260618_s3-b2_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>TiteSeq_03_bin3</td>
      <td>TiteSeq</td>
      <td>3</td>
      <td>3.0</td>
      <td>260618</td>
      <td>815421</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/ALT/2026/260723_illumina_SARS2-BA32/260618_s3-b3_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>TiteSeq_03_bin4</td>
      <td>TiteSeq</td>
      <td>4</td>
      <td>3.0</td>
      <td>260618</td>
      <td>2031835</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/ALT/2026/260723_illumina_SARS2-BA32/260618_s3-b4_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>TiteSeq_04_bin1</td>
      <td>TiteSeq</td>
      <td>1</td>
      <td>4.0</td>
      <td>260618</td>
      <td>2829658</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/ALT/2026/260723_illumina_SARS2-BA32/260618_s4-b1_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>TiteSeq_04_bin2</td>
      <td>TiteSeq</td>
      <td>2</td>
      <td>4.0</td>
      <td>260618</td>
      <td>993553</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/ALT/2026/260723_illumina_SARS2-BA32/260618_s4-b2_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>TiteSeq_04_bin3</td>
      <td>TiteSeq</td>
      <td>3</td>
      <td>4.0</td>
      <td>260618</td>
      <td>1022976</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/ALT/2026/260723_illumina_SARS2-BA32/260618_s4-b3_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>TiteSeq_04_bin4</td>
      <td>TiteSeq</td>
      <td>4</td>
      <td>4.0</td>
      <td>260618</td>
      <td>715195</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/ALT/2026/260723_illumina_SARS2-BA32/260618_s4-b4_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>TiteSeq_05_bin1</td>
      <td>TiteSeq</td>
      <td>1</td>
      <td>5.0</td>
      <td>260618</td>
      <td>3185438</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/ALT/2026/260723_illumina_SARS2-BA32/260618_s5-b1_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>TiteSeq_05_bin2</td>
      <td>TiteSeq</td>
      <td>2</td>
      <td>5.0</td>
      <td>260618</td>
      <td>1050530</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/ALT/2026/260723_illumina_SARS2-BA32/260618_s5-b2_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>TiteSeq_05_bin3</td>
      <td>TiteSeq</td>
      <td>3</td>
      <td>5.0</td>
      <td>260618</td>
      <td>774453</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/ALT/2026/260723_illumina_SARS2-BA32/260618_s5-b3_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>TiteSeq_05_bin4</td>
      <td>TiteSeq</td>
      <td>4</td>
      <td>5.0</td>
      <td>260618</td>
      <td>138623</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/ALT/2026/260723_illumina_SARS2-BA32/260618_s5-b4_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>TiteSeq_06_bin1</td>
      <td>TiteSeq</td>
      <td>1</td>
      <td>6.0</td>
      <td>260618</td>
      <td>5203442</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/ALT/2026/260723_illumina_SARS2-BA32/260618_s6-b1_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>TiteSeq_06_bin2</td>
      <td>TiteSeq</td>
      <td>2</td>
      <td>6.0</td>
      <td>260618</td>
      <td>610454</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/ALT/2026/260723_illumina_SARS2-BA32/260618_s6-b2_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>TiteSeq_06_bin3</td>
      <td>TiteSeq</td>
      <td>3</td>
      <td>6.0</td>
      <td>260618</td>
      <td>35989</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/ALT/2026/260723_illumina_SARS2-BA32/260618_s6-b3_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>TiteSeq_06_bin4</td>
      <td>TiteSeq</td>
      <td>4</td>
      <td>6.0</td>
      <td>260618</td>
      <td>3295</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/ALT/2026/260723_illumina_SARS2-BA32/260618_s6-b4_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>TiteSeq_07_bin1</td>
      <td>TiteSeq</td>
      <td>1</td>
      <td>7.0</td>
      <td>260618</td>
      <td>5627537</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/ALT/2026/260723_illumina_SARS2-BA32/260618_s7-b1_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>TiteSeq_07_bin2</td>
      <td>TiteSeq</td>
      <td>2</td>
      <td>7.0</td>
      <td>260618</td>
      <td>341866</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/ALT/2026/260723_illumina_SARS2-BA32/260618_s7-b2_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>TiteSeq_07_bin3</td>
      <td>TiteSeq</td>
      <td>3</td>
      <td>7.0</td>
      <td>260618</td>
      <td>11695</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/ALT/2026/260723_illumina_SARS2-BA32/260618_s7-b3_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>TiteSeq_07_bin4</td>
      <td>TiteSeq</td>
      <td>4</td>
      <td>7.0</td>
      <td>260618</td>
      <td>2784</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/ALT/2026/260723_illumina_SARS2-BA32/260618_s7-b4_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>TiteSeq_08_bin1</td>
      <td>TiteSeq</td>
      <td>1</td>
      <td>8.0</td>
      <td>260618</td>
      <td>5034985</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/ALT/2026/260723_illumina_SARS2-BA32/260618_s8-b1_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>TiteSeq_08_bin2</td>
      <td>TiteSeq</td>
      <td>2</td>
      <td>8.0</td>
      <td>260618</td>
      <td>324193</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/ALT/2026/260723_illumina_SARS2-BA32/260618_s8-b2_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>TiteSeq_08_bin3</td>
      <td>TiteSeq</td>
      <td>3</td>
      <td>8.0</td>
      <td>260618</td>
      <td>11598</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/ALT/2026/260723_illumina_SARS2-BA32/260618_s8-b3_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>TiteSeq_08_bin4</td>
      <td>TiteSeq</td>
      <td>4</td>
      <td>8.0</td>
      <td>260618</td>
      <td>2618</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/ALT/2026/260723_illumina_SARS2-BA32/260618_s8-b4_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>TiteSeq_09_bin1</td>
      <td>TiteSeq</td>
      <td>1</td>
      <td>9.0</td>
      <td>260618</td>
      <td>5299400</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/ALT/2026/260723_illumina_SARS2-BA32/260618_s9-b1_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>TiteSeq_09_bin2</td>
      <td>TiteSeq</td>
      <td>2</td>
      <td>9.0</td>
      <td>260618</td>
      <td>266926</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/ALT/2026/260723_illumina_SARS2-BA32/260618_s9-b2_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>TiteSeq_09_bin3</td>
      <td>TiteSeq</td>
      <td>3</td>
      <td>9.0</td>
      <td>260618</td>
      <td>5513</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/ALT/2026/260723_illumina_SARS2-BA32/260618_s9-b3_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>TiteSeq_09_bin4</td>
      <td>TiteSeq</td>
      <td>4</td>
      <td>9.0</td>
      <td>260618</td>
      <td>2298</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/ALT/2026/260723_illumina_SARS2-BA32/260618_s9-b4_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool2</td>
      <td>TiteSeq_01_bin1</td>
      <td>TiteSeq</td>
      <td>1</td>
      <td>1.0</td>
      <td>260624</td>
      <td>1014758</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/ALT/2026/260723_illumina_SARS2-BA32/260624_s11-b1_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool2</td>
      <td>TiteSeq_01_bin2</td>
      <td>TiteSeq</td>
      <td>2</td>
      <td>1.0</td>
      <td>260624</td>
      <td>852513</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/ALT/2026/260723_illumina_SARS2-BA32/260624_s11-b2_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool2</td>
      <td>TiteSeq_01_bin3</td>
      <td>TiteSeq</td>
      <td>3</td>
      <td>1.0</td>
      <td>260624</td>
      <td>720293</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/ALT/2026/260723_illumina_SARS2-BA32/260624_s11-b3_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool2</td>
      <td>TiteSeq_01_bin4</td>
      <td>TiteSeq</td>
      <td>4</td>
      <td>1.0</td>
      <td>260624</td>
      <td>3068144</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/ALT/2026/260723_illumina_SARS2-BA32/260624_s11-b4_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool2</td>
      <td>TiteSeq_02_bin1</td>
      <td>TiteSeq</td>
      <td>1</td>
      <td>2.0</td>
      <td>260624</td>
      <td>947108</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/ALT/2026/260723_illumina_SARS2-BA32/260624_s12-b1_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool2</td>
      <td>TiteSeq_02_bin2</td>
      <td>TiteSeq</td>
      <td>2</td>
      <td>2.0</td>
      <td>260624</td>
      <td>877246</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/ALT/2026/260723_illumina_SARS2-BA32/260624_s12-b2_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool2</td>
      <td>TiteSeq_02_bin3</td>
      <td>TiteSeq</td>
      <td>3</td>
      <td>2.0</td>
      <td>260624</td>
      <td>681516</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/ALT/2026/260723_illumina_SARS2-BA32/260624_s12-b3_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool2</td>
      <td>TiteSeq_02_bin4</td>
      <td>TiteSeq</td>
      <td>4</td>
      <td>2.0</td>
      <td>260624</td>
      <td>2735404</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/ALT/2026/260723_illumina_SARS2-BA32/260624_s12-b4_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool2</td>
      <td>TiteSeq_03_bin1</td>
      <td>TiteSeq</td>
      <td>1</td>
      <td>3.0</td>
      <td>260624</td>
      <td>1407844</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/ALT/2026/260723_illumina_SARS2-BA32/260624_s13-b1_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool2</td>
      <td>TiteSeq_03_bin2</td>
      <td>TiteSeq</td>
      <td>2</td>
      <td>3.0</td>
      <td>260624</td>
      <td>915716</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/ALT/2026/260723_illumina_SARS2-BA32/260624_s13-b2_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool2</td>
      <td>TiteSeq_03_bin3</td>
      <td>TiteSeq</td>
      <td>3</td>
      <td>3.0</td>
      <td>260624</td>
      <td>890870</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/ALT/2026/260723_illumina_SARS2-BA32/260624_s13-b3_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool2</td>
      <td>TiteSeq_03_bin4</td>
      <td>TiteSeq</td>
      <td>4</td>
      <td>3.0</td>
      <td>260624</td>
      <td>2398672</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/ALT/2026/260723_illumina_SARS2-BA32/260624_s13-b4_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool2</td>
      <td>TiteSeq_04_bin1</td>
      <td>TiteSeq</td>
      <td>1</td>
      <td>4.0</td>
      <td>260624</td>
      <td>1534482</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/ALT/2026/260723_illumina_SARS2-BA32/260624_s14-b1_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool2</td>
      <td>TiteSeq_04_bin2</td>
      <td>TiteSeq</td>
      <td>2</td>
      <td>4.0</td>
      <td>260624</td>
      <td>932797</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/ALT/2026/260723_illumina_SARS2-BA32/260624_s14-b2_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool2</td>
      <td>TiteSeq_04_bin3</td>
      <td>TiteSeq</td>
      <td>3</td>
      <td>4.0</td>
      <td>260624</td>
      <td>1154320</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/ALT/2026/260723_illumina_SARS2-BA32/260624_s14-b3_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool2</td>
      <td>TiteSeq_04_bin4</td>
      <td>TiteSeq</td>
      <td>4</td>
      <td>4.0</td>
      <td>260624</td>
      <td>1612220</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/ALT/2026/260723_illumina_SARS2-BA32/260624_s14-b4_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool2</td>
      <td>TiteSeq_05_bin1</td>
      <td>TiteSeq</td>
      <td>1</td>
      <td>5.0</td>
      <td>260624</td>
      <td>2827911</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/ALT/2026/260723_illumina_SARS2-BA32/260624_s15-b1_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool2</td>
      <td>TiteSeq_05_bin2</td>
      <td>TiteSeq</td>
      <td>2</td>
      <td>5.0</td>
      <td>260624</td>
      <td>1385571</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/ALT/2026/260723_illumina_SARS2-BA32/260624_s15-b2_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool2</td>
      <td>TiteSeq_05_bin3</td>
      <td>TiteSeq</td>
      <td>3</td>
      <td>5.0</td>
      <td>260624</td>
      <td>1080614</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/ALT/2026/260723_illumina_SARS2-BA32/260624_s15-b3_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool2</td>
      <td>TiteSeq_05_bin4</td>
      <td>TiteSeq</td>
      <td>4</td>
      <td>5.0</td>
      <td>260624</td>
      <td>166080</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/ALT/2026/260723_illumina_SARS2-BA32/260624_s15-b4_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool2</td>
      <td>TiteSeq_06_bin1</td>
      <td>TiteSeq</td>
      <td>1</td>
      <td>6.0</td>
      <td>260624</td>
      <td>4569238</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/ALT/2026/260723_illumina_SARS2-BA32/260624_s16-b1_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool2</td>
      <td>TiteSeq_06_bin2</td>
      <td>TiteSeq</td>
      <td>2</td>
      <td>6.0</td>
      <td>260624</td>
      <td>527668</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/ALT/2026/260723_illumina_SARS2-BA32/260624_s16-b2_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool2</td>
      <td>TiteSeq_06_bin3</td>
      <td>TiteSeq</td>
      <td>3</td>
      <td>6.0</td>
      <td>260624</td>
      <td>9381</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/ALT/2026/260723_illumina_SARS2-BA32/260624_s16-b3_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool2</td>
      <td>TiteSeq_06_bin4</td>
      <td>TiteSeq</td>
      <td>4</td>
      <td>6.0</td>
      <td>260624</td>
      <td>1812</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/ALT/2026/260723_illumina_SARS2-BA32/260624_s16-b4_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool2</td>
      <td>TiteSeq_07_bin1</td>
      <td>TiteSeq</td>
      <td>1</td>
      <td>7.0</td>
      <td>260624</td>
      <td>4977097</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/ALT/2026/260723_illumina_SARS2-BA32/260624_s17-b1_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool2</td>
      <td>TiteSeq_07_bin2</td>
      <td>TiteSeq</td>
      <td>2</td>
      <td>7.0</td>
      <td>260624</td>
      <td>343102</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/ALT/2026/260723_illumina_SARS2-BA32/260624_s17-b2_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool2</td>
      <td>TiteSeq_07_bin3</td>
      <td>TiteSeq</td>
      <td>3</td>
      <td>7.0</td>
      <td>260624</td>
      <td>4041</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/ALT/2026/260723_illumina_SARS2-BA32/260624_s17-b3_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool2</td>
      <td>TiteSeq_07_bin4</td>
      <td>TiteSeq</td>
      <td>4</td>
      <td>7.0</td>
      <td>260624</td>
      <td>1912</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/ALT/2026/260723_illumina_SARS2-BA32/260624_s17-b4_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool2</td>
      <td>TiteSeq_08_bin1</td>
      <td>TiteSeq</td>
      <td>1</td>
      <td>8.0</td>
      <td>260624</td>
      <td>5264860</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/ALT/2026/260723_illumina_SARS2-BA32/260624_s18-b1_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool2</td>
      <td>TiteSeq_08_bin2</td>
      <td>TiteSeq</td>
      <td>2</td>
      <td>8.0</td>
      <td>260624</td>
      <td>349821</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/ALT/2026/260723_illumina_SARS2-BA32/260624_s18-b2_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool2</td>
      <td>TiteSeq_08_bin3</td>
      <td>TiteSeq</td>
      <td>3</td>
      <td>8.0</td>
      <td>260624</td>
      <td>4001</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/ALT/2026/260723_illumina_SARS2-BA32/260624_s18-b3_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool2</td>
      <td>TiteSeq_08_bin4</td>
      <td>TiteSeq</td>
      <td>4</td>
      <td>8.0</td>
      <td>260624</td>
      <td>1759</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/ALT/2026/260723_illumina_SARS2-BA32/260624_s18-b4_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool2</td>
      <td>TiteSeq_09_bin1</td>
      <td>TiteSeq</td>
      <td>1</td>
      <td>9.0</td>
      <td>260624</td>
      <td>4991289</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/ALT/2026/260723_illumina_SARS2-BA32/260624_s19-b1_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool2</td>
      <td>TiteSeq_09_bin2</td>
      <td>TiteSeq</td>
      <td>2</td>
      <td>9.0</td>
      <td>260624</td>
      <td>223371</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/ALT/2026/260723_illumina_SARS2-BA32/260624_s19-b2_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool2</td>
      <td>TiteSeq_09_bin3</td>
      <td>TiteSeq</td>
      <td>3</td>
      <td>9.0</td>
      <td>260624</td>
      <td>2202</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/ALT/2026/260723_illumina_SARS2-BA32/260624_s19-b3_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool2</td>
      <td>TiteSeq_09_bin4</td>
      <td>TiteSeq</td>
      <td>4</td>
      <td>9.0</td>
      <td>260624</td>
      <td>1688</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/ALT/2026/260723_illumina_SARS2-BA32/260624_s19-b4_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>SortSeq_bin1</td>
      <td>SortSeq</td>
      <td>1</td>
      <td>NaN</td>
      <td>260618</td>
      <td>840000</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/ALT/2026/260723_illumina_SARS2-BA32/260618_exp1_RBD_bin1_1_R1_001.fastq.gz, /uufs/chpc.utah.edu/common/home/starr-group1/sequencing/ALT/2026/260723_illumina_SARS2-BA32/260618_exp1_RBD_bin1_2_R1_001.fastq.gz, /uufs/chpc.utah.edu/common/home/starr-group1/sequencing/ALT/2026/260723_illumina_SARS2-BA32/260618_exp1_RBD_bin1_3_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>SortSeq_bin2</td>
      <td>SortSeq</td>
      <td>2</td>
      <td>NaN</td>
      <td>260618</td>
      <td>2419472</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/ALT/2026/260723_illumina_SARS2-BA32/260618_exp1_RBD_bin2_1_R1_001.fastq.gz, /uufs/chpc.utah.edu/common/home/starr-group1/sequencing/ALT/2026/260723_illumina_SARS2-BA32/260618_exp1_RBD_bin2_2_R1_001.fastq.gz, /uufs/chpc.utah.edu/common/home/starr-group1/sequencing/ALT/2026/260723_illumina_SARS2-BA32/260618_exp1_RBD_bin2_3_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>SortSeq_bin3</td>
      <td>SortSeq</td>
      <td>3</td>
      <td>NaN</td>
      <td>260618</td>
      <td>2805239</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/ALT/2026/260723_illumina_SARS2-BA32/260618_exp1_RBD_bin3_1_R1_001.fastq.gz, /uufs/chpc.utah.edu/common/home/starr-group1/sequencing/ALT/2026/260723_illumina_SARS2-BA32/260618_exp1_RBD_bin3_2_R1_001.fastq.gz, /uufs/chpc.utah.edu/common/home/starr-group1/sequencing/ALT/2026/260723_illumina_SARS2-BA32/260618_exp1_RBD_bin3_3_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>SortSeq_bin4</td>
      <td>SortSeq</td>
      <td>4</td>
      <td>NaN</td>
      <td>260618</td>
      <td>3683689</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/ALT/2026/260723_illumina_SARS2-BA32/260618_exp1_RBD_bin4_1_R1_001.fastq.gz, /uufs/chpc.utah.edu/common/home/starr-group1/sequencing/ALT/2026/260723_illumina_SARS2-BA32/260618_exp1_RBD_bin4_2_R1_001.fastq.gz, /uufs/chpc.utah.edu/common/home/starr-group1/sequencing/ALT/2026/260723_illumina_SARS2-BA32/260618_exp1_RBD_bin4_3_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool2</td>
      <td>SortSeq_bin1</td>
      <td>SortSeq</td>
      <td>1</td>
      <td>NaN</td>
      <td>260624</td>
      <td>630000</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/ALT/2026/260723_illumina_SARS2-BA32/260624_exp2_RBD_bin1_1_R1_001.fastq.gz, /uufs/chpc.utah.edu/common/home/starr-group1/sequencing/ALT/2026/260723_illumina_SARS2-BA32/260624_exp2_RBD_bin1_2_R1_001.fastq.gz, /uufs/chpc.utah.edu/common/home/starr-group1/sequencing/ALT/2026/260723_illumina_SARS2-BA32/260624_exp2_RBD_bin1_3_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool2</td>
      <td>SortSeq_bin2</td>
      <td>SortSeq</td>
      <td>2</td>
      <td>NaN</td>
      <td>260624</td>
      <td>1947407</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/ALT/2026/260723_illumina_SARS2-BA32/260624_exp2_RBD_bin2_1_R1_001.fastq.gz, /uufs/chpc.utah.edu/common/home/starr-group1/sequencing/ALT/2026/260723_illumina_SARS2-BA32/260624_exp2_RBD_bin2_2_R1_001.fastq.gz, /uufs/chpc.utah.edu/common/home/starr-group1/sequencing/ALT/2026/260723_illumina_SARS2-BA32/260624_exp2_RBD_bin2_3_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool2</td>
      <td>SortSeq_bin3</td>
      <td>SortSeq</td>
      <td>3</td>
      <td>NaN</td>
      <td>260624</td>
      <td>2016169</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/ALT/2026/260723_illumina_SARS2-BA32/260624_exp2_RBD_bin3_1_R1_001.fastq.gz, /uufs/chpc.utah.edu/common/home/starr-group1/sequencing/ALT/2026/260723_illumina_SARS2-BA32/260624_exp2_RBD_bin3_2_R1_001.fastq.gz, /uufs/chpc.utah.edu/common/home/starr-group1/sequencing/ALT/2026/260723_illumina_SARS2-BA32/260624_exp2_RBD_bin3_3_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool2</td>
      <td>SortSeq_bin4</td>
      <td>SortSeq</td>
      <td>4</td>
      <td>NaN</td>
      <td>260624</td>
      <td>2102706</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/ALT/2026/260723_illumina_SARS2-BA32/260624_exp2_RBD_bin4_1_R1_001.fastq.gz, /uufs/chpc.utah.edu/common/home/starr-group1/sequencing/ALT/2026/260723_illumina_SARS2-BA32/260624_exp2_RBD_bin4_2_R1_001.fastq.gz, /uufs/chpc.utah.edu/common/home/starr-group1/sequencing/ALT/2026/260723_illumina_SARS2-BA32/260624_exp2_RBD_bin4_3_R1_001.fastq.gz]</td>
    </tr>
  </tbody>
</table>


Make sure library / sample combinations are unique:


```python
assert len(barcode_runs) == len(barcode_runs.groupby(['library', 'sample']))
```

Make sure the the libraries for which we have barcode runs are all in our variant table:


```python
unknown_libs = set(barcode_runs['library']) - set(variants['library'])
if unknown_libs:
    raise ValueError(f"Libraries with barcode runs not in variant table: {unknown_libs}")
```

Now we initialize an [IlluminaBarcodeParser](https://jbloomlab.github.io/dms_variants/dms_variants.illuminabarcodeparser.html#dms_variants.illuminabarcodeparser.IlluminaBarcodeParser) for each library.

First, get the length of the barcode from the alignment target after making sure the same length for all targets:


```python
bclen = len(targets.targets[0].get_feature('barcode').seq)

assert (bclen == len(target.get_feature('barcode').seq) for target in targets.targets)

print(f"Barcodes of length {bclen}")
```

    Barcodes of length 16


The other barcode parsing params come from the config file:


```python
parser_params = config['illumina_barcode_parser_params']

display(HTML(
    pd.Series(parser_params, name='value')
    .rename_axis(index='parameter')
    .reset_index()
    .to_html(index=False)
    ))
```


<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th>parameter</th>
      <th>value</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>upstream</td>
      <td>GGCCGC</td>
    </tr>
    <tr>
      <td>downstream</td>
      <td></td>
    </tr>
    <tr>
      <td>minq</td>
      <td>20</td>
    </tr>
    <tr>
      <td>upstream_mismatch</td>
      <td>1</td>
    </tr>
    <tr>
      <td>downstream_mismatch</td>
      <td>0</td>
    </tr>
  </tbody>
</table>


The parser needs to know the set of valid barcodes, which are stored in the variant table and are different for each library.
So we create a different parser for each library using these valid barcode sets:


```python
# create dict keyed by library, value is parser for library
parsers = {lib: dms_variants.illuminabarcodeparser.IlluminaBarcodeParser(
                    bclen=bclen,
                    valid_barcodes=variants.loc[variants['library']==lib]['barcode'],
                    **parser_params)
           for lib in set(variants['library'])}

print('Number of valid barcodes searched for by each parser:')
display(HTML(
    pd.DataFrame([(lib, len(p.valid_barcodes)) for lib, p in parsers.items()],
                 columns=['library', 'number of valid barcodes'])
    .to_html(index=False)
    ))
```

    Number of valid barcodes searched for by each parser:



<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th>library</th>
      <th>number of valid barcodes</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>pool1</td>
      <td>101800</td>
    </tr>
    <tr>
      <td>pool2</td>
      <td>92080</td>
    </tr>
  </tbody>
</table>


## Parse barcodes
We now parse the barcodes.
Since this will take a while, we utilize multiple CPUs via the Python [multiprocessing](https://docs.python.org/3.6/library/multiprocessing.html) module.
First, determine how many CPUs to use.
We use the minimum of the user-specified number hardcoded below and the number actually available.
(If you are running *interactively* on the Hutch cluster, you may need to reduce the number below in order to avoid an error as there is an enforced CPU limit on the home `rhino` nodes):


```python
ncpus = min(config['max_cpus'], multiprocessing.cpu_count())
print(f"Using {ncpus} CPUs")
```

    Using 8 CPUs


Parse the barcodes in parallel via a [multiprocessing.Pool](https://docs.python.org/3.6/library/multiprocessing.html#multiprocessing.pool.Pool) using all the available CPUs to get a list of the data frames with barcode counts / fates for each sample:


```python
def process_func(parser, r1files, library, sample):
    """Convenience function to be starmapped to multiprocessing pool."""
    return parser.parse(r1files, add_cols={'library': library, 'sample': sample})

# parallel computation of list of data frames
with multiprocessing.pool.Pool(processes=ncpus) as pool:
    bclist = pool.starmap(
                process_func,
                [(parsers[run.library], run.R1, run.library, run.sample)
                  for run in barcode_runs.itertuples()],
                )
```

Now concatenate the list into data frames of barcode counts and barcode fates:


```python
counts = pd.concat([samplecounts for samplecounts, _ in bclist],
                   sort=False,
                   ignore_index=True)

print('First few lines of counts data frame:')
display(HTML(counts.head().to_html(index=False)))

fates = pd.concat([samplefates for _, samplefates in bclist],
                  sort=False,
                  ignore_index=True)

print('First few lines of fates data frame:')
display(HTML(fates.head().to_html(index=False)))
```

    First few lines of counts data frame:



<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th>barcode</th>
      <th>count</th>
      <th>library</th>
      <th>sample</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>TAGGTAATGCCATACA</td>
      <td>1614</td>
      <td>pool1</td>
      <td>TiteSeq_01_bin1</td>
    </tr>
    <tr>
      <td>CTTCGCACAAACAACC</td>
      <td>1577</td>
      <td>pool1</td>
      <td>TiteSeq_01_bin1</td>
    </tr>
    <tr>
      <td>TTAGTCCTACTCCTGT</td>
      <td>1525</td>
      <td>pool1</td>
      <td>TiteSeq_01_bin1</td>
    </tr>
    <tr>
      <td>CTAGACTTCCGTCTAC</td>
      <td>1463</td>
      <td>pool1</td>
      <td>TiteSeq_01_bin1</td>
    </tr>
    <tr>
      <td>AGGAAACTCCAATGAA</td>
      <td>1386</td>
      <td>pool1</td>
      <td>TiteSeq_01_bin1</td>
    </tr>
  </tbody>
</table>


    First few lines of fates data frame:



<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th>fate</th>
      <th>count</th>
      <th>library</th>
      <th>sample</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>valid barcode</td>
      <td>2121228</td>
      <td>pool1</td>
      <td>TiteSeq_01_bin1</td>
    </tr>
    <tr>
      <td>invalid barcode</td>
      <td>848779</td>
      <td>pool1</td>
      <td>TiteSeq_01_bin1</td>
    </tr>
    <tr>
      <td>low quality barcode</td>
      <td>280924</td>
      <td>pool1</td>
      <td>TiteSeq_01_bin1</td>
    </tr>
    <tr>
      <td>unparseable barcode</td>
      <td>56802</td>
      <td>pool1</td>
      <td>TiteSeq_01_bin1</td>
    </tr>
    <tr>
      <td>failed chastity filter</td>
      <td>0</td>
      <td>pool1</td>
      <td>TiteSeq_01_bin1</td>
    </tr>
  </tbody>
</table>


## Examine fates of parsed barcodes
First, we'll analyze the "fates" of the parsed barcodes.
These fates represent what happened to each Illumina read we parsed:
 - Did the barcode read fail the Illumina chastity filter?
 - Was the barcode *unparseable* (i.e., the read didn't appear to be a valid barcode based on flanking regions)?
 - Was the barcode sequence too *low quality* based on the Illumina quality scores?
 - Was the barcode parseable but *invalid* (i.e., not in our list of variant-associated barcodes in the codon variant table)?
 - Was the barcode *valid*, and so will be added to variant counts.
 
First, we just write a CSV file with all the barcode fates:


```python
fatesfile = os.path.join(config['counts_dir'], 'barcode_fates.csv')
print(f"Writing barcode fates to {fatesfile}")
fates.to_csv(fatesfile, index=False)
```

    Writing barcode fates to results/counts/barcode_fates.csv


Next, we tabulate the barcode fates in wide format:


```python
display(HTML(fates
             .pivot_table(columns='fate',
                          values='count',
                          index=['library', 'sample'])
             .to_html()
             ))
```


<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>fate</th>
      <th>failed chastity filter</th>
      <th>invalid barcode</th>
      <th>low quality barcode</th>
      <th>unparseable barcode</th>
      <th>valid barcode</th>
    </tr>
    <tr>
      <th>library</th>
      <th>sample</th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th rowspan="40" valign="top">pool1</th>
      <th>SortSeq_bin1</th>
      <td>0</td>
      <td>590412</td>
      <td>192477</td>
      <td>327913</td>
      <td>1347260</td>
    </tr>
    <tr>
      <th>SortSeq_bin2</th>
      <td>0</td>
      <td>1639742</td>
      <td>578982</td>
      <td>161426</td>
      <td>4004398</td>
    </tr>
    <tr>
      <th>SortSeq_bin3</th>
      <td>0</td>
      <td>1992983</td>
      <td>700617</td>
      <td>180389</td>
      <td>5220095</td>
    </tr>
    <tr>
      <th>SortSeq_bin4</th>
      <td>0</td>
      <td>2630292</td>
      <td>856873</td>
      <td>150628</td>
      <td>7191005</td>
    </tr>
    <tr>
      <th>TiteSeq_01_bin1</th>
      <td>0</td>
      <td>848779</td>
      <td>280924</td>
      <td>56802</td>
      <td>2121228</td>
    </tr>
    <tr>
      <th>TiteSeq_01_bin2</th>
      <td>0</td>
      <td>476258</td>
      <td>162546</td>
      <td>41541</td>
      <td>1194318</td>
    </tr>
    <tr>
      <th>TiteSeq_01_bin3</th>
      <td>0</td>
      <td>385511</td>
      <td>125980</td>
      <td>19675</td>
      <td>1022183</td>
    </tr>
    <tr>
      <th>TiteSeq_01_bin4</th>
      <td>0</td>
      <td>2118762</td>
      <td>744744</td>
      <td>106028</td>
      <td>5852147</td>
    </tr>
    <tr>
      <th>TiteSeq_02_bin1</th>
      <td>0</td>
      <td>1205044</td>
      <td>405502</td>
      <td>120175</td>
      <td>3042368</td>
    </tr>
    <tr>
      <th>TiteSeq_02_bin2</th>
      <td>0</td>
      <td>528377</td>
      <td>184986</td>
      <td>35340</td>
      <td>1324339</td>
    </tr>
    <tr>
      <th>TiteSeq_02_bin3</th>
      <td>0</td>
      <td>559759</td>
      <td>164039</td>
      <td>35756</td>
      <td>1475791</td>
    </tr>
    <tr>
      <th>TiteSeq_02_bin4</th>
      <td>0</td>
      <td>2301453</td>
      <td>775505</td>
      <td>117205</td>
      <td>6308978</td>
    </tr>
    <tr>
      <th>TiteSeq_03_bin1</th>
      <td>0</td>
      <td>1529207</td>
      <td>515863</td>
      <td>127783</td>
      <td>3861283</td>
    </tr>
    <tr>
      <th>TiteSeq_03_bin2</th>
      <td>0</td>
      <td>604490</td>
      <td>213980</td>
      <td>180247</td>
      <td>1549142</td>
    </tr>
    <tr>
      <th>TiteSeq_03_bin3</th>
      <td>0</td>
      <td>641542</td>
      <td>235727</td>
      <td>37654</td>
      <td>1725966</td>
    </tr>
    <tr>
      <th>TiteSeq_03_bin4</th>
      <td>0</td>
      <td>1713285</td>
      <td>617964</td>
      <td>92102</td>
      <td>4707552</td>
    </tr>
    <tr>
      <th>TiteSeq_04_bin1</th>
      <td>0</td>
      <td>1997576</td>
      <td>689075</td>
      <td>234683</td>
      <td>4960056</td>
    </tr>
    <tr>
      <th>TiteSeq_04_bin2</th>
      <td>0</td>
      <td>705702</td>
      <td>256524</td>
      <td>48127</td>
      <td>1880894</td>
    </tr>
    <tr>
      <th>TiteSeq_04_bin3</th>
      <td>0</td>
      <td>731528</td>
      <td>249357</td>
      <td>36183</td>
      <td>2010967</td>
    </tr>
    <tr>
      <th>TiteSeq_04_bin4</th>
      <td>0</td>
      <td>537882</td>
      <td>209701</td>
      <td>35259</td>
      <td>1488459</td>
    </tr>
    <tr>
      <th>TiteSeq_05_bin1</th>
      <td>0</td>
      <td>2373777</td>
      <td>847247</td>
      <td>144256</td>
      <td>6226001</td>
    </tr>
    <tr>
      <th>TiteSeq_05_bin2</th>
      <td>0</td>
      <td>618332</td>
      <td>227490</td>
      <td>94918</td>
      <td>1672222</td>
    </tr>
    <tr>
      <th>TiteSeq_05_bin3</th>
      <td>0</td>
      <td>538595</td>
      <td>203916</td>
      <td>38692</td>
      <td>1485574</td>
    </tr>
    <tr>
      <th>TiteSeq_05_bin4</th>
      <td>0</td>
      <td>134241</td>
      <td>42359</td>
      <td>6937</td>
      <td>366640</td>
    </tr>
    <tr>
      <th>TiteSeq_06_bin1</th>
      <td>0</td>
      <td>3479332</td>
      <td>1303007</td>
      <td>207287</td>
      <td>9271543</td>
    </tr>
    <tr>
      <th>TiteSeq_06_bin2</th>
      <td>0</td>
      <td>384253</td>
      <td>131079</td>
      <td>28721</td>
      <td>1031676</td>
    </tr>
    <tr>
      <th>TiteSeq_06_bin3</th>
      <td>0</td>
      <td>29218</td>
      <td>10917</td>
      <td>5462</td>
      <td>70783</td>
    </tr>
    <tr>
      <th>TiteSeq_06_bin4</th>
      <td>0</td>
      <td>12953</td>
      <td>4576</td>
      <td>32759</td>
      <td>34822</td>
    </tr>
    <tr>
      <th>TiteSeq_07_bin1</th>
      <td>0</td>
      <td>3870106</td>
      <td>1456204</td>
      <td>231675</td>
      <td>10318502</td>
    </tr>
    <tr>
      <th>TiteSeq_07_bin2</th>
      <td>0</td>
      <td>271014</td>
      <td>98465</td>
      <td>16479</td>
      <td>719536</td>
    </tr>
    <tr>
      <th>TiteSeq_07_bin3</th>
      <td>0</td>
      <td>1298</td>
      <td>529</td>
      <td>175</td>
      <td>3019</td>
    </tr>
    <tr>
      <th>TiteSeq_07_bin4</th>
      <td>0</td>
      <td>4399</td>
      <td>2257</td>
      <td>16351</td>
      <td>10218</td>
    </tr>
    <tr>
      <th>TiteSeq_08_bin1</th>
      <td>0</td>
      <td>3271770</td>
      <td>1159440</td>
      <td>194134</td>
      <td>8799622</td>
    </tr>
    <tr>
      <th>TiteSeq_08_bin2</th>
      <td>0</td>
      <td>190275</td>
      <td>73940</td>
      <td>59446</td>
      <td>507840</td>
    </tr>
    <tr>
      <th>TiteSeq_08_bin3</th>
      <td>0</td>
      <td>608</td>
      <td>280</td>
      <td>161</td>
      <td>1252</td>
    </tr>
    <tr>
      <th>TiteSeq_08_bin4</th>
      <td>0</td>
      <td>4461</td>
      <td>1838</td>
      <td>3474</td>
      <td>10785</td>
    </tr>
    <tr>
      <th>TiteSeq_09_bin1</th>
      <td>0</td>
      <td>2854147</td>
      <td>983008</td>
      <td>156045</td>
      <td>7578083</td>
    </tr>
    <tr>
      <th>TiteSeq_09_bin2</th>
      <td>0</td>
      <td>175141</td>
      <td>58574</td>
      <td>9631</td>
      <td>469493</td>
    </tr>
    <tr>
      <th>TiteSeq_09_bin3</th>
      <td>0</td>
      <td>1734</td>
      <td>2556</td>
      <td>13889</td>
      <td>3391</td>
    </tr>
    <tr>
      <th>TiteSeq_09_bin4</th>
      <td>0</td>
      <td>169</td>
      <td>91</td>
      <td>4150</td>
      <td>237</td>
    </tr>
    <tr>
      <th rowspan="40" valign="top">pool2</th>
      <th>SortSeq_bin1</th>
      <td>0</td>
      <td>590841</td>
      <td>159192</td>
      <td>115664</td>
      <td>1022407</td>
    </tr>
    <tr>
      <th>SortSeq_bin2</th>
      <td>0</td>
      <td>1736341</td>
      <td>426127</td>
      <td>94885</td>
      <td>3094458</td>
    </tr>
    <tr>
      <th>SortSeq_bin3</th>
      <td>0</td>
      <td>2270324</td>
      <td>637943</td>
      <td>113379</td>
      <td>4157056</td>
    </tr>
    <tr>
      <th>SortSeq_bin4</th>
      <td>0</td>
      <td>2267981</td>
      <td>600680</td>
      <td>100225</td>
      <td>4227183</td>
    </tr>
    <tr>
      <th>TiteSeq_01_bin1</th>
      <td>0</td>
      <td>1080669</td>
      <td>293834</td>
      <td>48016</td>
      <td>1907042</td>
    </tr>
    <tr>
      <th>TiteSeq_01_bin2</th>
      <td>0</td>
      <td>833218</td>
      <td>203377</td>
      <td>37472</td>
      <td>1492518</td>
    </tr>
    <tr>
      <th>TiteSeq_01_bin3</th>
      <td>0</td>
      <td>867069</td>
      <td>213588</td>
      <td>36038</td>
      <td>1585093</td>
    </tr>
    <tr>
      <th>TiteSeq_01_bin4</th>
      <td>0</td>
      <td>3113174</td>
      <td>821089</td>
      <td>129088</td>
      <td>5840188</td>
    </tr>
    <tr>
      <th>TiteSeq_02_bin1</th>
      <td>0</td>
      <td>1093912</td>
      <td>299953</td>
      <td>49547</td>
      <td>1944138</td>
    </tr>
    <tr>
      <th>TiteSeq_02_bin2</th>
      <td>0</td>
      <td>827353</td>
      <td>227915</td>
      <td>52712</td>
      <td>1487709</td>
    </tr>
    <tr>
      <th>TiteSeq_02_bin3</th>
      <td>0</td>
      <td>735229</td>
      <td>206131</td>
      <td>34299</td>
      <td>1333993</td>
    </tr>
    <tr>
      <th>TiteSeq_02_bin4</th>
      <td>0</td>
      <td>2765434</td>
      <td>825149</td>
      <td>126871</td>
      <td>5155431</td>
    </tr>
    <tr>
      <th>TiteSeq_03_bin1</th>
      <td>0</td>
      <td>1573650</td>
      <td>435712</td>
      <td>73732</td>
      <td>2789067</td>
    </tr>
    <tr>
      <th>TiteSeq_03_bin2</th>
      <td>0</td>
      <td>876856</td>
      <td>231356</td>
      <td>38026</td>
      <td>1579333</td>
    </tr>
    <tr>
      <th>TiteSeq_03_bin3</th>
      <td>0</td>
      <td>1015315</td>
      <td>277987</td>
      <td>53578</td>
      <td>1842963</td>
    </tr>
    <tr>
      <th>TiteSeq_03_bin4</th>
      <td>0</td>
      <td>2595224</td>
      <td>740667</td>
      <td>113130</td>
      <td>4931865</td>
    </tr>
    <tr>
      <th>TiteSeq_04_bin1</th>
      <td>0</td>
      <td>1507049</td>
      <td>417025</td>
      <td>69883</td>
      <td>2661093</td>
    </tr>
    <tr>
      <th>TiteSeq_04_bin2</th>
      <td>0</td>
      <td>993925</td>
      <td>275530</td>
      <td>51017</td>
      <td>1794164</td>
    </tr>
    <tr>
      <th>TiteSeq_04_bin3</th>
      <td>0</td>
      <td>987336</td>
      <td>267068</td>
      <td>41909</td>
      <td>1815085</td>
    </tr>
    <tr>
      <th>TiteSeq_04_bin4</th>
      <td>0</td>
      <td>1679425</td>
      <td>478450</td>
      <td>73503</td>
      <td>3214329</td>
    </tr>
    <tr>
      <th>TiteSeq_05_bin1</th>
      <td>0</td>
      <td>2702476</td>
      <td>755250</td>
      <td>140394</td>
      <td>4831368</td>
    </tr>
    <tr>
      <th>TiteSeq_05_bin2</th>
      <td>0</td>
      <td>1252503</td>
      <td>345519</td>
      <td>54810</td>
      <td>2310120</td>
    </tr>
    <tr>
      <th>TiteSeq_05_bin3</th>
      <td>0</td>
      <td>1372838</td>
      <td>392208</td>
      <td>64972</td>
      <td>2625581</td>
    </tr>
    <tr>
      <th>TiteSeq_05_bin4</th>
      <td>0</td>
      <td>173288</td>
      <td>48838</td>
      <td>7251</td>
      <td>340057</td>
    </tr>
    <tr>
      <th>TiteSeq_06_bin1</th>
      <td>0</td>
      <td>4470002</td>
      <td>1284217</td>
      <td>229689</td>
      <td>8135217</td>
    </tr>
    <tr>
      <th>TiteSeq_06_bin2</th>
      <td>0</td>
      <td>593479</td>
      <td>165184</td>
      <td>26061</td>
      <td>1113917</td>
    </tr>
    <tr>
      <th>TiteSeq_06_bin3</th>
      <td>0</td>
      <td>33750</td>
      <td>7915</td>
      <td>4638</td>
      <td>46025</td>
    </tr>
    <tr>
      <th>TiteSeq_06_bin4</th>
      <td>0</td>
      <td>422</td>
      <td>96</td>
      <td>20</td>
      <td>631</td>
    </tr>
    <tr>
      <th>TiteSeq_07_bin1</th>
      <td>0</td>
      <td>6158500</td>
      <td>1737196</td>
      <td>280244</td>
      <td>11308562</td>
    </tr>
    <tr>
      <th>TiteSeq_07_bin2</th>
      <td>0</td>
      <td>350074</td>
      <td>101419</td>
      <td>44817</td>
      <td>648204</td>
    </tr>
    <tr>
      <th>TiteSeq_07_bin3</th>
      <td>0</td>
      <td>2837</td>
      <td>898</td>
      <td>55380</td>
      <td>3553</td>
    </tr>
    <tr>
      <th>TiteSeq_07_bin4</th>
      <td>0</td>
      <td>2187</td>
      <td>539</td>
      <td>46</td>
      <td>3653</td>
    </tr>
    <tr>
      <th>TiteSeq_08_bin1</th>
      <td>0</td>
      <td>4777252</td>
      <td>1248894</td>
      <td>253824</td>
      <td>8811087</td>
    </tr>
    <tr>
      <th>TiteSeq_08_bin2</th>
      <td>0</td>
      <td>306272</td>
      <td>86506</td>
      <td>17028</td>
      <td>573219</td>
    </tr>
    <tr>
      <th>TiteSeq_08_bin3</th>
      <td>0</td>
      <td>2310</td>
      <td>660</td>
      <td>1069</td>
      <td>3391</td>
    </tr>
    <tr>
      <th>TiteSeq_08_bin4</th>
      <td>0</td>
      <td>372</td>
      <td>277</td>
      <td>1224</td>
      <td>605</td>
    </tr>
    <tr>
      <th>TiteSeq_09_bin1</th>
      <td>0</td>
      <td>4531087</td>
      <td>1172431</td>
      <td>193816</td>
      <td>8287281</td>
    </tr>
    <tr>
      <th>TiteSeq_09_bin2</th>
      <td>0</td>
      <td>299531</td>
      <td>84429</td>
      <td>109719</td>
      <td>549223</td>
    </tr>
    <tr>
      <th>TiteSeq_09_bin3</th>
      <td>0</td>
      <td>2043</td>
      <td>608</td>
      <td>3955</td>
      <td>3293</td>
    </tr>
    <tr>
      <th>TiteSeq_09_bin4</th>
      <td>0</td>
      <td>1051</td>
      <td>2712</td>
      <td>21616</td>
      <td>594</td>
    </tr>
  </tbody>
</table>


Now we plot the barcode-read fates for each library / sample, showing the bars for valid barcodes in orange and the others in gray.
We see that the largest fraction of barcode reads correspond to valid barcodes, and most of the others are invalid barcodes (probably because the map to variants that aren't present in our variant table since we didn't associate all variants with barcodes). The exception to this is lib2 Titeseq_03_bin3; the PCR for this sample in the original sequencing run failed, so we followed it up with a single MiSeq lane. We did not filter out the PhiX reads from this data before parsing, so these PhiX reads will deflate the fraction of valid barcode reads as expected, but does not indicate any problems.


```python
barcode_fate_plot = (
    ggplot(
        fates
        .assign(sample=lambda x: pd.Categorical(x['sample'],
                                                x['sample'].unique(),
                                                ordered=True),
                fate=lambda x: pd.Categorical(x['fate'],
                                              x['fate'].unique(),
                                              ordered=True),
                is_valid=lambda x: x['fate'] == 'valid barcode'
                ), 
        aes('fate', 'count', fill='is_valid')) +
    geom_bar(stat='identity') +
    facet_grid('sample ~ library') +
    facet_grid('sample ~ library') +
    scale_fill_manual(CBPALETTE, guide=False) +
    theme(figure_size=(1.4 * (1 + fates['library'].nunique()),
                       1.7 * (1.2 + fates['sample'].nunique())),
          axis_text_x=element_text(angle=90),
          panel_grid_major_x=element_blank()
          ) +
    scale_y_continuous(labels=dms_variants.utils.latex_sci_not,
                       name='number of reads')
    )

_ = barcode_fate_plot.draw()
```


    
![png](count_variants_files/count_variants_42_0.png)
    


## Output csv of barcode counts in variant-barcode lookup table


```python
print(f"Writing variant counts to {config['variant_counts_file']}")
counts.to_csv(config['variant_counts_file'], index=False)
```

    Writing variant counts to results/counts/variant_counts.csv


The [CodonVariantTable](https://jbloomlab.github.io/dms_variants/dms_variants.codonvarianttable.html#dms_variants.codonvarianttable.CodonVariantTable) has lots of nice functions that can be used to analyze the counts it contains.
However, we do that in the next notebook so we don't have to re-run this entire (rather computationally intensive) notebook every time we want to analyze a new aspect of the counts.


```python

```
