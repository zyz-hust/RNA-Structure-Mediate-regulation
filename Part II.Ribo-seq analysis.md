# 0）数据说明
Ribo-seq测序数据位于`/data/TA_QUIZ_RNA_regulation/data/riboshape_liulab_batch4/Ribo-seq` 下,共12套。

# 1）Mapping
Ribo-seq测序数据的`fastqc`、`trimmed`、`remove_rRNA`和`mapping`流程与`differential expression`分析的流程大体一致。

这里我们采用了**Arabidopsis_thaliana.TAIR10.34.gtf**的注释文件。

|files name| path|
|:------:|:------:|
|rRNA_bowtie_index|/data/TA_QUIZ_RNA_regulation/data/ATH/index/bowtie1/rRNA/Arabidopsis_thaliana.TAIR10.34.rRNA|
|GTF|/data/TA_QUIZ_RNA_regulation/data/ATH/GTF/Arabidopsis_thaliana.TAIR10.34.gtf|
|STAR_genome_index|/data/TA_QUIZ_RNA_regulation/data/ATH/index/STAR/genome|

## 1a) QC-Trim-QC
### 1a.1) Step one - QC of raw data
**Input：**
`/data/TA_QUIZ_RNA_regulation/data/riboshape_liulab_batch4/Ribo-seq`下
|data type|path|
|:------:|:------:|
|raw data|\~/uvr8\*.fq.gz;\~/wt\*.fq.gz|

**Software/Parameters：**
`fastqc`

|options|function|
|:------:|:------:|
|-q --quiet|禁止在标准输出上显示所有进度，仅报告错误；安静模式|
|-o --outdir|指定输出目录|
|-h --help|选项的详细介绍|
|-t 4|设置线程数|
|--noextract|指定结果文件压缩|

**Output**
`\*.html`：质控报告解读的网页文件。
`\*.zip`：网页文件中画图的数值及图片。
因此`fastqc`结果仅看`html`即可。

### 1a.2) Step two -cut adaptor & trim long read
**input：**
`/data/TA_QUIZ_RNA_regulation/data/riboshape_liulab_batch4/Ribo-seq`下
|data type|path|
|:------:|:------:|
|raw data|\~/uvr8\*.fq.gz;\~/wt\*.fq.gz|

**Software/Parameters：**

`fastp`

|options|function|
|:------:|:------:|
|-i|read1 输入文件名|
|-o|read1 输出文件名|
|-I|read2 输入文件名|
|-O|read2 输出文件名|
|--thread=4|指定线程数为4|
|-l --length_required 15|短于指定长度的reads将被丢弃，默认为15，即长度<15的read被去掉|
|--length_limit|设置read最大长度，默认为0，即没有最大长度限制。|
|-j --json filename|设置输出的json格式的质控结果文件名，不设置则默认json文件名为fastp.json|
|-h --html filename|设置输出html格式的质控结果文件名，不设置则默认html 文件名为fastp.html|

`usage: fastp -i <in1> -o <out1> [-I <in1> -O <out2>] [options...]`

**Output :**
eg:`uvr8no1.clean.fq.gz` fastp质控后的数据。
`uvr8no1.html`/`uvr8no1.json`:质控结果报告文件。

## 1.b) Clean rRNA reads
使用`bowtie`将上一步中得到的`fastp`质控后的`*.clean.[1/2].fastq.gz`文件比对到`rRNA index`上，除去这部分比对到`rRNA index`的部分，从而得到不含`rRNA reads`的文件`.fastq`文件

**Input :**
1.a.2)操作结束后的`\*.clean.fastq.gz`

**Software/Parameters :**

`bowtie` 可以Clean rRNA reads得到不含rRNA reads的`.fastq`文件

|options with Parameter Setting|function|
|:------:|:------:|
|-n N|代表在高保真区内的错配不能超过N个，这里设置为0，即不允许任何错配|
|-a|在允许的错配碱基数量下的全部可能的比对情况都输出|
|-m1 --best -strata |只报告最好的比对情况的那一个输出|
|-p 4|指定使用4个线程|
|-l L|代表序列高保真区的长度，最短不能少于5，这里我们设置为L=15|
|--un \*rm_rRNA.fq|输出不能map到指定基因组上的reads，fasta格式。即我们所需要的去除rRNA后的文件|
|-S|输出文件名，格式为.sam|
|-1/-2|对于PE测序数据，-1指定输入的read1文件，-2指定输入的read2文件|

**reference：**
```bash
bowtie -n 0 -y -a --norc --best --strata -S -p 4 -l 15 \
--un ~/*.rm_rRNA.fq \
~/Arabidopsis_thaliana.TAIR10.34.rRNA \
-1 ~/*.clean.fastq.gz \
-2 ~/*.clean.fastq.gz \
~/*.aligned_rRNA.txt
```


**Output**
不含rRNA reads的`.fastq`文件`*rm_rRNA.fq`,位于fastq文件夹，map到rRNA index上的`*aligned_rRNA.txt`。
`gzip`压缩`*rm_rRNA.fq`文件，得到`*rm_rRNA.fq.gz`文件。

## 1.c) mapping
我们使用`STAR`软件对1.b)操作后得到的`*rm_rRNA.fq.gz`进行比对。
### 1.c.1) Generating genome index
由于PartI.RNA-seq analysis中我们已经创建了参考基因组目录。这里直接使用已经创建好的`genome directory`即可，不做详细描述。

### 1.c.2) Running mapping jobs
**input**
`trimmed`&`remove rRNA`后的`*rm_rRNA.fq.gz`

**Software/Parameters :**

`STAR`

|options with Parameter Setting|function|
|:------:|:------:|
|--runThreadN Number_of_Threads|set number of threads|
|--limitBAMsortRAM 20000000000|maximum available RAM for sorting BAM|
|--outFilterType BySJout|reduces the number of "spurious" junctions|
|--outFilterMismatchNmax 10|alignment will be output only if it has no more mismatches than this value|
|--genomeDir /path/to/genomeDir|path to genom directory|
|--readFilesIn /path/to/read1 /path/to/read2|path to input files|
|--readFilesCommand 'zcat'|If the read files are compressed, use the --readFilesCommand UncompressionCommand option,for gzipped files (.gz) use --readFilesCommand zcat |
|--outFileNamePrefix|output files name prefix|
|--outSAMtype BAM Unsorted|output unsorted Aligned.out.bam file|
|--quantMode TranscriptomeSAM GeneCounts|With --quantMode TranscriptomeSAM option STAR will output alignments translated into transcript coordinates in the Aligned.toTranscriptome.out.bam file|
|--outSAMattributes All|The SAM attributes can be specified by the user using --outSAMattributes|
|--outSAMstrandField intronMotif|For unstranded RNA-seq data, Cufflinks/Cuffdiff require spliced alignments with XS strand attribute, which STAR will generate with --outSAMstrandField intronMotif option|
|--outBAMcompression 6|int:-1 to 10 BAM compression level|
|--outReadsUnmapped Fastx| output of unmapped and partially mapped  reads in separate file;Fastx:output in separate fasta/fastq files, Unmapped.out.mate1/2|

**reference：**
以`wtnouvb1`样本为例：
```bash
STAR \
--runThreadN 4 \
--outFilterType BySJout \
--outFilterMismatchNmax 2 \
--outFilterMultimapNmax 1 \
--genomeDir genomedir \
--readFilesIn ~/wtnouvb1.rm_rRNA.fq.gz \
--readFilesCommand 'zcat' \
--outFileNamePrefix wtnouvb1. \
--outSAMtype BAM SortedByCoordinate \
--quantMode TranscriptomeSAM GeneCounts \
--outSAMattributes All \
--outSAMattrRGline ID:1 LB:ribo_seq PL:ILLUMINA SM:wtnouvb1 \
--outBAMcompression 6 \
--outReadsUnmapped Fastx
```

- outFilterMultimapNmax=1，是因为不考虑Multimap reads
注意部分参数的选择与`differential expression`的选择不同，为什么不同，请大家自行思考。

**output：**
`.Aligned.out.sorted` :比对结果文件
`.toTranscriptome.out.sorted`：比对到转录本上的reads组成的文件

### 1.c.2) samtools sort &index
这里我们使用`samtools sort -T`按TAG值排序，随后使用`samtools index`对排序后的序列建立索引。输出为`bai`文件。

**output :**
- `.Aligned.sortedByCoord.out.bam`
- `.Aligned.toTranscriptome.out.sorted.bam`
- `.Aligned.sortedByCoord.out.bam.bai`
- `.Aligned.toTranscriptome.out.sorted.bam.bai`

## 1.d) Reads counts
这里我们使用`mapping`后并且经过`samtools`排序索引后的比对结果文件作为输入。使用`htseq-count`软件来计算reads数。
**input：**
`.Aligned.sortedByCoord.out.bam`

**Software/Parameters :**
`htseq-count`

**Software/Parameters :**

`htseq-count`

|options with Parameter Setting|function|
|:------:|:------:|
|-f bam|设置输入文件的格式|
|-s no|设置是否是链特异性测序|
|-i gene_id|设置feature ID是由gtf/gff文件第9列哪一个标签决定的|
|-t CDS|设置指定的feature进行表达量计算|
|-m union|设置表达量计算模式，可有的参数值有union, intersection-strict and intersection-nonempty。真核生物推荐使用union模式|

**reference：**

``` bash
htseq-count \
-f bam \
-s no \
-i gene_id \
-t CDS \
-m union \
~/*.Aligned.sortedByCoord.out.bam \
~/Arabidopsis_thaliana.TAIR10.34.gtf > ~/.read_count.HTSeq.txt
```
注意，由于Ribo-seq测序数据，仅测到被核糖体保护的片段，因此这里我们仅设置`CDS`模式。

得到的`read_count.HTSeq.txt`中包含着以`__`开头的注释信息。因此需要进行以下处理。
```
# 提取其中的__开头的信息。
grep __ ~/*.read_count.HTSeq.txt > ~/*.read_count.HTSeq.txt.summary

# 除去其中的__开头的信息
sed -i '/^__/d' ~/*.read_count.HTSeq.txt
```

**output：**
`read_count.HTSeq.txt`：包含gene_id与reads counts的信息。
`read_count.HTSeq.txt.summary`：汇总信息。

# 2）周期性和ORF分析
## 2. a) 获得周期性分析报告
Ribosome profiling data with a good quality tend to have a good 3-nt periodicity.
这里我们使用`RiboCode`中的`metaplots`来获得`3-nt periodicity`分析报告。

**input：**
`.Aligned.toTranscriptome.out.bam`：1.c.2)中`mapping`后并且经过`samtools`排序索引后的比对结果文件
`/data/TA_QUIZ_RNA_regulation/data/Ribocode/RiboCode_annot`:由`RiboCode`产生的一些注释文件目录。

**Software/Parameters：**
`metaplots`

|options with Parameter Setting|function|
|:------:|:------:|
|-a|设置RiboCode注释文件夹路径|
|-r|设置输入的sam文件路径|
|-o|设置输出路径|
|-m 26|设置最小值|
|-M 34|设置最大值|

**reference：**

```
metaplots \
-a /data/TA_QUIZ_RNA_regulation/data/Ribocode/RiboCode_annot \
-r ~/*.Aligned.toTranscriptome.out.sorted.bam \
-o ~/RiboCode/metaplot \
-m 26 \
-M 34 \
-s no \
-pv1 0.05 \
-pv2 0.05
```

**output：**
`~/RiboCode/metaplot`：输出目录，目录中包含`3-nt periodicity`分析的pdf文件。
`*.Aligned.toTranscriptome.out.sorted.pdf`：pdf文件
`*`

## 2. b) 获得ORF区域结果
**input：**
`.Aligned.toTranscriptome.out.bam`：1.c.2)中`mapping`后并且经过`samtools`排序索引后的比对结果文件
`/data/TA_QUIZ_RNA_regulation/data/Ribocode/RiboCode_annot`:由`RiboCode`产生的一些注释文件目录。
`*._pre_config.txt`:产生的配置文件

**Software/Parameters：**
`RiboCode`
`RiboCode -a <RiboCode_annot> -c <config.txt> -l no -g -o <RiboCode_ORFs_result>`

|options with Parameter Setting|function|
|:------:|:------:|
|-a|设置RiboCode注释文件夹路径|
|-c|设置配置文件夹路径|
|-g|获得预测ORF的"gtf"格式文件|
|-o|设置输出文件夹路径|

**reference：**
以`wtnouvb1`样本为例子
```bash
RiboCode \
-a /data/TA_QUIZ_RNA_regulation/data/Ribocode/RiboCode_annot \
-c ~/RiboCode/wtnouvb1._pre_config.txt \
-l no \
-g -o ~/RiboCode/wtnouvb1
```

**output：**
`~/RiboCode/*_collapsed.txt`：ORF结果

# 3）Differential translation efficiency
我们利用`Xtail`基于`Riboseq`的`count matrix`和`RNAseq`的`count matrix`计算差异翻译效率。为了方便对比`RNAseq`只统计`CDS`区域的`reads`。
## 3.a) RNA-seq 汇总
`RNA-seq`

**input：**
`differential expression`的2.a)`Reads counts`步骤产生的统计CDS文件，`.featurecounts.txt`

同计算`counts matrix`的方法相同，合并['CD1_1','CD1_2','CD1_3','CD0_1','CD0_2','CD0_3','UD1_1','UD1_2','UD1_3','UD0_1','UD0_2','UD0_3']12个样本的`Gene id`和`counts`信息。得到`count_CDS.txt`。

这里建议使用 `python`的`pandas`包来提取文件信息。
`count_CDS.txt`的每行为一个基因，每列为一个样本，矩阵中间的数据为表达值。
示例：
![count_CDS.png](https://github.com/zyz-hust/RNA-Structure-Mediate-regulation/blob/11f6ed84ec5d823fbb6a6dcab14d5f87fe9eb707/Images/count_CDS.png)

## 3.b) Ribo-seq 汇总
`Ribo-seq`

**input：**
`Ribo-seq`的1.d)的`Reads counts`步骤产生的统计CDS文件，`read_count.HTSeq.txt`

与上述相同,合并['wtnouvb1','wtnouvb2','wtnouvb3','wtuvb1','wtuvb2','wtuvb3','uvr8no1','uvr8no2','uvr8no3','uvr8uvb1','uvr8uvb2','uvr8uvb3']12个样本的`Gene id`和`counts`信息。得到`WT_count.txt`、`UVR_count.txt`。

示例：

![WT_count.png](https://github.com/zyz-hust/RNA-Structure-Mediate-regulation/blob/11f6ed84ec5d823fbb6a6dcab14d5f87fe9eb707/Images/WT_count.png)

## 3.c) Xtail 分析
使用R包`xtail`进行差异翻译分析。
### 3.c.1) WT样本分析
**input：**
`WT_count.txt`、`count_CDS.txt`

**reference：**
```
# 读取riboseq和mrna的reads counts matrix。
ribo <- read.table('~/wt_count.txt',header=T,quote='',check.names=F, sep='\t',row.names=1)
mrna<- read.table('~/count_CDS.txt',header=T,quote='',check.names=F, sep='\t',row.names=1)

# 样本信息
condition <- c("control","control","treat","treat","treat")

# 获得xtail分析结果，按照q-value排序，选择输出log2FCs和log2R
results <- xtail(mrna,ribo,condition,minMeanCount=1,bins=10000)
results_tab <-resultsTable(results,sort.by="pvalue.adjust",log2FCs=TRUE,log2Rs=TRUE)

# 输出xtail分析结果
write.table(results_tab,"~/wt.uvb-vs-nouvb.TE_new.xls",quote=F,sep="\t")

```

**output：**

![wt.uvb-vs-nouvb.TE_new.png](https://github.com/zyz-hust/RNA-Structure-Mediate-regulation/blob/11f6ed84ec5d823fbb6a6dcab14d5f87fe9eb707/Images/wt.uvb-vs-nouvb.TE_new.png)

### 3.c.2) UVR8样本分析
**input：**
`UVR_count.txt`、`count_CDS.txt`

**reference：**
```
# 读取riboseq和mrna的reads counts matrix。
ribo <- read.table('~/uvr8_count.txt',header=T,quote='',check.names=F, sep='\t',row.names=1)
mrna<- read.table('~/count_CDS.txt',header=T,quote='',check.names=F, sep='\t',row.names=1)

# 样本信息
condition <- c("control","control","treat","treat","treat")

# 获得xtail分析结果，按照q-value排序，选择输出log2FCs和log2R
results <- xtail(mrna,ribo,condition,minMeanCount=1,bins=10000)
results_tab <-resultsTable(results,sort.by="pvalue.adjust",log2FCs=TRUE,log2Rs=TRUE)

# 输出xtail分析结果
write.table(results_tab,"~/uvr8.uvb-vs-nouvb.TE_new.xls",quote=F,sep="\t")

```

**output：**

![uvr8.uvb-vs-nouvb_new.TE.png](https://github.com/zyz-hust/RNA-Structure-Mediate-regulation/blob/11f6ed84ec5d823fbb6a6dcab14d5f87fe9eb707/Images/uvr8.uvb-vs-nouvb_new.TE.png)

包含所有基因结果，未进行差异过滤。



