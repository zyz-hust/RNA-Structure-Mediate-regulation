# 0）数据说明
SHAPE-MaP测序结果中不加NAI的样本可以作为RNA-seq样本，即`/data/zhaoyizi/TA_QUIZ_RNA_regulation/data/riboshape_liulab_batch4`路径下，`WT`/`UVR8`文件夹下`control`下的以`UD`/`CD`开头的`.fastq.gz`文件。共12套。

`P-cluster`上使用队列提交任务的方法可参考supplement中的2）部分
# 1）Mapping

## 1a) QC-Trim-QC
这步操作主要目的有两个，一个是检查数据的质量，另一个是减掉接头序列。

### 1a.1) Step one - QC of raw data
**Input :**

|data type|path|
|:------:|:------:|
|raw data|\~/UD\*.fastq.gz;\~/CD\*.fastq.gz|

**Software/Parameters：**

`fastqc`

|options|function|
|:------:|:------:|
|-q --quiet|禁止在标准输出上显示所有进度，仅报告错误；安静模式|
|-o --outdir|指定输出目录|
|-h --help|选项的详细介绍|
|--noextract|指定结果文件压缩|



**Output :**
`\*.html`：质控报告解读的网页文件。
`\*.zip`：网页文件中画图的数值及图片。
因此`fastqc`结果仅看`html`即可。

### 1a.2) Step two -cut adaptor & trim long read
这里使用`fastp`默认参数，对数据自动进行全方位质控。
**Input :**

|data type|path|
|:------:|:------:|
|raw data|\~/UD\*.fastq.gz;\~/CD\*.fastq.gz|

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
eg:`CD1_1.clean.1.fastq.gz`:fastp质控后数据
`CD1_1.html/CD1_1.json`：质控结果报告文件

## 1.b) Clean rRNA reads
使用`bowtie`将上一步中得到的`fastp`质控后的`*.clean.[1/2].fastq.gz`文件比对到`rRNA index`上，除去这部分比对到`rRNA index`的部分，从而得到不含`rRNA reads`的文件`.fastq`文件

**Input :**
1.a.2)操作结束后的`\*.clean.[1/2].fastq.gz`

**Software/Parameters :**
`bowtie`可以Clean rRNA reads得到不含rRNA reads的`.fastq`文件

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

**reference :**
```bash
bowtie -n 0 -y -a --norc --best --strata -S -p 4 -l 15 \
--un  ~/*.rm_rRNA.fq \
 ~/Arabidopsis_thaliana.TAIR10.34.rRNA \
-1  ~/*.clean.1.fastq.gz \
-2  ~/*.clean.2.fastq.gz \
 ~/*.aligned_rRNA.txt
```

**Output**
不含rRNA reads的`.fastq`文件`*rm_rRNA_[1/2].fq`,位于fastq文件夹，map到rRNA index上的`*aligned_rRNA.txt`

从`fastqc`、`trimmed`到`remove_rRNA`,差异表达分析与差异剪接分析使用的数据与方法完全一致。在`Mapping`过程中有部分参数选择不同，差异表达分析比差异剪接分析多调用一个`--outFilterType BySJout`参数。下面仅介绍`differential expression`的mapping过程。
## 1.c) mapping 
我们使用`STAR`软件对1.b)操作后得到的`*rm_rRNA_[1/2].fq`进行比对。
### 1.c.1) Generating genome index
这里我们采用的参考基因组序列文件`Arabidopsis_thaliana.TAIR10.dna.toplevel.fa`与参考基因组注释文件`Arabidopsis_thaliana.TAIR10.34.gtf`（文件路径见数据介绍部分）来构建`genome index`。

**Software/Parameters：**

`STAR`

|options with Parameter Setting|function|
|:------:|:------:|
|--runThreadN Number_of_Threads|设置调用的线程数|
|--runMode genomeGenerate|指示STAR运行构建genome index的功能|
|--genomeDir /path/to/genomeDir|指定genome index输出目录|
|--genomeFastaFiles /path/to/genome/fasta|指定参考基因组序列文件|
|-sjdbGTFfile /path/to/annotations.gtf|指定参考基因组注释文件|

**reference :**
```bash
STAR \
--runMode genomeGenerate \
--runThreadN 4 \
--genomeDir genomedir \
--genomeFastaFiles Arabidopsis_thaliana.TAIR10.dna.toplevel.fa \
--sjdbGTFfile Arabidopsis_thaliana.TAIR10.34.gtf
```

**output：**
`/path/to/genomeDir`—`genome index`

### 1.c.2) Running mapping jobs
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
```bash
STAR \
--runThreadN 8 \
--limitBAMsortRAM 20000000000 \
--outFilterType BySJout \
--outFilterMismatchNmax 10  \
--genomeDir genomedir \
--readFilesIn .rm_rRNA_1.fq.gz \
.rm_rRNA_2.fq.gz \
--readFilesCommand 'zcat' \
--outFileNamePrefix  name_prefix \
--outSAMtype BAM Unsorted \
--quantMode TranscriptomeSAM GeneCounts \
--outSAMattributes All  \
--outSAMstrandField intronMotif \
--outBAMcompression 6 \
--outReadsUnmapped Fastx
```
**output :**
`.Aligned.out.sorted` :比对结果文件
`.toTranscriptome.out.sorted`：比对到转录本上的reads组成的文件

### 1.c.2)samtools sort &index
为了防止`STAR`在`sort`过程中，内存溢出错误，这里我们设置`limitBAMsortRAM 20000000000;outSAMtype BAM Unsorted`。不在`STAR`过程中进行排序，而是利用samtools单独进行排序。
这里我们使用`samtools sort -T`按TAG值排序，随后使用`samtools index`对排序后的序列建立索引。输出为`bai`文件。

**output :**
- `.Aligned.sortedByCoord.out.bam`
- `.Aligned.toTranscriptome.out.sorted.bam`
- `.Aligned.sortedByCoord.out.bam.bai`
- `.Aligned.toTranscriptome.out.sorted.bam.bai`

注：`differential_expression`与`differential_splicing`的`mapping`过程仅在`--outFilterType BySJout`命令选择上不同。


# 2）Differential expression
## 2.a) Reads counts
这里我们使用`mapping`后并且经过`samtools`排序索引后的比对结果文件作为输入。使用`featureCounts`软件来计算reads数。

**input :**
`.Aligned.sortedByCoord.out.bam`

**Software/Parameters :**
`featureCounts`

|options with Parameter Setting|function|
|:------:|:------:|
|-T 8|指定线程数，这里设置为8|
|-s 0|链特异性的read counting。默认值为0，0代表(非链)，1代表(正向链)，2代表(反向链)。|
|-p |针对PE数据，会统计fragement而不统计read|
|-t CDS|设置feature-type ，-t指定的必须是gtf中有的feature，同时read只有落到这些feature上才会被统计到，默认是“exon”|
|-g gene_id|指定注释文件中的属性信息|
|-a|注释文件，GTF或GFF文件，通过这个文件区分哪些区域为外显子。|
|-o|指定输出|

由于数据采用的建库方法是非链特异性建库方法，因此这里使用`-s 0`模式。读者也可以使用`RSeQC`的`infer_experiment.py`自行判断一下数据是否是非链特性数据。
我们分别选用不同`feature-type`模式`CDS`/`exon`来分别计算`.featurecounts.txt`和`.featurecounts.all.txt`

** output : **
`.featurecounts.all.txt`：统计exon，用于差异表达分析
`.featurecounts.txt`：统计CDS，用于差异分析分析(与Riboseq数据统计方法统一)。


## 2.b) counts matrix
这里使用`python`脚本，将2.a)中的计算的`.featurecounts.all.txt`整理成`counts matrix`。提取其中的`Geneid`(第一行)、`counts`(最后一行)等信息，将['CD1_1','CD1_2','CD1_3','CD0_1','CD0_2','CD0_3','UD1_1','UD1_2','UD1_3','UD0_1','UD0_2','UD0_3']12个样本的信息整合在一起，得到`count_all.txt`。
`count_all.txt`的每行为一个基因，每列为一个样本，矩阵中间的数据为表达值。

这里建议使用 `python`的`pandas`包来提取文件信息。
示例：
![10.2.1.count_matrix.png](../../.gitbook/assets/10.2.1.count_matrix.png)

## 2.c) Differential expression analysis
这里分别使用`edgeR`和`DESeq2`来进行差异表达分析。

### 2.c.1) DESeq2 进行差异表达分析
**input:**
`count_all.txt`：2.b)中产生的表达矩阵。
分别提取其中的`WT`和`UVR8`分别进行分析。

**总的来说分为三步：构建dds矩阵，标准化，以及进行差异分析**
#### (1) 构建dds矩阵
```R
 dds <- DESeqDataSetFromMatrix(countData = cts, colData = coldata, design= ~ batch + condition) 
 #~在R里面用于构建公式对象，~左边为因变量，右边为自变量。
```
> countData ：表达矩阵，通过read count计算后并融合生成的矩阵，行为各个基因，列为各个样本，中间为计算reads或fragment得到的整数。
> colData ：样本信息矩阵，一个数据框，第一列是样本名称，第二列是样本的处理情况，即condition，condition的类型是一个factor。
>  design：差异比较矩阵，差异比较矩阵就是告诉差异分析函数是要从要分析哪些变量间的差异，简单来说就是说明哪些是对照，哪些是处理。

#### (2) 标准化DESeq()
对原始`dds`进行`normalize`，`DESeq`包含三步，`estimation of size factors（estimateSizeFactors)`，` estimation of dispersion（estimateDispersons)`， `Negative Binomial GLM fitting and Wald statistics（nbinomWaldTest）`，可以分布运行，也可用一步到位，最后返回 results可用的`DESeqDataSet`对象。

**分步:**
1）归一化系数sizeFactor
```R
dds <- estimateSizeFactors(dds) 
```
2) 估计基因的离散程度estimateDispersions
```R
dds<- estimateDispersions(dds)
```
3) 统计检验，差异分析
```R
dds <- nbinomWaldTest(dds) 
```

**或者一步就位：**
```R

dds <- DESeq(dds) 

```
DESeq2假定基因的表达量符合负二项分布，有两个关键参数，总体均值和离散程度α值。α值衡量均值和方差之间的关系。

### 3. 获得分析结果，差异分析
```R
# 获取分析结果
res<-results(dds)
# 对结果按照padj值大小排序
res <- res[order(res$padj),]
# 筛选出其中 padj<0.05 & log2FoldChange绝对值大于1的基因。
diff_gene_deseq_col0 <- subset(res,padj<0.05&(log2FoldChange>1|log2FoldChange< -1))

```

* 1） 默认使用样本信息的最后一个因子与第一个因子进行比较
* 2）返回一个数据库res，包含6列，baseMean、log2FC、lfcSE、stat、pvalue、padj
* 3）baseMean表示所有样本经过归一化系数矫正的read counts=(counts/sizeFactor)的均值。
* 4）log2Foldchange表示该基因的表达发生了多大的变化，是估计的效应大小effect size。变化倍数=2^log2Foldchange。==DESeq2在差异分析的过程中已经考虑了样本本身的差异，其最终提供的log2FC只包含了分组间的差异==
* 5) lfcSE 是对log2Foldchange估计的标准误差估计，效应大小估计具有不确定性。
* 6）stat是Wald统计量，是由log2Foldchange/标准差所得。
* 7）pvalue和padj分别代表原始p值及经过校正后的p值。

**output:**
`wt_gene_list.txt`:差异表达的基因列表
`wt_rawdata.csv`:DESeq2产生的原始结果文件。
**results:**
![10.2.2.DESeq2_result.png](![](../../.gitbook/assets/10.2.2.DESeq2_result.png))

### 2.c.2) edgeR 进行差异表达分析
**input :**
`count_all.txt`：2.b)中产生的表达矩阵。
分别提取其中的`WT`和`UVR8`分别进行分析。

同`DESeq2`分析一致，`edgeR`分析也分大致分三步

#### (1) 构建DGEList 对象
```R
dgListGroups <- c(rep("Control",3),rep("Treat",3))
dgList <- DGEList(counts=countData_col0, genes=rownames(countData_col0),group=factor(dgListGroups))
```
> countData ：表达矩阵
> group ：分组信息数据
> genes ：基因名称

#### (2) 过滤掉low counts数据
可以在构建DGEList之前就已经构建完毕，或者可以根据以下方法进行构建。
```R
# 筛掉raw_count在6个样本中总数少于100的基因。
col0_raw_count <- col0_raw_count[rowSums(col0_raw_count)>100,]

```

#### (3) 使用TMM算法对DGEList标准化
```R
dgList <- calcNormFactors(dgList, method="TMM")
countsPerMillion <- cpm(dgList, normalized.lib.sizes=TRUE)

design.mat <- model.matrix(~0 + dgList$sample$group)
colnames(design.mat) <- levels(dgList$sample$group)

# 估计common离散度
d2 <- estimateGLMCommonDisp(dgList, design=design.mat)
# 估计trended离散度
d2 <- estimateGLMTrendedDisp(d2, design=design.mat)
# 估计tagwise离散度
d2 <- estimateGLMTagwiseDisp(d2, design=design.mat)

# glmFit 和 glmLRT 函数是配对使用的，用于 likelihood ratio test (似然比检验)
fit <- glmFit(d2, design.mat)
lrt <- glmLRT(fit,contrast=c(-1,1))

# 计算差异基因得到结果
edgeR_result <- topTags(lrt,n=nrow(dgList))
```

**results：**
![10.2.3.edgeR_result.png](../../.gitbook/assets/10.2.3.edgeR_result.png)
每列依次是`genes`、`logFC`、`logCPM`、`LR`、`PValue`和`FDR`。可根据`logFC`和`FDR`的值对结果进行筛选。

最终的结果如下：
```
5.DESeq2/
├── uvr8
│   ├── uvr8_gene_list.txt
│   └── uvr8_rawdata.csv
└── wt
    ├── wt_gene_list.txt
    └── wt_rawdata.csv

6.edgeR/
├── uvr8
│   └── uvr8_rawdata.csv
└── wt
    └── wt_rawdata.csv
```

## 2.d) GO pathway analysis
通过差异表达分析得到的`gene_list`,可以在KEGG、GO等基因注释数据库中进行通路分析。这里不做详细讨论，大家自行分析。


# 3）Differential splicing
## 3.a) 利用rMATS计算splicing events 
这里我们使用`mapping`后并且经过`samtools`排序索引后的比对结果文件进行后续分析。注意差异表达和差异分析的`mapping`不同。使用`rMATS`来计算`splicing event`

**input：**
`wt.noUVB.txt`、`wt.UVB.txt`
`uvr8.noUVB.txt`、`uvr8.UVB.txt`
mapping后的bam文件路径，文件内以逗号分隔重复样本的bam文件名。
eg：wt.noUVB.txt:
`~/CD1_1Aligned.sortedByCoord.out.bam,~/CD1_2Aligned.sortedByCoord.out.bam,~/CD1_3Aligned.sortedByCoord.out.bam`

**Software/Parameters：**

`rmats.py`

|options with Parameter Setting|function|
|:------:|:------:|
|--b1 b1.txt|输入sample1的txt格式的文件。文件内以逗号分隔重复样本的bam文件名|
|--b2 b2.txt|同--b1|
|--gtf gtf_file|需要输入的gtf文件路径|
|--od outDir|输出文件的文件夹路径|
|--cstat |设置splicing difference的阈值，这里设置为0.0001|
|--readLength |给定读段长度，这里设置为151|
|--variable-read-length|允许除设置的readLength之外的其他长度的reads|
|--tmp tmpDir|设置临时文件夹路径|

**reference：**

```bash
python2 rmats.py \
--b1 ~/wt.noUVB.txt \
--b2 ~/wt.UVB.txt \
--gtf ~/Arabidopsis_thaliana.TAIR10.34.gtf \
--od ~/wt.UVB-vs-noUVB \
-t paired \
--readLength 151 \ 
--cstat 0.0001 \
--tmp wt.UVB-vs-noUVB/tmp \
--nthread 4 \
--variable-read-length
```

**output：**
rMATS的结果文件是以各个可变剪切事件的分布。
XX 指代SE\RI\MXE\A5SS\A3SS 五项可变剪切时间。
* 1. fromGTF.XX.txt 系列： 直接从GTF文件和数据文件读出的结果
* 2. fromGTF.novelEvents.XX.txt : 从数据文件发现的新的可变剪切事件

* 3. XX.MATS.JC.txt 和XX.MATS.JECE.txt,是JC.raw.input.XX.txt 和 JCEC.raw.input.XX.txt 经过统计模型分析后的结果，多了P值和FDR值。
* 4. JC和JCEC的区别在于前者考虑跨越剪切位点的reads，而后者不仅考虑前者的reads还考虑到只比对到第一张图的条纹区(没有跨越剪切位点的reads)

**XX.MATS.JC.txt中包含的信息比较多，以SE.MATS.JC.txt为例：**
* 1. 1-5列分别为：ID、GeneID、geneSymbol、chr、strand
* 2. 6-11列分别为外显子位置信息：分别为exonStart_0base，exonEnd，upstreamES，upstreamEE，downstreamES，downstreamEE。
![10.2.4.skipped_exon.png](../../.gitbook/assets/10.2.4.skipped_exon.png)
* 3. 13-16列  展示两组样品在IJC(inclusion junction) 和SJC(skipping junction counts)下的counts数，重复样本的结果以逗号分隔：列名分别为IJC_SAMPLE_1，SJC_SAMPLE_1，IJC_SAMPLE_2，SJC_SAMPLE_2。

![10.2.5.skipped_exon_AS.png](../../.gitbook/assets/10.2.5.skipped_exon_AS.png)

* 4. lncFormLen和SkipFormLen分别是inclusion form和skipping form的有效长度值,虽然有计算公式，还是要根据reads跨越时的具体情况来定。
* 5. IncLevel 可被认作为exon inclusion level(φ),是exon inclusion isoform在总(Exon inclusion isoform +exon skipping isoform 所占比例)
* 6. IncLevelDifference则是指两组样本IncLevel的差异，如果一组内多个样本，那么则是各自组的均值之间差值。

## 3.b) splicing events 结果处理
上述通过`rmats.py`计算出来的数据是未经过过滤筛选的。下面我们将对rMATS脚本进行过滤,使用位于`/data/TA_QUIZ_RNA_regulation/data/script/PartI.RNA-se_analysis/differential_splicing`路径下的`splice_sig_psi.py`脚本。
过滤指标`p value <0.05 & psi>0.1`

**input ：**
`~/wt.UVB-vs-noUVB`：rMATS结果的输出文件夹路径

**reference：**
```
# p value<0.05 & psi>0.1
python splice_sig_psi.py PValue 0.05 0.1 wt.UVB-vs-noUVB wt.UVB-vs-noUVB_filtered_psi_0.1
```

**output：**
`uvr8.UVB-vs-noUVB_filtered_psi_0.1` ,` wt.UVB-vs-noUVB_filtered_psi_0.1`

## 3.c) 对过滤后的differential splicing gene进行通路注释
从`uvr8.UVB-vs-noUVB_filtered_psi_0.1` ,` wt.UVB-vs-noUVB_filtered_psi_0.1`中分别提取过滤后的不同剪接事件的差异剪接基因的`xx_gene_list`。从服务器上下载下来。到[DAVID Functional Annotation Tool]
(https://david.ncifcrf.gov/summary.jsp)网站上进行分析
****

### 3.c.1) DAVID Functional Annotation Tool
点击`shortcut to DAVID Tool`->`Function annotation clustering` 左侧出现`Upload Gene List`界面。

![10.2.6.DAVID.png](../../.gitbook/assets/10.2.6.DAVID.png)

按照提示步骤提交文件即可，注：Step2选择`TAIR_ID`拟南芥基因ID，最后`Submit list`得到功能注释结果。


## 3.d) splicing events画图
使用`rmats2sashimiplot`软件对差异剪接事件进行可视化处理

**input：**
`~/wt.UVB-vs-noUVB_filtered_psi_0.1/${event}.MATS.JCEC.txt`：3.c)经过过滤的rMATS结果输出文件夹中各种剪接事件的结果文件。
`~/*Aligned.sortedByCoord.out.bam`：1.c.2)中sort后的mapped reads的bam文件

**Software/Parameters：**

`rmats2sashimiplot`

|options with Parameter Setting|function|
|:------:|:------:|
|-b1 B1.1,B1.2,B1.3|指定bam文件，需要是sort后的，输入的mapping files 能够通过逗号划分不同的group|
|--b2 B2.1,B2.2,B2.3|同上|
|--l1 L1 --l2 L2|指定label|
|-t {SE,A5SS,A3SS,MXE,RI} |指定rMATS事件类型|
|-e EVENTS_FILE	|rMATS输出文件|
|-o OUT_DIR|指定输出目录|
|-exon_s EXON_S|缩小多少倍外显子大小，Default:1|
|-intron_s INTRON_S|缩小多少倍内含子大小|
|-min-counts 0|read count 小于--min-counts 的individual junction 会被忽略。|

**reference：**
```linux
rmats2sashimiplot \
--b1 CD1_1.Aligned.sortedByCoord.out.bam,CD1_2.Aligned.sortedByCoord.out.bam,CD1_3.Aligned.sortedByCoord.out.bam \
--b2 CD0_1.Aligned.sortedByCoord.out.bam,CD0_2.Aligned.sortedByCoord.out.bam,CD0_3.Aligned.sortedByCoord.out.bam \
--l1 Control_wt_noUVB \
--l2 Treat_wt_UVB\
--exon_s 1 \
--intron_s 2 \
--min-counts 0 \
-t $event \
-e $inputDir/${event}.MATS.JCEC.txt \
-o wt.UVB-vs-noUVB_sashimiplot_all/$event
```
其中$event 代表着{SE,A5SS,A3SS,MXE,RI}不同剪接事件。

**output：**
`wt.UVB-vs-noUVB_sashimiplot_all`
`uvr8.UVB-vs-noUVB_sashimiplot_all`
其中每个文件为画图结果，开头ID与{event}.MATS.JCEC.all.sashimiplot.txt文件第一列ID相对应。







