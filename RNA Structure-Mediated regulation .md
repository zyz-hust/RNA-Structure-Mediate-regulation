# 2.RNA Structure-Mediated regulation 
## RNA regulation : RNA Structure-Mediated Post-transcriptional regulations affect translation efficiency

## 0) 背景介绍
RNA的转录后调控对于基因表达来说是至关重要的，单链RNA由于其不稳定性，总是倾向于形成有稳定的碱基互补配对的结构,其中最典型典型的就是茎环结构。不同的位置上的茎环结构的形成会直接影响到表达量，RNA分子的舒展程度可以很明显地影响核糖体沿RNA移动的过程。因此，RNA二级结构在不同条件下的转变最直接影响到的就是翻译过程。在以往的研究中发现，这种影响往往是介导转录调控(可变剪接AS、多聚腺苷酸加尾APA、嵌合体RNA等)进行的。

相较于使用RNA酶选择性切割单链或双链核苷酸的体外酶促RNA结构分析由于RNA酶体积过大且受限于体外反应不能准确的反映活细胞内的RNA真实的折叠状态，基于化学修饰的RNA结构探测方法由于相对体积更小能够更精细更敏锐地体内检测特定的RNA结构。通过引物延伸分析的选择性2’羟基酰化(SHAPE)与硫酸二甲酯化(DMS)是RNA修饰的两种主要的化学加合物。DMS方法的一个显著局限性就是缺失了将近一半的转录组RNA结构信息。因为DMS仅探测As(腺苷酸)和Cs(胞嘧啶)的结构信息。而缺少Us(尿嘧啶)和Gs(鸟苷酸)的配对状态信息。SHAPE方法由于其高分辨率与普适性，在近期研究中被广泛地用来探测mRNA二级结构调控元件以及其调控功能。

本Quiz依托于Lulab现有的一些研究方法，结合与上海生命科学研究院植物生理生态研究所的刘宏涛老师实验室合作产生的SHAPE-MaP测序数据及Ribo-seq测序数据，测定在UVB光照前后，拟南芥可能产生的各种转录后调控，结构变化以及最终导致的翻译效率的变化。          通过各种已有的工具和统计学分析方法，我们希望将两种条件下的差异剪接，差异表达，结构变化等信息与差异翻译联系起来，希望从中找寻到具有生物学意义的结果，例如建立与UVR8或光照通路密切相关的RNA二级结构动态模型，探究结构展开动力学与翻译效率之间的关系。

## 1) 总体流程图
![RNA regulation protocol2.png](https://github.com/zyz-hust/RNA-Structure-Mediate-regulation/blob/baa37b5eb26fd367ce321c457642f89b9a8cb97e/Images/RNA%20regulation%20protocol2.png)

### Part I. RNA-seq analysis
- 1 完成样本的**Reads Processing、Remove RNA and Mapping**工作，得到Mapped reads(bam)并绘制质量控制相关图。
- 2 完成**differential expression**分析
	- 2.1) 计算RNA-seq reads count matrix
	- 2.2) 利用DESeq2进行差异表达分析
	- 2.3) 利用edgeR进行差异表达分析
- 3 完成**differential splicing**分析(选做)
	- 3.1) 利用rMATS计算splicing events
	- 3.2) 对差异剪接基因进行通路注释并作图
### Part II. Ribo-seq analysis
- 1 完成样本的**Reads Processing、Remove RNA and Mapping**工作，得到Mapped reads(bam)并绘制质量控制相关图。
- 2 ORF分析，计算Ribo-seq的reads count matrix
- 3 差异翻译效率分析
	- 3.1) 利用Xtail基于Ribo-seq的count matrix 和RNA-seq的count matrix计算差异翻译效率。为了方便对比RNAseq只统计CDS区域的reads。
### Part III. SHAPE-seq analysis
- 1 利用Shapemapper 进行数据预处理并计算每个转录本的reactivity
- 2 计算得到结构改变区域
	- 2.1）计算hit level，根据hit level分布确定阈值
	- 2.2) 计算归一化因子，对Reactivity进行预处理和归一化。
	- 2.3）计算滑动窗口Gini index，利用滑动窗口判断每个转录本加光前后的结构改变区域，对连续的滑动窗口进行合并，汇总结构。
### Part IV. comprehensive analysis
- 1 转录本丰度变化和翻译效率(TE)变化程度之间关系
- 2 转录本结构改变程度和翻译效率TE变化程度之间关系
- 3 翻译效率变化基因的motif分析
	- 3.1) TE上/下调的3'UTR和5'UTR的motif富集。
	- 3.2) 富集的motif区域是否发生结构变化。 

## 2) 报告要求
**报告要求**：提交一份完整的工作报告，中英文不限(鼓励英文，可以参考一些发表的文献，如([Pervasive Regulatory Functions of mRNA Structure
Revealed by High-Resolution SHAPE Probing](https://www.cell.com/cell/fulltext/S0092-8674(18)30211-3)),同时提交源代码。请读者使用我们提供的数据，完成上述分析，并回答以下问题。

### Part I.RNA-seq analysis
1) WT/UVR8型中UVB光照前后发生差异表达(DESseq)的基因总共有多少个，其中多少个基因上调，多少个基因下调；其中光通路相关基因PIF3(AT1G09530)在UVB光照前后差异表达的log2FoldChange值为多少。
2) 绘制差异表达基因的GO/KEGG通路分析图(选做)

### Part II.Ribo-seq analysis
3) WT/UVR8型中UVB光照前后发生差异翻译的基因总共有多少，其中多少基因上调，多少基因下调；PIF3(AT1G09530)在UVB光照前后差异翻译log2FC_TE_final值为多少。
4) 绘制差异翻译基因的GO/KEGG通路分析图(选做)

### Part III.SHAPE-seq analysis
5) 绘制转录本随hit level变化趋势图，思考应该选择何值作为hit level阈值，为什么
6) 用于归一化的转录本为(),思考为什么要对样本做归一化
7) WT型UVB光照后PIF3(AT1G09530.3)的平均reactivity值为多少,hit level值为多少。

### Part IV. comprehensive analysis
8) 转录本丰度变化和翻译效率(TE)变化程度之间关系
- 绘制Log2FoldChange(TE)和Log2FoldChange(RNA-seq)的关系图，计算相关系数。
- 说明转录本丰度变化和TE是否存在相关关系，若存在，试解释产生这种相关性的原因。
9) 转录本结构改变程度和翻译效率TE变化程度之间关系
- 绘制TE变化与结构改变关系图
10) 翻译效率变化基因的motif分析
- TE上/下调的3'UTR和5'UTR的motif富集。
- 说明富集的motif区域是否发生结构变化。 

## 3）数据介绍
数据来自于我们合作的上海生命科学研究院生理生态研究所的刘宏涛研究员产生的部分SHAPE-MaP测序数据及Ribo-seq测序数据。实验设计如下
![Arabidopsis UVB protocal.png](https://github.com/zyz-hust/RNA-Structure-Mediate-regulation/blob/a0a18cc17de65e83b8181d09bd5311ffb290e21f/Images/Arabidopsis%20UVB%20protocal.png)

实验使用两种植株，分别是Columbia(Col-0)做为野生型，以及对uvr8基因进行转移DNA插入突变的uvr8突变株。分别设置实验组UVR8和对照组WT。实验组放在313-315nm的窄带UV-B紫外光下处理1h，而对照组不做紫外光处理，在白光下处理1h，共产生四种不同处理状态的植株——WT_UV-(野生型，白光处理)、UVR8_UV-(UVR8突变型，白光处理)、WT_UV+(野生型,UV-B处理)、UVR8_UV+(UVR8突变型，UV-B处理)。

随后对四种不同处理状态的植株进行SHAPE-MaP实验以及不进行化学修饰的对照实验和Ribo-seq实验。SHAPE-MaP 使用NAI作为化学修饰物，而对照实验加入等量的DMSO溶液处理。每组实验均进行三次生物重复。RNA提取使用的为天根试剂盒DP432，测序前的反转录使用takara试剂盒。SHAPE测序及Ribo-seq测序均由测序公司完成。共24套SHAPE-MaP数据以及12套Ribo-seq数据，实验数据如下。
<div align="center">
<table style="text-align:center">
	<caption>表1 SHAPE-MaP数据统计</caption>
   <tr>
      <th>样本名称</th>
      <th>植株</th>
      <th>光照处理</th>
      <th>RNA 标记处理</th>
      <th>样本数</th>
   </tr>
   <tr>
      <td>CD_0</td>
      <td>Col-0</td>
      <td>UV-B</td>
      <td>DMSO</td>
      <td>3</td>
   </tr>
   <tr>
      <td>CD_1</td>
      <td>Col-0</td>
      <td>White light</td>
      <td>DMSO</td>
      <td>3</td>
   </tr>
   <tr>
      <td>CN_0</td>
      <td>Col-0</td>
      <td>UV-B</td>
      <td>NAI</td>
      <td>3</td>
   </tr>
   <tr>
      <td>CN_1</td>
      <td>Col-0</td>
      <td>White light</td>
      <td>NAI</td>
      <td>3</td>
   </tr>
   <tr>
      <td>UD_0</td>
      <td>UVR8</td>
      <td>UV-B</td>
      <td>DMSO</td>
      <td>3</td>
   </tr>
   <tr>
      <td>UD_1</td>
      <td>UVR8</td>
      <td>White light</td>
      <td>DMSO</td>
      <td>3</td>
   </tr>
   <tr>
      <td>UN_0</td>
      <td>UVR8</td>
      <td>UV-B</td>
      <td>NAI</td>
      <td>3</td>
   </tr>
   <tr>
      <td>UN_1</td>
      <td>UVR8</td>
      <td>White light</td>
      <td>NAI</td>
      <td>3</td>
   </tr>
</table>
</div >

<div align="center">
<table style="text-align:center">
	<caption>表2 Ribo-seq数据统计</caption>
   <tr>
      <td>样本名称</td>
      <td>植株</td>
      <td>光照处理</td>
      <td>样本数</td>
   </tr>
   <tr>
      <td>wtuvb</td>
      <td>Col-0</td>
      <td>UV-B</td>
      <td>3</td>
   </tr>
   <tr>
      <td>wtnouvb</td>
      <td>Col-0</td>
      <td>White light</td>
      <td>3</td>
   </tr>
   <tr>
      <td>uvr8uvb</td>
      <td>UVR8</td>
      <td>UV-B</td>
      <td>3</td>
   </tr>
   <tr>
      <td>uvr8no</td>
      <td>UVR8</td>
      <td>White light</td>
      <td>3</td>
   </tr>
</table>
</div >

表1 样本名称中，C代表Col-0植株，U代表UVR8突变株；D代表SHAPE-MaP实验中使用DMSI对照处理，N代表NAI处理；0代表UV-B处理，1代表白光处理。

由于SHAPE-MaP实验中作为对照组的，未加入小分子修饰的组别，与正常进行的RNA-seq得到的测序数据没有明显差异，因此用来做转录后调控分析的RNA-seq数据即为SHAPE-MaP实验中的对照组数据CD/UD。

## 4）数据及软件存放地址
作业所用到的`SHAPE-MaP`数据及`Ribo-seq`数据均存放于`P-cluster`的`/data/TA_QUIZ_RNA_regulation/data/riboshape_liulab_batch4`下
### 4.1）SHAPE-MaP 数据

```shell
SHAPE-MaP/
├── UVR8
│   ├── control
│   │   ├── nouvb
│   │   │   ├── UD1_1.clean.1.fastq.gz
│   │   │   ├── UD1_1.clean.2.fastq.gz
│   │   │   ├── UD1_2.clean.1.fastq.gz
│   │   │   ├── UD1_2.clean.2.fastq.gz
│   │   │   ├── UD1_3.clean.1.fastq.gz
│   │   │   └── UD1_3.clean.2.fastq.gz
│   │   └── uvb
│   │       ├── UD0_1.clean.1.fastq.gz
│   │       ├── UD0_1.clean.2.fastq.gz
│   │       ├── UD0_2.clean.1.fastq.gz
│   │       ├── UD0_2.clean.2.fastq.gz
│   │       ├── UD0_3.clean.1.fastq.gz
│   │       └── UD0_3.clean.2.fastq.gz
│   └── modified
│       ├── nouvb
│       │   ├── UN1_1.clean.1.fastq.gz
│       │   ├── UN1_1.clean.2.fastq.gz
│       │   ├── UN1_2.clean.1.fastq.gz
│       │   ├── UN1_2.clean.2.fastq.gz
│       │   ├── UN1_3.clean.1.fastq.gz
│       │   └── UN1_3.clean.2.fastq.gz
│       └── uvb
│           ├── UN0_1.clean.1.fastq.gz
│           ├── UN0_1.clean.2.fastq.gz
│           ├── UN0_2.clean.1.fastq.gz
│           ├── UN0_2.clean.2.fastq.gz
│           ├── UN0_3.clean.1.fastq.gz
│           └── UN0_3.clean.2.fastq.gz
└── WT
    ├── control
    │   ├── nouvb
    │   │   ├── CD1_1.clean.1.fastq.gz
    │   │   ├── CD1_1.clean.2.fastq.gz
    │   │   ├── CD1_2.clean.1.fastq.gz
    │   │   ├── CD1_2.clean.2.fastq.gz
    │   │   ├── CD1_3.clean.1.fastq.gz
    │   │   └── CD1_3.clean.2.fastq.gz
    │   └── uvb
    │       ├── CD0_1.clean.1.fastq.gz
    │       ├── CD0_1.clean.2.fastq.gz
    │       ├── CD0_2.clean.1.fastq.gz
    │       ├── CD0_2.clean.2.fastq.gz
    │       ├── CD0_3.clean.1.fastq.gz
    │       └── CD0_3.clean.2.fastq.gz
    └── modified
        ├── nouvb
        │   ├── CN1_1.clean.1.fastq.gz
        │   ├── CN1_1.clean.2.fastq.gz
        │   ├── CN1_2.clean.1.fastq.gz
        │   ├── CN1_2.clean.2.fastq.gz
        │   ├── CN1_3.clean.1.fastq.gz
        │   └── CN1_3.clean.2.fastq.gz
        └── uvb
            ├── CN0_1.clean.1.fastq.gz
            ├── CN0_1.clean.2.fastq.gz
            ├── CN0_2.clean.1.fastq.gz
            ├── CN0_2.clean.2.fastq.gz
            ├── CN0_3.clean.1.fastq.gz
            └── CN0_3.clean.2.fastq.gz

```

### 4.2）Ribo-seq 数据
```shell
Ribo-seq/
├── uvr8no1.fq.gz
├── uvr8no2.fq.gz
├── uvr8no3.fq.gz
├── uvr8uvb1.fq.gz
├── uvr8uvb2.fq.gz
├── uvr8uvb3.fq.gz
├── wtnouvb1.fq.gz
├── wtnouvb2.fq.gz
├── wtnouvb3.fq.gz
├── wtuvb1.fq.gz
├── wtuvb2.fq.gz
└── wtuvb3.fq.gz
```

### 4.3）注释文件
注释文件均位于`/data/TA_QUIZ_RNA_regulation/data/ATH`下，以下文件路径中`~/`均代指`/data/TA_QUIZ_RNA_regulation/data/ATH`

|data|path|
|:------:|:------:|
|rRNA bowtie index|~/GTF/Arabidopsis_thaliana.TAIR10.34.gtf|
|GTF|~/index/bowtie1/rRNA/Arabidopsis_thaliana.TAIR10.34.rRNA|
|STAR genome index|~/index/STAR/genome|

所用到的软件可以自行下载，也可以使用已经配置好的软件，其位于`P-cluster`的`/data/zhaoyizi/software/anaconda3/envs/Riboshape/bin`
