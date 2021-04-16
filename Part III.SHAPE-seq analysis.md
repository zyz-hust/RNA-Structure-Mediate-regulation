# 0）数据说明
SHAPE-MaP测序数据共24套。位于`/data/zhaoyizi/TA_QUIZ_RNA_regulation/data/riboshape_liulab_batch4`路径下，`WT`/`UVR8`文件夹下,`control`/`motified`下的以`CN`/`CD`或`UN`/`UD`开头的`.fastq.gz`文件。
本章我们分三步进行，(1)利用shapemapper软件计算每个转录本的reactivity文件；(2) 计算结构改变区域，包括hit level 计算，数据预处理，滑动窗口结构计算和结构改变区域合并；(3) motif 分析，包括序列motif，已知结构motif，从头发现的结构motif。

# 1）Shapemapper

这里使用`shapemapper`软件对数据进行预处理操作。我们选择计算==每个基因最长的转录本==，
由于`shapemapper`的运行速度非常慢，这里我们仅选择部分example文件供大家练习`Shapemapper`的使用，本节最后会提供使用`Shapemapper`处理好的24套数据的每个转录本的`reactivity`文件。

**input：**
这里，我们将`SHAPE-MaP`测序数据中的`CD1_1.clean.1/2.fastq.gz`,`CN1_1.clean.1/2.fastq.gz`使用`SEQKIT`分割20份，各取其中一份得到的`CD1_1.clean.part1.R1.fastq.gz`/`CD1_1.clean.part1.R2.fastq.gz`和`CN1_1.clean.part1.R1.fastq.gz`/`CN1_1.clean.part1.R2.fastq.gz`作为`Shapemapper`的输入，文件位于`/data/TA_QUIZ_RNA_regulation/data/shapemapper_test`的`control`/`modified`文件夹下。

为了进一步降低运行时间，我们参考转录本文件也进行了简化，仅取三条转录本作练习使用，文件位于`/data/TA_QUIZ_RNA_regulation/data/shapemapper_test/Arabidopsis_thaliana.TAIR10.34.transcripts_new_2.fa`

**Software/Parameters：**

`shapemapper`

|options|function|
|:------:|:------:|
|--target /path/to/transcript_index|FASTA  file or list of files containing one or more target DNA sequences('T' not 'U')|
|--name "Ath_C1_1"|所有输出文件的前缀，在同一个文件下运行多个Shapemapper时十分推荐。|
|--out|输出文件夹、默认为"shapemapper_out"|
|--temp|临时文件夹，默认为"shapemapper_temp"|
|--min-depth 100| Minimum effective sequencing depth for including data|
|--min-qual-to-count 20|Only count mutations with all basecall quality scores meeting this minimum score|
|--overwrite|覆盖输出和临时文件夹中的现有警告。默认值=False。|
|--modified --folder /path/to/modified|设置modified文件夹路径|
|--untreated --folder /path/to/control|设置control文件夹路径|
|--star-aligner|使用STAR代替Bowtie2进行序列比对|
|--nproc 8|设置线程数|
|--verbose|显示每个执行过程的完整命令，并在发生错误时，显示更多过程输出信息。默认值为=False。|

**reference：**

```bash
shapemapper \
--target ~/Arabidopsis_thaliana.TAIR10.34.transcripts_new_1_1.fa \
--name "Ath_test_C1_1" \
--min-depth 100 \
--min-qual-to-count 20 \
--overwrite \
--modified --folder ~/data/modified/  \
--untreated --folder ~/data/control/ \
--star-aligner \
--nproc 8 \
--verbose
```

**result:**
结果位于输出文件夹`shapemapper_out`中，以`AT1G01010.1`这个转录本为例,其中的`Ath_test_C1_1_transcript_id_profile.txt`中的`Sequence`、`Modified_mutations`、`Modified_effective_depth`、`Untreated_mutations`、`Untreated_effective_depth`列。

Ath_test_C1_1_transcript_id_profile.txt：

![10.4.1.transcript_profile.png](../../.gitbook/assets/10.4.1.transcript_profile.png)


将C1_1的所有转录本的信息整合起来
得到的整合文件`C1_1/final.modified_unmodified_new`如图所示：

![10.4.2.final.modified_unmodified_new.png](../../.gitbook/assets/10.4.2.final.modified_unmodified_new.png)

文件的每列内容为`transcript_id`,`Nucleotide`,`Modified_mutations`,`Modified_effective_depth`,`Untreated_mutations`,`Untreated_effective_depth`。

将`C1_1`、`C1_2`、`C1_3`三个“野生型，非光照”的replicates样本的转录本信息整合到一起,得到`col_nouv/final.modified_unmodified`。同理我们整合得到`col_uv/final.modified_unmodified`;`uvr8_nouv/final.modified_unmodified`;`uvr8_uv/final.modified_unmodified`。

==整合部分不需要大家完成，只需了解是如何得到最终结果文件即可。这里将直接提供整合好的转录本信息文件用于后续分析使用。==文件位于`/data/TA_QUIZ_RNA_regulation/data/riboshape_liulab_batch4/final.modified_umodified`下,`col_nouv`/`col_uv`/`uvr8_nouv`/`uvr8_uv`文件夹下,`final.modified_unmodified`用于后续计算`hit level`。

# 2）Structure changes analysis
## 2.a) 数据及原理
为了方便理解，我们介绍`Reactivity`、`Hit level`、`Gini index` 等几个基础概念。
- Reactivity
R = mutation rate(S)- mutation rate(B)
S corresponds to a SHAPE modified sample, B to untreated.

- Hit level

![10.4.3.hit_level.png](../../.gitbook/assets/10.4.3.hit_level.png)

The hit level metric quantifies the total background-subtracted signal per nucleotide of transcript.

- Gini index

![10.4.4.Gini index.png](../../.gitbook/assets/10.4.4.Gini index.png)

We then applied these metrics to windows containing a total of 50 nucleotides Requiring a minimum number of 100 total reads, number of nan in reactivity less than 10%. 

`Reactivity`、`Hit level` refer to [rnA motif discovery by shAPe and mutational
profiling (shAPe-maP)](https://www.nature.com/articles/nmeth.3029);
`Gini index` refer to [A Stress Response that Monitors and Regulates
mRNA Structure Is Central to Cold Shock Adaptation](https://www.sciencedirect.com/science/article/pii/S1097276518301801?via%3Dihub)

我们以Shapemapper计算得到的整合好的`col_nouv`、`col_uv`、`uvr8_nouv`、`uvr8_uv`中所有转录本结果来计算得到结构改变区域。数据位于`/data/TA_QUIZ_RNA_regulation/data/riboshape_liulab_batch4/final.modified_umodified`下。

我们的计算过程主要包括以下几步：
- 计算hit level
- 查看hit level的分布，并确定阈值
- 计算归一化因子
- Reactivity预处理和归一化
- 合并结构变化区域
- 结果汇总

## 2.b) 计算 hit level
根据结果计算hit level 值，给出每条对应的hit level，按照深度进行区间划分的hit level平均值。
`hit_level.py` refer to `/data/TA_QUIZ_RNA_regulation/data/script/PartIII.SHAPE-seq_analysis/hit_level.py`

**reference:**

```linux
### hit level
### 1）calculate the average mutation events above background for each transcript

path0=/data/TA_QUIZ_RNA_regulation/data/script/PartIII.SHAPE-seq_analysis
path1=/data/TA_QUIZ_RNA_regulation/data/riboshape_liulab_batch4/final.modified_umodified/col_nouv

## 注意，这里的输入文件是final.modified_unmodified而非final.modified_unmodified.hit
python $path0/hit_level.py \
--data_path $path1/final.modified_unmodified \
--savepath_hit $path1/final.modified_unmodified.hit

echo -e "cutoff\ttranscript_id\tmodified_depth_median\tunmodified_depth_median\tmodified_depth_sum\tunmodified_depth_sum\thit"       >       $path1/cutoff.hit.group;
awk -F '\t' '$3>0 && $2<=25 && $2>0{print "0\t"$0}'     $path1/final.modified_unmodified.hit  >>        $path1/cutoff.hit.group;
awk -F '\t' '$3>0 && $2<=50 && $2 > 25{print "25\t"$0}' $path1/final.modified_unmodified.hit  >>        $path1/cutoff.hit.group;
awk -F '\t' '$3>0 && $2<=100 && $2 > 50{print "50\t"$0}' $path1/final.modified_unmodified.hit >>      $path1/cutoff.hit.group;
awk -F '\t' '$3>0 && $2<=200 && $2 > 100{print "100\t"$0}' $path1/final.modified_unmodified.hit >>      $path1/cutoff.hit.group;
awk -F '\t' '$3>0 && $2<=300 && $2 > 200{print "200\t"$0}' $path1/final.modified_unmodified.hit >>      $path1/cutoff.hit.group;
awk -F '\t' '$3>0 && $2<=500 && $2 >300{print "300\t"$0}' $path1/final.modified_unmodified.hit >>      $path1/cutoff.hit.group;
awk -F '\t' '$3>0 && $2<=750 && $2 >500{print "500\t"$0}' $path1/final.modified_unmodified.hit >>      $path1/cutoff.hit.group;
awk -F '\t' '$3>0 && $2<=1000 && $2 >750{print "750\t"$0}' $path1/final.modified_unmodified.hit >>      $path1/cutoff.hit.group;
awk -F '\t' '$3>0 && $2<=2000 && $2 >1000{print "1000\t"$0}' $path1/final.modified_unmodified.hit >>      $path1/cutoff.hit.group;
awk -F '\t' '$3>0 && $2<=5000 && $2 >2000{print "2000\t"$0}' $path1/final.modified_unmodified.hit >>      $path1/cutoff.hit.group;
awk -F '\t' '$3>0 && $2 > 5000{print "5000\t"$0}' $path1/final.modified_unmodified.hit >>      $path1/cutoff.hit.group;


```

**result:**
`final.modified_unmodified.hit`：其中第6列即为转录本对应的`hit level`的值。
`cutoff.hit.group`：是对`final.modified_unmodified.hit`文件的质控过滤后的结果。筛选，其中`unmodified.median>0,modified.median处于各种范围`的行。

## 2.c) 查看hit level 分布并确定阈值
脚本参考：`/data/TA_QUIZ_RNA_regulation/script/PartIII.SHAPE-seq_analysis/structure_changes_analysis/2.hit_level.py`

```python
#####1.统计不同hit level的转录本个数#######
import numpy as np
import pandas as pd
import re
import seaborn as sns
import matplotlib.pyplot as plt
from numpy import median
from numba import jit

#########转录本注释和数据读入#################
#### path #####
exon_gtf_path='/data/TA_QUIZ_RNA_regulation/data/ATH/GTF/shape_map/ath_exons.gtf'
col_uv_f_path='/data/TA_QUIZ_RNA_regulation/data/riboshape_liulab_batch4/final.modified_umodified/col_nouv/'
col_uv_z_path='/data/TA_QUIZ_RNA_regulation/data/riboshape_liulab_batch4/final.modified_umodified/col_uv/'
uvr_uv_f_path='/data/TA_QUIZ_RNA_regulation/data/riboshape_liulab_batch4/final.modified_umodified/uvr8_nouv/'
uvr_uv_z_path='/data/TA_QUIZ_RNA_regulation/data/riboshape_liulab_batch4/final.modified_umodified/uvr8_uv/'

### 读取注释文件，在exon_gtf文件的第8列中提取'transcript_id','gene_type','gene_name','chr'等信息
gtf_data=pd.read_csv(exon_gtf_path,sep='\t',header=None)
gtf_data_new=pd.DataFrame(columns={'transcript_id','gene_type','gene_name','chr'})
gtf_data_new['transcript_id']=gtf_data.iloc[:,8].apply(lambda x:x.split(';')[1].split('"')[1])
gtf_data_new['gene_type'] = gtf_data.iloc[:,8].apply(lambda x:re.findall('gene_biotype ".*?"',x)[0].split('"')[1])
gtf_data_new['gene_name'] = gtf_data.iloc[:,8].apply(lambda x:re.findall('gene_name ".*?"',x)[0].split('"')[1] if 'gene_name' in x else np.nan)
gtf_data_new['chr'] = gtf_data.iloc[:,0].apply(lambda x: 6 if x=='Mt' else 7 if x=='Pt' else x ).astype('int')
gtf_data_new = gtf_data_new.drop_duplicates()
gtf_data_new.index = range(len(gtf_data_new))

#  提取col_nouv、col_uv、uvr8_nouv、uvr8_uv中的'transcript_id','modified.median','unmodified.median','modified.sum','unmodified.sum','hit'列，并将ath_exons参考基因组中的'gene_type'，'gene_name','chr'的信息合并到原数据中

hit_level_col_uv_f = pd.read_csv(col_uv_f_path+'/final.modified_unmodified.hit',sep='\t',header=None)
hit_level_col_uv_f.columns =['transcript_id','modified.median','unmodified.median','modified.sum','unmodified.sum','hit']
hit_level_col_uv_f = pd.merge(hit_level_col_uv_f,gtf_data_new,on='transcript_id',how='left')
hit_level_col_uv_f['Type'] = 'WT_UV-'

hit_level_col_uv_z = pd.read_csv(col_uv_z_path+'/final.modified_unmodified.hit',sep='\t',header=None)
hit_level_col_uv_z.columns =['transcript_id','modified.median','unmodified.median','modified.sum','unmodified.sum','hit']
hit_level_col_uv_z = pd.merge(hit_level_col_uv_z,gtf_data_new,on='transcript_id',how='left')
hit_level_col_uv_z['Type'] = 'WT_UV+'

hit_level_uvr_uv_f = pd.read_csv(uvr_uv_f_path+'/final.modified_unmodified.hit',sep='\t',header=None)
hit_level_uvr_uv_f.columns =['transcript_id','modified.median','unmodified.median','modified.sum','unmodified.sum','hit']
hit_level_uvr_uv_f = pd.merge(hit_level_uvr_uv_f,gtf_data_new,on='transcript_id',how='left')
hit_level_uvr_uv_f['Type'] = 'UVR8_UV-'

hit_level_uvr_uv_z = pd.read_csv(uvr_uv_z_path+'/final.modified_unmodified.hit',sep='\t',header=None)
hit_level_uvr_uv_z.columns =['transcript_id','modified.median','unmodified.median','modified.sum','unmodified.sum','hit']
hit_level_uvr_uv_z = pd.merge(hit_level_uvr_uv_z,gtf_data_new,on='transcript_id',how='left')
hit_level_uvr_uv_z['Type'] = 'UVR8_UV+'

##计算col_nouv、col_uv、uvr8_nouv、uvr8_uv中hit_level>-1000,hit_level>0.hit_level>1,hit_level>2,hit_level>5,hit_level>10,hit_level>15的数目并整理成数据框格式。
data1=pd.DataFrame(columns={'Type','Number of transcripts','hit'})
for num in [-1000,0,1,2,5,10,15]:
    WT_UV_f = len(hit_level_col_uv_f.loc[(hit_level_col_uv_f['hit'] > num) & (hit_level_col_uv_f['modified.median'] > 100), :])
    WT_UV_z = len(hit_level_col_uv_z.loc[(hit_level_col_uv_z['hit'] > num) & (hit_level_col_uv_z['modified.median'] > 100)])
    Uvr8_UV_f = len(hit_level_uvr_uv_f.loc[(hit_level_uvr_uv_f['hit'] > num) & (hit_level_uvr_uv_f['modified.median'] > 100)])
    Uvr8_UV_z = len(hit_level_uvr_uv_z.loc[(hit_level_uvr_uv_z['hit'] > num) & (hit_level_uvr_uv_z['modified.median'] > 100)])

    data2 = pd.DataFrame(columns={'Type', 'Number of transcripts', 'hit'})
    data2['Type'] = ['WT_UV-', 'WT_UV+', 'UVR8_UV-', 'UVR8_UV+']
    data2['Number of transcripts'] = [WT_UV_f, WT_UV_z, Uvr8_UV_f, Uvr8_UV_z]

    if num==-1000:
        data2['hit'] ='ALL transcripts'
    else:
        data2['hit']='Hit level>'+str(num)
    data1 = pd.concat([data1,data2])

# 将hit level分布情况画图
#注意matplotlib.pyplot的默认画板在linux下无法被调用，需要用plt.switch_backend('agg')命令来设定画板。
plt.switch_backend('agg')
plt.figure(figsize=(8, 6))
sns.set(style="ticks", context="talk")
sns.axes_style({'font.family': ['sans-serif'],'font.sans-serif': ['Arial']})
g = sns.barplot(y='Number of transcripts',x='hit',hue='Type',data=data1,palette=["#909497", "#212F3C", "#3498DB", "#1F618D"])

sns.despine()
font1 = {'family' : 'Arial','weight' : 'roman','size': 22}
plt.xticks(rotation=60)

plt.legend(fontsize='small')
plt.tight_layout()
plt.savefig('/data/TA_QUIZ_RNA_regulation/result/PartIII.SHAPE-seq_analysis/plot/2a.transcript_num_hit.png')
plt.close()
```

注意matplotlib.pyplot的默认画板在linux下无法被调用，需要用plt.switch_backend('agg')命令来设定画板。

**result:**
`2a.transcript_num_hit.png`：转录本个数随hit level 变化趋势图。

![10.4.5.transcript_num_hit.png](../../.gitbook/assets/10.4.5.transcript_num_hit.png)

根据转录本个数变化趋势，我们取hit level大于2的转录本为可信转录本，进行后续分析。

## 2.d) 归一化因子
rRNAs,lncRNA and tRNA, which were sequenced to high depths, and showed large hit level (hit level >10) across experimental probing conditions.
找到用于归一化的rRNAs/lncRNA/tRNA，脚本参考:`/data/TA_QUIZ_RNA_regulation/script/PartIII.SHAPE-seq_analysis/structure_changes_analysis/3.normalization_factor.py`

```python
import numpy as np
import pandas as pd
import re
import seaborn as sns
import matplotlib.pyplot as plt
from numpy import median
from numba import jit

#########转录本注释和数据读入#################
#### path #####
exon_gtf_path='/data/TA_QUIZ_RNA_regulation/data/ATH/GTF/shape_map/ath_exons.gtf'
col_uv_f_path='/data/TA_QUIZ_RNA_regulation/data/riboshape_liulab_batch4/final.modified_umodified/col_nouv/'
col_uv_z_path='/data/TA_QUIZ_RNA_regulation/data/riboshape_liulab_batch4/final.modified_umodified/col_uv/'
uvr_uv_f_path='/data/TA_QUIZ_RNA_regulation/data/riboshape_liulab_batch4/final.modified_umodified/uvr8_nouv/'
uvr_uv_z_path='/data/TA_QUIZ_RNA_regulation/data/riboshape_liulab_batch4/final.modified_umodified/uvr8_uv/'

### 读取注释文件，在exon_gtf文件的第8列中提取'transcript_id','gene_type','gene_name','chr'等信息
gtf_data=pd.read_csv(exon_gtf_path,sep='\t',header=None)
gtf_data_new=pd.DataFrame(columns={'transcript_id','gene_type','gene_name','chr'})
gtf_data_new['transcript_id']=gtf_data.iloc[:,8].apply(lambda x:x.split(';')[1].split('"')[1])
gtf_data_new['gene_type'] = gtf_data.iloc[:,8].apply(lambda x:re.findall('gene_biotype ".*?"',x)[0].split('"')[1])
gtf_data_new['gene_name'] = gtf_data.iloc[:,8].apply(lambda x:re.findall('gene_name ".*?"',x)[0].split('"')[1] if 'gene_name' in x else np.nan)
gtf_data_new['chr'] = gtf_data.iloc[:,0].apply(lambda x: 6 if x=='Mt' else 7 if x=='Pt' else x ).astype('int')
gtf_data_new = gtf_data_new.drop_duplicates()
gtf_data_new.index = range(len(gtf_data_new))

#  提取col_nouv、col_uv、uvr8_nouv、uvr8_uv中'cutoff.hit.group'的'group','transcript_id','modified.median','unmodified.median','modified.sum','unmodified.sum','hit'列，并将ath_exons参考基因组中的'gene_type'，'gene_name','chr'等信息合并到原数据中

hit_level_col_uv_f = pd.read_csv(col_uv_f_path+'/cutoff.hit.group',sep='\t')
hit_level_col_uv_f.columns =['group','transcript_id','modified.median','unmodified.median','modified.sum','unmodified.sum','hit']
hit_level_col_uv_f = pd.merge(hit_level_col_uv_f,gtf_data_new,on='transcript_id',how='left')
hit_level_col_uv_f['spe'] = 'WT_UV-'

hit_level_col_uv_z = pd.read_csv(col_uv_z_path+'/cutoff.hit.group',sep='\t')
hit_level_col_uv_z.columns =['group','transcript_id','modified.median','unmodified.median','modified.sum','unmodified.sum','hit']
hit_level_col_uv_z = pd.merge(hit_level_col_uv_z,gtf_data_new,on='transcript_id',how='left')
hit_level_col_uv_z['spe'] = 'WT_UV+'

hit_level_uvr_uv_f = pd.read_csv(uvr_uv_f_path+'/cutoff.hit.group',sep='\t')
hit_level_uvr_uv_f.columns =['group','transcript_id','modified.median','unmodified.median','modified.sum','unmodified.sum','hit']
hit_level_uvr_uv_f = pd.merge(hit_level_uvr_uv_f,gtf_data_new,on='transcript_id',how='left')
hit_level_uvr_uv_f['spe'] = 'uvr8_UV-'

hit_level_uvr_uv_z = pd.read_csv(uvr_uv_z_path+'/cutoff.hit.group',sep='\t')
hit_level_uvr_uv_z.columns =['group','transcript_id','modified.median','unmodified.median','modified.sum','unmodified.sum','hit']
hit_level_uvr_uv_z = pd.merge(hit_level_uvr_uv_z,gtf_data_new,on='transcript_id',how='left')
hit_level_uvr_uv_z['spe'] = 'uvr8_UV+'

###从col_uv_f、col_uv_z、uvr_uv_f、uvr_uv_z中筛选出gene_type属于['lncRNA','rRNA','tRNA']，hit_level的值>10 && modified.median>5000的转录本，并输出它的转录本ID###
col_uv_f_id = hit_level_col_uv_f.loc[(hit_level_col_uv_f.gene_type.isin(['lncRNA','rRNA','tRNA']))&(hit_level_col_uv_f['modified.median']>5000)&(hit_level_col_uv_f['hit']>10),'transcript_id']
col_uv_z_id = hit_level_col_uv_z.loc[(hit_level_col_uv_z.gene_type.isin(['lncRNA','rRNA','tRNA']))&(hit_level_col_uv_z['modified.median']>5000)&(hit_level_col_uv_z['hit']>10),'transcript_id']
uvr_uv_f_id = hit_level_uvr_uv_f.loc[(hit_level_uvr_uv_f.gene_type.isin(['lncRNA','rRNA','tRNA']))&(hit_level_uvr_uv_f['modified.median']>5000)&(hit_level_uvr_uv_f['hit']>10),'transcript_id']
uvr_uv_z_id = hit_level_uvr_uv_z.loc[(hit_level_uvr_uv_z.gene_type.isin(['lncRNA','rRNA','tRNA']))&(hit_level_uvr_uv_z['modified.median']>5000)&(hit_level_uvr_uv_z['hit']>10),'transcript_id']

print(col_uv_f_id)
print(col_uv_z_id)
print(uvr_uv_f_id)
print(uvr_uv_z_id)
print(set(col_uv_f_id)&set(col_uv_z_id)&set(uvr_uv_f_id)&set(uvr_uv_z_id))
```
输出用于归一化的因子为{'AT3G41768.1', 'AT3G06355.1'}

## 2.e) Reactivity预处理和归一化
首先我们定义深度不够、不够可信位点的Reactivity为空值。然后我们利用归一化的因子(这里为AT3G41768.1，AT3G06355.1)进行Reactivity归一化。

- **1） Reactivity is NAN**
Nucleotides with apparent mutation rates above 0.02 in any untreated sample were excluded from analysis. 
These artifacts were identified as regions of at least 10 nucleotides in which three or more of the 10 nucleotides showed mutation rates above 0.03 in the absence of NAI treatment, or modified mutation rates above 0.1 in any condition, and were also excluded from analysis.
Read depth of nucleotides below 100 in any condition were also excluded from analysis.

- **2）Reactivity normalization**
Reactivities were normalized with in each probing condition to the mean of the 92-98th percentile reactivities of nucleotides from the rRNAs, which were sequenced to high depths, and showed large hit level (hit level >0) across experimental probing conditions.

脚本参考:`/data/TA_QUIZ_RNA_regulation/script/PartIII.SHAPE-seq_analysis/structure_changes_analysis/4.data_normalization.py`

```python
import numpy as np
import pandas as pd
import re
from numba import jit
import math

@jit(nopython=False)
def calc_quartile(x, q, qtype=7):
    # 计算给定的四分位数，x为所需计算的四分位数的列表，q代表是第25%百分位数还是75%百分位数。qtype代表所选用的方法。
    # source: http://adorio-research.org/wordpress/?p=125
    # x = array, q = quartile (in % as a decimal)
    y = np.copy(x)
    n = len(y)
    abcd = [(0, 0, 1, 0),  # inverse empirical distrib.function., R type 1
            (0.5, 0, 1, 0),  # similar to type 1, averaged, R type 2
            (0.5, 0, 0, 0),  # nearest order statistic,(SAS) R type 3
            (0, 0, 0, 1),  # California linear interpolation, R type 4
            (0.5, 0, 0, 1),  # hydrologists method, R type 5
            (0, 1, 0, 1),  # mean-based estimate(Weibull method), (SPSS,Minitab), type 6
            (1, -1, 0, 1),  # mode-based method,(S, S-Plus), R type 7
            (1.0 / 3, 1.0 / 3, 0, 1),  # median-unbiased ,  R type 8
            (3 / 8.0, 0.25, 0, 1)  # normal-unbiased, R type 9.
            ]
    a, b, c, d = abcd[qtype - 1]
    g, j = math.modf(a + (n + b) * q - 1)
    #第1四分位数的计算公式为 （n+3)/4,这里为(n-1)/4;第3四分位数的计算公式为(3n+1)/4,这里为(3n-3)/4。是由于数组中的元素是从x[0]开始，因此需要-1.
    if j < 0:
        return x[0]
    elif j >= n:
        return x[n - 1]
    j = int(math.floor(j))
    if g == 0:
        return x[j]
    else:
        return y[j] + (y[j + 1] - y[j]) * (c + d * g)

# @jit(nopython=False)
def find_boxplot_factor(array):
    x, o, a = [], [], 0
    # Following deprecated line is behavior that normalization and
    # structure modeling were optimized with, but this behavior
    # is probably not ideal. For RNAs with regions of poor sequencing
    # depth, treating those regions as unreactive could bias the
    # normalized reactivities. This is especially important for
    # larger RNAs.
    # x = np.fromiter((n if not isnan(n) else 0 for n in array))
    x = array[np.where(np.isfinite(array))]
    # 除去数组中的非有限数的值。
    x = x[np.nonzero(x)]
    # 输出不为0的元素的下标
    if x.shape[0] < 10:
        norm_factor = np.nan
    else:
        x.sort()
        ten_pct = len(x) // 10
        five_pct = len(x) // 20
        # calculate the interquartile range *1.5
        q_limit = 1.5 * abs(calc_quartile(x, 0.25) - calc_quartile(x, 0.75))
        ten_limit = x[x.shape[0] - 1 - ten_pct]
        five_limit = x[x.shape[0] - 1 - five_pct]
        # choose the cutoff that eliminates the fewest points
        limit = max(q_limit, ten_limit)
        if len(x) < 100:
            limit = max(q_limit, five_limit)
        # make new list without the outliers
        for i in range(len(x)):
            if x[i] < limit:
                o.append(x[i])
        # avg next ten percent
        try:
            for i in range(-ten_pct, 0):
                a = o[i] + a
            norm_factor = a / ten_pct
        except IndexError:
            norm_factor = np.nan
    return norm_factor

@jit(nopython=False)
def find_boxplot_factor_new(array):
    ####the mean of the 92-98th percentile reactivities#####
    # x, o, a = [], [], []
    x = array[np.where(np.isfinite(array))]
    x = x[np.nonzero(x)]
    if x.shape[0] < 10:
        norm_factor = np.nan
    else:

        # calculate the interquartile range *1.5
        o = x[(x>=np.quantile(x,0.92))&(x<=np.quantile(x,0.98))]
        norm_factor = np.mean(o)
    return norm_factor


@jit(nopython=False)
def filter(X,Mutru_cut_off=0.02,read_depth_cutoff=100,R_window=10,window_mutru_cutoff=0.03,window_mutru_num=3,window_mutrs_cutoff=0.1,window_mutrs_num=3):

    name = X[0]
    nucleotide = X[1]
    modified = np.array(X[2].split(',')).astype('float').astype('int')
    modified_depth = np.array(X[3].split(',')).astype('float').astype('int')
    unmodified = np.array(X[4].split(',')).astype('float').astype('int')
    unmodified_depth = np.array(X[5].split(',')).astype('float').astype('int')

    # print(nucleotide)
    Mutrs = modified/modified_depth
    Mutru = unmodified/unmodified_depth
    R = Mutrs - Mutru
    n = len(modified)
    ##### 校正reactivity###########
    R[Mutru>Mutru_cut_off] = np.nan
    R[(modified_depth <= read_depth_cutoff)|(unmodified_depth)<=read_depth_cutoff] = np.nan
    R[R<0]= np.nan
    for i in range(len(R)-R_window+1):
        data_Mutru = Mutru[i:i+R_window-1]
        data_Mutrs = Mutrs[i:i+R_window-1]
        if (len(data_Mutru[data_Mutru>window_mutru_cutoff])>=window_mutru_num)|(len(data_Mutrs[data_Mutrs>window_mutrs_cutoff])>=window_mutrs_num):
            R[i:i + R_window-1] = np.nan
    R_new = ",".join(list(R.astype('str'))).replace('nan','-999')
    X = X.append(pd.Series(R_new))
    return X
    #
@jit(nopython=False)
def get_fatcor(data):
    R_all = np.array(data.iloc[0,6].split(',')).astype('float')
    for i in range(1,len(data)):
        R_ = np.array(data.iloc[i,6].split(',')).astype('float')
        R_all = np.hstack((R_all,R_))
    R_all[R_all==-999]=np.nan
    factor = find_boxplot_factor_new(R_all)
    return factor

@jit(nopython=False)
def normalization_all(X,factor):
    #########读取数据#############
    R = np.array(X[6].split(',')).astype('float')
    R_nan = len(R[R == -999])/len(R)
    R_0 = len(R[R == 0])/len(R)
    R[R == -999] = np.nan
    R_new = R/factor
    R_new = ",".join(list(R_new.astype('str'))).replace('nan','-999')
    X[6] =R_new
    return X,R_nan,R_0
def main(data,transcript_list=[]):

    data_filter = pd.DataFrame(np.zeros([len(data),7]))
    data_filter.columns = ['transcript_id','Nucleotide' ,'Modified_mutations', 'Modified_effective_depth',
                             'Untreated_mutations', 'Untreated_effective_depth', 'reactivity']
    data_new = pd.DataFrame(np.zeros([len(data), 7]))
    data_new.columns = ['transcript_id','Nucleotide' ,'Modified_mutations', 'Modified_effective_depth',
                           'Untreated_mutations', 'Untreated_effective_depth', 'reactivity']
    print(len(data))
    for i in range(len(data)):
        if i %1000 == 0:
            print(i)
        X=filter(data.iloc[i,:])
        data_filter.iloc[i,:] = list(X)
    if len(transcript_list)> 0:
        factor = get_fatcor(data_filter.loc[data_filter.transcript_id.isin(transcript_list),:])
    else:
        factor = get_fatcor(data_filter)
    for i in range(len(data_filter)):
        if i %1000 == 0:
            print(i)
            # print data_
        data_new.iloc[i,:],R_nan,R_0 = normalization_all(data_filter.iloc[i,:],factor)
        # data_new.loc[i, 'ratio_nan']=R_nan
        # data_new.loc[i, 'ratio_0'] = R_0
    data_new = data_new[['transcript_id','Nucleotide', 'Modified_mutations', 'Modified_effective_depth',
                           'Untreated_mutations', 'Untreated_effective_depth', 'reactivity']]
    return data_new

def loctaion(X):
    X_1 = X.split('(')[1].split(')')[0].split(',')
    loctaion_list = []
    for i in range(len(X_1)):
        loctaion_list_ = [i for i in range(int(X_1[i].split(':')[0]),int(X_1[i].split(':')[1])+1)]
        loctaion_list.extend(loctaion_list_)
    loctaion_list =np.array(loctaion_list)
    loctaion_list = ",".join(list(loctaion_list.astype('str')))

    return loctaion_list



def mapping(exon_data):
    data_location = pd.DataFrame(np.zeros([len(exon_data),4]))
    data_location.columns = ['transcript_id','chr','strand','location']

    for i in range(len(exon_data)):
        if i%1000 ==0:
            print(i)
        data_location.loc[i, 'transcript_id'] = exon_data.loc[i, 'transcript_id']
        data_location.loc[i,'chr'] = exon_data.loc[i,'chr'].split('.')[0]
        data_location.loc[i, 'strand'] = exon_data.loc[i, 'chr'].split('.')[1]
        data_location.loc[i, 'location'] = loctaion(exon_data.loc[i, 'start_end'])
    return data_location

if __name__ == '__main__':
    col_uv_f = '/data/TA_QUIZ_RNA_regulation/data/riboshape_liulab_batch4/final.modified_umodified/col_nouv/'
    col_uv_z = '/data/TA_QUIZ_RNA_regulation/data/riboshape_liulab_batch4/final.modified_umodified/col_uv/'
    uvr_uv_f = '/data/TA_QUIZ_RNA_regulation/data/riboshape_liulab_batch4/final.modified_umodified/uvr8_nouv/'
    uvr_uv_z = '/data/TA_QUIZ_RNA_regulation/data/riboshape_liulab_batch4/final.modified_umodified/uvr8_uv/'
    path=[col_uv_f,col_uv_z,uvr_uv_f,uvr_uv_z]
    for i in range(len(path)):
        transcript_list = [ 'AT3G41768.1', 'AT3G06355.1']
        data_uv_z = pd.read_csv(path[i] + '/final.modified_unmodified', sep='\t')
        data_uv_z.columns = ['transcript_id', 'Nucleotide', 'Modified_mutations', 'Modified_effective_depth','Untreated_mutations', 'Untreated_effective_depth']
        data_uv_z_new=main(data_uv_z,transcript_list)
        print(data_uv_z_new)
        data_uv_z_new.to_csv(path[i]+'/final.modified_unmodified_new',sep='\t',header=True,index=False)
```

**output：**
`final.modified_unmodified_new`：相较于`final.modified_unmodified`多了一列`reactivity`信息。

## 2.f) 滑动窗口Gini index 计算
我们利用滑动窗口计算一个转录本的结构变化区域。需要注意以下几点：
- 每个滑动窗口非空值大于90%（ratio_nan=0.9）
- 滑动窗口大小为50nt，步长为1nt
- UV-和UV+数据相比，存在abs(delta_gini)>0.1的滑动窗口，该转录本为结构改变的转录本，该区域为结构改变区域，一个转录本可能有多个结构改变区域。
- 计算了整个转录本、转录本3'UTR、转录本5'UTR、转录本CDS区域的Gini index

脚本参考：`/data/TA_QUIZ_RNA_regulation/script/PartIII.SHAPE-seq_analysis/structure_changes_analysis/5.calculate_gini_index.py`
```python
########################计算gini_index#############################
import numpy as np
import pandas as pd
import re
from numba import jit
from scipy.stats import ks_2samp
import statsmodels.stats.multitest as multi

jit(nopython=False)
def get_gini(R_, nucleotide_, cut_AC=False, ratio_nan=0.9, ratio_nan_AC=0.8):
    ########
    if len(R_)==0:
        gini=np.nan
    else:
        R_nonull = R_[np.isnan(R_) == False]
        # R_nonull = R_[R_>0]
        ratio = len(R_nonull) / len(R_)
        if ratio <= ratio_nan:
            gini = np.nan
        else:
            if cut_AC:
                R_AC = R_[(nucleotide_ == b'A') | (nucleotide_ == b'C')]
                R_AC_nonull = R_AC[np.isnan(R_AC) == False]
                ratio_AC = len(R_AC_nonull) / len(R_AC)
                if (ratio_AC <= ratio_nan_AC) | len(R_AC) <= 1 | len(R_AC_nonull) <= 1:
                    gini = np.nan
                else:
                    sorted = np.sort(R_AC_nonull)
                    height, area = 0, 0
                    for i in range(0, len(sorted)):
                        height += sorted[i]
                        area += height - sorted[i] / 2.
                    fair_area = height * len(sorted) / 2.
                    if fair_area == 0:
                        gini = np.nan
                    else:
                        gini = (fair_area - area) / fair_area
            else:
                sorted = np.sort(R_nonull)
                height, area = 0, 0
                for i in range(0, len(sorted)):
                    height += sorted[i]
                    area += height - sorted[i] / 2.
                fair_area = height * len(sorted) / 2.
                if fair_area == 0:
                    gini = np.nan
                else:
                    gini = (fair_area - area) / fair_area
    return gini

@jit(nopython=False)
def get_p(R_1,R_2,ratio_nan=0.9):
    R_1_nonull = R_1[np.isnan(R_1) == False]
    ratio_1 = len(R_1_nonull) / len(R_1)
    R_2_nonull = R_2[np.isnan(R_2) == False]
    ratio_2 = len(R_2_nonull) / len(R_2)
    if (ratio_1 <= ratio_nan)|(ratio_2 <= ratio_nan):
        p = np.nan
    else:
        s,p=ks_2samp(R_1, R_2)
    return p

@jit(nopython=False)
def get_window(X, window=50, step=1, cut_AC=False, ratio_nan=0.9, ratio_nan_AC=0.8):
    #########读取数据#############
    name = X[0]

    nucleotide = np.array(list(X[1]))

    modified = np.array(X[2].split(',')).astype('float').astype('int')


    modified_depth = np.array(X[3].split(',')).astype('float').astype('int')

    #
    unmodified = np.array(X[4].split(',')).astype('float').astype('int')


    unmodified_depth = np.array(X[5].split(',')).astype('float').astype('int')
    Mutrs = modified / modified_depth
    Mutru = unmodified / unmodified_depth
    # R = Mutrs - Mutru
    R = np.array(X[6].split(',')).astype('float')
    R[R == -999] = np.nan
    n = len(nucleotide)

    #####滑动窗口计算change##################
    nucleotide_ = np.zeros([len(range(0, n - (window), step)), window], dtype=np.string_)
    R_ = np.zeros([len(range(0, n - (window), step)), window])
    gini = np.zeros([len(range(0, n - (window), step)), ])
    j = 0
    for i in range(0, n - (window), step):
        nucleotide_[j, :] = nucleotide[i:i + window]
        R_[j, :] = R[i:i + window]
        gini[j] = get_gini(R[i:i + window], nucleotide[i:i + window], cut_AC=cut_AC, ratio_nan=ratio_nan,
                           ratio_nan_AC=ratio_nan_AC)
        j = j + 1

    return nucleotide_, R_, gini

@jit(nopython=False)
def get_window_p(X1,X2, window=50, step=1,ratio_nan=0.9):
    #########读取数据#############
    name = X1[0]
    nucleotide = np.array(list(X1[1]))
    R_z = np.array(X1[6].split(',')).astype('float')
    R_z[R_z == -999] = np.nan
    R_f = np.array(X2[6].split(',')).astype('float')
    R_f[R_f == -999] = np.nan
    n = len(nucleotide)

    #####滑动窗口计算change##################
    nucleotide_ = np.zeros([len(range(0, n - (window), step)), window], dtype=np.string_)
    R_z_ = np.zeros([len(range(0, n - (window), step)), window])
    R_f_ = np.zeros([len(range(0, n - (window), step)), window])
    p = np.zeros([len(range(0, n - (window), step)), ])
    j = 0
    for i in range(0, n - (window), step):
        nucleotide_[j, :] = nucleotide[i:i + window]
        R_z_[j, :] = R_z[i:i + window]
        R_f_[j, :] = R_f[i:i + window]
        p[j] = get_p(R_z[i:i + window], R_f[i:i + window],ratio_nan=ratio_nan,)
        j = j + 1
    p_=p.copy()
    if len(p[~np.isnan(p)])>0:
        a,p_bh,b,c =multi.multipletests(p[~np.isnan(p)],method='fdr_bh')
        p_[~np.isnan(p_)]=p_bh
    return p_,p


@jit(nopython=False)
def calculate_delta_gini(R_1, R_2, gini_1, gini_2, ratio_nan=0.9):
    '''Calculates Standard RMSD on two vectors of numbers of the same length'''
    # Check to see the vectors are of equal length.
    #计算delta_gini,(1)当R_1,R_2长度不同时，输出nan值;(2)当滑动窗口非空值比率<90%，gini_index输出空值。(3)除(1)(2)情况外，输出delta_gini。
    if len(R_1) != len(R_2):
        return np.nan
    else:
        R = R_1 - R_2
        if len(R[np.isnan(R) == False]) / len(R) <= ratio_nan:
            return np.nan
        else:
            delta_gini = gini_1 - gini_2
    return delta_gini

@jit(nopython=False)
def get_all(X,cut_AC=False, ratio_nan=0, ratio_nan_AC=0.8):
    #########读取数据#############
    name = X[0]
    nucleotide = np.array(list(X[1]))
    R = np.array(X[6].split(',')).astype('float')
    R[R == -999] = np.nan
    n = len(nucleotide)
    gini= get_gini(R, nucleotide, cut_AC=cut_AC, ratio_nan=ratio_nan,ratio_nan_AC=ratio_nan_AC)
    return gini

@jit(nopython=False)
def get_se(X,location ,start,end,strand,cut_AC=False, ratio_nan=0, ratio_nan_AC=0.8):
    #########读取数据#############
    name = X[0]
    nucleotide = np.array(list(X[1]))
    R = np.array(X[6].split(',')).astype('float')
    R[R == -999] = np.nan
    n = len(nucleotide)
    if strand=='+':
        R_=R[location<start]
        nucleotide_=nucleotide[location<start]
        gini_5= get_gini(R_, nucleotide_, cut_AC=cut_AC, ratio_nan=ratio_nan,ratio_nan_AC=ratio_nan_AC)
        R_ = R[location > end]
        nucleotide_ = nucleotide[location > end]
        gini_3 = get_gini(R_, nucleotide_, cut_AC=cut_AC, ratio_nan=ratio_nan, ratio_nan_AC=ratio_nan_AC)
    else:
        R_ = R[location < start]
        nucleotide_ = nucleotide[location < start]
        gini_3 = get_gini(R_, nucleotide_, cut_AC=cut_AC, ratio_nan=ratio_nan, ratio_nan_AC=ratio_nan_AC)
        R_ = R[location > end]
        nucleotide_ = nucleotide[location > end]
        gini_5 = get_gini(R_, nucleotide_, cut_AC=cut_AC, ratio_nan=ratio_nan, ratio_nan_AC=ratio_nan_AC)
    R_ = R[(location<= end)&(location>=start)]
    nucleotide_ = nucleotide[(location<= end)&(location>=start)]
    gini_cds = get_gini(R_, nucleotide_, cut_AC=cut_AC, ratio_nan=ratio_nan, ratio_nan_AC=ratio_nan_AC)
    return gini_3,gini_cds,gini_5


@jit(nopython=False)
def get_statistics(data_z, data_f,data_gff_,location_list,strand):
    data_z = np.array(data_z).reshape([7, ])
    data_f = np.array(data_f).reshape([7, ])
    nucleotide_z, R_z, gini_z = get_window(data_z)
    nucleotide_f, R_f, gini_f = get_window(data_f)
    ###uv+###
    gini_z_=gini_z[pd.notnull(gini_z)]
    if len(gini_z_)==0:
        gini_max_z=np.nan
    else:
        gini_max_z=gini_z[pd.notnull(gini_z)].max()
    gini_all_z=get_all(data_z)
    if len(data_gff_.loc[data_gff_['location']=='CDS','strat'])==0:
        gini_3_z=np.nan
        gini_5_z=np.nan
        gini_cds_z=np.nan
    else:
        CDS_strat=list(data_gff_.loc[data_gff_['location']=='CDS','strat'])[0]
        CDS_end = list(data_gff_.loc[data_gff_['location']=='CDS','end'])[0]
        gini_3_z, gini_cds_z, gini_5_z = get_se(data_z,location_list,CDS_strat,CDS_end,strand)
    ###uv-###
    gini_f_ = gini_f[pd.notnull(gini_f)]
    if len(gini_f_) == 0:
        gini_max_f = np.nan
    else:
        gini_max_f = gini_f[pd.notnull(gini_f)].max()
    gini_all_f = get_all(data_f)
    if len(data_gff_.loc[data_gff_['location'] == 'CDS', 'strat']) == 0:
        gini_3_f = np.nan
        gini_5_f = np.nan
        gini_cds_f = np.nan
    else:
        CDS_strat = list(data_gff_.loc[data_gff_['location'] == 'CDS', 'strat'])[0]
        CDS_end = list(data_gff_.loc[data_gff_['location'] == 'CDS', 'end'])[0]
        gini_3_f, gini_cds_f, gini_5_f = get_se(data_f, location_list, CDS_strat, CDS_end, strand)

    if len(R_z) != len(R_f):
        print("error")
    else:
        delta_gini = np.zeros([len(R_z), ])
        for i in range(len(R_z)):
            delta_gini[i] = calculate_delta_gini(R_z[i, :], R_f[i, :], gini_z[i], gini_f[i])
    return delta_gini,gini_max_z,gini_all_z,gini_3_z, gini_cds_z, gini_5_z,gini_max_f,gini_all_f,gini_3_f, gini_cds_f, gini_5_f


def main_sum_gini(data_uv_z, data_uv_f, data_location, data_gff, transcript_id_list):
    statistics_sum = pd.DataFrame(columns={'transcript_id', 'num', 'num_0.1','delta_max', 'delta_min'})
    statistics_z = pd.DataFrame(columns={'transcript_id', 'gini_all', 'gini_max', 'gini_3UTR', 'gini_5UTR', 'gini_CDS'})
    statistics_f = pd.DataFrame(columns={'transcript_id', 'gini_all', 'gini_max', 'gini_3UTR', 'gini_5UTR', 'gini_CDS'})
    i = 0
    for transcript in transcript_id_list:
        data_z = data_uv_z.loc[data_uv_z['transcript_id'] == transcript, :]
        data_f = data_uv_f.loc[data_uv_f['transcript_id'] == transcript, :]
        data_location_ = data_location.loc[data_location['transcript_id'] == transcript,:]
        data_gff_ = data_gff.loc[data_gff['transcript_id'] == transcript, :]
        location_list=np.array(list(data_location_['location'])[0].split(',')).astype('int')
        chr = list(data_location_['chr'])[0]
        strand=list(data_location_['strand'])[0]
        delta_gini, gini_max_z, gini_all_z, gini_3_z, gini_cds_z, gini_5_z, gini_max_f, gini_all_f, gini_3_f, gini_cds_f, gini_5_f= get_statistics(data_z,data_f, data_gff_, location_list, strand)

        statistics_sum.loc[i, 'transcript_id'] = transcript
        num = sum(np.abs(delta_gini) >=0.2)
        statistics_sum.loc[i, 'num'] = num
        num_1 = sum(np.abs(delta_gini) >=0.1)
        statistics_sum.loc[i, 'num_0.1'] = num_1
        statistics_sum.loc[i, 'delta_gini_list'] = ','.join(list(delta_gini.astype('str')))

        if len(delta_gini[np.isnan(delta_gini) == False]) > 0:
            statistics_sum.loc[i, 'delta_max'] = np.max(delta_gini[np.isnan(delta_gini) == False])
            statistics_sum.loc[i, 'delta_min'] = np.min(delta_gini[np.isnan(delta_gini) == False])
        else:
            statistics_sum.loc[i, 'delta_max'] = np.nan
            statistics_sum.loc[i, 'delta_min'] = np.nan

        statistics_z.loc[i,'transcript_id']=transcript
        statistics_z.loc[i,'gini_all']=gini_all_z
        statistics_z.loc[i, 'gini_max'] = gini_max_z
        statistics_z.loc[i, 'gini_3UTR'] = gini_3_z
        statistics_z.loc[i, 'gini_5UTR'] = gini_5_z
        statistics_z.loc[i, 'gini_CDS'] = gini_cds_z

        statistics_f.loc[i,'transcript_id']=transcript
        statistics_f.loc[i,'gini_all']=gini_all_f
        statistics_f.loc[i, 'gini_max'] = gini_max_f
        statistics_f.loc[i, 'gini_3UTR'] = gini_3_f
        statistics_f.loc[i, 'gini_5UTR'] = gini_5_f
        statistics_f.loc[i, 'gini_CDS'] = gini_cds_f
        i = i + 1
        if i % 100 == 0:
            print(i)
    return statistics_sum,statistics_z,statistics_f

if __name__ == '__main__':
    col_uv_f = '/data/TA_QUIZ_RNA_regulation/data/riboshape_liulab_batch4/final.modified_umodified/col_nouv/'
    col_uv_z = '/data/TA_QUIZ_RNA_regulation/data/riboshape_liulab_batch4/final.modified_umodified/col_uv/'
    uvr_uv_f = '/data/TA_QUIZ_RNA_regulation/data/riboshape_liulab_batch4/final.modified_umodified/uvr8_nouv/'
    uvr_uv_z = '/data/TA_QUIZ_RNA_regulation/data/riboshape_liulab_batch4/final.modified_umodified/uvr8_uv/'

###col###
    result = '/data/TA_QUIZ_RNA_regulation/result/PartIII.SHAPE-seq_analysis/col'
    uv_z = col_uv_z
    uv_f = col_uv_f

    data_uv_z = pd.read_csv(uv_z + '/final.modified_unmodified_new', sep='\t')
    data_uv_z = data_uv_z[['transcript_id', 'Nucleotide', 'Modified_mutations', 'Modified_effective_depth',
                         'Untreated_mutations', 'Untreated_effective_depth', 'reactivity']]
    data_uv_z = data_uv_z.drop_duplicates()

    data_uv_f = pd.read_csv(uv_f + '/final.modified_unmodified_new', sep='\t')
    data_uv_f = data_uv_f[['transcript_id', 'Nucleotide', 'Modified_mutations', 'Modified_effective_depth',
                         'Untreated_mutations', 'Untreated_effective_depth', 'reactivity']]
    data_uv_f = data_uv_f.drop_duplicates()
    transcript_uv_z = pd.read_csv(uv_z + '/cutoff.hit.group', sep='\t')
    transcript_uv_z.columns = ['cut_off', 'transcript_id', 'modified_depth_median',
                               'unmodified_depth_median', 'modified_depth_sum', 'unmodified_depth_sum', 'hit_level']

    transcript_uv_f = pd.read_csv(uv_f + '/cutoff.hit.group', sep='\t')
    transcript_uv_f.columns = ['cut_off', 'transcript_id', 'modified_depth_median',
                               'unmodified_depth_median', 'modified_depth_sum', 'unmodified_depth_sum', 'hit_level']

    uv_z_transcript = list(transcript_uv_z.loc[(transcript_uv_z['modified_depth_median'] > 100) & (
            transcript_uv_z['hit_level'] > 0), 'transcript_id'])

    uv_f_transcript = list(transcript_uv_f.loc[(transcript_uv_f['modified_depth_median'] > 100) & (
            transcript_uv_f['hit_level'] > 0), 'transcript_id'])
    transcript_all = list(set(uv_z_transcript) & set(uv_f_transcript))
    print(len(transcript_all))
    gff_path = '/data/TA_QUIZ_RNA_regulation/data/ATH/GFF/Ath_genes.gff'
    data_gff = pd.read_csv(gff_path, sep='\t')
    data_location = pd.read_csv(
        '/data/TA_QUIZ_RNA_regulation/data/ATH/GTF/shape_map/result/transcript_exon_location.csv',
        sep='\t')
    statistics_sum, statistics_z, statistics_f = main_sum_gini(data_uv_z, data_uv_f, data_location, data_gff, transcript_all)
    pd.DataFrame(statistics_sum).to_csv(result + '/gini_summary_50_1.csv', sep='\t', index=False, header=True)
    pd.DataFrame(statistics_z).to_csv(result + '/gini_summary_UV+_50_1.csv', sep='\t', index=False, header=True)
    pd.DataFrame(statistics_f).to_csv(result + '/gini_summary_UV-_50_1.csv', sep='\t', index=False, header=True)

###uvr8###
    result = '/data/TA_QUIZ_RNA_regulation/result/PartIII.SHAPE-seq_analysis/uvr8'
    uv_z = uvr_uv_z
    uv_f = uvr_uv_f

    data_uv_z = pd.read_csv(uv_z + '/final.modified_unmodified_new', sep='\t')
    data_uv_z = data_uv_z[['transcript_id', 'Nucleotide', 'Modified_mutations', 'Modified_effective_depth',
                         'Untreated_mutations', 'Untreated_effective_depth', 'reactivity']]
    data_uv_z = data_uv_z.drop_duplicates()

    data_uv_f = pd.read_csv(uv_f + '/final.modified_unmodified_new', sep='\t')
    data_uv_f = data_uv_f[['transcript_id', 'Nucleotide', 'Modified_mutations', 'Modified_effective_depth',
                         'Untreated_mutations', 'Untreated_effective_depth', 'reactivity']]
    data_uv_f = data_uv_f.drop_duplicates()
    transcript_uv_z = pd.read_csv(uv_z + '/cutoff.hit.group', sep='\t')
    transcript_uv_z.columns = ['cut_off', 'transcript_id', 'modified_depth_median',
                               'unmodified_depth_median', 'modified_depth_sum', 'unmodified_depth_sum', 'hit_level']

    transcript_uv_f = pd.read_csv(uv_f + '/cutoff.hit.group', sep='\t')
    transcript_uv_f.columns = ['cut_off', 'transcript_id', 'modified_depth_median',
                               'unmodified_depth_median', 'modified_depth_sum', 'unmodified_depth_sum', 'hit_level']

    uv_z_transcript = list(transcript_uv_z.loc[(transcript_uv_z['modified_depth_median'] > 100) & (
            transcript_uv_z['hit_level'] > 0), 'transcript_id'])

    uv_f_transcript = list(transcript_uv_f.loc[(transcript_uv_f['modified_depth_median'] > 100) & (
            transcript_uv_f['hit_level'] > 0), 'transcript_id'])
    transcript_all = list(set(uv_z_transcript) & set(uv_f_transcript))
    print(len(transcript_all))
    gff_path = '/data/TA_QUIZ_RNA_regulation/data/ATH/GFF/Ath_genes.gff'
    data_gff = pd.read_csv(gff_path, sep='\t')
    data_location = pd.read_csv(
        '/data/TA_QUIZ_RNA_regulation/data/ATH/GTF/shape_map/result/transcript_exon_location.csv',
        sep='\t')
    statistics_sum, statistics_z, statistics_f = main_sum_gini(data_uv_z, data_uv_f, data_location, data_gff, transcript_all)
    pd.DataFrame(statistics_sum).to_csv(result + '/gini_summary_50_1.csv', sep='\t', index=False, header=True)
    pd.DataFrame(statistics_z).to_csv(result + '/gini_summary_UV+_50_1.csv', sep='\t', index=False, header=True)
    pd.DataFrame(statistics_f).to_csv(result + '/gini_summary_UV-_50_1.csv', sep='\t', index=False, header=True)
```

**output:**
`gini_summary_50_1.csv`为滑动窗口delta gini index结果
`gini_summary_UV+_50_1.csv`和`gini_summary_UV-_50_1.csv`为UV-和UV+下整个转录本、转录本3'UTR、转录本5'UTR、转录本CDS区域的Gini index。

## 2.g) 合并结构变化区域
我们利用滑动窗口判断每个转录本加光前后的结构变化区域，一个转录本可能有多个滑动窗口满足|delta gini index|>0.1。我们对连续的滑动窗口进行合并，我们采取逐步向前的策略，也就是如果新的滑动窗口加入后，重新计算整个区域|delta gini index|，如果整个区域|delta gini index|>0.1则合并，否则不合并。

连续的滑动窗口为滑动窗口起始点连续的，例如窗口[1,51],[2,52]。有重叠的滑动窗口为存在交叉重叠位点的滑动窗口，例如[1,51],[20,70]。我们这里合并连续的滑动窗口。

参考脚本:`/data/TA_QUIZ_RNA_regulation/script/PartIII.SHAPE-seq_analysis/structure_changes_analysis/6.gini_transcript_merge_windows_new.py`

```python
import numpy as np
import pandas as pd
import re
from numba import jit
from scipy.stats import ks_2samp
import statsmodels.stats.multitest as multi

####################函数定义部分####################
@jit(nopython=False)
def get_gini(R_, nucleotide_, cut_AC=False, ratio_nan=0.9, ratio_nan_AC=0.8):
    R_nonull = R_[np.isnan(R_)==False]
    ratio = len(R_nonull)/len(R_)
    if ratio <= ratio_nan:
        gini = np.nan
    else:
        if cut_AC:
            R_AC = R_[(nucleotide_==b'A')|(nucleotide_==b'C')]
            R_AC_nonull =  R_AC[np.isnan(R_AC)==False]
            ratio_AC = len(R_AC_nonull) / len(R_AC)
            if (ratio_AC <= ratio_nan_AC)|len(R_AC)<=1|len(R_AC_nonull)<=1:
                gini = np.nan
            else:
                sorted = np.sort(R_AC_nonull)
                height, area = 0, 0
                for i in range(0,len(sorted)):
                    height += sorted[i]
                    area += height - sorted[i] / 2.
                fair_area = height * len(sorted) / 2.
                if fair_area ==0:
                    gini=np.nan
                else:
                    gini =  (fair_area - area) / fair_area
        else:
            sorted = np.sort(R_nonull)
            height, area = 0, 0
            for i in range(0, len(sorted)):
                height += sorted[i]
                area += height - sorted[i] / 2.
            fair_area = height * len(sorted) / 2.
            if fair_area == 0:
                gini = np.nan
            else:
                gini = (fair_area - area) / fair_area
    return gini

@jit(nopython=False)
def get_p(R_1,R_2,ratio_nan=0.9):
    R_1_nonull = R_1[np.isnan(R_1) == False]
    ratio_1 = len(R_1_nonull) / len(R_1)
    R_2_nonull = R_2[np.isnan(R_2) == False]
    ratio_2 = len(R_2_nonull) / len(R_2)
    if (ratio_1 <= ratio_nan)|(ratio_2 <= ratio_nan):
        p = np.nan
    else:
        s,p=ks_2samp(R_1, R_2)
    return p

@jit(nopython=False)
def get_window(X,location,window=50,step=1,cut_AC=False,ratio_nan=0.9,ratio_nan_AC=0.8):
    #########读取数据#############
    name = X[0]
    nucleotide = np.array(list(X[1]))
    modified = np.array(X[2].split(',')).astype('int')
    modified_depth = np.array(X[3].split(',')).astype('int')
    unmodified = np.array(X[4].split(',')).astype('int')
    unmodified_depth = np.array(X[5].split(',')).astype('int')
    Mutrs = modified/modified_depth
    Mutru = unmodified/unmodified_depth
    # R = Mutrs - Mutru
    R = np.array(X[6].split(',')).astype('float')
    R[R == -999] = np.nan
    n = len(nucleotide)
    index = np.array(range(0, n))
    location_index = np.array(np.array(location['location'])[0].split(',')).astype('int')
    #####滑动窗口计算change##################
    nucleotide_ = np.zeros([len(range(0, n-(window), step)), window],dtype=np.string_)
    R_ = np.zeros([len(range(0, n-(window), step)), window])
    gini = np.zeros([len(range(0, n-(window), step)),])
    index_ = np.zeros([len(range(0, n-(window), step)), window])
    location_index_ = np.zeros([len(range(0, n - (window), step)), window])
    j =0
    for i in range(0, n-(window), step):
        index_[j,:]=index[i:i+window]
        location_index_[j,:] = location_index[i:i+window]
        nucleotide_[j,:] = nucleotide[i:i+window]
        R_[j,:] = R[i:i+window]
        gini[j] = get_gini(R[i:i+window],nucleotide[i:i+window],cut_AC=cut_AC,ratio_nan=ratio_nan,ratio_nan_AC=ratio_nan_AC)
        j=j+1

    return index_,location_index_,nucleotide_,R_,gini


@jit(nopython=False)
def get_window_p(X1,X2, window=50, step=1,ratio_nan=0.9):
    #########读取数据#############
    name = X1[0]
    nucleotide = np.array(list(X1[1]))
    R_z = np.array(X1[6].split(',')).astype('float')
    R_z[R_z == -999] = np.nan
    R_f = np.array(X2[6].split(',')).astype('float')
    R_f[R_f == -999] = np.nan
    n = len(nucleotide)

    #####滑动窗口计算change##################
    nucleotide_ = np.zeros([len(range(0, n - (window), step)), window], dtype=np.string_)
    R_z_ = np.zeros([len(range(0, n - (window), step)), window])
    R_f_ = np.zeros([len(range(0, n - (window), step)), window])
    p = np.zeros([len(range(0, n - (window), step)), ])
    j = 0
    for i in range(0, n - (window), step):
        nucleotide_[j, :] = nucleotide[i:i + window]
        R_z_[j, :] = R_z[i:i + window]
        R_f_[j, :] = R_f[i:i + window]
        p[j] = get_p(R_z[i:i + window], R_f[i:i + window],ratio_nan=ratio_nan,)
        j = j + 1
    p_=p.copy()
    if len(p[~np.isnan(p)])>0:
        a,p_bh,b,c =multi.multipletests(p[~np.isnan(p)],method='fdr_bh')
        p_[~np.isnan(p_)]=p_bh
    return p_,p


@jit(nopython=False)
def calculate_delta_gini(R_1, R_2, gini_1, gini_2, ratio_nan=0.9):
    '''Calculates Standard RMSD on two vectors of numbers of the same length'''
    # Check to see the vectors are of equal length.
    if len(R_1) != len(R_2):
        return np.nan
    else:
        R = R_1 - R_2
        if len(R[np.isnan(R) == False]) / len(R) <= ratio_nan:
            return np.nan
        else:
            delta_gini = gini_1 - gini_2
    return delta_gini

@jit(nopython=False)
def get_gff(X,chr,strand,data_gff_):
    gff=[]
    if len(data_gff_.loc[data_gff_['location'] == 'CDS', 'strat']) == 0:
        for i in range(len(X)):
            gff_='erro'
            gff.append(gff_)
    else:
        CDS_strat = list(data_gff_.loc[data_gff_['location'] == 'CDS', 'strat'])[0]
        CDS_end = list(data_gff_.loc[data_gff_['location'] == 'CDS', 'end'])[0]
        for i in range(len(X)):
            if strand == '+':
                if X[i]<CDS_strat:
                    gff_='five_prime_UTR'
                elif X[i]>CDS_end:
                    gff_ = 'three_prime_UTR'
                else:
                    gff_='CDS'
            else:
                if X[i]<CDS_strat:
                    gff_='three_prime_UTR'
                elif X[i]>CDS_end:
                    gff_ = 'five_prime_UTR'
                else:
                    gff_='CDS'
            gff.append(gff_)
    return gff


@jit(nopython=False)
def merge_gini(X1,X2,strat,end):
    #########读取数据#############
    name = X1[0]
    nucleotide = np.array(list(X1[1]))
    R_z = np.array(X1[6].split(',')).astype('float')
    R_z[R_z == -999] = np.nan
    R_f = np.array(X2[6].split(',')).astype('float')
    R_f[R_f == -999] = np.nan
    n = len(nucleotide)
    R_z_=R_z[strat:end+1]
    R_f_=R_f[strat:end+1]
    gini_z=get_gini(R_z_,nucleotide_=[],ratio_nan=0)
    gini_f = get_gini(R_f_, nucleotide_=[],ratio_nan=0)
    delta=gini_z-gini_f
    p=get_p(R_z_,R_f_,ratio_nan=0)
    R_f_mean = R_f_[R_f_>=0].mean()
    R_z_mean = R_z_[R_z_ >= 0].mean()
    return delta,p,R_f_mean,R_z_mean,gini_z,gini_f



@jit(nopython=False)
def get_statistics(data_z,data_f,data_location_,data_gff_,strand):
    data_z = np.array(data_z).reshape([7,])
    data_f = np.array(data_f).reshape([7,])
    chr=list(data_location_['chr'])[0]
    index_z,location_index_z,nucleotide_z, R_z, gini_z = get_window(data_z,data_location_)
    index_f,location_index_f,nucleotide_f, R_f, gini_f = get_window(data_f,data_location_)
    # print('ok')
    if len(R_z)!=len(R_f):
        print("error")
    else:
        location_exon = get_gff(location_index_z[:,0],chr,strand,data_gff_)
        # print('ok')
        location_exon = np.array(location_exon)
        p_bh, p = get_window_p(data_z, data_f)
        delta_gini = np.zeros([len(R_z), ])
        index_list = []
        nucleotide_list = []
        location_index_list = []
        for i in range(len(R_z)):
            index_list.append( ','.join(list(index_z[i,:].astype('int').astype('str'))))
            nucleotide_list.append(','.join(list(nucleotide_z[i, :].astype('str'))))
            location_index_list.append( ','.join(list(location_index_z[i,:].astype('int').astype('str'))))
            delta_gini[i] = calculate_delta_gini(R_z[i, :], R_f[i, :], gini_z[i], gini_f[i])
        index_list = np.array(index_list)
        location_index_list = np.array(location_index_list)
        nucleotide_list = np.array(nucleotide_list)
    return delta_gini,p_bh,p, gini_z, gini_f,index_list,location_index_list,nucleotide_list,location_exon


# @jit(nopython=False,error_model="numpy")
def merge_data(delta_gini,index_list,data_z,data_f,location_list,chr,strand,data_gff_,windows=50):
    data_z = np.array(data_z).reshape([7,])
    data_f = np.array(data_f).reshape([7,])
    nucleotide = np.array(list(data_z[1]))
    index = np.array([x.split(',')[0] for x in index_list]).astype('int')
    delta_gini_ =delta_gini.copy()
    delta_gini_[np.isnan(delta_gini_)]=0
    i=0
    windows_num=[]
    n=index[abs(delta_gini_)>0.1][0]+windows-1
    for j in range(len(index)-1):
        if abs(delta_gini_[j])<0.1:
            windows_num.append(np.nan)
        else:
            if index[j]>n+1:
                i=i+1
            n=index[j]+windows-1
            windows_num.append(i)
    if abs(delta_gini_[len(index)-1])<0.1:
        windows_num.append(np.nan)
    else:
        if index[j] > n + 1:
            i = i + 1
        windows_num.append(i)

    data_all= pd.DataFrame(np.zeros([len(set(windows_num) - set([np.nan])),9]))
    for i in list(set(windows_num) - set([np.nan])):
        windows_num_=np.array(windows_num)
        strat = index[windows_num_==i].astype('int')[0]
        end = index[windows_num_==i].astype('int')[-1]+windows-1
        strat_chr=location_list[strat]
        end_chr = location_list[end]
        delta,p=merge_gini(data_z,data_f,strat,end)
        location = get_gff(np.array([strat_chr,end_chr]),chr,strand,data_gff_)
        nucleotide_ =nucleotide[strat:end+1]
        nucleotide_str = ','.join(list(nucleotide_.astype('str')))
        data_all.iloc[i,:]=[strat,end,strat_chr,end_chr,location[0],location[1],delta,p,nucleotide_str]
    data_all.columns=['start','end','start_chr','end_chr','location_start','location_end','delta','p','nucleotide']
    return data_all


def merge_data_2(delta_gini,index_list,data_z,data_f,location_list,chr,strand,data_gff_,windows=50):
    data_z = np.array(data_z).reshape([7,])
    data_f = np.array(data_f).reshape([7,])
    nucleotide = np.array(list(data_z[1]))
    index = np.array([x.split(',')[0] for x in index_list]).astype('int')
    delta_gini_ =delta_gini.copy()
    delta_gini_[np.isnan(delta_gini_)]=0
    i=0
    windows_num=[]
    if abs(delta_gini_[0]) < 0.1:
        windows_num.append(np.nan)
    else:
        i = i + 1
        windows_num.append(i)

    for j in range(1,len(index)):
        if abs(delta_gini_[j])<0.1:
            windows_num.append(np.nan)
        else:
            if abs(delta_gini_[j-1])>0.1:
                windows_num.append(i)
            else:
                i=i+1
                windows_num.append(i)

    data_all= pd.DataFrame(np.zeros([len(set(windows_num) - set([np.nan])),13]))
    for i in list(set(windows_num) - set([np.nan])):
        windows_num_=np.array(windows_num)
        strat = index[windows_num_==i].astype('int')[0]
        end = index[windows_num_==i].astype('int')[-1]+windows-1
        delta, p, R_f_mean, R_z_mean,gini_z,gini_f=merge_gini(data_z,data_f,strat,end)
        if abs(delta)<=0.1:
            num = len(index[windows_num_==i].astype('int'))
            for n in range(1,num):
                end_=end-n
                delta, p, R_f_mean, R_z_mean,gini_z,gini_f = merge_gini(data_z, data_f, strat, end_)
                if abs(delta)>0.1:
                    end = end_
                    break
        strat_chr=location_list[strat]
        end_chr = location_list[end]
        location = get_gff(np.array([strat_chr,end_chr]),chr,strand,data_gff_)
        nucleotide_ =nucleotide[strat:end+1]
        nucleotide_str = ','.join(list(nucleotide_.astype('str')))
        data_all.iloc[i-1,:]=[strat,end,strat_chr,end_chr,location[0],location[1],delta,p,R_f_mean, R_z_mean,gini_z,gini_f,nucleotide_str]
    data_all.columns=['start','end','start_chr','end_chr','location_start','location_end','delta','p','R_f_mean','R_z_mean','gini_z','gini_f','nucleotide']
    return data_all



def main_sum(data_uv_z,data_uv_f,transcript_id_list,data_location,data_gff,save_path):
    i=0
    data_all = pd.DataFrame(columns=['transcript_id','start','end','start_chr','end_chr','location_start','location_end','delta','p','R_f_mean','R_z_mean','gini_z','gini_f','nucleotide'])
    for transcript in transcript_id_list:
        data_z = data_uv_z.loc[data_uv_z['transcript_id']==transcript,:]
        data_f = data_uv_f.loc[data_uv_f['transcript_id'] == transcript,:]
        data_location_ = data_location.loc[data_location['transcript_id'] == transcript,:]
        data_gff_ = data_gff.loc[data_gff['transcript_id']== transcript,:]
        location_list=np.array(list(data_location_['location'])[0].split(',')).astype('int')
        chr = list(data_location_['chr'])[0]
        strand = list(data_location_['strand'])[0]
        delta_gini, p_bh,p,gini_z, gini_f, index_list, location_index_list, nucleotide_list,location_exon = get_statistics(data_z,data_f,data_location_,data_gff_,strand)
        # print(delta_gini)

        data_new = pd.DataFrame(np.vstack([location_exon,index_list,location_index_list,nucleotide_list,gini_z, gini_f,delta_gini,p_bh,p])).T
        data_new.columns =['location_exon','index','location_index','nucleotide','gini_z', 'gini_f','delta_gini','p_bh','p']
        data_new['transcript']=transcript
        # pd.DataFrame(data_new).to_csv(save_path+'/'+transcript+'_gini.csv',sep='\t',header=True,index=False)
        data_all_=merge_data_2(delta_gini, index_list, data_z, data_f, location_list,chr,strand,data_gff_)
        data_all_['transcript_id']=transcript
        data_all = pd.concat([data_all,data_all_])
        i = i+1
        print(i)
    return data_all

####################主函数部分####################
if __name__ == '__main__':
    col_uv_f = '/data/TA_QUIZ_RNA_regulation/data/riboshape_liulab_batch4/final.modified_umodified/col_nouv/'
    col_uv_z = '/data/TA_QUIZ_RNA_regulation/data/riboshape_liulab_batch4/final.modified_umodified/col_uv/'
    uvr_uv_f = '/data/TA_QUIZ_RNA_regulation/data/riboshape_liulab_batch4/final.modified_umodified/uvr8_nouv/'
    uvr_uv_z = '/data/TA_QUIZ_RNA_regulation/data/riboshape_liulab_batch4/final.modified_umodified/uvr8_uv/'

    ###col###
    result = '/data/TA_QUIZ_RNA_regulation/result/PartIII.SHAPE-seq_analysis/col'
    uv_z = col_uv_z
    uv_f = col_uv_f

    data_uv_z = pd.read_csv(uv_z + '/final.modified_unmodified_new', sep='\t')
    data_uv_z.columns = ['transcript_id', 'Nucleotide', 'Modified_mutations', 'Modified_effective_depth',
                         'Untreated_mutations', 'Untreated_effective_depth', 'R1']
    data_uv_z = data_uv_z[['transcript_id', 'Nucleotide', 'Modified_mutations', 'Modified_effective_depth',
                         'Untreated_mutations', 'Untreated_effective_depth', 'R1']]
    data_uv_z = data_uv_z.drop_duplicates()

    data_uv_f = pd.read_csv(uv_f + '/final.modified_unmodified_new', sep='\t', header=None)
    data_uv_f.columns = ['transcript_id', 'Nucleotide', 'Modified_mutations', 'Modified_effective_depth',
                         'Untreated_mutations', 'Untreated_effective_depth', 'R1']
    data_uv_f = data_uv_f[['transcript_id', 'Nucleotide', 'Modified_mutations', 'Modified_effective_depth',
                         'Untreated_mutations', 'Untreated_effective_depth', 'R1']]
    data_uv_f = data_uv_f.drop_duplicates()
    statistics_sum=pd.read_csv(result + '/gini_summary_50_1.csv', sep='\t')
    # statistics_sum_new = pd.read_csv(result + '/summary_result_merge_50_1_0.csv',sep='\t')

    transcript_all = list(statistics_sum.loc[statistics_sum['num_0.1']>0,'transcript_id'])
    # transcript_all =transcript_all[:1000]
    # transcript_all=['AT5G26000.1','AT5G26000.2']
    # transcript_all = ['AT5G45260.2']
    print(len(transcript_all))
    gff_path = '/data/TA_QUIZ_RNA_regulation/data/ATH/GFF/Ath_genes.gff'
    data_gff = pd.read_csv(gff_path, sep='\t')
    # data_gff.columns = ['exosome','name','location','strat','end','.','+/-','num','id']
    data_location = pd.read_csv('/data/TA_QUIZ_RNA_regulation/data/ATH/GTF/shape_map/result/transcript_exon_location.csv', sep='\t')
    data_all = main_sum(data_uv_z,data_uv_f,transcript_all,data_location,data_gff,result+'/transcript_gini_merge')
    data_all.to_csv(result+'/summary_result_merge_50_1.csv',sep='\t',index=False)
    hit_level_coil_uv_f = pd.read_csv(uv_f+'cutoff.hit.group',sep='\t')
    #hit_level_coil_uv_f = hit_level_coil_uv_f.reset_index()
    hit_level_coil_uv_f.columns =['group','transcript_id','modified.median','unmodified.median','modified.sum','unmodified.sum','hit_f']
    hit_level_coil_uv_z = pd.read_csv(uv_z+'cutoff.hit.group',sep='\t')
    #hit_level_coil_uv_z = hit_level_coil_uv_z.reset_index()
    hit_level_coil_uv_z.columns =['group','transcript_id','modified.median','unmodified.median','modified.sum','unmodified.sum','hit_z']
    shape_data_coil =pd.merge(pd.merge(data_all,hit_level_coil_uv_f[['transcript_id','hit_f']],on='transcript_id',how='left'),hit_level_coil_uv_z[['transcript_id','hit_z']],on='transcript_id',how='left')
    shape_data_coil_2 =shape_data_coil.loc[(shape_data_coil['hit_f']>2)&(shape_data_coil['hit_z']>2),:]
    shape_data_coil.to_csv(result+'/summary_result_merge_50_1_new.csv',sep='\t',index=False)

    ###uvr8###
    result = '/data/TA_QUIZ_RNA_regulation/result/PartIII.SHAPE-seq_analysis/uvr8'
    uv_z = uvr_uv_z
    uv_f = uvr_uv_f

    data_uv_z = pd.read_csv(uv_z + '/final.modified_unmodified_new', sep='\t')
    data_uv_z.columns = ['transcript_id', 'Nucleotide', 'Modified_mutations', 'Modified_effective_depth',
                         'Untreated_mutations', 'Untreated_effective_depth', 'R1']
    data_uv_z = data_uv_z[['transcript_id', 'Nucleotide', 'Modified_mutations', 'Modified_effective_depth',
                           'Untreated_mutations', 'Untreated_effective_depth', 'R1']]
    data_uv_z = data_uv_z.drop_duplicates()

    data_uv_f = pd.read_csv(uv_f + '/final.modified_unmodified_new', sep='\t', header=None)
    data_uv_f.columns = ['transcript_id', 'Nucleotide', 'Modified_mutations', 'Modified_effective_depth',
                         'Untreated_mutations', 'Untreated_effective_depth', 'R1']
    data_uv_f = data_uv_f[['transcript_id', 'Nucleotide', 'Modified_mutations', 'Modified_effective_depth',
                           'Untreated_mutations', 'Untreated_effective_depth', 'R1']]
    data_uv_f = data_uv_f.drop_duplicates()
    statistics_sum = pd.read_csv(result + '/gini_summary_50_1.csv', sep='\t')
    # statistics_sum_new = pd.read_csv(result + '/summary_result_merge_50_1_0.csv',sep='\t')

    transcript_all = list(statistics_sum.loc[statistics_sum['num_0.1'] > 0, 'transcript_id'])
    # transcript_all =transcript_all[:1000]
    # transcript_all=['AT5G26000.1','AT5G26000.2']
    # transcript_all = ['AT5G45260.2']
    print(len(transcript_all))
    gff_path = '/data/TA_QUIZ_RNA_regulation/data/ATH/GFF/Ath_genes.gff'
    data_gff = pd.read_csv(gff_path, sep='\t')
    # data_gff.columns = ['exosome','name','location','strat','end','.','+/-','num','id']
    data_location = pd.read_csv(
        '/data/TA_QUIZ_RNA_regulation/data/ATH/GTF/shape_map/result/transcript_exon_location.csv', sep='\t')
    data_all = main_sum(data_uv_z, data_uv_f, transcript_all, data_location, data_gff,
                        result + '/transcript_gini_merge')
    data_all.to_csv(result + '/summary_result_merge_50_1.csv', sep='\t', index=False)
    hit_level_coil_uv_f = pd.read_csv(uv_f + 'cutoff.hit.group', sep='\t')
    #hit_level_coil_uv_f = hit_level_coil_uv_f.reset_index()
    hit_level_coil_uv_f.columns = ['group', 'transcript_id', 'modified.median', 'unmodified.median', 'modified.sum',
                                   'unmodified.sum', 'hit_f']
    hit_level_coil_uv_z = pd.read_csv(uv_z + 'cutoff.hit.group', sep='\t')
    #hit_level_coil_uv_z = hit_level_coil_uv_z.reset_index()
    hit_level_coil_uv_z.columns = ['group', 'transcript_id', 'modified.median', 'unmodified.median', 'modified.sum',
                                   'unmodified.sum', 'hit_z']
    shape_data_coil = pd.merge(
        pd.merge(data_all, hit_level_coil_uv_f[['transcript_id', 'hit_f']], on='transcript_id', how='left'),
        hit_level_coil_uv_z[['transcript_id', 'hit_z']], on='transcript_id', how='left')
    shape_data_coil_2 = shape_data_coil.loc[(shape_data_coil['hit_f'] > 2) & (shape_data_coil['hit_z'] > 2), :]
    shape_data_coil.to_csv(result + '/summary_result_merge_50_1_new.csv', sep='\t', index=False)
```

**output:**
`summary_result_merge_50_1.csv`：为hit level 大于0的转录本结果。
`summary_result_merge_50_1_new.csv`：为hit level 大于2的转录本结果。

## 2.h) 结果汇总
我们汇总结构改变区域的结构，包括该转录本是否差异表达，是否差异翻译，是否有剪切事件。
**input:**
`summary_result_merge_50_1_new.csv`：结构变化区域的汇总结果
`wt_rawdata.csv`：差异表达汇总结果
`wt.uvb-vs-nouvb.TE.csv`：差异翻译汇总结果
`result_splicing_WT.csv`：差异剪接事件汇总结果

脚本参考：`/data/TA_QUIZ_RNA_regulation/script/PartIII.SHAPE-seq_analysis/structure_changes_analysis/7.result_merge.py`

```python
import pandas as pd
import numpy as np
import re

exon_gtf_path ='/data/TA_QUIZ_RNA_regulation/data/ATH/GTF/shape_map/ath_exons.gtf'
gtf_data=pd.read_csv(exon_gtf_path,sep='\t',header=None)
gtf_data_new=pd.DataFrame(columns={'transcript_id','gene_type','gene_name','chr','+/-'})
gtf_data_new['+/-'] = gtf_data.iloc[:,6]
gtf_data_new['transcript_id'] = gtf_data.iloc[:,8].apply(lambda x:x.split(';')[1].split('"')[1])
gtf_data_new['gene_type'] = gtf_data.iloc[:,8].apply(lambda x:re.findall('gene_biotype ".*?"',x)[0].split('"')[1])
gtf_data_new['gene_name'] = gtf_data.iloc[:,8].apply(lambda x:re.findall('gene_name ".*?"',x)[0].split('"')[1] if 'gene_name' in x else np.nan)
gtf_data_new['chr'] = gtf_data.iloc[:,0].apply(lambda x: 6 if x=='Mt' else 7 if x=='Pt' else x ).astype('int')
gtf_data_new = gtf_data_new.drop_duplicates()
gtf_data_new.index = range(len(gtf_data_new))

#######################################WT########################################################
#######SHAPE############
shape_data_wt = pd.read_csv('/data/TA_QUIZ_RNA_regulation/result/PartIII.SHAPE-seq_analysis/col/summary_result_merge_50_1_new.csv',sep='\t')
shape_data_wt['up_down']= shape_data_wt['delta'].map(lambda x: 'up' if x>0 else 'down')
shape_data_wt['gene']=shape_data_wt['transcript_id'].map(lambda x:x.split('.')[0])

#####EXP############
exp_data_wt = pd.read_csv('/data/TA_QUIZ_RNA_regulation/result/PartI.RNA-seq_analysis/differential_expression/5.DESeq2/wt/wt_rawdata.csv',sep=',')
exp_data_wt =exp_data_wt.rename(columns={'Row.names':'gene','log2FoldChange':'log2FoldChange(EXP)','pvalue':'pvalue(EXP)','padj':'FDR(EXP)'})

#####riboshape############
ribo_data_wt = pd.read_csv('/data/TA_QUIZ_RNA_regulation/result/PartII.Ribo-seq_analysis/7.TE/wt.uvb-vs-nouvb.TE_new.csv',sep='\t')
ribo_data_wt = ribo_data_wt.reset_index()
ribo_data_wt =ribo_data_wt.rename(columns={'index':'gene','log2FC_TE_final':'log2FoldChange(TE)','pvalue_final':'pvalue(TE)','pvalue.adjust':'FDR(TE)'})

#####splicing############
splicing_data_wt = pd.read_csv('/data/TA_QUIZ_RNA_regulation/result/PartI.RNA-seq_analysis/differential_splicing/wt.UVB-vs-noUVB_filtered_psi_0.1/result_splicing_WT.csv',sep='\t')

####merge###########
data_wt = pd.merge(pd.merge(pd.merge(shape_data_wt[['gene','transcript_id','delta','start_chr','end_chr','start','end','location_start','location_end','hit_f', 'hit_z','R_f_mean', 'R_z_mean','gini_f', 'gini_z','up_down']]
                   ,exp_data_wt[['gene','log2FoldChange(EXP)','pvalue(EXP)', 'FDR(EXP)']],on='gene',how='left')
                   ,ribo_data_wt[['gene','log2FoldChange(TE)', 'pvalue(TE)', 'FDR(TE)']],on='gene',how='left')
                   ,splicing_data_wt[['gene','PSI(AS)', 'pvalue(AS)', 'FDR(AS)','AS_type']],on='gene',how='left')


#######################################UVR########################################################
#####SHAPE############
shape_data_uvr = pd.read_csv('/data/TA_QUIZ_RNA_regulation/result/PartIII.SHAPE-seq_analysis/uvr8/summary_result_merge_50_1_new.csv',sep='\t')
shape_data_uvr['up_down']= shape_data_uvr['delta'].map(lambda x: 'up' if x>0 else 'down')
shape_data_uvr['gene']=shape_data_uvr['transcript_id'].map(lambda x:x.split('.')[0])

#####EXP############
exp_data_uvr = pd.read_csv('/data/TA_QUIZ_RNA_regulation/result/PartI.RNA-seq_analysis/differential_expression/5.DESeq2/uvr8/uvr8_rawdata.csv',sep=',')
exp_data_uvr =exp_data_uvr.rename(columns={'Row.names':'gene','log2FoldChange':'log2FoldChange(EXP)','pvalue':'pvalue(EXP)','padj':'FDR(EXP)'})

#####riboshape############
ribo_data_uvr = pd.read_csv('/data/TA_QUIZ_RNA_regulation/result/PartII.Ribo-seq_analysis/7.TE/uvr8.uvb-vs-nouvb.TE.csv',sep='\t')
ribo_data_uvr = ribo_data_uvr.reset_index()
ribo_data_uvr =ribo_data_uvr.rename(columns={'index':'gene','log2FC_TE_final':'log2FoldChange(TE)','pvalue_final':'pvalue(TE)','pvalue.adjust':'FDR(TE)'})

#####splicing############
splicing_data_uvr = pd.read_csv('/data/TA_QUIZ_RNA_regulation/result/PartI.RNA-seq_analysis/differential_splicing/uvr8.UVB-vs-noUVB_filtered_psi_0.1/result_splicing_UVR8.csv',sep='\t')

####merge###########
data_uvr = pd.merge(pd.merge(pd.merge(shape_data_uvr[['gene','transcript_id','delta','start_chr','end_chr','start','end','location_start','location_end','hit_f', 'hit_z','R_f_mean', 'R_z_mean','gini_f', 'gini_z','up_down']]
                   ,exp_data_uvr[['gene','log2FoldChange(EXP)','pvalue(EXP)', 'FDR(EXP)']],on='gene',how='left')
                   ,ribo_data_uvr[['gene','log2FoldChange(TE)', 'pvalue(TE)', 'FDR(TE)']],on='gene',how='left')
                   ,splicing_data_uvr[['gene','PSI(AS)', 'pvalue(AS)', 'FDR(AS)','AS_type']],on='gene',how='left')

data_wt.to_csv('/data/TA_QUIZ_RNA_regulation/result/PartIII.SHAPE-seq_analysis/merge/merge_data_WT.csv',sep='\t',index=False)
data_uvr.to_csv('/data/TA_QUIZ_RNA_regulation/result/PartIII.SHAPE-seq_analysis/merge/merge_data_UVR8.csv',sep='\t',index=False)
data_wt.loc[(data_wt['hit_f']>2)&(data_wt['hit_z']>2),:].to_csv('/data/TA_QUIZ_RNA_regulation/result/PartIII.SHAPE-seq_analysis/merge/merge_data_WT_2.csv',sep='\t',index=False)
data_uvr.loc[(data_wt['hit_f']>2)&(data_wt['hit_z']>2),:].to_csv('/data/TA_QUIZ_RNA_regulation/result/PartIII.SHAPE-seq_analysis/merge/merge_data_UVR8_2.csv',sep='\t',index=False)
```

**output：**
`merge_data_WT.csv`：WT型，hit level>0的汇总结果
`merge_data_UVR8.csv`：UVR8型，hit level>0的汇总结果
`merge_data_WT_2.csv`：WT型，hit level>2的汇总结果
`merge_data_UVR8_2.csv`：UVR8型，hit level>2的汇总结果


# 3) motif analysis(选做/了解即可)
我们根据结构变化区域进行motif分析，主要分为三个部分：
- de novo sequence motif discovery
- known structure motif
- de novo structure motif

## 3.a) de novo sequence motif discovery
我们将结构变化区域分为上调(gini index of UV+>gini index of UV-)和下调(gini index of UV+<gini index of UV-)的两个部分进行motif富集。首先汇总结构变化区域的序列信息，然后利用MEME进行de novo sequence motif discovery。

### 3.a.1) 汇总结构变化区域的序列信息
这里我们将上调结构变化区域和下调结构变化区域分别汇总到一个文件，要求每个结构变化区域序列长度相同，我们取50nt固定长度，超过50nt的结构变化区域进行截断。

脚本参考：`/data/TA_QUIZ_RNA_regulation/script/PartIII.SHAPE-seq_analysis/motif/1.change_region_fasta.py`

```python
# de novo sequence motif discovery
# 将结构变化区域分为上调和下调的两个部分进行motif富集。首先汇总结构变化区域的序列信息，然后利用MEME进行de novo sequence motif discovery

import pandas as pd
path='/data/TA_QUIZ_RNA_regulation/result/PartIII.SHAPE-seq_analysis/3.3.motif/08.sequence/UVR8/'
## 相当于按delta的值，大于0的设置为up，小于0的设置为down。并筛选出其中hit_f,hit_z>2，传到data.change中
shape_data_coil = pd.read_csv('/data/TA_QUIZ_RNA_regulation/result/PartIII.SHAPE-seq_analysis/uvr8/summary_result_merge_50_1_new.csv',sep='\t')
shape_data_coil['up_down']= shape_data_coil['delta'].map(lambda x: 'up' if x>0 else 'down')
shape_data_coil['gene']= shape_data_coil['transcript_id'].map(lambda x:x.split('.')[0])
data_change =shape_data_coil.loc[(shape_data_coil['hit_f']>2)&(shape_data_coil['hit_z']>2),:]

## 将data_change按照up和down分开。
data_change_up = data_change.loc[data_change['up_down']=='up',:]
data_change_up.index=range(len(data_change_up))
data_change_down = data_change.loc[data_change['up_down']=='down',:]
data_change_down.index=range(len(data_change_down))

f = open(path + 'merge_up.fasta','w')
# 打开了一个文件 /data/TA_QUIZ_RNA_regulation/result/PartIII.SHAPE-seq_analysis/3.3.motif/08.sequence/UVR8/merge_up.fasta
# 'w' 打开一个文件只用于写入，如果该文件已存在则打开文件，并从开头开始编辑，即原有内容会被删除。如果该文件不存在，创建新文件。
for i in range(len(data_change_up)):
    transcript = data_change_up.loc[i,'transcript_id']
    # .loc[行名，列名]，行名为i，列名为transcript_i的项。
    location = data_change_up.loc[i, 'location_start']
    nucleotide = data_change_up.loc[i, 'nucleotide']
    nucleotide = nucleotide.replace(',', '')[:50]
    up_down = data_change_up.loc[i, 'up_down']
    start = data_change_up.loc[i, 'start']
    f.write('>' + transcript + '_' + location + '_' + str(int(start)) + '_' + up_down + '\n')
    #每两行为一端核苷酸长度为50nt的序列，第一行为转录本名称，所在部位及起始位置及上调还是下调
    #第二行为具体的核苷酸序列。
    f.write(nucleotide + '\n')
f.close()

# 同上，得到的是明显差异下调的序列
f = open(path + 'merge_down.fasta', 'w')
for i in range(len(data_change_down)):
    transcript = data_change_down.loc[i, 'transcript_id']
    location = data_change_down.loc[i, 'location_start']
    nucleotide = data_change_down.loc[i, 'nucleotide']
    nucleotide = nucleotide.replace(',', '')[:50]
    up_down = data_change_down.loc[i, 'up_down']
    start = data_change_down.loc[i, 'start']
    f.write('>' + transcript + '_' + location + '_' + str(int(start)) +'_'+ up_down+ '\n')
    f.write(nucleotide + '\n')

f.close()

# 同上，只不过现在计算的是WT而非UVR了
path = '/data/TA_QUIZ_RNA_regulation/result/PartIII.SHAPE-seq_analysis/3.3.motif/08.sequence/WT/'
shape_data_coil = pd.read_csv('/data/TA_QUIZ_RNA_regulation/result/PartIII.SHAPE-seq_analysis/col/summary_result_merge_50_1_new.csv',sep='\t')
shape_data_coil['up_down']= shape_data_coil['delta'].map(lambda x: 'up' if x>0 else 'down')
shape_data_coil['gene']= shape_data_coil['transcript_id'].map(lambda x:x.split('.')[0])
data_change =shape_data_coil.loc[(shape_data_coil['hit_f']>2)&(shape_data_coil['hit_z']>2),:]

data_change_up = data_change.loc[data_change['up_down']=='up',:]
data_change_up.index=range(len(data_change_up))
data_change_down = data_change.loc[data_change['up_down']=='down',:]
data_change_down.index=range(len(data_change_down))

f = open(path + 'merge_up.fasta', 'w')
for i in range(len(data_change_up)):
    transcript = data_change_up.loc[i, 'transcript_id']
    location = data_change_up.loc[i, 'location_start']
    nucleotide = data_change_up.loc[i, 'nucleotide']
    nucleotide = nucleotide.replace(',', '')[:50]
    up_down = data_change_up.loc[i, 'up_down']
    start = data_change_up.loc[i, 'start']
    f.write('>' + transcript + '_' + location + '_' + str(int(start)) +'_'+ up_down+ '\n')
    f.write(nucleotide + '\n')

f.close()

f = open(path + 'merge_down.fasta', 'w')
for i in range(len(data_change_down)):
    transcript = data_change_down.loc[i, 'transcript_id']
    location = data_change_down.loc[i, 'location_start']
    nucleotide = data_change_down.loc[i, 'nucleotide']
    nucleotide = nucleotide.replace(',', '')[:50]
    up_down = data_change_down.loc[i, 'up_down']
    start = data_change_down.loc[i, 'start']
    f.write('>' + transcript + '_' + location + '_' + str(int(start)) +'_'+ up_down+ '\n')
    f.write(nucleotide + '\n')

f.close()
```

**output:**
WT样本的:`merge_down.fasta`、`merge_up.fasta`
UVR8突变体的:`merge_down.fasta`、`merge_up.fasta`

![10.4.6.merge_down.fasta.png](../../.gitbook/assets/10.4.6.merge_down.fasta.png)

第一行为转录本名称，所在部位、起始位置、上调还是下调。
第二行为具体的核苷酸序列。

### 3.a.2）MEME在线motif分析
我们利用MEME在线版进行motif分析，地址为[http://meme-suite.org/tools/meme](http://meme-suite.org/tools/meme)。
使用MEME的Motif discovery工具，预测输入序列上的motif信息

![10.4.7.MEME.png](../../.gitbook/assets/10.4.7.MEME.png)

我们上传merge_down.fasta、merge_up.fasta后提交运算。

![10.4.8.discovery_motif.png](../../.gitbook/assets/10.4.8.discovery_motif.png)

在得到的结果中我们根据Evalue进行筛选，Sites为找到的motif在几个结构改变区域中出现过，由此可以计算出该motif的覆盖率。我们选择覆盖率较高的motif，这里取大于4%。

## 3.b) Know structure motif
我们将结构变化区域序列与Rfam中已知的motif进行对比，找到结构变化区域中包含的已知的结构motif。

### 3.b.1) 获取结构变化区域序列
我们将同一转录本结构改变区域的序列放在一个文件中，不需要对序列长度进行处理。
参考脚本：`/data/TA_QUIZ_RNA_regulation/script/PartIII.SHAPE-seq_analysis/motif/2.known_structure_motif.py`

```python
import pandas as pd

path='/data/TA_QUIZ_RNA_regulation/result/PartIII.SHAPE-seq_analysis/3.3.motif/01.fasta/UVR8/'
data = pd.read_csv('/data/TA_QUIZ_RNA_regulation/data/riboshape_liulab_batch4/final.modified_umodified/uvr8_nouv/final.modified_unmodified', sep='\t')
shape_data_coil = pd.read_csv('/data/TA_QUIZ_RNA_regulation/result/PartIII.SHAPE-seq_analysis/uvr8/summary_result_merge_50_1.csv',sep='\t')
shape_data_coil['up_down']= shape_data_coil['delta'].map(lambda x: 'up' if x>0 else 'down')
hit_level_coil_uv_f = pd.read_csv('/data/TA_QUIZ_RNA_regulation/data/riboshape_liulab_batch4/final.modified_umodified/uvr8_nouv/cutoff.hit.group',sep='\t')
# hit_level_coil_uv_f = hit_level_coil_uv_f.reset_index()
hit_level_coil_uv_f.columns =['group','transcript_id','modified.median','unmodified.median','modified.sum','unmodified.sum','hit_f']
hit_level_coil_uv_z = pd.read_csv('/data/TA_QUIZ_RNA_regulation/data/riboshape_liulab_batch4/final.modified_umodified/uvr8_uv/cutoff.hit.group',sep='\t')
# hit_level_coil_uv_z = hit_level_coil_uv_z.reset_index()
hit_level_coil_uv_z.columns =['group','transcript_id','modified.median','unmodified.median','modified.sum','unmodified.sum','hit_z']
shape_data_coil = pd.merge(pd.merge(shape_data_coil,hit_level_coil_uv_f[['transcript_id','hit_f']],on='transcript_id',how='left'),hit_level_coil_uv_z[['transcript_id','hit_z']],on='transcript_id',how='left')
shape_data_coil['gene']=shape_data_coil['transcript_id'].map(lambda x:x.split('.')[0])
data_change = shape_data_coil.loc[(shape_data_coil['hit_f']>2)&(shape_data_coil['hit_z']>2),:]
UVR_only_sample = pd.read_csv('/data/TA_QUIZ_RNA_regulation/data/shape_only/Uvr_shape_only.txt',sep='\t',header=None)
data_change = data_change.loc[data_change.gene.isin(list(UVR_only_sample[0])),:]
data_change.index=range(len(data_change))

for i in range(len(data_change)):
    if i%1000 == 0:
        print(i)
    transcript = data_change.loc[i,'transcript_id']
    Nucleotide = list(list(data.loc[data['transcript_id']==transcript,'Nucleotide'])[0])
    start = int(data_change.loc[i,'start'])
    end = int(data_change.loc[i,'end'])
    data_change.loc[i,'nucleotide']=''.join(Nucleotide[start:end+1])

data_5UTR=data_change
transcript_list=list(set(data_5UTR['transcript_id']))
gene = [i.split('.')[0] for i in transcript_list]

for i in range(len(transcript_list)):
    transcript=transcript_list[i]
    f = open(path + transcript + '.fasta', 'w')
    data_5UTR_ = data_5UTR.loc[data_5UTR['transcript_id'] == transcript, :]
    data_5UTR_.index = range(len(data_5UTR_))
    for j in range(len(data_5UTR_)):
        start = int(data_5UTR_.loc[j, 'start'])
        end = int(data_5UTR_.loc[j, 'end'])
        location = data_5UTR_.loc[j, 'location_start']
        up_down = data_5UTR_.loc[j, 'up_down']
        seq = data_5UTR_.loc[j, 'nucleotide']
        f.write('>' + transcript + '_' + location + '_' + up_down + '_' + str(start) + '_' + str(end) + '\n')
        f.write(seq + '\n')
    f.close


path = '/data/TA_QUIZ_RNA_regulation/result/PartIII.SHAPE-seq_analysis/3.3.motif/01.fasta/WT/'
data = pd.read_csv('/data/TA_QUIZ_RNA_regulation/data/riboshape_liulab_batch4/final.modified_umodified/col_nouv/final.modified_unmodified',sep='\t')
shape_data_coil = pd.read_csv('/data/TA_QUIZ_RNA_regulation/result/PartIII.SHAPE-seq_analysis/col/summary_result_merge_50_1.csv',sep='\t')
shape_data_coil['up_down']= shape_data_coil['delta'].map(lambda x: 'up' if x>0 else 'down')
hit_level_coil_uv_f = pd.read_csv('/data/TA_QUIZ_RNA_regulation/data/riboshape_liulab_batch4/final.modified_umodified/col_nouv/cutoff.hit.group',sep='\t')
# hit_level_coil_uv_f = hit_level_coil_uv_f.reset_index()
hit_level_coil_uv_f.columns =['group','transcript_id','modified.median','unmodified.median','modified.sum','unmodified.sum','hit_f']
hit_level_coil_uv_z = pd.read_csv('/data/TA_QUIZ_RNA_regulation/data/riboshape_liulab_batch4/final.modified_umodified/col_uv/cutoff.hit.group',sep='\t')
# hit_level_coil_uv_z = hit_level_coil_uv_z.reset_index()
hit_level_coil_uv_z.columns =['group','transcript_id','modified.median','unmodified.median','modified.sum','unmodified.sum','hit_z']
shape_data_coil =pd.merge(pd.merge(shape_data_coil,hit_level_coil_uv_f[['transcript_id','hit_f']],on='transcript_id',how='left'),hit_level_coil_uv_z[['transcript_id','hit_z']],on='transcript_id',how='left')
shape_data_coil['gene']= shape_data_coil['transcript_id'].map(lambda x:x.split('.')[0])
data_change =shape_data_coil.loc[(shape_data_coil['hit_f']>2)&(shape_data_coil['hit_z']>2),:]
UVR_only_sample = pd.read_csv('/data/TA_QUIZ_RNA_regulation/data/shape_only/WT_shape_only.txt',sep='\t',header=None)
data_change = data_change.loc[data_change.gene.isin(list(UVR_only_sample[0])),:]
data_change.index=range(len(data_change))

for i in range(len(data_change)):
    if i%1000==0:
        print(i)
    trasnscript = data_change.loc[i,'transcript_id']
    Nucleotide = list(list(data.loc[data['transcript_id'] == trasnscript, 'Nucleotide'])[0])
    start = int(data_change.loc[i,'start'])
    end = int(data_change.loc[i, 'end'])
    data_change.loc[i,'nucleotide']=''.join(Nucleotide[start:end+1])


data_5UTR=data_change
trasnscript_list=list(set(data_5UTR['transcript_id']))
gene = [i.split('.')[0] for i in trasnscript_list]

for i in range(len(trasnscript_list)):
    trasnscript=trasnscript_list[i]
    f = open(path+trasnscript+'.fasta','w')
    data_5UTR_ = data_5UTR.loc[data_5UTR['transcript_id']==trasnscript,:]
    data_5UTR_.index=range(len(data_5UTR_))
    for j in range(len(data_5UTR_)):
        start = int(data_5UTR_.loc[j, 'start'])
        end = int(data_5UTR_.loc[j, 'end'])
        location = data_5UTR_.loc[j, 'location_start']
        up_down =data_5UTR_.loc[j, 'up_down']
        seq = data_5UTR_.loc[j, 'nucleotide']
        f.write('>'+trasnscript+'_'+location+'_'+up_down+'_'+str(start)+'_'+str(end)+'\n')
        f.write(seq+'\n')
    f.close()
```

**output：**
输出目录`01.fasta/WT`,`01.fasta/UVR8`，文件中每个转录本为一个文件。

### 3.b.2) cmscan软件对比Rfam
每个转录本都需要进行Rfam,单个转录本进行对比的脚本如下：
以`AT3G16830.1`为例
```linux
#!/bin/sh
#SBATCH -J AT3G16830.1
#SBATCH -p CN_BIOT
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --output=AT3G16830.1.out
#SBATCH --error=AT3G16830.1.err
#source ~/.bashrc

export PATH=/data/zhaoyizi/software/anaconda3/envs/Riboshape/bin:/data/zhaoyizi/software/infernal_bin/bin:/data/zhaoyizi/software:$PATH

echo AT3G16830.1 start `date`
/data/zhaoyizi/software/infernal_bin/bin/cmscan \
/data/zhaoyizi/software/Rfam.cm \
/data/TA_QUIZ_RNA_regulation/result/PartIII.SHAPE-seq_analysis/3.3.motif/01.fasta/WT/AT1G23730.2.fasta > /data/TA_QUIZ_RNA_regulation/result/PartIII.SHAPE-seq_analysis/3.3.motif/05.cmscan/WT/AT1G23730.2.cmscan.txt
echo AT1G23730.2 end `date`
```
为了方便，我们写脚本批量生成所有转录本的比对脚本
脚本参考:`/data/TA_QUIZ_RNA_regulation/script/PartIII.SHAPE-seq_analysis/motif/3.cmscan_batch.sh`
```linux
#!/bin/sh
#SBATCH -J PartIII.SHAPE-seq_analysis
#SBATCH -p CN_BIOT
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --output=3.cmscan_batch.out
#SBATCH --error=3.cmscan_batch.err
#source ~/.bashrc


export PATH=/data/zhaoyizi/software/anaconda3/envs/Riboshape/bin:$PATH

mkdir -p /data/TA_QUIZ_RNA_regulation/result/PartIII.SHAPE-seq_analysis/3.3.motif/05.cmscan/WT
python /data/TA_QUIZ_RNA_regulation/data/script/PartIII.SHAPE-seq_analysis/3.cmscan.py \
/data/TA_QUIZ_RNA_regulation/result/PartIII.SHAPE-seq_analysis/3.3.motif/01.fasta/WT \
/data/TA_QUIZ_RNA_regulation/result/PartIII.SHAPE-seq_analysis/3.3.motif/05.cmscan/WT \
/data/TA_QUIZ_RNA_regulation/result/PartIII.SHAPE-seq_analysis/3.3.motif/05.cmscan/WT.cmscan.sh

mkdir -p /data/TA_QUIZ_RNA_regulation/result/PartIII.SHAPE-seq_analysis/3.3.motif/05.cmscan/UVR
python /data/TA_QUIZ_RNA_regulation/data/script/PartIII.SHAPE-seq_analysis/3.cmscan.py \
/data/TA_QUIZ_RNA_regulation/result/PartIII.SHAPE-seq_analysis/3.3.motif/01.fasta/UVR \
/data/TA_QUIZ_RNA_regulation/result/PartIII.SHAPE-seq_analysis/3.3.motif/05.cmscan/UVR8 \
/data/TA_QUIZ_RNA_regulation/result/PartIII.SHAPE-seq_analysis/3.3.motif/05.cmscan/UVR.cmscan.sh
```

其中`3.cmscan.py`为功能脚本，后面需要三个输入参数，分别是数据输入目录路径、输出目录路径和生成脚本的路径及名称。`3.cmscan.py`脚本如下。
```linux
import os,re,sys
import pandas as pd
import numpy as np

inputDir = sys.argv[1]
outputDir = sys.argv[2]
script = sys.argv[3]

SCRIPT = open(script,"w")
SCRIPT.write("#!/bin/sh\n")
SCRIPT.write("#SBATCH -J PartIII.SHAPE-seq_analysis\n")
SCRIPT.write("#SBATCH -p CN_BIOT\n")
SCRIPT.write("#SBATCH --nodes=1\n")
SCRIPT.write("#SBATCH --ntasks=1\n")
SCRIPT.write("#SBATCH --output=%j.out\n")
SCRIPT.write("#SBATCH --error=%j.err\n")
SCRIPT.write("#source ~/.bashrc\n")
SCRIPT.write("export PATH=/data/zhaoyizi/software/anaconda3/envs/Riboshape/bin:/data/zhaoyizi/software/infernal_bin/bin:/data/zhaoyizi/software:$PATH\n")

for root,dirs,files in os.walk(inputDir):
        for file_name in files:
                if re.search("(.*).fasta",file_name):
                        file = os.path.join(root,file_name)
                        transcript_ID = re.search("(.*).fasta",file_name).group(1)
                        output_file = os.path.join(outputDir,transcript_ID+".cmscan.txt")
                        SCRIPT.write("echo "+transcript_ID+" start `date`\n")
                        #SCRIPT.write("/BioII/lulab_b/chenyinghui/software/infernal/infernal-1.1.2/src/cmscan /BioII/lulab_b/chenyinghui/database/Rfam/14.1/Rfam.cm "+file+" > "+output_file+"\n")
                        SCRIPT.write("/data/zhaoyizi/software/infernal_bin/bin/cmscan /data/zhaoyizi/software/Rfam.cm " +file+" > "+output_file+"\n")
                        SCRIPT.write("echo "+transcript_ID+" end `date`\n")
SCRIPT.close()
```

**output：**
`WT.cmscan.sh`和`UVR.cmscan.sh`为所有转录本比对脚本，分别运行两个脚本即可完成与Rfam的比对,得到输出文件夹`~/05.cmscan/WT`和`~/05.cmscan/UVR8`

### 3.b.3) 过滤并汇总结果
过滤Evalue<1的motif，并进行汇总
脚本参考：`/data/TA_QUIZ_RNA_regulation/script/PartIII.SHAPE-seq_analysis/motif/4.motif_filter.sh`

```linux
#!/bin/bash
#SBATCH -J PartIII.SHAPE-seq_analysis
#SBATCH -p CN_BIOT
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --output=4.motif_filter.out
#SBATCH --error=4.motif_filter.err
#source ~/.bashrc

export PATH=/data/zhaoyizi/software/anaconda3/envs/Riboshape/bin:$PATH

echo start `date`

mkdir -p /data/zhaoyizi/3.3.motif/data/06.motif_alignment_info/WT
python /data/TA_QUIZ_RNA_regulation/data/script/PartIII.SHAPE-seq_analysis/4.cmscan_result_summary.sig.py \
/data/TA_QUIZ_RNA_regulation/result/PartIII.SHAPE-seq_analysis/3.3.motif/05.cmscan/WT \
/data/TA_QUIZ_RNA_regulation/result/PartIII.SHAPE-seq_analysis/3.3.motif/06.motif_alignment_info/WT \
/data/TA_QUIZ_RNA_regulation/result/PartIII.SHAPE-seq_analysis/3.3.motif/06.motif_alignment_info/WT.motif.stat.xls

mkdir -p /data/zhaoyizi/3.3.motif/data/06.motif_alignment_info/UVR8
python /data/TA_QUIZ_RNA_regulation/data/script/PartIII.SHAPE-seq_analysis/4.cmscan_result_summary.sig.py \
/data/TA_QUIZ_RNA_regulation/result/PartIII.SHAPE-seq_analysis/3.3.motif/05.cmscan/UVR8 \
/data/TA_QUIZ_RNA_regulation/result/PartIII.SHAPE-seq_analysis/3.3.motif/06.motif_alignment_info/UVR8 \
/data/TA_QUIZ_RNA_regulation/result/PartIII.SHAPE-seq_analysis/3.3.motif/06.motif_alignment_info/UVR8.motif.stat.xls

echo finish `date`
```

**output：**
输出文件`~/WT.motif.stat.xls`和`~/UVR.motif.stat.xls`包括Evalue>1的motif及其所在的转录本。

## 3.c) de novo structure motif
利用BEAR软件找到的de novo structure motif,将结构变化区域分为上调(gini index of UV+>gini index of UV-)和下调(gini index of UV+<gini index of UV-)的两个部分进行motif富集
主要步骤分为：
- 1. 获得上调结构区域和下调结构区域对应的序列、UV-下reactivity、UV+下reactivity。
- 2. 基于RNAstructure，根据每个结构变化区域的序列和reactivity进行RNA二级结构预测，得到结构的DB文件。
- 3. 将所有上调结构变化区域的序列和结构DB汇总到一个文件上，将所有下调结构变化区域的序列和结构DB汇总到一个文件上。
- 4. 基于BEAM进行motif 富集分析。
- 5. 对找到的结构motif进行过滤和汇总。

### 3.c.1) 结构区域对应的序列和reactivity
脚本参考：`/data/TA_QUIZ_RNA_regulation/script/PartIII.SHAPE-seq_analysis/motif/5.1.get_fasta_shape.py`

```python
import pandas as pd

col_uv_f='/data/TA_QUIZ_RNA_regulation/data/riboshape_liulab_batch4/final.modified_umodified/col_nouv/'
col_uv_z='/data/TA_QUIZ_RNA_regulation/data/riboshape_liulab_batch4/final.modified_umodified/col_uv/'
uvr_uv_f='/data/TA_QUIZ_RNA_regulation/data/riboshape_liulab_batch4/final.modified_umodified/uvr8_nouv/'
uvr_uv_z='/data/TA_QUIZ_RNA_regulation/data/riboshape_liulab_batch4/final.modified_umodified/uvr8_uv/'

data_1 = pd.read_csv(col_uv_z + '/final.modified_unmodified_new', sep='\t', header=None)
data_1.columns =['transcript_id', 'Nucleotide', 'Modified_mutations', 'Modified_effective_depth','Untreated_mutations', 'Untreated_effective_depth', 'R']
data_1 = data_1[['transcript_id', 'Nucleotide', 'Modified_mutations', 'Modified_effective_depth','Untreated_mutations', 'Untreated_effective_depth', 'R']]
data_1 = data_1.drop_duplicates()

data_2=pd.read_csv(col_uv_f + '/final.modified_unmodified_new', sep='\t', header=None)
data_2.columns = ['transcript_id', 'Nucleotide', 'Modified_mutations', 'Modified_effective_depth','Untreated_mutations', 'Untreated_effective_depth', 'R']
data_2 = data_2[['transcript_id', 'Nucleotide', 'Modified_mutations', 'Modified_effective_depth','Untreated_mutations', 'Untreated_effective_depth', 'R']]
data_2 = data_2.drop_duplicates()

data_3 = pd.read_csv(uvr_uv_z + '/final.modified_unmodified_new', sep='\t', header=None)
data_3.columns = ['transcript_id', 'Nucleotide', 'Modified_mutations', 'Modified_effective_depth','Untreated_mutations', 'Untreated_effective_depth', 'R']
data_3 = data_3[['transcript_id', 'Nucleotide', 'Modified_mutations', 'Modified_effective_depth','Untreated_mutations', 'Untreated_effective_depth', 'R']]
data_3 = data_3.drop_duplicates()

data_4 = pd.read_csv(uvr_uv_f + '/final.modified_unmodified_new', sep='\t', header=None)
data_4.columns = ['transcript_id', 'Nucleotide', 'Modified_mutations', 'Modified_effective_depth','Untreated_mutations', 'Untreated_effective_depth', 'R']
data_4 = data_4[['transcript_id', 'Nucleotide', 'Modified_mutations', 'Modified_effective_depth','Untreated_mutations', 'Untreated_effective_depth', 'R']]
data_4 = data_4.drop_duplicates()

##### WT #####
path='/data/TA_QUIZ_RNA_regulation/result/PartIII.SHAPE-seq_analysis/col/'
savepath='/data/TA_QUIZ_RNA_regulation/result/PartIII.SHAPE-seq_analysis/3.3.motif/07.BEAM/WT/'
shape_data_coil=pd.read_csv(path+'/summary_result_merge_50_1_new.csv',sep='\t')
shape_data_coil = shape_data_coil.loc[(shape_data_coil['hit_f']>2)&(shape_data_coil['hit_z']>2),:]
shape_data_coil.index=range(len(shape_data_coil))
shape_data_coil['up_down']= shape_data_coil['delta'].map(lambda x: 'up' if x>0 else 'down')
# shape_data_coil = shape_data_coil.reset_index()
transcript_list=list(shape_data_coil['transcript_id'])
shape_data_coil[['index','transcript_id','up_down','start']].to_csv(savepath+'index.csv',sep='\t',index=False)

for i in range(len(transcript_list)):
    print(i)
    transcript=transcript_list[i]
    f = open(savepath+'/fasta/'+str(i)+'.fasta','w')
    start = int(shape_data_coil.loc[i,'start'])
    end = int(shape_data_coil.loc[i,'end'])
    location = shape_data_coil.loc[i, 'location_start']
    up_down = shape_data_coil.loc[i, 'up_down']
    seq = shape_data_coil.loc[i, 'nucleotide'].replace(',', '')
    f.write('>' + transcript + '_' + location + '_' + up_down + '_' + str(start) + '_' + str(end) + '\n')
    f.write(seq + '\n')
    f.close()
    data_z_ = list(data_1.loc[data_1['transcript_id'] == transcript, 'R'])[0].split(',')[start:end]
    data_f_ = list(data_2.loc[data_2['transcript_id'] == transcript, 'R'])[0].split(',')[start:end]
    data_z = pd.DataFrame(columns=['Positions', 'Reacivity'])
    data_z['Positions'] = [j for j in range(len(data_z_))]
    data_z['Reactivity'] = data_z_
    data_z = data_z[['Positions', 'Reactivity']]
    data_f = pd.DataFrame(columns=['Positions', 'Reactivity'])
    data_f['Positions'] = [j for j in range(len(data_f_))]
    data_f['Reactivity'] = data_f_
    data_f = data_f[['Positions', 'Reactivity']]
    data_z.to_csv(savepath + '/shape/' + str(i) + '_z.shape', sep='\t', index=False, header=None)
    data_f.to_csv(savepath + '/shape/' + str(i) + '_f.shape', sep='\t', index=False, header=None)


##### UVR8 #####
path = '/data/TA_QUIZ_RNA_regulation/result/PartIII.SHAPE-seq_analysis/uvr8/'
savepath = '/data/TA_QUIZ_RNA_regulation/result/PartIII.SHAPE-seq_analysis/3.3.motif/07.BEAM/UVR8/'
shape_data_coil = pd.read_csv(path+'/summary_result_merge_50_1_new.csv',sep='\t')
shape_data_coil = shape_data_coil.loc[(shape_data_coil['hit_f']>2)&(shape_data_coil['hit_z']>2),:]
shape_data_coil.index=range(len(shape_data_coil))
shape_data_coil['up_down']= shape_data_coil['delta'].map(lambda x: 'up' if x>0 else 'down')
# shape_data_coil = shape_data_coil.reset_index()
transcript_list=list(shape_data_coil['transcript_id'])
shape_data_coil[['index','transcript_id','up_down','start']].to_csv(savepath+'index.csv',sep='\t',index=False)
print(shape_data_coil)

for i in range(len(transcript_list)):
    print(i)
    transcript=transcript_list[i]
    f = open(savepath+'/fasta/'+str(i)+'.fasta','w')
    start = int(shape_data_coil.loc[i,'start'])
    end = int(shape_data_coil.loc[i,'end'])
    location = shape_data_coil.loc[i, 'location_start']
    up_down = shape_data_coil.loc[i, 'up_down']
    seq = shape_data_coil.loc[i, 'nucleotide'].replace(',', '')
    f.write('>' + transcript + '_' + location + '_' + up_down + '_' + str(start) + '_' + str(end) + '\n')
    f.write(seq + '\n')
    f.close()
    data_z_ = list(data_3.loc[data_3['transcript_id'] == transcript, 'R'])[0].split(',')[start:end]
    data_f_ = list(data_4.loc[data_4['transcript_id'] == transcript, 'R'])[0].split(',')[start:end]
    data_z =pd.DataFrame(columns=['Positions','Reacivity'])
    data_z['Positions'] = [j for j in range(len(data_z_))]
    data_z['Reactivity'] = data_z_
    data_z = data_z[['Positions', 'Reactivity']]
    data_f = pd.DataFrame(columns=['Positions', 'Reactivity'])
    data_f['Positions'] = [j for j in range(len(data_f_))]
    data_f['Reactivity'] = data_f_
    data_f = data_f[['Positions', 'Reactivity']]
    data_z.to_csv(savepath + '/shape/' + str(i) + '_z.shape', sep='\t', index=False, header=None)
    data_f.to_csv(savepath + '/shape/' + str(i) + '_f.shape', sep='\t', index=False, header=None)
```

**output：**
得到`fasta`和`shape`两个输出目录以及`index.csv`输出文件。
`fasta`文件夹中，每一个文件为一个转录本。
第一行为转录本id、所在位置、上调还是下调、起始坐标。
第二行为具体核苷酸信息。

`shape`文件夹中，每一个文件为与fasta文件对应的结构信息。
`f`代表不加光照;`z`代表被光照。

`index.csv`为目录索引文件，方便查询。

### 3.c.2) 预测RNA structure
 我们利用`RNAstructure`软件中的`ShapeKnots`和`ct2dot`函数，基于结构改变区域的序列和reactivity文件进行二级结构预测。
 
 ```linux
 #!/bin/sh
#SBATCH -J PartIII.SHAPE-seq_analysis
#SBATCH -p CN_BIOT
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --output=5.2.RNA_structure_predicted.out
#SBATCH --error=5.2.RNA_structure_predicted.err
#source ~/.bashrc

export PATH=/data/zhaoyizi/software/anaconda3/envs/Riboshape/bin:$PATH
# 添加热力学参数的路径
export DATAPATH=/data/zhaoyizi/software/anaconda3/envs/Riboshape/share/rnastructure/data_tables/

#####WT#####
cd /data/TA_QUIZ_RNA_regulation/result/PartIII.SHAPE-seq_analysis/3.3.motif/07.BEAM/WT
mkdir -p dot
fasta=$"/data/TA_QUIZ_RNA_regulation/result/PartIII.SHAPE-seq_analysis/3.3.motif/07.BEAM/WT/fasta/"
shape=$"/data/TA_QUIZ_RNA_regulation/result/PartIII.SHAPE-seq_analysis/3.3.motif/07.BEAM/WT/shape/"
ct_path=$"/data/TA_QUIZ_RNA_regulation/result/PartIII.SHAPE-seq_analysis/3.3.motif/07.BEAM/WT/dot/"
# 获得每个结构变化区域在UV-和UV+下的RNA二级结构，输出文件
for i in `seq 0 1580`
# 这里的1580是由具体情况决定的，这里共产生了1580个fasta文件，因此是1580
do
echo $i
ShapeKnots $fasta/"$i".fasta $ct_path/"$i"_z.ct -sh $shape/"$i"_z.shape -sm 1.9 -si -0.7
ShapeKnots $fasta/"$i".fasta $ct_path/"$i"_f.ct -sh $shape/"$i"_f.shape -sm 1.9 -si -0.7
ct2dot $ct_path/"$i"_z.ct 1 $ct_path/"$i"_z.dt
ct2dot $ct_path/"$i"_f.ct 1 $ct_path/"$i"_f.dt
done

cd /data/TA_QUIZ_RNA_regulation/result/PartIII.SHAPE-seq_analysis/3.3.motif/07.BEAM/UVR8/
mkdir -p dot
fasta=$"/data/TA_QUIZ_RNA_regulation/result/PartIII.SHAPE-seq_analysis/3.3.motif/07.BEAM/UVR8/fasta/"
shape=$"/data/TA_QUIZ_RNA_regulation/result/PartIII.SHAPE-seq_analysis/3.3.motif/07.BEAM/UVR8/shape/"
ct_path=$"/data/TA_QUIZ_RNA_regulation/result/PartIII.SHAPE-seq_analysis/3.3.motif/07.BEAM/UVR8/dot/"

for i in `seq 0 1313`
# 同上
do
echo $i
ShapeKnots $fasta/"$i".fasta $ct_path/"$i"_z.ct -sh $shape/"$i"_z.shape -sm 1.9 -si -0.7
ShapeKnots $fasta/"$i".fasta $ct_path/"$i"_f.ct -sh $shape/"$i"_f.shape -sm 1.9 -si -0.7
ct2dot $ct_path/"$i"_z.ct 1 $ct_path/"$i"_z.dt
ct2dot $ct_path/"$i"_f.ct 1 $ct_path/"$i"_f.dt
done
 ```
 
我们可以获得每个结构变化区域在UV-和UV+下的RNA二级结构，输出文件夹为`/data/TA_QUIZ_RNA_regulation/result/PartIII.SHAPE-seq_analysis/3.3.motif/07.BEAM/WT/dot`。其中`_f.ct`为UV-下的RNA结构，`_z.ct`为UV+下的RNA结构。

### 3.c.3) 数据汇总
我们对UV-下上调的结构变化区域和下调结构变化区域进行motif富集分析，需要对UV-下上调/下调的所有结构变化区域进行汇总。

脚本参考：`/data/TA_QUIZ_RNA_regulation/script/PartIII.SHAPE-seq_analysis/motif/6.merge.py`
 ```
 import pandas as pd
path='/data/TA_QUIZ_RNA_regulation/result/PartIII.SHAPE-seq_analysis/3.3.motif/07.BEAM/WT/'
index = pd.read_csv(path+'index.csv',sep='\t')
## 相当于是提取了   index中，up_down的值为up的那一行的index值。
index_up = list(index.loc[index['up_down']=='up','index'])
index_down = list(index.loc[index['up_down']=='down','index'])

# 按结构变化区域的上调和下调以及是否加光，将文件合并。
f_new = open(path + 'merge_up_nouv.fa','w')
# 打开一个名叫 merge_up_nouv.fa 的可写文件，如果不存在就创一个
for i in range(len(index_up)):
    # 以上调的结构变化区域的数量进行循环
    try:
        f = open(path + '/dot/' + str(index_up[i]) + '_f.dt')
        # 打开一个已经存在_f.dt文件
        j = 0
        for line in f:
            # 对于f的每一行进行操作
            if j == 0:
                transcript = line.split(' ')[-1]
                # 第一行提取转录本信息
                f_new.write('>' + transcript)
            else:
                # 对除第一行外的其他行，全部保留
                f_new.write(line)
            j = j + 1
        f.close()
    except FileNotFoundError:
        # 如果文件不存在就 跳过
        continue
f_new.close()

f_new=open(path + 'merge_down_nouv.fa','w')
for i in range(len(index_down)):
    try:
        f = open(path + '/dot/' + str(index_down[i])+'_f.dt')
        j = 0
        for line in f:
            # print(line)
            if j == 0:
                transcript = line.split(' ')[-1]
                f_new.write('>' + transcript)
            else:
                f_new.write(line)
            j = j + 1
        f.close()
    except FileNotFoundError:
        continue
f_new.close()

path = '/data/TA_QUIZ_RNA_regulation/result/PartIII.SHAPE-seq_analysis/3.3.motif/07.BEAM/UVR8/'
index = pd.read_csv(path+'index.csv',sep='\t')
index_up = list(index.loc[index['up_down']=='up','index'])
index_down = list(index.loc[index['up_down']=='down','index'])

f_new = open(path + 'merge_up_nouv.fa', 'w')
for i in range(len(index_up)):
    try:
        f = open(path + '/dot/' + str(index_up[i]) + '_f.dt')
        j=0
        for line in f:
            # print(line)
            if j==0:
                transcript = line.split(' ')[-1]
                f_new.write('>'+transcript)
            else:
                f_new.write(line)
            j=j+1
        f.close()
    except FileNotFoundError:
        continue
f_new.close()

f_new = open(path + 'merge_down_nouv.fa', 'w')
for i in range(len(index_down)):
    try:
        f = open(path + '/dot/' + str(index_down[i]) + '_f.dt')
        j=0
        for line in f:
            # print(line)
            if j==0:
                transcript = line.split(' ')[-1]
                f_new.write('>'+transcript)
            else:
                f_new.write(line)
            j=j+1
        f.close()
    except FileNotFoundError:
        continue
f_new.close()
 ```

**output：**
WT输出文件：`merge_up_nouv.fa`和`merge_down_uv.fa`。
分别是上调和下调的转录本，其中文件的第一行为转录本信息，第二行为核苷酸信息，第三行为结构信息。

### 3.c.4) motif discovery
利用BEAM，对motif进行富集分析。
**input：**
`merge_up_nouv.fa`、`merge_down_uv.fa`

**Software/Parameter：**
`BEAR_encoder.jar`
`java -jar BEAR_encoder.jar -i <sequence.fasta> -o <output.files> `
|files name|path|
|:------:|:------:|
|-i|Sequence(s) filename in FASTA or FASTB format|
|-o|Path to output file name.|

`BEAM_release2.5.0.jar`
`java -jar BEAM_release2.5.0.jar -f <BEAM_ready> -w 10 -W 50 -M 100`

|files name|path|
|:------:|:------:|
|-w|Min motif Width (10)|
|-W|Max motif Width (50)|
|-M|number of masks(motifs) to be computed|

**reference：**
```linux
#!/bin/sh
#SBATCH -J PartIII.SHAPE-seq_analysis
#SBATCH -p CN_BIOT
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --output=7.motif_discovery.out
#SBATCH --error=7.motif_discovery.err
#source ~/.bashrc

export PATH=/data/zhaoyizi/software/anaconda3/envs/Riboshape/bin:$PATH

cd /data/TA_QUIZ_RNA_regulation/result/PartIII.SHAPE-seq_analysis/3.3.motif/07.BEAM/WT

java -jar /data/zhaoyizi/software/BEAR_encoder/jar/BEAR_encoder.jar -i merge_down_nouv.fa -o merge_down_nouv_ready.fa
java -jar /data/zhaoyizi/software/BEAR_encoder/jar/BEAR_encoder.jar -i merge_up_nouv.fa -o merge_up_nouv_ready.fa

java -jar /data/zhaoyizi/software/beam/BEAM_release2.5.0.jar -f merge_down_nouv_ready.fa -w 10 -W 50 -M 100
java -jar /data/zhaoyizi/software/beam/BEAM_release2.5.0.jar -f merge_up_nouv_ready.fa -w 10 -W 50 -M 100

cd /data/TA_QUIZ_RNA_regulation/result/PartIII.SHAPE-seq_analysis/3.3.motif/07.BEAM/UVR8

java -jar /data/zhaoyizi/software/BEAR_encoder/jar/BEAR_encoder.jar -i merge_down_nouv.fa -o merge_down_nouv_ready.fa
java -jar /data/zhaoyizi/software/BEAR_encoder/jar/BEAR_encoder.jar -i merge_up_nouv.fa -o merge_up_nouv_ready.fa

java -jar /data/zhaoyizi/software/beam/BEAM_release2.5.0.jar -f merge_down_nouv_ready.fa -w 10 -W 50 -M 100
java -jar /data/zhaoyizi/software/beam/BEAM_release2.5.0.jar -f merge_up_nouv_ready.fa -w 10 -W 50 -M 100
```

**output：**
`~/results/merge_down_nouv_ready`
`~/results/merge_up_nouv_ready`

### 3.c.5) motif 汇总
我们根据计算得到结构motif进行过滤、汇总和绘图。

**1. 过滤结构motif**
我们以WT样本上调结构变化区域找的结构motif为例，说明过滤规则。其中 上面步骤文件夹`~/results/merge_up_nouv_ready/merge_up_nouv_ready_summary.txt`文件为结构motif汇总。
我们的过滤规则如下：
- motif宽度大于5
- pvalueMW:pvalue of MannWhitney U Test against,the null hypothesis that the two samples of partial scores of the RNA with the motif and the partials scores of top scoring background with the model PFM，要求pvalue<0.05
- AUC:derived from the U statistic, it rates the classifying power of the motif，要求AUC>0.65
- coverage:fraction of RNA in the input containing the motif,要求coverage>20%

**2. 结构motif汇总**
我们根据`merge_up_nouv_ready_summary.txt`文件可以得到过滤后motif的`pvalueMW'，'AUC`，`qBEAR`。
另外我们在`~/results/merge_down_nouv_ready/benchmark/motifs/*_ready_m1_run1.txt`中找到每个motif对应的结果改变区域的具体信息。

![10.4.9.weblogo.png](../../.gitbook/assets/10.4.9.weblogo.png)

**3.绘图**
我们利用weblogo进行结构motif绘图。我们对每个过滤后的motif进行绘图。
```linux
weblogo -a 'ZAQXSWCDEVFRBGTNHY' -f merge_up_ck_ready_m1_run1_wl.fa  -D fasta \
-o m1.pdf -F pdf --composition="none" \
-C red ZAQ 'Stem' -C blue XSW 'Loop' -C forestgreen CDE 'InternalLoop' \
-C orange VFR 'StemBranch' -C DarkOrange B 'Bulge' \
-C lime G 'BulgeBranch' -C purple T 'Branching' \
-C limegreen NHY 'InternalLoopBranch'
```

![10.4.10.weblogo3.6.0.png](../../.gitbook/assets/10.4.10.weblogo3.6.0.png)

图中为每个位点qBEAR，不同字母的大小代表这个位点为该字母的可能性（因为motif允许一定错配），字母含义为motif结构编码，Z:big Stem,A:medial Stem,Q:small Stem,X:big Loop,S:medial Loop,W:small Loop,C:big InternalLoop,D:medial InternalLoop,E:small InternalLoop,V:big StemBranch,F:medial StemBranch,R:small StemBranch,B:Bulge,G:BulgeBranch,T:Branching,N:big InternalLoopBranch,H:medial InternalLoopBranch,Y:small InternalLoopBranch




