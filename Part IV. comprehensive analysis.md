# 1) 转录本丰度变化和翻译效率(TE)变化程度之间关系
根据Part I 中(DESeq)计算得到的`wt_rawdata.csv`差异表达矩阵与Part II 中计算得到的差异翻译矩阵`wt.uvb-vs-nouvb.TE_new.csv`。我们来绘制转录本丰度变化和翻译效率(TE)变化关系图。

```python
import numpy as np
import pandas as pd
import re
import seaborn as sns
import matplotlib.pyplot as plt

light_data=pd.read_csv('/data/TA_QUIZ_RNA_regulation/data/gene_list/light_gene/light_gene_list.txt',sep='\t')
light_list=list(light_data['gene'])

RF_data=pd.read_csv('/data/TA_QUIZ_RNA_regulation/result/PartII.Ribo-seq_analysis/7.TE/wt.uvb-vs-nouvb.TE_new.csv',sep='\t')
RF_data=RF_data.reset_index()
RF_data=RF_data.rename(columns={'index':'gene'})

RS_data=pd.read_csv('/data/TA_QUIZ_RNA_regulation/result/PartI.RNA-seq_analysis/differential_expression/5.DESeq2/wt/wt_rawdata.csv',sep=',')
RS_data=RS_data.rename(columns={'Row.names':'gene'})

data=pd.merge(RS_data[['gene','pvalue','padj','log2FoldChange']],RF_data[['gene','pvalue_final','pvalue.adjust','log2FC_TE_final']],on='gene',how='outer')
data.columns=['gene','pvalue(RNA-seq)','padj(RNA-seq)','log2FoldChange(RNA-seq)','pvalue(TE)','padj(TE)','log2FoldChange(TE)']
data['padj(RNA-seq)']=data['padj(RNA-seq)'].fillna(1)
data['pvalue(TE)']=data['pvalue(TE)'].fillna(1)

data['group'] = 'darkgray'
result_1=data.loc[(data['pvalue(TE)']<0.05)&(data['log2FoldChange(TE)']>0)&(data['padj(RNA-seq)']>0.05),:]
result_2=data.loc[(data['pvalue(TE)']<0.05)&(data['log2FoldChange(TE)']<0)&(data['padj(RNA-seq)']>0.05),:]
result_3=data.loc[data.gene.isin(light_list),:]

#确定坐标轴显示范围
xmin=-3
xmax=6
ymin=-8
ymax=8

#绘制散点图
fig = plt.figure(figsize=plt.figaspect(5/6)) #确定fig比例（h/w）
ax = fig.add_subplot()
ax.set(xlim=(xmin, xmax), ylim=(ymin, ymax), title='')
ax.scatter(data['log2FoldChange(RNA-seq)'], data['log2FoldChange(TE)'], s=15, c=data['group'])

ax.scatter(result_1['log2FoldChange(RNA-seq)'], result_1['log2FoldChange(TE)'], s=20, marker='.',c='#cc0000',label = '580 TE up mRNAs')
ax.scatter(result_2['log2FoldChange(RNA-seq)'], result_2['log2FoldChange(TE)'], s=20, marker='.',c='steelblue',label = '1109 TE down mRNAs')
# ax.scatter(result_3['log2FoldChange(RNA-seq)'], result_3['log2FoldChange(TE)'], s=20, marker='^',c='black',label = '126 light-related mRNAs')

ax.spines['right'].set_visible(False) #去掉右边框
ax.spines['top'].set_visible(False) #去掉上边框
ax.spines['bottom'].set_linewidth(2)###设置底部坐标轴的粗细
ax.spines['left'].set_linewidth(2)####设置左边坐标轴的粗细

plt.tick_params(labelsize=15)
ax.set_xticks(range(xmin,xmax,1)) #设置x轴刻度起点和步长
ax.set_yticks(range(ymin,ymax,2)) #设置y轴刻度起点和步长
# font2 = {'family': 'Times New Roman','weight': 'normal','size': 15}
font2 = {'weight': 'normal','size': 18}
plt.xlabel('log2FoldChange(RNA-seq)',font2)
plt.ylabel('log2FoldChange(TE)',font2)
plt.legend(fontsize=10)
font3 = {'weight': 'normal','size': 15}
plt.text(3,4, 'r = -0.270',font3)
plt.tight_layout()
# plt.show()
# 图片输出路径
plt.savefig('/data/TA_QUIZ_RNA_regulation/result/PartIV.comprehensive_analysis/RNA-seq_TE/riboseq/logFC_TE_corr.png')
plt.close()

data_= data[["log2FoldChange(RNA-seq)","log2FoldChange(TE)"]]
data_corr = data_.corr()
```

此外，这里，我们提供光通路相关的`gene list`，路径位于`/data/TA_QUIZ_RNA_regulation/data/gene_list/light_gene/light_gene_list.txt`

**Questions：**
1. 你得到的转录本丰度变化和翻译效率(TE)变化关系图是怎样的？
2. 二者是否存在相关关系，若存在，试解释产生这种相关性的原因。
3. 我们提供了光通路相关的`gene list`,可自行探索`light-related mRNAs`在图中的趋势。(选做)


# 2) 转录本结构改变程度和翻译效率TE变化程度之间关系

....

# 3) 翻译效率变化基因的motif分析
本节我们将分别对TE上调和TE下调的所有基因的3‘UTR和5'UTR富集序列motif。

## 提取所有3’UTR/5'UTR的fastq文件
这里以`WT型未加光`为例
**input：**
`/data/TA_QUIZ_RNA_regulation/data/ATH/GFF/Ath_genes.gff`：参考基因组注释文件
`/data/TA_QUIZ_RNA_regulation/data/riboshape_liulab_batch4/final.modified_umodified/col_nouv/`：包含转录本全长序列信息文件
`/data/TA_QUIZ_RNA_regulation/data/ATH/GTF/shape_map/result/transcript_exon_location.csv`：转录本位置信息文件

```python
import numpy as np
import pandas as pd

col_uv_f = '/data/TA_QUIZ_RNA_regulation/data/riboshape_liulab_batch4/final.modified_umodified/col_nouv/'

gff_path = '/data/TA_QUIZ_RNA_regulation/data/ATH/GFF/Ath_genes.gff'
data_gff = pd.read_csv(gff_path, sep='\t')
data_location = pd.read_csv( '/data/TA_QUIZ_RNA_regulation/data/ATH/GTF/shape_map/result/transcript_exon_location.csv',sep='\t')
data_gff=pd.merge(data_gff,data_location[['transcript_id','strand']],on='transcript_id',how='left')
gene=pd.read_csv('/data/TA_QUIZ_RNA_regulation/data/gene_list/wt/ribo_wt_gene_list.txt',sep='\t',header=None)

gene_list=list(gene[0])

data_1 = pd.read_csv(col_uv_f + '/final.modified_unmodified_new', sep='\t')
data_1.columns = ['transcript_id', 'Nucleotide', 'Modified_mutations', 'Modified_effective_depth',
                     'Untreated_mutations', 'Untreated_effective_depth', 'R1']
data_1 = data_1[['transcript_id', 'Nucleotide', 'Modified_mutations', 'Modified_effective_depth',
                       'Untreated_mutations', 'Untreated_effective_depth', 'R1']]
data_1 = data_1.drop_duplicates()
data_1['gene']=data_1['transcript_id'].map(lambda x:x.split('.')[0])


data=data_1
data_nuc=pd.DataFrame(columns={'gene','transcript_id','5UTR','CDS','3UTR'})
j=0

for i in gene_list:
    data_ = data.loc[data['gene'] == i, :]

    transcript_id = list(data_['transcript_id'])[0]

    data_location_ = data_location.loc[data_location['transcript_id'] == transcript_id, :]
    data_gff_ = data_gff.loc[data_gff['transcript_id'] == transcript_id, :]
    if len(data_gff_) == 0:
        continue
    else:
        if list(data_gff_['strand'])[0] == '+':
            five_UTR_CDS = list(data_gff_.loc[data_gff_['location'] == 'CDS', 'strat'])[0]
            three_UTR_CDS = list(data_gff_.loc[data_gff_['location'] == 'CDS', 'end'])[0]
        else:
            five_UTR_CDS = list(data_gff_.loc[data_gff_['location'] == 'CDS', 'end'])[0]
            three_UTR_CDS = list(data_gff_.loc[data_gff_['location'] == 'CDS', 'strat'])[0]
        Nuc = list(data_['Nucleotide'])[0]
        location = np.array(list(data_location_['location'])[0].split(',')).astype('int')
        five_UTR_CDS_ = np.where(location == five_UTR_CDS)[0][0]
        three_UTR_CDS_ = np.where(location == three_UTR_CDS)[0][0]

        data_5UTR_ = Nuc[:five_UTR_CDS_]
        data_CDS_ = Nuc[five_UTR_CDS_:three_UTR_CDS_ + 1]
        data_3UTR_ = Nuc[three_UTR_CDS_ + 1:]
    data_nuc.loc[j, 'gene'] = i
    data_nuc.loc[j, 'transcript_id'] = transcript_id
    data_nuc.loc[j, '5UTR'] = data_5UTR_
    data_nuc.loc[j, '5UTR_len'] = len(data_5UTR_)
    data_nuc.loc[j, 'CDS'] = data_CDS_
    data_nuc.loc[j, 'CDS_len'] = len(data_CDS_)
    data_nuc.loc[j, '3UTR'] = data_3UTR_
    data_nuc.loc[j, '3UTR_len'] = len(data_3UTR_)
    j = j + 1
    print(j)

data_nuc=data_nuc[['gene','transcript_id','5UTR_len','CDS_len','3UTR_len','5UTR','CDS','3UTR']]
print('ALL',len(data_nuc))

n=15
data_nuc_=data_nuc.loc[data_nuc['5UTR_len']>n]
data_nuc_.index=range(len(data_nuc_))
print('5UTR',len(data_nuc_))

f = open('/data/TA_QUIZ_RNA_regulation/data/gene_list/motif/merge_5UTR.fasta', 'w')
for i in range(len(data_nuc_)):
    transcript = data_nuc_.loc[i, 'transcript_id']
    nucleotide= data_nuc_.loc[i, '5UTR']
    f.write('>' + transcript+'\n')
    f.write(nucleotide + '\n')

f.close()

data_nuc_=data_nuc.loc[data_nuc['3UTR_len']>n]
data_nuc_.index=range(len(data_nuc_))
print('3UTR',len(data_nuc_))
f = open('/data/TA_QUIZ_RNA_regulation/data/gene_list/motif/merge_3UTR.fasta', 'w')
for i in range(len(data_nuc_)):
    transcript = data_nuc_.loc[i, 'transcript_id']
    nucleotide= data_nuc_.loc[i, '3UTR']
    f.write('>' + transcript+'\n')
    f.write(nucleotide + '\n')

f.close()
```

**output：**
`merge_5UTR.fasta`：汇总的所有5‘UTR的序列信息
`merge_3UTR.fasta`：汇总的所有3’UTR的序列信息

## 提取TE变化区域的3‘UTR、5’UTR序列信息
**方法一：**
根据Part II 得到的TE上/下调的`gene list`，替换上一步中我们所使用的未经筛选的所有WT的`gene list`—`ribo_wt_gene_list.txt`。最终得到的即为TE变化基因的5'UTR，3‘UTR序列信息。

**方法二：**
根据Part II 得到的TE上/下调的`gene list`，来筛选上一步中得到的所有5’UTR序列信息和3‘UTR序列信息。——`merge_5UTR.fasta`和`merge_3UTR.fasta`

**output：**
`merge_TE_up_5UTR.fasta`/`merge_TE_down_5UTR.fasta`:汇总的TE上/下调的所有基因5’UTR的fasta信息
`merge_TE_up_3UTR.fasta`/`merge_TE_down_3UTR.fasta`:汇总的TE上/下调的所有基因3’UTR的fasta信息

### MEME在线motif分析
我们利用MEME在线版进行motif分析，地址为[http://meme-suite.org/tools/meme](http://meme-suite.org/tools/meme)。
使用MEME的Motif discovery工具，预测输入序列上的motif信息

![10.4.7.MEME.png](../../.gitbook/assets/10.4.7.MEME.png)

我们上传`merge_TE_up_5UTR.fasta`/`merge_TE_down_5UTR.fasta`后提交运算。

![10.4.8.discovery_motif.png](../../.gitbook/assets/10.4.8.discovery_motif.png)

在得到的结果中我们根据Evalue进行筛选，Sites为找到的motif在几个结构改变区域中出现过，由此可以计算出该motif的覆盖率。我们选择覆盖率较高的motif，这里取大于4%。


**question：**
1. `WT/UVR8`型TE上/下调的基因的3‘UTR/5'UTR富集到的motif分别是什么？
2. motif区域是否发生结构变化
