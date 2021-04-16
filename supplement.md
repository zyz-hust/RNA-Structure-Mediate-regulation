# 0. 编程工具介绍
由于本次作业需要一定的计算资源支持，我们为各个小组提供了集群账户。大作业需要使用python完成，推荐读者使用python3。建议使用pycharm进行代码的编辑、运行和调试。可以参考以下的`PyCharm 连接远程服务器`，ssh远程连接P-cluster。

## 1) PyCharm 连接远程服务器
注：只有professional版才有ssh 远程连接的功能，社区版没有。
可参考[学生如何免费试用Pycharm专业版](https://blog.csdn.net/weixin_45459911/article/details/104767525)

### 1.1) 打开PyCharm 的Preferences
![Pycharm的Preferecence.png](https://github.com/zyz-hust/RNA-Structure-Mediate-regulation/blob/30a8f2cc44ea0a77c60e923505c5349cf27a5d62/Images/Pycharm%E7%9A%84Preferecence.png)

### 1.2) 设置远程服务器&设置账号密码
![设置远程服务器.png](https://github.com/zyz-hust/RNA-Structure-Mediate-regulation/blob/30a8f2cc44ea0a77c60e923505c5349cf27a5d62/Images/%E8%AE%BE%E7%BD%AE%E8%BF%9C%E7%A8%8B%E6%9C%8D%E5%8A%A1%E5%99%A8.png)
- Build,Execition,Deployment -> Deployment。
- 点击左上角的“+”键，给server取名。
- `type` 选择 `SFTP`。
- 点击`SSH configuration` 设置`Host`、`Username`、`password`。
![SSH Configuration.png](https://github.com/zyz-hust/RNA-Structure-Mediate-regulation/blob/30a8f2cc44ea0a77c60e923505c5349cf27a5d62/Images/SSH%20Configuration.png)
- `Test Connection`检查连接是否正常，如出现问题，可能是上步SSH设置问题。
- 设置`Root path`,点击上方的`mappings`,设置本地映射目录和远程映射目录。
- 设置好后`apply`

### 1.3) 设置python解释器
![设置python解释器.png](https://github.com/zyz-hust/RNA-Structure-Mediate-regulation/blob/30a8f2cc44ea0a77c60e923505c5349cf27a5d62/Images/%E8%AE%BE%E7%BD%AEpython%E8%A7%A3%E9%87%8A%E5%99%A8.png)
- Preferences -> Project:xxxx -> Python interpreter
- 点击设置，再点击`add`,在`Python Interpreter`中设置调用的python路径
![ssh Interpreter.png](https://github.com/zyz-hust/RNA-Structure-Mediate-regulation/blob/30a8f2cc44ea0a77c60e923505c5349cf27a5d62/Images/ssh%20Interpreter.png)
- 点击`SSH Interpreter`,选择`Existing server configuration`
- 选择我们刚刚设置好的`SSH configuration`
- `next`,实现PyCharm远程服务器的连接。

### 1.4) 文件的传输以及Terminal的打开
- 选择左上角的`tools`->`Deployment`
- `upload`和`download`能够实现文件在本机和服务器上的互传。

- `Tools`->`start SSH session`
- 选择需要连接的服务器
- 随后在下方的`Terminal`处，即可进行linux 操作。
----
## 2) 使用Cluster队列提交任务
具体可参考[1.2.Cluster](https://lulab2.gitbook.io/teaching/part-i.-basic-skills/1.setup/1.2-cluster)的`4) How to use cluster`部分。这里简单介绍一下：

我们的作业需要P-cluster上完成，cluster是为多用户同时使用而设计，包含多个运算节点。为了不发生多用户同时使用相同运算节点造成拥堵和死机，需要按照queue（队列）进行使用。简而言之，本次作用中所有需要集群后台处理或大规模运算的`Shell script(.sh)`都需要在脚本包含以下部分：

```Shell
### 1. 提交任务部分，任务的说明
#!/bin/bash
#SBATCH -J tophat_test
#SBATCH -p CN_BIOT
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --output=%j.out
#SBATCH --error=%j.err

### 2. 声明该脚本中所用到的所有软件的路径，并export到PATH变量中
export PATH=/WORK/teaching/bin:$PATH

### 3. 运行代码部分
```

### 1. 提交任务部分
|Name|mean|
|:-------:|:------:|
|#SBATCH -J |给该任务命名，可自行定义|
|#SBATCH -p CN_BIOT |使用CN_BIOT这个队列|
|#SBATCH --nodes=1|使用的节点数为1|
|#SBATCH --ntasks=8|使用的任务线程数为8|
|#SBATCH --output=%j.out|运行日志输出到当前目录中，以.out结尾，以JOBID为名|
|#SBATCH --error=%j.err|运行错误日志输出到当前目录中，以 .err 结尾，以JOBID为名|

### 2.声明所用软件路径部分
`export PATH=/WORK/teaching/bin:$PATH`的`/WORK/teaching/bin`为所用到的软件所在目录的路径，如有多个软件路径，路径之间用`:`相隔，例如`export PATH=/data/zhaoyizi/software/anaconda3/envs/Riboshape/bin:/data/zhaoyizi/software/shapemapper-2.1.5:$PATH`

### 3. 运行代码部分
该部分输入实际的运行代码。
**如果是`python脚本`需要提交的集群中运算**，可以在2部分，声明python(环境)所在路径。在此部分写入`python  /path_to_your_python_script/XXX.py`，R脚本同理。

### 4. 提交、查看、取消任务
在完成脚本的编辑后，使用`sbatch`提交任务。
- 使用`sbatch`命令提交任务
```shell
sbatch test1.sh
```
提交任务后，会在当前目录下，生成`.out`和`.err`文件，`less .out/.err`文件可以查看任务运行状况和是否出错。

- `squeue`查看队列信息 

|JOBID|PARTITION|NAME|USER|ST|TIME|NODES|NODELIST|
|:------:|:------:|:------:|:------:|:------:|:------:|:------:|
|任务ID|队列名称|任务名称|用户名|脚本类型|运行时间|节点数|节点名|

- `scancel `取消任务
```shell
scancel jobid
```


==Attention ：用shell脚本提交python脚本或R脚本任务到队列时，注意写绝对路径。相对路径可能会报错。==
