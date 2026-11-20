# SynCom-MetaFlow

SynCom-MetaFlow 接收双端 FASTQ、样本 TSV、SynComDB 目录和模块列表，生成可直接执行的 Bash 脚本。支持 M01–M06 独立运行、自动加入依赖、分批处理及合并。

## Python 主控与分析目录

Python 入口 `syncom-metaflow` 调用 `syncom/cli.py` 解析样本、校验参考和工具、选择依赖并生成样本级任务。执行模式共用同一份计划。运行 `./syncom-metaflow modules` 可列出六个模块。模块的工具命令定义在 `syncom/modules/`，包括功能定量 `m04.py`、SNV `m05.py` 和基因 PAV/CNV `m06.py`。

`--outdir` 接受用户指定的分析目录。每个目录下的 `shell/` 指向该次计划快照，含可直接执行的 `run_all.sh`、每样本任务、命令定义、配置和输入清单。单样本片段由总控脚本加载，共用锁、缓存和日志。实际分析结果与 `shell/` 同级。历史命令中的 `scripts/` 为兼容入口。

## 代码与数据范围

本仓库仅提供核心程序、数据库资产构建脚本、示例样本表和工具配置。数据库、参考序列、测序数据、分析结果、开发测试脚本以及开发提示词不随代码发布。下文的 `db/SynCom122/` 指开发环境中的本地数据库，使用前须自行准备兼容的数据库目录和软件。Python 主程序仅使用标准库；运行环境为 Linux、Python 3.10 或以上版本及 Bash。

## 运行

在项目目录中运行以下命令。将 `/path/to/samples.tsv` 替换为本批实验样本表。

```bash
./syncom-metaflow doctor
./syncom-metaflow run --samples /path/to/samples.tsv --db db/SynCom122 \
  --modules M01,M02,M03,M04,M05,M06 --threads 4 --outdir results/batch01
```

也可生成脚本后后台启动。

```bash
./syncom-metaflow plan --samples /path/to/samples.tsv --db db/SynCom122 \
  --modules M01,M02,M03,M04,M05,M06 --threads 4 --outdir results/batch01
nohup bash results/batch01/shell/run_all.sh > results/batch01/run.log 2>&1 &
```

`run` 与 `plan` 生成相同的脚本结构。输出目录里的 `shell/run_all.sh` 可以从任意工作目录调用。它读取同目录的样本路径、参数和模块脚本，在 Linux 上直接调用配置好的分析程序。HPC 作业可运行同一个入口，当前按样本串行执行。比对阶段同时运行 BWA 和排序，两个进程至少各占一个核。

原始 reads 默认经过 M01。只做质控可指定 `--modules M01`。只有一对已清洗 FASTQ 的样本可用 `--clean-reads --modules M02,M03` 直接分析。输入格式见 `examples/samples.tsv`。TSV 必须有 `sample_id`、`r1` 和 `r2`，路径相对于样本表所在目录解析。多 lane 样本须为每行填写唯一的 `lane_id`，按样本表行序合并后质控。分组、对照和接种组合等额外列保存在运行快照中。

## 软件与数据库

公开仓库中的 `config/tools.tsv` 使用 PATH 中的通用命令名。已在原开发环境验证的软件包括 fastp 0.23.1、BWA 0.7.17、samtools 1.20、featureCounts 1.6.3、bcftools 1.19 和已有的 mOTUs 环境。换机器时可用 `--tools /path/to/tools.tsv` 指定另一张两列 TSV，左列是工具名，右列是可执行路径或 PATH 中的命令。软件是否能被定位可由 `doctor` 检查，实际版本写入每次运行的 `provenance/versions.tsv`。

共享 mOTUs 环境用个性化库运行时报告 `3.0.1_NEWDB`，这是该运行输出的版本标记。库中包含 `ref_mOTU_v31_*` 及 `NEWDB_*` 条目。已经用现有库完成测试。更换 mOTUs 软件或数据库前应另做兼容性验证。

`db/SynCom122/manifest.tsv` 引用已有的 BWA FASTA、索引、contig 映射与 mOTUs 数据库。清单中的相对路径以清单目录为基准。菌株表 `strains.tsv` 来自用户提供的工作簿，工作表和源文件校验和记在 `strains.source.tsv`。122 个元数据 ID 与参考映射一一对应。运行前检查 contig 名称、长度及菌株元数据的连接。

现有 GFF、CDS 和蛋白 FASTA 的路径见 `db/SynCom122/sequence_sources.tsv`。从 GFF 取 CDS 坐标后，按链方向从联合参考取序列，并与原 FFN 逐条比较。122 个基因组的 569,717 条 CDS 均已通过。派生的 SAF、BED、基因映射和 GO、COG/KOG、KEGG KO、CAZy、Pfam 注释表位于 `db/SynCom122/gene_assets/`，包含 7,361,788 条基因–功能映射。每株的导入计数在 `import_report.tsv`，被排除的记录写入 `rejected_cds.tsv`。如需从新的源文件重建，可在新目录运行 `python3 scripts/build_gene_assets.py --help` 查看参数。

## 结果与续跑

矩阵以实体为行、样本为列。主要汇总文件都在 `results/<batch>/summary/`。

| 模块 | 汇总结果 | 单样本产物 |
|---|---|---|
| M01 | `qc.tsv`，含输入和清洗后的配对数 | `qc/<sample>/` 下的 clean FASTQ、fastp JSON 与 HTML |
| M02 | `taxonomy.tsv`、`taxonomy_matrix.tsv` | `taxonomy/<sample>/` 下的原始 mOTUs 表和标准长表 |
| M03 | `strains.tsv`、`strains_matrix.tsv`、菌株组长表及矩阵 | `strains/<sample>/` 下的 contig 覆盖及基因组统计 |
| M04 | `genes.tsv`、`functions.tsv`、`gene_fpkm_matrix.tsv`、`function_fpkm_matrix.tsv` | `functions/<sample>/` 下的原始 featureCounts 计数与汇总 |
| M05 | `snv.tsv`、`snv_eligibility.tsv` | `snv/<sample>/` 下的可调用区间 BED、SNV 和参考筛选状态 |
| M06 | `pav_cnv.tsv`、`gene_pav_matrix.tsv`、`gene_cnv_matrix.tsv`、`pav_cnv_eligibility.tsv` | `pav_cnv/<sample>/` 下的逐基因覆盖、状态和相对拷贝数 |

M02 使用 mOTUs 原始相对丰度。个性化库若重复输出相同的 `unassigned`，标准表只取一次，无标记信号时原工具可能输出全零。

M03、M04 共用 `work/<sample>/alignment/` 下的联合参考排序 BAM 和 BAI。M03 的 `reads_mapq20` 是 MAPQ ≥ 20 的读段数。平均深度按 contig 长度加权。`target_depth_fraction` 是单个参考的平均深度占所有目标参考平均深度之和的比例，全库无覆盖时为 `NA`。组级矩阵按菌株表中的 ANI ≥99.99% 分组汇总高 MAPQ 证据。多重比对 reads 不参与分配，组级数值只反映可用的高置信证据。

M04 在联合参考上用符合比对条件的双端片段计数。FPKM 分母为样本全部 clean read pairs，长度为参考 CDS 长度。功能表按注释源和条目把对应基因的 FPKM 相加。一条基因可进入多个功能条目，因此不同条目之和可能超过基因丰度总和。

M05、M06 自动加入 M03。默认要求参考基因组覆盖广度 ≥0.5、平均深度 ≥5，并排除元数据中属于多成员 ANI ≥99.99% 组的参考。可通过 `--min-genome-depth` 和 `--min-genome-breadth` 调整参考筛选门槛。自建库无 ANI 列时只根据覆盖筛选。

- M05 用 bcftools pileup 统计高质量碱基，不进行二倍体基因型推断。位点深度至少 5，ALT 至少 3 个碱基、正反链各至少 1 个，频率至少 5%。`snv.tsv` 使用 1 基坐标，`callable.bed` 使用 0 基半开区间。空变异表应与可调用区间和筛选状态一起查看。
- M06 使用基因与参考基因组的平均深度比估计相对拷贝数。基因广度 ≥80% 且深度 ≥1 判为 `present`。参考基因组平均深度 ≥10 且基因完全无覆盖时判为 `absent`。其余为 `uncertain`，拷贝数为 `NA`。这是参考基因层面的估计，覆盖不均、参考差异和背景近缘菌会影响判定。

使用原命令或重新执行 `shell/run_all.sh` 可续跑。脚本核对输入标识、数据库、参数、软件和产物后重跑受影响任务。每项任务有独立日志。以 `provenance/status.tsv` 为准读取完成批次，`provenance/summary_catalog.tsv` 列出该次运行发布的结果。

## 合并批次

```bash
./syncom-metaflow merge --outdir results/combined results/batch01 results/batch02
```

合并校验运行状态、数据库文件标识、软件版本、模块和参数。样本 ID 必须跨批次唯一，输出目录须为新目录。合并从批次长表重新生成矩阵，并记录样本来源。单样本 BAM、日志和可调用区间保留在原批次。

## 验证范围

开发环境中的本地测试已核对多 lane、带空格的路径、前台与后台一致性、BAI 缺失后的续跑、已知 SNV、功能计数、基因存在与缺失、低覆盖跳过和跨批合并。测试脚本不随公开代码发布。现有 SynCom122 全库和个性化 mOTUs 库也用模拟 reads 做过六模块联测。正式样本需要环境对照验证污染背景、菌株可区分性和检测限。
