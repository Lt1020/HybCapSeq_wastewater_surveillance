# ERCC 实验记录核查与校准统计

可在 Linux 安装的 Python 3.10+ 命令行项目，运行时不需要第三方 Python 库。

**当前功能范围：输入记录核查与人工合成 ERCC 内标的统计校准。此项目不是 FASTQ 到病原绝对定量的流程，不包含病原数据库、病原比对、病原定量或 EsViritu 的自动调用。**

版本：0.1.0。此代码作为独立子项目位于 `Lt1020/HybCapSeq_wastewater_surveillance` 的 `tools/ercc-record-audit/`，非 EsViritu 官方组件。尚未发布到 PyPI 或 Bioconda。

## 安装和演示

从项目分支获取源码后，在子项目目录运行：

```bash
git clone --branch ercc-calibration-20261016 --single-branch https://github.com/Lt1020/HybCapSeq_wastewater_surveillance.git
cd HybCapSeq_wastewater_surveillance/tools/ercc-record-audit
conda env create -f environment.yml
conda activate ercc-audit
python -m pip install --no-build-isolation --no-deps .
ercc-audit --help
ercc-audit fit-ercc \
  --samples examples/simulated/samples.tsv \
  --standards examples/simulated/standards.tsv \
  --observations examples/simulated/observations.tsv \
  --outdir results-demo
python -m unittest discover -s tests -v
```

已有 Python 3.10+ 和 setuptools 的环境可直接执行 `python -m pip install --no-build-isolation --no-deps .`。也可不安装运行 `PYTHONPATH=src python -m ercc_audit --help`。输出目录必须为新目录或空目录，防止覆盖以前的结果。

演示表全部为手工构造的数学模拟数据，`ERCC-SIM-*` 不对应真实 ERCC 产品或序列；不得作为本次实验结果、真实标准浓度或产品参数。示例故意采用不等长参考以检验长度归一化。

## 核对样本记录

将本地实验记录整理为 `samples.tsv`，即可运行：

```bash
ercc-audit audit-records --samples /path/to/samples.tsv --outdir results-records
```

RNA和DNA同时进入文库时，程序按两者质量浓度之和与实际投入体积核对总输入。空白字段表示未提供，不转换为0；未加入ERCC的对照不借用其他样本的曲线。这里仅发布模拟数据和空白模板，真实实验记录保留在本地。

## 三种输入表

运行 `ercc-audit init --outdir inputs` 创建空白表头。TSV为UTF-8、制表符分隔；使用Excel编辑时须保留列名和样本编号。单位已包含在字段名中。

### samples.tsv

| 字段 | 含义 |
| --- | --- |
| sample_id | 唯一文库编号 |
| biological_id | 原始生物样本编号 |
| replicate | 技术重复标识 |
| condition | 处理分组；加标和PBS保留为不同组 |
| sample_role | sample 或 negative_control |
| rna_ng_ul, dna_ng_ul | 同一提取液实测RNA和DNA质量浓度 |
| library_input_ul | 实际进入文库的提取液体积，不能用完整反应体积替代 |
| recorded_sample_mass_ng | 原始表记载的DNA+RNA总输入质量，不包含ERCC，可留空 |
| ercc_volume_ul | ERCC工作液加入体积；未加入时明确填0 |
| ercc_dilution_factor | 相对原液的绝对稀释因子；1∶1000填1000，不填1X |
| ercc_lot | 与standards表一致的真实批号 |
| total_effective_bases | 与ERCC比对统计口径一致的总有效碱基数；非reads数、非mapped reads数 |

记录核查允许缺批号或测序碱基数并显式标记；ERCC拟合必须补齐。`recorded_sample_mass_ng`与计算值不符时标记 `RECORDED_MASS_MISMATCH`，拟合依赖实际浓度和体积，不悄悄采用手工质量。

### standards.tsv

| 字段 | 含义 |
| --- | --- |
| lot_id | 实际ERCC产品批号 |
| ercc_id | 逐条ERCC标识，必须以ERCC-开头 |
| length_nt | 实际ERCC参考全长，整数 |
| stock_copies_per_ul | 该条ERCC在未稀释原液中的已知浓度 |

应提供该工作池包含的完整标准清单，不硬编码48或92条。同一批号不得重复；不同批号可出现相同ERCC标识。不能只填参加拟合的少数条目，否则ERCC总质量分母会不完整。输入的长度和浓度由使用者按当批资料确认，本程序不能证明其真实性。

### observations.tsv

| 字段 | 含义 |
| --- | --- |
| sample_id, ercc_id | 与样本和标准表严格匹配 |
| aligned_bases | 已完成既定上游筛选的ERCC比对碱基数；非read数或RPKMF |
| covered_bases | 按同一参考及筛选口径统计的覆盖位点数，整数 |
| quant_include | 1纳入拟合，0不纳入 |
| exclusion_reason | 不纳入时必填，保留预先设定的筛选原因 |

对每个有ERCC的样本，必须提供该批号全部标准的观测；未检出也应显式填0并排除，缺行不当作未检出。没有ERCC投入的NC不做拟合；本接口不接收该NC的背景统计。实际NC背景可在独立原始记录中保留，不在此接口里伪装成校准数据。

## 输出与公式

`sample_audit.tsv`输出RNA质量、DNA质量、两者合计、ERCC状态及质量差异标记。按已有`biological_id`和`condition`保留实验设计，不将a、b、c自动合并计算CV。

`ercc_detail.tsv`输出每条标准的已知投入量、等效提取液浓度、参考长度、比对碱基数、平均深度、覆盖广度及Y(seq)：

- 已知投入拷贝数 = 原液copies/µL × ERCC加入体积 / 绝对稀释因子。
- ERCC质量ng = 拷贝数 × 长度 × 340 / 阿伏伽德罗常数 × 10^9。
- 总输入质量 = (RNA浓度 + DNA浓度) × 提取液投入体积 + 全部ERCC质量。
- Y(seq) = (ERCC比对碱基数 / 总有效碱基数) / (该条ERCC质量 / 总输入质量)。
- 平均深度 = 比对碱基数 / 完整参考长度，包含零覆盖位置的长度。
- 归一化深度 = 平均深度 / 提取液投入体积。

340为原方案使用的RNA分子量近似，不等同于逐条产品的计量认证。软件不以RNA质量校准DNA绝对浓度。

`ercc_fits.tsv`输出ERCC投入质量分数与比对碱基分数的log10回归，以及ERCC归一化深度与已知等效copies/µL的log10回归，并报告合格ERCC的Y(seq)均值和中位数。至少4个不同浓度水平才能拟合；程序没有声称它们必然是相邻产品梯度，也不自动判定饱和或覆盖均一性。示例检查条件R²≥0.95和斜率0.8–1.2来自所附方案，仅作为探索性检查，非已验证的通用验收标准。

不会删除离群点提高R²，不生成LoD、LoQ、病原浓度或回收率。零值不会强行取对数或添加伪计数。`run_manifest.json`保留输入SHA256、软件和Python版本、运行时间及假设。

## 与EsViritu的接口事实

官方资料检索日期为2026-10-16：

- 仓库：https://github.com/cmmr/EsViritu
- 自定义库：https://esviritu.readthedocs.io/en/latest/custom-database/
- 输出解释：https://esviritu.readthedocs.io/en/latest/interpretting-outputs/

官方自定义库由FASTA、索引和元数据组成，参考标识需一致。ERCC在解释上仍是人工合成对照，不能因参与比对就认定为病毒。RPKMF不等于本程序接收的比对碱基数。

本软件不追加ERCC到病原库、不伪造NCBI TaxID、不调用EsViritu、不从其汇总丰度反推比对碱基数，也未验证EsViritu各阶段对合成标准的保留行为。需要的是已有ERCC统计表，而不是直接传入EsViritu病毒汇总表。没有上游真实统计时，只有软件演示，不是实样定量。

## 版本与发布状态

附带测试覆盖量纲、输入质量、缺失与零值、重复主键、非完整ERCC观测、批号匹配、无ERCC对照、非整数覆盖、输出覆盖保护及模拟曲线。本次未使用真实FASTQ进行验证。

子项目包含标准`pyproject.toml`、命令行入口、MIT许可证、测试及空白模板。仓库根目录的 `.github/workflows/ercc-calibration-tests.yml` 在GitHub运行Python 3.10、3.11和3.12测试；实际运行状态以Actions页面为准。

可以直接安装该分支中的子项目：

```bash
python -m pip install 'git+https://github.com/Lt1020/HybCapSeq_wastewater_surveillance.git@ercc-calibration-20261016#subdirectory=tools/ercc-record-audit'
```

