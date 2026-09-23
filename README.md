# 11B_p

使用 Geant4、ROOT 和 C++17 模拟质子入射硼靶的核反应、三 α 衰变和探测器响应。当前工程位于 `g4/`，由原 `xc_260921/` 整理而来；历史日期版本通过 Git 提交和标签保留。

## 目录

| 路径 | 内容 |
| --- | --- |
| `g4/HB.cc`、`g4/CMakeLists.txt` | 程序入口、线程与随机种子、构建配置 |
| `g4/include/`、`g4/src/` | 反应模型、衰变运动学、探测器与 ROOT 数据结构 |
| `g4/macros/`、`g4/batch.sh` | 运行配置、批量运行、线程 ROOT 文件合并 |
| `g4/analysis/` | α 能谱、通道、Dalitz 图、探测器响应、效率与拟合分析 |
| `g4/cs_model/` | 截面模型计算与绘图脚本 |
| [HISTORY.md](HISTORY.md) | 来源目录、日期标签、快照提交、版本差异与回溯方法 |

`jupyter/`、`papers/`、`presentation/`、所有 `.ipynb`、ROOT 数据、构建目录、生成图像、日志、导出表格及压缩包只在本地保存，排除规则见 [.gitignore](.gitignore)。克隆仓库得到源码、宏、分析脚本和说明文件，运行后会在本地产生数据和图像。

## 当前模型与输出

- `H11BCrossSection` 和 `H11BReaction` 实现 162/675 keV 共振分量、顺序三 α 衰变、直接三 α 分支及 γ 捕获；各通道是否启用由宏决定。
- 675 keV 分量采用含能量相关库仑穿透率的 Breit-Wigner 截面；三 α 衰变支持 `symmetrizedCoherentL1L3`，由 `/h11b/` 命令配置。
- 靶、束斑、Si 阵列和虚拟球响应可在宏中配置，另包含 LaBr3、HPGe 响应。几何相关设置应放在 `/run/initialize` 之前。
- ROOT 保存反应、运行信息和按配置启用的探测器响应；字段以 `g4/include/DataStructure.hh`、`g4/src/RootIO.cc` 为准。分析旧数据时核对对应版本及 ROOT 内的运行配置。
- `/run/beamOn` 表示入射事件数，实际反应数由 `reaction` 树统计。偏置采样下的原始计数不直接等于绝对物理产额。

## 编译与运行

需要 CMake、C++17 编译器、Geant4、ROOT，以及合并线程输出所用的 `hadd`。整理时本机环境为 Geant4 11.3.2、ROOT 6.34.10；应先加载安装环境，让 CMake 能找到二者。

在仓库根目录执行：

```bash
cmake -S g4 -B g4/build -DCMAKE_BUILD_TYPE=Release
cmake --build g4/build -j8
cd g4
./batch.sh macros/run.mac 8
```

脚本默认 8 个工作线程，允许 1–8；从仓库根目录也可运行 `./g4/batch.sh g4/macros/run.mac 8`。数据写入 `g4/data/`，线程输出成功合并后默认清理线程文件；`CLEAN_THREADS=0` 可保留它们。同一数据目录应按能量点顺序使用，避免同时运行多个 `batch.sh`。

旧构建缓存含原目录的绝对路径，迁移时保存在本地归档；应在 `g4/build/` 重新配置和编译。最早的历史标签 `g4_251202` 仍使用 `B2a/` 和 `exampleB2a`，不能直接套用当前 HB 命令。

## α 能谱工作流

以下命令从 `g4/` 执行。Python 分析使用 NumPy、Matplotlib、uproot；部分探测器分析还使用 Awkward Array。ROOT 数据不随仓库分发。

```bash
# 162 keV 单点：仅启用 162 分量，100 万个入射事件
./batch.sh macros/spectrum_162_lab_05_95_00.mac 8

# 150–300 keV、步长 10 keV，每点 100 万个入射事件
python3 analysis/alpha_spectrum/run_scan.py

# 50–200 keV，各点单独绘图与总览
python3 analysis/alpha_spectrum/run_scan.py --energies 50 200 10 --separate
python3 analysis/alpha_spectrum/plot_lab_total.py --energies 50 200 10 --separate --overview
```

`run_scan.py` 校验、复用已有单点数据，并按约定写入 `ScanInfo`。独立运行 `batch.sh` 的输出带时间戳；162 单点绘图默认读取 `data/spectrum_162_lab_05_95_00.root`，也可显式传入合并 ROOT 文件：

```bash
python3 analysis/alpha_spectrum/plot_162_lab_total.py data/your_merged_output.root --total-only
python3 analysis/alpha_spectrum/plot_lab_total.py --total-only
```

每次反应的三颗 α 刚产生时的 LAB 动能分别填入同一张直方图，全立体角、每颗计一次。默认 **25 keV/bin**，横轴 MeV，纵轴 Counts；162 keV 单点范围 0–6.5 MeV，扫描谱范围 0–7 MeV。多能量点合并采用计数直接相加，不额外乘截面或束流权重。

能谱宏使用 `0.12 um` 的 `Enriched_B11` 靶、无背衬、点束流和 `1e11` 截面偏置。162 分量的基态／2⁺／直接分支设为 5%／95%／0%；扫描宏同时启用 675 分量，该分量取 0%／100%／0%。这些是当前样本配置，具体以宏和 ROOT 保存的信息为准。

## 查看与回溯版本

```bash
git log --oneline --decorate --all
git tag --list 'g4_*'
git show g4_260921
git diff g4_260629 g4_260714 -- g4/src g4/include g4/macros

# 在独立工作区查看旧版本，保留当前目录的数据
git worktree add --detach ../11B_p_260629 g4_260629
```

11 个日期标签保存对应目录的源码快照。标签日期来自目录名称，部分目录在该日期之后继续修改过；提交说明记录原 Git HEAD 和是否存在本地修改。所有迁移提交使用实际整理时间。无日期的 `xc/` 单独保存在 `g4_xc_legacy` 标签，不假定其开发日期。

原目录、原 Git 元数据和历史构建产物保存在本机 `.local_archive/unify_20260923/`，不上传。当前本地 ROOT、图像、TXT 和 `Note.ipynb` 随最新工程保留在 `g4/`。标签用于还原公开源码；完整运行现场还需相应的本地归档和 ROOT 数据。

## 验证与已知限制

整理时逐文件核对了所有快照的内容与可执行权限，并检查提交历史中的排除规则。迁移适配只调整 `batch.sh` 的工程目录识别，使它在外层统一 Git 仓库中定位 `g4/`；历史标签保留原脚本。其余物理源码按来源保存。

2026-09-23 在新路径完成 CMake 配置和 `-j8` 编译；Bash 与当前全部 Python 文件语法检查通过，批处理目录定位在仓库根目录、`g4/`、`g4/build/` 三种启动位置通过。以 162 keV 能谱宏临时缩减为 20 个入射事件、2 线程运行，确认日志完成 20 个事件、ROOT 含 20 次反应和 60 颗 α，事件号覆盖 0–19，能量守恒残差的展宽为 `3.92e-12 MeV`。测试宏和数据写在 `/tmp/`，未重跑正式百万事件样本。

历史版本未逐一重新编译或重跑大规模模拟。当前 `HB.cc` 尚未将 Geant4 宏执行错误转换为非零进程退出码，`batch.sh` 在无新 ROOT 文件时也可能仅警告并返回成功；生产计算应检查 Geant4 日志与实际完成事件数。这是后续应优先修复的运行可靠性问题。

后续修改建议按物理模型、几何、ROOT 数据结构、分析分开提交；重要计算使用有说明的标签，在本地运行记录中保留提交号、宏、随机种子、入射事件数及数据文件名。

## 许可证

保留原仓库的 [MIT License](LICENSE)。源文件中已有的 Geant4 等组件版权与许可声明一并保留。
