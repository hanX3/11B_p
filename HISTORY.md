# Geant4 版本沿革与回溯

保留远端初始提交 `c4489869d4c1481e429f08bdee801402ef6ce1c0`，在其后追加公开范围配置及源码快照，不改写该历史。

日期版本依目录名按顺序导入 main，每版均使用 g4/ 路径并创建带说明的标签。快照取自迁移前的实际文件，包含尚未提交的源码修改。逐一核对文件内容与可执行权限，只剔除笔记、原始数据、构建产物、生成结果和本地工具配置。

目录日期仅作为版本标签，不能据此断言文件修改日期或模拟发生日期；提交时间为实际整理时间。“原 HEAD”指原目录自带仓库的提交，不能替代包含本地修改的快照。原 Git 仓库及详细开发历史在本地归档中保留。

## 日期版本

| 原目录 | 标签 | 快照提交 | 公开文件数 | 原 HEAD | 含未提交修改 |
| --- | --- | --- | ---: | --- | --- |
| `xb_251202` | `g4_251202` | `2e8c3358f255` | 49 | `无独立仓库` | 不适用 |
| `xb_251231` | `g4_251231` | `a2fdf09bf7dc` | 38 | `无独立仓库` | 不适用 |
| `xb_260614` | `g4_260614` | `513d0f8b0869` | 64 | `无独立仓库` | 不适用 |
| `xb_260616` | `g4_260616` | `0a3f7b57ad6e` | 65 | `无独立仓库` | 不适用 |
| `xb_260624` | `g4_260624` | `651efc895b9f` | 142 | `8d7e9a76ac5f` | 是 |
| `xb_260629` | `g4_260629` | `f9bc5c2f7335` | 161 | `8e3b6e8307eb` | 否 |
| `xb_260703` | `g4_260703` | `94fbc95ee9c9` | 133 | `8e3b6e8307eb` | 是 |
| `xb_260707` | `g4_260707` | `edcdaee04b54` | 107 | `8e3b6e8307eb` | 是 |
| `xb_260710` | `g4_260710` | `593e378f8a11` | 118 | `0f8a232c6b2b` | 否 |
| `xb_260714` | `g4_260714` | `4f38be0ebbc1` | 131 | `7b9c35071bda` | 否 |
| `xc_260921` | `g4_260921` | `b9cfe5378983` | 136 | `7b9c35071bda` | 是 |

## 各版本说明

### 2025-12-02 · `g4_251202`

快照：`2e8c3358f2550c67b6c0bca47bf508c261c5a646`。来源：`xb_251202/`。

- 保留 Geant4 B2a 示例结构、exampleB2a 入口和 JUNAReaction/JUNACrossSection 原型。
- 收录条形探测器 StripSD/StripHit、自定义物理过程、运行宏与 tryTGen.C。
- 顶层 CMake 引用 B2a；此历史布局按原样保存，未把它改写为现代 HB 布局。

源码清单 SHA-256：`bda6616359677d419e6c38a0fc98e8be5df9609ffd6b38e585321ff18c1c616b`。

### 2025-12-31 · `g4_251231`

快照：`a2fdf09bf7dc376f89462cd0d59a465ac61bef9a`。来源：`xb_251231/`。

- 从 B2a 目录结构转换为顶层 HB.cc、include/、src/ 和 macro/。
- 保留 p11BReaction/p11BCrossSection、ProtonLimiter 与 StripSD/StripHit 版本。
- 使用 C++17、Geant4 和 ROOT 构建，运行宏目录仍为单数 macro/。

源码清单 SHA-256：`c7bb5ac1062057e54b10dc07c35e264c797d925b33d9064af46190142b2433af`。

### 2026-06-14 · `g4_260614`

快照：`513d0f8b086969e53c2f7a3814183520a121f5c4`。来源：`xb_260614/`。

- 保存 H11BReaction/H11BCrossSection/H11BAngularDistribution 及库仑穿透率表。
- 采用 Si、LaBr3、HPGe 阵列与独立 RootIO/DataStructure，运行宏统一为 macros/。
- 收录 Dalitz 分析与库仑曲线检查脚本；保留历史 README 的原文供追溯。

源码清单 SHA-256：`5bc522bb98818aa0de9ba6ceb08a102265c5acf5d04266062b3adc3a7ac107fc`。

### 2026-06-16 · `g4_260616`

快照：`0a3f7b57ad6e7862a98291f6eda7a08a70c488ca`。来源：`xb_260616/`。

- 新增 batch.sh 并纳入 CMake 运行文件，保留线程输出合并流程。
- H11BReaction 在三 alpha 总四动量定义的质心系计算各 alpha 动能，并写入 ROOT。
- 同步保存 DataStructure、RootIO、Dalitz 分析及探测器常数调整。

源码清单 SHA-256：`d8d86e597461b2b12667944946ac47729c0f8b59c946b6896cc4527c0d959ab9`。

### 2026-06-24 · `g4_260624`

快照：`651efc895b9f3b6b662a893635766f557e30e404`。来源：`xb_260624/`。

- 保存 H11BConfig、OutputConfig/OutputPath 和运行时宏控制，以及反应和探测器诊断字段。
- 保留 165/675 keV 分量、gamma 捕获、历史背景/直接三 alpha 配置及角分布。
- 收录截面比较、三 alpha 运动学、角分布、gamma 和 Si 响应验证脚本。
- 原工作树中 DetectorConstruction 的两处未提交文件修改包含在快照内。

源码清单 SHA-256：`b450de914263678908e6f5bb65e905116612e8f782ccd7f07b59295b4ea23394`。

### 2026-06-29 · `g4_260629`

快照：`f9bc5c2f7335776ee956948844264403297e3416`。来源：`xb_260629/`。

- 保存 675 keV 对称化相干 L=1/L=3 三 alpha 生成器及严格采样诊断。
- 此版本的 675 分量使用 weighted_chebE16 拟合，直接过程按各共振的顺序衰变比例拆分。
- 验证脚本集中到 validation/ 的通道、角分布、运动学和探测器响应子目录。
- 保留 CMake 宏目录刷新逻辑与 HPGe、ROOT 诊断调整。

源码清单 SHA-256：`2f398c87ee3d48ea639a0b8acb8beec56f0ec45d61f9e1811fd5e9fecbda400b`。

### 2026-07-03 · `g4_260703`

快照：`94fbc95ee9c927959c2a283fc7fe1a9c59fc8731`。来源：`xb_260703/`。

- 新增 SiArrayConfig，保存运行时 Si 几何、条带读出及束流相关设置。
- 收录 three_alpha_coincidence、alpha0 角分布和主 alpha 能谱分析。
- 原工作树删除和新增的宏、验证脚本按实际状态纳入；不以旧 HEAD 代替当前源码。

源码清单 SHA-256：`23204d7598a6e28e17ac9e793c9a1841884ea0ad63b45715347a85ce7a2a65d9`。

### 2026-07-07 · `g4_260707`

快照：`edcdaee04b543ac013c7552f7d57116c972d461c`。来源：`xb_260707/`。

- 保存从 165 标记到 162 keV 模型与宏命令的更新，以及对应截面和 ROOT 字段调整。
- 收录 162 keV 主/次 alpha、直接过程、gamma 与 Dalitz 分析。
- 保存 675 严格模型检查脚本和该工作树对旧验证目录的精简。

源码清单 SHA-256：`27fc4ed7a89aca8755a64d68a8edde9b461052bf43ab9d0f227ede762a5fa350`。

### 2026-07-10 · `g4_260710`

快照：`593e378f8a1152bb1f22fa53bcf9d6462476e992`。来源：`xb_260710/`。

- 保存鼓形 Si 阵列配置、独立双面条带读出及几何初始化顺序。
- 保留重新组织的探测器事件 ROOT 输出，以及三 alpha 符合分析的配套修改。
- 收录 675 Wang 截面分析、675 通道验证和靶厚扫描宏。

源码清单 SHA-256：`877aa18c3623cbbbc595fdcca3fa4bc37e5b2babd4fc072e458238e83fb6022d`。

### 2026-07-14 · `g4_260714`

快照：`4f38be0ebbc10fb9519cbb9b513b640f05a0e6b6`。来源：`xb_260714/`。

- 新增 BeamConfig、VirtualSphereConfig/Hit/SD 与粒子标签/轨迹信息。
- 保存 Si 几何、靶、束斑、探测器事件和虚拟球输出的运行时配置。
- 675 截面转为带能量相关库仑穿透率的 Breit-Wigner；分析集中到通道、Dalitz、效率、响应和贝叶斯拟合目录。
- Note.ipynb 和所有图像、拟合输出保留本地，公开快照只含程序和宏。

源码清单 SHA-256：`7f33e8191605e92b9c6baa330b1acb1889643f6869cb583f92178688baf0eb29`。

### 2026-09-21 · `g4_260921`

快照：`b9cfe53789835d84ffc2a2f36563a08f469c9752`。来源：`xc_260921/`。

- 在 2026-07-14 工作树基础上新增 162 keV 单点与 162/675 分量共同参与的 LAB 能谱扫描宏。
- 新增 alpha_spectrum 绘图和扫描脚本：支持合并、逐能量点及总览；每颗 alpha 独立填谱，默认 25 keV/bin。
- 保存 150–300 keV 扫描与 50–200 keV 独立谱入口；现有 ROOT、PNG/PDF、TXT 和 notebook 不纳入 Git。
- HB.cc 和 batch.sh 默认 8 线程、允许 1–8；快照包含这些已有本地修改。

源码清单 SHA-256：`88cdef64e6c0d62d4d4d6c8486148ec3a3301732901d861588b07440ad7ac8fb`。

## 无日期 xc 版本

标签：`g4_xc_legacy`；快照：`e3d935cd47469880484db5f8cea5f5edd62757f6`；原 HEAD：`3b977234ed6a36861917561b5b75dd20f3d2ff21`。

- 保存无日期 xc 工作树：HB 入口、H11B 反应、Si/LaBr3/HPGe 阵列、RootIO 和角分布分析。
- 原 Git 来自 HFRS_gamma，当前工作树包含大量未提交修改；不将祖先提交时间当作此快照的开发日期。
- 该快照从迁移准备提交单独分出，仅由 g4_xc_legacy 标签保留，不插入日期版本的演进顺序。

该快照从迁移准备提交独立分出，未合并到日期快照主线。`git log --all --graph --oneline` 可查看它。

## 导入后的适配

- 当前 main 使用 g4_260921 的模型与分析程序；迁移适配提交单独修正 batch.sh 将 Git 顶层误当作 Geant4 工程目录的问题，并补充仓库说明。
- 历史标签保留原文件。需要在新布局运行历史 batch.sh 时，应核对其工程目录识别逻辑。
- 新布局使用 g4/build/ 和 g4/data/。原最新版本的构建目录与 Git 元数据分别保存在 `.local_archive/unify_20260923/xc_260921.build/`、`.local_archive/unify_20260923/xc_260921.git/`。
- 其他原目录整体移动到 `.local_archive/unify_20260923/<原目录名>/`，目录内旧数据、笔记、构建结果和 .git 保留；本地归档不上传。
- 最新版原 batch.sh 另存为归档中的 `xc_260921.batch.sh`，本地协作说明的原文另存为 `xc_260921.AGENTS.md`。

## 使用示例

```bash
git show g4_260629
git diff g4_260707 g4_260710 -- g4/include/DataStructure.hh g4/src/RootIO.cc
git worktree add --detach ../11B_p_260714 g4_260714
git archive --format=tar --output=/tmp/g4_260921.tar g4_260921 g4
```

使用旧版本时配套读取该版本的宏和 ROOT 字段定义。历史快照只保证源码回溯，不宣称已在当前环境通过编译或物理验证。结果再现还需要相应数据、随机种子和完整运行配置。
