"""One-off, repository-scoped issue maintenance. No product or main writes.

Run only from the dedicated maintenance workflow. The plan is based on the
read-only audit in Actions run 38104084787. Every issue is re-read before writes;
original bodies are preserved verbatim in a historical details section.
"""
import hashlib
import json
import os
import time
import urllib.error
import urllib.request

REPO = 'Aequiludium/cocycle-rs'
BASE = 'https://api.github.com/repos/' + REPO
BRANCH = 'refs/heads/maintenance/p4-dependencies-20261011'
MARKER = '<!-- p4-scope:20261011:interfaces-only:v2 -->'
HISTORY = '<details>\n<summary>本次范围修订前的完整记录（历史；旧算法实施要求及旧依赖不再生效）</summary>\n\n'
EXPECTED_MAIN = 'a7ee9e911dcdeaea137b682bf383f7c2e1d090d5'
OLD = {i: [] for i in range(57, 98)}
OLD.update({60:[59],61:[59],62:[60,61],63:[62],65:[59],66:[63,65],68:[60],69:[63,66,68],71:[60,61],72:[63,71],73:[63,66,71],74:[73],75:[63,66,71],76:[71],78:[61,71],79:[63,78,73],81:[61,65,71],82:[63,66,73,81],84:[61,71],85:[63,79,82,84],87:[69,72,79,82,85],88:[66,73,74,79,85],89:[69,73,79,85],90:[69,75,79,85],91:[69,73,79,82,85],92:[69,73,85],93:[69,71,84],94:[87,88,89,90,91,92],96:[94]})
NEW = {n: list(ds) for n, ds in OLD.items()}
NEW.update({69:[63,68],71:[63],72:[71,87],73:[60,61],74:[60,61],75:[60,61],76:[],78:[61],79:[71,78,87],81:[61],82:[71,81,87],84:[61],85:[71,84,87],87:[63],88:[73,74],89:[69,73],90:[69,75],91:[81],92:[84],93:[],94:[72,79,82,85],96:[]})
HASHES = {
57:'b3ce6396585fbe43aa2d77f84484e21b73df8f2534681242102b519cfc23f077',
64:'4db5ea0aad7a733d7d9cc864849fd66522645a5d8e82073e7b80e1f5f7050ad2',
67:'2b7abca22d5e0a0dddbcd4089a9057a0474ccf0b2f4e5263adcc900935cbd040',
69:'f274726bb141e78a1c81e9c546f5f55354c9b5af1614b2eae7ed2cb3bdd57c11',
70:'493daf2167ef7e372478b335dcf31caca3e50f131b2f5e936d1776dfccd488a0',
71:'2580e197000863fa4f7071d396d7b4a92cfd1492e5a9ad73e1d59ddc679dc77d',
72:'9404e3762084aeee4fc8527d0a9f0bc123ddc31a624ad620120fbcabb04ad5b8',
73:'ed016733f2065a11aaa9f85d402734b9792b3745813a43984062b793d8cd015f',
74:'d20cbb4aaf997b8d47b63271b39cc33c093e3373deef8ef76c32d5d8d550af20',
75:'a263e4f19d47c10da0f5fa8716eba992848095dd931947f55dcdef09a3a9babe',
76:'22b8dd597fd805edc1c6978dc43dcab0deee309c0f27eeb713068503d2b913e4',
77:'aec55fd84753fa0d8e0de9c8f1453522dc359a46873367e789b21bf3110d4dce',
78:'ed769c4528b025a045e9cfee073f2b978ee5e25b18ec3a51673fd090aa9c8874',
79:'d34bc69178c07e3aec6006fc0925dff7d5df757190e5849c28e296f8790dfff0',
80:'e35592e3e1eb264a93b704e8a218d365fcfaadda06f535d43feff04ebb2b26dd',
81:'6917a759a3e2eb8f4a7a1ab2532af2f2e5fdf237184a9d0b7971ecf6afd4c4ca',
82:'8230304d14959285ea4e21636b0809c761cbdb6a6347799803665db6677cbd9a',
83:'9aa8eec1c439c91066f74c23c27c0bb20d6b491a44475bebe8daddc06abd478a',
84:'d49eb6abc5f74dde1cdd67ee7ee8e7395b7a96687e1094ef0e2bfa26cb7f117b',
85:'df1f406692b6459aff94e853cf94ee75289244cd1b07bc9a921595b625831b4d',
86:'9b841799f54e33c19704037e5f1a533e87e56e95d4c5038279e9a11af22c38e7',
87:'193825b6c707fa10717879eb29a69892b49120305eaa2241c11559622bedfaac',
88:'eaa1a9bb170a24e8bb867eed9644ba783b7e7632a26ba8b17dad818ae9a4ab80',
89:'c12f68ff70505eace35b2c1efd68dbfa3d57840040c49349fe54ba031b64f22a',
90:'0d21b8c5865822cde9ef1939eda7f52ec9196d28849023d326fef6c8771664c2',
91:'93291bfa4bdd7958cc4e823e695132e71c58f583c40e1fa31cfb7001a83a3a08',
92:'da043b6808bed0a2a4a4fdfcdde0ec9bd529cc4b31722ed71ee5a7df1352b4f0',
93:'2a86dbe51696953ddca6a6f6819cb6873329344ca18a1bae6c22940f77855246',
94:'e1f4ca3225bd3948e51a4b89a88f7441cc7c59429b180ddd9e3326c8fbc06583',
96:'23b9cbcd9a5cea9263bd422c67b1860c1408f24357dbda8ca187e313924f618c'
}
PLAN = {}
def plan(n, title, text):
    PLAN[n] = (title, text.strip())

plan(57, '[Phase 4] 重构既有代码与整理新算法接口（不实施新算法）', '''
## 本轮边界

2026-10-11 用户进一步明确：当前是旧代码重构阶段。拉普拉斯等新增算法尚未开始编写，本阶段只提供它们的接口信息，不引入算法实现。上一版“旧代码迁移完成后继续交付新算法”的阶段目标在此收回。

当前交付只有两类：① 已有代码的必要内部重构、调用者迁移和兼容回归；② 新算法未来接入所需的输入、输出、数学条件、所有权、错误和资源契约资料。后者可以使用手算预期、读取示例或明确标注的契约替身，不能被写成真实算法、谱求解或学习功能已经交付。

不加入 Laplacian、Dirac、持久拉普拉斯、同调投影求解器、新非图距离、特征算法或训练运行时；不为了占位公开空 trait、万能结果容器、返回伪数据的生产接口或未实现的公共方法。未来实现另开明确实施任务，不进入本轮重构硬前置。

## 依赖的含义

原生 blocked-by 只表达当前任务确实需要的交付。父子关系用于归类，参考资料、未来扩展、同一框架内的关联和最终交叉核对不自动成为阻塞。移除“只完成父 Issue 的某个 M/N/L 分段就算解除 block”的模糊规则；本轮任务已缩窄到可以独立完成的范围。

### 旧代码主线

| 任务 | 直接硬前置 | 本轮出口 |
| --- | --- | --- |
| #71 迁移方案 | #63 | 文件/符号/调用者的真实处置清单，不含五类新算法研究 |
| #87 旧行为基线 | #63 | 固定现有版本与回归样例；先于迁移完成 |
| #72 PH 重构 | #71、#87 | 已确认源码调整及旧入口回归 |
| #78 / #81 / #84 契约核对 | 各自仅 #61 | 结果、比较、分析的既有契约及未来信息说明 |
| #79 结果迁移 | #71、#78、#87 | 迁移已有结果与消费者，无谱结果实现要求 |
| #82 距离迁移 | #71、#81、#87 | 迁移已有图距离，无新增非图算法要求 |
| #85 分析迁移 | #71、#84、#87 | 迁移已有描述子，不等待距离或外部训练 |
| #94 重构总验收 | #72、#79、#82、#85 | 同版本集成回归与兼容/成本记录 |

#72/#79/#82/#85 可按实际修改文件并行。共享文件变化用明确的实施 PR 依赖和同 PR 测试处理，不把所有消费者永久串成一条等待链。逐 PR 回归仍是实施任务自己的验收，不再把 #87 或 #94 同时用作实施前置与实施后总验收。

### 新算法只整理接口资料

#73 Laplacian/谱、#74 Dirac、#75 F2 同调算子均只依赖已完成的 #60/#61，彼此独立。#74 不再等待 #73 的算法；共同内积、基与伴随约定是资料协同关系。#88 核对 #73/#74 的单尺度交接；#89 核对 #69/#73 的阶段/持久谱接口；#90 核对 #69/#75 的窗口/传输接口；#91 只依赖 #81；#92 只依赖 #84。上述交接均不调用未实现的算法，也不阻塞旧代码重构 #94。

#76 广义拓扑和 #93 学习反馈保留为后续候选，不设置本轮硬阻塞；不开工不等于失败或功能交付。现有父子关联保留为导航，候选未关闭不应被解读为旧代码不能验收。

## 已有交付与状态

#58–#63 的完成状态不重置。#66 的建模交付和 #69 的既有阶段实现仍沿用原 PR；#69 只保留 #63/#68 为前置，#66 的示例交付不再阻塞阶段代码。两份指南的组合回归由实际合并 PR 承担。#26/#49 的准备实现、依赖和性能/review 门槛保持原归属；#56 实验候选不借此次重构进入默认路径。

本次不更改 Issue 开闭、assignee、标签、Project 自定义字段或原生父子归属。已关闭 #95–#97 不重新打开；移除已关闭 #96 对未完成 #94 的现行硬依赖，本轮采用与退出由 #94 和此总计划记录，历史关闭不当作本轮通过证据。

## 验收

旧代码用真实测试和必要的同版本成本检查验收；接口资料按定义完整性和交接一致性验收，两者分开记账。原有快速路径、素数域、覆盖、维度、代表、错误/预算与所有权语义保留。新增算法是否完成不再是本轮任何任务的验收条件。总体记录分别列出“旧代码已迁移 / 接口资料已说明 / 算法尚未实现”，不得相互替代。

此次只维护执行计划和 GitHub 原生依赖，没有实施产品源码重构或运行新算法。下方旧记录完整保存，但旧实施要求和旧硬前置不再作为当前规则。
''')
plan(64, None, '''
## 本轮范围说明

当前阶段保留 #65/#66 已有的数据对象复用、读取示例与原 PR 交付。后续 Alpha、立方等构造算法仅保留输入/身份/扩展信息，本轮不新增其算法实现。

本父任务用于归类，无额外原生阻塞。研究完成、#66 的实现/CI/审阅状态继续以原记录和对应 PR 为准，本次不重置状态，也不把指南中的手算 Hodge 矩阵称为谱计算能力。
''')
plan(67, None, '''
## 本轮范围与依赖说明

保留 #68/#69 已经开展的借用阶段与同源包含交付，不重新实施或取消其已有代码。用户限制的是尚未开始的新算法；当前阶段不因此增加持久拉普拉斯、持久 Dirac、多参数求解或一般链映射算法。

#69 的直接硬前置调整为 #63 和 #68。#66 的建模指南不是其阶段源代码的技术前置；两份指南若有组合改动，在实际集成 PR 中检查，不再用 #66 的 Issue 关闭状态阻塞 #69。已有实现、CI、review 和合并记录保留，未由本次依赖维护改写。
''')
plan(69, None, '''
## 当前范围与硬前置修正

只保留本项已经实现并进入原 PR #100 的有限单参数阶段读取、同源 inclusion、来源限制及其原定回归。已有代码与验证记录继续有效；本次不新写过滤算法，也不重置本项进度。

直接硬前置只保留 #63、#68，删除对 #66 的原生 blocked-by。理由：#66 的实际交付是已有 API 上的建模/读取示例，阶段源代码无须等待该示例 Issue 关闭。两 PR 的指南交叉修改仍须在合并版本完成组合检查，不能用关系解除代替 review 或测试。

本项可向 #89/#90 提供真实阶段/编号/映射资料；这些下游只整理接口，不要求在本轮调用新的持久谱或算子求解器。
''')
plan(70, '[P4-B03] 重构旧计算代码，整理新算法接口信息', '''
## 当前目标

用户正在推进 B03。本轮只重构已有计算代码，并把未来计算需要的数据和结果契约说明清楚。Laplacian、Dirac、持久谱及新的 F2 算子集成均不在本轮实现范围。

#71 只交付已有代码的迁移方案，#72 据此实施 PH 整理，#87 在改动前固定回归基线。#71 不再夹带五类新算法的完整研究，避免该 Issue 的原生阻塞把旧代码工作拖住。

#73/#74/#75 是独立接口资料任务：输入/输出、基/系数/内积或代价、所有权、来源、错误、精度/证书和资源边界。它们不依赖新的算法代码，不阻塞 #72。#76 只是非阻塞的后续模型候选。

已有共享 BoundaryInput/归约、隐式快路径、高维接续、代表与临界集是需要保留和审查的旧资产；不得重新申报或机械重写。新算法只在接口资料中说明缺口，不公开空实现，不增加求解器或运行时依赖。

## 完成记录

分别记录 #71/#72 的迁移及回归证据与 #73–#75 的接口资料；接口说明通过不写成算法能力已具备。父子归类不作为额外 blocked-by，#76 的后续研究也不作为旧计算代码验收门槛。
''')
plan(71, '[P4-B03R] 确定既有代码迁移清单与重构方案', '''
## 本轮唯一交付

接续 #59–#63 的既有结论，针对当前源码差异形成可执行的迁移清单。删除上一版 N 段“五类新计算研究”的完成要求；相关未来接口信息分别由 #73–#75、#78/#81/#84 承担，不成为本项关闭条件。

- [ ] 固定实际 main/实施 base，逐文件和符号列出来源、调用者、现有测试、目标责任、唯一负责 Issue/PR。
- [ ] 动作明确选择原位保留、内部提取迁移、薄兼容转发、清理失效私有重复或有依据的延期。已共享边界与已有快路径不重复实现。
- [ ] 覆盖 PH、diagram、diagram_distances、descriptors 及共用输入/控制/错误/导出；用同一套原子重构清单交接 #72/#79/#82/#85。
- [ ] 核对待合并 PR 的 base/head 与文件重叠；#26/#49/#56 及已有阶段/建模 PR 仍有独立采用门槛。
- [ ] 为实际缺口指定小批 PR、兼容处理、相关回归/成本检查和回退办法。当前没有必要变化的边界给源码及消费者依据，不制造空框架。

## 出口

清单和处置方案可独立完成；不等待 #73/#74/#75 或任何新增算法。#87 的基线工作可与本项并行，二者完成后开展实际迁移。本文不是重构实现或测试已通过的声明。
''')
plan(72, '[P4-B03P] 重构既有 PH 实现并保持旧入口兼容', '''
## 本轮工作

只实现 #71 中已确认的旧 PH 内部重构；以 #87 的固定旧行为为基线。无新算法家族、公共后端注册器或求解器引入。

- [ ] 整理输入适配、请求编排、数学内核和结果装配，迁移真实调用者并保持旧公开路径；清理失效私有重复，避免双份同义实现。
- [ ] 沿用 BoundaryInput、boundary::diagram、代数归约，以及外部输入验证/可信冻结 owner 读取的不同职责。
- [ ] 保留隐式 Rips/flag 快路径、F2 H0/H1→H2+ 接续、完整 clearing/虚拟配对、奇素数与溢出回退、代表元及近似源限制；不强制显式展开或重复归约低维。
- [ ] 保持 computed-empty/未计算、finite/essential/censored、coverage、context、代表关联及每 terminal 预算/错误语义。
- [ ] 每个实施 PR 带相关测试和旧新对照；若 dispatch、验证、拷贝、分配、生命周期或泛型调用成本变化，检查实际 base/head 的相关性能，不能用“公式未变”豁免。
- [ ] 更新受影响示例、doctest、架构/贡献说明及 worker。已经符合目标的边界提供 no-op 依据，真实缺口需要源码交付。

前置是完整且已缩窄的 #71 和 #87，不能再依赖“#71-M 局部通过”这种原生关系无法表达的条件。#94 是本项完成后的总验收，不反向阻塞本项；逐 PR 回归由本项自己完成。

#73–#75 只作为未来接口资料，均不阻塞本项，未来实现不得混入本轮重构。已有 PR 的合并/review 门槛保持，计划修改不等于代码交付。
''')
plan(73, '[P4-B03L] 仅定义 Laplacian/谱接口：输入、输出与数值契约', '''
## 范围

用户明确尚未开始新增算法。本项从“交付拉普拉斯与谱计算”改为接口资料，不实施 Laplacian 构造、矩阵作用、特征分解或新的数值后端。

- [ ] 写清输入所需的带序基、真实有符号边界、实/复系数、Gram/内积、权重含义、单尺度选择和来源身份。
- [ ] 给出未来入口与返回信息的拟定签名/说明，列明算子作用、矩阵和谱结果的不同输出要求；候选名称不自动成为公共 API。
- [ ] 说明维度、空矩阵、误差/残差、部分谱、核维度、失败和预算的契约，不把 F2 边界证书提升为实数证书。
- [ ] 用已有拓扑只演示必要数据能被读取；允许手算矩阵和标明为预期的输出形状，不运行或伪装新的谱算法。

验收是接口资料完整、与 #60/#61 一致及未实现项明确。算法、性能实验和外部求解器验证放在未来独立实施任务中，不是本轮出口。#74 与本项并行，共用约定在 #88 做资料一致性检查；本项不阻塞旧结果、距离、分析或 PH 迁移。
''')
plan(74, '[P4-B03D] 仅定义 Dirac 接口：分次边界、伴随与结果约定', '''
## 范围

只说明未来 Dirac 接入需要什么，不写 Dirac 算子、谱求解或持久 Dirac 实现。与 #73 的共用约定是资料协作关系，删除对 #73 的原生阻塞。

- [ ] 说明参与的分次空间、相邻维度、基顺序、有符号边界、内积/Gram 和伴随约定。
- [ ] 列出未来算子作用/矩阵/谱结果接口需要的身份、维度、误差及失败信息。
- [ ] 把平方关系和谱对称性写成未来实现的数学条件及手算预期，不能当作本库已经提供的功能。
- [ ] 明确空维度、截断分次范围、借用/显式物化、控制与依赖边界。

验收为接口说明与 #60/#61 相容。#88 核对本项和 #73 的交接一致性，不要求先实现任一算法；旧代码重构不等待本项。
''')
plan(75, '[P4-B03H] 仅定义 F2 同调算子接入契约，不移植求解器', '''
## 范围

外部 homology-operator 的已有研究可作接口参考，但本轮不移植或调用新的投影构造、最优求解、类查询和条形码算法，也不引入 Python/外部求解运行时。

- [ ] 描述三组带序基、两级边界、F2 系数、正坐标代价、来源和权重身份，给出从既有 incidence 读取到窗口的数据映射说明。
- [ ] 列明未来投影/作用/代表/类查询的输入输出及错误，区分可行、最优认证、未计算和失败。
- [ ] 明确所有权、预算、显式转换及代价变化后的兼容条件；已有 CriticalSetWorkspace 不等于本接口所指同调投影。
- [ ] 允许仅读取窗口、手算不变量或明确的契约测试替身，不将替身标为生产算子，不开展新 solver/性能实验。

验收是接入信息和数学限制完整。实际求解器集成另开未来实施任务，本项与 #73/#74 并行；#90 只核对窗口/传输接口资料。
''')
plan(76, '[后续候选] 广义拓扑模型与接口差异（不阻塞重构）', '''
## 后续范围

sheaf、路径/超图和 Mayer 等保留为后续候选。当前重构无需选定并实现其中一种，也不把全面文献研究作为重构门槛。

需要时可补充输入模型、边界规则、系数及无法由现有接口表达的信息，明确尚未决定部分。当前仅保留资料索引；正式研究或算法实施另行安排。

移除本项所有本轮硬前置；本项不得成为 #72/#79/#82/#85/#94 的阻塞。保留 open 和既有父子导航，不以“移出重构”冒充研究完成或功能交付。
''')
plan(77, '[P4-B04] 迁移既有结果表示，新增结果只描述契约', '''
## 当前交付

#78 核对既有结果读取/转换契约并记录未来结果需要的信息；#79 只迁移已有 diagram/result 与其消费者。删除本阶段必须出现一个真实谱结果或同调算子结果的要求。

未来谱、投影的接口资料可以描述基、域/Gram、误差/证书和可用性；不新增空生产者、假的谱输出或万能结果容器。真实新结果的实现跟随未来算法任务，不阻塞当前结果迁移。

本父任务通过子项聚合记录进度，不额外设置框架间 block。
''')
plan(78, '[P4-B04R] 确认既有结果契约及未来结果所需信息', '''
## 当前交付

接续 #61，明确旧 PersistenceDiagram/PersistenceData/PersistenceResult、逻辑维度视图、代表关联及借用/移动转换的兼容目标。只核对实际缺口，不重做已完成的结果研究。

- [ ] 列出现有读取和消费者必须保持的重数、维度、coverage、端点、context、代表关联与显式信息损失。
- [ ] 给 #79 一份可独立实施的契约/测试清单，不要求先出现新谱或投影的实例。
- [ ] 对未来谱、算子仅记录输入输出信息和未知项，允许引用 #73–#75 草案，不等待其算法或将引用设置为硬前置。

验收为旧结果迁移契约清楚、未来结果状态诚实；不要求新算法输出、数值求解或公共占位类型。只保留 #61 为直接前置。
''')
plan(79, '[P4-B04M] 迁移既有结果表示与消费者，不引入新算法结果实现', '''
## 当前实施范围

只实施 #71/#78 确认的已有结果/读取边界调整，并以 #87 为旧行为基线。上一版 N 段的谱和算子结果交付要求移出本轮。

- [ ] 处置 diagram 中的已有类型、逻辑维度读取、result/context 移动及代表关联，明确保留、内部提取、兼容转发与失效私有重复清理。
- [ ] 迁移实际 PH 装配、距离、描述子和示例调用者；公开路径与原数学语义保留，不维护两套同义结果。
- [ ] 检查 computed-empty/未计算、多重性、finite/essential/censored、coverage/context、代表索引及 into_parts/into_data/into_diagram 的信息保留/丢失。
- [ ] 同 PR 完成受影响测试、doctest 和必要分配/生命周期核对；未来接口资料只能说明形状和未知保证，不制造结果实例。

#73/#74/#75 不再是本项任何分段或关闭条件；本项可完整结束，无需等待未来算法。与 #72/#82/#85 有共享文件时按实际 PR 协调，不默认互相阻塞。最终集成由 #94 汇总。
''')
plan(80, '[P4-B05] 整理既有距离代码，新增比较仅留接口信息', '''
## 当前交付

#81 核对已有图距离及未来比较接口边界；#82 只迁移既有图距离、准备/solver 边界和调用者。不在重构阶段实现非图距离、新核函数、批调度器或新的样本矩阵算法。

比较入口仍需数学上可比的输入。未来对象需要的维度、单位、基、来源和样本身份只作接口信息，不构成已支持比较的声明。

#26/#49 的既有实现归属及 review/性能门槛保持；它们未被此次修改自动采用。父任务不额外设置框架间 block。
''')
plan(81, '[P4-B05R] 确认既有比较接口与未来比较边界', '''
## 当前交付

接续 #61 核对已有 raw/contextual/Execution 图距离入口与数学保证，给 #82 明确的兼容目标。删除“选出并交付一条新非图比较”作为本轮必须成果的要求。

- [ ] 明确既有 bottleneck、W1、W2 的度量约定、维度、完整性、essential、数值失败与数学无穷边界。
- [ ] 记录 facade、准备、solver 与结果消费者的既有责任，不接管 #26/#49 的未完成工作。
- [ ] 未来谱/过滤/算子比较仅列所需输入与身份、单位、基对齐、返回类型及不支持情况；尚未选择的数学定义保持待定，不为定义新算法扩大本项。

仅 #61 为直接硬前置，不等待 #71 五类研究、谱实现或新数据构造。接口资料可供 #91 核对，本项通过不能被写成新增度量已实现。
''')
plan(82, '[P4-B05M] 迁移既有图距离及调用者，不新增非图算法', '''
## 当前实施范围

按照 #71/#81 的方案迁移已有 diagram_distances，使用 #87 的旧行为基线。取消 N 段新增非图比较、通用成对矩阵和新批处理能力的本轮交付要求。

- [ ] 整理已有公开 raw/contextual/Execution 入口、准备/solver/数值回退与其调用者，保留或提取的动作有具体文件依据。
- [ ] 保留既有度量、coverage、计算维度、多重性、essential、field/声明尺度、数值错误及预算作用域。
- [ ] 修改必要旧调用者、测试、示例、worker，保持公开路径；按受影响分配、复制、分派和生命周期检查真实版本成本。
- [ ] #26 的逻辑准备、#49 的重复操作数准备保持原唯一实现和 main←#26←#49 依赖，不能复制或绕过原有 review/性能门槛。

删除对 #66/#73 的 block，不等待新的谱算法或真实非图结果。#79 的具体接口变动以同 PR/叠加 PR 协调；只需 #71/#81/#87 为稳定任务前置，集成回归归 #94。
''')
plan(83, '[P4-B06] 整理既有分析代码与可选学习接口边界', '''
## 当前交付

#84 记录已有描述子契约与未来外部数据交接信息；#85 迁移已有描述子和分析入口。本轮不实现新特征算法、训练器、分类/聚类算法或外部学习运行时集成。

数据可以直接进入未来外部工具，不强制 PH 或距离。当前只说明这种接口关系及样本/基/尺度身份，不要求真实训练或模型效果。

已有可选 CriticalSetWorkspace 保留计算层责任，不能因“学习框架”机械搬迁。#93 反馈研究是非阻塞后续候选。
''')
plan(84, '[P4-B06R] 明确已有描述子契约与外部数据交接信息', '''
## 当前交付

以 #61 的结果语义为基础，给 #85 确定旧 betti_curve、finite_lifetime_summary 及其消费者的读取/错误/范围契约。

- [ ] 保留 computed-empty/缺失维度、网格、coverage、essential/censored 排除计数、空值及数值规则。
- [ ] 明确未来数据交换中的样本顺序、拓扑/基身份、特征形状、单位和来源；未决定的编码/绑定方案记录备选，不引入新的序列化或 ML 依赖。
- [ ] 说明直接数据/过滤到外部工具的输入信息和未支持点，给 #92 文档级交接案例。

不要求运行训练、实现新特征或选择并证明一个学习算法。删除对 #71 五类新计算调研的硬依赖，只保留 #61。
''')
plan(85, '[P4-B06M] 迁移既有描述子与分析入口，不引入学习运行时', '''
## 当前实施范围

只落实 #71/#84 确认的旧 descriptors/分析入口整理，并以 #87 的现有行为为基线。上一版 N 段的新分析算法、外部训练和学习接续运行要求移出当前重构。

- [ ] 迁移真实旧调用者，保留既有描述子实现及公开路径，不另造同义版本，不重新执行 PH 或读取归约临时数组。
- [ ] 验证 computed-empty/缺失维度、严格递增有限网格、coverage、essential/censored 排除计数、空输出和单位/来源解释。
- [ ] 更新旧示例、doctest 和受影响测试；读取、排序、分配或复制变化按同版本范围检查。
- [ ] 保留已有 optimization 的独立责任及其与代数/复形的关系，不搬入训练层，不把临界集当作同调投影。

只依赖 #71/#84/#87，不再整体依赖 #79/#82 或未实现谱。真实共享类型变动在对应 PR 协调并由 #94 集成检查；未来学习接口资料由 #84/#92 负责。
''')
plan(86, '[P4-C] 旧链路回归与未来接口交接核对', '''
## 两类检查分开

#87 只建立改动前可独立完成的旧链路基线；#72/#79/#82/#85 每个实施 PR 自带回归；#94 在实际旧代码迁移之后汇总最终集成。删除“#87 既是前置基线又等待所有迁移关闭”的混合角色。

#88–#92 只核对未来单尺度、多尺度、算子、比较和学习的数据交接资料。可使用已有读取能力、手算预期与标注为契约替身的样本；不要求跑通尚未实现的 Laplacian、Dirac、持久谱、算子求解或训练。

#94 不被 #88–#92 阻塞；旧代码回归与接口资料完成度分别汇报。#93 是非阻塞后续研究。父子关系只保留导航，不自动要求未来算法完成后当前重构才能验收。
''')
plan(87, '[P4-C1] 建立重构前旧主线回归基线（不等待迁移完成）', '''
## 可独立完成的基线任务

本项只承担重构前基线，不再同时承担所有迁移结束后的总验收。直接前置仅 #63；删除对 #69/#72/#79/#82/#85 的 block。

- [ ] 固定当前已交付 main SHA、既有测试/样例及命令，覆盖 Rips/flag/外部复形、代表元可选请求、图距离和已有描述子。
- [ ] 记录适用范围内的旧入口输出、维度/重数/端点/coverage/context、代表关联、错误和预算行为；保留失败与限制。
- [ ] 给出重构 PR 可重复使用的输入/断言与相关成本边界。优先复用已有 oracle/worker，不因建基线重复所有历史昂贵研究。
- [ ] 检查三种既有消费方式：图比较、图分析、距离后消费。不得借此增加新的分析/距离算法。

交付基线包或可复现记录即可完成本项。之后 #72/#79/#82/#85 在各自 PR 内对照该基线，最终集成由 #94 负责；未来算法和接口资料均不影响本项退出。
''')
plan(88, '[P4-C2] 核对单尺度谱与 Dirac 接口交接，不实现计算', '''
## 接口级核对

仅消费 #73/#74 的接口资料，检查同一已有拓扑的数据能否按一致的基、方向、域、Gram/权重和来源交给未来两类计算。

交付数据读取说明、输入/输出形状表、合法/非法与空维度案例，以及平方关系适用条件的手算预期。允许仅使用已有边界读取或明确标记的契约替身。

不写或运行 Laplacian/Dirac 求解器，不实现谱分析，不声称功能或数值验证已完成。删除对 #66/#79/#85 的运行交付依赖；本项不阻塞旧代码重构 #94。
''')
plan(89, '[P4-C3] 核对阶段映射与持久谱接口，不实现持久拉普拉斯', '''
## 接口级核对

基于 #69 的既有阶段/同源包含资料和 #73 的谱接口信息，描述未来持久谱所需的阶段身份、基、边界、包含映射、内积及返回信息。

给出有限阶段的交接表、映射/截断/单位不兼容案例，说明普通逐阶段谱序列不自动提供持久谱语义。既有阶段读取可用于确认数据可得，未实现的谱结果只能作为手算预期或明示契约样本。

取消本项实现持久构造、求特征值、外部谱对拍和接入分析算法的要求。只依赖 #69/#73，不再依赖 #79/#85 的新结果/新分析运行；本项不阻塞旧代码重构 #94。
''')
plan(90, '[P4-C4] 核对 F2 窗口与传输接口，不引入算子求解', '''
## 接口级核对

基于 #69 已有阶段/包含与 #75 的 F2 算子接口资料，描述局部窗口、三组基、F2 边界、正代价、来源和阶段映射怎样交接。

交付窗口/映射字段对照、边界关系的手算不变量、同源/跨源与权重混用的拒绝条件，区分 selected representative 信息和最短类/最优证书。

不移植投影构造、类查询、传输 reducer 或算子条形码；不为接口演示调用外部 Python solver。原先两条真实算法运行例移至未来实施范围。

只依赖 #69/#75，移除 #79/#85 的算法结果/分析要求；本项不阻塞旧代码重构 #94。
''')
plan(91, '[P4-C5] 核对比较入口与样本身份，不新增比较算法', '''
## 接口级核对

只依赖 #81，说明已有图距离和未来数据/过滤/计算结果比较各自需要的输入、单位、维度、基或样本身份以及输出类型。

交付允许/不支持/未定义表和行列样本对照示例。已有图距离可作为真实调用例；未来对象仅用明确的契约样本，未知数学定义保留未知。

不要求真实谱输出、非图距离实现、通用相似性矩阵或新分析算法。删除对 #69/#73/#79/#82/#85 的整体硬依赖，必要参考仅保留文档链接；本项不阻塞 #94。
''')
plan(92, '[P4-C6] 定义数据到外部学习的交接契约，不运行训练', '''
## 接口级核对

只依赖 #84。记录已有拓扑或过滤数据怎样提供原编号、带序基、边界、特征形状、样本顺序、尺度和来源给未来外部工具。

交付数据字段/所有权/编码候选、维度和身份不匹配案例、明确的外部运行时边界。可用已有读取结果或静态样本验证形状，不要求增加绑定、序列化框架或新的特征算法。

删除实际外部前向/小训练成功的本轮验收，也不依赖 #73 谱计算、#69 新阶段运行或 #85 的学习集成。直接学习路径不强制 PH/距离，本项完成仅代表接口信息完整，不代表训练能力已提供；不阻塞 #94。
''')
plan(93, '[后续候选] 学习反馈与可导性研究（不阻塞重构）', '''
## 后续范围

学习反馈、可微过滤和优化算法均不是当前旧代码重构的交付。保留原研究线索，正式选题和实施另行安排。

本轮至多记录参数更新、身份失效与核心不能反向依赖学习运行时的边界；不要求梯度/训练实验或新增计算。

删除本项本轮硬前置，保持 open 与原父子导航，不成为任何旧代码迁移或 #94 的阻塞。暂不开展不等于数学不可行，也不等于已经完成研究。
''')
plan(94, '[P4-C8] 汇总既有代码重构验收，不等待新算法实现', '''
## 明确的末端验收

本项只在 #72/#79/#82/#85 的旧代码迁移完成后汇总同版本集成。它不再同时承担实施前必须关闭的基线任务；改动前基线是 #87，逐 PR 检查是各实施任务自己的验收。

- [ ] 固定实际集成 SHA 和各迁移 PR，核对文件/符号处置、公开路径、旧调用者和失效私有重复清理。
- [ ] 同一版本复核既有直接/外部/隐式输入、PH、代表、图结果、距离、描述子、错误和预算作用域；旧源/维度/coverage 等含义保持。
- [ ] 按触及的分派、验证、拷贝、分配、所有权与 workspace 生命周期检查实际成本，区分测试、CI、外部核对、性能、合并与发布。
- [ ] 确认本次重构未夹带新的 Laplacian/Dirac/持久谱/投影求解器、非图算法或 ML 运行时；未来接口资料不冒充运行能力。
- [ ] 汇总已采用旧代码范围的通过项、失败/回退和限制，在 #57 记录退出。

删除对 #88–#92 的 block；其接口交接资料独立审阅，不要求真实新算法输出。#73–#76/#93 不得直接或经其他任务阻塞此重构出口。#66/#69 的既有独立 PR 仅在实际使用时做对应版本组合核对，不通过无关 Issue 关闭状态扩大本项。

已关闭的 #95–#97 不被本项重新打开或当作本次通过证据。总体报告将旧代码迁移和未来接口资料分栏记录，所有未实现算法继续明确为未实现。
''')
plan(96, None, '''
## 当前依赖整理说明

本项已在本轮维护前关闭，保持原状态与历史记录，不重新打开，也不据此宣布本次重构已验收。

移除本项对仍未完成 #94 的现行原生 blocked-by，避免已关闭的历史记录继续充当本轮待执行出口。当前旧代码重构的总验收和采用记录归 #94 与 #57；新算法实现不在当前阶段范围。
''')

HEADERS = {}
LAST_WRITE = 0.0

def api(path, method='GET', payload=None):
    global LAST_WRITE
    if not path.startswith('/') or path.startswith('//'):
        raise ValueError('Invalid repository-relative API path')
    if method != 'GET':
        time.sleep(max(0.0, 1.15 - (time.monotonic() - LAST_WRITE)))
        LAST_WRITE = time.monotonic()
    data = None if payload is None else json.dumps(payload, ensure_ascii=False).encode()
    req = urllib.request.Request(BASE + path, data=data, headers=HEADERS, method=method)
    try:
        with urllib.request.urlopen(req, timeout=45) as response:
            raw = response.read()
            return json.loads(raw) if raw else None
    except urllib.error.HTTPError as error:
        # Report API diagnostics, never request headers or credentials.
        print('API_ERROR', method, path, error.code, error.read().decode()[:1500], flush=True)
        raise

def list_api(path):
    result = []
    for page in range(1, 11):
        sep = '&' if '?' in path else '?'
        values = api(path + sep + 'per_page=100&page=' + str(page))
        if not isinstance(values, list):
            raise RuntimeError('Expected paginated list: ' + path)
        result.extend(values)
        if len(values) < 100:
            return result
    raise RuntimeError('Pagination bound exceeded: ' + path)

def dependencies(n):
    values = list_api('/issues/' + str(n) + '/dependencies/blocked_by')
    if any(d['repository_url'] != BASE for d in values):
        raise RuntimeError('Unexpected cross-repository dependency')
    return {d['number']: d for d in values}

def ancestors(n, graph, visiting=None):
    visiting = set() if visiting is None else visiting
    if n in visiting:
        raise RuntimeError('Dependency cycle at ' + str(n))
    out = set()
    for p in graph[n]:
        out.add(p)
        out.update(ancestors(p, graph, visiting | {n}))
    return out

def validate_graph(graph):
    if set(graph) != set(range(57, 98)):
        raise RuntimeError('Unexpected audit scope')
    for n in graph:
        if not set(graph[n]).issubset(graph) or n in graph[n]:
            raise RuntimeError('Unknown/self dependency')
        ancestors(n, graph)
    future = {73,74,75,76,88,89,90,91,92,93}
    for n in (71,72,78,79,81,82,84,85,87,94):
        if ancestors(n, graph) & future:
            raise RuntimeError('Future capability blocks migration: ' + str(n))
    for optional in (76, 93):
        if graph[optional] or any(optional in ps for ps in graph.values()):
            raise RuntimeError('Optional research has a hard dependency')

def digest(text):
    return hashlib.sha256(text.encode()).hexdigest()

def render(n, original, title):
    predecessor = '、'.join('#' + str(p) for p in NEW[n]) or '无；父子归类与资料引用不构成额外硬阻塞'
    run = 'https://github.com/' + REPO + '/actions/runs/' + os.environ.get('GITHUB_RUN_ID', 'unknown')
    header = MARKER + '\n## 2026-10-11 当前有效范围与依赖\n\n'
    header += '**直接硬前置：' + predecessor + '。**\n\n'
    header += PLAN[n][1] + '\n\n'
    header += '维护依据：用户本轮“重构阶段仅提供新算法接口信息”的要求；原生依赖只按本段前置维护。'
    header += '本次不更改本 Issue 的开闭状态、指派、标签或父子归属。原始正文逐字保存在下方；其中旧实施验收和旧依赖已被本段取代。\n\n'
    header += '[本次依赖维护运行](' + run + ')；[原生依赖 API 规范](https://docs.github.com/en/rest/issues/issue-dependencies)。运行是否完成以日志和回读记录为准。\n\n'
    return {'title': title, 'body': header + HISTORY + original + '\n\n</details>\n'}

def main():
    if os.environ.get('GITHUB_REPOSITORY') != REPO or os.environ.get('GITHUB_REF') != BRANCH:
        raise RuntimeError('Dedicated repository/branch guard failed')
    token = os.environ.get('GH_TOKEN')
    if not token:
        raise RuntimeError('Missing workflow token')
    HEADERS.update({'Authorization':'Bearer ' + token, 'Accept':'application/vnd.github+json', 'Content-Type':'application/json', 'X-GitHub-Api-Version':'2022-11-28', 'User-Agent':'p4-dependency-maintenance'})
    validate_graph(NEW)
    if set(PLAN) != set(HASHES):
        raise RuntimeError('Plan/hash allowlist mismatch')
    main_before = api('/branches/main')['commit']['sha']
    if main_before != EXPECTED_MAIN:
        raise RuntimeError('Main changed; review the baseline before writing')
    issues, current, patches = {}, {}, {}
    # Read and validate the entire scope before the first mutation.
    for n in range(57, 98):
        item = api('/issues/' + str(n))
        if item['number'] != n or 'pull_request' in item:
            raise RuntimeError('Unexpected issue identity')
        issues[n] = item
        deps = dependencies(n)
        observed = set(deps)
        if not observed.issubset(set(OLD[n]) | set(NEW[n])) or not (set(OLD[n]) & set(NEW[n])).issubset(observed):
            raise RuntimeError('Unreviewed dependency edits at #' + str(n))
        current[n] = observed
        print('BEFORE', n, sorted(observed), item['state'], flush=True)
        if n not in PLAN:
            continue
        original = item.get('body') or ''
        if MARKER in original:
            if not original.startswith(MARKER) or HISTORY not in original or not original.endswith('\n\n</details>\n'):
                raise RuntimeError('Cannot safely resume modified issue #' + str(n))
            original = original.split(HISTORY, 1)[1][:-len('\n\n</details>\n')]
        if digest(original) != HASHES[n]:
            print('BODY_UPDATED_SINCE_AUDIT', n, 'preserve_latest_text', flush=True)
            # Current text is preserved beneath the new scoped header; concurrent edits
            # are rejected by the fresh comparison immediately before each write.
        title = PLAN[n][0] or item['title']
        patches[n] = render(n, original, title)
        if len(patches[n]['body']) > 65000:
            raise RuntimeError('Issue body too large after history preservation')
    removals = [(n,p) for n in NEW for p in sorted(current[n] - set(NEW[n]))]
    additions = [(n,p) for n in NEW for p in sorted(set(NEW[n]) - current[n])]
    print('PLAN', json.dumps({'remove':removals,'add':additions,'body_targets':sorted(patches)}), flush=True)
    # First put the scoped task descriptions in place. Check for concurrent edits.
    for n, patch in patches.items():
        fresh = api('/issues/' + str(n))
        if fresh['title'] != issues[n]['title'] or fresh.get('body') != issues[n].get('body') or fresh['state'] != issues[n]['state']:
            raise RuntimeError('Concurrent issue edit at #' + str(n))
        if fresh['title'] != patch['title'] or fresh.get('body') != patch['body']:
            api('/issues/' + str(n), 'PATCH', patch)
        print('BODY_UPDATED', n, digest(patch['body']), flush=True)
    # Remove obsolete edges first so no transient cycle is introduced.
    for n, p in removals:
        found = dependencies(n)
        if p not in found:
            raise RuntimeError('Concurrent dependency removal')
        if found[p]['id'] != issues[p]['id']:
            raise RuntimeError('Dependency identity mismatch')
        api('/issues/' + str(n) + '/dependencies/blocked_by/' + str(issues[p]['id']), 'DELETE')
        print('REMOVED', p, 'blocks', n, flush=True)
    for n, p in additions:
        found = dependencies(n)
        if p in found:
            raise RuntimeError('Concurrent dependency addition')
        api('/issues/' + str(n) + '/dependencies/blocked_by', 'POST', {'issue_id':issues[p]['id']})
        print('ADDED', p, 'blocks', n, flush=True)
    # Fresh readback of every dependency and state; no optimistic success claims.
    actual = {}
    for n in range(57, 98):
        item = api('/issues/' + str(n))
        actual[n] = sorted(dependencies(n))
        if set(actual[n]) != set(NEW[n]):
            raise RuntimeError('Native graph verification failed at #' + str(n))
        for key in ('state','state_reason','closed_at','assignees','labels','milestone','parent_issue_url'):
            if item.get(key) != issues[n].get(key):
                raise RuntimeError('Non-scope state changed: #' + str(n) + ' ' + key)
        if n in patches:
            if item['title'] != patches[n]['title'] or item.get('body') != patches[n]['body']:
                raise RuntimeError('Body/title verification failed')
        elif item['title'] != issues[n]['title'] or item.get('body') != issues[n].get('body'):
            raise RuntimeError('Unaffected issue changed during maintenance')
        print('VERIFIED', n, json.dumps(actual[n]), item['state'], flush=True)
    validate_graph(actual)
    # Check both sides of representative native relationships.
    for n in (73,74,75,76,87,93):
        reverse = list_api('/issues/' + str(n) + '/dependencies/blocking')
        same_scope = sorted(x['number'] for x in reverse if x['repository_url'] == BASE and x['number'] in NEW)
        expected = sorted(k for k, ps in NEW.items() if n in ps)
        if same_scope != expected:
            raise RuntimeError('Reverse blocking verification failed')
        print('REVERSE_VERIFIED', n, json.dumps(same_scope), flush=True)
    main_after = api('/branches/main')['commit']['sha']
    if main_after != main_before:
        raise RuntimeError('Main moved concurrently; do not claim unchanged main')
    summary = {'issues_audited':41,'bodies_updated':len(patches),'dependency_subjects_changed':sum(set(OLD[n]) != set(NEW[n]) for n in NEW),'edges_before':sum(map(len,OLD.values())),'edges_after':sum(map(len,NEW.values())),'edges_removed_this_run':len(removals),'edges_added_this_run':len(additions),'acyclic':True,'future_interfaces_block_migration':False,'states_and_parents_preserved':True,'main_sha':main_after}
    print('VERIFIED_SUMMARY ' + json.dumps(summary), flush=True)
    note = '<!-- p4-native-dependencies-verified:20261011:v2 -->\n## 原生依赖与范围修订已回读验证\n\n'
    note += '已核对 #57–#97 共 41 个 Issue；更新 30 项正文/范围，调整 22 个 Issue 的原生 blocked-by。原图 78 条边，现图 45 条边；计划删除 54 条、增加 21 条，已逐项回读匹配，并验证无环。\n\n'
    note += '旧代码链路：#63 → #71 与 #87；#71/#87 → #72；#79 另依 #78，#82 另依 #81，#85 另依 #84；#72/#79/#82/#85 → #94。#73/#74/#75 仅为接口资料，均不直接或间接阻塞这些迁移任务；#76/#93 无本轮硬阻塞。#69 已删除对 #66 的依赖。\n\n'
    note += 'Laplacian、Dirac、持久谱、算子求解、新非图算法和学习运行均已移出本轮实现/验收要求。所有 Issue 开闭状态、指派、标签、父子归属保持不变；已关闭 #96 与待完成 #94 的硬依赖已移除。没有产品代码修改或 PR 合并，main 仍为 `' + main_after + '`。\n\n'
    note += '验证运行：https://github.com/' + REPO + '/actions/runs/' + os.environ['GITHUB_RUN_ID'] + '\n'
    existing = list_api('/issues/57/comments')
    if not any(c.get('body','').startswith('<!-- p4-native-dependencies-verified:20261011:v2 -->') for c in existing):
        comment = api('/issues/57/comments','POST',{'body':note})
        print('VERIFICATION_COMMENT', comment['html_url'], flush=True)
    with open(os.environ['GITHUB_STEP_SUMMARY'], 'a', encoding='utf-8') as f:
        f.write(note)
    print('MAINTENANCE_COMPLETE', flush=True)

if __name__ == '__main__':
    main()
