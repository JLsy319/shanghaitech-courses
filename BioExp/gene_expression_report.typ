#import "@preview/rubber-article:0.1.0": *
#import "@preview/whalogen:0.2.0": ce

#show: article.with()

#set page(
  paper: "a4",
  margin: (x: 2.54cm, y: 3.18cm),
)

#set par(
  leading: 0.5em,
  first-line-indent: (amount: 2em, all: true)
)

#set heading(numbering: "1.")

#maketitle(
  title: "中级生物学基因表达检测实验报告",
  authors: (
    "刘思昀",
  ),
  date: datetime.today().display("[day padding:none]. [month repr:long] [year]"),
)

= 知识掌握

== 简要介绍几种外源基因导入真核细胞的基本方法


+ 化学介导的转染: 利用化学物质辅助核酸分子穿过细胞膜

  - 磷酸钙共沉淀法: 其原理是将外源DNA与氯化钙溶液混合，然后缓慢加入磷酸盐缓冲液，形成磷酸钙-DNA的微小沉淀颗粒。这些沉淀颗粒被细胞通过内吞作用吸收，从而将DNA带入细胞内。

  - 脂质体介导的转染: 其原理是利用阳离子脂质体或非脂质聚合物包裹带负电荷的核酸分子，形成DNA-脂质复合物。由于细胞膜也带负电荷，带正电荷的复合物易于与细胞膜吸附并融合或通过内吞作用，从而将核酸释放到细胞质中。

  - 阳离子聚合物法: 利用聚乙烯亚胺或DEAE葡聚糖。这些带正电荷的聚合物能与带负电荷的DNA结合，形成复合物，然后通过内吞作用被细胞摄取。

+ 物理介导的转染: 利用物理手段在细胞膜上形成暂时的通道，以便核酸分子进入

  - 电穿孔法: 将细胞悬液置于高强度电场脉冲中，瞬间提高细胞膜的通透性，形成暂时的微孔，从而允许外源DNA分子直接进入细胞质甚至细胞核。

  - 显微注射法: 在显微镜下，利用极细的玻璃毛细管针直接将外源DNA注射到单个细胞的细胞质或细胞核内。

  - 基因枪法: 将包裹有外源DNA的微小金属颗粒 (通常是金或钨颗粒) 通过高压气体加速，直接射入靶细胞或组织中。这种方法可以克服细胞壁的障碍。

+ 生物介导的转染: 利用病毒具有高效感染细胞并将遗传物质导入宿主细胞的天然能力

  - 病毒转导: 利用经过基因工程改造的病毒作为载体，携带外源基因感染靶细胞。常用的病毒载体包括逆转录病毒、腺病毒、腺相关病毒等。这些病毒载体通常被改造成复制缺陷型，以确保生物安全性。

== 本次实验中所采用的脂质体介导的真核细胞转染技术的基本原理是什么？

脂质体介导的转染技术利用阳离子脂质体与核酸形成复合物，通过静电吸附和内吞的方式将核酸导入细胞，并最终通过内涵体逃逸将核酸释放到细胞质或细胞核中，从而实现外源基因的表达或功能调控。

== RNA抽提过程中最应注意的是尽量避免RNA酶的污染，实验过程中有哪些措施或方法可以实现这一目的？

由于环境中也存在大量的RNA酶，因此在提取RNA时，应该尽量创造一个无RNA酶的环境，包括去除外源性的RNA酶污染和抑制内源性RNA酶的活性。

主要采用 DEPC 除去外源性DNA酶，通过RNA酶的阻抑蛋白 Rnasin 和强力的蛋白质变性剂和异硫氰酸胍抑制内源性RNA酶。

= 实验概要

== 实验目的

+ 解真核细胞脂质体转染技术、RNA 抽提及 RT-PCR 技术的基本原理
+ 掌握转染技术、TRIzol 法提取哺乳细胞 RNA 及 RT-PCR 的实验方法
+ 了解流式细胞术 (FCM) 的基本原理

== 实验材料、试剂与仪器

=== 材料

+ 实验组质粒 pEGFP-PD1-N1, 浓度为 $536.5 " " mu l$
+ 阳性对照组质粒 pEGFP-N1, 由助教提供
+ Hela 细胞

=== 试剂

+ opti-MEM 培养基
+ lipofectamine 2000
+ DPBS
+ TRIzol
+ 氯仿
+ 异丙醇
+ 无水乙醇
+ DEPC $#ce("H2O")$
+ TAKARA PrimeScript#super[TM] RT Master Mix (Perfect Real Time) RR036A
+ Vazyme 2x Rapid Taq Master Mix P222-01-AA
+ baygene Regular Agarose BY-R0100
+ 1x TAE buffer
+ Tiangen 10000x GeneGreen Nucleic Acid Y2013
+ TransGen Biotech 100bp DNA Ladder BM301 

=== 仪器

+ 细胞培养箱
+ 荧光显微镜
+ eppendorf Centrifuge 5424R
+ Thermo Fisher Scientific Nanodrop One
+ Thermo Fisher Scientific ProFlex PCR System
+ Bio-rad 电泳槽
+ Azure Biosystems C150

== 实验方法

=== pEGFP-PD1-N1 及 pEGFP-N1 的转染

+ 转染前一天，将 Hela 细胞接种到6孔培养板中 ($4 times 10^5 "/ml"$); 转染时，细胞要达到 $80-95 #sym.percent$ 的融合，细胞培养于 $2 "ml"$ 含血清，不含抗生素的正常生长培养基中。 

+ 实验共包括三组: 实验组 (pEGFP-PD1-N1)，阳性对照 (pEGFP-N1) 和阴性对照组 (转染试剂 Lipo 2000); 阴性对照组只加入转染试剂 $490 " " mu l$ opti-MEM; 无血清培养基和 $10 " " mu l$ Lipo 2000 混合物，实验组及阳性对照组分别依照步骤 3 和4 进行。

+ 溶液 A: $240 " " mu l$ opti-MEM 培养基 + $10 " " mu l$ lipofectamine 2000 (总体积 $250 " " mu l$)，室温孵育 $5 "min"$. 如果采用 DMEM 做为 Lipofectamine 2000 的稀释液，必须在 $5 "min"$ 内同稀释的DNA混合。

+ 溶液 B: $(250 - x) " " mu l$ opti-MEM 培养基 + $4 " " mu g$ 质粒 $x " " mu l$ (总体积 $250 " " mu l$) (此为推荐剂量)。可以尝试不同质量质粒与不同体积 Lipo 的组合，找寻最优转染条件)。

+ 将溶液A 与溶液 B 混合，室温下静置 $5-20 "min"$.

+ 与此同时，将6孔板中的细胞用 DPBS 漂洗细胞两遍后 (每次 $1 "ml"$)，加入 $1.5 "ml"$ opti-MEM培养基。

+ 将溶液 A 与溶液 B 的混合液 (共 $0.5 "ml"$) 逐滴缓缓地加入孔中，摇动培养板，轻轻混匀。于 $37 #sym.degree"C", 5#sym.percent #ce("CO2")$ 培养箱中培养 $5#sym.tilde 6$ 小时。

+ 培养 24 小时后每孔加入 $1 "ml"$ TRIzol 收集细胞，一般在转染后 $48#sym.tilde 72$ 小时检测基因表达水平。 

+ 如果做稳定转染，更换完全生长培养基培养后 24 小时，即可以 $1: 10$ 或更高的稀释比例 (据细胞的生长情况) 接种到新的培养板，加相应的抗生素进行筛选。

+ 可采用 FCM 进行实验组转染效率的分析，则在进行转染实验时需要多设置一组实验组用于后续 FCM 分析GFP 阳性细胞比例。

=== TRIzol 法提取哺乳细胞 RNA (低温进行)

+ 样品处理: 吸弃旧培养液，每孔加入 $1 "ml"$ DPBS 洗涤 2 次。吸尽 DPBS 后，每孔加入 $1 "ml"$ TRIzol 反复吹打， 收获细胞裂解液于 $1.5 "ml"$ 离心管中，混匀，室温静置 $5 "min"$.

+ 加入 $0.2 "ml"$ 氯仿，振荡 $15 "s"$，静置 $2 "min"$.

+ $4#sym.degree"C"$ 离心 $12000 "g", 15 "min"$，取上清 (无色透明相约 $500 " " mu l$)  。

+ 加入 $0.5 "ml"$ 异丙醇，将管中液体轻轻混匀，室温静置 $10 "min"$.

+ $4#sym.degree"C"$ 离心 $12000 "g", 10 "min"$，弃上清。

+ 加入 $1 "ml"$ $75#sym.percent$ 乙醇，轻轻洗涤沉淀。$4#sym.degree"C"$ 离心 $7500 "g", 5 "min"$，弃上清。

+ 晾干，加入适量的 ($ 10-20 " " mu l$) DEPC $#ce("H2O")$ 溶解，取 $1 " " mu l$ RNA 稀释 10 倍后用 Nanodrop 测定浓度。

=== 逆转录反应 (RT)

使用TAKARA 公司的 PrimeScript#super[TM] RT Master Mix (Perfect Real Time) 进行逆转录反应，货为号RR036A. 每个样本取 $500 "ng"$ Total RNA，按照下列反应体系配制 RT 反应液 (配制过程请在冰上进行): 

#show table.cell.where(y: 0): strong
#set table(
  stroke: (x, y) => if y == 0 {
    (bottom: 0.7pt + black)
  },
  align: center
)

#align(center)[
  #table(
    columns: 3,
    table.header(
      [试剂],
      [使用量],
      [终浓度],
    ),
    [5x PrimeScript RT Master Mix],
    [$2 mu l$], [1x],
    [Total RNA],
    table.cell(colspan: 2)[$x mu l$ (每个样本$500 "ng"$)],
    [RNase Free d$#ce("H2O")$],
    table.cell(colspan: 2)[至$10 mu l$],
  )
]

$#sym.tilde 10 " " mu l$ 体系最大可使用 $500 "ng"$ 的 total RNA.

轻柔混匀 RT 反应液后，以下列条件进行逆转录反应: 

- $37#sym.degree"C", 15 "min"$ 逆转录反应
- $85#sym.degree"C", 5 "sec"$ 逆转录酶失活反应 
- $4#sym.degree"C"$ $#sym.tilde #sym.tilde$

=== 聚合酶链式反应 (PCR)

RT 样本包括: 阴性对照组 (NC), 阳性对照组 (PC, pEGFP-N1), 实验组 (pEGFP-N1-PD1)。每组扩增两个基因: β-actin (105bp); PD1 (171bp)，退火温度均采用 $57#sym.degree"C"$.

采用诺唯赞公司 2x Rapid Taq Master Mix (货号 P222) 进行 PCR 反应，体系总体积为 $20 " " mu l$: 

- cDNA $1 " " mu l$
- F-Primer $1 " " mu l$
- R-Primer $1 " " mu l$
- 2x Rapid Taq Master Mix $10 " " mu l$
- dd $#ce("H2O")$ $7 " " mu l$

β-actin (105bp) 反应条件为: 
- $95#sym.degree"C", 3 "min"$ 预变性后，进入循环 (循环 25 次):
  - $95#sym.degree"C", 15 "sec"$
  - $57#sym.degree"C", 15 "sec"$
  - $72#sym.degree"C", 5 "sec"$
- 然后 $72#sym.degree"C", 5 "min"$

$1 #sym.percent$ 琼脂糖凝胶电泳或置于 $4#sym.degree"C"$ 保存 (Marker: 100bp 或 DL2000).

== 原始数据

实验组质粒: $7.46 " " mu l$

阳性对照组质粒: $6 " " mu l$

逆转录反应 Total RNA: 实验组 $1.15 " " mu l$, 阳性对照组 $2.35 " " mu l$, 阴性对照组 $8.77 " " mu l$ 

== 实验结果

=== 转染后实验中三组细胞荧光表达情况

#figure(
  image("NC.jpg", width: 80%),
  caption: [
    阴性对照组细胞荧光表达情况，曝光时间: $9 "ms"$
  ],
) <NC>

#figure(
  image("PC.jpg", width: 80%),
  caption: [
    阳性对照组细胞荧光表达情况，曝光时间: $4.5 "ms"$
  ],
) <PC>

#figure(
  image("LSY.jpg", width: 80%),
  caption: [
    实验组细胞荧光表达情况，曝光时间: $4.5 "ms"$
  ],
) <LSY>

通过对于荧光结果的观察可以得到，质粒成功转染。

=== 三组样本的RNA抽提结果

#align(center)[
  #table(
    columns: 4,
    // caption: [RNA抽提后 Nanodrop 测定结果],
    table.header(
      [样本类型],
      [浓度 ($ "ng /" mu l$)],
      [OD260/OD280],
      [逆转录反应使用量 ($ " " mu l$)],
    ),
    [实验组 (稀释10倍)],
    [43.6], [1.34], [1.15],
    [阳性对照组 (稀释10倍)],
    [21.3], [1.28], [2.36],
    [阴性对照组],
    [57.0], [1.40], [8.78],
  )
]

=== 三组样本中 β-actin 和 PD1 基因的转录表达情况

#figure(
  image("gel1-report.png", width: 90%),
  caption: [
    第一次电泳结果, PC: 阳性对照组, NC: 阴性对照组, Exp: 实验组
  ],
) <gel1>

通过观察图@gel1，阴性对照组有正确的内参蛋白条带，同时不含有PD1基因表达条带，符合预期；阳性对照组不含有PD1基因表达条带，符合预期，然而预期结果需要含有内参蛋白的条带，结果中并没有观察到；实验组中含有PD1基因表达条带，符合预期，然而缺少内参蛋白的条带。

由于内参蛋白条带的缺少，因此使用RNA抽提样本再次进行 RT-PCR 实验，试剂添加量与第一次一致。(由于NC的RNA抽提样本没有剩余，因此没有重复阴性对照组的实验)

#figure(
  image("pd1-report.png", width: 60%),
  caption: [
    第二次电泳PD1基因表达结果, PC: 阳性对照组, NC: 阴性对照组, Exp: 实验组
  ],
) <pd1>

观察图@pd1，PD1基因表达情况为: 实验组有条带，阳性对照组无条带，符合预期。

#figure(
  image("actin-report.png", width: 60%),
  caption: [
    第二次电泳内参蛋白结果, PC: 阳性对照组, NC: 阴性对照组, Exp: 实验组
  ],
) <actin>

观察图@actin，内参蛋白表达情况为: 实验组有条带，而阳性对照组仍然没有条带。

这里借用了其他组的阳性对照组RNA样本进行 RT-PCR 实验，浓度为 $66.0 "ng/"mu l$，OD260/OD280 为 $1.60$，逆转录过程中添加量为 $0.8 " " mu l$.

#figure(
  image("gel3-report.png", width: 60%),
  caption: [
    第三次电泳内参蛋白结果, PC: 阳性对照组, NC: 阴性对照组
  ],
) <gel3>

观察图@gel3，结果显示阳性对照组的内参单板有条带，PD1基因表达没有条带，符合预期。

= 分析与讨论

+ 溶液 A 与溶液 B 混合后，可视情况决定静置时间，通常 $5 "min"$ 就可以保证转染效率。

+ 混匀细胞裂解液时，以液体重新回到较稀的流体状态即可。

+ 加入乙醇并离心后，弃上清剩余液体越少越好，以减少晾干的时间; 晾干时，开盖正立，观察到沉淀刚刚变为透明的时刻，加入 DEPC $#ce("H2O")$ 溶解 (最好时机)。

+ 对于实验结果中阳性对照组内参蛋白条带缺失的原因，分析如下: 

  - 收细胞时，观察到明显的荧光，因此转染和细胞培养的过程并没有问题，转染质粒也成功表达。
  - 由于 RT-PCR 是所有实验组和对照组同时出现的实验，因此可以排除这一部分的操作问题。
  - 考虑到RNA提取后，PC样本测定结果中OD260/OD280的比值为1.28，通常较纯的RNA样本比值在1.8-2.0之间，因此可能提取得到的RNA已经被污染或降解，导致其转录后无法扩增出所需的条带。(因为阳性对照组中不表达PD1基因，因此无法确认RNA提取的样本质量)
