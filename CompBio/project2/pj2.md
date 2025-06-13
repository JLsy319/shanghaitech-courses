# PurZ/PurZ0/PurA 同源序列搜索及 profile hmm 的建立

Z基因组，即ZTCG基因组，在基因组中Z完全取代了A。目前已知几百种噬菌体具有Z基因组。Z基因组的生物合成通路在2021年被两个研究团队发现（Science,  2021, 372(6541):512-516.  Science.  2021,  372(6541):516-520）。 在Z基因组生物合成通路中，PurZ是一个关键蛋白，它是嘌呤核苷合成通路中关键蛋白PurA的同源物。PurZ和PurA最关键的区别，在PurA中参与催化的Asp13 （按E  coli  PurA的序列编号）被Ser所取代。后续工作中，一种PurZ的变体PurZ0被发现 （Nat Microbiol. 2023 Jul;8(7):1330-1338.），它是PurA到PurZ进化的中间体。

1. 从文献/UniProt/PDB出发，找出五条PurZ序列，五条PurZ0序列， 五条PurA序列。分别提交名为PurZ.fasta，PurZ0.fasta和PurA.fasta的文件。
2. 请尝试在云桌面上本地安装BLAST。按顺序记录下来你安装过程中的每一个命令并提交这个命令集合文件。下载nr库，从PurZ0.fasta中选一条PurZ0序列作为查询序列，使用本地安装的BLASTp搜索nr库来获取库中所有PurZ0序列（注意不要在结果中混入PurA或PurZ，阐述你如何实现这一点）。你可以找到多少条？分别来自什么物种？提交一个Excel文件，表头为NCBI Accession，物种，coverage (alignment length/length of query)，  seqid (aka pident)，  E-value，bit score，fasta序列。  为了获取Excel的内容，你选择了BLAST的哪种输出格式？
3. 下载IMG/VR  v4，使用本地BLAST的makeblastdb命令为IMG/VR  v4建立索引。以（2）中用到的那条PurZ0序列作为查询序列，使用BLASTp搜索IMG/VR v4库，来获取库中所有的PurZ0序列（注意不要在结果中混入PurA或PurZ）。你可以找到多少条？分别来自什么物种？提交一个Excel文件，表头为IMG/VR Accession，物种，fasta序列，E-value，bit score。根据你找到的PurZ0的序列条数，简要评价NCBI nr 和IMG/VR v4。
4. 对于（2）和（3）中找到   的所有的PurZ0序列，使用mafft工具分别做多序列比对。在云桌面本地安装hmmer3，使用hmmer3中的hmmbuild工具，为两个多序列比对生成profile hmm并使用Skylign（http://skylign.org/）来获得这两个hmm文件的可视化结果。提交名为PurZ0_nr_hmm和PurZ0_imgvr4_hmm可视化结果图片。比较这两个图片，你可以得出什么结论？
5. 对于（2）中找到的所有的PurZ0序列，先使用cdhit工具以80% seqid进行压缩，使用压缩后的序列文件重复（4）的所有操作。提交名为PurZ0_nr_hmm_cdhit0.8可视化结果图片。比较PurZ0_nr_hmm和PurZ0_nr_hmm_cdhit0.8，你有什么结论？
6. 在PurA.fasta中选一条PurA作为查询序列，使用本地安装的BLASTp搜索nr库，使用输出结果的前300条，用cdhit以80%  seqid进行压缩，用压缩后的序列文件重复（4）的所有操作。提交名为PurA_nr_hmm_cdhit0.8的可视化结果图片。对PurZ0_nr_hmm_cdhit0.8和PurA_nr_hmm_cdhit0.8，分别标出底物结合口袋残基的位置，比较两者的异同。提交标记好的图片。
7. 对于IMG/VR v4数据库，使用jackhmmer进行搜索其中的PurZ0，和（3）中BLASTp的结果进行比较。提交搜索命令和搜索结果。
8. 还有哪些方法可以用来搜索IMG/VR v4数据库中的PurZ0？尝试选取一种，并提交搜索命令和搜索结果。这种方法和BLASTp以及jackhmmer有什么不同？
