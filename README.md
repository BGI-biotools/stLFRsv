## Content 
##### 	[Introuction](#jump1)
##### 	[Directory Structure](#jump2)
##### 	[Parameter Description](#jump3)
##### 	[Result file type](#jump4)
##### 	[Noun explanation](#jump5)
##### 	[Parameter details and Algorithm](#jump6)
##### 	[Result format explaination](#jump7)
##### 	[Tools](#jump8)
##### 	[Contact us](#jump9)

# stLFRsv
## <span id="jump1"> Introuction: </span>
Structure variation(SV) pipeline for stLFR co-barcode reads data.
This tool is applicable to stLFR technology and similar co-barcode data. Currently running on stLFR data, theoretically it is also applicable to other linked read data. You can analyze and test this tools by converting the data to the `read_id#XXX_XXX_XXX` format on read ID(`XXX_XXX_XXX` is the read barcode).

**Share barcode** information is used to detect breakpoint signals of structural variation (SV), such as: equilibrium translocation, inversion, deletion, duplication, and more complex structural breakpoints, which can be combined with CNV results and phase results. The accuracy of the structural variation is limited by the distribution density and length of linked reads on DNA molecules. For example, stLFR kits built with 1.5ng human DNA starting volume(about 30x significant depth) can produce data that can guarantee SV detection accuracy above 20K.

## <span id="jump2"> Directory Structure: </span>

* **bin**: necessary binarys to run the pipeline
* **data**: pre-build control list and black list
* **lib**: necessary libraris used in the pipeline
* **tools**: useful tools for the pipeline
* **src**: all source codes of the pipeline
* **LFR-sv**: main program of the pipeline 
* **run.sh**: example script

**Note**: This pipeline is build by `Perl` and `C++`. For the convenience of using (not need to install Perl and additional packages), the `Perl` scripts are packed into binarys. Please add the `lib` PATH to LD_LIBRARY_PATH when get a `error while loading shared libraries......` error. If you get the pipeline through `git clone`, give executive permission(`chmod +x`) to `LFR-sv bin/*` files or you can download the package form the [release page](https://github.com/BGI-biotools/stLFRsv/releases).
(Of course, you can replace the binarys by *.pl in the src with a little modification in LFR-sv.pl)


## <span id="jump3"> Parameter Description：</span>
###### Example:
```
./LFR-sv -bam example/L0.sort.rmdup.bam -out example/result -ncpu 20 -phase example/phase_out -bl data/bad_region_hg19_withchr.bllist -cl data/human_hg19_2000_20000_20000_0.9995_0.95_withchr.conlist -human Y
```
###### Parameters:
|  Parameter  |  Type | Description   |
| :------------ | :------------ | :------------ |
|-bam |\<string> |  Original sorted and markduped bam file,if the index dose not exist, will be created.\[necessary\]|
|-out |\<string> |  Output SV dir.\[necessary\](warning: If exists, the output dir will be cleaned first!!!!!)|
|-ncpu |\<int>  |   Thread number for running pipeline.[default 1]|
|-bar_th |\<int> |At least N read pairs in one barcode.[default 8]|
|-seg_th| \<int> |At least N read pairs in one segment.[default 4]|
|-gap |\<int> |Define the gap size which should be considered as different segment.|
|-size |\<int>| Output SV length.[default 20000]|
|-is |\<int> |Proper IS size for read pair library, read pairs with too large IS will be abandoned.[default 300]|
|-bin |\<int>| Bin size for cluster the segments.|
|-merge1 |\<int>| N continue outline bins could be considered as the same break point and will be merged into one evidence.|
|-merge2 |\<int>| SVs nearby under N bin size will be considered as one event.[default 5]|
|-mmax |\<int> |The max SVs allowed in one event.[default 4]|
|-low |\<int>|Lowest shared barcode counts threshold.[default 4]|
|-sd |\<int>| Break ends with a depth higher than avg_dep+N*sd will be considered as candidates.[default 3]|
|-p_th |\<float> |Break ends significantly high with P value lower than this threshold will be considered as candidates.[default 0.1]|
|-phase |\<string> |Formatted phase result directory including phased barcode and region by chromosome.[default NULL]|
|-bl| \<string> |Black list file(BED format).[default NULL]|
|-cl| \<string>| Sorted control list file(BEDPE format).\[default NULL\](Be sure the chromosome and position are sorted in one line!!!)|
|-sc |\<int>| Allow max sv counts for the same position in one direction.[default 4]|
|-human| \<Y/N>| For Homo sapiens,keep only [1234567890XYM] chromosome.[default N]|
|-qc1| \<float>| Valid read pair ratio for SV detection.[default 0.60]|
|-qc2 |\<int>| Average read pair count for one barcode.[default 30]|
|-qc3 |\<int>| Average segment end count for one bin.[default 15]|
|-sp |\<float>| Sample percentage for DNA fragment length statistic.[default 0.2]|
|-cn| \<int> |Sample count for read pair distance statistic.[default 20000000]|
|-rlen| \<int> |Read length of one read.[default 100]|
|-mlen |\<int>| Physical limit for the long DNA segment.[default 400000]|
|-help| |Show help message.|

## <span id="jump4"> Result file type（by generate order）： </span>
**sbf file**  
The binary segment file which generated form bam with the gap, bar_th, seg_th parameter.  
Segment binary format：  
`(bar:int64)8byte(index:int32)4byte(contig_name:char)32byte(start:int32)4byte(end:int32)4byte(pe_count:int32)4byte`  
The `index` means the Nth segment in this barcode，and the `bar` is combined barcode which constructed by `20bit_20bit_20bit` corresponding to `XXX_XXX_XXX`.(a simple read script at /src/bar-sort/read-sbf.pl)  
**bfi file**  
Binary index file for random access on sbf file.  
Index binary format:  
`(bar:int64)8byte(offset:int64)8bit`  
The `bar` is only valid in the lower 40bit corresponding to the 1st and 2st part of `XXX_XXX_XXX`.(a simple read script at /src/bar-sort/read-bfi.pl)  
**gap file**  
Include samples of gap between reads on segments of the largest contig, can be used to estimate parameters in `Auto Mode`.   
**all.gap file**  
Include all gaps form all segments(barcodes), can be used to estimate parameters in `Manual Mode`.   
**stat file**  
Include some statistical info from the bam.   
**HQ.seg file**  
Include samples of high quality segment size in sbf file for statistics.   
**seg file**  
Include all segment size for each barcode split by `4294967295` in bam file.   
**freq file**  
Include the possibilities that one barcode could cross a certain `N bp` gap without SVs for ALL and HQ barcodes.   
**sin file**  
Single end cluster file，the `sin.raw` is the original sin file, and `sin` file is the merged file by `merge1` parameter.  
**lnd.all file**  
all suspect breakpoints with details.  
**lnd file**  
all passed breakpoints in `lnd.all` file that will be send to downstream analysis.  
**lns file**  
segment link file, link the single end cluster to each segment.  
**sln file**  
split link file, split the passed breakpoints by co-barcode on the same haplotype.  
**judge file**  
passed breakpoints judged by several quality filter.   
**filter file**  
filter by additional LFR rules.   
**region file**  
mark breakpoints by black list file and control list file.  
**final file**  
the final PASS SVs, and the `final.NoRegionFilter` is another final file that not considering the region markers.  
**heatmap_plot folder**  
the heatmap of PASS SVs in final file. Or you can do it yourself using `/tools/plot_script` depend on which SVs you want.  

## <span id="jump5"> Noun explanation  </span>
![](https://github.com/BGI-biotools/stLFRsv/blob/master/graph/Fig1.png)
**segment**   
Composed by several continue read pairs within a `gap` size, can be regarded as a DNA fragment without SVs. Segment with read pairs more than `seg_th` is defined as `High Quality segment`.   
**barcode**   
Include all read pairs with the same barcode tag. One barcode may contains one or more segment. Barcode with read pairs more than `bar_th` is defined as `High Quality barcode`.   
**single end cluster**   
Or may be called "single end breakpoint". Several segments break at the same position and the same orientation, the cluster of these segments is defined as one `single end cluster`.   
**SV breakpoint**   
Or may be called "pair end breakpoint". When two `single end cluster` are linked by co-barcode, it is called one `SV breakpoint`.   
**SV Event**   
A general SV such as Deletion，Inversion，Duplication or Transaction etc. One `SV Event` maybe constructed by one or two (more than two sometimes) SV breakpoint.   

## <span id="jump6"> Parameter details and Algorithm  </span>
**bin, gap and merge1**   
These three parameters are the most important for the pipeline which define how to generate segment, how to cluster segment and judge the SV breakpoint.
![](https://github.com/BGI-biotools/stLFRsv/blob/master/graph/Fig2.png)
As shown in the pic above, the distance between read pairs follows a distribution which affected by the sequence depth, input DNA amount and DNA fragment length etc.
* bin: Usually use a 65% quantile of read-pair size. This means 65% read pairs should be in a distance of bin size. Meanwhile, 65% quantile of read-pair size keeps a reasonable precision of breakpoint position.
* merge1: Usually use a 93% quantile of read-pair size. When clustering segment ends using bin size, several clusters at the same orientation in nearby `merge1` bins should be form one real breakpoint and will be merged in `sin` file.
* gap: Usually use a 98% quantile of read-pair size. For several continue read pairs on the same chromosome with the same barcode, segments are generated by split these read pairs using gap size.

Pipeline offers `Auto Mode` and `Manual Mode` for these three parameters. 
If users understand their data well and know exactly the gaps distribution of their standard library, they could specify the parameters manually. This should be more efficient and  avoiding statistical errors.
If users didn't specify the parameters manually, `Auto Mode` will gather the gaps form the largest contig(chromosome) and generate these three parameters.   

**bar_th and seg_th**   
As described in Noun explanation, these two parameters are used to define `HQ` barcode and segment. Usually it is not recommended to set lower values than defaults unless your DNA fragment length is too short.

**SV detection size**   
`size` parameter set the SV size for the output. However, due to the `gap` and DNA fragment size, the pipeline will estimate the capability of SV detection size on specific data. Finally, pipeline will combine the `size` and the capability for output SV size(shown in the output message)

**Specificity and Sensitivity**   
Three parameters are used to adjust the specificity and sensitivity: low, sd and p_th.
![](https://github.com/BGI-biotools/stLFRsv/blob/master/graph/Fig3.png)
* low: The lowest co-barcode between two single end cluster bin which will be chosen as candidate.
* sd: As shown in pic above, one bin with segment single end depth higher than AVG+ sd*var is chosen as candidate.
* p_th: As shown in pic above, depth of one bin significantly higher than bins around with P-value lower than p_th is chosen as candidate.

Tips: Usually, when the data is poor(with QC warnings) or you just care SVs with large size or different contig(e.g. equilibrium translocation or other SV far greater than the capability for output SV size)，you could set a large `size`, lower `sd` and higher `p_th` to ensure the sensitivity with a little loss of specificity.

**black list and control list**   
* bl: The parameter which set the black list file. Black list region should contain gaps, low coverage, low complex and reference assembly issue regions which are likely to give false positive SVs in BED format.
* cl: The parameter which set the control list file. Control list region should contain segmental duplications, high population frequency and other systematical SV region which caused by aligner, sequencer and library method etc. The control list file should be provided in BEDPE format with inline sorted. The rule is:  

		First order: contig name. 
			When the contig names are different, remove `chr` form the name. If the remains are composed by numbers, then sorted by number. Otherwise, sorted by string. (from small to large) 
			When the contig names are same, go to the second order.
		Second order: position.
			When the left positions are different, sorted by the left position.(from small to large) 
			When the left positions are same, sorted by the right position.(from small to large）
		Examples:
			chr1 100000 300000 chr1 500000 800000
			chr4 100000 300000 chr5 500000 800000
			chr12 2200000 3000000 chr15 500000 800000
			12 100000 300000 19 500000 800000
			abc 100000 300000 cde 500000 800000
			m 1050000 3080000 n 500000 800000

Note: The pre-build control list file in `data` may be only applicable for the data which capability of SV detection size over: DEL(10k bp), other(30-40k bp) and bin size around 1500-2000 bp for humam beings. Other data may need to build a custom control list file using `tools/make_control`.

**Phase files**   
The formatted phase files could be generated by raw [HapCUT2](https://github.com/vibansal/HapCUT2) results with a tool in `tools/gen_phase`. The `phase` parameter set the path of phase files directory.

**Quantity control**   
There are four values in the `stat` file split by tab in the first line: Total valid read pair, HQ valid read pair, HQ barcode and HQ segment.
Three QC parameters was defined by:
* qc1: HQ read pair ratio = HQ valid read pair/Total valid read pair
* qc2: avg read pairs on one segment = HQ valid read pair/HQ segment
* qc3: avg segments end coverage in one bin = HQ segment*2/(ref_total_len/bin_size)

Note: The default values of the pipeline are tested and examined with Homo sapiens data of standard stLFRkit. If you use linked read data of other protocol or species, please make some test runs of standard samples to adjust the qc parameters.

**merge2 and mmax**   
This two parameters is about SV Event. SV breakpoints within a distance of merge2*bin are considered as one possible SV Event. However, one Event which the number of SV breakpoints over mmax may be a high frequency false positive region and will be abandoned.   

## <span id="jump7"> Result format explaination   </span>
There are two kinds of final result file, `final ` and `final.NoRegionFilter`, as described in Result file type above.
They have the same format below:   

|  Parameter  | Description   |
| :------------ | :------------ |
|EventID | The ID of a SV Event, may contains one or more SV breakpoints and start with 'S' character. And take notice that it's not a high confidence SV breakpoint but just a nearby SV breakpoint for indication only When start with '*S'.|
|SvID | The ID of each unique SV breakpoint which counld be tracked since the `judge` file.|
|BreakID1 | The ID of the first single end cluster which could be tracked since the `lnd.all` file.|
|BreakID2 | The ID of the second single end cluster which could be tracked since the `lnd.all` file.|
| ChrA | The contig name of the first single end cluster.|
| PosA | The contig position of the first single end cluster.|
| ChrB | The contig name of the second single end cluster.|
| PosB | The contig position of the second single end cluster.|
|ShareBarcode | The number of share barcodes between two single end clusters.|
|RealType | The link direction type of two single end clusters: R means Right and L means Left.|
|SimpleType | More human readable type. For example: RL for DEL, LR for DUP, LL or RR for INV.(Not applicable for complex SV cases)|
|QualityScore | Comprehensive SV breakpoint score for all QC filter.|
|ComprehensiveFilter | Simple judgment of SV breakpoint quality.|
|HeatmapFilter | The score of 2D share barcode matrix by Wilcoxon rank-sum test. "NULL" means there is not enough data for this judgment.|
|PhaseFilter | The judgment of barcode phase info.<br>The format is constructed with: `judge:(HP of BP1 - HP of BP2):phaseblock1:phaseblock2:P_value11,P_value12,P_value13,P_value21,P_value22,P_value23`<br>The P_values are P value of (1&#124;0), (0&#124;1) and (1&#124;1) for both SV breakpoint ends respectively by Fisher's exact test. "NULL:REGION" means this region is not phased, "NULL:COUNT" means there is not enough data and "NULL:UNPHASED" means that unphased barcodes are over 75%.
|MapQFilter | The high map quantity ratio(>10) of two single end clusters nearby.| 
| PairEndFilter | The pair-end support for this SV breakpoint.<br>The format is constructed with: `judge:+-,-+,--,++:PE support position`|
|LocalizationFilter | This filter judge the symmetry and localization of share barcodes base on DNA fragment length statistics. <br>The format is constructed with: `symmetry,localization1,localization2`<br> "SYM" means the share barcodes distribution is symmetrical on two single end cluster directions, otherwise it is unsymmetrical. "LONGER:xxxxx" means the SV end on its opposite direction is longer than the avg fragment length, which indicates a reliable SV breakpoint.|
|BlackList | Mark the black region on two ends of the SV breakpoint.|
|ControlList |Mark the control  region on the SV breakpoint.|
|SegmentCheck1 | Additional segment check by `sc` parameter.|
|SegmentCheck2 | Additional segment check by share barcode number with the same single end cluster.|
|SVchain | An experimental function to link the SV Events in one chain.|

## <span id="jump8"> Tools   </span>
**gen_phase**   
A tool that convert HapCUT2 result files to phase files which can be used by the pipeline.   
**make_control**   
A tool that help users to build their own control list file through a batch of `lnd.all` files form their own data.   
**merge_smoove**   
A tool that help users to combine the stLFRsv result and [smoove](https://github.com/brentp/smoove) result.   
**plot_script**   
A tool that help users to plot 2D share barcode heatmap from the `final` file or any other positions you want.   

## <span id="jump9"> Contact us  </span>
If you have any problem while using this pipeline, please add an issue at [issues](https://github.com/BGI-biotools/stLFRsv/issues) or send an email to guojunfu@genomics.cn.