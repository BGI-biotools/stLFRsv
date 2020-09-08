# stLFRsv
## Introuction:
Structure variation(SV) pipeline for stLFR co-barcode reads data.
This tool is applicable to stLFR technology and similar co-barcode data. Currently running on stLFR data, theoretically it is also applicable to other linked read data. You can analyze and test this tools by converting the data to the `read_id#XXX_XXX_XXX` format on read ID(`XXX_XXX_XXX` is the read barcode).

**Share barcode** information is used to detect breakpoint signals of structural variation (SV), such as: equilibrium translocation, inversion, deletion, duplication, and more complex structural breakpoints, which can be combined with CNV results and phase results. The accuracy of the structural variation is limited by the distribution density and length of linked reads on DNA molecules. For example, stLFR kits built with 1.5ng human DNA starting volume(about 30x significant depth) can produce data that can guarantee SV detection accuracy above 20K.

## Directory Structure:

* **bin**: necessary binarys to run the pipeline
* **data**: pre-build control list and black list
* **lib**: necessary libraris used in the pipeline
* **tools**: useful tools for the pipeline
* **src**: all source codes of the pipeline
* **LFR-sv**: main program of the pipeline 
* **run.sh**: example script

**Note**: this pipeline is build by `Perl` and `C++`. For the convenience of using (not need to install Perl and additional packages), the `Perl` scripts are packed into binarys. Please add the `lib` PATH to LD_LIBRARY_PATH when get a `error while loading shared libraries......` error. If you get the pipeline through `git clone`, give executive permission(`chmod +x`) to `LFR-sv bin/*` files or you can download the package form the [release page](https://github.com/BGI-biotools/stLFRsv/releases).
(Of course, you can replace the binarys by *.pl in the src with a little modification in LFR-sv.pl)


## Parameter Description：
###### Example:
```
./LFR-sv -bam example/L0.sort.rmdup.bam -out example/result -ncpu 20 -phase example/phase_out -bl data/bad_region_hg19_withchr.bllist -cl data/human_hg19_2000_20000_20000_0.9995_0.95_withchr.conlist -human Y
```
###### Parameters:
|  Parameter  |  Type | Description   |
| :------------ | :------------ | :------------ |
|-bam |<string> |  original sorted and markduped bam file,if the index dose not exist, will be created.\[necessary\]|
|-out |<string> |  output SV dir.\[necessary\](warning:if exists, the output dir will be cleaned first!!!!!)|
|-ncpu |<int>  |   thread number for running pipeline.[default 1]|
|-bar_th |<int> |at least N read pairs in one barcode.[default 8]|
|-seg_th| <int> |at least N read pairs in one segment.[default 4]|
|-gap |<int> |define the gap size which should be considered as different segment.|
|-size |<int>| output SV length.[default 20000]|
|-is |<int> |proper IS size for read pair library, read pairs with too large IS will be abandoned.[default 300]|
|-bin |<int>| bin size for cluster the segments.|
|-merge1 |<int>| N continue outline bins could be considered as the same break point and will be merged into one evidence.|
|-merge2 |<int>| SVs nearby under N binsize will be considered as one event.[default 5]|
|-mmax |<int> |the max SVs allowed in one event.[default 4]|
|-low |<int>|lowest shared barcode counts threshold.[default 4]|
|-sd |<int>| break ends with a depth higher than avg_dep+N*sd will be considered as candidates.[default 3]|
|-p_th |<float> |break ends significantly high with P value lower than this threshold will be considered as candidates.[default 0.1]|
|-phase |<string> |formatted phase result directory including phased barcode and region by chromosome.[default NULL]|
|-bl| <string> |black list file(BED format).[default NULL]|
|-cl| <string>| sorted control list file(BEDPE format).\[default NULL\](Be sure the chromosome and position are sorted in one line!!!)|
|-sc |<int>| allow max sv counts for the same position in one direction.[default 4]|
|-human| <Y/N>| for Homo sapiens,keep only [1234567890XYM] chromosome.[default N]|
|-qc1| <float>| valid read pair ratio for SV detection.[default 0.60]|
|-qc2 |<float>| average read pair count for one barcode.[default 30]|
|-qc3 |<float>| average segment end count for one bin.[default 15]|
|-sp |<float>| sample percentage for DNA fragment length statistic.[default 0.2]|
|-cn| <int> |sample count for read pair distance statistic.[default 20000000]|
|-rlen| <int> |read length of one read.[default 100]|
|-mlen |<int>| physical limit for the long DNA segment.[default 400000]|
|-help| |Show help message.|

## Result file type（by generate order）：
**sbf file**  
the binary segment file which generated form bam with the gap,bar_th,seg_th parameter  
segment binary format：  
`(bar:int64)8byte(index:int32)4byte(contig_name:char)32byte(start:int32)4byte(end:int32)4byte(pe_count:int32)4byte`  
the `index` means the Nth segment in this barcode，and the `bar` is combined barcode which constructed by `20bit_20bit_20bit` corresponding to `XXX_XXX_XXX`(a simple read script at /src/bar-sort/read-sbf.pl)  
**bfi file**  
binary index file for random access on sbf file  
index binary format:  
`(bar:int64)8byte(offset:int64)8bit`  
the `bar` is only valid in the lower 40bit corresponding to the 1st and 2st part of `XXX_XXX_XXX`(a simple read script at /src/bar-sort/read-bfi.pl)  
**gap file**  
include samples of gap between reads on segments of the largest contig, can be used to estimate parameters in `Auto Mode` 
**all.gap file**  
include all gaps form all segments(barcodes), can be used to estimate parameters in `Manual Mode` 
**stat file**  
include some statistical info from the bam   
**HQ.seg file**  
include samples of high quality segment size in sbf file for statistics   
**seg file**  
include all segment size for each barcode split by `4294967295` in bam file   
**freq file**  
include the possibilities that one barcode could cross a certain `N bp` gap without SVs for ALL and HQ barcodes   
**sin file**  
single end cluster file，the `sin.raw` is the original sin file, and `sin` file is the merged file by `merge1` parameter  
**lnd.all file**  
all suspect breakpoints with details.  
**lnd file**  
all passed breakpoints in `lnd.all` file that will be send to downstream analysis  
**lns file**  
segment link file, link the single end cluster to each segment  
**sln file**  
split link file, split the passed breakpoints by co-barcode on the same haplotype  
**judge file**  
passed breakpoints judged by several quality filter   
**filter file**  
filter by additional LFR rules   
**region file**  
mark breakpoints by black list file and control list file  
**final file**  
the final PASS SVs, and the `final.NoRegionFilter` is another final file that not considering the region markers  
**heatmap_plot folder**  
the heatmap of PASS SVs in final file. Or you can do it yourself using `/tools/plot_script` depend on which SVs you want  

## Noun explanation
![](https://github.com/BGI-biotools/stLFRsv/blob/master/graph/Fig1.png)
**segment** 
composed by several continue read pairs within a `gap` size, can be regarded as a DNA fragment without SVs. segment with read pairs more than `seg_th` is defined as `High Quality segment`   
**barcode**
include all read pairs with the same barcode tag. one barcode may contains one or more segment. barcode with read pairs more than `bar_th` is defined as `High Quality barcode`   
**single end cluster**
or may be called "single end breakpoint". several segments break at the same position and the same orientation, the cluster of these segments is defined as one `single end cluster`   
**SV breakpoint**
or may be called "pair end breakpoint". when two `single end cluster` are linked by co-barcode, it is called one `SV breakpoint`   
**SV Event**
a general SV such as Deletion，Inversion，Duplication or Transaction etc. one `SV Event` maybe constructed by one or two (more than two sometimes) SV breakpoint   