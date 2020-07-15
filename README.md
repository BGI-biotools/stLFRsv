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
* **Readme_CN.pdf**: manual(Chinese version)

**Note**: this pipeline is build by `Perl` and `C++`. For the convenience of using (not need to install Perl and additional packages), the `Perl` scripts are packed into binarys. If you get the pipeline through `git clone`, give executive permission(`chmod +x`) to `LFR-sv bin/*` files or you can download the package form the [release page](https://github.com/BGI-biotools/stLFRsv/releases).
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
|-gap |<int> |define the gap size which should be considered as different segment.[default 20000]|
|-size |<int>| output SV length.[default 20000]|
|-is |<int> |proper IS size for read pair library, read pairs with too large IS will be abandoned.[default 300]|
|-bin |<int>| bin size for cluster the segments.[default 2000]|
|-merge1 |<int>| N continue outline bins could be considered as the same break point and will be merged into one evidence.[default 5]|
|-merge2 |<int>| SVs nearby under N binsize will be considered as one event.[default 5]|
|-mmax |<int> |the max SVs allowed in one event.[default 4]|
|-low1 |<float/int>|single end barcode counts threshold, 0-1 float: higher than X percentage counts; 1> int: higher than X counts.[default 0.95]|
|-low2 |<float/int>| end to end barcode counts threshold, 0-1 float: higher than X percentage counts; 1> int: higher than X counts.[default 0.9995]|
|-ex1 |<float> |when low1 is a float of 0-1, exclude the bins which depth under ex1.[default 0.2]|
|-ex2 |<float>| when low2 is a float of 0-1, exclude the bins which depth under ex2.[default 0.2]|
|-phase |<string> |formatted phase result directory including phased barcode and region by chromosome.[default NULL]|
|-bl| <string> |black list file(BED format).[default NULL]|
|-cl| <string>| sorted control list file(BEDPE format).\[default NULL\](Be sure the chromosome and position are sorted in one line!!!)|
|-sc |<int>| allow max sv counts for the same position in one direction.[default 4]|
|-human| <Y/N>| for Homo sapiens,keep only [1234567890XYM] chromosome.[default N]|
|-qc1| <float>| valid read pair ratio for SV detection.[default 0.60]|
|-qc2 |<float>| average read pair count for one barcode.[default 30]|
|-qc3 |<float>| average segment end count for one bin.[default 8]|
|-rlen| <int> |read length of one read.[default 100]|
|-mlen |<int>| physical limit for the long DNA segment.[default 400000]|
|-help| |Show help message.|

###### Result file type（by generate order）：
**sbf file**

the binary segment file which generated form bam with the gap,bar_th,seg_th parameter.
segment binary format：
`(bar:int64)8byte(index:int32)4byte(contig_name:char)32byte(start:int32)4byte(end:int32)4byte(pe_count:int32)4byte`
the `index` means the Nth segment in this barcode，and the `bar` is combined barcode which constructed by `20bit_20bit_20bit` corresponding to `XXX_XXX_XXX`.(a simple read script at /src/bar-sort/read-sbf.pl)

**bfi file**

binary index file for random access on sbf file
index binary format:
`(bar:int64)8byte(offset:int64)8bit`
the `bar` is only valid in the lower 40bit corresponding to the 1st and 2st part of `XXX_XXX_XXX`.(a simple read script at /src/bar-sort/read-bfi.pl)

**gap file**

include all gaps form all segments(barcodes), can be used to gap statistic and parameter setting

**stat file**

include some statistical info from the bam 

**sin file**

single end cluster file，the `sin.raw` is the original sin file, and `sin` file is the merged file by `merge1` parameter

**lnd.all file**

all suspect breakpoints over low1 and low2 parameter.

**lnd file**

all passed breakpoints in lnd.all file that will be send to downstream analysis

**lns file**

segment link file, link the single end cluster to each segment

**judge file**

passed breakpoints judged by several quality filter 

**filter file**

filter by additional LFR rules like the `sc` parameter 

**region file**

mark breakpoints by black list file and control list file

**final file**

the final PASS SVs, and the `final.NoRegionFilter` is another final file that not considering the region markers

**heatmap_plot folder**

the heatmap of PASS SVs in final file. Or you can do it yourself using `/tools/plot_script` depend on which SVs you wan