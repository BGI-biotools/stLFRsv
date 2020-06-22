Usage: perl stLFRsv_merge_result.pl result1 result2 flank ratio break sample filt

result1: stLFRsv final result(.final/.final.NoRegionFilter)
result2: smoove final result(smoove.genotyped.vcf.gz)
flank: maxium length between two break points(only effect BND type)
ratio: minium ratio of(overlap_length/result1_length and overlap_length/result2_length) to merge(effect DEL/DUP/INV)
break: cutoff of length,only sv length >break in result and sv length <break in result2 will be merged(effect DEL/DUP/INV)
sample: sample name in final out
filt: 1-only include "PASS" in result1;2-inlcude "PASS" and "PASS|COMMON" in result


output:
merge.dels.vcf : vcf only include all dels from stLFRsv and smoove
merge.other.sv.vcf : vcf include other type svs(DUP,INV,BND)

