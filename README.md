# denovoSAVI

In-house python script for identifying de novo germline variants and new homozygous variants from patient-unaffected-parant trio data

# Status
Active Development

# Usage
* denovoSAVI takes four input files that can be retrieved from the output files of [SAVI](https://github.com/WangLabHKUST/SAVI).

* You will need mother-blood, father-blood, patient-blood, and patient-diseased-specimen bams (if you have) to run SAVI.<br />
* All SAVI should have been run with bams in exact order as shown below.<br />
* All input files for denovoSAVI listed below should have patient ID in the first column and include all chromosome information.<br />

input file 1 - PDfilter results from SAVI (i.e., report.coding.PDfilter.txt) ran with "merged-parant-bam, patient-blood-bam" <br />
input file 2 - PDfilter results from SAVI (i.e., report.coding.PDfilter.txt) ran with "merged-parant-bam, patient-blood-bam, patient-diseased-specimen" <br />
input file 3 - Coding variants from SAVI (i.e., report.coding.txt) ran with "merged-parant-bam, patient-blood-bam, father-bam, mother-bam" <br />
input file 4 - Coding variants from SAVI (i.e., report.coding.txt) ran with "merged-parant-bam, patient-blood-bam, patient-diseased-specimen, father-bam, mother-bam" <br />

* Note that input file 1 and 3 will be from the patients without patient-diseased-specimen data, while input file 2 and 4 are from those with patient-diseased-specimen.<br />

When the input files are ready, you can run denovoSAVI by:
```
python denovosavi.py input1 input2 input3 input4
```

denovoSAVI will output the number of de novo germline mutations (DNMs) within the cohort as well as a histogram showing the number of DNMs per patient.<br />
It will also create pie charts showing the distribution of variant types.


# Contact
For technical questions, contact Yoonhee via email: ynam@connect.ust.hk
