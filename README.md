# dragenRNAseq

A workflow for calling SNVs on tumor-only or tumor-normal inputs in somatic mode

## Overview

## Dependencies

* [gsi modules : gencode 31](https://gitlab.oicr.on.ca/ResearchIT/modulator)
* [gsi modules : dragen-scripts 0.3](https://gitlab.oicr.on.ca/ResearchIT/modulator)


## Usage

### Cromwell
```
java -jar cromwell.jar run dragenRNAseq.wdl --inputs inputs.json
```

### Inputs

#### Required workflow parameters:
Parameter|Value|Description
---|---|---
`fastqInputs`|Array[InputGroup]|Input structure with fastq files and read group strings
`outputFileNamePrefix`|String|Prefix for output files
`reference`|String|The genome reference build. For example: hg19, hg38, mm10


#### Optional workflow parameters:
Parameter|Value|Default|Description
---|---|---|---


#### Optional task parameters:
Parameter|Value|Default|Description
---|---|---|---
`extractInfoLine.parsingScript`|String|"$DRAGEN_SCRIPTS_ROOT/bin/composeList.py"|Script for parsing inputs into a line
`extractInfoLine.timeout`|Int|4|Timeout for the job
`extractInfoLine.jobMemory`|Int|4|Job allocated RAM
`extractInfoLine.modules`|String|"dragen-scripts/0.3"|dependency modules
`composeList.listWritingScript`|String|"$DRAGEN_SCRIPTS_ROOT/bin/writeFile.py"|Script for writing out list of inputs
`composeList.jobMemory`|Int|4|Job allocated RAM
`composeList.timeout`|Int|4|Timeout for the job
`composeList.modules`|String|"dragen-scripts/0.3"|dependency modules
`runDragenRNAseq.additionalParameters`|String?|None|Additional dragen parameters
`runDragenRNAseq.correctGcBias`|Boolean|false|Correct for GC bias in fragment counts
`runDragenRNAseq.enableRNAfusion`|Boolean|true|Enable detection of transcript fusions
`runDragenRNAseq.restrictToProteinCodingGenes`|Boolean|true|A flag to restrict the analysis to protein-coding genes, only
`runDragenRNAseq.timeout`|Int|96|Hours before task timeout


### Outputs

Output | Type | Description | Labels
---|---|---|---
`splicingJunctions`|File|output splice junctions|


## Commands
 This section lists command(s) run by dragenRNAseqWORKFLOW workflow
 
 * dragenRNAseq
 
 This workflow generates results for the following analyses -
 
 * Gene fusiion analysis
 * Gene expression quantification
 * Alignment QC and other metrics
 
 ### Extract info for a single input file 
 
 ```
     python3 ~{parsingScript} -i ~{write_json(fastqInput)}
 ```
 
 ### Compose a list of inputs
 
 ```
    python3 ~{listWritingScript} -o ~{outputFileName} -l "~{sep=';' inputLines}"
 ```
 
 ### Run RNAseq pipeline on the input fastq files
 
 ```
    dragen -f -r ~{refDir} \
    --enable-rna true \
    --annotation-file ~{annotationGTF} \
    --enable-rna-quantification true \
    --rna-quantification-gc-bias ~{correctGcBias} \
    --enable-rna-gene-fusion ~{enableRNAfusion} \
    --rna-gf-restrict-genes ~{restrictToProteinCodingGenes} \
    --fastq-list ~{fastqList} \
    --output-directory . \
    --output-file-prefix ~{outputFileNamePrefix} ~{additionalParameters}
 ```
 ## Support

For support, please file an issue on the [Github project](https://github.com/oicr-gsi) or send an email to gsi@oicr.on.ca .

_Generated with generate-markdown-readme (https://github.com/oicr-gsi/gsi-wdl-tools/)_
