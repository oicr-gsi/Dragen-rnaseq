version 1.0

struct InputGroup {
  File fastqR1
  File fastqR2
  String readGroup
}

struct GenomeResources {
    String referenceDirectory
    String annotationGTF
    String dragenVersion
}

workflow dragenRNAseq {
    input {
        Array[InputGroup] fastqInputs
        String outputFileNamePrefix
        String reference
    }

    parameter_meta {
        fastqInputs: "Input structure with fastq files and read group strings"
        outputFileNamePrefix: "Prefix for output files"
        reference: "The genome reference build. For example: hg19, hg38, mm10"
    }

    Map[String,GenomeResources] dragen_resources_by_genome = { 
    "hg38": {
      "annotationGTF": "/.mounts/labs/gsiprojects/gsi/Dragen/reference/gencode/gencode.v31.annotation.gtf",
      "referenceDirectory": "/.mounts/labs/gsiprojects/gsi/Dragen/reference/hg38fa.p12/",
      "dragenVersion": "4.2.4"
    }}

    String dragen_ref = dragen_resources_by_genome [ reference ].referenceDirectory
    String dragen_gtf = dragen_resources_by_genome [ reference ].annotationGTF
    String dragen_version = dragen_resources_by_genome [ reference ].dragenVersion

    meta {
        author: "Peter Ruzanov"
        email: "pruzanov@oicr.on.ca"
        description: "A workflow for calling SNVs on tumor-only or tumor-normal inputs in somatic mode"
        dependencies: [
        {
          name: "gsi modules : gencode/31",
          url: "https://gitlab.oicr.on.ca/ResearchIT/modulator"
        },
        { 
          name: "gsi modules : dragen-scripts/0.3",
          url: "https://gitlab.oicr.on.ca/ResearchIT/modulator"
        }]
        output_meta: {
          unfilteredVcf: {
            description: "SNV calls before applying any filters",
            vidarr_label: "unfilteredVcf"
          }
        }
    }

    # Compose input lines using read group data
    scatter(t in fastqInputs) {
      call extractInfoLine {
        input:
        fastqInput = object{fastqR1: t.fastqR1, fastqR2: t.fastqR2, readGroup: t.readGroup}
      }
    }

    call composeList {
      input:
        inputLines = extractInfoLine.outputLine,
        outputFileName = "fastq_inputs.csv"
    } 
    
    call runDragenRNAseq {
      input:
        fastqList = composeList.inputList,
        refDir = dragen_ref,
        annotationGTF = dragen_gtf,
        outputFileNamePrefix = outputFileNamePrefix,
        dragenVersion = dragen_version
    } 

    output {
      File splicingJunctions = runDragenRNAseq.splicingJunctions
      File transcriptQuantification = runDragenRNAseq.transcriptQuantification
      File geneQuantification = runDragenRNAseq.geneQuantification
      File chimericJunctions = runDragenRNAseq.chimericJunctions
      File fusionCandidates = runDragenRNAseq.fusionCandidates
    }
}

# =====================================================================
# A scripted extraction of info from RG line to dragen-compliant string
# =====================================================================
task extractInfoLine {
   input {
       InputGroup fastqInput
       String parsingScript = "$DRAGEN_SCRIPTS_ROOT/bin/composeList.py"
       Int timeout = 4
       Int jobMemory = 4
       String modules = "dragen-scripts/0.3"
   }

   parameter_meta {
     fastqInput: "InputGroup struct entry with fastq files"
     parsingScript: "Script for parsing inputs into a line"
     timeout: "Timeout for the job"
     jobMemory: "Job allocated RAM"
     modules: "dependency modules"
   }

   command <<<
    python3 ~{parsingScript} -i ~{write_json(fastqInput)}
   >>>

   runtime {
     timeout: "~{timeout}"
     modules: "~{modules}"
     memory:  "~{jobMemory} GB"
   }

   output {
     String outputLine = read_string(stdout())
   }

   meta {
     output_meta: {
       outputLine: "Output line to use in a list of fastq files in dragen-compliant format"
     }
   }
}

# =====================================================================
#  Compose a dragen-compliant list of inputs to use with snv caller
# =====================================================================
task composeList {
   input  {
      Array[String] inputLines
      String listWritingScript = "$DRAGEN_SCRIPTS_ROOT/bin/writeFile.py"
      String outputFileName
      Int jobMemory = 4
      Int timeout = 4
      String modules = "dragen-scripts/0.3"
   }

   parameter_meta {
     inputLines: "Array of input lines to print"
     listWritingScript: "Script for writing out list of inputs"
     outputFileName: "Name of an output file, list of inputs"
     jobMemory: "Job allocated RAM"
     timeout: "Timeout for the job"
     modules: "dependency modules"
   }

   command<<<
   python3 ~{listWritingScript} -o ~{outputFileName} -l "~{sep=';' inputLines}"
   >>>
   

   runtime {
      timeout: "~{timeout}"
      modules: "~{modules}"
      memory:  "~{jobMemory} GB"
   }

   output {
     File inputList = "~{outputFileName}"
   }

   meta {
     output_meta: {
       inputList: "Output file to use with dragen SNV caller"
     }
   }
}
# ================================================================
# Main task for generating SNV calls in somatic mode (DRAGEN mode)
#
# we need CSV files with a header and data lines organized as:
#
# RGID Read Group
# RGSM Sample ID
# RGLB Library
# Lane Flow cell lane
# Read1File - Full path to a valid FASTQ input file
# Read2File - Full path to a valid FASTQ input file. Required for paired-end input. If not using paired-end input, leave empty.
# Each FASTQ file can only be referenced once in the CSV list.
# All values in the Read2File column must be reference valid files or must all be empty.
# ================================================================
task runDragenRNAseq {
    input {
        File fastqList
        String refDir
        String annotationGTF
        String dragenVersion
        String? additionalParameters
        String outputFileNamePrefix
        Boolean correctGcBias = false
        Boolean enableRNAfusion = true
        Boolean restrictToProteinCodingGenes = true         
        Int timeout = 96
    }

    parameter_meta {
        fastqList: "List of fastq files, required input"
        refDir: "The reference genome directoty"
        annotationGTF: "Annotation GTF file for transcript quantification"
        dragenVersion: "Expected version of dragen software on the DRAGEN node"
        additionalParameters: "Additional dragen parameters"
        outputFileNamePrefix: "Output file name prefix"
        correctGcBias: "Correct for GC bias in fragment counts"
        enableRNAfusion: "Enable detection of transcript fusions"
        restrictToProteinCodingGenes: "A flag to restrict the analysis to protein-coding genes, only"
        timeout: "Hours before task timeout"
    }
    
    String resultVcf = "~{outputFileNamePrefix}.vcf.gz"
    String hardfilteredVcfName = "~{outputFileNamePrefix}.hard-filtered.vcf.gz"
    String ploidyVcfName = "~{outputFileNamePrefix}.ploidy.vcf.gz"

    command <<<
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
    >>>
    runtime {
        backend: "DRAGEN"
        dragen_version: "~{dragenVersion}"
        timeout: "~{timeout}"
    }
    
    output {
        File splicingJunctions = "~{outputFileNamePrefix}.SJ.out.tab"
        File transcriptQuantification = "~{outputFileNamePrefix}.quant.sf"
        File geneQuantification = "~{outputFileNamePrefix}.quant.genes.sf"
        File chimericJunctions = "~{outputFileNamePrefix}.Chimeric.out.junction"
        File fusionCandidates = "~{outputFileNamePrefix}.fusion_candidates.features.csv"
    }

    meta {
        output_meta: {
            splicingJunctions: "output splice junctions",
            transcriptQuantification: "transcript expression quantification",
            geneQuantification: "gene-level expression quantification",
            chimericJunctions: "Predicted chimeric junction",
            fusionCandidates: "Fusion candidates, CSV file"
        }
    }

}

# ========================================================
# For outputs we may want to provision:
# HG008.quant.sf  - transcript expression quantification
# HG008.fusion_candidates.features.csv - fusion detection results
# HG008.Chimeric.out.junction - 
# HG008.SJ.out.tab
#
# There is a number of metrics files, we can parse them into a single .json file
# if required...
