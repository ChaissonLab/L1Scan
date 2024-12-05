# Example workflow
# Declare WDL version 1.0 if working in Terra
version 1.0

workflow l1Scan {

    input {
        File BAM
        File BAM_INDEX
        File L1HS
	File REFERENCE
        String SAMPLE = basename(BAM, ".bam")
        Int cpu
        Int mem
        Int diskSizeGb
        Int maxCoverage
    }

    call l1scanAnnotation {
        input:
        bam = BAM,
        bam_index = BAM_INDEX,
	reference = REFERENCE,
        l1File = L1HS,
        sample = SAMPLE,
	taskCpu = cpu,
	taskMem = mem,
	taskDiskSizeGb = diskSizeGb,
	taskMaxCoverage = maxCoverage
    }

    output {
        File l1hsOut = l1scanAnnotation.outTab
    }
}

task l1scan {
    input {
        File bam
        File bam_index
	File reference
	File l1hs
        String sample
        Int taskCpu
        Int taskMem
        Int taskDiskSizeGb
        Int taskMaxCoverage
    }

    command <<<
        l1scan ~{bam} ~{reference} ~{l1hs} > ~{sample}.tab ; 
        gzip ~{sample}.tab
    >>>

    output {
        File outTab = "~{sample}.tab.gz"
    }

    runtime {
        docker: "mchaisso/l1scan:v1"
        cpu: taskCpu
        memory: taskMem+"GB"
        disks: "local-disk " + taskDiskSizeGb + " LOCAL"
    }
}