nextflow.enable.dsl=2

process FASTQC {
  tag "$sample"
  publishDir "${params.outdir}/qc", mode: 'copy'
  input:
    tuple val(sample), path(read1), path(read2)
  output:
    path("*_fastqc.html"), emit: html
    path("*_fastqc.zip"), emit: zip
  script:
    """
    fastqc --threads 2 --outdir . $read1 $read2
    """
  stub:
    """
    touch ${sample}_1_fastqc.html ${sample}_2_fastqc.html ${sample}_1_fastqc.zip ${sample}_2_fastqc.zip
    """
}

process FASTP {
  tag "$sample"
  publishDir "${params.outdir}/trimmed", mode: 'copy'
  input:
    tuple val(sample), val(patient), val(condition), path(read1), path(read2)
  output:
    tuple val(sample), val(patient), val(condition),
      path("*_1.trimmed.fastq.gz"), path("*_2.trimmed.fastq.gz"),
      path("*_fastp.html"), path("*_fastp.json"), emit: trimmed
  script:
    """
    fastp --thread 2 --in1 $read1 --in2 $read2 --out1 ${sample}_1.trimmed.fastq.gz --out2 ${sample}_2.trimmed.fastq.gz --html ${sample}_fastp.html --json ${sample}_fastp.json
    """
  stub:
    """
    touch ${sample}_1.trimmed.fastq.gz ${sample}_2.trimmed.fastq.gz ${sample}_fastp.html ${sample}_fastp.json
    """
}

process STAR_INDEX {
  cpus 4
  input:
    path(fasta)
    path(gtf)
  output:
    path("star_index"), emit: index
  script:
    """
    STAR --runMode genomeGenerate --runThreadN 4 --genomeDir star_index --genomeFastaFiles $fasta --sjdbGTFfile $gtf --sjdbOverhang 100
    """
  stub:
    """
    mkdir star_index && touch star_index/STUB_INDEX.txt
    """
}

process STAR_ALIGN {
  tag "$sample"
  cpus 4
  publishDir "${params.outdir}/alignment", mode: 'copy'
  input:
    tuple val(sample), val(patient), val(condition), path(read1), path(read2), path(index)
  output:
    tuple val(sample), val(patient), val(condition), path("*Aligned.sortedByCoord.out.bam"), path("*Log.final.out"), emit: aligned
  script:
    """
    STAR --runThreadN 4 --genomeDir $index --readFilesIn $read1 $read2 --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --outFileNamePrefix ${sample}.
    """
  stub:
    """
    touch ${sample}.Aligned.sortedByCoord.out.bam ${sample}.Log.final.out
    """
}

process FEATURECOUNTS {
  cpus 4
  publishDir "${params.outdir}/counts", mode: 'copy'
  input:
    path(bams)
    path(gtf)
  output:
    path("counts.txt"), emit: counts
  script:
    """
    featureCounts -T 4 -a $gtf -o counts.txt -p -B -C $bams
    """
  stub:
    """
    printf 'Geneid\\tChr\\tStart\\tEnd\\tStrand\\tLength\\t%s\\n' "\\$(echo $bams | sed 's/ /\\t/g')" > counts.txt
    printf 'ENSG000001\\tchr1\\t1\\t100\\t+\\t100\\t10\\t12\\t8\\t14\\n' >> counts.txt
    """
}

process MULTIQC {
  publishDir "${params.outdir}/multiqc", mode: 'copy'
  input:
    path(reports)
  output:
    path("multiqc_report.html")
  script:
    """
    multiqc $reports --filename multiqc_report.html --force
    """
  stub:
    """
    echo '<html><body>Stub report only</body></html>' > multiqc_report.html
    """
}

process DESEQ2 {
  publishDir "${params.outdir}/deseq2", mode: 'copy'
  input:
    path(counts)
    path(samplesheet)
    path(rscript)
    val(control)
  output:
    path("deseq2_results")
  script:
    """
    Rscript $rscript $counts $samplesheet deseq2_results $control
    """
  stub:
    """
    mkdir -p deseq2_results && echo 'Stub only; no biological analysis.' > deseq2_results/VALIDATION_STATUS.txt
    """
}

workflow {
  if (!params.input || !params.fasta || !params.gtf) error "Provide --input, --fasta, and --gtf."
  samples = Channel.fromPath(params.input).splitCsv(header: true).map { row ->
    tuple(row.sample.toString(), row.patient.toString(), row.condition.toString(), file(row.fastq_1), file(row.fastq_2))
  }
  raw = Channel.fromPath(params.input).splitCsv(header: true).map { row ->
    tuple(row.sample.toString(), file(row.fastq_1), file(row.fastq_2))
  }
  FASTQC(raw)
  FASTP(samples)
  STAR_INDEX(file(params.fasta), file(params.gtf))
  align_inputs = FASTP.out.trimmed.map { sample, patient, condition, r1, r2, html, json -> tuple(sample, patient, condition, r1, r2) }.combine(STAR_INDEX.out.index)
  STAR_ALIGN(align_inputs)
  bams = STAR_ALIGN.out.aligned.map { sample, patient, condition, bam, log -> bam }.collect()
  FEATURECOUNTS(bams, file(params.gtf))
  reports = FASTQC.out.html.mix(FASTQC.out.zip)
    .mix(FASTP.out.trimmed.map { sample, patient, condition, r1, r2, html, json -> html })
    .mix(FASTP.out.trimmed.map { sample, patient, condition, r1, r2, html, json -> json })
    .collect()
  MULTIQC(reports)
  DESEQ2(FEATURECOUNTS.out.counts, file(params.input), file("$projectDir/scripts/run_deseq2.R"), params.control_condition)
}
