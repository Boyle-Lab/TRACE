#CAPER docker quay.io/ouyang/trace:version_1
#CAPER singularity docker://quay.io/ouyang/trace:version_1

workflow TRACE {
  #input {
    Boolean skipTrain = false
    Int THREAD = 40
    Int ITER = 200
    Int model_size = 10
    File seq_file = "./hg38.fa"
    File bam_file = ".bam"
    File bam_index_file = ".bami"
    File peak_file = ".bed"
    File peak_motif_file = "CTCF.bed"
    String genome = "hg38"
    String prefix = "test"
    Array[File] model_file_list = []
    Array[String] motif_list = []
  #}
  if (skipTrain) {
    if (length(motif_list) != length(model_file_list)) {
      call raise_exception as error_file_num  { input:
        msg = 'number of TFs and TRACE models provided are not same.'
      }
    }
  }
  scatter (i in range(length(motif_list))){
    call data_process { input:
      motif = motif_list[i],
      model_size = model_size,
      bam_file = bam_file,
      bam_index_file = bam_index_file,
      seq_file = seq_file,
      peak_file = peak_file,
      prefix = prefix,
      genome = genome
    }
    if (skipTrain) {
      call Trace_viterbi {input:
        THREAD = THREAD,
        motif = motif_list[i],
        prefix = prefix,
        model_file = model_file_list[i],
        count_file = data_process.count_output,
        seq_file = data_process.seq_output,
        slope_file = data_process.slope_2_output,
        peak_motif_file = peak_motif_file,
        peak_file = peak_file
      }
    }
    if (!skipTrain)  {
      call Trace_train { input:
        THREAD = THREAD,
        ITER = ITER,
        motif = motif_list[i],
        prefix = prefix,
        model_file = data_process.model_file,
        count_file = data_process.count_output,
        seq_file = data_process.seq_output,
        slope_file = data_process.slope_2_output,
        peak_motif_file = peak_motif_file,
        peak_file = peak_file
      }
    }
  }
  output {
    Array[File] init_data = data_process.model_file
    Array[Array[File]?] footprints_file = Trace_train.out
    Array[Array[File]?] footprints_v_file = Trace_viterbi.out
  }
}


task data_process {
  #input {
    Int model_size
    String genome
    String motif
    String prefix
    File bam_file
    File bam_index_file
    File seq_file
    File peak_file
  #}
  command {
    source activate TRACE_env
    /footprinting/TRACE/scripts/init_hmm.py ${motif}  --prefix ${prefix}_ --motif-number ${model_size}
    sort -k1,1 -k2,2n ${peak_file} | uniq > test_peak_3.bed
    /footprinting/TRACE/scripts/dataProcessing.py test_peak_3.bed ${bam_file} ${seq_file} --genome ${genome} --prefix ${prefix}
    rm test_peak_3.bed
  }
  output {
    File model_file = "${prefix}_${motif}_init_model.txt"
    File count_output = glob("*_count.txt")[0]
    File seq_output = glob("*_seq.txt")[0]
    File slope_2_output = glob("*_slope_2.txt")[0]
  }
}

task Trace_train {
  #input {
    Int THREAD
    Int ITER
    String motif
    String prefix
    File model_file
    File count_file
    File seq_file
    File slope_file
    File peak_file
    File peak_motif_file
  #}
  command {
    source activate TRACE_env
    sort -k1,1 -k2,2n ${peak_file} | uniq > test_peak_3.bed
    sort -k1,1 -k2,2n ${peak_motif_file} > test_peak_7.bed
    /footprinting/TRACE/TRACE ${seq_file} ${count_file} ${slope_file} --initial-model ${model_file} --final-model ${motif}_hmm.txt --motif-file test_peak_7.bed --peak-file test_peak_3.bed --thread THREAD --max-inter ITER
    rm test_peak_3.bed
    rm test_peak_7.bed
  }
  output {
    Array[File] out = glob("*.txt")
  }
}

task Trace_viterbi {
  #input {
    Int THREAD
    String motif
    String prefix
    File model_file
    File count_file
    File seq_file
    File slope_file
    File peak_file
    File peak_motif_file
  #}
  command {
    source activate TRACE_env
    sort -k1,1 -k2,2n ${peak_file} | uniq > test_peak_3.bed
    sort -k1,1 -k2,2n ${peak_motif_file} > test_peak_7.bed
    /footprinting/TRACE/TRACE --viterbi ${seq_file} ${count_file} ${slope_file} --final-model ${model_file} --motif-file test_peak_7.bed --peak-file test_peak_3.bed --thread THREAD
    rm test_peak_3.bed
    rm test_peak_7.bed
  }
  output {
    Array[File] out = glob("*.txt")
  }
}

task raise_exception {
  String msg
  command {
    echo -e "\n* Error: ${msg}\n" >&2
    exit 2
  }
  output {
    String error_msg = '${msg}'
  }
  runtime {
    maxRetries : 0
  }
}
