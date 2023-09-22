version 1.0

import "https://raw.githubusercontent.com/PRIMED-Ancestry-Inference/PCA_projection/optional_steps/create_pca_projection.wdl" as tasks

workflow basic_admixture {
  input {
    File bed
    File bim
    File fam
    File? ref_pop
    Int n_ancestral_populations
    Boolean cross_validation = false
    Boolean remove_relateds = true
    Float? max_kinship_coefficient
    Boolean prune_variants = true
		Int? window_size
		Int? shift_size
		Int? r2_threshold
  }

  if (remove_relateds) {
    call tasks.removeRelateds {
      input:
        bed = bed,
        bim = bim,
        fam = fam,
        max_kinship_coefficient = max_kinship_coefficient
    }
  }

  if (prune_variants) {
    call tasks.pruneVars {
      input:
        bed = select_first([removeRelateds.out_bed, bed]),
        bim = select_first([removeRelateds.out_bim, fam]),
        fam = select_first([removeRelateds.out_fam, fam]),
        window_size = window_size,
        shift_size = shift_size,
        r2_threshold = r2_threshold
    }
  }

  # if(defined(thin_count)) {
  #   call ThinVariants {
  #     input:
  #       bed = bed,
  #       bim = bim,
  #       fam = fam,
  #       thin_count = select_first([thin_count])
  #   }
  # }

  call Admixture_t {
    input:
      #bed = select_first([ThinVariants.thinned_bed, bed]),
      #bim = select_first([ThinVariants.thinned_bim, bim]),
      #fam = select_first([ThinVariants.thinned_fam, fam]),
      bed = select_first([pruneVars.out_bed, removeRelateds.out_bed, bed]),
      bim = select_first([pruneVars.out_bim, removeRelateds.out_bim, bim]),
      fam = select_first([pruneVars.out_fam, removeRelateds.out_fam, fam]),
      n_ancestral_populations = n_ancestral_populations
  }

  output {
    File ancestry_fractions = Admixture_t.ancestry_fractions
    File allele_frequencies = Admixture_t.allele_frequencies
  }
}

task Admixture_t {
  input {
    File bed
    File bim
    File fam
    File? ref_pop
    Int n_ancestral_populations
    Boolean cross_validation = false
    Int mem = 16
    Int n_cpus = 4
  }

  Int disk_size = ceil(1.5*(size(bed, "GB") + size(bim, "GB") + size(fam, "GB")))
  String basename = basename(bed, ".bed")

  command <<<
    /admixture_linux-1.3.0/admixture ~{if (cross_validation) then "--cv" else ""} \
      ~{bed} ~{n_ancestral_populations} ~{if defined(ref_pop) then "--supervised ~{ref_pop}" else ""} \
      -j~{n_cpus}
  >>>

  runtime {
    docker: "us.gcr.io/broad-dsde-methods/admixture_docker:v1.0.0"
    disks: "local-disk " + disk_size + " HDD"
    memory: mem + " GB"
    cpu: n_cpus
  }

  output {
    File ancestry_fractions = "~{basename}.~{n_ancestral_populations}.Q"
    File allele_frequencies = "~{basename}.~{n_ancestral_populations}.P"
  }
}

task ThinVariants {
  input {
    File bed
    File bim
    File fam
    Int thin_count
    Int mem = 8
  }

  Int disk_size = ceil(1.5*(size(bed, "GB") + size(bim, "GB") + size(fam, "GB")))
  String basename = basename(bed, ".bed")

  command <<<
    ln -s ~{bim} input.bim
    ln -s ~{bed} input.bed
    ln -s ~{fam} input.fam

    /plink2 --bfile input \
    --thin-count ~{thin_count} \
    --make-bed \
    --out ~{basename}.thinned

  >>>

  output {
    File thinned_bed = "~{basename}.thinned.bed"
    File thinned_bim = "~{basename}.thinned.bim"
    File thinned_fam = "~{basename}.thinned.fam"
  }

  runtime {
    docker: "us.gcr.io/broad-dsde-methods/plink2_docker@sha256:4455bf22ada6769ef00ed0509b278130ed98b6172c91de69b5bc2045a60de124"
    disks: "local-disk " + disk_size + " HDD"
    memory: mem + " GB"
  }
}
