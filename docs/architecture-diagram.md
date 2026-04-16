# DeFrABB Pipeline Architecture Diagram

_Prepared: 2026-03-24_

## High-Level Pipeline Flow

```mermaid
flowchart TD
    subgraph CONFIG["Configuration Layer"]
        AT["analyses*.tsv\n(per-row evaluation definitions)"]
        RY["resources.yml\n(references, assemblies, exclusions,\ncomparisons, compute settings)"]
        SC["schema/*.yml\n(JSON Schema validation)"]
        AT --> CM
        RY --> CM
        SC --> CM
    end

    subgraph COMMON["rules/common.smk — Config Parsing & Helpers"]
        CM["load_df() / load_analyses()\nschema validation"]
        CM --> VCT["vc_tbl\n(variant calling jobs)"]
        CM --> BT["bench_tbl\n(benchmark jobs)"]
        CM --> AN["analyses\n(evaluation jobs)"]
        VCT --> DCT["dipcall_tbl"]
        AN --> HA["happy_analyses"]
        AN --> TA["truvari_analyses"]
        AN --> TRA["truvari_refine_analyses"]
        BT --> BET["bench_excluded_tbl"]
    end

    subgraph DOWNLOAD["rules/download_resources.smk"]
        GR["get_ref\n(reference FASTA)"]
        GA["get_assemblies\n(maternal/paternal FASTAs)"]
        GS["get_strats\n(stratification tarballs)"]
        GCV["get_comparison_vcf/bed\n(truth/query callsets)"]
        DBG["download_bed_gz\n(exclusion BEDs)"]
    end

    subgraph INDEX["rules/utils.smk — Indexing"]
        IR["index_ref (.fai)"]
        BWA["bwa_index (.amb/.ann/.bwt/.pac/.sa)"]
        MMI["index_ref_mmi (.mmi)"]
        SDF["index_ref_sdf (.sdf/)"]
        TAB["tabix (.tbi)"]
        IDB["index_dip_bam (.bai)"]
    end

    subgraph VARCALL["rules/asm-varcall.smk — Variant Calling"]
        DIP["run_dipcall\n(dipcall 0.3)\nhap1/hap2 BAMs + VCF"]
        PAV["run_pav\n(PAV container)\nVCF + callable regions"]
        REN["rename_dipcall_vcf_sample"]
        IPC["intersect_pav_callable_regions\n(diploid callable BED)"]
        STD["standardize_vcasm_output\n(unified VCF + baseline BED)"]
        DIP --> REN --> STD
        PAV --> IPC --> STD
    end

    subgraph VCFPROC["rules/bench_vcf_processing.smk — VCF Processing Chain"]
        COPY["copy_std_asm_vcf_to_annotations"]
        FIX["fix_XY_genotype"]
        SPLIT["split_multiallelic_sites"]
        NORM["normalize_vars / filter_lt19_and_norm"]
        END_["add_end_info_header"]
        TRF["run_truvari_anno_trf"]
        SVI["run_truvari_anno_svinfo"]
        RPM["run_truvari_anno_repmask"]
        RMP["run_truvari_anno_remap"]
        LCR["run_truvari_anno_lcr"]
        INV["remove_inv"]
        MOV["move_processed_draft_bench_vcf"]
        COPY --> FIX & SPLIT & NORM & END_
        FIX & SPLIT & NORM & END_ --> TRF & SVI & RPM & RMP & LCR
        TRF & SVI & RPM & RMP & LCR --> INV --> MOV
    end

    subgraph EXCL["rules/exclusions.smk — Exclusion Processing"]
        GAPS["make_gaps_bed\n(reference Ns)"]
        SLOP["add_slop\n(15kb buffer)"]
        SLPM["add_slop_and_merge\n(15kb slop + 10kb merge)"]
        SVW["get_sv_widen_coords\n(SV regions from VCF)"]
        ISRS["intersect_SVs_and_simple_repeats"]
        CSV["get_consecutive_svs\n(consecutive indels in BAM)"]
        SD["self_discrep_happy\n(self-comparison FP/FN)"]
        SDE["self_discrep_extract_fpfns"]
        SDI["self_discrep_intersect_slop"]
        PINV["exclude_pav_inversions"]
        ISE["intersect_start_and_end\n(assembly breaks)"]
        FLK["get_flanks\n(15kb flanks)"]
        SUB["subtract_exclusions\n(final benchmark BED)"]
        GIS["generate_intersection_summary"]
        GAPS & SLOP & SLPM --> SUB
        SVW --> ISRS --> SUB
        CSV --> SUB
        SD --> SDE --> SDI --> SUB
        PINV --> SUB
        ISE --> SUB
        FLK --> SUB
        SUB --> GIS
    end

    subgraph BENCH["Benchmark Generation"]
        GVBR["get_variants_in_benchmark_regions\n(filter VCF to benchmark BED)"]
    end

    subgraph EVAL["rules/evaluation.smk"]
        RI["region_intersection\n(benchmark ∩ comparison BEDs)"]
        HPY["run_happy\n(hap.py / vcfeval)\nsummary.csv, extended.csv"]
        TRU["run_truvari\n(truvari bench)\nsummary.json, TP/FP/FN VCFs"]
        TRR["truvari_refine\n(mafft realignment)\nrefine.variant_summary.json"]
        RI --> HPY & TRU
        TRU --> TRR
    end

    subgraph REPORT["rules/report.smk"]
        AST["run_assembly_stats"]
        BCS["get_bcftools_stats"]
        RTS["get_rtg_vcf_stats"]
        BDS["get_bed_stats"]
        WRP["write_report_params\n(analysis_params.yml)"]
        RND["render_report\n(Quarto → analysis.html)"]
        AST & BCS & RTS & BDS --> WRP --> RND
    end

    subgraph RELEASE["run_defrabb wrapper"]
        PIPE["execute_snakemake_pipeline()"]
        SRPT["generate_snakemake_report()"]
        ARCH["generate_snakemake_archive()"]
        REL["release_run()\n(S3 / local NAS)"]
        MAN["create_data_manifest()"]
        PIPE --> SRPT --> ARCH --> REL --> MAN
    end

    %% Main flow connections
    CONFIG --> DOWNLOAD
    DOWNLOAD --> INDEX
    INDEX --> VARCALL
    VARCALL --> VCFPROC
    VCFPROC --> EXCL
    VCFPROC --> BENCH
    EXCL --> BENCH
    DOWNLOAD --> EVAL
    BENCH --> EVAL
    EVAL --> REPORT
    REPORT --> RELEASE
```

## Configuration Data Model

```mermaid
erDiagram
    ANALYSES_TSV {
        string eval_id PK "unique evaluation ID"
        string bench_id "benchmark set ID"
        string eval_cmd "happy | truvari | unhappy"
        string eval_params "additional eval args"
        string eval_comp_id "comparison callset ID"
        boolean eval_comp_id_is_truth "swap truth/query roles"
        boolean eval_truth_regions "use truth BED in eval"
        boolean eval_target_regions "use query BED in eval"
        string vc_id "variant caller config ID"
        string bench_type "smvar | stvar"
        string bench_vcf_processing "dot-separated processing steps"
        string bench_bed_processing "none | exclude"
        string exclusion_set "named exclusion preset or none"
        string asm_id "assembly ID"
        string ref "reference ID"
        string vc_cmd "dipcall | pav"
        string vc_param_id "caller parameter set"
        string vc_params "parameter values"
    }

    RESOURCES_YML {
        map references "ref_id -> ref_url, par_bed, adotto_db, exclusions, strats"
        map assemblies "asm_id -> maternal, paternal, is_male, sample_id"
        map comparisons "ref_id -> comp_id -> vcf_url, bed_url, tbi_url"
        map exclusion_set "set_name -> [region_names]"
        map exclusion_slop_regions "set_name -> [regions for 15kb slop]"
        map exclusion_slopmerge_regions "set_name -> [regions for slop+merge]"
        map exclusion_asm_intersect "set_name -> [regions for asm breaks]"
        map exclusion_asm_agnostic "set_name -> [assembly-independent regions]"
        map exclusion_ref_agnostic "set_name -> [reference-independent regions]"
        int _dipcall_threads "5"
        int _dipcall_jobs "4"
        int _dipcall_mem "32000"
        int _pav_threads "24"
        int _happy_threads "12"
        int _happy_mem "64000"
        int _truvari_refine_threads "24"
    }

    ANALYSES_TSV ||--o{ RESOURCES_YML : "references by asm_id, ref, eval_comp_id"
```

## Conda Environment Dependency Map

```mermaid
graph LR
    subgraph ENVS["Conda Environments (envs/*.yml)"]
        dipcall["dipcall.yml\ndipcall==0.3"]
        bedtools["bedtools.yml\nbedtools=2.31.0\npybedtools=0.8.2"]
        bcftools["bcftools.yml\nbcftools==1.14"]
        bcf_bed["bcftools_and_bedtools.yml\nbcftools==1.17\nbedtools=2.30.0"]
        happy["happy.yml\nhap.py==0.3.14\nrtg-tools==3.12.1"]
        rtgtools["rtgtools.yml\nrtg-tools==3.10.1"]
        truvari["truvari.yml\ntruvari==5.3.0\ntrf==4.09.1"]
        truv_remap["truvari_remap.yml\nTruvari==4.2.2 (pip)\ntrf==4.09.1\nmafft, repeatmasker"]
        sv_widen["sv_widen_coords.yml\npython=3.8\npybedtools"]
        consec["consecutive_svs.yml\npybedtools=0.10.0\npysam=0.22.1"]
        seqtk["seqtk.yml\nseqtk==1.4"]
        download["download_remotes.yml\ncurl, htslib"]
        quarto["quarto.yml\nquarto, r-tidyverse"]
    end

    subgraph ISSUES["Version Inconsistencies"]
        I1["bcftools: 1.14 vs 1.17"]
        I2["rtg-tools: 3.10.1 vs 3.12.1"]
        I3["truvari: 5.3.0 vs 4.2.2"]
        I4["python: 3.8 vs 3.12+"]
        I5["pybedtools: 0.8.2 vs 0.10.0"]
    end

    bcftools -.->|mismatch| bcf_bed
    happy -.->|mismatch| rtgtools
    truvari -.->|mismatch| truv_remap

    style I1 fill:#fee,stroke:#c00
    style I2 fill:#fee,stroke:#c00
    style I3 fill:#fee,stroke:#c00
    style I4 fill:#fee,stroke:#c00
    style I5 fill:#fee,stroke:#c00
```

## Output Directory Structure

```text
results/
+-- asm_varcalls/{vc_id}/
|   +-- {ref}_{asm}_{vc_cmd}-{vc_param}.dip.vcf.gz[.tbi]      # dipcall output
|   +-- {ref}_{asm}_{vc_cmd}-{vc_param}.hap[12].bam[.bai]      # dipcall alignments
|   +-- {ref}_{asm}_{vc_cmd}-{vc_param}.baseline.bed            # callable regions
|   +-- annotations/                                             # processed VCFs
|       +-- {ref}_{asm}_{vc_cmd}-{vc_param}.{steps}.vcf[.gz]
|
+-- draft_benchmarksets/{bench_id}/
|   +-- {ref}_{asm}_{bench_type}_{vc_cmd}-{vc_param}.vcf.gz     # benchmark VCF
|   +-- {ref}_{asm}_{bench_type}_{vc_cmd}-{vc_param}.benchmark.bed
|   +-- {ref}_{asm}_{bench_type}_{vc_cmd}-{vc_param}.exclusion_stats.txt
|   +-- exclusions/                                              # per-exclusion BEDs
|
+-- evaluations/
|   +-- happy/{eval_id}_{bench_id}/                              # hap.py outputs
|   |   +-- {ref}_{comp}_{asm}_{bench_type}_{vc_cmd}-{vc_param}.summary.csv
|   |   +-- {ref}_{comp}_{asm}_{bench_type}_{vc_cmd}-{vc_param}.extended.csv
|   +-- truvari/{eval_id}_{bench_id}/                            # truvari outputs
|       +-- {ref}_{comp}_{asm}_{bench_type}_{vc_cmd}-{vc_param}/
|           +-- summary.json, tp-*.vcf.gz, fp.vcf.gz, fn.vcf.gz
|
+-- report/
|   +-- assemblies/{asm_id}_{haplotype}_stats.txt
+-- analysis_params.yml
+-- analysis.html
```
