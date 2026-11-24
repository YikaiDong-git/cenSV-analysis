import os, sys, time
import os as _censv_os
WORKSPACE_BUCKET = _censv_os.environ["WORKSPACE_BUCKET"]
os.environ["JAVA_HOME"] = "/usr/lib/jvm/temurin-11-jdk-amd64"
import hail as hl

CHROM = sys.argv[1]
CONN = "/usr/local/share/google/dataproc/lib/gcs-connector.jar"
VDS = "gs://vwb-aou-datasets-controlled/v8/wgs/short_read/snpindel/vds/hail.vds"
LOCI = "/tmp/vdsx/predictor_loci_all.tsv"
OUTB = (WORKSPACE_BUCKET + "/ydon268/Project1_Centromere/analyses/"
        "disease_replication/aou_replication/results_vds/vds_predictor_snps_full")
OUT = f"{OUTB}/{CHROM}/{CHROM}.predictor_snps"
RES = f"/tmp/vdsx/ex2_{CHROM}.txt"
def log(m):
    open(RES, "a").write(m + "\n"); print(m, flush=True)
open(RES, "w").close(); t0 = time.time()


def get_call_expr(mt):
    ef = set(mt.entry.dtype.fields)
    return mt.GT if "GT" in ef else mt.LGT


def biallelic_snp(rows):
    return (hl.len(rows.alleles) == 2) & (hl.len(rows.alleles[0]) == 1) & (hl.len(rows.alleles[1]) == 1)


pos = sorted({int(l.split("\t")[1]) for l in open(LOCI) if l.split("\t")[0] == CHROM})
iv = []; s = pos[0]; prev = pos[0]
for p in pos[1:]:
    if p - prev > 1_000_000:
        iv.append((max(1, s - 1000), prev + 1000)); s = p
    prev = p
iv.append((max(1, s - 1000), prev + 1000))
log(f"CHROM {CHROM}: {len(pos)} loci -> {len(iv)} intervals ({pos[0]/1e6:.1f}-{pos[-1]/1e6:.1f}Mb)")

try:
    hl.init(master="local[16]", default_reference="GRCh38", quiet=True,
            log=f"/tmp/vdsx/hail2_{CHROM}.log", tmp_dir="/tmp/vdsx/tmp",
            spark_conf={
                "spark.jars": CONN, "spark.driver.extraClassPath": CONN,
                "spark.hadoop.fs.gs.impl": "com.google.cloud.hadoop.fs.gcs.GoogleHadoopFileSystem",
                "spark.hadoop.fs.AbstractFileSystem.gs.impl": "com.google.cloud.hadoop.fs.gcs.GoogleHadoopFS",
                "spark.hadoop.google.cloud.auth.service.account.enable": "true",
                "spark.hadoop.fs.gs.requester.pays.mode": "CUSTOM",
                "spark.hadoop.fs.gs.requester.pays.project.id": "wb-quick-beet-1004",
                "spark.hadoop.fs.gs.requester.pays.buckets": "vwb-aou-datasets-controlled",
                "spark.driver.memory": "36g",
            })
    log("INIT_OK")
    intervals = [hl.parse_locus_interval(f"{CHROM}:{a}-{b}", reference_genome="GRCh38") for a, b in iv]
    vds = hl.vds.read_vds(VDS)
    vds = hl.vds.filter_intervals(vds, intervals, split_reference_blocks=False)
    vds = hl.vds.split_multi(vds, filter_changed_loci=True)
    posset = hl.literal(set(pos))
    rows = vds.variant_data
    rows = rows.filter_rows(biallelic_snp(rows) & posset.contains(rows.locus.position))
    variant_keys = rows.rows().key_by("locus", "alleles").select()
    vds = hl.vds.filter_variants(vds, variant_keys)
    mt = hl.vds.to_dense_mt(vds)
    mt = mt.annotate_rows(vid=mt.locus.contig + ":" + hl.str(mt.locus.position) + ":"
                          + mt.alleles[0] + ":" + mt.alleles[1])
    ce = get_call_expr(mt)
    ce = hl.bind(lambda gt: hl.if_else(hl.is_missing(gt), hl.missing(hl.tcall),
                 hl.if_else(gt.ploidy == 1, hl.call(gt[0], gt[0]),
                 hl.if_else(gt.ploidy == 2, gt, hl.missing(hl.tcall)))), ce)
    hl.export_plink(mt, OUT, call=ce, fam_id=mt.s, ind_id=mt.s, varid=mt.vid)
    log(f"EXPORT_OK -> {OUT}  ({time.time()-t0:.0f}s)")
    log("DONE_OK")
except Exception as e:
    log(f"ERROR {type(e).__name__}: {str(e)[:500]}"); raise
