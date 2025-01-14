let doc = """
Pileup from VCF using hileup and determine genotype probabilities based on read base quality.
Usage:
  annotate_w_pileup [options] <BAM> <FASTA> <VCF> <OUT>

Arguments:
  <BAM>     Input bam file
  <FASTA>   Reference fasta file
  <VCF>     Input vcf file
  <OUT>     Output file

Options:
  -h --help                 Show this screen.
  --mapq <mapq>             Minimum mapping quality [default: 20].
  --baseq <baseq>           Minimum base quality [default: 20].
  --depth <depth>           Minimum depth [default: 20].
  --snvs <snvs>             VCF of SNVs or Indels [default: true].
"""

import hts
import hts/private/hts_concat
import hile
import times
import strutils
import sequtils
import progress
import math
import docopt
import strformat
import tables
import sets
import terminal
import algorithm

let args = docopt(doc, version = "0.1.0")

const max_reads = 300
const default_mapq = 0

let
    time = now()
    bam_path = $args["<BAM>"]
    vcf_path = $args["<VCF>"]
    fai_path = $args["<FASTA>"]
    depth_thresh = parseInt($args["--depth"])
    snvs = parseBool($args["--snvs"])
    min_mapq = uint8(parseUInt($args["--mapq"]))
    min_baseq = uint8(parseUInt($args["--baseq"]))
var
    bam: Bam
    fai: Fai
    vcf: VCF
    output = open($args["<OUT>"], fmWrite)

setForegroundColor(fgCyan)
stdout.styledWriteLine(fgCyan,"Input Bam: " & $args["<BAM>"])
stdout.styledWriteLine(fgCyan,"Output: " & $args["<OUT>"])
stdout.styledWriteLine(fgCyan,"VCF: " & $args["<VCF>"])
stdout.styledWriteLine(fgCyan,"FASTA: " & $args["<FASTA>"])

if not open(bam, bam_path, index=true, fai=fai_path):
    quit "bad bam"

if not open(fai, fai_path):
    quit "bad fai"

discard open(vcf, vcf_path)
var total = toSeq(vcf.items()).len
close(vcf)

var cfg = Config(TrackBaseQualities: true, MinBaseQuality: uint8(parseUInt($args["--baseq"])), TrackMappingQualities: true, MinMappingQuality: uint8(parseUInt($args["--mapq"])), TrackReadNames: true, ExcludeFlags: BAM_FUNMAP or BAM_FSECONDARY or BAM_FQCFAIL or BAM_FDUP)

var bar = newProgressBar(total=total)
bar.start()
discard open(vcf, vcf_path)

proc phred_to_p(x: int): float = 
    return pow(10, -1*float(x)/10)

proc strand_info_snv(h: Hile, ref_base: string, alt_base: string) : (int, int, int, int) =
    var
        rf = 0
        af = 0
        rr = 0
        ar = 0
    for b in h.bases:
        case $b.reverse_strand
        of "0":
            if $b.base.char == $ref_base:
                rf += 1
            elif $b.base.char == $alt_base:
                af += 1
        of "1":
            if $b.base.char == $ref_base:
                rr += 1
            elif $b.base.char == $alt_base:
                ar += 1
        else: continue
           
    return (rf, af, rr, ar)

proc strand_info_indel(h: Hile, var_reads: HashSet[string]) : (int, int, int, int) =
    var
        rf = 0
        af = 0
        rr = 0
        ar = 0
    for pairs in zip(h.read_names, h.bases):
        let rname = pairs[0]
        let reverse_strand = bool(pairs[1].reverse_strand)
        if var_reads.contains(rname):
            if reverse_strand:
                ar += 1
            else:
                af += 1
        else:
            if reverse_strand:
                rr += 1
            else:
                rf += 1 
           
    return (rf, af, rr, ar)

proc get_avg_alt_qualities(h: Hile, var_reads: HashSet[string]): (float, float) =
    var 
        avgAltBq: float = 0
        avgAltMq: float = 0
        alt_bqs = newSeq[int](0)
        alt_mqs = newSeq[int](0)
        bqs = map(h.bqs, proc(x: uint8): int = int(x))
        mqs = map(h.mqs, proc(x: uint8): int = int(x))
    for idx, rname in h.read_names:
        let bq = bqs[idx]
        let mq = mqs[idx]
        if var_reads.contains($rname):
            alt_bqs.add(bq)
            alt_mqs.add(mq)
    if len(alt_bqs) > 0:
        avgAltBq = sum(alt_bqs) / int(alt_bqs.len)
    if len(alt_mqs) > 0:
        avgAltMq = sum(alt_mqs) / int(alt_mqs.len)
    return (avgAltBq, avgAltMq)

proc get_avg_alt_mate_mapq(bam: Bam, chrom: string, position: int64, var_reads: HashSet[string]): float =
    type
        read_info = tuple
            rname: string
            chrom: string
            pos: int64
    var
        seen = initHashSet[string]()
        to_search = newSeq[read_info]()
        mate_mapqs = newSeq[float](0)
        avgAltMateMq: float = 0

    for r in bam.query(chrom, int(position-1), int(position)):
        if not r.flag.qcfail and not r.flag.secondary and not r.flag.supplementary and not r.flag.dup:
            if not seen.contains(r.qname) and var_reads.contains(r.qname):
                seen.incl(r.qname)
                to_search.add((r.qname, r.mate_chrom, r.mate_pos))
    # now search through all mates
    for entry in to_search:
        var
            mqname = entry[0]
            mchrom = entry[1]
            mpos = entry[2]
            reads = 0
        for r in bam.query(mchrom, int(mpos-1), int(mpos)):
            reads+=1
            if reads > max_reads:
                mate_mapqs.add(default_mapq)
                break
            if not r.flag.qcfail and not r.flag.secondary and not r.flag.supplementary and not r.flag.dup:
                if r.qname == mqname:
                    mate_mapqs.add(float(r.mapping_quality))
                    break
    
    if len(mate_mapqs) > 0:
        avgAltMateMq = sum(mate_mapqs) / float(len(mate_mapqs))
    return avgAltMateMq

proc median(xs: seq[float]): float =
  var ys = xs
  sort(ys, system.cmp[float])
  return 0.5 * (ys[ys.high div 2] + ys[ys.len div 2])

proc get_avg_insert_sizes(bam: Bam, chrom: string, position: int64, var_reads: HashSet[string]): (float, float) =
    var
        seen = initHashSet[string]()
        isizes = newSeq[float](0)
        alt_isizes = newSeq[float](0)

    for r in bam.query(chrom, int(position-1), int(position)):
        if not r.flag.qcfail and not r.flag.secondary and not r.flag.supplementary and not r.flag.dup:
            if not seen.contains(r.qname):
                seen.incl(r.qname)
                isizes.add(float(abs(r.isize)))
                if var_reads.contains(r.qname):
                    alt_isizes.add(float(abs(r.isize)))
    var 
        avgIs = median(isizes)
        avgAltIs = median(alt_isizes)
    return (avgIs, avgAltIs)

proc get_alt_alleles(baseinfo: seq, ref_base: string): int =
    var
        bases: array[4, string] = ["A", "C", "G", "T"]
        tmpVaf: float
        num_alt = 0
    for alt in bases:
        if alt != ref_base:
            tmpVaf = (count(map(baseinfo, proc(x: basestrand): string = $(x.base.char)), alt)) / (len(baseinfo))
            if tmpVaf >= 0.02:
                num_alt += 1
    return num_alt

proc get_indel_pileup(bam: Bam, chrom: string, position: int64, fai: Fai, cfg: Config, REF: string, ALT: string): (Hile, string) =
    if len(REF) > len(ALT):
        let h = bam.hileup(chrom, int(position-1), fai, cfg, snv=false)
        return (h, "DEL")
    else:
        let h = bam.hileup(chrom, int(position-1), fai, cfg, snv=false)
        return (h, "INS")

proc get_var_reads(h: Hile, ALT: string, var_type: string): HashSet[string] =
    var
        var_reads = initHashSet[string]()
    if var_type == "INS":
        for ins in h.ins:
            var_reads.incl(h.read_names[ins.index])
    elif var_type == "DEL":
        for del in h.del:
            var_reads.incl(h.read_names[del.index])
    else:
        for pairs in zip(h.read_names, h.bases):
            if $pairs[1].base.char == ALT:
                var_reads.incl(pairs[0])

    return var_reads

proc get_indel_count(bam: Bam, chrom: string, position: int64, alt_count: int): (float) =
    var
        del_ct : int = 0
        ins_ct: int = 0
        indel_count: float
    for r in bam.query(chrom, int(position-1), int(position)):
        if not r.flag.qcfail and not r.flag.secondary and not r.flag.supplementary and not r.flag.dup:
            if r.mapping_quality >= min_mapq:
                for op in r.cigar:
                    if op.op == CigarOp.deletion:
                        del_ct += 1
                    elif op.op == CigarOp.insert:
                        ins_ct += 1
    if del_ct+ins_ct > 0:
        indel_count = alt_count/(del_ct+ins_ct)
    else:
        indel_count = 0
    return indel_count

proc get_avg_edit_dist(bam: Bam, chrom: string, position: int64, var_reads: HashSet[string]): (float) = 
    var
        avgEditDist: float = 0
        edit_dists = newSeq[int](0) 
    for r in bam.query(chrom, int(position-1), int(position)):
        if not r.flag.qcfail and not r.flag.secondary and not r.flag.supplementary and not r.flag.dup:
            if var_reads.contains(r.qname):
                var mismatches = tag[int](r, "NM")
                if not mismatches.isNone:
                    edit_dists.add(mismatches.get)
    if len(edit_dists) > 0:
        avgEditDist = sum(edit_dists) / int(edit_dists.len)
    return avgEditDist

proc get_avg_read_bal(bam: Bam, chrom: string, position: int64, var_reads: HashSet[string]): (float) = 
    var
        avgReadBal: float = 0
        read_bals = newSeq[float](0)
        strt_dist: int = 0
        stop_dist: int = 0

    for r in bam.query(chrom, int(position-1), int(position)):
        if not r.flag.qcfail and not r.flag.secondary and not r.flag.supplementary and not r.flag.dup:
            if var_reads.contains(r.qname):
                strt_dist = int(position-r.start+1)
                stop_dist = int(r.stop-position+1)
                read_bals.add(abs(ln(strt_dist/stop_dist)))
    if len(read_bals) > 0:
        avgReadBal = sum(read_bals) / float(read_bals.len)
    return avgReadBal

var
    chrom : string
    position : int64
    REF : string
    ALT : string
    var_type : string
    depth : int
    alt_count : int
    ref_count : int
    vaf: float
    avgBq : float
    avgMq : float
    alt_alleles: int
    avgAltBq : float
    avgAltMq : float
    avgAltMateMq : float
    avgIs : float
    avgAltIs : float
    avgEditDist : float
    avgReadBal : float
    indel_count : float
    rf : int
    af : int
    rr : int
    ar : int
    var_reads = initHashSet[string]()

# write header
output.write(fmt("CHR\tSTART\tEND\tREF\tALT\tDEPTH\tVAF\tAVG_BQ\tAVG_ALT_BQ\tAVG_MQ\tAVG_ALT_MQ\tAVG_ALT_MATE_MQ\tAVG_IS\tAVG_ALT_IS\tAVG_EDIT_DIST\tAVG_READ_BAL\tVARIANT_READS\tVARIANT_ALLELES\tFR\tFA\tRR\tRA"))

if snvs:
    output.write(fmt("\n"))
    var_type = "SNV"
    for variant in vcf:
        chrom = $variant.CHROM
        position = variant.POS
        REF = $variant.REF
        ALT = variant.ALT[0]

        let h = bam.hileup(chrom, int(position-1), fai, cfg, snv=true)
        depth = len(h.mqs)

        if depth >= depth_thresh:
            var_reads = get_var_reads(h, ALT, var_type)
            # VAF
            alt_count = count(map(h.bases, proc(x: basestrand): string = $(x.base.char)), ALT)
            ref_count = count(map(h.bases, proc(x: basestrand): string = $(x.base.char)), REF)
            vaf = alt_count/depth
            # STRAND INFO
            (rf, af, rr, ar) = strand_info_snv(h, REF, ALT)
            # QUALITIES
            let bqs = map(h.bqs, proc(x: uint8): int = int(x))
            let mqs = map(h.mqs, proc(x: uint8): int = int(x))
            avgBq = sum(bqs) / int(bqs.len)
            avgMq = sum(mqs) / int(mqs.len)
            (avgAltBq, avgAltMq) = get_avg_alt_qualities(h, var_reads)
            avgAltMateMq = get_avg_alt_mate_mapq(bam, chrom, position, var_reads)
            # INSERT SIZE
            (avgIs, avgAltIs) = get_avg_insert_sizes(bam, chrom, position, var_reads)
            avgEditDist = get_avg_edit_dist(bam, chrom, position, var_reads)
            avgReadBal = get_avg_read_bal(bam, chrom, position, var_reads)
            alt_alleles = get_alt_alleles(h.bases, REF)
            output.write(fmt("{chrom}\t{position}\t{position}\t{REF}\t{ALT}\t{depth}\t{vaf}\t{avgBq}\t{avgAltBq}\t{avgMq}\t{avgAltMq}\t{avgAltMateMq}\t{avgIs}\t{avgAltIs}\t{avgEditDist}\t{avgReadBal}\t{alt_count}\t{alt_alleles}\t{rf}\t{af}\t{rr}\t{ar}\n"))
        bar.increment()
    bar.finish()
    output.close()
    close(bam)
    close(vcf)
    echo now() - time

else:
    output.write(fmt("\tINDEL_COUNT\n"))
    for variant in vcf:
        chrom = $variant.CHROM
        position = variant.POS
        REF = $variant.REF
        ALT = variant.ALT[0]

        let pileup = get_indel_pileup(bam, chrom, position, fai, cfg, REF, ALT)
        let h = pileup[0]
        var_type = pileup[1]
        depth = len(h.mqs)

        if depth >= depth_thresh:
            var_reads = get_var_reads(h, ALT, var_type)
            # VAF
            if var_type == "DEL":
                alt_count = len(h.del)
                ref_count = depth - alt_count
            else:
                alt_count = len(h.ins)
                ref_count = depth - alt_count
            vaf = alt_count/depth
            # STRAND INFO
            (rf, af, rr, ar) = strand_info_indel(h, var_reads)
            # QUALITIES
            let bqs = map(h.bqs, proc(x: uint8): int = int(x))
            let mqs = map(h.mqs, proc(x: uint8): int = int(x))
            avgBq = sum(bqs) / int(bqs.len)
            avgMq = sum(mqs) / int(mqs.len)
            (avgAltBq, avgAltMq) = get_avg_alt_qualities(h, var_reads)
            avgAltMateMq = get_avg_alt_mate_mapq(bam, chrom, position, var_reads)
            # INSERT SIZE
            (avgIs, avgAltIs) = get_avg_insert_sizes(bam, chrom, position, var_reads)
            avgEditDist = get_avg_edit_dist(bam, chrom, position, var_reads)
            avgReadBal = get_avg_read_bal(bam, chrom, position, var_reads)
            alt_alleles = 1
            # INDEL COUNT
            indel_count = get_indel_count(bam, chrom, position, alt_count)
            output.write(fmt("{chrom}\t{position}\t{position+len(ALT)}\t{REF}\t{ALT}\t{depth}\t{vaf}\t{avgBq}\t{avgAltBq}\t{avgMq}\t{avgAltMq}\t{avgAltMateMq}\t{avgIs}\t{avgAltIs}\t{avgEditDist}\t{avgReadBal}\t{alt_count}\t{alt_alleles}\t{rf}\t{af}\t{rr}\t{ar}\t{indel_count}\n"))
        bar.increment()
    bar.finish()
    output.close()
    close(bam)
    close(vcf)
    echo now() - time