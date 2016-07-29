#!/usr/bin/env python

"""Create the merged VCF. Combine all callers' genotypes for each actual sample.
"""
from __future__ import print_function, division

import collections
import sys
from itertools import chain

import numpy as np
import vcf
from vcf.model import make_calldata_tuple, _Record as Record, _Call as Call
from vcf.parser import _Format as Format

CNT_TYPICAL = 0
CNT_ATYPICAL = 0

# _____________________________________________________________________________
# Load all VCF data

def cleanse_callers(caller_data):
    """If a sample was called twice by a caller, choose the paired one.

    E.g. FreeBayes runs in single-sample and paired modes; if a paired calling
    is seen, drop the single-sample calls.

    In multi-tumor (recurrence) cases there could be multiple paired callings
    where the normal sample is reused. Keep those, but drop the singles.
    """
    for caller, cdata in caller_data.iteritems():
        print("Cleansing calls by", caller, file=sys.stderr)
        if len(cdata) == 1:
            idx, samples = cdata[0]
            print("Single calling of:", samples, file=sys.stderr)
            yield idx, samples, caller
        else:
            # Check for samples called repeatedly
            paired_samples = collections.defaultdict(list)
            single_samples = collections.defaultdict(list)
            for idx, samples in cdata:
                # NB: If multiple paired or single callings exist for the same
                # samples, keep them all (we'll squash them later)
                if len(samples) > 1:
                    paired_samples[tuple(samples)].append(idx)
                else:
                    single_samples[samples[0]].append(idx)
            if paired_samples:
                seen_paired_samples = set()
                for samples, idxes in paired_samples.iteritems():
                    for idx in idxes:
                        yield idx, samples, caller
                    seen_paired_samples.update(samples)
                print("Emitting paired samples:", *sorted(seen_paired_samples),
                      file=sys.stderr)
                # Ensure no single samples were missed
                seen_single_samples = set(single_samples)
                print("Skipped the singles:", *sorted(seen_single_samples),
                      file=sys.stderr)
                for missed_single in seen_single_samples - seen_paired_samples:
                    print("Recovering missed single", missed_single,
                          file=sys.stderr)
                    # Missed some ...
                    for idx in single_samples[missed_single]:
                        yield idx, [missed_single], caller
            else:
                print("Emitting single sample(s):",
                      *sorted(single_samples.keys()), file=sys.stderr)
                for sample, idxes in single_samples.iteritems():
                    for idx in idxes:
                        yield idx, [sample], caller


def combine_index(vcf_fnames):
    """Organize all VCF records into a data structure.

    Returns a list of the items collected by `index_add`, and a list of unique sample
    names from the VCF header.
    """
    data_pile = collections.defaultdict(list)
    all_idx = collections.OrderedDict()
    all_samples = set()
    vcf_to_samples={}
    for vcf_fname in vcf_fnames:
        vr = vcf.Reader(filename=vcf_fname)
        caller = guess_caller(vr)
        samples = clean_samples(vr)
        vr.samples = samples
        idx = index_records(vr)
        print(vcf_fname, len(idx), caller, *samples, sep='\t', file=sys.stderr)
        vcf_to_samples[vcf_fname] = samples
        data_pile[caller].append((idx, samples))
        all_samples.update(samples)
    for idx, samples, caller in cleanse_callers(data_pile):
        index_add(all_idx, idx, samples, caller)
    return sorted(all_idx.iteritems(), key=sortkey), sorted(all_samples), vcf_to_samples


def guess_caller(vr):
    """Detect the caller that created the given VCF."""
    if "source" in vr.metadata and len(vr.metadata["source"]) == 1:
        # Callers that follow the VCF spec: FreeBayes, pindel
        caller = vr.metadata["source"][0].split(None, 1)[0]
    elif "GATKCommandLine.MuTect" in vr.metadata:
        # GATK/SATK 3.4+
        caller = "MuTect"
    elif "GATKCommandLine.HaplotypeCaller" in vr.metadata:
        caller = "HaplotypeCaller"
    elif "GATKCommandLine.UnifiedGenotyper" in vr.metadata:
        caller = "UnifiedGenotyper"
    elif "GATKCommandLine" not in vr.metadata:
        raise ValueError("Bad VCF header missing caller info:\n%s"
                         % vr.metadata)
    else:
        if len(vr.metadata["GATKCommandLine"]) == 2:
            # It's "private" to UG vs. HC, via vcf_comp
            caller = "UnifiedGenotyper"
        else:
            # GATK tools don't follow the spec
            gatk_info = vr.metadata["GATKCommandLine"]
            assert len(gatk_info) == 1
            ##GATKCommandLine=<ID=UnifiedGenotyper,CommandLineOptions="...
            caller = gatk_info[0]["ID"]
    return caller


def clean_samples(vr):
    """Drop any filename suffixes from each sample name in the VCF header.

    Example: "CGP-993.deduplicated.realign.recal" -> "CGP-993"
    """
    return [name.split('.', 1)[0] for name in vr.samples]



def index_records(vr):
    """Create an ordered lookup of variant keys to sample records.

    Exhausts the VCF Reader iterable.
    """
    return collections.OrderedDict((record2key(rec), clean_sample_index(rec))
                                   for rec in vr)


def record2key(record):
    """Construct a tuple key from a VCF variant record."""
    # VcfKey = collections.namedtuple("VcfKey", "chrom pos ref alt")
    return (record.CHROM,
            record.POS,
            record.REF,
            str(record.ALT[0]))


def clean_sample_index(record):
    """Like `clean_samples`, but for individual records.
    """
    record._sample_indexes = {name.split('.', 1)[0]: i
                              for name, i in record._sample_indexes.iteritems()}
    return record


def index_add(all_index, this_index, samples, caller):
    """Add newly indexed records to the shared collection, in-place.

    The shared index looks like:
        { record_key: {
            sample_id: {
                caller: [ records ... ] } } }
    """
    for key, record in this_index.iteritems():
        if key not in all_index:
            all_index[key] = {}
        for sample_id in samples:
            if sample_id not in all_index[key]:
                all_index[key][sample_id] = {caller: []}
            elif caller not in all_index[key][sample_id]:
                all_index[key][sample_id][caller] = []
            # NB: If caller was run twice, will have 2 records here
            all_index[key][sample_id][caller].append(record)


def sortkey(item):
    """Make a Python-sortable element from the variant index key.

    Autosomes are sorted first, as integers, not strings.
    """
    chrom, pos, ref, alt = item[0]
    if chrom.startswith('chr'):
        chrom = chrom[3:]
    if chrom.isdigit():
        chrom = int(chrom)
    return (chrom, pos, len(ref), len(alt))


# _____________________________________________________________________________
# Merge and write VCF records

CALLER_NAMES = collections.OrderedDict([
    ("HaplotypeCaller", "CallHC"),
    ("UnifiedGenotyper", "CallUG"),
    ("freeBayes", "CallFB"),
    ("pindel", "CallPI"),
    ("SomaticIndelDetector", "CallSID"),
    ("MuTect", "CallMU")])

CALLER_CODES = CALLER_NAMES.values()


def caller_code_formats():
    for caller_name, callcode in CALLER_NAMES.iteritems():
        yield callcode, Format(callcode,   # id
                               1,          # num
                               'Integer',  # type
                               "Variant was called by " + caller_name)  # desc


def pedigree(tumors, normal):
    """Generate a PEDIGREE tag for the VCF header.

    Like:
        PEDIGREE=<Derived=TumorSample1,Original=NormalSample>
        PEDIGREE=<Derived=TumorSample2,Original=NormalSample>
    """
    print("Tumor sample(s):", *tumors, file=sys.stderr)
    print("Normal sample:", normal, file=sys.stderr)
    return ["<Derived={},Original={}>".format(tum, normal) for tum in tumors]


def consense_values(values):
    """Take the appropriate "consensus" of a list of values of any type."""
    values = filter_none(values)
    if not values:
        return None
    if len(values) == 1:
        return values[0]
    if all(isinstance(v, bool) for v in values):
        # e.g. SOMATIC
        return any(values)
    elif all(isinstance(v, list) for v in values):
        # List of numbers, e.g. PL = 11395,783,0
        return consense_lists(values)

    values = flatten(values)
    if all(isinstance(v, (int, float)) for v in values):
        # Numeric, e.g. DP, GMICOV, GMIMAF
        val = np.median(values)
        if all(isinstance(v, int) for v in values):
            val = int(val)
    else:
        # Default: Use the most frequently occurring (modal) value
        ctr = collections.Counter(values)
        val = ctr.most_common()[0][0]
    return val


def consense_lists(lists):
    """Consensus of a list of lists."""
    clean_columns = [filter_none(column) for column in zip(*lists)]
    vals = [consense_values(col) for col in clean_columns]
    if all(isinstance(c, int) for c in chain(*clean_columns)):
        vals = [int(v) for v in vals]
    return vals


def consense_dicts(dicts):
    """Consensus of a list of dictionaries."""
    all_keys = sorted(set(chain(*(dc for dc in dicts))))
    return {key: consense_values([dc[key] for dc in dicts if key in dc])
            for key in all_keys}


def consensus_of_gtype(calls, fmt):
    values = [getattr(c.data, fmt, None) for c in calls]
    if fmt == 'GMIMAF':
        # Values used to be fraction [0,1], now are percentages [0,100]
        pct_values = [(int(round(v * 100)) if 0 < v < 1 else v)
                      for v in filter_none(flatten(values))]
        if pct_values:
            try:
                return int(np.median(pct_values))
            except TypeError:
                print("GMIMAF:", values, file=sys.stderr)
                print("pct_values:", pct_values, file=sys.stderr)
                return pct_values[0]
        else:
            return None
    else:
        return consense_values(values)


def filter_none(elems):
    """Remove list elements that are None."""
    return [x for x in elems if x is not None]


def flatten(values):
    """Unpack any stray multiallelic calls to avoid crashing.

    Normally this should not be triggered in our pipeline, but let's not
    crash on old or outside VCFs.
    """
    values = list(values)
    if any(isinstance(v, list) for v in values):
        flatvals = []
        for v in values:
            if isinstance(v, list):
                flatvals.extend(v)
            else:
                flatvals.append(v)
        return flatvals
    else:
        return values


def unique(arr):
    """Remove duplicate entries from a sequence."""
    seen = set()
    for elem in arr:
        if elem not in seen:
            yield elem
            seen.add(elem)

# Options for switching a GT field to not-called or called
SWITCH_GT = {
    # Diploid substitutions
    "0/0": ("0/0", "0/1"),
    "0/1": ("0/0", "0/1"),
    "1/1": ("0/0", "1/1"),
    # Inserts
    "./1": ("0/.", "./1"),
    "0/.": ("0/.", "./1"),
    # Deletions, hemizygous calls?
    "0":   ("0", "1"),
}

def make_call_data(formats, gtypes, caller_flags):
    """Create a pyvcf-compatible namedtuple of a sample's genotype values."""
    CallData = make_calldata_tuple(formats + CALLER_CODES)
    call_vals = [consensus_of_gtype(gtypes, fmt) for fmt in formats]
    # Set GT field by whether any caller called the variant
    idx = formats.index("GT")
    if idx != 0:
        print("*** oddly ordered FORMAT:", formats, file=sys.stderr)
    orig_gt = call_vals[idx] or "0/1"
    is_called = int(any(caller_flags.viewvalues()))
    call_vals[idx] = SWITCH_GT[orig_gt][is_called]
    # Rebuild the genotype data structure
    return CallData(*(call_vals + caller_flags.values()))


def squash_records(key, sample_lookup, sample_ids, sample_indexes, normal_id):
    """Create an aggregate record, merging samples across callers."""
    chrom, pos, ref, alt = key
    ids = []
    quals = []
    filters = []
    infos = []
    formats = []
    gtypes = {sid: [] for sid in sample_ids}
    caller_flags = {sid: collections.OrderedDict([(key, None)
                                                  for key in CALLER_CODES])
                    for sid in sample_ids}
    for sample_id, caller_lookup in sample_lookup.iteritems():
        caller_ids = []
        caller_quals = []
        caller_filters = []
        caller_infos = []
        for caller, records in caller_lookup.iteritems():
            # Squash repeat calls on a sample (e.g. FreeBayes single/somatic)
            caller_ids.append(consense_values([rec.ID for rec in records]))
            caller_quals.append(consense_values([rec.QUAL for rec in records]))
            caller_filters.extend(chain(*[rec.FILTER
                                          for rec in records if rec.FILTER]))
            caller_infos.append(consense_dicts([rec.INFO for rec in records]))
            formats.extend(chain(*[rec.FORMAT.split(':') for rec in records]))
            calls = [rec.genotype(sample_id) for rec in records]
            gtypes[sample_id].extend(calls)  # ENH: consense
            caller_code = CALLER_NAMES[caller]
            caller_val = max(call.is_variant for call in calls)
            if caller_val is not None:
                caller_val = int(caller_val)
            caller_flags[sample_id][caller_code] = caller_val
        ids.extend(caller_ids)
        quals.extend(caller_quals)
        filters.extend(caller_filters)
        infos.extend(caller_infos)
    formats = list(unique(formats))
    out_record = Record(chrom, pos, consense_values(ids),
                        ref, [alt], consense_values(quals),
                        list(unique(filters)),
                        consense_dicts(infos),
                        ":".join(formats + CALLER_CODES),
                        sample_indexes)
    out_record.samples = [Call(out_record, sid,
                               make_call_data(formats, gtypes[sid],
                                              caller_flags[sid]))
                          for sid in sample_ids]

    # Post-processing
    # Fix AN and AC, which are defined specially in the VCF spec
    all_alleles = list(chain(*[gt.gt_alleles for gt in out_record.samples
                               if gt.gt_nums]))
    # AN (1): "Total number of alleles in called genotypes"
    out_record.INFO["AN"] = sum(a != '.' for a in all_alleles)
    # AC (A): "Total number of alternate alleles in called genotypes"
    # ENH: don't crash on multiallelics?
    out_record.INFO["AC"] = sum(sum(gta not in ('.', '0') for gta in gt.gt_alleles)
                                for gt in out_record.samples
                                if gt.gt_nums)

    # NB: Last 2 validation errors:
    # CGP-1240 FORMAT tag [GL] expected different number of values (expected 2, found 3)
    # chr2	209113113	.	G	A	2714.98	PASS	SOMATIC;NS=1;DP=653;DPB=716.75;AC=1;AN=3;AF=0.375;RO=590;AO=123;PRO=0.0;PAO=0.0;QR=21870;QA=4741;PQR=0.0;PQA=0.0;SRF=386;SRR=203;SAF=90;SAR=33;SRP=127.1;SAP=60.3689;AB=0.233397;ABP=328.364;RUN=1;RPP=24.6368;RPPR=13.5915625;RPL=44.0;RPR=79.0;EPP=4.44029;EPPR=8.549435;DPRA=1.56225;ODDS=272.57225;GTI=0;TYPE=snp;CIGAR=1X;NUMALT=1;MEANALT=2.0;LEN=1;MQM=60.0;MQMR=60.0177;PAIRED=0.95122;PAIREDR=0.958916;technology.ILLUMINA=1.0;MQ0=0;VT=SNP;BaseQRankSum=3.369;Dels=0.0;FS=2.338;HaplotypeScore=31.1133;MLEAC=1;MLEAF=0.5;MQ=60.02;MQRankSum=0.461;QD=7.05;ReadPosRankSum=-0.95;SOR=0.953;SF=0	GT:AD:BQ:DP:FA:SS:GMIMUT:GMIMAF:GMICOV:GQ:RO:QR:AO:QA:GL:PL:CallHC:CallUG:CallFB:CallPI:CallSID:CallMU:LR	0/1:387,117:38.0:506:0.239:2:117:23.0:504:120.456:401:14917:123:4741:-285.03,0.0,-1199.66:3594,0,14408:.:1:1:.:.:1:0.0723	0:200,0:.:227:0.0:0:0:0.0:226:141.912:252:9271:0:0:0.0,-66.9517,-832.671:.:.:.:0:.:.:0:0.0607
    # chr5	176715884	.	T	G	2652.62	PASS	SOMATIC;NS=1;DP=718;DPB=806.5;AC=1;AN=3;AF=0.375;RO=684;AO=122;PRO=0.0;PAO=0.0;QR=26587;QA=4642;PQR=0.0;PQA=0.0;SRF=334;SRR=350;SAF=62;SAR=60;SRP=3.803265;SAP=3.0815;AB=0.225508;ABP=357.065;RUN=1;RPP=5.57335;RPPR=3.3676575;RPL=67.0;RPR=55.0;EPP=4.14943;EPPR=9.662735;DPRA=1.1461875;ODDS=318.98275;GTI=0;TYPE=snp;CIGAR=1X;NUMALT=1;MEANALT=1.0;LEN=1;MQM=60.0;MQMR=60.0;PAIRED=0.991803;PAIREDR=0.97694275;technology.ILLUMINA=1.0;MQ0=0;VT=SNP;BaseQRankSum=-2.525;Dels=0.0;FS=0.0;HaplotypeScore=15.7349;MLEAC=1;MLEAF=0.5;MQ=60.0;MQRankSum=2.816;QD=6.74;ReadPosRankSum=1.308;SOR=0.691;SF=0	GT:AD:BQ:DP:FA:SS:GMIMUT:GMIMAF:GMICOV:GQ:RO:QR:AO:QA:GL:PL:CallHC:CallUG:CallFB:CallPI:CallSID:CallMU:LR	0/1:408,121:37.0:529:0.22:2:121:23.0:529:117.5855:419:16330:122:4642:-270.938,0.0,-1322.0:3592,0,16251:.:1:1:.:.:1:0.022	0:284,0:.:319:0.0:0:0:0.0:319:136.171:354:13676:0:0:0.0,-95.905,-1230.06:.:.:.:0:.:.:0:0.0006
    # * FreeBayes somatic produces the 3 GL values like this, w/e
    # chr2	209113113	.	G	A	2714.98	.	AB=0.233397;ABP=328.364;AC=1;AF=0.25;AN=4;AO=123;CIGAR=1X;DP=780;DPB=780;DPRA=2.083;EPP=4.44029;EPPR=9.74419;GTI=0;LEN=1;MEANALT=2;MQM=60;MQMR=60.0153;NS=2;NUMALT=1;ODDS=155.143;PAIRED=0.95122;PAIREDR=0.960184;PAO=0;PQA=0;PQR=0;PRO=0;QA=4741;QR=24188;RO=653;RPL=44;RPP=24.6368;RPPR=15.384;RPR=79;RUN=1;SAF=90;SAP=60.3689;SAR=33;SRF=424;SRP=129.458;SRR=229;TYPE=snp;technology.ILLUMINA=1;SOMATIC	GT:GQ:DP:RO:QR:AO:QA:GL:GMIMUT:GMIMAF:GMICOV	0/1:141.912:527:401:14917:123:4741:-285.03,0,-1199.66:123:23:524	0/0:141.912:253:252:9271:0:0:0,-66.9517,-832.671:0:0:252
    # chr5	176715884	.	T	G	2652.62	.	AB=0.225508;ABP=357.065;AC=1;AF=0.25;AN=4;AO=122;CIGAR=1X;DP=895;DPB=895;DPRA=1.52825;EPP=4.14943;EPPR=7.73248;GTI=0;LEN=1;MEANALT=1;MQM=60;MQMR=60;NS=2;NUMALT=1;ODDS=221.81;PAIRED=0.991803;PAIREDR=0.978008;PAO=0;PQA=0;PQR=0;PRO=0;QA=4642;QR=30006;RO=773;RPL=67;RPP=5.57335;RPPR=3.48505;RPR=55;RUN=1;SAF=62;SAP=3.0815;SAR=60;SRF=377;SRP=4.0244;SRR=396;TYPE=snp;technology.ILLUMINA=1;SOMATIC	GT:GQ:DP:RO:QR:AO:QA:GL:GMIMUT:GMIMAF:GMICOV	0/1:136.171:541:419:16330:122:4642:-270.938,0,-1322:122:23:541	0/0:136.171:354:354:13676:0:0:0,-95.905,-1230.06:0:0:354

    # Set SOMATIC flag (only add, don't remove existing flags)
    global CNT_ATYPICAL
    global CNT_TYPICAL
    if normal_id:
        normal_gt = out_record.genotype(normal_id)
        if ((not normal_gt.is_variant) and
            (not out_record.INFO.get("SOMATIC")) and
            any(gt.is_variant for gt in out_record.samples)):
            out_record.INFO["SOMATIC"] = True
            if normal_gt.gt_nums is None:
                CNT_TYPICAL += 1
            else:
                CNT_ATYPICAL += 1

    return out_record


def write_combined_vcf(outstream, example_vcf, record_tree, sample_ids,
                       tumors, normal):
    """Create the merged output VCF and write to the given handle."""
    template = vcf.Reader(filename=example_vcf)
    # Update the VCF metadata
    template.samples = sample_ids
    template.formats.update(caller_code_formats())
    template.metadata["source"] = ["SNV-Unifier"]
    if tumors and normal:
        print("Adding pedigree tag(s)", file=sys.stderr)
        template.metadata["PEDIGREE"] = pedigree(tumors, normal)

    # Merge & write variant records
    sample_indexes = {sid: i for i, sid in enumerate(sample_ids)}
    vw = vcf.Writer(outstream, template)
    prev_chrom = None
    for key, sample_lookup in record_tree:
        if key[0] != prev_chrom:
            print("Merging variants on chromosome", key[0], file=sys.stderr)
            prev_chrom = key[0]
        out_record = squash_records(key, sample_lookup, sample_ids,
                                    sample_indexes, normal)
        vw.write_record(out_record)


# _____________________________________________________________________________
# Extra stuff added by Courtney 10/16/15 to deal with multiple normal samples.

def find_good_normal(normals, vcf_to_samples):
    """Select one normal sample ID from several.

    Preferably one that was used in tumor-normal pairing; ideally there will be
    only one of those, but we're not picky here.
    """
    used_in_tn_pairing = {normal: any(len(vcf_ids) == 2 and normal in vcf_ids
                                      for vcf_ids in vcf_to_samples.values())
                          for normal in normals}
    for normal_id, is_paired in used_in_tn_pairing.viewitems():
        if is_paired:
            return normal_id
    return normals[0]


def remove_bad_vcfs(vcf_to_samples, bad_normals):
    """Keep the VCF objects that don't contain any bad-normal IDs."""
    return [vcf
            for vcf, vcf_sample_ids in vcf_to_samples.viewitems()
            if not any(bn_id in vcf_sample_ids for bn_id in bad_normals)]


# _____________________________________________________________________________

def update_ids(sids, exids):
    """Check if each specified ID is a substring of an expected ID.

    If so, return the corresponding expected IDs (found and remainder);
    otherwise err.
    """
    out = []
    seen_exids = dict()
    for sid in sids:
        for eid in exids:
            if sid in eid:
                if eid in seen_exids:
                    raise RuntimeError(
                        "Ambiguous VCF sample ID %s matches both specified "
                        "sample IDs %s and %s" % (eid, sid, seen_exids[eid]))
                print("Matched specified sample ID", sid, "to VCF's", eid,
                      file=sys.stderr)
                out.append(eid)
                seen_exids[eid] = sid
                break
        else:
            raise RuntimeError("Specified sample %s not in VCF header" % sid)
    return out, exids - set(seen_exids)


def check_sample_ids(tumor_ids, normal_ids, expected_ids):
    """Ensure input sample IDs match the VCF header / filenames.

    If the specified IDs are uniquely truncated versions of the expected IDs
    found in the VCF header (e.g. "CGP-1620" vs. "CGP-1620+CGP-1620"), then
    update the specified IDs to match the expected IDs.
    """
    if tumor_ids and normal_ids:
        arg_ids = set(tumor_ids + normal_ids)
        if arg_ids.symmetric_difference(expected_ids):
            print("Specified sample IDs:",
                  ", ".join(tumor_ids), "(tumor) /",
                  ", ".join(normal_ids), "(normal);"
                  "but found", ", ".join(expected_ids),
                  "in VCF header", file=sys.stderr)
            tumor_ids, expect_nids = update_ids(tumor_ids, expected_ids)
            normal_ids, remaining_ids = update_ids(normal_ids, expect_nids)
            if remaining_ids:
                raise RuntimeError("Didn't specify sample type of: "
                                   + ", ".join(remaining_ids))
        return tumor_ids, normal_ids

    if tumor_ids or normal_ids:
        print("Ignoring pedigree; didn't specify both tumor and normal",
              file=sys.stderr)
    return [], []


def main(args):
    # Load
    record_tree, samples, vcf_to_samples = combine_index(args.vcfs)
    print("Indexed", len(record_tree), "variants from", len(args.vcfs), "VCFs",
          file=sys.stderr)

    # Validate
    tumor_ids, normal_ids = check_sample_ids(args.tumor, args.normal,
                                             set(samples))
    if tumor_ids and normal_ids:
        # We can only handle one normal in the output, but to avoid crashing we
        # accept multiple normal IDs as input, then select one.
        if len(normal_ids) == 1:
            normal_id = normal_ids[0]
        else:
            # Get rid of the extra normal sample(s)
            print("Too many normal samples:", ", ".join(normal_ids),
                  file=sys.stderr)
            bad_normals = normal_ids[:]
            normal_id = find_good_normal(normal_ids, vcf_to_samples)
            bad_normals.remove(normal_id)
            print("Dropping bad/extra normal sample(s)", ", ".join(bad_normals),
                  "from input VCF list", file=sys.stderr)
            args.vcfs = remove_bad_vcfs(vcf_to_samples, bad_normals)
            record_tree, samples, vcf_to_samples = combine_index(args.vcfs)
    else:
        normal_id = None

    # Merge & write
    write_combined_vcf(args.output, args.vcf_header or args.vcfs[0],
                       record_tree, samples, tumor_ids, normal_id)

    print("Flagged", CNT_TYPICAL + CNT_ATYPICAL, "additional somatic variants,",
          "of which", CNT_ATYPICAL, "were genotyped in normal and",
          CNT_TYPICAL, "were not.", file=sys.stderr)



if __name__ == '__main__':
    import argparse
    AP = argparse.ArgumentParser(description=__doc__)
    AP.add_argument("vcfs", nargs='+', help="VCF files to merge.")
    AP.add_argument("-t", "--tumor", action='append', default=[],
                    help="Tumor sample ID (can be used >1x)")
    AP.add_argument("-n", "--normal", action='append', default=[],
                    help="Normal sample ID (should not be used >1x)")
    AP.add_argument("-v", "--vcf-header", # default="header.vcf",
                    help="VCF header file to use in output")
    AP.add_argument("-o", "--output", help="Output VCF filename",
                    default=sys.stdout, type=argparse.FileType("w"))
    main(AP.parse_args())
