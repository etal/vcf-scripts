#!/usr/bin/env python

"""Add CNV log2 ratio values to each SNV in a VCF.

"""
from __future__ import division, print_function

import sys
from os.path import basename

import pandas as pd
import vcf
from vcf.model import make_calldata_tuple, _Call as Call
from vcf.parser import _Format as Format


def logratio_at_site(record, segs, seen_skipped_chroms=set()):
    """Get the segment log2 ratio at the site of the SNV record."""
    chrom_segs = segs[segs['chromosome'] == record.CHROM]

    if len(chrom_segs) == 0:
        if record.CHROM not in seen_skipped_chroms:
            print("Skipping chromosome", record.CHROM, file=sys.stderr)
            seen_skipped_chroms.add(record.CHROM)
        return

    seg_idx = max(0, chrom_segs.start.searchsorted(record.start)[0] - 1)
    segment = chrom_segs.iloc[seg_idx]
    if (0 < seg_idx < len(chrom_segs) - 1) and not (
        (segment.start < record.start < segment.end) and
        segment.chromosome == record.CHROM):
        # SNV is between segments, e.g. in/near centromere
        print("** WARNING **",
              "VCF record {}:{}".format(record.CHROM, record.POS - 1),
              "not in segment {}:{}-{}".format(segment.chromosome,
                                               segment.start,
                                               segment.end),
              "at chrom. segment index", seg_idx,
              "of", len(chrom_segs),
              file=sys.stderr)
        return

    return segment.log2


def add_log2ratios_to_records(vcf_reader, segments_dict):
    """Add the CNV log2 ratio to each record's genotype for a sample."""
    for record in vcf_reader:
        orig_formats = record.FORMAT.split(':')
        if orig_formats[-1] == 'LR':
            # log2 ratios were already added to this file -- replace them
            CallData = make_calldata_tuple(orig_formats)
            orig_formats = orig_formats[:-1]
        else:
            # No existing log2 ratios (the common case)
            CallData = make_calldata_tuple(orig_formats + ['LR'])
            record.add_format("LR")
        out_samples = []
        for genotype in record.samples:
            if genotype.sample in segments_dict:
                segs = segments_dict[genotype.sample]
                site_log2 = logratio_at_site(record, segs)
            else:
                site_log2 = None
            # Create a new genotype field with LR added for this sample
            call_vals = [getattr(genotype.data, fmt, None)
                         for fmt in orig_formats]
            calldata = CallData(*(call_vals + [site_log2]))
            out_samples.append(Call(record, genotype.sample, calldata))
        record.samples = out_samples
        yield record


def get_segment_ids(segments, seg_ids, vcf_ids, sample_ids):
    """Map CNV segments to VCF sample IDs.

    Ensure every sample in the VCF header has a corresponding segment file.

    Return a dict of: { VCF sample IDs: DataFrame of CNV segments }
    """
    if sample_ids:
        if not set(sample_ids).issubset(vcf_ids):
            raise RuntimeError("Specified sample IDs %s are not all in VCF: %s"
                               % (", ".join(sample_ids), " ".join(vcf_ids)))
    else:
        if len(vcf_ids) == 1:
            # Use the only sample ID available in the VCF; ignore seg. filename
            sample_ids = vcf_ids
        else:
            # Otherwise, try to match the segments filename base
            vcf_ids = set(vcf_ids)
            if set(seg_ids).issubset(vcf_ids):
                sample_ids = seg_ids
            else:
                # Try to match segment filenames to VCF header
                sample_ids, extra_vcf_ids = update_ids(seg_ids, vcf_ids)
                if extra_vcf_ids:
                    raise ValueError("Segment filenames don't match VCF samples:"
                                     + "\n%s\nvs.\n%s" % (seg_ids, vcf_ids))

                # XXX Is this ever reached?
                extra_segments = sorted(set(sample_ids) - vcf_ids)
                if extra_segments:
                    print("Found extra segments which will not be used:",
                          " ".join(extra_segments), file=sys.stderr)
                    # Drop the extra sample_ids and segments
                    sample_ids = list(vcf_ids)
                    segments = [segments[sample_ids.index(s)] for s in vcf_ids]

        print("Selecting sample ID(s)", sample_ids, file=sys.stderr)
    assert len(sample_ids) == len(segments)
    return dict(zip(sample_ids, segments))


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
                        "Ambiguous VCF sample ID %s matches both segment "
                        "sample IDs %s and %s" % (eid, sid, seen_exids[eid]))
                print("Matched segment sample ID", sid, "to VCF's", eid)
                out.append(eid)
                seen_exids[eid] = sid
                break
        else:
            print("Segment sample %s not in VCF header; skipping" % sid)
    return out, exids - set(seen_exids)


def main(args):
    """."""
    vr = vcf.Reader(filename=args.snvs)
    segments = [pd.read_table(f, na_filter=False,
                              dtype={'chromosome': 'string'})
                for f in args.cnvs]
    seg_ids = [basename(f).split('.', 1)[0] for f in args.cnvs]
    seg_dict = get_segment_ids(segments, seg_ids, vr.samples, args.sample_ids)
    if not args.output:
        args.output = args.snvs.rsplit('.', 1)[0] + "log2ratios.vcf"
    idx = 0
    # Add log2 ratio FORMAT to header
    vr.formats["LR"] = Format("LR", 1, "Float", "CNV log2 ratio")
    with open(args.output, 'w') as outstream:
        vw = vcf.Writer(outstream, vr)
        for idx, record in enumerate(add_log2ratios_to_records(vr, seg_dict)):
            vw.write_record(record)
    print("Wrote", args.output, "with", idx + 1, "records", file=sys.stderr)


if __name__ == '__main__':
    import argparse
    AP = argparse.ArgumentParser(description=__doc__)
    AP.add_argument("snvs", help="VCF file containing SNVs.")
    AP.add_argument("cnvs", nargs='+', help="CNVkit segments.")
    AP.add_argument("-i", "--sample-ids", action="append", default=[],
                    help="Sample ID used in the VCF.")
    AP.add_argument("-o", "--output", help="Output VCF filename.")
    main(AP.parse_args())
