#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Convert a BLAST XML file to tabular output.
https://github.com/peterjc/galaxy_blast/blob/master/tools/ncbi_blast_plus/blastxml_to_tabular.py
"""
from __future__ import print_function
import re
import sys


if sys.version_info[:2] >= (2, 5):
    try:
        from xml.etree import cElementTree as ElementTree
    except ImportError:
        from xml.etree import ElementTree as ElementTree
else:
    from galaxy import eggs  # noqa - ignore flake8 F401
    import pkg_resources

    pkg_resources.require("elementtree")
    from elementtree import ElementTree


re_default_query_id = re.compile(r"^Query_\d+$")
assert re_default_query_id.match(r"Query_101")
assert not re_default_query_id.match(r"Query_101a")
assert not re_default_query_id.match(r"MyQuery_101")
re_default_subject_id = re.compile(r"^Subject_\d+$")
assert re_default_subject_id.match(r"Subject_1")
assert not re_default_subject_id.match(r"Subject_")
assert not re_default_subject_id.match(r"Subject_12a")
assert not re_default_subject_id.match(r"TheSubject_1")


blastxml_filename = r"E:\Desktop\3.xml"
output_handle = open("E:\Desktop\out.txt", 'w')

blast_program = None
# get an iterable
try:
    context = ElementTree.iterparse(blastxml_filename, events=("start", "end"))
except Exception:
    sys.exit("Invalid data format.")
# turn it into an iterator
context = iter(context)
# get the root element
try:
    event, root = next(context)
except Exception:
    sys.exit("Invalid data format.")
for event, elem in context:
    if event == "end" and elem.tag == "BlastOutput_program":
        blast_program = elem.text
    # for every <Iteration> tag
    if event == "end" and elem.tag == "Iteration":
        # Expecting either this, from BLAST 2.2.25+ using FASTA vs FASTA
        # <Iteration_query-ID>sp|Q9BS26|ERP44_HUMAN</Iteration_query-ID>
        # <Iteration_query-def>Endoplasmic reticulum resident protein 44
        # OS=Homo sapiens GN=ERP44 PE=1 SV=1</Iteration_query-def>
        # <Iteration_query-len>406</Iteration_query-len>
        # <Iteration_hits></Iteration_hits>
        #
        # Or, from BLAST 2.2.24+ run online
        # <Iteration_query-ID>Query_1</Iteration_query-ID>
        # <Iteration_query-def>Sample</Iteration_query-def>
        # <Iteration_query-len>516</Iteration_query-len>
        # <Iteration_hits>...
        qseqid = elem.findtext("Iteration_query-ID")
        if re_default_query_id.match(qseqid):
            # Place holder ID, take the first word of the query definition
            qseqid = elem.findtext("Iteration_query-def").split(None, 1)[0]
        qlen = int(elem.findtext("Iteration_query-len"))

        # for every <Hit> within <Iteration>
        for hit in elem.findall("Iteration_hits/Hit"):
            # Expecting either this,
            # <Hit_id>gi|3024260|sp|P56514.1|OPSD_BUFBU</Hit_id>
            # <Hit_def>RecName: Full=Rhodopsin</Hit_def>
            # <Hit_accession>P56514</Hit_accession>
            # or,
            # <Hit_id>Subject_1</Hit_id>
            # <Hit_def>gi|57163783|ref|NP_001009242.1|
            # rhodopsin [Felis catus]</Hit_def>
            # <Hit_accession>Subject_1</Hit_accession>
            #
            # apparently depending on the parse_deflines switch
            #
            # Or, with a local database not using -parse_seqids can get this,
            # <Hit_id>gnl|BL_ORD_ID|2</Hit_id>
            # <Hit_def>chrIII gi|240255695|ref|NC_003074.8| Arabidopsis
            # thaliana chromosome 3, complete sequence</Hit_def>
            # <Hit_accession>2</Hit_accession>
            sseqid = hit.findtext("Hit_id").split(None, 1)[0]
            hit_def = sseqid + " " + hit.findtext("Hit_def")
            if re_default_subject_id.match(sseqid) and sseqid == hit.findtext(
                "Hit_accession"
            ):
                # Place holder ID, take the first word of the subject definition
                hit_def = hit.findtext("Hit_def")
                sseqid = hit_def.split(None, 1)[0]
            if sseqid.startswith(
                "gnl|BL_ORD_ID|"
            ) and sseqid == "gnl|BL_ORD_ID|" + hit.findtext("Hit_accession"):
                # Alternative place holder ID, again take the first word of hit_def
                hit_def = hit.findtext("Hit_def")
                sseqid = hit_def.split(None, 1)[0]
            # for every <Hsp> within <Hit>
            for hsp in hit.findall("Hit_hsps/Hsp"):
                nident = hsp.findtext("Hsp_identity")
                length = hsp.findtext("Hsp_align-len")
                # As of NCBI BLAST+ 2.4.0 this is given to 3dp (not 2dp)
                pident = "%0.3f" % (100 * float(nident) / float(length))

                q_seq = hsp.findtext("Hsp_qseq")
                h_seq = hsp.findtext("Hsp_hseq")
                m_seq = hsp.findtext("Hsp_midline")
                assert len(q_seq) == len(h_seq) == len(m_seq) == int(length)
                gapopen = str(
                    len(q_seq.replace("-", " ").split())
                    - 1
                    + len(h_seq.replace("-", " ").split())
                    - 1
                )

                mismatch = (
                    m_seq.count(" ")
                    + m_seq.count("+")
                    - q_seq.count("-")
                    - h_seq.count("-")
                )
                # TODO - Remove this alternative mismatch calculation and test
                # once satisifed there are no problems
                expected_mismatch = len(q_seq) - sum(
                    1
                    for q, h in zip(q_seq, h_seq)
                    if q == h or q == "-" or h == "-"
                )
                xx = sum(1 for q, h in zip(q_seq, h_seq) if q == "X" and h == "X")
                if not (
                    expected_mismatch - q_seq.count("X")
                    <= int(mismatch)
                    <= expected_mismatch + xx
                ):
                    sys.exit(
                        "%s vs %s mismatches, expected %i <= %i <= %i"
                        % (
                            qseqid,
                            sseqid,
                            expected_mismatch - q_seq.count("X"),
                            int(mismatch),
                            expected_mismatch,
                        )
                    )

                # TODO - Remove this alternative identity calculation and test
                # once satisifed there are no problems
                expected_identity = sum(1 for q, h in zip(q_seq, h_seq) if q == h)
                if not (
                    expected_identity - xx
                    <= int(nident)
                    <= expected_identity + q_seq.count("X")
                ):
                    sys.exit(
                        "%s vs %s identities, expected %i <= %i <= %i"
                        % (
                            qseqid,
                            sseqid,
                            expected_identity,
                            int(nident),
                            expected_identity + q_seq.count("X"),
                        )
                    )

                evalue = hsp.findtext("Hsp_evalue")
                if evalue == "0":
                    evalue = "0.0"
                else:
                    evalue = "%0.0e" % float(evalue)

                bitscore = float(hsp.findtext("Hsp_bit-score"))
                if bitscore < 100:
                    # Seems to show one decimal place for lower scores
                    bitscore = "%0.1f" % bitscore
                else:
                    # Note BLAST does not round to nearest int, it truncates
                    bitscore = "%i" % bitscore

                values = [
                    qseqid,
                    sseqid,
                    pident,
                    length,  # hsp.findtext("Hsp_align-len")
                    str(mismatch),
                    gapopen,
                    hsp.findtext("Hsp_query-from"),  # qstart,
                    hsp.findtext("Hsp_query-to"),  # qend,
                    hsp.findtext("Hsp_hit-from"),  # sstart,
                    hsp.findtext("Hsp_hit-to"),  # send,
                    evalue,  # hsp.findtext("Hsp_evalue") in scientific notation
                    bitscore,  # hsp.findtext("Hsp_bit-score") rounded
                ]


                try:
                    sallseqid = ";".join(
                        name.split(None, 1)[0] for name in hit_def.split(" >")
                    )
                    salltitles = "<>".join(
                        name.split(None, 1)[1] for name in hit_def.split(" >")
                    )
                except IndexError as e:
                    sys.exit(
                        "Problem splitting multuple hits?\n%r\n--> %s"
                        % (hit_def, e)
                    )
                # print(hit_def, "-->", sallseqid)
                positive = hsp.findtext("Hsp_positive")
                ppos = "%0.2f" % (100 * float(positive) / float(length))
                qframe = hsp.findtext("Hsp_query-frame")
                sframe = hsp.findtext("Hsp_hit-frame")
                if blast_program == "blastp":
                    # Probably a bug in BLASTP that they use 0 or 1
                    # depending on format
                    if qframe == "0":
                        qframe = "1"
                    if sframe == "0":
                        sframe = "1"
                slen = int(hit.findtext("Hit_len"))
                values.extend(
                    [
                        sallseqid,
                        hsp.findtext("Hsp_score"),  # score,
                        nident,
                        positive,
                        hsp.findtext("Hsp_gaps"),  # gaps,
                        ppos,
                        qframe,
                        sframe,
                        # NOTE - for blastp, XML shows original seq,
                        # tabular uses XXX masking
                        q_seq,
                        h_seq,
                        str(qlen),
                        str(slen),
                        salltitles,
                    ]
                )


                output_handle.write("\t".join(values) + "\n")
        # prevents ElementTree from growing large datastructure
        root.clear()
        elem.clear()




