"""To define micro-homology or repeat-mediated deletions in Indel list."""
#!/usr/bin/python

# pylint: skip-file

# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Authors                : Manasa Ramakrishna
# Date created           : 18/01/2013
# Aim                    : To define micro-homology or repeat-mediated deletions in Indel list
# Last modification      : 4/02/2013, Introduced elif in line 213; introduced extra test between lines 144-155; Call accuracy up to 100% from 94% (based on 350 calls)
# 			   15/5/2014, Changing code such that if there is 1 repeat and 1 microhomology, then overall call is repeat-mediated, mr9
# Structures/loops used  : Functions from http://stackoverflow.com/questions/3551423/python-comparing-two-strings
# Make the python script executable: chmod +x microrep.py
# Usage                  : ./microrep.py -l debug -w /nfs/users/nfs_m/mr9/Unix/08_Python/02_Del-annot -i test_with-context.txt -o test-with-mh.txt
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------

# no longer annotates complex inderls

import argparse
import csv

# import difflib
import logging
import os
import re
import sys

s = re.compile(r"(.+)\1+")


LOGGING_LEVELS = {
    "critical": logging.CRITICAL,
    "error": logging.ERROR,
    "warning": logging.WARNING,
    "debug": logging.DEBUG,
    "info": logging.INFO,
}


# =============================================
# This is the part that runs the code
# =============================================

# Function used to call microhomology or repeat-mediated calls


def mh(a, b):  # noqa
    data = csv.DictReader(
        open(a, "r", encoding="utf-8"), delimiter="\t", quoting=csv.QUOTE_NONE
    )
    header = "\t".join(data.fieldnames)
    # print header

    of = open(b, "w")
    of.write(header + "\tMH-count\tMicrohomology\tRep-count\tRepeat\tClassification\n")

    for line in data:

        info = (
            line["VARIANT_ID"]
            + "\t"
            + line["VARIANT_TYPE"]
            + "\t"
            + line["ACCESSION"]
            + "\t"
            + line["LENGTH"]
            + "\t"
            + line["REPEATS"]
            + "\t"
            + line["MIN_POSITION"]
            + "\t"
            + line["MAX_POSITION"]
            + "\t"
            + line["SUM_MS"]
            + "\t"
            + line["CHANGE"]
            + "\t"
            + line["5'-context"]
            + "\t"
            + line["Actual_Change"]
            + "\t"
            + line["3'-context"]
            + "\t"
        )

        if line["VARIANT_TYPE"] == "D":
            a = line["Actual_Change"]  # The actual deletion
            b = line["3'-context"]  # Sequence 3' to deletion
            c = line["5'-context"]  # Sequence 5' to deletion

            # Look for Microhomology first and then for repeats - tandem/normal

            (mhcount, mh) = mhcaller(a, b)
            (repcount, repeat) = repcaller(a, b, c, line["LENGTH"])

            finalcall = finalcaller(mhcount, repcount * len(repeat), repeat)
            # print a, b, mhcount, mh, repcount*len(repeat), repeat, finalcall

            output = "\t".join([str(mhcount), mh, str(repcount), repeat, finalcall])
            # print output

            of.write(info + output + "\n")

        else:
            of.write(info + "-\t-\t-\t-\t-\n")
            # print a, b, line['VARIANT_TYPE'], line['LENGTH']


# ------------------------------------------------
# This is the Microhomology caller
# ------------------------------------------------

# Taken in the deletion sequence, the 3' context


def mhcaller(d, prime3):  # noqa
    countmh = 0
    seq = "-"
    # Dealing with microhomology or lack of microhmomology in the first position

    for i in range(len(d)):
        if d[i] == prime3[i]:
            # print d[i],prime3[i]
            countmh = countmh + 1
            seq = d[0:countmh]
        else:
            break

    return (countmh, seq)


# ------------------------------------------------
# This is the Repeat caller
# ------------------------------------------------

# Taken in the deletion sequence, the 3' context and the length of the deletion
# I have joined deletion and 3' context. To remove this, the rules for microhom/rep calling need to change


def repcaller(d, prime3, prime5, l):  # noqa
    countrep = 0
    # c = "".join([d,prime3])

    # This is for counting single base repeats
    if int(l) == 1:
        # print d
        for i in range(len(prime3)):
            if d[0] == prime3[i]:
                # print d[i],prime3[i]
                countrep = countrep + 1
                continue
            else:
                break
        return (countrep, d)

    # This is for counting whole deletion/DI repeats that are present in the flanking region
    elif prime3.find(d) == 0:
        countrep = tandemcount(d, prime3)
        rep = findsmallestrep(d)
        # print countrep, rep, tandemcount(rep,d), tandemcount(rep,prime3), "Whole indel repeat"
        countrep = max(countrep, tandemcount(rep, prime3))
        return (countrep, rep)

    else:
        # This is for counting anything in between 1bp and whole del repetition

        rep = "-"

        # Look for repeats of 2bp to n-1 bp of the length of the indel
        for t in range(len(d))[1:][::-1]:
            # print d[:t+1]
            if prime3.find(d[: t + 1]) == 0:

                countrep = tandemcount(d[: t + 1], prime3)
                rep = findsmallestrep(d[: t + 1])
                unit = tandemcount(rep, d) * len(rep)

                # The false calls arise in examples such as these : del = AACCCCATCTCTACT; 3' = AAAATTACAAACAAAT; rep = 'A'; repcount in 3' = 4 which is greater than MH = 2; Therefore, it is called Repeat-mediated
                # In fact, it should be repeat count = 0; Therefore, call should be microhomology mediated.
                # To do this, compare check how far the repeat stretched into the indel itself. Eg: 'A' is counted twice in the deletion. So compare del[:2] to del[2:4]. If they are the same,then keep it, else false

                if d[:unit] == d[unit : unit * 2]:
                    # print countrep, rep, tandemcount(rep,d), tandemcount(rep,prime3), "Repeat", unit, d[:unit], d[unit:unit*2]
                    countrep = max(countrep, tandemcount(rep, prime3))
                    return countrep, rep

                else:
                    return 0, "-"

            else:
                continue

        return countrep, rep


# -----------------------------------------------------------------------------------------------------
# Given a pattern and a string, counts how many times the pattern occurs in the string in tandem
# Eg: pattern = TG ; string = TGTGTGTGTGTGTTTT; result = 6
# -----------------------------------------------------------------------------------------------------


def tandemcount(pat, string):  # noqa
    sum = 0
    for k in range(0, len(string), len(pat)):
        if string[k : k + len(pat)] == pat:
            sum = sum + 1
        else:
            break
    return sum


# ----------------------------------------------------------
# This finds the smallest repeating subunit of the deletion
# ----------------------------------------------------------


def findsmallestrep(d):  # noqa
    # print d
    starts = []
    count = 1
    rep = ""
    #   s = re.compile(r"(.+)\1+")
    # Checking if deletions are repeat-mediated
    t = s.findall(d)

    # print t, len(t)

    if len(t) > 0:

        for i in range(len(t)):
            count = d.count(t[i])
            e = re.finditer(t[i], d)

            for m in e:
                starts.append(m.start())

            # If the start position of the repeat is the same as range of increasing units of repeat length, keep it
            if starts == range(0, len(d), len(t[i])):
                # print starts, range(0,len(d),len(t[i]))
                # print "===================="
                rep = t[i]
                rep = findsmallestrep(rep)
                # print rep

            elif len(t[i]) > 1:
                rep = findsmallestrep(t[i])
                # rep = d
                # print "===================="

            else:
                # print d
                rep = d
                break

        return rep

    else:
        return d


# ------------------------------------------
# This makes the final call on Deletion type
# ------------------------------------------
def finalcaller(mhcount, repcount, repeat):  # noqa

    # Given microhomology is in bases, repeat should be in bases too
    # If there is a single repeat of the indel 3' of it, then it should be labelled as Repeat-mediated not MH. Eg : TTTA	TTTATTATTAAGATTTTTAAATTTTAATT has 4bp MH and 1 repeat of TTTA. Counting itself, this is a repeat of 2, so it is repeat-mediated.15.05.14
    # Except if it is a single base from a longer indel that is repeating, then it is treated as MH

    if repcount >= mhcount:
        if repcount / len(repeat) >= 1:
            return "Repeat-mediated"
        elif mhcount > 0:
            return "Microhomology-mediated"
        else:
            return "None"
    else:
        if mhcount > 0:
            return "Microhomology-mediated"
        else:
            return "None"


# -----------------------------
# This is the help guide
# -----------------------------


def main():  # noqa
    aparser = argparse.ArgumentParser(
        description="This code calls a deletion as a microhomology-mediated or repeated mediated deletion."
    )
    aparser.add_argument(
        "-i",
        dest="infile",
        metavar="FILE",
        default=False,
        help="File that contains context added to indels",
    )
    aparser.add_argument(
        "-o",
        dest="outfile",
        metavar="FILE",
        default=False,
        help="Output file to write results into.",
    )
    aparser.add_argument(
        "-w", dest="wdir", metavar="PATH", default=".", help="Working directory"
    )
    aparser.add_argument(
        "-l",
        dest="loglevel",
        metavar="STRING",
        default="info",
        help="Logging level: info/warning/debug/error/critical",
    )

    args = aparser.parse_args()
    input = os.path.abspath(args.infile)
    wdir = args.wdir
    output = os.path.join(wdir, args.outfile)

    # Set logging
    logging_level = LOGGING_LEVELS.get(args.loglevel)
    logging.basicConfig(
        level=logging_level,
        format="%(asctime)s %(levelname)s: %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )

    mh(input, output)


if __name__ == "__main__":
    main()
