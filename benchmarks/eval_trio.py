import sys
import argparse
import logging
import time


def pase_info(seq):
    info = {'SVLEN': 0, 'END': 0, "SVTYPE": '', "RE": 0, "CHR2": ''}
    for i in seq.split(';'):
        if i.split('=')[0] in ["SVLEN", "END", "RE"]:
            try:
                info[i.split('=')[0]] = abs(int(float(i.split('=')[1])))
            except:
                pass
        if i.split('=')[0] in ["CHR2"]:
            info[i.split('=')[0]] = i.split('=')[1]
        if i.split('=')[0] in ["SVTYPE"]:
            info[i.split('=')[0]] = i.split('=')[1][0:3]

    return info


def phase_GT(seq):
    i = seq.split(':')
    if i[0] in ['0/1', '1/0']:
        return 'het'
    elif i[0] == '1/1':
        return 'hom'
    else:
        return 'unknown'


def load_callset(path):
    callset = dict()
    file = open(path, 'r')
    for line in file:
        seq = line.strip('\n').split('\t')
        if seq[0][0] == '#':
            continue

        chr = seq[0]
        pos = int(seq[1])
        info = pase_info(seq[7])

        if info['SVTYPE'] in ['DEL', 'INS']:
            if info['SVTYPE'] not in callset:
                callset[info['SVTYPE']] = list()
            if info['SVLEN'] == 0:
                info['SVLEN'] = info['END'] - pos + 1
            if info['SVLEN'] < 50:
                continue
            callset[info['SVTYPE']].append([chr, pos, info['END'], info['SVLEN'], phase_GT(seq[9]), 0])

        # if info['SVTYPE'] == "BND":
        #     # continue
        #     if seq[4][0] == ']':
        #         form = ']]N'
        #         chr2 = seq[4].split(':')[0][1:]
        #         pos2 = int(seq[4].split(':')[1][:-2])
        #     elif seq[4][0] == '[':
        #         form = '[[N'
        #         chr2 = seq[4].split(':')[0][1:]
        #         pos2 = int(seq[4].split(':')[1][:-2])
        #     else:
        #         if seq[4][1] == ']':
        #             form = 'N]]'
        #             chr2 = seq[4].split(':')[0][2:]
        #             pos2 = int(seq[4].split(':')[1][:-1])
        #         else:
        #             form = 'N[['
        #             chr2 = seq[4].split(':')[0][2:]
        #             pos2 = int(seq[4].split(':')[1][:-1])
        #     if info['SVTYPE'] not in callset:
        #         callset[info['SVTYPE']] = list()
        #     if info['END'] == 0:
        #         info['CHR2'] = chr2
        #         info['END'] = pos2
        #
        #     callset[info['SVTYPE']].append([chr, pos, info['CHR2'], info['END'], form, phase_GT(seq[9]), 0])

    file.close()
    return callset


def eva_record(call_A, call_B, bias, offect, gt):
    # call_A 0/1
    # call_B 1/1
    for svtype in call_A:
        if svtype not in call_B:
            continue
        else:
            for i in call_B[svtype]:
                if i[-2] not in gt:
                    continue
                for j in call_A[svtype]:
                    if i[0] != j[0]:
                        continue
                    if svtype == 'INS':
                        if abs(i[1] - j[1]) <= offect and float(min(i[3], j[3]) / max(i[3], j[3])) >= bias:
                            i[-1] = 1
                    # elif svtype == 'BND':
                    #     if i[2] == j[2] and i[4] == j[4]:
                    #         # chromosome + strand
                    #         if abs(i[1] - j[1]) <= offect and abs(i[3] - j[3]) <= offect:
                    #             i[-1] = 1
                    elif svtype == 'DEL':
                        if max(i[1] - offect, j[1]) <= min(i[2] + offect, j[2]) and float(
                                min(i[3], j[3]) / max(i[3], j[3])) >= bias:
                            # callset[info['SVTYPE']].append([chr, pos, info['END'], info['SVLEN'], phase_GT(seq[9]), 0])
                            i[-1] = 1


def statistics_true_possitive(callset, SVTYPE, gt):
    record = 0
    true_record = 0
    denovo_rec = open('denovo_sniffles.bed', 'w')
    if SVTYPE == "ALL":
        for svtype in callset:
            for i in callset[svtype]:
                if i[-2] in gt:
                    record += 1
                    if i[-1] == 1:
                        true_record += 1
                    else:
                        # if i[0] in ref_list:
                        denovo_rec.write("%s\t%d\t%d\t%d\t%s\t%s\n" % (i[0], i[1], i[2], i[3], svtype, i[-2]))
                # else:
                #     if i[0] in ref_list:
                #         denovo_rec.write("%s\t%d\t%d\t%d\t%s\t%s\n" % (i[0], i[1], i[2], i[3], svtype, i[-2]))

    else:
        if SVTYPE not in callset:
            return record, true_record
        for i in callset[SVTYPE]:
            if i[-2] in gt:
                record += 1
                if i[-1] == 1:
                    true_record += 1

    return record, true_record


def main_ctrl(args):
    logging.info("Load SV callset of selected caller.")

    call_child = load_callset(args.F1)
    call_father = load_callset(args.MP)
    call_mother = load_callset(args.FP)

    logging.info("Evaluate accuracy and sensitivity.")
    eva_record(call_child, call_father, args.bias, args.offect, ['hom'])
    eva_record(call_child, call_mother, args.bias, args.offect, ['hom'])
    eva_record(call_father, call_child, args.bias, args.offect, ['hom', 'het'])
    eva_record(call_mother, call_child, args.bias, args.offect, ['hom', 'het'])
    svtype = ["DEL", "INS", 'ALL']
    for i in svtype:
        record, true_record = statistics_true_possitive(call_child, i, ['hom', 'het'])
        if record == 0:
            logging.info("%s-%s: %d\t%d\t%.2f." % ('F1', i, record, record-true_record, 0))
            continue
        logging.info("%s-%s: %d\t%d\t%.2f." % ('F1', i, record, record-true_record, 100 * float((record-true_record)/ record)))
        # record, true_record = statistics_true_possitive(call_father, i, ['hom'])
        # if record == 0:
        #     logging.info("%s-%s: %d\t%d\t%.2f." % ('MP', i, record, record-true_record, 0))
        #     continue
        # logging.info("%s-%s: %d\t%d\t%.2f." % ('MP', i, record, true_record, 100 * float(true_record / record)))
        # record, true_record = statistics_true_possitive(call_mother, i, ['hom'])
        # if record == 0:
        #     logging.info("%s-%s: %d\t%d\t%.2f." % ('FP', i, record, record-true_record, 0))
        #     continue
        # logging.info("%s-%s: %d\t%d\t%.2f." % ('FP', i, record, true_record, 100 * float(true_record / record)))


def main(argv):
    args = parseArgs(argv)
    setupLogging(False)
    # print args
    starttime = time.time()
    main_ctrl(args)
    logging.info("Finished in %0.2f seconds." % (time.time() - starttime))


USAGE = """\
	Evaluate SV callset generated by cuteSV
"""


def parseArgs(argv):
    parser = argparse.ArgumentParser(prog="Trio_eval", description=USAGE,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    # parser.add_argument("CALLER", type=str, help="Choose a caller")
    parser.add_argument("MP", type=str, help="Male parent callsets")
    parser.add_argument('FP', type=str, help="Female parent callsets")
    parser.add_argument('F1', type=str, help="Offspring callsets")
    # parser.add_argument("MP_c", type=int, help="Coverage of male parent callsets")
    # parser.add_argument('FP_c', type=int, help = "Coverage of female parent callsets")
    # parser.add_argument('F1_c', type=int, help = "Coverage of offspring callsets")
    parser.add_argument('-b', '--bias', help="Bias of overlaping.[%(default)s]", default=0.5, type=float)
    parser.add_argument('-o', '--offect', help="Offect of translocation overlaping.[%(default)s]", default=1000,
                        type=int)
    args = parser.parse_args(argv)
    return args


def setupLogging(debug=False):
    logLevel = logging.DEBUG if debug else logging.INFO
    logFormat = "%(asctime)s [%(levelname)s] %(message)s"
    logging.basicConfig(stream=sys.stderr, level=logLevel, format=logFormat)
    logging.info("Running %s" % " ".join(sys.argv))


if __name__ == '__main__':
    main(sys.argv[1:])
