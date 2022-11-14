# -*- coding: UTF-8 -*-
"""
__project_ = 'GenePredict'
__file_name__ = 'create_cnn_model'
__author__ = 'Dyh'
__time__ = '2022/2/15 10:07'
__product_name = PyCharm
# code is far away from bugs with the god animal protecting
    I love animals. They taste delicious.
              ┏┓      ┏┓
            ┏┛┻━━━┛┻┓
            ┃        ┃
            ┃  ┳┛  ┗┳  ┃
            ┃      ┻      ┃
            ┗━┓      ┏━┛
                ┃      ┗━━━┓
                ┃  神兽保佑    ┣┓
                ┃　永无BUG！   ┏┛
                ┗┓┓┏━┳┓┏┛
                  ┃┫┫  ┃┫┫
                  ┗┻┛  ┗┻┛
"""
import re
import random
from Bio import SeqIO
from .CNN_train import CNN_model


def get_all_mrna_information(gff_file):
    # 打开基因注释文件 逐行阅读
    print('正在打开基因注释文件')
    f = open(gff_file, 'r')
    line = f.readline()
    print('打开完毕\n正在提取信息')
    # 将每行信息进行切割
    # NC_047034.1	Gnomon	mRNA	234056	239154	.	+	.	ID=rna-XM_019945566.2;
    # NC_047034.1	Gnomon	exon	234056	234153	.	+	.	ID=exon-XM_019945566.2-1;
    each_line_infor = line.split('\t')

    all_mRNA = []  # 用于保存全部mRNA的信息
    while line != '':
        # 每条mRNA的信息 {'local': 'NW_022983470.1', 'direction': '+', 'ID': 'ID=rna-XM_033850496.1', 'exon': [[106415, 106458], [107000, 107804]]}
        single_mRNA = {}
        if len(each_line_infor) > 2 and each_line_infor[2] == 'mRNA':
            single_mRNA['local'] = each_line_infor[0]  # 位于基因组的位置 NW_022983470.1
            single_mRNA['direction'] = each_line_infor[6]  # mRNA的转录方向 + or -
            single_mRNA['ID'] = each_line_infor[8].split(';')[0]  # mRNA的登录号
            line = f.readline()
            each_line_infor = line.split('\t')

            # 保存mRNA每个外显子的具体位置 'exon': [[106415, 106458], [107000, 107804]]
            single_mRNA['exon'] = []
            while len(each_line_infor) > 3:
                if each_line_infor[2] != 'exon':
                    break
                single_mRNA['exon'].append([int(each_line_infor[3]), int(each_line_infor[4])])
                line = f.readline()
                each_line_infor = line.split('\t')
            # 剔除只含有一个外显子的mRNA
            if len(single_mRNA['exon']) != 1:
                all_mRNA.append(single_mRNA)
        else:
            line = f.readline()
            each_line_infor = line.split('\t')
    f.close()
    return all_mRNA


def get_data_seqs(all_mRNA, genome_file: str, train_num):
    # 打开数据
    all_mRNA_infor = all_mRNA
    print(len(all_mRNA_infor))
    genome = SeqIO.parse(genome_file, 'fasta')
    all_donor_seq = []  # 供体剪切位点序列
    all_recep_seq = []  # 受体剪切位点序列
    all_non_donor_seq = []  # 非供体剪切位点序列
    all_non_recep_seq = []  # 非受体剪切位点序列
    seq_len = 50  # 序列长度
    j = 1  # 基因组中的序列数 1，2，3...n
    for RefSeq in genome:
        print('正在提取第%d条染色体信息。。。' % j)
        j = j + 1
        # 基因组文件中每条序列的序列号
        ID = RefSeq.id
        for single_mRNA in all_mRNA_infor:
            # 当gff文件中的序列号与之相匹配时 'local': 'NW_022983470.1' == 'NW_022983470.1'
            # ’+‘ 正向转录的序列
            if single_mRNA['local'] == ID and single_mRNA['direction'] == '+':
                # print(single_mRNA)
                # all_mRNA_information中的exon数据
                for i in range(len(single_mRNA['exon'])):
                    if (len(all_non_recep_seq) < train_num * 3) or (len(all_non_donor_seq) < train_num * 3):
                        # -----------外显子中的非剪切位点的AG GT 数据集------------
                        exon_seq = str(RefSeq.seq[single_mRNA['exon'][i][0] - 1:single_mRNA['exon'][i][1]]).upper()
                        AG_site = [i.start() for i in re.finditer('AG', exon_seq)]
                        GT_site = [i.start() for i in re.finditer('GT', exon_seq)]
                        for site in AG_site:
                            # print(site)
                            non_receptor = str(RefSeq.seq[single_mRNA['exon'][i][0] + site - (seq_len // 2):
                                                          single_mRNA['exon'][i][0] + site + (seq_len // 2)]).upper()
                            all_non_recep_seq.append(non_receptor)
                            # print('non_receptor:' + non_receptor)
                        for site in GT_site:
                            # print(site)
                            non_donor = str(RefSeq.seq[single_mRNA['exon'][i][0] + site - (seq_len // 2):
                                                       single_mRNA['exon'][i][0] + site + (seq_len // 2)]).upper()
                            all_non_donor_seq.append(non_donor)
                            # print('non_donor:' + non_donor)
                        # print('exon %d:' % (i+1) + RefSeq.seq[single_mRNA['exon'][i][0]-1:single_mRNA['exon'][i][1]])
                        # -----------外显子中的非剪切位点的AG GT 数据集------------

                    if i == 0:
                        # -----------外显子中的剪切位点数据集------------
                        donor = RefSeq.seq[single_mRNA['exon'][0][1] - ((seq_len // 2) - 1):
                                           single_mRNA['exon'][0][1] + ((seq_len // 2) + 1)]
                        all_donor_seq.append(str(donor).upper())
                        # print('donor\t' + donor)
                        # -----------外显子中的剪切位点数据集------------

                        # # -----------内含子中的非剪切位点的AG GT 数据集------------
                        if (len(all_non_recep_seq) < train_num * 3) or (
                                len(all_non_donor_seq) < train_num * 3):
                            intron = str(
                                RefSeq.seq[single_mRNA['exon'][i][1]:single_mRNA['exon'][i + 1][0] - 1]).upper()
                            # print('intron:' + intron[:50] + '...' + intron[len(intron) - 50:])
                            AG_site = [i.start() for i in re.finditer('AG', intron)]
                            GT_site = [i.start() for i in re.finditer('GT', intron)]
                            for site in AG_site:
                                # print(site)
                                non_receptor = str(RefSeq.seq[single_mRNA['exon'][i][1] + 1 + site - (seq_len // 2):
                                                              single_mRNA['exon'][i][1] + 1 + site + (
                                                                      seq_len // 2)]).upper()
                                all_non_recep_seq.append(non_receptor)
                                # print('non_receptor:' + non_receptor)
                            for site in GT_site:
                                # print(site)
                                non_donor = str(RefSeq.seq[single_mRNA['exon'][i][1] + 1 + site - (seq_len // 2):
                                                           single_mRNA['exon'][i][1] + 1 + site + (
                                                                   seq_len // 2)]).upper()
                                all_non_donor_seq.append(non_donor)
                                # print('non_donor:' + non_donor)
                        # # -----------内含子中的非剪切位点的AG GT 数据集------------

                    elif i < len(single_mRNA['exon']) - 1:
                        # # -----------内含子中的非剪切位点的AG GT 数据集------------
                        if (len(all_non_recep_seq) < train_num * 3) or (
                                len(all_non_donor_seq) < train_num * 3):
                            intron = str(
                                RefSeq.seq[single_mRNA['exon'][i][1]:single_mRNA['exon'][i + 1][0] - 1]).upper()
                            # print('intron:' + intron[:50] + '...' + intron[len(intron) - 50:])
                            AG_site = [i.start() for i in re.finditer('AG', intron)]
                            GT_site = [i.start() for i in re.finditer('GT', intron)]
                            for site in AG_site:
                                # print(site)
                                non_receptor = str(RefSeq.seq[single_mRNA['exon'][i][1] + 1 + site - (seq_len // 2):
                                                              single_mRNA['exon'][i][1] + 1 + site + (
                                                                      seq_len // 2)]).upper()
                                all_non_recep_seq.append(non_receptor)
                                # print('non_receptor:' + non_receptor)
                            for site in GT_site:
                                # print(site)
                                non_donor = str(RefSeq.seq[single_mRNA['exon'][i][1] + 1 + site - (seq_len // 2):
                                                           single_mRNA['exon'][i][1] + 1 + site + (
                                                                   seq_len // 2)]).upper()
                                all_non_donor_seq.append(non_donor)
                                # print('non_donor:' + non_donor)
                        # # -----------内含子中的非剪切位点的AG GT 数据集------------

                        # -----------外显子中的剪切位点数据集------------
                        receptor = RefSeq.seq[
                                   single_mRNA['exon'][i][0] - ((seq_len // 2) + 2):
                                   single_mRNA['exon'][i][0] + ((seq_len // 2) - 2)]
                        all_recep_seq.append(str(receptor).upper())
                        # print('recep\t' + receptor)
                        donor = RefSeq.seq[
                                single_mRNA['exon'][i][1] - ((seq_len // 2) - 1):
                                single_mRNA['exon'][i][1] + ((seq_len // 2) + 1)]
                        all_donor_seq.append(str(donor).upper())
                        # print('donor\t' + donor)
                        # -----------外显子中的剪切位点数据集------------
                    else:
                        # -----------外显子中的剪切位点数据集------------
                        receptor = RefSeq.seq[
                                   single_mRNA['exon'][-1][0] - ((seq_len // 2) + 2):
                                   single_mRNA['exon'][-1][0] + ((seq_len // 2) - 2)]
                        all_recep_seq.append(str(receptor).upper())
                        # print('recep\t' + receptor)
                        # -----------外显子中的剪切位点数据集------------
                # break

            # 反向转录的序列
            elif single_mRNA['local'] == ID and single_mRNA['direction'] == '-':
                # print(single_mRNA)
                # all_mRNA_information中的exon数据
                for i in range(len(single_mRNA['exon'])):
                    # # -----------外显子中的非剪切位点的AG GT 数据集------------
                    if (len(all_non_recep_seq) < train_num * 3) or (len(all_non_donor_seq) < train_num * 3):
                        exon_seq = str(RefSeq.seq[single_mRNA['exon'][i][0]:
                                                  single_mRNA['exon'][i][1] + 1].reverse_complement()).upper()
                        # print('exon %d:' % (i + 1) + exon_seq)
                        AG_site = [i.start() for i in re.finditer('AG', exon_seq)]
                        GT_site = [i.start() for i in re.finditer('GT', exon_seq)]
                        for site in AG_site:
                            # print(site)
                            non_receptor = str(RefSeq.seq[
                                               single_mRNA['exon'][i][1] - (site + (seq_len // 2)):
                                               single_mRNA['exon'][i][1] - (
                                                       site - (seq_len // 2))].reverse_complement()).upper()
                            all_non_recep_seq.append(non_receptor)
                            # print('non_receptor:' + non_receptor)
                        for site in GT_site:
                            # print(site)
                            non_donor = str(RefSeq.seq[
                                            single_mRNA['exon'][i][1] - (site + (seq_len // 2)):
                                            single_mRNA['exon'][i][1] - (
                                                    site - (seq_len // 2))].reverse_complement()).upper()
                            all_non_donor_seq.append(non_donor)
                            # print('non_donor:' + non_donor)
                    # # -----------外显子中的非剪切位点的AG GT 数据集------------

                    if i == 0:

                        # -----------外显子中的剪切位点数据集------------
                        donor = RefSeq.seq[
                                single_mRNA['exon'][0][0] - ((seq_len // 2) + 2):
                                single_mRNA['exon'][0][0] + ((seq_len // 2) - 2)].reverse_complement()
                        all_donor_seq.append(str(donor).upper())
                        # print('donor\t' + donor)
                        # -----------外显子中的剪切位点数据集------------

                        # # -----------内含子中的非剪切位点的AG GT 数据集------------
                        if (len(all_non_recep_seq) < train_num * 3) or (
                                len(all_non_donor_seq) < train_num * 3):
                            intron = str(RefSeq.seq[single_mRNA['exon'][i + 1][1]:
                                                    single_mRNA['exon'][i][0] - 1].reverse_complement()).upper()
                            if len(intron) < 20:
                                continue
                            # print('intron:' + intron[:50] + '...' + intron[len(intron) - 50:])
                            AG_site = [i.start() for i in re.finditer('AG', intron)]
                            GT_site = [i.start() for i in re.finditer('GT', intron)]
                            for site in AG_site:
                                # print(site)
                                non_receptor = str(RefSeq.seq[
                                                   single_mRNA['exon'][i][0] - 2 - (site + (seq_len // 2)):
                                                   single_mRNA['exon'][i][0] - 2 - (
                                                           site - (seq_len // 2))].reverse_complement()).upper()
                                all_non_recep_seq.append(non_receptor)
                                # print('non_receptor:' + non_receptor)
                            for site in GT_site:
                                # print(site)
                                non_donor = str(RefSeq.seq[
                                                single_mRNA['exon'][i][0] - 2 - (site + (seq_len // 2)):
                                                single_mRNA['exon'][i][0] - 2 - (
                                                        site - (seq_len // 2))].reverse_complement()).upper()
                                all_non_donor_seq.append(non_donor)
                                # print('non_donor:' + non_donor)
                        # # -----------内含子中的非剪切位点的AG GT 数据集------------
                        # break
                    elif i < len(single_mRNA['exon']) - 1:
                        # -----------外显子中的剪切位点数据集------------
                        receptor = RefSeq.seq[
                                   single_mRNA['exon'][i][1] - ((seq_len // 2) - 1):
                                   single_mRNA['exon'][i][1] + ((seq_len // 2) + 1)].reverse_complement()
                        # print('recep\t' + str(receptor).upper())
                        all_recep_seq.append(str(receptor).upper())
                        # print('recep\t' + receptor)
                        donor = RefSeq.seq[
                                single_mRNA['exon'][i][0] - ((seq_len // 2) + 2):
                                single_mRNA['exon'][i][0] + ((seq_len // 2) - 2)].reverse_complement()
                        all_donor_seq.append(str(donor).upper())
                        # print('donor\t' + donor)
                        # -----------外显子中的剪切位点数据集------------

                        # # -----------内含子中的非剪切位点的AG GT 数据集------------
                        if (len(all_non_recep_seq) < train_num * 3) or (
                                len(all_non_donor_seq) < train_num * 3):
                            # print(len(all_non_recep_seq), len(all_non_donor_seq))
                            intron = str(RefSeq.seq[single_mRNA['exon'][i + 1][1]:
                                                    single_mRNA['exon'][i][0] - 1].reverse_complement()).upper()
                            if len(intron) < 20:
                                continue
                            # print('intron:' + intron[:50] + '...' + intron[len(intron) - 50:])
                            AG_site = [i.start() for i in re.finditer('AG', intron)]
                            GT_site = [i.start() for i in re.finditer('GT', intron)]
                            for site in AG_site:
                                # print(site)
                                non_receptor = str(RefSeq.seq[
                                                   single_mRNA['exon'][i][0] - 2 - (site + (seq_len // 2)):
                                                   single_mRNA['exon'][i][0] - 2 - (
                                                           site - (seq_len // 2))].reverse_complement()).upper()
                                all_non_recep_seq.append(non_receptor)
                                # print('non_receptor:' + non_receptor)
                            for site in GT_site:
                                # print(site)
                                non_donor = str(RefSeq.seq[
                                                single_mRNA['exon'][i][0] - 2 - (site + (seq_len // 2)):
                                                single_mRNA['exon'][i][0] - 2 - (
                                                        site - (seq_len // 2))].reverse_complement()).upper()
                                all_non_donor_seq.append(non_donor)
                                # print('non_donor:' + non_donor)
                        # # -----------内含子中的非剪切位点的AG GT 数据集------------

                    else:
                        # -----------外显子中的剪切位点数据集------------
                        receptor = RefSeq.seq[
                                   single_mRNA['exon'][-1][1] - ((seq_len // 2) - 1):
                                   single_mRNA['exon'][-1][1] + ((seq_len // 2) + 1)].reverse_complement()
                        all_recep_seq.append(str(receptor).upper())
                        # print('recep\t' + receptor)
                        # -----------外显子中的剪切位点数据集------------
                # break

        # break
    print('共提取出%d条receptor_seqs' % len(all_recep_seq))
    print('共提取出%d条donor_seqs' % len(all_donor_seq))
    print('共提取出%d条non_recep_seqs' % len(all_non_recep_seq))
    print('共提取出%d条non_donor_seqs' % len(all_non_donor_seq))
    try:
        all_non_recep_seq = list(dict.fromkeys(all_non_recep_seq))
        all_non_donor_seq = list(dict.fromkeys(all_non_donor_seq))
        all_non_donor_seq = random.sample(all_non_donor_seq, train_num)
        all_non_recep_seq = random.sample(all_non_recep_seq, train_num)
        all_recep_seq = list(dict.fromkeys(all_recep_seq))
        all_donor_seq = list(dict.fromkeys(all_donor_seq))
        all_recep_seq = random.sample(all_recep_seq, train_num)
        all_donor_seq = random.sample(all_donor_seq, train_num)
    except Exception as e:
        pass
    print('共提取出%d条receptor_seqs' % len(all_recep_seq), all_recep_seq[:10])
    print('共提取出%d条donor_seqs' % len(all_donor_seq), all_donor_seq[:10])
    print('共提取出%d条non_recep_seqs' % len(all_non_recep_seq), all_non_recep_seq[:10])
    print('共提取出%d条non_donor_seqs' % len(all_non_donor_seq), all_non_donor_seq[:10])
    AG_seq_label = []
    GT_seq_label = []
    for seq in all_recep_seq:
        if ('N' not in seq) and ('W' not in seq) and ('Q' not in seq):
            AG_seq_label.append([seq, 1])
    for seq in all_non_recep_seq:
        if ('N' not in seq) and ('W' not in seq) and ('Q' not in seq):
            AG_seq_label.append([seq, 0])
    random.shuffle(AG_seq_label)
    AG_train_data = AG_seq_label[:int(len(AG_seq_label) * 0.8)]
    AG_test_data = AG_seq_label[int(len(AG_seq_label) * 0.8):]

    for seq in all_donor_seq:
        if ('N' not in seq) and ('W' not in seq) and ('Q' not in seq):
            GT_seq_label.append([seq, 1])
    for seq in all_non_donor_seq:
        if ('N' not in seq) and ('W' not in seq) and ('Q' not in seq):
            GT_seq_label.append([seq, 0])
    random.shuffle(GT_seq_label)
    GT_train_data = GT_seq_label[:int(len(GT_seq_label) * 0.8)]
    GT_test_data = GT_seq_label[int(len(GT_seq_label) * 0.8):]
    print('完成')
    return AG_train_data, AG_test_data, GT_train_data, GT_test_data



