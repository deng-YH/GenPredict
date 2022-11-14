# -*- coding: UTF-8 -*-
"""
__project_ = '卷积神经网络构建'
__file_name__ = 'test'
__author__ = '鲨鱼辣椒'
__time__ = '2021/5/24 13:18'
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
import json
import os
import re

from Bio import SeqIO
from CNN_VGG_16 import Net
import torch
import torch.nn.functional as F
import torch.nn as nn
import CNN_get_data

os.environ['KMP_DUPLICATE_LIB_OK'] = 'TRUE'
dataset_dir = './data/test/'  # 数据集路径
model_file = './model/model.pth'  # 模型保存路径


def get_AG_GT_site(seq: str, seq_len: int):
    all_AG_seq = []
    all_GT_seq = []
    # -----------seq的AG GT 数据集------------
    seq = seq.upper()
    AG_site = [i.start() for i in re.finditer('AG', seq)]
    GT_site = [i.start() for i in re.finditer('GT', seq)]
    for site in AG_site:
        # print(site)
        AG_seq = seq[site - (seq_len // 2) + 1: site + (seq_len // 2) + 1].upper()
        if len(AG_seq) == seq_len:
            all_AG_seq.append(AG_seq)

    for site in GT_site:
        # print(site)
        GT_seq = seq[site - (seq_len // 2) + 1: site + (seq_len // 2) + 1].upper()
        if len(GT_seq) == seq_len:
            all_GT_seq.append(GT_seq)

    return all_AG_seq, all_GT_seq
    # -----------seq的AG GT 数据集------------


# new version
def test():
    # setting model
    model = Net()  # 实例化一个网络
    # model.cuda()                                        # 送入GPU，利用GPU计算
    model = nn.DataParallel(model)
    model.load_state_dict(torch.load(model_file))  # 加载训练好的模型参数
    model.eval()  # 设定为评估模式，即计算过程中不要dropout

    # get data
    seqs = []  # img
    labels = []
    seqs_data = []  # img data
    dir = './data/test/AG.json'  # 测试集路径
    with open(dir) as f:
        seq_label = json.load(f)  # 打开测试集
    # handle = SeqIO.read(r'./1.fa', 'fasta')
    # seq = handle.seq
    # AG, GT = get_AG_GT_site(str(seq), 50)
    # for seq in AG:
    #     seqs.append(seq)
    #     # labels.append(seq[1])
    #     seq = CNN_get_data.seq_to_OneHot(seq)
    #     seqs_data.append(torch.unsqueeze(torch.tensor(seq, dtype=torch.float), dim=0))

    for seq in seq_label[:10000]:
        seqs.append(seq[0])
        labels.append(seq[1])
        seq = CNN_get_data.seq_to_OneHot(seq[0])
        seqs_data.append(torch.unsqueeze(torch.tensor(seq, dtype=torch.float), dim=0))

    seqs_data = torch.stack(seqs_data)  # tensor list合成一个4D tensor

    # calculation
    out = model(seqs_data)  # 对每个图像进行网络计算
    out = F.softmax(out, dim=1)  # 输出概率化
    out = out.data.cpu().numpy()  # 转成numpy数据
    print(len(out), len(labels))
    # pring results         显示结果
    right_num = 0
    for idx in range(len(out)):
        # plt.figure()
        if out[idx, 0] > out[idx, 1]:
            if labels[idx] == 0:
                right_num += 1
            print(seqs[idx] + '否:{:.1%},是:{:.1%}'.format(out[idx, 0], out[idx, 1]))
            print('', end='')
        else:
            if labels[idx] == 1:
                right_num += 1

            print(seqs[idx] + '是:{:.1%},否:{:.1%}'.format(out[idx, 1], out[idx, 0]))
        # plt.imshow(imgs[idx])
    print(right_num/len(out))


if __name__ == '__main__':
    test()
