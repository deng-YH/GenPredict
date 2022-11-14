# -*- coding: UTF-8 -*-
"""
__project_ = 'GenePredict'
__file_name__ = 'cnn'
__author__ = 'Dyh'
__time__ = '2021/11/13 13:52'
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

import torch
import torch.nn.functional as F
import torch.nn as nn
import torch.utils.data as data


# import numpy as np
# from torch.utils.data import DataLoader as DataLoader
#
# from Bio import SeqIO


# ---------------- k-mer one-hot编码 ---------------- #
def seq_to_OneHot(seq: str):
    # ---------得到ATGC的所有的3联组合列表----------
    labels_K_mer = []
    ATGC = ['A', 'T', 'G', 'C']
    for i in ATGC:
        for j in ATGC:
            for k in ATGC:
                labels_K_mer.append(str(i) + str(j) + str(k))

    # ----------将序列转换成K_mer编码 K = 3---------
    k = 3
    list_comb = []
    for index in range(len(seq)):
        t = seq[index:index + k]
        if (len(t)) == k:
            list_comb.append(t)

    # -------K_mer编码转换成one—hot编码------------
    data = list_comb
    ind_to_char = labels_K_mer
    # 定义字符到整数的映射
    char_to_int = dict((c, i) for i, c in enumerate(ind_to_char))  # 枚举
    int_to_char = dict((i, c) for i, c in enumerate(ind_to_char))
    # 整数编码
    integer_encoded = [char_to_int[char] for char in data]
    # one-hot编码
    one_hot_encoded = list()
    for value in integer_encoded:
        letter = [0 for _ in range(len(ind_to_char))]  # one-hot编码长度
        letter[value] = 1
        one_hot_encoded.append(letter)
    return one_hot_encoded


# ---------------- 获取剪切位点序列 ---------------- #
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


# ---------------- 网络类 ---------------- #
class Net(nn.Module):
    """
    网络类，继承PyTorch的nn.Module父类
    """

    def __init__(self):  # 构造函数，用于设定网络层
        super(Net, self).__init__()  # 标准语句
        self.Conv = torch.nn.Sequential(
            torch.nn.Conv2d(1, 64, kernel_size=3, stride=1, padding=1),
            torch.nn.ReLU(),
            torch.nn.Conv2d(64, 64, kernel_size=3, stride=1, padding=1),
            torch.nn.ReLU(),
            torch.nn.MaxPool2d(kernel_size=2, stride=2),

            torch.nn.Conv2d(64, 128, kernel_size=3, stride=1, padding=1),
            torch.nn.ReLU(),
            torch.nn.Conv2d(128, 128, kernel_size=3, stride=1, padding=1),
            torch.nn.ReLU(),
            torch.nn.MaxPool2d(kernel_size=2, stride=2),

            torch.nn.Conv2d(128, 256, kernel_size=3, stride=1, padding=1),
            torch.nn.ReLU(),
            torch.nn.Conv2d(256, 256, kernel_size=3, stride=1, padding=1),
            torch.nn.ReLU(),
            torch.nn.Conv2d(256, 256, kernel_size=3, stride=1, padding=1),
            torch.nn.ReLU(),
            torch.nn.MaxPool2d(kernel_size=2, stride=2),

            torch.nn.Conv2d(256, 512, kernel_size=3, stride=1, padding=1),
            torch.nn.ReLU(),
            torch.nn.Conv2d(512, 512, kernel_size=3, stride=1, padding=1),
            torch.nn.ReLU(),
            torch.nn.Conv2d(512, 512, kernel_size=3, stride=1, padding=1),
            torch.nn.ReLU(),
            torch.nn.MaxPool2d(kernel_size=2, stride=2),

        )

        self.Classes = torch.nn.Sequential(
            torch.nn.Linear(3 * 4 * 512, 1024),
            torch.nn.ReLU(),
            torch.nn.Dropout(p=0.5),
            torch.nn.Linear(1024, 1024),
            torch.nn.ReLU(),
            torch.nn.Dropout(p=0.5),
            torch.nn.Linear(1024, 2)

        )

    def forward(self, x):  # 重写父类forward方法，即前向计算，通过该方法获取网络输入数据后的输出值
        x = self.Conv(x)
        x = x.view(-1, 3 * 4 * 512)
        y = self.Classes(x)
        return y


# ---------------- 数据集类 ---------------- #
class SeqsData(data.Dataset):  # 新建一个数据集类，并且需要继承PyTorch中的data.Dataset父类
    def __init__(self, mode, dataset):  # 默认构造函数，传入数据集类别（训练或测试），以及数据集路径
        self.mode = mode
        self.list_seq = []  # 新建一个seq list，用于存放所有的seq文本 例：'GTCTCGTGTCCTGTGACTTCCTCAGGCCTCCACCAACA'
        self.list_label = []  # 新建一个label list，用于存放seq对应是否为真剪切位点的标签，其中数值0表示假，1表示真
        self.data_size = 0  # 记录数据集大小
        # self.transform = dataTransform      # 转换关系

        if self.mode == 'train':  # 训练集模式下，需要提取图片的路径和标签
            dataset_dir = r'./data/train/' + dataset  # 训练集路径在"dir"/train/
            with open(dataset_dir) as f:
                seq_label = json.load(f)  # 打开训练集
            for seq in seq_label:
                self.list_seq.append(seq[0])  # 序列读取
                self.data_size += 1  # 数据集增1
                self.list_label.append(seq[1])  # label读取 与序列一一对应
        elif self.mode == 'test':  # 测试集模式下，只需要提取路径就行
            dataset_dir = r'./data/test/' + dataset  # 测试集路径在"dir"/test/
            with open(dataset_dir) as f:
                seq_label = json.load(f)  # 打开测试集
            for seq in seq_label:
                self.list_seq.append(seq[0])  # 序列读取
                self.data_size += 1  # 数据集增1
                self.list_label.append(2)  # label读取 与序列一一对应
        else:
            print('Undefined Dataset!')

    def __getitem__(self, item):  # 重载data.Dataset父类方法，获取数据集中数据内容
        if self.mode == 'train':  # 训练集模式下需要读取数据集的seq和label
            seq = seq_to_OneHot(self.list_seq[item])  # 加载序列为one-hot编码
            label = self.list_label[item]  # 获取seq对应的label
            seq = torch.tensor(seq, dtype=torch.float)
            return torch.unsqueeze(seq, dim=0).float(), torch.LongTensor([label])  # 将image和label转换成PyTorch形式并返回
        elif self.mode == 'test':  # 测试集只需读取image
            seq = seq_to_OneHot(self.list_seq[item])  # 加载序列
            seq = torch.tensor(seq, dtype=torch.float)
            return torch.unsqueeze(seq, dim=0).float()  # 将image转换成PyTorch形式并返回
        else:
            print('None')

    def __len__(self):
        return self.data_size  # 返回数据集大小


class CNN_model(object):
    os.environ['KMP_DUPLICATE_LIB_OK'] = 'TRUE'

    def __init__(self):
        self.network = Net()  # 实例化网络
        self.network = nn.DataParallel(self.network)
        # self.device = torch.device('cuda:0' if torch.cuda.is_available() else 'cpu')
        # self.net = Net()

    def ini(self, flag, model_file):
        # device = torch.device('cuda:0' if torch.cuda.is_available() else 'cpu')
        # self.network.load_state_dict(torch.load(model_file))
        # self.network = self.network.to(device)
        # # self.network.eval()
        if flag == 1:
            self.network.load_state_dict(
                torch.load(model_file, map_location=lambda storage, loc: storage))  # 加载训练好的模型参数
        else:
            self.network.load_state_dict(torch.load(model_file))
        self.network.eval()  # 设定为评估模式，即计算过程中不要dropout
        # self.network.eval()  # 设定为评估模式，即计算过程中不要dropout

    def test1(self, seq: str, AG_MODEL, GT_MODEL):
        device = torch.device('cuda:0' if torch.cuda.is_available() else 'cpu')
        net = Net()
        net.load_state_dict(torch.load(AG_MODEL))
        net = net.to(device)
        seqs = []  # seqs
        AG_seqs_data = []  # seqs data
        AG, GT = get_AG_GT_site(seq, 50)
        for seq in AG:
            seqs.append(seq)
            # labels.append(seq[1])
            seq = seq_to_OneHot(seq)
            AG_seqs_data.append(torch.unsqueeze(torch.tensor(seq, dtype=torch.float), dim=0))
        seqs_data = torch.stack(AG_seqs_data)  # tensor list合成一个4D tensor

        # calculation
        AG_sites = []
        out = net(seqs_data)  # 对每个图像进行网络计算
        out = F.softmax(out, dim=1)  # 输出概率化
        out = out.data.cpu().numpy()  # 转成numpy数据

        for idx in range(len(out)):
            if out[idx, 1] > 0.9:
                # print(seqs[idx] + '是:{:.1%},否:{:.1%}'.format(out[idx, 1], out[idx, 0]))
                AG_sites.append(seqs[idx])

        # GT
        GT_seqs = []
        GT_seqs_data = []  # seqs data
        device = torch.device('cuda:0' if torch.cuda.is_available() else 'cpu')
        net = Net()
        net.load_state_dict(torch.load(GT_MODEL))
        net = net.to(device)
        for seq in GT:
            GT_seqs.append(seq)
            # labels.append(seq[1])
            seq = seq_to_OneHot(seq)
            GT_seqs_data.append(torch.unsqueeze(torch.tensor(seq, dtype=torch.float), dim=0))
        seqs_data = torch.stack(GT_seqs_data)  # tensor list合成一个4D tensor

        # calculation

        GT_sites = []
        out = net(seqs_data)  # 对每个图像进行网络计算
        out = F.softmax(out, dim=1)  # 输出概率化
        out = out.data.cpu().numpy()  # 转成numpy数据

        for idx in range(len(out)):
            if out[idx, 1] > 0.9:
                # print(seqs[idx] + '是:{:.1%},否:{:.1%}'.format(out[idx, 1], out[idx, 0]))
                GT_sites.append(GT_seqs[idx])
        return AG_sites, GT_sites

    def test(self, seq: str, AG_MODEL, GT_MODEL):
        # self.network.cuda()          # 送入GPU，利用GPU计算
        seqs = []  # seqs
        AG_seqs_data = []  # seqs data
        AG, GT = get_AG_GT_site(seq, 50)

        # AG
        self.ini(0, AG_MODEL)
        for seq in AG:
            seqs.append(seq)
            # labels.append(seq[1])
            seq = seq_to_OneHot(seq)
            AG_seqs_data.append(torch.unsqueeze(torch.tensor(seq, dtype=torch.float), dim=0))
        seqs_data = torch.stack(AG_seqs_data)  # tensor list合成一个4D tensor

        # calculation
        AG_sites = []
        out = self.network(seqs_data)  # 对每个图像进行网络计算
        out = F.softmax(out, dim=1)  # 输出概率化
        out = out.data.cpu().numpy()  # 转成numpy数据

        for idx in range(len(out)):
            if out[idx, 1] > 0.9:
                # print(seqs[idx] + '是:{:.1%},否:{:.1%}'.format(out[idx, 1], out[idx, 0]))
                AG_sites.append(seqs[idx])

        # GT
        GT_seqs = []
        GT_seqs_data = []  # seqs data
        self.ini(1, GT_MODEL)
        for seq in GT:
            GT_seqs.append(seq)
            # labels.append(seq[1])
            seq = seq_to_OneHot(seq)
            GT_seqs_data.append(torch.unsqueeze(torch.tensor(seq, dtype=torch.float), dim=0))
        seqs_data = torch.stack(GT_seqs_data)  # tensor list合成一个4D tensor

        # calculation

        GT_sites = []
        out = self.network(seqs_data)  # 对每个图像进行网络计算
        out = F.softmax(out, dim=1)  # 输出概率化
        out = out.data.cpu().numpy()  # 转成numpy数据

        for idx in range(len(out)):
            if out[idx, 1] > 0.9:
                # print(seqs[idx] + '是:{:.1%},否:{:.1%}'.format(out[idx, 1], out[idx, 0]))
                GT_sites.append(GT_seqs[idx])
        return AG_sites, GT_sites


if __name__ == '__main__':
    cnn = CNN_model()
