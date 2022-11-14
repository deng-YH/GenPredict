# -*- coding: UTF-8 -*-
"""
__project_ = '卷积神经网络构建'
__file_name__ = 'get_data2'
__author__ = '鲨鱼辣椒'
__time__ = '2021/5/24 19:42'
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
import os

import torch.utils.data as data
import json
import numpy as np
import torch
from torch.utils.data import DataLoader as DataLoader


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


# seq = torch.tensor(seq_to_OneHot('CTTTATGAGTAGTCAGGGTTTCTCAGGACCCTGTGGTAAGTAAGGCCTGC'), dtype=torch.float)
# seq = torch.unsqueeze(seq, dim=0).float()
# print(seq.shape)

class SeqsData(data.Dataset):  # 新建一个数据集类，并且需要继承PyTorch中的data.Dataset父类
    script_path = os.path.split(os.path.realpath(__file__))[0]

    def __init__(self, mode, dataset):  # 默认构造函数，传入数据集类别（训练或测试），以及数据集路径
        self.mode = mode
        self.list_seq = []  # 新建一个seq list，用于存放所有的seq文本 例：'GTCTCGTGTCCTGTGACTTCCTCAGGCCTCCACCAACA'
        self.list_label = []  # 新建一个label list，用于存放seq对应是否为真剪切位点的标签，其中数值0表示假，1表示真
        self.data_size = 0  # 记录数据集大小
        # self.transform = dataTransform      # 转换关系
        data_path = self.script_path.replace(r'src\CNN', '')
        if self.mode == 'train':  # 训练集模式下，需要提取图片的路径和标签
            dataset_dir = os.path.join(data_path, r'data/train/' + dataset)  # 训练集路径在"dir"/train/
            with open(dataset_dir) as f:
                seq_label = json.load(f)  # 打开训练集
            for seq in seq_label:
                self.list_seq.append(seq[0])  # 序列读取
                self.data_size += 1  # 数据集增1
                self.list_label.append(seq[1])  # label读取 与序列一一对应
        elif self.mode == 'test':  # 测试集模式下，只需要提取图片路径就行
            dataset_dir = os.path.join(data_path, r'data/test/' + dataset)  # 测试集路径在"dir"/test/
            with open(dataset_dir) as f:
                seq_label = json.load(f)  # 打开测试集
            for seq in seq_label:
                self.list_seq.append(seq[0])  # 序列读取
                self.data_size += 1  # 数据集增1
                self.list_label.append(seq[1])  # label读取 与序列一一对应
        elif self.mode == 'predict':  # 测试集模式下，只需要提取图片路径就行
            dataset_dir = os.path.join(data_path, r'data/test/' + dataset)  # 测试集路径在"dir"/test/
            with open(dataset_dir) as f:
                seq_label = json.load(f)  # 打开测试集
            for seq in seq_label:
                self.list_seq.append(seq[0])  # 序列读取
                self.data_size += 1  # 数据集增1
                # self.list_label.append(seq[1])  # label读取 与序列一一对应
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
            label = self.list_label[item]
            seq = torch.tensor(seq, dtype=torch.float)
            return torch.unsqueeze(seq, dim=0).float(), torch.LongTensor([label])  # 将image转换成PyTorch形式并返回
        elif self.mode == 'predict':  # 测试集只需读取image
            seq = seq_to_OneHot(self.list_seq[item])  # 加载序列
            label = self.list_seq[item]
            seq = torch.tensor(seq, dtype=torch.float)
            return torch.unsqueeze(seq, dim=0).float(), label  # 将image转换成PyTorch形式并返回
        else:
            print('None')

    def __len__(self):
        return self.data_size  # 返回数据集大小

# dataset_dir = './data/'
# datafile = SeqsData('train', dataset_dir)  # 实例化一个数据集
# dataloader = DataLoader(datafile, batch_size=512, shuffle=True, num_workers=4,
#                         drop_last=True)  # 用PyTorch的DataLoader类封装，实现数据集顺序打乱，多线程读取，一次取多个数据等效果
#
# print('Dataset loaded! length of train set is {0}'.format(len(datafile)))
