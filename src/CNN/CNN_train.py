# -*- coding: UTF-8 -*-
"""
__project_ = '卷积神经网络构建'
__file_name__ = 'train'
__author__ = '鲨鱼辣椒'
__time__ = '2021/5/24 10:33'
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

from PyQt5.QtCore import pyqtSignal, QObject
from .CNN_get_data import SeqsData
from torch.utils.data import DataLoader as DataLoader
from .CNN_VGG_16 import VGG16, Net
import torch


# os.environ['CUDA_VISIBLE_DEVICES'] = '1'


class CNN_model(QObject):
    single_out = pyqtSignal(str)
    single_finish = pyqtSignal()
    script_path = os.path.split(os.path.realpath(__file__))[0]
    print(script_path)

    def __init__(self):
        super(CNN_model, self).__init__()
        self.datasets = ['AG.json', 'GT.json']  # 数据集路径
        self.model_path = ''

    def train2(self, batch_size=64, workers=4, lr=0.0001, nepoch=10, save_name='model'):
        # CPU 或者 GPU
        device = torch.device('cuda:0' if torch.cuda.is_available() else 'cpu')
        print(device)
        model = Net().to(device)  # 实例化一个网络

        # 优化器设置
        criterion = torch.nn.CrossEntropyLoss()  # 定义loss计算方法，cross entropy，交叉熵，可以理解为两者数值越接近其值越小
        optimizer = torch.optim.Adam(model.parameters(), lr=lr)  # 实例化一个优化器，即调整网络参数，优化方式为adam方法

        for dataset in self.datasets:
            # 加载数据集
            train_dataset = SeqsData('train', dataset)
            test_dataset = SeqsData('test', dataset)
            self.single_out.emit('Dataset loaded! length of train set is {0}'.format(len(train_dataset)))
            self.single_out.emit('Dataset loaded! length of test set is {0}'.format(len(test_dataset)))
            train_loader = DataLoader(train_dataset, batch_size=batch_size, shuffle=True, num_workers=workers)

            test_loader = DataLoader(test_dataset, batch_size=batch_size, shuffle=True, num_workers=workers)
            save_path = os.path.join(self.model_path, '{}.pth'.format(dataset.replace('.json', '')))

            for epoch in range(nepoch):
                total_step = len(train_loader)
                train_epoch_loss = 0
                model.train()

                for i, (seqs, labels) in enumerate(train_loader):
                    # 梯度清零
                    optimizer.zero_grad()
                    # 加载序列和标签
                    seqs = seqs.to(device)
                    labels = labels.to(device)
                    # 前向计算
                    output = model(seqs)
                    loss = criterion(output, labels.squeeze())
                    # 反向传播与优化
                    loss.backward()
                    optimizer.step()
                    # 累加每代中所有步数的loss
                    train_epoch_loss += loss.item()
                    # 打印部分结果
                    if (i + 1) % 2 == 0:
                        self.single_out.emit('Epoch [{}/{}], Step [{}/{}], Loss: {:.5f}'
                                             .format(epoch + 1, nepoch, i + 1, total_step, loss.item()))
                    if (i + 1) == total_step:
                        epoch_eva_loss = train_epoch_loss / total_step
                        # evaloss.append(epoch_eva_loss)
                        self.single_out.emit('Epoch_eva loss is : {:.5f}'.format(epoch_eva_loss))
                # test过程

                model.eval()
                with torch.no_grad():
                    correct = 0
                    total = 0
                    # print(test_loader)
                    # help(test_loader)
                    for seqs, labels in test_loader:
                        seqs = seqs.to(device)
                        labels = labels.to(device)
                        output = model(seqs)
                        _, predicted = torch.max(output.data, 1)

                        predicted = predicted[:, None]
                        # print(predicted, labels)
                        total += labels.size(0)
                        correct += (predicted == labels).sum().item()
                    # acc.append(100 * (correct / total))
                    # print(total, correct)
                    self.single_out.emit('Test Accuracy  {} %'.format(100 * (correct / total)))
            torch.save(obj=model.state_dict(), f=save_path)
        self.single_out.emit(
            '\n*************************\n'
            'Model training completed!\n'
            '*************************'
        )
        self.single_finish.emit()

    def test(self):
        # CPU 或者 GPU
        device = torch.device('cuda:0' if torch.cuda.is_available() else 'cpu')
        net = Net()
        net.load_state_dict(torch.load(r'D:\ddd\GenePredict\model\fish\AG.pth'))
        net = net.to(device)
        test_dataset = SeqsData('predict', 'AG.json')
        test_loader = DataLoader(test_dataset, batch_size=64, shuffle=True, num_workers=4)
        with torch.no_grad():
            correct = 0
            total = 0
            # print(test_loader)
            # help(test_loader)
            for seqs, labels in test_loader:
                seqs = seqs.to(device)
                # labels = labels.to(device)
                output = net(seqs)
                _, predicted = torch.max(output.data, 1)
                print(predicted)
                print(labels)
                # predicted = predicted[:, None]
                # print(predicted, labels)
                # total += labels.size(0)
                # correct += (predicted == labels).sum().item()
            # acc.append(100 * (correct / total))
            # print(total, correct)
            # self.single_out.emit('Test Accuracy  {} %'.format(100 * (correct / total)))
            # print('Test Accuracy  {} %'.format(100 * (correct / total)))


if __name__ == '__main__':
    cnn_model = CNN_model()
    cnn_model.test()
