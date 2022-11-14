# -*- coding: UTF-8 -*-
"""
__project_ = '剪切位点识别'
__file_name__ = 'VG_16'
__author__ = '鲨鱼辣椒'
__time__ = '2021/5/25 16:36'
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
import torch
import torch.nn as nn


class Net(nn.Module):  # 新建一个网络类，就是需要搭建的网络，必须继承PyTorch的nn.Module父类
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


class VGG16(nn.Module):
    def __init__(self):
        super(VGG16, self).__init__()
        # the vgg's layers
        # self.features = features
        # 13 conv + 3 FC
        cfg = [64, 64, 'M', 128, 128, 'M', 256, 256, 256, 'M', 512, 512, 512, 'M', 512, 512, 512, 'M']
        layers = []
        batch_norm = False
        in_channels = 1  # The input channel: 1
        for v in cfg:
            if v == 'M':
                layers += [nn.MaxPool2d(kernel_size=2, stride=2)]
            else:
                conv2d = nn.Conv2d(in_channels, v, kernel_size=3, padding=1)
                if batch_norm:
                    layers += [conv2d, nn.Batchnorm2d(v), nn.ReLU(inplace=True)]
                else:
                    layers += [conv2d, nn.ReLU(inplace=True)]
                in_channels = v
        # use the vgg layers to get the feature
        self.features = nn.Sequential(*layers)
        # 全局池化
        self.avgpool = nn.AdaptiveAvgPool2d((7, 7))
        # 决策层：分类层
        self.classifier = nn.Sequential(
            nn.Linear(512 * 7 * 7, 4096),
            nn.ReLU(True),
            nn.Dropout(),
            nn.Linear(4096, 4096),
            nn.ReLU(True),
            nn.Dropout(),
            nn.Linear(4096, 1000),
        )

        for m in self.modules():
            if isinstance(m, nn.Conv2d):
                nn.init.kaiming_normal_(m.weight, mode='fan_out', nonlinearity='relu')
                if m.bias is not None:
                    nn.init.constant_(m.bias, 0)
            elif isinstance(m, nn.BatchNorm2d):
                nn.init.constant_(m.weight, 1)
                nn.init.constant_(m.bias, 1)
            elif isinstance(m, nn.Linear):
                nn.init.normal_(m.weight, 0, 0.01)
                nn.init.constant_(m.bias, 0)

    def forward(self, x):
        x = self.features(x)
        x_fea = x
        x = self.avgpool(x)
        x_avg = x
        # batch*pixel
        x = x.view(x.size(0), -1)
        x = self.classifier(x)
        return x


# # ------------------------------------------------------------------------------
#
# class VGG(nn.Module):
#     """
#     VGG通用网络模型
#     输入features为网络的特征提取部分网络层列表
#     分类数为 1000
#     """
#
#     def __init__(self, features, num_classes=1000, init_weights=True):
#         super(VGG, self).__init__()
#
#         # 特征提取部分
#         self.features = features
#
#         # 自适应平均池化，特征图池化到 7×7 大小
#         self.avgpool = nn.AdaptiveAvgPool2d((7, 7))
#
#         # 分类部分
#         self.classifier = nn.Sequential(
#             nn.Linear(512 * 7 * 7, 4096),  # 512*7*7 --> 4096
#             nn.ReLU(True),
#             nn.Dropout(),
#             nn.Linear(4096, 4096),  # 4096 --> 4096
#             nn.ReLU(True),
#             nn.Dropout(),
#             nn.Linear(4096, num_classes),  # 4096 --> 1000
#         )
#
#         # 权重初始化
#         if init_weights:
#             self._initialize_weights()
#
#     def forward(self, x):
#
#         # 特征提取
#         x = self.features(x)
#         # 自适应平均池化
#         x = self.avgpool(x)
#         # 特征图展平成向量
#         x = torch.flatten(x, 1)
#         # 分类器分类输出
#         x = self.classifier(x)
#         return x
#
#     def _initialize_weights(self):
#         '''
#         权重初始化
#         '''
#         for m in self.modules():
#             if isinstance(m, nn.Conv2d):
#                 # 卷积层使用 kaimming 初始化
#                 nn.init.kaiming_normal_(
#                     m.weight, mode='fan_out', nonlinearity='relu')
#                 # 偏置初始化为0
#                 if m.bias is not None:
#                     nn.init.constant_(m.bias, 0)
#             # 批归一化层权重初始化为1
#             elif isinstance(m, nn.BatchNorm2d):
#                 nn.init.constant_(m.weight, 1)
#                 nn.init.constant_(m.bias, 0)
#             # 全连接层权重初始化
#             elif isinstance(m, nn.Linear):
#                 nn.init.normal_(m.weight, 0, 0.01)
#                 nn.init.constant_(m.bias, 0)
#
#
# # ------------------------------------------------------------------------------
# def make_layers(cfg, batch_norm=False):
#     '''
#     根据配置表，返回模型层列表
#     '''
#     layers = []  # 层列表初始化
#
#     in_channels = 1  # 输入1通道
#
#     # 遍历配置列表
#     for v in cfg:
#         if v == 'M':  # 添加池化层
#             layers += [nn.MaxPool2d(kernel_size=2, stride=2)]
#         else:  # 添加卷积层
#
#             # 3×3 卷积
#             conv2d = nn.Conv2d(in_channels, v, kernel_size=3, padding=1)
#
#             # 卷积-->批归一化（可选）--> ReLU激活
#             if batch_norm:
#                 layers += [conv2d, nn.BatchNorm2d(v), nn.ReLU(inplace=True)]
#             else:
#                 layers += [conv2d, nn.ReLU(inplace=True)]
#
#             # 通道数方面，下一层输入即为本层输出
#             in_channels = v
#
#     # 以sequencial类型返回模型层列表
#     return nn.Sequential(*layers)
#
#
# # 网络参数配置表
# '''
# 数字代表通道数，如 64 表示输出 64 通道特征图，对应于论文中的 Conv3-64;
# M 代表最大池化操作，对应于论文中的 maxpool
# A-LRN使用了局部归一化响应，C网络存在1×1卷积，这两个网络比较特殊，所以排除在配置表中
# '''
# cfgs = {
#     'vgg11': [64, 'M', 128, 'M', 256, 256, 'M', 512, 512, 'M', 512, 512, 'M'],
#     'vgg13': [64, 64, 'M', 128, 128, 'M', 256, 256, 'M', 512, 512, 'M', 512, 512, 'M'],
#     'vgg16': [64, 64, 'M', 128, 128, 'M', 256, 256, 256, 'M', 512, 512, 512, 'M', 512, 512, 512, 'M'],
#     'vgg19': [64, 64, 'M', 128, 128, 'M', 256, 256, 256, 256, 'M', 512, 512, 512, 512, 'M', 512, 512, 512, 512, 'M'],
# }
#
#
# # ------------------------------------------------------------------------------
#
# def vgg(model_name="vgg16", **kwargs):  # 双星号(**)将参数以字典的形式导入
#     try:
#         cfg = cfgs[model_name]
#     except:
#         print("Warning: model number {} not in cfgs dict!".format(model_name))
#         exit(-1)
#     model = VGG(make_layers(cfg), **kwargs)
#     return model
