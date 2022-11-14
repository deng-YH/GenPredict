# -*- coding: UTF-8 -*-
"""
__project_ = 'blast_X32_GUI'
__file_name__ = 'main_window'
__author__ = '鲨鱼辣椒'
__time__ = '2021/5/29 11:59'
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
import sys
import subprocess
from PyQt5.QtCore import Qt, QCoreApplication, pyqtSlot
from PyQt5.QtWidgets import QMainWindow, QApplication, QSplitter, QHBoxLayout, QFileSystemModel

import GenePredict
from uifiles import MainWindow
from src import diamond, local_blast, blast_make_db, blastdbcmd, translate, nbci_entry


class MainWindowSRC(QMainWindow, MainWindow.Ui_MainWindow):
    root = os.path.split(os.path.realpath(__file__))[0]

    def __init__(self):
        super(MainWindowSRC, self).__init__()
        # 加载界面
        self.setupUi(self)  # 加载父类UI

        self.GenePredict = GenePredict.GenePredict()  # 实例GenePredict对象
        self.groupBox_layout = QHBoxLayout()  # groupBox中的layout

        self.model = QFileSystemModel()  # 这里得到目录结构
        self.tree_view_get_file()  # treeView 获取当前目录下所有文件

        self.splitter1 = QSplitter(Qt.Horizontal)  # QSplitter允许用户通过拖动子窗口之间的边界来控制子窗口小部件的大小
        self.splitter1.addWidget(self.treeView)  # splitter1添加treeView控件
        self.splitter1.addWidget(self.GenePredict)  # 将GenePredict添加到splitter1中

        self.windowLayout = QHBoxLayout()  # 水平布局
        self.windowLayout.addWidget(self.splitter1)  # 添加splitter1 控件
        self.centralwidget.setLayout(self.windowLayout)  # 设置主窗体为水平布局

        # 实例对象
        self.frame_diamond = diamond.Diamond()
        self.frame_local_blast = local_blast.LocalBlast()
        self.frame_make_blast_db = blast_make_db.MakeBlastDB()
        self.frame_blastdbcmd = blastdbcmd.Blastdbcmd()
        self.frame_translate = translate.C_translate()
        self.frame_ncbi_entry = nbci_entry.NCBIEntry()

        # 信号和槽
        self.treeView.doubleClicked.connect(self.treeView_doubleClicked)  # treeView连接槽函数 鼠标单击事件
        self.actionDiamond.triggered.connect(lambda: self.show_frame(self.frame_diamond))
        self.actionBlast.triggered.connect(lambda: self.show_frame(self.frame_local_blast))
        self.action_makeblastdb.triggered.connect(lambda: self.show_frame(self.frame_make_blast_db))
        self.action_Blastdbcmd.triggered.connect(lambda: self.show_frame(self.frame_blastdbcmd))
        self.actionTranslate_tool.triggered.connect(lambda: self.show_frame(self.frame_translate))
        self.actionEntry.triggered.connect(lambda: self.show_frame(self.frame_ncbi_entry))

    def tree_view_get_file(self):
        self.model.setRootPath(self.root)  # 设置文件结构根目录
        self.treeView.setModel(self.model)  # treeView设置目录模型
        self.treeView.setRootIndex(self.model.index(self.root))  # 从设置目录root加载目录信息

        self.treeView.setColumnHidden(1, True)  # tree_view隐藏SIZE
        self.treeView.setColumnHidden(2, True)  # tree_view隐藏file type
        self.treeView.setColumnHidden(3, True)  # # tree_view隐藏time

    # 槽函数
    def treeView_doubleClicked(self, Qmodelidx):
        """
        单击返回路径
        :param Qmodelidx: 选中的目录idx
        :return:  路径
        """
        filePath = self.model.filePath(Qmodelidx)  # 返回目录信息
        os.system(f'start {filePath}')  # 执行目录

    def show_frame(self, frame):
        frame.show()

    def closeEvent(self, event):
        """
        删除tempfile目录下的所有文件或文件夹
        """
        del_list = os.listdir('./tempfile')
        for f in del_list:
            file_path = os.path.join(r'./tempfile', f)
            if os.path.isfile(file_path):
                os.remove(file_path)


if __name__ == '__main__':
    QCoreApplication.setAttribute(Qt.AA_EnableHighDpiScaling)
    app = QApplication(sys.argv)
    MainWindow = MainWindowSRC()
    MainWindow.show()
    sys.exit(app.exec_())
