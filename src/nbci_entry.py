# -*- coding: UTF-8 -*-
"""
__project_ = 'blast_X32_GUI'
__file_name__ = 'entyr库'
__author__ = '鲨鱼辣椒'
__time__ = '2021/6/4 10:30'
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
import re
import sys

from Bio import Entrez
from PyQt5.QtCore import pyqtSignal
from threading import Thread
from PyQt5 import QtCore
from PyQt5.QtGui import QIntValidator
from uifiles import ui_ncbi_entry
from PyQt5.QtWidgets import QApplication, QFileDialog, QMessageBox, QWidget

Accession = re.compile('([a-zA-z]+[\D\W]?[0-9]{3,}\.?[0-9]+)')     # re匹配NCBI数据库中的Accession登录号


class NCBIEntry(QWidget, ui_ncbi_entry.Ui_NCBI_entry):
    single_efetch = pyqtSignal([str, int, int], [str])  # 定义完成信号
    single_warning = pyqtSignal(str)  # error信号

    def __init__(self, parent=None):
        super(NCBIEntry, self).__init__(parent)
        # 加载界面
        self.setupUi(self)
        self.title = 'NCBI Entry'
        self.email = ''
        self.batch_size = 5
        self.id_list = []
        Intvalidator = QIntValidator(self)  # 实例化整型验证器
        self.lineEdit_2.setValidator(Intvalidator)  # 设置验证器

        # 信号和槽
        self.pushButton_get.clicked.connect(self.pushButton_get_clicked)
        self.single_efetch[str, int, int].connect(self.label_2_update)
        self.single_efetch[str].connect(self.plainTextEdit_2_update)
        self.single_warning.connect(self.warning_message)
        self.pushButton_save.clicked.connect(self.pushButton_save_clicked)

    def get_fasta_by_entry(self, id_list):
        """
        通过biopython Entrez.efetch下载序列
        :param id_list: 登录号列表
        :return: none
        """
        Entrez.email = self.email  # 访问entry库必须的参数 邮箱信息 如果遇到什么问题，NCBI可以通过邮件联系到你。
        Entrez.tool = 'Biopython'  # 访问entry库使用的工具

        for db in ["nucleotide", "protein"]:    # 核算和蛋白数据库都进行搜索
            batch_size = int(self.batch_size)  # 步长 一次访问下载量
            for start in range(0, len(id_list), batch_size):
                end = min(len(id_list), start + batch_size)  # 一次下载的最后一个记录位置
                # print("Going to download %s %i to %i" % (db, start + 1, end))
                self.single_efetch[str, int, int].emit(db, int(start + 1), int(end))    # 下载进度 label更新
                try:
                    # 进行下载操作
                    hd_fetch = Entrez.efetch(db=db, rettype="fasta", retmode="text", id=','.join(id_list[start:end]))
                    records = hd_fetch.read()
                    self.single_efetch[str].emit(records)

                except Exception:
                    # 如果在数据库中没有搜索到会出现HTTP ERROR：bad request
                    continue
        self.single_efetch[str, int, int].emit('', 1, 1)    # 下载全部完成信号

    # 槽函数
    def pushButton_get_clicked(self):
        """
        get_fasta button 点击事件函数 下载左侧编辑框中的对应登录号的序列
        :return:
        """
        self.email = self.lineEdit.text()   # 读取邮箱
        if self.email == '':
            return QMessageBox.warning(self, 'warning', 'The email is not input!')
        self.batch_size = self.lineEdit_2.text()    # 步长
        if self.batch_size == '':
            return QMessageBox.warning(self, 'warning', 'The batch_size is not input!')
        self.plainTextEdit_2.clear()    # 每次get都清空显示结果的编辑框文本

        id_list = re.findall(Accession, self.plainTextEdit.toPlainText())   # 正则搜索登录号
        self.id_list = id_list
        # print(id_list)
        self.plainTextEdit.setPlainText('\n'.join(id_list))     # 显示提取出的登录号
        if len(id_list) == 0:
            return QMessageBox.warning(self, 'warning', 'The ID list is none or wrong format!')
        # 另起一线程进行下载 防止界面崩溃
        t = Thread(target=self.get_fasta_by_entry, args=(id_list,))
        t.setDaemon(True)
        t.start()

    def plainTextEdit_2_update(self, text: str):
        self.plainTextEdit_2.appendPlainText(text)

    def label_2_update(self, db: str, start: int, end: int):
        if db == '':
            self.label_2.setText("Download Complete!")
            all_download_seqs = self.plainTextEdit_2.toPlainText()
            not_found = [y for y in self.id_list if y not in all_download_seqs]  # 在list2列表中而不在list1列表中
            # print('noy_found:', not_found)  # 没有结果的登录号列表
        else:
            self.label_2.setText("Going to download %s %i to %i" % (db, start, end))

    def warning_message(self, warning: str):
        QMessageBox.warning(self, 'warning', warning)

    def pushButton_save_clicked(self):
        file_path, _ = QFileDialog.getSaveFileName(self, 'Save date', os.getcwd(), 'Sequences(*.fasta *.fa *.txt)')
        if len(file_path) != 0:
            with open(file_path, 'w') as f:
                f.write(self.plainTextEdit_2.toPlainText())


if __name__ == '__main__':
    QtCore.QCoreApplication.setAttribute(QtCore.Qt.AA_EnableHighDpiScaling)
    app = QApplication(sys.argv)
    MainWindow = NCBIEntry()
    MainWindow.show()
    sys.exit(app.exec_())
