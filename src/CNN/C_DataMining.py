# -*- coding: UTF-8 -*-
"""
__project_ = 'blast_X32_GUI'
__file_name__ = 'C_clastalW'
__author__ = '鲨鱼辣椒'
__time__ = '2021/6/9 14:51'
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
import sys
import DataProcessing
import CNN_train
from PyQt5.QtCore import pyqtSignal
from threading import Thread
from PyQt5 import QtCore
from uifiles import DataMining
from PyQt5.QtWidgets import QApplication, QFileDialog, QMessageBox, QWidget

CNN_file = os.getcwd()


class C_DataMining(QWidget, DataMining.Ui_DM):
    single_out = pyqtSignal(str)  # 输出信号

    def __init__(self, parent=None):
        super(C_DataMining, self).__init__(parent)
        # 加载界面
        self.setupUi(self)
        self.sub_data_dic = ''  # 子数据集文件夹
        self.input_directory = os.getcwd()  # 起始文件目录
        self.genome_file = ''  # 基因组文件路径
        self.gff_file = ''  # 基因注释文件路径
        self.donor_num = 0  # 供体剪切位点数据量
        self.acceptor_num = 0  # 受体剪切位点数据量
        self.non_donor_num = 0  # 非供体剪切位点数据量
        self.non_acceptor_num = 0  # 非受体剪切位点数据量
        self.k_mer_size = 0  # K-MER 编码长度
        self.CNN = CNN_train.CNN_model()
        # 信号和槽
        self.pushButton_get_data.clicked.connect(self.pushButton_get_data_clicked)
        self.toolButton.clicked.connect(lambda: self.open_file_clicked(1))  # 读取genome文件路径
        self.toolButton_2.clicked.connect(lambda: self.open_file_clicked(0))  # 读取gff文件路径
        self.single_out.connect(self.plainTextEdit_out_update)  # 更新文本框
        self.CNN.single_out.connect(self.plainTextEdit_out_update)

    def get_all_mrna_information(self):
        # 打开基因注释文件 逐行阅读
        self.single_out.emit('Extracting information from the GFF file......')
        all_mrna_information = DataProcessing.get_all_mrna_information(self.gff_file)
        self.single_out.emit('After the extraction of information, %d mRNAs with exons and introns were found.\n'
                             % len(all_mrna_information))
        # 根据提取出的基因注释信息 在基因组数据中提取供体和受体序列 并处理为训练数据和测试数据 并保存
        self.single_out.emit('According to the extracted gene annotation information, the donor and recipient '
                             'sequences were extracted from the genome data.')
        AG_train_data, AG_test_data, GT_train_data, GT_test_data = DataProcessing.get_data_seqs(all_mrna_information,
                                                                                                self.genome_file,
                                                                                                self.donor_num,
                                                                                                self.acceptor_num,
                                                                                                self.non_donor_num,
                                                                                                self.non_acceptor_num)
        with open(CNN_file + r'\data\train\AG.json', 'w') as f:
            json.dump(AG_train_data, f)
        with open(CNN_file + r'\data\test\AG.json', 'w') as f:
            json.dump(AG_test_data, f)
        with open(CNN_file + r'\data\train\GT.json', 'w') as f:
            json.dump(GT_train_data, f)
        with open(CNN_file + r'\data\test\GT.json', 'w') as f:
            json.dump(GT_test_data, f)
        self.single_out.emit('Training data and test data have been extracted!')
        print('建立模型')
        self.CNN.dataset = 'GT.json'
        self.CNN.train()
        print('结束')

    def pushButton_get_data_clicked(self):
        self.plainTextEdit_out.setPlainText('')
        self.genome_file = self.lineEdit_genome_file.text()
        if not os.path.exists(self.genome_file):
            return QMessageBox.warning(self, 'warning', 'The genome file is not input!')
        self.sub_data_dic = os.path.split(self.genome_file)[1].split('.')[0]
        self.gff_file = self.lineEdit_gff_file.text()
        if not os.path.exists(self.gff_file):
            return QMessageBox.warning(self, 'warning', 'The gff file is not input!')
        if self.checkBox_donor.isChecked():
            if self.lineEdit_donor.text() != '':
                self.donor_num = int(self.lineEdit_donor.text())
        if self.checkBox_acceptor.isChecked():
            if self.lineEdit_acceptor.text() != '':
                self.acceptor_num = int(self.lineEdit_acceptor.text())
        if self.checkBox_non_donor.isChecked():
            if self.lineEdit_non_donor.text() != '':
                self.non_donor_num = int(self.lineEdit_non_donor.text()) * 5
        if self.checkBox_non_acceptor.isChecked():
            if self.lineEdit_non_acceptor.text() != '':
                self.non_acceptor_num = int(self.lineEdit_non_acceptor.text()) * 5
        print(self.genome_file)
        t = Thread(target=self.get_all_mrna_information)
        t.setDaemon(True)
        t.start()

    def plainTextEdit_out_update(self, out: str):
        self.plainTextEdit_out.appendPlainText(out)

    def open_file_clicked(self, q_s):
        """
        读取genome_file or gff_file文件路径
        :param q_s: 1：genome_file or 0：gff_file
        :return:None
        """
        file_path, file_type = QFileDialog.getOpenFileName(self, 'open file', self.input_directory,
                                                           "All(*)")
        if file_path != '':
            if q_s == 1:
                self.genome_file = file_path
                # path = os.path.split(self.query_file)[1]
                self.input_directory = os.path.split(self.genome_file)[0]
                self.lineEdit_genome_file.setText(self.genome_file)
            elif q_s == 0:
                self.gff_file = file_path
                # path = os.path.split(self.gff_file)[1]
                self.input_directory = os.path.split(self.gff_file)[0]
                self.lineEdit_gff_file.setText(self.gff_file)


if __name__ == '__main__':
    QtCore.QCoreApplication.setAttribute(QtCore.Qt.AA_EnableHighDpiScaling)
    app = QApplication(sys.argv)
    MainWindow = C_DataMining()
    MainWindow.show()
    sys.exit(app.exec_())
