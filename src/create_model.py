# -*- coding: UTF-8 -*-
"""
__project_ = 'GenePredict'
__file_name__ = 'create_model'
__author__ = 'Dyh'
__time__ = '2022/2/15 9:45'
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
import ctypes
import inspect
import json
import subprocess
import sys

from PyQt5 import QtCore
from PyQt5.QtCore import pyqtSignal
from PyQt5.QtWidgets import QFileDialog, QWidget, QMessageBox, QApplication
import os
from uifiles import ui_CNN
from src.CNN import create_cnn_model, CNN_train
from threading import Thread


def _async_raise(tid, exctype):
    """raises the exception, performs cleanup if needed"""
    try:
        tid = ctypes.c_long(tid)
        if not inspect.isclass(exctype):
            exctype = type(exctype)
        res = ctypes.pythonapi.PyThreadState_SetAsyncExc(tid, ctypes.py_object(exctype))
        if res == 0:
            # pass
            raise ValueError("invalid thread id")
        elif res != 1:
            # """if it returns a number greater than one, you're in trouble,
            # and you should call it again with exc=NULL to revert the effect"""
            ctypes.pythonapi.PyThreadState_SetAsyncExc(tid, None)
            raise SystemError("PyThreadState_SetAsyncExc failed")
    except Exception as err:
        print(err)


def stop_thread(thread):
    """终止线程"""
    _async_raise(thread.ident, SystemExit)


class CreateModel(QWidget, ui_CNN.Ui_Form_CNN):
    single_out = pyqtSignal(str)
    single_get_data = pyqtSignal()
    script_path = os.path.split(os.path.realpath(__file__))[0]

    def __init__(self, parent=None):
        super(CreateModel, self).__init__(parent)
        # 加载界面
        self.setupUi(self)

        self.model_path = ''
        self.CNN_model = CNN_train.CNN_model()
        self.t = None

        self.CNN_model.single_finish.connect(self.train_finished)
        self.CNN_model.single_out.connect(self.update_out_textEdit)
        self.toolButton_genome_file.clicked.connect(lambda: self.get_file_name(0))
        self.toolButton_GFF_file.clicked.connect(lambda: self.get_file_name(1))
        self.pushButton_start.clicked.connect(self.start)
        self.single_out.connect(self.update_out_textEdit)
        self.single_get_data.connect(self.start_btn_clicked)

    def get_file_name(self, flag):
        """
        打开文件
        :return: open file name
        """
        file_path, file_type = QFileDialog.getOpenFileName(self, 'open file', os.getcwd(),
                                                           "Sequences (*.txt *.fna *.fa *.gff "
                                                           "*.fas *.gff3);;All(*)")
        if flag == 0:
            if file_path != '':
                if os.path.exists(file_path):
                    self.lineEdit_genome_file.setText(file_path)

        elif flag == 1:
            if file_path != '':
                if os.path.exists(file_path):
                    self.lineEdit_GFF_file.setText(file_path)

        return None

    def update_out_textEdit(self, info):
        self.textEdit.append(info)

    def make_model_dir(self):
        dir_name = self.lineEdit_model_name.text()
        model_path = self.script_path.replace(r'\src', r'\model')

        if not dir_name:
            QMessageBox.information(self, 'information', 'Please fill in the model name.')
            return False
        path = os.path.join(model_path, dir_name)
        print(path)
        if not os.path.exists(path):
            os.makedirs(path)
            print('The model directory was created successfully!')
        self.model_path = path
        return True

    def get_AG_GT_data(self):
        # 打开基因注释文件 逐行阅读
        self.single_out.emit('Extracting information from the GFF file......')
        gff_file = self.lineEdit_GFF_file.text()
        if os.path.exists(gff_file):
            all_mrna_information = create_cnn_model.get_all_mrna_information(gff_file)
            self.single_out.emit('After the extraction of information, %d mRNAs with exons and introns were found.\n'
                                 % len(all_mrna_information))
        else:
            QMessageBox.information(self, 'information', 'Please fill in the gff file path.')
            return False
        # 根据提取出的基因注释信息 在基因组数据中提取供体和受体序列 并处理为训练数据和测试数据 并保存
        self.single_out.emit('According to the extracted gene annotation information, the donor and recipient '
                             'sequences were extracted from the genome data.')
        train_num = int(self.lineEdit_train_num.text())
        genome_file = self.lineEdit_genome_file.text()
        if os.path.exists(genome_file):
            AG_train_data, AG_test_data, GT_train_data, GT_test_data = create_cnn_model.get_data_seqs(
                all_mrna_information,
                genome_file,
                train_num)
            self.script_path = self.script_path.replace(r'\src', '')
            print(os.path.join(self.script_path, r'data\train\AG.json'))
            with open(os.path.join(self.script_path, r'data\train\AG.json'), 'w') as f:
                json.dump(AG_train_data, f)
            with open(os.path.join(self.script_path, r'data\test\AG.json'), 'w') as f:
                json.dump(AG_test_data, f)
            with open(os.path.join(self.script_path, r'data\train\GT.json'), 'w') as f:
                json.dump(GT_train_data, f)
            with open(os.path.join(self.script_path, r'data\test\GT.json'), 'w') as f:
                json.dump(GT_test_data, f)
            self.single_out.emit('Training data and test data have been extracted!')
        else:
            QMessageBox.information(self, 'information', 'Please fill in the genome file path.')
            return False
        self.single_get_data.emit()
        return True

    def start(self):
        flag = self.pushButton_start.text()
        if flag == 'Start':
            if not self.make_model_dir():
                return False
            self.t = Thread(target=self.get_AG_GT_data)
            self.t.setDaemon(True)
            self.t.start()
            self.pushButton_start.setText('Stop')
        elif flag == 'Stop':
            self.stop_btn_clicked()

    def start_btn_clicked(self):
        self.CNN_model.model_path = self.model_path
        self.t = Thread(target=self.CNN_model.train2)
        self.t.setDaemon(True)
        self.t.start()

    def stop_btn_clicked(self):
        # _async_raise(self.t.ident, SystemExit)
        stop_thread(self.t)  # Signal termination
        print('sssss')
        self.pushButton_start.setText('Start')
        # self.pushButton_start.clicked.connect(self.start())
        # self.pushButton_start.clicked.disconnect(self.stop_btn_clicked)

    def train_finished(self):
        self.pushButton_start.setText('Start')
        # self.pushButton_start.clicked.connect(self.start())
        # self.pushButton_start.clicked.disconnect(self.stop_btn_clicked)


if __name__ == '__main__':
    QtCore.QCoreApplication.setAttribute(QtCore.Qt.AA_EnableHighDpiScaling)
    app = QApplication(sys.argv)
    MainWindow = CreateModel()
    MainWindow.show()
    sys.exit(app.exec_())
