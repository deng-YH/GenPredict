# -*- coding: UTF-8 -*-
"""
__project_ = 'GenePredict'
__file_name__ = 'genscan_out'
__author__ = 'Dyh'
__time__ = '2021/11/11 21:08'
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
import subprocess
import sys
# from Bio import SeqIO, AlignIO
from PyQt5 import QtCore
from PyQt5.QtCore import pyqtSignal
from PyQt5.QtWidgets import QWidget, QApplication
from uifiles import genscan_out
from src import DiamondMakeDB
from threading import Thread


class GENSCAN_OUT(QWidget, genscan_out.Ui_Frame_GENSCAN_OUT):
    single_diamond = pyqtSignal()  # 定义完成信号
    single_finish = pyqtSignal()

    def __init__(self, genscan_prot_cds_list=None, parent=None):
        super(GENSCAN_OUT, self).__init__(parent)
        # 加载界面
        self.setupUi(self)

        self.genscan_prot_cds_list = genscan_prot_cds_list
        self.db_path = os.path.abspath(r'./DB/diamonddb')
        self.app_path = os.path.abspath(r'./tools/diamond')
        self.get_all_db_file()  # 初始化DB列表
        self.protein_file = os.path.abspath(r'./tempfile/protein.fasta')
        self.cds_file = os.path.abspath(r'./tempfile/cds.fasta')
        self.out_file = os.path.abspath(r'./tempfile/out.txt')
        self.diamond_output_information_file = os.path.abspath(r'./doc/diamond_output_information.txt')
        self.diamond_makedb = DiamondMakeDB.MakeDB()
        self.out_flag = 0

        self.pushButton_makedb.clicked.connect(self.make_db_clicked)
        self.single_diamond.connect(self.plainTextEdit_out_update)
        self.diamond_makedb.single_finish.connect(self.get_all_db_file)
        self.pushButton_start.clicked.connect(self.btn_start_clicked)
        self.pushButton_help.clicked.connect(self.btn_help_clicked)
        self.single_finish.connect(self.accept_single_finished)

    def get_all_db_file(self):
        """
        获取数据库文件目录路径list
        :return: 所有目录下的全部数据库文件路径list
        """
        F = []  # 储存所有文件夹
        for _, dirs, files in os.walk(self.db_path):
            for file in files:
                if os.path.splitext(file)[1] == '.dmnd':
                    F.append(os.path.splitext(file)[0])
        self.comboBox_db_files.addItems(F)

    def plainTextEdit_load_seqs(self):
        for protein_cds in self.genscan_prot_cds_list:
            text = f'{protein_cds[0].format("fasta")}\n' \
                   f'{protein_cds[1].format("fasta")}\n' \
                   f'{"-" * 50}\n'
            self.plainTextEdit_out.appendPlainText(text)

    def s_popen(self, command_line):
        """
        独立线程执行命令行程序
        :param command_line: cmd 命令行
        :return: 执行结果 成功则返回打印信息 失败则弹出信息框
        """
        p = subprocess.Popen(command_line,
                             shell=True,
                             stdout=subprocess.PIPE,
                             stdin=subprocess.PIPE,
                             stderr=subprocess.STDOUT,
                             universal_newlines=True,
                             encoding='utf8',
                             )
        p.wait()  # 等待运行
        # return_code = str(p.returncode)  # 获取cmd执行结果代码
        # out = p.stdout.read()  # 返回cmd执行结果
        # print(return_code)

        self.single_diamond.emit()

    def make_db_clicked(self):
        self.diamond_makedb.show()

    def get_diamond_command(self, query_file):
        diamond_app = f'"{self.app_path}" {self.comboBox_app.currentText()}'
        db = os.path.join(self.db_path, self.comboBox_db_files.currentText())
        k = self.lineEdit_k.text()
        e_value = self.lineEdit_e.text()
        output = self.lineEdit_output.text()
        command = f'{diamond_app} --query "{query_file}" --db "{db}" -k {k} -e {e_value} --out "{self.out_file}" --outfmt 6 ' \
                  f'{output}'
        return command

    def show_result(self):
        self.plainTextEdit_out.clear()
        for protein_cds in self.genscan_prot_cds_list:
            if self.comboBox_app.currentIndex() == 0:
                # print(self.protein_file)
                with open(self.protein_file, 'w') as f:
                    f.write(str(protein_cds[0].format('fasta')))
                command = self.get_diamond_command(self.protein_file)
            else:
                # print(self.cds_file)
                with open(self.cds_file, 'w') as f:
                    f.write(str(protein_cds[1].format('fasta')))
                command = self.get_diamond_command(self.cds_file)
            self.s_popen(command)
        self.single_finish.emit()

    def plainTextEdit_out_update(self):
        with open(self.out_file, 'r') as f:
            out = f.read()
        text = f'{self.genscan_prot_cds_list[self.out_flag][0].format("fasta")}\n' \
               f'{self.genscan_prot_cds_list[self.out_flag][1].format("fasta")}\n' \
               f'diamond by {self.comboBox_app.currentText()}:\n' \
               f'{out}\n{"-" * 50}\n'
        self.plainTextEdit_out.appendPlainText(text)
        self.out_flag = self.out_flag + 1
        if self.out_flag == len(self.genscan_prot_cds_list):
            self.out_flag = 0

    def accept_single_finished(self):
        self.pushButton_start.setText('start')
        self.pushButton_start.setEnabled(True)

    def btn_help_clicked(self):
        os.system(f'start notepad {self.diamond_output_information_file}')

    def btn_start_clicked(self):
        self.pushButton_start.setText('please wait...')
        self.pushButton_start.setDisabled(True)
        t = Thread(target=self.show_result)
        t.setDaemon(True)
        t.start()


if __name__ == '__main__':
    QtCore.QCoreApplication.setAttribute(QtCore.Qt.AA_EnableHighDpiScaling)
    app = QApplication(sys.argv)
    MainWindow = GENSCAN_OUT()
    MainWindow.show()
    sys.exit(app.exec_())
