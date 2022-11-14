# -*- coding: UTF-8 -*-
"""
__project_ = 'GenePredict'
__file_name__ = 'cnn_rec'
__author__ = 'Dyh'
__time__ = '2021/11/13 10:40'
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
# import subprocess
import sys
from Bio import SeqIO
from PyQt5 import QtCore
from PyQt5.QtCore import pyqtSignal
from PyQt5.QtWidgets import QWidget, QApplication
from uifiles import cnn_rec
from threading import Thread
from src import cnn, translate
from io import StringIO


class CNN_REC(QWidget, cnn_rec.Ui_Form_CNN_REC):
    script_path = os.path.split(os.path.realpath(__file__))[0]
    single_calculated_result = pyqtSignal(str)

    def __init__(self, parent=None):
        super(CNN_REC, self).__init__(parent)
        # 加载界面
        self.setupUi(self)

        self.model_path = os.path.abspath(os.path.join(self.script_path, '../model'))
        self.seqs = []
        self.CNN_MODEL_AG = cnn.CNN_model()
        # self.CNN_MODEL_GT = cnn.CNN_model()
        self.combobox_model_load_model_path()
        self.translate = translate.TRANSLATE()

        self.textEdit.selectionChanged.connect(self.textedit_selection_changed)
        self.pushButton_calculate.clicked.connect(self.text_edit_load_seqs)
        self.single_calculated_result.connect(self.accept_single_calculated_result)

    def text_edit_load_seqs(self):
        self.seqs.clear()
        textEdit_txt = self.textEdit.toPlainText()
        # print(textEdit_txt)
        handle = SeqIO.parse(StringIO(textEdit_txt), 'fasta')
        for fasta_seq in handle:
            self.seqs.append(fasta_seq)
        self.textEdit.setText('')

        all_seq_str = ''
        for fasta_seq in self.seqs:
            AG_model_file = os.path.join(self.comboBox_model.currentText(), 'AG.pth')
            GT_model_file = os.path.join(self.comboBox_model.currentText(), 'GT.pth')

            AG_sites, GT_sites = self.CNN_MODEL_AG.test1(str(fasta_seq.seq), AG_model_file, GT_model_file)

            flag_sites_num = {}

            seq_str = f'>{fasta_seq.description}<br/>{fasta_seq.seq}<br/><br/>'
            for AG_site in AG_sites:
                num = seq_str.find(AG_site)

                flag_sites_num[num] = 'AG'

            for GT_site in GT_sites:
                num = seq_str.find(GT_site)

                flag_sites_num[num] = 'GT'

            html_seq = ''
            num1 = 0
            num2 = 0
            i = 0
            len_site = len(AG_sites[0])
            for num in sorted(flag_sites_num):
                i = i + 1

                num2 = num
                flag = flag_sites_num[num]
                if num1 == 0:
                    html_seq = html_seq + seq_str[:num2 + len_site // 2 - 1]
                    if flag == 'AG':
                        html_seq = html_seq + f'<font style="font-weight:bold;background-color:red">{seq_str[num2 + len_site // 2 - 1:num2 + len_site // 2 + 1]}</font>'

                    else:

                        html_seq = html_seq + f'<font style="font-weight:bold;background-color:yellow">{seq_str[num2 + len_site // 2 - 1:num2 + len_site // 2 + 1]}</font>'

                else:
                    html_seq = html_seq + seq_str[num1 + len_site // 2 + 1:num2 + len_site // 2 - 1]
                    if flag == 'AG':
                        html_seq = html_seq + f'<font style="font-weight:bold;background-color:red">{seq_str[num2 + len_site // 2 - 1:num2 + len_site // 2 + 1]}</font>'

                    else:
                        html_seq = html_seq + f'<font style="font-weight:bold;background-color:yellow">{seq_str[num2 + len_site // 2 - 1:num2 + len_site // 2 + 1]}</font>'

                num1 = num
                if i == len(flag_sites_num):
                    html_seq = html_seq + seq_str[num2 + len_site // 2 + 1:]
            all_seq_str = all_seq_str + html_seq + '<br/><br/>'

        self.single_calculated_result.emit(all_seq_str)

    def combobox_model_load_model_path(self):
        for root, dirs, files in os.walk(self.model_path):
            if root != self.model_path:
                self.comboBox_model.addItem(root)

    def textedit_selection_changed(self):
        select_txt = self.textEdit.createMimeDataFromSelection().text()
        num = len(select_txt)
        self.label_selection_changed.setText(f'Number of selected characters：{num}')
        Codon_table = self.comboBox_genetic_codes.currentText()
        translate_result = self.translate.translate(Codon_table, select_txt)
        if translate_result is None:
            return
        # 将结果加入到编辑框中
        self.textEdit_tanslate.clear()
        for i in range(len(translate_result[0])):
            self.textEdit_tanslate.append(translate_result[1].format(translate_result[0][i]))
            self.textEdit_tanslate.append(translate_result[2][i])

    def accept_single_calculated_result(self, seq_str):
        self.textEdit.setHtml(seq_str)

    def btn_calculate_clicked(self):
        t = Thread(target=self.text_edit_load_seqs)
        t.setDaemon(True)
        t.start()


if __name__ == '__main__':
    QtCore.QCoreApplication.setAttribute(QtCore.Qt.AA_EnableHighDpiScaling)
    app = QApplication(sys.argv)
    MainWindow = CNN_REC()
    MainWindow.show()
    sys.exit(app.exec_())
