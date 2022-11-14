# -*- coding: UTF-8 -*-
"""
__project_ = 'GenePredict'
__file_name__ = 'translate'
__author__ = 'Dyh'
__time__ = '2021/11/13 19:19'
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
import re
import sys

from Bio import SeqIO
from io import StringIO

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from PyQt5.QtCore import QCoreApplication, Qt, pyqtSignal
from PyQt5.QtWidgets import QApplication, QMessageBox, QWidget
from uifiles import ui_translate

genetic_codes = {
    'Standard': 1,
    'Vertebrate Mitochondrial': 2,
    'Yeast Mitochondrial': 3,
    'Mold, Protozoan, and Coelenterate Mitochondrial Code and the Mycoplasma/Spiroplasma': 4,
    'Invertebrate Mitochondrial': 5,
    'Ciliate, Dasycladacean and Hexamita Nuclear': 6,
    'Echinoderm and Flatworm Mitochondrial': 9,
    'Euplotid Nuclear': 10,
    'Bacterial, Archaeal and Plant Plastid': 11,
    'Alternative Yeast Nuclear': 12,
    'Ascidian Mitochondrial': 13,
    'Alternative Flatworm Mitochondrial': 14,
    'Chlorophycean Mitochondrial': 16,
    'Trematode Mitochondrial': 21,
    'Scenedesmus obliquus Mitochondrial': 22,
    'Thraustochytrium Mitochondrial': 23,
    'Rhabdopleuridae Mitochondrial': 24,
    'Candidate Division SR1 and Gracilibacteria': 25,
    'Pachysolen tannophilus Nuclear': 26,
    'Karyorelict Nuclear': 27,
    'Condylostoma Nuclear': 28,
    'Mesodinium Nuclear': 29,
    'Peritrich Nuclear': 30,
    'Blastocrithidia Nuclear': 31,
    'Cephalodiscidae Mitochondrial UAA-Tyr': 32
}


class TRANSLATE(object):

    def __init__(self):
        # 初始化参数
        self.title = 'Translate'
        self.seq = None  # 翻译序列
        self.table = None  # 密码子表
        self.protein_list = None  # 蛋白序列列表

    def translate_seq(self):
        """
        读取序列，正向3种组合 反向3种组合 分别进行翻译
        :return:蛋白序列列表
        """
        seq = self.seq

        def three_comb_trans(seqs):
            """
            所有三联组合进行翻译
            """
            n = len(seqs) % 3
            three_comb = seqs[:len(seqs) - n]
            # print(str(three_comb.translate(stop_symbol='*', table=self.table)))
            return str(three_comb.translate(stop_symbol='*', table=self.table))

        seq_list = ['', '', '', '', '', '']
        self.protein_list = []
        for i in range(3):
            seq_list[i] = (seq[i:])
            seq_list[i + 3] = (seq.reverse_complement()[i:])
        for seq in seq_list:
            # print(seq)
            self.protein_list.append(three_comb_trans(seq))

    def prot_to_xml(self):
        """
        转换为html格式
        :return:
        """
        protein_xml = []
        # 给起始M和ORF加标签
        M_xml = '<a style="color:red;font-weight:bold">{0}</a>'
        orf_xml = '<a style="background-color:#F9B1B6">{0}</a>'
        for protein_seq in self.protein_list:
            orfs = re.findall(r'(M[A-Z]*\*?)', protein_seq)
            for orf in orfs:
                M_orf = str(orf).replace('M', M_xml.format('M'))
                protein_seq = protein_seq.replace(orf, orf_xml.format(M_orf))
            protein_xml.append(protein_seq)
        # print(protein_xml)
        return protein_xml

    # 槽函数
    def translate(self, Codon_table, txt):
        """
        读取序列并进行正反向三三组合  翻译 然后添加标签修改字体 ORF为红色背景 起始M加粗标红
        :return:
        """
        # 单序列翻译
        txt = str(txt).strip().replace(' ', '')
        nucl_code = [65, 71, 84, 67, 97, 103, 116, 99, 88, 120]
        nucl = ''
        for i in txt:
            if ord(i) in nucl_code:
                nucl += i
        if nucl == '':
            return None
        # self.seq = SeqIO.read(StringIO(txt), 'fasta').seq
        fasta_seq = SeqRecord(
            seq=Seq(nucl),
            description='seq_selected'
        )
        self.seq = fasta_seq.seq

        self.table = genetic_codes[Codon_table]  # 读取编码方式

        self.translate_seq()

        head = ["5'3' Frame 1", "5'3' Frame 2", "5'3' Frame 3", "3'5' Frame 1", "3'5' Frame 2", "3'5' Frame 3"]
        title = '<a style="background-color: #FFE6FF">>{0}</a>'

        trans_xml_list = self.prot_to_xml()
        return head, title, trans_xml_list
        # # 将结果加入到编辑框中
        # for i in range(len(head)):
        #     self.textEdit.append(title.format(head[i]))
        #     self.textEdit.append(trans_xml[i])


class C_translate(QWidget, ui_translate.Ui_Translate):
    single_finish = pyqtSignal(str)  # 定义完成信号

    def __init__(self, parent=None):
        super(C_translate, self).__init__(parent)
        # 加载界面
        self.setupUi(self)
        # 初始化参数
        self.title = 'Translate'
        self.seq = None  # 翻译序列
        self.table = None  # 密码子表
        self.protein_list = None  # 蛋白序列列表
        # 信号和槽
        self.pushButton_translate.clicked.connect(self.translate)

    def translate_seq(self):
        """
        读取序列，正向3种组合 反向3种组合 分别进行翻译
        :return:蛋白序列列表
        """
        seq = self.seq

        def three_comb_trans(seqs):
            """
            所有三联组合进行翻译
            """
            n = len(seqs) % 3
            three_comb = seqs[:len(seqs) - n]
            # print(str(three_comb.translate(stop_symbol='*', table=self.table)))
            return str(three_comb.translate(stop_symbol='*', table=self.table))

        seq_list = ['', '', '', '', '', '']
        self.protein_list = []
        for i in range(3):
            seq_list[i] = (seq[i:])
            seq_list[i + 3] = (seq.reverse_complement()[i:])
        for seq in seq_list:
            # print(seq)
            self.protein_list.append(three_comb_trans(seq))
        # print(self.protein_list)

    def prot_to_xml(self):
        """
        转换为html格式
        :return:
        """
        protein_xml = []
        # 给起始M和ORF加标签
        M_xml = '<a style="color:red;font-weight:bold">{0}</a>'
        orf_xml = '<a style="background-color:#F9B1B6">{0}</a>'
        for protein_seq in self.protein_list:
            orfs = re.findall(r'(M[A-Z]*\*?)', protein_seq)
            for orf in orfs:
                M_orf = str(orf).replace('M', M_xml.format('M'))
                protein_seq = protein_seq.replace(orf, orf_xml.format(M_orf))
            protein_xml.append(protein_seq)
        # print(protein_xml)
        return protein_xml

    # 槽函数
    def translate(self):
        """
        读取序列并进行正反向三三组合  翻译 然后添加标签修改字体 ORF为红色背景 起始M加粗标红
        :return:
        """
        text = self.plainTextEdit.toPlainText().strip().upper()  # 读取文本框中的序列文本
        if text == '':
            return QMessageBox.warning(self, 'warning', 'The sequence is not input!')
        # 单序列翻译
        if '>' in text:  # 如果是标准fasta序列则用SeqIO进行读取
            try:
                self.seq = SeqIO.read(StringIO(text), 'fasta').seq
            except Exception:
                QMessageBox.warning(self, 'warning', 'Please enter one sequence!')
        elif text != '':
            text = f'>seq\n{text}'
            self.seq = SeqIO.read(StringIO(text), 'fasta').seq

        self.table = genetic_codes[self.comboBox_genetic_codes.currentText()]  # 读取编码方式

        self.translate_seq()
        self.textEdit.setText('')  # 清空之前的数据
        head = ["5'3' Frame 1", "5'3' Frame 2", "5'3' Frame 3", "3'5' Frame 1", "3'5' Frame 2", "3'5' Frame 3"]
        title = '<a style="background-color: #FFE6FF">>{0}</a>'

        trans_xml = self.prot_to_xml()
        # 将结果加入到编辑框中
        for i in range(len(head)):
            self.textEdit.append(title.format(head[i]))
            self.textEdit.append(trans_xml[i])


if __name__ == '__main__':
    QCoreApplication.setAttribute(Qt.AA_EnableHighDpiScaling)
    app = QApplication(sys.argv)
    MainWindow = C_translate()
    MainWindow.show()
    sys.exit(app.exec_())
