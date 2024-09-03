# -*- coding : utf-8-*-

"""
__project_ = 'GenePredict'
__file_name__ = 'GenePredict'
__author__ = 'Dyh'
__time__ = '2021/11/10 19:09'
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
import subprocess
import sys
from io import StringIO

from PyQt5 import QtCore
from PyQt5.QtCore import pyqtSignal, pyqtSlot, QMetaObject, Qt, QThread
from PyQt5.QtGui import QIntValidator
from PyQt5.QtWidgets import *
import os
from uifiles import gene_predict
from threading import Thread
import requests
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from src import genscan_out, cnn_rec


def genscan_web(seq: str):
    genscan_result = re.compile(r'(>/tmp/[^<]*)', re.S)
    all_protein_cds = ''
    url = 'http://hollywood.mit.edu/cgi-bin/genscanw_py.cgi'
    header = {
        'User-Agent': 'Mozilla / 5.0(Windows NT 10.0; Win64;x64) AppleWebKit / 537.36(KHTML, likeGecko) Chrome / 84.0.4147.135 Safari / 537.36',
    }
    # POST请求构建DATA表单
    DATA = {
        '-o': 'Vertebrate',
        '-e': '1.00',
        '-p': 'Predicted CDS and peptides',

        '-s': seq
    }
    # # 字符串转二进制编码

    # DATA = urllib.parse.urlencode(DATA).encode(encoding='UTF8')
    # # 构建请求
    # request = urllib.request.Request(url=url, headers=header, data=DATA, method='POST')
    # # 防止网络错误
    try:
        response = requests.post(url, headers=header, timeout=100, data=DATA)
        html = response.text
    except Exception as e:
        print(e)
        return
    except requests.exceptions.Timeout as e:
        print(e)
        return
    except requests.exceptions.HTTPError as e:
        print(e)
        return

    results = re.findall(genscan_result, html)
    protein_cds = ''
    for result in results:
        protein_cds = protein_cds + result
        all_protein_cds = all_protein_cds + protein_cds + '\n' * 2

    return all_protein_cds


class GenePredict(QWidget, gene_predict.Ui_Form_GenePredict):
    single_cmd_out = pyqtSignal(str, str)  # 定义blast cmd完成信号
    single_blastdbcmd_out = pyqtSignal(list)
    single_blastdbcmd_out2 = pyqtSignal(str, str)
    single_genscan_out = pyqtSignal()
    single_blast_formatter_finish = pyqtSignal()
    single_cnn_show = pyqtSignal()

    def __init__(self, parent=None):
        super(GenePredict, self).__init__(parent)
        # 加载界面

        self.setupUi(self)

        # 初始化
        self.s_popen_pid = None
        self.db_path = os.path.abspath(r'./DB/blastdb')  # 本地blast数据库目录
        self.comboBox_db_files_addItems()  # 获取本地数据库 添加到组合框comboBox_db_file中
        self.blast_tool_bin_path = os.path.abspath(r'./tools/blast+/bin')
        self.tempfile_path = os.path.abspath(r'./tempfile')
        self.query_seq_path = os.path.abspath(r'./tempfile/query.fasta')
        self.out_file = os.path.abspath(r'./tempfile/out.txt')
        self.blastdbcmd_out_file = None
        self.genscan_prot_cds_list = []  # 用于保存genscan的结果 [[protein, cds],]
        self.blastdbcmd_out_seqs = []
        self.cnn_rec = cnn_rec.CNN_REC()
        self.genscan_out = genscan_out.GENSCAN_OUT()
        self.create_validator()

        # 信号 槽连接
        self.comboBox_app.currentIndexChanged.connect(self.comboBox_app_currentIndexChanged)
        self.pushButton_start.clicked.connect(self.pushButton_start_clicked)
        self.single_cmd_out.connect(self.table_view_get_single)
        self.pushButton_GENSCAN.clicked.connect(self.pushButton_GENSCAN_clicked)
        self.comboBox_db_file.currentIndexChanged.connect(self.combobox_db_currentIndexChanged)

        # self.comboBox_db_file.highlighted.connect(self.comboBox_db_files_addItems)
        self.pushButton_create_file.clicked.connect(self.btn_blastdbcmd_create_file_clicked)
        self.single_blastdbcmd_out2.connect(self.blastdbcmd_get_single)
        self.pushButton_GENSCAN_2.clicked.connect(self.pushButton_GENSCAN_2_clicked)
        self.single_genscan_out.connect(self.accept_single_genscan_out)
        self.pushButton_view.clicked.connect(self.btn_view_clicked)
        self.pushButton_cnn.clicked.connect(self.btn_cnn_clicked)
        self.single_cnn_show.connect(self.accept_single_cnn_show)

    def create_validator(self):
        validator = QIntValidator(self)
        self.lineEdit_begin.setValidator(validator)
        self.lineEdit_end.setValidator(validator)
        self.lineEdit_out_size.setValidator(validator)

    def comboBox_db_files_addItems(self):
        """
        获取数据库文件目录路径list
        :return: 所有目录下的全部数据库文件路径list
        """
        self.comboBox_db_file.clear()
        F = []  # 储存所有文件夹
        db_list = []
        for _, dirs, files in os.walk(self.db_path):
            for file in files:
                if self.comboBox_app.currentIndex() == 1 or self.comboBox_app.currentIndex() == 2:
                    if os.path.splitext(file)[1] == '.psq':
                        F.append(os.path.splitext(file)[0])
                    # 去除.nal中的DBLIST项
                    elif os.path.splitext(file)[1] == '.pal':
                        F.append(os.path.splitext(file)[0])
                        path = os.path.join(_, file)
                        with open(path, "r") as f:
                            i = 0
                            while 1:
                                line = f.readline()
                                if "DBLIST" in line or i > 10:
                                    break
                                i += 1
                        line = line.replace('DBLIST ', '').strip()
                        DBLIST = line.split(' ')
                        DBLIST = list(filter(None, DBLIST))
                        db_list = db_list + DBLIST
                else:
                    if os.path.splitext(file)[1] == '.nsq':
                        F.append(os.path.splitext(file)[0])
                    # 去除.nal中的DBLIST项
                    elif os.path.splitext(file)[1] == '.nal':
                        F.append(os.path.splitext(file)[0])
                        path = os.path.join(_, file)
                        with open(path, "r") as f:
                            i = 0
                            while 1:
                                line = f.readline()
                                if "DBLIST" in line or i > 10:
                                    break
                                i += 1
                        line = line.replace('DBLIST ', '').strip()
                        DBLIST = line.split(' ')
                        DBLIST = list(filter(None, DBLIST))
                        db_list = db_list + DBLIST
        # 去除.nal中的DBLIST项
        for db_name in db_list:
            if db_name in F:
                F.remove(db_name)
        # 添加db到组合框中
        self.comboBox_db_file.addItems(F)
        self.comboBox_cmd_db.addItems(F)

    def execute_command(self, tool_code, command):
        p = subprocess.Popen(command,
                             shell=True,
                             stdout=subprocess.PIPE,
                             stdin=subprocess.PIPE,
                             stderr=subprocess.STDOUT,
                             universal_newlines=True,
                             encoding='utf8',
                             )
        p.wait()  # 等待运行
        return_code = p.returncode  # 获取cmd执行结果代码
        # out = p.stdout.read()  # 返回cmd执行结果
        if return_code != 0:
            return QMessageBox.warning(self, 'warning', f'error code:{return_code}', QMessageBox.Yes, QMessageBox.Yes)
        if tool_code == 'blast_formatter':
            self.single_blast_formatter_finish.emit()

    def s_popen(self, tool_app, command_line):
        """
        独立线程执行命令行程序
        :param tool_app: 执行命令的工具 程序
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
        self.s_popen_pid = p.pid
        p.wait()  # 等待运行
        return_code = str(p.returncode)  # 获取cmd执行结果代码
        out = p.stdout.read()  # 返回cmd执行结果

        if tool_app == 'blast' and return_code == "0":
            # with open(self.out_file, 'r') as f:
            #     out = f.read()
            self.single_cmd_out.emit(return_code, out)  # 发出信号
            self.pushButton_start.setText('Start')
            self.pushButton_start.setEnabled(True)
        elif tool_app == 'blastdbcmd' and return_code == "0":
            with open(self.blastdbcmd_out_file, 'r') as f:
                out = f.read()
            self.single_blastdbcmd_out2.emit(return_code, out)
            self.pushButton_create_file.setText('Create File')
            self.pushButton_create_file.setEnabled(True)
        else:
            return QMessageBox.warning(self, 'warning', f'error code:{return_code}', QMessageBox.Yes, QMessageBox.Yes)

    def get_blast_command(self):
        """
        将blast cmd命令行组合起来
        :return: blast cmd命令
        """
        app_name = os.path.join(self.blast_tool_bin_path, self.comboBox_app.currentText())  # blast程序
        query = self.query_seq_path  # query序列文件目录
        db = os.path.join(self.db_path, self.comboBox_db_file.currentText())  # 数据库文件目录
        outfmt = '11'  # blast输出格式
        e_value = self.lineEdit_e_value.text()  # e值
        task = self.comboBox_task.currentText()  # blast工作模式
        more_option = self.lineEdit_more_option.text()  # 更多设置选项
        if self.comboBox_app.currentText() != 'tblastx':
            command = f'{app_name} -query {query} -db {db} -outfmt {outfmt} -evalue {e_value} -out {self.out_file} -task {task} {more_option}'
        else:
            command = f'{app_name} -query {query} -db {db} -outfmt {outfmt} -evalue {e_value} -out {self.out_file} {more_option}'
        return command

    def pushButton_start_clicked(self):
        """
        star 按钮点击事件
        :return:
        """
        sequence = self.plainTextEdit.toPlainText()
        if sequence == '':
            return QMessageBox.warning(self, 'warning', 'Please input a sequence.', QMessageBox.Yes, QMessageBox.Yes)
        else:
            with open(self.query_seq_path, 'w') as f:
                f.write(sequence)
        command = self.get_blast_command()
        # self.s_popen('blast', command)
        t = Thread(target=self.s_popen, args=('blast', command))
        t.setDaemon(True)
        t.start()
        self.pushButton_start.setText('please wait...')
        self.pushButton_start.setDisabled(True)

    def combobox_db_currentIndexChanged(self):
        Index = self.comboBox_db_file.currentIndex()
        self.comboBox_cmd_db.setCurrentIndex(Index)

    def comboBox_app_currentIndexChanged(self):
        """
        组合框conBox_app选中项改变事件
        将task加入到组合框comboBox_task中，并与conBox_app相对应
        """
        self.comboBox_db_files_addItems()
        self.comboBox_task.clear()
        # blasn
        if self.comboBox_app.currentIndex() == 0:
            self.comboBox_task.addItems(['megablast', 'dc-megablast', 'blastn', 'blastn-short', 'rmblastn'])
        # blastp
        elif self.comboBox_app.currentIndex() == 1:
            self.comboBox_task.addItems(['blastp', 'blastp-fast', 'blastp-short'])
        # blastx
        elif self.comboBox_app.currentIndex() == 2:
            self.comboBox_task.addItems(['blastx', 'blastx-fast'])
        # tblastn
        elif self.comboBox_app.currentIndex() == 3:
            self.comboBox_task.addItems(['tblastn', 'tblastn-fast'])
        # tblastx
        elif self.comboBox_app.currentIndex() == 4:
            self.comboBox_task.clear()

    def outfmt_blast_formatter(self, out_file):
        app_name = os.path.join(self.blast_tool_bin_path, 'blast_formatter')
        blast_out_0_file = os.path.abspath(r'./tempfile/blast_out_0.txt')
        # blast_out_6_file = os.path.abspath(r'./tempfile/blast_out_6.txt')
        command_0 = f'{app_name} -archive {out_file} -outfmt 0 -out {blast_out_0_file}'
        command_6 = f'{app_name} -archive {out_file} -outfmt 6 -out {self.out_file}'
        return command_0, command_6

    def table_view_get_single(self, return_code, out_str):
        """
        blast结果返回后，将结果加载进tableWidget中
        :param return_code: blast cmd返回代码
        :param out_str:blast 结果
        :return:
        """
        command0, command6 = self.outfmt_blast_formatter(self.out_file)
        os.system(command0)
        os.system(command6)
        with open(self.out_file, 'r') as f:
            out_str = f.read()
        if return_code != '0':
            return QMessageBox.warning(self, 'warning', f'error code{return_code}', QMessageBox.Yes, QMessageBox.Yes)
        elif return_code == '0' and out_str == '':
            return QMessageBox.information(self, 'No hit!', 'No matching sequence, try adjusting parameters.',
                                           QMessageBox.Yes, QMessageBox.Yes)

        out_lines = out_str.split('\n')  # blast fmt 6格式文本 按行分割
        row_num = len(out_lines) - 1  # 行数 最后有一空行 所以 -1
        self.tableWidget.setRowCount(row_num)  # tableWidget设置行数
        for i in range(row_num):
            out_line = out_lines[i]
            #  eg：NM_213948.1:25-525    NW_020185339.1	90.355	197	18	1	182	377	14243113	14243309	6.94e-69	267
            line_parts = out_line.split('\t')

            for j in range(len(line_parts)):
                if len(line_parts) == 12:
                    if str(line_parts[j]).isnumeric():
                        data2 = QTableWidgetItem()
                        data2.setData(Qt.EditRole, int(line_parts[j]))  # 将字符串转换成整数
                    else:
                        data2 = QTableWidgetItem(str(line_parts[j]))
                    # self.tableWidget.setColumnWidth(j, 60)  # 设置行宽
                    self.tableWidget.setItem(i, j, data2)  # 载入数据到tableWidget

    def blastdbcmd_get_single(self, return_code, out_str):
        if return_code != '0':
            return QMessageBox.warning(self, 'warning', f'error code{return_code}', QMessageBox.Yes, QMessageBox.Yes)
        if self.radioButton_reverse_complement.isChecked():
            fasta_sequence = SeqIO.read(StringIO(out_str), 'fasta')
            record = SeqRecord(
                fasta_sequence.seq.reverse_complement(),
                id=fasta_sequence.id,
                name=fasta_sequence.name,
                description=fasta_sequence.description)
            # record = fasta_sequence.reverse_complement
            out_str = record.format('fasta')
        with open(self.blastdbcmd_out_file, 'w') as f:
            f.write(out_str)
        os.system(f'start notepad {self.blastdbcmd_out_file}')

    def tableWidget_get_lines(self):
        """
        提取tableWidget选中的数据，如果未选，则返回全部内容
        :return: 返回tableWidget选中的数据
        """
        col_num = self.tableWidget.columnCount()  # 表头列数
        row_num = self.tableWidget.rowCount()  # 表行数
        selected_Indexes = self.tableWidget.selectedIndexes()  # 返回选中项
        if len(selected_Indexes) == 0:
            selected_row_num = [x for x in range(0, row_num)]
            # print(selected_row_num)
        else:
            selected_rows = []
            for index in selected_Indexes:
                selected_rows.append(index.row())  # 保存选中项row到列表中
            selected_row_num = list(set(selected_rows))  # 去重复项

        lines = []  # 用于保存所有行数据
        for row in selected_row_num:
            line = []  # 用于保存每行数据
            for col in range(col_num):
                line.append(self.tableWidget.item(row, col).text())
            lines.append(line)

        return lines

    def blastdbcmd_command(self, db, entry, out_file, seq_range):
        """
        从本地blastdb中提取相应序列的command命令
        :param out_file: 结果文件输出位置
        :param entry: 查找序列号
        :param seq_range: 提取范围
        :return: command命令
        """
        app_cmd = os.path.join(self.blast_tool_bin_path, 'blastdbcmd')  # blast程序
        # db = os.path.join(self.db_path, self.comboBox_db_file.currentText())  # 数据库文件目录
        command = f'{app_cmd} -entry {entry} -db "{db}" -out {out_file} -range {seq_range}'
        return command

    def get_blastdbcmd_commands(self, lines: list):
        out_size = int(self.lineEdit_out_size.text())
        db = os.path.join(self.db_path, self.comboBox_db_file.currentText())  # 数据库文件目录
        blastdbcmd_commands = []
        for line in lines:
            entry = line[1]
            sstart = int(line[8])
            send = int(line[9])
            if sstart > send:  # 反向序列
                start = (send - out_size) if (send - out_size) > 0 else 1
                end = sstart + out_size
                seq_range = f'{start}-{end}'
                strand = '-'
            else:  # 正向序列
                start = (sstart - out_size) if (sstart - out_size) > 0 else 1
                end = send + out_size
                seq_range = f'{start}-{end}'
                strand = '+'
            command = self.blastdbcmd_command(db, entry, self.out_file, seq_range)
            blastdbcmd_commands.append(command)
        return blastdbcmd_commands

    def seqio_genscan_result(self, genscan_result, seq_des: str):
        handle = SeqIO.parse(StringIO(genscan_result), 'fasta')
        fasta_sequences = []
        for fasta_sequence in handle:
            fasta_sequences.append(fasta_sequence)

        j = 1
        for i in range(0, len(fasta_sequences), 2):
            protein = SeqRecord(
                seq=fasta_sequences[i].seq,
                id=seq_des + str(j),
                name="protein",
                description=f'{seq_des + str(j)} protein|{fasta_sequences[i].description.split("|")[-1]}'
            )
            cds = SeqRecord(
                seq=fasta_sequences[i + 1].seq,
                id=seq_des + str(j),
                name="CDS",
                description=f'{seq_des + str(j)} CDS|{fasta_sequences[i + 1].description.split("|")[-1]}'
            )
            j = j + 1
            self.genscan_prot_cds_list.append([protein, cds])

        # print(self.genscan_prot_cds_list)

    def execute_blastdbcmd_commands(self, blastdbcmd_commands):

        i = 1
        for blastdbcmd_command in blastdbcmd_commands:
            p = subprocess.Popen(blastdbcmd_command,
                                 shell=True,
                                 stdout=subprocess.PIPE,
                                 stdin=subprocess.PIPE,
                                 stderr=subprocess.STDOUT,
                                 universal_newlines=True,
                                 encoding='utf8',
                                 )
            self.s_popen_pid = p.pid
            p.wait()  # 等待运行
            return_code = str(p.returncode)  # 获取cmd执行结果代码
            if return_code != '0':
                return QMessageBox.warning(self, 'Warning', f'blastdbcmd return code:{return_code}', QMessageBox.Yes,
                                           QMessageBox.Yes)
            with open(self.out_file, 'r') as f:
                out = f.read()
            genscan_result = genscan_web(out)
            if len(genscan_result) == 0:
                pass
            else:
                seq_des = f'hit{str(i)}-'
                self.seqio_genscan_result(genscan_result, seq_des)
            i = i + 1

        self.single_genscan_out.emit()

    def accept_single_genscan_out(self):
        self.genscan_out.genscan_prot_cds_list = self.genscan_prot_cds_list
        self.genscan_out.plainTextEdit_out.clear()
        self.genscan_out.show()
        self.genscan_out.plainTextEdit_load_seqs()
        self.pushButton_GENSCAN.setEnabled(True)
        self.pushButton_GENSCAN.setText('→GENSCAN Web Server')
        self.pushButton_GENSCAN_2.setEnabled(True)
        self.pushButton_GENSCAN_2.setText('→GENSCAN Web Server')

    def pushButton_GENSCAN_clicked(self):
        self.genscan_prot_cds_list.clear()
        selected_lines = self.tableWidget_get_lines()
        blastdbcmd_commands = self.get_blastdbcmd_commands(selected_lines)  # blastdbcmd 命令 list
        t = Thread(target=self.execute_blastdbcmd_commands, args=(blastdbcmd_commands,))
        t.setDaemon(True)
        t.start()
        self.pushButton_GENSCAN.setText('please wait...')
        self.pushButton_GENSCAN.setDisabled(True)
        # self.execute_blastdbcmd_commands(blastdbcmd_commands)

    def btn_blastdbcmd_create_file_clicked(self):
        db = os.path.join(self.db_path, self.comboBox_cmd_db.currentText())
        entry = self.lineEdit_entry.text()
        if entry == '':
            return QMessageBox.warning(self, 'warning', 'Please input a sequence.', QMessageBox.Yes, QMessageBox.Yes)
        self.blastdbcmd_out_file, _ = QFileDialog.getSaveFileName(self, 'save file', os.getcwd(),
                                                                  'Sequnces (*.fasta *.fa *.txt)')

        if self.radioButton_whole.isChecked():
            range_seq = f'1-'
        else:
            begin = self.lineEdit_begin.text()
            end = self.lineEdit_end.text()
            # print(end)
            if begin > end != '':
                range_seq = f'{end}-{begin}'
            else:
                range_seq = f'{begin}-{end}'
        command = self.blastdbcmd_command(db, entry, self.blastdbcmd_out_file, range_seq)
        t = Thread(target=self.s_popen, args=('blastdbcmd', command))
        t.setDaemon(True)
        t.start()
        self.pushButton_create_file.setText('please wait...')
        self.pushButton_create_file.setDisabled(True)
        # self.s_popen('blastdbcmd', command)

    def pushButton_GENSCAN_2_clicked(self):
        db = os.path.join(self.db_path, self.comboBox_cmd_db.currentText())
        entry = self.lineEdit_entry.text()
        command_list = []
        if entry == '':
            return QMessageBox.warning(self, 'warning', 'Please input a sequence.', QMessageBox.Yes, QMessageBox.Yes)
        if self.radioButton_whole.isChecked():
            range_seq = f'1-'
        else:
            begin = self.lineEdit_begin.text()
            end = self.lineEdit_end.text()
            if begin > end != '':
                range_seq = f'{end}-{begin}'
            else:
                range_seq = f'{begin}-{end}'
        command = self.blastdbcmd_command(db, entry, self.out_file, range_seq)
        command_list.append(command)
        # self.execute_blastdbcmd_commands(command_list)
        t = Thread(target=self.execute_blastdbcmd_commands, args=(command_list,))
        t.setDaemon(True)
        t.start()
        self.pushButton_GENSCAN_2.setText('please wait...')
        self.pushButton_GENSCAN_2.setDisabled(True)

    def btn_view_clicked(self):
        blast_out_0_file = os.path.abspath(r'./tempfile/blast_out_0.txt')
        if os.path.exists(blast_out_0_file):
            os.system(f'start notepad {blast_out_0_file}')

    def btn_cnn_clicked_fun(self):
        self.blastdbcmd_out_seqs.clear()
        lines = self.tableWidget_get_lines()
        out_size = int(self.lineEdit_out_size.text())
        db = os.path.join(self.db_path, self.comboBox_db_file.currentText())  # 数据库文件目录
        for line in lines:
            entry = line[1]
            sstart = int(line[8])
            send = int(line[9])
            if sstart > send:  # 反向序列
                start = (send - out_size) if (send - out_size) > 0 else 1
                end = sstart + out_size
                seq_range = f'{start}-{end}'
                strand = '-'
            else:  # 正向序列
                start = (sstart - out_size) if (sstart - out_size) > 0 else 1
                end = send + out_size
                seq_range = f'{start}-{end}'
                strand = '+'
            command = self.blastdbcmd_command(db, entry, self.out_file, seq_range)
            p = subprocess.Popen(command)
            p.wait()
            with open(self.out_file, 'r') as f:
                txt = f.read()
            fasta_seq = SeqIO.read(StringIO(txt), 'fasta')
            if strand == '-':
                fasta_seq = SeqRecord(
                    fasta_seq.seq.reverse_complement(),
                    id=fasta_seq.id,
                    name=fasta_seq.name,
                    description=fasta_seq.description)
            self.blastdbcmd_out_seqs.append(fasta_seq)
        # commands = self.get_blastdbcmd_commands(lines)
        # blastdbcmd_out_text = ''
        # for command in commands:
        #     # print(command)
        #     p = subprocess.Popen(command)
        #     p.wait()
        #     with open(self.out_file, 'r') as f:
        #         txt = f.read()
        #     blastdbcmd_out_text += f'{txt}\n'
        # handle = SeqIO.parse(StringIO(blastdbcmd_out_text), 'fasta')
        # for fasta_seq in handle:
        #     self.blastdbcmd_out_seqs.append(fasta_seq)

        # self.cnn_rec.seqs = self.blastdbcmd_out_seqs
        # self.cnn_rec.text_edit_load_seqs()
        self.single_cnn_show.emit()

    def accept_single_cnn_show(self):
        self.cnn_rec.show()
        self.cnn_rec.textEdit.clear()
        for blastdbcmd_out_seq in self.blastdbcmd_out_seqs:
            self.cnn_rec.textEdit.append(blastdbcmd_out_seq.format('fasta') + '\n')

    def btn_cnn_clicked(self):
        t = Thread(target=self.btn_cnn_clicked_fun)
        t.setDaemon(True)
        t.start()

    def __del__(self):
        """
        删除tempfile目录下的所有文件或文件夹
        """
        del_list = os.listdir('./tempfile')
        for f in del_list:
            file_path = os.path.join(r'./tempfile', f)
            if os.path.isfile(file_path):
                os.remove(file_path)


if __name__ == '__main__':
    QtCore.QCoreApplication.setAttribute(QtCore.Qt.AA_EnableHighDpiScaling)
    app = QApplication(sys.argv)
    MainWindow = GenePredict()
    MainWindow.show()
    sys.exit(app.exec_())
