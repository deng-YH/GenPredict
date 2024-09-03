# -*- coding: UTF-8 -*-
"""
__project_ = 'blast_X32_GUI'
__file_name__ = 'C_blastdbcmd'
__author__ = '鲨鱼辣椒'
__time__ = '2021/5/29 15:05'
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
import os
import subprocess
import sys
from threading import Thread

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from PyQt5.QtCore import QCoreApplication, Qt, pyqtSignal
from PyQt5.QtGui import QIntValidator
from PyQt5.QtWidgets import QWidget, QApplication, QFileDialog, QMessageBox
from uifiles import ui_blastdbcmd


class Blastdbcmd(QWidget, ui_blastdbcmd.Ui_Frame_Blastdbcmd):
    single_finish = pyqtSignal(str)  # 定义完成信号
    single_error = pyqtSignal(int, str)
    script_path = os.getcwd().replace('src', '')

    def __init__(self, parent=None):
        super(Blastdbcmd, self).__init__(parent)
        # 加载界面
        self.setupUi(self)
        # 初始化参数
        self.load_path = os.getcwd()

        # 限制输入
        self.create_validator()
        self.db_path = os.path.join(self.script_path, 'DB', 'blastdb')
        self.blastdbcmd = os.path.join(self.script_path, 'tools', 'blast+', 'bin', 'blastdbcmd')
        # 加载本地数据库列表
        self.comboBox_cmd_db_addItems()


        # 信号和槽
        self.single_finish.connect(self.accept_single_finish)
        self.single_error.connect(self.accept_single_error)
        self.pushButton_create_file.clicked.connect(self.btn_create_clicked)
        self.comboBox_cmd_db.highlighted.connect(self.comboBox_cmd_db_addItems)
    def create_validator(self):
        validator = QIntValidator(self)
        self.lineEdit_begin.setValidator(validator)
        self.lineEdit_end.setValidator(validator)

    def comboBox_cmd_db_addItems(self):
        """
        获取数据库文件目录路径list
        :return: 所有目录下的全部数据库文件路径list
        """
        self.comboBox_cmd_db.clear()
        F = []  # 储存所有文件夹
        db_list = []
        for _, dirs, files in os.walk(self.db_path):
            for file in files:
                if os.path.splitext(file)[1] == '.nsq':
                    F.append(os.path.splitext(file)[0])
                elif os.path.splitext(file)[1] == '.psq':
                    F.append(os.path.splitext(file)[0])
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
                    print(DBLIST)
                    DBLIST = list(filter(None, DBLIST))
                    db_list = db_list + DBLIST
        # 去除.nal中的DBLIST项
        for db_name in db_list:
            if db_name in F:
                F.remove(db_name)
        # 添加db到组合框中
        self.comboBox_cmd_db.addItems(F)

    def get_file_name(self):
        """
        QFileDialog.get???FileName()

        :return: open file name or save file name
        """
        file_path, file_type = QFileDialog.getSaveFileName(self, 'save file', self.load_path,
                                                           "Sequences (*.fasta *.txt);;All(*)")
        if file_path != '':
            self.load_path = os.path.split(file_path)[0]
            return file_path
        return None

    def _s_popen(self, command, out_file):
        """
        执行命令行命令
        :param command: cmd 命令
        :return: 执行结果
        """
        print(command)
        p = subprocess.Popen(command,
                             shell=True,
                             stdout=subprocess.PIPE,
                             stdin=subprocess.PIPE,
                             stderr=subprocess.STDOUT,
                             universal_newlines=True,
                             encoding='gbk',
                             )
        self.s_popen_pid = p.pid
        p.wait()  # 等待运行
        return_code = str(p.returncode)  # 获取cmd执行结果代码
        out = p.stdout.read()  # 返回cmd执行结果
        # print('out:', out)
        if 'error' not in out:  # 成功
            self.single_finish.emit(out_file)
        elif 'error' in out:  # 错误
            self.single_error.emit(return_code, out)

    def accept_single_finish(self, out_file):

        fasta_sequence = SeqIO.read(out_file, 'fasta')

        if self.radioButton_reverse_complement.isChecked():
            # 取反向互补序列
            record = SeqRecord(
                fasta_sequence.seq.reverse_complement(),
                id=fasta_sequence.id,
                name=fasta_sequence.name,
                description=fasta_sequence.description)
            # record = fasta_sequence.reverse_complement
            out_str = record.format('fasta')
        else:
            out_str = fasta_sequence.format('fasta')
        with open(out_file, 'w') as f:
            f.write(out_str)

        # 还原按钮状态
        self.pushButton_create_file.setText('Start')
        self.pushButton_create_file.setEnabled(True)

        # 打开
        os.system(f'start notepad "{out_file}"')

    def accept_single_error(self, return_code, out):
        QMessageBox.warning(self, 'error', f'error code:{return_code}\n{out}')
        self.pushButton_create_file.setText('Start')
        self.pushButton_create_file.setEnabled(True)

    def blastdbcmd_command(self, db, entry, out_file, seq_range):
        """
        从本地blastdb中提取相应序列的command命令
        :param db: 本地数据库
        :param out_file: 结果文件输出位置
        :param entry: 查找序列号
        :param seq_range: 提取范围
        :return: command命令
        """
        app_cmd = self.blastdbcmd
        # db = os.path.join(self.db_path, self.comboBox_db_file.currentText())  # 数据库文件目录
        command = f'{app_cmd} -entry {entry} -db "{db}" -out {out_file} -range {seq_range}'
        return command

    def blastdbcmd_start(self, command, blastdbcmd_out_file):
        t = Thread(target=self._s_popen, args=(command, blastdbcmd_out_file))
        t.setDaemon(True)
        t.start()

    def btn_create_clicked(self):
        entry = self.lineEdit_entry.text()
        if entry == '':
            return QMessageBox.information(self, 'information', 'Please input a sequence.', QMessageBox.Yes,
                                           QMessageBox.Yes)

        blastdbcmd_out_file = self.get_file_name()  # 输出文件目录

        if blastdbcmd_out_file is not None:
            db = os.path.join(self.db_path, self.comboBox_cmd_db.currentText())
            if self.radioButton_whole.isChecked():
                range_seq = f'1-'
            else:
                begin = int(self.lineEdit_begin.text())
                end = int(self.lineEdit_end.text())

                if begin > end != '':
                    range_seq = f'{end}-{begin}'
                else:
                    begin = 1 if (begin <= 0) else begin
                    range_seq = f'{begin}-{end}'
            command = self.blastdbcmd_command(db, entry, blastdbcmd_out_file, range_seq)
            self.pushButton_create_file.setText('please wait...')
            self.pushButton_create_file.setDisabled(True)
            self.blastdbcmd_start(command, blastdbcmd_out_file)
            pass


if __name__ == '__main__':
    QCoreApplication.setAttribute(Qt.AA_EnableHighDpiScaling)
    app = QApplication(sys.argv)
    MainWindow = Blastdbcmd()
    MainWindow.show()
    sys.exit(app.exec_())
