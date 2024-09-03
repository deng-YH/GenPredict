# -*- coding: UTF-8 -*-
"""
__project_ = 'GenePredict'
__file_name__ = 'local_blast'
__author__ = 'Dyh'
__time__ = '2021/11/14 15:29'
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
from threading import Thread

from PyQt5 import QtCore
from PyQt5.QtCore import pyqtSignal
from PyQt5.QtGui import QIntValidator
from PyQt5.QtWidgets import QApplication, QFileDialog, QMessageBox, QWidget

from uifiles import ui_local_blast

from src import blast_make_db


class LocalBlast(QWidget, ui_local_blast.Ui_Form_BLAST):
    script_path = os.getcwd().replace('src', '')
    single_finish = pyqtSignal(str)  # 定义完成信号
    single_error = pyqtSignal(int, str)

    def __init__(self, parent=None):
        super(LocalBlast, self).__init__(parent)
        # 加载界面
        self.setupUi(self)
        self.groupBox_subject.hide()
        # self.comboBox_db_files.setEditable(True)  # 下拉框是否允许编辑

        # 初始化参数
        self.load_path = os.getcwd()
        self.db_path = os.path.join(self.script_path, 'DB', 'blastdb')
        self.blast_path = os.path.join(self.script_path, 'tools', 'blast+', 'bin')
        self.query_file = os.path.join(self.script_path, 'tempfile', 'blast_query.fasta')
        self.subject_file = os.path.join(self.script_path, 'tempfile', 'blast_subject.fasta')
        self.out_file = os.path.join(self.script_path, 'tempfile', 'blast_out.fasta')
        print(self.blast_path)
        # 加载diamond_db_files
        self.comboBox_db_files_addItems()

        # 实例一个makeblastdb对象
        self.MakeBlastDB = blast_make_db.MakeBlastDB()
        # 限制输入
        self.create_validator()

        # 信号和槽
        self.open_file_query.clicked.connect(self.btn_open_file_query_clicked)
        self.open_file_subject.clicked.connect(self.btn_open_file_subject_clicked)
        self.pushButton_clear.clicked.connect(lambda: self.btn_clear_clicked(self.label_choose_file_q))
        self.pushButton_clear2.clicked.connect(lambda: self.btn_clear_clicked(self.label_choose_file_s))
        self.CheckBox_two_seq.stateChanged.connect(self.check_box_is_checked)
        self.comboBox_app.currentIndexChanged.connect(self.comboBox_app_currentIndexChanged)
        self.PushButton_make_db.clicked.connect(self.btn_make_blast_db_clicked)
        self.pushButton_BLAST.clicked.connect(self.btn_blast_clicked)
        # self.comboBox_db_files.activated.connect(self.comboBox_db_files_addItems)
        self.single_finish.connect(self.accept_single_finish)
        self.single_error.connect(self.accept_single_error)
        self.MakeBlastDB.single_finish.connect(self.comboBox_db_files_addItems)

    def create_validator(self):
        validator = QIntValidator(0, 18, self)
        self.lineEdit_outfmt.setValidator(validator)

    def comboBox_db_files_addItems(self):
        """
        获取数据库文件目录路径list
        :return: 所有目录下的全部数据库文件路径list
        """
        self.comboBox_db_files.clear()
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
        self.comboBox_db_files.addItems(F)

    def get_file_name(self, option):
        """
        QFileDialog.get???FileName()
        :param option: query or subject
        :return: query file name or subject file name
        """
        if option == 'query':
            file_path, file_type = QFileDialog.getOpenFileName(self, 'open file', self.load_path,
                                                               "Sequences (*.txt *.fna *.fa *.fasta *.nex "
                                                               "*.fas);;All(*)")
            if file_path != '':
                if os.path.exists(file_path):
                    self.load_path = os.path.split(file_path)[0]
                    return file_path
            return None
        elif option == 'subject':
            file_path, file_type = QFileDialog.getOpenFileName(self, 'open file', self.load_path,
                                                               "Sequences (*.txt *.fna *.fa *.fasta *.nex "
                                                               "*.fas);;All(*)")
            if file_path != '':
                self.load_path = os.path.split(file_path)[0]
                return file_path
            return None

    def btn_open_file_query_clicked(self):
        query_file = self.get_file_name('query')
        self.label_choose_file_q.setText(query_file)

    def btn_open_file_subject_clicked(self):
        query_file = self.get_file_name('subject')
        self.label_choose_file_s.setText(query_file)

    def btn_clear_clicked(self, label):
        label.clear()

    def btn_make_blast_db_clicked(self):
        """
        创建blast数据库
        """
        self.MakeBlastDB.lineEdit_input.setText('')
        self.MakeBlastDB.lineEdit_out.setText('')
        self.MakeBlastDB.show()

    def check_box_is_checked(self):
        if self.CheckBox_two_seq.isChecked():
            self.groupBox_subject.show()
        else:
            self.groupBox_subject.hide()

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

    def _s_popen(self, command, out_file):
        """
        执行命令行命令
        :param command: cmd 命令
        :return: 执行结果
        """
        p = subprocess.Popen(command,
                             shell=True,
                             stdout=subprocess.PIPE,
                             stdin=subprocess.PIPE,
                             stderr=subprocess.STDOUT,
                             universal_newlines=True,
                             encoding='gbk',
                             )
        print(command)
        self.s_popen_pid = p.pid
        p.wait()  # 等待运行
        return_code = str(p.returncode)  # 获取cmd执行结果代码
        out = p.stdout.read()  # 返回cmd执行结果
        print(return_code)
        if 'Error' not in out:  # 成功
            self.single_finish.emit(out_file)
        elif 'Error' in out:  # 错误
            self.single_error.emit(return_code, out)

    def accept_single_finish(self, out_file):
        # 还原按钮状态
        os.system(f'start notepad "{out_file}"')
        self.pushButton_BLAST.setText('Start')
        self.pushButton_BLAST.setEnabled(True)

    def accept_single_error(self, return_code, out):
        QMessageBox.warning(self, 'error', f'error code:{return_code}\n{out}')
        self.pushButton_BLAST.setText('Start')
        self.pushButton_BLAST.setEnabled(True)

    def get_blast_command(self, flag, query, db_subject):
        """
        将blast cmd命令行组合起来
        :return: blast cmd命令
        """
        app_name = os.path.join(self.blast_path, self.comboBox_app.currentText())  # blast程序
        db_file = os.path.join(self.db_path, self.comboBox_db_files.currentText())
        outfmt = self.lineEdit_outfmt.text()  # blast输出格式
        e_value = self.lineEdit_e_value.text()  # e值
        task = self.comboBox_task.currentText()  # blast工作模式
        more_option = self.lineEdit_more_option.text()  # 更多设置选项
        if flag == 0:
            command = f'{app_name} -query "{query}" -subject "{db_subject}" -outfmt {outfmt} -evalue {e_value} -out "{self.out_file}" -task {task} {more_option}'
            if self.comboBox_app.currentText() == 'tblastx':
                command = f'{app_name} -query "{query}" -subject "{db_subject}" -outfmt {outfmt} -evalue {e_value} -out "{self.out_file}" {more_option}'
        else:
            command = f'{app_name} -query "{query}" -db "{db_file}" -outfmt {outfmt} -evalue {e_value} -out "{self.out_file}" -task {task} {more_option}'
            if self.comboBox_app.currentText() == 'tblastx':
                command = f'{app_name} -query "{query}" -db "{db_file}" -outfmt {outfmt} -evalue {e_value} -out "{self.out_file}" {more_option}'
        return command

    def blast_start(self, command):
        # print(command)
        t = Thread(target=self._s_popen, args=(command, self.out_file))
        t.setDaemon(True)
        t.start()

    def btn_blast_clicked(self):
        # query file name
        query_seq = self.input_query_edit.toPlainText()
        if query_seq == '':
            query_file = self.label_choose_file_q.text()
            if not os.path.exists(query_file):
                return QMessageBox.information(self, 'information', 'Please paste a query sequence or upload a file.')
        else:
            with open(self.query_file, 'w') as f:
                f.write(query_seq)
            query_file = self.query_file
        # subject file name  0: 双序列比对
        if self.CheckBox_two_seq.isChecked():
            subject_seq = self.input_subject_edit.toPlainText()
            if subject_seq == '':
                subject_file = self.label_choose_file_s.text()
                if not os.path.exists(subject_file):
                    return QMessageBox.information(self, 'information',
                                                   'Please paste a subject sequence or upload a file.')
            else:
                with open(self.subject_file, 'w') as f:
                    f.write(subject_seq)
                subject_file = self.subject_file
            command = self.get_blast_command(0, query_file, subject_file)

        else:   # 与本地数据库进行比对 1：db
            db = self.comboBox_db_files.currentText()
            command = self.get_blast_command(1, query_file, db)
        self.pushButton_BLAST.setText('please wait...')
        self.pushButton_BLAST.setDisabled(True)
        self.blast_start(command)


if __name__ == '__main__':
    QtCore.QCoreApplication.setAttribute(QtCore.Qt.AA_EnableHighDpiScaling)
    app = QApplication(sys.argv)
    MainWindow = LocalBlast()
    MainWindow.show()
    sys.exit(app.exec_())
