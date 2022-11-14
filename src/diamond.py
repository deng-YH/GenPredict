# -*- coding: UTF-8 -*-
"""
__project_ = 'GenePredict'
__file_name__ = 'diamond'
__author__ = 'Dyh'
__time__ = '2021/11/14 12:58'
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
import subprocess
import sys

from PyQt5 import QtCore
from PyQt5.QtCore import pyqtSignal
from PyQt5.QtWidgets import QFileDialog, QWidget, QMessageBox, QApplication
import os
from uifiles import ui_diamond
from threading import Thread
from src import DiamondMakeDB


class Diamond(QWidget, ui_diamond.Ui_Frame_Diamond):
    single_finish = pyqtSignal(str)  # 定义完成信号
    single_error = pyqtSignal(int, str)
    script_path = os.getcwd().replace('src', '')

    def __init__(self, parent=None):
        super(Diamond, self).__init__(parent)
        # 加载界面
        self.setupUi(self)

        # 初始化
        self.diamond_tool = f'{self.script_path}\\tools\diamond'
        self.db_path = f'{self.script_path}\\DB\diamonddb'
        self.input_file = f'{self.script_path}\\tempfile\diamond_input_seq.fasta'
        self.input_path = os.getcwd()

        # 存储变量
        self.out_path = ''
        self.s_popen_pid = None
        self.query_id_list = []
        self.subject_id_list = []
        self.subject_fasta = ''
        self.query_subject_des = []

        # 实例对象
        self.DiamondMakeDB = DiamondMakeDB.MakeDB()

        # 加载diamond_db_files
        self.comboBox_db_files_addItems()

        # 信号和槽
        self.toolButton_input_file.clicked.connect(self.tool_btn_input_clicked)
        self.toolButton_out_file.clicked.connect(self.tool_btn_out_clicked)
        self.pushButton_makedb.clicked.connect(self.btn_make_db_clicked)
        self.pushButton_start.clicked.connect(self.btn_start_clicked)
        self.single_finish.connect(self.accept_single_finish)
        self.single_error.connect(self.accept_single_error)

    def comboBox_db_files_addItems(self):
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

    def get_file_name(self, option):
        """
        QFileDialog.get???FileName()
        :param option: open or save
        :return: open file name or save file name
        """
        if option == 'open':
            file_path, file_type = QFileDialog.getOpenFileName(self, 'open file', self.input_path,
                                                               "Sequences (*.txt *.fna *.fa *.fasta *.nex "
                                                               "*.fas);;All(*)")
            if file_path != '':
                if os.path.exists(file_path):
                    self.input_path = os.path.split(file_path)[0]
                    return file_path
            return None
        elif option == 'save':
            file_path, file_type = QFileDialog.getSaveFileName(self, 'save file', self.input_path, "*.txt")
            # print(file_path)
            if file_path != '':
                self.input_path = os.path.split(file_path)[0]
                return file_path
            return None

    def tool_btn_input_clicked(self):
        input_file = self.get_file_name('open')
        if input_file is not None:
            self.label_file_path.setText(input_file)

    def tool_btn_out_clicked(self):
        out_file = self.get_file_name('save')
        # print(out_file)
        if out_file is not None:
            self.lineEdit_out_file.setText(out_file)

    def btn_make_db_clicked(self):
        self.DiamondMakeDB.show()
        self.DiamondMakeDB.single_make_db_finished.connect(self.comboBox_db_files_addItems)

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
        self.s_popen_pid = p.pid
        p.wait()  # 等待运行
        return_code = str(p.returncode)  # 获取cmd执行结果代码
        out = p.stdout.read()  # 返回cmd执行结果
        # print('out:', out)
        if 'Error' not in out:  # 成功
            self.single_finish.emit(out_file)
        elif 'Error' in out:  # 错误
            self.single_error.emit(return_code, out)

    def accept_single_finish(self, out_file):
        # 还原按钮状态
        os.system(f'start notepad "{out_file}"')
        self.pushButton_start.setText('Start')
        self.pushButton_start.setEnabled(True)

    def accept_single_error(self, return_code, out):
        QMessageBox.warning(self, 'error', f'error code:{return_code}\n{out}')
        self.pushButton_start.setText('Start')
        self.pushButton_start.setEnabled(True)

    def _get_diamond_options(self, diamond_app, input_file, out_file, db_file, k, e_value, out_information):
        """
        构造diamond 命令行
        :return: diamond 命令行
        """
        option_app = f'{self.diamond_tool} {diamond_app}'
        command = f'{option_app} --query "{input_file}" --db "{db_file}" -k {k} -e {e_value} --out "{out_file}" ' \
                  f'--outfmt 6 {out_information}'
        return command

    def diamond_start(self, diamond_app, input_file, out_file, db_file, k, e_value, out_information):
        command = self._get_diamond_options(diamond_app, input_file, out_file, db_file, k, e_value, out_information)
        # self._s_popen(command, out_file)
        # 另开线程执行
        t = Thread(target=self._s_popen, args=(command, out_file))
        t.setDaemon(True)
        t.start()

    def btn_start_clicked(self):
        # input输入判断
        input_seq = self.plainTextEdit.toPlainText()
        if input_seq == '':
            input_file = self.label_file_path.text()
            if not os.path.exists(input_file):
                return QMessageBox.information(self, 'information', 'Please paste a sequence or upload a file.')
        else:
            with open(self.input_file, 'w') as f:
                f.write(input_seq)
            input_file = self.input_file
        out_file = self.lineEdit_out_file.text()
        # 输出文件
        if out_file == '':
            return QMessageBox.information(self, 'information', 'Please enter the correct save file.')
        # 读取参数
        diamond_app = self.comboBox_app.currentText()
        db_file = os.path.join(self.db_path, self.comboBox_db_files.currentText())
        # print(db_file)
        k = self.lineEdit_k.text()
        e_value = self.lineEdit_e.text()
        out_information = self.lineEdit_output.text()
        # 执行cmd命令
        self.pushButton_start.setText('please wait...')
        self.pushButton_start.setDisabled(True)
        self.diamond_start(diamond_app, input_file, out_file, db_file, k, e_value, out_information)


if __name__ == '__main__':
    QtCore.QCoreApplication.setAttribute(QtCore.Qt.AA_EnableHighDpiScaling)
    app = QApplication(sys.argv)
    MainWindow = Diamond()
    MainWindow.show()
    sys.exit(app.exec_())
