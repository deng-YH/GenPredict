# -*- coding: utf-8 -*-
# @Time : 2021/11/1 0:26
# @Author : DengYuhang
# @Email : dengyuhang19@mails.ucas.ac.cn
# @File : DiamondMakeDB.py
# @Project : GenePredict
import subprocess
import sys

from PyQt5 import QtCore
from PyQt5.QtCore import pyqtSignal
from PyQt5.QtWidgets import QFileDialog, QWidget, QMessageBox, QApplication
import os
from uifiles import DiamondMakeDB
from threading import Thread


class MakeDB(QWidget, DiamondMakeDB.Ui_FormDiamondMakeDB):
    single_finish = pyqtSignal(str)  # 定义完成信号
    single_make_db_finished = pyqtSignal()
    script_path = os.getcwd().replace('src', '')

    def __init__(self, parent=None):
        super(MakeDB, self).__init__(parent)
        # 加载界面
        self.setupUi(self)
        # 初始化参数
        self.title = 'Make blast db'
        self.app_path = f'{self.script_path}\\tools\diamond'
        self.input_directory = os.getcwd()
        self.dbpath = f'{self.script_path}\\DB\diamonddb'
        self.s_popen_pid = None
        # 信号和槽
        self.toolButton_input.clicked.connect(self.open_files)
        self.pushButton_build.clicked.connect(self.build_db)
        self.single_finish.connect(self.single_finish_fun)

    # 重写closeEvent方法，实现窗体关闭时执行一些代码
    def closeEvent(self, event):
        """
        重写closeEvent方法，实现窗体关闭时执行一些代码
        :param event: close()触发的事件
        :return: None
        """
        if self.s_popen_pid is not None:
            reply = QMessageBox.question(self,
                                         'Local BLAST',
                                         "The program is running, do you want to close it?",
                                         QMessageBox.Yes | QMessageBox.No,
                                         QMessageBox.No)
            if reply == QMessageBox.Yes:
                event.accept()  # 接受确认关闭
                # 关闭shell任务
                try:
                    os.system("taskkill /t /f /pid %s" % self.s_popen_pid)
                except Exception as e:
                    pass
            else:
                event.ignore()

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
        self.s_popen_pid = p.pid
        while subprocess.Popen.poll(p) is None:
            stream = p.stdout.readline()
            if stream != '':
                self.single_finish.emit(stream)
                self.s_popen_pid = None
        self.single_finish.emit(f'Database construction completed!\nDB file:{self.dbpath}')
        self.single_make_db_finished.emit()

    # 槽函数
    def open_files(self):
        """open_file按钮打开文件对话框"""
        file_path, file_type = QFileDialog.getOpenFileName(self, 'open file', self.input_directory,
                                                           "All Files (*);;Sequnces (*.fasta *.fa *.txt)")
        if file_path != '':
            self.input_directory = os.path.split(file_path)[0]
            self.lineEdit_input.setText(file_path)

    def build_db(self):
        # 读取输入 输出文件
        input_file = self.lineEdit_input.text()
        out_file = os.path.join(self.dbpath, self.lineEdit_out.text())

        # 构造cmd命令
        commandline = self.app_path + f' makedb --in "{input_file}" --db "{out_file}"'

        # 另起线程执行命令
        self.plainTextEdit_out.setPlainText('The program is running, please wait...')
        self.pushButton_build.setDisabled(True)
        t = Thread(target=self.s_popen, args=(commandline,))
        t.setDaemon(True)
        t.start()

    def single_finish_fun(self, out):
        """
        更新文本框文本
        :param out: s_popen返回的结果
        :return:
        """
        self.plainTextEdit_out.appendPlainText(out)
        self.pushButton_build.setEnabled(True)


if __name__ == '__main__':
    QtCore.QCoreApplication.setAttribute(QtCore.Qt.AA_EnableHighDpiScaling)
    app = QApplication(sys.argv)
    MainWindow = MakeDB()
    MainWindow.show()
    sys.exit(app.exec_())
