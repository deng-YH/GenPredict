# -*- coding: UTF-8 -*-
"""
__project_ = 'blast_X32_GUI'
__file_name__ = 'C_blast_out'
__author__ = '鲨鱼辣椒'
__time__ = '2021/5/31 0:23'
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
from Bio import SearchIO

from PyQt5.QtCore import pyqtSignal, Qt, QCoreApplication
from uifiles import ui_blast_out
from PyQt5.QtWidgets import QMessageBox, QWidget, QSplitter, QHBoxLayout, QFileDialog, QApplication
from threading import Thread
from  blastdbcmd import Blastdbcmd

class Blast_Out(QWidget, ui_blast_out.Ui_BLASTOUT):
    single_finish = pyqtSignal(str)  # 定义完成信号

    def __init__(self,  task_item=0, parent=None):
        super(Blast_Out, self).__init__(parent)
        # 加载界面
        self.setupUi(self)
        self.setWindowTitle('BlastOut task %s' % str(task_item))
        self.splitter1 = QSplitter(Qt.Horizontal)
        # 初始化参数
        # self.commandline = commandline
        self.s_popen_pid = None
        self.blastdbcmd_widget = None
        self.save_directory = os.getcwd()
        # 信号和槽
        self.pushButton_blastdbcmd.clicked.connect(self.pushButton_blastdbcmd_clicked)
        self.single_finish.connect(self.single_finish_fun)
        self.pushButton_save.clicked.connect(self.pushButton_save_clicked)

        # self.star()

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
        # print(command_line)
        p = subprocess.Popen(command_line, shell=True, stdout=subprocess.PIPE)
        self.s_popen_pid = p.pid
        stdout, stderr = p.communicate()
        out = str(stdout, encoding='gbk')
        # print(out)
        # print(p.returncode)
        if p.returncode == 0:
            # print('发出信号')
            self.single_finish.emit(out)
            self.s_popen_pid = None
        else:
            QMessageBox.warning(self, 'error', 'The program error! Please check your parameter Settings!')

    def star(self):
        self.plainTextEdit.setPlainText('The program is running, please wait...')
        t = Thread(target=self.s_popen, args=(self.commandline,))
        t.setDaemon(True)
        t.start()

    # 槽函数
    def pushButton_blastdbcmd_clicked(self):
        """
        创建Blastdbcmd对象，加入滑动条，并显示
        隐藏按钮pushButton_blastdbcmd
        """
        self.blastdbcmd_widget = Blastdbcmd()
        # 控制左右滑动

        self.splitter1.addWidget(self.plainTextEdit)
        self.splitter1.addWidget(self.blastdbcmd_widget)
        # print('111')
        windowLayout = QHBoxLayout()
        windowLayout.addWidget(self.splitter1)
        self.verticalLayout.addLayout(windowLayout)
        self.pushButton_blastdbcmd.hide()

    def single_finish_fun(self, out):
        """
        更新文本框文本
        :param out: s_popen返回的结果
        :return:
        """
        self.plainTextEdit.appendPlainText(out)

    def pushButton_save_clicked(self):
        plain_text = self.plainTextEdit.toPlainText()
        save_path, save_type = QFileDialog.getSaveFileName(self, 'Save File', self.save_directory,
                                                           "File (*.txt)")
        if save_path != '':
            self.save_directory = os.path.split(save_path)[0]
            with open(save_path, 'w') as f:
                f.write(plain_text)
if __name__ == '__main__':
    QCoreApplication.setAttribute(Qt.AA_EnableHighDpiScaling)
    app = QApplication(sys.argv)
    MainWindow = Blast_Out()
    MainWindow.show()
    sys.exit(app.exec_())











