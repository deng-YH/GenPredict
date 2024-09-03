# -*- coding: UTF-8 -*-
"""
__project_ = 'blast_X32_GUI'
__file_name__ = 'make_db'
__author__ = '鲨鱼辣椒'
__time__ = '2021/5/26 15:28'
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
from uifiles import ui_blast_make_db
from threading import Thread


class MakeBlastDB(QWidget, ui_blast_make_db.Ui_Form_make_db):
    single_finish = pyqtSignal(str)  # 定义完成信号
    script_path = os.getcwd().replace('src', '')

    def __init__(self, parent=None):
        super(MakeBlastDB, self).__init__(parent)
        # 加载界面
        self.setupUi(self)
        # 初始化参数
        self.title = 'Make blast db'
        self.makeblastdb = os.path.join(self.script_path, 'tools', 'blast+', 'bin', 'makeblastdb')
        self.db_path = os.path.join(self.script_path, 'DB', 'blastdb')
        # print(self.db_path)
        self.load_path = os.getcwd()
        self.s_popen_pid = None

        # 信号和槽
        self.toolButton_input.clicked.connect(self.tool_btn_input_clicked)
        self.pushButton_build.clicked.connect(self.btn_build_clicked)
        self.single_finish.connect(self.accept_single_finish)

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
                    print(e)
            else:
                event.ignore()

    def get_file_name(self):
        """
        打开文件
        :return: open file name
        """
        file_path, file_type = QFileDialog.getOpenFileName(self, 'open file', self.load_path,
                                                           "Sequences (*.txt *.fna *.fa *.fasta *.nex "
                                                           "*.fas);;All(*)")
        if file_path != '':
            if os.path.exists(file_path):
                self.load_path = os.path.split(file_path)[0]
                return file_path
        return None

    def accept_single_finish(self, out):
        """
        更新文本框文本
        :param out: s_popen返回的结果
        :return:
        """
        self.plainTextEdit_out.appendPlainText(out)
        if 'Adding' in out:
            self.pushButton_build.setEnabled(True)
            self.pushButton_build.setText('Build')

    def s_popen(self, command):
        """
        执行命令行命令
        :param command: cmd 命令
        :return: 执行结果
        """

        # 使用Popen创建进程，并与进程进行复杂的交互
        proc = subprocess.Popen(
            command,  # cmd特定的查询空间的命令
            stdin=None,  # 标准输入 键盘
            stdout=subprocess.PIPE,  # -1 标准输出（演示器、终端) 保存到管道中以便进行操作
            stderr=subprocess.PIPE,  # 标准错误，保存到管道
            shell=True)
        self.single_finish.emit("Please wait a moment, the program is running.")
        outinfo, errinfo = proc.communicate()  # 获取输出和错误信息
        print(outinfo.decode('gbk'))  # 外部程序 (windows系统)决定编码格式
        print(errinfo.decode('gbk'))

        if not errinfo.decode('gbk'):
            self.single_finish.emit(outinfo.decode('gbk'))
        else:
            self.single_finish.emit(errinfo.decode('gbk'))

    def get_makeblastdb_command(self, in_file, dbtype, db_name):
        out = os.path.join(self.db_path,db_name)
        command = f'"{self.makeblastdb}" -in "{in_file}" -dbtype {dbtype} -out "{out}" -parse_seqids'
        print(command)
        return command

    def makeblastdb_start(self, in_file, dbtype, db_name):
        command = self.get_makeblastdb_command(in_file, dbtype, db_name)
        # 另开线程执行
        t = Thread(target=self.s_popen, args=(command,))
        t.setDaemon(True)
        t.start()

    def tool_btn_input_clicked(self):
        input_file = self.get_file_name()
        if input_file is None:
            return
        self.lineEdit_input.setText(input_file)

    def btn_build_clicked(self):
        self.plainTextEdit_out.clear()
        input_file = self.lineEdit_input.text()
        if not os.path.exists(input_file):
            return QMessageBox.information(self, 'information', 'File does not exist, please check.')
        db_name = self.lineEdit_out.text()
        if db_name == '':
            return QMessageBox.information(self, 'information', 'Please enter a database name.')
        db_type = self.comboBox.currentText()

        self.pushButton_build.setText('please wait...')
        self.pushButton_build.setDisabled(True)
        self.makeblastdb_start(input_file, db_type, db_name)


if __name__ == '__main__':
    QtCore.QCoreApplication.setAttribute(QtCore.Qt.AA_EnableHighDpiScaling)
    app = QApplication(sys.argv)
    MainWindow = MakeBlastDB()
    MainWindow.show()
    sys.exit(app.exec_())
