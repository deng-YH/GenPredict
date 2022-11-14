# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'gene_predict.ui'
#
# Created by: PyQt5 UI code generator 5.15.4
#
# WARNING: Any manual changes made to this file will be lost when pyuic5 is
# run again.  Do not edit this file unless you know what you are doing.


from PyQt5 import QtCore, QtGui, QtWidgets


class Ui_Form_GenePredict(object):
    def setupUi(self, Form_GenePredict):
        Form_GenePredict.setObjectName("Form_GenePredict")
        Form_GenePredict.resize(635, 555)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Maximum, QtWidgets.QSizePolicy.Maximum)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(Form_GenePredict.sizePolicy().hasHeightForWidth())
        Form_GenePredict.setSizePolicy(sizePolicy)
        self.verticalLayout_8 = QtWidgets.QVBoxLayout(Form_GenePredict)
        self.verticalLayout_8.setObjectName("verticalLayout_8")
        self.verticalLayout_7 = QtWidgets.QVBoxLayout()
        self.verticalLayout_7.setObjectName("verticalLayout_7")
        self.horizontalLayout_12 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_12.setObjectName("horizontalLayout_12")
        self.plainTextEdit = QtWidgets.QPlainTextEdit(Form_GenePredict)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.plainTextEdit.sizePolicy().hasHeightForWidth())
        self.plainTextEdit.setSizePolicy(sizePolicy)
        self.plainTextEdit.setStatusTip("")
        self.plainTextEdit.setObjectName("plainTextEdit")
        self.horizontalLayout_12.addWidget(self.plainTextEdit)
        self.verticalLayout_9 = QtWidgets.QVBoxLayout()
        self.verticalLayout_9.setObjectName("verticalLayout_9")
        self.tabWidget = QtWidgets.QTabWidget(Form_GenePredict)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Maximum, QtWidgets.QSizePolicy.Ignored)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.tabWidget.sizePolicy().hasHeightForWidth())
        self.tabWidget.setSizePolicy(sizePolicy)
        self.tabWidget.setMinimumSize(QtCore.QSize(200, 0))
        self.tabWidget.setMaximumSize(QtCore.QSize(300, 400))
        self.tabWidget.setObjectName("tabWidget")
        self.tab = QtWidgets.QWidget()
        self.tab.setObjectName("tab")
        self.verticalLayout_2 = QtWidgets.QVBoxLayout(self.tab)
        self.verticalLayout_2.setObjectName("verticalLayout_2")
        self.horizontalLayout = QtWidgets.QHBoxLayout()
        self.horizontalLayout.setObjectName("horizontalLayout")
        self.label_app = QtWidgets.QLabel(self.tab)
        self.label_app.setObjectName("label_app")
        self.horizontalLayout.addWidget(self.label_app)
        self.comboBox_app = QtWidgets.QComboBox(self.tab)
        self.comboBox_app.setObjectName("comboBox_app")
        self.comboBox_app.addItem("")
        self.comboBox_app.addItem("")
        self.comboBox_app.addItem("")
        self.comboBox_app.addItem("")
        self.comboBox_app.addItem("")
        self.horizontalLayout.addWidget(self.comboBox_app)
        spacerItem = QtWidgets.QSpacerItem(40, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        self.horizontalLayout.addItem(spacerItem)
        self.label_task = QtWidgets.QLabel(self.tab)
        self.label_task.setObjectName("label_task")
        self.horizontalLayout.addWidget(self.label_task)
        self.comboBox_task = QtWidgets.QComboBox(self.tab)
        self.comboBox_task.setObjectName("comboBox_task")
        self.comboBox_task.addItem("")
        self.comboBox_task.addItem("")
        self.comboBox_task.addItem("")
        self.comboBox_task.addItem("")
        self.comboBox_task.addItem("")
        self.horizontalLayout.addWidget(self.comboBox_task)
        self.verticalLayout_2.addLayout(self.horizontalLayout)
        self.horizontalLayout_2 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_2.setObjectName("horizontalLayout_2")
        self.label_db = QtWidgets.QLabel(self.tab)
        self.label_db.setObjectName("label_db")
        self.horizontalLayout_2.addWidget(self.label_db)
        self.comboBox_db_file = QtWidgets.QComboBox(self.tab)
        self.comboBox_db_file.setObjectName("comboBox_db_file")
        self.horizontalLayout_2.addWidget(self.comboBox_db_file)
        spacerItem1 = QtWidgets.QSpacerItem(40, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        self.horizontalLayout_2.addItem(spacerItem1)
        self.verticalLayout_2.addLayout(self.horizontalLayout_2)
        self.horizontalLayout_3 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_3.setObjectName("horizontalLayout_3")
        self.label_e_value = QtWidgets.QLabel(self.tab)
        self.label_e_value.setObjectName("label_e_value")
        self.horizontalLayout_3.addWidget(self.label_e_value)
        self.lineEdit_e_value = QtWidgets.QLineEdit(self.tab)
        font = QtGui.QFont()
        font.setFamily("Times New Roman")
        self.lineEdit_e_value.setFont(font)
        self.lineEdit_e_value.setFrame(True)
        self.lineEdit_e_value.setObjectName("lineEdit_e_value")
        self.horizontalLayout_3.addWidget(self.lineEdit_e_value)
        spacerItem2 = QtWidgets.QSpacerItem(40, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        self.horizontalLayout_3.addItem(spacerItem2)
        self.verticalLayout_2.addLayout(self.horizontalLayout_3)
        self.verticalLayout = QtWidgets.QVBoxLayout()
        self.verticalLayout.setObjectName("verticalLayout")
        self.label_more_option = QtWidgets.QLabel(self.tab)
        self.label_more_option.setObjectName("label_more_option")
        self.verticalLayout.addWidget(self.label_more_option)
        self.lineEdit_more_option = QtWidgets.QLineEdit(self.tab)
        self.lineEdit_more_option.setObjectName("lineEdit_more_option")
        self.verticalLayout.addWidget(self.lineEdit_more_option)
        self.verticalLayout_2.addLayout(self.verticalLayout)
        self.pushButton_start = QtWidgets.QPushButton(self.tab)
        self.pushButton_start.setObjectName("pushButton_start")
        self.verticalLayout_2.addWidget(self.pushButton_start)
        spacerItem3 = QtWidgets.QSpacerItem(20, 40, QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Expanding)
        self.verticalLayout_2.addItem(spacerItem3)
        self.pushButton_view = QtWidgets.QPushButton(self.tab)
        self.pushButton_view.setObjectName("pushButton_view")
        self.verticalLayout_2.addWidget(self.pushButton_view)
        self.tabWidget.addTab(self.tab, "")
        self.tab_2 = QtWidgets.QWidget()
        self.tab_2.setObjectName("tab_2")
        self.verticalLayout_5 = QtWidgets.QVBoxLayout(self.tab_2)
        self.verticalLayout_5.setObjectName("verticalLayout_5")
        self.horizontalLayout_5 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_5.setObjectName("horizontalLayout_5")
        self.label_out_size = QtWidgets.QLabel(self.tab_2)
        self.label_out_size.setObjectName("label_out_size")
        self.horizontalLayout_5.addWidget(self.label_out_size)
        self.lineEdit_out_size = QtWidgets.QLineEdit(self.tab_2)
        self.lineEdit_out_size.setObjectName("lineEdit_out_size")
        self.horizontalLayout_5.addWidget(self.lineEdit_out_size)
        spacerItem4 = QtWidgets.QSpacerItem(40, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        self.horizontalLayout_5.addItem(spacerItem4)
        self.verticalLayout_5.addLayout(self.horizontalLayout_5)
        self.line = QtWidgets.QFrame(self.tab_2)
        self.line.setFrameShape(QtWidgets.QFrame.HLine)
        self.line.setFrameShadow(QtWidgets.QFrame.Sunken)
        self.line.setObjectName("line")
        self.verticalLayout_5.addWidget(self.line)
        self.horizontalLayout_6 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_6.setObjectName("horizontalLayout_6")
        self.label_Organism = QtWidgets.QLabel(self.tab_2)
        self.label_Organism.setObjectName("label_Organism")
        self.horizontalLayout_6.addWidget(self.label_Organism)
        self.comboBox_Organism = QtWidgets.QComboBox(self.tab_2)
        self.comboBox_Organism.setObjectName("comboBox_Organism")
        self.comboBox_Organism.addItem("")
        self.comboBox_Organism.addItem("")
        self.comboBox_Organism.addItem("")
        self.horizontalLayout_6.addWidget(self.comboBox_Organism)
        spacerItem5 = QtWidgets.QSpacerItem(40, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        self.horizontalLayout_6.addItem(spacerItem5)
        self.verticalLayout_5.addLayout(self.horizontalLayout_6)
        self.horizontalLayout_7 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_7.setObjectName("horizontalLayout_7")
        self.label_cutoff = QtWidgets.QLabel(self.tab_2)
        self.label_cutoff.setObjectName("label_cutoff")
        self.horizontalLayout_7.addWidget(self.label_cutoff)
        self.comboBox_cutoff = QtWidgets.QComboBox(self.tab_2)
        self.comboBox_cutoff.setObjectName("comboBox_cutoff")
        self.comboBox_cutoff.addItem("")
        self.comboBox_cutoff.addItem("")
        self.comboBox_cutoff.addItem("")
        self.comboBox_cutoff.addItem("")
        self.comboBox_cutoff.addItem("")
        self.comboBox_cutoff.addItem("")
        self.comboBox_cutoff.addItem("")
        self.horizontalLayout_7.addWidget(self.comboBox_cutoff)
        self.verticalLayout_5.addLayout(self.horizontalLayout_7)
        self.pushButton_GENSCAN = QtWidgets.QPushButton(self.tab_2)
        self.pushButton_GENSCAN.setObjectName("pushButton_GENSCAN")
        self.verticalLayout_5.addWidget(self.pushButton_GENSCAN)
        spacerItem6 = QtWidgets.QSpacerItem(20, 40, QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Expanding)
        self.verticalLayout_5.addItem(spacerItem6)
        self.pushButton_cnn = QtWidgets.QPushButton(self.tab_2)
        self.pushButton_cnn.setObjectName("pushButton_cnn")
        self.verticalLayout_5.addWidget(self.pushButton_cnn)
        self.tabWidget.addTab(self.tab_2, "")
        self.tab_3 = QtWidgets.QWidget()
        self.tab_3.setObjectName("tab_3")
        self.verticalLayout_6 = QtWidgets.QVBoxLayout(self.tab_3)
        self.verticalLayout_6.setObjectName("verticalLayout_6")
        self.horizontalLayout_10 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_10.setObjectName("horizontalLayout_10")
        self.label = QtWidgets.QLabel(self.tab_3)
        self.label.setObjectName("label")
        self.horizontalLayout_10.addWidget(self.label)
        self.comboBox_cmd_db = QtWidgets.QComboBox(self.tab_3)
        self.comboBox_cmd_db.setObjectName("comboBox_cmd_db")
        self.horizontalLayout_10.addWidget(self.comboBox_cmd_db)
        spacerItem7 = QtWidgets.QSpacerItem(40, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        self.horizontalLayout_10.addItem(spacerItem7)
        self.verticalLayout_6.addLayout(self.horizontalLayout_10)
        self.horizontalLayout_9 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_9.setObjectName("horizontalLayout_9")
        self.label_entry = QtWidgets.QLabel(self.tab_3)
        self.label_entry.setObjectName("label_entry")
        self.horizontalLayout_9.addWidget(self.label_entry)
        self.lineEdit_entry = QtWidgets.QLineEdit(self.tab_3)
        self.lineEdit_entry.setObjectName("lineEdit_entry")
        self.horizontalLayout_9.addWidget(self.lineEdit_entry)
        spacerItem8 = QtWidgets.QSpacerItem(40, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        self.horizontalLayout_9.addItem(spacerItem8)
        self.verticalLayout_6.addLayout(self.horizontalLayout_9)
        self.groupBox = QtWidgets.QGroupBox(self.tab_3)
        self.groupBox.setObjectName("groupBox")
        self.verticalLayout_4 = QtWidgets.QVBoxLayout(self.groupBox)
        self.verticalLayout_4.setObjectName("verticalLayout_4")
        self.verticalLayout_3 = QtWidgets.QVBoxLayout()
        self.verticalLayout_3.setObjectName("verticalLayout_3")
        self.radioButton_whole = QtWidgets.QRadioButton(self.groupBox)
        self.radioButton_whole.setChecked(True)
        self.radioButton_whole.setObjectName("radioButton_whole")
        self.verticalLayout_3.addWidget(self.radioButton_whole)
        self.radioButton_selfected = QtWidgets.QRadioButton(self.groupBox)
        self.radioButton_selfected.setObjectName("radioButton_selfected")
        self.verticalLayout_3.addWidget(self.radioButton_selfected)
        self.verticalLayout_4.addLayout(self.verticalLayout_3)
        self.horizontalLayout_8 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_8.setObjectName("horizontalLayout_8")
        self.label_3 = QtWidgets.QLabel(self.groupBox)
        self.label_3.setObjectName("label_3")
        self.horizontalLayout_8.addWidget(self.label_3)
        self.lineEdit_begin = QtWidgets.QLineEdit(self.groupBox)
        self.lineEdit_begin.setObjectName("lineEdit_begin")
        self.horizontalLayout_8.addWidget(self.lineEdit_begin)
        self.label_4 = QtWidgets.QLabel(self.groupBox)
        self.label_4.setObjectName("label_4")
        self.horizontalLayout_8.addWidget(self.label_4)
        self.lineEdit_end = QtWidgets.QLineEdit(self.groupBox)
        self.lineEdit_end.setObjectName("lineEdit_end")
        self.horizontalLayout_8.addWidget(self.lineEdit_end)
        self.verticalLayout_4.addLayout(self.horizontalLayout_8)
        self.verticalLayout_6.addWidget(self.groupBox)
        self.groupBox_2 = QtWidgets.QGroupBox(self.tab_3)
        self.groupBox_2.setTitle("")
        self.groupBox_2.setObjectName("groupBox_2")
        self.horizontalLayout_4 = QtWidgets.QHBoxLayout(self.groupBox_2)
        self.horizontalLayout_4.setObjectName("horizontalLayout_4")
        self.radioButton_seq = QtWidgets.QRadioButton(self.groupBox_2)
        self.radioButton_seq.setChecked(True)
        self.radioButton_seq.setObjectName("radioButton_seq")
        self.horizontalLayout_4.addWidget(self.radioButton_seq)
        self.radioButton_reverse_complement = QtWidgets.QRadioButton(self.groupBox_2)
        self.radioButton_reverse_complement.setObjectName("radioButton_reverse_complement")
        self.horizontalLayout_4.addWidget(self.radioButton_reverse_complement)
        self.verticalLayout_6.addWidget(self.groupBox_2)
        self.horizontalLayout_11 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_11.setObjectName("horizontalLayout_11")
        spacerItem9 = QtWidgets.QSpacerItem(40, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        self.horizontalLayout_11.addItem(spacerItem9)
        self.pushButton_create_file = QtWidgets.QPushButton(self.tab_3)
        self.pushButton_create_file.setObjectName("pushButton_create_file")
        self.horizontalLayout_11.addWidget(self.pushButton_create_file)
        self.label_5 = QtWidgets.QLabel(self.tab_3)
        self.label_5.setObjectName("label_5")
        self.horizontalLayout_11.addWidget(self.label_5)
        self.pushButton_GENSCAN_2 = QtWidgets.QPushButton(self.tab_3)
        self.pushButton_GENSCAN_2.setObjectName("pushButton_GENSCAN_2")
        self.horizontalLayout_11.addWidget(self.pushButton_GENSCAN_2)
        spacerItem10 = QtWidgets.QSpacerItem(40, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        self.horizontalLayout_11.addItem(spacerItem10)
        self.verticalLayout_6.addLayout(self.horizontalLayout_11)
        self.tabWidget.addTab(self.tab_3, "")
        self.verticalLayout_9.addWidget(self.tabWidget)
        self.horizontalLayout_12.addLayout(self.verticalLayout_9)
        self.verticalLayout_7.addLayout(self.horizontalLayout_12)
        self.tableWidget = QtWidgets.QTableWidget(Form_GenePredict)
        self.tableWidget.setObjectName("tableWidget")
        self.tableWidget.setColumnCount(12)
        self.tableWidget.setRowCount(0)
        item = QtWidgets.QTableWidgetItem()
        self.tableWidget.setHorizontalHeaderItem(0, item)
        item = QtWidgets.QTableWidgetItem()
        self.tableWidget.setHorizontalHeaderItem(1, item)
        item = QtWidgets.QTableWidgetItem()
        self.tableWidget.setHorizontalHeaderItem(2, item)
        item = QtWidgets.QTableWidgetItem()
        self.tableWidget.setHorizontalHeaderItem(3, item)
        item = QtWidgets.QTableWidgetItem()
        self.tableWidget.setHorizontalHeaderItem(4, item)
        item = QtWidgets.QTableWidgetItem()
        self.tableWidget.setHorizontalHeaderItem(5, item)
        item = QtWidgets.QTableWidgetItem()
        self.tableWidget.setHorizontalHeaderItem(6, item)
        item = QtWidgets.QTableWidgetItem()
        self.tableWidget.setHorizontalHeaderItem(7, item)
        item = QtWidgets.QTableWidgetItem()
        self.tableWidget.setHorizontalHeaderItem(8, item)
        item = QtWidgets.QTableWidgetItem()
        self.tableWidget.setHorizontalHeaderItem(9, item)
        item = QtWidgets.QTableWidgetItem()
        self.tableWidget.setHorizontalHeaderItem(10, item)
        item = QtWidgets.QTableWidgetItem()
        self.tableWidget.setHorizontalHeaderItem(11, item)
        self.verticalLayout_7.addWidget(self.tableWidget)
        self.verticalLayout_8.addLayout(self.verticalLayout_7)

        self.retranslateUi(Form_GenePredict)
        self.tabWidget.setCurrentIndex(0)
        self.lineEdit_begin.textChanged['QString'].connect(self.radioButton_selfected.click)
        self.lineEdit_end.textChanged['QString'].connect(self.radioButton_selfected.click)
        QtCore.QMetaObject.connectSlotsByName(Form_GenePredict)

    def retranslateUi(self, Form_GenePredict):
        _translate = QtCore.QCoreApplication.translate
        Form_GenePredict.setWindowTitle(_translate("Form_GenePredict", "GenePredict"))
        self.plainTextEdit.setPlaceholderText(_translate("Form_GenePredict", "Enter the homologous sequence here"))
        self.label_app.setText(_translate("Form_GenePredict", "APP:"))
        self.comboBox_app.setItemText(0, _translate("Form_GenePredict", "blastn"))
        self.comboBox_app.setItemText(1, _translate("Form_GenePredict", "blastp"))
        self.comboBox_app.setItemText(2, _translate("Form_GenePredict", "blastx"))
        self.comboBox_app.setItemText(3, _translate("Form_GenePredict", "tblastn"))
        self.comboBox_app.setItemText(4, _translate("Form_GenePredict", "tblastx"))
        self.label_task.setText(_translate("Form_GenePredict", "-task:"))
        self.comboBox_task.setItemText(0, _translate("Form_GenePredict", "megablast"))
        self.comboBox_task.setItemText(1, _translate("Form_GenePredict", "dc-megablast"))
        self.comboBox_task.setItemText(2, _translate("Form_GenePredict", "blastn"))
        self.comboBox_task.setItemText(3, _translate("Form_GenePredict", "blastn-short"))
        self.comboBox_task.setItemText(4, _translate("Form_GenePredict", "rmblastn"))
        self.label_db.setText(_translate("Form_GenePredict", "-db:"))
        self.label_e_value.setText(_translate("Form_GenePredict", "e-value:"))
        self.lineEdit_e_value.setText(_translate("Form_GenePredict", "0.01"))
        self.label_more_option.setText(_translate("Form_GenePredict", "more optin:"))
        self.pushButton_start.setText(_translate("Form_GenePredict", "Start"))
        self.pushButton_view.setText(_translate("Form_GenePredict", "View"))
        self.tabWidget.setTabText(self.tabWidget.indexOf(self.tab), _translate("Form_GenePredict", "blast"))
        self.label_out_size.setText(_translate("Form_GenePredict", "out size(bp):"))
        self.lineEdit_out_size.setText(_translate("Form_GenePredict", "500"))
        self.label_Organism.setText(_translate("Form_GenePredict", "Organism:"))
        self.comboBox_Organism.setItemText(0, _translate("Form_GenePredict", "Vertebrate"))
        self.comboBox_Organism.setItemText(1, _translate("Form_GenePredict", "Arabidopsis"))
        self.comboBox_Organism.setItemText(2, _translate("Form_GenePredict", "Maize"))
        self.label_cutoff.setText(_translate("Form_GenePredict", "Suboptimal exon cutoff (optional):"))
        self.comboBox_cutoff.setItemText(0, _translate("Form_GenePredict", "1.00"))
        self.comboBox_cutoff.setItemText(1, _translate("Form_GenePredict", "0.50"))
        self.comboBox_cutoff.setItemText(2, _translate("Form_GenePredict", "0.25"))
        self.comboBox_cutoff.setItemText(3, _translate("Form_GenePredict", "0.10"))
        self.comboBox_cutoff.setItemText(4, _translate("Form_GenePredict", "0.05"))
        self.comboBox_cutoff.setItemText(5, _translate("Form_GenePredict", "0.02"))
        self.comboBox_cutoff.setItemText(6, _translate("Form_GenePredict", "0.01"))
        self.pushButton_GENSCAN.setText(_translate("Form_GenePredict", "→GENSCAN Web Server"))
        self.pushButton_cnn.setText(_translate("Form_GenePredict", "CNN"))
        self.tabWidget.setTabText(self.tabWidget.indexOf(self.tab_2), _translate("Form_GenePredict", "genscan"))
        self.label.setText(_translate("Form_GenePredict", "   db:"))
        self.label_entry.setText(_translate("Form_GenePredict", "entry:"))
        self.groupBox.setTitle(_translate("Form_GenePredict", "Change region shown"))
        self.radioButton_whole.setText(_translate("Form_GenePredict", "Whole sequence"))
        self.radioButton_selfected.setText(_translate("Form_GenePredict", "Selected region"))
        self.label_3.setText(_translate("Form_GenePredict", "from"))
        self.lineEdit_begin.setPlaceholderText(_translate("Form_GenePredict", "begin"))
        self.label_4.setText(_translate("Form_GenePredict", "to"))
        self.lineEdit_end.setPlaceholderText(_translate("Form_GenePredict", "end"))
        self.radioButton_seq.setText(_translate("Form_GenePredict", "sequence"))
        self.radioButton_reverse_complement.setText(_translate("Form_GenePredict", "reverse complement"))
        self.pushButton_create_file.setText(_translate("Form_GenePredict", "Create File"))
        self.label_5.setText(_translate("Form_GenePredict", "or"))
        self.pushButton_GENSCAN_2.setText(_translate("Form_GenePredict", "→GENSCAN Web Server"))
        self.tabWidget.setTabText(self.tabWidget.indexOf(self.tab_3), _translate("Form_GenePredict", "blastdbcmd"))
        self.tableWidget.setSortingEnabled(True)
        item = self.tableWidget.horizontalHeaderItem(0)
        item.setText(_translate("Form_GenePredict", "qseqid"))
        item = self.tableWidget.horizontalHeaderItem(1)
        item.setText(_translate("Form_GenePredict", "sseqid"))
        item = self.tableWidget.horizontalHeaderItem(2)
        item.setText(_translate("Form_GenePredict", "pident"))
        item = self.tableWidget.horizontalHeaderItem(3)
        item.setText(_translate("Form_GenePredict", "length"))
        item = self.tableWidget.horizontalHeaderItem(4)
        item.setText(_translate("Form_GenePredict", "mismatch"))
        item = self.tableWidget.horizontalHeaderItem(5)
        item.setText(_translate("Form_GenePredict", "gapopen"))
        item = self.tableWidget.horizontalHeaderItem(6)
        item.setText(_translate("Form_GenePredict", "qstart"))
        item = self.tableWidget.horizontalHeaderItem(7)
        item.setText(_translate("Form_GenePredict", "qend"))
        item = self.tableWidget.horizontalHeaderItem(8)
        item.setText(_translate("Form_GenePredict", "sstart"))
        item = self.tableWidget.horizontalHeaderItem(9)
        item.setText(_translate("Form_GenePredict", "send"))
        item = self.tableWidget.horizontalHeaderItem(10)
        item.setText(_translate("Form_GenePredict", "evalue"))
        item = self.tableWidget.horizontalHeaderItem(11)
        item.setText(_translate("Form_GenePredict", "bitscore"))
