# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'ui_blastdbcmd.ui'
#
# Created by: PyQt5 UI code generator 5.15.4
#
# WARNING: Any manual changes made to this file will be lost when pyuic5 is
# run again.  Do not edit this file unless you know what you are doing.


from PyQt5 import QtCore, QtGui, QtWidgets


class Ui_Frame_Blastdbcmd(object):
    def setupUi(self, Frame_Blastdbcmd):
        Frame_Blastdbcmd.setObjectName("Frame_Blastdbcmd")
        Frame_Blastdbcmd.resize(313, 246)
        self.verticalLayout = QtWidgets.QVBoxLayout(Frame_Blastdbcmd)
        self.verticalLayout.setObjectName("verticalLayout")
        self.horizontalLayout_10 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_10.setObjectName("horizontalLayout_10")
        self.label = QtWidgets.QLabel(Frame_Blastdbcmd)
        self.label.setObjectName("label")
        self.horizontalLayout_10.addWidget(self.label)
        self.comboBox_cmd_db = QtWidgets.QComboBox(Frame_Blastdbcmd)
        self.comboBox_cmd_db.setObjectName("comboBox_cmd_db")
        self.horizontalLayout_10.addWidget(self.comboBox_cmd_db)
        spacerItem = QtWidgets.QSpacerItem(40, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        self.horizontalLayout_10.addItem(spacerItem)
        self.verticalLayout.addLayout(self.horizontalLayout_10)
        self.horizontalLayout_9 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_9.setObjectName("horizontalLayout_9")
        self.label_entry = QtWidgets.QLabel(Frame_Blastdbcmd)
        self.label_entry.setObjectName("label_entry")
        self.horizontalLayout_9.addWidget(self.label_entry)
        self.lineEdit_entry = QtWidgets.QLineEdit(Frame_Blastdbcmd)
        self.lineEdit_entry.setObjectName("lineEdit_entry")
        self.horizontalLayout_9.addWidget(self.lineEdit_entry)
        spacerItem1 = QtWidgets.QSpacerItem(40, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        self.horizontalLayout_9.addItem(spacerItem1)
        self.verticalLayout.addLayout(self.horizontalLayout_9)
        self.groupBox = QtWidgets.QGroupBox(Frame_Blastdbcmd)
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
        self.verticalLayout.addWidget(self.groupBox)
        self.groupBox_2 = QtWidgets.QGroupBox(Frame_Blastdbcmd)
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
        self.verticalLayout.addWidget(self.groupBox_2)
        self.horizontalLayout_11 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_11.setObjectName("horizontalLayout_11")
        spacerItem2 = QtWidgets.QSpacerItem(40, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        self.horizontalLayout_11.addItem(spacerItem2)
        self.pushButton_create_file = QtWidgets.QPushButton(Frame_Blastdbcmd)
        self.pushButton_create_file.setObjectName("pushButton_create_file")
        self.horizontalLayout_11.addWidget(self.pushButton_create_file)
        spacerItem3 = QtWidgets.QSpacerItem(40, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        self.horizontalLayout_11.addItem(spacerItem3)
        self.verticalLayout.addLayout(self.horizontalLayout_11)

        self.retranslateUi(Frame_Blastdbcmd)
        QtCore.QMetaObject.connectSlotsByName(Frame_Blastdbcmd)

    def retranslateUi(self, Frame_Blastdbcmd):
        _translate = QtCore.QCoreApplication.translate
        Frame_Blastdbcmd.setWindowTitle(_translate("Frame_Blastdbcmd", "Blastdbcmd"))
        self.label.setText(_translate("Frame_Blastdbcmd", "   db:"))
        self.label_entry.setText(_translate("Frame_Blastdbcmd", "entry:"))
        self.groupBox.setTitle(_translate("Frame_Blastdbcmd", "Change region shown"))
        self.radioButton_whole.setText(_translate("Frame_Blastdbcmd", "Whole sequence"))
        self.radioButton_selfected.setText(_translate("Frame_Blastdbcmd", "Selected region"))
        self.label_3.setText(_translate("Frame_Blastdbcmd", "from"))
        self.lineEdit_begin.setPlaceholderText(_translate("Frame_Blastdbcmd", "begin"))
        self.label_4.setText(_translate("Frame_Blastdbcmd", "to"))
        self.lineEdit_end.setPlaceholderText(_translate("Frame_Blastdbcmd", "end"))
        self.radioButton_seq.setText(_translate("Frame_Blastdbcmd", "sequence"))
        self.radioButton_reverse_complement.setText(_translate("Frame_Blastdbcmd", "reverse complement"))
        self.pushButton_create_file.setText(_translate("Frame_Blastdbcmd", "Create File"))
