# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'BlastX.ui'
#
# Created by: PyQt5 UI code generator 5.15.4
#
# WARNING: Any manual changes made to this file will be lost when pyuic5 is
# run again.  Do not edit this file unless you know what you are doing.


from PyQt5 import QtCore, QtGui, QtWidgets


class Ui_FormBlastX(object):
    def setupUi(self, FormBlastX):
        FormBlastX.setObjectName("FormBlastX")
        FormBlastX.resize(764, 485)
        self.verticalLayout = QtWidgets.QVBoxLayout(FormBlastX)
        self.verticalLayout.setObjectName("verticalLayout")
        self.horizontalLayout = QtWidgets.QHBoxLayout()
        self.horizontalLayout.setObjectName("horizontalLayout")
        self.label_file = QtWidgets.QLabel(FormBlastX)
        self.label_file.setObjectName("label_file")
        self.horizontalLayout.addWidget(self.label_file)
        self.lineEdit_input = QtWidgets.QLineEdit(FormBlastX)
        self.lineEdit_input.setObjectName("lineEdit_input")
        self.horizontalLayout.addWidget(self.lineEdit_input)
        self.toolButton_input = QtWidgets.QToolButton(FormBlastX)
        self.toolButton_input.setObjectName("toolButton_input")
        self.horizontalLayout.addWidget(self.toolButton_input)
        self.verticalLayout.addLayout(self.horizontalLayout)
        self.horizontalLayout_2 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_2.setObjectName("horizontalLayout_2")
        self.label_dbfile = QtWidgets.QLabel(FormBlastX)
        self.label_dbfile.setObjectName("label_dbfile")
        self.horizontalLayout_2.addWidget(self.label_dbfile)
        self.comboBox = QtWidgets.QComboBox(FormBlastX)
        self.comboBox.setObjectName("comboBox")
        self.horizontalLayout_2.addWidget(self.comboBox)
        spacerItem = QtWidgets.QSpacerItem(40, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        self.horizontalLayout_2.addItem(spacerItem)
        self.Button_makedb = QtWidgets.QPushButton(FormBlastX)
        self.Button_makedb.setObjectName("Button_makedb")
        self.horizontalLayout_2.addWidget(self.Button_makedb)
        self.verticalLayout.addLayout(self.horizontalLayout_2)
        self.plainTextEdit = QtWidgets.QPlainTextEdit(FormBlastX)
        self.plainTextEdit.setObjectName("plainTextEdit")
        self.verticalLayout.addWidget(self.plainTextEdit)
        self.ButtonStart = QtWidgets.QPushButton(FormBlastX)
        self.ButtonStart.setObjectName("ButtonStart")
        self.verticalLayout.addWidget(self.ButtonStart)

        self.retranslateUi(FormBlastX)
        QtCore.QMetaObject.connectSlotsByName(FormBlastX)

    def retranslateUi(self, FormBlastX):
        _translate = QtCore.QCoreApplication.translate
        FormBlastX.setWindowTitle(_translate("FormBlastX", "Form"))
        self.label_file.setText(_translate("FormBlastX", "Input file:"))
        self.toolButton_input.setText(_translate("FormBlastX", "..."))
        self.label_dbfile.setText(_translate("FormBlastX", "DB file:"))
        self.Button_makedb.setText(_translate("FormBlastX", "MakeDB"))
        self.ButtonStart.setText(_translate("FormBlastX", "Start"))
