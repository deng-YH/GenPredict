# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'Translate.ui'
#
# Created by: PyQt5 UI code generator 5.15.4
#
# WARNING: Any manual changes made to this file will be lost when pyuic5 is
# run again.  Do not edit this file unless you know what you are doing.


from PyQt5 import QtCore, QtGui, QtWidgets


class Ui_Translate(object):
    def setupUi(self, Translate):
        Translate.setObjectName("Translate")
        Translate.resize(615, 468)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Maximum, QtWidgets.QSizePolicy.Preferred)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(Translate.sizePolicy().hasHeightForWidth())
        Translate.setSizePolicy(sizePolicy)
        Translate.setMinimumSize(QtCore.QSize(0, 0))
        Translate.setMaximumSize(QtCore.QSize(435345, 142432))
        self.verticalLayout_3 = QtWidgets.QVBoxLayout(Translate)
        self.verticalLayout_3.setObjectName("verticalLayout_3")
        self.widget = QtWidgets.QWidget(Translate)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Preferred)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.widget.sizePolicy().hasHeightForWidth())
        self.widget.setSizePolicy(sizePolicy)
        self.widget.setObjectName("widget")
        self.verticalLayout_2 = QtWidgets.QVBoxLayout(self.widget)
        self.verticalLayout_2.setObjectName("verticalLayout_2")
        self.verticalLayout = QtWidgets.QVBoxLayout()
        self.verticalLayout.setSizeConstraint(QtWidgets.QLayout.SetFixedSize)
        self.verticalLayout.setObjectName("verticalLayout")
        self.label = QtWidgets.QLabel(self.widget)
        self.label.setFrameShape(QtWidgets.QFrame.NoFrame)
        self.label.setObjectName("label")
        self.verticalLayout.addWidget(self.label)
        self.plainTextEdit = QtWidgets.QPlainTextEdit(self.widget)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Preferred)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.plainTextEdit.sizePolicy().hasHeightForWidth())
        self.plainTextEdit.setSizePolicy(sizePolicy)
        self.plainTextEdit.setMaximumSize(QtCore.QSize(16777215, 500))
        self.plainTextEdit.setObjectName("plainTextEdit")
        self.verticalLayout.addWidget(self.plainTextEdit)
        self.horizontalLayout = QtWidgets.QHBoxLayout()
        self.horizontalLayout.setObjectName("horizontalLayout")
        self.comboBox_genetic_codes = QtWidgets.QComboBox(self.widget)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.MinimumExpanding, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.comboBox_genetic_codes.sizePolicy().hasHeightForWidth())
        self.comboBox_genetic_codes.setSizePolicy(sizePolicy)
        self.comboBox_genetic_codes.setMinimumSize(QtCore.QSize(350, 0))
        self.comboBox_genetic_codes.setMaximumSize(QtCore.QSize(300, 16777215))
        self.comboBox_genetic_codes.setObjectName("comboBox_genetic_codes")
        self.comboBox_genetic_codes.addItem("")
        self.comboBox_genetic_codes.addItem("")
        self.comboBox_genetic_codes.addItem("")
        self.comboBox_genetic_codes.addItem("")
        self.comboBox_genetic_codes.addItem("")
        self.comboBox_genetic_codes.addItem("")
        self.comboBox_genetic_codes.addItem("")
        self.comboBox_genetic_codes.addItem("")
        self.comboBox_genetic_codes.addItem("")
        self.comboBox_genetic_codes.addItem("")
        self.comboBox_genetic_codes.addItem("")
        self.comboBox_genetic_codes.addItem("")
        self.comboBox_genetic_codes.addItem("")
        self.comboBox_genetic_codes.addItem("")
        self.comboBox_genetic_codes.addItem("")
        self.comboBox_genetic_codes.addItem("")
        self.comboBox_genetic_codes.addItem("")
        self.comboBox_genetic_codes.addItem("")
        self.comboBox_genetic_codes.addItem("")
        self.comboBox_genetic_codes.addItem("")
        self.comboBox_genetic_codes.addItem("")
        self.comboBox_genetic_codes.addItem("")
        self.comboBox_genetic_codes.addItem("")
        self.comboBox_genetic_codes.addItem("")
        self.comboBox_genetic_codes.addItem("")
        self.horizontalLayout.addWidget(self.comboBox_genetic_codes)
        spacerItem = QtWidgets.QSpacerItem(40, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        self.horizontalLayout.addItem(spacerItem)
        self.pushButton_reset = QtWidgets.QPushButton(self.widget)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.pushButton_reset.sizePolicy().hasHeightForWidth())
        self.pushButton_reset.setSizePolicy(sizePolicy)
        self.pushButton_reset.setObjectName("pushButton_reset")
        self.horizontalLayout.addWidget(self.pushButton_reset)
        self.pushButton_translate = QtWidgets.QPushButton(self.widget)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.pushButton_translate.sizePolicy().hasHeightForWidth())
        self.pushButton_translate.setSizePolicy(sizePolicy)
        self.pushButton_translate.setObjectName("pushButton_translate")
        self.horizontalLayout.addWidget(self.pushButton_translate)
        self.verticalLayout.addLayout(self.horizontalLayout)
        self.label_2 = QtWidgets.QLabel(self.widget)
        self.label_2.setObjectName("label_2")
        self.verticalLayout.addWidget(self.label_2)
        self.textEdit = QtWidgets.QTextEdit(self.widget)
        self.textEdit.setObjectName("textEdit")
        self.verticalLayout.addWidget(self.textEdit)
        self.verticalLayout_2.addLayout(self.verticalLayout)
        self.verticalLayout_3.addWidget(self.widget)

        self.retranslateUi(Translate)
        self.pushButton_reset.clicked.connect(self.plainTextEdit.clear)
        QtCore.QMetaObject.connectSlotsByName(Translate)

    def retranslateUi(self, Translate):
        _translate = QtCore.QCoreApplication.translate
        Translate.setWindowTitle(_translate("Translate", "Translate"))
        self.label.setText(_translate("Translate", "Enter Sequence:"))
        self.comboBox_genetic_codes.setItemText(0, _translate("Translate", "Standard"))
        self.comboBox_genetic_codes.setItemText(1, _translate("Translate", "Vertebrate Mitochondrial"))
        self.comboBox_genetic_codes.setItemText(2, _translate("Translate", "Yeast Mitochondrial"))
        self.comboBox_genetic_codes.setItemText(3, _translate("Translate", "Mold, Protozoan, and Coelenterate Mitochondrial Code and the Mycoplasma/Spiroplasma"))
        self.comboBox_genetic_codes.setItemText(4, _translate("Translate", "Invertebrate Mitochondrial"))
        self.comboBox_genetic_codes.setItemText(5, _translate("Translate", "Ciliate, Dasycladacean and Hexamita Nuclear"))
        self.comboBox_genetic_codes.setItemText(6, _translate("Translate", "Echinoderm and Flatworm Mitochondrial"))
        self.comboBox_genetic_codes.setItemText(7, _translate("Translate", "Euplotid Nuclear"))
        self.comboBox_genetic_codes.setItemText(8, _translate("Translate", "Bacterial, Archaeal and Plant Plastid"))
        self.comboBox_genetic_codes.setItemText(9, _translate("Translate", "Alternative Yeast Nuclear"))
        self.comboBox_genetic_codes.setItemText(10, _translate("Translate", "Ascidian Mitochondrial"))
        self.comboBox_genetic_codes.setItemText(11, _translate("Translate", "Alternative Flatworm Mitochondrial"))
        self.comboBox_genetic_codes.setItemText(12, _translate("Translate", "Chlorophycean Mitochondrial"))
        self.comboBox_genetic_codes.setItemText(13, _translate("Translate", "Trematode Mitochondrial"))
        self.comboBox_genetic_codes.setItemText(14, _translate("Translate", "Scenedesmus obliquus Mitochondrial"))
        self.comboBox_genetic_codes.setItemText(15, _translate("Translate", "Thraustochytrium Mitochondrial"))
        self.comboBox_genetic_codes.setItemText(16, _translate("Translate", "Rhabdopleuridae Mitochondrial"))
        self.comboBox_genetic_codes.setItemText(17, _translate("Translate", "Candidate Division SR1 and Gracilibacteria"))
        self.comboBox_genetic_codes.setItemText(18, _translate("Translate", "Pachysolen tannophilus Nuclear"))
        self.comboBox_genetic_codes.setItemText(19, _translate("Translate", "Karyorelict Nuclear"))
        self.comboBox_genetic_codes.setItemText(20, _translate("Translate", "Condylostoma Nuclear"))
        self.comboBox_genetic_codes.setItemText(21, _translate("Translate", "Mesodinium Nuclear"))
        self.comboBox_genetic_codes.setItemText(22, _translate("Translate", "Peritrich Nuclear"))
        self.comboBox_genetic_codes.setItemText(23, _translate("Translate", "Blastocrithidia Nuclear"))
        self.comboBox_genetic_codes.setItemText(24, _translate("Translate", "Cephalodiscidae Mitochondrial UAA-Tyr"))
        self.pushButton_reset.setText(_translate("Translate", "reset"))
        self.pushButton_translate.setText(_translate("Translate", "Translate"))
        self.label_2.setText(_translate("Translate", "Results of translation:"))
        self.textEdit.setHtml(_translate("Translate", "<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.0//EN\" \"http://www.w3.org/TR/REC-html40/strict.dtd\">\n"
"<html><head><meta name=\"qrichtext\" content=\"1\" /><style type=\"text/css\">\n"
"p, li { white-space: pre-wrap; }\n"
"</style></head><body style=\" font-family:\'SimSun\'; font-size:9pt; font-weight:400; font-style:normal;\">\n"
"<p style=\"-qt-paragraph-type:empty; margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><br /></p></body></html>"))