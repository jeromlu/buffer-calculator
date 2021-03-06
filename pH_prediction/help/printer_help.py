# -*- coding: utf-8 -*-
"""
Created on Wed Apr  3 08:23:32 2019

@author: jeromlu2
"""

import sys, os
from PyQt5 import QtCore, QtWidgets, QtPrintSupport

class Window(QtWidgets.QWidget):
    def __init__(self):
        super(Window, self).__init__()
        self.setWindowTitle('Document Printer')
        self.editor = QtWidgets.QTextEdit(self)
        self.editor.textChanged.connect(self.handleTextChanged)
        self.buttonOpen = QtWidgets.QPushButton('Open', self)
        self.buttonOpen.clicked.connect(self.handleOpen)
        self.buttonPrint = QtWidgets.QPushButton('Print', self)
        self.buttonPrint.clicked.connect(self.handlePrint)
        self.buttonPreview = QtWidgets.QPushButton('Preview', self)
        self.buttonPreview.clicked.connect(self.handlePreview)
        layout = QtWidgets.QGridLayout(self)
        layout.addWidget(self.editor, 0, 0, 1, 3)
        layout.addWidget(self.buttonOpen, 1, 0)
        layout.addWidget(self.buttonPrint, 1, 1)
        layout.addWidget(self.buttonPreview, 1, 2)
        self.handleTextChanged()

    def handleOpen(self):
        path = QtWidgets.QFileDialog.getOpenFileName(
            self, 'Open file', '',
            'HTML files (*.html);;Text files (*.txt)')[0]
        if path:
            file = QtCore.QFile(path)
            if file.open(QtCore.QIODevice.ReadOnly):
                stream = QtCore.QTextStream(file)
                text = stream.readAll()
                info = QtCore.QFileInfo(path)
                if info.completeSuffix() == 'html':
                    self.editor.setHtml(text)
                else:
                    self.editor.setPlainText(text)
                file.close()

    def handlePrint(self):
        dialog = QtPrintSupport.QPrintDialog()
        print(type(self.editor.document()))
        if dialog.exec_() == QtWidgets.QDialog.Accepted:
            self.editor.document().print_(dialog.printer())

    def handlePreview(self):
        dialog = QtPrintSupport.QPrintPreviewDialog()
        dialog.paintRequested.connect(self.editor.print_)
        dialog.exec_()

    def handleTextChanged(self):
        enable = not self.editor.document().isEmpty()
        self.buttonPrint.setEnabled(enable)
        self.buttonPreview.setEnabled(enable)
        
    def closeEvent(self, evt):
        QtWidgets.QApplication.quit()

if __name__ == '__main__':

    def run_app():
        app = QtWidgets.QApplication(sys.argv)
        window = Window()
        window.resize(640, 480)
        window.show()
        sys.exit(app.exec_())
    run_app()